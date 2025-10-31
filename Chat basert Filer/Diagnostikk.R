# ================== Rough Bergomi ??? Diagnostics ==================
# Forutsetter at simulate_Y_paths_hybrid_k1_v3() og build_variance_paths() er lastet.

rb_diag <- function(H = NULL, eta = NULL,
                    N = NULL, M = NULL,
                    T_max = NULL,
                    xi0_curve = NULL, xi0_level = NULL,
                    seed = NULL, antithetic = TRUE, ncores = 1,
                    skip_small_t = 1L) {
  # ---- Plukk defaults fra milj??et ditt dersom ikke gitt ----
  if (is.null(H)) {
    if (exists("H")) H <- get("H")
    else if (exists("alpha")) H <- get("alpha") + 0.5
    else stop("Gi H eller ha 'H'/'alpha' i environment.")
  }
  if (is.null(eta))  eta  <- if (exists("eta"))  get("eta")  else 1.9
  if (is.null(N))    N    <- if (exists("NCal")) get("NCal") else 500L
  if (is.null(M))    M    <- if (exists("MCal")) get("MCal") else 2000L
  if (is.null(seed)) seed <- if (exists("seed")) get("seed") else 1L
  
  # Finn T_max: bruk xi0_curve, ellers options_df, ellers 1
  if (is.null(T_max)) {
    if (!is.null(xi0_curve) && is.data.frame(xi0_curve)) {
      stopifnot(all(c("T_start","T_end","forward_var") %in% names(xi0_curve)))
      T_max <- max(xi0_curve$T_end)
    } else if (exists("options_df")) {
      T_max <- max(as.numeric(get("options_df")$T), na.rm = TRUE)
    } else {
      T_max <- 1.0
    }
  }
  
  # ---- Simuler Y ----
  sim <- simulate_Y_paths_hybrid_k1_v3(
    M = M, N = N, T = T_max, H = H,
    seed = seed, antithetic = antithetic, ncores = ncores
  )
  tgrid <- sim$t
  dt    <- sim$dt
  
  # ---- Bygg xi0(t_i) p?? samme grid ----
  make_xi0_vec <- function(sim, xi0_curve, xi0_level) {
    if (!is.null(xi0_curve)) {
      if (is.data.frame(xi0_curve) && all(c("T_start","T_end","forward_var") %in% names(xi0_curve))) {
        xi <- numeric(length(sim$t))
        for (i in seq_along(sim$t)) {
          # bruk venstrelukket intervall [T_start, T_end)
          t_start <- (i - 1) * sim$dt
          idx <- which(xi0_curve$T_start <= t_start & t_start < xi0_curve$T_end)
          if (length(idx) == 0) idx <- nrow(xi0_curve)
          xi[i] <- xi0_curve$forward_var[idx]
        }
        return(xi)
      } else if (is.numeric(xi0_curve)) {
        return(if (length(xi0_curve) == 1L) rep(xi0_curve, length(sim$t)) else xi0_curve)
      } else {
        stop("xi0_curve m?? v??re data.frame(T_start,T_end,forward_var) eller numerisk.")
      }
    } else if (!is.null(xi0_level)) {
      return(rep(xi0_level, length(sim$t)))
    } else {
      stop("Gi enten xi0_curve eller xi0_level.")
    }
  }
  xi0_vec <- make_xi0_vec(sim, xi0_curve, xi0_level)
  
  # ---- 1) Var(Y_t) ??? t^(2H) ----
  # Skipp de f??rste 'skip_small_t' punktene for ?? unng?? numerisk st??y n??r t=0
  i0 <- max(1L, as.integer(skip_small_t))
  varY_emp <- apply(sim$Y, 1, var)
  varY_emp <- varY_emp[i0:length(tgrid)]
  varY_the <- (tgrid[i0:length(tgrid)])^(2*H)
  
  varY_rel_err <- (varY_emp - varY_the) / pmax(varY_the, 1e-16)
  
  # ---- 2) E[exp(eta Y_t - 0.5 eta^2 t^{2H})] ??? 1 ----
  centering <- 0.5 * eta^2 * (tgrid^(2*H))
  Eexp_emp  <- rowMeans(exp(eta * sim$Y - matrix(centering, nrow = length(tgrid), ncol = ncol(sim$Y), byrow = FALSE)))
  Eexp_emp  <- Eexp_emp[i0:length(tgrid)]
  Eexp_err  <- Eexp_emp - 1
  
  # ---- 3) v_t = xi0(t) * exp(eta Y_t - 0.5 eta^2 t^{2H}) ??? E[v_t] ??? xi0(t) ----
  v_mat <- build_variance_paths(sim, eta = eta, xi0_curve = xi0_vec)  # N x M
  Ev_emp <- rowMeans(v_mat)[i0:length(tgrid)]
  Ev_the <- xi0_vec[i0:length(tgrid)]
  Ev_rel_err <- (Ev_emp - Ev_the) / pmax(Ev_the, 1e-16)
  
  # ---- 4) E[???_0^t v_s ds] ??? ???_0^t xi0(s) ds (venstre Riemann-sum som i pricer) ----
  if (nrow(v_mat) == 1L) {
    v_prev <- matrix(rep(xi0_vec[1], ncol(v_mat)), nrow = 1L)
  } else {
    v_prev <- rbind(rep(xi0_vec[1], ncol(v_mat)),
                    v_mat[seq_len(nrow(v_mat) - 1L), , drop = FALSE])
  }
  I_paths <- apply(v_prev, 2, function(col) cumsum(col) * dt)  # N x M
  EI_emp  <- rowMeans(I_paths)[i0:length(tgrid)]
  
  # venstresum av xi0 p?? samme skjema som v_prev:
  if (length(xi0_vec) == 1L) {
    xi_left <- xi0_vec
  } else {
    xi_left <- c(xi0_vec[1], xi0_vec[seq_len(length(xi0_vec) - 1L)])
  }
  I_the   <- cumsum(xi_left) * dt
  I_the   <- I_the[i0:length(tgrid)]
  EI_rel_err <- (EI_emp - I_the) / pmax(I_the, 1e-16)
  
  # ---- Sammendrag ----
  summ <- data.frame(
    t = tgrid[i0:length(tgrid)],
    varY_emp = varY_emp,
    varY_the = varY_the,
    varY_rel_err = varY_rel_err,
    Eexp_emp = Eexp_emp,
    Eexp_err = Eexp_err,
    Ev_emp = Ev_emp,
    Ev_the = Ev_the,
    Ev_rel_err = Ev_rel_err,
    EI_emp = EI_emp,
    EI_the = I_the,
    EI_rel_err = EI_rel_err
  )
  
  cat("\n== Diagnostics (med ", M, " paths, N=", N, ", H=", H, ", eta=", eta, ", T_max=", T_max, ") ==\n", sep="")
  cat(sprintf("Var(Y_t) ~ t^{2H}:  median rel.err = %.3g,  95%%-range = [%.3g, %.3g]\n",
              median(summ$varY_rel_err, na.rm=TRUE),
              quantile(summ$varY_rel_err, 0.025, na.rm=TRUE),
              quantile(summ$varY_rel_err, 0.975, na.rm=TRUE)))
  cat(sprintf("E[exp(eta Y - .5 eta^2 t^{2H})] ~ 1:  median abs.err = %.3g,  max abs.err = %.3g\n",
              median(abs(summ$Eexp_err), na.rm=TRUE), max(abs(summ$Eexp_err), na.rm=TRUE)))
  cat(sprintf("E[v_t] ~ xi0(t):  median rel.err = %.3g,  95%%-range = [%.3g, %.3g]\n",
              median(summ$Ev_rel_err, na.rm=TRUE),
              quantile(summ$Ev_rel_err, 0.025, na.rm=TRUE),
              quantile(summ$Ev_rel_err, 0.975, na.rm=TRUE)))
  cat(sprintf("E[??? v ds] ~ ??? xi0 ds: median rel.err = %.3g,  95%%-range = [%.3g, %.3g]\n\n",
              median(summ$EI_rel_err, na.rm=TRUE),
              quantile(summ$EI_rel_err, 0.025, na.rm=TRUE),
              quantile(summ$EI_rel_err, 0.975, na.rm=TRUE)))
  
  invisible(list(summary = summ,
                 sim = sim,
                 xi0_vec = xi0_vec,
                 v_mat = v_mat))
}

# ---------------- Example usage ----------------
rb_out <- rb_diag(H = H, eta = eta, N = 500, M = 10000,
                   xi0_curve = xi0_curve, seed = seed, antithetic = TRUE)
head(rb_out$summary, 10)


xi0_curve_sanity <- function(xi0_curve, T_max = NULL, tol = 1e-12) {
  stopifnot(is.data.frame(xi0_curve),
            all(c("T_start","T_end","forward_var") %in% names(xi0_curve)))
  cur <- xi0_curve
  cur <- cur[order(cur$T_start, cur$T_end), , drop = FALSE]
  cur$len <- pmax(cur$T_end - cur$T_start, 0)
  
  # Varsler
  if (any(cur$len < -tol)) warning("Negativ segmentlengde funnet.")
  if (any(abs(cur$len) < tol)) message("Merk: segment(er) med ~0 lengde (kan droppes).")
  if (abs(cur$T_start[1]) > tol) warning("Kurven starter ikke p?? 0.")
  gaps <- cur$T_start[-1] - cur$T_end[-nrow(cur)]
  max_gap <- max(abs(gaps), na.rm = TRUE)
  if (max_gap > tol) warning(sprintf("Ujevn dekning: maks gap = %.3g", max_gap))
  
  if (is.null(T_max)) T_max <- max(cur$T_end)
  cover_ok <- max(cur$T_end) + tol >= T_max
  
  int_xi <- sum(cur$forward_var * cur$len)
  ave_xi <- int_xi / T_max
  
  out <- list(
    T_max = T_max,
    covers_Tmax = cover_ok,
    segments = cur,
    integral_xi0 = int_xi,
    average_xi0 = ave_xi
  )
  print(data.frame(
    T_max = out$T_max,
    covers_Tmax = out$covers_Tmax,
    n_segments = nrow(cur),
    max_adjacent_gap = max_gap,
    integral_xi0 = out$integral_xi0,
    average_xi0 = out$average_xi0
  ))
  invisible(out)
}

xi_info <- xi0_curve_sanity(xi0_curve)

# Bruker dine sim-/bygg-funksjoner:
# simulate_Y_paths_hybrid_k1_v3(), build_variance_paths()

idx_from_T <- function(T_vec, T_max, N) {
  dt <- T_max / N
  pmax(1L, pmin(N, as.integer(round(T_vec / dt))))
}

rb_breakpoint_report <- function(H, eta, xi0_curve,
                                 N = 1000, M = 2000, seed = 1L,
                                 antithetic = TRUE, skip_small_t = 3L) {
  # Simuler til maks T_end i kurven
  T_max <- max(xi0_curve$T_end)
  sim <- simulate_Y_paths_hybrid_k1_v3(M = M, N = N, T = T_max, H = H,
                                       seed = seed, antithetic = antithetic)
  tgrid <- sim$t; dt <- sim$dt
  
  # Lag xi0 p?? gridet (samme venstresum-logikk som i prisingen)
  xi_vec <- numeric(length(tgrid))
  for (i in seq_along(tgrid)) {
    t_start <- (i - 1) * dt
    k <- which(xi0_curve$T_start <= t_start & t_start < xi0_curve$T_end)
    if (!length(k)) k <- nrow(xi0_curve)
    xi_vec[i] <- xi0_curve$forward_var[k]
  }
  
  # Varians-styring
  var_the <- (tgrid)^(2*H)
  v_mat   <- build_variance_paths(sim, eta = eta, xi0_curve = xi_vec) # N x M
  
  # For integraler bruker vi venstresum (som i turbo-priseren)
  if (nrow(v_mat) == 1L) {
    v_prev <- matrix(rep(xi_vec[1], ncol(v_mat)), nrow = 1L)
  } else {
    v_prev <- rbind(rep(xi_vec[1], ncol(v_mat)),
                    v_mat[seq_len(nrow(v_mat) - 1L), , drop = FALSE])
  }
  I_paths <- apply(v_prev, 2, function(col) cumsum(col) * dt)
  EI_emp  <- rowMeans(I_paths)
  
  if (length(xi_vec) == 1L) {
    xi_left <- xi_vec
  } else {
    xi_left <- c(xi_vec[1], xi_vec[seq_len(length(xi_vec) - 1L)])
  }
  I_the   <- cumsum(xi_left) * dt
  
  # Breakpoints = T_end (dropp duplikater og T=0)
  Tbp <- sort(unique(xi0_curve$T_end))
  Tbp <- Tbp[Tbp > 0]
  j   <- idx_from_T(Tbp, T_max, N)
  j   <- j[j > skip_small_t]  # unng?? f??rste trinn som er mest st??yende
  Tbp <- tgrid[j]
  
  # Diagnoser
  Eexp_emp <- rowMeans(exp(eta * sim$Y - matrix(0.5 * eta^2 * var_the,
                                                nrow = length(tgrid), ncol = ncol(sim$Y), byrow = FALSE)))
  out <- data.frame(
    t          = Tbp,
    Eexp_err   = Eexp_emp[j] - 1,
    Ev_emp     = rowMeans(v_mat)[j],
    Ev_the     = xi_vec[j],
    Ev_rel_err = (rowMeans(v_mat)[j] - xi_vec[j]) / pmax(xi_vec[j], 1e-16),
    EI_emp     = EI_emp[j],
    EI_the     = I_the[j],
    EI_rel_err = (EI_emp[j] - I_the[j]) / pmax(I_the[j], 1e-16)
  )
  rownames(out) <- NULL
  out
}

# Eksempel:
bp <- rb_breakpoint_report(H = H, eta = eta, xi0_curve = xi0_curve,
                            N = 1000, M = 20000, seed = seed)
print(bp, digits = 4)