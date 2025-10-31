# ========================= rBergomi_turbo.R =========================
# Turbocharged Monte Carlo pricing for the rough Bergomi model using conditional Black???Scholes.
# McCrickerd & Pakkanen (2018): integrates out the independent price noise analytically.
# This uses the hybrid scheme to simulate the Y process and then prices via conditional BS.

# ---- Simulation of fractional noise (Y paths) ----
simulate_Y_paths_hybrid_k1_v3 <- function(
    M, N, T, H, seed = NULL, antithetic = FALSE, ncores = 1
) {
  if (!is.null(seed)) set.seed(seed)
  dt <- T / N
  a  <- H - 0.5
  t  <- (1:N) * dt
  
  # Compute tail weights for the Riemann???Liouville fractional integral (for convolution)
  if (abs(a) < 1e-12) {
    w_tail <- rep(sqrt(2*H) * dt^a, N - 1)   # = 1 when H=0.5
  } else {
    k <- 2:N
    bk_star <- ((k^(a+1) - (k-1)^(a+1)) / (a+1))^(1/a)
    w_tail  <- sqrt(2*H) * ((bk_star * dt)^a)
  }
  w_all <- c(0, w_tail)
  
  # Coefficients for the Y process
  c1 <- sqrt(2*H) * dt^a / (a + 1)
  c2 <- dt^(a + 0.5) * sqrt(2*H * (1/(2*a + 1) - 1/(a + 1)^2))
  
  # --- RNG: generate all noise first (deterministic w.rt. seed) ---
  if (antithetic) {
    M2 <- ceiling(M/2)
    dW_half <- matrix(rnorm(N * M2, sd = sqrt(dt)), nrow = N)
    Z_half  <- matrix(rnorm(N * M2), nrow = N)
    dW_mat  <- cbind(dW_half, -dW_half)[, 1:M, drop = FALSE]
    Z_mat   <- cbind(Z_half,  -Z_half)[, 1:M, drop = FALSE]
  } else {
    dW_mat <- matrix(rnorm(N * M, sd = sqrt(dt)), nrow = N)
    Z_mat  <- matrix(rnorm(N * M), nrow = N)
  }
  
  # --- Parallel convolution per path ---
  conv_one <- function(m) {
    dW <- dW_mat[, m]
    Z  <- Z_mat[, m]
    tail_conv <- as.numeric(stats::convolve(dW, rev(w_all), type = "open"))[1:N]
    return(c1 * dW + c2 * Z + tail_conv)
  }
  
  if (ncores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply on Unix/macOS
    Y_list <- parallel::mclapply(seq_len(ncol(dW_mat)), conv_one, mc.cores = ncores)
    Y <- do.call(cbind, Y_list)
  } else {
    # On Windows mclapply is unavailable; fall back to the sequential lapply path
    Y <- do.call(cbind, lapply(seq_len(ncol(dW_mat)), conv_one))
  }
  
  return(list(t = t, Y = Y, dt = dt, H = H, dW = dW_mat))
}

# ---------- 1) Build variance paths given Y-paths ----------
build_variance_paths <- function(sim, eta, xi0_vec) {
  t <- sim$t; H <- sim$H
  N <- length(t); M <- ncol(sim$Y)
  stopifnot(length(xi0_vec) == N)
  drift  <- -0.5 * eta^2 * (t^(2*H))
  xi_mat <- matrix(xi0_vec, nrow = N, ncol = M)
  xi_mat * exp(eta * sim$Y + matrix(drift, nrow = N, ncol = M, byrow = FALSE))
}


# ---------- 2) Map maturity T to time-step index ----------
idx_from_T <- function(T_vec, T_max, N) {
  dt <- T_max / N
  # Map each maturity T to the closest index on [1, N]
  return(pmax(1L, pmin(N, as.integer(round(T_vec / dt)))))
}

make_xi0_vec <- function(t_grid, dt, xi0_curve, xi0_level) {
  xi <- NULL
  if (!is.null(xi0_curve)) {
    if (is.data.frame(xi0_curve) && all(c("T_start", "T_end", "forward_var") %in% names(xi0_curve))) {
      xi <- numeric(length(t_grid))
      for (i in seq_along(t_grid)) {
        t_start <- (i - 1) * dt
        idx <- which(xi0_curve$T_start <= t_start & t_start < xi0_curve$T_end)
        if (length(idx) == 0) idx <- nrow(xi0_curve)
        xi[i] <- xi0_curve$forward_var[idx]
      }
    } else if (is.numeric(xi0_curve)) {
      xi <- if (length(xi0_curve) == 1L) rep(xi0_curve, length(t_grid)) else xi0_curve
    } else if (is.function(xi0_curve)) {
      xi <- xi0_curve(t_grid)
    } else {
      stop("Unsupported xi0_curve type. Use data.frame(T_start,T_end,forward_var), numeric, or function.")
    }
  } else if (!is.null(xi0_level)) {
    xi <- rep(xi0_level, length(t_grid))
  }
  if (is.null(xi)) stop("Provide either xi0_curve or xi0_level.")
  if (length(xi) != length(t_grid)) {
    stop("xi0 specification does not match simulation grid length.")
  }
  xi
}


# ---------- 3) Option pricing via turbo rough Bergomi ----------
price_many_rB_turbo <- function(options_df, S0, r, q,
                                alpha, eta, rho,
                                xi0_level = NULL, xi0_curve = NULL,
                                N, M, antithetic = TRUE, seed = NULL,
                                ncores = 1, clip_iv = TRUE) {
  stopifnot(all(c("type","K","T") %in% names(options_df)))
  # Clean and validate inputs:
  df <- options_df
  df$type <- tolower(trimws(df$type))
  df <- df[df$type %in% c("call","put"), , drop = FALSE]
  df <- df[is.finite(df$K) & is.finite(df$T) & df$T > 0, , drop = FALSE]
  df$K <- as.numeric(df$K)
  df$T <- as.numeric(df$T)
  # Stabilize T values to avoid floating-point grouping issues
  df$T <- round(df$T, 10)
  if (nrow(df) == 0L) stop("No valid options after cleaning (check type/K/T).")
  
  T_vec <- sort(unique(df$T))
  T_max <- max(T_vec)
  H     <- alpha + 0.5
  
  # Simulate fractional noise Y up to the longest maturity
  sim <- simulate_Y_paths_hybrid_k1_v3(M = M, N = N, T = T_max, H = H,
                                       seed = seed, antithetic = antithetic, ncores = ncores)
  # Build per-step variance process from xi0_curve or xi0_level
  xi0_vec <- make_xi0_vec(sim$t, sim$dt, xi0_curve = xi0_curve, xi0_level = xi0_level)
  
  # Simulate variance paths given xi0_vec
  v_mat <- build_variance_paths(sim, eta, xi0_vec)   # dimensions: N x M
  # v_prev[i,] = variance at *start* of step i (for i from 1 to N)
  if (nrow(v_mat) == 1L) {
    v_prev <- matrix(rep(xi0_vec[1], ncol(v_mat)), nrow = 1L)
  } else {
    v_prev <- rbind(rep(xi0_vec[1], ncol(v_mat)),
                    v_mat[seq_len(nrow(v_mat) - 1L), , drop = FALSE])
  }
  v_prev <- pmax(v_prev, 0)  # ensure non-negative variances
  
  # Precompute maturity indices and grouping of options by (type, T_id)
  T_idx_all <- idx_from_T(T_vec, T_max, N)
  T_id <- match(df$T, T_vec)                     # gives an integer ID for each unique T
  key  <- paste(df$type, T_id)
  groups <- split(seq_len(nrow(df)), key)
  # Map each maturity ID to its group keys (call/put groups)
  T_to_groups <- split(names(groups),
                       as.integer(sapply(names(groups), function(k) strsplit(k, " ")[[1]][2])))
  
  # Map grid step index -> all maturity IDs (T_ids) that snap to that step
  idx_by_Tid <- T_idx_all                                  # length = length(T_vec)
  ids_by_idx <- split(seq_along(T_vec), idx_by_Tid)        # e.g. "17" -> c(3,4)
  
  
  # --- Traverse the time grid once; price options at maturity steps ---
  I_acc <- rep(0.0, M)   # ??? V_u du
  J_acc <- rep(0.0, M)   # ??? sqrt(V_u) dW_u
  records <- list()
  
  for (i in seq_len(nrow(v_prev))) {
    # Accumulate integrals pathwise
    vi <- v_prev[i, ]
    I_acc <- I_acc + vi * sim$dt
    J_acc <- J_acc + sqrt(vi) * sim$dW[i, ]
    
    # All maturities that snap to this grid step
    t_ids_here <- ids_by_idx[[as.character(i)]]
    if (!length(t_ids_here)) next
    
    # Conditional BS inputs common to all maturities hitting this step
    A_base <- S0 * exp(-0.5 * I_acc + rho * J_acc)  # no dividend yet
    s2     <- (1 - rho^2) * I_acc
    
    for (tid in t_ids_here) {
      Tm     <- T_vec[tid]
      gnames <- T_to_groups[[as.character(tid)]]
      if (is.null(gnames)) next
      
      # Dividend discount uses the *actual* maturity Tm
      A_T <- A_base * exp(-q * Tm)
      
      for (gn in gnames) {
        ix   <- groups[[gn]]
        typ  <- unique(df$type[ix])
        Kvec <- df$K[ix]
        
        price_vec <- numeric(length(Kvec))
        se_vec    <- numeric(length(Kvec))
        
        for (j in seq_along(Kvec)) {
          Kd <- Kvec[j] * exp(-r * Tm)
          payoff <- .cond_bs_pathwise(A_T, s2, Kd, type = typ)
          price_vec[j] <- mean(payoff)
          se_vec[j]    <- stats::sd(payoff) / sqrt(length(payoff))
        }
        ci95  <- 1.96 * se_vec
        ci_lo <- price_vec - ci95
        ci_hi <- price_vec + ci95
        
        records[[length(records) + 1L]] <- data.frame(
          type = typ, K = Kvec, T = Tm,
          price = price_vec, se = se_vec, ci_lo = ci_lo, ci_hi = ci_hi,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  
  if (length(records) == 0L) stop("No prices were produced (check T grid mapping).")
  res <- do.call(rbind, records)
  res <- res[order(res$T, match(res$type, c("call","put")), res$K), , drop = FALSE]
  rownames(res) <- NULL
  
  if (clip_iv) {
    # Enforce no-arbitrage bounds (price between intrinsic lower bound and forward price upper bound)
    res <- within(res, {
      disc_r <- exp(-r * T); disc_q <- exp(-q * T)
      lower  <- ifelse(type == "call", pmax(0, S0*disc_q - K*disc_r),
                       pmax(0, K*disc_r - S0*disc_q))
      upper  <- ifelse(type == "call", S0*disc_q, K*disc_r)
      price  <- pmin(pmax(price, lower + 1e-10), upper - 1e-10)
    })
  }
  return(res)
}

.cond_bs_pathwise <- function(A, s2, Kd, type = c("call", "put"), tiny = 1e-14) {
  type <- match.arg(type)
  s2 <- pmax(s2, 0)
  sd <- sqrt(s2)
  out <- numeric(length(A))
  use_bs <- (sd > 0)
  
  z <- numeric(length(A))
  z[use_bs] <- log(pmax(A[use_bs], tiny) / pmax(Kd, tiny))
  
  d1 <- d2 <- numeric(length(A))
  d1[use_bs] <- (z[use_bs] + 0.5 * s2[use_bs]) / sd[use_bs]
  d2[use_bs] <- d1[use_bs] - sd[use_bs]
  
  if (type == "call") {
    out[use_bs]  <- A[use_bs] * pnorm(d1[use_bs]) - Kd * pnorm(d2[use_bs])
    out[!use_bs] <- pmax(A[!use_bs] - Kd, 0)
  } else {
    out[use_bs]  <- Kd * pnorm(-d2[use_bs]) - A[use_bs] * pnorm(-d1[use_bs])
    out[!use_bs] <- pmax(Kd - A[!use_bs], 0)
  }
  out
}
