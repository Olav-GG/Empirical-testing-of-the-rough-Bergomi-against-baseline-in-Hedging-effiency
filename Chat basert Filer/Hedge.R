#### ============================================================
####  Two-Asset Hedging (Paper-style) for rough Bergomi vs Black???Scholes
####  - Tenor matching: hedge factor U_t matches remaining T_t each day
####  - Correct mapping: bump integrated variance on [0, T_t]
####  - U_t is either VIX level ("VIXF") or forward variance Xi_t ("VAR")
#### ============================================================

suppressPackageStartupMessages({ library(dplyr) })

## ---------- Small utils ----------
T_round <- function(x) round(as.numeric(x), 10)
K_round <- function(x) as.numeric(x)
bs_d1 <- function(S,K,T,r,q,sigma) (log(S/K) + (r - q + 0.5*sigma^2)*T)/(sigma*sqrt(T))
bs_vega_mkt <- function(S,K,T,r,q,sigma) S*exp(-q*T)*dnorm(bs_d1(S,K,T,r,q,sigma))*sqrt(T)

## ---------- 1) Build OOS xi0-curve provider (per day) ----------
# If build_xi0_from_ATM is available we use it; else we fall back to flat xi0 at Xi=(VIX/100)^2
xi0_provider_default <- function(date_key, chain_df, r=0, q=0, VIX_level = NULL, extend_to = NULL) {
  has_builder <- exists("build_xi0_from_ATM") && exists("build_calibration_df")
  if (has_builder) {
    # sparse grid around each expiry, estimate ATM IV -> iso TV -> xi0
    exps <- unique(chain_df$Expiry)
    grid <- lapply(exps, function(e) list(expiration = e, ratios = c(0.8,0.9,1.0,1.1,1.2)))
    S0T <- {
      T_min <- suppressWarnings(min(as.numeric(chain_df$Texp), na.rm=TRUE))
      near  <- subset(chain_df, is.finite(Fwd) & is.finite(Texp) & abs(Texp - T_min) <= 1e-6)
      if (nrow(near)==0) near <- subset(chain_df, is.finite(Fwd) & is.finite(Texp))
      median(near$Fwd * exp(-(r - q) * near$Texp), na.rm = TRUE)
    }
    calib <- build_calibration_df(
      spxIvolList = list(tmp=chain_df), date_key="tmp",
      grid_list = grid, S0=S0T, r=r, q=q, use_chain_fwd=TRUE
    )
    xi0 <- build_xi0_from_ATM(calib, extend_to = extend_to %||% max(as.numeric(chain_df$Texp), na.rm=TRUE))
    return(xi0)
  } else {
    # flat curve at Xi=(VIX/100)^2
    if (!is.finite(VIX_level)) stop("Need VIX_level for flat xi0 fallback.")
    Xi <- (VIX_level/100)^2
    Tmax <- extend_to %||% max(as.numeric(chain_df$Texp), na.rm=TRUE)
    data.frame(T_start=0, T_end=Tmax, forward_var=Xi)
  }
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

## ---------- 2) Curve scaling on [0, Delta] ----------
# Multiply forward_var on [0, Delta] by a factor mu
scale_xi0_on_interval <- function(xi0_df, Delta, mu) {
  out <- xi0_df
  if (Delta <= 0) return(out)
  sel <- out$T_start < Delta  # left-closed scheme consistent with turbo-pricer
  out$forward_var[sel] <- out$forward_var[sel] * mu
  out
}

# VIX bump: VIX = 100*sqrt( avg_xi(0,Delta) ), change VIX by +/- eps_vix points
make_vix_bumped_curves <- function(xi0_df, VIX_level, Delta, eps_vix = 0.5) {
  sigma <- pmax(VIX_level/100, 1e-8)
  mu_up <- ((sigma + eps_vix/100)^2)/(sigma^2)   # ??? 1 + 2*(eps_vix/100)/sigma
  sigma_dn <- pmax(sigma - eps_vix/100, 1e-8)
  mu_dn <- (sigma_dn^2)/(sigma^2)
  list(
    up = scale_xi0_on_interval(xi0_df, Delta, mu_up),
    dn = scale_xi0_on_interval(xi0_df, Delta, mu_dn),
    eps = eps_vix
  )
}

# VAR bump: U = Xi = sigma^2, change average Xi on [0,Delta] by +/- dXi
make_var_bumped_curves <- function(xi0_df, VIX_level, Delta, dXi = 1e-4) {
  sigma2 <- pmax((VIX_level/100)^2, 1e-12)
  mu_up  <- (sigma2 + dXi)/sigma2        # ??? 1 + dXi/sigma^2
  mu_dn  <- pmax((sigma2 - dXi), 1e-12)/sigma2
  list(
    up = scale_xi0_on_interval(xi0_df, Delta, mu_up),
    dn = scale_xi0_on_interval(xi0_df, Delta, mu_dn),
    dXi = dXi
  )
}

## ---------- 3) rB sensitivities: Delta and dP/dU (paper-style) ----------
# Prices the target option under (i) base curve, (ii) bumped curves on [0, T_t]
rb_delta_dPdU <- function(
    S0, K, T, type = c("call","put"), r=0, q=0,
    HCal, EtaCal, RhoCal,
    xi0_curve_base,
    hedge_asset = c("VIXF","VAR"),
    vix_level,
    dS_rel = 1e-4,
    eps_vix = 0.5,        # VIX point bump (for U = VIX)
    dXi = 1e-4,           # Xi bump (for U = Xi), in variance units
    N, M, seed = 1, antithetic = TRUE, ncores = 1
) {
  type <- match.arg(type); hedge_asset <- match.arg(hedge_asset)
  price_with_curve <- function(xi_df, S_override = NULL) {
    S_use <- S_override %||% S0
    odf <- data.frame(type=type, K=K, T=T)
    out <- price_many_rB_turbo(
      options_df = odf,
      S0 = S_use, r=r, q=q,
      alpha = HCal - 0.5, eta = EtaCal, rho = RhoCal,
      xi0_curve = xi_df,
      N = N, M = M, antithetic = antithetic, seed = seed,
      ncores = ncores, clip_iv = FALSE
    )
    as.numeric(out$price[1])
  }
  
  # Delta via central diff in S
  dS <- max(1e-8, abs(dS_rel)*S0)
  P_upS <- price_with_curve(xi0_curve_base, S_override = S0 + dS)
  P_dnS <- price_with_curve(xi0_curve_base, S_override = S0 - dS)
  delta <- (P_upS - P_dnS)/(2*dS)
  
  # dP/dU via central diff by bumping integrated variance on [0,T]
  if (hedge_asset == "VIXF") {
    bumps <- make_vix_bumped_curves(xi0_curve_base, vix_level, Delta = T, eps_vix = eps_vix)
    P_upU <- price_with_curve(bumps$up)
    P_dnU <- price_with_curve(bumps$dn)
    dPdU  <- (P_upU - P_dnU)/(2*bumps$eps)     # units: $ per VIX point
  } else {
    bumps <- make_var_bumped_curves(xi0_curve_base, vix_level, Delta = T, dXi = dXi)
    P_upU <- price_with_curve(bumps$up)
    P_dnU <- price_with_curve(bumps$dn)
    dPdU  <- (P_upU - P_dnU)/(2*bumps$dXi)     # units: $ per variance unit
  }
  
  list(delta = unname(delta), dPdU = unname(dPdU))
}

## ---------- 4) Main hedging loop (paper-style) ----------
fht_two_asset_hedge_rb_paper <- function(
    BT_NDAYS = 30,
    hedge_asset = c("VIXF","VAR"),
    r = 0, q = 0,
    # transaction costs (kept simple; set to 0 for methodology tests)
    tc_bp_S = 0.0, tc_tick_U = 0.0,
    # calibrated rB params
    HCal, EtaCal, RhoCal,
    # rB grid
    N_rB, M_rB, seed = 42, antithetic = TRUE,
    # xi0 provider: function(date_key, chain_df, r, q, VIX_level, extend_to) -> xi0_df
    xi0_provider = xi0_provider_default,
    verbose_every = 5
) {
  hedge_asset <- match.arg(hedge_asset)
  
  # Requirements
  reqs <- c("spxIvolList","VIX30","robust_mid_call","pick_expiry_near_T","pick_row_nearest_strike",
            "price_many_rB_turbo","bs_implied_vol","bs_greeks")
  stopifnot(all(vapply(reqs, exists, logical(1))))
  
  # VIX 30D path aligned with spxIvolList keys (simple tail match)
  to_date_chr <- function(x) {
    if (inherits(x,"Date")) return(format(x,"%Y%m%d"))
    if (inherits(x,"POSIXt")) return(format(as.Date(x),"%Y%m%d"))
    if (is.character(x)) return(gsub("-","", substr(x,1,10)))
    format(as.Date(x), "%Y%m%d")
  }
  VIX30_df <- data.frame(date_chr = to_date_chr(VIX30$Date),
                         cm30 = as.numeric(VIX30$PX_LAST))
  VIX30_df <- VIX30_df[is.finite(VIX30_df$cm30), ]
  VIX30_df <- VIX30_df[order(VIX30_df$date_chr), ]
  
  keys_spx <- sort(names(spxIvolList))
  stopifnot(length(keys_spx) >= BT_NDAYS, nrow(VIX30_df) >= BT_NDAYS)
  keys_path <- tail(keys_spx, BT_NDAYS)
  vix_path  <- tail(VIX30_df$cm30, BT_NDAYS)
  
  ## ---- Select target option on day 1 (ATM-fwd ~30D) ----
  date0 <- keys_path[1]
  df0   <- spxIvolList[[date0]]
  cm0   <- robust_mid_call(df0)
  keep0 <- is.finite(df0$Strike) & df0$Strike>0 & is.finite(df0$Texp) & df0$Texp>0 &
    is.finite(df0$Fwd) & df0$Fwd>0 & is.finite(cm0) & cm0>0
  df0 <- transform(df0[keep0, , drop=FALSE], CallMid = cm0[keep0])
  S0_hat <- median(df0$Fwd * exp(-(r - q)*df0$Texp))
  expiry_target <- pick_expiry_near_T(df0, T_target_days = 30)
  sub_exp <- df0[df0$Expiry == expiry_target, , drop=FALSE]
  i_atm   <- which.min(abs(sub_exp$Strike - sub_exp$Fwd))
  K_tgt   <- as.numeric(sub_exp$Strike[i_atm])
  T0_tgt  <- as.numeric(sub_exp$Texp[i_atm])
  type_tgt <- "call"
  cat(sprintf("[PAPER] Target: %s  Exp=%s  K=%.2f  T0=%.6f\n", date0,
              as.character(expiry_target), K_tgt, T0_tgt))
  
  ## ---- Init portfolios ----
  state <- list(BS=list(x=NA_real_, y=NA_real_, cash=NA_real_),
                rB=list(x=NA_real_, y=NA_real_, cash=NA_real_))
  first <- TRUE; out_rows <- list()
  
  ## ---- Daily hedging ----
  for (i in seq_len(BT_NDAYS)) {
    k <- keys_path[i]
    VIX_t <- vix_path[i]; if (!is.finite(VIX_t)) next
    df <- spxIvolList[[k]]
    cm <- robust_mid_call(df)
    keep <- is.finite(df$Strike) & df$Strike>0 & is.finite(df$Texp) & df$Texp>0 &
      is.finite(df$Fwd) & df$Fwd>0 & is.finite(cm) & cm>0
    df <- transform(df[keep, , drop=FALSE], CallMid = cm[keep])
    row_t <- pick_row_nearest_strike(df, expiry = expiry_target, K_target = K_tgt)
    if (is.null(row_t)) next
    
    T_t   <- as.numeric(row_t$Texp)
    F_t   <- as.numeric(row_t$Fwd)
    S_t   <- F_t * exp(-(r - q)*T_t)
    P_tgt <- as.numeric(row_t$CallMid)
    Sigma_t <- VIX_t/100
    Xi_t    <- Sigma_t^2
    U_t     <- if (hedge_asset=="VIXF") VIX_t else Xi_t
    
    ## ---- BS greeks (simple mapping) ----
    iv_t <- bs_implied_vol(P_tgt, S=S_t, K=K_tgt, T=T_t, r=r, q=q, type=type_tgt)
    if (!is.finite(iv_t) || iv_t<=0) next
    g_bs <- bs_greeks(S=S_t, K=K_tgt, T=T_t, r=r, q=q, sigma=iv_t, type=type_tgt)
    x_BS_new <- g_bs$delta
    y_BS_new <- if (hedge_asset=="VIXF") g_bs$vega/100 else g_bs$vega/(2*Sigma_t)  # $/VIXpt or $/variance
    
    ## ---- rB greeks (paper-style dP/dU with tenor-matched bump on [0, T_t]) ----
    xi0_base <- tryCatch(
      xi0_provider(date_key = k, chain_df = df, r=r, q=q, VIX_level = VIX_t, extend_to = T_t),
      error = function(e) data.frame(T_start=0, T_end=T_t, forward_var=Xi_t) # flat fallback
    )
    rb_sens <- rb_delta_dPdU(
      S0=S_t, K=K_tgt, T=T_t, type=type_tgt, r=r, q=q,
      HCal=HCal, EtaCal=EtaCal, RhoCal=RhoCal,
      xi0_curve_base = xi0_base,
      hedge_asset = hedge_asset, vix_level = VIX_t,
      dS_rel = 1e-4, eps_vix = 0.5, dXi = 1e-4,
      N=N_rB, M=M_rB, seed=seed, antithetic=antithetic, ncores = max(1L, parallel::detectCores()-1L)
    )
    x_rB_new <- rb_sens$delta
    y_rB_new <- rb_sens$dPdU   # already $/VIXpt or $/variance depending on hedge_asset
    
    ## ---- Self-financing update + TC ----
    if (first) {
      state$BS$cash <- P_tgt - x_BS_new*S_t - y_BS_new*U_t
      state$BS$x <- x_BS_new; state$BS$y <- y_BS_new
      state$rB$cash <- P_tgt - x_rB_new*S_t - y_rB_new*U_t
      state$rB$x <- x_rB_new; state$rB$y <- y_rB_new
      first <- FALSE
    } else {
      dx <- x_BS_new - state$BS$x; dy <- y_BS_new - state$BS$y
      state$BS$cash <- state$BS$cash - dx*S_t - dy*U_t - (abs(dx)*abs(S_t)*tc_bp_S + abs(dy)*tc_tick_U)
      state$BS$x <- x_BS_new; state$BS$y <- y_BS_new
      
      dx <- x_rB_new - state$rB$x; dy <- y_rB_new - state$rB$y
      state$rB$cash <- state$rB$cash - dx*S_t - dy*U_t - (abs(dx)*abs(S_t)*tc_bp_S + abs(dy)*tc_tick_U)
      state$rB$x <- x_rB_new; state$rB$y <- y_rB_new
    }
    
    HE_BS <- state$BS$cash + state$BS$x*S_t + state$BS$y*U_t - P_tgt
    HE_rB <- state$rB$cash + state$rB$x*S_t + state$rB$y*U_t - P_tgt
    
    out_rows[[length(out_rows)+1L]] <- data.frame(
      day=i, date=k, S_t=S_t, U_t=U_t, T_t=T_t, P_tgt=P_tgt,
      VIX_t=VIX_t, Xi_t=Xi_t, iv_t=iv_t,
      x_BS=state$BS$x, y_BS=state$BS$y, HE_BS=HE_BS,
      x_rB=state$rB$x, y_rB=state$rB$y, HE_rB=HE_rB,
      stringsAsFactors = FALSE
    )
    
    if (verbose_every>0 && (i %% verbose_every == 1))
      cat(sprintf("[Day %02d | %s]  ?? (BS=%.3e,rB=%.3e)   dP/dU (BS=%.3e,rB=%.3e)\n",
                  i, hedge_asset, state$BS$x, state$rB$x, y_BS_new, y_rB_new))
  }
  
  path <- do.call(rbind, out_rows)
  stat <- function(x) c(n=length(x), mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE),
                        rmse=sqrt(mean(x^2,na.rm=TRUE)),
                        q05=unname(quantile(x,0.05,na.rm=TRUE)),
                        q95=unname(quantile(x,0.95,na.rm=TRUE)))
  cat("\n=== Hedging Error Summary (paper-style) ===\n")
  print(rbind(BS = stat(path$HE_BS), rB = stat(path$HE_rB)), digits=4)
  
  invisible(list(path=path,
                 summary=list(BS=stat(path$HE_BS), rB=stat(path$HE_rB))))
}

## ---------- Example run (uses your existing globals if present) ----------
set.seed(42)
res_vixf <- fht_two_asset_hedge_rb_paper(
  BT_NDAYS = 30,
  hedge_asset = "VIXF",
  r = get0("r", ifnotfound=0), q = get0("q", ifnotfound=0),
  tc_bp_S = 0.0, tc_tick_U = 0.0,
  HCal = get0("HCal", ifnotfound=0.12),
  EtaCal = get0("EtaCal", ifnotfound=2.0),
  RhoCal = get0("RhoCal", ifnotfound=-0.8),
  N_rB = get0("NHedge", ifnotfound=300),
  M_rB = get0("MHedge", ifnotfound=800),
  xi0_provider = xi0_provider_default
)

res_var <- fht_two_asset_hedge_rb_paper(
  BT_NDAYS = 30,
  hedge_asset = "VAR",
  r = get0("r", ifnotfound=0), q = get0("q", ifnotfound=0),
  tc_bp_S = 0.0, tc_tick_U = 0.0,
  HCal = get0("HCal", ifnotfound=0.12),
  EtaCal = get0("EtaCal", ifnotfound=2.0),
  RhoCal = get0("RhoCal", ifnotfound=-0.8),
  N_rB = get0("NHedge", ifnotfound=300),
  M_rB = get0("MHedge", ifnotfound=800),
  xi0_provider = xi0_provider_default
)
