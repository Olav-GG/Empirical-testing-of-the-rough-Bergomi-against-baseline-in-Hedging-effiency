M=10000
N=1000
alpha=-0.43
H=alpha+1/2
seed=1
#xi0_level=(tail(VIX$PX_LAST, 1)/100)^2; xi0_level
eta=1.9
rho=-0.9

simulate_Y_paths_hybrid_k1_v3 <- function(
    M, N, T, H, b = 0.5, seed = NULL, antithetic = FALSE,
    ncores = 1
) {
  if (!is.null(seed)) set.seed(seed)
  dt <- T / N
  a  <- H - 0.5
  t  <- (1:N) * dt
  
  # tail-vekter
  if (abs(a) < 1e-12) {
    w_tail <- rep(sqrt(2*H) * dt^a, N - 1)   # = 1 n??r H=0.5
  } else {
    k <- 2:N
    bk_star <- ((k^(a+1) - (k-1)^(a+1))/(a+1))^(1/a)
    w_tail  <- sqrt(2*H) * ((bk_star * dt)^a)
  }
  w_all <- c(0, w_tail)
  
  # koeffisienter
  c1 <- sqrt(2*H) * dt^a / (a + 1)
  c2 <- dt^(a + 0.5) * sqrt( 2*H * ( 1/(2*a + 1) - 1/(a + 1)^2 ) )
  
  # --- RNG: lag all st??y f??rst (deterministisk mht. seed) ---
  if (antithetic) {
    M2 <- ceiling(M/2)
    dW_half <- matrix(rnorm(N * M2, sd = sqrt(dt)), nrow = N)
    Z_half  <- matrix(rnorm(N * M2),                  nrow = N)
    dW_mat  <- cbind(dW_half, -dW_half)[, 1:M, drop = FALSE]
    Z_mat   <- cbind(Z_half,  -Z_half )[ , 1:M, drop = FALSE]
  } else {
    dW_mat <- matrix(rnorm(N * M, sd = sqrt(dt)), nrow = N)
    Z_mat  <- matrix(rnorm(N * M),                nrow = N)
  }
  
  # --- Parallell konvolusjon per-bane ---
  conv_one <- function(m) {
    dW <- dW_mat[, m]
    Z  <- Z_mat[, m]
    tail_conv <- as.numeric(stats::convolve(dW, rev(w_all), type = "open"))[1:N]
    c1 * dW + c2 * Z + tail_conv
  }
  
  if (ncores > 1 && .Platform$OS.type != "windows") {
    # Unix/macOS: mclapply
    Y_list <- parallel::mclapply(seq_len(ncol(dW_mat)), conv_one, mc.cores = ncores)
    Y <- do.call(cbind, Y_list)
  } else if (ncores > 1 && .Platform$OS.type == "windows") {
    
    Y <- do.call(cbind, lapply(seq_len(ncol(dW_mat)), conv_one))
  } else {
    # Sekvensielt
    Y <- do.call(cbind, lapply(seq_len(ncol(dW_mat)), conv_one))
  }
  
  list(t = t, Y = Y, dt = dt, H = H, dW = dW_mat)
}

ncores <- max(1L, parallel::detectCores() - 1L)

sim_hyb <- simulate_Y_paths_hybrid_k1_v3(
  M = M, N = N, T = T, H = H,
  seed = seed, antithetic = TRUE, ncores = ncores
)

#empirisk check
check_moments <- function(sim) {
  stopifnot(is.list(sim), !is.null(sim$t), !is.null(sim$Y), !is.null(sim$H))
  t <- sim$t; Y <- sim$Y; H <- sim$H
  mc_mean <- rowMeans(Y)
  mc_var  <- apply(Y, 1, var)     # sample variance over baner
  th_var  <- t^(2 * H)            # teoretisk varians
  data.frame(t = t,
             mc_mean = mc_mean,
             mc_var  = mc_var,
             th_var  = th_var,
             var_ratio = mc_var / th_var)
}

#volatilitets delen
t <- sim_hyb$t; H <- sim_hyb$H
drift_row <- -0.5 * eta^2 * t^(2 * H)
v <- xi0_level * exp(eta * sim_hyb$Y +
                       matrix(drift_row, nrow = length(t), ncol = ncol(sim_hyb$Y), byrow = FALSE))

#Empirisk sjekk av ulike momenter:
mom_hyb <- check_moments(sim_hyb)
idx <- unique(pmax(1, round(c(0.01,0.1,0.2, 0.4,0.5,0.6,0.7,0.8,
                              0.9, 1.00) * N)))
mv <- rowMeans(v)
idxv <- unique(pmax(1, round(c(0.01,0.1,0.2, 0.4,0.5,0.6,0.7,0.8,
                               0.9, 1.00) * length(t))))



mom_hyb[idxv, c("t", "mc_mean", "mc_var", "th_var", "var_ratio")]
cbind(t = t[idx], mc_mean_v = mv[idx], target = xi0_level)

mean(mv)
sqrt(mean(mv))

# Hent grid-indekser for hver T
idx_from_T <- function(T_vec, T_max, N) {
  dt <- T_max / N
  pmax(1L, pmin(N, round(T_vec / dt)))
}

# Bygg v(t) = xi0 * exp(eta*Y + drift)  (flat forward-niv??)
#build_variance_paths <- function(sim, eta, xi0_level){
#  t <- sim$t; H <- sim$H; N <- length(t); M <- ncol(sim$Y)
#  drift <- -0.5 * eta^2 * (t^(2*H))
#  xi_mat <- matrix(xi0_level, nrow=N, ncol=M)
#  xi_mat * exp( eta * sim$Y + matrix(drift, nrow=N, ncol=M, byrow=FALSE) )
#}


build_variance_paths_curve <- function(sim, eta, xi0_curve) {
  stopifnot(length(xi0_curve) == length(sim$t))
  t <- sim$t; H <- sim$H; N <- length(t); M <- ncol(sim$Y)
  drift <- -0.5 * eta^2 * (t^(2*H))
  xi_mat <- matrix(xi0_curve, nrow = N, ncol = M)  # repliker kurven nedover kolonner
  xi_mat * exp(eta * sim$Y + matrix(drift, nrow = N, ncol = M, byrow = FALSE))
}

# Simuler S og tapp ut S_T for mange modenheter i samme l??p
simulate_ST_paths <- function(S0, r, q, rho, sim, v_mat, xi0_0, T_idx, antithetic=TRUE){
  N <- nrow(v_mat); M <- ncol(v_mat); dt <- sim$dt; sqrt_dt <- sqrt(dt)
  v_prev <- rbind(rep(xi0_0, M), v_mat[1:(N-1), , drop = FALSE])
  logS <- rep(log(S0), M)
  
  if (antithetic) {
    M2 <- ceiling(M/2); Zp <- matrix(rnorm(N * M2), nrow = N)
    Z  <- cbind(Zp, -Zp)[, 1:M, drop = FALSE]
  }
  
  ST_list <- vector("list", length(T_idx)); k_take <- 1L
  for (i in 1:N) {
    dW_perp <- if (antithetic) sqrt_dt * Z[i, ] else sqrt_dt * rnorm(M)
    dB_i    <- rho * sim$dW[i, ] + sqrt(1 - rho^2) * dW_perp
    vi <- pmax(v_prev[i, ], 0)
    logS <- logS + (r - q - 0.5 * vi) * dt + sqrt(vi) * dB_i
    if (k_take <= length(T_idx) && i == T_idx[k_take]) {
      ST_list[[k_take]] <- exp(logS)
      k_take <- k_take + 1L
    }
  }
  ST_list
}

# options_df m?? ha kolonnene: type ("call"/"put"), K, T (??r)
price_many_rB <- function(options_df, S0, r, q,
                          alpha, eta, rho, xi0_level,
                          N=N, M=M, antithetic=TRUE, seed=seed,
                          clip_iv = TRUE) {
  stopifnot(all(c("type","K","T") %in% names(options_df)))
  df <- options_df
  df$type <- tolower(trimws(df$type))
  if (!all(df$type %in% c("call","put"))) stop("type m?? v??re 'call' eller 'put'.")
  
  # 1) Simuler til maks forfall (??N gang)
  T_vec <- sort(unique(df$T))
  T_max <- max(T_vec)
  H <- alpha + 0.5
  sim <- simulate_Y_paths_hybrid_k1_v3(M=M, N=N, T=T_max, H=H, seed=seed, antithetic=antithetic)
  v_mat <- build_variance_paths_curve(sim, eta, xi0_curve)  # <- bruk kurven
  xi0_0 <- xi0_curve[1]       
  
  # 2) Hent S_T for alle T i ett pass
  T_idx <- idx_from_T(T_vec, T_max, N)
  STs   <- simulate_ST_paths(S0, r, q, rho, sim, v_mat, xi0_0, T_idx, antithetic=antithetic)
  names(STs) <- as.character(T_vec)
  
  # 3) Vektoriser prising per (type,T) over alle strikes
  split_keys <- paste(df$type, df$T)
  groups <- split(seq_len(nrow(df)), split_keys)
  
  out <- lapply(groups, function(ix){
    typ <- unique(df$type[ix]); Tm <- unique(df$T[ix]); K_vec <- df$K[ix]
    ST  <- STs[[as.character(Tm)]]
    Mmc <- length(ST)
    
    # Payoff-matrise (M x nK)
    if (typ == "call") {
      P <- pmax(matrix(ST, nrow=Mmc, ncol=length(K_vec)) -
                  matrix(K_vec, nrow=Mmc, ncol=length(K_vec), byrow=TRUE), 0)
    } else {
      P <- pmax(matrix(K_vec, nrow=Mmc, ncol=length(K_vec), byrow=TRUE) -
                  matrix(ST, nrow=Mmc, ncol=length(K_vec)), 0)
    }
    disc <- exp(-r * Tm)
    price <- disc * colMeans(P)
    se    <- disc * apply(P, 2, function(x) sd(x)/sqrt(Mmc))
    ci_lo <- price - 1.96*se
    ci_hi <- price + 1.96*se
    
    # BS-iv fra dine funksjoner; optional clipping til arb-grenser for robusthet
    disc_r <- exp(-r*Tm); disc_q <- exp(-q*Tm)
    lower <- if (typ=="call") pmax(0, S0*disc_q - K_vec*disc_r) else pmax(0, K_vec*disc_r - S0*disc_q)
    upper <- if (typ=="call") rep(S0*disc_q, length(K_vec)) else K_vec*disc_r
    P_iv  <- if (clip_iv) pmin(pmax(price, lower + 1e-10), upper - 1e-10) else price
    
    iv <- vapply(seq_along(K_vec),
                 function(j) bs_implied_vol(P_iv[j], S=S0, K=K_vec[j], T=Tm, r=r, q=q, type=typ),
                 numeric(1))
    
    data.frame(type=typ, K=K_vec, T=Tm, price=price, se=se, ci_lo=ci_lo, 
               ci_hi=ci_hi, iv=iv, ST=ST)
  })
  
  res <- do.call(rbind, out)
  res <- res[order(res$T, match(res$type, c("call","put")), res$K), ]
  rownames(res) <- NULL
  res
}

## Prising + IV for hele options_df
rb_res <- price_many_rB(options_df, S0, r, q,
                        alpha=alpha, eta=eta, rho=rho, xi0_level=xi0_level,
                        N=N, M=M, antithetic=TRUE, seed=seed)

mean(rb_res$ST)

