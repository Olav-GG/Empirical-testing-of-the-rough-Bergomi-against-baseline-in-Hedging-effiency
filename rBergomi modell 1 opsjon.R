simulate_Y_paths_hybrid_k1 <- function(M, N, T, H, b = 0.5, seed = NULL, antithetic = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  dt <- T / N; a <- H - 0.5; t <- (1:N) * dt
  k <- 2:N
  bk_star <- ((k^(a+1) - (k-1)^(a+1))/(a+1))^(1/a)
  w_tail  <- sqrt(2*H) * ( (bk_star * dt)^a )
  w_all   <- c(0, w_tail)
  
  c1 <- sqrt(2*H) * dt^a / (a + 1)
  c2 <- dt^(a + 0.5) * sqrt( 2*H * ( 1/(2*a + 1) - 1/(a + 1)^2 ) )
  
  if (antithetic) {
    M2 <- ceiling(M/2)
    dW_half <- matrix(rnorm(N * M2, sd = sqrt(dt)), nrow = N)
    Z_half  <- matrix(rnorm(N * M2),                  nrow = N)
    dW_all  <- cbind(dW_half, -dW_half)[, 1:M, drop = FALSE]
    Z_all   <- cbind(Z_half,  -Z_half )[ , 1:M, drop = FALSE]
  }
  
  Y <- matrix(NA_real_, N, M)
  dW_mat <- matrix(NA_real_, N, M)
  
  for (m in 1:M) {
    dW <- if (antithetic) dW_all[, m] else sqrt(dt) * rnorm(N)
    Z  <- if (antithetic) Z_all[, m]  else rnorm(N)
    dW_mat[, m] <- dW
    tail_conv <- as.numeric(stats::convolve(dW, rev(w_all), type = "open"))[1:N]
    Y[, m] <- c1 * dW + c2 * Z + tail_conv
  }
  list(t = t, Y = Y, dt = dt, H = H, dW = dW_mat)
}

sim_hyb <- simulate_Y_paths_hybrid_k1(M = M, N = N, T = T, H = H,
                                      seed = seed, antithetic = TRUE)
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


#Prisings del:
rb_price_call_mc <- function(S0, K, r, q, rho, sim, v_mat, xi0_0, 
                             antithetic = FALSE, seed) {
  t  <- sim$t; dt <- sim$dt
  N  <- length(t)
  stopifnot(nrow(v_mat) == N, ncol(v_mat) == ncol(sim$dW))
  M  <- ncol(v_mat)
  sqrt_dt <- sqrt(dt)
  
  # Venstre-endepunkt: v_prev[i,] brukes for intervallet (t_{i-1}, t_i]
  # F??rste skritt bruker v(0) = xi0(0). For flat xi0 er det bare xi0_0.
  v_prev <- rbind(rep(xi0_0, M), v_mat[1:(N-1), , drop = FALSE])
  
  logS <- rep(log(S0), M)
  
  # Valgfritt: antithetisk for den uavhengige delen av dB
  if (antithetic) {
    M2 <- ceiling(M/2)
    Zp <- matrix(rnorm(N * M2), nrow = N)
    Z  <- cbind(Zp, -Zp)[, 1:M, drop = FALSE]
  }
  
  for (i in 1:N) {
    dW_perp <- if (antithetic) sqrt_dt * Z[i, ] else sqrt_dt * rnorm(M)
    dB_i    <- rho * sim$dW[i, ] + sqrt(1 - rho^2) * dW_perp
    
    vi <- pmax(v_prev[i, ], 0)  # numerisk trygghet
    logS <- logS + (r - q - 0.5 * vi) * dt + sqrt(vi) * dB_i
  }
  
  ST <- exp(logS)
  disc <- exp(-r * t[N])
  price <- disc * mean(pmax(ST - K, 0))
  list(price = price, ST = ST, mean_ST = mean(ST))
}

res <- rb_price_call_mc(S0 = S0, K = K, r = r, q = q,
                        rho = Rho, sim = sim_hyb, v_mat = v,
                        xi0_0 = xi0_level, antithetic = TRUE, seed = seed)

res$price

c(n = length(res$ST), positives = sum(res$ST > K), mean_ST = mean(res$ST))
summary(res$ST)
mean(pmax(res$ST - K, 0))  # skal v??re > 0


disc <- exp(-r * tail(sim_hyb$t, 1))
payoff <- disc * pmax(res$ST - K, 0)
mc_se <- sd(payoff) / sqrt(length(payoff))
c(price = res$price,
  se = mc_se,
  ci_lo = res$price - 1.96 * mc_se,
  ci_hi = res$price + 1.96 * mc_se)

