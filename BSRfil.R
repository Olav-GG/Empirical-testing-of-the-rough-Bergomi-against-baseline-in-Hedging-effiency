
bs_price <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) return(if (type == "call") pmax(S - K, 0) else pmax(K - S, 0))
  if (sigma <= 0) {
    F <- S * exp((r - q) * T); disc <- exp(-r * T)
    return(if (type == "call") disc * pmax(F - K, 0) else disc * pmax(K - F, 0))
  }
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  if (type == "call") {
    S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  } else {
    K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)
  }
}

bs_call_price <- bs_price(S, option$K, option$T, option$r, option$q, 
                          Sigma, type = option$type)
bs_call_price


.bs_d1d2 <- function(S, K, T, r, q, sigma) {
  sqrtT <- sqrt(T)
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrtT)
  d2 <- d1 - sigma * sqrtT
  list(d1 = d1, d2 = d2)
}


bs_greeks <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) {
    # Ved forfall: alt er stykkevis, sett konservativt til grenser
    intrinsic <- if (type == "call") pmax(S - K, 0) else pmax(K - S, 0)
    return(list(delta = if (type == "call") as.numeric(S > K) 
                else -as.numeric(S < K),
                gamma = 0, vega = 0, theta = 0, rho = 0, price = intrinsic))
  }
  if (sigma <= 0) sigma <- .Machine$double.eps
  
  dd <- .bs_d1d2(S, K, T, r, q, sigma)
  d1 <- dd$d1; d2 <- dd$d2
  Nd1 <- pnorm(d1); Nd2 <- pnorm(d2)
  nd1 <- dnorm(d1)
  
  disc_r <- exp(-r * T); disc_q <- exp(-q * T)
  sqrtT <- sqrt(T)
  
  # Pris (nyttig ?? f?? ut samtidig)
  price <- if (type == "call") {
    S * disc_q * Nd1 - K * disc_r * Nd2
  } else {
    K * disc_r * pnorm(-d2) - S * disc_q * pnorm(-d1)
  }
  
  delta <- if (type == "call") disc_q * Nd1 else disc_q * (Nd1 - 1)
  gamma <- disc_q * nd1 / (S * sigma * sqrtT)
  vega  <- S * disc_q * nd1 * sqrtT                 # per 1.00 i sigma
  theta_year <- -(S * disc_q * nd1 * sigma) / (2 * sqrtT) +
    if (type == "call") (q * S * disc_q * Nd1 - r * K * disc_r * Nd2)
  else                (-q * S * disc_q * pnorm(-d1) + r * K * disc_r * pnorm(-d2))
  theta <- theta_year/365
  rho   <- if (type == "call") (K * T * disc_r * Nd2) 
  else (-K * T * disc_r * pnorm(-d2))
  
  list(delta = delta, gamma = gamma, vega = vega, theta = theta, rho = rho, price = price)
}

# Implisitt volatilitet via robust uniroot (Brent)
bs_implied_vol <- function(mkt_price, S, K, T, r, q, type = c("call", "put"),
                           sigma_lo = 1e-8, sigma_hi = 5, tol = 1e-8, maxiter = 100) {
  type <- match.arg(type)
  if (T <= 0) return(NA_real_)
  
  disc_r <- exp(-r * T); disc_q <- exp(-q * T)
  
  # Enkle (arbitrage) grenser
  lower <- if (type == "call") max(0, S * disc_q - K * disc_r) else max(0, K * disc_r - S * disc_q)
  upper <- if (type == "call") S * disc_q else K * disc_r
  
  if (mkt_price < lower - 1e-12 || mkt_price > upper + 1e-12) return(NA_real_)
  
  f <- function(s) {
    # BS pris minus markedspris
    d1 <- (log(S / K) + (r - q + 0.5 * s^2) * T) / (s * sqrt(T))
    d2 <- d1 - s * sqrt(T)
    price <- if (type == "call") {
      S * disc_q * pnorm(d1) - K * disc_r * pnorm(d2)
    } else {
      K * disc_r * pnorm(-d2) - S * disc_q * pnorm(-d1)
    }
    price - mkt_price
  }
  
  # S??rg for at intervallet flankerer roten; utvid om n??dvendig
  flo <- f(sigma_lo); fhi <- f(sigma_hi)
  tries <- 0
  while (flo * fhi > 0 && tries < 5) {
    sigma_hi <- sigma_hi * 2
    fhi <- f(sigma_hi)
    tries <- tries + 1
  }
  if (flo * fhi > 0) return(NA_real_)
  
  uniroot(f, lower = sigma_lo, upper = sigma_hi, tol = tol, maxiter = maxiter)$root
}

# Greeks for opsjonen din (ATM call)
g <- bs_greeks(S, option$K, option$T, option$r, option$q, Sigma, option$type)
g$price                      # skal matche BS-prisen (??? 22.8715)
g$delta; g$gamma; g$vega; g$theta; g$rho

# Implisitt vol fra en markedspris (bruk BS-prisen som test -> b??r gi ~Sigma)
iv <- bs_implied_vol(mkt_price = g$price, S = S, K = option$K, T = option$T,
                     r = option$r, q = option$q, type = option$type)
iv
