# ===== Black-Scholes core functions =====
bs_price <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) {
    # Option at expiry: payoff only
    return(if (type == "call") max(S - K, 0) else max(K - S, 0))
  }
  if (sigma <= 0) {
    # Zero vol: option price = discounted intrinsic value
    F <- S * exp((r - q) * T); disc <- exp(-r * T)
    return(if (type == "call") disc * max(F - K, 0) else disc * max(K - F, 0))
  }
  sqrtT <- sqrt(T)
  d1 <- (log(S/K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrtT)
  d2 <- d1 - sigma * sqrtT
  if (type == "call") {
    return(S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2))
  } else {
    return(K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1))
  }
}

.bs_d1d2 <- function(S, K, T, r, q, sigma) {
  sqrtT <- sqrt(T)
  d1 <- (log(S/K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrtT)
  d2 <- d1 - sigma * sqrtT
  return(list(d1 = d1, d2 = d2))
}

bs_greeks <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) {
    # At expiry, delta = 1 (call if S>K) or -1 (put if S<K), other Greeks 0
    intrinsic <- if (type == "call") max(S - K, 0) else max(K - S, 0)
    return(list(delta = if (type == "call") as.numeric(S > K) else -as.numeric(S < K),
                gamma = 0, vega = 0, theta = 0, rho = 0, price = intrinsic))
  }
  if (sigma <= 0) sigma <- .Machine$double.eps  # avoid division by zero
  dd <- .bs_d1d2(S, K, T, r, q, sigma)
  d1 <- dd$d1; d2 <- dd$d2
  Nd1 <- pnorm(d1); Nd2 <- pnorm(d2); nd1 <- dnorm(d1)
  disc_r <- exp(-r * T); disc_q <- exp(-q * T); sqrtT <- sqrt(T)
  price <- if (type == "call") S * disc_q * Nd1 - K * disc_r * Nd2
  else                K * disc_r * pnorm(-d2) - S * disc_q * pnorm(-d1)
  delta <- if (type == "call") disc_q * Nd1 else disc_q * (Nd1 - 1)
  gamma <- disc_q * nd1 / (S * sigma * sqrtT)
  vega  <- S * disc_q * nd1 * sqrtT
  theta_year <- -(S * disc_q * nd1 * sigma)/(2 * sqrtT) +
    if (type == "call") (q * S * disc_q * Nd1 - r * K * disc_r * Nd2)
  else                (-q * S * disc_q * pnorm(-d1) + r * K * disc_r * pnorm(-d2))
  theta <- theta_year / 365
  rho   <- if (type == "call") (K * T * disc_r * Nd2) else (-K * T * disc_r * pnorm(-d2))
  return(list(delta=delta, gamma=gamma, vega=vega, theta=theta, rho=rho, price=price))
}

bs_implied_vol <- function(mkt_price, S, K, T, r, q, type = c("call", "put"),
                           sigma_lo = 1e-8, sigma_hi = 5, tol = 1e-8, maxiter = 100) {
  # Numerically invert the BS formula to find implied vol
  type <- match.arg(type)
  if (T <= 0) return(NA_real_)
  disc_r <- exp(-r * T); disc_q <- exp(-q * T)
  lower_bound <- if (type == "call") max(0, S * disc_q - K * disc_r)
  else                max(0, K * disc_r - S * disc_q)
  upper_bound <- if (type == "call") S * disc_q else K * disc_r
  if (mkt_price < lower_bound - 1e-12 || mkt_price > upper_bound + 1e-12) return(NA_real_)
  
  if (abs(mkt_price - lower_bound) <= 1e-12) return(0)
  if (abs(mkt_price - upper_bound) <= 1e-12) return(sigma_hi)
  f <- function(sigma) {
    if (sigma <= 0) return(lower_bound - mkt_price)  # force function negative if sigma=0
    sqrtT <- sqrt(T)
    d1 <- (log(S/K) + (r - q + 0.5*sigma^2)*T) / (sigma*sqrtT)
    d2 <- d1 - sigma*sqrtT
    price <- if (type == "call") S*disc_q*pnorm(d1) - K*disc_r*pnorm(d2)
    else                K*disc_r*pnorm(-d2) - S*disc_q*pnorm(-d1)
    return(price - mkt_price)
  }
  # Expand sigma_hi until root is bracketed
  flo <- f(sigma_lo); fhi <- f(sigma_hi); tries <- 0
  while (flo * fhi > 0 && tries < 5) {
    sigma_hi <- sigma_hi * 2; fhi <- f(sigma_hi); tries <- tries + 1
  }
  if (flo * fhi > 0) return(NA_real_)
  return(uniroot(f, lower = sigma_lo, upper = sigma_hi, tol = tol, maxiter = maxiter)$root)
}
# ===== end BS core =====

## --- Helper: calculate BS values for one row ---
.bs_row <- function(typ, K, T, S0, r, q, Sigma) {
  disc_r <- exp(-r*T); disc_q <- exp(-q*T); sqrtT <- sqrt(T)
  d1 <- (log(S0/K) + (r - q + 0.5*Sigma^2)*T) / (Sigma*sqrtT)
  d2 <- d1 - Sigma*sqrtT
  Nd1 <- pnorm(d1); Nd2 <- pnorm(d2); nd1 <- dnorm(d1)
  
  price <- if (typ=="call") S0*disc_q*Nd1 - K*disc_r*Nd2
  else               K*disc_r*pnorm(-d2) - S0*disc_q*pnorm(-d1)
  delta <- if (typ=="call") disc_q*Nd1 else disc_q*(Nd1 - 1)
  gamma <- disc_q*nd1/(S0*Sigma*sqrtT)
  vega  <- S0*disc_q*nd1*sqrtT
  theta_year <- -(S0*disc_q*nd1*Sigma)/(2*sqrtT) +
    if (typ=="call") (q*S0*disc_q*Nd1 - r*K*disc_r*Nd2)
  else              (-q*S0*disc_q*pnorm(-d1) + r*K*disc_r*pnorm(-d2))
  theta <- theta_year/365
  rho   <- if (typ=="call") (K*T*disc_r*Nd2) else (-K*T*disc_r*pnorm(-d2))
  
  list(price=as.numeric(price), delta=as.numeric(delta), gamma=as.numeric(gamma),
       vega=as.numeric(vega), theta=as.numeric(theta), rho=as.numeric(rho),
       d1=as.numeric(d1), d2=as.numeric(d2), Nd1=Nd1, Nd2=Nd2,
       lower = if (typ=="call") max(0, S0*disc_q - K*disc_r) else max(0, K*disc_r - S0*disc_q),
       upper = if (typ=="call") S0*disc_q else K*disc_r)
}


## --- Main: price all options in a DataFrame ---
bs_price_many <- function(options_df, S0, r, q, Sigma, use_mkt_price_col = NULL) {
  stopifnot(all(c("type","K","T") %in% names(options_df)))
  df <- options_df
  df$type <- tolower(trimws(df$type))
  
  # Apply .bs_row on each row using mapply; combine results into a data frame
  rows <- mapply(
    FUN = .bs_row,
    typ = df$type, K = df$K, T = df$T,
    MoreArgs = list(S0=S0, r=r, q=q, Sigma=Sigma),
    SIMPLIFY = FALSE
  )
  aux <- do.call(rbind, lapply(rows, as.data.frame))
  out <- cbind(df, aux)
  
  # Fill in some derived fields if not already present
  if (!"F0"        %in% names(out)) out$F0        <- S0 * exp((r - q) * out$T)
  if (!"disc"      %in% names(out)) out$disc      <- exp(-r * out$T)
  if (!"div"       %in% names(out)) out$div       <- exp(-q * out$T)
  if (!"moneyness" %in% names(out)) out$moneyness <- S0 / out$K
  if (!"option_id" %in% names(out)) out$option_id <- sprintf("%s_K%.0f_T%.4f", out$type, out$K, out$T)
  
  # Compute implied vol for each option:
  # - If a specific market price column is given (e.g. "mkt_price"), use that for IV calculation.
  # - Otherwise, use the model price (for sanity check, should equal Sigma).
  iv_input <- if (!is.null(use_mkt_price_col) && use_mkt_price_col %in% names(out)) {
    out[[use_mkt_price_col]]
  } else out$price
  out$iv <- mapply(function(typ, K, T, P) {
    bs_implied_vol(mkt_price = P, S = S0, K = K, T = T, r = r, q = q, type = typ)
  },
  out$type, out$K, out$T, iv_input)
  
  # Sanity check: ensure model price is within arbitrage bounds
  out$in_bounds <- (out$price >= out$lower - 1e-10 & out$price <= out$upper + 1e-10)
  
  # Sort output by maturity, type, and strike
  out <- out[order(out$T, out$type, out$K), ]
  rownames(out) <- NULL
  return(out)
}



