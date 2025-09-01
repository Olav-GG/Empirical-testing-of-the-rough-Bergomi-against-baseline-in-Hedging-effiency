# ===== Black???Scholes: kjerne =====
bs_price <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) return(if (type == "call") pmax(S - K, 0) else pmax(K - S, 0))
  if (sigma <= 0) {
    F <- S * exp((r - q) * T); disc <- exp(-r * T)
    return(if (type == "call") disc * pmax(F - K, 0) else disc * pmax(K - F, 0))
  }
  sqrtT <- sqrt(T)
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrtT)
  d2 <- d1 - sigma * sqrtT
  if (type == "call") S * exp(-q * T) * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
  else                K * exp(-r * T) * pnorm(-d2) - S * exp(-q * T) * pnorm(-d1)
}

.bs_d1d2 <- function(S, K, T, r, q, sigma) {
  sqrtT <- sqrt(T)
  d1 <- (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrtT)
  d2 <- d1 - sigma * sqrtT
  list(d1 = d1, d2 = d2)
}

bs_greeks <- function(S, K, T, r, q, sigma, type = c("call", "put")) {
  type <- match.arg(type)
  if (T <= 0) {
    intrinsic <- if (type == "call") pmax(S - K, 0) else pmax(K - S, 0)
    return(list(delta = if (type == "call") as.numeric(S > K) else -as.numeric(S < K),
                gamma = 0, vega = 0, theta = 0, rho = 0, price = intrinsic))
  }
  if (sigma <= 0) sigma <- .Machine$double.eps
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
  list(delta=delta, gamma=gamma, vega=vega, theta=theta, rho=rho, price=price)
}

bs_implied_vol <- function(mkt_price, S, K, T, r, q, type = c("call", "put"),
                           sigma_lo = 1e-8, sigma_hi = 5, tol = 1e-8, maxiter = 100) {
  type <- match.arg(type)
  if (T <= 0) return(NA_real_)
  disc_r <- exp(-r * T); disc_q <- exp(-q * T)
  lower <- if (type == "call") max(0, S * disc_q - K * disc_r) else max(0, K * disc_r - S * disc_q)
  upper <- if (type == "call") S * disc_q else K * disc_r
  if (mkt_price < lower - 1e-12 || mkt_price > upper + 1e-12) return(NA_real_)
  f <- function(s) {
    if (s <= 0) return(lower - mkt_price)  # tving funksjonen under
    sqrtT <- sqrt(T)
    d1 <- (log(S/K) + (r - q + 0.5*s^2)*T)/(s*sqrtT)
    d2 <- d1 - s*sqrtT
    price <- if (type == "call") S*disc_q*pnorm(d1) - K*disc_r*pnorm(d2)
    else                K*disc_r*pnorm(-d2) - S*disc_q*pnorm(-d1)
    price - mkt_price
  }
  flo <- f(sigma_lo); fhi <- f(sigma_hi); tries <- 0
  while (flo * fhi > 0 && tries < 5) { sigma_hi <- sigma_hi * 2; fhi <- f(sigma_hi); tries <- tries + 1 }
  if (flo * fhi > 0) return(NA_real_)
  uniroot(f, lower = sigma_lo, upper = sigma_hi, tol = tol, maxiter = maxiter)$root
}
# ===== slutt BS-kjerne =====


## --- hjelpe: beregn BS-data for ??n rad ---
.bs_row <- function(typ, K, T, S0, r, q, Sigma) {
  # Pris
  P <- bs_price(S=S0, K=K, T=T, r=r, q=q, sigma=Sigma, type=typ)
  # Greeks (fra din funksjon)
  g <- bs_greeks(S=S0, K=K, T=T, r=r, q=q, sigma=Sigma, type=typ)
  # d1/d2
  dd <- .bs_d1d2(S=S0, K=K, T=T, r=r, q=q, sigma=Sigma)
  # Arbitrasjegrenser
  disc_r <- exp(-r*T); disc_q <- exp(-q*T)
  lower <- if (typ=="call") max(0, S0*disc_q - K*disc_r) else max(0, K*disc_r - S0*disc_q)
  upper <- if (typ=="call") S0*disc_q else K*disc_r
  # Returner ??n rad som named vector/list
  list(
    price = as.numeric(P),
    delta = as.numeric(g$delta),
    gamma = as.numeric(g$gamma),
    vega  = as.numeric(g$vega),
    theta = as.numeric(g$theta),
    rho   = as.numeric(g$rho),
    d1    = as.numeric(dd$d1),
    d2    = as.numeric(dd$d2),
    Nd1   = pnorm(dd$d1),
    Nd2   = pnorm(dd$d2),
    lower = lower,
    upper = upper
  )
}

## --- hoved: pris alle rader i options_df ---
bs_price_many <- function(options_df, S0, r, q, Sigma, use_mkt_price_col = NULL) {
  stopifnot(all(c("type","K","T") %in% names(options_df)))
  df <- options_df
  df$type <- tolower(trimws(df$type))
  
  # Kj??r radvis med mapply; hver rad gir en liten liste -> bind sammen
  rows <- mapply(
    FUN = .bs_row,
    typ = df$type, K = df$K, T = df$T,
    MoreArgs = list(S0=S0, r=r, q=q, Sigma=Sigma),
    SIMPLIFY = FALSE
  )
  aux <- do.call(rbind, lapply(rows, as.data.frame))
  out <- cbind(df, aux)
  
  # Fyll ut noen avledete felt (om de ikke finnes fra f??r)
  if (!"F0"        %in% names(out)) out$F0        <- S0 * exp((r - q) * out$T)
  if (!"disc"      %in% names(out)) out$disc      <- exp(-r * out$T)
  if (!"div"       %in% names(out)) out$div       <- exp(-q * out$T)
  if (!"moneyness" %in% names(out)) out$moneyness <- S0 / out$K
  if (!"option_id" %in% names(out)) out$option_id <- sprintf("%s_K%d_T%.2f", out$type, out$K, out$T)
  
  # Implied vol:
  # - Hvis du har en kolonne med markedspris (f.eks. "mkt_price"), bruk den
  # - Ellers bruk modellprisen som test: IV skal da bli ~Sigma
  iv_base <- if (!is.null(use_mkt_price_col) && use_mkt_price_col %in% names(out)) out[[use_mkt_price_col]] else out$price
  out$iv <- mapply(function(typ, K, T, P)
    bs_implied_vol(mkt_price = P, S = S0, K = K, T = T, r = r, q = q, type = typ),
    out$type, out$K, out$T, iv_base)
  
  # Liten sanity: pris innen arbitrasjegrenser
  out$in_bounds <- (out$price >= out$lower - 1e-10 & out$price <= out$upper + 1e-10)
  
  # Sort??r pent
  out <- out[order(out$T, out$type, out$K), ]
  rownames(out) <- NULL
  out
}

bs_res <- bs_price_many(options_df, S0=S0, r=r, q=q, Sigma=Sigma)

# Se resultat
print(bs_res)

# Filtrer p?? f.eks. T=1.0
subset(bs_res, T == 1.0)

# Eksempel: ATM call (S0=K=1000, T=1, r=q=0) ??? fasit pris ??? 93.536155956
subset(bs_res, abs(T-1)<1e-12 & K==1000)[, c("type","price")]
subset(rb_res, abs(T - 1.0) < 1e-12)[, c("type","price")]
