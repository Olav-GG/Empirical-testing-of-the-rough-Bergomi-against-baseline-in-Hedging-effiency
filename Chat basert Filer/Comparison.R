# ====================== MULTI-DATE OOS COMPARISON with rB T-cache ======================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(digest)
})

# --- OOS date selection (assumes you have date_key from calibration) ---
all_dates     <- sort(names(spxIvolList))
date_keys_oos <- all_dates[all_dates > date_key]

# -------- sanity: required core objects/functions --------
stopifnot(exists("spxIvolList"), exists("price_many_rB_turbo"), exists("build_calibration_df"))

# -------- helpers: metrics, buckets, HTML --------
RMSE <- function(x) { y <- x[is.finite(x)]; if (!length(y)) NA_real_ else sqrt(mean(y^2)) }
MAE  <- function(x) { y <- x[is.finite(x)]; if (!length(y)) NA_real_ else mean(abs(y)) }
fmt  <- function(x, d=6) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "")

# Accept both 'left' and legacy 'left_align'
to_html <- function(df, left=c(), left_align=NULL) {
  if (!is.null(left_align)) left <- left_align
  heads <- vapply(names(df), function(nm) {
    cls <- if (nm %in% left) "left" else "right"
    sprintf("<th class='%s'>%s</th>", cls, nm)
  }, character(1))
  rows <- apply(df, 1, function(row) {
    paste0("<tr>", paste0(
      mapply(function(val, j) {
        cls <- if (names(df)[j] %in% left) "left" else "right"
        sprintf("<td class='%s'>%s</td>", cls, ifelse(is.na(val) || val=="NA","",as.character(val)))
      }, row, seq_along(row)), collapse=""), "</tr>")
  })
  paste0("<table><thead><tr>", paste0(heads, collapse=""),
         "</tr></thead><tbody>", paste0(rows, collapse=""), "</tbody></table>")
}

T_bucket <- function(T) {
  cut(T,
      breaks = c(0, 1/12, 0.25, 0.5, 1, 2, 5, Inf),
      labels = c("???1M","1???3M","3???6M","6???12M","1???2Y","2???5Y",">5Y"),
      right = TRUE, include.lowest = TRUE)
}

# -------- Black???Scholes helpers --------
bs_price <- function(S, K, T, r, q, sigma, type=c("call","put")) {
  type <- match.arg(type)
  if (!is.finite(T) || T<=0 || !is.finite(sigma) || sigma<=0) return(NA_real_)
  d1 <- (log(S/K) + (r - q + 0.5*sigma^2)*T) / (sigma*sqrt(T))
  d2 <- d1 - sigma*sqrt(T)
  if (type=="call") S*exp(-q*T)*pnorm(d1) - K*exp(-r*T)*pnorm(d2)
  else              K*exp(-r*T)*pnorm(-d2) - S*exp(-q*T)*pnorm(-d1)
}

bs_price_many <- function(options_df, S0, r, q, Sigma) {
  stopifnot(all(c("type","K","T") %in% names(options_df)))
  df <- options_df
  df$type <- tolower(trimws(df$type))
  df$price <- mapply(function(tp,K,T) bs_price(S0, K, T, r, q, Sigma, type=tp),
                     df$type, df$K, df$T)
  df
}

# Robust implied vol if none present
if (!exists("bs_implied_vol")) {
  bs_implied_vol <- function(mkt_price, S, K, T, r, q, type=c("call","put"), lo=1e-6, hi=5) {
    type <- match.arg(type)
    if (!is.finite(mkt_price) || !is.finite(S) || !is.finite(K) || !is.finite(T) || T<=0) return(NA_real_)
    disc_r <- exp(-r*T); disc_q <- exp(-q*T)
    lower  <- if (type=="call") max(0, S*disc_q - K*disc_r) else max(0, K*disc_r - S*disc_q)
    upper  <- if (type=="call") S*disc_q else K*disc_r
    P <- mkt_price
    if (P <= lower + 1e-12 || P >= upper - 1e-12) return(NA_real_)
    f <- function(sig) bs_price(S,K,T,r,q,sig,type) - P
    f_lo <- f(lo); f_hi <- f(hi)
    if (!is.finite(f_lo) || !is.finite(f_hi) || f_lo*f_hi > 0) {
      for (h in c(7,10)) { f_hi <- f(h); if (is.finite(f_hi) && f_lo*f_hi <= 0) { hi <- h; break } }
      if (!is.finite(f_hi) || f_lo*f_hi > 0) return(NA_real_)
    }
    out <- try(uniroot(f, c(lo,hi), tol=1e-10, maxiter=200), silent=TRUE)
    if (inherits(out,"try-error")) NA_real_ else out$root
  }
}

safe_bs_iv_one <- function(type, price, S, K, T, r, q) {
  if (!is.finite(price) || !is.finite(K) || !is.finite(T) || T <= 0) return(NA_real_)
  disc_r <- exp(-r*T); disc_q <- exp(-q*T)
  lower  <- if (type=="call") max(0, S*disc_q - K*disc_r) else max(0, K*disc_r - S*disc_q)
  upper  <- if (type=="call") S*disc_q else K*disc_r
  P <- min(max(price, lower+1e-12), upper-1e-12)
  out <- try(bs_implied_vol(mkt_price=P, S=S, K=K, T=T, r=r, q=q, type=type), silent=TRUE)
  if (inherits(out,"try-error")) NA_real_ else as.numeric(out)
}
safe_bs_iv_vec <- function(types, prices, S, Ks, Ts, r, q)
  mapply(safe_bs_iv_one, type=types, price=prices, S=S, K=Ks, T=Ts, MoreArgs=list(r=r, q=q))

# -------- build xi0 from ATM IVs (isotonic total variance) --------
build_xi0_from_ATM <- function(calib_df, extend_to = NULL) {
  req <- c("T","K","market_iv","fwd")
  stopifnot(all(req %in% names(calib_df)))
  df <- calib_df %>%
    transmute(T = as.numeric(.data$T),
              K = as.numeric(.data$K),
              market_iv = as.numeric(.data$market_iv),
              fwd = as.numeric(.data$fwd)) %>%
    filter(is.finite(T), T>0, is.finite(K), is.finite(market_iv), is.finite(fwd))
  stopifnot(nrow(df) > 0)
  
  atm_by_T <- df %>%
    group_by(T) %>%
    slice_min(order_by = abs(K - median(fwd)), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(T) %>%
    mutate(total_var = (market_iv^2) * T)
  stopifnot(nrow(atm_by_T) >= 1, min(atm_by_T$T) > 0)
  
  tv_iso <- stats::isoreg(x = atm_by_T$T, y = atm_by_T$total_var)$yf
  eps <- 1e-12
  T_nodes   <- atm_by_T$T
  fwd_first <- tv_iso[1] / T_nodes[1]
  fwd_rest  <- diff(tv_iso) / pmax(diff(T_nodes), eps)
  fwd_all   <- pmax(c(fwd_first, fwd_rest), 0)
  
  xi0 <- data.frame(
    T_start     = c(0, head(T_nodes, -1)),
    T_end       = T_nodes,
    forward_var = fwd_all
  )
  if (!is.null(extend_to) && is.finite(extend_to) && extend_to > max(xi0$T_end)) {
    xi0 <- dplyr::bind_rows(
      xi0,
      data.frame(T_start = max(xi0$T_end), T_end = extend_to,
                 forward_var = tail(xi0$forward_var, 1))
    )
  }
  xi0
}

# ============================ rB PRICE CACHE =============================
T_round <- function(x) round(as.numeric(x), 10)

# in-memory cache: name -> list(prices=data.frame, meta=list(...))
.rB_cache <- new.env(parent = emptyenv())

build_S0_from_chain <- function(chain, r, q) {
  T_min <- min(chain$T[is.finite(chain$T)], na.rm=TRUE)
  near  <- chain %>% dplyr::filter(abs(T - T_min) < 1e-6 & is.finite(F))
  if (nrow(near) == 0) near <- chain %>% dplyr::filter(is.finite(F))
  median(near$F * exp(-(r - q)*near$T), na.rm=TRUE)
}

global_T_max <- function(keys, spxIvolList) {
  max(vapply(keys, function(d) {
    ch <- spxIvolList[[d]]
    if (!is.data.frame(ch) || !nrow(ch)) return(NA_real_)
    suppressWarnings(max(as.numeric(ch$Texp), na.rm=TRUE))
  }, numeric(1)), na.rm=TRUE)
}

# Build a stable key that includes params + xi0 digest
.cache_key <- function(d, H, Rho, Eta, N, M, seed, S0, r, q, xi0_df) {
  paste(
    d, round(H,6), round(Rho,6), round(Eta,6), N, M, seed,
    round(S0,6), r, q, digest(xi0_df), sep="|"
  )
}

# Save/Load cache to disk (optional)
save_rb_cache <- function(path="rb_cache.rds") {
  lst <- as.list.environment(.rB_cache, all.names = TRUE)
  saveRDS(lst, path)
}
load_rb_cache <- function(path="rb_cache.rds") {
  if (!file.exists(path)) return(invisible(FALSE))
  lst <- readRDS(path)
  rm(list = ls(.rB_cache, all.names=TRUE), envir = .rB_cache)
  for (nm in names(lst)) assign(nm, lst[[nm]], envir = .rB_cache)
  invisible(TRUE)
}

# Pre-price whole date d (vectorized) and store
preprice_rb_for_date <- function(d, chain_raw, r, q, H, Rho, Eta, N=500L, M=10000L,
                                 seed_base=1L, extend_to=NULL) {
  stopifnot(is.data.frame(chain_raw), nrow(chain_raw) > 0)
  
  chain <- chain_raw %>%
    mutate(T = as.numeric(Texp), K = as.numeric(Strike), F = as.numeric(Fwd),
           CallMid = suppressWarnings(as.numeric(CallMid)))
  
  # S0 and option table (calls + puts via parity)
  S0_today <- build_S0_from_chain(chain, r, q)
  
  oos <- chain %>%
    filter(is.finite(T), T > 0, is.finite(K), is.finite(CallMid), is.finite(F)) %>%
    transmute(type="call", K, T=T_round(T), mkt_price = CallMid, date=d) %>%
    bind_rows(
      chain %>%
        filter(is.finite(T), T > 0, is.finite(K), is.finite(CallMid), is.finite(F)) %>%
        transmute(type="put",  K, T=T_round(T),
                  mkt_price = CallMid - exp(-r*T)*(F - K), date=d)
    )
  if (!nrow(oos)) return(NULL)
  
  # xi0 curve up to extend_to (global T_max recommended)
  calib_df_today <- build_calibration_df(
    spxIvolList, date_key = d,
    grid_list = lapply(unique(chain$Expiry),
                       function(e) list(expiration=e, ratios=c(0.8,0.9,1.0,1.1,1.2))),
    S0 = S0_today, r = r, q = q, use_chain_fwd = TRUE
  )
  xi0_curve_today <- build_xi0_from_ATM(
    calib_df_today,
    extend_to = if (is.null(extend_to)) max(oos$T, na.rm=TRUE) else extend_to
  )
  
  # Vectorized rB pricing for all (type,K,T) that day
  seed_i <- as.integer(abs((seed_base + sum(utf8ToInt(d))*131) %% .Machine$integer.max))
  rB_tbl <- price_many_rB_turbo(
    options_df = oos[, c("type","K","T")],
    S0 = S0_today, r = r, q = q,
    alpha = H - 0.5, eta = Eta, rho = Rho,
    xi0_curve = xi0_curve_today,
    N = N, M = M, antithetic = TRUE, seed = seed_i, clip_iv = TRUE
  ) %>%
    mutate(type = tolower(trimws(type)),
           K = as.numeric(K), T = T_round(as.numeric(T)),
           price_rB = price, .keep = "unused")
  
  meta <- list(
    date=d, H=H, Rho=Rho, Eta=Eta, N=N, M=M, seed=seed_i,
    S0=S0_today, r=r, q=q, xi0_digest = digest(xi0_curve_today)
  )
  
  key <- .cache_key(d, H, Rho, Eta, N, M, seed_i, S0_today, r, q, xi0_curve_today)
  assign(key, list(prices = rB_tbl, meta = meta), envir = .rB_cache)
  
  invisible(list(prices=rB_tbl, meta=meta, key=key))
}

# Retrieve cached prices for wanted (type,K,T); recompute if missing or params differ
get_rb_prices <- function(d, wanted_df, chain_raw, r, q, H, Rho, Eta, N=500L, M=10000L,
                          seed_base=1L, extend_to=NULL) {
  # Try to find any cache entry for this date d
  keys <- ls(.rB_cache, all.names = TRUE)
  cand <- grep(paste0("^", d, "\\|"), keys, value = TRUE)
  
  pick <- function() {
    if (!length(cand)) return(NULL)
    # pick the last inserted candidate
    get(tail(cand, 1), envir = .rB_cache)
  }
  
  entry <- pick()
  
  # If not found or params mismatch, build afresh
  needs_rebuild <- is.null(entry) ||
    is.null(entry$meta) ||
    any(abs(c(entry$meta$H, entry$meta$Rho, entry$meta$Eta) - c(H,Rho,Eta)) > 1e-12) ||
    entry$meta$N != N || entry$meta$M != M || entry$meta$r != r || entry$meta$q != q
  
  if (needs_rebuild) {
    pre <- preprice_rb_for_date(
      d, chain_raw, r, q, H, Rho, Eta, N, M, seed_base, extend_to
    )
    if (is.null(pre)) return(NULL)
    entry <- get(pre$key, envir = .rB_cache)
  }
  
  # Join cached prices
  wanted_df <- wanted_df %>% mutate(type=tolower(trimws(type)), K=as.numeric(K), T=T_round(T))
  out <- wanted_df %>% left_join(entry$prices, by=c("type","K","T"))
  
  # If any missing (new strikes/tenors), price just those and update the cache
  if (anyNA(out$price_rB)) {
    missing <- out %>% filter(!is.finite(price_rB)) %>% distinct(type,K,T)
    
    # Recreate S0 and xi0 for the day using the same meta seed and params
    chain <- chain_raw %>%
      mutate(T = as.numeric(Texp), K = as.numeric(Strike), F = as.numeric(Fwd),
             CallMid = suppressWarnings(as.numeric(CallMid)))
    S0_today <- build_S0_from_chain(chain, r, q)
    calib_df_today <- build_calibration_df(
      spxIvolList, date_key = d,
      grid_list = lapply(unique(chain$Expiry),
                         function(e) list(expiration=e, ratios=c(0.8,0.9,1.0,1.1,1.2))),
      S0 = S0_today, r = r, q = q, use_chain_fwd = TRUE
    )
    xi0_curve_today <- build_xi0_from_ATM(
      calib_df_today,
      extend_to = if (is.null(extend_to)) max(chain$T, na.rm=TRUE) else extend_to
    )
    
    extra <- price_many_rB_turbo(
      options_df = missing,
      S0 = S0_today, r = r, q = q,
      alpha = H - 0.5, eta = Eta, rho = Rho,
      xi0_curve = xi0_curve_today,
      N = N, M = M, antithetic = TRUE, seed = entry$meta$seed, clip_iv = TRUE
    ) %>% transmute(type=tolower(trimws(type)), K=as.numeric(K), T=T_round(as.numeric(T)),
                    price_rB = price)
    
    # update cache and output
    updated <- bind_rows(entry$prices, extra) %>% distinct(type,K,T,.keep_all=TRUE)
    entry$prices <- updated
    assign(ls(envir=.rB_cache, pattern=paste0("^", d, "\\|")) |> tail(1), entry, envir=.rB_cache)
    
    out <- wanted_df %>% left_join(updated, by=c("type","K","T"))
  }
  
  out
}

# ============================ MAIN OOS COMPARISON =============================
oos_compare_many <- function(date_keys, spxIvolList, r, q,
                             H_star, Rho_star, Eta_star,
                             N = 500L, M = 10000L,
                             seed_base = 1L,
                             outfile = "oos_multi_summary.html",
                             calib_date = get0("date_key", ifnotfound="(unknown)"),
                             cache_path = NULL # set e.g. "rb_cache.rds" to persist
) {
  # Optional: load existing cache
  if (!is.null(cache_path)) load_rb_cache(cache_path)
  
  # Global T_max for xi0 extension (covers the longest maturity in the set)
  Tmax_global <- global_T_max(date_keys, spxIvolList)
  if (!is.finite(Tmax_global)) Tmax_global <- NA_real_
  
  per_date <- function(d) {
    chain <- spxIvolList[[d]]
    if (!is.data.frame(chain) || !nrow(chain)) return(NULL)
    
    # Numerics
    chain <- chain %>%
      mutate(T = as.numeric(Texp),
             K = as.numeric(Strike),
             F = as.numeric(Fwd),
             CallMid = suppressWarnings(as.numeric(CallMid)))
    
    # Build OOS option table (calls + puts via parity)
    df_oos <- chain %>%
      filter(is.finite(T), T > 0, is.finite(K), is.finite(CallMid), is.finite(F)) %>%
      mutate(call_px = CallMid,
             put_px  = call_px - exp(-r * T) * (F - K))
    if (!nrow(df_oos)) return(NULL)
    
    calls_df <- df_oos %>% transmute(type="call", K, T=T_round(T), mkt_price=call_px, date=d)
    puts_df  <- df_oos %>% transmute(type="put",  K, T=T_round(T), mkt_price=put_px,  date=d)
    options_df_oos <- bind_rows(calls_df, puts_df)
    if (!nrow(options_df_oos)) return(NULL)
    
    # -------- rB prices via cache (preprices whole date once) --------
    rB_prices <- get_rb_prices(
      d = d, wanted_df = options_df_oos[, c("type","K","T")],
      chain_raw = chain,
      r = r, q = q, H = H_star, Rho = Rho_star, Eta = Eta_star,
      N = N, M = M, seed_base = seed_base,
      extend_to = if (is.finite(Tmax_global)) Tmax_global else max(options_df_oos$T, na.rm=TRUE)
    )
    if (is.null(rB_prices)) return(NULL)
    
    # -------- Black???Scholes with shortest-bucket ATM vol --------
    S0_today <- build_S0_from_chain(chain, r, q)
    calib_df_today <- build_calibration_df(
      spxIvolList, date_key = d,
      grid_list = lapply(unique(chain$Expiry),
                         function(e) list(expiration=e, ratios=c(0.8,0.9,1.0,1.1,1.2))),
      S0 = S0_today, r = r, q = q, use_chain_fwd = TRUE
    )
    xi0_curve_today <- build_xi0_from_ATM(
      calib_df_today, extend_to = if (is.finite(Tmax_global)) Tmax_global else max(options_df_oos$T, na.rm=TRUE)
    )
    Sigma_BS <- sqrt(max(xi0_curve_today$forward_var[1], 1e-12))
    
    BS_prices <- bs_price_many(
      options_df = options_df_oos[, c("type","K","T")],
      S0 = S0_today, r = r, q = q, Sigma = Sigma_BS
    ) %>% mutate(type=tolower(trimws(type)), K=as.numeric(K), T=T_round(as.numeric(T)),
                 price_BS = price, .keep="unused")
    
    # Market IV + errors
    mkt_iv <- safe_bs_iv_vec(options_df_oos$type, options_df_oos$mkt_price,
                             S=S0_today, Ks=options_df_oos$K, Ts=options_df_oos$T, r=r, q=q)
    
    cmp <- options_df_oos %>%
      mutate(mkt_iv = mkt_iv,
             type   = tolower(trimws(type)),
             K      = as.numeric(K),
             T      = T_round(T)) %>%
      left_join(rB_prices, by=c("type","K","T")) %>%
      left_join(BS_prices, by=c("type","K","T")) %>%
      mutate(
        iv_rB = safe_bs_iv_vec(type, price_rB, S=S0_today, Ks=K, Ts=T, r=r, q=q),
        iv_BS = safe_bs_iv_vec(type, price_BS, S=S0_today, Ks=K, Ts=T, r=r, q=q),
        errP_rB  = price_rB - mkt_price,
        errP_BS  = price_BS - mkt_price,
        errIV_rB = iv_rB - mkt_iv,
        errIV_BS = iv_BS - mkt_iv
      )
    
    daily <- cmp %>% summarise(
      date = first(date),
      N_options      = sum(is.finite(mkt_price)),
      IV_RMSE_rB     = RMSE(errIV_rB),
      IV_RMSE_BS     = RMSE(errIV_BS),
      Price_RMSE_rB  = RMSE(errP_rB),
      Price_RMSE_BS  = RMSE(errP_BS),
      IV_MSE_rB      = mean(errIV_rB[is.finite(errIV_rB)]^2),
      IV_MSE_BS      = mean(errIV_BS[is.finite(errIV_BS)]^2),
      Price_MSE_rB   = mean(errP_rB[is.finite(errP_rB)]^2),
      Price_MSE_BS   = mean(errP_BS[is.finite(errP_BS)]^2),
      WinRate_IV     = mean(abs(errIV_rB) < abs(errIV_BS), na.rm=TRUE),
      dIV_RMSE       = IV_RMSE_rB - IV_RMSE_BS,
      dPrice_RMSE    = Price_RMSE_rB  - Price_RMSE_BS
    )
    
    list(cmp=cmp, daily=daily)
  }
  
  message("Processing ", length(date_keys), " date(s) ...")
  results_list <- lapply(date_keys, per_date)
  results_list <- Filter(Negate(is.null), results_list)
  if (!length(results_list)) stop("No valid out-of-sample dates produced results.")
  
  daily_all <- bind_rows(lapply(results_list, `[[`, "daily"))
  cmp_all   <- bind_rows(lapply(results_list, `[[`, "cmp"))
  
  # Pooled metrics
  pooled <- cmp_all %>% summarise(
    N_options      = sum(is.finite(mkt_price)),
    IV_RMSE_rB     = RMSE(errIV_rB),   IV_RMSE_BS     = RMSE(errIV_BS),
    Price_RMSE_rB  = RMSE(errP_rB),    Price_RMSE_BS  = RMSE(errP_BS),
    IV_MAE_rB      = MAE(errIV_rB),    IV_MAE_BS      = MAE(errIV_BS),
    Price_MAE_rB   = MAE(errP_rB),     Price_MAE_BS   = MAE(errP_BS),
    WinRate_IV     = mean(abs(cmp_all$errIV_rB) < abs(cmp_all$errIV_BS), na.rm=TRUE)
  ) %>%
    mutate(
      dIV_RMSE      = IV_RMSE_rB - IV_RMSE_BS,
      dPrice_RMSE   = Price_RMSE_rB - Price_RMSE_BS,
      IV_impr_pct   = 100 * (1 - IV_RMSE_rB / IV_RMSE_BS),
      Price_impr_pct= 100 * (1 - Price_RMSE_rB / Price_RMSE_BS)
    )
  
  # Stats across days
  iv_pair_idx <- with(daily_all, is.finite(IV_MSE_rB) & is.finite(IV_MSE_BS))
  pr_pair_idx <- with(daily_all, is.finite(Price_MSE_rB) & is.finite(Price_MSE_BS))
  t_iv <- if (sum(iv_pair_idx) >= 2) t.test(daily_all$IV_MSE_rB[iv_pair_idx],
                                            daily_all$IV_MSE_BS[iv_pair_idx],
                                            paired=TRUE, alternative="less") else NULL
  t_pr <- if (sum(pr_pair_idx) >= 2) t.test(daily_all$Price_MSE_rB[pr_pair_idx],
                                            daily_all$Price_MSE_BS[pr_pair_idx],
                                            paired=TRUE, alternative="less") else NULL
  valid_days <- daily_all$dIV_RMSE[is.finite(daily_all$dIV_RMSE)]
  n_days <- length(valid_days); n_wins <- sum(valid_days < 0)
  p_binom <- if (n_days > 0) pbinom(q = n_wins - 1, size = n_days, prob = 0.5, lower.tail = FALSE) else NA_real_
  
  # HTML report
  score_df <- tibble(
    Metric            = c("IV RMSE", "Price RMSE", "IV MAE", "Price MAE",
                          "IV WinRate (|err|)", "IV Improvement %", "Price Improvement %"),
    rB                = c(pooled$IV_RMSE_rB, pooled$Price_RMSE_rB,
                          pooled$IV_MAE_rB, pooled$Price_MAE_rB,
                          pooled$WinRate_IV, pooled$IV_impr_pct/100, pooled$Price_impr_pct/100),
    BS                = c(pooled$IV_RMSE_BS, pooled$Price_RMSE_BS,
                          pooled$IV_MAE_BS, pooled$Price_MAE_BS,
                          1 - pooled$WinRate_IV, NA, NA),
    "?? (rB???BS)"       = c(pooled$dIV_RMSE, pooled$dPrice_RMSE,
                          pooled$IV_MAE_rB - pooled$IV_MAE_BS,
                          pooled$Price_MAE_rB - pooled$Price_MAE_BS,
                          NA, NA, NA)
  )
  score_html <- to_html(score_df %>% mutate(across(where(is.numeric), ~fmt(.))), left="Metric")
  
  daily_tbl <- daily_all %>% mutate(
    IV_RMSE_rB   = fmt(IV_RMSE_rB),  IV_RMSE_BS   = fmt(IV_RMSE_BS),  dIV_RMSE   = fmt(dIV_RMSE),
    Price_RMSE_rB= fmt(Price_RMSE_rB), Price_RMSE_BS= fmt(Price_RMSE_BS), dPrice_RMSE = fmt(dPrice_RMSE),
    WinRate_IV   = fmt(WinRate_IV, 4)
  )
  daily_html <- to_html(daily_tbl, left="date")
  
  byT_df <- if (nrow(cmp_all)) {
    cmp_all %>% mutate(T_bucket = T_bucket(T)) %>%
      group_by(date, T_bucket) %>%
      summarise(dIV_RMSE = RMSE(errIV_rB) - RMSE(errIV_BS), .groups="drop") %>%
      group_by(T_bucket) %>%
      summarise(n_dates = n(),
                median_dIV_RMSE = median(dIV_RMSE, na.rm=TRUE),
                mean_dIV_RMSE   = mean(dIV_RMSE, na.rm=TRUE),
                winrate_days    = mean(dIV_RMSE < 0, na.rm=TRUE),
                .groups="drop")
  } else {
    tibble(T_bucket=factor(levels=T_bucket(1)), n_dates=0, median_dIV_RMSE=NA_real_,
           mean_dIV_RMSE=NA_real_, winrate_days=NA_real_)
  }
  byT_html <- to_html(byT_df %>% mutate(across(where(is.numeric), ~fmt(.))), left="T_bucket")
  
  concl  <- sprintf(
    "Across %d out-of-sample dates, pooled IV RMSE: rB = %.4f vs BS = %.4f (%.1f%% lower with rB).",
    n_days, pooled$IV_RMSE_rB, pooled$IV_RMSE_BS, pooled$IV_impr_pct
  )
  concl2 <- sprintf("Paired t-test on IV MSE differences gives %s; sign test: %s.",
                    if(!is.null(t_iv)) sprintf("p=%.3g", t_iv$p.value) else "p=N/A",
                    if(is.finite(p_binom)) sprintf("%d/%d days with rB better (p=%.3g)", n_wins, n_days, p_binom)
                    else "insufficient data")
  concl3 <- sprintf("Price RMSE: rB = %.4f vs BS = %.4f (%s).",
                    pooled$Price_RMSE_rB, pooled$Price_RMSE_BS,
                    if(!is.null(t_pr)) sprintf("p=%.3g", t_pr$p.value) else "p=N/A")
  
  html_body <- paste0(
    "<!doctype html><meta charset='utf-8'><title>OOS rB vs BS Comparison</title>",
    "<style>body{font:14px sans-serif; margin:20px;} h1{margin-bottom:0.5em;} ",
    ".note{color:#555;} table{border-collapse:collapse; width:100%; margin:1em 0;} ",
    "th, td{border:1px solid #ccc; padding:4px 8px;} th{background:#f7f7f7;} ",
    "th.left, td.left{text-align:left;} th.right, td.right{text-align:right;} ",
    ".badge{display:inline-block; padding:3px 7px; border-radius:4px; background:#def; margin-right:5px;} </style>",
    sprintf("<h1>Out-of-sample Comparison ??? rB vs BS (calibrated on %s)</h1>", calib_date),
    sprintf("<p class='note'><span class='badge'>Conclusion:</span>%s %s %s</p>", concl, concl2, concl3),
    "<h2>Pooled Performance</h2>", score_html,
    "<h2>By Maturity Bucket (median ??IV RMSE per bucket)</h2>", byT_html,
    "<h2>Per-Date Performance</h2>", daily_html
  )
  writeLines(html_body, outfile)
  message("OOS comparison complete. Results written to ", normalizePath(outfile))
  
  # Optional: persist the cache to disk
  if (!is.null(cache_path)) save_rb_cache(cache_path)
  
  invisible(list(
    pooled = pooled, daily = daily_all, by_T = byT_df,
    tests = list(t_iv=t_iv, t_pr=t_pr, sign_test_p=p_binom, wins=n_wins, n_days=n_days),
    params = c(H_star=H_star, Rho_star=Rho_star, Eta_star=Eta_star)
  ))
}

# --- Run (adjust N/M as you like) ---
results <- oos_compare_many(
  date_keys = date_keys_oos, spxIvolList = spxIvolList,
  r = r, q = q,
  H_star = HCal, Rho_star = RhoCal, Eta_star = EtaCal,
  N = 500, M = 10000, seed_base = 1,
  outfile = "oos_multi_summary.html",
  calib_date = date_key,
  cache_path = NULL  # set e.g. "rb_cache.rds" to persist across sessions
)
