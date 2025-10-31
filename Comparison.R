# ====================== MULTI-DATE OOS COMPARISON (Tidy Tables) ======================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(digest)
})

# --- OOS date selection (assumes 'date_key' from calibration step exists) ---
all_dates     <- sort(names(spxIvolList))
date_keys_oos <- all_dates[all_dates > date_key]

# -------- sanity: required core objects/functions --------
stopifnot(exists("spxIvolList"),
          exists("price_many_rB_turbo"),
          exists("build_calibration_df"),
          exists("bs_implied_vol"))

# -------- helpers: metrics, buckets --------
RMSE <- function(x) { y <- x[is.finite(x)]; if (!length(y)) NA_real_ else sqrt(mean(y^2)) }
MAE  <- function(x) { y <- x[is.finite(x)]; if (!length(y)) NA_real_ else mean(abs(y)) }

T_bucket <- function(T) {
  cut(T,
      breaks = c(0, 1/12, 0.25, 0.5, 1, 2, 5, Inf),
      labels = c("???1M","1???3M","3???6M","6???12M","1???2Y","2???5Y",">5Y"),
      right = TRUE, include.lowest = TRUE)
}

# -------- optional: safe BS IV wrapper (uses your bs_implied_vol) --------
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
# Same logic you used earlier, factored as a helper.
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
  T_min <- suppressWarnings(min(chain$T[is.finite(chain$T)], na.rm=TRUE))
  near  <- chain %>% dplyr::filter(abs(T - T_min) < 1e-6 & is.finite(F))
  if (nrow(near) == 0) near <- chain %>% dplyr::filter(is.finite(F))
  median(near$F * exp(-(r - q)*near$T), na.rm=TRUE)
}

global_T_max <- function(keys, spxIvolList) {
  suppressWarnings(max(vapply(keys, function(d) {
    ch <- spxIvolList[[d]]
    if (!is.data.frame(ch) || !nrow(ch)) return(NA_real_)
    max(as.numeric(ch$Texp), na.rm=TRUE)
  }, numeric(1)), na.rm=TRUE))
}

# Build a stable key that includes params + xi0 digest
.cache_key <- function(d, H, Rho, Eta, N, M, seed, S0, r, q, xi0_df) {
  paste(
    d, round(H,6), round(Rho,6), round(Eta,6), N, M, seed,
    round(S0,6), r, q, digest(xi0_df), sep="|"
  )
}

# Pre-price whole date d (vectorized) and store
preprice_rb_for_date <- function(d, chain_raw, r, q,
                                 H, Rho, Eta,
                                 N=500L, M=10000L,
                                 seed_base=1L,
                                 extend_to=NULL,
                                 antithetic=TRUE,
                                 ncores=4) {
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
  
  seed_i <- as.integer(abs((seed_base + sum(utf8ToInt(d))*131) %% .Machine$integer.max))
  rB_tbl <- price_many_rB_turbo(
    options_df = oos[, c("type","K","T")],
    S0 = S0_today, r = r, q = q,
    alpha = H - 0.5, eta = Eta, rho = Rho,
    xi0_curve = xi0_curve_today,
    N = N, M = M, antithetic = antithetic, seed = seed_i,
    ncores = ncores, clip_iv = TRUE
  ) %>%
    mutate(type = tolower(trimws(type)),
           K = as.numeric(K), T = T_round(as.numeric(T)),
           price_rB = price, .keep = "unused")
  
  meta <- list(
    date=d, H=H, Rho=Rho, Eta=Eta, N=N, M=M, seed=seed_i,
    S0=S0_today, r=r, q=q, xi0_digest = digest(xi0_curve_today),
    antithetic=antithetic, ncores=ncores
  )
  
  key <- .cache_key(d, H, Rho, Eta, N, M, seed_i, S0_today, r, q, xi0_curve_today)
  assign(key, list(prices = rB_tbl, meta = meta), envir = .rB_cache)
  invisible(list(prices=rB_tbl, meta=meta, key=key))
}

# Retrieve cached prices; recompute if missing or params differ
get_rb_prices <- function(d, wanted_df, chain_raw, r, q,
                          H, Rho, Eta,
                          N=500L, M=10000L, seed_base=1L,
                          extend_to=NULL, antithetic=TRUE, ncores=1) {
  keys <- ls(.rB_cache, all.names = TRUE)
  cand <- grep(paste0("^", d, "\\|"), keys, value = TRUE)
  entry <- if (length(cand)) get(tail(cand,1), envir = .rB_cache) else NULL
  
  needs_rebuild <- is.null(entry) ||
    is.null(entry$meta) ||
    any(abs(c(entry$meta$H, entry$meta$Rho, entry$meta$Eta) - c(H,Rho,Eta)) > 1e-12) ||
    entry$meta$N != N || entry$meta$M != M || entry$meta$r != r || entry$meta$q != q ||
    isTRUE(entry$meta$antithetic) != isTRUE(antithetic) || entry$meta$ncores != ncores
  
  if (needs_rebuild) {
    pre <- preprice_rb_for_date(
      d, chain_raw, r, q, H, Rho, Eta, N, M, seed_base, extend_to, antithetic, ncores
    )
    if (is.null(pre)) return(NULL)
    entry <- get(pre$key, envir = .rB_cache)
  }
  
  wanted_df <- wanted_df %>% mutate(type=tolower(trimws(type)), K=as.numeric(K), T=T_round(T))
  out <- wanted_df %>% left_join(entry$prices, by=c("type","K","T"))
  
  if (anyNA(out$price_rB)) {
    # Price only the missing ones, reuse same meta seed/params
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
    missing <- out %>% filter(!is.finite(price_rB)) %>% distinct(type,K,T)
    extra <- price_many_rB_turbo(
      options_df = missing,
      S0 = S0_today, r = r, q = q,
      alpha = H - 0.5, eta = Eta, rho = Rho,
      xi0_curve = xi0_curve_today,
      N = entry$meta$N, M = entry$meta$M,
      antithetic = entry$meta$antithetic, seed = entry$meta$seed,
      ncores = entry$meta$ncores, clip_iv = TRUE
    ) %>%
      transmute(type=tolower(trimws(type)), K=as.numeric(K), T=T_round(as.numeric(T)),
                price_rB = price)
    updated <- bind_rows(entry$prices, extra) %>% distinct(type,K,T,.keep_all=TRUE)
    entry$prices <- updated
    assign(cand[length(cand)], entry, envir=.rB_cache)
    out <- wanted_df %>% left_join(updated, by=c("type","K","T"))
  }
  out
}

# ============================ MAIN OOS COMPARISON =============================
oos_compare_many_tbl <- function(date_keys, spxIvolList, r, q,
                                 H_star, Rho_star, Eta_star,
                                 N = 500L, M = 10000L,
                                 seed_base = 1L,
                                 antithetic = TRUE,
                                 ncores = 1,
                                 calib_date = get0("date_key", ifnotfound="(unknown)")
) {
  # Global T_max to extend xi0 per date
  Tmax_global <- global_T_max(date_keys, spxIvolList)
  if (!is.finite(Tmax_global)) Tmax_global <- NA_real_
  
  per_date <- function(d) {
    chain <- spxIvolList[[d]]
    if (!is.data.frame(chain) || !nrow(chain)) return(NULL)
    
    chain <- chain %>%
      mutate(T = as.numeric(Texp),
             K = as.numeric(Strike),
             F = as.numeric(Fwd),
             CallMid = suppressWarnings(as.numeric(CallMid)))
    
    # OOS option table (calls + puts via parity)
    df_oos <- chain %>%
      filter(is.finite(T), T > 0, is.finite(K), is.finite(CallMid), is.finite(F)) %>%
      mutate(call_px = CallMid,
             put_px  = call_px - exp(-r * T) * (F - K))
    if (!nrow(df_oos)) return(NULL)
    
    options_df_oos <- bind_rows(
      df_oos %>% transmute(type="call", K, T=T_round(T), mkt_price=call_px, date=d),
      df_oos %>% transmute(type="put",  K, T=T_round(T), mkt_price=put_px,  date=d)
    )
    
    # rB prices via cache
    rB_prices <- get_rb_prices(
      d = d, wanted_df = options_df_oos[, c("type","K","T")],
      chain_raw = chain, r = r, q = q,
      H = H_star, Rho = Rho_star, Eta = Eta_star,
      N = N, M = M, seed_base = seed_base,
      extend_to = if (is.finite(Tmax_global)) Tmax_global else max(options_df_oos$T, na.rm=TRUE),
      antithetic = antithetic, ncores = ncores
    )
    if (is.null(rB_prices)) return(NULL)
    
    # BS baseline: use shortest-bucket ATM vol (from xi0)
    S0_today <- build_S0_from_chain(chain, r, q)
    calib_df_today <- build_calibration_df(
      spxIvolList, date_key = d,
      grid_list = lapply(unique(chain$Expiry),
                         function(e) list(expiration=e, ratios=c(0.8,0.9,1.0,1.1,1.2))),
      S0 = S0_today, r = r, q = q, use_chain_fwd = TRUE
    )
    xi0_curve_today <- build_xi0_from_ATM(
      calib_df_today,
      extend_to = if (is.finite(Tmax_global)) Tmax_global else max(options_df_oos$T, na.rm=TRUE)
    )
    Sigma_BS <- sqrt(max(xi0_curve_today$forward_var[1], 1e-12))
    
    # Market IV + errors
    mkt_iv <- safe_bs_iv_vec(options_df_oos$type, options_df_oos$mkt_price,
                             S=S0_today, Ks=options_df_oos$K, Ts=options_df_oos$T, r=r, q=q)
    
    cmp <- options_df_oos %>%
      mutate(mkt_iv = mkt_iv,
             type   = tolower(trimws(type)),
             K      = as.numeric(K),
             T      = T_round(T)) %>%
      left_join(rB_prices, by=c("type","K","T")) %>%
      mutate(
        # BS price using your bs_price (loaded from core file)
        price_BS = mapply(function(tp, K_, T_) {
          bs_price(S=S0_today, K=K_, T=T_, r=r, q=q, sigma=Sigma_BS, type=tp)
        }, type, K, T),
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
      Price_MSE_BS   = mean(errP_BS[is.finite(errP_BS)]^2)
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
  
  # By maturity bucket (median ??RMSE per bucket)
  byT_df <- if (nrow(cmp_all)) {
    cmp_all %>% mutate(T_bucket = T_bucket(T)) %>%
      group_by(date, T_bucket) %>%
      summarise(RMSE_rB = RMSE(errIV_rB), RMSE_BS = RMSE(errIV_BS), .groups="drop") %>%
      mutate(dIV_RMSE = RMSE_rB - RMSE_BS) %>%
      group_by(T_bucket) %>%
      summarise(
        n_dates = n(),
        median_dIV_RMSE = median(dIV_RMSE, na.rm=TRUE),
        mean_dIV_RMSE   = mean(dIV_RMSE, na.rm=TRUE),
        winrate_days    = mean(dIV_RMSE < 0, na.rm=TRUE),
        .groups="drop"
      )
  } else {
    tibble(T_bucket=factor(levels=T_bucket(1)), n_dates=0, median_dIV_RMSE=NA_real_,
           mean_dIV_RMSE=NA_real_, winrate_days=NA_real_)
  }
  
  # Paired t-tests across days (IV/Price MSE)
  iv_pair_idx <- with(daily_all, is.finite(IV_MSE_rB) & is.finite(IV_MSE_BS))
  pr_pair_idx <- with(daily_all, is.finite(Price_MSE_rB) & is.finite(Price_MSE_BS))
  t_iv <- if (sum(iv_pair_idx) >= 2) t.test(daily_all$IV_MSE_rB[iv_pair_idx],
                                            daily_all$IV_MSE_BS[iv_pair_idx],
                                            paired=TRUE, alternative="less") else NULL
  t_pr <- if (sum(pr_pair_idx) >= 2) t.test(daily_all$Price_MSE_rB[pr_pair_idx],
                                            daily_all$Price_MSE_BS[pr_pair_idx],
                                            paired=TRUE, alternative="less") else NULL
  
  valid_days <- daily_all %>% mutate(dIV_RMSE = IV_RMSE_rB - IV_RMSE_BS) %>%
    filter(is.finite(dIV_RMSE)) %>% pull(dIV_RMSE)
  n_days <- length(valid_days); n_wins <- sum(valid_days < 0)
  p_binom <- if (n_days > 0) pbinom(q = n_wins - 1, size = n_days, prob = 0.5, lower.tail = FALSE) else NA_real_
  
  # Nicely formatted tables for RStudio printing
  pooled_fmt <- pooled %>%
    mutate(across(where(is.numeric), ~round(., 6))) %>%
    mutate(IV_impr_pct = round(IV_impr_pct, 2),
           Price_impr_pct = round(Price_impr_pct, 2))
  
  daily_fmt <- daily_all %>%
    mutate(across(where(is.numeric), ~round(., 6)))
  
  byT_fmt <- byT_df %>%
    mutate(across(where(is.numeric), ~round(., 6)))
  
  # Return tidy tables (no HTML)
  invisible(list(
    pooled = pooled_fmt,
    daily  = daily_fmt,
    by_T   = byT_fmt,
    params = c(H_star=H_star, Rho_star=Rho_star, Eta_star=Eta_star,
               N=N, M=M, antithetic=antithetic, ncores=ncores,
               calib_date=calib_date),
    tests  = list(t_iv=t_iv, t_pr=t_pr, sign_test=list(n_days=n_days, n_wins=n_wins, p=p_binom))
  ))
}

# --- Example run (adjust N/M and cores as you like) ---
results <- oos_compare_many_tbl(
  date_keys = date_keys_oos, spxIvolList = spxIvolList,
  r = r, q = q,
  H_star = HCal, Rho_star = RhoCal, Eta_star = EtaCal,
  N = 500, M = 2000,       # reasonable calibration-size grid for speed
  seed_base = 1,
  antithetic = TRUE,
  ncores = 4,              # propagate to rB pricer
  calib_date = date_key
)

# ---- One-table RMSE/MAE comparison ----
rmse_mae_table <- function(results, digits = 6, perc_digits = 2) {
  stopifnot(is.list(results), !is.null(results$pooled))
  p <- results$pooled
  stopifnot(all(c("IV_RMSE_rB","IV_RMSE_BS","Price_RMSE_rB","Price_RMSE_BS",
                  "IV_MAE_rB","IV_MAE_BS","Price_MAE_rB","Price_MAE_BS") %in% names(p)))
  out <- tibble::tibble(
    Metric = c("IV RMSE","Price RMSE","IV MAE","Price MAE"),
    rB     = c(p$IV_RMSE_rB,   p$Price_RMSE_rB,   p$IV_MAE_rB,   p$Price_MAE_rB),
    BS     = c(p$IV_RMSE_BS,   p$Price_RMSE_BS,   p$IV_MAE_BS,   p$Price_MAE_BS)
  ) |>
    dplyr::mutate(
      `?? (rB - BS)`        = rB - BS,
      `Improvement %`      = 100 * (1 - rB / BS)
    ) |>
    dplyr::mutate(
      dplyr::across(c(rB, BS, `?? (rB - BS)`), ~round(., digits)),
      `Improvement %` = round(`Improvement %`, perc_digits)
    )
  out
}

# ---- Print nicely (easy to copy from RStudio) ----
tab <- rmse_mae_table(results, digits = 6, perc_digits = 2)
print(tab, n = Inf)      # copy from the Console straight into Word


# Inspect in RStudio:
print(results$pooled)
print(results$by_T)
print(results$daily)



