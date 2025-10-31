#---------------- Load and clean the dataset ----------------
library(readxl)
library(janitor)
library(readr)
library(lubridate)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(nloptr)
library(numDeriv)

# ============ Core inputs ============
# You can override S0/r/q below. If use_chain_S0=TRUE, we infer S0 from the chain's Fwd & Texp.
r  <- 0.0
q  <- 0.0
day_count <- 365

# Monte Carlo sizes (pricing vs. calibration vs. hedging)
M    <- 200000;  N    <- 500     # for main pricing
MCal <- 10000;   NCal <- 500      # for calibration
MHedge <- 10000; NHedge <- 5000    # for hedging simulation

# ============ Load market data ============
# Expecting this .RData to provide `spxIvolList` as a named list of per-date data frames.
load("spxOptionMetricsIVols.rdata")
stopifnot(exists("spxIvolList"), is.list(spxIvolList), length(spxIvolList) > 0)

# (Optional) Read VIX term structure; we derive xi0 from ATM IVs anyway.
VIX30 <- tryCatch(read_xlsx("VIX 30.xlsx"), error = function(e) NULL)

# ----- Choose a trading date key from the panel -----
# Set this explicitly, e.g. date_key <- "19960104". If NULL, we'll take the first available.
date_key <- "20140805"
if (is.null(date_key)) date_key <- names(spxIvolList)[1]
stopifnot(date_key %in% names(spxIvolList))
chain <- spxIvolList[[date_key]]
stopifnot(is.data.frame(chain), nrow(chain) > 0)

# ============ Helper definitions (panel-friendly) ============

.parse_date_loose <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(as.character(x), format = "%Y%m%d"))
  try_formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d.%m.%Y", "%Y%m%d")
  for (fmt in try_formats) {
    d <- suppressWarnings(as.Date(x, format = fmt))
    if (!is.na(d)) return(d)
  }
  as.Date(NA)
}

# Build a simple ratio grid per expiry if you don't supply one:
make_grid_from_chain <- function(chain_df, ratios = c(0.8, 0.9, 1.0, 1.1, 1.2)) {
  exps <- unique(chain_df$Expiry)
  lapply(exps, function(e) list(expiration = e, ratios = ratios))
}

# Safe IV inversion using the BS helper in your BS file
.safe_bs_iv <- function(opt_type, S0, K, T, r, q, price) {
  out <- try(bs_implied_vol(mkt_price = price, S = S0, K = K, T = T, r = r, q = q, type = opt_type),
             silent = TRUE)
  if (inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}

# Build calibration table from one date in spxIvolList
build_calibration_df <- function(spxIvolList, date_key, grid_list,
                                 S0, r, q,
                                 use_chain_fwd = TRUE,
                                 call_col = "CallMid", put_col = "PutMid",
                                 date_tol_days = 7L) {
  stopifnot(date_key %in% names(spxIvolList))
  chain <- spxIvolList[[date_key]]
  stopifnot(all(c("Expiry","Texp","Strike", call_col) %in% names(chain)))
  chain$expiry_date <- .parse_date_loose(chain$Expiry)
  if (anyNA(chain$expiry_date)) stop("Could not parse some Expiry values into Dates.")
  chain$Strike <- as.numeric(chain$Strike)
  chain$Texp   <- as.numeric(chain$Texp)
  
  split_by_exp <- split(chain, chain$expiry_date)
  avail_exps   <- sort(unique(chain$expiry_date))
  
  out <- list(); k <- 0L
  
  for (entry in grid_list) {
    target_exp <- .parse_date_loose(entry$expiration)
    if (is.na(target_exp)) next
    # nearest match within tolerance
    diffs <- abs(as.integer(avail_exps - target_exp))
    i     <- if (length(diffs)) which.min(diffs) else NA_integer_
    if (!length(i) || is.na(i) || diffs[i] > date_tol_days) next
    
    exp_date <- avail_exps[i]
    sub_df   <- split_by_exp[[as.character(exp_date)]]
    if (is.null(sub_df) || !nrow(sub_df)) next
    
    T_val <- median(sub_df$Texp, na.rm = TRUE)
    disc_r <- exp(-r * T_val)
    
    # Forward for this expiry
    if (use_chain_fwd && ("Fwd" %in% names(sub_df)) && any(is.finite(sub_df$Fwd))) {
      Fwd_e <- median(sub_df$Fwd[is.finite(sub_df$Fwd)], na.rm = TRUE)
    } else {
      Fwd_e <- S0 * exp((r - q) * T_val)
    }
    
    # Prepare puts (via parity if PutMid missing)
    have_put_col <- (put_col %in% names(sub_df))
    if (!have_put_col) {
      rhs <- function(C, K) C - (Fwd_e - K) * disc_r                      # P = C - (F-K)e^{-rT}
      sub_df$PutMid_from_parity <- rhs(sub_df[[call_col]], sub_df$Strike)
    }
    
    for (ratio in entry$ratios) {
      target_K <- Fwd_e * ratio
      j <- which.min(abs(sub_df$Strike - target_K))
      if (!length(j) || !is.finite(j)) next
      
      K_j <- sub_df$Strike[j]
      opt_type <- if (ratio >= 1) "call" else "put"
      mid_price <- if (opt_type == "call") {
        sub_df[[call_col]][j]
      } else {
        if (have_put_col) sub_df[[put_col]][j] else sub_df$PutMid_from_parity[j]
      }
      if (!is.finite(mid_price) || mid_price <= 0) next
      
      iv <- .safe_bs_iv(opt_type, S0, K_j, T_val, r, q, mid_price)
      if (!is.finite(iv) || iv <= 0) next
      
      k <- k + 1L
      out[[k]] <- data.frame(
        date_key    = date_key,
        expiry_date = exp_date,
        type        = opt_type,
        K           = K_j,
        T           = T_val,
        fwd         = Fwd_e,
        ratio       = ratio,
        mid_price   = mid_price,
        market_iv   = iv,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (k == 0L) {
    warning("No calibration rows built; check grid_list vs. chain expiries.")
    return(data.frame())
  }
  
  df <- do.call(rbind, out)
  df$type <- tolower(trimws(df$type))
  df <- df[is.finite(df$market_iv) & is.finite(df$K) & is.finite(df$T) & df$T > 0, , drop = FALSE]
  df <- df[order(df$T, match(df$type, c("call","put")), df$K), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# ============ Prepare inputs from the chosen panel date ============
# Optionally infer S0 from chain if forwards are provided (recommended for 1996???2014 panel)
use_chain_S0 <- TRUE
if (use_chain_S0 && "Fwd" %in% names(chain)) {
  # S0 ??? median(Fwd * e^{-(r-q)T}) across the smallest Texp bucket
  chain <- chain %>% mutate(Texp = as.numeric(Texp))
  T_min <- min(chain$Texp[is.finite(chain$Texp)], na.rm = TRUE)
  near  <- chain %>% filter(is.finite(Fwd), is.finite(Texp), abs(Texp - T_min) <= 1e-6)
  if (nrow(near) == 0) near <- chain %>% filter(is.finite(Fwd), is.finite(Texp))
  if (nrow(near) > 0) {
    S0 <- median(near$Fwd * exp(-(r - q) * near$Texp), na.rm = TRUE)
  }
}

# Create a default strike-ratio grid per expiry
grid_list <- make_grid_from_chain(chain, ratios = c(0.8, 0.9, 1.0, 1.1, 1.2))

# ============ Build calibration DataFrame ============
options_df_calib <- build_calibration_df(
  spxIvolList = spxIvolList,
  date_key    = date_key,
  grid_list   = grid_list,
  S0 = S0, r = r, q = q,
  use_chain_fwd = TRUE
)
stopifnot(is.data.frame(options_df_calib), nrow(options_df_calib) > 0)

# ============ Build full options table (calls and puts) ============
# Use every (Strike, Texp) pair in the chain; duplicate into call and put rows.
stopifnot(all(c("Strike","Texp") %in% names(chain)))
K_vec <- as.numeric(chain$Strike)
T_vec <- as.numeric(chain$Texp)
# Stabilize tiny float noise in maturities
T_vec <- round(T_vec, 10)

options_df <- rbind(
  data.frame(type = "call", K = K_vec, T = T_vec),
  data.frame(type = "put",  K = K_vec, T = T_vec)
)
options_df$type <- tolower(trimws(options_df$type))
options_df$K    <- as.numeric(options_df$K)
options_df$T    <- as.numeric(options_df$T)
options_df <- options_df[is.finite(options_df$K) & is.finite(options_df$T) & options_df$T > 0, , drop = FALSE]

# ============ Compute forward variance curve from ATM IVs ============
# Build xi0 using the full chain (one true-ATM per maturity from *all* strikes)
has_put <- "PutMid" %in% names(chain)

atm_by_T_full <- chain %>%
  transmute(
    T = as.numeric(Texp),
    K = as.numeric(Strike),
    F_raw   = as.numeric(Fwd),
    CallMid = as.numeric(CallMid),
    PutMid  = if (has_put) as.numeric(.data$PutMid) else NA_real_
  ) %>%
  mutate(F = ifelse(is.finite(F_raw), F_raw, S0 * exp((r - q) * T))) %>%
  select(-F_raw) %>%
  filter(is.finite(T), T > 0, is.finite(K), is.finite(F), is.finite(CallMid)) %>%
  group_by(T) %>%
  slice_min(order_by = abs(K - median(F, na.rm = TRUE)), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(disc_r = exp(-r * T),
         P_mid  = ifelse(is.finite(PutMid), PutMid, CallMid - disc_r * (F - K))) %>%
  rowwise() %>%  # <-- make bs_implied_vol scalar-safe
  mutate(
    market_iv = bs_implied_vol(
      mkt_price = CallMid, S = S0, K = K, T = T, r = r, q = q, type = "call"
    )
  ) %>%
  ungroup() %>%
  transmute(T, K, market_iv, fwd = F)

xi0_curve <- build_xi0_from_ATM(atm_by_T_full, 
            extend_to = max(options_df$T, 
            na.rm = TRUE))



# Extend flat beyond last calibration maturity if pricing longer-dated options
T_max_calib <- max(xi0_curve$T_end)
T_max       <- max(options_df$T)
if (T_max > T_max_calib) {
  xi0_curve <- dplyr::bind_rows(
    xi0_curve,
    data.frame(T_start = T_max_calib, T_end = T_max, forward_var = tail(xi0_curve$forward_var, 1))
  )
}

# Initial instantaneous forward variance & BS baseline
xi0_level <- xi0_curve$forward_var[1]
Sigma     <- sqrt(xi0_level)