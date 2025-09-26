# --- 0) Packages --------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

#read and clean the file for everything we don't need

csv_name <- read.csv("sample_data.csv")
spx_clean <- csv_name %>%
  mutate(
    type = case_when(
      str_to_upper(option_type) %in% c("C", "CALL") ~ "C",
      str_to_upper(option_type) %in% c("P", "PUT") ~ "P",
      TRUE ~ NA_character_
    ),
    bid_use = coalesce(bid_1545, bid_eod),
    ask_use = coalesce(ask_1545, ask_eod),
    mid_use = (bid_use + ask_use)/2
  ) %>%
  filter(
    underlying_symbol %in% c("^SPX", "SPX"),
    root %in% c("SPX", "SPXW"),
    type %in% c("C", "P"),
    is.finite(strike), strike > 0
  ) %>%
  transmute(
    quote_date = as.Date(quote_date),
    expiration = as.Date(expiration),
    root,
    strike     = as.numeric(strike),
    type,
    bid_1545   = as.numeric(bid_1545),
    ask_1545   = as.numeric(ask_1545),
    bid_eod    = as.numeric(bid_eod),
    ask_eod    = as.numeric(ask_eod),
    bid        = as.numeric(bid_use),
    ask        = as.numeric(ask_use),
    mid        = as.numeric(mid_use)
  ) %>%
  arrange(expiration, strike, type) %>%
  distinct()

writexl::write_xlsx(spx_clean, "spx_clean.xlsx")

#small check to see if you have SPX and SPXW and check for amount of strikes and rows
spx_clean %>% count(root)
spx_clean %>% group_by(expiration) %>% summarise(n_rows=n(), n_strikes=n_distinct(strike)) %>% arrange(expiration) %>% print(n=10)

#This is a demo for the forward variance curve
#We are using a flat risk free curve, when doing this again later this need to be a curve
#We are working on a 365 days a year basis

r_flat <- 0.05
year_basis <- 365

#creating a mid_ function, this takes the average
mid_ <- function(bid, ask) (bid + ask)/2

deltaK_ <- function(K) {
  K <- sort(K)
  n <- length(K)
  dK <- numeric(n)
  if (n >= 2) {
    dK[1] <- K[2] - K[1]
    dK[n] <- K[n] - K[n-1]
  }
  if (n >= 3) dK[2:(n-1)] <- (K[3:n] - K[1:(n-2)])/2
  dK
}

infer_forward_from_parity <- function(wide, D_T) {
  ok <- is.finite(wide$mid_C) & is.finite(wide$mid_P) & wide$mid_C > 0 & wide$mid_P > 0
  stopifnot(any(ok))
  F_cands <- wide$K[ok] + (wide$mid_C[ok] - wide$mid_P[ok]) / D_T
  stats::median(F_cands, na.rm = TRUE)
}

trim_two_zerobids <- function(df, side = c("below", "above")) {
  side <- match.arg(side)
  if (side == "below") {
    df_scan <- df %>% arrange(desc(K))
  } else {
    df_scan <- df %>% arrange(K)
  }
  z <- (df_scan$bid_side <= 0) | !is.finite(df_scan$bid_side)
  cut_idx <- which(z & dplyr::lag(z, default = FALSE))
  if (length(cut_idx)) df_scan <- df_scan[seq_len(min(cut_idx) - 1), , drop = FALSE]
  df_scan %>% arrange(K)
}


varswap_one_expiry <- function(chain_one_expiry, valuation_date, exp_dt, year_basis = 365, r_flat = 0.0) {
  #define the expiration date as a date 
  #then check that there is only one expiration date
  exp_dt = as.Date(exp_dt)
  stopifnot(length(exp_dt) == 1)
  
  T_days <- as.numeric(difftime(exp_dt, as.Date(valuation_date), units = "days"))
  T_yrs <- T_days/year_basis
  if (T_yrs <= 0) return(tibble())
  
  D_T <- exp(-r_flat * T_yrs)
  
  wide <- chain_one_expiry %>%
    transmute(K = as.numeric(strike),
              type,
              bid = as.numeric(bid),
              ask = as.numeric(ask),
              mid = as.numeric(mid)) %>%
    group_by(K, type) %>%
    summarise(bid = mean(bid, na.rm = TRUE),
              mid = mean(mid, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = c(bid, mid), names_sep = "_") %>%
    arrange(K)
  
  for (nm in c("mid_C", "mid_P", "bid_C", "bid_P")) if (!nm %in% names(wide)) wide[[nm]] <- NA_real_
  
  F_T <- infer_forward_from_parity(wide, D_T)
  
  K0 <- max(wide$K[wide$K <= F_T])
  if (!is.finite(K0)) return(tibble())
  
  otm <- wide %>%
    transmute(
      K,
      Q = dplyr::case_when(
        K < K0 ~ mid_P,
        K > K0 ~ mid_C,
        TRUE   ~ 0.5 * (mid_C + mid_P)
      ),
      bid_side = dplyr::case_when(
        K < K0 ~ bid_P,
        K > K0 ~ bid_C,
        TRUE   ~ 0.5 * (bid_C + bid_P)
      )
    ) %>%
    filter(is.finite(Q))
  
  below <- otm %>% filter(K <= K0) %>% trim_two_zerobids("below")
  above <- otm %>% filter(K >= K0) %>% trim_two_zerobids("above")
  use   <- bind_rows(below, above) %>% distinct(K, .keep_all = TRUE) %>% arrange(K)
  if (nrow(use) < 3) return(tibble())
  
  dK <- deltaK_(use$K)
  
  term_sum <- sum( (dK/(use$K^2)) * (use$Q/D_T), na.rm = TRUE )
  sigma2   <- (2/T_yrs) * term_sum - (1/T_yrs) * ((F_T/K0 - 1)^2)
  sigma2   <- max(sigma2, 0)
  
  tibble(
    expiry    = exp_dt,
    T_years   = T_yrs,
    discount  = D_T,
    forward   = F_T,
    K0        = K0,
    n_strikes = nrow(use),
    K_var     = sigma2,
    TV        = T_yrs * sigma2
  )
}

#pick date
valuation_date <- spx_clean %>% summarise(d = first(quote_date)) %>% pull(d)


varswap_pts <- spx_clean %>%
  group_by(expiration) %>%
  group_modify(~ varswap_one_expiry(
    .x,
    valuation_date,
    exp_dt          = .y$expiration,
    year_basis = year_basis,
    r_flat = r_flat
    )) %>%
  ungroup() %>%
  arrange(T_years)



# 3) Interpolate in *total variance* and compute piecewise-constant forwards
tv <- varswap_pts %>% select(T = T_years, TV) %>% arrange(T)
# ensure non-decreasing TV
tv <- tv %>% mutate(TV = cummax(TV))

fwd_tbl <- tv %>%
  mutate(T_prev  = lag(T, default = 0),
         TV_prev = lag(TV, default = 0)) %>%
  filter(T > 0) %>%
  transmute(T_start = T_prev,
            T_end   = T,
            fwd_var = (TV - TV_prev) / (T - T_prev))

# 4) Sample xi0 on any simulation grid (here: daily trading grid to 5y)
t_grid <- seq(0, 5, by = 1/252)  # years
xi0_on_grid <- function(fwd_tbl, t_grid) {
  # piecewise-constant lookup
  sapply(t_grid, function(ti) {
    i <- max(which(fwd_tbl$T_start <= ti & ti < fwd_tbl$T_end), na.rm = TRUE)
    if (!is.finite(i) || i == -Inf) fwd_tbl$fwd_var[1] else fwd_tbl$fwd_var[i]
  })
}
xi_vec <- xi0_on_grid(fwd_tbl, t_grid)

# Optional: mass-correct per interval so the discrete sum matches TV exactly
fwd_tbl <- fwd_tbl %>%
  rowwise() %>%
  mutate(
    # indices in grid within this interval
    idx = list(which(t_grid >= T_start & t_grid < T_end)),
    # scale so sum(xi*dt) equals area
    scale = (fwd_var * (T_end - T_start)) / (sum(rep(1/252, length(idx))) * fwd_var + 1e-16)
  ) %>%
  ungroup()

# Apply scaling
for (r in seq_len(nrow(fwd_tbl))) {
  idx <- fwd_tbl$idx[[r]]
  if (length(idx)) xi_vec[idx] <- xi_vec[idx] * fwd_tbl$scale[r]
}

# 5) Quick plots (optional)
ggplot(tv, aes(T, TV)) + geom_line() +
  labs(title = "Total implied variance TV(T)", x = "T (years)", y = "TV")

fv_plot <- tibble( 
    Tp = c(rbind(fwd_tbl$T_start, fwd_tbl$T_end)),
    FVp = c(rbind(fwd_tbl$fwd_var, fwd_tbl$fwd_var))
  )     
ggplot(fv_plot, aes(Tp, FVp)) + geom_step(direction = "hv") +
  labs(title = "Forward variance (piecewise-constant)", x = "T (years)", y = "FV")

