# ============================================================
# ROUGH BERGOMI KALIBRERING (SLSQP) ??? EN DAG (date_key)
#   - Grid G1 (30 punkter), OTM-only
#   - Kost: sum_{i=1..30} (IV_model - IV_mkt)^2
#   - Parametre: (H, rho, eta) med SLSQP og boksgrenser
# Forutsetter at xi0_curve, price_many_rB_turbo, bs_implied_vol er lastet.
# ============================================================

# ---------- Avhengigheter ----------
suppressPackageStartupMessages({
  library(dplyr)
  library(nloptr)
  library(numDeriv)
})

stopifnot(exists("spxIvolList"), exists("date_key"), exists("xi0_curve"))
options_raw <- spxIvolList[[date_key]]
stopifnot(is.data.frame(options_raw), nrow(options_raw) > 0)
stopifnot(all(c("Texp","Strike","Fwd") %in% names(options_raw)))
stopifnot("CallMid" %in% names(options_raw))   # vi bruker CallMid + paritet for puts

# ---------- Basis inputs (bruker eksisterende om de finnes) ----------
r     <- if (exists("r")) get("r") else 0.0
q     <- if (exists("q")) get("q") else 0.0
MCal  <- if (exists("MCal")) get("MCal") else 1000L
NCal  <- if (exists("NCal")) get("NCal") else 500L
seed  <- if (exists("seed")) get("seed") else 1L

# ---------- Ankre S0 til forwards p?? korteste T ----------
options_raw$Texp   <- as.numeric(options_raw$Texp)
options_raw$Strike <- as.numeric(options_raw$Strike)
T_min <- min(options_raw$Texp[is.finite(options_raw$Texp)], na.rm = TRUE)
nearT <- subset(options_raw, is.finite(Fwd) & is.finite(Texp) & abs(Texp - T_min) <= 1e-6)
if (nrow(nearT) == 0) nearT <- subset(options_raw, is.finite(Fwd) & is.finite(Texp))
S0 <- median(nearT$Fwd * exp(-(r - q) * nearT$Texp), na.rm = TRUE)

# ---------- Hjelpere ----------
T_round <- function(x) round(as.numeric(x), 10)
K_round <- function(x) as.numeric(x)

# Bygg G1-punkter: 6 l??p. med ratio-settene som gir totalt 30 punkter
.build_G1_df <- function(df) {
  df <- df %>%
    mutate(T = as.numeric(Texp), K = as.numeric(Strike), F = as.numeric(Fwd))
  target_Ts <- c("1M"=1/12, "3M"=3/12, "6M"=0.5, "1Y"=1, "2Y"=2, "5Y"=5)
  ratio_map <- list(
    "1M" = c(1.00),
    "3M" = c(0.90, 1.00),
    "6M" = c(0.80, 0.90, 1.00, 1.10, 1.20),
    "1Y" = c(0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40),
    "2Y" = c(0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40),
    "5Y" = c(0.80, 0.90, 1.00, 1.10, 1.20)
  )
  all_T  <- sort(unique(df$T))
  pick_T <- function(Tt) all_T[ which.min(abs(all_T - Tt)) ]
  
  nearest_idx <- function(x, target, used = integer(0)) {
    ord <- order(abs(x - target))
    if (!length(used)) return(ord[1])
    for (i in ord) if (!(i %in% used)) return(i)
    ord[1]
  }
  
  rows <- list()
  for (lab in names(target_Ts)) {
    Tsel <- pick_T(target_Ts[[lab]])
    dfT  <- df[abs(df$T - Tsel) < 1e-10 | df$T == Tsel, , drop = FALSE]
    if (!nrow(dfT)) next
    Fsel <- median(dfT$F, na.rm = TRUE)
    ok   <- is.finite(dfT$CallMid) & is.finite(dfT$K)
    if (!any(ok)) next
    used <- integer(0)
    
    for (rat in ratio_map[[lab]]) {
      Ktgt <- rat * Fsel
      idx  <- which(ok)[ nearest_idx(dfT$K[ok], Ktgt, used = match(used, which(ok))) ]
      used <- c(used, idx)
      
      # OTM-regel for typevalg
      typ <- if (rat < 1) "put" else "call"
      
      # Markedspris: call fra CallMid; put via paritet P = C - e^{-rT}(F - K)
      C_mid <- as.numeric(dfT$CallMid[idx])
      Tm    <- as.numeric(Tsel)
      Km    <- as.numeric(dfT$K[idx])
      disc  <- exp(-r * Tm)
      P_mid <- C_mid - disc * (Fsel - Km)
      mkt_price <- if (identical(typ,"call")) C_mid else P_mid
      
      rows[[length(rows)+1]] <- data.frame(
        type = typ,
        K    = Km,
        T    = Tm,
        ratio= rat,
        maturity_bucket = lab,
        fwd  = Fsel,
        market_price = mkt_price,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  # Unike punkter (type,K,T)
  out <- out[!duplicated(out[, c("type","K","T")]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

options_df_calib <- .build_G1_df(options_raw)
stopifnot(nrow(options_df_calib) == 30L)

# ---------- Markeds-ivs (BS) ----------
options_df_calib$market_iv <- mapply(function(typ, K, T, P) {
  bs_implied_vol(mkt_price = P, S = S0, K = K, T = T, r = r, q = q, type = typ)
}, options_df_calib$type, options_df_calib$K, options_df_calib$T, options_df_calib$market_price)

# Rydd og sorter
options_df_calib <- options_df_calib %>%
  filter(is.finite(market_iv), is.finite(K), is.finite(T), T > 0) %>%
  mutate(type = tolower(trimws(type)),
         K = K_round(K), T = T_round(T)) %>%
  arrange(T, match(type, c("call","put")), K)
stopifnot(nrow(options_df_calib) == 30L)

cat("Kalibreringssett (G1) by bucket:\n")
print(options_df_calib %>%
        group_by(maturity_bucket, T) %>%
        summarise(n = dplyr::n(),
                  ratios = paste(sort(unique(ratio)), collapse=", "),
                  .groups="drop"))

# ---------- Kostnadsfunksjon (likt vektet SSE p?? IV) ----------
.cost_fun <- function(param) {
  H   <- param[1]; rho <- param[2]; eta <- param[3]
  alpha <- H - 0.5
  
  # Pris??r alle 30 OTM-opsjoner med din rB-turbo pricer
  pr_res <- price_many_rB_turbo(
    options_df = options_df_calib[, c("type","K","T")],
    S0 = S0, r = r, q = q,
    alpha = alpha, eta = eta, rho = rho,
    xi0_curve = xi0_curve,
    N = NCal, M = MCal, antithetic = TRUE, seed = seed,
    clip_iv = TRUE
  )
  
  # Normaliser n??kler og join
  pr_res <- pr_res %>%
    mutate(type = tolower(trimws(as.character(type))),
           K = K_round(K),
           T = T_round(T))
  joined <- merge(
    pr_res[, c("type","K","T","price")],
    options_df_calib[, c("type","K","T","market_iv")],
    by = c("type","K","T"), all = FALSE
  )
  
  # Modell-ivs fra prisene
  iv_model <- mapply(function(typ, K, T, P) {
    bs_implied_vol(mkt_price = P, S = S0, K = K, T = T, r = r, q = q, type = typ)
  }, joined$type, joined$K, joined$T, joined$price)
  
  # Likt vektet SSE (kap. 10)
  diff <- iv_model - joined$market_iv
  sum(diff^2)
}

# ---------- SLSQP-setup ----------
lower_bounds <- c(H = 0.01, rho = -0.99, eta = 0.50)
upper_bounds <- c(H = 0.49, rho = -0.01, eta = 5.50)
init_param   <- c(H = 0.1, rho = -0.50, eta = 2.00)

# Kj??r SLSQP (med numerisk gradient fra numDeriv)
result <- nloptr::nloptr(
  x0   = init_param,
  eval_f = .cost_fun,
  eval_grad_f = function(p) numDeriv::grad(.cost_fun, p, method = "simple"),
  lb = lower_bounds,
  ub = upper_bounds,
  opts = list(
    algorithm  = "NLOPT_LD_SLSQP",
    xtol_rel   = 1e-6,
    ftol_rel   = 1e-8,
    maxeval    = 200,
    print_level= 5
  )
)

# ---------- Utskrift ----------
calibrated_params <- result$solution
cat(sprintf("\nKalibrerte parametre (H*, rho*, eta*): %.6f, %.6f, %.6f\n",
            calibrated_params[1], calibrated_params[2], calibrated_params[3]))
cat("Initial kost (x0) =", .cost_fun(init_param), "\n")
str(list(status = result$status, message = result$message,
         f_best = result$objective, x_best = result$solution), max.level = 1)


HCal <- result$solution[1]; HCal
RhoCal <- result$solution[2]; RhoCal
EtaCal <- result$solution[3]; EtaCal

HCal; RhoCal; EtaCal
