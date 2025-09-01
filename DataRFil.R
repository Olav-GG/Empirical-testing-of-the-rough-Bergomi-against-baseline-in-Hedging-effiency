## ==== Grunnparametre (input) ===========================================
S0    <- 1000      # underliggende i dag
r     <- 0.00      # risikofri rente (kontinuerlig)
q     <- 0.00     # utbytte-rate (kontinuerlig)
Sigma <- 0.235     # BS-referansevol (til null-test / sanity)
day_count <- 252   # handelsdager pr ??r (bruk 365 hvis du heller vil ha ACT/365)

## ??n referanseopsjon (samme som du hadde)
K      <- 1000
T_days <- 252
T      <- T_days / day_count
T
T_max <- 1.0 
## For klarhet
S <- S0
moneyness <- S / K

# Opsjonsbeskrivelse (helt enkel list)
option <- list(
  type = "call",
  K = K,
  T = T,
  r = r,
  q = q,
  name = sprintf("Call_K%d_%dd", K, T_days)
)

S0 <- S0; r <- r; q <- q
alpha <- -0.43; eta <- 1.9; Rho <- -0.9
xi0_level <- 0.235^2; H <- alpha + 0.5
N <- min(500, ceiling(day_count * T_max))  # ??? ??n tidssteg per handelsdag
M <- 1000; seed <- 1

## ==== Opsjons-datasett (grid) ==========================================
# Behold din grid ??? vi standardiserer bare litt rundt typer og numerikk
options_df <- expand.grid(
  type = c("call", "put"),
  K    = c(400, 600, 800, 900, 1000, 1100, 1200, 1400, 1600),
  T    = c(0.1, 0.2, 0.25, 0.4, 0.5, 0.75, 1.0),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
options_df$type <- tolower(trimws(options_df$type))
options_df$K    <- as.numeric(options_df$K)
options_df$T    <- as.numeric(options_df$T)

## (Valgfritt) Avledete felt som er nyttige senere (plott, sanity, m.m.)
options_df$F0         <- S0 * exp((r - q) * options_df$T)
options_df$k          <- log(options_df$K / S0)     # log-strike
options_df$moneyness  <- S0 / options_df$K
options_df$disc       <- exp(-r * options_df$T)
options_df$div        <- exp(-q * options_df$T)
options_df$option_id  <- sprintf("%s_K%d_T%.2f", options_df$type, options_df$K, options_df$T)

## ==== Datavalidering ===================================================
stopifnot(all(options_df$type %in% c("call","put")))
stopifnot(all(is.finite(options_df$K) & options_df$K > 0))
stopifnot(all(is.finite(options_df$T) & options_df$T > 0))
options_df <- unique(options_df)  # fjern ev. duplikater

## ==== Modenheter til simulering (brukes av rB-priseren) ================
T_vec <- sort(unique(options_df$T))    # alle forfall vi skal prise p??
T_max <- max(T_vec)                    # simuler opp til lengste forfall

## ==== Kort oppsummering i konsoll ======================================
print(head(options_df, 10))
cat("Unique maturities (years):", paste(T_vec, collapse=", "), "\n")
cat("T_max =", T_max, "years;  N =", N, "time steps;  M =", M, "paths\n")
cat("S0 =", S0, " r =", r, " q =", q, " Sigma(BS ref.) =", Sigma, "\n")

