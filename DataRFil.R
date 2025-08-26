S0 <- 1000 # underliggende pris i dag
K <- 1500 # strike
T_days <- 10 # dager til forfall
r <- 0.00 # risikofri rente (kontinuerlig)
q <- 0.00 # utbytte-rate (sett 0 om ikke relevant)
Sigma <- 0.235
# ===========================


# Avledet: tid i ??r
T <- T_days/252


# Opsjonsbeskrivelse (helt enkel list)
option <- list(
  type = "call",
  K = K,
  T = T,
  r = r,
  q = q,
  name = sprintf("Call_K%d_%dd", K, T_days)
)


# For klarhet
S <- S0

moneyness <- S / K

