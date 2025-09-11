#----------------Load and clean the dataset----------

library(readxl)
library(janitor)
library(readr)   # parse_number()
library(lubridate)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

Raw_df <- read_xlsx("SPX Opsjonsdata.xlsx")
VIX30 <- read_xlsx("VIX 30.xlsx")
VIX60 <- read_xlsx("VIX 60.xlsx")
head(Raw_df)
head(VIX30)
head(VIX60)

df_clean <- Raw_df[-1,]

# Anta at df_clean er som i head(df_clean) du viste.
nm <- names(df_clean)

# 1) Finn date-kolonner: de heter typisk "Security...N" i Bloomberg-eksporter
date_idx <- which(grepl("^Security\\.\\.\\.", nm))
# 2) Tilh??rende "Last"-kolonne ligger rett til h??yre
last_idx <- date_idx + 1

# Sikkerhet: behold bare par som faktisk finnes
ok <- last_idx <= length(nm)
date_idx <- date_idx[ok]
last_idx <- last_idx[ok]

# 3) Gi date-kolonnene navn "Date_<instrumentnavn>"
new_names <- nm
new_names[date_idx] <- paste0("Date_", nm[last_idx])
names(df_clean) <- new_names

# 4) Konverter typer effektivt:
date_cols <- new_names[date_idx]
last_cols <- new_names[last_idx]

# Datoer: Excel-serial -> Date (origin 1899-12-30)
df_clean[date_cols] <- lapply(df_clean[date_cols],
                              function(x) as.Date(as.numeric(x), origin = "1899-12-30"))
# Priser: numeric
df_clean[last_cols] <- lapply(df_clean[last_cols], parse_number)

# Sjekk:
head(df_clean)

tail(df_clean)

#=======Gather data for input=========

# ---- 1) Finn opsjonskolonnene (ikke Date_-kolonner, og m?? ha C/P+strike i navnet)
last_cols <- setdiff(names(df_clean), grep("^Date_", names(df_clean), value = TRUE))
opt_cols  <- last_cols[str_detect(last_cols, "\\s[CP]\\d")]

# ---- 2) Parse metadata fra navnene
# M??nster: "<underlying> US <mm/dd/yy> <C|P><strike> Index"
opt_meta <- tibble(security = opt_cols) %>%
  extract(
    col   = security,
    into  = c("underlying", "market", "maturity_chr", "type_code", "strike_chr"),
    regex = "^(.+?)\\s+(\\w{2})\\s+(\\d{2}/\\d{2}/\\d{2})\\s+([CP])(\\d+(?:\\.?\\d+)?)\\s+Index$",
    remove = FALSE
  ) %>%
  mutate(
    maturity = mdy(maturity_chr),                    # 07/18/25 -> Date (2025-07-18)
    type     = recode(type_code, C = "Call", P = "Put"),
    strike   = parse_number(strike_chr),
    .keep    = "unused"
  )
# opt_meta har n??: security, underlying, market, maturity(Date), type, strike(numeric)

# ---- 3) (Valgfritt) Bygg et langt datasett med dato+pris for hver opsjon
# Knytter hver opsjonskolonne til riktig Date_-kolonne ved navn
date_cols <- paste0("Date_", opt_cols)
have <- date_cols %in% names(df_clean)
opt_cols  <- opt_cols[have]
date_cols <- date_cols[have]

# Effektivt og korrekt parvis "stack" av (dato, pris) per security:
opt_long <- purrr::map2_dfr(opt_cols, date_cols, function(lc, dc) {
  tibble(
    security = lc,
    date     = as.Date(df_clean[[dc]]),         # skal allerede v??re Date; as.Date() er no-op
    last     = as.numeric(df_clean[[lc]])
  )
})

# Sl?? p?? metadata
opt_long <- opt_long %>% left_join(opt_meta, by = "security")

# Eksempelsp??rringer:
# - Finn P6295 som forfaller 2025-07-18:
opt_meta %>% filter(type == "Put", maturity == as.Date("2025-07-18"), strike == 6295)

# - Sist observerte pris for en gitt opsjon:
opt_long %>%
  filter(type == "Put", maturity == as.Date("2025-07-18"), strike == 6295) %>%
  arrange(date) %>%
  summarise(last_obs = dplyr::last(na.omit(last))) %>%
  pull(last_obs)

head(opt_long)


#=============input=================
S0 <- tail(df_clean$`SPX Index`, 1); S0
r <- 0.0425
q <- 0.0156
Sigma <- tail(VIX$PX_LAST, 1)/100; Sigma
day_count <- 365

last_cols <- setdiff(nm, grep("^Date_", nm, value = TRUE))
opt_cols  <- last_cols[str_detect(last_cols, "\\s[CP]\\d")]  # har C/P + strike i navnet

opt_meta <- tibble(security = opt_cols) %>%
  extract(
    col   = security,
    into  = c("underlying", "market", "maturity_chr", "type_code", "strike_chr"),
    regex = "^(.+?)\\s+(\\w{2})\\s+(\\d{2}/\\d{2}/\\d{2})\\s+([CP])(\\d+(?:\\.?\\d+)?)\\s+Index$",
    remove = FALSE
  ) %>%
  mutate(
    maturity = mdy(maturity_chr),                 # 07/18/25 -> 2025-07-18
    type     = recode(type_code, C = "call", P = "put"),
    K        = parse_number(strike_chr)
  ) %>%
  select(security, underlying, market, maturity, type, K)

start_date <- tail(df_clean$`Date_SPX Index`,1); start_date

options_df <- opt_meta %>%
  transmute(
    type,
    K,
    T_days = as.integer(maturity - start_date),
    T      = pmax(0, T_days / day_count)
  ) %>%
  filter(T > 0) %>%            # dropp utl??pte/ikke-positive T
  arrange(T, type, K)

T_max <- max(options_df$T)


start_date <- as.Date("2025-05-30")
end_date  <- max(opt_meta$maturity)

T <- as.integer(end_date-start_date); T


## ------- 1) Velg riktige kolonner (enkelt) -------
# Sett disse til de faktiske navnene i filene dine:
date_col_30  <- grep("Date|Dato", names(VIX30), ignore.case = TRUE, value = TRUE)[1]
price_col_30 <- if ("PX_LAST" %in% names(VIX30)) "PX_LAST" else names(VIX30)[2]

date_col_60  <- grep("Date|Dato", names(VIX60), ignore.case = TRUE, value = TRUE)[1]
price_col_60 <- if ("PX_LAST" %in% names(VIX60)) "PX_LAST" else names(VIX60)[2]

## ------- 2) Gj??r om til (date, vol) serier -------
toDate <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(as.numeric(x), origin = "1899-12-30")  # funker for Excel-serial og chr-num
}

v30 <- tibble::tibble(
  date = toDate(VIX30[[date_col_30]]),
  vol  = as.numeric(VIX30[[price_col_30]])/100  # 20 -> 0.20
) |> dplyr::filter(!is.na(date), !is.na(vol)) |> dplyr::arrange(date)

v60 <- tibble::tibble(
  date = toDate(VIX60[[date_col_60]]),
  vol  = as.numeric(VIX60[[price_col_60]])/100
) |> dplyr::filter(!is.na(date), !is.na(vol)) |> dplyr::arrange(date)

## ------- 3) Ta verdien p??/like f??r start_date -------
pick_on_or_before <- function(df, start_date) {
  cand <- dplyr::filter(df, date <= start_date)
  if (nrow(cand)) dplyr::slice_tail(cand, n = 1) else dplyr::slice_head(df, n = 1)
}

sigma30 <- pick_on_or_before(v30, start_date)$vol[1]  # desimal
sigma60 <- pick_on_or_before(v60, start_date)$vol[1]

## ------- 4) Bygg enkel piecewise-constant xi0(t) -------
T1 <- 30 / day_count
T2 <- 60 / day_count
w30 <- T1 * sigma30^2
w60 <- T2 * sigma60^2


xi1 <- w30 / T1                              # 0???30d
xi2 <- (w60 - w30) / (T2 - T1)               # 30???60d (forlenges flatt >60d)

dt      <- T_max / N
t_grid  <- (1:N) * dt
xi0_curve <- ifelse(t_grid <= T1, xi1, xi2)  # vektor lengde N
xi0_curve <- pmax(xi0_curve, 1e-10)          # numerisk gulv

# (valgfritt) niv??et ved t -> 0:
xi0_level <- xi0_curve[1]; xi0_level
xi0_curve

## ------- 5) Kjapp sanity (valgfritt) -------
ggplot2::qplot(t_grid, sqrt(xi0_curve), geom = "line",
               xlab = "T (??r)", ylab = "Forward vol (desimal)")


