# 1) N??kkel for (type, K, T, maturity) fra options_df
options_key <- options_df %>%
  mutate(maturity = as.Date(start_date) + as.integer(round(T * day_count))) %>%
  distinct(type, K, T, maturity)

# 2) Siste markedspris per opsjon (bruker opt_long/opt_meta du allerede lagde)
market_last <- opt_long %>%
  filter(!is.na(last)) %>%
  group_by(security) %>%
  slice_max(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    type     = tolower(type),                          # ensartet sm?? bokstaver
    K        = as.numeric(dplyr::coalesce(K, strike)), # K hvis finnes, ellers strike
    maturity = maturity,
    mkt_last = last
  )

# 3) Sl?? sammen alt og regn ut feil + markeds-IV
cmp <- options_key %>%
  left_join(rb_res %>% select(type, K, T, price_rB = price, iv_rB = iv),
            by = c("type","K","T")) %>%
  left_join(bs_res %>% select(type, K, T, price_BS = price, iv_BS = iv),
            by = c("type","K","T")) %>%
  left_join(market_last, by = c("type","K","maturity")) %>%
  mutate(
    iv_mkt = ifelse(
      is.na(mkt_last), NA_real_,
      mapply(function(typ, KK, TT, P) bs_implied_vol(P, S = S0, K = KK, T = TT, r = r, q = q, type = typ),
             type, K, T, mkt_last)
    ),
    err_rB_mkt = price_rB - mkt_last,
    err_BS_mkt = price_BS - mkt_last,
    rel_rB_mkt = (price_rB - mkt_last) / mkt_last,
    rel_BS_mkt = (price_BS - mkt_last) / mkt_last
  ) %>%
  arrange(maturity, type, K)

# 4) Oppsummering per forfall/type
by_maturity_type <- cmp %>%
  filter(!is.na(mkt_last)) %>%
  group_by(maturity, T, type) %>%
  summarise(
    n        = n(),
    MAE_rB   = mean(abs(err_rB_mkt)),
    RMSE_rB  = sqrt(mean(err_rB_mkt^2)),
    MAPE_rB  = mean(abs(rel_rB_mkt), na.rm = TRUE),
    MAE_BS   = mean(abs(err_BS_mkt)),
    RMSE_BS  = sqrt(mean(err_BS_mkt^2)),
    MAPE_BS  = mean(abs(rel_BS_mkt), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(maturity, type)

# 5) Visning (bruker show_pretty hvis den finnes)
if (exists("show_pretty")) {
  show_pretty(cmp %>% select(maturity, T, type, K, mkt_last, price_rB, price_BS, iv_mkt, iv_rB, iv_BS,
                             err_rB_mkt, err_BS_mkt, rel_rB_mkt, rel_BS_mkt),
              title = "Market vs rB vs Black???Scholes", digits = 4)
  show_pretty(by_maturity_type, title = "Feilm??l per forfall og type", digits = 4)
} else {
  print(head(cmp, 20))
  print(by_maturity_type)
}


library(DT)

## --- 1) Market vs rB vs BS: pen tabell ---
cmp_view <- cmp %>%
  dplyr::select(maturity, T, type, K,
                mkt_last, price_rB, price_BS,
                iv_mkt, iv_rB, iv_BS,
                err_rB_mkt, err_BS_mkt, rel_rB_mkt, rel_BS_mkt) %>%
  dplyr::arrange(maturity, type, K)

num_cols <- names(cmp_view)[vapply(cmp_view, is.numeric, logical(1))]

cmp_dt <- datatable(
  cmp_view,
  rownames = FALSE,
  caption = htmltools::tags$caption(style="caption-side: top; text-align:left;",
                                    "Market vs rB vs Black???Scholes"),
  options = list(
    pageLength = 15,
    lengthMenu = c(10,15,25,50,100),
    scrollX = TRUE,
    autoWidth = TRUE,
    dom = "tip"
  )
) %>%
  DT::formatRound(columns = num_cols, digits = 4) %>%
  # farge p?? signerte feil (gr??nn = modell under markedet, r??d = over)
  DT::formatStyle(c("err_rB_mkt","err_BS_mkt","rel_rB_mkt","rel_BS_mkt"),
                  color = DT::styleInterval(0, c("#2e7d32", "#c62828"))) %>%
  # bar-heat p?? absolutt feil (bruker absoluttverdi for bredde)
  DT::formatStyle("err_rB_mkt",
                  background = DT::styleColorBar(abs(cmp_view$err_rB_mkt), "lightsteelblue"),
                  backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center") %>%
  DT::formatStyle("err_BS_mkt",
                  background = DT::styleColorBar(abs(cmp_view$err_BS_mkt), "lightsteelblue"),
                  backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")

cmp_dt

## --- 2) Oppsummering per forfall/type: pen tabell ---
by_view <- by_maturity_type %>%
  dplyr::arrange(maturity, type)

by_num <- names(by_view)[vapply(by_view, is.numeric, logical(1))]

by_dt <- datatable(
  by_view,
  rownames = FALSE,
  caption = htmltools::tags$caption(style="caption-side: top; text-align:left;",
                                    "Feilm??l per forfall og type"),
  options = list(
    pageLength = 10,
    lengthMenu = c(10,20,50),
    scrollX = TRUE,
    autoWidth = TRUE,
    dom = "tip"
  )
) %>%
  DT::formatRound(columns = by_num, digits = 4) %>%
  # varm farge for h??y feilkost
  DT::formatStyle(c("MAE_rB","RMSE_rB","MAPE_rB","MAE_BS","RMSE_BS","MAPE_BS"),
                  background = DT::styleColorBar(by_view$RMSE_rB, "lavender"),
                  backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")

by_dt
