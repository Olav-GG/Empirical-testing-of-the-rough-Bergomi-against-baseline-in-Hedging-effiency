# install.packages("DT")  # <- kj??r ??n gang hvis du ikke har DT
library(DT)

# 1) Hjelper: legg kolonner i en ryddig rekkef??lge hvis de finnes
.reorder_cols <- function(df, priority) {
  keep <- intersect(priority, names(df))
  rest <- setdiff(names(df), keep)
  df[, c(keep, rest), drop = FALSE]
}

# 2) Hjelper: rund av tall kolonnevis (default = 6 desimaler)
.round_numeric <- function(df, digits = 6) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) round(x, digits))
  df
}

# 3) Hoved: vis pent i RStudio med DT (interaktivt), fallback til base::print
show_pretty <- function(df, title = NULL, digits = 6) {
  df <- .round_numeric(df, digits)
  
  # Fors??k ?? velge en god kolonnerekkef??lge ut fra hva som faktisk finnes
  bs_order <- c("type","K","T","price","iv","delta","gamma","vega","theta","rho",
                "d1","d2","Nd1","Nd2","lower","upper","in_bounds",
                "se","ci_lo","ci_hi", "F0","disc","div","moneyness","option_id")
  df <- .reorder_cols(df, bs_order)
  
  if ("datatable" %in% getNamespaceExports("DT")) {
    dt <- datatable(
      df,
      rownames = FALSE,
      caption = if (!is.null(title)) htmltools::tags$caption(style="caption-side: top; text-align:left;", title),
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50, 100),
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "tip"  # table + info + pagination (rent og pent)
      )
    )
    
    # Litt formatering hvis kolonner finnes
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    if (length(num_cols)) dt <- DT::formatRound(dt, columns = num_cols, digits = digits)
    if ("in_bounds" %in% names(df)) {
      dt <- DT::formatStyle(
        dt, "in_bounds",
        backgroundColor = DT::styleEqual(c(TRUE, FALSE), c(NA, "#ffe6e6"))
      )
    }
    return(dt)
  }
  
  # Fallback uten DT
  print(df)
  invisible(df)
}

# ---- Bruk: pene tabeller for BS og rB ----
show_pretty(bs_res, title = "Black???Scholes (BS) ??? priser, greeks og IV")
show_pretty(rb_res, title = "rBergomi (rB) ??? MC-priser, SE/CI og BS-implied vol")



## ---- sammenlign T ~ 1 ----
tol <- 1e-8
bs_T1 <- subset(bs_res, abs(T - 1) < tol, select = c(type, K, T, price))
rb_T1 <- subset(rb_res, abs(T - 1) < tol, select = c(type, K, T, price))

cmp_T1 <- merge(bs_T1, rb_T1, by = c("type","K","T"), suffixes = c("_BS","_rB"))
cmp_T1$diff     <- cmp_T1$price_rB - cmp_T1$price_BS
cmp_T1$pct_diff <- 100 * cmp_T1$diff / cmp_T1$price_BS
cmp_T1 <- cmp_T1[order(cmp_T1$type, cmp_T1$K), ]

# Vis pent (bruk show_pretty hvis du la den inn, ellers print)
if (exists("show_pretty")) show_pretty(cmp_T1, title = "BS vs rB (T ??? 1 ??r)", digits = 4) else print(cmp_T1)


## ---- generell sammenligner for hele griden ----
compare_all <- function(bs_res, rb_res) {
  b <- bs_res[, c("type","K","T","price")]
  r <- rb_res[, c("type","K","T","price")]
  out <- merge(b, r, by = c("type","K","T"), suffixes = c("_BS","_rB"))
  out$diff     <- out$price_rB - out$price_BS
  out$pct_diff <- 100 * out$diff / out$price_BS
  out[order(out$T, match(out$type, c("call","put")), out$K), ]
}

cmp_all <- compare_all(bs_res, rb_res)
if (exists("show_pretty")) show_pretty(cmp_all, title = "BS vs rB (alle forfall)", digits = 4) else print(cmp_all)

