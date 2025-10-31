
# ?????? One-date builder (same logic as your OOS per-date) ???????????????????????????????????????????????????????????????????????????
build_cmp_for_date <- function(d, spxIvolList, r, q,
                               H_star, Rho_star, Eta_star,
                               N=500L, M=2000L, seed_base=1L,
                               antithetic=TRUE, ncores=1) {
  stopifnot(d %in% names(spxIvolList))
  chain <- spxIvolList[[d]]
  stopifnot(is.data.frame(chain), nrow(chain) > 0)
  
  chain <- chain %>%
    dplyr::mutate(
      T = as.numeric(Texp),
      K = as.numeric(Strike),
      F = as.numeric(Fwd),
      CallMid = suppressWarnings(as.numeric(CallMid))
    )
  
  # Market table (calls + puts via parity)
  df_oos <- chain %>%
    dplyr::filter(is.finite(T), T>0, is.finite(K), is.finite(CallMid), is.finite(F)) %>%
    dplyr::mutate(call_px = CallMid,
                  put_px  = call_px - exp(-r*T)*(F - K))
  stopifnot(nrow(df_oos) > 0)
  
  options_df <- dplyr::bind_rows(
    df_oos %>% dplyr::transmute(type="call", K, T=round(T,10), mkt_price=call_px, date=d),
    df_oos %>% dplyr::transmute(type="put",  K, T=round(T,10), mkt_price=put_px,  date=d)
  )
  
  # S0 + ??0(T) from your calibration helper
  S0_today <- build_S0_from_chain(chain, r, q)
  calib_df_today <- build_calibration_df(
    spxIvolList, date_key = d,
    grid_list = lapply(unique(chain$Expiry),
                       function(e) list(expiration=e, ratios=c(0.8,0.9,1.0,1.1,1.2))),
    S0 = S0_today, r = r, q = q, use_chain_fwd = TRUE
  )
  xi0_curve_today <- build_xi0_from_ATM(
    calib_df_today, extend_to = max(options_df$T, na.rm=TRUE)
  )
  Sigma_BS <- sqrt(max(xi0_curve_today$forward_var[1], 1e-12))
  
  # Ensure even M when antithetic
  if (antithetic) M <- 2L * ceiling(M/2L)
  seed_i <- as.integer(abs((seed_base + sum(utf8ToInt(d))*131) %% .Machine$integer.max))
  
  # rB prices
  pr_rB <- price_many_rB_turbo(
    options_df = options_df[, c("type","K","T")],
    S0 = S0_today, r = r, q = q,
    alpha = H_star - 0.5, eta = Eta_star, rho = Rho_star,
    xi0_curve = xi0_curve_today,
    N = N, M = M, antithetic = antithetic, seed = seed_i, ncores = ncores, clip_iv = TRUE
  ) %>%
    dplyr::transmute(type = tolower(trimws(type)),
                     K = as.numeric(K), T = round(as.numeric(T),10),
                     price_rB = price)
  
  # Join and compute IVs
  cmp <- options_df %>%
    dplyr::mutate(type = tolower(trimws(type)),
                  K = as.numeric(K), T = round(as.numeric(T),10)) %>%
    dplyr::left_join(pr_rB, by=c("type","K","T")) %>%
    dplyr::mutate(
      price_BS = mapply(function(tp, K_, T_) {
        bs_price(S=S0_today, K=K_, T=T_, r=r, q=q, sigma=Sigma_BS, type=tp)
      }, type, K, T),
      mkt_iv = safe_bs_iv_vec(type, mkt_price, S=S0_today, Ks=K, Ts=T, r=r, q=q),
      iv_rB  = safe_bs_iv_vec(type, price_rB,  S=S0_today, Ks=K, Ts=T, r=r, q=q),
      iv_BS  = safe_bs_iv_vec(type, price_BS,  S=S0_today, Ks=K, Ts=T, r=r, q=q)
    )
  
  # Attach useful extras for plotting
  cmp$S0 <- S0_today
  cmp$r  <- r; cmp$q <- q
  cmp$F  <- cmp$S0 * exp((cmp$r - cmp$q) * cmp$T)
  cmp$m  <- cmp$K / cmp$F                 # moneyness
  cmp$k  <- log(pmax(cmp$m, 1e-12))       # log-moneyness
  attr(cmp, "xi0_curve") <- xi0_curve_today
  cmp
}

library(ggplot2)

iv_surface_plot <- function(cmp, which = c("market","rB","BS"),
                            x = c("moneyness","logm","strike"),
                            geom = c("tile","raster")) {
  which <- match.arg(which)
  x     <- match.arg(x)
  geom  <- match.arg(geom)
  
  iv_col <- switch(which,
                   market = "mkt_iv",
                   rB     = "iv_rB",
                   BS     = "iv_BS")
  stopifnot(iv_col %in% names(cmp))
  
  # choose x-variable
  cmp$X <- switch(x,
                  moneyness = cmp$m,
                  logm      = cmp$k,
                  strike    = cmp$K)
  x_lab <- switch(x,
                  moneyness = "K / F(T)",
                  logm      = "log(K / F(T))",
                  strike    = "Strike K")
  
  df <- cmp %>% dplyr::filter(is.finite(.data[[iv_col]]), is.finite(X), is.finite(T))
  if (!nrow(df)) stop("No finite points to plot.")
  
  p <- ggplot(df, aes(X, T, fill = .data[[iv_col]]))
  p <- if (geom == "tile") p + geom_tile() else p + geom_raster(interpolate = TRUE)
  p +
    scale_fill_viridis_c(option = "C", name = "IV") +
    labs(x = x_lab, y = "Maturity T (years)",
         title = sprintf("Implied Volatility Surface (%s)", which)) +
    theme_minimal(base_size = 12)
}

iv_residual_heatmap <- function(cmp, model = c("rB","BS"),
                                x = c("moneyness","logm","strike"),
                                limits = NULL) {
  model <- match.arg(model); x <- match.arg(x)
  col <- if (model == "rB") "iv_rB" else "iv_BS"
  stopifnot(col %in% names(cmp))
  
  cmp$X <- switch(x, moneyness = cmp$m, logm = cmp$k, strike = cmp$K)
  df <- cmp %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(mkt_iv), is.finite(X), is.finite(T)) %>%
    dplyr::mutate(resid = .data[[col]] - mkt_iv)
  
  gg <- ggplot(df, aes(X, T, fill = resid)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(name = "??IV", low = "#3B82F6", mid = "white", high = "#EF4444",
                         midpoint = 0, limits = limits) +
    labs(x = switch(x, moneyness="K / F(T)", logm="log(K/F(T))", strike="Strike K"),
         y = "Maturity T", title = sprintf("Residual IV heatmap (%s ??? market)", model)) +
    theme_minimal(base_size = 12)
  gg
}

plot_smiles <- function(cmp, which = c("market","rB","BS"),
                        maturities = NULL, x = c("moneyness","logm","strike")) {
  which <- match.arg(which); x <- match.arg(x)
  col <- switch(which, market="mkt_iv", rB="iv_rB", BS="iv_BS")
  stopifnot(col %in% names(cmp))
  
  cmp$X <- switch(x, moneyness=cmp$m, logm=cmp$k, strike=cmp$K)
  df <- cmp %>% dplyr::filter(is.finite(.data[[col]]), is.finite(T), is.finite(X))
  if (!is.null(maturities)) {
    pick <- sapply(maturities, function(Tt) {
      idx <- which.min(abs(df$T - Tt)); df$T[idx]
    })
    df <- df %>% dplyr::filter(T %in% pick)
  }
  ggplot(df, aes(X, .data[[col]], color = factor(round(T,3)))) +
    geom_line(alpha = 0.8) +
    labs(x = switch(x, moneyness="K / F(T)", logm="log(K/F(T))", strike="Strike K"),
         y = "Implied vol",
         color = "T (years)",
         title = sprintf("Smiles by maturity ??? %s", which)) +
    theme_minimal(base_size = 12)
}


plot_atm_term_structure <- function(cmp, xi0_curve = attr(cmp, "xi0_curve")) {
  # ATM proxy: pick closest-to-forward per T
  atm <- cmp %>%
    dplyr::group_by(T) %>%
    dplyr::slice_min(order_by = abs(K - median(F, na.rm=TRUE)), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  p1 <- ggplot(atm, aes(T, mkt_iv)) +
    geom_point() + geom_line() +
    labs(title="ATM implied vol (market)", x="T (years)", y="IV") +
    theme_minimal(base_size = 12)
  
  if (is.null(xi0_curve) || !all(c("T_start","T_end","forward_var") %in% names(xi0_curve))) return(p1)
  
  xi_plot <- xi0_curve %>%
    dplyr::mutate(T_mid = 0.5*(T_start + T_end), inst_vol = sqrt(pmax(forward_var,0)))
  
  p2 <- ggplot(xi_plot, aes(T_mid, inst_vol)) +
    geom_step(direction = "mid") +
    labs(title="????? instantaneous vol (sqrt(forward_var))", x="T (years)", y="??_inst(T)") +
    theme_minimal(base_size = 12)
  
  list(atm_iv = p1, xi0 = p2)
}


# Pick an OOS date:
d <- date_keys_oos[1]

# Build comparison table for that date (respects antithetics + cores)
cmp <- build_cmp_for_date(
  d, spxIvolList, r, q,
  H_star = HCal, Rho_star = RhoCal, Eta_star = EtaCal,
  N = 500, M = 2000, antithetic = TRUE, ncores = 4
)

# Surfaces
iv_surface_plot(cmp, which="market", x="moneyness")
iv_surface_plot(cmp, which="rB",     x="moneyness")
iv_surface_plot(cmp, which="BS",     x="moneyness")

# Residuals
iv_residual_heatmap(cmp, model="rB", x="moneyness")

# Smiles
plot_smiles(cmp, which="market", maturities=c(1/12, 0.5, 1, 2))

# ATM term structure + ?????(T)
plots <- plot_atm_term_structure(cmp)
plots$atm_iv; plots$xi0




library(dplyr)

.grid_to_surface <- function(x, y, z, nx=60, ny=60, method=c("mba","akima"), extend=FALSE) {
  method <- match.arg(method)
  stopifnot(length(x)==length(y), length(y)==length(z))
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  x <- x[ok]; y <- y[ok]; z <- z[ok]
  if (!length(z)) stop("No finite points to grid.")
  
  # Regular target grid
  xo <- seq(min(x), max(x), length.out = nx)
  yo <- seq(min(y), max(y), length.out = ny)
  
  if (method == "mba" && requireNamespace("MBA", quietly=TRUE)) {
    inp <- cbind(x, y, z)
    est <- MBA::mba.surf(inp, no.X = nx, no.Y = ny, extend = extend)$xyz.est
    list(x = est$x, y = est$y, z = est$z)
  } else {
    if (!requireNamespace("akima", quietly=TRUE))
      stop("Neither MBA nor akima available. Please install one of them.")
    itp <- akima::interp(x, y, z, xo = xo, yo = yo, duplicate = "mean", linear = TRUE)
    list(x = itp$x, y = itp$y, z = itp$z)
  }
}


library(plotly)

iv_surface_3d <- function(cmp,
                          which = c("market","rB","BS"),
                          xaxis = c("moneyness","logm","strike"),
                          nx = 60, ny = 60,
                          method = c("mba","akima")) {
  which  <- match.arg(which)
  xaxis  <- match.arg(xaxis)
  method <- match.arg(method)
  
  iv_col <- switch(which, market = "mkt_iv", rB = "iv_rB", BS = "iv_BS")
  if (!iv_col %in% names(cmp)) stop("Column ", iv_col, " not found in cmp.")
  
  # Choose surface X variable
  X <- switch(xaxis,
              moneyness = cmp$m,
              logm      = cmp$k,
              strike    = cmp$K)
  x_lab <- switch(xaxis,
                  moneyness = "K / F(T)",
                  logm      = "log(K / F(T))",
                  strike    = "Strike K")
  
  df <- cmp %>% filter(is.finite(.data[[iv_col]]), is.finite(T), is.finite(X))
  if (!nrow(df)) stop("No finite points to plot.")
  
  # Deduplicate close points to make gridding robust
  df <- df %>%
    mutate(Xr = round(X, 6), Tr = round(T, 6)) %>%
    group_by(Xr, Tr) %>%
    summarise(IV = mean(.data[[iv_col]], na.rm=TRUE), .groups="drop")
  
  S <- .grid_to_surface(x = df$Xr, y = df$Tr, z = df$IV,
                        nx = nx, ny = ny, method = method)
  
  plot_ly(x = S$x, y = S$y, z = S$z, type = "surface",
          colorscale = "Viridis", showscale = TRUE) %>%
    layout(
      title = paste("Implied Volatility Surface ???", which),
      scene = list(
        xaxis = list(title = x_lab),
        yaxis = list(title = "Maturity T (years)"),
        zaxis = list(title = "IV")
      )
    )
}

residual_surface_3d <- function(
  cmp,
  model  = c("rB","BS"),
  xaxis  = c("moneyness","logm","strike"),
  nx = 60, ny = 60,
  method = c("mba","akima"),
  S0 = NULL, r = 0, q = 0
) {
  model  <- match.arg(model)
  xaxis  <- match.arg(xaxis)
  method <- match.arg(method)

  # Choose model IV column
  col <- if (model == "rB") "iv_rB" else "iv_BS"
  if (!all(c(col, "mkt_iv", "T", "K") %in% names(cmp))) {
    stop("cmp must contain columns: '", col, "', 'mkt_iv', 'T', 'K'.")
  }

  # Forward vector: use cmp$F if present; else derive from S0,r,q if provided
  if ("F" %in% names(cmp)) {
    F_vec <- as.numeric(cmp$F)
  } else if (!is.null(S0)) {
    F_vec <- as.numeric(S0) * exp((r - q) * as.numeric(cmp$T))
  } else {
    F_vec <- rep(NA_real_, nrow(cmp))
  }

  # X coordinate
  Xraw <- switch(
    xaxis,
    moneyness = if ("m" %in% names(cmp)) as.numeric(cmp$m)
                else if (all(is.finite(F_vec))) as.numeric(cmp$K) / F_vec
                else stop("Need 'F' in cmp or pass S0,r,q to derive moneyness."),
    logm      = if ("k" %in% names(cmp)) as.numeric(cmp$k)
                else if (all(is.finite(F_vec))) log(as.numeric(cmp$K) / F_vec)
                else stop("Need 'F' in cmp or pass S0,r,q to derive log-moneyness."),
    strike    = as.numeric(cmp$K)
  )
  x_lab <- switch(xaxis,
                  moneyness = "K / F(T)",
                  logm      = "log(K / F(T))",
                  strike    = "Strike K")

  # Build unique (X,T) grid with mean residuals
  df <- data.frame(
    X = as.numeric(Xraw),
    T = as.numeric(cmp$T),
    resid = as.numeric(cmp[[col]] - cmp$mkt_iv)
  )
  df <- df[is.finite(df$X) & is.finite(df$T) & is.finite(df$resid), , drop = FALSE]

  # round to merge duplicates safely before gridding
  df$Xr <- round(df$X, 6)
  df$Tr <- round(df$T, 6)

  df2 <- df |>
    dplyr::group_by(.data$Xr, .data$Tr) |>
    dplyr::summarise(dIV = mean(.data$resid), .groups = "drop")

  if (nrow(df2) < 4L) stop("Not enough points to build a surface.")

  # Interpolate onto a regular grid
  if (method == "mba") {
    if (!requireNamespace("MBA", quietly = TRUE)) {
      stop("Package 'MBA' is required for method='mba'. Install it via install.packages('MBA').")
    }
    xyz <- data.frame(x = df2$Xr, y = df2$Tr, z = df2$dIV)
    surf <- MBA::mba.surf(xyz, no.X = nx, no.Y = ny, extend = TRUE)
    gx <- surf$xyz.est$x
    gy <- surf$xyz.est$y
    gz <- surf$xyz.est$z  # matrix [length(x) x length(y)]
  } else {
    if (!requireNamespace("akima", quietly = TRUE)) {
      stop("Package 'akima' is required for method='akima'. Install it via install.packages('akima').")
    }
    itp <- akima::interp(x = df2$Xr, y = df2$Tr, z = df2$dIV,
                         duplicate = "mean", nx = nx, ny = ny, extrap = FALSE)
    gx <- itp$x
    gy <- itp$y
    gz <- itp$z  # matrix [length(x) x length(y)]
  }

  # 3D plotly surface (plotly expects z with dims [length(y) x length(x)])
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for 3D plotting. Install it via install.packages('plotly').")
  }
  p <- plotly::plot_ly(
    x = gx, y = gy, z = ~t(gz), type = "surface"
  ) |>
    plotly::layout(
      title = paste0("Residual IV surface (", model, " ??? market)"),
      scene = list(
        xaxis = list(title = x_lab),
        yaxis = list(title = "T (years)"),
        zaxis = list(title = "dIV (model ??? mkt)")
      )
    )

  # return plot and grid
  attr(p, "grid") <- list(x = gx, y = gy, z = gz)
  return(p)
}



cmp <- cmp %>% mutate(m = K / F, k = log(m))

iv_surface_3d(cmp, which="market", xaxis="moneyness", nx=80, ny=60, method="mba")
iv_surface_3d(cmp, which="rB",     xaxis="moneyness", nx=80, ny=60, method="mba")
residual_surface_3d(cmp, model="rB", xaxis="moneyness", nx=80, ny=60, method="mba")

