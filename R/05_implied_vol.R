library(rugarch)
library(quantmod)
library(yaml)
library(xts)

cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

W   <- cfg$rolling_window 
L   <- cfg$forecast_length
B   <- cfg$refit_every
H   <- 21
ann <- sqrt(252) * 100         # daily vol -> annualized, in %

#implied volatility
v2tx_raw <- suppressWarnings(
  read.table("data/v2tx_history.csv", sep = ";", header = TRUE,
             dec = ".", stringsAsFactors = FALSE)
)
v2tx <- xts(v2tx_raw$Indexvalue,
            as.Date(v2tx_raw$Date, format = "%d.%m.%Y"))
stopifnot(!any(is.na(index(v2tx))))
cat("V2TX:", as.character(start(v2tx)), "->", as.character(end(v2tx)), "\n")

make_spec <- function(model) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = cfg$distribution
  )
}

roll_avg_vol <- function(ret, model) {
  x <- as.numeric(ret); n <- length(x)
  starts <- seq(n - L, n, by = B)
  res <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    s <- starts[k]
    e <- min(s + B - 1, n)
    m <- e - s                              # out-of-sample obs in this block
    dat <- x[(s - W + 1):e]
    fit <- ugarchfit(make_spec(model), data = dat,
                     out.sample = m, solver = "hybrid")
    if (convergence(fit) != 0)
      warning(sprintf("%s: no convergence at origin index %d", model, s))
    fc  <- ugarchforecast(fit, n.ahead = H, n.roll = m)
    sg  <- as.matrix(sigma(fc))             # H x (m + 1)
    res[[k]] <- data.frame(idx = s:e,
                           fc  = sqrt(colMeans(sg^2)) * ann)
  }
  do.call(rbind, res)
}

#Newey-West standard error
nw_se <- function(x, lag = H - 1) {
  x <- x - mean(x); n <- length(x)
  s <- sum(x^2) / n
  for (k in 1:lag)
    s <- s + 2 * (1 - k / (lag + 1)) * sum(x[(k + 1):n] * x[1:(n - k)]) / n
  sqrt(s / n)
}

analyse <- function(ret, label) {
  cat(sprintf("\n==================== %s ====================\n", label))
  x <- as.numeric(ret); d <- index(ret); n <- length(x)
  
  cat("Rolling GARCH...\n");  g <- roll_avg_vol(ret, "sGARCH")
  cat("Rolling EGARCH...\n"); e <- roll_avg_vol(ret, "eGARCH")
  stopifnot(identical(g$idx, e$idx))
  
  df <- data.frame(date = d[g$idx], garch = g$fc, egarch = e$fc)
  # realized vol over the NEXT H trading days (uncentered, annualized %)
  df$rv_fwd <- sapply(g$idx, function(i)
    if (i + H <= n) sqrt(mean(x[(i + 1):(i + H)]^2)) * ann else NA_real_)
  df$implied <- as.numeric(v2tx)[match(df$date, index(v2tx))]
  
  cat(sprintf("Origins: %d (%s -> %s) | with implied vol: %d | with forward realized: %d\n",
              nrow(df), df$date[1], df$date[nrow(df)],
              sum(!is.na(df$implied)), sum(!is.na(df$rv_fwd))))
  
  ok <- complete.cases(df)
  d2 <- df[ok, ]
  
  cat("\nImplied (V2TX) vs model forecast:\n")
  for (m in c("garch", "egarch")) {
    cat(sprintf("  %-6s: corr = %.3f | mean(implied - forecast) = %+.2f vol pts\n",
                m, cor(d2$implied, d2[[m]]), mean(d2$implied - d2[[m]])))
  }
  
  vrp <- d2$implied - d2$rv_fwd
  cat(sprintf("\nVariance risk premium proxy (implied subsequent realized vol):\n"))
  cat(sprintf("  mean = %+.2f vol pts (Newey-West s.e. %.2f) | implied > realized on %.0f%% of days\n",
              mean(vrp), nw_se(vrp), 100 * mean(vrp > 0)))
  
  rmse <- function(a, b) sqrt(mean((a - b)^2))
  cat("\nForecasting the next-21-day realized volatility (RMSE, vol pts; mean error = forecast - realized):\n")
  cat(sprintf("  GARCH   : RMSE = %.2f | mean error = %+.2f\n",
              rmse(d2$garch, d2$rv_fwd), mean(d2$garch - d2$rv_fwd)))
  cat(sprintf("  EGARCH  : RMSE = %.2f | mean error = %+.2f\n",
              rmse(d2$egarch, d2$rv_fwd), mean(d2$egarch - d2$rv_fwd)))
  cat(sprintf("  Implied : RMSE = %.2f | mean error = %+.2f\n",
              rmse(d2$implied, d2$rv_fwd), mean(d2$implied - d2$rv_fwd)))
  cat("  (obs. overlap by 20 days: treat differences as descriptive)\n")
  
  cat("\nLatest origin:\n"); print(round(df[nrow(df), -1], 2))
  
  # ---- Plot ----
  ragg::agg_png(sprintf("figures/iv_vs_garch_%s.png", label),
                width = 1400, height = 650, units = "px", res = 110,
                background = "white")
  yl <- range(df[, c("garch", "egarch", "rv_fwd", "implied")], na.rm = TRUE)
  plot(df$date, df$implied, type = "l", col = "black", lwd = 2, ylim = yl,
       main = sprintf("%s: 21-day volatility forecasts vs VSTOXX (annualized, %%)", label),
       xlab = "", ylab = "%")
  lines(df$date, df$garch,  col = "steelblue", lwd = 1.5)
  lines(df$date, df$egarch, col = "darkred",   lwd = 1.5)
  lines(df$date, df$rv_fwd, col = "grey55",    lwd = 2, lty = 2)
  legend("topright",
         c("VSTOXX (implied)", "GARCH forecast", "EGARCH forecast",
           "Subsequent realized (21d)"),
         col = c("black", "steelblue", "darkred", "grey55"),
         lty = c(1, 1, 1, 2), lwd = c(2, 1.5, 1.5, 2), bty = "n")
  dev.off()
  
  write.csv(df, sprintf("outputs/iv_comparison_%s.csv", label), row.names = FALSE)
  invisible(df)
}

returns_cac <- readRDS("data/returns.rds")
res_cac <- analyse(returns_cac, "CAC40")
returns_stoxx <- tryCatch({
  px <- getSymbols("^STOXX50E", from = cfg$start_date, to = cfg$end_date,
                   auto.assign = FALSE)
  px <- Ad(px); px <- px[!is.na(px)]
  r  <- na.omit(diff(log(px))); colnames(r) <- "ret"; r
}, error = function(err) { message("STOXX50E download failed: ", conditionMessage(err)); NULL })

if (!is.null(returns_stoxx)) {
  saveRDS(returns_stoxx, "data/returns_stoxx50.rds")
  res_stoxx <- analyse(returns_stoxx, "STOXX50E")
}