library(rugarch)
library(quantmod)
library(yaml)
library(xts)

source("R/stats.R")
cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

W   <- cfg$rolling_window 
B   <- cfg$refit_every                    
H   <- 21                    
ann <- sqrt(252) * 100
data_start <- if (is.null(cfg$data_start)) "2000-01-01" else cfg$data_start
eval_start <- as.Date(if (is.null(cfg$eval_start)) "2008-01-01" else cfg$eval_start)

#implied vol 
read_stoxx_file <- function(path) {
  raw <- suppressWarnings(read.table(path, sep = ";", header = TRUE,
                                     dec = ".", stringsAsFactors = FALSE))
  if (!all(c("Date", "Indexvalue") %in% names(raw)))
    stop(path, ": expected columns Date;Symbol;Indexvalue, found: ",
         paste(names(raw), collapse = ", "))
  xts(raw$Indexvalue, as.Date(raw$Date, format = "%d.%m.%Y"))
}
v2tx <- read_stoxx_file("data/v2tx_history.csv")
cat("V2TX:", as.character(start(v2tx)), "->", as.character(end(v2tx)), "\n")

#returns
returns_yahoo <- function(symbol) {
  px <- getSymbols(symbol, from = data_start, to = cfg$end_date, auto.assign = FALSE)
  px <- Ad(px); px <- px[!is.na(px)]
  r  <- na.omit(diff(log(px))); colnames(r) <- "ret"; r
}
returns_file <- function(path) {
  px <- read_stoxx_file(path)
  px <- px[!is.na(px) & as.numeric(px) > 0]
  px <- px[paste0(data_start, "/")]
  r  <- na.omit(diff(log(px))); colnames(r) <- "ret"; r
}

make_spec <- function(model) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = cfg$distribution
  )
}

# Fit with fallback solvers; NULL if every solver fails
fit_safe <- function(model, dat, m) {
  for (sv in c("hybrid", "nlminb", "solnp")) {
    fit <- tryCatch(suppressWarnings(
      ugarchfit(make_spec(model), data = dat, out.sample = m, solver = sv)),
      error = function(e) NULL)
    if (!is.null(fit) && convergence(fit) == 0) return(fit)
  }
  NULL
}

#i0  to n
roll_avg_vol <- function(x, i0, model) {
  n <- length(x)
  starts <- seq(i0, n, by = B)
  res <- vector("list", length(starts)); n_fail <- 0
  for (k in seq_along(starts)) {
    s <- starts[k]; e <- min(s + B - 1, n); m <- e - s
    fit <- fit_safe(model, x[(s - W + 1):e], m)
    if (is.null(fit)) {
      n_fail <- n_fail + 1
      res[[k]] <- data.frame(idx = s:e, fc = NA_real_)
    } else {
      sg <- as.matrix(sigma(ugarchforecast(fit, n.ahead = H, n.roll = m)))
      res[[k]] <- data.frame(idx = s:e, fc = sqrt(colMeans(sg^2)) * ann)
    }
    if (k %% 10 == 0) cat(sprintf("  %s: block %d / %d\n", model, k, length(starts)))
  }
  out <- do.call(rbind, res); attr(out, "n_fail") <- n_fail; out
}

# ---- HAC helpers (overlapping 21-day windows) ----
nw_se <- function(x, lag = H - 1) {
  x <- x - mean(x); n <- length(x)
  s <- sum(x^2) / n
  for (k in 1:lag)
    s <- s + 2 * (1 - k / (lag + 1)) * sum(x[(k + 1):n] * x[1:(n - k)]) / n
  sqrt(s / n)
}
dm_hac <- function(d) { z <- mean(d) / nw_se(d); c(stat = z, p = 2 * pnorm(-abs(z))) }

periods <- list(
  "2008-2009 (GFC)"         = c("2008-01-01", "2009-12-31"),
  "2010-2012 (euro crisis)" = c("2010-01-01", "2012-12-31"),
  "2013-2019"               = c("2013-01-01", "2019-12-31"),
  "2020 (COVID)"            = c("2020-01-01", "2020-12-31"),
  "2021-2026"               = c("2021-01-01", "2026-12-31")
)

analyse <- function(ret, label) {
  cat(sprintf("\n==================== %s ====================\n", label))
  x <- as.numeric(ret); d <- index(ret); n <- length(x)
  cat(sprintf("Returns: %s -> %s (N = %d, zero returns: %d)\n",
              d[1], d[n], n, sum(x == 0)))
  
  i0 <- max(which(d >= eval_start)[1], W)
  cat(sprintf("First forecast origin: %s (requested >= %s, window = %d days)\n",
              d[i0], eval_start, W))
  if (d[i0] > eval_start + 30)
    cat("NOTE: first origin is later than requested because history is too short.\n")
  
  cat("garch. \n");  g <- roll_avg_vol(x, i0, "sGARCH")
  cat("egarch.\n"); e <- roll_avg_vol(x, i0, "eGARCH")
  cat(sprintf("Blocks with no converged fit (forecasts set to NA): GARCH = %d | EGARCH = %d\n",
              attr(g, "n_fail"), attr(e, "n_fail")))
  
  df <- data.frame(date = d[g$idx], garch = g$fc, egarch = e$fc)
  df$rv_fwd  <- sapply(g$idx, function(i)
    if (i + H <= n) sqrt(mean(x[(i + 1):(i + H)]^2)) * ann else NA_real_)
  df$implied <- as.numeric(v2tx)[match(df$date, index(v2tx))]
  
  d2 <- df[complete.cases(df), ]
  cat(sprintf("Origins: %d | complete (forecasts, implied, forward realized): %d\n",
              nrow(df), nrow(d2)))

  rmse <- function(a, b) sqrt(mean((a - b)^2))
  dem  <- function(v) v - mean(v)
  fs <- list(GARCH = d2$garch, EGARCH = d2$egarch, VSTOXX = d2$implied)
  tab <- data.frame(
    corr_vstoxx   = c(cor(d2$garch, d2$implied), cor(d2$egarch, d2$implied), NA),
    mean_error    = sapply(fs, function(f) mean(f - d2$rv_fwd)),
    rmse          = sapply(fs, function(f) rmse(f, d2$rv_fwd)),
    corr_realized = sapply(fs, function(f) cor(f, d2$rv_fwd)),
    rmse_demeaned = sapply(fs, function(f) rmse(dem(f), dem(d2$rv_fwd))))
  cat("\nFull sample (vol points, annualized; mean error = forecast - realized):\n")
  print(round(tab, 3))
  
  vrp <- d2$implied - d2$rv_fwd
  cat(sprintf("\nVariance risk premium (VSTOXX - subsequent realized): mean = %+.2f (NW s.e. %.2f) | positive on %.0f%% of origins\n",
              mean(vrp), nw_se(vrp), 100 * mean(vrp > 0)))
#DB mariano with HAC
  loss  <- function(f) (f - d2$rv_fwd)^2
  lossd <- function(f) (dem(f) - dem(d2$rv_fwd))^2
  cat("\nDiebold-Mariano, squared error vs realized, HAC lag 20 (stat > 0 => second model better):\n")
  r1 <- dm_hac(loss(d2$garch) - loss(d2$egarch))
  cat(sprintf("  raw      GARCH vs EGARCH : stat = %+.2f, p = %.3f\n", r1[1], r1[2]))
  cat("Same after removing each series' mean level (isolates tracking from bias):\n")
  for (pr in list(c("garch", "egarch"), c("egarch", "implied"), c("garch", "implied"))) {
    r2 <- dm_hac(lossd(d2[[pr[1]]]) - lossd(d2[[pr[2]]]))
    cat(sprintf("  demeaned %-6s vs %-7s: stat = %+.2f, p = %.3f\n",
                pr[1], pr[2], r2[1], r2[2]))
  }
  
  #by regime
  per <- do.call(rbind, lapply(names(periods), function(nm) {
    p <- periods[[nm]]
    s <- d2[d2$date >= as.Date(p[1]) & d2$date <= as.Date(p[2]), ]
    if (nrow(s) < 30) return(NULL)
    data.frame(period = nm, n = nrow(s),
               vrp = mean(s$implied - s$rv_fwd),
               rmse_garch = rmse(s$garch, s$rv_fwd),
               rmse_egarch = rmse(s$egarch, s$rv_fwd),
               rmse_vstoxx = rmse(s$implied, s$rv_fwd),
               cor_garch = cor(s$garch, s$rv_fwd),
               cor_egarch = cor(s$egarch, s$rv_fwd),
               cor_vstoxx = cor(s$implied, s$rv_fwd))
  }))
  cat("\nBy regime (vol points):\n")
  print(per, row.names = FALSE, digits = 3)
  ragg::agg_png(sprintf("figures/ext_iv_vs_garch_%s.png", label),
                width = 1600, height = 700, units = "px", res = 110,
                background = "white")
  yl <- range(df[, c("garch", "egarch", "rv_fwd", "implied")], na.rm = TRUE)
  plot(df$date, df$implied, type = "l", col = "black", lwd = 1.6, ylim = yl,
       main = sprintf("%s: 21-day volatility forecasts vs VSTOXX (annualized, %%)", label),
       xlab = "", ylab = "%")
  lines(df$date, df$rv_fwd, col = "grey60",    lwd = 1.6)
  lines(df$date, df$garch,  col = "steelblue", lwd = 1.2)
  lines(df$date, df$egarch, col = "darkred",   lwd = 1.2)
  legend("topright",
         c("VSTOXX (implied)", "Subsequent realized (21d)", "GARCH", "EGARCH"),
         col = c("black", "grey60", "steelblue", "darkred"),
         lty = 1, lwd = c(1.6, 1.6, 1.2, 1.2), bty = "n")
  dev.off()
  
  write.csv(df,  sprintf("outputs/ext_iv_comparison_%s.csv", label), row.names = FALSE)
  write.csv(per, sprintf("outputs/ext_iv_by_regime_%s.csv",  label), row.names = FALSE)
  invisible(df)
}
res_cac <- analyse(returns_yahoo("^FCHI"), "CAC40")
r_stoxx <- tryCatch({
  if (file.exists("data/sx5e_history.csv")) {
    cat("\nUsing data/sx5e_history.csv for the EURO STOXX 50\n")
    returns_file("data/sx5e_history.csv")
  } else {
    returns_yahoo("^STOXX50E")
  }
}, error = function(err) { message("EURO STOXX 50 load failed: ", conditionMessage(err)); NULL })
if (!is.null(r_stoxx)) res_stoxx <- analyse(r_stoxx, "STOXX50E")