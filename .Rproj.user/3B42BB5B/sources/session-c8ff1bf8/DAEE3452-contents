library(rugarch)
library(yaml)
library(xts)
library(zoo)

source("R/utils.R")
cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

returns <- readRDS("data/returns.rds")
n <- length(returns)

make_spec <- function(model) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = cfg$distribution
  )
}
roll <- function(model) {
  ugarchroll(
    make_spec(model), data = returns,
    n.ahead         = 1,
    forecast.length = cfg$forecast_length,
    refit.every     = cfg$refit_every,
    refit.window    = "moving",
    window.size     = cfg$rolling_window,
    solver          = "hybrid",
    calculate.VaR   = FALSE,
    keep.coef       = TRUE
  )
}

cat("Rolling GARCH...\n");  roll_g <- roll("sGARCH")
cat("Rolling EGARCH...\n"); roll_e <- roll("eGARCH")

cat("Non-converged refits: GARCH =", length(roll_g@model$noncidx),
    "| EGARCH =", length(roll_e@model$noncidx), "\n")

fc_g <- as.data.frame(roll_g, which = "density")
fc_e <- as.data.frame(roll_e, which = "density")
stopifnot(nrow(fc_g) == cfg$forecast_length, nrow(fc_e) == cfg$forecast_length)

sig2_g <- fc_g$Sigma^2
sig2_e <- fc_e$Sigma^2
real   <- pmax((fc_g$Realized - fc_g$Mu)^2, 1e-10)  # noisy proxy; floor avoids log(0)

#loss
qlike <- function(real, fc) real / fc - log(real / fc) - 1
mse   <- function(real, fc) (real - fc)^2

L_g <- qlike(real, sig2_g)
L_e <- qlike(real, sig2_e)

cat(sprintf("\nMean QLIKE: GARCH = %.4f | EGARCH = %.4f\n", mean(L_g), mean(L_e)))
cat(sprintf("Mean MSE  : GARCH = %.3e | EGARCH = %.3e\n",
            mean(mse(real, sig2_g)), mean(mse(real, sig2_e))))

#db
dm_test <- function(l1, l2) {
  d    <- l1 - l2
  stat <- mean(d) / sqrt(var(d) / length(d))
  list(statistic = stat, p.value = 2 * pnorm(-abs(stat)))
}
dm_q <- dm_test(L_g, L_e)
dm_m <- dm_test(mse(real, sig2_g), mse(real, sig2_e))
cat(sprintf("Diebold-Mariano QLIKE: stat = %.3f, p = %.4f\n", dm_q$statistic, dm_q$p.value))
cat(sprintf("Diebold-Mariano MSE  : stat = %.3f, p = %.4f\n", dm_m$statistic, dm_m$p.value))
cat("Positive stat => GARCH has higher loss => EGARCH forecasts better.\n")
oos <- data.frame(date = as.Date(rownames(fc_g)), real = real,
                  garch = sig2_g, egarch = sig2_e)
saveRDS(oos, "outputs/oos_forecasts.rds")
write.csv(oos, "outputs/oos_forecasts.csv", row.names = FALSE)
ret_xts  <- returns
rv22     <- sqrt(rollapply(ret_xts^2, 22, mean, align = "right", fill = NA)) * sqrt(252) * 100
rv_oos   <- as.numeric(rv22[oos$date])

ragg::agg_png("figures/oos_forecast_vs_realized.png",
              width = 1400, height = 600, units = "px", res = 110,
              background = "white")
yl <- range(c(sqrt(oos$garch), sqrt(oos$egarch)) * sqrt(252) * 100, rv_oos, na.rm = TRUE)
plot(oos$date, rv_oos, type = "l", col = "grey60", lwd = 2, ylim = yl,
     main = "CAC 40: 1-day-ahead volatility forecasts vs 22-day realized (annualized, %)",
     xlab = "", ylab = "%")
lines(oos$date, sqrt(oos$garch)  * sqrt(252) * 100, col = "steelblue", lwd = 1.2)
lines(oos$date, sqrt(oos$egarch) * sqrt(252) * 100, col = "darkred",   lwd = 1.2)
legend("topright", c("Realized (22d)", "GARCH", "EGARCH"),
       col = c("grey60", "steelblue", "darkred"), lty = 1, lwd = c(2, 1.2, 1.2), bty = "n")
dev.off()
fits <- readRDS("outputs/fits.rds")
H <- cfg$forecast_horizon

fc30_g <- ugarchforecast(fits$garch,  n.ahead = H)
fc30_e <- ugarchforecast(fits$egarch, n.ahead = H)

path_g <- as.numeric(sigma(fc30_g)) * sqrt(252) * 100
path_e <- as.numeric(sigma(fc30_e)) * sqrt(252) * 100

#30day avergae vol
avg30 <- function(sig) sqrt(mean(sig^2)) * sqrt(252) * 100
avg_g <- avg30(as.numeric(sigma(fc30_g)))
avg_e <- avg30(as.numeric(sigma(fc30_e)))

cat(sprintf("\nForecast origin: %s\n", as.character(end(returns))))
cat(sprintf("Average 30-day vol (annualized): GARCH = %.2f%% | EGARCH = %.2f%%\n",
            avg_g, avg_e))
cat(sprintf("Day 1 vol: GARCH = %.2f%% | EGARCH = %.2f%%  ||  Day %d: GARCH = %.2f%% | EGARCH = %.2f%%\n",
            path_g[1], path_e[1], H, path_g[H], path_e[H]))

fc30 <- data.frame(horizon = 1:H, garch_vol = path_g, egarch_vol = path_e)
saveRDS(list(path = fc30, avg_garch = avg_g, avg_egarch = avg_e,
             origin = end(returns)), "outputs/forecast_30d.rds")
write.csv(fc30, "outputs/forecast_30d.csv", row.names = FALSE)

#last 120 days of fitted vol+30-day forecast
hist_n <- 120
hist_g <- tail(as.numeric(sigma(fits$garch)),  hist_n) * sqrt(252) * 100
hist_e <- tail(as.numeric(sigma(fits$egarch)), hist_n) * sqrt(252) * 100
x_hist <- -(hist_n - 1):0
x_fc   <- 1:H

ragg::agg_png("figures/forecast_30d.png",
              width = 1200, height = 600, units = "px", res = 110,
              background = "white")
yl <- range(hist_g, hist_e, path_g, path_e)
plot(x_hist, hist_g, type = "l", col = "steelblue", ylim = yl,
     xlim = c(-(hist_n - 1), H),
     main = "CAC 40: conditional volatility and 30-day forecast (annualized, %)",
     xlab = "trading days from forecast origin", ylab = "%")
lines(x_hist, hist_e, col = "darkred")
lines(x_fc, path_g, col = "steelblue", lty = 2, lwd = 2)
lines(x_fc, path_e, col = "darkred",   lty = 2, lwd = 2)
abline(v = 0, col = "grey50", lty = 3)
legend("topleft", c("GARCH", "EGARCH"), col = c("steelblue", "darkred"),
       lty = 1, bty = "n")
dev.off()