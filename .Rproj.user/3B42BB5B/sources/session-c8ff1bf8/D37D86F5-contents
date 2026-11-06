library(rugarch)
library(yaml)
library(xts)

source("R/utils.R")
cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

returns <- readRDS("data/returns.rds")

make_spec <- function(model) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = cfg$distribution
  )
}

fit_garch  <- ugarchfit(make_spec("sGARCH"), data = returns, solver = "hybrid")
fit_egarch <- ugarchfit(make_spec("eGARCH"), data = returns, solver = "hybrid")

stopifnot(convergence(fit_garch) == 0, convergence(fit_egarch) == 0)

param_table <- function(fit) {
  m <- fit@fit$robust.matcoef
  data.frame(param = rownames(m), round(m, 5), row.names = NULL,
             check.names = FALSE)
}
tab_garch  <- param_table(fit_garch)
tab_egarch <- param_table(fit_egarch)

cat("\n=== GARCH(1,1) - Student-t ===\n");  print(tab_garch)
cat("\n=== EGARCH(1,1) - Student-t ===\n"); print(tab_egarch)

halflife <- function(p) log(0.5) / log(p)

cg <- coef(fit_garch)
pers_g <- unname(cg["alpha1"] + cg["beta1"])
lr_vol_g <- sqrt(unname(cg["omega"]) / (1 - pers_g)) * sqrt(252)

pers_e <- unname(coef(fit_egarch)["beta1"])

cat(sprintf("\nGARCH : alpha+beta = %.4f | half-life = %.1f days | long-run vol = %.1f%% annualized\n",
            pers_g, halflife(pers_g), 100 * lr_vol_g))
cat(sprintf("EGARCH: beta = %.4f | half-life = %.1f days\n",
            pers_e, halflife(pers_e)))

# info criteria
ic <- cbind(GARCH = infocriteria(fit_garch)[, 1],
            EGARCH = infocriteria(fit_egarch)[, 1])
cat("\nInformation criteria (per observation):\n"); print(round(ic, 5))
cat(sprintf("Log-likelihood: GARCH = %.1f | EGARCH = %.1f\n",
            likelihood(fit_garch), likelihood(fit_egarch)))

# residuals
resid_diag <- function(fit, name) {
  z <- as.numeric(residuals(fit, standardize = TRUE))
  cat(sprintf("\n--- %s standardized residuals ---\n", name))
  for (k in c(10, 20)) {
    lb1 <- Box.test(z,   lag = k, type = "Ljung-Box")
    lb2 <- Box.test(z^2, lag = k, type = "Ljung-Box")
    cat(sprintf("  lag %2d | LB(z): p = %.3f | LB(z^2): p = %.3f\n",
                k, lb1$p.value, lb2$p.value))
  }
  for (q in c(5, 10)) {
    a <- arch_lm_test(z, q)
    cat(sprintf("  ARCH-LM q = %2d: p = %.3f\n", q, a$p.value))
  }
  invisible(z)
}
z_g <- resid_diag(fit_garch,  "GARCH")
z_e <- resid_diag(fit_egarch, "EGARCH")
cat("\nGood result = p-values ABOVE 0.05 (no structure left).\n")
dates <- index(returns)
vol_g <- as.numeric(sigma(fit_garch))  * sqrt(252) * 100
vol_e <- as.numeric(sigma(fit_egarch)) * sqrt(252) * 100

ragg::agg_png("figures/cond_vol_garch_egarch.png",
              width = 1400, height = 600, units = "px", res = 110,
              background = "white")
plot(dates, vol_g, type = "l", col = "steelblue", lwd = 1.2,
     main = "CAC 40 conditional volatility (annualized, %)",
     xlab = "", ylab = "%")
lines(dates, vol_e, col = "darkred", lwd = 1.2)
legend("topright", c("GARCH(1,1)", "EGARCH(1,1)"),
       col = c("steelblue", "darkred"), lty = 1, bty = "n")
dev.off()

ragg::agg_png("figures/news_impact_egarch.png",
              width = 900, height = 600, units = "px", res = 110,
              background = "white")
ni <- newsimpact(fit_egarch)
plot(ni$zx, ni$zy, type = "l", lwd = 2, col = "darkred",
     main = "EGARCH news impact curve",
     xlab = "standardized shock z(t-1)", ylab = "variance response")
abline(v = 0, lty = 2, col = "grey50")
dev.off()

saveRDS(list(garch = fit_garch, egarch = fit_egarch), "outputs/fits.rds")
write.csv(tab_garch,  "outputs/params_garch.csv",  row.names = FALSE)
write.csv(tab_egarch, "outputs/params_egarch.csv", row.names = FALSE)