library(yaml)
library(xts)
library(tseries)

cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

returns <- readRDS("data/returns.rds")
ret_vec <- as.numeric(returns)

#stats

s <- c(
  n          = length(ret_vec),
  mean       = mean(ret_vec),
  sd         = sd(ret_vec),
  skewness   = mean((ret_vec - mean(ret_vec))^3) / sd(ret_vec)^3,
  excess_kurt= mean((ret_vec - mean(ret_vec))^4) / sd(ret_vec)^4 - 3,
  min        = min(ret_vec),
  max        = max(ret_vec)
)
cat("Summary statistics \n")
print(round(s, 6))

#JB test
jb <- jarque.bera.test(ret_vec)
cat("\nJarque-Bera p-value:", format.pval(jb$p.value), "\n")
cat("Interpretation: p < 0.05 => reject normality => fat tails\n\n")


png("figures/returns_and_squared.png", width = 1400, height = 700)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
plot(returns, main = "CAC 40 log returns", ylab = "return", col = "steelblue")
plot(returns^2, main = "Squared returns (proxy for variance)",
     ylab = "squared return", col = "darkred")
dev.off()

# ljung box 
cat("Ljung-Box on squared returns \n")
cat("H0: no autocorrelation in squared returns\n")
for (k in c(10, 20, 50)) {
  lb <- Box.test(ret_vec^2, lag = k, type = "Ljung-Box")
  cat(sprintf("  lag %2d: Q = %8.2f, p = %s\n",
              k, lb$statistic, format.pval(lb$p.value)))
}

# ARCH LM test

arch_lm_test <- function(x, q) {
  # Regress x^2 on q lags of x^2, get LM = n*R^2
  y  <- x^2
  n  <- length(y)
  X  <- embed(y, q + 1)                     # first col = y_t, rest = lags
  yt <- X[, 1]
  Xr <- cbind(1, X[, -1])
  fit <- lm(yt ~ Xr - 1)
  r2  <- summary(fit)$r.squared
  lm_stat <- (n - q) * r2
  pval <- pchisq(lm_stat, df = q, lower.tail = FALSE)
  list(statistic = lm_stat, df = q, p.value = pval)
}

cat("\n ARCH-LM test\n")
cat("H0: no ARCH effects up to lag q\n")
for (q in c(5, 10)) {
  lm_q <- arch_lm_test(ret_vec, q)
  cat(sprintf("  q = %2d: LM = %8.2f, df = %d, p = %s\n",
              q, lm_q$statistic, lm_q$df, format.pval(lm_q$p.value)))
}

png("figures/qq_returns.png", width = 700, height = 700)
qqnorm(ret_vec, main = "QQ plot: returns vs Normal",
       pch = 20, cex = 0.5, col = "steelblue")
qqline(ret_vec, col = "red", lwd = 2)
dev.off()

diag <- list(
  summary       = s,
  jb_pvalue     = jb$p.value,
  lb_lags       = c(10, 20, 50),
  arch_lags     = c(5, 10)
)
saveRDS(diag, "outputs/diagnostics.rds")