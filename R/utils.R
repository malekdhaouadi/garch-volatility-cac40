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