arch_lm_test <- function(x, q) {
  y   <- x^2
  X   <- embed(y, q + 1)
  fit <- lm(X[, 1] ~ X[, -1])
  lm_stat <- nrow(X) * summary(fit)$r.squared
  list(statistic = lm_stat, df = q,
       p.value = pchisq(lm_stat, df = q, lower.tail = FALSE))
}