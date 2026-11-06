test_that("arch_lm_test detects ARCH effects and not iid noise", {
  set.seed(1)
  iid <- rnorm(3000)
  expect_gt(arch_lm_test(iid, 5)$p.value, 0.001)

  e <- rnorm(3000); x <- numeric(3000)
  for (t in 2:3000) x[t] <- sqrt(0.1 + 0.7 * x[t - 1]^2) * e[t]
  expect_lt(arch_lm_test(x, 5)$p.value, 1e-6)
})

test_that("arch_lm_test matches FinTS::ArchTest", {
  skip_if_not_installed("FinTS")
  set.seed(3)
  x <- rnorm(1500) * (1 + 0.5 * sin(seq_len(1500) / 40))
  ours <- arch_lm_test(x, 5)
  ref  <- FinTS::ArchTest(x, lags = 5)
  expect_equal(unname(ours$statistic), unname(ref$statistic), tolerance = 1e-6)
  expect_equal(ours$p.value, unname(ref$p.value), tolerance = 1e-6)
})

test_that("kupiec: zero statistic at the nominal rate, rejection when far off", {
  v <- rep(FALSE, 1000); v[1:10] <- TRUE
  k <- kupiec(v, 0.01)
  expect_equal(unname(k["LR"]), 0, tolerance = 1e-10)
  expect_equal(unname(k["p"]), 1, tolerance = 1e-10)

  v2 <- rep(FALSE, 1000); v2[1:30] <- TRUE
  expect_lt(unname(kupiec(v2, 0.01)["p"]), 1e-4)

  expect_true(is.finite(kupiec(rep(FALSE, 500), 0.01)["LR"]))   # no exceedances at all
})

test_that("christoffersen_ind: independent hits pass, clustered hits fail", {
  set.seed(2)
  expect_gt(unname(christoffersen_ind(runif(5000) < 0.05)["p"]), 0.001)
  clustered <- rep(c(rep(TRUE, 10), rep(FALSE, 990)), 5)
  expect_lt(unname(christoffersen_ind(clustered)["p"]), 1e-6)
})

test_that("dq_test: independent hits pass, autocorrelated hits fail", {
  set.seed(4)
  VaR <- rnorm(4000)
  expect_gt(dq_test(runif(4000) < 0.05, VaR, 0.05), 0.001)
  clustered <- rep(c(rep(TRUE, 12), rep(FALSE, 228)), 16)
  expect_lt(dq_test(clustered, VaR[seq_along(clustered)], 0.05), 1e-6)
  expect_true(is.na(dq_test(rep(FALSE, 10), rnorm(10), 0.05)))      # too short
})

test_that("es_test: correct ES passes, ES that ignores the tail fails", {
  set.seed(5)
  r <- rnorm(60000); a <- 0.05
  VaR <- rep(qnorm(a), length(r)); ES <- rep(-dnorm(qnorm(a)) / a, length(r)); sg <- rep(1, length(r))
  expect_gt(es_test(r, VaR, ES, sg)["p_less"], 0.001)
  expect_lt(es_test(r, VaR, VaR, sg)["p_less"], 1e-6)                # ES set equal to VaR
})

test_that("tick and qlike behave as defined", {
  expect_equal(tick(-0.03, -0.02, 0.01), 0.0099, tolerance = 1e-12)   # exceedance
  expect_equal(tick(0.01, -0.02, 0.01), 0.0003, tolerance = 1e-12)    # no exceedance
  expect_equal(qlike(2, 2), 0)
  expect_gt(qlike(2, 1), 0); expect_gt(qlike(2, 4), 0)
})

test_that("nw_se exceeds the naive standard error for autocorrelated data", {
  set.seed(6)
  x <- as.numeric(arima.sim(list(ar = 0.9), 5000))
  expect_gt(nw_se(x, 20), 2 * sd(x) / sqrt(length(x)))
  y <- rnorm(5000)
  expect_equal(nw_se(y, 5), sd(y) / sqrt(5000), tolerance = 0.2)
  expect_error(nw_se(1:10, 20))
})

test_that("ols_hac reproduces lm coefficients and gives sensible errors", {
  set.seed(7)
  X <- cbind(rnorm(2000), rnorm(2000)); y <- 1 + 2 * X[, 1] - X[, 2] + rnorm(2000)
  f <- ols_hac(y, X, lag = 5); l <- lm(y ~ X)
  expect_equal(f$b, unname(coef(l)), tolerance = 1e-8)
  expect_equal(f$se, unname(summary(l)$coefficients[, 2]), tolerance = 0.2)
  expect_equal(f$r2, summary(l)$r.squared, tolerance = 1e-8)
})

test_that("dm_hac sign convention: positive when the second model has lower loss", {
  set.seed(8)
  d <- rnorm(2000, mean = 0.3)            # L1 - L2 > 0 on average
  expect_gt(dm_hac(d)["stat"], 0)
  expect_lt(dm_hac(d)["p"], 0.001)
})
