test_that("std_t_var is the unit-variance quantile used by rugarch", {
  skip_if_not(has_rugarch)
  for (nu in c(4.5, 5, 8)) for (a in c(0.01, 0.025, 0.05))
    expect_equal(std_t_var(a, nu),
                 qdist(distribution = "std", p = a, mu = 0, sigma = 1, shape = nu),
                 tolerance = 1e-8)
})

test_that("std_t_es equals the integral of the quantile function", {
  skip_if_not(has_rugarch)
  nu <- 5; a <- 0.01
  num <- integrate(function(u) qdist(distribution = "std", p = u, mu = 0, sigma = 1, shape = nu),
                   0, a)$value / a
  expect_equal(std_t_es(a, nu), num, tolerance = 1e-4)
  expect_lt(std_t_es(a, nu), std_t_var(a, nu))
  expect_lt(std_t_var(a, nu), 0)
})

test_that("dist_q_es: skewed-t with symmetry parameter 1 equals the Student-t", {
  skip_if_not(has_rugarch)
  a <- 0.025
  sym <- dist_q_es("sstd", a, c(shape = 6, skew = 1))
  std <- dist_q_es("std",  a, c(shape = 6))
  expect_equal(unname(sym), unname(std), tolerance = 1e-4)
})

test_that("hs_risk uses only past returns", {
  set.seed(9)
  x <- rnorm(600, 0, 0.01); i0 <- 300
  a <- hs_risk(x, i0)
  x2 <- x; x2[600] <- -0.5                      # change the last observation
  b <- hs_risk(x2, i0)
  expect_equal(a$v010, b$v010)                  # VaR for every target day, including the last, unchanged
  expect_true(all(a$e010 <= a$v010 + 1e-12))    # ES is below VaR
  expect_true(all(a$v010 <= a$v050 + 1e-12))    # 1% VaR is below 5% VaR
})

test_that("fit_safe recovers persistence on simulated GARCH data", {
  skip_if_not(has_rugarch)
  x <- sim_garch(2000)
  fit <- fit_safe("sGARCH", x)
  expect_false(is.null(fit))
  pers <- sum(coef(fit)[c("alpha1", "beta1")])
  expect_gt(pers, 0.85); expect_lt(pers, 1)
})

test_that("roll_risk has no look-ahead: the last return cannot affect any forecast", {
  skip_if_not(has_rugarch)
  x <- sim_garch(700, seed = 11)
  a <- roll_risk(x, i0 = 520, "sGARCH", window = 500, refit_every = 60, alphas = c(0.01, 0.05))
  x2 <- x; x2[700] <- x2[700] * 10
  b <- roll_risk(x2, i0 = 520, "sGARCH", window = 500, refit_every = 60, alphas = c(0.01, 0.05))
  expect_equal(nrow(a), 700 - 520 + 1)
  expect_false(anyNA(a$sigma))
  expect_equal(a$pv010, b$pv010, tolerance = 1e-8)
  expect_equal(a$fv050, b$fv050, tolerance = 1e-8)
  expect_true(all(a$pv010 < a$pv050))
})

test_that("roll_avg_vol has no look-ahead and returns plausible annualized volatility", {
  skip_if_not(has_rugarch)
  x <- sim_garch(700, seed = 12)
  a <- roll_avg_vol(x, i0 = 520, "sGARCH", window = 500, refit_every = 60, horizon = 5)
  x2 <- x; x2[700] <- x2[700] * 10
  b <- roll_avg_vol(x2, i0 = 520, "sGARCH", window = 500, refit_every = 60, horizon = 5)
  expect_equal(a$fc[a$idx < 700], b$fc[b$idx < 700], tolerance = 1e-8)   # origin 700 may use x[700]
  expect_true(all(a$fc > 1 & a$fc < 200))
  expect_equal(attr(a, "n_fail"), 0)
})
