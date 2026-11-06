# These tests recompute published numbers from the saved outputs. If an analysis function
# changes behaviour, they fail. They are skipped when the output files are absent (e.g. in CI).

test_that("backtest functions reproduce the saved VaR backtest (07)", {
  f1 <- out_path("var_es_forecasts_CAC40.csv"); f2 <- out_path("var_backtest_CAC40.csv")
  skip_if_not(file.exists(f1) && file.exists(f2), "saved VaR outputs not found")
  fc <- read.csv(f1); bt <- read.csv(f2)
  for (i in seq_len(nrow(bt))) {
    m <- bt$model[i]; a <- bt$alpha[i]; col <- paste0(m, "_var", sprintf("%02d", round(a * 100)))
    v <- fc$ret < fc[[col]]
    ku <- kupiec(v, a); ci <- christoffersen_ind(v)
    expect_equal(sum(v), bt$exceed[i], info = paste(m, a))
    expect_equal(unname(ku["p"]), bt$kupiec_p[i], tolerance = 1e-6, info = paste(m, a))
    expect_equal(unname(ci["p"]), bt$indep_p[i], tolerance = 1e-6, info = paste(m, a))
    expect_equal(pchisq(unname(ku["LR"] + ci["LR"]), 2, lower.tail = FALSE), bt$cc_p[i],
                 tolerance = 1e-6, info = paste(m, a))
    expect_equal(1e4 * mean(tick(fc$ret, fc[[col]], a)), bt$tick_x1e4[i],
                 tolerance = 1e-6, info = paste(m, a))
  }
})

test_that("ols_hac reproduces the saved encompassing regression (09)", {
  f1 <- out_path("ext_iv_comparison_STOXX50E.csv"); f2 <- out_path("encompassing_regressions.csv")
  skip_if_not(file.exists(f1) && file.exists(f2), "saved encompassing outputs not found")
  d <- read.csv(f1); d <- d[complete.cases(d), ]; enc <- read.csv(f2)
  fit <- ols_hac(d$rv_fwd, as.matrix(d[, c("implied", "garch")]), lag = 20)
  row <- enc[enc$dataset == "STOXX50E" & enc$scale == "vol" &
             enc$spec == "4: implied + GARCH" & enc$term == "garch", ]
  expect_equal(fit$b[3], row$est, tolerance = 1e-6)
  expect_equal(fit$b[3] / fit$se[3], row$t20, tolerance = 1e-6)
})

test_that("variance risk premium matches the reported numbers (EURO STOXX 50)", {
  f1 <- out_path("ext_iv_comparison_STOXX50E.csv")
  skip_if_not(file.exists(f1), "saved comparison not found")
  d <- read.csv(f1); d <- d[complete.cases(d), ]
  v <- d$implied - d$rv_fwd
  expect_lt(abs(mean(v) - 3.23), 0.01)
  expect_lt(abs(nw_se(v, 20) - 0.40), 0.01)
  expect_lt(abs(100 * mean(v > 0) - 78), 1)
})
