test_that("compute_log_returns is exact and keeps column name", {
  p <- xts(c(100, 110, 99), as.Date("2024-01-01") + 0:2)
  r <- compute_log_returns(p)
  expect_equal(as.numeric(r), c(log(110 / 100), log(99 / 110)))
  expect_equal(colnames(r), "ret")
})

test_that("drop_na_prices removes NA rows and never interpolates", {
  p <- xts(c(100, 101, NA, 103, 104), as.Date("2024-01-01") + 0:4)
  q <- drop_na_prices(p)
  expect_equal(length(q), 4)
  expect_false(anyNA(q))
  expect_equal(attr(q, "dropped"), "2024-01-03")
  expect_equal(as.numeric(q), c(100, 101, 103, 104))
})

test_that("validate_prices passes clean data", {
  res <- validate_prices(make_prices())
  expect_true(all(res$ok))
  expect_silent(assert_valid_prices(make_prices()))
})

test_that("validate_prices catches each kind of problem", {
  p <- make_prices()
  ok_of <- function(res, name) res$ok[res$check == name]

  bad <- p; bad[10] <- NA
  expect_false(ok_of(validate_prices(bad), "no missing prices"))

  bad <- p; bad[10] <- -5
  expect_false(ok_of(validate_prices(bad), "prices strictly positive"))

  d <- as.Date("2020-01-01") + c(0:299, 299, 300:599)
  dup <- suppressWarnings(xts(rep(100, 601) + seq_len(601) / 100, d))
  expect_false(ok_of(validate_prices(dup), "dates strictly increasing"))

  gap <- p[-(200:240)]
  expect_false(ok_of(validate_prices(gap), "no long gaps"))

  jump <- p; jump[300] <- as.numeric(jump[299]) * 2
  expect_false(ok_of(validate_prices(jump), "no extreme returns"))

  expect_false(ok_of(validate_prices(p[1:100]), "enough observations"))

  flat <- xts(rep(100, 600), as.Date("2020-01-01") + 0:599)
  expect_false(ok_of(validate_prices(flat), "few zero returns"))
})

test_that("staleness check uses the supplied date", {
  p <- make_prices()
  last <- end(p)
  expect_true(all(validate_prices(p, max_staleness_days = 5, today = last + 2)$ok))
  expect_false(all(validate_prices(p, max_staleness_days = 5, today = last + 30)$ok))
})

test_that("assert_valid_prices stops with a readable message", {
  bad <- make_prices(); bad[10] <- NA
  expect_error(assert_valid_prices(bad), "Price validation failed")
})
