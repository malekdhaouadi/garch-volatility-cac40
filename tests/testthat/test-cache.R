# Proves that roll_risk() (shared library) reproduces the forecasts cached by the ORIGINAL
# implementation of script 08. Needs outputs/var_cache/ and data/returns_CAC40_long.rds
# (written by analysis/06 or 08); skipped when they are absent.

test_that("roll_risk reproduces the cached forecasts of the original implementation", {
  files <- list.files(file.path(project_root, "outputs", "var_cache"),
                      pattern = "^sGARCH_mC_std_W1000_i[0-9]+\\.rds$", full.names = TRUE)
  f_ret <- file.path(project_root, "data", "returns_CAC40_long.rds")
  skip_if_not(length(files) == 1 && file.exists(f_ret) && has_rugarch,
              "cache or long return series not found")
  old <- readRDS(files)
  i0  <- as.integer(sub("^.*_i([0-9]+)\\.rds$", "\\1", basename(files)))
  x   <- as.numeric(readRDS(f_ret))
  skip_if_not(length(x) >= i0 + 99, "return series shorter than the cache")

  # two refit blocks (100 target days) are enough to prove equivalence
  new <- roll_risk(x[seq_len(i0 + 99)], i0, "sGARCH", window = 1000, refit_every = 50,
                   alphas = c(0.01, 0.025, 0.05))
  expect_equal(unname(as.matrix(new)), unname(as.matrix(old[1:100, ])), tolerance = 1e-6)
})
