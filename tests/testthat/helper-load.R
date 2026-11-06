# Loaded automatically by testthat before the tests run.
suppressPackageStartupMessages({ library(testthat); library(xts) })

# testthat runs tests with the working directory set to tests/testthat
project_root <- normalizePath(file.path(getwd(), "..", ".."))
out_path <- function(f) file.path(project_root, "outputs", f)

for (f in list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

has_rugarch <- requireNamespace("rugarch", quietly = TRUE)
if (has_rugarch) suppressPackageStartupMessages(library(rugarch))

# Synthetic helpers
make_prices <- function(n = 600, start = as.Date("2020-01-01"), seed = 1) {
  set.seed(seed)
  xts(100 * exp(cumsum(rnorm(n, 0, 0.01))), start + seq_len(n) - 1)
}

sim_garch <- function(n, omega = 2e-6, alpha = 0.08, beta = 0.90, nu = 6, seed = 1) {
  set.seed(seed)
  N <- n + 500
  z <- rt(N, nu) * sqrt((nu - 2) / nu)
  h <- numeric(N); r <- numeric(N)
  h[1] <- omega / (1 - alpha - beta); r[1] <- sqrt(h[1]) * z[1]
  for (t in 2:N) { h[t] <- omega + alpha * r[t - 1]^2 + beta * h[t - 1]; r[t] <- sqrt(h[t]) * z[t] }
  tail(r, n)
}
