# R/01_data.R
library(quantmod)
library(yaml)
library(xts)

cfg <- read_yaml("config.yml")

fetch_prices <- function(symbol, from, to) {
  stopifnot(is.character(symbol), length(symbol) == 1)
  px <- getSymbols(symbol, from = from, to = to, auto.assign = FALSE)
  px <- Ad(px)                       # adjusted close only
  colnames(px) <- "adj_close"
  px
}

compute_log_returns <- function(prices) {
  r <- na.omit(diff(log(prices)))
  colnames(r) <- "ret"
  r
}

prices  <- fetch_prices(cfg$index, cfg$start_date, cfg$end_date)
returns <- compute_log_returns(prices)

saveRDS(prices,  "data/prices.rds")
saveRDS(returns, "data/returns.rds")

# Quick sanity check
cat("Observations:", nrow(returns), "\n")
cat("Date range:  ", as.character(start(returns)), "->",
    as.character(end(returns)), "\n")
cat("Mean return: ", round(mean(returns), 6), "\n")
cat("SD return:   ", round(sd(returns),   6), "\n")

# First plot: returns over time
png("figures/returns_series.png", width = 1200, height = 500)
plot(returns, main = "CAC 40 daily log returns",
     ylab = "log return", col = "steelblue")
dev.off()