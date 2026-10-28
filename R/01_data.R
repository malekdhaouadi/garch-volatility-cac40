library(quantmod)
library(yaml)
library(xts)

cfg <- read_yaml("config.yml")

fetch_prices <- function(symbol, from, to) {
  stopifnot(is.character(symbol), length(symbol) == 1)
  px <- getSymbols(symbol, from = from, to = to, auto.assign = FALSE)
  px <- Ad(px)
  colnames(px) <- "adj_close"
  px
}

compute_log_returns <- function(prices) {
  r <- na.omit(diff(log(prices)))
  colnames(r) <- "ret"
  r
}

prices <- fetch_prices(cfg$index, cfg$start_date, cfg$end_date)

# Drop holiday NAs — do NOT interpolate: a non-trading day is not a 0% return
n_na_prices <- sum(is.na(prices))
cat("NA prices from Yahoo:", n_na_prices, "\n")
if (n_na_prices > 0) {
  cat("Dropping NA dates:", as.character(index(prices)[is.na(prices)]), "\n")
  prices <- prices[!is.na(prices)]
}

returns <- compute_log_returns(prices)

saveRDS(prices,  "data/prices.rds")
saveRDS(returns, "data/returns.rds")

cat("Observations:", nrow(returns), "\n")
cat("Date range:  ", as.character(start(returns)), "->",
    as.character(end(returns)), "\n")
cat("Mean return: ", round(mean(returns), 6), "\n")
cat("SD return:   ", round(sd(returns),   6), "\n")

stopifnot(
  nrow(returns) > 2000,
  sd(returns)  > 0.005,
  sd(returns)  < 0.03,
  !any(is.na(returns))
)

png("figures/returns_series.png", width = 1200, height = 500)
plot(returns, main = "CAC 40 daily log returns",
     ylab = "log return", col = "steelblue")
dev.off()