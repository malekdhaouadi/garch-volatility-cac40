library(yaml); library(xts); library(quantmod)
source("R/stats.R"); source("R/data.R"); source("R/risk.R")

cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

alphas     <- cfg$var_alphas
eval_start <- as.Date(cfg$eval_start)
split_date <- as.Date("2016-12-31")
k_grid     <- seq(cfg$iv_var_scale_range[1], cfg$iv_var_scale_range[2], by = 0.01)
ret <- get_returns("^FCHI", cfg$data_start, cfg$end_date)
v2x <- read_stoxx_file("data/v2tx_history.csv")          
idx <- index(ret)
v   <- as.numeric(v2x)[match(idx, index(v2x))]            
v   <- na.locf(v, na.rm = FALSE)                        

keep <- !is.na(v) & as.Date(idx) >= eval_start
r  <- as.numeric(ret)[keep]
v  <- v[keep]
d  <- idx[keep]
n  <- length(r)

cat(sprintf("Days: %s -> %s (N = %d)\n", d[1], d[n], n))
fits <- readRDS("outputs/fits.rds")
nu   <- as.numeric(coef(fits$garch)["shape"])

var_iv   <- function(k, alpha) {
  # per-day 1-day VaR at level alpha for a position of 1 unit, mean 0
  sigma_d <- k * v / 100 / sqrt(252)        
  mu_d    <- 0
  mu_d + sigma_d * std_t_var(alpha, nu)
}

sel <- d <= split_date
score <- sapply(k_grid, function(k) {
  mean(sapply(c(0.01, 0.05), function(a) {
    VaR <- var_iv(k, a)[sel]
    rate <- max(mean(r[sel] < VaR), 1 / sum(sel))
    abs(log(rate / a))
  }))
})
best_k <- k_grid[which.min(score)]
cat(sprintf("Selected k on 2008-2016: %.2f (score %.4f)\n", best_k, min(score)))
cat(sprintf("Encompassing-regression k: %.2f\n", cfg$iv_var_scaling))

eval_k <- function(k) {
  out <- list()
  for (a in alphas) {
    VaR <- var_iv(k, a)
    v_hit <- r < VaR
    ku <- kupiec(v_hit, a)
    ci <- christoffersen_ind(v_hit)
    dq <- dq_test(v_hit, VaR, a)
    out[[as.character(a)]] <- c(
      rate_pct  = 100 * mean(v_hit),
      kupiec_p  = unname(ku["p"]),
      indep_p   = unname(ci["p"]),
      dq_p      = dq,
      tick_x1e4 = 1e4 * mean(tick(r, VaR, a))
    )
  }
  out
}

cat("\n Hold-out 2017-2026 \n")
for (a in alphas) {
  cat(sprintf("\n  alpha = %.1f%%\n", 100 * a))
  cat(sprintf("    %-20s %8s %10s %10s %10s %12s\n",
              "model", "rate%", "kupiec_p", "indep_p", "dq_p", "tick_x1e4"))
  for (nm in c(sprintf("IV VaR k=%.2f (sel)", best_k),
               sprintf("IV VaR k=%.2f (enf)", cfg$iv_var_scaling),
               sprintf("IV VaR k=1.00 (raw)", 1.0))) {
    k <- switch(nm,
                sprintf("IV VaR k=%.2f (sel)", best_k) = best_k,
                sprintf("IV VaR k=%.2f (enf)", cfg$iv_var_scaling) = cfg$iv_var_scaling,
                sprintf("IV VaR k=1.00 (raw)", 1.0) = 1.0)
    res <- eval_k(k)[[as.character(a)]]
    cat(sprintf("    %-20s %8.2f %10.3g %10.3g %10.3g %12.2f\n",
                nm, res["rate_pct"], res["kupiec_p"], res["indep_p"],
                res["dq_p"], res["tick_x1e4"]))
  }
}

iv_var_sel <- lapply(alphas, function(a) var_iv(best_k, a))
names(iv_var_sel) <- paste0("var", sprintf("%02d", round(100 * alphas)))
iv_es_sel  <- lapply(alphas, function(a) {
  sigma_d <- best_k * v / 100 / sqrt(252)
  mu_d <- 0
  mu_d + sigma_d * std_t_es(a, nu)
})
names(iv_es_sel) <- paste0("es", sprintf("%02d", round(100 * alphas)))

out <- data.frame(date = d, ret = r, vstoxx = v,
                  var01_iv = iv_var_sel$var01,
                  var05_iv = iv_var_sel$var05,
                  es01_iv  = iv_es_sel$es01,
                  es05_iv  = iv_es_sel$es05)
saveRDS(list(k_selected = best_k, k_encompassing = cfg$iv_var_scaling,
             nu = nu, daily = out), "outputs/iv_var.rds")
write.csv(out, "outputs/iv_var_daily.csv", row.names = FALSE)

cat(sprintf("\nSaved outputs/iv_var.rds (k_selected = %.2f)\n", best_k))