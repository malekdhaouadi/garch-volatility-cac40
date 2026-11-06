library(rugarch)
library(quantmod)
library(yaml)
library(xts)

source("R/data.R"); source("R/models.R"); source("R/forecast.R"); source("R/stats.R")

cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

W          <- cfg$rolling_window
B          <- cfg$refit_every
data_start <- if (is.null(cfg$data_start)) "2000-01-01" else cfg$data_start
eval_start <- as.Date(if (is.null(cfg$eval_start)) "2008-01-01" else cfg$eval_start)
alphas     <- c(0.01, 0.05)
hs_win     <- 250

# Two-digit tags keep the column names of the saved CSVs (GARCH_var01, ...)
tag2 <- function(a) sprintf("%02d", round(a * 100))

# Views of the library output with the column names used below
param_view <- function(df) {
  out <- data.frame(idx = df$idx, sigma = df$sigma)
  for (a in alphas) {
    out[[paste0("var", tag2(a))]] <- df[[paste0("pv", tag(a))]]
    out[[paste0("es",  tag2(a))]] <- df[[paste0("pe", tag(a))]]
  }
  out
}
hs_view <- function(df) {
  out <- data.frame(idx = df$idx, sigma = df$sigma)
  for (a in alphas) {
    out[[paste0("var", tag2(a))]] <- df[[paste0("v", tag(a))]]
    out[[paste0("es",  tag2(a))]] <- df[[paste0("e", tag(a))]]
  }
  out
}

run_backtest <- function(ret, label) {
  cat(sprintf("\n==================== %s ====================\n", label))
  x <- as.numeric(ret); d <- index(ret); n <- length(x)
  i0 <- max(which(d >= eval_start)[1], W + 1, hs_win + 1)
  cat(sprintf("Returns %s -> %s (N = %d). First target day: %s, %d out-of-sample days.\n",
              d[1], d[n], n, d[i0], n - i0 + 1))
  
  cat("garch\n")
  g <- param_view(roll_risk(x, i0, "sGARCH", W, B, alphas, dist = cfg$distribution, verbose = TRUE))
  cat("egarch\n")
  e <- param_view(roll_risk(x, i0, "eGARCH", W, B, alphas, dist = cfg$distribution, verbose = TRUE))
  h <- hs_view(hs_risk(x, i0, alphas, hs_win))
  cat(sprintf("Days without a converged fit: GARCH = %d | EGARCH = %d\n",
              sum(is.na(g$sigma)), sum(is.na(e$sigma))))
  
  mods <- list(GARCH = g, EGARCH = e, HS250 = h)
  ok <- complete.cases(g[, -1], e[, -1], h[, -1])
  mods <- lapply(mods, function(m) m[ok, ])
  r  <- x[g$idx[ok]]; dt <- d[g$idx[ok]]
  cat(sprintf("Days used (all three models available): %d\n", length(r)))
  
  all_rows <- list()
  for (a in alphas) {
    cat(sprintf("\n----- %.0f%% VaR / ES (expected exceedances: %.1f of %d) -----\n",
                100 * a, a * length(r), length(r)))
    rows <- lapply(names(mods), function(nm) {
      m <- mods[[nm]]; VaR <- m[[paste0("var", tag2(a))]]; ES <- m[[paste0("es", tag2(a))]]
      v <- r < VaR
      ku <- kupiec(v, a); ci <- christoffersen_ind(v)
      cc_p <- pchisq(ku["LR"] + ci["LR"], 2, lower.tail = FALSE)
      et <- es_test(r, VaR, ES, m$sigma)
      data.frame(model = nm, exceed = sum(v), rate_pct = 100 * mean(v),
                 kupiec_p = ku["p"], indep_p = ci["p"], cc_p = unname(cc_p),
                 es_mean_resid = et["mean_resid"], es_p_two = et["p_two"],
                 es_p_less = et["p_less"],
                 tick_x1e4 = 1e4 * mean(tick(r, VaR, a)), row.names = NULL)
    })
    tab <- do.call(rbind, rows); tab$alpha <- a
    print(format(tab[, -ncol(tab)], digits = 3), row.names = FALSE)
    all_rows[[as.character(a)]] <- tab
    
    cat("Diebold-Mariano on tick loss (stat > 0 => second model has lower loss):\n")
    L <- lapply(mods, function(m) tick(r, m[[paste0("var", tag2(a))]], a))
    for (pr in list(c("GARCH", "EGARCH"), c("GARCH", "HS250"), c("EGARCH", "HS250"))) {
      z <- dm_plain(L[[pr[1]]] - L[[pr[2]]])
      cat(sprintf("  %-6s vs %-6s: stat = %+.2f, p = %.3f\n", pr[1], pr[2], z[1], z[2]))
    }
  }
  
  yr <- format(dt, "%Y")
  ex <- sapply(names(mods), function(nm) tapply(r < mods[[nm]]$var01, yr, sum))
  ex <- cbind(days = as.integer(table(yr)), ex)
  cat("\nExceedances of the 1% VaR by calendar year (expected about 2.5 per 250 days):\n")
  print(ex)
  
  ragg::agg_png(sprintf("figures/var_exceedances_%s.png", label),
                width = 1600, height = 700, units = "px", res = 110, background = "white")
  v_e <- r < mods$EGARCH$var01
  plot(dt, 100 * r, type = "l", col = "grey65", xlab = "", ylab = "daily return (%)",
       main = sprintf("%s: 1%% 1-day VaR (EGARCH) and exceedances", label))
  lines(dt, 100 * mods$EGARCH$var01, col = "darkred", lwd = 1.2)
  lines(dt, 100 * mods$GARCH$var01,  col = "steelblue", lwd = 1.2)
  points(dt[v_e], 100 * r[v_e], pch = 19, col = "darkred", cex = 0.6)
  legend("bottomright", c("Return", "VaR 1% EGARCH", "VaR 1% GARCH", "EGARCH exceedance"),
         col = c("grey65", "darkred", "steelblue", "darkred"),
         lty = c(1, 1, 1, NA), pch = c(NA, NA, NA, 19), bty = "n")
  dev.off()
  
  ragg::agg_png(sprintf("figures/var_cum_exceedances_%s.png", label),
                width = 1200, height = 650, units = "px", res = 110, background = "white")
  cum <- sapply(names(mods), function(nm) cumsum(r < mods[[nm]]$var01))
  matplot(dt, cbind(cum, 0.01 * seq_along(r)), type = "l", lty = c(1, 1, 1, 2),
          col = c("steelblue", "darkred", "darkgreen", "black"), lwd = c(1.5, 1.5, 1.5, 2),
          xlab = "", ylab = "cumulative exceedances",
          main = sprintf("%s: cumulative 1%% VaR exceedances vs expected", label))
  legend("topleft", c(names(mods), "Expected (1%)"),
         col = c("steelblue", "darkred", "darkgreen", "black"),
         lty = c(1, 1, 1, 2), lwd = c(1.5, 1.5, 1.5, 2), bty = "n")
  dev.off()
  
  res <- do.call(rbind, all_rows)
  write.csv(res, sprintf("outputs/var_backtest_%s.csv", label), row.names = FALSE)
  out <- data.frame(date = dt, ret = r)
  for (nm in names(mods)) for (a in alphas) {
    out[[paste0(nm, "_var", tag2(a))]] <- mods[[nm]][[paste0("var", tag2(a))]]
    out[[paste0(nm, "_es",  tag2(a))]] <- mods[[nm]][[paste0("es",  tag2(a))]]
  }
  write.csv(out, sprintf("outputs/var_es_forecasts_%s.csv", label), row.names = FALSE)
  invisible(list(table = res, forecasts = out))
}

ret_cac <- get_returns("^FCHI", data_start, cfg$end_date)
res_cac <- run_backtest(ret_cac, "CAC40")