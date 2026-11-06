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
alphas     <- c(0.01, 0.025, 0.05)
hs_win     <- 250
split_date <- as.Date("2016-12-31")     # selection: <= this date, hold-out: after

# ---- Experiment grid (fixed in advance; every cell is reported) ----
grid <- expand.grid(model = c("sGARCH", "eGARCH"),
                    mean0 = c(FALSE, TRUE),          # TRUE = zero-mean model
                    dist  = c("std", "sstd"),        # Student-t, skewed Student-t
                    Wc    = c(W, 2 * W),             # estimation window
                    stringsAsFactors = FALSE)
# Selection rule (fixed in advance): smallest mean |log(exceedance rate / alpha)|
# at the 1% and 5% levels on the SELECTION period. The hold-out period is not used.
sel_alphas <- c(0.01, 0.05)

dir.create("outputs/var_cache", recursive = TRUE, showWarnings = FALSE)

eval_cand <- function(cd, a, r, mask) {
  t <- tag(a); VaR <- cd$var[[t]][mask]; ES <- cd$es[[t]][mask]
  sg <- cd$sigma[mask]; rr <- r[mask]; v <- rr < VaR
  ku <- kupiec(v, a); ci <- christoffersen_ind(v); et <- es_test(rr, VaR, ES, sg)
  data.frame(alpha = a, n = length(rr), exceed = sum(v), rate_pct = 100 * mean(v),
             kupiec_p = unname(ku["p"]), indep_p = unname(ci["p"]),
             cc_p = pchisq(unname(ku["LR"] + ci["LR"]), 2, lower.tail = FALSE),
             dq_p = dq_test(v, VaR, a),
             es_resid = unname(et["mean_resid"]), es_p_less = unname(et["p_less"]),
             tick_x1e4 = 1e4 * mean(tick(rr, VaR, a)), row.names = NULL)
}

# =========================== Main ===========================
ret <- get_returns("^FCHI", data_start, cfg$end_date)
saveRDS(ret, "data/returns_CAC40_long.rds")
x <- as.numeric(ret); d <- index(ret); n <- length(x)
i0 <- max(which(d >= eval_start)[1], 2 * W + 1, hs_win + 1)     # same days for every candidate
cat(sprintf("Returns %s -> %s (N = %d). First target: %s | %d out-of-sample days\n",
            d[1], d[n], n, d[i0], n - i0 + 1))
cat(sprintf("Grid: %d configurations x 2 tail methods (+ HS%d benchmark)\n", nrow(grid), hs_win))
cat("Forecasts are cached in outputs/var_cache/; an interrupted run resumes where it stopped.\n")

cfg_dfs <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  id <- sprintf("%s_%s_%s_W%d", g$model, ifelse(g$mean0, "m0", "mC"), g$dist, g$Wc)
  f  <- sprintf("outputs/var_cache/%s_i%d.rds", id, i0)
  if (file.exists(f)) {
    cfg_dfs[[id]] <- readRDS(f); cat(sprintf("[%2d/%d] %s (cached)\n", i, nrow(grid), id))
  } else {
    t0 <- Sys.time()
    cfg_dfs[[id]] <- roll_risk(x, i0, g$model, window = g$Wc, refit_every = B,
                               alphas = alphas, mean0 = g$mean0, dist = g$dist)
    saveRDS(cfg_dfs[[id]], f)
    cat(sprintf("[%2d/%d] %s done in %.1f min\n", i, nrow(grid), id,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
hs <- hs_risk(x, i0, alphas, hs_win)

# Common days: drop any day where some configuration has no converged fit
ok <- Reduce(`&`, lapply(cfg_dfs, complete.cases))
cat(sprintf("Days used: %d of %d (days lost to non-converged blocks: %d)\n",
            sum(ok), length(ok), sum(!ok)))
r  <- x[i0:n][ok]; dt <- d[i0:n][ok]

mk <- function(df, pv, pe) list(
  sigma = df$sigma[ok],
  var = setNames(lapply(alphas, function(a) df[[paste0(pv, tag(a))]][ok]), tag(alphas)),
  es  = setNames(lapply(alphas, function(a) df[[paste0(pe, tag(a))]][ok]), tag(alphas)))
cands <- list()
for (id in names(cfg_dfs)) {
  cands[[paste0(id, "|param")]] <- mk(cfg_dfs[[id]], "pv", "pe")
  cands[[paste0(id, "|FHS")]]   <- mk(cfg_dfs[[id]], "fv", "fe")
}
cands[["HS250"]] <- mk(hs, "v", "e")
baseline <- "sGARCH_mC_std_W1000|param"          # the 07_var_backtest.R specification

sel  <- dt <= split_date
hold <- !sel
cat(sprintf("Selection period: %s -> %s (%d days) | hold-out: %s -> %s (%d days)\n",
            dt[1], max(dt[sel]), sum(sel), min(dt[hold]), dt[length(dt)], sum(hold)))

# ---- Selection on 2008-2016 only ----
score <- sapply(names(cands), function(nm) {
  mean(sapply(sel_alphas, function(a) {
    rate <- max(mean(r[sel] < cands[[nm]]$var[[tag(a)]][sel]), 1 / sum(sel))
    abs(log(rate / a)) }))
})
sel_tab <- do.call(rbind, lapply(names(cands), function(nm) {
  e1 <- eval_cand(cands[[nm]], 0.01, r, sel); e5 <- eval_cand(cands[[nm]], 0.05, r, sel)
  data.frame(candidate = nm, score = score[nm], rate1 = e1$rate_pct, kupiec1 = e1$kupiec_p,
             rate5 = e5$rate_pct, kupiec5 = e5$kupiec_p, tick1 = e1$tick_x1e4,
             row.names = NULL)
}))
sel_tab <- sel_tab[order(sel_tab$score), ]
sel_tab$rank <- seq_len(nrow(sel_tab))
cat("\n=== Selection period (2008-2016): top 12 candidates by calibration score ===\n")
print(format(head(sel_tab, 12), digits = 3), row.names = FALSE)
cat("\nBaseline and benchmark on the selection period:\n")
print(format(sel_tab[sel_tab$candidate %in% c(baseline, "HS250"), ], digits = 3), row.names = FALSE)
best <- sel_tab$candidate[1]
cat(sprintf("\nSelected on 2008-2016: %s\n", best))

# ---- Hold-out evaluation (2017-2026): selected, runner-ups, baseline, benchmark ----
show <- unique(c(best, sel_tab$candidate[2:3], baseline, "HS250"))
cat("\n=== HOLD-OUT (2017-2026): not used for selection ===\n")
for (a in alphas) {
  cat(sprintf("\n--- %.1f%% VaR / ES (expected exceedances: %.1f of %d) ---\n",
              100 * a, a * sum(hold), sum(hold)))
  tab <- do.call(rbind, lapply(show, function(nm) cbind(candidate = nm, eval_cand(cands[[nm]], a, r, hold))))
  print(format(tab[, -2], digits = 3), row.names = FALSE)
}

# ---- Everything, both periods (for the README and for transparency) ----
all_res <- do.call(rbind, lapply(names(cands), function(nm) rbind(
  cbind(candidate = nm, period = "selection 2008-2016", do.call(rbind, lapply(alphas, function(a) eval_cand(cands[[nm]], a, r, sel)))),
  cbind(candidate = nm, period = "hold-out 2017-2026", do.call(rbind, lapply(alphas, function(a) eval_cand(cands[[nm]], a, r, hold)))))))
write.csv(all_res, "outputs/var_experiments_all.csv", row.names = FALSE)
write.csv(sel_tab, "outputs/var_experiments_selection.csv", row.names = FALSE)

# ---- Exceedances per calendar year (1% VaR) ----
yr <- format(dt, "%Y")
ex <- sapply(show, function(nm) tapply(r < cands[[nm]]$var[["010"]], yr, sum))
ex <- cbind(days = as.integer(table(yr)), ex)
cat("\nExceedances of the 1% VaR by calendar year (expected about 2.5 per 250 days):\n")
print(ex)

# ---- Figure: cumulative 1% exceedances ----
ragg::agg_png("figures/var_experiments_cum_exceedances.png",
              width = 1300, height = 650, units = "px", res = 110, background = "white")
cols <- c("black", "darkred", "steelblue", "darkgreen")
sh <- c(best, baseline, "HS250")
cum <- sapply(sh, function(nm) cumsum(r < cands[[nm]]$var[["010"]]))
matplot(dt, cbind(0.01 * seq_along(r), cum), type = "l", lty = c(2, 1, 1, 1), lwd = 1.8,
        col = cols, xlab = "", ylab = "cumulative exceedances",
        main = "CAC 40: cumulative 1% VaR exceedances (dashed line = expected)")
abline(v = split_date, col = "grey50", lty = 3)
legend("topleft", c("Expected (1%)", paste("Selected:", best), paste("Baseline:", baseline), "HS250"),
       col = cols, lty = c(2, 1, 1, 1), lwd = 1.8, bty = "n", cex = 0.8)
dev.off()