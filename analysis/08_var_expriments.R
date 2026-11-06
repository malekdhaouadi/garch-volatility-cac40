library(rugarch)
library(quantmod)
library(yaml)
library(xts)

cfg <- read_yaml("config.yml")
set.seed(cfg$seed)

W          <- cfg$rolling_window
B          <- cfg$refit_every
data_start <- if (is.null(cfg$data_start)) "2000-01-01" else cfg$data_start
eval_start <- as.Date(if (is.null(cfg$eval_start)) "2008-01-01" else cfg$eval_start)
alphas     <- c(0.01, 0.025, 0.05)
hs_win     <- 250
split_date <- as.Date("2016-12-31")
grid <- expand.grid(model = c("sGARCH", "eGARCH"),
                    mean0 = c(FALSE, TRUE),         
                    dist  = c("std", "sstd"),        
                    Wc    = c(W, 2 * W),             
                    stringsAsFactors = FALSE)
sel_alphas <- c(0.01, 0.05)

dir.create("outputs/var_cache", recursive = TRUE, showWarnings = FALSE)

returns_yahoo <- function(symbol) {
  px <- getSymbols(symbol, from = data_start, to = cfg$end_date, auto.assign = FALSE)
  px <- Ad(px); px <- px[!is.na(px)]
  r  <- na.omit(diff(log(px))); colnames(r) <- "ret"; r
}

tag <- function(a) sprintf("%03d", round(a * 1000))

make_spec <- function(model, mean0, dist) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = !mean0),
    distribution.model = dist
  )
}

fit_safe <- function(model, mean0, dist, dat, m) {
  for (sv in c("hybrid", "nlminb", "solnp")) {
    fit <- tryCatch(suppressWarnings(
      ugarchfit(make_spec(model, mean0, dist), data = dat, out.sample = m, solver = sv)),
      error = function(e) NULL)
    if (!is.null(fit) && convergence(fit) == 0) return(fit)
  }
  NULL
}

#param var
dist_q_es <- function(dist, a, cf) {
  nu <- unname(cf["shape"])
  if (dist == "std") {
    ta <- qt(a, nu); sc <- sqrt((nu - 2) / nu)
    return(c(q = ta * sc, es = -dt(ta, nu) / a * (nu + ta^2) / (nu - 1) * sc))
  }
  xi <- unname(cf["skew"])
  qf <- function(u) qdist("sstd", u, mu = 0, sigma = 1, skew = xi, shape = nu)
  es <- tryCatch(integrate(qf, 0, a, rel.tol = 1e-6)$value / a,
                 error = function(e) mean(qf((seq_len(2000) - 0.5) / 2000 * a)))
  c(q = qf(a), es = es)
}
risk_config <- function(x, i0, model, mean0, dist, Wc) {
  n <- length(x); starts <- seq(i0, n, by = B)
  out <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    s <- starts[k]; e <- min(s + B - 1, n); m <- e - s
    fit <- fit_safe(model, mean0, dist, x[(s - Wc):(e - 1)], m)
    row <- data.frame(idx = s:e, mu = NA_real_, sigma = NA_real_)
    for (a in alphas) for (p in c("pv", "pe", "fv", "fe")) row[[paste0(p, tag(a))]] <- NA_real_
    if (!is.null(fit)) {
      fc <- ugarchforecast(fit, n.ahead = 1, n.roll = m)
      row$mu <- as.numeric(fitted(fc)); row$sigma <- as.numeric(sigma(fc))
      cf <- coef(fit); z <- as.numeric(residuals(fit, standardize = TRUE))
      for (a in alphas) {
        qe <- dist_q_es(dist, a, cf)
        qz <- quantile(z, a, type = 7, names = FALSE); ez <- mean(z[z <= qz])
        row[[paste0("pv", tag(a))]] <- row$mu + row$sigma * qe["q"]
        row[[paste0("pe", tag(a))]] <- row$mu + row$sigma * qe["es"]
        row[[paste0("fv", tag(a))]] <- row$mu + row$sigma * qz
        row[[paste0("fe", tag(a))]] <- row$mu + row$sigma * ez
      }
    }
    out[[k]] <- row
  }
  do.call(rbind, out)
}

hs_risk <- function(x, i0) {
  n <- length(x); idx <- i0:n
  df <- data.frame(idx = idx, sigma = NA_real_)
  for (a in alphas) { df[[paste0("v", tag(a))]] <- NA_real_; df[[paste0("e", tag(a))]] <- NA_real_ }
  for (j in seq_along(idx)) {
    win <- x[(idx[j] - hs_win):(idx[j] - 1)]
    df$sigma[j] <- sd(win)
    for (a in alphas) {
      q <- quantile(win, a, type = 1, names = FALSE)
      df[[paste0("v", tag(a))]][j] <- q
      df[[paste0("e", tag(a))]][j] <- mean(win[win <= q])
    }
  }
  df
}
xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))
kupiec <- function(v, alpha) {
  n <- length(v); x <- sum(v); p <- x / n
  ll <- function(q) xlogy(n - x, 1 - q) + xlogy(x, q)
  LR <- -2 * (ll(alpha) - ll(p)); c(LR = LR, p = pchisq(LR, 1, lower.tail = FALSE))
}
christoffersen_ind <- function(v) {
  v <- as.integer(v); a <- v[-length(v)]; b <- v[-1]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1)
  n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  p01 <- n01 / (n00 + n01); p11 <- ifelse(n10 + n11 == 0, 0, n11 / (n10 + n11))
  p   <- (n01 + n11) / (n00 + n01 + n10 + n11)
  ll0 <- xlogy(n00 + n10, 1 - p) + xlogy(n01 + n11, p)
  ll1 <- xlogy(n00, 1 - p01) + xlogy(n01, p01) + xlogy(n10, 1 - p11) + xlogy(n11, p11)
  LR <- -2 * (ll0 - ll1); c(LR = LR, p = pchisq(LR, 1, lower.tail = FALSE))
}
#engle manganelli
dq_test <- function(v, VaR, a, lags = 4) {
  h <- as.numeric(v) - a; n <- length(h)
  if (n < lags + 20) return(NA_real_)
  idx <- (lags + 1):n
  X <- cbind(1, sapply(1:lags, function(l) h[idx - l]), VaR[idx]); y <- h[idx]
  b <- tryCatch(solve(crossprod(X), crossprod(X, y)), error = function(e) NULL)
  if (is.null(b)) return(NA_real_)
  stat <- as.numeric(t(b) %*% crossprod(X) %*% b) / (a * (1 - a))
  pchisq(stat, df = ncol(X), lower.tail = FALSE)
}
es_test <- function(r, VaR, ES, sig) {
  v <- r < VaR
  if (sum(v) < 5) return(c(mean_resid = NA, p_less = NA))
  d <- (r[v] - ES[v]) / sig[v]; tt <- t.test(d)
  c(mean_resid = unname(tt$estimate), p_less = pt(unname(tt$statistic), tt$parameter))
}
tick <- function(r, q, alpha) (alpha - (r < q)) * (r - q)

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

ret <- returns_yahoo("^FCHI")
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
    cfg_dfs[[id]] <- risk_config(x, i0, g$model, g$mean0, g$dist, g$Wc)
    saveRDS(cfg_dfs[[id]], f)
    cat(sprintf("[%2d/%d] %s done in %.1f min\n", i, nrow(grid), id,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
hs <- hs_risk(x, i0)

#where config has no cvg fit
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
baseline <- "sGARCH_mC_std_W1000|param" 

sel  <- dt <= split_date
hold <- !sel
cat(sprintf("Selection period: %s -> %s (%d days) | hold-out: %s -> %s (%d days)\n",
            dt[1], max(dt[sel]), sum(sel), min(dt[hold]), dt[length(dt)], sum(hold)))
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
show <- unique(c(best, sel_tab$candidate[2:3], baseline, "HS250"))
cat("\n=== HOLD-OUT (2017-2026): not used for selection ===\n")
for (a in alphas) {
  cat(sprintf("\n--- %.1f%% VaR / ES (expected exceedances: %.1f of %d) ---\n",
              100 * a, a * sum(hold), sum(hold)))
  tab <- do.call(rbind, lapply(show, function(nm) cbind(candidate = nm, eval_cand(cands[[nm]], a, r, hold))))
  print(format(tab[, -2], digits = 3), row.names = FALSE)
}
all_res <- do.call(rbind, lapply(names(cands), function(nm) rbind(
  cbind(candidate = nm, period = "selection 2008-2016", do.call(rbind, lapply(alphas, function(a) eval_cand(cands[[nm]], a, r, sel)))),
  cbind(candidate = nm, period = "hold-out 2017-2026", do.call(rbind, lapply(alphas, function(a) eval_cand(cands[[nm]], a, r, hold)))))))
write.csv(all_res, "outputs/var_experiments_all.csv", row.names = FALSE)
write.csv(sel_tab, "outputs/var_experiments_selection.csv", row.names = FALSE)
yr <- format(dt, "%Y")
ex <- sapply(show, function(nm) tapply(r < cands[[nm]]$var[["010"]], yr, sum))
ex <- cbind(days = as.integer(table(yr)), ex)
cat("\nExceedances of the 1% VaR by calendar year (expected about 2.5 per 250 days):\n")
print(ex)

#cumulative 1 per cent
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