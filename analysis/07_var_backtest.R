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
alphas     <- c(0.01, 0.05)     
hs_win     <- 250            

returns_yahoo <- function(symbol) {
  px <- getSymbols(symbol, from = data_start, to = cfg$end_date, auto.assign = FALSE)
  px <- Ad(px); px <- px[!is.na(px)]
  r  <- na.omit(diff(log(px))); colnames(r) <- "ret"; r
}

make_spec <- function(model) {
  ugarchspec(
    variance.model     = list(model = model, garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = cfg$distribution
  )
}

fit_safe <- function(model, dat, m) {
  for (sv in c("hybrid", "nlminb", "solnp")) {
    fit <- tryCatch(suppressWarnings(
      ugarchfit(make_spec(model), data = dat, out.sample = m, solver = sv)),
      error = function(e) NULL)
    if (!is.null(fit) && convergence(fit) == 0) return(fit)
  }
  NULL
}

#rugarch (analytical var student-t)
# z = T * sqrt((nu-2)/nu) avec T ~ t(nu)
std_t_var <- function(alpha, nu) qt(alpha, nu) * sqrt((nu - 2) / nu)
std_t_es  <- function(alpha, nu) {
  ta <- qt(alpha, nu)
  -dt(ta, nu) / alpha * (nu + ta^2) / (nu - 1) * sqrt((nu - 2) / nu)
}
tag <- function(a) sprintf("%02d", round(a * 100))
garch_risk <- function(x, i0, model) {
  n <- length(x)
  starts <- seq(i0, n, by = B)
  out <- vector("list", length(starts)); n_fail <- 0
  for (k in seq_along(starts)) {
    s <- starts[k]; e <- min(s + B - 1, n); m <- e - s
    fit <- fit_safe(model, x[(s - W):(e - 1)], m)
    if (is.null(fit)) {
      n_fail <- n_fail + 1
      out[[k]] <- data.frame(idx = s:e, mu = NA_real_, sigma = NA_real_, nu = NA_real_)
    } else {
      fc <- ugarchforecast(fit, n.ahead = 1, n.roll = m)
      out[[k]] <- data.frame(idx = s:e, mu = as.numeric(fitted(fc)),
                             sigma = as.numeric(sigma(fc)),
                             nu = unname(coef(fit)["shape"]))
    }
    if (k %% 10 == 0) cat(sprintf("  %s: block %d / %d\n", model, k, length(starts)))
  }
  df <- do.call(rbind, out)
  for (a in alphas) {
    df[[paste0("var", tag(a))]] <- df$mu + df$sigma * std_t_var(a, df$nu)
    df[[paste0("es",  tag(a))]] <- df$mu + df$sigma * std_t_es(a, df$nu)
  }
  attr(df, "n_fail") <- n_fail
  df
}

hs_risk <- function(x, i0) {
  n <- length(x); idx <- i0:n
  df <- data.frame(idx = idx, sigma = NA_real_)
  for (a in alphas) { df[[paste0("var", tag(a))]] <- NA_real_
  df[[paste0("es",  tag(a))]] <- NA_real_ }
  for (j in seq_along(idx)) {
    win <- x[(idx[j] - hs_win):(idx[j] - 1)]
    df$sigma[j] <- sd(win)
    for (a in alphas) {
      q <- quantile(win, a, type = 1, names = FALSE)
      df[[paste0("var", tag(a))]][j] <- q
      df[[paste0("es",  tag(a))]][j]  <- mean(win[win <= q])
    }
  }
  df
}

xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))

kupiec <- function(v, alpha) {                
  n <- length(v); x <- sum(v); p <- x / n
  ll <- function(q) xlogy(n - x, 1 - q) + xlogy(x, q)
  LR <- -2 * (ll(alpha) - ll(p))
  c(LR = LR, p = pchisq(LR, 1, lower.tail = FALSE))
}
christoffersen_ind <- function(v) {   
  v <- as.integer(v); a <- v[-length(v)]; b <- v[-1]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1)
  n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  p01 <- n01 / (n00 + n01); p11 <- ifelse(n10 + n11 == 0, 0, n11 / (n10 + n11))
  p   <- (n01 + n11) / (n00 + n01 + n10 + n11)
  ll0 <- xlogy(n00 + n10, 1 - p) + xlogy(n01 + n11, p)
  ll1 <- xlogy(n00, 1 - p01) + xlogy(n01, p01) + xlogy(n10, 1 - p11) + xlogy(n11, p11)
  LR <- -2 * (ll0 - ll1)
  c(LR = LR, p = pchisq(LR, 1, lower.tail = FALSE))
}
#mcneil frey
es_test <- function(r, VaR, ES, sig) {
  v <- r < VaR
  if (sum(v) < 5) return(c(mean_resid = NA, p_two = NA, p_less = NA))
  d  <- (r[v] - ES[v]) / sig[v]
  tt <- t.test(d)
  c(mean_resid = unname(tt$estimate), p_two = tt$p.value,
    p_less = pt(unname(tt$statistic), tt$parameter))   # H1: losses worse than ES
}
tick <- function(r, q, alpha) (alpha - (r < q)) * (r - q)    # quantile (pinball) loss
dm_plain <- function(d) { z <- mean(d) / (sd(d) / sqrt(length(d)))
c(stat = z, p = 2 * pnorm(-abs(z))) }

run_backtest <- function(ret, label) {
  cat(sprintf("\n==================== %s ====================\n", label))
  x <- as.numeric(ret); d <- index(ret); n <- length(x)
  i0 <- max(which(d >= eval_start)[1], W + 1, hs_win + 1)
  cat(sprintf("Returns %s -> %s (N = %d). First target day: %s, %d out-of-sample days.\n",
              d[1], d[n], n, d[i0], n - i0 + 1))
  
  cat("garch\n");  g <- garch_risk(x, i0, "sGARCH")
  cat("egarch\n"); e <- garch_risk(x, i0, "eGARCH")
  h <- hs_risk(x, i0)
  cat(sprintf("Blocks without a converged fit: GARCH = %d | EGARCH = %d\n",
              attr(g, "n_fail"), attr(e, "n_fail")))
  
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
      m <- mods[[nm]]; VaR <- m[[paste0("var", tag(a))]]; ES <- m[[paste0("es", tag(a))]]
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
    L <- lapply(mods, function(m) tick(r, m[[paste0("var", tag(a))]], a))
    for (pr in list(c("GARCH", "EGARCH"), c("GARCH", "HS250"), c("EGARCH", "HS250"))) {
      z <- dm_plain(L[[pr[1]]] - L[[pr[2]]])
      cat(sprintf("  %-6s vs %-6s: stat = %+.2f, p = %.3f\n", pr[1], pr[2], z[1], z[2]))
    }
  }
  
  yr <- format(dt, "%Y")
  ex <- sapply(names(mods), function(nm) tapply(r < mods[[nm]]$var01, yr, sum))
  ex <- cbind(days = as.integer(table(yr)), ex)
  cat("\n Exceedances of the 1% VaR by calendar year:\n")
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
    out[[paste0(nm, "_var", tag(a))]] <- mods[[nm]][[paste0("var", tag(a))]]
    out[[paste0(nm, "_es",  tag(a))]] <- mods[[nm]][[paste0("es",  tag(a))]]
  }
  write.csv(out, sprintf("outputs/var_es_forecasts_%s.csv", label), row.names = FALSE)
  invisible(list(table = res, forecasts = out))
}

res_cac <- run_backtest(returns_yahoo("^FCHI"), "CAC40")