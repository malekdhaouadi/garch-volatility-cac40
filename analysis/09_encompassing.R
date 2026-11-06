library(ragg)

# Uses the out-of-sample forecasts saved by R/06_extended_eval.R==

H     <- 21
lags  <- c(20, 40)
files <- c(CAC40    = "outputs/ext_iv_comparison_CAC40.csv",
           STOXX50E = "outputs/ext_iv_comparison_STOXX50E.csv")

#OLS with newy west
ols_hac <- function(y, X, lag) {
  X <- cbind(1, X); n <- nrow(X)
  XtXi <- solve(crossprod(X))
  b <- as.numeric(XtXi %*% crossprod(X, y))
  e <- as.numeric(y - X %*% b)
  u <- X * e
  S <- crossprod(u) / n
  for (l in 1:lag) {
    G <- crossprod(u[(l + 1):n, , drop = FALSE], u[1:(n - l), , drop = FALSE]) / n
    S <- S + (1 - l / (lag + 1)) * (G + t(G))
  }
  V <- n * XtXi %*% S %*% XtXi
  list(b = b, se = sqrt(diag(V)), n = n,
       r2 = 1 - sum(e^2) / sum((y - mean(y))^2))
}

specs <- list(
  "1: implied"          = c("implied"),
  "2: GARCH"            = c("garch"),
  "3: EGARCH"           = c("egarch"),
  "4: implied + GARCH"  = c("implied", "garch"),
  "5: implied + EGARCH" = c("implied", "egarch")
)

periods <- list(
  "2008-2012"  = c("2008-01-01", "2012-12-31"),
  "2013-2019"  = c("2013-01-01", "2019-12-31"),
  "2020-2026"  = c("2020-01-01", "2026-12-31")
)

all_rows <- list(); coef_plot <- list()

for (lab in names(files)) {
  if (!file.exists(files[lab])) { cat(sprintf("\n%s: %s not found, skipped\n", lab, files[lab])); next }
  d <- read.csv(files[lab]); d$date <- as.Date(d$date)
  d <- d[complete.cases(d), ]
  cat(sprintf("\n==================== %s ====================\n", lab))
  cat(sprintf("n = %d origins (%s -> %s), about %d independent 21-day windows\n",
              nrow(d), min(d$date), max(d$date), round(nrow(d) / H)))
  cat(sprintf("Correlation of regressors: implied-GARCH %.3f | implied-EGARCH %.3f | GARCH-EGARCH %.3f\n",
              cor(d$implied, d$garch), cor(d$implied, d$egarch), cor(d$garch, d$egarch)))
  
  for (scale in c("vol", "log")) {
    tr <- if (scale == "log") log else identity
    cat(sprintf("\n--- %s: %s regressions, dependent variable = realized vol over next 21 days ---\n",
                lab, ifelse(scale == "vol", "level (vol points)", "log")))
    cat(sprintf("(t20, t40 = t-statistics with Newey-West lag %d and %d)\n", lags[1], lags[2]))
    rows <- list()
    for (sn in names(specs)) {
      v <- specs[[sn]]
      y <- tr(d$rv_fwd); X <- as.matrix(tr(d[, v, drop = FALSE]))
      f20 <- ols_hac(y, X, lags[1]); f40 <- ols_hac(y, X, lags[2])
      rows[[sn]] <- data.frame(dataset = lab, scale = scale, spec = sn,
                               term = c("(Intercept)", v), est = f20$b, se20 = f20$se,
                               t20 = f20$b / f20$se, t40 = f40$b / f40$se,
                               r2 = f20$r2, n = f20$n, row.names = NULL)
    }
    tab <- do.call(rbind, rows)
    print(format(tab[, c("spec", "term", "est", "se20", "t20", "t40", "r2")], digits = 3),
          row.names = FALSE)
    all_rows[[paste(lab, scale)]] <- tab
    
    if (scale == "vol") {
      g <- function(sn, tm, col) tab[tab$spec == sn & tab$term == tm, col]
      cat("\nKey tests (level regressions, HAC lag 20):\n")
      cat(sprintf("  Slope on VSTOXX alone = %.3f; test slope = 1: t = %+.2f\n",
                  g("1: implied", "implied", "est"),
                  (g("1: implied", "implied", "est") - 1) / g("1: implied", "implied", "se20")))
      cat(sprintf("  Does GARCH add to VSTOXX?   coef on GARCH  = %+.3f (t = %+.2f) | coef on VSTOXX = %+.3f (t = %+.2f)\n",
                  g("4: implied + GARCH", "garch", "est"), g("4: implied + GARCH", "garch", "t20"),
                  g("4: implied + GARCH", "implied", "est"), g("4: implied + GARCH", "implied", "t20")))
      cat(sprintf("  Does EGARCH add to VSTOXX?  coef on EGARCH = %+.3f (t = %+.2f) | coef on VSTOXX = %+.3f (t = %+.2f)\n",
                  g("5: implied + EGARCH", "egarch", "est"), g("5: implied + EGARCH", "egarch", "t20"),
                  g("5: implied + EGARCH", "implied", "est"), g("5: implied + EGARCH", "implied", "t20")))
      cat(sprintf("  R2: VSTOXX %.3f | GARCH %.3f | EGARCH %.3f | VSTOXX+GARCH %.3f | VSTOXX+EGARCH %.3f\n",
                  g("1: implied", "implied", "r2"), g("2: GARCH", "garch", "r2"),
                  g("3: EGARCH", "egarch", "r2"), g("4: implied + GARCH", "garch", "r2"),
                  g("5: implied + EGARCH", "egarch", "r2")))
      for (pr in list(c("4: implied + GARCH", "garch", "GARCH"), c("5: implied + EGARCH", "egarch", "EGARCH"))) {
        for (tm in c("implied", pr[2])) {
          est <- g(pr[1], tm, "est"); se <- g(pr[1], tm, "se20")
          coef_plot[[length(coef_plot) + 1]] <- data.frame(
            label = sprintf("%s: %s (with %s)", lab, ifelse(tm == "implied", "VSTOXX", pr[3]), pr[3]),
            est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
        }
      }
    }
  }
  
  cat(sprintf("\n--- %s by sub-period (level regressions, HAC lag 20; descriptive only) ---\n", lab))
  for (pn in names(periods)) {
    s <- d[d$date >= as.Date(periods[[pn]][1]) & d$date <= as.Date(periods[[pn]][2]), ]
    if (nrow(s) < 250) next
    for (sp in c("4: implied + GARCH", "5: implied + EGARCH")) {
      v <- specs[[sp]]
      f <- ols_hac(s$rv_fwd, as.matrix(s[, v, drop = FALSE]), lags[1])
      cat(sprintf("  %-10s n = %4d  %-20s b_implied = %+.3f (t = %+.2f) | b_%s = %+.3f (t = %+.2f)\n",
                  pn, nrow(s), sp, f$b[2], f$b[2] / f$se[2], v[2], f$b[3], f$b[3] / f$se[3]))
    }
  }
}
if (length(all_rows)) {
  write.csv(do.call(rbind, all_rows), "outputs/encompassing_regressions.csv", row.names = FALSE)
}
if (length(coef_plot)) {
  cp <- do.call(rbind, coef_plot)
  agg_png("figures/encompassing_coefficients.png", width = 1200, height = 650,
          units = "px", res = 110, background = "white")
  par(mar = c(4, 17, 3, 1))
  yy <- rev(seq_len(nrow(cp)))
  plot(cp$est, yy, xlim = range(cp$lo, cp$hi, 0), yaxt = "n", pch = 19, col = "darkred",
       xlab = "coefficient (95% Newey-West interval)", ylab = "",
       main = "Encompassing regressions: forward realized volatility")
  segments(cp$lo, yy, cp$hi, yy, col = "darkred", lwd = 2)
  axis(2, at = yy, labels = cp$label, las = 1, cex.axis = 0.8)
  abline(v = 0, lty = 2, col = "grey50")
  dev.off()
}