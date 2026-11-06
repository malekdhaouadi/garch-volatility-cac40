suppressPackageStartupMessages({
  library(yaml); library(xts); library(rugarch); library(quantmod)
})
source("R/data.R"); source("R/models.R"); source("R/forecast.R"); source("R/stats.R")

args     <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) > 0) args[1] else "config.yml"
if (!file.exists(cfg_path)) stop("Config file not found: ", cfg_path)
cfg <- read_yaml(cfg_path)
set.seed(cfg$seed)
dir.create(cfg$pipeline$outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$pipeline$models_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$pipeline$log_dir,     recursive = TRUE, showWarnings = FALSE)

is_live <- grepl("live", basename(cfg_path), ignore.case = TRUE)

ts       <- format(Sys.time(), "%Y%m%dT%H%M%SZ")
log_path <- file.path(cfg$pipeline$log_dir, sprintf("run_%s.log", ts))
log_con  <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)

log_line <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n"); cat(msg, "\n", file = log_con)
}
fail <- function(msg) { log_line("FAIL: ", msg); quit(status = 1) }

cache_path <- file.path(cfg$pipeline$outputs_dir, "..", "data", "prices.rds")
if (!is_live && file.exists(cache_path)) {
  prices <- readRDS(cache_path)
  log_line("Loaded cached prices: ", nrow(prices), " rows, ",
           as.character(start(prices)), " -> ", as.character(end(prices)),
           " (from ", cache_path, ")")
} else {
  log_line("Config: ", cfg_path, " | fetch ", cfg$index,
           " ", cfg$start_date, " -> ", cfg$end_date)
  prices <- tryCatch(
    fetch_prices(cfg$index, cfg$start_date, cfg$end_date),
    error = function(e) fail(paste("fetch_prices:", conditionMessage(e)))
  )
  prices <- prices[index(prices) <= as.Date(cfg$end_date)]
}
n_na <- sum(is.na(prices))
if (n_na > 0) {
  log_line("Dropping ", n_na, " NA price row(s): ",
           paste(as.character(index(prices)[is.na(prices)]), collapse = ", "))
  prices <- prices[!is.na(prices)]
}

v <- if (is_live) {
  validate_prices(prices, max_staleness_days = cfg$api$max_staleness_days)
} else {
  validate_prices(prices, max_staleness_days = cfg$api$max_staleness_days,
                  today = end(prices))
}
if (!all(v$ok)) {
  fail(paste0("Price validation failed: ",
              paste(v$check[!v$ok], collapse = "; "),
              " | last price: ", as.character(end(prices)),
              " | today: ", Sys.Date(),
              " | mode: ", if (is_live) "live" else "reproduction"))
}
log_line("Validation passed (", nrow(prices), " prices, ",
         as.character(end(prices)), ", mode = ",
         if (is_live) "live" else "reproduction", ")")

returns <- compute_log_returns(prices)
if (anyNA(returns) || nrow(returns) < 500)
  fail("returns series invalid or too short")
log_line("Returns: ", nrow(returns), " obs, sd = ", round(sd(returns), 6))

log_line("Fitting GARCH(1,1) + EGARCH(1,1) with ", cfg$distribution)
fit_g <- fit_safe("sGARCH", as.numeric(returns), dist = cfg$distribution)
fit_e <- fit_safe("eGARCH", as.numeric(returns), dist = cfg$distribution)
if (is.null(fit_g) || is.null(fit_e)) fail("model fit failed")

log_line("GARCH  persistence = ",
         round(sum(coef(fit_g)[c("alpha1","beta1")]), 4),
         " | omega = ", signif(coef(fit_g)["omega"], 3),
         " | loglik = ", round(likelihood(fit_g), 1))
log_line("EGARCH beta = ", round(coef(fit_e)["beta1"], 4),
         " | alpha1 (leverage) = ", round(coef(fit_e)["alpha1"], 4),
         " | loglik = ", round(likelihood(fit_e), 1))

H <- cfg$forecast_horizon
fc_g <- ugarchforecast(fit_g, n.ahead = H)
fc_e <- ugarchforecast(fit_e, n.ahead = H)
path_g <- as.numeric(sigma(fc_g)) * sqrt(252) * 100
path_e <- as.numeric(sigma(fc_e)) * sqrt(252) * 100
avg_g  <- sqrt(mean(path_g^2))
avg_e  <- sqrt(mean(path_e^2))

log_line("Forecast origin: ", as.character(end(returns)),
         " | day1: GARCH = ", round(path_g[1], 2), "% | EGARCH = ", round(path_e[1], 2), "%",
         " | 30d avg: GARCH = ", round(avg_g, 2), "% | EGARCH = ", round(avg_e, 2), "%")
if (abs(avg_g - avg_e) > 5) {
  log_line("WARN: 30d forecast divergence GARCH vs EGARCH = ",
           round(abs(avg_g - avg_e), 2), " pp (weak identification when omega ~ 0)")
}

iv_info <- tryCatch({
  if (file.exists("outputs/iv_var.rds")) {
    iv <- readRDS("outputs/iv_var.rds")
    list(k_selected = iv$k_selected, k_encompassing = iv$k_encompassing,
         nu = iv$nu, last_vstoxx = tail(iv$daily$vstoxx, 1))
  } else {
    log_line("NOTE: outputs/iv_var.rds not found — skipping IV-VaR section")
    NULL
  }
}, error = function(e) {
  log_line("WARN: IV-VaR read failed: ", conditionMessage(e)); NULL
})

latest <- list(
  meta = list(
    run_utc       = Sys.time(),
    config        = cfg_path,
    mode          = if (is_live) "live" else "reproduction",
    index         = cfg$index,
    end_date      = as.character(end(returns)),
    n_returns     = nrow(returns),
    distribution  = cfg$distribution,
    forecast_horizon = H,
    git_commit    = tryCatch(system("git rev-parse --short HEAD", intern = TRUE),
                             error = function(e) NA_character_)
  ),
  fit = list(
    garch = list(
      coef = as.list(coef(fit_g)),
      loglik = as.numeric(likelihood(fit_g)),
      aic = infocriteria(fit_g)["Akaike", 1],
      persistence = unname(sum(coef(fit_g)[c("alpha1","beta1")]))
    ),
    egarch = list(
      coef = as.list(coef(fit_e)),
      loglik = as.numeric(likelihood(fit_e)),
      aic = infocriteria(fit_e)["Akaike", 1],
      persistence = unname(coef(fit_e)["beta1"])
    )
  ),
  forecast = list(
    horizon = seq_len(H),
    garch_vol = path_g,
    egarch_vol = path_e,
    day1_garch = path_g[1],
    day1_egarch = path_e[1],
    avg30_garch = avg_g,
    avg30_egarch = avg_e
  ),
  iv_var = iv_info,
  diagnostics = list(
    skewness = mean((as.numeric(returns) - mean(returns))^3) / sd(returns)^3,
    excess_kurt = mean((as.numeric(returns) - mean(returns))^4) / sd(returns)^4 - 3,
    jb_pvalue = tryCatch(tseries::jarque.bera.test(as.numeric(returns))$p.value,
                         error = function(e) NA_real_)
  )
)
saveRDS(latest, file.path(cfg$pipeline$outputs_dir, "latest.rds"))
saveRDS(list(garch = fit_g, egarch = fit_e),
        file.path(cfg$pipeline$models_dir, "latest_fits.rds"))

log_line("Wrote outputs/latest.rds and outputs/models/latest_fits.rds")
log_line("PIPELINE OK")
quit(status = 0)