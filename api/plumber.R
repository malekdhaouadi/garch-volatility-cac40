suppressPackageStartupMessages({
  library(plumber); library(jsonlite)
})

if (!exists(".GARCH_CFG", envir = globalenv(), inherits = FALSE))
  stop("api/plumber.R must be launched via api/run.R (it needs .GARCH_CFG)")
if (!exists(".GARCH_LATEST_FP", envir = globalenv(), inherits = FALSE))
  stop("api/plumber.R must be launched via api/run.R (it needs .GARCH_LATEST_FP)")

GARCH_CFG       <- get(".GARCH_CFG",       envir = globalenv())
GARCH_LATEST_FP <- get(".GARCH_LATEST_FP", envir = globalenv())
latest_fp       <- GARCH_LATEST_FP

state <- new.env(parent = emptyenv())
state$loaded_utc <- NULL
state$latest     <- NULL
state$mtime      <- NULL

load_if_needed <- function() {
  if (!file.exists(latest_fp)) {
    stop(jsonlite::toJSON(list(status = "unavailable",
                               error  = paste("missing", latest_fp)),
                          auto_unbox = TRUE))
  }
  mt <- file.info(latest_fp)$mtime
  if (!is.null(state$mtime) && !is.na(mt) && identical(mt, state$mtime)) return()
  state$latest     <- readRDS(latest_fp)
  state$mtime      <- mt
  state$loaded_utc <- Sys.time()
}

#* @apiTitle GARCH forecast API
#* @apiDescription Serves the latest GARCH/EGARCH fit and forecast.

#* @get /health
#* @serializer unboxedJSON
function() {
  load_if_needed()
  age_days <- as.numeric(difftime(Sys.time(), state$mtime, units = "days"))
  stale    <- age_days > GARCH_CFG$api$max_staleness_days
  list(
    status          = if (stale) "stale" else "ok",
    mode            = state$latest$meta$mode,
    index           = state$latest$meta$index,
    end_date        = state$latest$meta$end_date,
    run_utc         = format(state$latest$meta$run_utc, "%Y-%m-%dT%H:%M:%SZ"),
    loaded_utc      = format(state$loaded_utc, "%Y-%m-%dT%H:%M:%SZ"),
    cache_age_days  = round(age_days, 2),
    max_staleness_days = GARCH_CFG$api$max_staleness_days,
    git_commit      = state$latest$meta$git_commit
  )
}

#* @get /model-info
#* @serializer unboxedJSON
function(req, res) {
  tryCatch({
    load_if_needed()
    list(meta = state$latest$meta,
         fit  = state$latest$fit,
         diagnostics = state$latest$diagnostics,
         iv_var = state$latest$iv_var)
  }, error = function(e) {
    res$status <- 503L
    list(status = "unavailable", error = conditionMessage(e))
  })
}

#* @param horizon Truncate the forecast path to this many days (default 30, max 30)
#* @param annualize "true" for annualized percent (default), "false" for daily
#* @get /forecast
#* @serializer unboxedJSON
function(req, res, horizon = 30L, annualize = "true") {
  tryCatch({
    load_if_needed()
    H <- min(as.integer(horizon), GARCH_CFG$forecast_horizon)
    if (is.na(H) || H < 1L) {
      res$status <- 400L
      return(list(status = "bad_request", error = "horizon must be >= 1"))
    }
    ann <- !identical(tolower(annualize), "false")
    f <- state$latest$forecast
    scale <- if (ann) 1 else 1 / sqrt(252)
    
    list(
      origin_date  = state$latest$meta$end_date,
      mode         = state$latest$meta$mode,
      units        = if (ann) "annualized_percent" else "daily_fraction",
      horizon      = H,
      day1_garch   = f$day1_garch   * scale,
      day1_egarch  = f$day1_egarch  * scale,
      avg30_garch  = f$avg30_garch  * scale,
      avg30_egarch = f$avg30_egarch * scale,
      garch_vol    = as.numeric(f$garch_vol[seq_len(H)])  * scale,
      egarch_vol   = as.numeric(f$egarch_vol[seq_len(H)]) * scale
    )
  }, error = function(e) {
    res$status <- 503L
    list(status = "unavailable", error = conditionMessage(e))
  })
}