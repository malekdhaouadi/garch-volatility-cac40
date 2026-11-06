suppressPackageStartupMessages({
  library(rugarch); library(yaml)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path  <- normalizePath(sub("^--file=", "", file_arg))
  project_root <- normalizePath(file.path(dirname(script_path), ".."))
  setwd(project_root)
}

cfg <- yaml::read_yaml("config.yml")

baseline_fp <- "outputs/baseline/paper_fits.rds"
current_fp  <- file.path(cfg$pipeline$models_dir, "latest_fits.rds")
current_meta_fp <- file.path(cfg$pipeline$outputs_dir, "latest.rds")

if (!file.exists(baseline_fp))     stop("Baseline not found: ", baseline_fp)
if (!file.exists(current_fp))      stop("Current fits not found: ", current_fp,
                                        "\nRun `Rscript run_pipeline.R config.live.yml` first")
if (!file.exists(current_meta_fp)) stop("Current metadata not found: ", current_meta_fp)

baseline <- readRDS(baseline_fp)
current  <- readRDS(current_fp)
current_meta <- readRDS(current_meta_fp)$meta

# The paper baseline is pinned to 2026-10-26 by config.yml at the time it was
# saved. Recorded here so the report is self-describing. If you regenerate the
# baseline, update this line to the new end_date.
baseline_end <- "2026-10-26"
current_end  <- current_meta$end_date

# --- Coefficients to monitor ---
garch_b <- coef(baseline$garch)
garch_c <- coef(current$garch)
eg_b    <- coef(baseline$egarch)
eg_c    <- coef(current$egarch)

# --- Thresholds (diagnostic, not automated alerts: see README) ---
th <- list(
  garch_persistence = 0.02,
  garch_alpha       = 0.02,
  garch_beta        = 0.02,
  garch_shape       = 1.0,
  egarch_alpha      = 0.03,
  egarch_beta       = 0.02,
  egarch_gamma      = 0.02,
  egarch_shape      = 1.0
)

pers_b <- sum(garch_b[c("alpha1","beta1")])
pers_c <- sum(garch_c[c("alpha1","beta1")])

checks <- data.frame(
  metric = c("GARCH persistence (a+b)",
             "GARCH alpha1",
             "GARCH beta1",
             "GARCH nu (Student-t)",
             "EGARCH alpha1 (leverage)",
             "EGARCH beta1",
             "EGARCH gamma1",
             "EGARCH nu (Student-t)"),
  baseline  = c(pers_b, garch_b["alpha1"], garch_b["beta1"], garch_b["shape"],
                eg_b["alpha1"], eg_b["beta1"], eg_b["gamma1"], eg_b["shape"]),
  current   = c(pers_c, garch_c["alpha1"], garch_c["beta1"], garch_c["shape"],
                eg_c["alpha1"], eg_c["beta1"], eg_c["gamma1"], eg_c["shape"]),
  threshold = unlist(th),
  row.names = NULL
)
checks$delta <- checks$current - checks$baseline
checks$drift <- abs(checks$delta) > checks$threshold

cat("\n===== Drift report: current vs paper baseline =====\n")
cat("Baseline: outputs/baseline/paper_fits.rds (paper model, end_date ",
    baseline_end, ")\n", sep = "")
cat("Current:  outputs/models/latest_fits.rds (end_date ",
    current_end, ", mode ", current_meta$mode, ")\n\n", sep = "")
print(format(checks, digits = 4), row.names = FALSE)

drifted <- checks[checks$drift, ]
if (nrow(drifted) == 0) {
  cat("\nSTATUS: OK -- no metric exceeds its threshold\n")
  quit(status = 0)
} else {
  cat("\nSTATUS: DRIFT DETECTED\n")
  cat("Affected metrics:\n")
  print(drifted[, c("metric","delta","threshold")], row.names = FALSE)
  quit(status = 1)
}