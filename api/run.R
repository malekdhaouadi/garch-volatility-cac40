# api/run.R
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) stop("Run via: Rscript api/run.R")
script_path  <- normalizePath(sub("^--file=", "", file_arg))
project_root <- normalizePath(file.path(dirname(script_path), ".."))
setwd(project_root)

suppressPackageStartupMessages({
  library(plumber); library(yaml)
})

cfg <- read_yaml("config.yml")

# Absolute paths — plumber re-parses api/plumber.R in its own working directory,
# so anything plumber.R opens must be an absolute path.
latest_fp <- normalizePath(file.path(cfg$pipeline$outputs_dir, "latest.rds"),
                           mustWork = FALSE)

if (!file.exists(latest_fp))
  stop("latest.rds not found at ", latest_fp, " — run `Rscript run_pipeline.R config.live.yml` first")

.GARCH_CFG      <- cfg
.GARCH_LATEST_FP <- latest_fp

cat("[run.R] project_root =", project_root, "\n")
cat("[run.R] latest_fp    =", .GARCH_LATEST_FP, "\n")
cat("[run.R] port         =", .GARCH_CFG$api$port, "\n")

plumber::pr("api/plumber.R") |>
  plumber::pr_set_debug(FALSE) |>
  plumber::pr_run(host = "0.0.0.0", port = .GARCH_CFG$api$port)