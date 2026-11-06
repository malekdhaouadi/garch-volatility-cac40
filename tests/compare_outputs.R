# Compare the CSVs in outputs/ with copies saved BEFORE the refactor.
# Before re-running scripts:
#   dir.create("archive/outputs_before_refactor", recursive = TRUE)
#   file.copy(list.files("outputs", pattern = "\\.csv$", full.names = TRUE),
#             "archive/outputs_before_refactor")
# After re-running:   source("tests/compare_outputs.R")

old_dir <- "archive/outputs_before_refactor"
stopifnot(dir.exists(old_dir))
res <- do.call(rbind, lapply(list.files(old_dir, pattern = "\\.csv$"), function(f) {
  new_f <- file.path("outputs", f)
  if (!file.exists(new_f)) return(data.frame(file = f, same = NA, note = "missing in outputs/"))
  a <- read.csv(file.path(old_dir, f)); b <- read.csv(new_f)
  cmp <- all.equal(a, b, tolerance = 1e-8)
  data.frame(file = f, same = isTRUE(cmp), note = if (isTRUE(cmp)) "" else paste(cmp[1:min(2, length(cmp))], collapse = " | "))
}))
print(res, right = FALSE, row.names = FALSE)
cat(sprintf("\n%d of %d files identical\n", sum(res$same, na.rm = TRUE), nrow(res)))
