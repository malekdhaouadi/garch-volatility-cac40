# Run from the project root:  source("tests/run_tests.R")
testthat::test_dir("tests/testthat", reporter = "progress", stop_on_failure = TRUE)
