# run-tests.R: the acceptance test on your func_general_nipals.R. Run from this folder: Rscript run-tests.R
source("load-candidate.R"); load_candidate()
library(testthat)
r <- as.data.frame(test_file("test-identity.R", reporter = "summary", stop_on_failure = FALSE))
cat(sprintf("\nRESULT tests %d expectations %d failed %d error %d skipped %d\n", nrow(r), sum(r$nb), sum(r$failed), sum(r$error), sum(r$skipped)))
