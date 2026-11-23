# KIT-172 section 7(b): re-dump the normalisation tab's frozen interface values through the test file's own dump path
# (MULTISCHOLAR_KIT105_W1_DUMP, `omicNormDump`, tests/testthat/test-omic-norm-module.R), as the per-file runner loads
# the package (pkgload::load_all, export_all and helpers). The dump is written to <out dir>/omic-norm-dump.R; the
# test file's own results go to <out dir>/omic-norm-results.tsv.
#
# Usage (run-tests scripts lane, from the package root): Rscript docs/qa/workflow-execution/kit/kit-172/redump-norm-ui.R <out dir>
args <- commandArgs(trailingOnly = TRUE)
dir.create(args[[1L]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[1L]])
dump <- file.path(out, "omic-norm-dump.R")
unlink(dump)
Sys.setenv(MULTISCHOLAR_KIT105_W1_DUMP = dump, NOT_CRAN = "true")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
reporter <- testthat::ListReporter$new()
testthat::test_file("tests/testthat/test-omic-norm-module.R", reporter = reporter, stop_on_failure = FALSE)
results <- as.data.frame(reporter$get_results())
utils::write.table(results[, intersect(c("test", "nb", "failed", "skipped", "error", "warning"), names(results))],
  file.path(out, "omic-norm-results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("dump %s: %s\n", dump, file.exists(dump)))
