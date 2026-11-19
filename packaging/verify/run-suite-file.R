# KIT-134 (branch run-k134-rc): one test file in its own R process against the bundle's INSTALLED MultiScholaR
# (testthat::test_file with package = "MultiScholaR", load_package = "installed": the test environment is a child of the
# installed namespace, the helpers are sourced from tests/testthat, as R CMD check runs tests). Adapted from the
# run-tests workflow's run-file.R (branch run-tests), whose only difference is pkgload::load_all() in place of the
# installed package. Started by run-suite.R with the bundle's R.
#
# Usage: Rscript --vanilla run-suite-file.R <checkout> <test-tool library> <out dir> <test file>
# Writes <out dir>/rows/<file>.tsv and <out dir>/results/<file>.rds (per test counts and failing messages). A process
# that dies or is killed writes no row; run-suite.R records it as crashed or timeout.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) stop("usage: Rscript run-suite-file.R <checkout> <test-tool library> <out dir> <test file>")
tools_lib <- normalizePath(arguments[[2L]])
out_dir <- normalizePath(arguments[[3L]])
path <- normalizePath(arguments[[4L]])
name <- basename(path)
started <- Sys.time()
.libPaths(c(.libPaths(), tools_lib))

suppressPackageStartupMessages(library(testthat))
list_reporter <- testthat::ListReporter$new()
reporter <- testthat::MultiReporter$new(list(testthat::SummaryReporter$new(), list_reporter))
outcome <- tryCatch({
    testthat::test_file(path, reporter = reporter, stop_on_failure = FALSE, package = "MultiScholaR",
        load_package = "installed")
    NULL
}, error = function(error) conditionMessage(error))
cat(sprintf("\nINSTALLED MultiScholaR from %s\n", find.package("MultiScholaR")))

results <- list_reporter$get_results()
df <- as.data.frame(results)
per_test <- data.frame(test = character(), passed = integer(), failed = integer(), errors = integer(),
    skipped = integer(), stringsAsFactors = FALSE)
messages <- data.frame(test = character(), kind = character(), message = character(), stringsAsFactors = FALSE)
for (test in results) {
    name_of <- if (is.null(test$test)) "<outside test_that>" else test$test
    kinds <- vapply(test$results, function(r) class(r)[[1L]], "")
    per_test[nrow(per_test) + 1L, ] <- list(name_of, sum(kinds == "expectation_success"), sum(kinds == "expectation_failure"),
        sum(kinds == "expectation_error"), sum(kinds == "expectation_skip"))
    for (result in test$results) {
        if (!inherits(result, c("expectation_failure", "expectation_error"))) next
        messages[nrow(messages) + 1L, ] <- list(name_of, if (inherits(result, "expectation_error")) "error" else "failure",
            conditionMessage(result))
    }
}
if (!is.null(outcome)) messages[nrow(messages) + 1L, ] <- list("<test_file itself>", "error", outcome)

failed <- if (nrow(df)) sum(df$failed) else 0L
errors <- (if (nrow(df)) sum(df$error) else 0L) + as.integer(!is.null(outcome))
row <- data.frame(file = name, passed = if (nrow(df)) sum(df$passed) else 0L, failed = failed, errors = errors,
    skipped = if (nrow(df)) sum(df$skipped) else 0L, warnings = if (nrow(df)) sum(df$warning) else 0L,
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    status = if (failed + errors > 0L) "failed" else "ok", stringsAsFactors = FALSE)
stem <- sub("[.]R$", "", name)
saveRDS(list(tests = per_test, messages = messages), file.path(out_dir, "results", paste0(stem, ".rds")))
utils::write.table(row, file.path(out_dir, "rows", paste0(stem, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nCOUNTS %s pass %d fail %d error %d skip %d warn %d\n", name, row$passed, row$failed, row$errors,
    row$skipped, row$warnings))
