# run-tests: one test file in its own R process, as the laptop's per-file runs (tools/testing/test-baseline.R, load
# mode load_all): pkgload::load_all(<package>, export_all = TRUE, helpers = TRUE), then testthat::test_file() with the
# summary reporter (printed to this process's log) and a list reporter (the counts).
#
# Usage: Rscript run-file.R <package dir> <test file> <out dir>
# Writes <out dir>/rows/<file>.tsv (file, passed, failed, errors, skipped, warnings, seconds, status ok|failed) and
# <out dir>/failures/<file>.txt (one block per failing or erroring expectation: file, test, message, first trace lines).
# A process that dies or is killed writes no row; the shard runner records it as crashed or timeout.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("usage: Rscript run-file.R <package dir> <test file> <out dir>")
package <- normalizePath(arguments[[1L]])
path <- normalizePath(arguments[[2L]])
out_dir <- normalizePath(arguments[[3L]])
name <- basename(path)
started <- Sys.time()

suppressMessages(pkgload::load_all(package, quiet = TRUE, export_all = TRUE, helpers = TRUE))
list_reporter <- testthat::ListReporter$new()
reporter <- testthat::MultiReporter$new(list(testthat::SummaryReporter$new(), list_reporter))
outcome <- tryCatch({
    testthat::test_file(path, reporter = reporter, stop_on_failure = FALSE)
    NULL
}, error = function(error) conditionMessage(error))

results <- list_reporter$get_results()
df <- as.data.frame(results)
oneLine <- function(text) gsub("[\t\r\n]+", " ", text)
blocks <- character()
# Per test (in run order) its counts, and per failing or erroring expectation its full message: the base-versus-head
# comparison (assemble.R) reads them from <out dir>/results/<file>.rds.
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
        trace <- if (!is.null(result$trace)) tryCatch(format(result$trace), error = function(e) character()) else character()
        blocks <- c(blocks, paste(c(
            sprintf("file: %s", name),
            sprintf("test: %s", oneLine(if (is.null(test$test)) "<outside test_that>" else test$test)),
            sprintf("kind: %s", if (inherits(result, "expectation_error")) "error" else "failure"),
            "message:", paste0("  ", strsplit(conditionMessage(result), "\n", fixed = TRUE)[[1L]]),
            if (length(trace)) c("trace (first lines):", paste0("  ", utils::head(trace, 12L)))), collapse = "\n"))
    }
}
if (!is.null(outcome)) {
    blocks <- c(blocks, paste(c(sprintf("file: %s", name), "test: <test_file itself>", "kind: error",
        "message:", paste0("  ", strsplit(outcome, "\n", fixed = TRUE)[[1L]])), collapse = "\n"))
}

failed <- if (nrow(df)) sum(df$failed) else 0L
errors <- (if (nrow(df)) sum(df$error) else 0L) + as.integer(!is.null(outcome))
row <- data.frame(file = name, passed = if (nrow(df)) sum(df$passed) else 0L, failed = failed, errors = errors,
    skipped = if (nrow(df)) sum(df$skipped) else 0L, warnings = if (nrow(df)) sum(df$warning) else 0L,
    seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
    status = if (failed + errors > 0L) "failed" else "ok", stringsAsFactors = FALSE)
stem <- sub("[.]R$", "", name)
if (!is.null(outcome)) messages[nrow(messages) + 1L, ] <- list("<test_file itself>", "error", outcome)
dir.create(file.path(out_dir, "results"), showWarnings = FALSE)
saveRDS(list(tests = per_test, messages = messages), file.path(out_dir, "results", paste0(stem, ".rds")))
if (length(blocks)) writeLines(paste0(blocks, "\n"), file.path(out_dir, "failures", paste0(stem, ".txt")))
utils::write.table(row, file.path(out_dir, "rows", paste0(stem, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nCOUNTS %s pass %d fail %d error %d skip %d warn %d\n", name, row$passed, row$failed, row$errors,
    row$skipped, row$warnings))
