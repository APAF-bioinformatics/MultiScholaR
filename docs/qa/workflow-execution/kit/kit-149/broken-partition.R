# KIT-149, audit ENG-06 closure: the identity test of the partitioned alignment must fail when partitioning is done
# wrong. The head's test file runs once as it is and once per deliberately broken mechanism, each in its own R process
# (pkgload::load_all, export_all, helpers, as the runner's per-file runs load the package):
#   - none: the implementation as committed (every test passes);
#   - ties_reversed: the most-variable selection breaks ties by the later feature first;
#   - offset_lost: every page is handed over as if it were the first (offset 0), so the internal standards of the
#     second page are dropped by the wrong rows;
#   - pages_reversed: the pages are handed over last page first (with their own offsets), so the selection sees the
#     features in another order and ties across the page boundary are taken from the wrong page.
# Expected: `none` passes; each broken mechanism fails the ENG-06 test.
#
# Usage (from the package root): Rscript broken-partition.R <out dir> [<test file>]
# Writes <out dir>/broken-partition.tsv (mode, test, passed, failed, errors, skipped, first message).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 1L)
out_dir <- arguments[[1L]]
file <- if (length(arguments) >= 2L) arguments[[2L]] else "test-integration-partitioned-alignment.R"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
package <- normalizePath(".")

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, export_all = TRUE, helpers = TRUE))
mode <- arguments[[3L]]
ns <- asNamespace("MultiScholaR")
rows_of <- get("partitionByRows", ns)
mocks <- switch(mode,
    none = list(),
    ties_reversed = list(integrationSelectTopFeatures = function(statistic, n) {
        sort(order(-statistic, -seq_along(statistic))[seq_len(n)])
    }),
    offset_lost = list(partitionByRows = function(scope, columns, kernel, ...) {
        rows_of(scope, columns, function(page, index, offset) kernel(page, index, 0L), ...)
    }),
    pages_reversed = list(partitionByRows = function(scope, columns, kernel, consume, ...) {
        held <- list()
        result <- rows_of(scope, columns, function(page, index, offset) {
            held[[length(held) + 1L]] <<- list(page = page, index = index, offset = offset)
            NULL
        }, function(value, index) NULL, ...)
        for (h in rev(held)) consume(kernel(h$page, h$index, h$offset), h$index)
        result
    }))
reporter <- testthat::ListReporter$new()
run <- function() testthat::test_file(arguments[[2L]], reporter = reporter, stop_on_failure = FALSE)
outcome <- tryCatch({
    if (length(mocks)) do.call(testthat::with_mocked_bindings, c(list(code = quote(run())), mocks, list(.package = "MultiScholaR"))) else run()
    ""
}, error = function(e) conditionMessage(e))
one <- function(text) substr(gsub("[\\t\\r\\n]+", " ", text), 1L, 400L)
rows <- lapply(reporter$get_results(), function(test) {
    kinds <- vapply(test$results, function(r) class(r)[[1L]], "")
    bad <- Filter(function(r) inherits(r, c("expectation_failure", "expectation_error")), test$results)
    data.frame(test = one(test$test), passed = sum(kinds == "expectation_success"),
        failed = sum(kinds == "expectation_failure"), errors = sum(kinds == "expectation_error"),
        skipped = sum(kinds == "expectation_skip"), first_message = if (length(bad)) one(conditionMessage(bad[[1L]])) else "",
        stringsAsFactors = FALSE)
})
rows <- do.call(rbind, rows)
if (nzchar(outcome)) rows <- rbind(rows, data.frame(test = "<test_file itself>", passed = 0L, failed = 0L, errors = 1L,
    skipped = 0L, first_message = one(outcome), stringsAsFactors = FALSE))
utils::write.table(rows, arguments[[4L]], sep = "\\t", quote = FALSE, row.names = FALSE)
'
script <- file.path(tempdir(), "broken-partition-child.R")
writeLines(child, script)

results <- list()
for (mode in c("none", "ties_reversed", "offset_lost", "pages_reversed")) {
    target <- file.path(tempdir(), sprintf("%s.tsv", mode))
    log <- file.path(out_dir, sprintf("%s.log", mode))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, package,
        file.path(package, "tests", "testthat", file), mode, target), stdout = log, stderr = log)
    rows <- if (file.exists(target)) utils::read.delim(target, stringsAsFactors = FALSE, quote = "") else
        data.frame(test = "<process>", passed = 0L, failed = 0L, errors = 1L, skipped = 0L,
            first_message = sprintf("the process exited with status %d and wrote no results", status), stringsAsFactors = FALSE)
    results[[mode]] <- cbind(mode = mode, rows, stringsAsFactors = FALSE)
}
table <- do.call(rbind, results)
utils::write.table(table, file.path(out_dir, "broken-partition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table[c("mode", "test", "passed", "failed", "errors")], right = FALSE)
eng06 <- grepl("ENG-06", table$test)
verdict <- all(table$failed[table$mode == "none"] == 0L & table$errors[table$mode == "none"] == 0L) &&
    all(vapply(c("ties_reversed", "offset_lost", "pages_reversed"), function(mode) {
        rows <- table[table$mode == mode & eng06, , drop = FALSE]
        nrow(rows) > 0L && any(rows$failed + rows$errors > 0L)
    }, logical(1)))
cat(sprintf("KIT-149 broken-partition verdict: %s\n", if (verdict) "the identity test fails on every broken mechanism and passes on the implementation" else "NOT AS EXPECTED"))
quit(save = "no", status = if (verdict) 0L else 1L)
