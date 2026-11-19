# KIT-129 integration, second pass: the failing-first evidence of the head's tests. The head's versions of the named
# test files run against the base commit's package (its R/, src/ and helpers) and against the head's, one R process
# per file and side, as the runner's per-file runs load the package (pkgload::load_all, export_all, helpers). A test
# written for a change fails on base and passes on head; every other test passes on both.
#
# Usage (from the package root; full git history): Rscript failing-first.R <base commit> <out dir> <test file> ...
# Writes <out dir>/failing-first-tests.tsv (side, file, test, passed, failed, errors, skipped, first message) and
# <out dir>/failing-first-summary.tsv (file, test, base and head counts, verdict).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 3L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
files <- arguments[-(1:2)]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")

tree <- file.path(tempdir(), "failing-first-base")
stopifnot(system2("git", c("worktree", "add", "--detach", tree, base)) == 0L)
for (file in files) {
    target <- file.path(tree, "tests", "testthat", file)
    stopifnot(system2("git", c("show", paste0("HEAD:tests/testthat/", file)), stdout = target) == 0L)
}

child <- '
arguments <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(arguments[[1L]], quiet = TRUE, export_all = TRUE, helpers = TRUE))
reporter <- testthat::ListReporter$new()
outcome <- tryCatch({ testthat::test_file(arguments[[2L]], reporter = reporter, stop_on_failure = FALSE); "" },
    error = function(e) conditionMessage(e))
one <- function(text) substr(gsub("[\\t\\r\\n]+", " ", text), 1L, 400L)
rows <- lapply(reporter$get_results(), function(test) {
    kinds <- vapply(test$results, function(r) class(r)[[1L]], "")
    bad <- Filter(function(r) inherits(r, c("expectation_failure", "expectation_error")), test$results)
    data.frame(test = one(test$test %||% "<outside test_that>"), passed = sum(kinds == "expectation_success"),
        failed = sum(kinds == "expectation_failure"), errors = sum(kinds == "expectation_error"),
        skipped = sum(kinds == "expectation_skip"), first_message = if (length(bad)) one(conditionMessage(bad[[1L]])) else "",
        stringsAsFactors = FALSE)
})
rows <- do.call(rbind, rows)
if (nzchar(outcome)) rows <- rbind(rows, data.frame(test = "<test_file itself>", passed = 0L, failed = 0L, errors = 1L, skipped = 0L,
    first_message = one(outcome), stringsAsFactors = FALSE))
utils::write.table(rows, arguments[[3L]], sep = "\\t", quote = FALSE, row.names = FALSE)
'
script <- file.path(tempdir(), "failing-first-child.R")
writeLines(child, script)

results <- list()
for (side in c("base", "head")) {
    package <- if (side == "base") tree else head
    for (file in files) {
        target <- file.path(tempdir(), sprintf("%s-%s.tsv", side, file))
        status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", script, package,
            file.path(package, "tests", "testthat", file), target), stdout = file.path(out_dir, sprintf("%s-%s.log", side, file)),
            stderr = file.path(out_dir, sprintf("%s-%s.log", side, file)))
        rows <- if (file.exists(target)) utils::read.delim(target, stringsAsFactors = FALSE, quote = "") else
            data.frame(test = "<process>", passed = 0L, failed = 0L, errors = 1L, skipped = 0L,
                first_message = sprintf("the process exited with status %d and wrote no results", status), stringsAsFactors = FALSE)
        results[[length(results) + 1L]] <- cbind(side = side, file = file, rows, stringsAsFactors = FALSE)
    }
}
tests <- do.call(rbind, results)
utils::write.table(tests, file.path(out_dir, "failing-first-tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

key <- function(t) paste(t$file, t$test, sep = "\r")
base_rows <- tests[tests$side == "base", ]
head_rows <- tests[tests$side == "head", ]
keys <- unique(c(key(head_rows), key(base_rows)))
bad <- function(t) t$failed + t$errors > 0L
summary <- do.call(rbind, lapply(keys, function(k) {
    b <- base_rows[key(base_rows) == k, ][1L, ]
    h <- head_rows[key(head_rows) == k, ][1L, ]
    verdict <- if (is.na(h$file)) "absent on head" else if (is.na(b$file)) "absent on base" else
        if (bad(b) && !bad(h)) "fails on base, passes on head" else if (!bad(b) && !bad(h)) "passes on both" else
        if (bad(b) && bad(h)) "fails on both" else "passes on base, fails on head"
    data.frame(file = strsplit(k, "\r")[[1L]][[1L]], test = strsplit(k, "\r")[[1L]][[2L]],
        base_passed = b$passed, base_failed = b$failed + b$errors, head_passed = h$passed, head_failed = h$failed + h$errors,
        head_skipped = h$skipped, verdict = verdict, base_first_message = b$first_message, stringsAsFactors = FALSE)
}))
utils::write.table(summary, file.path(out_dir, "failing-first-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("base:", base), paste("head:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("files:", paste(files, collapse = ", ")), R.version.string), file.path(out_dir, "failing-first-meta.txt"))
print(table(summary$verdict))
system2("git", c("worktree", "remove", "--force", tree))
