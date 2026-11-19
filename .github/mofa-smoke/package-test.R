# KIT-062 (branch run-mofa-smoke-k062): run MultiScholaR test files against the installed package, with MOFA2 and the
# bundle's Python (started through the bundle's launcher, which sets RETICULATE_PYTHON). Each file runs in the
# package's namespace (testthat::test_file, load_package = "installed") from the package checkout's tests/testthat.
#
# Usage: Rscript --vanilla package-test.R <package checkout> <out dir> <test file> [<test file> ...]
# Writes <out dir>/tests.tsv (one row per test) and exits 1 when an expectation failed or errored, or when every test
# of a file skipped (on this lane a skip means MOFA2 or Python was not reached).

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 3L)
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
out_dir <- arguments[[2L]]
files <- arguments[-(1:2)]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("MultiScholaR %s from %s; RETICULATE_PYTHON %s\n", utils::packageVersion("MultiScholaR"),
    find.package("MultiScholaR"), Sys.getenv("RETICULATE_PYTHON")))
rows <- list()
status <- 0L
for (file in files) {
    path <- file.path(checkout, "tests", "testthat", file)
    started <- Sys.time()
    result <- tryCatch(testthat::test_file(path, package = "MultiScholaR", load_package = "installed",
        reporter = "summary", stop_on_failure = FALSE), error = function(e) {
        message("FILE ERROR ", file, ": ", conditionMessage(e))
        NULL
    })
    if (is.null(result)) {
        status <- 1L
        rows[[length(rows) + 1L]] <- data.frame(file = file, test = NA_character_, passed = 0L, failed = 0L,
            error = TRUE, skipped = FALSE, seconds = NA_real_, stringsAsFactors = FALSE)
        next
    }
    table <- as.data.frame(result)
    table <- data.frame(file = file, test = table$test, passed = table$nb - table$failed - table$skipped,
        failed = table$failed, error = table$error, skipped = table$skipped, seconds = round(table$real, 2),
        stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- table
    if (any(table$failed > 0L) || any(table$error) || all(table$skipped)) status <- 1L
    cat(sprintf("%s: %d tests, %d failed, %d errors, %d skipped, %.1f s\n", file, nrow(table), sum(table$failed > 0L),
        sum(table$error), sum(table$skipped), as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
utils::write.table(do.call(rbind, rows), file.path(out_dir, "tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
quit(save = "no", status = status)
