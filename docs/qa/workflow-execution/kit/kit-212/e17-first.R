# KIT-212 E17, failing first (implementation audit K212-I02): the row E17 of tests/testthat/test-kit212-da-current.R run on
# the tree of the tests-first commit (before the repair), in its own process. Evidence code under docs/qa: no package code.
# Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the commit under test.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-212/e17-first.R <out dir> <tests-first commit>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
src <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k212-e17-src")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
stopifnot(system(sprintf("git archive %s | tar -x -C %s", shQuote(args[[2L]]), shQuote(src))) == 0L)
desc <- "KIT-212 E17: lipid, the Summary tab mounted: after Copy, a second DA, then Report asks for Copy again"
code <- paste0("suppressMessages(pkgload::load_all(", deparse(src), ", export_all = TRUE, helpers = TRUE, quiet = TRUE)); ",
    "r <- testthat::test_file(file.path(", deparse(src), ", 'tests', 'testthat', 'test-kit212-da-current.R'), desc = ",
    deparse(desc), ", reporter = 'summary', stop_on_failure = FALSE); ",
    "d <- as.data.frame(r); write.table(d[, c('test', 'nb', 'failed', 'skipped', 'error')], ", deparse(file.path(out, "e17-first.tsv")),
    ", sep = '\\t', quote = FALSE, row.names = FALSE)")
log <- file.path(out, "e17-first.log.txt")
exit <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)), stdout = log, stderr = log, env = "NOT_CRAN=true")
writeLines(sprintf("exit=%s", exit), file.path(out, "e17-first-exit.txt"))
