# KIT-173 T7 (a) and (b) in one runner job: a record made by the base code replays on this commit.
#
# Usage, from the package root on a run-tests runner, after the commit under test (head) is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t7-replay.R <out dir> [<base commit>]
#
# R1 ran capture-record.R on the tests-only commit and its upload withheld the project store, which the replay needs
# (the replay checks the record's import and design states by entity id). Here the base commit's tree (default
# 62a41ba1, whose R/ equals a60c4e9a) is installed into a private library in the job, capture-record.R runs on it, and
# test-go-propagation.R runs on the installed head with KIT173_REPLAY_BASE pointing at that capture, so T7(b) replays
# the base record on its own store. The store stays in the job. Writes to <out dir>: t7-base-check.tsv (which code made
# the capture), t7-capture.tsv (the capture's values), t7-tests.tsv (per test: expectations, failures, skips,
# errors) and t7-failures.tsv (each failing expectation: test, class, first 12 message lines; added under VAL-54). Install and test consoles stay in the job; their exit statuses are in t7-base-check.tsv.

args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(args) >= 1L) args[[1L]] else "k173-t7-out", mustWork = FALSE)
base <- if (length(args) >= 2L) args[[2L]] else "62a41ba1ed35cecc8562c2ad9bd603b8bf9151fb"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k173-t7")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
capture <- file.path(work, "capture")
for (d in c(src, lib, capture)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
checks <- list(base_commit = base)

checks$archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
checks$install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(work, "install.log"),
    stderr = file.path(work, "install.log"))
probe <- system2(rscript, c("-e", shQuote(paste("suppressPackageStartupMessages(library(MultiScholaR));",
    "cat(find.package('MultiScholaR'), exists('propagateGoAnnotations', asNamespace('MultiScholaR')), sep = '\\t')"))),
    stdout = TRUE, env = paste0("R_LIBS=", lib))
checks$base_package_path_is_private_library <- startsWith(probe[[length(probe)]], lib)
checks$base_has_propagateGoAnnotations <- sub("^.*\t", "", probe[[length(probe)]])
checks$capture_exit <- system2(rscript, c(file.path("docs", "qa", "workflow-execution", "kit", "kit-173",
    "capture-record.R"), shQuote(capture)), stdout = file.path(work, "capture.log"), stderr = file.path(work, "capture.log"),
    env = paste0("R_LIBS=", lib))
file.copy(file.path(capture, "capture.tsv"), file.path(out, "t7-capture.tsv"))

code <- paste0("res <- testthat::test_file('tests/testthat/test-go-propagation.R', package = 'MultiScholaR', ",
    "load_package = 'installed', reporter = 'silent', stop_on_failure = FALSE); ",
    "df <- as.data.frame(res); utils::write.table(df[, c('test', 'nb', 'failed', 'skipped', 'error')], ",
    "file = '", file.path(out, "t7-tests.tsv"), "', sep = '\\t', quote = FALSE, row.names = FALSE); ",
    # VAL-54 instrument fix: the failing expectations' test name, class and the first 12 lines of their message
    "bad <- do.call(rbind, lapply(res, function(t) do.call(rbind, lapply(t$results, function(e) ",
    "if (inherits(e, 'expectation_success') || inherits(e, 'expectation_skip')) NULL else data.frame(test = t$test, ",
    "class = class(e)[[1L]], message = gsub('\\t', ' ', paste(utils::head(strsplit(conditionMessage(e), '\\n')[[1L]], 12L), ",
    "collapse = ' | ')), stringsAsFactors = FALSE))))); ",
    "if (is.null(bad)) bad <- data.frame(test = character(), class = character(), message = character()); ",
    "utils::write.table(bad, file = '", file.path(out, "t7-failures.tsv"), "', sep = '\\t', quote = FALSE, row.names = FALSE)")
checks$tests_exit <- system2(rscript, c("-e", shQuote(code)), stdout = file.path(work, "tests.log"),
    stderr = file.path(work, "tests.log"), env = c(paste0("KIT173_REPLAY_BASE=", capture), "NOT_CRAN=true"))
utils::write.table(data.frame(measure = names(checks), value = vapply(checks, as.character, "")),
    file.path(out, "t7-base-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
if (file.exists(file.path(out, "t7-tests.tsv"))) print(utils::read.delim(file.path(out, "t7-tests.tsv")))
