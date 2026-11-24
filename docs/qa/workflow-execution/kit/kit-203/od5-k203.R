# KIT-203 OD5: the control for re-freezing the enrichment tab's stored UI capture `omicEnrichW1FrozenUi`
# (tests/testthat/test-omic-enrich-module.R). Evidence code under docs/qa: no package code. Runs on a GitHub-hosted
# runner (run-tests `scripts`) from the package root of the commit under test (head).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-203/od5-k203.R <out dir> <base commit>
#
# Three separate processes, each `pkgload::load_all()` of a source tree, then test (a) of that tree's
# test-omic-enrich-module.R alone (`testthat::test_file(desc = )`) with MULTISCHOLAR_KIT107_W1_DUMP set, so test (a)
# writes its capture as R code: once on the base commit's tree (git archive), twice on head. Precondition (K203-P12):
# the base dump equals the base tree's stored block. Writes to <out dir>: dump-base.R, dump-head-1.R, dump-head-2.R,
# od5.tsv (check, pass, detail). Exit status 0.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
out <- normalizePath(args[[1L]], mustWork = FALSE)
base <- args[[2L]]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k203-od5")
src <- file.path(work, "base-src")
dir.create(src, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
DESC <- "(a) the enrichment UI keeps its HTML, ids, test ids, tab titles, labels, help texts, choices and defaults"
archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))

dump <- function(tree, file) {
    code <- paste0("suppressMessages(pkgload::load_all(", deparse(tree), ", export_all = TRUE, helpers = TRUE, ",
        "quiet = TRUE)); res <- testthat::test_file(file.path(", deparse(tree), ", 'tests', 'testthat', ",
        "'test-omic-enrich-module.R'), desc = ", deparse(DESC), ", reporter = 'summary', stop_on_failure = FALSE); ",
        "df <- as.data.frame(res); cat('TESTS', nrow(df), 'FAILED', sum(df$failed), 'ERROR', sum(df$error), '\\n')")
    log <- file.path(out, paste0(basename(file), ".log.txt"))
    exit <- system2(rscript, c("-e", shQuote(code)), stdout = log, stderr = log,
        env = c(paste0("MULTISCHOLAR_KIT107_W1_DUMP=", file), "NOT_CRAN=true"))
    list(exit = exit, tail = paste(utils::tail(readLines(log, warn = FALSE), 3L), collapse = " | "))
}
runs <- list(base = dump(src, file.path(out, "dump-base.R")),
    head_1 = dump(normalizePath("."), file.path(out, "dump-head-1.R")),
    head_2 = dump(normalizePath("."), file.path(out, "dump-head-2.R")))

read_bytes <- function(f) if (file.exists(f)) readBin(f, "raw", file.info(f)$size) else raw()
stored_block <- function(tree) {
    lines <- readLines(file.path(tree, "tests", "testthat", "test-omic-enrich-module.R"), warn = FALSE)
    at <- c(grep("^# OMIC_ENRICH_FROZEN_UI_START$", lines), grep("^# OMIC_ENRICH_FROZEN_UI_END$", lines))
    paste(lines[(at[[1L]] + 1L):(at[[2L]] - 1L)], collapse = "\n")
}
dump_text <- function(f) if (file.exists(f)) sub("\\s+$", "", paste(readLines(f, warn = FALSE), collapse = "\n")) else ""
checks <- data.frame(check = c("base archive", paste("dump", names(runs)), "head dumps byte-identical",
    "base dump equals the base tree's stored block (precondition)", "head dump differs from the stored block"),
    pass = c(archive_exit == 0L, vapply(runs, function(r) r$exit == 0L, logical(1)),
        length(read_bytes(file.path(out, "dump-head-1.R"))) > 0L &&
            identical(read_bytes(file.path(out, "dump-head-1.R")), read_bytes(file.path(out, "dump-head-2.R"))),
        nzchar(dump_text(file.path(out, "dump-base.R"))) && identical(dump_text(file.path(out, "dump-base.R")),
            stored_block(src)),
        !identical(dump_text(file.path(out, "dump-head-1.R")), stored_block(src))),
    detail = c(paste("exit", archive_exit), vapply(runs, function(r) paste("exit", r$exit, ";", r$tail), ""), "", "",
        ""), stringsAsFactors = FALSE)
utils::write.table(checks, file.path(out, "od5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
