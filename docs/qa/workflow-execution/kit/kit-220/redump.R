# KIT-220 control C and run D (ticket 10; owner 2026-11-25). Evidence code under docs/qa: no package code. Runs on a
# GitHub-hosted runner (run-tests `scripts`) from the package root: base (staging) once, the merged head twice.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-220/redump.R <out dir> <test file>...
#
# Each test file runs in its own process as run-tests does (pkgload::load_all with helpers, testthat::test_file), with
# the dump switches set to files under <out dir>:
#   MULTISCHOLAR_KIT112_W1_DUMP (test-workflow-saved-da.R), MULTISCHOLAR_KIT106_W1_DUMP (test-omic-da-module.R),
#   MULTISCHOLAR_KIT107_W1_DUMP (test-omic-enrich-module.R), MULTISCHOLAR_GUI_ENTRY_RECORD_DIR (each GUI-entry lane's
#   record as JSON, compared afterwards with the stored oracle by key). The console of each file is kept.
# The comparison of the dumps (base, head, head) and the stored values is made from these files; nothing here judges.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
out_dir <- normalizePath(args[[1L]], mustWork = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
files <- args[-1L]
status <- data.frame(file = character(), exit = integer(), seconds = numeric(), stringsAsFactors = FALSE)
for (file in files) {
    stem <- sub("[.]R$", "", file)
    env <- c(sprintf("MULTISCHOLAR_KIT112_W1_DUMP=%s", file.path(out_dir, paste0(stem, "-kit112.R"))),
        sprintf("MULTISCHOLAR_KIT106_W1_DUMP=%s", file.path(out_dir, paste0(stem, "-kit106.R"))),
        sprintf("MULTISCHOLAR_KIT107_W1_DUMP=%s", file.path(out_dir, paste0(stem, "-kit107.R"))),
        sprintf("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=%s", file.path(out_dir, "gui-entry")),
        "MULTISCHOLAR_TRX_MANIFEST_NOTE=1", "NOT_CRAN=true")
    expr <- sprintf(paste0("pkgload::load_all('.', export_all = TRUE, helpers = TRUE, quiet = TRUE); ",
        "r <- testthat::test_file(file.path('tests', 'testthat', '%s'), stop_on_failure = FALSE, reporter = 'summary'); ",
        "d <- as.data.frame(r); utils::write.table(d[, intersect(c('test', 'nb', 'failed', 'skipped', 'error'), names(d))], ",
        "'%s', sep = '\\t', quote = FALSE, row.names = FALSE)"), file, file.path(out_dir, paste0(stem, "-results.tsv")))
    started <- Sys.time()
    exit <- system2("Rscript", c("--vanilla", "-e", shQuote(expr)), env = env,
        stdout = file.path(out_dir, paste0(stem, "-console.log")), stderr = file.path(out_dir, paste0(stem, "-console.log")))
    status[nrow(status) + 1L, ] <- list(file, as.integer(exit), round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1))
    print(status[nrow(status), ])
}
utils::write.table(status, file.path(out_dir, "status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
