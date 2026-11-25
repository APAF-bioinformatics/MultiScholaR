# KIT-226 control C (ticket section 5; VAL-196, VAL-208; audit K226-I07). Evidence code under docs/qa: no package code.
# Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root: base (staging 16112917 plus this script)
# once, the merged head twice, each run its own job. The form is KIT-220's redump.R.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-226/redump.R <out dir> <test file>...
#
# Each test file runs in its own process as run-tests does (pkgload::load_all with helpers, testthat::test_file), with
# the dump switches set to files under <out dir>: MULTISCHOLAR_KIT107_W1_DUMP (test-omic-enrich-module.R's W1
# captures as R code) and MULTISCHOLAR_GUI_ENTRY_RECORD_DIR (each GUI-entry lane's record as JSON). The console of each
# file is kept. The comparison of the dumps (base, head, head) with the stored values is made from these files; nothing
# here judges.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(args[[1L]], mustWork = TRUE)
files <- args[-1L]
status <- data.frame(file = character(), exit = integer(), seconds = numeric(), stringsAsFactors = FALSE)
for (file in files) {
    stem <- sub("[.]R$", "", file)
    env <- c(sprintf("MULTISCHOLAR_KIT107_W1_DUMP=%s", file.path(out_dir, paste0(stem, "-kit107.R"))),
        sprintf("MULTISCHOLAR_GUI_ENTRY_RECORD_DIR=%s", file.path(out_dir, "gui-entry")), "NOT_CRAN=true")
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
