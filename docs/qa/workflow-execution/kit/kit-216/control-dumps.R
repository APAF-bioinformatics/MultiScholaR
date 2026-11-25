# KIT-216 controls (ticket 6; VAL-174's form): the dumps of the three re-frozen files on this commit.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-216/control-dumps.R <label> <repeats>
#   Each repeat runs test-omic-enrich-module.R (MULTISCHOLAR_KIT107_W1_DUMP, MULTISCHOLAR_KIT113_W1_DUMP),
#   test-omic-summary-module.R (MULTISCHOLAR_KIT108_W1_DUMP) and test-gui-entry-proteomics.R
#   (MULTISCHOLAR_GUI_ENTRY_RECORD_DIR), each in its own R process, into
#   docs/qa/workflow-execution/kit/kit-216/control/<label>-<repeat>/.
arguments <- commandArgs(trailingOnly = TRUE)
label <- arguments[[1L]]
repeats <- as.integer(arguments[[2L]])
root <- file.path("docs", "qa", "workflow-execution", "kit", "kit-216", "control")
runner <- file.path(root, "run-one.R")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "suppressMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = TRUE))",
    "testthat::test_file(file.path('tests', 'testthat', args[[1L]]), reporter = 'summary', stop_on_failure = FALSE)"),
    runner)
for (r in seq_len(repeats)) {
    out <- normalizePath(file.path(root, sprintf("%s-%d", label, r)), mustWork = FALSE)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    jobs <- list(
        list(file = "test-omic-enrich-module.R", env = c(MULTISCHOLAR_KIT107_W1_DUMP = file.path(out, "kit107-w1.R"),
            MULTISCHOLAR_KIT113_W1_DUMP = file.path(out, "kit113-w1.R"))),
        list(file = "test-omic-summary-module.R", env = c(MULTISCHOLAR_KIT108_W1_DUMP = file.path(out, "kit108-w1.R"))),
        list(file = "test-gui-entry-proteomics.R", env = c(MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = file.path(out, "gui-entry"))))
    for (job in jobs) {
        log <- file.path(out, paste0(sub("[.]R$", "", job$file), ".log"))
        env <- paste0(names(job$env), "=", job$env)
        status <- system2("env", c(env, "Rscript", "--vanilla", runner, job$file), stdout = log, stderr = log)
        cat(sprintf("KIT-216 control %s-%d %s exit %d\n", label, r, job$file, status))
    }
}
