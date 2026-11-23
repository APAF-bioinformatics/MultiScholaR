# KIT-173 (VAL-61): after the geneID fix, a record captured on the head and replayed on the head in two separate
# Rscript processes (t7-cross-capture.R, then t7-cross-replay.R with variant "headhead"), on the installed head, in
# the runner's locale; the same columns as the declared controls. Changes no test and no product code.
#
# Usage, from the package root on a run-tests runner, after the commit under test is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t7-headhead.R <out dir>

args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(args) >= 1L) args[[1L]] else "k173-headhead-out", mustWork = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k173-headhead")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
here <- file.path("docs", "qa", "workflow-execution", "kit", "kit-173")
step <- function(name, script, script_args) {
    log <- file.path(work, paste0(name, ".log"))
    status <- system2(rscript, c(file.path(here, script), script_args), stdout = log, stderr = log)
    if (!identical(status, 0L)) {
        writeLines(utils::tail(readLines(log, warn = FALSE), 15L), file.path(out, paste0(name, "-error-tail.txt")))
    }
    status
}
dir <- file.path(work, "pair")
checks <- list(capture_exit = step("capture", "t7-cross-capture.R", shQuote(dir)))
checks$replay_exit <- step("replay", "t7-cross-replay.R", c("headhead", shQuote(dir), shQuote(out)))
utils::write.table(data.frame(measure = names(checks), value = vapply(checks, as.character, "")),
    file.path(out, "headhead-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
