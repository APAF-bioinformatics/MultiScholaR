# KIT-173 declared control 3 (outcome.md): a record captured on the base and replayed on the base, in two separate
# Rscript processes (t7-cross-capture.R, then t7-cross-replay.R with variant "basebase"), both on the base installed in
# a private library, in the runner's locale. Changes no test and no product code.
#
# Usage, from the package root on a run-tests runner:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t7-basebase.R <out dir> [<base commit>]

args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(args) >= 1L) args[[1L]] else "k173-basebase-out", mustWork = FALSE)
base <- if (length(args) >= 2L) args[[2L]] else "62a41ba1ed35cecc8562c2ad9bd603b8bf9151fb"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k173-basebase")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
for (d in c(src, lib)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
here <- file.path("docs", "qa", "workflow-execution", "kit", "kit-173")
checks <- list(base_commit = base)
checks$archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
checks$install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(work, "install.log"),
    stderr = file.path(work, "install.log"))
step <- function(name, script, script_args) {
    log <- file.path(work, paste0(name, ".log"))
    status <- system2(rscript, c(file.path(here, script), script_args), stdout = log, stderr = log,
        env = paste0("R_LIBS=", lib))
    if (!identical(status, 0L)) {
        writeLines(utils::tail(readLines(log, warn = FALSE), 15L), file.path(out, paste0(name, "-error-tail.txt")))
    }
    status
}
dir <- file.path(work, "pair")
checks$capture_exit <- step("capture", "t7-cross-capture.R", shQuote(dir))
checks$replay_exit <- step("replay", "t7-cross-replay.R", c("basebase", shQuote(dir), shQuote(out)))
utils::write.table(data.frame(measure = names(checks), value = vapply(checks, as.character, "")),
    file.path(out, "basebase-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
