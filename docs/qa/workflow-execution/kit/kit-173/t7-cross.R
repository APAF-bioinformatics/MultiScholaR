# KIT-173 declared control 2 (outcome.md): the cross pair (base capture, head replay; two processes, as in T7i1) run
# with both processes in C collation (variant a) and in the runner's locale (variant b), plus the collation of the
# T7i1 arrangement.
#
# Usage, from the package root on a run-tests runner, after the head commit is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t7-cross.R <out dir> [<base commit>]
#
# Writes the variants' tables (t7-cross-replay.R), t7i1-collation.tsv and cross-check.tsv (exit statuses; a failed
# process's last 15 console lines as <step>-error-tail.txt). Changes no test and no product code.

args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(args) >= 1L) args[[1L]] else "k173-cross-out", mustWork = FALSE)
base <- if (length(args) >= 2L) args[[2L]] else "62a41ba1ed35cecc8562c2ad9bd603b8bf9151fb"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k173-cross")
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

step <- function(name, script, script_args, env) {
    log <- file.path(work, paste0(name, ".log"))
    status <- system2(rscript, c(file.path(here, script), script_args), stdout = log, stderr = log, env = env)
    if (!identical(status, 0L)) {
        writeLines(utils::tail(readLines(log, warn = FALSE), 15L), file.path(out, paste0(name, "-error-tail.txt")))
    }
    status
}
variants <- list(a = "LC_COLLATE=C", b = character())
for (variant in names(variants)) {
    dir <- file.path(work, paste0("variant-", variant))
    checks[[paste0(variant, "_capture_exit")]] <- step(paste0(variant, "-capture"), "t7-cross-capture.R", shQuote(dir),
        c(paste0("R_LIBS=", lib), variants[[variant]]))
    checks[[paste0(variant, "_replay_exit")]] <- step(paste0(variant, "-replay"), "t7-cross-replay.R",
        c(variant, shQuote(dir), shQuote(out)), variants[[variant]])
}

# The T7i1 arrangement's collation: the capture as an Rscript process in the runner's locale; the replay inside a
# testthat test in a second process (test_file runs each test so).
runner <- system2(rscript, c("-e", shQuote("cat(Sys.getlocale('LC_COLLATE'))")), stdout = TRUE)
in_test <- system2(rscript, c("-e", shQuote(paste0("testthat::test_that('collation', { cat(Sys.getlocale('LC_COLLATE'))",
    "; testthat::succeed() })"))), stdout = TRUE)
utils::write.table(data.frame(variant = "T7i1", process = c("capture (base): Rscript", "replay (head): testthat test"),
    collation = c(paste(runner, collapse = " "), paste(in_test, collapse = " "))), file.path(out, "t7i1-collation.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(data.frame(measure = names(checks), value = vapply(checks, as.character, "")),
    file.path(out, "cross-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
