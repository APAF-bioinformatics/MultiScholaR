# KIT-173 declared control (outcome.md, "Declared control"): is T7(b)'s geneID-order difference KIT-173's, or replay's?
#
# Usage, from the package root on a run-tests runner, after the head commit is installed:
#   Rscript docs/qa/workflow-execution/kit/kit-173/t7-control.R <out dir> [<base commit>]
#
# Installs the base commit (default 62a41ba1, R/ equal to a60c4e9a) into a private library as t7-replay.R does, then
# runs t7-control-pair.R twice: pair "base" (capture and replay on the base) and pair "head" (capture and replay on the
# installed head). Changes no test and no product code. Writes the pairs' tables and control-check.tsv.

args <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(args) >= 1L) args[[1L]] else "k173-control-out", mustWork = FALSE)
base <- if (length(args) >= 2L) args[[2L]] else "62a41ba1ed35cecc8562c2ad9bd603b8bf9151fb"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k173-control")
src <- file.path(work, "base-src")
lib <- file.path(work, "base-lib")
for (d in c(src, lib)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
pair_script <- file.path("docs", "qa", "workflow-execution", "kit", "kit-173", "t7-control-pair.R")
checks <- list(base_commit = base)
checks$archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(src)))
checks$install_exit <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load",
    paste0("--library=", shQuote(lib)), shQuote(src)), stdout = file.path(work, "install.log"),
    stderr = file.path(work, "install.log"))
checks$pair_base_exit <- system2(rscript, c(pair_script, "base", shQuote(out)), stdout = file.path(work, "base.log"),
    stderr = file.path(work, "base.log"), env = paste0("R_LIBS=", lib))
checks$pair_head_exit <- system2(rscript, c(pair_script, "head", shQuote(out)), stdout = file.path(work, "head.log"),
    stderr = file.path(work, "head.log"))
# A pair that failed returns the last 15 lines of its console (an error message), nothing else.
for (side in c("base", "head")) if (!identical(checks[[paste0("pair_", side, "_exit")]], 0L))
    writeLines(utils::tail(readLines(file.path(work, paste0(side, ".log")), warn = FALSE), 15L),
        file.path(out, paste0(side, "-error-tail.txt")))
utils::write.table(data.frame(measure = names(checks), value = vapply(checks, as.character, "")),
    file.path(out, "control-check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(checks)
for (f in list.files(out, "\\.tsv$", full.names = TRUE)) print(utils::read.delim(f)[, 1:4])
