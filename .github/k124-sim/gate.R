# KIT-124 confirmation round on hosted runners: the gate (carried over, decision M19) and the classic legacy
# comparison, against the committed reference rows (scored/confirm-reference/; this script decides nothing).
#
# Usage: Rscript --vanilla gate.R <sim: s1..s4> <reference dir> <tables dir> <label> <compare-results.R> <out dir>
#
# The run's rows: s<i>-runs.tsv of every artifact k124-<label>-<sim>-* under <tables dir>. The reference rows:
# <reference dir>/s<i>-runs.tsv.gz (and, for s4, s4-runs.laptop-extra.tsv.gz). Two readings, each through
# compare-laptop.R (the same rows of both sides, keyed by cell, run, arm and subset; a run key the reference does not
# hold stops it) and the fixed rule compare-results.R:
#   gate     arms 1 and 2: before a cell is scored they equal the earlier rows of the same seeds
#   legacy   arm 3 (s1 only): the installed package's legacy design against the earlier runner rows
# A cell for which the reference holds no row at all is listed as "no reference" (a warning, not a failure). Writes
# <out dir>/<sim>-<reading>-summary.tsv (one row per cell and column: verdict, largest relative difference; the
# rule-bearing columns of M19 marked) and exits 1 when a rule-bearing value of a referenced cell is outside the rule or
# differs, or when a reference row is missing.

RULE_BEARING <- c("status", "k", "n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) stop("usage: Rscript --vanilla gate.R <sim> <reference dir> <tables dir> <label> <compare-results.R> <out dir>")
sim <- arguments[[1L]]; reference_dir <- arguments[[2L]]; tables_dir <- arguments[[3L]]; label <- arguments[[4L]]
rule <- arguments[[5L]]; out_dir <- arguments[[6L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
read <- function(path) utils::read.delim(path, colClasses = "character", quote = "", na.strings = "NA")

parts <- list.files(tables_dir, pattern = sprintf("^k124-%s-%s-", label, sim), full.names = TRUE)
files <- file.path(parts, sprintf("%s-runs.tsv", sim))
files <- files[file.exists(files)]
if (!length(files)) {
    cat(sprintf("::warning::%s: no result table among the artifacts (a failed job uploads nothing); no gate\n", sim))
    quit(save = "no", status = 0L)
}
run_rows <- do.call(rbind, lapply(files, read))
reference_files <- file.path(reference_dir, c(sprintf("%s-runs.tsv.gz", sim), if (sim == "s4") "s4-runs.laptop-extra.tsv.gz"))
reference_files <- reference_files[file.exists(reference_files)]
if (!length(reference_files)) stop(sprintf("::error::%s: no committed reference table under %s", sim, reference_dir))
reference_rows <- do.call(rbind, lapply(reference_files, \(path) read(gzfile(path))))

failed <- FALSE
reading <- function(name, arms) {
    mine <- run_rows[run_rows$arm %in% arms, , drop = FALSE]
    if (!nrow(mine)) return(invisible(NULL))
    referenced <- unique(reference_rows$cell[reference_rows$arm %in% arms])
    unreferenced <- setdiff(unique(mine$cell), referenced)
    rows <- list()
    for (cell in unreferenced) {
        cat(sprintf("::warning::%s %s: cell %s has no reference rows (reported, not gated)\n", sim, name, cell))
        rows[[length(rows) + 1L]] <- data.frame(simulator = toupper(sim), reading = name, cell = cell, column = NA,
            verdict = "no reference", max_rel = NA, rule_bearing = NA)
    }
    mine <- mine[mine$cell %in% referenced, , drop = FALSE]
    if (nrow(mine)) {
        work <- file.path(out_dir, sprintf("%s-%s", sim, name))
        dir.create(file.path(work, "reference"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(work, "run"), showWarnings = FALSE)
        write <- function(x, dir) utils::write.table(x, file.path(work, dir, sprintf("%s-runs.tsv", sim)), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
        write(reference_rows, "reference")
        write(mine, "run")
        stem <- file.path(out_dir, sprintf("%s-%s", sim, name))
        status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", file.path(here, "compare-laptop.R"), sim,
            file.path(work, "reference"), file.path(work, "run"), paste(arms, collapse = ","), paste0(stem, "-a.rds"),
            paste0(stem, "-b.rds")))
        if (status != 0L) {
            cat(sprintf("::error::%s %s: the reference does not hold every row of this run (compare-laptop.R stopped)\n", sim, name))
            failed <<- TRUE
            return(invisible(NULL))
        }
        status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", rule, paste0(stem, "-a.rds"), paste0(stem, "-b.rds"),
            paste0(stem, "-compare.tsv")), stdout = FALSE)
        if (status != 0L) stop("compare-results.R failed for ", basename(stem))
        unlink(work, recursive = TRUE)
        compared <- utils::read.delim(paste0(stem, "-compare.tsv"), stringsAsFactors = FALSE, na.strings = "")
        path_name <- sub("\\$.*$", "", compared$path)
        column <- sub("^[^$]*\\$", "", compared$path)
        cell <- sub("^\\S+ ", "", path_name)
        for (key in unique(paste(cell, column))) {
            at <- paste(cell, column) == key
            verdicts <- compared$verdict[at]
            verdict <- if (all(verdicts == "agree")) "agree" else if (any(verdicts == "differ")) "differ" else "outside"
            rows[[length(rows) + 1L]] <- data.frame(simulator = toupper(sim), reading = name, cell = cell[at][[1L]],
                column = column[at][[1L]], verdict = verdict,
                max_rel = if (any(is.finite(compared$max_rel[at]))) max(compared$max_rel[at], na.rm = TRUE) else NA_real_,
                rule_bearing = column[at][[1L]] %in% RULE_BEARING)
        }
    }
    summary <- do.call(rbind, rows)
    utils::write.table(summary, file.path(out_dir, sprintf("%s-%s-summary.tsv", sim, name)), sep = "\t", quote = FALSE,
        row.names = FALSE, na = "")
    bad <- summary[summary$verdict %in% c("outside", "differ") & summary$rule_bearing %in% TRUE, , drop = FALSE]
    for (i in seq_len(nrow(bad))) cat(sprintf("::error::%s %s: cell %s, %s %s (max rel %s)\n", sim, name, bad$cell[[i]],
        bad$column[[i]], bad$verdict[[i]], format(bad$max_rel[[i]])))
    if (nrow(bad)) failed <<- TRUE
    cat(sprintf("%s %s: %d cell(s) gated, %d without reference; rule-bearing values outside or differing: %d\n", sim, name,
        length(unique(summary$cell[summary$verdict != "no reference"])), length(unreferenced), nrow(bad)))
}
reading("gate", c("1", "2"))
if (sim == "s1") reading("legacy", "3")
if (failed) quit(save = "no", status = 1L)
