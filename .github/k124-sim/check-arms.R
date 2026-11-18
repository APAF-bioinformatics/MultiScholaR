# KIT-124 confirmation round on hosted runners: the shipped-kernel arms of a simulator's out dir, checked after the
# upload (so that the rows of a failing job are kept).
#
# Usage: Rscript --vanilla check-arms.R <sim: s1..s4> <out dir> <cell> <arms requested (K124_ARMS)>
#
# The shipped-kernel arms are 34 to 41. Of those requested, the ones the simulator forms in <cell> (the rules below,
# as each simulator's arm list has them) must have rows in every run of the table; a requested arm the simulator does
# not form in the cell is listed as such. Exit status 1 when a formed arm is missing from a run or has a crashed row,
# so that a crash of the package under test fails the job instead of being scored silently. Prints rows per arm and
# status and the seconds per run of each shipped arm.

FORMED <- list(
    s1 = function(cell) c("34", "35", if (cell != "C4") "36", "37", if (cell %in% c("C6", "C8")) "38", "39", "40", "41"),
    s2 = function(cell) c("34", "35", if (cell != "C4") "36", "37", "39", "40", "41"),
    s3 = function(cell) c("34", "35", "36", "37", "39", "40", "41"),
    s4 = function(cell) c("34", "35", "36", "37", if (cell %in% c("C6", "C8")) "38", "39", "40", "41")
)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) stop("usage: Rscript --vanilla check-arms.R <sim> <out dir> <cell> <arms>")
sim <- arguments[[1L]]; out_dir <- arguments[[2L]]; cell <- arguments[[3L]]
requested <- intersect(strsplit(arguments[[4L]], ",", fixed = TRUE)[[1L]], as.character(34:41))
formed <- intersect(requested, FORMED[[sim]](cell))
problems <- character()
not_formed <- setdiff(requested, formed)
if (length(not_formed)) cat(sprintf("requested, not formed by %s in %s (by the simulator's rule): %s\n", sim, cell, paste(not_formed, collapse = ", ")))
if (!length(formed)) {
    cat("no shipped-kernel arm is formed here: nothing to check\n")
    quit(save = "no", status = 0L)
}
rows <- utils::read.delim(file.path(out_dir, sprintf("%s-runs.tsv", sim)), colClasses = "character", quote = "")
runs <- sort(unique(as.integer(rows$run)))
for (arm in formed) {
    mine <- rows[rows$arm == arm, , drop = FALSE]
    missing_runs <- setdiff(runs, as.integer(mine$run))
    if (length(missing_runs)) problems <- c(problems, sprintf("arm %s has no row in run(s) %s", arm, paste(missing_runs, collapse = ", ")))
    crashed <- mine[mine$status == "crashed", , drop = FALSE]
    if (nrow(crashed)) problems <- c(problems, sprintf("arm %s crashed in %d row(s) (runs %s)", arm, nrow(crashed),
        paste(sort(unique(as.integer(crashed$run))), collapse = ", ")))
    first_subset <- !duplicated(mine[c("run", "arm")])
    seconds <- as.numeric(mine$seconds[first_subset])
    cat(sprintf("arm %s: %d rows; status %s; seconds per run mean %.2f, max %.2f\n", arm, nrow(mine),
        paste(sprintf("%s %d", names(table(mine$status)), as.integer(table(mine$status))), collapse = ", "),
        mean(seconds, na.rm = TRUE), suppressWarnings(max(seconds, na.rm = TRUE))))
}
if (length(problems)) {
    cat(sprintf("::error::%s %s: %s\n", sim, cell, problems), sep = "")
    quit(save = "no", status = 1L)
}
cat("shipped-kernel arms: every formed arm has a row in every run and none crashed\n")
