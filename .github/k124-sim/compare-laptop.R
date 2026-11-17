# KIT-124 on hosted runners: the same rows from the laptop's and a runner's out dir, saved for the fixed comparison
# rule (packaging/verify/compare-results.R on performance-v2; this script decides nothing).
#
# Usage: Rscript --vanilla compare-laptop.R <sim: s1..s4> <laptop out dir> <runner out dir> <arms> <a.rds> <b.rds>
#
# Reads s<i>-runs.tsv of both. The rows compared are the runner's cells and runs, of the arms given (a comma list),
# keyed by cell, run, arm and subset. <a.rds> (laptop) and <b.rds> (runner) each hold one data frame per cell
# ("S1 C0"), rows in key order, with the key, the status and every measure column; `seconds` (timing) is left out.
# Columns whose values are whole numbers in both tables are stored as integers, so the rule compares them exactly;
# the others as doubles. A key present on one side only stops the script.

DISCRETE <- c("seed", "k", "n_features", "n_true", "n_calls")
MEASURES <- c("null_p05", "null_p01", "fdp", "power", "bias", "rmse")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) stop("usage: Rscript --vanilla compare-laptop.R <sim> <laptop dir> <runner dir> <arms> <a.rds> <b.rds>")
sim <- arguments[[1L]]
arms <- strsplit(arguments[[4L]], ",", fixed = TRUE)[[1L]]
readRuns <- function(dir) {
    utils::read.delim(file.path(dir, sprintf("%s-runs.tsv", sim)), colClasses = "character", quote = "", fill = TRUE,
        na.strings = "NA")
}
laptop <- readRuns(arguments[[2L]])
runner <- readRuns(arguments[[3L]])
stopifnot(identical(names(laptop), names(runner)), all(c(DISCRETE, MEASURES, "seconds") %in% names(runner)))
key <- function(x) paste(x$cell, x$run, x$arm, x$subset, sep = "\t")
runner <- runner[runner$arm %in% arms, ]
laptop <- laptop[key(laptop) %in% key(runner), ]
if (!setequal(key(laptop), key(runner)) || anyDuplicated(key(laptop)) || anyDuplicated(key(runner))) {
    stop(sprintf("the keys differ: laptop %d rows, runner %d rows, %d of the runner's not on the laptop", nrow(laptop),
        nrow(runner), sum(!key(runner) %in% key(laptop))))
}
ordered <- function(x) x[order(x$cell, as.integer(x$run), x$arm, x$subset), ]
laptop <- ordered(laptop); runner <- ordered(runner)
whole <- function(x) { v <- suppressWarnings(as.numeric(x)); all(is.na(v) == is.na(x)) && all(v[!is.na(v)] == floor(v[!is.na(v)])) }
typed <- function(x, as_integer) {
    out <- x[c("simulator", "cell", "cell_class", "run", "arm", "subset", "status")]
    out$run <- as.integer(out$run)
    for (column in DISCRETE) out[[column]] <- if (as_integer[[column]]) as.integer(x[[column]]) else as.numeric(x[[column]])
    for (column in MEASURES) out[[column]] <- as.numeric(x[[column]])
    rownames(out) <- NULL
    out
}
as_integer <- vapply(DISCRETE, \(column) whole(laptop[[column]]) && whole(runner[[column]]), logical(1L))
if (!all(as_integer)) cat("not whole numbers, compared as doubles:", paste(DISCRETE[!as_integer], collapse = ", "), "\n")
byCell <- function(x) {
    parts <- split(typed(x, as_integer), x$cell)
    parts <- lapply(parts, \(part) { rownames(part) <- NULL; part })
    stats::setNames(parts, paste(toupper(sim), names(parts)))
}
saveRDS(byCell(laptop), arguments[[5L]])
saveRDS(byCell(runner), arguments[[6L]])
cat(sprintf("%s: %d rows per side (cells %s; runs 1 to %d; arms %s)\n", toupper(sim), nrow(runner),
    paste(unique(runner$cell), collapse = ", "), max(as.integer(runner$run)), paste(arms, collapse = ", ")))
