# KIT-124 confirmation round on hosted runners: a shipped-kernel arm against the simulator arm it replaces, on the same
# data sets, under the fixed rule (packaging/verify/compare-results.R of the package under test; this script decides
# nothing).
#
# Usage: Rscript --vanilla compare-arms.R <sim: s1..s4> <tables dir> <label> <pairs "new:base,..."> <compare-results.R> <out dir>
#
# Reads s<i>-runs.tsv from every artifact directory k124-<label>-<sim>-* under <tables dir>. Per pair, the rows of the
# two arms are keyed by cell, run and subset (a key on one side only is listed and left out); <out dir>/<sim>-<new>-<base>-a.rds
# (the simulator arm) and -b.rds (the shipped arm) hold one data frame per cell ("S1 C6 34~32") with the key, the
# status and every measure column (`seconds`, timing, left out); whole-number columns as integers so the rule compares
# them exactly, as compare-laptop.R stores them. Then compare-results.R writes <out dir>/<sim>-<new>-<base>-compare.tsv,
# and <out dir>/<sim>-summary.tsv has one row per pair and column: the worst verdict over the cells and the largest
# relative difference. The rule-bearing columns (decision M19) are status, k, the counts, null_p05, null_p01, fdp and
# power; bias and rmse are reported.

DISCRETE <- c("seed", "k", "n_features", "n_true", "n_calls")
MEASURES <- c("null_p05", "null_p01", "fdp", "power", "bias", "rmse")
RULE_BEARING <- c("status", "k", "n_features", "n_true", "n_calls", "null_p05", "null_p01", "fdp", "power")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 6L) stop("usage: Rscript --vanilla compare-arms.R <sim> <tables dir> <label> <pairs> <compare-results.R> <out dir>")
sim <- arguments[[1L]]; tables_dir <- arguments[[2L]]; label <- arguments[[3L]]
pairs <- strsplit(strsplit(arguments[[4L]], ",", fixed = TRUE)[[1L]], ":", fixed = TRUE)
rule <- arguments[[5L]]; out_dir <- arguments[[6L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parts <- list.files(tables_dir, pattern = sprintf("^k124-%s-%s-", label, sim), full.names = TRUE)
files <- file.path(parts, sprintf("%s-runs.tsv", sim))
files <- files[file.exists(files)]
summary_path <- file.path(out_dir, sprintf("%s-summary.tsv", sim))
if (!length(files)) {
    utils::write.table(data.frame(simulator = toupper(sim), pair = NA, column = NA, cells = 0L, verdict = "no table",
        max_rel = NA, rule_bearing = NA), summary_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    cat(sprintf("%s: no result table among the artifacts (a failed job uploads nothing)\n", sim))
    quit(save = "no", status = 0L)
}
rows <- do.call(rbind, lapply(files, utils::read.delim, colClasses = "character", quote = "", na.strings = "NA"))
key <- function(x) paste(x$cell, x$run, x$subset, sep = "\t")
whole <- function(x) { v <- suppressWarnings(as.numeric(x)); all(is.na(v) == is.na(x)) && all(v[!is.na(v)] == floor(v[!is.na(v)])) }

summary <- list()
for (pair in pairs) {
    new <- pair[[1L]]; base <- pair[[2L]]
    a <- rows[rows$arm == base, , drop = FALSE]
    b <- rows[rows$arm == new, , drop = FALSE]
    common <- intersect(key(a), key(b))
    only <- c(setdiff(key(a), common), setdiff(key(b), common))
    if (length(only)) cat(sprintf("%s %s~%s: %d key(s) on one side only, left out\n", sim, new, base, length(only)))
    a <- a[match(common, key(a)), , drop = FALSE]
    b <- b[match(common, key(b)), , drop = FALSE]
    as_integer <- vapply(DISCRETE, \(column) whole(a[[column]]) && whole(b[[column]]), logical(1L))
    typed <- function(x) {
        out <- x[c("cell", "run", "subset", "status")]
        out$run <- as.integer(out$run)
        for (column in DISCRETE) out[[column]] <- if (as_integer[[column]]) as.integer(x[[column]]) else as.numeric(x[[column]])
        for (column in MEASURES) out[[column]] <- as.numeric(x[[column]])
        out <- out[order(out$cell, out$run, out$subset), , drop = FALSE]
        rownames(out) <- NULL
        out
    }
    byCell <- function(x) {
        t <- typed(x)
        cells <- split(t, t$cell)
        cells <- lapply(cells, \(part) { rownames(part) <- NULL; part })
        stats::setNames(cells, sprintf("%s %s %s~%s", toupper(sim), names(cells), new, base))
    }
    stem <- file.path(out_dir, sprintf("%s-%s-%s", sim, new, base))
    saveRDS(byCell(a), paste0(stem, "-a.rds"))
    saveRDS(byCell(b), paste0(stem, "-b.rds"))
    status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", rule, paste0(stem, "-a.rds"), paste0(stem, "-b.rds"),
        paste0(stem, "-compare.tsv")), stdout = FALSE)
    if (status != 0L || !file.exists(paste0(stem, "-compare.tsv"))) stop("compare-results.R failed for ", basename(stem))
    compared <- utils::read.delim(paste0(stem, "-compare.tsv"), stringsAsFactors = FALSE, na.strings = "")
    column <- sub("^[^$]*\\$", "", compared$path)
    compared <- compared[column %in% c(RULE_BEARING, "bias", "rmse", "seed", "cell", "run", "subset"), , drop = FALSE]
    column <- sub("^[^$]*\\$", "", compared$path)
    for (name in unique(column)) {
        mine <- compared[column == name, , drop = FALSE]
        verdict <- if (all(mine$verdict == "agree")) "agree" else if (any(mine$verdict == "differ")) "differ" else "outside"
        summary[[length(summary) + 1L]] <- data.frame(simulator = toupper(sim), pair = paste0(new, "~", base), column = name,
            cells = nrow(mine), rows = length(common), verdict = verdict,
            max_rel = if (any(is.finite(mine$max_rel))) max(mine$max_rel, na.rm = TRUE) else NA_real_,
            outside = sum(mine$outside, na.rm = TRUE), rule_bearing = name %in% RULE_BEARING, stringsAsFactors = FALSE)
    }
}
summary <- do.call(rbind, summary)
utils::write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("%s: %d pair(s) compared; rule-bearing columns outside the rule or differing: %d\n", sim, length(pairs),
    sum(summary$rule_bearing & summary$verdict != "agree")))
