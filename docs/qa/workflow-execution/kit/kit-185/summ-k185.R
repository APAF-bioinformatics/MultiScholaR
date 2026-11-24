# KIT-185: REPORTED-ONLY summary of the estimands the pinned scoring does not summarise (implementation audit I03;
# owner VAL-98 (2)). Reads the committed part-A rows; no classification, no band, nothing here decides. Per cell, mask
# and quantity: replicates, mean, sd, Monte Carlo error and a two-sided 95 % t interval.
#
# Usage (from the package root): Rscript docs/qa/workflow-execution/kit/kit-185/summ-k185.R <records dir> <out dir>
# Writes <out dir>/reported-extra.tsv and <out dir>/qvalue-methods.tsv.

args <- commandArgs(TRUE)
if (length(args) != 2L) stop("usage: Rscript summ-k185.R <records dir> <out dir>")
records <- args[[1L]]
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(file.path(records, "rows"), pattern = "^a-.*[0-9]\\.tsv\\.gz$", full.names = TRUE)
A <- do.call(rbind, lapply(files, function(f) utils::read.delim(gzfile(f), stringsAsFactors = FALSE, na.strings = "NA")))
A <- A[A$variant == "base", , drop = FALSE]
keys <- c("cell", "rep", "mask")

pick <- function(corr, stratum, column) {
    x <- A[A$corr == corr & A$stratum == stratum, c(keys, column), drop = FALSE]
    names(x)[4L] <- "v"
    x
}
joinDelta <- function(a, b, sign = 1) {
    x <- merge(a, b, by = keys)
    data.frame(x[keys], v = sign * (x$v.x - x$v.y), stringsAsFactors = FALSE)
}
summ <- function(x, quantity) {
    if (is.null(x) || !nrow(x)) return(NULL)
    groups <- unique(x[c("cell", "mask")])
    do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
        v <- x$v[x$cell == groups$cell[[i]] & x$mask == groups$mask[[i]]]
        v <- v[is.finite(v)]
        n <- length(v)
        m <- if (n) mean(v) else NA_real_
        s <- if (n > 1L) stats::sd(v) else NA_real_
        se <- s / sqrt(n)
        t <- if (n > 1L) stats::qt(0.975, n - 1L) else NA_real_
        data.frame(cell = groups$cell[[i]], mask = groups$mask[[i]], quantity = quantity, n = n, mean = m, sd = s, mcse = se,
            lo95 = m - t * se, hi95 = m + t * se, stringsAsFactors = FALSE)
    }))
}
out <- list()
add <- function(x, quantity) out[[length(out) + 1L]] <<- summ(x, quantity)

# (f) features lost, by residual-dimension stratum and planted status (S and C)
for (corr in c("S", "C")) for (stratum in c("filt", "rd_lt", "rd_eq", "rd_gt")) for (column in c("lost", "lost_planted", "lost_null")) {
    add(pick(corr, stratum, column), sprintf("%s %s %s", column, corr, stratum))
}
# B2 (residual variance relative to the oracle's), the contrast mu(m) - S0(m)
for (corr in c("S", "C")) add(joinDelta(pick(corr, "all", "B2"), pick("S0", "all", "B2")), sprintf("Delta B2 %s", corr))
# the incomplete-feature recovery (median share of the removed vector in the planted span), S and C on the mask
for (corr in c("S", "C")) add(pick(corr, "inc", "inc_share"), sprintf("inc_share %s", corr))
# the canonical correlations: every rho (S, F^comp), and the rank of the removed part
rho <- A[A$corr == "S" & A$stratum == "comp" & !is.na(A$rho), c(keys, "rho"), drop = FALSE]
if (nrow(rho)) {
    parts <- strsplit(rho$rho, ";", fixed = TRUE)
    width <- max(lengths(parts))
    for (j in seq_len(width)) {
        add(data.frame(rho[keys], v = vapply(parts, function(p) if (length(p) >= j) as.numeric(p[[j]]) else NA_real_, numeric(1L))),
            sprintf("rho_%d S comp", j))
    }
}
add(pick("S", "comp", "rank"), "rank S comp")
# the comp and inc strata of B and C*: mu(m) - S0(m)
for (corr in c("S", "C")) for (stratum in c("comp", "inc")) for (column in c("B", "Cstar")) {
    add(joinDelta(pick(corr, stratum, column), pick("S0", stratum, column)), sprintf("Delta %s %s %s", column, corr, stratum))
}
# the residual-dimension strata of the DA measures, oracle-referenced as the deciding contrast
for (corr in c("S", "C")) for (stratum in c("rd_lt", "rd_eq", "rd_gt")) for (column in c("fdp", "power", "null_rate")) {
    sign <- if (identical(column, "power")) -1 else 1
    a <- joinDelta(pick(corr, stratum, column), pick("O", stratum, column))
    b <- joinDelta(pick("S0", stratum, column), pick("O0", stratum, column))
    add(joinDelta(a, b, sign), sprintf("Delta %s %s %s", column, corr, stratum))
}
# the oracle's own change, e(O, m) - e(O, none) (on the same rows), and the uncorrected arm's departure e(N, m) - e(O, m)
for (column in c("fdp", "power", "null_rate")) {
    add(joinDelta(pick("O", "filt", column), pick("O0", "filt", column)), sprintf("oracle change %s", column))
    add(joinDelta(pick("N", "filt", column), pick("O", "filt", column)), sprintf("uncorrected departure %s", column))
}
# the arms' mechanics (5.2)
for (column in c("removed", "miss_share", "miss_run_min", "miss_run_med", "miss_run_max", "n_rd_lt", "n_rd_eq", "n_rd_gt", "block_share")) {
    add(pick("mask", "all", column), sprintf("arm %s", column))
}
result <- do.call(rbind, out)
result$reported_only <- "yes: not classified, decides nothing"
utils::write.table(result, file.path(out_dir, "reported-extra.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# the q-value methods (5.6): per cell and mask, the DA calls by method, and the replicates whose S(m) and S0(m) calls
# used different methods (the deciding FDP contrast's pair)
q <- A[A$stratum == "filt" & !is.na(A$qmethod), c(keys, "corr", "qmethod"), drop = FALSE]
methods <- stats::aggregate(list(calls = q$qmethod), by = list(cell = q$cell, mask = q$mask, corr = q$corr, qmethod = q$qmethod), FUN = length)
s <- merge(q[q$corr == "S", c(keys, "qmethod")], q[q$corr == "S0", c(keys, "qmethod")], by = keys)
mismatch <- stats::aggregate(list(mismatched = s$qmethod.x != s$qmethod.y), by = list(cell = s$cell, mask = s$mask), FUN = sum)
methods <- merge(methods, mismatch, by = c("cell", "mask"), all.x = TRUE)
utils::write.table(methods, file.path(out_dir, "qvalue-methods.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("KIT-185 reported-only summary: %d rows\n", nrow(result)))
