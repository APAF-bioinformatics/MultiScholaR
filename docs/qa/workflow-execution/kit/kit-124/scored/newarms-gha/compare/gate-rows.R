# Per-row reading of the equivalence comparison, and the gate as declared and as read under decision M19 (same rule as packaging/verify/compare-results.R: rel 1e-6, abs 1e-9).
args <- commandArgs(trailingOnly = TRUE); dir <- args[[1L]]
out <- list(); summ <- list()
for (s in c("s1", "s2", "s3", "s4")) {
    a <- readRDS(file.path(dir, sprintf("%s-laptop.rds", s))); b <- readRDS(file.path(dir, sprintf("%s-runner.rds", s)))
    stopifnot(identical(names(a), names(b)))
    for (cell in names(a)) {
        x <- a[[cell]]; y <- b[[cell]]
        keys <- c("simulator", "cell", "cell_class", "run", "arm", "subset", "status", "seed", "k", "n_features", "n_true", "n_calls")
        discrete <- sum(vapply(keys, \(col) sum(!(x[[col]] == y[[col]] | (is.na(x[[col]]) & is.na(y[[col]]))) | xor(is.na(x[[col]]), is.na(y[[col]])), na.rm = TRUE), numeric(1L)))
        bad <- rep(FALSE, nrow(x))
        for (m in c("null_p05", "null_p01", "fdp", "power", "bias", "rmse")) {
            u <- x[[m]]; v <- y[[m]]; fin <- is.finite(u) & is.finite(v)
            d <- abs(u - v); mx <- pmax(abs(u), abs(v)); ok <- !fin | d <= 1e-6 * mx | d <= 1e-9
            nak <- identical(is.na(u), is.na(v))
            if (any(!ok)) out[[length(out) + 1L]] <- data.frame(cell = cell, run = x$run[!ok], arm = x$arm[!ok], subset = x$subset[!ok],
                measure = m, laptop = u[!ok], runner = v[!ok], abs_diff = d[!ok], rel_diff = (d / mx)[!ok])
            bad <- bad | !ok
            if (!nak) stop("NA positions differ: ", cell, " ", m)
        }
        summ[[length(summ) + 1L]] <- data.frame(simulator = toupper(s), cell = sub("^.* ", "", cell), rows_compared = nrow(x),
            runs = length(unique(x$run)), arms = paste(sort(unique(x$arm)), collapse = ","), subsets = paste(sort(unique(x$subset)), collapse = ","),
            discrete_differences = discrete, rows_outside_rule = sum(bad),
            gate_as_declared = if (discrete == 0 && !any(bad)) "pass" else "fail",
            measures_outside = paste(unique(unlist(lapply(out, \(o) if (identical(o$cell[[1L]], cell)) o$measure[[1L]]))), collapse = ","),
            stringsAsFactors = FALSE)
    }
}
summ <- do.call(rbind, summ); out <- if (length(out)) do.call(rbind, out) else data.frame()
# Decision M19 (user, 2026-11-18, design-review-20261118.md section 6): the gate is read on the values the rules use
# (null rates, false discovery proportion, power, k, counts); rows outside the tolerance on bias and rmse only are
# accepted and reported. The gate as first declared (every value) stays in gate_as_declared.
summ$gate_as_read_M19 <- ifelse(summ$gate_as_declared == "pass", "pass", ifelse(summ$discrete_differences == 0 &
    vapply(strsplit(summ$measures_outside, ","), \(m) all(m %in% c("bias", "rmse")), logical(1L)),
    "pass: outside on bias/rmse only, accepted (M19); rows in gate-rows-outside.tsv", "fail"))
write.table(summ, file.path(dir, "gate-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out, file.path(dir, "gate-rows-outside.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE); options(width = 200, digits = 10); print(out, row.names = FALSE)
