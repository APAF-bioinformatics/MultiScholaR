# KIT-124 exactness round: per-row reading of the equivalence comparison (rule of packaging/verify/compare-results.R:
# relative 1e-6, absolute floor 1e-9, discrete values exact), the gate as declared and as read under decision M19.
# Same reading as scored/newarms-gha/compare/gate-rows.R, generalised to named pairs of rds files written by
# .github/k124-sim/compare-laptop.R (sim-k124-exact). Reads tables only; decides nothing beyond the fixed rule.
# Usage: Rscript --vanilla gate-rows.R <dir> <tag>=<reference.rds>,<runner.rds>[,<reference name>] ...
args <- commandArgs(trailingOnly = TRUE); dir <- args[[1L]]; pairs <- args[-1L]
keys <- c("simulator", "cell", "cell_class", "run", "arm", "subset", "status", "seed", "k", "n_features", "n_true", "n_calls")
rule_measures <- c("null_p05", "null_p01", "fdp", "power")
out <- list(); summ <- list()
for (p in pairs) {
    tag <- sub("=.*$", "", p); f <- strsplit(sub("^[^=]*=", "", p), ",")[[1L]]
    a <- readRDS(file.path(dir, f[[1L]])); b <- readRDS(file.path(dir, f[[2L]])); ref <- if (length(f) > 2L) f[[3L]] else "laptop"
    stopifnot(identical(names(a), names(b)))
    for (cell in names(a)) for (arm in sort(unique(a[[cell]]$arm))) {
        x <- a[[cell]][a[[cell]]$arm == arm, ]; y <- b[[cell]][b[[cell]]$arm == arm, ]
        stopifnot(identical(x$run, y$run), identical(x$subset, y$subset))
        same <- function(u, v) (u == v | (is.na(u) & is.na(v))) & !xor(is.na(u), is.na(v))
        dk <- vapply(keys, \(col) sum(!same(x[[col]], y[[col]]), na.rm = TRUE), numeric(1L))
        bad <- rep(FALSE, nrow(x)); badm <- character()
        maxrel <- 0
        for (m in c(rule_measures, "bias", "rmse")) {
            u <- x[[m]]; v <- y[[m]]; fin <- is.finite(u) & is.finite(v)
            d <- abs(u - v); mx <- pmax(abs(u), abs(v)); ok <- !fin | d <= 1e-6 * mx | d <= 1e-9
            if (!identical(is.na(u), is.na(v))) { ok <- ok & (is.na(u) == is.na(v)); }
            rel <- ifelse(fin & mx > 0, d / mx, 0); maxrel <- max(maxrel, rel, na.rm = TRUE)
            if (any(!ok)) {
                out[[length(out) + 1L]] <- data.frame(comparison = tag, reference = ref, cell = sub("^.* ", "", cell),
                    run = x$run[!ok], arm = arm, subset = x$subset[!ok], measure = m, reference_value = u[!ok],
                    runner_value = v[!ok], abs_diff = d[!ok], rel_diff = (d / mx)[!ok])
                badm <- c(badm, m)
            }
            bad <- bad | !ok
        }
        rule_bad <- any(badm %in% rule_measures)
        summ[[length(summ) + 1L]] <- data.frame(comparison = tag, reference = ref, simulator = sub(" .*$", "", cell),
            cell = sub("^.* ", "", cell), arm = arm, rows_compared = nrow(x), runs = length(unique(x$run)),
            subsets = paste(sort(unique(x$subset)), collapse = ","), discrete_differences = sum(dk),
            discrete_columns = paste(names(dk)[dk > 0], collapse = ","), rows_outside_rule = sum(bad),
            measures_outside = paste(unique(badm), collapse = ","), largest_rel_diff = signif(maxrel, 3),
            gate_as_declared = if (sum(dk) == 0 && !any(bad)) "pass" else "fail",
            gate_as_read_M19 = if (sum(dk) == 0 && !rule_bad) {
                if (any(bad)) "pass: outside on bias/rmse only, accepted (M19)" else "pass" } else "fail",
            stringsAsFactors = FALSE)
    }
}
summ <- do.call(rbind, summ); out <- if (length(out)) do.call(rbind, out) else data.frame()
write.table(summ, file.path(dir, "gate-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(out, file.path(dir, "gate-rows-outside.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250); print(summ, row.names = FALSE)
cat("\nrows outside the rule:", nrow(out), "\n"); if (nrow(out)) print(utils::head(out, 60), row.names = FALSE)
