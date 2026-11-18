# KIT-125: the per-omic tables and the verdict of ticket revision 1.2, from the run folders.
#   Rscript k125_aggregate.R <folder holding <omic>/runs.tsv> <kit-124/scored folder> <out tsv>
# Per omic, split type and design: splits, statuses, the mean null rate at 0.05 with its range and its 5th and 95th
# percentiles over the splits (no confidence bound: the splits share samples, revision 1.1), the share of splits
# with any call, and beside them the same design in KIT-124's matching simulation cells (C4, the global null; C5, the
# factor tied to the groups at 0.3, for the batch-aligned splits).
# Verdict (revision 1.2; the "passed"/"failed" table is the declaration's, fixed before the runs):
#   within = mean null rate <= 0.075 and (any-call share <= 2 x no-RUV's share on the same splits, or <= 0.10)
#   a passed design not within on some split type          -> contradiction
#   a failed design within on every split type of its data -> contradiction
#   fewer than 100 scored splits                           -> not informative
#   otherwise                                              -> consistent
args <- commandArgs(trailingOnly = TRUE)
runs_dir <- args[[1L]]; scored_dir <- args[[2L]]; out <- args[[3L]]
SIMULATOR <- c(classic = "s1", dia_limpa = "s2", trx = "s3", metab = "s4")
PASSED <- rbind(
    data.frame(omic = "classic", arm = c("1", "3", "7"), kit124 = c("passed", "failed", "passed")),
    data.frame(omic = "dia_limpa", arm = c("1", "3", "4", "7"), kit124 = c("failed", "failed", "failed", "passed")),
    data.frame(omic = "trx", arm = c("1", "3", "7"), kit124 = c("passed", "failed", "passed")),
    data.frame(omic = "metab", arm = c("1", "3", "7"), kit124 = c("passed", "failed", "passed")))
rows <- list()
for (omic in names(SIMULATOR)) {
    file <- file.path(runs_dir, omic, "runs.tsv")
    if (!file.exists(file)) next
    runs <- utils::read.delim(file, stringsAsFactors = FALSE, quote = "", colClasses = c(arm = "character"))
    cells <- utils::read.delim(file.path(scored_dir, paste0(SIMULATOR[[omic]], "-core"), "k124-cells.tsv"), stringsAsFactors = FALSE,
        colClasses = c(arm = "character"))
    cells <- cells[cells$subset == "full", ]
    sim <- function(cell, arm, column) { v <- cells[[column]][cells$cell == cell & cells$arm == arm]; if (length(v) == 1L) v else NA_real_ }
    for (type in unique(runs$split_type)) {
        of_type <- runs[runs$split_type == type, ]
        base <- of_type[of_type$arm == "1", ]
        base_any <- mean(base$n_calls[is.finite(base$n_calls)] > 0)
        for (arm in unique(of_type$arm)) {
            r <- of_type[of_type$arm == arm, ]
            scored <- r[is.finite(r$null_p05), ]
            any_share <- mean(scored$n_calls > 0)
            mean_null <- mean(scored$null_p05)
            rows[[length(rows) + 1L]] <- data.frame(omic = omic, split_type = type, arm = arm, design = r$design[[1L]],
                splits = nrow(r), ok = sum(r$status == "ok"), refused = sum(r$status == "refused"), crashed = sum(r$status == "crashed"),
                scored = nrow(scored), mean_k = mean(scored$k, na.rm = TRUE), k_zero_share = mean(scored$k == 0, na.rm = TRUE),
                mean_features = mean(scored$n_features), null05_mean = mean_null, null05_min = min(scored$null_p05),
                null05_p05 = stats::quantile(scored$null_p05, 0.05, names = FALSE), null05_median = stats::median(scored$null_p05),
                null05_p95 = stats::quantile(scored$null_p05, 0.95, names = FALSE), null05_max = max(scored$null_p05),
                null01_mean = mean(scored$null_p01), anycall_splits = sum(scored$n_calls > 0), anycall_share = any_share,
                no_ruv_anycall_share = base_any, median_calls_when_any = if (any(scored$n_calls > 0)) stats::median(scored$n_calls[scored$n_calls > 0]) else NA_real_,
                within_null = mean_null <= 0.075, within_anycall = any_share <= 2 * base_any || any_share <= 0.10,
                sim_C4_null05 = sim("C4", arm, "null05"), sim_C4_anycall = sim("C4", arm, "anycall_share"),
                sim_C5_null05 = sim("C5", arm, "null05"), stringsAsFactors = FALSE)
        }
    }
}
table <- do.call(rbind, rows)
table$within <- table$within_null & table$within_anycall
table$kit124 <- PASSED$kit124[match(paste(table$omic, table$arm), paste(PASSED$omic, PASSED$arm))]
table$verdict <- NA_character_
for (key in unique(paste(table$omic, table$arm))) {
    at <- which(paste(table$omic, table$arm) == key)
    verdict <- if (any(table$scored[at] < 100L)) "not informative" else if (table$kit124[at[[1L]]] == "passed") {
        if (all(table$within[at])) "consistent" else "contradiction"
    } else if (all(table$within[at])) "contradiction" else "consistent"
    table$verdict[at] <- verdict
}
numeric <- vapply(table, is.double, logical(1L))
table[numeric] <- lapply(table[numeric], signif, digits = 4L)
utils::write.table(table, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
print(table[c("omic", "split_type", "arm", "scored", "null05_mean", "null05_min", "null05_max", "anycall_share", "within", "kit124", "verdict")], row.names = FALSE)
