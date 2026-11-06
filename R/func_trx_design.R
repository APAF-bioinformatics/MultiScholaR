# KIT-093: transcriptomics design validation. Labels attach by run ID, never by position. The validated design keeps the
# design-file order; its sample column is `Run` and its group column `group` (a factor with the declared levels), whatever
# the file calls them, because the notebook's contrast strings and model frames depend on those names (rmd:529-533).

#' Validate a transcriptomics design against the count sample columns.
#' @param sample_columns Character vector of the count table's sample columns, in their file order.
#' @param design Data frame, one row per run.
#' @param group_levels The declared group levels, in order: two or more, each non-empty (KIT-093 revision 1.2, the user's
#'   decision of 2026-11-06: any number of groups, DA as pairwise contrasts, as the other omics).
#' @param sample_id,group_id Column names of the run ID and the group in `design`.
#' @return A list: `design` (design-file order, columns renamed to `Run` and `group`, `group` a factor), `group` (that
#'   factor) and `count_index` (so that `counts[, count_index]` is in design order, as `count_matrix[, design_matrix$Run]`
#'   at rmd:152). Refusals inherit `trx_design_invalid`.
#' @noRd
validateTrxDesign <- function(sample_columns, design, group_levels, sample_id = "Run", group_id = "group") {
    refuse <- function(message, fault) trxAbort(message, c(fault, "trx_design_invalid", "trx_error"))
    if (!is.data.frame(design) || !nrow(design) || !all(c(sample_id, group_id) %in% names(design))) {
        refuse(sprintf("the design is missing, empty or lacks the columns '%s' and '%s'", sample_id, group_id), "trx_design_missing")
    }
    runs <- as.character(design[[sample_id]])
    if (anyDuplicated(runs) || anyDuplicated(sample_columns)) {
        refuse(sprintf("duplicated run(s): %s", paste(unique(c(runs[duplicated(runs)], sample_columns[duplicated(sample_columns)])),
            collapse = ", ")), "trx_design_runs_duplicated")
    }
    extra <- setdiff(sample_columns, runs)
    if (length(extra)) refuse(sprintf("count column(s) not in the design: %s", paste(extra, collapse = ", ")), "trx_design_runs_extra")
    missing <- setdiff(runs, sample_columns)
    if (length(missing)) refuse(sprintf("design run(s) not in the counts: %s", paste(missing, collapse = ", ")), "trx_design_runs_missing")
    if (length(group_levels) < 2L || anyDuplicated(group_levels)) {
        trxAbort(sprintf("at least two distinct declared groups are needed, not %d", length(unique(group_levels))), c("trx_design_invalid", "trx_error"))
    }
    labels <- as.character(design[[group_id]])
    unknown <- is.na(labels) | !labels %in% group_levels
    if (any(unknown)) refuse(sprintf("group label(s) outside the declared levels: %s", paste(unique(labels[unknown]), collapse = ", ")),
        "trx_design_level_unknown")
    empty <- setdiff(group_levels, labels)
    if (length(empty)) refuse(sprintf("empty group(s): %s", paste(empty, collapse = ", ")), "trx_design_group_empty")
    names(design)[match(c(sample_id, group_id), names(design))] <- c("Run", "group")
    design[["Run"]] <- runs
    design[["group"]] <- factor(labels, levels = group_levels)
    list(design = design, group = design[["group"]], count_index = match(runs, sample_columns))
}
