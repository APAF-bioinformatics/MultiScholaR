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

# KIT-103 W2 (C7): the design's level order, KIT-093's validation at save, the refusal text and the design state's
# object, moved with their names from R/mod_transcriptomics.R (70af7bb4) for the shared design core and the trx
# configuration (`omicGuiConfig("transcriptomics")$design`). KIT-103 W5 deleted the R/mod_transcriptomics.R copies, so
# these are the only definitions.

#' The group level order the user entered: the assigned groups ordered by the positions of their factor levels in the
#' builder's factor list (the order they were added); a group of several factors by each factor in turn.
#' @noRd
trxDesignGroupLevels <- function(design, factors) {
    assigned <- !is.na(design$group) & design$group != ""
    groups <- unique(design$group[assigned])
    first <- design[match(groups, design$group), , drop = FALSE]
    position <- function(column) if (column %in% names(first)) match(first[[column]], factors) else rep(NA_integer_, length(groups))
    groups[order(position("factor1"), position("factor2"), position("factor3"), seq_along(groups))]
}

#' A refusal's text as the tab shows it: the message and its class.
#' @noRd
trxConditionText <- function(error) {
    text <- conditionMessage(error)
    class <- class(error)[[1L]]
    if (grepl(class, text, fixed = TRUE)) text else sprintf("%s (%s)", text, class)
}

#' KIT-093's design validation at save, by name only, against every imported count column except the runs the user
#' removed (an unassigned run is a count column outside the design), with at least two runs per group (RUVs and RUVIII
#' need replicates). Returns the lipid preflight's shape.
#' @noRd
trxDesignCheck <- function(results, columnMapping) {
    runs <- setdiff(columnMapping$sample_columns, results$excluded_runs)
    tryCatch({
        levels <- results$group_levels %||% trxDesignGroupLevels(results$design_matrix, character())
        checked <- validateTrxDesign(runs, results$design_matrix, levels)
        # KIT-128: a group's samples are its biological units whenever the declared column repeats a unit: technical
        # replicates are summed per unit at the design step (`trxSumTechnicalReplicates`), and a design whose every group
        # is one unit (`spanned_by_design`) has one sample per group, which is refused below.
        replicate <- resolveProtTechnicalReplicateBlock(as.data.frame(checked$design), "group", "tech_rep_group")
        units <- replicate$reason %in% c("technical_replicates", "spanned_by_design")
        sizes <- table(if (units) checked$group[!duplicated(replicate$block)] else checked$group)
        small <- names(sizes)[sizes < 2L]
        if (length(small)) {
            trxAbort(sprintf("group(s) with fewer than two samples: %s", paste(small, collapse = ", ")),
                c("trx_design_group_singleton", "trx_design_invalid", "trx_error"))
        }
        list(valid = TRUE, errors = character())
    }, trx_error = function(error) list(valid = FALSE, errors = trxConditionText(error)))
}

#' KIT-104 D10. The QC finalise's alignment for trx: KIT-093's design check (`validateTrxDesign`) on the count table's
#' sample columns, taken from the shell's zero-row prototype as `trxAssayInputs` takes them, and the recorded level
#' order (`args$group_levels`). It refuses what every other trx step refuses.
#' @param shell The head state's shell (one count table, possibly with no rows).
#' @return `list(valid, errors)`, the errors as the tab shows them (`trxConditionText`).
#' @noRd
trxDesignAlignment <- function(shell) {
    tryCatch({
        tables <- shell@count_data
        if (length(tables) != 1L || !is.data.frame(tables[[1L]])) {
            trxAbort("count_data must hold exactly one count table", c("trx_counts_invalid", "trx_error"))
        }
        ids <- c(shell@gene_id_column, shell@annotation_id_column)
        if (!all(ids %in% names(tables[[1L]]))) {
            trxAbort(sprintf("the count table lacks the column(s) %s", paste(setdiff(ids, names(tables[[1L]])), collapse = ", ")),
                c("trx_counts_invalid", "trx_error"))
        }
        validateTrxDesign(setdiff(names(tables[[1L]]), ids), shell@design_matrix, shell@args$group_levels, shell@sample_id,
            shell@group_id)
        list(valid = TRUE, errors = character())
    }, trx_error = function(error) list(valid = FALSE, errors = trxConditionText(error)))
}

#' Sum the technical replicates of a design per biological unit (KIT-128, decisions D1 and D2): the documented
#' treatment of library-level technical replicates for the edgeR quasi-likelihood test, which takes no block (edgeR
#' User's Guide 4.9.1, limma User's Guide 18.2.5: `sumTechReps` before filtering, normalisation and the fit). The rule
#' is proteomics' (`resolveProtTechnicalReplicateBlock`): a design it does not call technical is returned unchanged.
#' @param count_data The builder's count tables (one table: id columns and one column per run).
#' @param design The builder's design (`Run`, `group` and the technical replicate column).
#' @param technical_replicate_id The design column naming each run's biological unit.
#' @return A list: `count_data` (the run columns replaced by one summed column per unit, in first-occurrence order),
#'   `design` (one row per unit, the unit's first run's row with `Run` the unit's name) and `summed` (unit to its runs,
#'   or NULL when nothing was summed).
#' @noRd
trxSumTechnicalReplicates <- function(count_data, design, technical_replicate_id = "tech_rep_group") {
    replicate <- resolveProtTechnicalReplicateBlock(as.data.frame(design), "group", technical_replicate_id)
    if (!isTRUE(replicate$technical)) return(list(count_data = count_data, design = design, summed = NULL))
    requireOptionalPackage("edgeR", "summing the technical replicates")
    runs <- as.character(design$Run)
    units <- unique(replicate$block)
    count_data <- lapply(count_data, function(table) {
        summed <- edgeR::sumTechReps(as.matrix(table[, runs, drop = FALSE]), ID = replicate$block)
        cbind(table[, setdiff(names(table), runs), drop = FALSE], as.data.frame(summed[, units, drop = FALSE], check.names = FALSE))
    })
    design <- design[!duplicated(replicate$block), , drop = FALSE]
    design$Run <- units
    rownames(design) <- NULL
    logger::log_info(sprintf("Technical replicates: the counts of %d runs summed into %d biological samples", length(runs), length(units)))
    list(count_data = count_data, design = design, summed = split(runs, factor(replicate$block, levels = units)))
}

#' The design state's object: the builder's count table and design, the recorded level order, the chosen id and name
#' columns, and the builder's `tech_rep_group` as the technical replicates. KIT-128: declared technical replicates are
#' summed per biological unit here (`trxSumTechnicalReplicates`), so every later step works on the units; the runs that
#' were summed are recorded in `args$technical_replicates_summed`.
#' @noRd
trxDesignObject <- function(results, columnMapping) {
    summed <- trxSumTechnicalReplicates(results$data_cln, results$design_matrix)
    createTranscriptomicsAssayData(summed$count_data, summed$design, results$group_levels,
        roles = list(gene_id_column = columnMapping$lipid_id_col, annotation_id_column = columnMapping$annotation_col,
            technical_replicate_id = "tech_rep_group"),
        args = if (is.null(summed$summed)) list() else list(technical_replicates_summed = summed$summed))
}
