# KIT-093: the transcriptomics count class, in the lipid/metab shape (R/func_lipid_s4_assay_data_class.R). Recorded
# divergence from conventions/R-bioinformatics.md:28-33: like the lipid and metab classes it does not inherit
# SummarizedExperiment; its container contract is the kit's class declaration (TRX-005 decision 1). Validity runs the
# count policy (validateTrxCounts) and the design checks (validateTrxDesign), which refuse by class.

#' Transcriptomics count data
#'
#' @slot count_data Named list holding one wide genes-by-samples count table with the gene ID and gene name columns.
#' @slot gene_id_column,annotation_id_column The gene ID and gene name (annotation) columns of the count table.
#' @slot design_matrix The design, one row per run.
#' @slot sample_id,group_id,technical_replicate_id Design columns: run ID, group, technical-replicate set (or NA).
#' @slot args Parameters; `args$group_levels` holds the declared group level order every kernel builds its factor from.
#' @noRd
setClass("TranscriptomicsAssayData",
    slots = c(count_data = "list", gene_id_column = "character", annotation_id_column = "character",
        design_matrix = "data.frame", sample_id = "character", group_id = "character",
        technical_replicate_id = "character", args = "list"),
    prototype = list(count_data = list(), gene_id_column = "GeneID", annotation_id_column = "GeneName",
        design_matrix = data.frame(), sample_id = "Run", group_id = "group", technical_replicate_id = NA_character_,
        args = list()),
    validity = function(object) {
        trxAssayInputs(object)
        TRUE
    })

#' The design-ordered count matrix and the validated design of a TranscriptomicsAssayData object.
#' Runs the count policy and the design checks (both refuse by class).
#' @param object A TranscriptomicsAssayData object.
#' @return A list: `counts` (numeric genes-by-runs matrix, rows named by gene ID, columns in design order, built as
#'   rmd:150-152), `design`, `group` and `count_index` (from validateTrxDesign).
#' @noRd
trxAssayInputs <- function(object) {
    tables <- object@count_data
    if (length(tables) != 1L || !is.data.frame(tables[[1L]])) {
        trxAbort("count_data must hold exactly one count table", c("trx_counts_invalid", "trx_error"))
    }
    table <- tables[[1L]]
    ids <- c(object@gene_id_column, object@annotation_id_column)
    if (!all(ids %in% names(table))) {
        trxAbort(sprintf("the count table lacks the column(s) %s", paste(setdiff(ids, names(table)), collapse = ", ")), c("trx_counts_invalid", "trx_error"))
    }
    sample_columns <- setdiff(names(table), ids)
    checked <- validateTrxDesign(sample_columns, object@design_matrix, object@args$group_levels, object@sample_id, object@group_id)
    counts <- as.matrix(table[, sample_columns])
    rownames(counts) <- table[[object@gene_id_column]]
    validateTrxCounts(counts)
    c(list(counts = counts[, checked$count_index]), checked)
}

#' Create a TranscriptomicsAssayData object.
#' @param count_data Named list of one count table (gene ID, gene name and one numeric column per run).
#' @param design_matrix The design data frame.
#' @param group_levels The declared group levels, in order (stored in `args$group_levels`).
#' @param roles Optional named list overriding `gene_id_column` ("GeneID"), `annotation_id_column` ("GeneName"),
#'   `sample_id` ("Run"), `group_id` ("group") and `technical_replicate_id` (NA).
#' @param args Further parameters.
#' @return A valid TranscriptomicsAssayData object; refusals are classed (`trx_counts_*`, `trx_design_*`).
#' @noRd
createTranscriptomicsAssayData <- function(count_data, design_matrix, group_levels, roles = list(), args = list()) {
    defaults <- list(gene_id_column = "GeneID", annotation_id_column = "GeneName", sample_id = "Run", group_id = "group",
        technical_replicate_id = NA_character_)
    unknown <- setdiff(names(roles), names(defaults))
    if (length(unknown)) trxAbort(sprintf("unknown role(s): %s", paste(unknown, collapse = ", ")), c("trx_design_invalid", "trx_error"))
    roles <- utils::modifyList(defaults, roles)
    args$group_levels <- as.character(group_levels)
    methods::new("TranscriptomicsAssayData", count_data = count_data, gene_id_column = roles$gene_id_column,
        annotation_id_column = roles$annotation_id_column, design_matrix = design_matrix, sample_id = roles$sample_id,
        group_id = roles$group_id, technical_replicate_id = as.character(roles$technical_replicate_id), args = args)
}
