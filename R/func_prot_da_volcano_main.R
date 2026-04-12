# ----------------------------------------------------------------------------
# writeInteractiveVolcanoPlotProteomicsMain
# ----------------------------------------------------------------------------
#' @export
writeInteractiveVolcanoPlotProteomicsMain <- function(
  da_analysis_results_list,
  theObject,
  uniprot_tbl,
  de_analysis_results_list = NULL,
  publication_graphs_dir = NULL,
  file_prefix = NULL,
  plots_format = NULL,
  args_row_id = NULL,
  da_q_val_thresh = NULL,
  de_q_val_thresh = NULL,
  gene_names_column = NULL,
  fdr_column = NULL,
  raw_p_value_column = NULL,
  log2fc_column = NULL,
  uniprot_id_column = NULL,
  display_columns = NULL
) {
  # Accept legacy de_* aliases but normalize to the current da_* path internally.
  results_list <- da_analysis_results_list
  if (is.null(results_list) && !is.null(de_analysis_results_list)) {
    results_list <- de_analysis_results_list
  }
  if (is.null(results_list)) {
    stop("A DA/DE analysis results list must be supplied.", call. = FALSE)
  }

  safe_get_slot <- function(res_list, primary_name, legacy_name = NULL) {
    if (!is.null(res_list[[primary_name]])) {
      return(res_list[[primary_name]])
    }
    if (!is.null(legacy_name) && !is.null(res_list[[legacy_name]])) {
      return(res_list[[legacy_name]])
    }
    NULL
  }

  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", NULL)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", NULL)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", "uniprot_acc")
  da_q_val_thresh <- if (!is.null(da_q_val_thresh)) {
    da_q_val_thresh
  } else if (!is.null(de_q_val_thresh)) {
    de_q_val_thresh
  } else {
    checkParamsObjectFunctionSimplify(theObject, "da_q_val_thresh", 0.05)
  }
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", "gene_names")
  fdr_column <- checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
  raw_p_value_column <- checkParamsObjectFunctionSimplify(theObject, "raw_p_value_column", "raw_pvalue")
  log2fc_column <- checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", "Entry")
  display_columns <- checkParamsObjectFunctionSimplify(theObject, "display_columns", c("best_uniprot_acc"))

  theObject <- updateParamInObject(theObject, "uniprot_tbl")
  theObject <- updateParamInObject(theObject, "publication_graphs_dir")
  theObject <- updateParamInObject(theObject, "args_row_id")
  theObject <- updateParamInObject(theObject, "da_q_val_thresh")
  theObject <- updateParamInObject(theObject, "gene_names_column")
  theObject <- updateParamInObject(theObject, "fdr_column")
  theObject <- updateParamInObject(theObject, "raw_p_value_column")
  theObject <- updateParamInObject(theObject, "log2fc_column")
  theObject <- updateParamInObject(theObject, "uniprot_id_column")
  theObject <- updateParamInObject(theObject, "display_columns")


  ## Write interactive volcano plot

  da_proteins_long <- safe_get_slot(results_list, "da_proteins_long", "de_proteins_long")
  contrasts_results <- safe_get_slot(results_list, "contrasts_results")

  if (is.null(da_proteins_long) || is.null(contrasts_results)) {
    stop("Results list does not contain the expected DA/DE volcano inputs.", call. = FALSE)
  }

  # Use helper to extract counts table instead of direct slot access
  result_object <- safe_get_slot(results_list, "theObject")
  if (is.null(result_object)) {
    result_object <- theObject
  }

  counts_mat <- getCountsTable(result_object) |>
    column_to_rownames(args_row_id) |>
    as.matrix()

  this_design_matrix <- result_object@design_matrix

  rownames(this_design_matrix) <- this_design_matrix[, result_object@sample_id]

  this_groups <- this_design_matrix[colnames(counts_mat), "group"]

  writeInteractiveVolcanoPlotProteomics(da_proteins_long,
    uniprot_tbl = uniprot_tbl,
    fit.eb = contrasts_results$fit.eb,
    publication_graphs_dir = publication_graphs_dir,
    args_row_id = args_row_id,
    fdr_column = fdr_column,
    raw_p_value_column = raw_p_value_column,
    log2fc_column = log2fc_column,
    da_q_val_thresh = da_q_val_thresh,
    counts_tbl = counts_mat,
    groups = this_groups,
    uniprot_id_column = uniprot_id_column,
    gene_names_column = gene_names_column,
    display_columns = display_columns
  )
}
