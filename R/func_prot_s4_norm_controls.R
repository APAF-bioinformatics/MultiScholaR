## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Get Low Coefficient of Variation Proteins
#'
#' Identifies proteins with the lowest coefficient of variation for use as
#' negative controls in RUV normalization.
#'
#' @param theObject A ProteinQuantitativeData object
#' @param percentage_as_neg_ctrl Percentage of proteins to use as negative controls (default: 10)
#' @param num_neg_ctrl Number of negative control proteins (overrides percentage if set)
#' @return A logical vector indicating which proteins are selected as negative controls
#' @export
setMethod(
  f = "getLowCoefficientOfVariationProteins",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL
  ) {
    percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify(theObject, "percentage_as_neg_ctrl", 10)
    num_neg_ctrl <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_neg_ctrl",
      round(nrow(theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0)
    )

    theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
    theObject <- updateParamInObject(theObject, "num_neg_ctrl")

    list_of_control_genes <- theObject@protein_quant_table |>
      column_to_rownames(theObject@protein_id_column) |>
      t() |>
      as.data.frame() |>
      summarise(across(everything(), ~ sd(.) / mean(.))) |>
      t() |>
      as.data.frame() |>
      dplyr::rename(coefficient_of_variation = "V1") |>
      tibble::rownames_to_column(theObject@protein_id_column) |>
      arrange(coefficient_of_variation) |>
      head(num_neg_ctrl)

    control_gene_index_helper <- theObject@protein_quant_table |>
      dplyr::select(theObject@protein_id_column) |>
      mutate(index = row_number()) |>
      left_join(list_of_control_genes, by = theObject@protein_id_column) |>
      mutate(is_selected = case_when(
        is.na(coefficient_of_variation) ~ FALSE,
        TRUE ~ TRUE
      )) |>
      arrange(index) |>
      dplyr::select(!!sym(theObject@protein_id_column), is_selected) |>
      column_to_rownames(theObject@protein_id_column) |>
      t()

    control_gene_index <- control_gene_index_helper[1, ]

    control_gene_index
  }
)

#' @export
setMethod(
  f = "getNegCtrlProtAnova",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    ruv_grouping_variable = NULL,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL,
    ruv_qval_cutoff = NULL,
    ruv_fdr_method = NULL,
    exclude_pool_samples = TRUE
  ) {
    message("--- DEBUG66: Entering getNegCtrlProtAnova (S4) ---")

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    group_id <- theObject@group_id
    sample_id <- theObject@sample_id

    normalised_frozen_protein_matrix_filt <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", "replicates")
    percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify(theObject, "percentage_as_neg_ctrl", 10)
    num_neg_ctrl <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_neg_ctrl",
      round(nrow(theObject@protein_quant_table) * percentage_as_neg_ctrl / 100, 0)
    )
    ruv_qval_cutoff <- checkParamsObjectFunctionSimplify(theObject, "ruv_qval_cutoff", 0.05)
    ruv_fdr_method <- checkParamsObjectFunctionSimplify(theObject, "ruv_fdr_method", "BH")
    exclude_pool_samples <- checkParamsObjectFunctionSimplify(theObject, "exclude_pool_samples", TRUE)

    message(sprintf("   DEBUG66 S4 Param: ruv_grouping_variable = %s", ruv_grouping_variable))
    message(sprintf("   DEBUG66 S4 Param: percentage_as_neg_ctrl = %s", percentage_as_neg_ctrl))
    message(sprintf("   DEBUG66 S4 Param: num_neg_ctrl = %s", num_neg_ctrl))
    message(sprintf("   DEBUG66 S4 Param: exclude_pool_samples = %s", exclude_pool_samples))

    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
    theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
    theObject <- updateParamInObject(theObject, "num_neg_ctrl")
    theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")
    theObject <- updateParamInObject(theObject, "ruv_fdr_method")
    theObject <- updateParamInObject(theObject, "exclude_pool_samples")

    # Prepare design matrix: keep only ruv_grouping_variable column, set sample IDs as rownames
    # This ensures we have the grouping information needed for Pool/QC detection
    design_matrix_for_anova <- design_matrix |>
      column_to_rownames(sample_id) |>
      dplyr::select(!!sym(ruv_grouping_variable))

    message("   DEBUG66 S4: Calling getNegCtrlProtAnovaHelper...")
    control_genes_index <- getNegCtrlProtAnovaHelper(normalised_frozen_protein_matrix_filt[, design_matrix |> dplyr::pull(!!sym(sample_id))],
      design_matrix = design_matrix_for_anova,
      grouping_variable = ruv_grouping_variable,
      percentage_as_neg_ctrl = percentage_as_neg_ctrl,
      num_neg_ctrl = num_neg_ctrl,
      ruv_qval_cutoff = ruv_qval_cutoff,
      ruv_fdr_method = ruv_fdr_method,
      exclude_pool_samples = exclude_pool_samples
    )

    message(sprintf("   DEBUG66 S4: Helper returned %d control genes", sum(control_genes_index)))
    return(control_genes_index)
  }
)
