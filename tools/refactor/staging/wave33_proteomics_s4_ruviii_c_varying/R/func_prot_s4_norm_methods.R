## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
setMethod(
  f = "ruvIII_C_Varying",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    group_id <- theObject@group_id
    sample_id <- theObject@sample_id
    replicate_group_column <- theObject@technical_replicate_id


    ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", NULL)
    k <- checkParamsObjectFunctionSimplify(theObject, "ruv_number_k", NULL)
    ctrl <- checkParamsObjectFunctionSimplify(theObject, "ctrl", NULL)

    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
    theObject <- updateParamInObject(theObject, "ruv_number_k")
    theObject <- updateParamInObject(theObject, "ctrl")

    normalised_frozen_protein_matrix_filt <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    Y <- t(normalised_frozen_protein_matrix_filt[, design_matrix |> dplyr::pull(!!sym(sample_id))])

    M <- getRuvIIIReplicateMatrixHelper(
      design_matrix,
      !!sym(sample_id),
      !!sym(ruv_grouping_variable)
    )

    cln_mat <- RUVIII_C_Varying(
      k = ruv_number_k,
      Y = Y,
      M = M,
      toCorrect = colnames(Y),
      potentialControls = names(ctrl[which(ctrl)])
    )

    # Remove samples with no values
    cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat), ]

    # Remove proteins with no values
    cln_mat_3 <- t(cln_mat_2)
    cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3), ]

    ruv_normalised_results_cln <- cln_mat_4 |>
      as.data.frame() |>
      rownames_to_column(protein_id_column)

    theObject@protein_quant_table <- ruv_normalised_results_cln

    theObject <- cleanDesignMatrix(theObject)

    return(theObject)
  }
)

