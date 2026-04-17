## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Normalise between Arrays
#' @export
#' @param theObject Object of class ProteinQuantitativeData
#' @param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
setMethod(
  f = "normaliseBetweenSamples",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, normalisation_method = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    normalisation_method <- checkParamsObjectFunctionSimplify(
      theObject,
      "normalisation_method",
      "cyclicloess"
    )

    theObject <- updateParamInObject(theObject, "normalisation_method")

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix[!is.finite(frozen_protein_matrix)] <- NA

    normalised_frozen_protein_matrix <- frozen_protein_matrix

    print(paste0("normalisation_method = ", normalisation_method))

    switch(normalisation_method,
      cyclicloess = {
        normalised_frozen_protein_matrix <- normalizeCyclicLoess(frozen_protein_matrix)
      },
      quantile = {
        normalised_frozen_protein_matrix <- normalizeQuantiles(frozen_protein_matrix)
      },
      scale = {
        normalised_frozen_protein_matrix <- normalizeMedianAbsValues(frozen_protein_matrix)
      },
      none = {
        normalised_frozen_protein_matrix <- frozen_protein_matrix
      }
    )

    normalised_frozen_protein_matrix[!is.finite(normalised_frozen_protein_matrix)] <- NA

    # normalised_frozen_protein_matrix_filt <- as.data.frame( normalised_frozen_protein_matrix ) |>
    #   dplyr::filter( if_all( everything(), \(x) { !is.na(x) } ) ) |>
    #   as.matrix()

    theObject@protein_quant_table <- normalised_frozen_protein_matrix |>
      as.data.frame() |>
      rownames_to_column(protein_id_column)

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

#' @export
setMethod(
  f = "ruvCancor",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, ctrl = NULL, num_components_to_impute = NULL, ruv_grouping_variable = NULL) {
    message("--- DEBUG66: Entering ruvCancor (S4) ---")

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    group_id <- theObject@group_id
    sample_id <- theObject@sample_id

    ctrl <- checkParamsObjectFunctionSimplify(theObject, "ctrl", NULL)
    num_components_to_impute <- checkParamsObjectFunctionSimplify(theObject, "num_components_to_impute", 2)
    ruv_grouping_variable <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", NULL)

    message(sprintf("   DEBUG66 S4 Param: num_components_to_impute = %d", num_components_to_impute))
    message(sprintf("   DEBUG66 S4 Param: ruv_grouping_variable = %s", ruv_grouping_variable))
    message(sprintf("   DEBUG66 S4 Param: ctrl length = %d (TRUE count: %d)", length(ctrl), sum(ctrl, na.rm = TRUE)))

    theObject <- updateParamInObject(theObject, "ctrl")
    theObject <- updateParamInObject(theObject, "num_components_to_impute")
    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

    if (!ruv_grouping_variable %in% colnames(design_matrix)) {
      stop(paste0(
        "The 'ruv_grouping_variable = ",
        ruv_grouping_variable,
        "' is not a column in the design matrix."
      ))
    }

    if (is.na(num_components_to_impute) || num_components_to_impute < 1) {
      stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
    }

    if (length(ctrl) < 5) {
      stop(paste0("The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
    }

    normalised_frozen_protein_matrix_filt <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    Y <- t(normalised_frozen_protein_matrix_filt[, design_matrix |> dplyr::pull(!!sym(sample_id))])
    if (length(which(is.na(normalised_frozen_protein_matrix_filt))) > 0) {
      message("   DEBUG66 S4: Performing imputation (NIPALS)...")
      Y <- impute.nipals(t(normalised_frozen_protein_matrix_filt[, design_matrix |> dplyr::pull(!!sym(sample_id))]),
        ncomp = num_components_to_impute
      )
    }

    message("   DEBUG66 S4: Calling ruv_cancorplot...")
    cancorplot_r2 <- ruv_cancorplot(Y,
      X = design_matrix |>
        dplyr::pull(!!sym(ruv_grouping_variable)),
      ctl = ctrl
    )
    message("   DEBUG66 S4: ruv_cancorplot returned.")
    cancorplot_r2
  }
)

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
