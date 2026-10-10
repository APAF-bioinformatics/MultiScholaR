## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
resolveProtLimmaNormalizer <- function(functionName, envir = parent.frame()) {
  normalizer <- get0(functionName, envir = envir, mode = "function", inherits = TRUE)
  if (is.function(normalizer)) {
    return(normalizer)
  }
  getExportedValue("limma", functionName)
}

orientProtRuvCorrectedMatrix <- function(correctedMatrix, featureIds, sampleIds) {
  if (!is.matrix(correctedMatrix)) {
    stop("RUV correction returned a non-matrix result.")
  }

  rowNames <- rownames(correctedMatrix)
  colNames <- colnames(correctedMatrix)

  if (!is.null(rowNames) && !is.null(colNames)) {
    if (all(featureIds %in% rowNames) && all(sampleIds %in% colNames)) {
      return(correctedMatrix[featureIds, sampleIds, drop = FALSE])
    }

    if (all(sampleIds %in% rowNames) && all(featureIds %in% colNames)) {
      return(t(correctedMatrix[sampleIds, featureIds, drop = FALSE]))
    }
  }

  if (nrow(correctedMatrix) == length(featureIds) && ncol(correctedMatrix) == length(sampleIds)) {
    rownames(correctedMatrix) <- featureIds
    colnames(correctedMatrix) <- sampleIds
    return(correctedMatrix)
  }

  if (nrow(correctedMatrix) == length(sampleIds) && ncol(correctedMatrix) == length(featureIds)) {
    correctedMatrix <- t(correctedMatrix)
    rownames(correctedMatrix) <- featureIds
    colnames(correctedMatrix) <- sampleIds
    return(correctedMatrix)
  }

  stop("RUV correction returned a matrix with incompatible dimensions.")
}

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
        normalised_frozen_protein_matrix <- resolveProtLimmaNormalizer("normalizeCyclicLoess")(frozen_protein_matrix)
      },
      quantile = {
        normalised_frozen_protein_matrix <- resolveProtLimmaNormalizer("normalizeQuantiles")(frozen_protein_matrix)
      },
      scale = {
        normalised_frozen_protein_matrix <- resolveProtLimmaNormalizer("normalizeMedianAbsValues")(frozen_protein_matrix)
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

    Y <- prepareRuvCancorInput(theObject, num_components_to_impute)
    ruvCancorFromInput(Y, theObject, ctrl, ruv_grouping_variable)
  }
)

#' The sample-by-protein input of the canonical-correlation plot: the normalised matrix, NIPALS-imputed
#' when it has missing values. It does not depend on the control set, so a percentage sweep prepares it
#' once (KIT-005, user decision 2026-10-10) and `ruvCancor` prepares it for a single call; both paths
#' compute the same matrix.
#' @noRd
prepareRuvCancorInput <- function(theObject, num_components_to_impute) {
  design_matrix <- theObject@design_matrix
  sample_columns <- design_matrix |> dplyr::pull(!!sym(theObject@sample_id))
  normalised_frozen_protein_matrix_filt <- theObject@protein_quant_table |>
    column_to_rownames(theObject@protein_id_column) |>
    as.matrix()
  Y <- t(normalised_frozen_protein_matrix_filt[, sample_columns])
  if (length(which(is.na(normalised_frozen_protein_matrix_filt))) > 0) {
    message("   DEBUG66 S4: Performing imputation (NIPALS)...")
    Y <- mixOmics::impute.nipals(t(normalised_frozen_protein_matrix_filt[, sample_columns]),
      ncomp = num_components_to_impute
    )
  }
  Y
}

#' The canonical-correlation plot for one control set on a prepared input (see `prepareRuvCancorInput`).
#' @noRd
ruvCancorFromInput <- function(Y, theObject, ctrl, ruv_grouping_variable) {
  message("   DEBUG66 S4: Calling ruv_cancorplot...")
  cancorplot_r2 <- ruv_cancorplot(Y,
    X = theObject@design_matrix |>
      dplyr::pull(!!sym(ruv_grouping_variable)),
    ctl = ctrl
  )
  message("   DEBUG66 S4: ruv_cancorplot returned.")
  cancorplot_r2
}

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
    k <- validateProtNormRuvK(k, "ruv_number_k")

    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
    theObject <- updateParamInObject(theObject, "ruv_number_k")
    theObject <- updateParamInObject(theObject, "ctrl")

    normalised_frozen_protein_matrix_filt <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    sample_ids <- design_matrix |>
      dplyr::pull(!!sym(sample_id)) |>
      as.character()

    feature_ids <- rownames(normalised_frozen_protein_matrix_filt)

    Y <- t(normalised_frozen_protein_matrix_filt[, sample_ids, drop = FALSE])
    # RUV-III cannot estimate more factors than it has control proteins observed in every sample; with
    # fewer it returns no finite values and the corrected table would be empty (KIT-005, 2026-10-10).
    complete_controls <- sum(ruvCompleteControlRows(theObject)[which(ctrl)], na.rm = TRUE)
    if (k > complete_controls) {
      cli::cli_abort(c("RUV-III cannot estimate {k} factor{?s} from {complete_controls} control protein{?s} observed in every sample",
        i = "choose a smaller k or a larger control set"), class = "prot_ruv_insufficient_controls")
    }

    # Replicate sets: technical replicates when the design has them, else the grouping variable.
    replicate_blocks <- resolveProtRuvReplicateBlocks(design_matrix, sample_id, group_id,
      replicate_group_column, ruv_grouping_variable)
    ruv_replicate_structure <- replicate_blocks$structure
    theObject <- updateParamInObject(theObject, "ruv_replicate_structure")
    M <- getRuvIIIReplicateMatrixHelper(
      replicate_blocks$design,
      !!sym(sample_id),
      .ruv_replicate_block
    )

    cln_mat <- RUVIII_C_Varying(
      k = k,
      Y = Y,
      M = M,
      toCorrect = colnames(Y),
      potentialControls = names(ctrl[which(ctrl)])
    )

    cln_mat_4 <- orientProtRuvCorrectedMatrix(
      correctedMatrix = cln_mat,
      featureIds = feature_ids,
      sampleIds = sample_ids
    )
    cln_mat_4 <- cln_mat_4[rowSums(is.finite(cln_mat_4)) > 0, , drop = FALSE]
    cln_mat_4 <- cln_mat_4[, colSums(is.finite(cln_mat_4)) > 0, drop = FALSE]
    if (nrow(cln_mat_4) == 0L || ncol(cln_mat_4) == 0L) {
      cli::cli_abort("RUV-III produced no finite corrected values (k = {k}, {sum(ctrl)} controls)",
        class = "prot_ruv_empty_result")
    }

    ruv_normalised_results_cln <- cln_mat_4 |>
      as.data.frame() |>
      rownames_to_column(protein_id_column)

    theObject@protein_quant_table <- ruv_normalised_results_cln

    theObject <- cleanDesignMatrix(theObject)

    return(theObject)
  }
)
