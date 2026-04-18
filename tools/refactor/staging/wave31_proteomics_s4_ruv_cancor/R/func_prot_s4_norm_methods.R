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

