#' @export
setMethod(
  f = "removeRowsWithMissingValuesPercent",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    ruv_grouping_variable = NULL,
    groupwise_percentage_cutoff = NULL,
    max_groups_percentage_cutoff = NULL,
    proteins_intensity_cutoff_percentile = NULL
  ) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering removeRowsWithMissingValuesPercent S4 Method          |")
    message("+===========================================================================+")
    flush.console()

    # Memory tracking - Entry checkpoint
    entry_mem <- checkMemoryBoth("Entry", context = "removeRowsWithMissingValuesPercent")

    message("   DEBUG66 STEP 1: Extracting slots from S4 object...")
    flush.console()

    message("      DEBUG66: Extracting protein_quant_table...")
    flush.console()
    protein_quant_table <- theObject@protein_quant_table
    message(sprintf("      DEBUG66: protein_quant_table extracted. Class: %s", class(protein_quant_table)[1]))
    flush.console()

    message("      DEBUG66: Extracting protein_id_column...")
    flush.console()
    protein_id_column <- theObject@protein_id_column
    message(sprintf("      DEBUG66: protein_id_column = %s", protein_id_column))
    flush.console()

    message("      DEBUG66: Extracting design_matrix...")
    flush.console()
    design_matrix <- theObject@design_matrix
    message(sprintf("      DEBUG66: design_matrix extracted. Class: %s", class(design_matrix)[1]))
    flush.console()

    message("      DEBUG66: Extracting group_id...")
    flush.console()
    group_id <- theObject@group_id
    message(sprintf("      DEBUG66: group_id = %s", group_id))
    flush.console()

    message("      DEBUG66: Extracting sample_id...")
    flush.console()
    sample_id <- theObject@sample_id
    message(sprintf("      DEBUG66: sample_id = %s", sample_id))
    flush.console()

    message("      DEBUG66: Extracting technical_replicate_id...")
    flush.console()
    replicate_group_column <- theObject@technical_replicate_id
    message(sprintf("      DEBUG66: replicate_group_column = %s", ifelse(is.null(replicate_group_column), "NULL", replicate_group_column)))
    flush.console()

    message("   DEBUG66 STEP 2: Logging data dimensions...")
    flush.console()
    message(sprintf("      Data State (protein_quant_table): Dims = %d rows, %d cols", nrow(protein_quant_table), ncol(protein_quant_table)))
    flush.console()
    message(sprintf("      Data State (design_matrix): Dims = %d rows, %d cols", nrow(design_matrix), ncol(design_matrix)))
    flush.console()
    message("      DEBUG66: design_matrix columns:")
    message(paste("        ", paste(names(design_matrix), collapse = ", ")))
    flush.console()
    message("      DEBUG66: About to call head(design_matrix)...")
    flush.console()
    dm_head <- head(design_matrix)
    message("      DEBUG66: head() completed successfully")
    flush.console()
    message("      DEBUG66: Skipping print() - known to cause hangs in Shiny context")
    flush.console()
    # Use message instead of print to avoid Shiny reactive context issues
    message("      DEBUG66: First row sample_id: ", dm_head[[1]][1])
    flush.console()

    # Check for any issues with the protein_quant_table data
    message("      DEBUG66: Checking protein_quant_table for NaN/Inf values...")
    flush.console()
    pqt_matrix <- protein_quant_table[, -1] # Exclude protein ID column
    nan_count <- sum(is.nan(as.matrix(pqt_matrix)), na.rm = TRUE)
    inf_count <- sum(is.infinite(as.matrix(pqt_matrix)), na.rm = TRUE)
    na_count <- sum(is.na(as.matrix(pqt_matrix)))
    message(sprintf("      DEBUG66: protein_quant_table has %d NaN, %d Inf, %d NA values", nan_count, inf_count, na_count))
    flush.console()

    message("   DEBUG66 STEP 3: Resolving ruv_grouping_variable...")
    flush.console()

    # DEBUG66: Quick check of @args (simplified to avoid potential issues)
    message("      DEBUG66: Checking theObject@args...")
    flush.console()
    if (!is.null(theObject@args)) {
      args_names <- names(theObject@args)
      message(sprintf("      DEBUG66: theObject@args has %d entries", length(args_names)))
      flush.console()

      # Only check specific known function args, don't iterate all
      ruv_var_found <- NULL
      if ("ruvIII_C_Varying" %in% args_names &&
        !is.null(theObject@args$ruvIII_C_Varying$ruv_grouping_variable)) {
        ruv_var_found <- theObject@args$ruvIII_C_Varying$ruv_grouping_variable
        message(sprintf("      DEBUG66: Found ruv_grouping_variable in ruvIII_C_Varying: '%s'", ruv_var_found))
        flush.console()
      }
      if ("getNegCtrlProtAnova" %in% args_names &&
        !is.null(theObject@args$getNegCtrlProtAnova$ruv_grouping_variable)) {
        ruv_var_found <- theObject@args$getNegCtrlProtAnova$ruv_grouping_variable
        message(sprintf("      DEBUG66: Found ruv_grouping_variable in getNegCtrlProtAnova: '%s'", ruv_var_found))
        flush.console()
      }
    } else {
      message("      DEBUG66 WARNING: theObject@args is NULL!")
      flush.console()
    }

    message("      DEBUG66: About to call checkParamsObjectFunctionSimplify for ruv_grouping_variable...")
    flush.console()
    ruv_grouping_variable <- checkParamsObjectFunctionSimplify(
      theObject,
      "ruv_grouping_variable",
      NULL
    )
    message(sprintf("      DEBUG66: ruv_grouping_variable resolved = %s", ifelse(is.null(ruv_grouping_variable), "NULL *** POTENTIAL PROBLEM ***", ruv_grouping_variable)))
    flush.console()

    # CRITICAL CHECK: If ruv_grouping_variable is NULL, we need to get it from another source
    if (is.null(ruv_grouping_variable)) {
      message("      DEBUG66 WARNING: ruv_grouping_variable is NULL!")
      message("      DEBUG66: Attempting to find it from ruvIII_C_Varying or getNegCtrlProtAnova args...")
      flush.console()

      # Try to get from ruvIII_C_Varying first
      if (!is.null(theObject@args$ruvIII_C_Varying$ruv_grouping_variable)) {
        ruv_grouping_variable <- theObject@args$ruvIII_C_Varying$ruv_grouping_variable
        message(sprintf("      DEBUG66: Found in @args$ruvIII_C_Varying: '%s'", ruv_grouping_variable))
        flush.console()
      } else if (!is.null(theObject@args$getNegCtrlProtAnova$ruv_grouping_variable)) {
        ruv_grouping_variable <- theObject@args$getNegCtrlProtAnova$ruv_grouping_variable
        message(sprintf("      DEBUG66: Found in @args$getNegCtrlProtAnova: '%s'", ruv_grouping_variable))
        flush.console()
      } else {
        # Fallback to group_id from the S4 object
        ruv_grouping_variable <- theObject@group_id
        message(sprintf("      DEBUG66: FALLBACK to theObject@group_id: '%s'", ruv_grouping_variable))
        flush.console()
      }
    }

    # Validate that the column exists in design_matrix
    if (!ruv_grouping_variable %in% names(design_matrix)) {
      message(sprintf("      DEBUG66 CRITICAL ERROR: ruv_grouping_variable '%s' NOT FOUND in design_matrix!", ruv_grouping_variable))
      message(sprintf("      DEBUG66: design_matrix columns are: %s", paste(names(design_matrix), collapse = ", ")))
      flush.console()
      stop(sprintf(
        "ruv_grouping_variable '%s' not found in design_matrix columns: %s",
        ruv_grouping_variable, paste(names(design_matrix), collapse = ", ")
      ))
    }
    message(sprintf("      DEBUG66: Confirmed '%s' exists in design_matrix", ruv_grouping_variable))
    flush.console()

    message("   DEBUG66 STEP 4: Resolving groupwise_percentage_cutoff...")
    flush.console()
    groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify(
      theObject,
      "groupwise_percentage_cutoff",
      50
    )
    message(sprintf("      DEBUG66: groupwise_percentage_cutoff resolved = %g", groupwise_percentage_cutoff))
    flush.console()

    message("   DEBUG66 STEP 5: Resolving max_groups_percentage_cutoff...")
    flush.console()
    max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify(
      theObject,
      "max_groups_percentage_cutoff",
      50
    )
    message(sprintf("      DEBUG66: max_groups_percentage_cutoff resolved = %g", max_groups_percentage_cutoff))
    flush.console()

    message("   DEBUG66 STEP 6: Resolving proteins_intensity_cutoff_percentile...")
    flush.console()
    proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify(
      theObject,
      "proteins_intensity_cutoff_percentile",
      1
    )
    message(sprintf("      DEBUG66: proteins_intensity_cutoff_percentile resolved = %g", proteins_intensity_cutoff_percentile))
    flush.console()

    message("   DEBUG66 STEP 7: Updating parameters in S4 object...")
    flush.console()

    message("      DEBUG66: Updating ruv_grouping_variable...")
    flush.console()
    theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
    message("      DEBUG66: ruv_grouping_variable updated")
    flush.console()

    message("      DEBUG66: Updating groupwise_percentage_cutoff...")
    flush.console()
    theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
    message("      DEBUG66: groupwise_percentage_cutoff updated")
    flush.console()

    message("      DEBUG66: Updating max_groups_percentage_cutoff...")
    flush.console()
    theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
    message("      DEBUG66: max_groups_percentage_cutoff updated")
    flush.console()

    message("      DEBUG66: Updating proteins_intensity_cutoff_percentile...")
    flush.console()
    theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")
    message("      DEBUG66: proteins_intensity_cutoff_percentile updated")
    flush.console()

    message("   DEBUG66 STEP 8: Preparing to call helper function...")
    flush.console()
    message(sprintf("      Helper Args: cols = %s", protein_id_column))
    message(sprintf("      Helper Args: sample_id = %s", sample_id))
    message(sprintf("      Helper Args: row_id = %s", protein_id_column))
    message(sprintf("      Helper Args: grouping_variable = %s", ruv_grouping_variable))
    message(sprintf("      Helper Args: groupwise_percentage_cutoff = %g", groupwise_percentage_cutoff))
    message(sprintf("      Helper Args: max_groups_percentage_cutoff = %g", max_groups_percentage_cutoff))
    message(sprintf("      Helper Args: proteins_intensity_cutoff_percentile = %g", proteins_intensity_cutoff_percentile))
    flush.console()

    message("   DEBUG66 STEP 9: CALLING removeRowsWithMissingValuesPercentHelper NOW...")
    flush.console()

    # Memory tracking - Before helper
    pre_helper_mem <- checkMemoryBoth("Before helper", context = "removeRowsWithMissingValuesPercent")

    theObject@protein_quant_table <- removeRowsWithMissingValuesPercentHelper(protein_quant_table,
      cols = protein_id_column,
      design_matrix = design_matrix,
      sample_id = !!sym(sample_id),
      row_id = !!sym(protein_id_column),
      grouping_variable = !!sym(ruv_grouping_variable),
      groupwise_percentage_cutoff = groupwise_percentage_cutoff,
      max_groups_percentage_cutoff = max_groups_percentage_cutoff,
      proteins_intensity_cutoff_percentile = proteins_intensity_cutoff_percentile,
      temporary_abundance_column = "Log_Abundance"
    )

    message("   DEBUG66 STEP 10: Helper function returned successfully!")
    flush.console()
    message(sprintf(
      "      DEBUG66: New protein_quant_table dims = %d rows, %d cols",
      nrow(theObject@protein_quant_table), ncol(theObject@protein_quant_table)
    ))
    flush.console()

    # Memory tracking - After helper
    reportMemoryDelta(pre_helper_mem, "helper function", context = "removeRowsWithMissingValuesPercent")

    message("   DEBUG66 STEP 11: Cleaning design matrix...")
    flush.console()

    # Memory tracking - Before cleanDesignMatrix
    pre_clean_mem <- checkMemoryBoth("Before cleanDesignMatrix", context = "removeRowsWithMissingValuesPercent")

    theObject <- cleanDesignMatrix(theObject)

    # Memory tracking - After cleanDesignMatrix (THE SUSPECTED CULPRIT)
    reportMemoryDelta(pre_clean_mem, "cleanDesignMatrix", context = "removeRowsWithMissingValuesPercent")

    message("   DEBUG66: cleanDesignMatrix completed")
    flush.console()

    # Memory tracking - Exit checkpoint
    reportMemoryDelta(entry_mem, "TOTAL removeRowsWithMissingValuesPercent", context = "removeRowsWithMissingValuesPercent")

    message("+===========================================================================+")
    message("|  DEBUG66: Exiting removeRowsWithMissingValuesPercent S4 Method           |")
    message("+===========================================================================+")
    flush.console()
    return(theObject)
  }
)

