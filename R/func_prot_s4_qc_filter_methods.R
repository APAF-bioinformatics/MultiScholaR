# Resolve S4 QC method parameters without relying on call-stack names. Several
# module apply-step tests inject S4 generics via function arguments, which makes
# stack-introspective helpers see the injected argument name instead of the S4
# method section in @args.
.resolveProteinQcMethodParam <- function(theObject,
                                         section_name,
                                         param_name,
                                         explicit_value = NULL,
                                         default_value = NULL,
                                         accept_null = FALSE,
                                         check_fn = checkParamsObjectFunctionSimplify) {
  if (!is.null(explicit_value)) {
    return(explicit_value)
  }

  section <- theObject@args[[section_name]]
  has_object_value <- is.list(section) && param_name %in% names(section)
  object_value <- if (has_object_value) {
    section[[param_name]]
  } else {
    NULL
  }
  if (has_object_value && (!is.null(object_value) || accept_null)) {
    return(object_value)
  }

  checked_value <- tryCatch(
    suppressWarnings(check_fn(theObject, param_name, default_value)),
    error = function(e) NULL
  )
  if (!is.null(checked_value) || accept_null) {
    return(checked_value)
  }

  default_value
}

.updateProteinQcMethodParam <- function(theObject,
                                        section_name,
                                        param_name,
                                        value) {
  if (is.null(theObject@args[[section_name]]) ||
      !is.list(theObject@args[[section_name]])) {
    theObject@args[[section_name]] <- list()
  }
  theObject@args[[section_name]][[param_name]] <- value
  theObject
}

#' @export
setMethod(
  f = "removeProteinsWithOnlyOneReplicate",
  signature = c(
    theObject = "ProteinQuantitativeData",
    core_utilisation = "ANY",
    grouping_variable = "ANY"
  ),
  definition = function(theObject, core_utilisation = NULL, grouping_variable = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    samples_id_tbl <- theObject@design_matrix
    sample_id_tbl_sample_id_column <- theObject@sample_id
    # replicate_group_column <- theObject@technical_replicate_id
    protein_id_column <- theObject@protein_id_column

    input_table_sample_id_column <- theObject@sample_id
    quantity_column <- "log_values"

    grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull(
      theObject,
      "grouping_variable",
      NULL
    )

    core_utilisation <- checkParamsObjectFunctionSimplify(
      theObject,
      "core_utilisation",
      NA
    )

    theObject <- updateParamInObject(theObject, "grouping_variable")
    theObject <- updateParamInObject(theObject, "core_utilisation")

    data_long_cln <- protein_quant_table |>
      pivot_longer(
        cols = !matches(protein_id_column),
        names_to = input_table_sample_id_column,
        values_to = quantity_column
      )

    protein_quant_table <- removeProteinsWithOnlyOneReplicateHelper(
      input_table = data_long_cln,
      samples_id_tbl = samples_id_tbl,
      input_table_sample_id_column = !!sym(input_table_sample_id_column),
      sample_id_tbl_sample_id_column = !!sym(sample_id_tbl_sample_id_column),
      replicate_group_column = !!sym(grouping_variable),
      protein_id_column = !!sym(protein_id_column),
      quantity_column = !!sym(quantity_column),
      core_utilisation = core_utilisation
    )


    theObject@protein_quant_table <- protein_quant_table |>
      pivot_wider(
        id_cols = !!sym(protein_id_column),
        names_from = !!sym(input_table_sample_id_column),
        values_from = !!sym(quantity_column)
      )

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

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

    message("      DEBUG66: About to resolve ruv_grouping_variable...")
    flush.console()
    ruv_grouping_variable <- .resolveProteinQcMethodParam(
      theObject = theObject,
      section_name = "removeRowsWithMissingValuesPercent",
      param_name = "ruv_grouping_variable",
      explicit_value = ruv_grouping_variable,
      default_value = NULL,
      accept_null = TRUE
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
    groupwise_percentage_cutoff <- .resolveProteinQcMethodParam(
      theObject = theObject,
      section_name = "removeRowsWithMissingValuesPercent",
      param_name = "groupwise_percentage_cutoff",
      explicit_value = groupwise_percentage_cutoff,
      default_value = 50
    )
    message(sprintf("      DEBUG66: groupwise_percentage_cutoff resolved = %g", groupwise_percentage_cutoff))
    flush.console()

    message("   DEBUG66 STEP 5: Resolving max_groups_percentage_cutoff...")
    flush.console()
    max_groups_percentage_cutoff <- .resolveProteinQcMethodParam(
      theObject = theObject,
      section_name = "removeRowsWithMissingValuesPercent",
      param_name = "max_groups_percentage_cutoff",
      explicit_value = max_groups_percentage_cutoff,
      default_value = 50
    )
    message(sprintf("      DEBUG66: max_groups_percentage_cutoff resolved = %g", max_groups_percentage_cutoff))
    flush.console()

    message("   DEBUG66 STEP 6: Resolving proteins_intensity_cutoff_percentile...")
    flush.console()
    proteins_intensity_cutoff_percentile <- .resolveProteinQcMethodParam(
      theObject = theObject,
      section_name = "removeRowsWithMissingValuesPercent",
      param_name = "proteins_intensity_cutoff_percentile",
      explicit_value = proteins_intensity_cutoff_percentile,
      default_value = 1
    )
    message(sprintf("      DEBUG66: proteins_intensity_cutoff_percentile resolved = %g", proteins_intensity_cutoff_percentile))
    flush.console()

    message("   DEBUG66 STEP 7: Updating parameters in S4 object...")
    flush.console()

    message("      DEBUG66: Updating ruv_grouping_variable...")
    flush.console()
    theObject <- .updateProteinQcMethodParam(
      theObject,
      "removeRowsWithMissingValuesPercent",
      "ruv_grouping_variable",
      ruv_grouping_variable
    )
    message("      DEBUG66: ruv_grouping_variable updated")
    flush.console()

    message("      DEBUG66: Updating groupwise_percentage_cutoff...")
    flush.console()
    theObject <- .updateProteinQcMethodParam(
      theObject,
      "removeRowsWithMissingValuesPercent",
      "groupwise_percentage_cutoff",
      groupwise_percentage_cutoff
    )
    message("      DEBUG66: groupwise_percentage_cutoff updated")
    flush.console()

    message("      DEBUG66: Updating max_groups_percentage_cutoff...")
    flush.console()
    theObject <- .updateProteinQcMethodParam(
      theObject,
      "removeRowsWithMissingValuesPercent",
      "max_groups_percentage_cutoff",
      max_groups_percentage_cutoff
    )
    message("      DEBUG66: max_groups_percentage_cutoff updated")
    flush.console()

    message("      DEBUG66: Updating proteins_intensity_cutoff_percentile...")
    flush.console()
    theObject <- .updateProteinQcMethodParam(
      theObject,
      "removeRowsWithMissingValuesPercent",
      "proteins_intensity_cutoff_percentile",
      proteins_intensity_cutoff_percentile
    )
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

#' @export
setMethod(
  f = "filterSamplesByProteinCorrelationThreshold",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering filterSamplesByProteinCorrelationThreshold             |")
    message("+===========================================================================+")

    # Memory tracking - Entry
    entry_mem <- checkMemoryBoth("Entry", context = "filterSamplesByProteinCorrelationThreshold")

    pearson_correlation_per_pair <- checkParamsObjectFunctionSimplify(theObject,
      "pearson_correlation_per_pair",
      default_value = NULL
    )
    min_pearson_correlation_threshold <- checkParamsObjectFunctionSimplify(theObject,
      "min_pearson_correlation_threshold",
      default_value = 0.75
    )

    theObject <- updateParamInObject(theObject, "pearson_correlation_per_pair")
    theObject <- updateParamInObject(theObject, "min_pearson_correlation_threshold")

    # Memory tracking - Before helper
    pre_helper_mem <- checkMemoryBoth("Before helper", context = "filterSamplesByProteinCorrelationThreshold")

    filtered_table <- filterSamplesByProteinCorrelationThresholdHelper(
      pearson_correlation_per_pair,
      protein_intensity_table = theObject@protein_quant_table,
      min_pearson_correlation_threshold = min_pearson_correlation_threshold,
      filename_column_x = !!sym(paste0(theObject@sample_id, ".x")),
      filename_column_y = !!sym(paste0(theObject@sample_id, ".y")),
      protein_id_column = theObject@protein_id_column,
      correlation_column = pearson_correlation
    )

    # Memory tracking - After helper
    reportMemoryDelta(pre_helper_mem, "helper function", context = "filterSamplesByProteinCorrelationThreshold")

    theObject@protein_quant_table <- filtered_table

    # Memory tracking - Before cleanDesignMatrix
    pre_clean_mem <- checkMemoryBoth("Before cleanDesignMatrix", context = "filterSamplesByProteinCorrelationThreshold")

    theObject <- cleanDesignMatrix(theObject)

    # Memory tracking - After cleanDesignMatrix
    reportMemoryDelta(pre_clean_mem, "cleanDesignMatrix", context = "filterSamplesByProteinCorrelationThreshold")

    # Memory tracking - Exit
    reportMemoryDelta(entry_mem, "TOTAL filterSamplesByProteinCorrelationThreshold", context = "filterSamplesByProteinCorrelationThreshold")

    message("+===========================================================================+")
    message("|  DEBUG66: Exiting filterSamplesByProteinCorrelationThreshold              |")
    message("+===========================================================================+")

    theObject
  }
)

#' @export
setMethod(
  f = "cleanDesignMatrix",
  signature = "ProteinQuantitativeData",
  definition = function(theObject) {
    # --- MEMORY OPTIMIZED: Using base R to avoid tidyverse environment capture ---

    # Memory tracking - Entry
    entry_mem <- checkMemoryBoth("Entry", context = "cleanDesignMatrix")

    # Get sample IDs from protein_quant_table columns (excluding the protein ID column)
    samples_id_vector <- setdiff(colnames(theObject@protein_quant_table), theObject@protein_id_column)

    # --- Validate Design Matrix --- #
    design_samples <- tryCatch(
      as.character(theObject@design_matrix[[theObject@sample_id]]),
      error = function(e) {
        character(0)
      }
    )
    if (length(design_samples) == 0) {
      warning(sprintf("cleanDesignMatrix: Could not extract valid sample IDs from design matrix column '%s'. Returning object unchanged.", theObject@sample_id), immediate. = TRUE)
      return(theObject)
    }

    # Find samples that exist in both protein_quant_table and design_matrix
    # This handles cases where pool samples might be in protein_quant_table but not in design_matrix
    sample_cols_in_design <- intersect(samples_id_vector, design_samples)
    if (length(sample_cols_in_design) == 0) {
      warning("cleanDesignMatrix: No sample columns identified in protein_quant_table matching design matrix sample IDs. Returning object unchanged.")
      return(theObject)
    }

    # Ensure columns are treated as character for matching consistency
    samples_id_vector_char <- as.character(sample_cols_in_design)

    # --- Filter and Reorder Design Matrix (Base R) --- #
    # Create a working copy of design matrix to avoid modifying original
    design_matrix_copy <- theObject@design_matrix

    # Ensure the sample ID column is character for matching
    design_matrix_copy[[theObject@sample_id]] <- as.character(design_matrix_copy[[theObject@sample_id]])

    # Find matching rows in the order of samples_id_vector_char
    matched_rows <- match(samples_id_vector_char, design_matrix_copy[[theObject@sample_id]])

    # Filter out NA matches (samples in data but not in design matrix)
    valid_matches <- matched_rows[!is.na(matched_rows)]

    if (length(valid_matches) == 0) {
      warning("cleanDesignMatrix: No matching samples found after filtering. Returning object unchanged.")
      return(theObject)
    }

    # Memory tracking - Before subsetting
    checkMemoryBoth("Before subset", context = "cleanDesignMatrix")

    # Subset and reorder design matrix using base R indexing
    cleaned_design_matrix <- design_matrix_copy[valid_matches, , drop = FALSE]

    # Reset row names to avoid confusion
    rownames(cleaned_design_matrix) <- NULL

    # Memory tracking - Before slot assignment (CRITICAL - This can trigger S4 copy)
    pre_assign_mem <- checkMemoryBoth("Before slot assignment", context = "cleanDesignMatrix")

    # Assign back to object
    theObject@design_matrix <- as.data.frame(cleaned_design_matrix)

    # Memory tracking - After slot assignment
    reportMemoryDelta(pre_assign_mem, "slot assignment (S4 copy trigger)", context = "cleanDesignMatrix")

    # Clean up intermediate objects to free memory
    rm(design_matrix_copy, matched_rows, valid_matches, cleaned_design_matrix)
    gc()

    # Memory tracking - Exit
    reportMemoryDelta(entry_mem, "TOTAL cleanDesignMatrix", context = "cleanDesignMatrix")

    return(theObject)
  }
)

#' @export
setMethod(
  f = "proteinIntensityFiltering",
  signature = "ProteinQuantitativeData",
  definition = function(
    theObject,
    proteins_intensity_cutoff_percentile = NULL,
    proteins_proportion_of_samples_below_cutoff = NULL,
    core_utilisation = NULL
  ) {
    protein_quant_table <- theObject@protein_quant_table

    proteins_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify(
      theObject,
      "proteins_intensity_cutoff_percentile",
      NULL
    )
    proteins_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify(
      theObject,
      "proteins_proportion_of_samples_below_cutoff",
      NULL
    )
    core_utilisation <- checkParamsObjectFunctionSimplify(
      theObject,
      "core_utilisation",
      NA
    )

    theObject <- updateParamInObject(theObject, "proteins_intensity_cutoff_percentile")
    theObject <- updateParamInObject(theObject, "proteins_proportion_of_samples_below_cutoff")
    theObject <- updateParamInObject(theObject, "core_utilisation")


    data_long_cln <- protein_quant_table |>
      pivot_longer(
        cols = !matches(theObject@protein_id_column),
        names_to = theObject@sample_id,
        values_to = "log_values"
      ) |>
      mutate(temp = "")

    min_peptide_intensity_threshold <- ceiling(quantile(data_long_cln$log_values, na.rm = TRUE, probs = c(proteins_intensity_cutoff_percentile)))[1]

    peptide_normalised_pif_cln <- peptideIntensityFilteringHelper(data_long_cln,
      min_peptide_intensity_threshold = min_peptide_intensity_threshold,
      proteins_proportion_of_samples_below_cutoff = proteins_proportion_of_samples_below_cutoff,
      protein_id_column = !!sym(theObject@protein_id_column),
      peptide_sequence_column = temp,
      peptide_quantity_column = log_values,
      core_utilisation = core_utilisation
    )


    theObject@protein_quant_table <- peptide_normalised_pif_cln |>
      dplyr::select(-temp) |>
      pivot_wider(id_cols = theObject@protein_id_column, names_from = !!sym(theObject@sample_id), values_from = log_values)

    theObject <- cleanDesignMatrix(theObject)

    updated_object <- theObject

    return(updated_object)
  }
)

#' @export
setMethod(
  f = "filterMinNumPeptidesPerProtein",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, ...) {
    # Extract specific parameters from ...
    args <- list(...)
    num_peptides_per_protein_thresh <- args$num_peptides_per_protein_thresh
    num_peptidoforms_per_protein_thresh <- args$num_peptidoforms_per_protein_thresh
    verbose <- args$verbose

    # --- Parameter validation and defaults ---
    num_peptides_per_protein_thresh <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_peptides_per_protein_thresh",
      1
    )

    num_peptidoforms_per_protein_thresh <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_peptidoforms_per_protein_thresh",
      2
    )

    verbose <- checkParamsObjectFunctionSimplify(theObject, "verbose", TRUE)

    # Update parameters in object
    theObject <- updateParamInObject(theObject, "num_peptides_per_protein_thresh")
    theObject <- updateParamInObject(theObject, "num_peptidoforms_per_protein_thresh")
    theObject <- updateParamInObject(theObject, "verbose")

    if (verbose) {
      log_info("Starting protein filtering based on peptide and peptidoform evidence...")
      log_info("Minimum unique peptides per protein: {num_peptides_per_protein_thresh}")
      log_info("Minimum total peptidoforms per protein: {num_peptidoforms_per_protein_thresh}")
    }

    # Get the peptide summary table (which is now guaranteed to be in sync)
    peptide_summary <- theObject@args$limpa_dpc_quant_results$peptide_counts_per_protein
    if (is.null(peptide_summary)) {
      stop("Could not find the peptide summary table. Please run chooseBestProteinAccession first.")
    }

    # --- Perform the filtering ---
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    proteins_before <- nrow(protein_quant_table)

    protein_ids_to_keep <- peptide_summary |>
      dplyr::filter(peptide_count >= num_peptides_per_protein_thresh & peptidoform_count >= num_peptidoforms_per_protein_thresh) |>
      dplyr::pull(!!sym(protein_id_column))

    filtered_protein_table <- protein_quant_table |>
      dplyr::filter(!!sym(protein_id_column) %in% protein_ids_to_keep)

    proteins_after <- nrow(filtered_protein_table)

    if (verbose) {
      log_info("Proteins before filtering: {proteins_before}")
      log_info("Proteins after filtering: {proteins_after}")
      log_info("Proteins removed: {proteins_before - proteins_after}")
      if (proteins_before > 0) {
        log_info("Retention rate: {round(100 * proteins_after / proteins_before, 1)}%")
      }
    }

    # Update the main data table
    theObject@protein_quant_table <- filtered_protein_table

    # Also filter the EList for consistency
    if (!is.null(theObject@args$limpa_dpc_quant_results$quantified_elist)) {
      original_elist <- theObject@args$limpa_dpc_quant_results$quantified_elist
      indices_to_keep <- which(original_elist$genes$protein.id %in% protein_ids_to_keep)
      theObject@args$limpa_dpc_quant_results$quantified_elist <- original_elist[indices_to_keep, ]
    }

    return(theObject)
  }
)

