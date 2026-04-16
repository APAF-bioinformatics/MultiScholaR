#' @export
setMethod(
  f = "removeProteinsWithOnlyOneReplicate",
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
  f = "plotRle",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)

    if (!is.null(sample_label)) {
      if (sample_label %in% colnames(design_matrix)) {
        rownames(design_matrix) <- design_matrix[, sample_label]
        colnames(frozen_protein_matrix) <- design_matrix[, sample_label]
      }
    } else {
      rownames(design_matrix) <- design_matrix[, sample_id]
    }

    # print( design_matrix)

    rowinfo_vector <- NA
    if (!is.na(grouping_variable)) {
      rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), grouping_variable]
    }

    print(rownames(design_matrix))
    print(colnames(frozen_protein_matrix))
    print(rowinfo_vector)
    # Handle missing/non-finite values
    working_matrix <- frozen_protein_matrix
    working_matrix[!is.finite(working_matrix)] <- NA

    rle_plot_before_cyclic_loess <- plotRleHelper(t(working_matrix),
      rowinfo = rowinfo_vector,
      yaxis_limit = yaxis_limit
    )

    return(rle_plot_before_cyclic_loess)
  }
)

#' @export
setMethod(
  f = "plotRleList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, list_of_columns, yaxis_limit = c()) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)
    rownames(design_matrix) <- design_matrix[, sample_id]

    # print( design_matrix)

    runOneRle <- function(column_name) {
      rowinfo_vector <- NA

      if (column_name %in% colnames(design_matrix)) {
        rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), column_name]
      }

      rle_plot_before_cyclic_loess <- plotRleHelper(t(frozen_protein_matrix),
        rowinfo = rowinfo_vector,
        yaxis_limit = yaxis_limit
      )

      return(rle_plot_before_cyclic_loess)
    }

    list_of_rle_plots <- purrr::map(list_of_columns, runOneRle)

    names(list_of_rle_plots) <- list_of_columns

    return(list_of_rle_plots)
  }
)

#' @export
savePlotRleList <- function(input_list, prefix = "RLE", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0("RLE", "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )


  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y){
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

#' @export
setMethod(
  f = "plotPca",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size = 8, cv_percentile = 0.90) {
    # Defensive checks
    if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
      stop("grouping_variable must be a single character string")
    }

    if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
      stop("shape_variable must be NULL or a single character string")
    }

    if (!grouping_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
    }

    if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
    }

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.na(label_column) || label_column == "") {
      label_column <- ""
    }

    required_cols <- c(sample_id, grouping_variable)
    if (!is.null(shape_variable)) {
      required_cols <- c(required_cols, shape_variable)
    }
    missing_cols <- setdiff(required_cols, colnames(design_matrix))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing columns in design matrix: %s", paste(missing_cols, collapse = ", ")))
    }

    tryCatch(
      {
        pca_plot <- plotPcaHelper(frozen_protein_matrix_pca,
          design_matrix,
          sample_id_column = sample_id,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable,
          label_column = label_column,
          title = title,
          geom.text.size = font_size
        )
        return(pca_plot)
      },
      error = function(e) {
        stop(sprintf("Error in plotPcaHelper: %s", e$message))
      }
    )
  }
)

#' @export
setMethod(
  f = "plotPcaList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, label_column, title, font_size = 8, cv_percentile = 0.90) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.na(label_column) || label_column == "") {
      label_column <- ""
    }

    pca_plots_list <- plotPcaListHelper(frozen_protein_matrix_pca,
      design_matrix,
      sample_id_column = sample_id,
      grouping_variables_list = grouping_variables_list,
      label_column = label_column,
      title = title,
      geom.text.size = font_size
    )

    return(pca_plots_list)
  }
)

#' @export
savePlotPcaList <- function(input_list, prefix = "PCA", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0("RLE", "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )


  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y){
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

#' @export
setMethod(
  f = "getPcaMatrix",
  signature = "ProteinQuantitativeData",
  definition = function(theObject) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id


    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA


    pca_mixomics_before_cyclic_loess <- mixOmics::pca(t(as.matrix(frozen_protein_matrix_pca)))$variates$X |>
      as.data.frame() |>
      rownames_to_column(sample_id) |>
      left_join(design_matrix, by = sample_id)


    return(pca_mixomics_before_cyclic_loess)
  }
)

#' @export
setMethod(
  f = "proteinTechRepCorrelation",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id
    tech_rep_column <- theObject@technical_replicate_id

    tech_rep_num_column <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_num_column", NULL)
    tech_rep_remove_regex <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", NULL)

    theObject <- updateParamInObject(theObject, "tech_rep_num_column")
    theObject <- updateParamInObject(theObject, "tech_rep_remove_regex")

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    protein_matrix_tech_rep <- proteinTechRepCorrelationHelper(design_matrix, frozen_protein_matrix_pca,
      protein_id_column = protein_id_column,
      sample_id_column = sample_id,
      tech_rep_column = tech_rep_column,
      tech_rep_num_column = tech_rep_num_column,
      tech_rep_remove_regex = tech_rep_remove_regex
    )

    return(protein_matrix_tech_rep)
  }
)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Plot Pearson Correlation
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex DEPRECATED - use exclude_pool_samples instead
#' @param correlation_group is the group where every pair of samples are compared
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" or "QC" in their name from correlation analysis.
#' @export
setMethod(
  f = "plotPearson",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    correlation_group_to_use <- correlation_group

    if (is.na(correlation_group)) {
      correlation_group_to_use <- theObject@technical_replicate_id
    }

    # Handle deprecated parameter (backward compatibility)
    if (!is.null(tech_rep_remove_regex)) {
      message("*** plotPearson: WARNING - tech_rep_remove_regex is deprecated, use exclude_pool_samples instead ***")
    }

    correlation_vec <- pearsonCorForSamplePairs(theObject,
      tech_rep_remove_regex = tech_rep_remove_regex,
      correlation_group = correlation_group_to_use,
      exclude_pool_samples = exclude_pool_samples
    )

    pearson_plot <- correlation_vec |>
      ggplot(aes(pearson_correlation)) +
      geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
      scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
      xlab("Pearson Correlation") +
      ylab("Counts") +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    return(pearson_plot)
  }
)

#' @title Pearson Correlation for Sample Pairs
#' @param theObject is an object of the type ProteinQuantitativeData
#' @param tech_rep_remove_regex DEPRECATED - use exclude_pool_samples instead
#' @param correlation_group is the group where every pair of samples are compared
#' @param exclude_pool_samples Logical. If TRUE (default), automatically exclude samples from groups containing "Pool" or "QC" in their name from correlation analysis. Pool/QC samples are excluded from within-group correlation calculations but remain in RUV-III correction.
#' @export
setMethod(
  f = "pearsonCorForSamplePairs",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering pearsonCorForSamplePairs                               |")
    message("+===========================================================================+")

    # Memory tracking - Entry
    entry_mem <- checkMemoryBoth("Entry", context = "pearsonCorForSamplePairs")

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    replicate_group_column <- theObject@technical_replicate_id
    if (!is.na(correlation_group)) {
      replicate_group_column <- correlation_group
    }

    # Handle deprecated parameter (backward compatibility)
    if (!is.null(tech_rep_remove_regex)) {
      message("*** PEARSON: WARNING - tech_rep_remove_regex is deprecated, use exclude_pool_samples instead ***")
    }

    exclude_pool_samples <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "exclude_pool_samples", TRUE)
    theObject <- updateParamInObject(theObject, "exclude_pool_samples")

    # --- OPTIMIZED MATRIX APPROACH ---
    message("--- DEBUG66 [pearsonCorForSamplePairs]: Using optimized matrix correlation ---")

    # 1. Extract numeric matrix
    # Get sample columns (exclude ID)
    sample_cols <- setdiff(colnames(protein_quant_table), protein_id_column)

    # Ensure they match design matrix samples
    # Note: protein_quant_table might contain pools/controls not in design if not cleaned, but usually they should align.
    valid_samples <- intersect(sample_cols, design_matrix[[sample_id]])

    if (length(valid_samples) == 0) {
      stop("No matching samples found between protein data and design matrix")
    }

    # Subset and convert to matrix
    mat <- as.matrix(protein_quant_table[, valid_samples, drop = FALSE])

    message(sprintf("   [pearsonCorForSamplePairs] Matrix dimensions: %d proteins x %d samples", nrow(mat), ncol(mat)))

    # 2. Calculate Correlation Matrix
    # fast vectorized correlation
    start_time <- Sys.time()
    cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
    end_time <- Sys.time()
    message(sprintf("   [pearsonCorForSamplePairs] Correlation matrix computed in %.2f seconds", as.numeric(difftime(end_time, start_time, units = "secs"))))

    # 3. Convert to Long Format (Pairwise) and Filter

    # Map samples to groups using design matrix
    # Create a named vector for fast lookup
    sample_to_group <- setNames(as.character(design_matrix[[replicate_group_column]]), as.character(design_matrix[[sample_id]]))

    # Get upper triangle indices to avoid duplicates and self-correlation
    upper_tri_indices <- which(upper.tri(cor_mat), arr.ind = TRUE)

    # Extract names and values
    samples1 <- rownames(cor_mat)[upper_tri_indices[, 1]]
    samples2 <- colnames(cor_mat)[upper_tri_indices[, 2]]
    cor_values <- cor_mat[upper_tri_indices]

    # Define column names expected by filterSamplesByProteinCorrelationThreshold
    col_x <- paste0(sample_id, ".x")
    col_y <- paste0(sample_id, ".y")

    # Create dataframe
    # Pre-allocate list for speed then convert
    result_list <- list()
    result_list[[col_x]] <- samples1
    result_list[[col_y]] <- samples2
    result_list[["pearson_correlation"]] <- cor_values

    result_df <- as.data.frame(result_list, stringsAsFactors = FALSE)

    # Add groups
    result_df$Group_1 <- sample_to_group[result_df[[col_x]]]
    result_df$Group_2 <- sample_to_group[result_df[[col_y]]]

    # Filter for within-group correlations only (Group_1 == Group_2)
    # Handle NA groups just in case
    result_df <- result_df[!is.na(result_df$Group_1) & !is.na(result_df$Group_2) & result_df$Group_1 == result_df$Group_2, ]

    # Add the group column required by the S4 object / downstream logic
    # The helper seems to expect the group info might be joined or passed, but let's see.
    # The return value is just the correlation table.
    # However, for `exclude_pool_samples`, we need the group.

    result_df[[replicate_group_column]] <- result_df$Group_1

    # Clean up temporary group columns
    result_df$Group_1 <- NULL
    result_df$Group_2 <- NULL

    # 4. Exclude Pool Samples Logic
    if (exclude_pool_samples) {
      # Detect Pool/QC groups (case-insensitive detection for "pool" or "qc")
      is_pool_qc_group <- grepl("pool|qc", result_df[[replicate_group_column]], ignore.case = TRUE)
      pool_qc_count <- sum(is_pool_qc_group)

      if (pool_qc_count > 0) {
        message(sprintf("*** PEARSON: Excluded %d sample pairs from Pool/QC groups ***", pool_qc_count))
        result_df <- result_df[!is_pool_qc_group, ]
      } else {
        message("*** PEARSON: No Pool/QC groups detected to exclude ***")
      }
    }

    message(sprintf("   [pearsonCorForSamplePairs] Final pair count: %d", nrow(result_df)))

    # Memory tracking - Exit
    reportMemoryDelta(entry_mem, "TOTAL pearsonCorForSamplePairs", context = "pearsonCorForSamplePairs")

    message("+===========================================================================+")
    message("|  DEBUG66: Exiting pearsonCorForSamplePairs                                |")
    message("+===========================================================================+")

    # Ensure the returned dataframe has the expected columns for the next step
    # filterSamplesByProteinCorrelationThreshold expects:
    # - {sample_id}.x
    # - {sample_id}.y
    # - pearson_correlation
    # And it might carry over replicate_group_column (it was in the original helper output)

    return(result_df)
  }
)

#' Create PC1/PC2 Boxplots from PCA ggplot Object
#'
#' @description Extracts PCA data from a ggplot object and creates boxplots
#' for PC1 and PC2 grouped by a specified variable. Works with both classic
#' S3 ggplot objects and new S7-based ggplot2 (v3.5+) objects.
#'
#' @param theObject A ggplot object containing PCA data with PC1 and PC2 columns
#' @param grouping_variable Character string specifying the grouping variable
#' @param title Plot title (default: "")
#' @param font_size Font size for plot text (default: 8)
#'
#' @return A patchwork combined plot with PC1 and PC2 boxplots
#' @export
setMethod(
  f = "plotPcaBox",
  signature = "ANY",
  definition = function(theObject, grouping_variable, title = "", font_size = 8, show_legend = FALSE) {
    # Validate input is a ggplot-like object (works with both S3 and S7 ggplot)
    if (!inherits(theObject, c("gg", "ggplot"))) {
      stop("theObject must be a ggplot object. Got class: ", paste(class(theObject), collapse = ", "))
    }

    # Extract data directly from the ggplot object
    if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
      pca_data <- as_tibble(theObject$data)
    } else {
      # Fall back to other extraction methods
      pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])

      # If the data doesn't have PC1/PC2, try to extract from the plot's environment
      if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
        # Try to get the data from the plot's environment
        if (exists("data", envir = environment(theObject$mapping$x))) {
          pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
        } else {
          stop("Could not extract PCA data from the ggplot object")
        }
      }
    }

    # Check if grouping variable exists in the data
    if (!grouping_variable %in% colnames(pca_data)) {
      stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
    }

    # Determine legend position
    legend_pos <- if (show_legend) "right" else "none"

    # Create PC1 boxplot
    pc1_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC1, fill = !!sym(grouping_variable))) +
      geom_boxplot(notch = TRUE) +
      theme_bw() +
      labs(
        title = title,
        x = "",
        y = "PC1"
      ) +
      theme(
        legend.position = legend_pos,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = font_size),
        plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    # Add explicit fill scale to support >6 discrete levels
    categorical_colors <- getCategoricalColourPalette()
    pc1_box <- pc1_box + scale_fill_manual(values = categorical_colors)

    # Create PC2 boxplot
    pc2_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC2, fill = !!sym(grouping_variable))) +
      geom_boxplot(notch = TRUE) +
      theme_bw() +
      labs(
        x = "",
        y = "PC2"
      ) +
      theme(
        legend.position = legend_pos,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = font_size),
        plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    # Add explicit fill scale to support >6 discrete levels
    pc2_box <- pc2_box + scale_fill_manual(values = categorical_colors)

    # Combine plots with minimal spacing
    # If legend is enabled, we might want to collect guides to avoid duplication
    # but for now let's keep it simple as patchwork/cowplot handles it
    combined_plot <- pc1_box / pc2_box +
      plot_layout(heights = c(1, 1), guides = if (show_legend) "collect" else NULL) +
      plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

    return(combined_plot)
  }
)

#' @export
setMethod(
  f = "plotDensityList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, title = "", font_size = 8) {
    # Create a list of density plots for each grouping variable
    density_plots_list <- purrr::map(grouping_variables_list, function(group_var) {
      tryCatch(
        {
          plotDensity(theObject,
            grouping_variable = group_var,
            title = title,
            font_size = font_size
          )
        },
        error = function(e) {
          warning(sprintf("Error creating density plot for %s: %s", group_var, e$message))
          return(NULL)
        }
      )
    })

    # Name the list elements with the grouping variables
    names(density_plots_list) <- grouping_variables_list

    # Remove any NULL elements (failed plots)
    density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

    return(density_plots_list)
  }
)

#' @export
savePlotDensityList <- function(input_list, prefix = "Density", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0(prefix, "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )

  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y) {
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

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
