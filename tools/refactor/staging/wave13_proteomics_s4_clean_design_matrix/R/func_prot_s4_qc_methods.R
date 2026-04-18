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

