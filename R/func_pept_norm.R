# ============================================================================
# func_pept_norm.R
# ============================================================================
# Purpose: Peptide normalization functions
# 
# This file contains functions for peptide-level normalization, including
# RUV-III-C normalization for peptides, negative control selection, and
# log2 transformation. Functions in this file are used by mod_prot_norm.R
# and related normalization modules.
#
# Functions to extract here:
# - ruvIII_C_Varying(): S4 method for RUV-III-C normalization (peptide)
# - findBestNegCtrlPercentagePeptides(): Find optimal negative control for peptides
# - getNegCtrlProtAnovaPeptides(): Get negative control peptides using ANOVA
# - log2TransformPeptideMatrix(): S4 method for log2 transformation
# - Additional peptide normalization helper functions
#
# Dependencies:
# - limma, RUVSeq (or custom RUV implementation)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: ruvIII_C_Varying() (peptide method)
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Applies RUV-III-C normalization to peptides with varying k
# setMethod(f = "ruvIII_C_Varying", signature = "PeptideQuantitativeData", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 2: findBestNegCtrlPercentagePeptides()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Finds optimal percentage of negative control peptides
# setMethod(f = "findBestNegCtrlPercentagePeptides", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: getNegCtrlProtAnovaPeptides()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets negative control peptides using ANOVA
# setMethod(f = "getNegCtrlProtAnovaPeptides", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 4: log2TransformPeptideMatrix()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Applies log2 transformation to peptide matrix
# setMethod(f = "log2TransformPeptideMatrix", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 5: log2Transformation()
# Current location: R/helper_functions.R
# Description: General log2 transformation function
# log2Transformation <- function(...) {
#   # Extract from R/helper_functions.R
# }


# ----------------------------------------------------------------------------
# log2Transformation
# ----------------------------------------------------------------------------
#' @export
#' @title Log2 Transformation with Pseudo-count
#' @description Log 2 transformation with pseudo count
log2Transformation <- function(input_matrix) {

  pseudo_count <- min( input_matrix[input_matrix> 0] , na.rm=TRUE)/100
  input_matrix[input_matrix> 0 & !is.na(input_matrix)] <- input_matrix[input_matrix> 0 & !is.na(input_matrix)] + pseudo_count
  input_matrix <- log2(input_matrix)

  return(input_matrix )
}


# ----------------------------------------------------------------------------
# findBestNegCtrlPercentagePeptides
# ----------------------------------------------------------------------------
#' @export
setMethod(f="findBestNegCtrlPercentagePeptides"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject,
                                percentage_range = seq(1, 20, by = 1),
                                num_components_to_impute = 5,
                                ruv_grouping_variable = "group",
                                ruv_qval_cutoff = 0.05,
                                ruv_fdr_method = "qvalue",
                                separation_metric = "max_difference",
                                k_penalty_weight = 0.5,
                                max_acceptable_k = 3,
                                adaptive_k_penalty = TRUE,
                                verbose = TRUE,
                                ensure_matrix = TRUE) {
  
  # Input validation
  if (!inherits(theObject, "PeptideQuantitativeData")) {
    stop("theObject must be a PeptideQuantitativeData object")
  }
  
  # Ensure peptide matrix is calculated if requested
  if (ensure_matrix && (!"peptide_matrix" %in% slotNames(theObject) || is.null(theObject@peptide_matrix) || length(theObject@peptide_matrix) == 0)) {
    if (verbose) {
      log_info("Peptide matrix not found. Calculating peptide matrix...")
    }
    theObject <- calcPeptideMatrix(theObject)
  }
  
  if (length(percentage_range) == 0 || any(percentage_range <= 0) || any(percentage_range > 100)) {
    stop("percentage_range must contain values between 0 and 100")
  }
  
  if (!separation_metric %in% c("max_difference", "mean_difference", "auc", "weighted_difference")) {
    stop("separation_metric must be one of: 'max_difference', 'mean_difference', 'auc', 'weighted_difference'")
  }
  
  if (k_penalty_weight < 0 || k_penalty_weight > 1) {
    stop("k_penalty_weight must be between 0 and 1")
  }
  
  if (max_acceptable_k < 1 || !is.numeric(max_acceptable_k)) {
    stop("max_acceptable_k must be a positive number >= 1")
  }
  
  if (adaptive_k_penalty && max_acceptable_k < 3) {
    stop("max_acceptable_k must be at least 3 when adaptive_k_penalty is TRUE")
  }
  
  # Calculate adaptive max_acceptable_k if requested
  if (adaptive_k_penalty) {
    # Get sample size from the peptide matrix
    sample_size <- ncol(theObject@peptide_matrix)
    
    # Calculate adaptive max_acceptable_k based on sample size
    adaptive_max_k <- .peptide_calculateAdaptiveMaxK(sample_size)
    
    if (verbose) {
      log_info("Adaptive penalty enabled: Sample size = {sample_size}, Adaptive max_acceptable_k = {adaptive_max_k} (original = {max_acceptable_k})")
    }
    
    max_acceptable_k <- adaptive_max_k
  }
  
  # Detect small datasets and warn/adjust percentage range if needed
  sample_size <- ncol(theObject@peptide_matrix)
  if (sample_size < 15 && max(percentage_range) < 30) {
    if (verbose) {
      log_warn("Small dataset detected (n={sample_size}). Consider testing higher percentages (up to 30-50%) for better negative control identification.")
      log_warn("Current range: {paste(range(percentage_range), collapse = '-')}%. May need wider range for optimal results.")
    }
  }
  
  if (verbose) {
    log_info("Starting optimization of negative control percentage for PEPTIDE data with k value consideration...")
    log_info("Testing {length(percentage_range)} different percentages: {paste(range(percentage_range), collapse = '-')}%")
    log_info("K penalty weight: {k_penalty_weight}, Max acceptable k: {max_acceptable_k}")
    if (adaptive_k_penalty) {
      log_info("Using adaptive k penalty based on sample size")
    }
  }
  
  # Process all percentages using functional programming
  if (verbose) {
    log_info("Processing {length(percentage_range)} percentages using vectorized operations...")
  }
  
  # Create a function to process a single percentage
  process_percentage <- function(current_percentage, index) {
    if (verbose && index %% 5 == 0) {
      log_info("Testing percentage {index}/{length(percentage_range)}: {current_percentage}%")
    }
    
    tryCatch({
      # Get negative control peptides for current percentage
      control_genes_index <- getNegCtrlProtAnovaPeptides(
        theObject,
        ruv_grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = current_percentage,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
      
      # Check if we have enough control peptides
      num_controls <- sum(control_genes_index, na.rm = TRUE)
      if (num_controls < 5) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Only {num_controls} control peptides found (minimum 5 required). Skipping.")
        }
        return(list(
          percentage = current_percentage,
          separation_score = NA_real_,
          best_k = NA_real_,
          composite_score = NA_real_,
          num_controls = num_controls,
          valid_plot = FALSE,
          control_genes_index = NULL,
          cancor_plot = NULL
        ))
      }
      
      # Generate canonical correlation plot using FAST version (skips expensive DPC imputation)
      cancorplot <- ruvCancorFast(
        theObject,
        ctrl = control_genes_index,
        num_components_to_impute = num_components_to_impute,
        ruv_grouping_variable = ruv_grouping_variable,
        simple_imputation_method = "mean"  # Use simple mean imputation for speed
      )
      
      # Calculate separation score
      separation_score <- .peptide_calculateSeparationScore(cancorplot, separation_metric)
      
      # Calculate the best k using the existing findBestK function
      best_k <- tryCatch({
        findBestK(cancorplot)
      }, error = function(e) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Error calculating best k: {e$message}")
        }
        return(NA_real_)
      })
      
      # Calculate composite score that considers both separation and k value
      composite_score <- .peptide_calculateCompositeScore(
        separation_score, 
        best_k, 
        k_penalty_weight, 
        max_acceptable_k
      )
      
      return(list(
        percentage = current_percentage,
        separation_score = separation_score,
        best_k = best_k,
        composite_score = composite_score,
        num_controls = num_controls,
        valid_plot = TRUE,
        control_genes_index = control_genes_index,
        cancor_plot = cancorplot
      ))
      
    }, error = function(e) {
      if (verbose) {
        log_warn("Percentage {current_percentage}%: Error occurred - {e$message}")
      }
      return(list(
        percentage = current_percentage,
        separation_score = NA_real_,
        best_k = NA_real_,
        composite_score = NA_real_,
        num_controls = NA_integer_,
        valid_plot = FALSE,
        control_genes_index = NULL,
        cancor_plot = NULL
      ))
    })
  }
  
  # Use purrr::imap() for functional processing
  all_results <- percentage_range |>
    purrr::imap(process_percentage)
  
  # Extract results into proper data frame
  results <- all_results |>
    purrr::map_dfr(~ data.frame(
      percentage = .x$percentage,
      separation_score = .x$separation_score,
      best_k = .x$best_k,
      composite_score = .x$composite_score,
      num_controls = .x$num_controls,
      valid_plot = .x$valid_plot
    ))
  
  # Find the best result using composite score (considers both separation and k value)
  valid_results <- all_results[!is.na(purrr::map_dbl(all_results, "composite_score"))]
  
  if (length(valid_results) == 0) {
    stop("No valid percentage found. Please check your data and parameters.")
  }
  
  best_index <- which.max(purrr::map_dbl(valid_results, "composite_score"))
  best_result <- valid_results[[best_index]]
  
  best_percentage <- best_result$percentage
  best_control_genes_index <- best_result$control_genes_index
  best_cancor_plot <- best_result$cancor_plot
  best_separation_score <- best_result$separation_score
  best_composite_score <- best_result$composite_score
  best_k <- best_result$best_k
  
  # Final validation and logging
  if (verbose) {
    log_info("Optimization complete!")
    log_info("Best percentage: {best_percentage}% (composite score: {round(best_composite_score, 4)})")
    log_info("  - Separation score: {round(best_separation_score, 4)}")
    log_info("  - Best k value: {best_k}")
    log_info("  - Number of control peptides: {sum(best_control_genes_index, na.rm = TRUE)}")
  }
  
  # Return comprehensive results
  return(list(
    best_percentage = best_percentage,
    best_k = best_k,
    best_control_genes_index = best_control_genes_index,
    best_separation_score = best_separation_score,
    best_composite_score = best_composite_score,
    optimization_results = results,
    best_cancor_plot = best_cancor_plot,
    separation_metric_used = separation_metric,
    k_penalty_weight = k_penalty_weight,
    max_acceptable_k = max_acceptable_k,
    adaptive_k_penalty_used = adaptive_k_penalty,
    sample_size = if(adaptive_k_penalty) ncol(theObject@peptide_matrix) else NA
  ))
})


# ----------------------------------------------------------------------------
# getNegCtrlProtAnovaPeptides
# ----------------------------------------------------------------------------
#'@export
setMethod(f="getNegCtrlProtAnovaPeptides"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject
                                 , ruv_grouping_variable = NULL
                                 , percentage_as_neg_ctrl = NULL
                                 , num_neg_ctrl = NULL
                                 , ruv_qval_cutoff = NULL
                                 , ruv_fdr_method = NULL ) {
            
            peptide_data <- theObject@peptide_data
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id
            design_matrix <- theObject@design_matrix
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column
            group_id <- theObject@group_id
            
            normalised_frozen_peptide_matrix_filt <- theObject@peptide_matrix
            
            ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", "replicates")
            percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
            num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                               , "num_neg_ctrl"
                                                               , round(nrow( normalised_frozen_peptide_matrix_filt) * percentage_as_neg_ctrl / 100, 0))
            ruv_qval_cutoff <- checkParamsObjectFunctionSimplify( theObject, "ruv_qval_cutoff", 0.05)
            ruv_fdr_method <- checkParamsObjectFunctionSimplify( theObject, "ruv_fdr_method", "BH")
            
            theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
            theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
            theObject <- updateParamInObject(theObject, "num_neg_ctrl")
            theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")
            theObject <- updateParamInObject(theObject, "ruv_fdr_method")
            
            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id_column)) ]
                                                              , design_matrix = design_matrix |>
                                                                column_to_rownames(sample_id_column) |>
                                                                dplyr::select( -!!sym(group_id))
                                                              , grouping_variable = ruv_grouping_variable
                                                              , percentage_as_neg_ctrl = percentage_as_neg_ctrl
                                                              , num_neg_ctrl = num_neg_ctrl
                                                              , ruv_qval_cutoff = ruv_qval_cutoff
                                                              , ruv_fdr_method = ruv_fdr_method )
            
            return(control_genes_index)
          })


# ----------------------------------------------------------------------------
# log2TransformPeptideMatrix
# ----------------------------------------------------------------------------
#' Log2 Transform Peptide Matrix
#'
#' Transforms raw peptide intensity values to log2 scale for downstream normalization and RUV analysis.
#' This should be called after calcPeptideMatrix() and before normaliseBetweenSamples().
#' @export
#'
#' @param theObject A PeptideQuantitativeData object
#' @return PeptideQuantitativeData object with log2 transformed data
#' @export
setMethod(f="log2TransformPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject) {
            
            if (theObject@is_logged_data) {
              warning("Data appears to already be log-transformed (is_logged_data = TRUE). Skipping transformation.")
              return(theObject)
            }
            
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            current_quant_column <- theObject@norm_quantity_column
            
            # Log2 transform the matrix
            log2_peptide_matrix <- peptide_matrix
            
            # Handle zeros and negative values
            log2_peptide_matrix[log2_peptide_matrix <= 0] <- NA
            log2_peptide_matrix <- log2(log2_peptide_matrix)
            
            # Update matrix
            theObject@peptide_matrix <- log2_peptide_matrix
            
            # Also update peptide_data to maintain consistency
            # Create long format from log2 matrix
            log2_long <- log2_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_row_id") |>
              pivot_longer(cols = -peptide_row_id, 
                          names_to = theObject@sample_id, 
                          values_to = "log2_value")

            # Update peptide_data: match by protein%peptide ID and sample
            updated_peptide_data <- peptide_data |>
              mutate(peptide_row_id = paste(!!sym(theObject@protein_id_column), 
                                           !!sym(theObject@peptide_sequence_column), 
                                           sep = "%")) |>
              left_join(log2_long, by = c("peptide_row_id", theObject@sample_id)) |>
              mutate(!!sym(current_quant_column) := ifelse(!is.na(log2_value), 
                                                           log2_value, 
                                                           !!sym(current_quant_column))) |>
              select(-peptide_row_id, -log2_value)

            theObject@peptide_data <- updated_peptide_data
            
            # Mark as logged
            theObject@is_logged_data <- TRUE
            
            theObject <- cleanDesignMatrixPeptide(theObject)
            
            message("Peptide data successfully log2 transformed. Raw intensities converted to log2 scale.")
            message(paste("is_logged_data flag set to:", theObject@is_logged_data))
            
            return(theObject)
          })

