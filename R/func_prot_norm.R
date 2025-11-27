# ============================================================================
# func_prot_norm.R
# ============================================================================
# Purpose: Protein normalization functions
# 
# This file contains functions for protein-level normalization, including
# RUV-III-C normalization, negative control selection, and normalization
# parameter optimization. Functions in this file are used by mod_prot_norm.R
# and related normalization modules.
#
# Functions to extract here:
# - normaliseBetweenSamples(): S4 method for between-sample normalization
# - normaliseUntransformedData(): S4 method for untransformed data normalization
# - ruvIII_C_Varying(): S4 method for RUV-III-C normalization (protein)
# - ruvCancor(): S4 method for RUV canonical correlation analysis
# - ruvCancorFast(): S4 method for fast RUV canonical correlation
# - findBestNegCtrlPercentage(): Find optimal negative control percentage
# - findBestK(): Find optimal k value for RUV
# - getNegCtrlProtAnova(): Get negative control proteins using ANOVA
# - Additional normalization helper functions
#
# Dependencies:
# - limma, RUVSeq (or custom RUV implementation)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: normaliseBetweenSamples()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Normalizes data between samples using various methods
# setMethod(f = "normaliseBetweenSamples", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 2: normaliseUntransformedData()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Normalizes untransformed data
# setMethod(f = "normaliseUntransformedData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 3: ruvIII_C_Varying() (protein method)
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Applies RUV-III-C normalization with varying k
# setMethod(f = "ruvIII_C_Varying", signature = "ProteinQuantitativeData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 4: ruvCancor()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Performs RUV canonical correlation analysis
# setMethod(f = "ruvCancor", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 5: ruvCancorFast()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Fast version of RUV canonical correlation analysis
# setMethod(f = "ruvCancorFast", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 6: findBestNegCtrlPercentage()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds optimal percentage of negative control proteins
# findBestNegCtrlPercentage <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 7: findBestK()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds optimal k value for RUV normalization
# findBestK <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 8: findBestKForAssayList()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds best k for a list of assays
# findBestKForAssayList <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 9: getNegCtrlProtAnova()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets negative control proteins using ANOVA
# setMethod(f = "getNegCtrlProtAnova", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 10: getNegCtrlProtAnovaHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for getting negative control proteins
# getNegCtrlProtAnovaHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 11: getRuvIIIReplicateMatrixHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for creating RUV-III replicate matrix
# getRuvIIIReplicateMatrixHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 12: extractRuvResults()
# Current location: R/helper_functions.R
# Description: Extracts RUV normalization results
# extractRuvResults <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 13: updateRuvParameters()
# Current location: R/helper_functions.R
# Description: Updates RUV parameters in config
# updateRuvParameters <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 14: scaleCenterAndFillMissing()
# Current location: R/helper_functions.R
# Description: Scales, centers, and fills missing values
# scaleCenterAndFillMissing <- function(...) {
#   # Extract from R/helper_functions.R
# }


# ----------------------------------------------------------------------------
# updateRuvParameters
# ----------------------------------------------------------------------------
##################################################################################################################
#' @export
updateRuvParameters <- function(config_list, best_k, control_genes_index, percentage_as_neg_ctrl) {
  config_list$ruvParameters$best_k <- best_k
  config_list$ruvParameters$num_neg_ctrl <- length(control_genes_index)
  config_list$ruvParameters$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
  
  # Print the number of negative controls (as in the original code)
  config_list$ruvParameters$num_neg_ctrl
  
  # Return the updated config list
  return(config_list)
}


# ----------------------------------------------------------------------------
# getRuvIIIReplicateMatrixHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'Converts a design matrix to a biological replicate matrix for use with ruvIII.
#'@param design_matrix The design matrix with the sample ID in one column and the experimental group in another column
#'@param sample_id_column The name of the column with the sample ID, tidyverse style input.
#'@param grouping_variable The name of the column with the experimental group, tidyverse style input.
#'@param temp_column The name of the temporary column that indicates which samples are biological replicates of the same experimental group.
#'@return A numeric matrix with rows as samples, columns as experimental group, and a value of 1 for samples within the same experimental group represented by the same column, and a value of zero otherwise.
#'@export
getRuvIIIReplicateMatrixHelper <- function(design_matrix, sample_id_column, grouping_variable, temp_column = is_replicate_temp) {

  ruvIII_replicates_matrix <- design_matrix |>
    dplyr::select({ { sample_id_column } }, { { grouping_variable } }) |>
    mutate({ { temp_column } } := 1) |>
    pivot_wider(id_cols = as_string( as_name(enquo(sample_id_column ))),
                names_from = { { grouping_variable } },
                values_from = { { temp_column } },
                values_fill = 0) |>
    column_to_rownames(as_string(as_name(enquo(sample_id_column)))) |>
    as.matrix()

  ruvIII_replicates_matrix
}


# ----------------------------------------------------------------------------
# getNegCtrlProtAnovaHelper
# ----------------------------------------------------------------------------
#' Identify negative control proteins for use in removal of unwanted variation, using an ANOVA test.
#' @param data_matrix A matrix containing the log (base 2) protein abundance values where each column represents a sample and each row represents a protein group, and proteins as rows. The row ID are the protein accessions. The data is preferably median-scaled with missing values imputed.
#' @param design_matrix A data frame with the design matrix. Matches sample IDs to group IDs.
#' @param grouping_variable The name of the column with the experimental group, as a string.
#' @param num_neg_ctrl The number of negative control genes to select. Typically the number of genes with the highest q-value (e.g. least statistically significant). Default is 100
#' @param ruv_qval_cutoff The FDR threshold. No proteins with q-values lower than this value are included in the list of negative control proteins. This means the number of negative control proteins could be less than the number specified in \code{num_neg_ctrl} when some were excluded by this threshold.
#' @param ruv_fdr_method The FDR calculation method, default is "qvalue". The other option is "BH"
#' @return A boolean vector which indicates which row in the input data matrix is a control gene. The row is included if the value is TRUE. The names of each element is the row ID / protein accessions of the input data matrix.
#'@export
getNegCtrlProtAnovaHelper <- function(data_matrix
                                , design_matrix
                                , grouping_variable = "group"
                                , percentage_as_neg_ctrl = 10
                                , num_neg_ctrl = round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
                                , ruv_qval_cutoff = 0.05
                                , ruv_fdr_method = "qvalue") {

  ## Both percentage_as_neg_ctrl and num_neg_ctrl is missing, and number of proteins >= 50 use only 10 percent of the proteins as negative control by default
  if((is.null(percentage_as_neg_ctrl) ||
     is.na(percentage_as_neg_ctrl) ) &&
     (is.null(num_neg_ctrl) ||
      is.na(num_neg_ctrl) ) &&
     nrow(data_matrix) >= 50 ) {
    num_neg_ctrl <- round( nrow(data_matrix)*10/100, 0)
    warnings( paste0( getFunctionName(), ": Using 10% of proteins from the input matrix as negative controls by default.\n"))
  } else if (!is.null(percentage_as_neg_ctrl) &
             !is.na(percentage_as_neg_ctrl)) {
    num_neg_ctrl <- round( nrow(data_matrix)*percentage_as_neg_ctrl/100, 0)
  } else if(!is.null(num_neg_ctrl) &
       !is.na(num_neg_ctrl)) {
    num_neg_ctrl <- as.integer(num_neg_ctrl)
  } else {
    stop(paste0( getFunctionName(), ": Please provide either percentage_as_neg_ctrl or num_neg_ctrl.\n"))
  }

  ## Inspired by matANOVA function from PhosR package: http://www.bioconductor.org/packages/release/bioc/html/PhosR.html

  grps <- design_matrix[colnames(data_matrix), grouping_variable]

  ps <- apply(data_matrix, 1, function(x) {
       if( length( unique( grps[!is.na(x)] )  ) > 1 ) {
         summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
       } else {
          return(NA_real_)
       }
    })

  ps[is.na(ps)] <- 1

  aov <- c()

  if ( ruv_fdr_method == "qvalue") {
    aov <- qvalue(unlist(ps))$qvalues
  } else if ( ruv_fdr_method == "BH") {
    aov <- qvalue(unlist(ps), pi0=1)$qvalues
  } else {
    error( paste( "Input FDR method", ruv_fdr_method, "not valid") )
  }

  filtered_list <- aov[aov > ruv_qval_cutoff]

  list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)

  control_genes <- names(sort(filtered_list, decreasing = TRUE)[1:list_size])

  #nrow(data_matrix) - length(control_genes)
  control_genes_index <- rownames(data_matrix) %in% control_genes
  names(control_genes_index) <- rownames(data_matrix)

  return(control_genes_index)

}


# ----------------------------------------------------------------------------
# findBestK
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
findBestK <- function( cancorplot_r1) {
  controls_idx <- which(cancorplot_r1$data$featureset == "Control")
  all_idx <- which( cancorplot_r1$data$featureset == "All")
  difference_between_all_ctrl <- cancorplot_r1$data$cc[all_idx] - cancorplot_r1$data$cc[controls_idx]
  max_difference <- max(difference_between_all_ctrl, na.rm=TRUE)
  best_idx <- which( difference_between_all_ctrl == max_difference)
  best_k <- (cancorplot_r1$data$K[controls_idx] )[best_idx]
  return( best_k)
}


# ----------------------------------------------------------------------------
# findBestKForAssayList
# ----------------------------------------------------------------------------
#' Find the Best K Value for RUV from a List of Canonical Correlation Plots
#'
#' This function iterates over a list of canonical correlation plots (typically
#' generated by `ruvCancor` for multi-assay objects like `MetaboliteAssayData`)
#' and applies the `findBestK` logic to each plot to determine the optimal
#' number of unwanted variation factors (k) for each assay.
#'
#' @param cancor_plots_list A named list where each element is a ggplot object
#'   returned by `ruv_cancorplot` (the output of `ruvCancor` for a multi-assay
#'   object).
#'
#' @return A named list where keys are the assay names (from the input list)
#'   and values are the determined best K for each assay. Returns `NA_integer_`
#'   for assays where `findBestK` fails or returns an invalid result.
#'
#' @importFrom purrr map set_names
#' @importFrom logger log_warn
#' @export
findBestKForAssayList <- function(cancor_plots_list) {

    # --- Input Validation ---
    if (!is.list(cancor_plots_list)) {
        stop("`cancor_plots_list` must be a list.")
    }
    if (length(cancor_plots_list) == 0) {
        log_warn("`cancor_plots_list` is empty. Returning an empty list.")
        return(list())
    }
    if (is.null(names(cancor_plots_list)) || any(names(cancor_plots_list) == "")) {
        log_warn("Input list is not fully named. Assigning default names (Plot_1, Plot_2, ...)")
        names(cancor_plots_list) <- paste0("Plot_", seq_along(cancor_plots_list))
    }

    # --- Iterate and Apply findBestK ---
    best_k_list <- purrr::map(seq_along(cancor_plots_list), function(i) {
        assay_name <- names(cancor_plots_list)[i]
        current_plot <- cancor_plots_list[[i]]

        # Check if it's a ggplot object
        if (!inherits(current_plot, "ggplot")) {
            log_warn("Element '{assay_name}' in the list is not a ggplot object. Skipping.", .logr = TRUE)
            return(NA_integer_)
        }

        # Call the original findBestK function
        best_k_assay <- tryCatch({
            k <- findBestK(current_plot)
            # Add a check for Inf or non-numeric result which might cause issues later
            if (is.infinite(k) || !is.numeric(k) || length(k) != 1) {
                 log_warn("Assay '{assay_name}': findBestK returned an invalid value ({k}). Returning NA.", .logr = TRUE)
                 NA_integer_
            } else {
                 # Ensure it's an integer
                 as.integer(k)
            }
        }, warning = function(w) {
            # Catch the specific warning from findBestK if max returns -Inf
            if (grepl("no non-missing arguments to max", w$message, ignore.case = TRUE)) {
                log_warn("Assay '{assay_name}': Warning in findBestK (likely no valid data in plot): {w$message}. Returning NA.", .logr = TRUE)
            } else {
                # Log other warnings but still try to return NA
                log_warn("Assay '{assay_name}': Warning during findBestK: {w$message}. Returning NA.", .logr = TRUE)
            }
            NA_integer_
        }, error = function(e) {
            log_warn("Assay '{assay_name}': Error calling findBestK: {e$message}. Returning NA.", .logr = TRUE)
            NA_integer_
        })

        return(best_k_assay)
    })

    # Set names for the result list
    names(best_k_list) <- names(cancor_plots_list)

    return(best_k_list)
}


# ----------------------------------------------------------------------------
# findBestNegCtrlPercentage
# ----------------------------------------------------------------------------
#' Find the Best Negative Control Percentage for RUV-III Analysis
#'
#' This function automatically determines the optimal percentage of proteins to use
#' as negative controls for RUV-III analysis by testing different percentages and
#' evaluating the separation quality between "All" and "Control" groups in canonical
#' correlation plots.
#'
#' @param normalised_protein_matrix_obj A ProteinQuantitativeData object containing
#'   the normalized protein quantification data
#' @param percentage_range A numeric vector specifying the range of percentages to test.
#'   Default is seq(1, 20, by = 1) for testing 1% to 20% in 1% increments
#' @param num_components_to_impute Number of components to use for imputation in ruvCancor.
#'   Default is 5
#' @param ruv_grouping_variable The grouping variable to use for RUV analysis.
#'   Default is "group"
#' @param ruv_qval_cutoff The FDR threshold for negative control selection.
#'   Default is 0.05
#' @param ruv_fdr_method The FDR calculation method. Default is "qvalue"
#' @param separation_metric The metric to use for evaluating separation quality.
#'   Options: "max_difference" (default), "mean_difference", "auc", "weighted_difference"
#' @param k_penalty_weight Weight for penalizing high k values in composite score.
#'   Default is 0.5. Higher values penalize high k more strongly
#' @param max_acceptable_k Maximum acceptable k value. k values above this get heavy penalty.
#'   Default is 3
#' @param adaptive_k_penalty Whether to automatically adjust max_acceptable_k based on sample size.
#'   Default is TRUE (recommended). Set to FALSE only if you need exact reproducibility with previous results
#' @param verbose Whether to print progress messages. Default is TRUE
#'
#' @return A list containing:
#'   \itemize{
#'     \item best_percentage: The optimal percentage as a numeric value
#'     \item best_k: The optimal k value from findBestK() for the best percentage
#'     \item best_control_genes_index: The control genes index for the best percentage
#'     \item best_separation_score: The separation score for the best percentage
#'     \item best_composite_score: The composite score (separation penalized by k value)
#'     \item optimization_results: A data frame with all tested percentages and their scores
#'     \item best_cancor_plot: The canonical correlation plot for the best percentage
#'     \item separation_metric_used: The separation metric that was used
#'     \item k_penalty_weight: The k penalty weight that was used
#'     \item max_acceptable_k: The maximum acceptable k value that was used
#'   }
#'
#' @importFrom logger log_info log_warn
#' @importFrom purrr imap map_dfr map_dbl
#' @export
findBestNegCtrlPercentage <- function(normalised_protein_matrix_obj,
                                      percentage_range = seq(1, 20, by = 1),
                                      num_components_to_impute = 5,
                                      ruv_grouping_variable = "group",
                                      ruv_qval_cutoff = 0.05,
                                      ruv_fdr_method = "qvalue",
                                      separation_metric = "max_difference",
                                      k_penalty_weight = 0.5,
                                      max_acceptable_k = 3,
                                      adaptive_k_penalty = TRUE,
                                      verbose = TRUE) {
  
  # Input validation
  if (!inherits(normalised_protein_matrix_obj, "ProteinQuantitativeData")) {
    stop("normalised_protein_matrix_obj must be a ProteinQuantitativeData object")
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
    # Get sample size from the data object
    sample_size <- ncol(normalised_protein_matrix_obj@protein_quant_table) - 1  # Subtract 1 for protein ID column
    
    # Calculate adaptive max_acceptable_k based on sample size
    # Conservative approach: roughly 1 k per 10-15 samples, but with reasonable bounds
    adaptive_max_k <- calculateAdaptiveMaxK(sample_size)
    
    if (verbose) {
      log_info("Adaptive penalty enabled: Sample size = {sample_size}, Adaptive max_acceptable_k = {adaptive_max_k} (original = {max_acceptable_k})")
    }
    
    max_acceptable_k <- adaptive_max_k
  }
  
  # Detect small datasets and warn/adjust percentage range if needed
  sample_size <- ncol(normalised_protein_matrix_obj@protein_quant_table) - 1
  if (sample_size < 15 && max(percentage_range) < 30) {
    if (verbose) {
      log_warn("Small dataset detected (n={sample_size}). Consider testing higher percentages (up to 30-50%) for better negative control identification.")
      log_warn("Current range: {paste(range(percentage_range), collapse = '-')}%. May need wider range for optimal results.")
    }
  }
  
  if (verbose) {
    log_info("Starting optimization of negative control percentage with k value consideration...")
    log_info("Testing {length(percentage_range)} different percentages: {paste(range(percentage_range), collapse = '-')}%")
    log_info("K penalty weight: {k_penalty_weight}, Max acceptable k: {max_acceptable_k}")
    if (adaptive_k_penalty) {
      log_info("Using adaptive k penalty based on sample size")
    }
  }
  
  # Process all percentages using functional programming (proper R way)
  if (verbose) {
    log_info("Processing {length(percentage_range)} percentages using vectorized operations...")
  }
  
  # Create a function to process a single percentage
  process_percentage <- function(current_percentage, index) {
    if (verbose && index %% 5 == 0) {
      log_info("Testing percentage {index}/{length(percentage_range)}: {current_percentage}%")
    }
    
    tryCatch({
      # Get negative control genes for current percentage
      control_genes_index <- getNegCtrlProtAnova(
        normalised_protein_matrix_obj,
        ruv_grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = current_percentage,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
      
      # Check if we have enough control genes
      num_controls <- sum(control_genes_index, na.rm = TRUE)
      if (num_controls < 5) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Only {num_controls} control genes found (minimum 5 required). Skipping.")
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
      
      # Generate canonical correlation plot
      cancorplot <- ruvCancor(
        normalised_protein_matrix_obj,
        ctrl = control_genes_index,
        num_components_to_impute = num_components_to_impute,
        ruv_grouping_variable = ruv_grouping_variable
      )
      
      # Calculate separation score
      separation_score <- calculateSeparationScore(cancorplot, separation_metric)
      
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
      composite_score <- calculateCompositeScore(
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
  
  # Use purrr::imap() for functional processing (the R way)
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
    log_info("  - Number of control genes: {sum(best_control_genes_index, na.rm = TRUE)}")
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
    sample_size = if(adaptive_k_penalty) ncol(normalised_protein_matrix_obj@protein_quant_table) - 1 else NA
  ))
}


# ----------------------------------------------------------------------------
# extractRuvResults
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
extractRuvResults <- function(results_list) {

  extracted <- purrr::map(results_list, \(x){ x$results })

  names(extracted) <- names(results_list)

  return(extracted)
}


# ----------------------------------------------------------------------------
# scaleCenterAndFillMissing
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
#' @title Scale, Center and Fill Missing Values
#'@param input_matrix Samples are columns, rows are proteins or peptides
scaleCenterAndFillMissing <- function( input_matrix) {

  input_matrix_scaled <- scale(input_matrix, center = TRUE, scale = TRUE)


  min_data_point <- min(input_matrix_scaled, na.rm=TRUE)


  input_matrix_scaled_fill_missing <-  input_matrix_scaled
  input_matrix_scaled_fill_missing[is.na(input_matrix_scaled_fill_missing)] <- min_data_point*2

  input_matrix_scaled_fill_missing
}


# ----------------------------------------------------------------------------
# calculateSeparationScore
# ----------------------------------------------------------------------------
#' Calculate Separation Score for Canonical Correlation Plot
#'
#' Internal helper function to calculate separation quality between "All" and "Control"
#' groups in a canonical correlation plot.
#'
#' @param cancorplot A ggplot object from ruvCancor
#' @param metric The separation metric to calculate
#'
#' @return A numeric separation score (higher is better)
#'
#' @keywords internal
calculateSeparationScore <- function(cancorplot, metric = "max_difference") {
  
  # Extract data from the plot
  if (!inherits(cancorplot, "ggplot") || is.null(cancorplot$data)) {
    return(NA_real_)
  }
  
  plot_data <- cancorplot$data
  
  # Check required columns exist
  if (!all(c("featureset", "cc", "K") %in% colnames(plot_data))) {
    return(NA_real_)
  }
  
  # Get indices for Control and All groups
  controls_idx <- which(plot_data$featureset == "Control")
  all_idx <- which(plot_data$featureset == "All")
  
  if (length(controls_idx) == 0 || length(all_idx) == 0) {
    return(NA_real_)
  }
  
  # Calculate differences between All and Control
  difference_between_all_ctrl <- plot_data$cc[all_idx] - plot_data$cc[controls_idx]
  
  # Remove any NA or infinite values
  valid_diffs <- difference_between_all_ctrl[is.finite(difference_between_all_ctrl)]
  
  if (length(valid_diffs) == 0) {
    return(NA_real_)
  }
  
  # Calculate score based on specified metric
  score <- switch(metric,
    "max_difference" = max(valid_diffs, na.rm = TRUE),
    "mean_difference" = mean(valid_diffs, na.rm = TRUE),
    "auc" = {
      # Area under the curve (trapezoidal rule approximation)
      k_values <- plot_data$K[all_idx][is.finite(difference_between_all_ctrl)]
      if (length(k_values) < 2) return(NA_real_)
      
      # Sort by K value
      sorted_idx <- order(k_values)
      k_sorted <- k_values[sorted_idx]
      diff_sorted <- valid_diffs[sorted_idx]
      
      # Calculate AUC using trapezoidal rule
      sum(diff(k_sorted) * (head(diff_sorted, -1) + tail(diff_sorted, -1)) / 2)
    },
    "weighted_difference" = {
      # Weight differences by their K value (higher K gets more weight)
      k_values <- plot_data$K[all_idx][is.finite(difference_between_all_ctrl)]
      if (length(k_values) == 0) return(NA_real_)
      
      weights <- k_values / max(k_values, na.rm = TRUE)
      sum(valid_diffs * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
    },
    NA_real_
  )
  
  return(as.numeric(score))
}


# ----------------------------------------------------------------------------
# calculateCompositeScore
# ----------------------------------------------------------------------------
#' Calculate Composite Score for Percentage Optimization
#'
#' Internal helper function to calculate a composite score that considers both
#' separation quality and the resulting k value from findBestK(). This prevents
#' over-optimization towards percentages that give good separation but unreasonably
#' high k values that would remove biological signal.
#'
#' @param separation_score The separation score from calculateSeparationScore()
#' @param best_k The best k value from findBestK()
#' @param k_penalty_weight Weight for k penalty (0-1). Higher values penalize high k more
#' @param max_acceptable_k Maximum acceptable k value. k values above this get heavy penalty
#'
#' @return A numeric composite score (higher is better)
#'
#' @keywords internal
calculateCompositeScore <- function(separation_score, best_k, k_penalty_weight, max_acceptable_k) {
  
  # Handle NA cases
  if (is.na(separation_score) || is.na(best_k)) {
    return(NA_real_)
  }
  
  # Ensure positive values
  if (separation_score <= 0) {
    return(0)
  }
  
  # Calculate k penalty
  if (best_k <= max_acceptable_k) {
    # Linear penalty within acceptable range: penalty = 0 at k=1, penalty = k_penalty_weight at k=max_acceptable_k
    k_penalty <- k_penalty_weight * (best_k - 1) / (max_acceptable_k - 1)
  } else {
    # Heavy penalty for k values above max_acceptable_k
    # Exponential penalty: starts at k_penalty_weight and increases rapidly
    excess_k <- best_k - max_acceptable_k
    k_penalty <- k_penalty_weight + (1 - k_penalty_weight) * (1 - exp(-excess_k))
  }
  
  # Ensure k_penalty is between 0 and 1
  k_penalty <- pmax(0, pmin(1, k_penalty))
  
  # Calculate composite score: separation_score * (1 - k_penalty)
  # This means:
  # - k=1: no penalty (multiply by 1)
  # - k=max_acceptable_k: penalty = k_penalty_weight (multiply by 1-k_penalty_weight)
  # - k>max_acceptable_k: heavy penalty (multiply by value approaching 0)
  composite_score <- separation_score * (1 - k_penalty)
  
  return(as.numeric(composite_score))
}


# ----------------------------------------------------------------------------
# calculateAdaptiveMaxK
# ----------------------------------------------------------------------------
#' Calculate Adaptive Maximum Acceptable K Based on Sample Size
#'
#' Internal helper function to determine an appropriate max_acceptable_k value
#' based on the number of samples in the dataset. This prevents over-correction
#' in small datasets and allows more flexibility in large datasets.
#'
#' @param sample_size Number of samples in the dataset
#'
#' @return An integer representing the adaptive max_acceptable_k value
#'
#' @details
#' The adaptive calculation follows these principles:
#' - Small datasets (n < 15): Conservative approach, max_k = 2
#' - Medium datasets (n = 15-40): Standard approach, max_k = 3  
#' - Large datasets (n = 40-80): Moderate approach, max_k = 4
#' - Very large datasets (n > 80): Permissive approach, max_k = 5
#' 
#' The rationale is that with more samples, you have more degrees of freedom
#' and statistical power, making higher k values less problematic.
#'
#' @keywords internal
calculateAdaptiveMaxK <- function(sample_size) {
  
  if (sample_size < 15) {
    # Small datasets: be very conservative
    # Each k factor consumes significant degrees of freedom
    return(2L)
  } else if (sample_size < 40) {
    # Medium datasets: standard approach
    # This is the typical proteomics experiment size
    return(3L)
  } else if (sample_size < 80) {
    # Large datasets: can afford one extra k factor
    # Sufficient statistical power to handle k=4
    return(4L)
  } else {
    # Very large datasets: most permissive
    # Abundant statistical power allows k=5 if separation justifies it
    return(5L)
  }
}

