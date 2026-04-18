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
  
  message("--- DEBUG66: Entering findBestNegCtrlPercentage ---")
  message(sprintf("   DEBUG66: percentage_range: %s", paste(range(percentage_range), collapse="-")))
  message(sprintf("   DEBUG66: adaptive_k_penalty: %s", adaptive_k_penalty))
  message(sprintf("   DEBUG66: separation_metric: %s", separation_metric))
  
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
    message(sprintf("   DEBUG66 Loop: Processing percentage %d%%", current_percentage))
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
      message(sprintf("      DEBUG66: num_controls = %d", num_controls))
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
      message(sprintf("      DEBUG66: separation_score = %.4f", separation_score))
      
      # Calculate the best k using the existing findBestK function
      best_k <- tryCatch({
        findBestK(cancorplot)
      }, error = function(e) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Error calculating best k: {e$message}")
        }
        return(NA_real_)
      })
      message(sprintf("      DEBUG66: best_k = %s", as.character(best_k)))
      
      # Calculate composite score that considers both separation and k value
      composite_score <- calculateCompositeScore(
        separation_score, 
        best_k, 
        k_penalty_weight, 
        max_acceptable_k
      )
      message(sprintf("      DEBUG66: composite_score = %.4f", composite_score))
      
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
      message(sprintf("      DEBUG66 ERROR in process_percentage: %s", e$message))
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
  
  message("--- DEBUG66: Exiting findBestNegCtrlPercentage ---")

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

