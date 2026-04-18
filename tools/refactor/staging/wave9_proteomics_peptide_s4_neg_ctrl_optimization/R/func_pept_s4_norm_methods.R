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

#' Calculate Separation Score for Canonical Correlation Plot (Peptide version)
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
.peptide_calculateSeparationScore <- function(cancorplot, metric = "max_difference") {
  
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

#' Calculate Composite Score for Percentage Optimization (Peptide version)
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
.peptide_calculateCompositeScore <- function(separation_score, best_k, k_penalty_weight, max_acceptable_k) {
  
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

#' Calculate Adaptive Maximum Acceptable K Based on Sample Size (Peptide version)
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
.peptide_calculateAdaptiveMaxK <- function(sample_size) {
  
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

