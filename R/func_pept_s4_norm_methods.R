# Internal wrapper for ruvCancor S4 dispatch (mockable in tests without corrupting S4 method table)
#' @keywords internal
.call_ruvCancor <- function(obj, ctrl, num_comp, grp_var) {
  ruvCancor(obj, ctrl = ctrl, num_components_to_impute = num_comp, ruv_grouping_variable = grp_var)
}

## normalise between Arrays
#'@export
#'@param theObject Object of class PeptideQuantitativeData
#'@param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
#' @title Normalize Between Samples for PeptideQuantitativeData
#' @name normaliseBetweenSamples,PeptideQuantitativeData-method
#' @export
setMethod(f="normaliseBetweenSamples"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, normalisation_method = NULL) {
            peptide_data <- theObject@peptide_data
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            current_quant_column <- theObject@norm_quantity_column

            normalisation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                      , "normalisation_method"
                                                                      , "cyclicloess")

            theObject <- updateParamInObject(theObject, "normalisation_method")

            # Create matrix from peptide_data (like protein version from protein_quant_table)
            # Create unique peptide identifier for matrix rownames
            peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]],
                                      peptide_data[[theObject@peptide_sequence_column]],
                                      sep = "%")

            # Create wide format data frame then convert to matrix (exactly like protein version)
            temp_peptide_wide <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
              pivot_wider(names_from = !!sym(sample_id),
                         values_from = !!sym(current_quant_column),
                         values_fn = mean)

            # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
            frozen_peptide_matrix <- temp_peptide_wide |>
              column_to_rownames("peptide_id") |>
              as.matrix()

            frozen_peptide_matrix[!is.finite(frozen_peptide_matrix)] <- NA

            normalised_frozen_peptide_matrix <- frozen_peptide_matrix

            print(paste0("normalisation_method = ", normalisation_method))

            switch(normalisation_method
                   , cyclicloess = {
                     normalised_frozen_peptide_matrix <- limma::normalizeCyclicLoess(frozen_peptide_matrix)
                   }
                   , quantile = {
                     normalised_frozen_peptide_matrix <- limma::normalizeQuantiles(frozen_peptide_matrix)
                   }
                   , scale = {
                     normalised_frozen_peptide_matrix <- limma::normalizeMedianAbsValues(frozen_peptide_matrix)
                   }
                   , none = {
                     normalised_frozen_peptide_matrix <- frozen_peptide_matrix
                   }
            )

            normalised_frozen_peptide_matrix[!is.finite(normalised_frozen_peptide_matrix)] <- NA

            # Convert back to data frame (exactly like protein version)
            normalised_peptide_table <- normalised_frozen_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_id")

            # Update peptide_data by joining with normalized values and replacing the quant column
            updated_peptide_data <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(-!!sym(current_quant_column)) |>
              left_join(normalised_peptide_table |>
                       pivot_longer(cols = -peptide_id,
                                   names_to = sample_id,
                                   values_to = current_quant_column),
                       by = c("peptide_id", sample_id)) |>
              select(-peptide_id)

            # Update both slots
            theObject@peptide_data <- updated_peptide_data
            theObject@peptide_matrix <- normalised_frozen_peptide_matrix

            theObject <- cleanDesignMatrixPeptide(theObject)

            return(theObject)
          })

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
              mutate(!!sym(current_quant_column) := log2_value) |>
              select(-peptide_row_id, -log2_value)

            theObject@peptide_data <- updated_peptide_data

            # Mark as logged
            theObject@is_logged_data <- TRUE

            theObject <- cleanDesignMatrixPeptide(theObject)

            message("Peptide data successfully log2 transformed. Raw intensities converted to log2 scale.")
            message(paste("is_logged_data flag set to:", theObject@is_logged_data))

            return(theObject)
          })

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

            design_for_anova <- design_matrix |>
              tibble::column_to_rownames(sample_id_column)

            if (group_id != ruv_grouping_variable) {
              design_for_anova <- design_for_anova |>
                dplyr::select(-dplyr::all_of(group_id))
            }

            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id_column)) ]
                                                              , design_matrix = design_for_anova
                                                              , grouping_variable = ruv_grouping_variable
                                                              , percentage_as_neg_ctrl = percentage_as_neg_ctrl
                                                              , num_neg_ctrl = num_neg_ctrl
                                                              , ruv_qval_cutoff = ruv_qval_cutoff
                                                              , ruv_fdr_method = ruv_fdr_method )

            return(control_genes_index)
          })

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

  # Emit weighted_difference deprecation warning once per public call
  if (identical(separation_metric, "weighted_difference")) {
    warning(
      paste0(
        "'weighted_difference' is deprecated because it up-weights high-K values ",
        "while the composite score penalizes high K. ",
        "Use 'max_difference' or 'mean_difference' instead."
      ),
      call. = FALSE
    )
  }

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

  # Pre-compute quantities constant across all percentages
  candidate_feature_count <- nrow(theObject@peptide_matrix)
  sample_size <- ncol(theObject@peptide_matrix)

  # Calculate adaptive max_acceptable_k if requested
  if (adaptive_k_penalty) {
    adaptive_max_k <- .peptide_calculateAdaptiveMaxK(sample_size)

    if (verbose) {
      log_info("Adaptive penalty enabled: Sample size = {sample_size}, Adaptive max_acceptable_k = {adaptive_max_k} (original = {max_acceptable_k})")
    }

    max_acceptable_k <- adaptive_max_k
  }

  # Detect small datasets and warn if needed
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
    log_info("Processing {length(percentage_range)} percentages using vectorized operations...")
  }

  # Helper: build a trace row for a failed or skipped percentage
  .make_failed_row <- function(current_percentage, status, error_reason = "",
                                realized_num_controls = NA_integer_,
                                realized_percentage = NA_real_,
                                separation_score = NA_real_,
                                best_k = NA_integer_,
                                composite_score = NA_real_,
                                control_genes_index = NULL,
                                cancor_plot = NULL) {
    list(
      percentage_requested    = current_percentage,
      candidate_feature_count = candidate_feature_count,
      realized_num_controls   = realized_num_controls,
      realized_percentage     = realized_percentage,
      sample_size             = sample_size,
      best_k                  = best_k,
      separation_score        = separation_score,
      composite_score         = composite_score,
      status                  = status,
      error_reason            = error_reason,
      control_genes_index     = control_genes_index,
      cancor_plot             = cancor_plot
    )
  }

  # Process a single percentage, returning a fully annotated trace record
  process_percentage <- function(current_percentage, index) {
    if (verbose && index %% 5 == 0) {
      log_info("Testing percentage {index}/{length(percentage_range)}: {current_percentage}%")
    }

    # Step 1: Select negative controls
    ctrl_error <- NULL
    control_genes_index <- tryCatch({
      getNegCtrlProtAnovaPeptides(
        theObject,
        ruv_grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = current_percentage,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
    }, error = function(e) {
      ctrl_error <<- e$message
      NULL
    })

    if (is.null(control_genes_index)) {
      if (verbose) {
        log_warn(paste("Percentage", current_percentage, "% : Error in getNegCtrlProtAnovaPeptides:", ctrl_error))
      }
      return(.make_failed_row(current_percentage, "error_neg_ctrl_selection",
        error_reason = if (is.null(ctrl_error)) "" else ctrl_error))
    }

    # Step 2: Check minimum control count
    realized_num_controls <- sum(control_genes_index, na.rm = TRUE)
    realized_percentage <- if (candidate_feature_count > 0L) {
      100 * realized_num_controls / candidate_feature_count
    } else {
      NA_real_
    }

    if (realized_num_controls < 5L) {
      if (verbose) {
        log_warn("Percentage {current_percentage}%: Only {realized_num_controls} control peptides found (minimum 5 required). Skipping.")
      }
      return(.make_failed_row(current_percentage, "skipped_insufficient_controls",
        realized_num_controls = realized_num_controls,
        realized_percentage = realized_percentage))
    }

    # Step 3: Canonical correlation (full objective, not surrogate)
    cancor_error <- NULL
    cancorplot <- tryCatch({
      .call_ruvCancor(
        theObject,
        ctrl = control_genes_index,
        num_comp = num_components_to_impute,
        grp_var = ruv_grouping_variable
      )
    }, error = function(e) {
      cancor_error <<- e$message
      NULL
    })

    if (is.null(cancorplot)) {
      if (verbose) {
        log_warn(paste("Percentage", current_percentage, "% : Error in ruvCancor:", cancor_error))
      }
      return(.make_failed_row(current_percentage, "error_cancor",
        realized_num_controls = realized_num_controls,
        realized_percentage = realized_percentage,
        error_reason = if (is.null(cancor_error)) "" else cancor_error,
        control_genes_index = control_genes_index))
    }

    # Step 4: Score the result
    separation_score <- .peptide_calculateSeparationScore(cancorplot, separation_metric)

    best_k <- tryCatch(
      findBestKElbow(cancorplot),
      error = function(e) NA_integer_
    )

    if (is.na(best_k)) {
      return(.make_failed_row(current_percentage, "invalid_cancor_plot",
        realized_num_controls = realized_num_controls,
        realized_percentage = realized_percentage,
        separation_score = separation_score,
        error_reason = "findBestKElbow returned NA",
        control_genes_index = control_genes_index,
        cancor_plot = cancorplot))
    }

    best_k <- as.integer(best_k)
    composite_score <- .peptide_calculateCompositeScore(
      separation_score,
      best_k,
      k_penalty_weight,
      max_acceptable_k
    )

    if (!is.finite(composite_score)) {
      return(.make_failed_row(current_percentage, "error_scoring",
        realized_num_controls = realized_num_controls,
        realized_percentage = realized_percentage,
        separation_score = separation_score,
        best_k = best_k,
        composite_score = composite_score,
        error_reason = "composite score is non-finite",
        control_genes_index = control_genes_index,
        cancor_plot = cancorplot))
    }

    list(
      percentage_requested    = current_percentage,
      candidate_feature_count = candidate_feature_count,
      realized_num_controls   = realized_num_controls,
      realized_percentage     = realized_percentage,
      sample_size             = sample_size,
      best_k                  = best_k,
      separation_score        = separation_score,
      composite_score         = composite_score,
      status                  = "ok",
      error_reason            = "",
      control_genes_index     = control_genes_index,
      cancor_plot             = cancorplot
    )
  }

  all_results <- percentage_range |>
    purrr::imap(process_percentage)

  # Build optimization trace: one row per tested percentage
  # New schema primary; deprecated aliases (percentage, num_controls, valid_plot) kept for one release cycle
  results <- all_results |>
    purrr::map_dfr(~ data.frame(
      percentage_requested    = .x$percentage_requested,
      candidate_feature_count = .x$candidate_feature_count,
      realized_num_controls   = .x$realized_num_controls,
      realized_percentage     = .x$realized_percentage,
      sample_size             = .x$sample_size,
      best_k                  = .x$best_k,
      separation_score        = .x$separation_score,
      composite_score         = .x$composite_score,
      status                  = .x$status,
      error_reason            = .x$error_reason,
      percentage              = .x$percentage_requested,
      num_controls            = .x$realized_num_controls,
      valid_plot              = identical(.x$status, "ok"),
      stringsAsFactors = FALSE
    ))

  # Deterministic winner selection — stable ordering replaces which.max()
  valid_df <- results |>
    dplyr::filter(is.finite(composite_score)) |>
    dplyr::arrange(
      dplyr::desc(composite_score),
      best_k,
      dplyr::desc(separation_score),
      percentage_requested
    )

  if (nrow(valid_df) == 0L) {
    stop("No valid percentage found. Please check your data and parameters.")
  }

  best_percentage <- valid_df$percentage_requested[[1L]]

  # Retrieve the full record (with control_genes_index and cancor_plot) for the winner
  winner_idx <- which(
    vapply(all_results, function(r) r$percentage_requested == best_percentage, logical(1L))
  )[[1L]]
  best_result <- all_results[[winner_idx]]

  best_k                   <- as.integer(best_result$best_k)
  best_control_genes_index <- best_result$control_genes_index
  best_cancor_plot         <- best_result$cancor_plot
  best_separation_score    <- best_result$separation_score
  best_composite_score     <- best_result$composite_score

  if (verbose) {
    log_info("Optimization complete!")
    log_info("Best percentage: {best_percentage}% (composite score: {round(best_composite_score, 4)})")
    log_info("  - Separation score: {round(best_separation_score, 4)}")
    log_info("  - Best k value: {best_k}")
    log_info("  - Number of control peptides: {sum(best_control_genes_index, na.rm = TRUE)}")
  }

  list(
    best_percentage          = best_percentage,
    best_k                   = best_k,
    best_control_genes_index = best_control_genes_index,
    best_separation_score    = best_separation_score,
    best_composite_score     = best_composite_score,
    optimization_results     = results,
    best_cancor_plot         = best_cancor_plot,
    separation_metric_used   = separation_metric,
    k_penalty_weight         = k_penalty_weight,
    max_acceptable_k         = max_acceptable_k,
    adaptive_k_penalty_used  = adaptive_k_penalty,
    sample_size              = sample_size
  )
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
  curve <- .extractCancorCurveData(cancorplot)
  if (is.null(curve)) {
    return(NA_real_)
  }

  valid_rows <- curve[is.finite(curve$delta), , drop = FALSE]
  if (nrow(valid_rows) == 0) {
    return(NA_real_)
  }

  score <- switch(metric,
    "max_difference" = max(valid_rows$delta),
    "mean_difference" = mean(valid_rows$delta),
    "auc" = {
      if (nrow(valid_rows) < 2) return(NA_real_)
      sorted <- valid_rows[order(valid_rows$K), ]
      sum(diff(sorted$K) * (head(sorted$delta, -1) + tail(sorted$delta, -1)) / 2)
    },
    "weighted_difference" = {
      weights <- valid_rows$K / max(valid_rows$K)
      sum(valid_rows$delta * weights) / sum(weights)
    },
    NA_real_
  )

  as.numeric(score)
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
  if (!is.finite(separation_score) || !is.finite(best_k)) {
    return(NA_real_)
  }

  if (!is.finite(max_acceptable_k) || max_acceptable_k < 1) {
    return(NA_real_)
  }

  if (separation_score <= 0) {
    return(0)
  }

  if (max_acceptable_k == 1) {
    k_penalty <- if (best_k <= 1) 0 else 1
  } else if (best_k <= max_acceptable_k) {
    k_penalty <- k_penalty_weight * (best_k - 1) / (max_acceptable_k - 1)
  } else {
    excess_k <- best_k - max_acceptable_k
    k_penalty <- k_penalty_weight + (1 - k_penalty_weight) * (1 - exp(-excess_k))
  }

  k_penalty <- pmax(0, pmin(1, k_penalty))

  as.numeric(separation_score * (1 - k_penalty))
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

#' @title Fast version of ruvCancor for optimization
#' @name ruvCancorFast,PeptideQuantitativeData-method
#' @description This function provides a lightweight alternative to ruvCancor that skips
#' the expensive DPC imputation step, making it suitable for use during
#' optimization loops where speed is critical.
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Control features for RUV
#' @param num_components_to_impute Number of components to impute
#' @param ruv_grouping_variable Grouping variable for RUV
#' @param simple_imputation_method Method for simple missing value handling.
#'   Options: "none" (default), "mean", "median", "min"
#' @export
setMethod( f = "ruvCancorFast"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL,
                                  ruv_grouping_variable = NULL, simple_imputation_method = "none") {

             .Deprecated("ruvCancor")

             peptide_matrix <- theObject@peptide_matrix
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             sample_id <- theObject@sample_id

             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }

             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }

             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             # Skip expensive DPC imputation - use simple approach for optimization
             peptide_matrix_working <- log2(peptide_matrix + 1)  # Add small constant to avoid log(0)

             # Handle missing values with simple method if requested
             if (simple_imputation_method != "none" && anyNA(peptide_matrix_working)) {
               if (simple_imputation_method == "mean") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- mean(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "median") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- median(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "min") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- min(x, na.rm = TRUE)
                   return(x)
                 })
               }
             }

             # Remove or winsorize extreme values (lightweight version)
             peptide_matrix_clean <- apply(peptide_matrix_working, 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })

             # Use the matrix directly without expensive DPC imputation
             Y <- t(peptide_matrix_clean[, design_matrix |> dplyr::pull(!!sym(sample_id))])

             # Generate canonical correlation plot (this is what we actually need for optimization)
             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)

             return(cancorplot_r2)
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' RUV Canonical Correlation Analysis for Peptide Data
#'
#' Performs Remove Unwanted Variation (RUV) using canonical correlation analysis
#' with DPC imputation for missing values.
#'
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Vector of control gene indices for RUV
#' @param num_components_to_impute Number of principal components for imputation
#' @param ruv_grouping_variable Column name in design matrix for RUV grouping
#'
#' @export
#' @title RUV Canonical Correlation for PeptideQuantitativeData
#' @name ruvCancor,PeptideQuantitativeData-method
#' @export
setMethod( f = "ruvCancor"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {


             peptide_matrix <- theObject@peptide_matrix
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             sample_id <- theObject@sample_id

             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)

             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }

             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }

             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             # Remove or winsorize extreme values
             peptide_matrix_clean <- apply(log2(peptide_matrix), 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })


             dpcfit <- dpc(peptide_matrix_clean)

             peptide_matrix_complete <- dpcImpute(peptide_matrix_clean, dpc=dpcfit)$E


             Y <-  t( peptide_matrix_complete[,design_matrix |> dplyr::pull(!!sym(sample_id))])

             print("steps")

             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2


           })

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {

             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id

             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             k <- checkParamsObjectFunctionSimplify( theObject, "ruv_number_k", NULL)
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)

             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "ruv_number_k")
             theObject <- updateParamInObject(theObject, "ctrl")

             # Create matrix from peptide_data (exactly like protein version)
             current_quant_column <- theObject@norm_quantity_column

             # Create unique peptide identifier for matrix rownames
             peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]],
                                       peptide_data[[theObject@peptide_sequence_column]],
                                       sep = "%")

             # Create wide format data frame then convert to matrix (exactly like protein version)
             temp_peptide_wide <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
               pivot_wider(names_from = !!sym(sample_id),
                          values_from = !!sym(current_quant_column),
                          values_fn = mean)

             # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
             normalised_frozen_peptide_matrix_filt <- temp_peptide_wide |>
               column_to_rownames("peptide_id") |>
               as.matrix()

             Y <-  t( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])

             M <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                  , !!sym(sample_id)
                                                  , !!sym(ruv_grouping_variable))

             cln_mat <- RUVIII_C_Varying( k = ruv_number_k
                                          , Y = Y
                                          , M = M
                                          , toCorrect = colnames(Y)
                                          , potentialControls = names( ctrl[which(ctrl)] ) )

             # Remove samples with no values
             cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat),]

             # Remove proteins with no values
             cln_mat_3 <- t(cln_mat_2)
             cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3),]

             # Convert back to data frame (exactly like protein version)
             ruv_normalised_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column("peptide_id")

             # Update peptide_data by joining with RUV-corrected values and replacing the quant column
             updated_peptide_data <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(-!!sym(current_quant_column)) |>
               left_join(ruv_normalised_results_cln |>
                        pivot_longer(cols = -peptide_id,
                                    names_to = sample_id,
                                    values_to = current_quant_column),
                        by = c("peptide_id", sample_id)) |>
               select(-peptide_id)

             # Update both slots
             theObject@peptide_data <- updated_peptide_data
             theObject@peptide_matrix <- cln_mat_4

             theObject <- cleanDesignMatrixPeptide(theObject)

             return( theObject )
           })
