# Vendored verbatim from APAF-bioinformatics/e_35036_MorpheusTex_GarryNunn_20260329_b (private), commit 4aa97fc3b85e897ee457052682d8bba2a7afbbc4,
# file scripts/proteomics_isoseq/DIA_workflow_limpa_isoseq_tissue.Rmd, lines 1086-1172 (fetched 2026-10-21 for KIT-041;
# the reference the package optimiser findBestNegCtrlPercentageCapped is checked against). Not package code.
findBestNegCtrlPercentageCustom <- function(normalised_protein_matrix_obj,
                                            percentage_range = seq(1, 30, by = 1),
                                            num_components_to_impute = 5,
                                            ruv_grouping_variable = "group",
                                            ruv_qval_cutoff = 0.05,
                                            ruv_fdr_method = "BH",
                                            separation_metric = "max_difference",
                                            k_penalty_weight = 0.5,
                                            max_acceptable_k = 3,
                                            adaptive_k_penalty = TRUE,
                                            verbose = TRUE) {
  sample_size <- ncol(normalised_protein_matrix_obj@protein_quant_table) - 1
  if (adaptive_k_penalty) {
    max_acceptable_k <- MultiScholaR:::calculateAdaptiveMaxK(sample_size)
  }
  
  process_percentage <- function(current_percentage, index) {
    if (verbose && index %% 5 == 0) {
      message(paste0("Testing percentage ", index, "/", length(percentage_range), ": ", current_percentage, "%"))
    }
    tryCatch({
      control_genes_index <- getNegCtrlProtAnova(
        normalised_protein_matrix_obj,
        ruv_grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = current_percentage,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
      num_controls <- sum(control_genes_index, na.rm = TRUE)
      if (num_controls < 5) {
        return(list(percentage = current_percentage, separation_score = NA_real_, best_k = NA_real_,
                    composite_score = NA_real_, num_controls = num_controls, valid_plot = FALSE,
                    control_genes_index = NULL, cancor_plot = NULL))
      }
      cancorplot <- ruvCancor(
        normalised_protein_matrix_obj,
        ctrl = control_genes_index,
        num_components_to_impute = num_components_to_impute,
        ruv_grouping_variable = ruv_grouping_variable
      )
      
      # Capping K to be strictly less than num_controls AND at most max_acceptable_k
      limit_k <- min(max_acceptable_k, num_controls - 1)
      cancorplot$data <- cancorplot$data |> dplyr::filter(K <= limit_k)
      
      separation_score <- MultiScholaR:::calculateSeparationScore(cancorplot, separation_metric)
      
      # Find best k within the capped range
      controls_idx <- which(cancorplot$data$featureset == "Control")
      all_idx <- which(cancorplot$data$featureset == "All")
      difference_between_all_ctrl <- cancorplot$data$cc[all_idx] - cancorplot$data$cc[controls_idx]
      max_difference <- max(difference_between_all_ctrl, na.rm = TRUE)
      best_idx <- which(difference_between_all_ctrl == max_difference)[1]
      best_k <- (cancorplot$data$K[controls_idx])[best_idx]
      
      composite_score <- MultiScholaR:::calculateCompositeScore(
        separation_score, best_k, k_penalty_weight, max_acceptable_k
      )
      return(list(percentage = current_percentage, separation_score = separation_score,
                  best_k = best_k, composite_score = composite_score,
                  num_controls = num_controls, valid_plot = TRUE,
                  control_genes_index = control_genes_index, cancor_plot = cancorplot))
    }, error = function(e) {
      if (verbose) message(paste0("Error at ", current_percentage, "%: ", e$message))
      return(list(percentage = current_percentage, separation_score = NA_real_, best_k = NA_real_,
                  composite_score = NA_real_, num_controls = NA_integer_, valid_plot = FALSE,
                  control_genes_index = NULL, cancor_plot = NULL))
    })
  }
  
  all_results <- purrr::imap(percentage_range, process_percentage)
  results <- purrr::map_dfr(all_results, ~data.frame(percentage = .x$percentage,
      separation_score = .x$separation_score, best_k = .x$best_k,
      composite_score = .x$composite_score, num_controls = .x$num_controls,
      valid_plot = .x$valid_plot))
  valid_results <- all_results[!is.na(purrr::map_dbl(all_results, "composite_score"))]
  if (length(valid_results) == 0) {
    stop("No valid percentage found. Please check your data and parameters.")
  }
  best_index <- which.max(purrr::map_dbl(valid_results, "composite_score"))
  best_result <- valid_results[[best_index]]
  return(list(best_percentage = best_result$percentage, best_k = best_result$best_k,
      best_control_genes_index = best_result$control_genes_index,
      best_separation_score = best_result$separation_score, best_composite_score = best_result$best_composite_score,
      optimization_results = results, best_cancor_plot = best_result$best_cancor_plot))
}

