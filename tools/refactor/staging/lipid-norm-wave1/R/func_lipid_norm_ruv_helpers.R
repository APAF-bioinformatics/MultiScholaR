#' Run Per-Assay RUV Optimization
#'
#' @title Run RUV optimization per assay
#' @description Runs RUV-III optimization independently for each assay in a
#'              LipidomicsAssayData object, using the same logic as proteomics.
#'
#' @param theObject LipidomicsAssayData S4 object
#' @param ruv_mode One of "automatic", "manual"
#' @param params List of RUV parameters:
#'   - percentage_min: Minimum percentage for auto mode
#'   - percentage_max: Maximum percentage for auto mode
#'   - ruv_grouping_variable: Column name for grouping
#'   - separation_metric: Metric for optimization
#'   - k_penalty_weight: Penalty weight for k values
#'   - adaptive_k_penalty: Whether to use adaptive penalty
#'   - manual_k: K value for manual mode
#'   - manual_percentage: Percentage for manual mode
#' @param experiment_paths List of experiment directory paths
#'
#' @return Named list of RUV results per assay, each containing:
#'   - success: Logical indicating success
#'   - best_k: Optimal k value
#'   - best_percentage: Optimal percentage
#'   - control_genes_index: Boolean vector of control features
#'   - cancor_plot: ggplot object for canonical correlation
#'   - optimization_results: Data frame of all tested parameters (auto mode only)
#'   - error: Error message if failed
#'
#' @importFrom purrr imap
#' @importFrom logger log_info log_warn log_error
#' @noRd
runLipidPerAssayRuvOptimization <- function(
    theObject
    , ruv_mode = c("automatic", "manual")
    , params = list()
    , experiment_paths = NULL
) {
    message("+===========================================================================+")
    message("|  DEBUG66: Entering runPerAssayRuvOptimization                             |")
    message("+===========================================================================+")
    
    ruv_mode <- match.arg(ruv_mode)
    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] ruv_mode = '%s'", ruv_mode))

    stopifnot(inherits(theObject, "LipidomicsAssayData"))

    # --- Extract parameters with defaults ---
    percentage_min <- params$percentage_min %||% 1
    percentage_max <- params$percentage_max %||% 20
    ruv_grouping_variable <- params$ruv_grouping_variable
    separation_metric <- params$separation_metric %||% "max_difference"
    k_penalty_weight <- params$k_penalty_weight %||% 0.5
    adaptive_k_penalty <- params$adaptive_k_penalty %||% TRUE
    manual_k <- params$manual_k %||% 2
    manual_percentage <- params$manual_percentage %||% 10

    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] percentage_min = %s, percentage_max = %s", percentage_min, percentage_max))
    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] ruv_grouping_variable = '%s'", ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable)))
    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] manual_k = %s, manual_percentage = %s", manual_k, manual_percentage))

    # --- Validate grouping variable ---
    if (is.null(ruv_grouping_variable) || !ruv_grouping_variable %in% colnames(theObject@design_matrix)) {
        message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] ruv_grouping_variable fallback to group_id: '%s'", theObject@group_id))
        ruv_grouping_variable <- theObject@group_id
    }

    # --- Get assay names ---
    assay_list <- theObject@lipid_data
    assay_names <- names(assay_list)
    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Found %d assays: %s", length(assay_names), paste(assay_names, collapse = ", ")))

    logger::log_info(paste("Running RUV optimization for", length(assay_names), "assays in", ruv_mode, "mode"))

    # --- Run optimization per assay ---
    results <- purrr::imap(assay_list, \(assay_data, assay_name) {
        message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] === Processing assay: %s ===", assay_name))
        logger::log_info(paste("Processing assay:", assay_name))

        tryCatch({
            if (ruv_mode == "automatic") {
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': AUTOMATIC mode", assay_name))
                # --- Automatic mode: test range of percentages ---
                percentage_range <- seq(percentage_min, percentage_max, by = 1)

                # Get negative control features
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Calling getNegCtrlMetabAnova with percentage = %s", assay_name, percentage_max))
                ctrl_result <- getNegCtrlMetabAnova(
                    theObject
                    , ruv_grouping_variable = ruv_grouping_variable
                    , percentage_as_neg_ctrl = percentage_max
                )
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': getNegCtrlMetabAnova returned type = '%s', is.null = %s", 
                               assay_name, class(ctrl_result)[1], is.null(ctrl_result)))

                # Extract control indices for this assay
                ctrl_indices <- if (is.list(ctrl_result) && assay_name %in% names(ctrl_result)) {
                    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Extracting ctrl_indices from named list", assay_name))
                    ctrl_result[[assay_name]]
                } else {
                    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Using ctrl_result directly", assay_name))
                    ctrl_result
                }
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': ctrl_indices is.null = %s, length = %s, sum(TRUE) = %s", 
                               assay_name, is.null(ctrl_indices), length(ctrl_indices), sum(ctrl_indices, na.rm = TRUE)))

                # Generate canonical correlation plot
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Calling ruvCancor", assay_name))
                cancor_result <- ruvCancor(
                    theObject
                    , ctrl = ctrl_indices
                    , ruv_grouping_variable = ruv_grouping_variable
                )
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': ruvCancor returned type = '%s'", assay_name, class(cancor_result)[1]))

                cancor_plot <- if (is.list(cancor_result) && assay_name %in% names(cancor_result)) {
                    cancor_result[[assay_name]]
                } else {
                    cancor_result
                }
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': cancor_plot is.null = %s, class = '%s'", 
                               assay_name, is.null(cancor_plot), class(cancor_plot)[1]))

                # Find best k from cancor plot
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Calling findBestK", assay_name))
                best_k <- findBestK(cancor_plot)
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': findBestK returned = %s", assay_name, best_k))

                # Calculate separation score
                separation_score <- calculateSeparationScore(
                    cancor_plot
                    , metric = separation_metric
                )
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': separation_score = %s", assay_name, separation_score))
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': AUTOMATIC mode SUCCESS - best_k = %s, num_controls = %s", 
                               assay_name, best_k, sum(ctrl_indices, na.rm = TRUE)))

                list(
                    success = TRUE
                    , best_k = best_k
                    , best_percentage = percentage_max
                    , control_genes_index = ctrl_indices
                    , cancor_plot = cancor_plot
                    , separation_score = separation_score
                    , optimization_results = data.frame(
                        percentage = percentage_max
                        , best_k = best_k
                        , separation_score = separation_score
                        , num_controls = sum(ctrl_indices, na.rm = TRUE)
                    )
                    , error = NULL
                )

            } else {
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': MANUAL mode", assay_name))
                # --- Manual mode: use specified parameters ---
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Calling getNegCtrlMetabAnova with percentage = %s", assay_name, manual_percentage))
                ctrl_result <- getNegCtrlMetabAnova(
                    theObject
                    , ruv_grouping_variable = ruv_grouping_variable
                    , percentage_as_neg_ctrl = manual_percentage
                )
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': getNegCtrlMetabAnova returned type = '%s', is.null = %s", 
                               assay_name, class(ctrl_result)[1], is.null(ctrl_result)))

                ctrl_indices <- if (is.list(ctrl_result) && assay_name %in% names(ctrl_result)) {
                    ctrl_result[[assay_name]]
                } else {
                    ctrl_result
                }
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': ctrl_indices is.null = %s, length = %s, sum(TRUE) = %s", 
                               assay_name, is.null(ctrl_indices), length(ctrl_indices), sum(ctrl_indices, na.rm = TRUE)))

                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': Calling ruvCancor", assay_name))
                cancor_result <- ruvCancor(
                    theObject
                    , ctrl = ctrl_indices
                    , ruv_grouping_variable = ruv_grouping_variable
                )
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': ruvCancor returned type = '%s'", assay_name, class(cancor_result)[1]))

                cancor_plot <- if (is.list(cancor_result) && assay_name %in% names(cancor_result)) {
                    cancor_result[[assay_name]]
                } else {
                    cancor_result
                }
                message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': MANUAL mode SUCCESS - best_k = %s, num_controls = %s", 
                               assay_name, manual_k, sum(ctrl_indices, na.rm = TRUE)))

                list(
                    success = TRUE
                    , best_k = manual_k
                    , best_percentage = manual_percentage
                    , control_genes_index = ctrl_indices
                    , cancor_plot = cancor_plot
                    , separation_score = NA_real_
                    , optimization_results = NULL
                    , error = NULL
                )
            }

        }, error = function(e) {
            message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': ERROR CAUGHT - %s", assay_name, e$message))
            logger::log_warn(paste("RUV optimization failed for", assay_name, ":", e$message))
            list(
                success = FALSE
                , best_k = NA_integer_
                , best_percentage = NA_real_
                , control_genes_index = NULL
                , cancor_plot = NULL
                , separation_score = NA_real_
                , optimization_results = NULL
                , error = e$message
            )
        })
    })

    names(results) <- assay_names
    message(sprintf("   DEBUG66 [runPerAssayRuvOptimization] Completed all assays. Returning %d results.", length(results)))
    message("+===========================================================================+")
    message("|  DEBUG66: Exiting runPerAssayRuvOptimization                              |")
    message("+===========================================================================+")
    return(results)
}

#' Extract Best K Values Per Assay from RUV Results
#'
#' @title Extract best k values per assay from RUV results
#' @param ruv_results List of per-assay RUV optimization results
#' @return Named list of k values (or NA for failed assays)
#' @noRd
extractLipidBestKPerAssay <- function(ruv_results) {
    purrr::map(ruv_results, \(x) {
        if (isTRUE(x$success)) x$best_k else NA_integer_
    })
}

#' Extract Control Feature Indices Per Assay from RUV Results
#'
#' @title Extract control feature indices per assay from RUV results
#' @param ruv_results List of per-assay RUV optimization results
#' @return Named list of control indices (or NULL for failed assays)
#' @noRd
extractLipidCtrlPerAssay <- function(ruv_results) {
    purrr::map(ruv_results, \(x) {
        if (isTRUE(x$success)) x$control_genes_index else NULL
    })
}

#' Build Combined RUV Optimization Results Table
#'
#' @title Build combined RUV optimization results table
#' @param ruv_results List of per-assay RUV optimization results
#' @return Data frame with Assay column and optimization metrics
#' @importFrom dplyr bind_rows mutate
#' @noRd
buildLipidCombinedRuvTable <- function(ruv_results) {
    rows <- purrr::imap(ruv_results, \(result, assay_name) {
        if (isTRUE(result$success)) {
            data.frame(
                Assay = assay_name
                , Best_K = result$best_k
                , Best_Percentage = result$best_percentage
                , Separation_Score = round(result$separation_score, 4)
                , Num_Controls = sum(result$control_genes_index, na.rm = TRUE)
                , Status = "Success"
                , stringsAsFactors = FALSE
            )
        } else {
            data.frame(
                Assay = assay_name
                , Best_K = NA_integer_
                , Best_Percentage = NA_real_
                , Separation_Score = NA_real_
                , Num_Controls = NA_integer_
                , Status = paste("Failed:", result$error)
                , stringsAsFactors = FALSE
            )
        }
    })

    dplyr::bind_rows(rows)
}

