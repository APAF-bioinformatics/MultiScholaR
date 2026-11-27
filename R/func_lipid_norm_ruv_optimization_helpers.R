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
#' @importFrom purrr map map_dfr imap
#' @importFrom dplyr filter arrange desc
#' @importFrom logger log_info log_warn
#' @noRd
runLipidPerAssayRuvOptimization <- function(
    theObject
    , ruv_mode = c("automatic", "manual")
    , params = list()
    , experiment_paths = NULL
) {

    ruv_mode <- match.arg(ruv_mode)

    stopifnot(inherits(theObject, "LipidomicsAssayData"))

    # --- Extract parameters with defaults ---
    percentage_min        <- params$percentage_min        %||% 1
    percentage_max        <- params$percentage_max        %||% 20
    max_acceptable_k      <- params$max_acceptable_k      %||% 3L
    ruv_grouping_variable <- params$ruv_grouping_variable
    separation_metric     <- params$separation_metric     %||% "max_difference"
    k_penalty_weight      <- params$k_penalty_weight      %||% 0.5
    adaptive_k_penalty    <- params$adaptive_k_penalty    %||% TRUE
    manual_k              <- params$manual_k              %||% 2
    manual_percentage     <- params$manual_percentage     %||% 10

    # Emit deprecation warning once at public entry boundary
    if (ruv_mode == "automatic" && identical(separation_metric, "weighted_difference")) {
        warning(
            paste0(
                "'weighted_difference' is deprecated because it up-weights high-K values "
                , "while the composite score penalizes high K. "
                , "Use 'max_difference' or 'mean_difference' instead."
            )
            , call. = FALSE
        )
    }


    # --- Validate grouping variable ---
    if (is.null(ruv_grouping_variable) || !ruv_grouping_variable %in% colnames(theObject@design_matrix)) {
        ruv_grouping_variable <- theObject@group_id
    }

    assay_list  <- theObject@lipid_data
    assay_names <- names(assay_list)
    logger::log_info(paste("Running RUV optimization for", length(assay_names), "assays in", ruv_mode, "mode"))

    # =========================================================================
    # AUTOMATIC MODE
    # =========================================================================
    if (ruv_mode == "automatic") {

        percentage_range <- seq(percentage_min, percentage_max, by = 1)

        # --- Percentage-first whole-object loop ---
        # getNegCtrlMetabAnova() and ruvCancor() are called ONCE PER PERCENTAGE,
        # not once per assay.  Both return whole-object named lists.
        cancor_object <- tryCatch(
            prepareOmicsRuvCancorObject(theObject, "lipid_data", "lipid_id_column")
            , error = function(e) theObject
        )
        percentage_results <- purrr::map(percentage_range, function(pct) {

            ctrl_list <- tryCatch(
                getNegCtrlMetabAnova(
                    theObject
                    , ruv_grouping_variable = ruv_grouping_variable
                    , percentage_as_neg_ctrl = pct
                )
                , error = function(e) e
            )

            if (inherits(ctrl_list, "error")) {
                return(list(
                    percentage_requested = pct
                    , object_error = paste0("error_neg_ctrl_selection: ", ctrl_list$message)
                    , ctrl_list   = NULL
                    , cancor_list = NULL
                    , refusal     = if (inherits(ctrl_list, "multischolar_ruv_grouping_unusable")) ctrl_list
                ))
            }

            cancor_list <- tryCatch(
                ruvCancor(
                    cancor_object
                    , ctrl = ctrl_list
                    , ruv_grouping_variable = ruv_grouping_variable
                )
                , error = function(e) e
            )

            if (inherits(cancor_list, "error")) {
                return(list(
                    percentage_requested = pct
                    , object_error = paste0("error_cancor: ", cancor_list$message)
                    , ctrl_list   = ctrl_list
                    , cancor_list = NULL
                ))
            }

            list(
                percentage_requested = pct
                , object_error = NULL
                , ctrl_list   = ctrl_list
                , cancor_list = cancor_list
            )
        })

        # A grouping that cannot be used is refused to the caller, not recorded as failed rows (KIT-190). Signalled
        # here, outside purrr::map, so the refusal reaches the caller with its own class and message.
        for (pct_res in percentage_results) {
            if (!is.null(pct_res$refusal)) {
                stop(pct_res$refusal)
            }
        }

        # --- Per-assay derivation from shared percentage results ---
        design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])

        results <- purrr::imap(assay_list, function(assay_data, assay_name) {

            # Sample size: assay columns that appear in the design matrix
            sample_size <- length(intersect(colnames(assay_data), design_samples))
            candidate_feature_count <- nrow(assay_data)

            effective_max_k <- if (adaptive_k_penalty) {
                calculateAdaptiveMaxK(sample_size)
            } else {
                as.integer(max_acceptable_k)
            }

            # Build one result row per tested percentage for this assay
            complete_rows <- omicsRuvCompleteRows(assay_data, theObject@lipid_id_column, design_samples)

            assay_pct_rows <- purrr::map(percentage_results, function(pct_res) {
                pct <- pct_res$percentage_requested

                # Whole-object failure → failed row for this assay, run continues
                if (!is.null(pct_res$object_error)) {
                    err_status <- if (startsWith(pct_res$object_error, "error_neg_ctrl_selection")) {
                        "error_neg_ctrl_selection"
                    } else {
                        "error_cancor"
                    }
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = NA_integer_
                        , realized_percentage     = NA_real_
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = err_status
                        , error_reason            = pct_res$object_error
                        , ctrl_indices            = NULL
                        , cancor_plot             = NULL
                    ))
                }

                # KIT-242 (VAL-373 Q-A): under the current rule an assay the control selection skipped has no control set
                # of its own (the earlier rule took the whole list, on which sum() below fails); the reason travels on.
                if (omicsRuvFallbackCurrent() && is.list(pct_res$ctrl_list) && !is.null(names(pct_res$ctrl_list)) &&
                    !assay_name %in% names(pct_res$ctrl_list)) {
                    skip <- attr(pct_res$ctrl_list, "skipped")[[assay_name]] %||% list(reason = "not_selected")
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = 0L
                        , realized_percentage     = 0
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "skipped_no_controls_selected"
                        , error_reason            = skip$reason
                        , ctrl_indices            = NULL
                        , cancor_plot             = NULL
                        , skip_reason             = skip
                    ))
                }

                # Extract per-assay slices from whole-object results
                ctrl_indices <- if (is.list(pct_res$ctrl_list) && assay_name %in% names(pct_res$ctrl_list)) {
                    pct_res$ctrl_list[[assay_name]]
                } else {
                    pct_res$ctrl_list
                }

                cancor_plot <- if (is.list(pct_res$cancor_list) && assay_name %in% names(pct_res$cancor_list)) {
                    pct_res$cancor_list[[assay_name]]
                } else {
                    pct_res$cancor_list
                }

                realized_num_controls <- sum(ctrl_indices, na.rm = TRUE)
                realized_percentage <- if (candidate_feature_count > 0L) {
                    100 * realized_num_controls / candidate_feature_count
                } else {
                    NA_real_
                }

                if (realized_num_controls < 5L) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "skipped_insufficient_controls"
                        , error_reason            = ""
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = NULL
                    ))
                }

                # KIT-242 (VAL-242 Q6.2 (c)): under the current rule an assay the canonical-correlation step returned no
                # plot for (it could not be computed) is a failed computation, not a curve without an elbow.
                if (omicsRuvFallbackCurrent() && is.list(pct_res$cancor_list) && !assay_name %in% names(pct_res$cancor_list)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = NA_real_
                        , composite_score         = NA_real_
                        , status                  = "error_cancor"
                        , error_reason            = "error_cancor: the canonical correlations could not be computed for this assay"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = NULL
                    ))
                }

                separation_score <- tryCatch(
                    calculateSeparationScore(cancor_plot, metric = separation_metric)
                    , error = function(e) NA_real_
                )

                best_k <- tryCatch(
                    findBestKElbow(cancor_plot)
                    , error = function(e) NA_integer_
                )

                if (is.na(best_k)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = NA_integer_
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "invalid_cancor_plot"
                        , error_reason            = "findBestKElbow returned NA"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }

                best_k <- as.integer(best_k)
                complete_controls <- omicsRuvCompleteControlCount(complete_rows, ctrl_indices)
                # RUV-III cannot estimate more factors than it has controls observed in every sample (KIT-005).
                if (isTRUE(best_k > complete_controls)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = best_k
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "skipped_k_exceeds_complete_controls"
                        , error_reason            = sprintf("k = %d exceeds %d controls complete in every sample", best_k, complete_controls)
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }
                composite_score <- tryCatch(
                    calculateCompositeScore(separation_score, best_k, k_penalty_weight, effective_max_k)
                    , error = function(e) NA_real_
                )

                if (!is.finite(composite_score)) {
                    return(list(
                        percentage_requested    = pct
                        , candidate_feature_count = candidate_feature_count
                        , realized_num_controls   = realized_num_controls
                        , realized_percentage     = realized_percentage
                        , sample_size             = sample_size
                        , best_k                  = best_k
                        , separation_score        = separation_score
                        , composite_score         = NA_real_
                        , status                  = "error_scoring"
                        , error_reason            = "composite score is non-finite"
                        , ctrl_indices            = ctrl_indices
                        , cancor_plot             = cancor_plot
                    ))
                }

                list(
                    percentage_requested    = pct
                    , candidate_feature_count = candidate_feature_count
                    , realized_num_controls   = realized_num_controls
                    , realized_percentage     = realized_percentage
                    , sample_size             = sample_size
                    , best_k                  = best_k
                    , separation_score        = separation_score
                    , composite_score         = composite_score
                    , status                  = "ok"
                    , error_reason            = ""
                    , ctrl_indices            = ctrl_indices
                    , cancor_plot             = cancor_plot
                )
            })

            # Build optimization_results trace (one row per tested percentage)
            optimization_results <- purrr::map_dfr(assay_pct_rows, function(row) {
                data.frame(
                    percentage_requested    = row$percentage_requested
                    , candidate_feature_count = row$candidate_feature_count
                    , realized_num_controls   = row$realized_num_controls %||% NA_integer_
                    , realized_percentage     = row$realized_percentage   %||% NA_real_
                    , sample_size             = row$sample_size
                    , best_k                  = row$best_k                %||% NA_integer_
                    , separation_score        = row$separation_score      %||% NA_real_
                    , composite_score         = row$composite_score       %||% NA_real_
                    , status                  = row$status
                    , error_reason            = row$error_reason
                    , stringsAsFactors = FALSE
                )
            })

            # Deterministic winner selection: highest composite → lowest k →
            # highest separation → lowest percentage_requested
            valid_df <- optimization_results |>
                dplyr::filter(is.finite(composite_score)) |>
                dplyr::arrange(
                    dplyr::desc(composite_score)
                    , best_k
                    , dplyr::desc(separation_score)
                    , percentage_requested
                )

            if (nrow(valid_df) == 0L) {
                logger::log_warn(paste("RUV: no valid percentage for assay:", assay_name))
                refusal <- ruvNoValidPercentageMessage(optimization_results$percentage_requested
                    , optimization_results$status
                    , optimization_results$best_k
                    , vapply(assay_pct_rows, function(row) {
                        if (identical(row$status, "skipped_k_exceeds_complete_controls")) omicsRuvCompleteControlCount(complete_rows, row$ctrl_indices) else NA_integer_
                    }, integer(1L))
                    , "lipids")
                failed <- list(
                    success                      = FALSE
                    , best_k                     = NA_integer_
                    , best_percentage            = NA_real_
                    , best_realized_num_controls = NA_integer_
                    , best_realized_percentage   = NA_real_
                    , control_genes_index        = NULL
                    , cancor_plot                = NULL
                    , separation_score           = NA_real_
                    , composite_score            = NA_real_
                    , optimization_results       = optimization_results
                    , error                      = refusal %||% "No valid percentage found"
                )
                # KIT-242 (VAL-373 Q-A): why the control selection skipped the assay, when it did (else no field).
                skip <- Find(Negate(is.null), lapply(assay_pct_rows, `[[`, "skip_reason"))
                if (!is.null(skip)) failed$skip_reason <- skip
                return(failed)
            }

            best_pct <- valid_df$percentage_requested[[1L]]
            winner_idx <- which(
                vapply(assay_pct_rows, function(r) r$percentage_requested == best_pct, logical(1L))
            )[[1L]]
            winner <- assay_pct_rows[[winner_idx]]


            list(
                success                      = TRUE
                , best_k                     = as.integer(winner$best_k)
                , best_percentage            = winner$percentage_requested
                , best_realized_num_controls = as.integer(winner$realized_num_controls)
                , best_realized_percentage   = winner$realized_percentage
                , control_genes_index        = winner$ctrl_indices
                , cancor_plot                = winner$cancor_plot
                , separation_score           = winner$separation_score
                , composite_score            = winner$composite_score
                , optimization_results       = optimization_results
                , error                      = NULL
            )
        })

        names(results) <- assay_names
        return(results)
    }

    # =========================================================================
    # MANUAL MODE (unchanged semantics)
    # =========================================================================
    results <- purrr::imap(assay_list, function(assay_data, assay_name) {
        logger::log_info(paste("Processing assay:", assay_name))

        tryCatch({
            ctrl_result <- getNegCtrlMetabAnova(
                theObject
                , ruv_grouping_variable = ruv_grouping_variable
                , percentage_as_neg_ctrl = manual_percentage
            )

            # KIT-242 (VAL-373 Q-A): under the current rule an assay the control selection skipped fails the Manual
            # search with its reason (the earlier rule took the whole list as that assay's controls).
            if (omicsRuvFallbackCurrent() && is.list(ctrl_result) && !is.null(names(ctrl_result)) &&
                !assay_name %in% names(ctrl_result)) {
                skip <- attr(ctrl_result, "skipped")[[assay_name]] %||% list(reason = "not_selected")
                return(list(
                    success               = FALSE
                    , best_k              = NA_integer_
                    , best_percentage     = NA_real_
                    , control_genes_index = NULL
                    , cancor_plot         = NULL
                    , separation_score    = NA_real_
                    , optimization_results = NULL
                    , error               = paste("the control selection skipped the assay:", skip$reason)
                    , skip_reason         = skip
                ))
            }

            ctrl_indices <- if (is.list(ctrl_result) && assay_name %in% names(ctrl_result)) {
                ctrl_result[[assay_name]]
            } else {
                ctrl_result
            }

            cancor_result <- ruvCancor(
                theObject
                , ctrl = ctrl_indices
                , ruv_grouping_variable = ruv_grouping_variable
            )

            cancor_plot <- if (is.list(cancor_result) && assay_name %in% names(cancor_result)) {
                cancor_result[[assay_name]]
            } else {
                cancor_result
            }


            list(
                success               = TRUE
                , best_k              = manual_k
                , best_percentage     = manual_percentage
                , control_genes_index = ctrl_indices
                , cancor_plot         = cancor_plot
                , separation_score    = NA_real_
                , optimization_results = NULL
                , error               = NULL
            )

        }, error = function(e) {
            if (inherits(e, "multischolar_ruv_grouping_unusable")) {
                return(list(refusal = e))
            }
            logger::log_warn(paste("RUV manual mode failed for", assay_name, ":", e$message))
            list(
                success               = FALSE
                , best_k              = NA_integer_
                , best_percentage     = NA_real_
                , control_genes_index = NULL
                , cancor_plot         = NULL
                , separation_score    = NA_real_
                , optimization_results = NULL
                , error               = e$message
            )
        })
    })

    # A grouping that cannot be used is refused to the caller (KIT-190); signalled outside purrr::imap so the
    # refusal keeps its own class and message.
    for (res in results) {
        if (!is.null(res$refusal)) {
            stop(res$refusal)
        }
    }

    names(results) <- assay_names
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

