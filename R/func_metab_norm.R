# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization standalone functions (non-S4)
# 
# NOTE: S4 methods for metabolomics normalization (setMethod calls) remain in
# their source files (metabolite_normalization.R, metaboliteVsSamplesS4Objects.R)
# to maintain proper class-method coupling.
#
# This file is a placeholder for any future standalone (non-S4) metabolomics
# normalization helper functions that may be extracted.
#
# S4 methods for metabolomics normalization:
# - logTransformAssays() -> metabolite_normalization.R
# - normaliseUntransformedData() -> metabolite_normalization.R
# - normaliseBetweenSamples() -> metaboliteVsSamplesS4Objects.R
# - getNegCtrlMetabAnova() -> metaboliteVsSamplesS4Objects.R
#
# Dependencies:
# - limma
# - func_general_helpers.R (for utility functions)
# ============================================================================

# ============================================================================
# Helper Functions for Metabolomics Normalization Module
# ============================================================================



#' Generate QC Plots for Metabolomics Data and Save to Disk
#'
#' @title Generate QC plots for metabolomics data and save to disk
#' @description Generates PCA, RLE, Density, and Correlation plots for each assay
#'              in a MetaboliteAssayData object and saves them as PNG files.
#'
#' @param theObject MetaboliteAssayData S4 object
#' @param experiment_paths List of experiment directory paths (must contain metabolite_qc_dir)
#' @param stage One of "post_filter", "post_norm", "ruv_corrected"
#' @param grouping_variable Column name from design matrix for plot coloring
#' @param shape_variable Column name from design matrix for plot shapes (optional)
#' @param font_size Font size for plots
#'
#' @return Invisible list of file paths saved (named by assay and plot type)
#'
#' @importFrom ggplot2 ggsave
#' @importFrom purrr imap set_names
#' @importFrom logger log_info log_warn log_error
#' @noRd
generateMetabQcPlots <- function(
    theObject
    , experiment_paths
    , stage = c("post_filter", "post_norm", "ruv_corrected")
    , grouping_variable = NULL
    , shape_variable = NULL
    , font_size = 8
) {
    # --- Input validation ---
    stage <- match.arg(stage)

    stopifnot(
        inherits(theObject, "MetaboliteAssayData")
        , !is.null(experiment_paths$metabolite_qc_dir)
    )

    qc_dir <- experiment_paths$metabolite_qc_dir
    if (!dir.exists(qc_dir)) {
        dir.create(qc_dir, recursive = TRUE)
    }

    # --- Get assay names ---
    assay_list <- theObject@metabolite_data
    assay_names <- names(assay_list)

    if (length(assay_names) == 0) {
        logger::log_warn("No assays found in MetaboliteAssayData object")
        return(invisible(list()))
    }

    # --- Determine grouping variable ---
    if (is.null(grouping_variable)) {
        grouping_variable <- theObject@group_id
    }

    # --- Stage prefix mapping ---
    stage_prefix <- switch(
        stage
        , "post_filter" = "pre_norm"
        , "post_norm" = "post_norm"
        , "ruv_corrected" = "ruv_corrected"
    )

    # --- Generate plots for each assay ---
    saved_paths <- list()

    for (assay_name in assay_names) {
        safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
        logger::log_info(paste("Generating QC plots for assay:", assay_name))

        # --- PCA Plot ---
        tryCatch({
            pca_plots <- plotPca(
                theObject
                , grouping_variable = grouping_variable
                , label_column = ""
                , shape_variable = shape_variable
                , title = ""
                , font_size = font_size
            )

            # Select plot for this assay
            pca_plot <- if (is.list(pca_plots) && assay_name %in% names(pca_plots)) {
                pca_plots[[assay_name]]
            } else if (inherits(pca_plots, "ggplot")) {
                pca_plots
            } else if (is.list(pca_plots) && length(pca_plots) > 0) {
                pca_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(pca_plot)) {
                pca_path <- file.path(qc_dir, sprintf("%s_%s_pca.png", safe_name, stage_prefix))
                ggplot2::ggsave(pca_path, pca_plot, width = 8, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_pca")]] <- pca_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating PCA for", assay_name, ":", e$message))
        })

        # --- RLE Plot ---
        tryCatch({
            rle_plots <- plotRle(
                theObject
                , group = grouping_variable
                , yaxis_limit = c(-6, 6)
            )

            rle_plot <- if (is.list(rle_plots) && assay_name %in% names(rle_plots)) {
                rle_plots[[assay_name]]
            } else if (inherits(rle_plots, "ggplot")) {
                rle_plots
            } else if (is.list(rle_plots) && length(rle_plots) > 0) {
                rle_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(rle_plot)) {
                rle_path <- file.path(qc_dir, sprintf("%s_%s_rle.png", safe_name, stage_prefix))
                ggplot2::ggsave(rle_path, rle_plot, width = 10, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_rle")]] <- rle_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating RLE for", assay_name, ":", e$message))
        })

        # --- Density Plot ---
        tryCatch({
            # Generate PCA first for density plot
            pca_for_density <- plotPca(
                theObject
                , grouping_variable = grouping_variable
                , label_column = ""
                , shape_variable = shape_variable
                , title = ""
                , font_size = font_size
            )

            density_plots <- plotDensity(
                theObject
                , grouping_variable = grouping_variable
                , title = ""
                , font_size = font_size
            )

            density_plot <- if (is.list(density_plots) && assay_name %in% names(density_plots)) {
                density_plots[[assay_name]]
            } else if (inherits(density_plots, "ggplot")) {
                density_plots
            } else if (is.list(density_plots) && length(density_plots) > 0) {
                density_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(density_plot)) {
                density_path <- file.path(qc_dir, sprintf("%s_%s_density.png", safe_name, stage_prefix))
                ggplot2::ggsave(density_path, density_plot, width = 8, height = 6, dpi = 150)
                saved_paths[[paste0(assay_name, "_density")]] <- density_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating Density for", assay_name, ":", e$message))
        })

        # --- Pearson Correlation Plot ---
        tryCatch({
            corr_plots <- plotPearson(
                theObject
                , correlation_group = grouping_variable
            )

            corr_plot <- if (is.list(corr_plots) && assay_name %in% names(corr_plots)) {
                corr_plots[[assay_name]]
            } else if (inherits(corr_plots, "ggplot")) {
                corr_plots
            } else if (is.list(corr_plots) && length(corr_plots) > 0) {
                corr_plots[[1]]
            } else {
                NULL
            }

            if (!is.null(corr_plot)) {
                corr_path <- file.path(qc_dir, sprintf("%s_%s_correlation.png", safe_name, stage_prefix))
                ggplot2::ggsave(corr_path, corr_plot, width = 10, height = 8, dpi = 150)
                saved_paths[[paste0(assay_name, "_correlation")]] <- corr_path
            }
        }, error = function(e) {
            logger::log_warn(paste("Error generating Correlation for", assay_name, ":", e$message))
        })
    }

    logger::log_info(paste("Generated", length(saved_paths), "QC plots for stage:", stage))
    return(invisible(saved_paths))
}


#' Run Per-Assay RUV Optimization
#'
#' @title Run RUV optimization per assay
#' @description Automatic mode uses a percentage-first whole-object search loop:
#'   \code{getNegCtrlMetabAnova()} and \code{ruvCancor()} are called once per
#'   tested percentage (not once per assay), and per-assay result rows are then
#'   derived from the shared whole-object results.  Manual mode calls each helper
#'   once and uses the caller-supplied k and percentage directly.
#'
#' @param theObject MetaboliteAssayData S4 object
#' @param ruv_mode One of "automatic", "manual"
#' @param params List of RUV parameters:
#'   - percentage_min: Minimum percentage for auto mode (default 1)
#'   - percentage_max: Maximum percentage for auto mode (default 20)
#'   - max_acceptable_k: Hard cap on acceptable k (default 3, overridden by
#'     adaptive_k_penalty)
#'   - ruv_grouping_variable: Column name for grouping
#'   - separation_metric: One of "max_difference", "mean_difference", "auc",
#'     "weighted_difference" (deprecated)
#'   - k_penalty_weight: Penalty weight in [0,1] (default 0.5)
#'   - adaptive_k_penalty: Use sample-size-adaptive max_acceptable_k (default TRUE)
#'   - manual_k: k value for manual mode (default 2)
#'   - manual_percentage: Percentage for manual mode (default 10)
#' @param experiment_paths List of experiment directory paths (currently unused,
#'   kept for interface stability)
#'
#' @return Named list (one entry per assay).  Automatic-mode entries include:
#'   \itemize{
#'     \item success, best_k, best_percentage, best_realized_num_controls,
#'       best_realized_percentage, control_genes_index, cancor_plot,
#'       separation_score, composite_score, optimization_results, error
#'   }
#'   Manual-mode entries include: success, best_k, best_percentage,
#'   control_genes_index, cancor_plot, separation_score, optimization_results, error.
#'
#' @importFrom purrr map imap map_dfr
#' @importFrom dplyr filter arrange desc
#' @importFrom logger log_info log_warn
#' @noRd
runPerAssayRuvOptimization <- function(
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

    stopifnot(inherits(theObject, "MetaboliteAssayData"))

    # --- Extract parameters ---
    percentage_min       <- params$percentage_min       %||% 1
    percentage_max       <- params$percentage_max       %||% 20
    max_acceptable_k     <- params$max_acceptable_k     %||% 3L
    ruv_grouping_variable <- params$ruv_grouping_variable
    separation_metric    <- params$separation_metric    %||% "max_difference"
    k_penalty_weight     <- params$k_penalty_weight     %||% 0.5
    adaptive_k_penalty   <- params$adaptive_k_penalty   %||% TRUE
    manual_k             <- params$manual_k             %||% 2
    manual_percentage    <- params$manual_percentage    %||% 10

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

    message(sprintf(
        "   DEBUG66 [runPerAssayRuvOptimization] percentage_min=%s, percentage_max=%s"
        , percentage_min, percentage_max
    ))
    message(sprintf(
        "   DEBUG66 [runPerAssayRuvOptimization] ruv_grouping_variable='%s'"
        , ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable)
    ))

    # --- Validate grouping variable ---
    if (is.null(ruv_grouping_variable) || !ruv_grouping_variable %in% colnames(theObject@design_matrix)) {
        message(sprintf(
            "   DEBUG66 [runPerAssayRuvOptimization] fallback ruv_grouping_variable -> '%s'"
            , theObject@group_id
        ))
        ruv_grouping_variable <- theObject@group_id
    }

    assay_list  <- theObject@metabolite_data
    assay_names <- names(assay_list)
    message(sprintf(
        "   DEBUG66 [runPerAssayRuvOptimization] Found %d assays: %s"
        , length(assay_names), paste(assay_names, collapse = ", ")
    ))
    logger::log_info(paste("Running RUV optimization for", length(assay_names), "assays in", ruv_mode, "mode"))

    # =========================================================================
    # AUTOMATIC MODE
    # =========================================================================
    if (ruv_mode == "automatic") {

        percentage_range <- seq(percentage_min, percentage_max, by = 1)
        message(sprintf(
            "   DEBUG66 [runPerAssayRuvOptimization] Testing %d percentages: %s-%s"
            , length(percentage_range), percentage_min, percentage_max
        ))

        # --- Percentage-first whole-object loop ---
        # getNegCtrlMetabAnova() and ruvCancor() are called ONCE PER PERCENTAGE,
        # not once per assay.  Both return whole-object named lists.
        percentage_results <- purrr::map(percentage_range, function(pct) {
            message(sprintf(
                "   DEBUG66 [runPerAssayRuvOptimization] Whole-object pass: pct=%d%%", pct
            ))

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
                ))
            }

            cancor_list <- tryCatch(
                ruvCancor(
                    theObject
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

        # --- Per-assay derivation from shared percentage results ---
        design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])

        results <- purrr::imap(assay_list, function(assay_data, assay_name) {
            message(sprintf(
                "   DEBUG66 [runPerAssayRuvOptimization] Deriving per-assay results: '%s'"
                , assay_name
            ))

            # Sample size: assay columns that appear in the design matrix
            sample_size <- length(intersect(colnames(assay_data), design_samples))
            candidate_feature_count <- nrow(assay_data)

            effective_max_k <- if (adaptive_k_penalty) {
                calculateAdaptiveMaxK(sample_size)
            } else {
                as.integer(max_acceptable_k)
            }

            # Build one result row per tested percentage for this assay
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
                return(list(
                    success                  = FALSE
                    , best_k                 = NA_integer_
                    , best_percentage        = NA_real_
                    , best_realized_num_controls = NA_integer_
                    , best_realized_percentage   = NA_real_
                    , control_genes_index    = NULL
                    , cancor_plot            = NULL
                    , separation_score       = NA_real_
                    , composite_score        = NA_real_
                    , optimization_results   = optimization_results
                    , error                  = "No valid percentage found"
                ))
            }

            best_pct <- valid_df$percentage_requested[[1L]]
            winner_idx <- which(
                vapply(assay_pct_rows, function(r) r$percentage_requested == best_pct, logical(1L))
            )[[1L]]
            winner <- assay_pct_rows[[winner_idx]]

            message(sprintf(
                "   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': winner pct=%d%%, k=%d, composite=%.4f"
                , assay_name
                , winner$percentage_requested
                , winner$best_k
                , winner$composite_score
            ))

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
        message(sprintf(
            "   DEBUG66 [runPerAssayRuvOptimization] Automatic mode complete. %d assay results."
            , length(results)
        ))
        message("+===========================================================================+")
        message("|  DEBUG66: Exiting runPerAssayRuvOptimization                              |")
        message("+===========================================================================+")
        return(results)
    }

    # =========================================================================
    # MANUAL MODE (unchanged semantics)
    # =========================================================================
    results <- purrr::imap(assay_list, function(assay_data, assay_name) {
        message(sprintf(
            "   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': MANUAL mode", assay_name
        ))
        logger::log_info(paste("Processing assay:", assay_name))

        tryCatch({
            ctrl_result <- getNegCtrlMetabAnova(
                theObject
                , ruv_grouping_variable = ruv_grouping_variable
                , percentage_as_neg_ctrl = manual_percentage
            )

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

            message(sprintf(
                "   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': MANUAL mode SUCCESS - k=%s"
                , assay_name, manual_k
            ))

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
            message(sprintf(
                "   DEBUG66 [runPerAssayRuvOptimization] Assay '%s': MANUAL ERROR - %s"
                , assay_name, e$message
            ))
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

    names(results) <- assay_names
    message(sprintf(
        "   DEBUG66 [runPerAssayRuvOptimization] Manual mode complete. %d assay results."
        , length(results)
    ))
    message("+===========================================================================+")
    message("|  DEBUG66: Exiting runPerAssayRuvOptimization                              |")
    message("+===========================================================================+")
    return(results)
}








