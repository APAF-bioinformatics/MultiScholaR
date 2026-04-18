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
#' @description Runs RUV-III optimization independently for each assay in a
#'              MetaboliteAssayData object, using the same logic as proteomics.
#'
#' @param theObject MetaboliteAssayData S4 object
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
    assay_list <- theObject@metabolite_data
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








