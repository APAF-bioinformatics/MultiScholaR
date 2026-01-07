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

#' Build ITSD Selection Table for DT Display
#'
#' @title Build ITSD selection table for DT display
#' @description Creates a data frame suitable for DT::datatable display with
#'              pre-selection of candidate internal standards based on pattern matching.
#'
#' @param assay_data Data frame for a single assay
#' @param metabolite_id_col Column name for metabolite IDs
#' @param annotation_cols Character vector of columns to search for IS patterns
#' @param is_patterns Named list of regex patterns for IS detection
#'
#' @return Data frame with columns: feature_id, annotation, mean_intensity, cv_percent, is_candidate
#'
#' @importFrom dplyr select mutate across all_of any_of filter arrange desc
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect str_to_lower
#' @importFrom stats sd
#' @noRd
buildItsdSelectionTable <- function(
    assay_data
    , metabolite_id_col
    , annotation_cols = NULL
    , is_patterns = list(
        itsd = "(?i)itsd|istd|internal.?standard"
        , rt = "(?i)^rt_|_rt$|retention.?time"
        , isd = "(?i)^is_|_is$|^isd_|_isd$"
        , labeled = "(?i)_d[0-9]+$|_13c|_15n|deuterated|labeled"
    )
) {
    # --- Input validation ---
    stopifnot(
        is.data.frame(assay_data)
        , is.character(metabolite_id_col)
        , metabolite_id_col %in% colnames(assay_data)
    )

    # --- Identify sample columns (numeric) vs metadata columns ---
    numeric_cols <- names(assay_data)[sapply(assay_data, is.numeric)]
    sample_cols <- setdiff(numeric_cols, metabolite_id_col)

    if (length(sample_cols) == 0) {
        warning("No sample columns found in assay_data")
        return(data.frame(
            feature_id = character(0)
            , annotation = character(0)
            , mean_intensity = numeric(0)
            , cv_percent = numeric(0)
            , is_candidate = logical(0)
        ))
    }

    # --- Determine annotation columns to search ---
    if (is.null(annotation_cols)) {
        # Default: search metabolite ID column and common annotation columns
        potential_anno_cols <- c(
            metabolite_id_col
            , "metabolite"
            , "metabolite_identification"
            , "annotation"
            , "compound_name"
            , "name"
        )
        annotation_cols <- intersect(potential_anno_cols, colnames(assay_data))
    }

    if (length(annotation_cols) == 0) {
        annotation_cols <- metabolite_id_col
    }

    # --- Calculate summary statistics per feature ---
    feature_stats <- assay_data |>
        dplyr::mutate(
            feature_id = as.character(.data[[metabolite_id_col]])
            , mean_intensity = rowMeans(
                dplyr::across(dplyr::all_of(sample_cols))
                , na.rm = TRUE
            )
            , sd_intensity = apply(
                dplyr::across(dplyr::all_of(sample_cols))
                , 1
                , sd
                , na.rm = TRUE
            )
        ) |>
        dplyr::mutate(
            cv_percent = ifelse(
                mean_intensity > 0
                , (sd_intensity / mean_intensity) * 100
                , NA_real_
            )
        )

    # --- Create annotation string for pattern matching ---
    if (length(annotation_cols) > 1) {
        feature_stats$annotation <- apply(
            feature_stats[, annotation_cols, drop = FALSE]
            , 1
            , \(x) paste(na.omit(x), collapse = " | ")
        )
    } else {
        feature_stats$annotation <- as.character(feature_stats[[annotation_cols[1]]])
    }

    # --- Detect IS candidates using patterns ---
    combined_pattern <- paste(unlist(is_patterns), collapse = "|")

    feature_stats$is_candidate <- stringr::str_detect(
        stringr::str_to_lower(feature_stats$annotation)
        , combined_pattern
    )

    # --- Select and arrange output columns ---
    result <- feature_stats |>
        dplyr::select(
            feature_id
            , annotation
            , mean_intensity
            , cv_percent
            , is_candidate
        ) |>
        dplyr::arrange(dplyr::desc(is_candidate), cv_percent)

    return(result)
}


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
    ruv_mode <- match.arg(ruv_mode)

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

    # --- Validate grouping variable ---
    if (is.null(ruv_grouping_variable) || !ruv_grouping_variable %in% colnames(theObject@design_matrix)) {
        ruv_grouping_variable <- theObject@group_id
    }

    # --- Get assay names ---
    assay_list <- theObject@metabolite_data
    assay_names <- names(assay_list)

    logger::log_info(paste("Running RUV optimization for", length(assay_names), "assays in", ruv_mode, "mode"))

    # --- Run optimization per assay ---
    results <- purrr::imap(assay_list, \(assay_data, assay_name) {
        logger::log_info(paste("Processing assay:", assay_name))

        tryCatch({
            if (ruv_mode == "automatic") {
                # --- Automatic mode: test range of percentages ---
                percentage_range <- seq(percentage_min, percentage_max, by = 1)

                # Get negative control features
                ctrl_result <- getNegCtrlProtAnova(
                    theObject
                    , ruv_grouping_variable = ruv_grouping_variable
                    , percentage_as_neg_ctrl = percentage_max
                )

                # Extract control indices for this assay
                ctrl_indices <- if (is.list(ctrl_result) && assay_name %in% names(ctrl_result)) {
                    ctrl_result[[assay_name]]
                } else {
                    ctrl_result
                }

                # Generate canonical correlation plot
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

                # Find best k from cancor plot
                best_k <- findBestK(cancor_plot)

                # Calculate separation score
                separation_score <- calculateSeparationScore(
                    cancor_plot
                    , metric = separation_metric
                )

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
                # --- Manual mode: use specified parameters ---
                ctrl_result <- getNegCtrlProtAnova(
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
    return(results)
}


#' Extract Best K Values Per Assay from RUV Results
#'
#' @title Extract best k values per assay from RUV results
#' @param ruv_results List of per-assay RUV optimization results
#' @return Named list of k values (or NA for failed assays)
#' @noRd
extractBestKPerAssay <- function(ruv_results) {
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
extractCtrlPerAssay <- function(ruv_results) {
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
buildCombinedRuvTable <- function(ruv_results) {
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


#' Build Normalization Configuration Object
#'
#' @title Build normalization configuration object
#' @param input Shiny input object
#' @return List of configuration parameters for state saving
#' @noRd
buildNormConfig <- function(input) {
    list(
        itsd = list(
            applied = input$apply_itsd %||% FALSE
            , method = input$itsd_method %||% "median"
        )
        , log2 = list(
            offset = input$log_offset %||% 1
        )
        , normalization = list(
            method = input$norm_method %||% "cyclicloess"
        )
        , ruv = list(
            mode = input$ruv_mode %||% "skip"
            , grouping_variable = input$ruv_grouping_variable
            , auto_percentage_min = input$auto_percentage_min
            , auto_percentage_max = input$auto_percentage_max
            , separation_metric = input$separation_metric
            , k_penalty_weight = input$k_penalty_weight
            , adaptive_k_penalty = input$adaptive_k_penalty
            , manual_k = input$ruv_k
            , manual_percentage = input$ruv_percentage
        )
        , timestamp = Sys.time()
    )
}
