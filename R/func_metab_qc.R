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
# func_metab_qc.R
# ============================================================================
# Purpose: Metabolomics quality control and filtering functions
#
# This file contains functions for metabolomics QC filtering, including
# intensity filtering, missing value analysis, CV calculations, and
# internal standard metrics. Functions in this file are used by metabolomics
# QC modules and related workflows.
#
# Ownership: metabolomics QC orchestration retained here. Correlation, filtering,
# metrics, plotting, progress, and variability live in func_metab_qc_* owners.
# ============================================================================

#' @title Update and Visualize Metabolomics Filtering Progress
#' @description Tracks and visualizes the impact of filtering steps on metabolomics
#'              data, updating a global `FilteringProgressMetabolomics` object.
#'              Generates QC plots summarizing the changes.
#'
#' @details
#' This function serves as the primary interface for tracking metabolomics QC.
#' It performs the following:
#' \itemize{
#'   \item Initializes or retrieves the global `FilteringProgressMetabolomics` object.
#'   \item Takes a `MetaboliteAssayData` object (or similar S4 containing assays list
#'         and design matrix) as input for the current processing step.
#'   \item Extracts the list of assay data frames/tibbles and the design matrix.
#'   \item For each assay, calls helper functions (`countUniqueMetabolites`,
#'         `countMetabolitesPerSample`, `calculateMissingness`,
#'         `calculateSumIntensityPerSample`, `calculateMetaboliteCVs`,
#'         `getInternalStandardMetrics`) to calculate QC metrics.
#'   \item Calls `calculateTotalUniqueMetabolitesAcrossAssays` for an overall count.
#'   \item Updates the global `FilteringProgressMetabolomics` object with the
#'         calculated metrics for the given `step_name`.
#'   \item Generates summary plots visualizing the tracked metrics across steps.
#'   \item Optionally saves plots to disk if `publication_graphs_dir` and a time
#'         directory are available.
#'   \item Returns either a combined grid plot or an invisible list of plots.
#' }
#'
#' By default, this preserves the legacy global filtering-progress behavior.
#' Artifact workflows supply `progress_owner` and explicit output paths so progress
#' and plot ownership remain within the workflow session.
#'
#' @param theObject A S4 object (e.g., `MetaboliteAssayData`, `SummarizedExperiment`,
#'                  `MultiAssayExperiment`) containing metabolomics data. Must provide
#'                  access to a list of assays (data frames/tibbles with metabolite rows
#'                  and sample columns) and a colData/design matrix linking samples to groups.
#' @param step_name Character string uniquely identifying the current filtering step.
#' @param publication_graphs_dir Optional path for saving plots.
#' @param omics_type Optional character string specifying the omics type (e.g., "metabolomics").
#'                  If provided and project_dirs exists in the global environment, will use
#'                  `project_dirs[[omics_type]]$time_dir` for plot saving.
#' @param time_dir Optional explicit path to the time directory. If provided, this overrides
#'                other methods of finding the time directory.
#' @param overwrite Logical, whether to overwrite existing data for `step_name`.
#' @param return_grid Logical, whether to return a `gridExtra` combined plot.
#' @param group_id_col Character, name of the column in `colData(theObject)` specifying groups.
#' @param sample_id_col Character, name of the sample ID column in `colData(theObject)`.
#' @param metabolite_id_col Character, name of the metabolite ID column in assay data.
#' @param is_pattern Character, regex for identifying internal standards. If not provided,
#'                   attempts to get from `theObject@internal_standard_regex` if slot exists.
#' @param progress_owner Optional workflow/session environment that owns progress state.
#'
#' @return If `return_grid` is `TRUE`, a `grob` object. Otherwise, an invisible list
#'         containing individual `ggplot` objects.
#'
#' @importFrom methods slotNames is
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom purrr map map_dfr imap imap_dfr walk iwalk
#' @export
updateMetaboliteFiltering <- function(theObject,
                                      step_name,
                                      publication_graphs_dir = NULL,
                                      omics_type = NULL,
                                      time_dir = NULL,
                                      overwrite = FALSE,
                                      return_grid = FALSE,
                                      group_id_col = NULL,
                                      sample_id_col = NULL,
                                      metabolite_id_col = NULL,
                                      is_pattern = NULL,
                                      progress_owner = NULL) {
    prog_met <- if (is.null(progress_owner)) {
        getFilteringProgressMetabolomics()
    } else {
        getFilteringProgressMetabolomics(progress_owner)
    }


    if (!isS4(theObject)) {
        stop("`theObject` must be an S4 object.")
    }

    # Check for specific class before checking generic slots
    if (inherits(theObject, "MetaboliteAssayData")) {
        # Specific handling for MetaboliteAssayData
        design_matrix <- theObject@design_matrix
        assay_list <- theObject@metabolite_data # Directly access the slot
        assay_names <- names(assay_list)
        if (is.null(assay_names)) assay_names <- paste0("Assay_", seq_along(assay_list))
        names(assay_list) <- assay_names
    } else {
        # Basic checks for required methods/slots for non-MetaboliteAssayData objects
        if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
            stop("Package 'SummarizedExperiment' needed for this function to work.")
        }
        if (!("assays" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
            stop("`theObject` must have an accessible `assays` method or slot.")
        }
        if (!("colData" %in% methods::slotNames(theObject) || canCoerce(theObject, "SummarizedExperiment"))) {
            stop("`theObject` must have an accessible `colData` method or slot.")
        }

        # Use SE accessors if available
        design_matrix <- SummarizedExperiment::colData(theObject)
        assay_list <- SummarizedExperiment::assays(theObject)
        # Ensure assay_list elements are data frames/tibbles if needed by helpers
        assay_list <- lapply(assay_list, as.data.frame)
        if (is.null(names(assay_list))) assay_names <- paste0("Assay_", seq_along(assay_list)) else assay_names <- names(assay_list)
        names(assay_list) <- assay_names
    }

    # Attempt to get parameters from object slots if not provided
    if (is.null(group_id_col) && "group_id" %in% slotNames(theObject)) group_id_col <- theObject@group_id
    if (is.null(sample_id_col) && "sample_id" %in% slotNames(theObject)) sample_id_col <- theObject@sample_id
    if (is.null(metabolite_id_col) && "metabolite_id_column" %in% slotNames(theObject)) metabolite_id_col <- theObject@metabolite_id_column
    if (is.null(is_pattern) && "internal_standard_regex" %in% slotNames(theObject)) is_pattern <- theObject@internal_standard_regex

    # Check if essential parameters are now available
    if (is.null(group_id_col)) stop("group_id_col must be provided or accessible via theObject@group_id")
    if (is.null(sample_id_col)) stop("sample_id_col must be provided or accessible via theObject@sample_id")
    if (is.null(metabolite_id_col)) stop("metabolite_id_col must be provided or accessible via theObject@metabolite_id_column")

    # Convert design matrix rownames to column if needed
    if (!sample_id_col %in% colnames(design_matrix)) {
        # If not, check if the rownames seem to match the sample IDs
        if (identical(rownames(design_matrix), as.character(design_matrix[[sample_id_col]]))) {
            # This case is unlikely if sample_id_col isn't a column name
            warning("Sample ID column '", sample_id_col, "' not found, but rownames might match? Check object structure.")
        } else if (!is.null(rownames(design_matrix)) && sample_id_col == "Run") { # Heuristic: If rownames exist and user expects 'Run'
            message("Moving rownames of design matrix to column: ", sample_id_col)
            design_matrix <- as.data.frame(design_matrix)
            design_matrix[[sample_id_col]] <- rownames(design_matrix)
            rownames(design_matrix) <- NULL # Remove rownames after moving
        } else {
            warning("Sample ID column '", sample_id_col, "' not found in design matrix and rownames don't seem to match or weren't checked.")
        }
    }

    design_matrix <- as.data.frame(design_matrix) # Ensure it's a data frame

    # Extract actual sample column names from design matrix
    # This is CRITICAL - sample_columns are the values in the sample_id_col (e.g., Run column)
    # which correspond to column names in the assay data
    sample_columns <- as.character(design_matrix[[sample_id_col]])

    metrics_list_this_step <- list()
    if (length(assay_list) > 0) {
        metrics_list_this_step <- purrr::imap(assay_list, function(current_assay_data, current_assay_name) {
            # Ensure current_assay_data is a data frame/tibble
            if (!is.data.frame(current_assay_data)) {
                warning("Assay ", current_assay_name, " is not a data frame. Skipping metrics calculation.")
                # Return placeholder for non-data frames
                return(list(
                    n_metabolites = 0,
                    detected_per_sample = data.frame(Run = character(), n_detected = integer()),
                    missingness = NA_real_,
                    sum_intensity_per_sample = data.frame(Run = character(), sum_intensity = numeric()),
                    cv_distribution = data.frame(metabolite_id = character(), group = character(), cv = numeric()),
                    is_metrics = data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
                ))
            }

            n_met <- tryCatch(
                {
                    countUniqueMetabolites(current_assay_data, metabolite_id_col)
                },
                error = function(e) {
                    stop(e)
                }
            )

            det_per_sample <- tryCatch(
                {
                    countMetabolitesPerSample(current_assay_data, sample_id_col, metabolite_id_col, sample_columns = sample_columns)
                },
                error = function(e) {
                    stop(e)
                }
            )

            miss <- tryCatch(
                {
                    calculateMissingness(current_assay_data, sample_id_col, sample_columns = sample_columns)
                },
                error = function(e) {
                    stop(e)
                }
            )

            sum_int <- tryCatch(
                {
                    calculateSumIntensityPerSample(current_assay_data, sample_id_col, sample_columns = sample_columns)
                },
                error = function(e) {
                    stop(e)
                }
            )

            cv_dist <- tryCatch(
                {
                    calculateMetaboliteCVs(current_assay_data, design_matrix, group_id_col, NULL, sample_id_col, metabolite_id_col, sample_columns = sample_columns)
                },
                error = function(e) {
                    stop(e)
                }
            )

            is_met <- tryCatch(
                {
                    getInternalStandardMetrics(current_assay_data, is_pattern, metabolite_id_col, sample_id_col, sample_columns = sample_columns)
                },
                error = function(e) {
                    stop(e)
                }
            )


            list(
                n_metabolites = n_met,
                detected_per_sample = det_per_sample,
                missingness = miss,
                sum_intensity_per_sample = sum_int,
                cv_distribution = cv_dist,
                is_metrics = is_met
            )
        })
    } else {
        warning("No assays found in theObject.")
        # Handle empty assay list case - maybe stop or proceed with empty metrics
        return(invisible(NULL)) # Or handle appropriately
    }


    total_metabolites <- tryCatch(
        {
            calculateTotalUniqueMetabolitesAcrossAssays(assay_list, metabolite_id_col)
        },
        error = function(e) {
            stop(e)
        }
    )

    tryCatch(
        {
            args <- list(
                prog_met,
                step_name,
                assay_names,
                metrics_list_this_step,
                total_metabolites,
                overwrite
            )
            if (!is.null(progress_owner)) args$owner <- progress_owner
            do.call(updateFilteringProgressMetabolomics, args)
        },
        error = function(e) {
            stop(e)
        }
    )

    plot_list <- tryCatch(
        {
            progress <- if (is.null(progress_owner)) {
                getFilteringProgressMetabolomics()
            } else {
                getFilteringProgressMetabolomics(progress_owner)
            }
            generateMetaboliteFilteringPlots(progress)
        },
        error = function(e) {
            stop(e)
        }
    )

    actual_save_dir <- resolveMetabQcPlotDirectory(
        publication_graphs_dir,
        time_dir
    )

    if (!is.null(actual_save_dir)) {
        if (!dir.exists(actual_save_dir)) {
            dir.create(actual_save_dir, recursive = TRUE)
        }

        # Save individual plots directly into actual_save_dir (which is the time_dir)
        purrr::iwalk(plot_list, function(plot, plot_name) {
            filename <- file.path(actual_save_dir, sprintf("%s_%s.png", step_name, plot_name))
            ggsave(filename,
                plot = plot,
                width = 10,
                height = 8,
                dpi = 300
            )
        })

        # Save combined grid if return_grid is TRUE and plots exist
        if (return_grid && length(plot_list) > 0L &&
            !is.null(plot_list[[1L]]) &&
            inherits(plot_list[[1L]], "ggplot")) {
            # Use arrangeGrob (not grid.arrange) to create grob without drawing
            # Wrap in pdf(NULL)/dev.off() to prevent Rplots.pdf error
            pdf(NULL)
            grid_plot_obj <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = 2))
            invisible(dev.off())
            filename_grid <- file.path(actual_save_dir, sprintf("%s_combined_plots.png", step_name))
            ggsave(filename_grid, plot = grid_plot_obj, width = 15, height = 15, dpi = 300)
        }
    } else {
        message("No valid metabolomics QC plot directory; plots were not saved.")
    }

    if (return_grid) {
        # Check conditions one by one
        cond1 <- length(plot_list) > 0
        cond2 <- !is.null(plot_list[[1]])
        cond3 <- if (cond2) inherits(plot_list[[1]], "ggplot") else FALSE


        if (cond1 && cond2 && cond3) {
            # Use arrangeGrob (not grid.arrange) to create grob without drawing
            # This allows the grob to be stored in a reactiveVal and rendered later by Shiny
            # CRITICAL: Wrap in pdf(NULL)/dev.off() to prevent "cannot open file 'Rplots.pdf'" error
            # in Shiny's sandboxed environment where arrangeGrob tries to create a temp PDF device
            grid_plot_obj <- tryCatch(
                {
                    pdf(NULL)
                    result <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = 2))
                    invisible(dev.off())
                    result
                },
                error = function(e) {
                    # Make sure to close device even on error
                    tryCatch(invisible(dev.off()), error = function(e2) NULL)
                    NULL
                }
            )


            if (is.null(grid_plot_obj)) {
                return(NULL)
            }

            return(grid_plot_obj)
        } else {
            return(NULL)
        }
    } else {
        # Print each plot individually
        if (length(plot_list) > 0) {
            message("Printing plots individually as return_grid is FALSE or grid could not be formed.")
            purrr::walk(plot_list, function(plot) {
                if (inherits(plot, "ggplot")) {
                    print(plot)
                } else {
                    message("Encountered a non-ggplot object in plot_list when trying to print individually.")
                }
            })
        } else {
            message("Plot list is empty, nothing to print.")
        }
        # Return the list invisibly
        invisible(plot_list)
    }
}
