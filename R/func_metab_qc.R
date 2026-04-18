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
# Functions to extract here:
# - metaboliteIntensityFiltering(): S4 method for metabolite intensity filtering
# - metaboliteIntensityFilteringHelper(): Helper for intensity filtering
# - updateMetaboliteFiltering(): Updates filtering progress tracking
# - getFilteringProgressMetabolomics(): Gets filtering progress object
# - countUniqueMetabolites(): Counts unique metabolites
# - countMetabolitesPerSample(): Counts metabolites per sample
# - calculateMissingness(): Calculates missing value percentage
# - calculateSumIntensityPerSample(): Calculates sum intensity per sample
# - calculateMetaboliteCVs(): Calculates coefficient of variation
# - getInternalStandardMetrics(): Gets internal standard metrics
# - Additional metabolomics QC helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: metaboliteIntensityFiltering()
# Current location: R/metabolite_qc.R
# Type: S4 method (exportMethods)
# Description: Filters metabolites based on intensity thresholds
# setMethod(f = "metaboliteIntensityFiltering", ...) {
#   # Extract from R/metabolite_qc.R
# }

# Function 2: metaboliteIntensityFilteringHelper()
# Current location: R/metabolite_qc.R
# Description: Helper function for metabolite intensity filtering
# metaboliteIntensityFilteringHelper <- function(...) {
#   # Extract from R/metabolite_qc.R
# }

# Function 3: updateMetaboliteFiltering()
# Current location: R/QC_visualisation.R
# Description: Updates and visualizes metabolomics filtering progress
# updateMetaboliteFiltering <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 4: getFilteringProgressMetabolomics()
# Current location: R/QC_visualisation.R
# Description: Gets or initializes filtering progress object
# getFilteringProgressMetabolomics <- function() {
#   # Extract from R/QC_visualisation.R
# }

# Function 5: updateFilteringProgressMetabolomics()
# Current location: R/QC_visualisation.R
# Description: Updates the filtering progress object
# updateFilteringProgressMetabolomics <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 6: countUniqueMetabolites()
# Current location: R/QC_visualisation.R
# Description: Counts unique metabolites in an assay
# countUniqueMetabolites <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 7: countMetabolitesPerSample()
# Current location: R/QC_visualisation.R
# Description: Counts detected metabolites per sample
# countMetabolitesPerSample <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 8: calculateMissingness()
# Current location: R/QC_visualisation.R
# Description: Calculates overall missingness percentage
# calculateMissingness <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 9: calculateSumIntensityPerSample()
# Current location: R/QC_visualisation.R
# Description: Calculates sum intensity per sample (TIC proxy)
# calculateSumIntensityPerSample <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 10: calculateMetaboliteCVs()
# Current location: R/QC_visualisation.R
# Description: Calculates within-group metabolite CVs
# calculateMetaboliteCVs <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 11: calculateMetabolitePairCorrelation()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Description: Calculates correlation between metabolite pairs
# calculateMetabolitePairCorrelation <- function(...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 12: getInternalStandardMetrics()
# Current location: R/QC_visualisation.R
# Description: Calculates internal standard metrics
# getInternalStandardMetrics <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 13: calculateTotalUniqueMetabolitesAcrossAssays()
# Current location: R/QC_visualisation.R
# Description: Calculates total unique metabolites across assays
# calculateTotalUniqueMetabolitesAcrossAssays <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 14: findDuplicateFeatureIDs()
# Current location: R/metabolite_qc.R
# Description: Finds duplicate feature IDs
# findDuplicateFeatureIDs <- function(...) {
#   # Extract from R/metabolite_qc.R
# }

# Function 15: resolveDuplicateFeatures()
# Current location: R/metabolite_qc.R
# Type: S4 method (exportMethods)
# Description: Resolves duplicate features
# setMethod(f = "resolveDuplicateFeatures", ...) {
#   # Extract from R/metabolite_qc.R
# }

# Function 16: resolveDuplicateFeaturesByIntensity()
# Current location: R/metabolite_qc.R
# Description: Resolves duplicates by intensity
# resolveDuplicateFeaturesByIntensity <- function(...) {
#   # Extract from R/metabolite_qc.R
# }

# Function 17: generateMetaboliteFilteringPlots()
# Current location: R/QC_visualisation.R
# Description: Generates filtering progress plots
# generateMetaboliteFilteringPlots <- function(...) {
#   # Extract from R/QC_visualisation.R
# }






# ----------------------------------------------------------------------------
# updateMetaboliteFiltering
# ----------------------------------------------------------------------------
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
#' **Important:** Relies on and modifies the global `filtering_progress_metabolomics` object.
#' Requires helper functions (defined previously in this file) to be available.
#' For plot saving, requires either `time_dir` in the global environment or access to
#' the appropriate directory via `project_dirs$omics_type$time_dir`.
#'
#' @param theObject A S4 object (e.g., `MetaboliteAssayData`, `SummarizedExperiment`,
#'                  `MultiAssayExperiment`) containing metabolomics data. Must provide
#'                  access to a list of assays (data frames/tibbles with metabolite rows
#'                  and sample columns) and a colData/design matrix linking samples to groups.
#' @param step_name Character string uniquely identifying the current filtering step.
#' @param publication_graphs_dir Optional path for saving plots. If provided, the function
#'                              will try to find the corresponding time_dir.
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
                                      is_pattern = NULL) {
    prog_met <- getFilteringProgressMetabolomics()


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
            updateFilteringProgressMetabolomics(prog_met, step_name, assay_names, metrics_list_this_step, total_metabolites, overwrite)
        },
        error = function(e) {
            stop(e)
        }
    )

    plot_list <- tryCatch(
        {
            generateMetaboliteFilteringPlots(getFilteringProgressMetabolomics())
        },
        error = function(e) {
            stop(e)
        }
    )

    # --- 9. Directory handling and plot saving --- #
    actual_save_dir <- NULL # Will hold the final directory path for saving

    message("--- Plot Saving Diagnostics ---")
    message(sprintf("Value of publication_graphs_dir argument in function call: %s", ifelse(is.null(publication_graphs_dir), "NULL", publication_graphs_dir)))
    # The omics_type and time_dir arguments to this function are not directly used for path construction here;
    # we rely on global omic_type, experiment_label, and project_dirs.
    message(sprintf("Value of omics_type argument in function call: %s", ifelse(is.null(omics_type), "NULL", omics_type)))
    message(sprintf("Value of time_dir argument in function call: %s", ifelse(is.null(time_dir), "NULL", time_dir)))

    message("Attempting to determine save directory using global project_dirs, omic_type, and experiment_label...")

    if (exists("project_dirs", envir = .GlobalEnv) &&
        exists("omic_type", envir = .GlobalEnv) &&
        exists("experiment_label", envir = .GlobalEnv)) {
        message("Global variables project_dirs, omic_type, and experiment_label found.")
        local_project_dirs <- get("project_dirs", envir = .GlobalEnv)
        local_omic_type <- get("omic_type", envir = .GlobalEnv)
        local_experiment_label <- get("experiment_label", envir = .GlobalEnv)
        message(sprintf("Global omic_type value: '%s', Global experiment_label value: '%s'", local_omic_type, local_experiment_label))

        omics_key <- paste0(local_omic_type, "_", local_experiment_label)
        message(sprintf("Constructed omics_key for project_dirs: '%s'", omics_key))

        if (omics_key %in% names(local_project_dirs) &&
            !is.null(local_project_dirs[[omics_key]]) &&
            "time_dir" %in% names(local_project_dirs[[omics_key]])) {
            message(sprintf("omics_key '%s' found in project_dirs and has a 'time_dir' entry.", omics_key))
            retrieved_time_dir <- local_project_dirs[[omics_key]]$time_dir
            message(sprintf("Retrieved time_dir from project_dirs: %s", ifelse(is.null(retrieved_time_dir), "NULL", retrieved_time_dir)))

            if (is.null(retrieved_time_dir) || !is.character(retrieved_time_dir) || !nzchar(retrieved_time_dir)) {
                warning(sprintf("project_dirs[['%s']]$time_dir is NULL, not a character string, or empty. Plots will not be saved.", omics_key))
                message("Reason: retrieved_time_dir is invalid.")
            } else {
                actual_save_dir <- retrieved_time_dir
                message(sprintf("Successfully set actual_save_dir for plot saving to: '%s'", actual_save_dir))
            }
        } else {
            warning(sprintf("Could not find omics_key '%s' in project_dirs, or it lacks a 'time_dir' entry. Plots will not be saved.", omics_key))
            message(sprintf(
                "Details: omics_key '%s' in names(project_dirs): %s. project_dirs[['%s']] is NULL: %s. 'time_dir' in names(project_dirs[['%s']]): %s.",
                omics_key, omics_key %in% names(local_project_dirs),
                omics_key, is.null(local_project_dirs[[omics_key]]),
                omics_key, if (omics_key %in% names(local_project_dirs) && !is.null(local_project_dirs[[omics_key]])) "time_dir" %in% names(local_project_dirs[[omics_key]]) else NA
            ))
            message("Available keys in project_dirs: ", paste(names(local_project_dirs), collapse = ", "))
        }
    } else {
        warning("One or more global variables ('project_dirs', 'omic_type', 'experiment_label') not found. Plots will not be saved.")
        message(sprintf(
            "Exists project_dirs: %s, Exists omic_type: %s, Exists experiment_label: %s",
            exists("project_dirs", envir = .GlobalEnv),
            exists("omic_type", envir = .GlobalEnv),
            exists("experiment_label", envir = .GlobalEnv)
        ))
    }

    # Proceed with saving ONLY if actual_save_dir was successfully determined from global project_dirs
    if (!is.null(actual_save_dir)) {
        message(sprintf("Proceeding to save plots to derived directory: %s", actual_save_dir))
        if (!dir.exists(actual_save_dir)) {
            dir.create(actual_save_dir, recursive = TRUE)
            message("Created directory for QC plots: ", actual_save_dir)
        }

        # Save individual plots directly into actual_save_dir (which is the time_dir)
        purrr::iwalk(plot_list, function(plot, plot_name) {
            filename <- file.path(actual_save_dir, sprintf("%s_%s.png", step_name, plot_name))
            message(sprintf("Saving plot: %s", filename))
            ggsave(filename,
                plot = plot,
                width = 10,
                height = 8,
                dpi = 300
            )
        })

        # Save combined grid if return_grid is TRUE and plots exist
        if (return_grid && length(plot_list) > 0 && !is.null(plot_list[[1]]) && inherits(plot_list[[1]], "ggplot")) {
            # Use arrangeGrob (not grid.arrange) to create grob without drawing
            # Wrap in pdf(NULL)/dev.off() to prevent Rplots.pdf error
            pdf(NULL)
            grid_plot_obj <- do.call(gridExtra::arrangeGrob, c(plot_list, ncol = 2))
            invisible(dev.off())
            filename_grid <- file.path(actual_save_dir, sprintf("%s_combined_plots.png", step_name))
            message(sprintf("Saving combined grid plot: %s", filename_grid))
            ggsave(filename_grid, plot = grid_plot_obj, width = 15, height = 15, dpi = 300)
        }
        message("Metabolomics QC plots saved to: ", actual_save_dir)
    } else {
        # This block means actual_save_dir is still NULL.
        # This implies either globals were missing or project_dirs structure was invalid for deriving time_dir.
        message("No valid save directory determined from global project_dirs. Plots will not be saved.")
        # If publication_graphs_dir was provided in the call, and we still ended up here, it means the global lookup failed.
        if (!is.null(publication_graphs_dir)) {
            warning("Function was called with a publication_graphs_dir path, but plot saving still failed because a valid time_dir could not be derived from global project_dirs.")
        }
    }
    message("--- End of Plot Saving Diagnostics ---")


    if (length(plot_list) > 0) {
        if (!is.null(plot_list[[1]])) {}
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
























# ----------------------------------------------------------------------------
# metaboliteIntensityFiltering
# ----------------------------------------------------------------------------
#' @title Metabolite Intensity Filtering Method for MetaboliteAssayData
#'
#' @description
#' Filters metabolites in *all* assays of a MetaboliteAssayData object.
#' It removes metabolites that have intensities below a certain percentile threshold
#' in a proportion of samples exceeding a defined cutoff. The threshold is calculated
#' independently for each assay.
#'
#' @describeIn metaboliteIntensityFiltering Method for MetaboliteAssayData
#'
#' @param theObject A MetaboliteAssayData object.
#' @param metabolites_intensity_cutoff_percentile See generic definition.
#' @param metabolites_proportion_of_samples_below_cutoff See generic definition.
#'
#' @importFrom dplyr pull select all_of across
#' @importFrom rlang sym
#' @importFrom stats quantile
#'
#' @return An updated MetaboliteAssayData object.
#' @export
setMethod(
    f = "metaboliteIntensityFiltering",
    signature = "MetaboliteAssayData",
    definition = function(theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {
        # --- Parameter Resolution (Done once) ---
        config_intensity_percentile <- "metabolites_intensity_cutoff_percentile"
        raw_intensity_percentile <- checkParamsObjectFunctionSimplify(
            theObject,
            config_intensity_percentile,
            metabolites_intensity_cutoff_percentile
        )
        message("Raw intensity percentile from config/param: ", raw_intensity_percentile)
        cleaned_intensity_percentile <- trimws(sub("#.*$", "", raw_intensity_percentile))
        intensity_cutoff_percentile_final <- as.numeric(cleaned_intensity_percentile)

        config_proportion_cutoff <- "metabolites_proportion_of_samples_below_cutoff"
        raw_proportion_cutoff <- checkParamsObjectFunctionSimplify(
            theObject,
            config_proportion_cutoff,
            metabolites_proportion_of_samples_below_cutoff
        )
        message("Raw proportion cutoff from config/param: ", raw_proportion_cutoff)
        cleaned_proportion_cutoff <- trimws(sub("#.*$", "", raw_proportion_cutoff))
        proportion_of_samples_below_cutoff_final <- as.numeric(cleaned_proportion_cutoff)

        if (is.na(intensity_cutoff_percentile_final)) {
            stop("Failed to convert cleaned metabolites_intensity_cutoff_percentile ('", cleaned_intensity_percentile, "' from raw '", raw_intensity_percentile, "') to numeric. Check config.ini or parameter value.")
        }
        if (is.na(proportion_of_samples_below_cutoff_final)) {
            stop("Failed to convert cleaned metabolites_proportion_of_samples_below_cutoff ('", cleaned_proportion_cutoff, "' from raw '", raw_proportion_cutoff, "') to numeric. Check config.ini or parameter value.")
        }

        # --- Update Object Parameters (Done once) ---
        theObject <- updateParamInObject(theObject, config_intensity_percentile)
        theObject <- updateParamInObject(theObject, config_proportion_cutoff)

        # --- Process Each Assay in the List ---
        metabolite_id_col <- theObject@metabolite_id_column
        original_assay_list <- theObject@metabolite_data
        original_assay_names <- names(original_assay_list)

        if (length(original_assay_list) == 0) {
            warning("MetaboliteAssayData object has no assays in 'metabolite_data' slot. No filtering performed.")
            return(theObject)
        }

        # Iterate using indices
        filtered_assay_list <- lapply(seq_along(original_assay_list), function(i) {
            assay_table <- original_assay_list[[i]]
            # Determine assay name for messages (use index if no name)
            assay_name_for_msg <- if (!is.null(original_assay_names) && nzchar(original_assay_names[i])) {
                original_assay_names[i]
            } else {
                as.character(i) # Use index as fallback name
            }
            message("\nProcessing assay: ", assay_name_for_msg)

            if (!(metabolite_id_col %in% names(assay_table))) {
                warning("Metabolite ID column '", metabolite_id_col, "' not found in assay '", assay_name_for_msg, "'. Skipping this assay.")
                return(assay_table) # Return the original table if ID is missing
            }

            # Identify numeric sample columns for this assay
            sample_cols <- names(assay_table)[sapply(assay_table, is.numeric)]

            if (length(sample_cols) == 0) {
                warning("No numeric sample columns found in assay '", assay_name_for_msg, "'. Skipping filtering for this assay.")
                return(assay_table)
            }

            # Extract intensity values for this assay
            all_intensity_values <- assay_table |>
                dplyr::select(all_of(sample_cols)) |>
                unlist()

            if (length(all_intensity_values) == 0 || all(is.na(all_intensity_values))) {
                warning("No valid intensity values found in assay '", assay_name_for_msg, "' to calculate threshold. Skipping filtering for this assay.")
                return(assay_table)
            }

            # Calculate threshold specifically for this assay
            min_metabolite_intensity_threshold <- ceiling(quantile(all_intensity_values,
                na.rm = TRUE,
                probs = c(intensity_cutoff_percentile_final / 100)
            ))[1]

            message("Calculated minimum intensity threshold for assay '", assay_name_for_msg, "': ", min_metabolite_intensity_threshold)

            # Filter using Helper
            filtered_assay <- metaboliteIntensityFilteringHelper(
                assay_table = assay_table,
                min_metabolite_intensity_threshold = min_metabolite_intensity_threshold,
                metabolites_proportion_of_samples_below_cutoff = proportion_of_samples_below_cutoff_final,
                metabolite_id_column = metabolite_id_col
            )

            message("Filtered assay '", assay_name_for_msg, "'. Original rows: ", nrow(assay_table), ", Filtered rows: ", nrow(filtered_assay))
            return(filtered_assay)
        })

        # Restore original names if they existed
        if (!is.null(original_assay_names)) {
            names(filtered_assay_list) <- original_assay_names
        }

        # Assign the list of filtered assays back to the object
        theObject@metabolite_data <- filtered_assay_list

        # Optional: Call a generic cleanup/design matrix function if applicable
        # theObject <- cleanDesignMatrix(theObject) # If a generic method exists

        return(theObject)
    }
)
