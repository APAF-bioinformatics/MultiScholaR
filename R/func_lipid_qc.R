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
# func_lipid_qc.R
# ============================================================================
# Purpose: Lipidomics quality control and filtering functions
#
# This file contains functions for lipidomics QC filtering, including
# intensity filtering, missing value analysis, CV calculations, and
# internal standard metrics. Functions in this file are used by lipidomics
# QC modules and related workflows.
#
# Functions to extract here:
# - lipidIntensityFiltering(): S4 method for lipid intensity filtering
# - lipidIntensityFilteringHelper(): Helper for intensity filtering
# - updateLipidFiltering(): Updates filtering progress tracking
# - getFilteringProgressLipidomics(): Gets filtering progress object
# - countUniqueLipids(): Counts unique lipids
# - countLipidsPerSample(): Counts lipids per sample
# - calculateLipidMissingness(): Calculates missing value percentage
# - calculateLipidSumIntensityPerSample(): Calculates sum intensity per sample
# - calculateLipidCVs(): Calculates coefficient of variation
# - getLipidInternalStandardMetrics(): Gets internal standard metrics
# - Additional lipidomics QC helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: lipidIntensityFiltering()
# Current location: R/lipid_qc.R
# Type: S4 method (exportMethods)
# Description: Filters lipids based on intensity thresholds
# setMethod(f = "lipidIntensityFiltering", ...) {
#   # Extract from R/lipid_qc.R
# }

# Function 2: lipidIntensityFilteringHelper()
# Current location: R/lipid_qc.R
# Description: Helper function for lipid intensity filtering
# lipidIntensityFilteringHelper <- function(...) {
#   # Extract from R/lipid_qc.R
# }

# Function 3: updateLipidFiltering()
# Current location: R/QC_visualisation.R
# Description: Updates and visualizes lipidomics filtering progress
# updateLipidFiltering <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 4: getFilteringProgressLipidomics()
# Current location: R/QC_visualisation.R
# Description: Gets or initializes filtering progress object
# getFilteringProgressLipidomics <- function() {
#   # Extract from R/QC_visualisation.R
# }

# Function 5: updateFilteringProgressLipidomics()
# Current location: R/QC_visualisation.R
# Description: Updates the filtering progress object
# updateFilteringProgressLipidomics <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 6: countUniqueLipids()
# Current location: R/QC_visualisation.R
# Description: Counts unique lipids in an assay
# countUniqueLipids <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 7: countLipidsPerSample()
# Current location: R/QC_visualisation.R
# Description: Counts detected lipids per sample
# countLipidsPerSample <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 8: calculateLipidMissingness()
# Current location: R/QC_visualisation.R
# Description: Calculates overall missingness percentage
# calculateLipidMissingness <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 9: calculateLipidSumIntensityPerSample()
# Current location: R/QC_visualisation.R
# Description: Calculates sum intensity per sample (TIC proxy)
# calculateLipidSumIntensityPerSample <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 10: calculateLipidCVs()
# Current location: R/QC_visualisation.R
# Description: Calculates within-group lipid CVs
# calculateLipidCVs <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 11: calculateLipidPairCorrelation()
# Current location: R/lipidVsSamplesS4Objects.R
# Description: Calculates correlation between lipid pairs
# calculateLipidPairCorrelation <- function(...) {
#   # Extract from R/lipidVsSamplesS4Objects.R
# }

# Function 12: getLipidInternalStandardMetrics()
# Current location: R/QC_visualisation.R
# Description: Calculates internal standard metrics
# getLipidInternalStandardMetrics <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 13: calculateTotalUniqueLipidsAcrossAssays()
# Current location: R/QC_visualisation.R
# Description: Calculates total unique lipids across assays
# calculateTotalUniqueLipidsAcrossAssays <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 14: findDuplicateFeatureIDs()
# Current location: R/lipid_qc.R
# Description: Finds duplicate feature IDs
# findDuplicateFeatureIDs <- function(...) {
#   # Extract from R/lipid_qc.R
# }

# Function 15: resolveDuplicateFeatures()
# Current location: R/lipid_qc.R
# Type: S4 method (exportMethods)
# Description: Resolves duplicate features
# setMethod(f = "resolveDuplicateFeatures", ...) {
#   # Extract from R/lipid_qc.R
# }

# Function 16: resolveLipidDuplicateFeaturesByIntensity()
# Current location: R/lipid_qc.R
# Description: Resolves duplicates by intensity
# resolveLipidDuplicateFeaturesByIntensity <- function(...) {
#   # Extract from R/lipid_qc.R
# }

# Function 17: generateLipidFilteringPlots()
# Current location: R/QC_visualisation.R
# Description: Generates filtering progress plots
# generateLipidFilteringPlots <- function(...) {
#   # Extract from R/QC_visualisation.R
# }


# ----------------------------------------------------------------------------
# lipidIntensityFilteringHelper
# ----------------------------------------------------------------------------
#' @title Helper function for lipid intensity filtering
#' @name lipidIntensityFilteringHelper
#' @description Filter lipids based on an intensity threshold and the proportion of samples below that threshold in a wide-format table.
#' @param assay_table A wide data frame where rows are lipids and columns include a lipid identifier and numeric sample intensities.
#' @param min_lipid_intensity_threshold The calculated minimum intensity value. Lipids in samples below this threshold are considered 'below threshold'.
#' @param lipids_proportion_of_samples_below_cutoff The maximum allowed proportion (0 to 1) of samples where a lipid can be below the threshold. If a lipid exceeds this proportion, it's removed.
#' @param lipid_id_column A string specifying the name of the column containing the unique lipid identifiers.
#' @return A filtered wide data frame containing only the lipids that pass the filter.
#' @export
lipidIntensityFilteringHelper <- function(
  assay_table,
  min_lipid_intensity_threshold,
  lipids_proportion_of_samples_below_cutoff,
  lipid_id_column
) {
    # Identify numeric columns representing sample intensities
    sample_cols <- names(assay_table)[sapply(assay_table, is.numeric)]
    num_samples <- length(sample_cols)

    if (num_samples == 0) {
        warning("No numeric sample columns found in the assay table. Returning original table.")
        return(assay_table)
    }

    # Calculate the number of samples below threshold for each lipid
    lipids_below_threshold <- assay_table |>
        # Ensure id column is character for safe rowwise operations if needed
        # mutate({{lipid_id_column}} := as.character({{lipid_id_column}})) |>
        rowwise() |>
        mutate(
            num_below_threshold = sum(c_across(all_of(sample_cols)) < min_lipid_intensity_threshold, na.rm = TRUE),
            proportion_below_threshold = num_below_threshold / num_samples
        ) |>
        ungroup()

    # Filter lipids based on the proportion cutoff
    filtered_assay_table <- lipids_below_threshold |>
        dplyr::filter(proportion_below_threshold < lipids_proportion_of_samples_below_cutoff) |>
        # Remove the temporary calculation columns
        dplyr::select(-num_below_threshold, -proportion_below_threshold)

    return(filtered_assay_table)
}
















#' @importFrom methods slotNames is
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggsave
#' @importFrom purrr map map_dfr imap imap_dfr walk iwalk
#' @export
updateLipidFiltering <- function(theObject,
                                 step_name,
                                 publication_graphs_dir = NULL,
                                 omics_type = NULL,
                                 time_dir = NULL,
                                 overwrite = FALSE,
                                 return_grid = FALSE,
                                 group_id_col = NULL,
                                 sample_id_col = NULL,
                                 lipid_id_col = NULL,
                                 is_pattern = NULL) {
    prog_met <- getFilteringProgressLipidomics()
    filtering_context <- prepareLipidFilteringContext(
        theObject = theObject,
        group_id_col = group_id_col,
        sample_id_col = sample_id_col,
        lipid_id_col = lipid_id_col,
        is_pattern = is_pattern
    )
    design_matrix <- filtering_context$design_matrix
    assay_list <- filtering_context$assay_list
    assay_names <- filtering_context$assay_names
    group_id_col <- filtering_context$group_id_col
    sample_id_col <- filtering_context$sample_id_col
    lipid_id_col <- filtering_context$lipid_id_col
    is_pattern <- filtering_context$is_pattern
    sample_columns <- filtering_context$sample_columns

    metrics_list_this_step <- list()
    if (length(assay_list) > 0) {
        metrics_list_this_step <- purrr::imap(assay_list, function(current_assay_data, current_assay_name) {
            calculateLipidFilteringAssayMetrics(
                current_assay_data = current_assay_data,
                current_assay_name = current_assay_name,
                design_matrix = design_matrix,
                group_id_col = group_id_col,
                sample_id_col = sample_id_col,
                lipid_id_col = lipid_id_col,
                is_pattern = is_pattern,
                sample_columns = sample_columns
            )
        })
    } else {
        warning("No assays found in theObject.")
        # Handle empty assay list case - maybe stop or proceed with empty metrics
        return(invisible(NULL)) # Or handle appropriately
    }

    filtering_step_outputs <- finalizeLipidFilteringStep(
        prog_met = prog_met,
        step_name = step_name,
        assay_names = assay_names,
        metrics_list_this_step = metrics_list_this_step,
        assay_list = assay_list,
        lipid_id_col = lipid_id_col,
        overwrite = overwrite
    )
    plot_list <- filtering_step_outputs$plot_list

    # --- 9. Directory handling and plot saving --- #
    actual_save_dir <- resolveLipidFilteringPlotSaveDir(
        publication_graphs_dir = publication_graphs_dir,
        omics_type = omics_type,
        time_dir = time_dir
    )
    saveLipidFilteringPlots(
        plot_list = plot_list,
        step_name = step_name,
        actual_save_dir = actual_save_dir,
        return_grid = return_grid,
        publication_graphs_dir = publication_graphs_dir
    )

    returnLipidFilteringPlots(
        plot_list = plot_list,
        return_grid = return_grid
    )
}
























# ----------------------------------------------------------------------------
# lipidIntensityFiltering
# ----------------------------------------------------------------------------
#' @title Lipid Intensity Filtering Method for LipidomicsAssayData
#'
#' @description
#' Filters lipids in *all* assays of a LipidomicsAssayData object.
#' It removes lipids that have intensities below a certain percentile threshold
#' in a proportion of samples exceeding a defined cutoff. The threshold is calculated
#' independently for each assay.
#'
#' @describeIn lipidIntensityFiltering Method for LipidomicsAssayData
#'
#' @param theObject A LipidomicsAssayData object.
#' @param lipids_intensity_cutoff_percentile See generic definition.
#' @param lipids_proportion_of_samples_below_cutoff See generic definition.
#'
#' @importFrom dplyr pull select all_of across
#' @importFrom rlang sym
#' @importFrom stats quantile
#'
#' @return An updated LipidomicsAssayData object.
#' @export
setMethod(
    f = "lipidIntensityFiltering",
    signature = "LipidomicsAssayData",
    definition = function(theObject, lipids_intensity_cutoff_percentile = NULL, lipids_proportion_of_samples_below_cutoff = NULL) {
        # --- Parameter Resolution (Done once) ---
        config_intensity_percentile <- "lipids_intensity_cutoff_percentile"
        raw_intensity_percentile <- checkParamsObjectFunctionSimplify(
            theObject,
            config_intensity_percentile,
            lipids_intensity_cutoff_percentile
        )
        message("Raw intensity percentile from config/param: ", raw_intensity_percentile)
        cleaned_intensity_percentile <- trimws(sub("#.*$", "", raw_intensity_percentile))
        intensity_cutoff_percentile_final <- as.numeric(cleaned_intensity_percentile)

        config_proportion_cutoff <- "lipids_proportion_of_samples_below_cutoff"
        raw_proportion_cutoff <- checkParamsObjectFunctionSimplify(
            theObject,
            config_proportion_cutoff,
            lipids_proportion_of_samples_below_cutoff
        )
        message("Raw proportion cutoff from config/param: ", raw_proportion_cutoff)
        cleaned_proportion_cutoff <- trimws(sub("#.*$", "", raw_proportion_cutoff))
        proportion_of_samples_below_cutoff_final <- as.numeric(cleaned_proportion_cutoff)

        if (is.na(intensity_cutoff_percentile_final)) {
            stop("Failed to convert cleaned lipids_intensity_cutoff_percentile ('", cleaned_intensity_percentile, "' from raw '", raw_intensity_percentile, "') to numeric. Check config.ini or parameter value.")
        }
        if (is.na(proportion_of_samples_below_cutoff_final)) {
            stop("Failed to convert cleaned lipids_proportion_of_samples_below_cutoff ('", cleaned_proportion_cutoff, "' from raw '", raw_proportion_cutoff, "') to numeric. Check config.ini or parameter value.")
        }

        # --- Update Object Parameters (Done once) ---
        theObject <- updateParamInObject(theObject, config_intensity_percentile)
        theObject <- updateParamInObject(theObject, config_proportion_cutoff)

        # --- Process Each Assay in the List ---
        lipid_id_col <- theObject@lipid_id_column
        original_assay_list <- theObject@lipid_data
        original_assay_names <- names(original_assay_list)

        if (length(original_assay_list) == 0) {
            warning("LipidomicsAssayData object has no assays in 'lipid_data' slot. No filtering performed.")
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

            if (!(lipid_id_col %in% names(assay_table))) {
                warning("Lipid ID column '", lipid_id_col, "' not found in assay '", assay_name_for_msg, "'. Skipping this assay.")
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
            min_lipid_intensity_threshold <- ceiling(quantile(all_intensity_values,
                na.rm = TRUE,
                probs = c(intensity_cutoff_percentile_final / 100)
            ))[1]

            message("Calculated minimum intensity threshold for assay '", assay_name_for_msg, "': ", min_lipid_intensity_threshold)

            # Filter using Helper
            filtered_assay <- lipidIntensityFilteringHelper(
                assay_table = assay_table,
                min_lipid_intensity_threshold = min_lipid_intensity_threshold,
                lipids_proportion_of_samples_below_cutoff = proportion_of_samples_below_cutoff_final,
                lipid_id_column = lipid_id_col
            )

            message("Filtered assay '", assay_name_for_msg, "'. Original rows: ", nrow(assay_table), ", Filtered rows: ", nrow(filtered_assay))
            return(filtered_assay)
        })

        # Restore original names if they existed
        if (!is.null(original_assay_names)) {
            names(filtered_assay_list) <- original_assay_names
        }

        # Assign the list of filtered assays back to the object
        theObject@lipid_data <- filtered_assay_list

        # Optional: Call a generic cleanup/design matrix function if applicable
        # theObject <- cleanDesignMatrix(theObject) # If a generic method exists

        return(theObject)
    }
)
