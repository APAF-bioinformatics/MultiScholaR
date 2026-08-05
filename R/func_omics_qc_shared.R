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
# func_omics_qc_shared.R
# ============================================================================
# Purpose: Shared Quality Control functions for multiple omics types
#
# This file contains QC functions that are common across different omics
# workflows (e.g., lipidomics, metabolomics) or have been unified to handle
# multiple object types to avoid code duplication and namespace conflicts.
#
# Dependencies:
# - dplyr
# - purrr
# - methods
# ============================================================================

# ----------------------------------------------------------------------------
# findDuplicateFeatureIDs (Unified)
# ----------------------------------------------------------------------------
#' @title Find Duplicate Feature IDs in Omics Assay Data
#'
#' @description
#' Identifies duplicate feature IDs within each assay of an omics data object.
#' Supports both LipidomicsAssayData and MetaboliteAssayData objects.
#'
#' @param theObject A LipidomicsAssayData or MetaboliteAssayData object.
#'
#' @return A named list where each element corresponds to an assay and contains
#'         a tibble of duplicate feature IDs (if any), or NULL.
#'
#' @importFrom dplyr count filter pull %>%
#' @importFrom purrr map set_names
#' @importFrom methods slot
#' @importFrom rlang sym !!
#'
#' @examples
#' \dontrun{
#' # For Lipidomics
#' duplicates_lipid <- findDuplicateFeatureIDs(lipid_assay_obj)
#'
#' # For Metabolomics
#' duplicates_metab <- findDuplicateFeatureIDs(metab_assay_obj)
#' }
#' @export
findDuplicateFeatureIDs <- function(theObject) {
    if (inherits(theObject, "LipidomicsAssayData")) {
        return(findLipidDuplicateFeatureIDs(theObject))
    } else if (inherits(theObject, "MetaboliteAssayData")) {
        return(findMetabDuplicateFeatureIDs(theObject))
    } else {
        stop("findDuplicateFeatureIDs: Unsupported object type: ", class(theObject))
    }
}


# ----------------------------------------------------------------------------
# validateColumnMapping (Unified)
# ----------------------------------------------------------------------------
#' @title Validate Column Mapping for Omics Data
#' @description Checks that required columns exist in the data and returns validation status.
#'              Supports both lipidomics and metabolomics validation.
#'
#' @param data Data frame to validate.
#' @param id_column Name of the primary feature ID column (e.g., lipid_id or metabolite_id).
#' @param sample_columns Character vector of sample column names.
#' @param omics_type Character specifying "lipidomics" or "metabolomics". 
#'                   If NULL, attempts to guess based on context.
#'
#' @return A list containing validation status and summary statistics.
#' @keywords internal
#' @noRd
validateOmicsColumnMapping <- function(data, id_column, sample_columns, omics_type = NULL) {
    # If omics_type is not provided, we check which function to call.
    # We default to metabolomics if unknown, but prefer explicit type.
    if (!is.null(omics_type) && tolower(omics_type) == "lipidomics") {
        return(validateLipidColumnMapping(data, id_column, sample_columns))
    } else {
        # Default or explicit metabolomics
        return(validateMetabColumnMapping(data, id_column, sample_columns))
    }
}
# ----------------------------------------------------------------------------
# resolveDuplicateFeaturesByIntensity
# ----------------------------------------------------------------------------
#' Resolve Duplicate Features by Keeping Highest Average Intensity
#'
#' Within an assay tibble, identifies features with duplicate IDs and keeps only
#' the one with the highest average intensity across sample columns.
#'
#' @param assay_tibble A data frame or tibble representing one assay.
#' @param id_col Character string. The name of the column containing the feature IDs.
#' @param sample_cols Character vector. The names of the columns containing quantitative sample data.
#'
#' @return A tibble with duplicate features resolved based on highest average intensity.
#' @keywords internal
#' @importFrom dplyr group_by summarise ungroup filter slice_max select rowwise mutate c_across any_of
#' @importFrom rlang sym !!
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
resolveDuplicateFeaturesByIntensity <- function(assay_tibble, id_col, sample_cols) {
    if (!id_col %in% colnames(assay_tibble)) {
        warning(sprintf("ID column '%s' not found in assay tibble. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    if (length(sample_cols) == 0) {
        warning("No sample columns provided. Returning original tibble.")
        return(assay_tibble)
    }

    # Check for duplicates first
    id_counts <- assay_tibble %>% dplyr::count(!!rlang::sym(id_col), name = "feature_count")
    duplicates_exist <- any(id_counts$feature_count > 1)

    if (!duplicates_exist) {
        # message(sprintf("No duplicates found in ID column '%s'. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    message(sprintf("Resolving duplicates in ID column '%s' by keeping highest average intensity feature...", id_col))

    # Ensure sample columns are numeric for mean calculation
    assay_tibble_numeric <- assay_tibble %>%
        dplyr::mutate(dplyr::across(dplyr::any_of(sample_cols), as.numeric))

    # Calculate average intensity (handle NAs)
    # Using rowwise is more robust to non-numeric columns than converting to matrix first
    resolved_tibble <- assay_tibble_numeric %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
        ) %>%
        dplyr::ungroup() %>%
        # Handle cases where avg_intensity might be NaN (if all samples are NA)
        dplyr::mutate(avg_intensity = ifelse(is.nan(avg_intensity), -Inf, avg_intensity)) %>%
        # Group by the ID and keep the one with the highest average intensity
        dplyr::group_by(!!rlang::sym(id_col)) %>%
        # slice_max keeps ties by default; with_ties = FALSE ensures only one row per ID
        dplyr::slice_max(order_by = avg_intensity, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        # Remove the temporary average intensity column
        dplyr::select(-avg_intensity)

    # Report how many rows were removed
    rows_removed <- nrow(assay_tibble) - nrow(resolved_tibble)
    if (rows_removed > 0) {
        message(sprintf("Removed %d lower-intensity duplicate feature row(s).", rows_removed))
    }

    return(resolved_tibble)
}
