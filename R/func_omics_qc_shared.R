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
#' @export
validateColumnMapping <- function(data, id_column, sample_columns, omics_type = NULL) {
    # If omics_type is not provided, we check which function to call.
    # We default to metabolomics if unknown, but prefer explicit type.
    if (!is.null(omics_type) && tolower(omics_type) == "lipidomics") {
        return(validateLipidColumnMapping(data, id_column, sample_columns))
    } else {
        # Default or explicit metabolomics
        return(validateMetabColumnMapping(data, id_column, sample_columns))
    }
}
