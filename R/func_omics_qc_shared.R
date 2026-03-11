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
# findDuplicateFeatureIDs <- function(theObject) {
#     # Deprecated in favor of findMetabDuplicateFeatureIDs and findLipidDuplicateFeatureIDs
#     # to resolve namespace collisions in the Shiny apps.
#     stop("findDuplicateFeatureIDs is deprecated. Use findMetabDuplicateFeatureIDs or findLipidDuplicateFeatureIDs.")
# }
