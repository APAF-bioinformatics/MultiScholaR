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
    # Determine object type and extract relevant slots
    if (inherits(theObject, "LipidomicsAssayData")) {
        assay_list <- methods::slot(theObject, "lipid_data")
        feature_id_col <- methods::slot(theObject, "lipid_id_column")
        obj_type <- "Lipidomics"
    } else if (inherits(theObject, "MetaboliteAssayData")) {
        assay_list <- methods::slot(theObject, "metabolite_data")
        feature_id_col <- methods::slot(theObject, "metabolite_id_column")
        obj_type <- "Metabolomics"
    } else {
        stop("Input must be a LipidomicsAssayData or MetaboliteAssayData object.")
    }

    if (length(assay_list) == 0) {
        warning(sprintf("No assays found in %s data slot.", obj_type))
        return(list())
    }

    # Ensure list is named
    assay_names <- names(assay_list)
    if (is.null(assay_names)) {
        assay_names <- paste0("Assay_", seq_along(assay_list))
        names(assay_list) <- assay_names
        warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
    } else if (any(assay_names == "")) {
        needs_name <- which(assay_names == "")
        assay_names[needs_name] <- paste0("Assay_", needs_name)
        names(assay_list) <- assay_names
        warning("Some assays were unnamed. Using default names for them.")
    }

    duplicate_list <- purrr::map(assay_names, function(assay_name) {
        current_assay_data <- assay_list[[assay_name]]

        # Check if the feature ID column exists
        if (!feature_id_col %in% colnames(current_assay_data)) {
            warning(sprintf(
                "Assay '%s': Feature ID column '%s' not found. Skipping duplicate check.",
                assay_name, feature_id_col
            ))
            return(NULL)
        }

        # Find duplicates
        id_counts <- current_assay_data %>%
            dplyr::count(!!rlang::sym(feature_id_col), name = "count")

        duplicates_found <- id_counts %>%
            dplyr::filter(.data$count > 1)

        if (nrow(duplicates_found) > 0) {
            message(sprintf("Duplicates found in Assay: '%s' (Column: '%s')", assay_name, feature_id_col))
            return(duplicates_found)
        } else {
            # message(sprintf("No duplicates found in Assay: '%s' (Column: '%s')", assay_name, feature_id_col))
            return(NULL) # Return NULL if no duplicates
        }
    }) %>%
        purrr::set_names(assay_names) # Set names for the final list

    if (all(sapply(duplicate_list, is.null))) {
        message("No duplicate feature IDs found in any assay.")
    }

    return(duplicate_list)
}
