#' @importFrom dplyr count filter pull %>%
#' @importFrom purrr map set_names
#' @importFrom methods slot
#'
#' @examples
#' \dontrun{
#' # Assuming 'met_assay_obj' is your LipidomicsAssayData object
#' duplicate_ids_list <- findDuplicateFeatureIDs(met_assay_obj)
#' print(duplicate_ids_list)
#'
#' # To get duplicates from the first assay (if any)
#' duplicates_assay1 <- duplicate_ids_list[[1]]
#' print(duplicates_assay1)
#' }
#' @export
findLipidDuplicateFeatureIDs <- function(theObject) {
    if (!inherits(theObject, "LipidomicsAssayData")) {
        stop("Input must be a LipidomicsAssayData object.")
    }

    assay_list <- methods::slot(theObject, "lipid_data")
    feature_id_col <- methods::slot(theObject, "lipid_id_column")

    if (length(assay_list) == 0) {
        warning("No assays found in `lipid_data` slot.")
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

    # Filter out assays with no duplicates for a cleaner output,
    # or keep them as NULL to indicate they were checked.
    # For clarity, let's keep the NULLs
    # duplicate_list <- duplicate_list[!sapply(duplicate_list, is.null)]

    if (all(sapply(duplicate_list, is.null))) {
        message("No duplicate feature IDs found in any assay.")
    }

    return(duplicate_list)
}

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

