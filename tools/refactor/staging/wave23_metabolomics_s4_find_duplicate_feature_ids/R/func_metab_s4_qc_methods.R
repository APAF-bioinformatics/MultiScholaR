#' @importFrom dplyr count filter pull %>%
#' @importFrom purrr map set_names
#' @importFrom methods slot
#'
#' @examples
#' \dontrun{
#' # Assuming 'met_assay_obj' is your MetaboliteAssayData object
#' duplicate_ids_list <- findDuplicateFeatureIDs(met_assay_obj)
#' print(duplicate_ids_list)
#'
#' # To get duplicates from the first assay (if any)
#' duplicates_assay1 <- duplicate_ids_list[[1]]
#' print(duplicates_assay1)
#' }
#' @export
findMetabDuplicateFeatureIDs <- function(theObject) {
    if (!inherits(theObject, "MetaboliteAssayData")) {
        stop("Input must be a MetaboliteAssayData object.")
    }

    assay_list <- methods::slot(theObject, "metabolite_data")
    feature_id_col <- methods::slot(theObject, "metabolite_id_column")

    if (length(assay_list) == 0) {
        warning("No assays found in `metabolite_data` slot.")
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

