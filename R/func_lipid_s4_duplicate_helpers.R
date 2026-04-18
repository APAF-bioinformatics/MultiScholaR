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

# Keep the public S4 method shell in func_lipid_s4_objects.R and route the
# structural implementation through this helper file.
resolveDuplicateFeaturesForLipidObject <- function(theObject, itsd_pattern_columns = NULL) {
    assay_list <- methods::slot(theObject, "lipid_data")
    id_col <- methods::slot(theObject, "lipid_id_column")
    annot_col_slot_name <- methods::slot(theObject, "annotation_id_column")
    specific_name_col <- "lipid"
    itsd_regex <- methods::slot(theObject, "internal_standard_regex")
    sample_id_col <- methods::slot(theObject, "sample_id")
    design_matrix <- methods::slot(theObject, "design_matrix")

    columns_to_check_for_itsd <- itsd_pattern_columns
    if (is.null(columns_to_check_for_itsd)) {
        columns_to_check_for_itsd <- annot_col_slot_name
        message(sprintf(
            "Parameter `itsd_pattern_columns` is NULL. Defaulting to check column specified in slot `annotation_id_column`: '%s'",
            annot_col_slot_name
        ))
    } else {
        if (!is.character(columns_to_check_for_itsd)) {
            stop("`itsd_pattern_columns` must be NULL or a character vector.")
        }
        message(sprintf(
            "Checking for ITSD regex in user-specified columns: %s",
            paste(columns_to_check_for_itsd, collapse = ", ")
        ))
    }

    if (length(assay_list) == 0) {
        warning("No assays found in `lipid_data` slot. Returning object unchanged.")
        return(theObject)
    }
    if (is.na(itsd_regex) || itsd_regex == "") {
        message("No internal standard regex provided. No special handling for ITSDs will occur.")
    }

    assay_names <- names(assay_list)
    if (is.null(assay_names)) {
        assay_names <- paste0("Assay_", seq_along(assay_list))
        names(assay_list) <- assay_names
        warning("Assay list was unnamed. Using default names...", immediate. = TRUE)
    } else if (any(assay_names == "")) {
        needs_name <- which(assay_names == "")
        assay_names[needs_name] <- paste0("Assay_", needs_name)
        names(assay_list) <- assay_names
        warning("Some assays were unnamed. Using default names...", immediate. = TRUE)
    }

    resolved_assay_list <- purrr::map(assay_names, function(assay_name) {
        message(sprintf("-- Processing Assay: %s --", assay_name))
        assay_tibble <- assay_list[[assay_name]]

        required_cols <- c(id_col, specific_name_col)
        itsd_check_cols_exist <- columns_to_check_for_itsd %in% colnames(assay_tibble)
        missing_itsd_check_cols <- columns_to_check_for_itsd[!itsd_check_cols_exist]

        assay_columns_to_check_for_itsd <- columns_to_check_for_itsd
        if (length(missing_itsd_check_cols) > 0) {
            warning(sprintf(
                "Specified ITSD pattern column(s) missing from assay '%s': %s. Cannot check these columns for ITSDs.",
                assay_name, paste(missing_itsd_check_cols, collapse = ", ")
            ), immediate. = TRUE)
            assay_columns_to_check_for_itsd <- assay_columns_to_check_for_itsd[itsd_check_cols_exist]
            if (length(assay_columns_to_check_for_itsd) == 0) {
                message("No valid columns remain for ITSD pattern checking.")
            }
        }

        missing_req_cols <- required_cols[!required_cols %in% colnames(assay_tibble)]
        if (length(missing_req_cols) > 0) {
            warning(sprintf(
                "Required column(s) missing from assay '%s': %s. Skipping assay.",
                assay_name, paste(missing_req_cols, collapse = ", ")
            ), immediate. = TRUE)
            return(assay_tibble)
        }

        design_samples <- as.character(design_matrix[[sample_id_col]])
        all_assay_cols <- colnames(assay_tibble)
        sample_cols <- intersect(all_assay_cols, design_samples)
        if (length(sample_cols) == 0) {
            warning("No sample columns found matching design matrix. Skipping assay.", immediate. = TRUE)
            return(assay_tibble)
        }

        assay_tibble_numeric <- assay_tibble %>%
            dplyr::mutate(dplyr::across(dplyr::any_of(sample_cols), as.numeric)) %>%
            dplyr::mutate(.original_row_id = dplyr::row_number())

        is_itsd_present <- FALSE
        itsd_rows <- NULL
        non_itsd_rows <- assay_tibble_numeric

        if (!is.na(itsd_regex) && itsd_regex != "" && length(assay_columns_to_check_for_itsd) > 0) {
            list_of_matches <- purrr::map(assay_columns_to_check_for_itsd, function(col_name) {
                values <- assay_tibble_numeric[[col_name]]
                if (!is.character(values)) {
                    warning(sprintf(
                        "Column '%s' used for ITSD check is not character. Coercing.",
                        col_name
                    ), immediate. = TRUE)
                    values <- as.character(values)
                }
                matches <- stringr::str_detect(values, itsd_regex)
                matches[is.na(matches)] <- FALSE
                matches
            })

            itsd_indices <- purrr::reduce(list_of_matches, `|`)

            if (any(itsd_indices)) {
                message(sprintf(
                    "Identified ITSDs using regex '%s' on column(s): %s.",
                    itsd_regex,
                    paste(assay_columns_to_check_for_itsd, collapse = ", ")
                ))
                is_itsd_present <- TRUE
                itsd_rows <- assay_tibble_numeric %>% dplyr::filter(itsd_indices)
                non_itsd_rows <- assay_tibble_numeric %>% dplyr::filter(!itsd_indices)
                message(sprintf(
                    "Separated %d ITSD row(s) and %d non-ITSD row(s).",
                    nrow(itsd_rows),
                    nrow(non_itsd_rows)
                ))
            } else {
                message(sprintf(
                    "No ITSDs found using regex '%s' on column(s): %s.",
                    itsd_regex,
                    paste(assay_columns_to_check_for_itsd, collapse = ", ")
                ))
            }
        }

        processed_itsd_rows <- NULL
        if (is_itsd_present && !is.null(itsd_rows) && nrow(itsd_rows) > 0) {
            message("Processing ITSD rows...")
            itsd_rows_intensity <- itsd_rows %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
                ) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(avg_intensity = ifelse(is.nan(.data$avg_intensity), -Inf, .data$avg_intensity))

            itsd_specific_name_counts <- itsd_rows_intensity %>%
                dplyr::count(!!rlang::sym(specific_name_col), name = "itsd_specific_count")
            duplicate_itsd_specific_names <- itsd_specific_name_counts %>%
                dplyr::filter(.data$itsd_specific_count > 1) %>%
                dplyr::pull(!!rlang::sym(specific_name_col))

            if (length(duplicate_itsd_specific_names) > 0) {
                message(sprintf(
                    "Found %d duplicated descriptive ITSD name(s) in column '%s': %s",
                    length(duplicate_itsd_specific_names),
                    specific_name_col,
                    paste(duplicate_itsd_specific_names, collapse = ", ")
                ))
                message("Appending intensity rank suffix (_1, _2, ...) to descriptive name before updating ID column.")

                duplicates_to_rename <- itsd_rows_intensity %>%
                    dplyr::filter(!!rlang::sym(specific_name_col) %in% duplicate_itsd_specific_names) %>%
                    dplyr::group_by(!!rlang::sym(specific_name_col)) %>%
                    dplyr::arrange(dplyr::desc(.data$avg_intensity), .by_group = TRUE) %>%
                    dplyr::mutate(
                        rank_suffix = paste0("_", dplyr::row_number()),
                        .new_itsd_id = paste0(!!rlang::sym(specific_name_col), .data$rank_suffix)
                    ) %>%
                    dplyr::ungroup()

                non_duplicates_itsd <- itsd_rows_intensity %>%
                    dplyr::filter(!(!!rlang::sym(specific_name_col) %in% duplicate_itsd_specific_names)) %>%
                    dplyr::mutate(.new_itsd_id = !!rlang::sym(specific_name_col))

                processed_itsd_rows_temp <- dplyr::bind_rows(non_duplicates_itsd, duplicates_to_rename) %>%
                    dplyr::mutate(!!rlang::sym(id_col) := .data$.new_itsd_id) %>%
                    dplyr::select(-.new_itsd_id, -rank_suffix, -avg_intensity)
            } else {
                message("No duplicate descriptive ITSD names found.")
                processed_itsd_rows_temp <- itsd_rows_intensity %>%
                    dplyr::mutate(!!rlang::sym(id_col) := !!rlang::sym(specific_name_col)) %>%
                    dplyr::select(-avg_intensity)
            }

            processed_itsd_rows <- processed_itsd_rows_temp
            message(sprintf(
                "Finished processing ITSD rows. Resulting ITSD rows: %d",
                nrow(processed_itsd_rows)
            ))
        }

        processed_non_itsd_rows <- NULL
        rows_removed_non_itsd <- 0
        if (!is.null(non_itsd_rows) && nrow(non_itsd_rows) > 0) {
            message("Processing non-ITSD rows...")
            non_itsd_id_counts <- non_itsd_rows %>%
                dplyr::count(!!rlang::sym(id_col), name = "feature_count")
            duplicates_exist_non_itsd <- any(non_itsd_id_counts$feature_count > 1)

            if (!duplicates_exist_non_itsd) {
                message(sprintf("No duplicates found in non-ITSD ID column '%s'.", id_col))
                processed_non_itsd_rows <- non_itsd_rows
            } else {
                message(sprintf(
                    "Resolving duplicates in non-ITSD ID column '%s' by keeping highest average intensity feature...",
                    id_col
                ))
                processed_non_itsd_rows <- non_itsd_rows %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(
                        avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(avg_intensity = ifelse(is.nan(.data$avg_intensity), -Inf, .data$avg_intensity)) %>%
                    dplyr::group_by(!!rlang::sym(id_col)) %>%
                    dplyr::slice_max(order_by = .data$avg_intensity, n = 1, with_ties = FALSE) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(-avg_intensity)

                rows_removed_non_itsd <- nrow(non_itsd_rows) - nrow(processed_non_itsd_rows)
                if (rows_removed_non_itsd > 0) {
                    message(sprintf(
                        "Removed %d lower-intensity duplicate non-ITSD feature row(s).",
                        rows_removed_non_itsd
                    ))
                }
            }
            message(sprintf(
                "Finished processing non-ITSD rows. Resulting non-ITSD rows: %d",
                nrow(processed_non_itsd_rows)
            ))
        }

        final_resolved_tibble <- dplyr::bind_rows(processed_itsd_rows, processed_non_itsd_rows) %>%
            dplyr::select(-.original_row_id)
        message(sprintf(
            "Total rows after combining ITSD and non-ITSD: %d",
            nrow(final_resolved_tibble)
        ))
        message("-----------------------------")
        final_resolved_tibble
    }) %>%
        purrr::set_names(assay_names)

    methods::slot(theObject, "lipid_data") <- resolved_assay_list
    theObject
}
