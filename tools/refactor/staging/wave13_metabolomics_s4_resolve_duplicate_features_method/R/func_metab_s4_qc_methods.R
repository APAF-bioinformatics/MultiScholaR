#' @title Resolve Duplicate Features for MetaboliteAssayData
#' @name resolveDuplicateFeatures,MetaboliteAssayData-method
#' @export
setMethod("resolveDuplicateFeatures",
    signature = "MetaboliteAssayData",
    definition = function(theObject, itsd_pattern_columns = NULL) {
        assay_list <- methods::slot(theObject, "metabolite_data")
        id_col <- methods::slot(theObject, "metabolite_id_column") # e.g., database_identifier
        annot_col_slot_name <- methods::slot(theObject, "annotation_id_column") # Name stored in slot
        specific_name_col <- "metabolite" # *** Still assuming this holds the descriptive name ***
        itsd_regex <- methods::slot(theObject, "internal_standard_regex")
        sample_id_col <- methods::slot(theObject, "sample_id")
        design_matrix <- methods::slot(theObject, "design_matrix")

        # --- Determine which columns to check for ITSD pattern --- #
        columns_to_check_for_itsd <- itsd_pattern_columns
        if (is.null(columns_to_check_for_itsd)) {
            columns_to_check_for_itsd <- annot_col_slot_name
            message(sprintf("Parameter `itsd_pattern_columns` is NULL. Defaulting to check column specified in slot `annotation_id_column`: '%s'", annot_col_slot_name))
        } else {
            if (!is.character(columns_to_check_for_itsd)) {
                stop("`itsd_pattern_columns` must be NULL or a character vector.")
            }
            message(sprintf("Checking for ITSD regex in user-specified columns: %s", paste(columns_to_check_for_itsd, collapse = ", ")))
        }
        # ---------------------------------------------------- #

        if (length(assay_list) == 0) {
            warning("No assays found in `metabolite_data` slot. Returning object unchanged.")
            return(theObject)
        }
        if (is.na(itsd_regex) || itsd_regex == "") {
            message("No internal standard regex provided. No special handling for ITSDs will occur.")
        }

        # Ensure list is named
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

            # --- Check required columns (primary ID and specific name col are always needed) ---
            required_cols <- c(id_col, specific_name_col)
            # Also check if the columns targeted for ITSD pattern matching exist
            itsd_check_cols_exist <- columns_to_check_for_itsd %in% colnames(assay_tibble)
            missing_itsd_check_cols <- columns_to_check_for_itsd[!itsd_check_cols_exist]

            if (length(missing_itsd_check_cols) > 0) {
                warning(sprintf(
                    "Specified ITSD pattern column(s) missing from assay '%s': %s. Cannot check these columns for ITSDs.",
                    assay_name, paste(missing_itsd_check_cols, collapse = ", ")
                ), immediate. = TRUE)
                # Adjust columns to check to only those that exist
                columns_to_check_for_itsd <- columns_to_check_for_itsd[itsd_check_cols_exist]
                if (length(columns_to_check_for_itsd) == 0) {
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
            # -------------------------------------------------------------------------------------- #

            # --- Identify Sample Columns --- #
            design_samples <- as.character(design_matrix[[sample_id_col]])
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) == 0) {
                warning(sprintf("No sample columns found matching design matrix. Skipping assay."), immediate. = TRUE)
                return(assay_tibble)
            }
            # -------------------------------- #

            # --- Ensure sample columns are numeric --- #
            assay_tibble_numeric <- assay_tibble %>%
                dplyr::mutate(dplyr::across(dplyr::any_of(sample_cols), as.numeric))
            # --------------------------------------- #

            # --- Add unique row identifier before splitting --- #
            assay_tibble_numeric <- assay_tibble_numeric %>%
                dplyr::mutate(.original_row_id = dplyr::row_number())
            # --------------------------------------------- #

            # --- Separate ITSDs --- #
            is_itsd_present <- FALSE
            itsd_rows <- NULL
            non_itsd_rows <- assay_tibble_numeric # Start assuming all are non-ITSD

            if (!is.na(itsd_regex) && itsd_regex != "" && length(columns_to_check_for_itsd) > 0) {
                # Apply regex to each specified column and combine results
                list_of_matches <- purrr::map(columns_to_check_for_itsd, function(col_name) {
                    values <- assay_tibble_numeric[[col_name]]
                    # Handle potential non-character columns gracefully
                    if (!is.character(values)) {
                        warning(sprintf("Column '%s' used for ITSD check is not character. Coercing.", col_name), immediate. = TRUE)
                        values <- as.character(values)
                    }
                    matches <- stringr::str_detect(values, itsd_regex)
                    matches[is.na(matches)] <- FALSE # Treat NA as no match
                    return(matches)
                })

                # Combine matches with OR logic: TRUE if it matches in *any* column
                itsd_indices <- purrr::reduce(list_of_matches, `|`)

                if (any(itsd_indices)) {
                    message(sprintf("Identified ITSDs using regex '%s' on column(s): %s.", itsd_regex, paste(columns_to_check_for_itsd, collapse = ", ")))
                    is_itsd_present <- TRUE
                    itsd_rows <- assay_tibble_numeric %>% dplyr::filter(itsd_indices)
                    non_itsd_rows <- assay_tibble_numeric %>% dplyr::filter(!itsd_indices)
                    message(sprintf("Separated %d ITSD row(s) and %d non-ITSD row(s).", nrow(itsd_rows), nrow(non_itsd_rows)))
                } else {
                    message(sprintf("No ITSDs found using regex '%s' on column(s): %s.", itsd_regex, paste(columns_to_check_for_itsd, collapse = ", ")))
                }
            }
            # ---------------------- #

            # --- Process ITSDs (if separated) --- #
            processed_itsd_rows <- NULL
            if (is_itsd_present && !is.null(itsd_rows) && nrow(itsd_rows) > 0) {
                message("Processing ITSD rows...")
                # Calculate average intensity first
                itsd_rows_intensity <- itsd_rows %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(
                        avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(avg_intensity = ifelse(is.nan(.data$avg_intensity), -Inf, .data$avg_intensity))

                # Use the specific_name_col (e.g., 'metabolite') as the basis for checking duplicates among ITSDs
                itsd_specific_name_counts <- itsd_rows_intensity %>% dplyr::count(!!rlang::sym(specific_name_col), name = "itsd_specific_count")
                duplicate_itsd_specific_names <- itsd_specific_name_counts %>%
                    dplyr::filter(.data$itsd_specific_count > 1) %>%
                    dplyr::pull(!!rlang::sym(specific_name_col))

                if (length(duplicate_itsd_specific_names) > 0) {
                    message(sprintf(
                        "Found %d duplicated descriptive ITSD name(s) in column '%s': %s",
                        length(duplicate_itsd_specific_names), specific_name_col, paste(duplicate_itsd_specific_names, collapse = ", ")
                    ))
                    message("Appending intensity rank suffix (_1, _2, ...) to descriptive name before updating ID column.")

                    # Process duplicates: Rank within each descriptive name group and create suffixed name
                    duplicates_to_rename <- itsd_rows_intensity %>%
                        dplyr::filter(!!rlang::sym(specific_name_col) %in% duplicate_itsd_specific_names) %>%
                        dplyr::group_by(!!rlang::sym(specific_name_col)) %>%
                        dplyr::arrange(dplyr::desc(.data$avg_intensity), .by_group = TRUE) %>%
                        dplyr::mutate(
                            rank_suffix = paste0("_", dplyr::row_number()),
                            # Create the new, unique ID based on the specific name + suffix
                            .new_itsd_id = paste0(!!rlang::sym(specific_name_col), .data$rank_suffix)
                        ) %>%
                        dplyr::ungroup()

                    # Get non-duplicates
                    non_duplicates_itsd <- itsd_rows_intensity %>%
                        dplyr::filter(!(!!rlang::sym(specific_name_col) %in% duplicate_itsd_specific_names)) %>% # Should not have .new_itsd_id yet
                        dplyr::mutate(.new_itsd_id = !!rlang::sym(specific_name_col)) # Use specific name directly as new ID

                    # Combine and overwrite the main ID column (id_col)
                    processed_itsd_rows_temp <- dplyr::bind_rows(non_duplicates_itsd, duplicates_to_rename) %>%
                        dplyr::mutate(!!rlang::sym(id_col) := .data$.new_itsd_id) %>% # Overwrite id_col (e.g. database_identifier)
                        dplyr::select(-.data$.new_itsd_id, -.data$rank_suffix, -.data$avg_intensity) # Clean up temp columns
                } else {
                    message("No duplicate descriptive ITSD names found.")
                    # Just overwrite id_col with the specific name if no duplicates
                    processed_itsd_rows_temp <- itsd_rows_intensity %>% # Use the one with intensity calculated
                        dplyr::mutate(!!rlang::sym(id_col) := !!rlang::sym(specific_name_col)) %>% # Overwrite id_col
                        dplyr::select(-.data$avg_intensity) # Clean up temp col
                }

                # Final assignment after processing ITSDs
                processed_itsd_rows <- processed_itsd_rows_temp
                message(sprintf("Finished processing ITSD rows. Resulting ITSD rows: %d", nrow(processed_itsd_rows)))
            }
            # --------------------------------- #

            # --- Process Non-ITSDs --- #
            processed_non_itsd_rows <- NULL
            rows_removed_non_itsd <- 0
            if (!is.null(non_itsd_rows) && nrow(non_itsd_rows) > 0) {
                message("Processing non-ITSD rows...")
                # *** Use the correct ID column (id_col) for non-ITSD duplicate checks ***
                non_itsd_id_counts <- non_itsd_rows %>% dplyr::count(!!rlang::sym(id_col), name = "feature_count")
                duplicates_exist_non_itsd <- any(non_itsd_id_counts$feature_count > 1)

                if (!duplicates_exist_non_itsd) {
                    # *** Refer to the correct ID column in the message ***
                    message(sprintf("No duplicates found in non-ITSD ID column '%s'.", id_col))
                    processed_non_itsd_rows <- non_itsd_rows # Already has .original_row_id
                } else {
                    # *** Refer to the correct ID column in the message ***
                    message(sprintf("Resolving duplicates in non-ITSD ID column '%s' by keeping highest average intensity feature...", id_col))
                    # Resolve duplicates using the original intensity logic based on id_col
                    processed_non_itsd_rows <- non_itsd_rows %>%
                        dplyr::rowwise() %>%
                        dplyr::mutate(
                            avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
                        ) %>%
                        dplyr::ungroup() %>%
                        dplyr::mutate(avg_intensity = ifelse(is.nan(.data$avg_intensity), -Inf, .data$avg_intensity)) %>%
                        # *** Group by the correct ID column ***
                        dplyr::group_by(!!rlang::sym(id_col)) %>%
                        dplyr::slice_max(order_by = .data$avg_intensity, n = 1, with_ties = FALSE) %>%
                        dplyr::ungroup() %>%
                        dplyr::select(-.data$avg_intensity)

                    rows_removed_non_itsd <- nrow(non_itsd_rows) - nrow(processed_non_itsd_rows)
                    if (rows_removed_non_itsd > 0) {
                        message(sprintf("Removed %d lower-intensity duplicate non-ITSD feature row(s).", rows_removed_non_itsd))
                    }
                }
                message(sprintf("Finished processing non-ITSD rows. Resulting non-ITSD rows: %d", nrow(processed_non_itsd_rows)))
            }
            # -------------------------- #

            # --- Combine Results and remove helper column --- #
            final_resolved_tibble <- dplyr::bind_rows(processed_itsd_rows, processed_non_itsd_rows) %>%
                dplyr::select(-.original_row_id)
            message(sprintf("Total rows after combining ITSD and non-ITSD: %d", nrow(final_resolved_tibble)))
            message("-----------------------------")
            return(final_resolved_tibble)
        }) %>% purrr::set_names(assay_names)

        # Update the slot in the object
        methods::slot(theObject, "metabolite_data") <- resolved_assay_list

        return(theObject)
    }
)

