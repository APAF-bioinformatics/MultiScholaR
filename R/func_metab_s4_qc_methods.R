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

#' @title Filter Samples by Metabolite Correlation Threshold
#' @name filterSamplesByMetaboliteCorrelationThreshold,MetaboliteAssayData-method
#' @description Removes samples from a MetaboliteAssayData object based on
#'   Pearson correlation thresholds. Samples that do not have at least one
#'   replicate pair with correlation above the threshold are removed.
#' @param theObject A MetaboliteAssayData object
#' @param pearson_correlation_per_pair A list of data frames (one per assay)
#'   containing pair-wise correlation results from \code{pearsonCorForSamplePairs}.
#' @param min_pearson_correlation_threshold A numeric value (0-1). Samples with
#'   correlation below this threshold in their replicate group are removed.
#' @return An updated MetaboliteAssayData object with poorly correlated samples removed.
#' @importFrom dplyr filter select distinct pull
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @export
setMethod(
    f = "filterSamplesByMetaboliteCorrelationThreshold",
    signature = "MetaboliteAssayData",
    definition = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = 0.5) {
        message("+===========================================================================+")
        message("|        Metabolite Sample Filtering by Correlation Threshold (S4)          |")
        message("+===========================================================================+")

        if (is.null(pearson_correlation_per_pair) || !is.list(pearson_correlation_per_pair)) {
            stop("`pearson_correlation_per_pair` must be a list of correlation data frames (one per assay).")
        }

        if (is.null(min_pearson_correlation_threshold) || !is.numeric(min_pearson_correlation_threshold)) {
            stop("`min_pearson_correlation_threshold` must be a numeric value.")
        }

        design_matrix <- theObject@design_matrix
        sample_id_col_name <- theObject@sample_id
        assay_list <- theObject@metabolite_data

        if (length(assay_list) == 0) {
            warning("No assays found in MetaboliteAssayData object.")
            return(theObject)
        }

        samples_to_remove_total <- character()

        filtered_assay_list <- purrr::map2(assay_list, pearson_correlation_per_pair, function(current_assay_data, correlation_results) {
            if (is.null(correlation_results) || nrow(correlation_results) == 0) {
                warning("No correlation results provided for assay. Skipping filtering.")
                return(current_assay_data)
            }

            run_id_col_x <- paste0(sample_id_col_name, ".x")
            run_id_col_y <- paste0(sample_id_col_name, ".y")

            if (!all(c(run_id_col_x, run_id_col_y, "pearson_correlation") %in% colnames(correlation_results))) {
                warning("Correlation results table missing expected columns. Skipping filtering.")
                return(current_assay_data)
            }

            all_samples_in_analysis <- correlation_results |>
                tidyr::pivot_longer(cols = c(!!rlang::sym(run_id_col_x), !!rlang::sym(run_id_col_y)), values_to = "sample_id") |>
                dplyr::distinct(sample_id) |>
                dplyr::pull(sample_id)

            passing_pairs <- correlation_results |>
                dplyr::filter(pearson_correlation >= min_pearson_correlation_threshold)

            samples_to_keep <- passing_pairs |>
                tidyr::pivot_longer(cols = c(!!rlang::sym(run_id_col_x), !!rlang::sym(run_id_col_y)), values_to = "sample_id") |>
                dplyr::distinct(sample_id) |>
                dplyr::pull(sample_id)

            samples_to_remove <- setdiff(all_samples_in_analysis, samples_to_keep)
            samples_to_remove_total <<- c(samples_to_remove_total, samples_to_remove)

            if (length(samples_to_remove) > 0) {
                message(sprintf(
                    "  Removing %d samples below correlation threshold: %s",
                    length(samples_to_remove), paste(samples_to_remove, collapse = ", ")
                ))
                cols_to_keep <- setdiff(colnames(current_assay_data), samples_to_remove)
                current_assay_data <- current_assay_data[, cols_to_keep, drop = FALSE]
            } else {
                message("  No samples below correlation threshold.")
            }

            return(current_assay_data)
        })

        names(filtered_assay_list) <- names(assay_list)
        theObject@metabolite_data <- filtered_assay_list

        if (length(samples_to_remove_total) > 0) {
            samples_to_remove_unique <- unique(samples_to_remove_total)
            theObject@design_matrix <- design_matrix |>
                dplyr::filter(!(!!rlang::sym(sample_id_col_name) %in% samples_to_remove_unique))
            message(sprintf("Total samples removed across all assays: %d", length(samples_to_remove_unique)))
        }

        message("+===========================================================================+")

        return(theObject)
    }
)

#' @title Calculate Pearson Correlation for Sample Pairs
#' @name pearsonCorForSamplePairs,MetaboliteAssayData-method
#' @importFrom purrr map set_names map_df
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join select mutate filter distinct arrange group_by summarise ungroup pull if_else case_when row_number n rename relocate
#' @importFrom tibble add_column
#' @importFrom stringr str_detect
#' @importFrom rlang sym !! :=
#' @export
#' @export
setMethod(
    f = "pearsonCorForSamplePairs",
    signature = "MetaboliteAssayData",
    definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
        message("+===========================================================================+")
        message("|  DEBUG66: Entering pearsonCorForSamplePairs (MetaboliteAssayData)        |")
        message("+===========================================================================+")

        # --- Input Validation ---
        # tech_rep_remove_regex can be NULL, checked inside helper/later use
        # correlation_group can be NA, checked below

        design_matrix <- theObject@design_matrix
        sample_id_col_name <- theObject@sample_id
        metabolite_id_col_name <- theObject@metabolite_id_column
        tech_rep_col_name <- theObject@technical_replicate_id # Default grouping if correlation_group is NA
        assay_list <- theObject@metabolite_data

        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] sample_id_col = '%s', metabolite_id_col = '%s'", sample_id_col_name, metabolite_id_col_name))
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] tech_rep_col_name = '%s', correlation_group = '%s'", tech_rep_col_name, correlation_group))
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Number of assays: %d", length(assay_list)))

        if (length(assay_list) == 0) {
            warning("No assays found in `metabolite_data` slot. Returning empty list.")
            return(list())
        }

        # Ensure list is named
        if (is.null(names(assay_list))) {
            names(assay_list) <- paste0("Assay_", seq_along(assay_list))
            warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
        }
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay names: %s", paste(names(assay_list), collapse = ", ")))

        # Determine the actual grouping column to use for pairing samples
        replicate_group_column_name <- correlation_group
        if (is.na(correlation_group)) {
            message("   DEBUG66 [pearsonCorForSamplePairs] correlation_group is NA, falling back to tech_rep_col")
            if (is.na(tech_rep_col_name) || !tech_rep_col_name %in% colnames(design_matrix)) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] FAIL - tech_rep_col '%s' not valid", tech_rep_col_name))
                stop("`correlation_group` is NA and `technical_replicate_id` ('", tech_rep_col_name, "') is NA or not found in design_matrix. Cannot determine sample pairing.")
            }
            replicate_group_column_name <- tech_rep_col_name
        } else {
            if (!correlation_group %in% colnames(design_matrix)) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] FAIL - correlation_group '%s' not in design_matrix", correlation_group))
                stop(sprintf("Specified `correlation_group` ('%s') not found in design_matrix.", correlation_group))
            }
        }
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Using replicate_group_column_name = '%s'", replicate_group_column_name))

        # Resolve tech_rep_remove_regex from config if needed (or use default)
        # Assuming the helper function or subsequent filtering handles NULL regex gracefully (meaning no filtering)
        tech_rep_remove_regex_final <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", tech_rep_remove_regex) # Allow override
        # theObject <- updateParamInObject(theObject, "tech_rep_remove_regex") # Update object if needed


        # --- Correlation Logic per Assay ---
        correlation_results_list <- purrr::map(seq_along(assay_list), function(i) {
            assay_name <- names(assay_list)[i]
            current_assay_data <- assay_list[[i]]
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] === Processing assay: %s ===", assay_name))

            # --- Correctly identify sample columns based on design matrix ---
            design_samples <- as.character(design_matrix[[sample_id_col_name]]) # Get sample IDs from design matrix
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_samples count = %d", assay_name, length(design_samples)))
            all_assay_cols <- colnames(current_assay_data)
            sample_cols <- intersect(all_assay_cols, design_samples) # Find which design samples are columns in the assay
            metadata_cols <- setdiff(all_assay_cols, sample_cols) # All other columns are metadata/ID
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': sample_cols count = %d, metadata_cols count = %d", assay_name, length(sample_cols), length(metadata_cols)))

            # Ensure the primary metabolite ID column is considered metadata
            if (metabolite_id_col_name %in% sample_cols) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' is also listed as a sample ID. Check configuration.", assay_name, metabolite_id_col_name))
            }
            metadata_cols <- union(metadata_cols, metabolite_id_col_name) # Ensure metabolite ID is not treated as a sample column
            sample_cols <- setdiff(all_assay_cols, metadata_cols) # Final list of sample columns
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': final sample_cols count = %d", assay_name, length(sample_cols)))

            if (length(sample_cols) < 2) { # Need at least 2 samples for correlation
                warning(sprintf("Assay '%s': Fewer than 2 sample columns found matching design matrix sample IDs. Skipping Pearson correlation.", assay_name))
                return(NULL)
            }
            # --- End Correction ---


            # Check sample consistency (now based on correctly identified sample_cols)
            design_samples_check <- design_matrix[[sample_id_col_name]] # Use original type for check
            missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
            if (length(missing_samples_in_design) > 0) {
                # This condition should theoretically not be met if sample_cols were derived from design_samples,
                # but keeping as a safeguard against type issues or unexpected data.
                warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping Pearson correlation.", assay_name, paste(missing_samples_in_design, collapse = ", ")))
                return(NULL)
            }


            # Filter design matrix to match assay samples
            design_matrix_filtered <- design_matrix[design_matrix[[sample_id_col_name]] %in% sample_cols, ]
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_matrix_filtered rows = %d", assay_name, nrow(design_matrix_filtered)))

            # Ensure metabolite ID column exists
            if (!metabolite_id_col_name %in% colnames(current_assay_data)) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping Pearson correlation.", assay_name, metabolite_id_col_name))
                return(NULL)
            }


            # Prepare long data for helper
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Pivoting to long format...", assay_name))
            assay_long <- current_assay_data |>
                tidyr::pivot_longer(
                    cols = all_of(sample_cols),
                    names_to = sample_id_col_name,
                    values_to = "abundance"
                )
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': assay_long rows = %d", assay_name, nrow(assay_long)))


            # Prepare the design matrix subset for the helper
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Creating design_subset with cols '%s' and '%s'", assay_name, sample_id_col_name, replicate_group_column_name))
            design_subset <- design_matrix_filtered |>
                dplyr::select(!!rlang::sym(sample_id_col_name), !!rlang::sym(replicate_group_column_name))
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_subset rows = %d, cols = %s", assay_name, nrow(design_subset), paste(colnames(design_subset), collapse = ", ")))


            # --- Ensure consistent Sample ID type (character) --- #
            # Convert the sample ID column in the design subset to character
            # to match the type expected from pivot_longer names_to
            design_subset <- design_subset |>
                dplyr::mutate(!!rlang::sym(sample_id_col_name) := as.character(!!rlang::sym(sample_id_col_name)))

            # Also ensure the assay_long sample ID column is character (pivot_longer usually does this)
            assay_long <- assay_long |>
                dplyr::mutate(!!rlang::sym(sample_id_col_name) := as.character(!!rlang::sym(sample_id_col_name)))
            # ---------------------------------------------------- #


            # --- Calculate Correlations Directly --- #

            # 1. Get pairs of samples to compare based on the replicate grouping column
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Calling getPairsOfSamplesTable...", assay_name))
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_subset unique groups = %s", assay_name, paste(unique(design_subset[[replicate_group_column_name]]), collapse = ", ")))
            pairs_for_comparison <- tryCatch(
                {
                    getPairsOfSamplesTable(design_subset, # Contains sample_id and replicate_group_column
                        run_id_column = sample_id_col_name,
                        replicate_group_column = replicate_group_column_name
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': ERROR in getPairsOfSamplesTable: %s", assay_name, e$message))
                    warning(sprintf("Assay '%s': Error getting sample pairs: %s. Skipping correlation.", assay_name, e$message))
                    return(NULL)
                }
            )

            if (is.null(pairs_for_comparison) || nrow(pairs_for_comparison) == 0) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': No valid sample pairs found. Skipping.", assay_name))
                warning(sprintf("Assay '%s': No valid sample pairs found for correlation. Skipping.", assay_name))
                return(NULL)
            }
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': pairs_for_comparison rows = %d", assay_name, nrow(pairs_for_comparison)))

            # Get the names of the columns containing paired sample IDs (e.g., "Run.x", "Run.y")
            run_id_col_x <- paste0(sample_id_col_name, ".x")
            run_id_col_y <- paste0(sample_id_col_name, ".y")

            # Check if these columns exist in the pairs table
            if (!all(c(run_id_col_x, run_id_col_y) %in% colnames(pairs_for_comparison))) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Missing expected columns. pairs_for_comparison cols = %s", assay_name, paste(colnames(pairs_for_comparison), collapse = ", ")))
                warning(sprintf("Assay '%s': Expected paired sample columns ('%s', '%s') not found in pairs table. Skipping correlation.", assay_name, run_id_col_x, run_id_col_y))
                return(NULL)
            }

            # Calculate correlations as a separate vector first
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Calculating correlations for %d pairs...", assay_name, nrow(pairs_for_comparison)))
            calculated_correlations <- tryCatch(
                {
                    purrr::map2_dbl(
                        .x = pairs_for_comparison[[run_id_col_x]], # Directly access columns
                        .y = pairs_for_comparison[[run_id_col_y]], # Directly access columns
                        .f = ~ {
                            # Filter the long assay data for the current pair
                            assay_pair_filtered <- assay_long |>
                                dplyr::filter(!!rlang::sym(sample_id_col_name) %in% c(.x, .y))

                            # Call the new metabolite-specific helper
                            correlation_val <- calculateMetabolitePairCorrelation(
                                input_pair_table = assay_pair_filtered,
                                feature_id_column = metabolite_id_col_name,
                                sample_id_column = sample_id_col_name,
                                value_column = "abundance"
                            )
                            return(correlation_val) # Explicit return
                        }
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': ERROR in map2_dbl: %s", assay_name, e$message))
                    warning(sprintf("Assay '%s': Error during map2_dbl correlation calculation: %s. Returning NULL results.", assay_name, e$message))
                    return(NULL) # Return NULL if map2_dbl fails
                }
            )

            # Check if calculation succeeded and add the column
            if (is.null(calculated_correlations)) {
                correlation_results_raw <- NULL # Propagate failure
            } else if (length(calculated_correlations) != nrow(pairs_for_comparison)) {
                warning(sprintf("Assay '%s': Number of calculated correlations (%d) does not match number of pairs (%d). Skipping.", assay_name, length(calculated_correlations), nrow(pairs_for_comparison)))
                correlation_results_raw <- NULL
            } else {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Successfully calculated %d correlations", assay_name, length(calculated_correlations)))
                correlation_results_raw <- pairs_for_comparison |>
                    dplyr::mutate(pearson_correlation = calculated_correlations)
            }
            # ---------------------------------------------------------- #

            # Initialize correlation_results_filtered with raw results (will be filtered if regex provided)
            correlation_results_filtered <- correlation_results_raw

            if (is.null(correlation_results_raw)) {
                # If calculation failed earlier, correlation_results_raw is NULL
                return(NULL)
            } else if (!is.null(tech_rep_remove_regex_final) && tech_rep_remove_regex_final != "") {
                # Ensure the replicate group column exists in the result before filtering
                if (replicate_group_column_name %in% colnames(correlation_results_raw)) {
                    correlation_results_filtered <- correlation_results_raw |>
                        dplyr::filter(!stringr::str_detect(!!rlang::sym(replicate_group_column_name), tech_rep_remove_regex_final))
                } else {
                    warning(sprintf("Assay '%s': Replicate group column '%s' not found in correlation results. Cannot apply `tech_rep_remove_regex`. Returning unfiltered results.", assay_name, replicate_group_column_name))
                    correlation_results_filtered <- correlation_results_raw
                }
            }
            # If no regex, correlation_results_filtered already holds correlation_results_raw

            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Returning %d correlation results", assay_name, nrow(correlation_results_filtered)))
            return(correlation_results_filtered)
        })

        # Set names for the list of results
        names(correlation_results_list) <- names(assay_list)

        # Remove NULL elements (skipped assays)
        non_null_count_before <- sum(!sapply(correlation_results_list, is.null))
        correlation_results_list <- correlation_results_list[!sapply(correlation_results_list, is.null)]

        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Finished. Returning %d assay results (removed %d NULL)", length(correlation_results_list), length(assay_list) - non_null_count_before))
        message("+===========================================================================+")
        message("|  DEBUG66: Exiting pearsonCorForSamplePairs                               |")
        message("+===========================================================================+")

        return(correlation_results_list)
    }
)

#' @title Plot Pearson Correlation
#' @name plotPearson,MetaboliteAssayData-method
#' @importFrom purrr map set_names
#' @importFrom ggplot2 ggplot aes geom_histogram scale_y_continuous xlab ylab theme element_blank
#' @export
setMethod(
    f = "plotPearson",
    signature = "MetaboliteAssayData",
    definition = function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {
        # Get the list of correlation tibbles (one per assay)
        # tech_rep_remove_regex and correlation_group are passed down
        correlation_list <- pearsonCorForSamplePairs(theObject,
            tech_rep_remove_regex = tech_rep_remove_regex,
            correlation_group = correlation_group
        )

        if (length(correlation_list) == 0) {
            warning("No correlation results generated (likely no valid assays). Returning empty list.")
            return(list())
        }

        # Ensure list is named (pearsonCorForSamplePairs should have handled this, but double-check)
        if (is.null(names(correlation_list))) {
            names(correlation_list) <- paste0("Assay_", seq_along(correlation_list))
        }


        # --- Plotting Logic per Assay's Correlation Results ---
        pearson_plots_list <- purrr::map(seq_along(correlation_list), function(i) {
            assay_name <- names(correlation_list)[i]
            correlation_vec <- correlation_list[[i]]

            # Check if the correlation data is valid
            if (is.null(correlation_vec) || nrow(correlation_vec) == 0 || !"pearson_correlation" %in% colnames(correlation_vec)) {
                warning(sprintf("Assay '%s': Invalid or empty correlation data provided. Skipping Pearson plot.", assay_name))
                return(NULL)
            }

            # Check for all NA values
            if (all(is.na(correlation_vec$pearson_correlation))) {
                warning(sprintf("Assay '%s': All Pearson correlation values are NA. Skipping plot.", assay_name))
                return(NULL)
            }

            # Calculate breaks carefully, handling potential NAs and edge cases
            min_cor <- min(correlation_vec$pearson_correlation, na.rm = TRUE)
            # Ensure min_cor is finite; default if not
            if (!is.finite(min_cor)) min_cor <- 0

            # Use finer breaks, similar to protein version, clamped to [0, 1]
            # Note: Protein version uses 0.001 step, using 0.01 here for potentially better visibility first.
            hist_breaks <- seq(0, 1, 0.01)

            # --- Create Plot ---
            tryCatch(
                {
                    pearson_plot <- correlation_vec |>
                        ggplot(aes(pearson_correlation)) +
                        geom_histogram(breaks = hist_breaks, na.rm = TRUE) +
                        # Set x-axis limits and breaks
                        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
                        # # Set fixed y-axis scale, similar to protein version
                        # scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
                        xlab("Pearson Correlation") +
                        ylab("Counts") +
                        # ggtitle(paste(assay_name)) +
                        theme_bw() +
                        theme(
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()
                        )

                    return(pearson_plot)
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error creating Pearson histogram: %s. Skipping.", assay_name, e$message))
                    return(NULL)
                }
            )
        })

        # Set names for the list of plots
        names(pearson_plots_list) <- names(correlation_list)

        # Remove NULL elements (skipped assays)
        pearson_plots_list <- pearson_plots_list[!sapply(pearson_plots_list, is.null)]

        return(pearson_plots_list)
    }
)

# --- Internal Helper for Metabolite Pair Correlation --- #
#' Calculate Pearson correlation for a pair of samples from long-format data
#'
#' Internal helper function specifically for metabolomics data structure.
#' Assumes input table is already filtered for the two relevant samples.
#'
#' @param input_pair_table Tibble in long format with columns for feature ID,
#'   sample ID, and abundance values. Must contain exactly two unique sample IDs.
#' @param feature_id_column String name of the column containing feature IDs.
#' @param sample_id_column String name of the column containing sample IDs.
#' @param value_column String name of the column containing abundance values.
#'
#' @return Numeric Pearson correlation value, or NA_real_ on error or insufficient data.
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym !!
#' @importFrom stats cor
#' @keywords internal
#' @export
calculateMetabolitePairCorrelation <- function(input_pair_table, feature_id_column, sample_id_column, value_column) {
    # Get the two unique sample IDs from the input table
    sample_ids <- unique(input_pair_table[[sample_id_column]])
    if (length(sample_ids) != 2) {
        warning("calculateMetabolitePairCorrelation: Input table does not contain exactly two samples.")
        return(NA_real_)
    }
    sample_x_id <- sample_ids[1]
    sample_y_id <- sample_ids[2]

    # Pivot wider to get features as rows and the two samples as columns
    wide_pair_table <- tryCatch(
        {
            input_pair_table |>
                dplyr::select(!!rlang::sym(feature_id_column), !!rlang::sym(sample_id_column), !!rlang::sym(value_column)) |>
                tidyr::pivot_wider(
                    names_from = !!rlang::sym(sample_id_column),
                    values_from = !!rlang::sym(value_column)
                )
        },
        error = function(e) {
            warning(sprintf("Error pivoting data wider for correlation between %s and %s: %s", sample_x_id, sample_y_id, e$message))
            return(NULL)
        }
    )

    if (is.null(wide_pair_table) || nrow(wide_pair_table) < 2) {
        # Need at least 2 features for correlation
        return(NA_real_)
    }

    # --- Added Check ---
    # Check if expected columns exist after pivot (using the character IDs)
    expected_colnames <- as.character(sample_ids)
    if (!all(expected_colnames %in% colnames(wide_pair_table))) {
        warning(sprintf("Expected sample columns %s or %s not found after pivoting.", expected_colnames[1], expected_colnames[2]))
        return(NA_real_)
    }
    # --- End Added Check ---

    # Extract the value vectors for the two samples
    # Column names will be the actual sample IDs (e.g., "51581", "51582")
    values_x <- wide_pair_table[[expected_colnames[1]]] # Use verified name
    values_y <- wide_pair_table[[expected_colnames[2]]] # Use verified name

    # Calculate correlation
    cor_result <- tryCatch(
        {
            stats::cor(values_x, values_y, use = "pairwise.complete.obs")
        },
        error = function(e) {
            warning(sprintf("calculateMetabolitePairCorrelation: Error in stats::cor for samples %s and %s: %s", sample_x_id, sample_y_id, e$message))
            return(NA_real_) # Returns NA_real_ on cor error
        }
    )

    # --- Modified Check ---
    # Ensure the result is a single, finite numeric value
    if (length(cor_result) != 1 || !is.numeric(cor_result) || !is.finite(cor_result)) {
        return(NA_real_)
    }
    # --- End Modified Check ---

    return(cor_result)
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

# Optional: Add other normalization methods (e.g., PQN, Median) here later
# setMethod(f = "normaliseUntransformedData",
#           signature = signature(theObject = "MetaboliteAssayData", method = "character"),
#           definition = function(theObject, method = "PQN", ...) { ... }
# )
# ==========================================
# Content from metabolite_qc.R
# ==========================================
#' @title Helper function for metabolite intensity filtering
#' @name metaboliteIntensityFilteringHelper
#' @description Filter metabolites based on an intensity threshold and the proportion of samples below that threshold in a wide-format table.
#' @param assay_table A wide data frame where rows are metabolites and columns include a metabolite identifier and numeric sample intensities.
#' @param min_metabolite_intensity_threshold The calculated minimum intensity value. Metabolites in samples below this threshold are considered 'below threshold'.
#' @param metabolites_proportion_of_samples_below_cutoff The maximum allowed proportion (0 to 1) of samples where a metabolite can be below the threshold. If a metabolite exceeds this proportion, it's removed.
#' @param metabolite_id_column A string specifying the name of the column containing the unique metabolite identifiers.
#' @return A filtered wide data frame containing only the metabolites that pass the filter.
#' @export
metaboliteIntensityFilteringHelper <- function(
  assay_table,
  min_metabolite_intensity_threshold,
  metabolites_proportion_of_samples_below_cutoff,
  metabolite_id_column
) {
    # Identify numeric columns representing sample intensities
    sample_cols <- names(assay_table)[sapply(assay_table, is.numeric)]
    num_samples <- length(sample_cols)

    if (num_samples == 0) {
        warning("No numeric sample columns found in the assay table. Returning original table.")
        return(assay_table)
    }

    # Calculate the number of samples below threshold for each metabolite
    metabolites_below_threshold <- assay_table |>
        # Ensure id column is character for safe rowwise operations if needed
        # mutate({{metabolite_id_column}} := as.character({{metabolite_id_column}})) |>
        rowwise() |>
        mutate(
            num_below_threshold = sum(c_across(all_of(sample_cols)) < min_metabolite_intensity_threshold, na.rm = TRUE),
            proportion_below_threshold = num_below_threshold / num_samples
        ) |>
        ungroup()

    # Filter metabolites based on the proportion cutoff
    filtered_assay_table <- metabolites_below_threshold |>
        dplyr::filter(proportion_below_threshold < metabolites_proportion_of_samples_below_cutoff) |>
        # Remove the temporary calculation columns
        dplyr::select(-num_below_threshold, -proportion_below_threshold)

    return(filtered_assay_table)
}

