#' Filter Metabolites by Intensity and Sample Proportion (Helper)
#'
#' @description
#' A helper function that filters rows (metabolites) from a wide-format assay table.
#' The filtering is based on a minimum intensity threshold and the maximum allowed
#' proportion of samples that can fall below this threshold for any given metabolite.
#'
#' @details
#' For each metabolite (row), this function calculates the number and proportion of
#' samples that have an intensity value less than `min_metabolite_intensity_threshold`.
#' If this proportion is greater than or equal to `metabolites_proportion_of_samples_below_cutoff`,
#' the metabolite is removed from the table.
#'
#' This function is typically called by the `metaboliteIntensityFiltering` S4 method.
#'
#' @param assay_table A wide data frame where rows are metabolites and columns
#'   include a metabolite identifier and numeric sample intensities.
#' @param min_metabolite_intensity_threshold A numeric value. Intensities below
#'   this threshold are considered 'below threshold'.
#' @param metabolites_proportion_of_samples_below_cutoff A numeric value between 0 and 1.
#'   The maximum allowed proportion of samples where a metabolite can be below the
#'   threshold. If a metabolite's proportion of low-intensity samples is greater than
#'   or equal to this value, it is removed.
#' @param metabolite_id_column A character string specifying the name of the column
#'   containing unique metabolite identifiers. This column is ignored during calculations.
#'
#' @return A filtered wide data frame containing only the metabolites that pass the filter.
#'   The returned data frame has the same columns as the input `assay_table`.
#'
#' @importFrom dplyr rowwise mutate c_across all_of ungroup filter select
#' @export
#' @keywords internal
metaboliteIntensityFilteringHelper <- function(assay_table
                                               , min_metabolite_intensity_threshold
                                               , metabolites_proportion_of_samples_below_cutoff
                                               , metabolite_id_column) {

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
      num_below_threshold = sum(c_across(all_of(sample_cols)) < min_metabolite_intensity_threshold, na.rm = TRUE)
      , proportion_below_threshold = num_below_threshold / num_samples
    ) |>
    ungroup()

  # Filter metabolites based on the proportion cutoff
  filtered_assay_table <- metabolites_below_threshold |>
    dplyr::filter(proportion_below_threshold < metabolites_proportion_of_samples_below_cutoff) |>
    # Remove the temporary calculation columns
    dplyr::select(-num_below_threshold, -proportion_below_threshold)

  return(filtered_assay_table)
}

#-------------------------------------------------------------------------------

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
setMethod( f="metaboliteIntensityFiltering"
           , signature="MetaboliteAssayData"
           , definition = function( theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {

             # --- Parameter Resolution (Done once) ---
             config_intensity_percentile <- "metabolites_intensity_cutoff_percentile"
             raw_intensity_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                           , config_intensity_percentile
                                                                           , metabolites_intensity_cutoff_percentile)
             message("Raw intensity percentile from config/param: ", raw_intensity_percentile)
             cleaned_intensity_percentile <- trimws(sub("#.*$", "", raw_intensity_percentile))
             intensity_cutoff_percentile_final <- as.numeric(cleaned_intensity_percentile)

             config_proportion_cutoff <- "metabolites_proportion_of_samples_below_cutoff"
             raw_proportion_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                        , config_proportion_cutoff
                                                                        , metabolites_proportion_of_samples_below_cutoff)
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
                 min_metabolite_intensity_threshold <- ceiling( quantile( all_intensity_values
                                                                          , na.rm=TRUE
                                                                          , probs = c(intensity_cutoff_percentile_final/100) ))[1]

                 message("Calculated minimum intensity threshold for assay '", assay_name_for_msg, "': ", min_metabolite_intensity_threshold)

                 # Filter using Helper
                 filtered_assay <- metaboliteIntensityFilteringHelper(assay_table = assay_table
                                                                      , min_metabolite_intensity_threshold = min_metabolite_intensity_threshold
                                                                      , metabolites_proportion_of_samples_below_cutoff = proportion_of_samples_below_cutoff_final
                                                                      , metabolite_id_column = metabolite_id_col
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
           }) 

#' Find Duplicate Feature IDs in MetaboliteAssayData Assays
#'
#' @description
#' Iterates through each assay tibble in a `MetaboliteAssayData` object
#' and identifies any duplicated values in the specified feature ID column. This
#' is useful for diagnosing data quality issues before attempting to resolve them.
#'
#' @param theObject A `MetaboliteAssayData` object.
#'
#' @return A named list. Each element corresponds to an assay in the input object.
#'   If duplicates are found in an assay, the element will be a tibble showing
#'   the duplicated identifiers and their counts. If no duplicates are found for an
#'   assay, the corresponding element in the list will be `NULL`.
#'
#' @importFrom dplyr count filter pull %>%
#' @importFrom purrr map set_names
#' @importFrom methods slot
#' @importFrom rlang sym .data
#'
#' @examples
#' \dontrun{
#' # Assuming 'met_assay_obj' is your MetaboliteAssayData object
#' duplicate_ids_list <- findDuplicateFeatureIDs(met_assay_obj)
#'
#' # The list contains results for all assays (NULL for those without duplicates)
#' print(duplicate_ids_list)
#'
#' # To get duplicates from the first assay (if any)
#' duplicates_assay1 <- duplicate_ids_list[[1]]
#' if (!is.null(duplicates_assay1)) {
#'   print(duplicates_assay1)
#' }
#' }
#' @export
findDuplicateFeatureIDs <- function(theObject) {
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
      warning(sprintf("Assay '%s': Feature ID column '%s' not found. Skipping duplicate check.",
                      assay_name, feature_id_col))
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

  if(all(sapply(duplicate_list, is.null))) {
      message("No duplicate feature IDs found in any assay.")
  }

  return(duplicate_list)
}



##-----------------------------------------------------------------------------
## Helper Function to Resolve Duplicate Features by Intensity
##-----------------------------------------------------------------------------

#' Resolve Duplicate Features by Highest Average Intensity (Helper)
#'
#' @description
#' A helper function that resolves duplicated feature IDs within a single assay tibble.
#' For each set of duplicate IDs, it keeps only the feature (row) with the highest
#' average intensity across all specified sample columns.
#'
#' @details
#' The function first calculates the mean intensity for each row across the provided
#' `sample_cols`. For rows that share the same identifier in `id_col`, it compares
#' their mean intensities and discards all but the one with the maximum value.
#' If there's a tie in average intensity, `slice_max` with `with_ties = FALSE`
#' ensures that only one row is arbitrarily kept.
#'
#' This is a common strategy to consolidate redundant feature measurements into a
#' single, most representative data point.
#'
#' @param assay_tibble A data frame or tibble representing one assay, with features
#'   as rows and samples as columns.
#' @param id_col A character string specifying the name of the column containing
#'   the feature identifiers to check for duplicates.
#' @param sample_cols A character vector specifying the names of the numeric columns
#'   that contain sample intensity data. These are used to calculate the average intensity.
#'
#' @return A tibble with duplicate features resolved. The returned tibble has the
#'   same structure as the input but with fewer rows if duplicates were found and removed.
#'
#' @keywords internal
#' @importFrom dplyr group_by summarise ungroup filter slice_max select rowwise mutate c_across any_of count
#' @importFrom rlang sym !! .data
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



#' Resolve Duplicate Feature IDs in MetaboliteAssayData
#'
#' @description
#' This method resolves duplicate feature identifiers within each assay of a
#' `MetaboliteAssayData` object. It uses different logic for Internal Standards (ITSDs)
#' versus regular (non-ITSD) features.
#'
#' @details
#' The method iterates through each assay and performs the following steps:
#'
#' 1.  **Separation of ITSDs:** It first identifies ITSD features using the regex
#'     pattern stored in the `internal_standard_regex` slot. The pattern is matched
#'     against the columns specified in `itsd_pattern_columns` (which defaults to the
#'     `annotation_id_column` slot).
#'
#' 2.  **ITSD Duplicate Resolution:**
#'     - For ITSDs, duplicates are identified based on a **descriptive name column**
#'       (assumed to be named `"metabolite"`).
#'     - If duplicate descriptive names are found (e.g., two ITSDs both named "Valine-d8"),
#'       they are ranked by their average intensity across all samples.
#'     - The primary identifier (in the `metabolite_id_column`) for these duplicated ITSDs
#'       is then **rewritten** to be unique by appending a rank suffix (e.g.,
#'       "Valine-d8_1", "Valine-d8_2"). This ensures that all ITSDs are retained but
#'       with unique primary IDs.
#'
#' 3.  **Non-ITSD Duplicate Resolution:**
#'     - For all other features, duplicates are identified based on the **primary ID column**
#'       (specified in the `metabolite_id_column` slot).
#'     - The function `resolveDuplicateFeaturesByIntensity` is used to resolve these
#'       duplicates: for each set of duplicate IDs, only the feature with the
#'       **highest average intensity** across all samples is kept. The others are discarded.
#'
#' 4.  **Recombination:** The processed ITSD and non-ITSD features are combined back
#'     into a single, resolved assay tibble.
#'
#' This dual-logic approach ensures that critical internal standards are not lost
#' due to ID clashes, while still cleaning up redundant measurements for other metabolites.
#'
#' @param theObject A `MetaboliteAssayData` object.
#' @param itsd_pattern_columns A character vector specifying which columns in the
#'   assay tibbles to search for the ITSD regex pattern. If `NULL` (default),
#'   it uses the column name stored in the `annotation_id_column` slot.
#'
#' @return An updated `MetaboliteAssayData` object with duplicate features resolved
#'   in its `metabolite_data` slot.
#'
#' @importFrom dplyr count filter pull %>% bind_rows group_by arrange desc mutate row_number ungroup slice_max select any_of rowwise c_across
#' @importFrom purrr map set_names reduce
#' @importFrom methods slot
#' @importFrom rlang sym !! .data :=
#' @importFrom stringr str_detect
#'
#' @describeIn resolveDuplicateFeatures Method for MetaboliteAssayData
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
                 message(sprintf("Checking for ITSD regex in user-specified columns: %s", paste(columns_to_check_for_itsd, collapse=", ")))
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
                      warning(sprintf("Specified ITSD pattern column(s) missing from assay '%s': %s. Cannot check these columns for ITSDs.",
                                      assay_name, paste(missing_itsd_check_cols, collapse=", ")), immediate. = TRUE)
                      # Adjust columns to check to only those that exist
                      columns_to_check_for_itsd <- columns_to_check_for_itsd[itsd_check_cols_exist]
                      if(length(columns_to_check_for_itsd) == 0) {
                          message("No valid columns remain for ITSD pattern checking.")
                      }
                 }

                missing_req_cols <- required_cols[!required_cols %in% colnames(assay_tibble)]
                if(length(missing_req_cols) > 0) {
                     warning(sprintf("Required column(s) missing from assay '%s': %s. Skipping assay.",
                                     assay_name, paste(missing_req_cols, collapse=", ")), immediate. = TRUE)
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
                         message(sprintf("Identified ITSDs using regex '%s' on column(s): %s.", itsd_regex, paste(columns_to_check_for_itsd, collapse=", ")))
                         is_itsd_present <- TRUE
                         itsd_rows <- assay_tibble_numeric %>% dplyr::filter(itsd_indices)
                         non_itsd_rows <- assay_tibble_numeric %>% dplyr::filter(!itsd_indices)
                         message(sprintf("Separated %d ITSD row(s) and %d non-ITSD row(s).", nrow(itsd_rows), nrow(non_itsd_rows)))
                    } else {
                         message(sprintf("No ITSDs found using regex '%s' on column(s): %s.", itsd_regex, paste(columns_to_check_for_itsd, collapse=", ")))
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
                    duplicate_itsd_specific_names <- itsd_specific_name_counts %>% dplyr::filter(.data$itsd_specific_count > 1) %>% dplyr::pull(!!rlang::sym(specific_name_col))

                    if (length(duplicate_itsd_specific_names) > 0) {
                        message(sprintf("Found %d duplicated descriptive ITSD name(s) in column '%s': %s",
                                        length(duplicate_itsd_specific_names), specific_name_col, paste(duplicate_itsd_specific_names, collapse=", ")))
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

