
' @title Log2 Transform Assay Data for MetaboliteAssayData
#'
#' @description
#' Applies a log2 transformation (log2(x + offset)) to the numeric sample columns
#' in all assay tibbles stored within the `metabolite_data` slot of a
#' `MetaboliteAssayData` object.
#'
#' @details
#' This transformation is often applied to reduce skewness and stabilize variance
#' in quantitative omics data, especially before visualization or statistical modeling.
#' Zeros in the data are handled by adding a small `offset` before transformation.
#' Non-numeric columns and the metabolite ID column are preserved.
#'
#' @describeIn logTransformAssays Method for MetaboliteAssayData
#'
#' @param theObject A `MetaboliteAssayData` object.
#' @param offset A small positive number to add before log transformation (default: 1).
#' @param ... Currently unused.
#'
#' @return An updated `MetaboliteAssayData` object with log2 transformed values in the `metabolite_data` slot.
#'
#' @importFrom dplyr mutate across all_of select filter intersect union setdiff
#' @importFrom rlang sym !!
#' @importFrom purrr map set_names
#' @importFrom methods slot slot<- is
#' @importFrom tibble as_tibble is_tibble
#'
#' @export
setMethod(f = "logTransformAssays",
          signature = "MetaboliteAssayData",
          definition = function(theObject, offset = 1, ...) {

            # --- Input Validation ---
            if (!is.numeric(offset) || length(offset) != 1 || offset <= 0) {
              stop("`offset` must be a single positive numeric value.")
            }
            message(sprintf("Applying log2(x + %s) transformation to assays.", as.character(offset)))

            # --- Get Object Slots ---
            assay_list <- methods::slot(theObject, "metabolite_data")
            metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
            design_matrix <- methods::slot(theObject, "design_matrix")
            sample_id_col_name <- methods::slot(theObject, "sample_id")

            if (length(assay_list) == 0) {
              warning("MetaboliteAssayData object has no assays in 'metabolite_data' slot. No transformation performed.")
              return(theObject)
            }

            # Ensure list is named
            original_assay_names <- names(assay_list)
             if (is.null(original_assay_names)) {
                 names(assay_list) <- paste0("Assay_", seq_along(assay_list))
                 warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).", immediate. = TRUE)
             } else if (any(original_assay_names == "")) {
                 needs_name <- which(original_assay_names == "")
                 original_assay_names[needs_name] <- paste0("Assay_", needs_name)
                 names(assay_list) <- original_assay_names
                 warning("Some assays were unnamed. Using default names for them.", immediate. = TRUE)
             }
             assay_names <- names(assay_list) # Use the potentially corrected names


            # --- Process Each Assay ---
            transformed_assay_list <- lapply(seq_along(assay_list), function(i) {
                assay_index_name <- assay_names[i] # Use the name from the corrected list
                assay_tibble <- assay_list[[i]]
                 message(sprintf("-- Processing assay: %s", assay_index_name))

                 if (!tibble::is_tibble(assay_tibble)) {
                     warning(sprintf("Assay '%s' is not a tibble. Attempting to coerce.", assay_index_name), immediate. = TRUE)
                     assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                         warning(sprintf("Failed to coerce assay '%s' to tibble: %s. Skipping transformation.", assay_index_name, e$message), immediate. = TRUE)
                         return(NULL)
                     })
                     if (is.null(assay_tibble)) return(assay_list[[i]]) # Return original if coercion failed
                 }

                 # Check for metabolite ID column
                 if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                     warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping transformation.", assay_index_name, metabolite_id_col_name), immediate. = TRUE)
                     return(assay_tibble) # Return original
                 }

                 # --- Identify Sample Columns ---
                 if (!methods::is(design_matrix, "data.frame")) {
                    stop("Slot 'design_matrix' is not a data.frame.")
                 }
                 if (!sample_id_col_name %in% colnames(design_matrix)) {
                     stop(sprintf("Sample ID column '%s' not found in design_matrix. Cannot identify sample columns for transformation.", sample_id_col_name))
                 }
                 design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e){
                     stop(sprintf("Could not extract sample IDs from design_matrix column '%s': %s", sample_id_col_name, e$message))
                 })
                 all_assay_cols <- colnames(assay_tibble)
                 sample_cols <- intersect(all_assay_cols, design_samples)
                 metadata_cols <- setdiff(all_assay_cols, sample_cols)
                 # Ensure metabolite ID is not treated as a sample column
                 metadata_cols <- union(metadata_cols, metabolite_id_col_name)
                 sample_cols <- setdiff(all_assay_cols, metadata_cols)
                 # ----------------------------- #

                 if (length(sample_cols) == 0) {
                     warning(sprintf("Assay '%s': No numeric sample columns identified matching design matrix sample IDs. Skipping transformation.", assay_index_name), immediate. = TRUE)
                     return(assay_tibble) # Return original
                 }

                 # --- Apply Log2 Transformation ---
                 transformed_tibble <- tryCatch({
                     assay_tibble %>%
                         # Replace negative values with 0 before log transformation
                         dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), ~ ifelse(!is.na(.x) & .x < 0, 0, .x))) %>%
                         # Ensure target columns are numeric before transformation
                         dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), as.numeric)) %>%
                         dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), ~ log2(.x + offset)))
                 }, error = function(e) {
                      warning(sprintf("Assay '%s': Error during log2 transformation: %s. Returning original assay data.", assay_index_name, e$message), immediate. = TRUE)
                      return(assay_tibble) # Return original on error
                 })

                 # Check if transformation actually happened (e.g., if error occurred)
                 if (identical(transformed_tibble, assay_tibble)) {
                    message(sprintf("Assay '%s': Transformation skipped or failed.", assay_index_name))
                 } else {
                    message(sprintf("Assay '%s': Successfully applied log2 transformation to %d sample column(s).", assay_index_name, length(sample_cols)))
                 }

                 return(transformed_tibble)
            })

            # Restore original names if they existed
            names(transformed_assay_list) <- assay_names

            # Assign the list of transformed assays back to the object
            methods::slot(theObject, "metabolite_data") <- transformed_assay_list

            # Optional: Add log transformation status to args (consider structure)
            # Ensure args is a list
            if (!is.list(theObject@args)) {
                warning("Slot 'args' is not a list. Cannot record log transformation status.", immediate. = TRUE)
            } else {
                theObject@args$log_transformed <- TRUE
                theObject@args$log_transform_offset <- offset
            }


            message("Log2 transformation complete for all assays.")
            return(theObject)
          }
)


#' @title Normalize by Internal Standard (ITSD) for MetaboliteAssayData
#'
#' @description
#' Normalizes the untransformed metabolite intensity data in each assay of a
#' `MetaboliteAssayData` object using Internal Standards (ITSDs).
#'
#' @details
#' This method corrects for systematic variations between samples (e.g., differences
#' in sample loading or instrument sensitivity) based on the signal of Internal
#' Standards (ITSDs). The process involves:
#' 1. Identifying ITSD features using a regex pattern (`internal_standard_regex` slot)
#'    applied to specified annotation columns (`itsd_pattern_columns`).
#' 2. Calculating an aggregate ITSD signal for *each sample* by combining the
#'    intensities of all identified ITSDs within that sample. The aggregation method
#'    is chosen via `itsd_aggregation` ("sum", "mean", or "median").
#' 3. Calculating the *average* of these aggregate ITSD signals across *all samples*.
#' 4. For each sample, calculating a normalization factor:
#'    `Factor = Average_Aggregate_ITSD / Sample_Aggregate_ITSD`.
#' 5. Multiplying the intensity of every non-ITSD feature in that sample by this
#'    sample-specific `Factor`.
#'
#' This method operates on the raw (or previously normalized, but typically not
#' log-transformed) data in the `metabolite_data` slot. It aims to preserve the
#' overall scale of the data while adjusting relative differences between samples,
#' centering the normalized data around the average ITSD response observed.
#'
#' Assays where no ITSDs are found will issue a warning and remain unchanged. Samples
#' where the aggregate ITSD signal is zero or NA will also result in warnings,
#' and the normalized values for those samples will likely become NA.
#'
#' @describeIn normaliseUntransformedData Method for MetaboliteAssayData using ITSD
#'
#' @param theObject A `MetaboliteAssayData` object.
#' @param method Must be "ITSD" for this method.
#' @param itsd_pattern_columns A character vector specifying the column names within
#'   each assay tibble where the ITSD pattern should be searched. If `NULL` (default),
#'   it uses the column name stored in the `annotation_id_column` slot of the object.
#' @param itsd_aggregation The method used to aggregate intensities of multiple ITSDs
#'   within each sample. Options are "sum" (default),
#'   "mean", or "median".
#' @param remove_itsd_after_norm Logical (default: TRUE). If TRUE, rows identified as
#'   ITSDs are removed from the assay tibbles after normalization is complete.
#' @param ... Currently unused.
#'
#' @return An updated `MetaboliteAssayData` object with ITSD-normalized values in the
#'   `metabolite_data` slot. Normalization status is recorded in the `args` slot.
#'
#' @importFrom dplyr mutate across all_of filter select matches summarise group_by pull if_else
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom rlang sym !! := .data
#' @importFrom purrr map set_names reduce
#' @importFrom methods slot slot<- is
#' @importFrom tibble as_tibble is_tibble tibble
#' @importFrom stringr str_detect
#'
#' @export
setMethod(f = "normaliseUntransformedData",
          signature = signature(theObject = "MetaboliteAssayData", method = "character"),
          definition = function(theObject,
                                method = "ITSD",
                                itsd_pattern_columns = NULL, # Default set later
                                itsd_aggregation = "sum",
                                remove_itsd_after_norm = TRUE,
                                ...) {

            # --- Input Validation ---
            if (tolower(method) != "itsd") {
                stop("This method currently only supports method = 'ITSD'.")
            }
            valid_aggregations <- c("sum", "mean", "median")
            if (!tolower(itsd_aggregation) %in% valid_aggregations) {
                stop(sprintf("`itsd_aggregation` must be one of: %s", paste(valid_aggregations, collapse = ", ")))
            }
            itsd_aggregation_func <- switch(tolower(itsd_aggregation),
                                            "sum" = sum,
                                            "mean" = mean,
                                            "median" = stats::median)
            message(sprintf("Applying Internal Standard (ITSD) normalization using '%s' aggregation and scaling to average ITSD response.", tolower(itsd_aggregation)))

            # --- Get Object Slots ---
            assay_list <- methods::slot(theObject, "metabolite_data")
            metabolite_id_col <- methods::slot(theObject, "metabolite_id_column")
            design_matrix <- methods::slot(theObject, "design_matrix")
            sample_id_col <- methods::slot(theObject, "sample_id")
            itsd_regex <- methods::slot(theObject, "internal_standard_regex")
            annotation_col <- methods::slot(theObject, "annotation_id_column") # Default ITSD column

            # Set default for itsd_pattern_columns if NULL
            if (is.null(itsd_pattern_columns)) {
                itsd_pattern_columns <- annotation_col
                message(sprintf("Using default annotation column '%s' to identify ITSDs.", annotation_col))
            } else {
                 message(sprintf("Using column(s) '%s' to identify ITSDs.", paste(itsd_pattern_columns, collapse="', '")))
            }

            if (length(assay_list) == 0) {
              warning("MetaboliteAssayData object has no assays in 'metabolite_data' slot. No normalization performed.")
              return(theObject)
            }
             if (is.null(itsd_regex) || itsd_regex == "") {
                 stop("The `internal_standard_regex` slot is empty. Cannot identify ITSDs.")
             }


            # Ensure list is named
            original_assay_names <- names(assay_list)
             if (is.null(original_assay_names)) {
                 names(assay_list) <- paste0("Assay_", seq_along(assay_list))
                 warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).", immediate. = TRUE)
             } else if (any(original_assay_names == "")) {
                 needs_name <- which(original_assay_names == "")
                 original_assay_names[needs_name] <- paste0("Assay_", needs_name)
                 names(assay_list) <- original_assay_names
                 warning("Some assays were unnamed. Using default names for them.", immediate. = TRUE)
             }
             assay_names <- names(assay_list) # Use the potentially corrected names

            # --- Identify Sample Columns ---
            if (!methods::is(design_matrix, "data.frame")) {
              stop("Slot 'design_matrix' is not a data.frame.")
            }
            if (!sample_id_col %in% colnames(design_matrix)) {
               stop(sprintf("Sample ID column '%s' not found in design_matrix. Cannot identify sample columns for normalization.", sample_id_col))
            }
            design_samples <- tryCatch(as.character(design_matrix[[sample_id_col]]), error = function(e){
               stop(sprintf("Could not extract sample IDs from design_matrix column '%s': %s", sample_id_col, e$message))
            })
            if (length(design_samples) == 0) {
                 stop("No sample IDs found in the design matrix.")
            }


            # --- Process Each Assay ---
            normalized_assay_list <- lapply(seq_along(assay_list), function(i) {
                assay_index_name <- assay_names[i] # Use the name from the corrected list
                assay_tibble <- assay_list[[i]]
                 message(sprintf("-- Processing assay: %s", assay_index_name))

                 # --- Basic Checks ---
                 if (!tibble::is_tibble(assay_tibble)) {
                     warning(sprintf("Assay '%s' is not a tibble. Attempting to coerce.", assay_index_name), immediate. = TRUE)
                     assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                         warning(sprintf("Failed to coerce assay '%s' to tibble: %s. Skipping normalization.", assay_index_name, e$message), immediate. = TRUE)
                         return(NULL) # Signal to skip this assay
                     })
                     if (is.null(assay_tibble)) return(assay_list[[i]]) # Return original if coercion failed
                 }
                if (!metabolite_id_col %in% colnames(assay_tibble)) {
                    warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping normalization.", assay_index_name, metabolite_id_col), immediate. = TRUE)
                    return(assay_tibble)
                }
                 # Check if *any* of the itsd_pattern_columns exist
                 if (!any(itsd_pattern_columns %in% colnames(assay_tibble))) {
                    warning(sprintf("Assay '%s': None of the specified ITSD pattern columns ('%s') found. Cannot identify ITSDs. Skipping normalization.",
                                    assay_index_name, paste(itsd_pattern_columns, collapse="', '")), immediate. = TRUE)
                    return(assay_tibble)
                 }
                 # Filter itsd_pattern_columns to only those present in the current assay
                 actual_itsd_cols <- intersect(itsd_pattern_columns, colnames(assay_tibble))
                  if (length(actual_itsd_cols) == 0) { # Should be caught above, but double check
                      warning(sprintf("Assay '%s': No ITSD identification columns found after checking existence. Skipping normalization.", assay_index_name), immediate. = TRUE)
                      return(assay_tibble)
                  }


                # --- Identify Sample Columns in this Assay ---
                 all_assay_cols <- colnames(assay_tibble)
                 sample_cols <- intersect(all_assay_cols, design_samples)
                 if (length(sample_cols) == 0) {
                     warning(sprintf("Assay '%s': No sample columns identified matching design matrix sample IDs. Skipping normalization.", assay_index_name), immediate. = TRUE)
                     return(assay_tibble)
                 }
                 # Ensure sample columns are numeric
                 non_numeric_samples <- sample_cols[!sapply(assay_tibble[sample_cols], is.numeric)]
                 if (length(non_numeric_samples) > 0) {
                    warning(sprintf("Assay '%s': Non-numeric sample columns found: %s. Attempting coercion, but this may indicate upstream issues.",
                                    assay_index_name, paste(non_numeric_samples, collapse=", ")), immediate. = TRUE)
                    assay_tibble <- assay_tibble |>
                        dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
                 }

                 # --- Identify ITSD Rows ---
                 itsd_rows_logical <- assay_tibble |>
                    dplyr::select(dplyr::all_of(actual_itsd_cols)) |>
                    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringr::str_detect(as.character(.), itsd_regex))) |>
                    # Row is ITSD if pattern matches in *any* of the specified columns
                    purrr::reduce(`|`)

                 if (!any(itsd_rows_logical, na.rm = TRUE)) {
                     warning(sprintf("Assay '%s': No ITSD features identified using regex '%s' in columns '%s'. Skipping normalization.",
                                     assay_index_name, itsd_regex, paste(actual_itsd_cols, collapse="', '")), immediate. = TRUE)
                     return(assay_tibble)
                 }
                 itsd_data <- assay_tibble |> dplyr::filter(itsd_rows_logical)
                 non_itsd_data <- assay_tibble |> dplyr::filter(!itsd_rows_logical)
                 n_itsd <- nrow(itsd_data)
                 message(sprintf("   Identified %d ITSD features.", n_itsd))

                # --- Calculate Normalization Factors ---
                # Step 1: Calculate aggregate ITSD per sample
                norm_factors_long <- tryCatch({
                    itsd_data |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col, sample_cols))) |>
                        tidyr::pivot_longer(cols = dplyr::all_of(sample_cols), names_to = "Sample", values_to = "Intensity") |>
                        dplyr::group_by(.data$Sample) |>
                        # Aggregate ITSD intensities per sample
                        # Use SampleNormFactor to be clear about the value's meaning
                        dplyr::summarise(SampleNormFactor = itsd_aggregation_func(.data$Intensity, na.rm = TRUE), .groups = "drop")
                }, error = function(e) {
                     warning(sprintf("Assay '%s': Error calculating normalization factors: %s. Skipping normalization.", assay_index_name, e$message), immediate. = TRUE)
                     return(NULL)
                 })
                if (is.null(norm_factors_long)) return(assay_tibble) # Return original on error

                # Check for zero or NA factors
                problematic_factors <- norm_factors_long |>
                    dplyr::filter(is.na(.data$SampleNormFactor) | .data$SampleNormFactor == 0)
                if (nrow(problematic_factors) > 0) {
                    warning(sprintf("Assay '%s': Aggregate ITSD signal (SampleNormFactor) is NA or zero for samples: %s. Normalized values will be NA for these samples.",
                                    assay_index_name, paste(problematic_factors$Sample, collapse=", ")), immediate. = TRUE)
                    # Set factor to NA to ensure division results in NA
                    norm_factors_long <- norm_factors_long |>
                        dplyr::mutate(SampleNormFactor = dplyr::if_else(is.na(.data$SampleNormFactor) | .data$SampleNormFactor == 0, NA_real_, .data$SampleNormFactor))
                }

                # Step 2: Calculate the average aggregate ITSD signal across all samples
                average_norm_factor <- mean(norm_factors_long$SampleNormFactor, na.rm = TRUE)
                 if (is.na(average_norm_factor) || average_norm_factor == 0) {
                     warning(sprintf("Assay '%s': Average aggregate ITSD signal is NA or zero. Cannot perform average-centered normalization. Skipping.", assay_index_name), immediate. = TRUE)
                     return(assay_tibble)
                 }

                # Create a named vector for easy lookup (Sample -> NormFactor)
                sample_norm_factors_vec <- stats::setNames(norm_factors_long$SampleNormFactor, norm_factors_long$Sample)

                # --- Apply Normalization ---
                 message(sprintf("   Applying normalization to %d non-ITSD features...", nrow(non_itsd_data)))
                 normalized_non_itsd_data <- tryCatch({
                    non_itsd_data |>
                        dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols),
                                             # Multiply intensity by (Average Factor / Sample Factor)
                                             ~ .x * (average_norm_factor / sample_norm_factors_vec[dplyr::cur_column()])))
                 }, error = function(e) {
                     warning(sprintf("Assay '%s': Error applying normalization factors: %s. Returning unnormalized data for non-ITSD features.", assay_index_name, e$message), immediate. = TRUE)
                     return(non_itsd_data) # Return unnormalized if error
                 })

                # --- Reconstruct Assay Tibble ---
                if (remove_itsd_after_norm) {
                    final_assay_tibble <- normalized_non_itsd_data
                    message(sprintf("   Removed %d ITSD features after normalization.", n_itsd))
                } else {
                    # If keeping ITSDs, they remain unnormalized (or could be normalized like others)
                    # Here, we keep them unnormalized as their purpose was the factor calculation.
                    final_assay_tibble <- dplyr::bind_rows(normalized_non_itsd_data, itsd_data) |>
                        # Optional: arrange back by original ID or similar if needed
                         dplyr::arrange(!!rlang::sym(metabolite_id_col))
                    message(sprintf("   Kept %d ITSD features (unnormalized) in the output.", n_itsd))
                }

                 message(sprintf("   Assay '%s' normalization complete.", assay_index_name))
                 return(final_assay_tibble)
            }) # End lapply

            # Filter out any assays that failed (returned NULL - though currently returning original)
            # This check might be redundant if errors return original, but good practice.
            successful_assays <- !sapply(normalized_assay_list, is.null)
             if (!all(successful_assays)) {
                 warning("Normalization failed or was skipped for some assays.", immediate. = TRUE)
                 # Keep original data for failed assays
                 normalized_assay_list[!successful_assays] <- assay_list[!successful_assays]
             }


            # Restore original names if they existed
            names(normalized_assay_list) <- assay_names

            # Assign the list of normalized assays back to the object
            methods::slot(theObject, "metabolite_data") <- normalized_assay_list

            # --- Update Args Slot ---
            if (!is.list(theObject@args)) {
                warning("Slot 'args' is not a list. Cannot record ITSD normalization status.", immediate. = TRUE)
            } else {
                 # Ensure the specific list exists
                if (!"ITSDNormalization" %in% names(theObject@args)) {
                    theObject@args$ITSDNormalization <- list()
                }
                theObject@args$ITSDNormalization$applied <- TRUE
                theObject@args$ITSDNormalization$method_type <- "average_centered"
                theObject@args$ITSDNormalization$itsd_aggregation <- tolower(itsd_aggregation)
                theObject@args$ITSDNormalization$itsd_pattern_columns <- itsd_pattern_columns # Record actual columns used (potentially default)
                theObject@args$ITSDNormalization$removed_itsd <- remove_itsd_after_norm
                theObject@args$ITSDNormalization$timestamp <- Sys.time()
            }

            message("ITSD normalization process complete for all applicable assays.")
            return(theObject)
          }
)

# Optional: Add other normalization methods (e.g., PQN, Median) here later
# setMethod(f = "normaliseUntransformedData",
#           signature = signature(theObject = "MetaboliteAssayData", method = "character"),
#           definition = function(theObject, method = "PQN", ...) { ... }
# ) 