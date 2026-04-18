#' @title Log2 Transform Assay Data for MetaboliteAssayData
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
#' @title Log Transform Assays for MetaboliteAssayData
#' @name logTransformAssays,MetaboliteAssayData-method
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
#' @export
setMethod(
    f = "logTransformAssays",
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
                if (is.null(assay_tibble)) {
                    return(assay_list[[i]])
                } # Return original if coercion failed
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
            design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) {
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
            transformed_tibble <- tryCatch(
                {
                    assay_tibble %>%
                        # Replace negative values with 0 before log transformation
                        dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), ~ ifelse(!is.na(.x) & .x < 0, 0, .x))) %>%
                        # Ensure target columns are numeric before transformation
                        dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), as.numeric)) %>%
                        dplyr::mutate(dplyr::across(dplyr::all_of(sample_cols), ~ log2(.x + offset)))
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error during log2 transformation: %s. Returning original assay data.", assay_index_name, e$message), immediate. = TRUE)
                    return(assay_tibble) # Return original on error
                }
            )

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
#' @param itsd_feature_ids Optional named list of feature ID vectors for manual ITSD
#'   selection per assay. Names must match assay names. When provided, overrides
#'   regex-based ITSD detection. Example: `list(LCMS_Pos = c("ITSD_1", "ITSD_2"), LCMS_Neg = c("ITSD_3"))`.
#' @param itsd_pattern_columns A character vector specifying the column names within
#'   each assay tibble where the ITSD pattern should be searched. If `NULL` (default),
#'   it uses the column name stored in the `annotation_id_column` slot of the object.
#'   Ignored if `itsd_feature_ids` is provided.
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
setMethod(
    f = "normaliseUntransformedData",
    signature = signature(theObject = "MetaboliteAssayData", method = "character"),
    definition = function(theObject,
                          method = "ITSD",
                          itsd_feature_ids = NULL,
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
            "median" = stats::median
        )

        # Determine ITSD selection mode
        use_manual_itsd <- !is.null(itsd_feature_ids) && is.list(itsd_feature_ids) && length(itsd_feature_ids) > 0
        if (use_manual_itsd) {
            message(sprintf("Applying Internal Standard (ITSD) normalization using MANUAL feature selection and '%s' aggregation.", tolower(itsd_aggregation)))
        } else {
            message(sprintf("Applying Internal Standard (ITSD) normalization using REGEX pattern and '%s' aggregation.", tolower(itsd_aggregation)))
        }

        # --- Get Object Slots ---
        assay_list <- methods::slot(theObject, "metabolite_data")
        metabolite_id_col <- methods::slot(theObject, "metabolite_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col <- methods::slot(theObject, "sample_id")
        itsd_regex <- methods::slot(theObject, "internal_standard_regex")
        annotation_col <- methods::slot(theObject, "annotation_id_column") # Default ITSD column

        # Set default for itsd_pattern_columns if NULL (only relevant for regex mode)
        if (!use_manual_itsd) {
            if (is.null(itsd_pattern_columns)) {
                itsd_pattern_columns <- annotation_col
                message(sprintf("Using default annotation column '%s' to identify ITSDs.", annotation_col))
            } else {
                message(sprintf("Using column(s) '%s' to identify ITSDs.", paste(itsd_pattern_columns, collapse = "', '")))
            }

            # Validate regex is available for regex mode
            if (is.null(itsd_regex) || itsd_regex == "") {
                stop("The `internal_standard_regex` slot is empty and no manual itsd_feature_ids provided. Cannot identify ITSDs.")
            }
        }

        if (length(assay_list) == 0) {
            warning("MetaboliteAssayData object has no assays in 'metabolite_data' slot. No normalization performed.")
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

        # --- Identify Sample Columns ---
        if (!methods::is(design_matrix, "data.frame")) {
            stop("Slot 'design_matrix' is not a data.frame.")
        }
        if (!sample_id_col %in% colnames(design_matrix)) {
            stop(sprintf("Sample ID column '%s' not found in design_matrix. Cannot identify sample columns for normalization.", sample_id_col))
        }
        design_samples <- tryCatch(as.character(design_matrix[[sample_id_col]]), error = function(e) {
            stop(sprintf("Could not extract sample IDs from design_matrix column '%s': %s", sample_id_col, e$message))
        })
        if (length(design_samples) == 0) {
            stop("No sample IDs found in the design matrix.")
        }

        # --- Initialize ITSD Feature Collector ---
        # Use environment for side-effect collection during lapply
        itsd_collector <- new.env(parent = emptyenv())
        itsd_collector$features_per_assay <- list()
        itsd_collector$counts_per_assay <- list()

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
                if (is.null(assay_tibble)) {
                    return(assay_list[[i]])
                } # Return original if coercion failed
            }
            if (!metabolite_id_col %in% colnames(assay_tibble)) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping normalization.", assay_index_name, metabolite_id_col), immediate. = TRUE)
                return(assay_tibble)
            }
            # --- Validate ITSD Identification Prerequisites ---
            # Only check pattern columns if using regex mode
            actual_itsd_cols <- NULL
            if (!use_manual_itsd) {
                # Check if *any* of the itsd_pattern_columns exist
                if (!any(itsd_pattern_columns %in% colnames(assay_tibble))) {
                    warning(sprintf(
                        "Assay '%s': None of the specified ITSD pattern columns ('%s') found. Cannot identify ITSDs. Skipping normalization.",
                        assay_index_name, paste(itsd_pattern_columns, collapse = "', '")
                    ), immediate. = TRUE)
                    return(assay_tibble)
                }
                # Filter itsd_pattern_columns to only those present in the current assay
                actual_itsd_cols <- intersect(itsd_pattern_columns, colnames(assay_tibble))
                if (length(actual_itsd_cols) == 0) { # Should be caught above, but double check
                    warning(sprintf("Assay '%s': No ITSD identification columns found after checking existence. Skipping normalization.", assay_index_name), immediate. = TRUE)
                    return(assay_tibble)
                }
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
                warning(sprintf(
                    "Assay '%s': Non-numeric sample columns found: %s. Attempting coercion, but this may indicate upstream issues.",
                    assay_index_name, paste(non_numeric_samples, collapse = ", ")
                ), immediate. = TRUE)
                assay_tibble <- assay_tibble |>
                    dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # --- Identify ITSD Rows ---
            if (use_manual_itsd) {
                # Manual mode: use provided feature IDs for this assay
                assay_itsd_ids <- itsd_feature_ids[[assay_index_name]]
                if (is.null(assay_itsd_ids) || length(assay_itsd_ids) == 0) {
                    warning(sprintf("Assay '%s': No manual ITSD features provided in itsd_feature_ids. Skipping normalization.", assay_index_name), immediate. = TRUE)
                    return(assay_tibble)
                }
                # Match feature IDs against metabolite_id_col
                feature_ids <- as.character(assay_tibble[[metabolite_id_col]])
                itsd_rows_logical <- feature_ids %in% as.character(assay_itsd_ids)

                if (!any(itsd_rows_logical, na.rm = TRUE)) {
                    warning(sprintf(
                        "Assay '%s': None of the %d manual ITSD feature IDs matched rows in the data. Skipping normalization.",
                        assay_index_name, length(assay_itsd_ids)
                    ), immediate. = TRUE)
                    return(assay_tibble)
                }
                message(sprintf("   Using %d manually selected ITSD features.", sum(itsd_rows_logical)))
            } else {
                # Regex mode: identify ITSDs by pattern matching
                itsd_rows_logical <- assay_tibble |>
                    dplyr::select(dplyr::all_of(actual_itsd_cols)) |>
                    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringr::str_detect(as.character(.), itsd_regex))) |>
                    # Row is ITSD if pattern matches in *any* of the specified columns
                    purrr::reduce(`|`)

                if (!any(itsd_rows_logical, na.rm = TRUE)) {
                    warning(sprintf(
                        "Assay '%s': No ITSD features identified using regex '%s' in columns '%s'. Skipping normalization.",
                        assay_index_name, itsd_regex, paste(actual_itsd_cols, collapse = "', '")
                    ), immediate. = TRUE)
                    return(assay_tibble)
                }
                message(sprintf("   Identified %d ITSD features via regex.", sum(itsd_rows_logical)))
            }

            itsd_data <- assay_tibble |> dplyr::filter(itsd_rows_logical)
            non_itsd_data <- assay_tibble |> dplyr::filter(!itsd_rows_logical)
            n_itsd <- nrow(itsd_data)
            message(sprintf("   Using %d ITSD features for normalization.", n_itsd))

            # --- Capture ITSD feature names for reporting ---
            itsd_feature_names <- as.character(itsd_data[[metabolite_id_col]])
            itsd_collector$features_per_assay[[assay_index_name]] <- itsd_feature_names
            itsd_collector$counts_per_assay[[assay_index_name]] <- n_itsd

            # --- Calculate Normalization Factors ---
            # Step 1: Calculate aggregate ITSD per sample
            norm_factors_long <- tryCatch(
                {
                    itsd_data |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col, sample_cols))) |>
                        tidyr::pivot_longer(cols = dplyr::all_of(sample_cols), names_to = "Sample", values_to = "Intensity") |>
                        dplyr::group_by(.data$Sample) |>
                        # Aggregate ITSD intensities per sample
                        # Use SampleNormFactor to be clear about the value's meaning
                        dplyr::summarise(SampleNormFactor = itsd_aggregation_func(.data$Intensity, na.rm = TRUE), .groups = "drop")
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error calculating normalization factors: %s. Skipping normalization.", assay_index_name, e$message), immediate. = TRUE)
                    return(NULL)
                }
            )
            if (is.null(norm_factors_long)) {
                return(assay_tibble)
            } # Return original on error

            # Check for zero or NA factors
            problematic_factors <- norm_factors_long |>
                dplyr::filter(is.na(.data$SampleNormFactor) | .data$SampleNormFactor == 0)
            if (nrow(problematic_factors) > 0) {
                warning(sprintf(
                    "Assay '%s': Aggregate ITSD signal (SampleNormFactor) is NA or zero for samples: %s. Normalized values will be NA for these samples.",
                    assay_index_name, paste(problematic_factors$Sample, collapse = ", ")
                ), immediate. = TRUE)
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
            normalized_non_itsd_data <- tryCatch(
                {
                    non_itsd_data |>
                        dplyr::mutate(dplyr::across(
                            dplyr::all_of(sample_cols),
                            # Multiply intensity by (Average Factor / Sample Factor)
                            ~ .x * (average_norm_factor / sample_norm_factors_vec[dplyr::cur_column()])
                        ))
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error applying normalization factors: %s. Returning unnormalized data for non-ITSD features.", assay_index_name, e$message), immediate. = TRUE)
                    return(non_itsd_data) # Return unnormalized if error
                }
            )

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

            # Store per-assay ITSD feature information for reporting
            if (length(itsd_collector$features_per_assay) > 0) {
                theObject@args$ITSDNormalization$itsd_features_per_assay <- as.list(itsd_collector$features_per_assay)
                theObject@args$ITSDNormalization$itsd_counts_per_assay <- as.list(itsd_collector$counts_per_assay)
                message(sprintf("   Stored ITSD feature names for %d assays", length(itsd_collector$features_per_assay)))
            }
        }

        message("ITSD normalization process complete for all applicable assays.")
        return(theObject)
    }
)

#' @title Normalize Between Samples for MetaboliteAssayData
#' @name normaliseBetweenSamples,MetaboliteAssayData-method
#' @param theObject Object of class MetaboliteAssayData
#' @param normalisation_method Method to use for normalization. Options are
#'   "cyclicloess", "quantile", "scale", "none". If NULL, the value is retrieved
#'   from the object's configuration arguments (looking for "normalisation_method").
#'
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom dplyr select all_of left_join relocate any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom methods slot slot<- is
#' @importFrom logger log_info
#' @export
setMethod(
    f = "normaliseBetweenSamples",
    signature = "MetaboliteAssayData",
    definition = function(theObject, normalisation_method = NULL) {
        assay_list <- methods::slot(theObject, "metabolite_data")
        metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col_name <- methods::slot(theObject, "sample_id")

        # --- Get Normalization Method ---
        # Use the general parameter name as in protein version for consistency
        normalisation_method_final <- checkParamsObjectFunctionSimplify(
            theObject,
            "normalisation_method",
            default = "cyclicloess" # Default if not found in args or user override
        )
        # Store the *actually used* method back into args
        theObject@args$normalisation_method <- normalisation_method_final
        log_info("Applying between-sample normalization method: {normalisation_method_final}")


        if (length(assay_list) == 0) {
            warning("No assays found in `metabolite_data` slot. Skipping normalization.")
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
        assay_names <- names(assay_list) # Use potentially corrected names


        # --- Process Each Assay ---
        normalized_assay_list <- lapply(seq_along(assay_list), function(i) {
            assay_name <- assay_names[i]
            assay_tibble <- assay_list[[i]]
            message(sprintf("-- Processing assay for normalization: %s", assay_name))

            # --- Basic Checks ---
            if (!tibble::is_tibble(assay_tibble)) {
                warning(sprintf("Assay '%s' is not a tibble. Attempting coercion.", assay_name), immediate. = TRUE)
                assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                    warning(sprintf("Failed to coerce assay '%s' to tibble: %s. Skipping normalization.", assay_name, e$message), immediate. = TRUE)
                    return(NULL) # Signal to skip
                })
                if (is.null(assay_tibble)) {
                    return(assay_list[[i]])
                } # Return original if coercion failed
            }
            if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping normalization.", assay_name, metabolite_id_col_name), immediate. = TRUE)
                return(assay_tibble)
            }

            # --- Identify Sample Columns ---
            design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) {
                character(0)
            })
            if (length(design_samples) == 0) {
                warning(sprintf("Assay '%s': Could not extract valid sample IDs from design matrix column '%s'. Skipping normalization.", assay_name, sample_id_col_name), immediate. = TRUE)
                return(assay_tibble)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) == 0) {
                warning(sprintf("Assay '%s': No sample columns identified matching design matrix sample IDs. Skipping normalization.", assay_name), immediate. = TRUE)
                return(assay_tibble)
            }
            # Ensure sample columns are numeric
            non_numeric_samples <- sample_cols[!sapply(assay_tibble[sample_cols], is.numeric)]
            if (length(non_numeric_samples) > 0) {
                warning(sprintf("Assay '%s': Non-numeric sample columns found: %s. Attempting coercion, but this may indicate upstream issues.", assay_name, paste(non_numeric_samples, collapse = ", ")), immediate. = TRUE)
                assay_tibble <- assay_tibble |>
                    dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # --- Prepare Matrix for Normalization ---
            assay_matrix <- tryCatch(
                {
                    assay_tibble |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |> # Select ID + Samples
                        tibble::column_to_rownames(var = metabolite_id_col_name) |>
                        as.matrix()
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error converting tibble to matrix: %s. Skipping normalization.", assay_name, e$message), immediate. = TRUE)
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(assay_tibble)
            } # Return original if matrix conversion failed

            assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf

            # Check if matrix is valid for normalization
            if (nrow(assay_matrix) < 1 || ncol(assay_matrix) < 1) {
                warning(sprintf("Assay '%s': Matrix is empty after preparation. Skipping normalization.", assay_name), immediate. = TRUE)
                return(assay_tibble)
            }
            # Check if all values are NA in any column (causes issues for some methods)
            if (any(colSums(!is.na(assay_matrix)) == 0)) {
                warning(sprintf("Assay '%s': At least one sample column contains only NA values. Skipping normalization.", assay_name), immediate. = TRUE)
                return(assay_tibble)
            }

            # --- Apply Normalization Method ---
            normalized_matrix <- assay_matrix # Default to original if method is 'none' or fails

            if (normalisation_method_final != "none") {
                normalized_matrix <- tryCatch(
                    {
                        switch(normalisation_method_final,
                            cyclicloess = {
                                message("   Applying cyclic loess normalization...")
                                limma::normalizeCyclicLoess(assay_matrix)
                            },
                            quantile = {
                                message("   Applying quantile normalization...")
                                limma::normalizeQuantiles(assay_matrix)
                            },
                            scale = {
                                message("   Applying scale (median absolute deviation) normalization...")
                                limma::normalizeMedianAbsValues(assay_matrix)
                            },
                            { # Default case for switch if none of the above match (should not happen due to checkParams)
                                warning(sprintf("Assay '%s': Unknown normalization method '%s'. Skipping normalization.", assay_name, normalisation_method_final), immediate. = TRUE)
                                assay_matrix # Return original matrix
                            }
                        )
                    },
                    error = function(e) {
                        warning(sprintf("Assay '%s': Error during '%s' normalization: %s. Returning unnormalized data for this assay.", assay_name, normalisation_method_final, e$message), immediate. = TRUE)
                        return(assay_matrix) # Return original matrix on error
                    }
                )
            } else {
                message("   Normalization method is 'none'. Skipping application.")
            }

            normalized_matrix[!is.finite(normalized_matrix)] <- NA # Ensure NAs remain NAs

            # --- Reconstruct Tibble ---
            # Get original metadata columns
            metadata_cols <- setdiff(colnames(assay_tibble), sample_cols)

            reconstructed_tibble <- tryCatch(
                {
                    normalized_data_tibble <- normalized_matrix |>
                        as.data.frame() |> # Convert matrix to data frame
                        tibble::rownames_to_column(var = metabolite_id_col_name) |>
                        tibble::as_tibble()

                    original_metadata_tibble <- assay_tibble |>
                        dplyr::select(dplyr::any_of(metadata_cols)) # Use any_of in case some metadata cols were dynamic

                    # Ensure join column types match (rownames_to_column creates character)
                    original_metadata_tibble_char <- original_metadata_tibble |>
                        dplyr::mutate(!!rlang::sym(metabolite_id_col_name) := as.character(!!rlang::sym(metabolite_id_col_name)))
                    normalized_data_tibble_char <- normalized_data_tibble |>
                        dplyr::mutate(!!rlang::sym(metabolite_id_col_name) := as.character(!!rlang::sym(metabolite_id_col_name)))

                    # Join normalized data with original metadata
                    dplyr::left_join(original_metadata_tibble_char, normalized_data_tibble_char, by = metabolite_id_col_name) |>
                        # Ensure original column order (metadata first, then samples in original order)
                        dplyr::relocate(dplyr::all_of(metadata_cols), dplyr::all_of(sample_cols))
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error reconstructing tibble after normalization: %s. Returning original data.", assay_name, e$message), immediate. = TRUE)
                    return(assay_tibble) # Return original on error
                }
            )

            message(sprintf("   Assay '%s' normalization complete.", assay_name))
            return(reconstructed_tibble)
        })

        # Restore original names
        names(normalized_assay_list) <- assay_names

        # Update the slot in the object
        methods::slot(theObject, "metabolite_data") <- normalized_assay_list

        # --- Clean Design Matrix (as done in protein version) ---
        # Ensure the cleanDesignMatrix method exists for MetaboliteAssayData
        # (From handover.md, this should exist)
        theObject <- tryCatch(
            {
                cleanDesignMatrix(theObject)
            },
            error = function(e) {
                warning(sprintf("Error running cleanDesignMatrix after normalization: %s. Design matrix might not be fully synchronized.", e$message))
                return(theObject) # Return object even if cleaning fails
            }
        )

        log_info("Between-sample normalization process finished for all assays.")
        return(theObject)
    }
)

#' @title Clean Design Matrix for MetaboliteAssayData
#' @name cleanDesignMatrix,MetaboliteAssayData-method
#' @importFrom dplyr inner_join select rename filter all_of any_of
#' @importFrom rlang sym !!
#' @importFrom methods slot
#' @importFrom tibble tibble
#' @export
setMethod(
    f = "cleanDesignMatrix",
    signature = "MetaboliteAssayData",
    definition = function(theObject) {
        assay_list <- methods::slot(theObject, "metabolite_data")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col_name <- methods::slot(theObject, "sample_id")
        metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column") # Needed to exclude from sample cols

        if (length(assay_list) == 0) {
            warning("cleanDesignMatrix: No assays found in `metabolite_data`. Returning object unchanged.")
            return(theObject)
        }

        # Assume sample columns are consistent across assays (enforced by validity)
        # Get sample columns from the first assay
        first_assay <- assay_list[[1]]

        # --- Identify Sample Columns in the Assay --- #
        design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) {
            character(0)
        })
        if (length(design_samples) == 0) {
            warning(sprintf("cleanDesignMatrix: Could not extract valid sample IDs from design matrix column '%s'. Returning object unchanged.", sample_id_col_name), immediate. = TRUE)
            return(theObject)
        }
        all_assay_cols <- colnames(first_assay)
        sample_cols_in_assay <- intersect(all_assay_cols, design_samples)
        if (length(sample_cols_in_assay) == 0) {
            warning("cleanDesignMatrix: No sample columns identified in the first assay matching design matrix sample IDs. Returning object unchanged.")
            return(theObject)
        }
        # Ensure columns are treated as character for join consistency
        sample_cols_vector <- as.character(sample_cols_in_assay)

        # --- Filter and Reorder Design Matrix --- #
        # Ensure the sample ID column in the original design matrix is character for join
        design_matrix_char_id <- design_matrix |>
            dplyr::mutate(!!rlang::sym(sample_id_col_name) := as.character(!!rlang::sym(sample_id_col_name)))

        cleaned_design_matrix <- tryCatch(
            {
                # Create a tibble with just the sample IDs in the order they appear in the data
                sample_order_tibble <- tibble::tibble(temp_sample_id = sample_cols_vector)

                # Join with the design matrix to filter and reorder
                sample_order_tibble |>
                    dplyr::inner_join(design_matrix_char_id,
                        by = c("temp_sample_id" = sample_id_col_name)
                    )
            },
            error = function(e) {
                warning(sprintf("cleanDesignMatrix: Error during inner_join: %s. Returning object unchanged.", e$message))
                return(NULL) # Signal error
            }
        )

        if (is.null(cleaned_design_matrix)) {
            return(theObject) # Return original if join failed
        }

        # Rename the temporary column back to the original sample ID column name
        cleaned_design_matrix <- cleaned_design_matrix |>
            dplyr::rename(!!rlang::sym(sample_id_col_name) := "temp_sample_id")

        # Final check to ensure only expected samples remain (redundant but safe)
        final_cleaned_design <- cleaned_design_matrix |>
            dplyr::filter(!!rlang::sym(sample_id_col_name) %in% sample_cols_vector)

        theObject@design_matrix <- as.data.frame(final_cleaned_design) # Ensure it's stored as data.frame

        return(theObject)
    }
)

#' Get Negative Control Features using ANOVA (Metabolites)
#'
#' Identifies potential negative control features (metabolites) for RUV correction
#' based on ANOVA across a specified grouping variable. Features with the least
#' significant variation across the groups are selected.
#'
#' This method iterates through each assay in the `MetaboliteAssayData` object.
#'
#' @param theObject A `MetaboliteAssayData` object.
#' @param ruv_grouping_variable Character string. The column name in the
#'   `design_matrix` to use for grouping in the ANOVA (e.g., biological replicate
#'   group, batch). Defaults are looked up via `checkParamsObjectFunctionSimplify`
#'   using the key `"ruv_grouping_variable"`.
#' @param percentage_as_neg_ctrl Numeric (0-100). The percentage of total features
#'   to select as negative controls based on ANOVA p-value ranking. Overridden by
#'   `num_neg_ctrl` if provided. Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"metabolites_percentage_as_neg_ctrl"`.
#' @param num_neg_ctrl Integer. The absolute number of features to select as
#'   negative controls. Overrides `percentage_as_neg_ctrl`. Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"metabolites_num_neg_ctrl"`.
#' @param ruv_qval_cutoff Numeric. The q-value (adjusted p-value) threshold used
#'   internally by the ANOVA helper function (typically for filtering before ranking,
#'   though ranking is the primary selection method here). Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"ruv_qval_cutoff"`.
#' @param ruv_fdr_method Character string. The method used for p-value adjustment
#'   (e.g., "BH", "fdr"). Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"ruv_fdr_method"`.
#'
#' @return A named list, where each element corresponds to an assay in the input
#'   object. Each element contains a logical vector indicating which features
#'   (metabolites) in that assay were selected as negative controls. The vector
#'   is named with the feature IDs.
#'
#' @importFrom methods slot
#' @importFrom purrr map set_names map_lgl
#' @importFrom tibble column_to_rownames as_tibble is_tibble
#' @importFrom dplyr pull select filter all_of any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom logger log_info log_warn
#' @describeIn getNegCtrlMetabAnova Method for MetaboliteAssayData
#' @export
setMethod(
    f = "getNegCtrlMetabAnova",
    signature = "MetaboliteAssayData",
    definition = function(theObject,
                          ruv_grouping_variable = NULL, # These args are now effectively ignored for default resolution
                          percentage_as_neg_ctrl = NULL,
                          num_neg_ctrl = NULL,
                          ruv_qval_cutoff = NULL,
                          ruv_fdr_method = NULL) {
        message("+===========================================================================+")
        message("|  DEBUG66: Entering getNegCtrlMetabAnova (MetaboliteAssayData)             |")
        message("+===========================================================================+")

        assay_list <- methods::slot(theObject, "metabolite_data")
        metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        group_id <- methods::slot(theObject, "group_id") # Needed for helper
        sample_id <- methods::slot(theObject, "sample_id")

        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] metabolite_id_col = '%s', group_id = '%s', sample_id = '%s'", metabolite_id_col_name, group_id, sample_id))
        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Number of assays: %d, design_matrix rows: %d", length(assay_list), nrow(design_matrix)))
        message(sprintf(
            "   DEBUG66 [getNegCtrlMetabAnova] Function args: ruv_grouping_variable = %s, percentage_as_neg_ctrl = %s",
            ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable),
            ifelse(is.null(percentage_as_neg_ctrl), "NULL", as.character(percentage_as_neg_ctrl))
        ))

        # --- Resolve Global Parameters (Mimicking Protein version exactly) ---
        # Get values from object args slot, falling back to hardcoded defaults.
        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = "group") # Hardcoded default
        ruv_qval_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_qval_cutoff", default_value = 0.05) # Hardcoded default
        ruv_fdr_method_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_fdr_method", default_value = "BH") # Hardcoded default

        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Resolved: ruv_grouping_variable_final = '%s'", ruv_grouping_variable_final))
        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Resolved: ruv_qval_cutoff_final = %s, ruv_fdr_method_final = '%s'", ruv_qval_cutoff_final, ruv_fdr_method_final))

        # Update object args with resolved values (Mimicking Protein version)
        theObject <- updateParamInObject(theObject, "ruv_grouping_variable") # Correct: Only 2 args
        theObject <- updateParamInObject(theObject, "ruv_qval_cutoff") # Correct: Only 2 args
        theObject <- updateParamInObject(theObject, "ruv_fdr_method") # Correct: Only 2 args

        log_info("Starting Negative Control selection using ANOVA for metabolites.")
        log_info("Parameters (Resolved):")
        log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
        log_info("  - RUV Q-value Cutoff: {ruv_qval_cutoff_final}")
        log_info("  - RUV FDR Method: {ruv_fdr_method_final}")
        # Percentage/Num are resolved per assay

        if (length(assay_list) == 0) {
            message("   DEBUG66 [getNegCtrlMetabAnova] WARNING: No assays found! Returning empty list.")
            log_warn("No assays found in `metabolite_data` slot. Returning empty list.")
            return(list())
        }
        # Ensure list is named
        assay_names <- names(assay_list)
        if (is.null(assay_names)) {
            assay_names <- paste0("Assay_", seq_along(assay_list))
            message("   DEBUG66 [getNegCtrlMetabAnova] Assay list was unnamed. Using default names.")
            log_warn("Assay list was unnamed. Using default names.")
        }
        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay names: %s", paste(assay_names, collapse = ", ")))

        # --- Process Each Assay ---
        control_features_list <- lapply(seq_along(assay_list), function(i) {
            assay_name <- assay_names[i]
            assay_tibble <- assay_list[[i]]
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] === Processing assay: %s ===", assay_name))
            message(sprintf("-- Processing assay for NegCtrl ANOVA: %s", assay_name))

            # --- Basic Checks ---
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': tibble rows = %d, cols = %d", assay_name, nrow(assay_tibble), ncol(assay_tibble)))
            if (!tibble::is_tibble(assay_tibble)) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Not a tibble, attempting coercion", assay_name))
                log_warn("Assay '{assay_name}' is not a tibble. Attempting coercion.", .logr = TRUE)
                assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Coercion FAILED - %s", assay_name, e$message))
                    log_warn("Failed to coerce assay '{assay_name}' to tibble: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                })
                if (is.null(assay_tibble)) {
                    return(NULL)
                } # Skip assay
            }
            if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - metabolite ID column '%s' not found", assay_name, metabolite_id_col_name))
                log_warn("Assay '{assay_name}': Metabolite ID column '{metabolite_id_col_name}' not found. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_tibble) == 0) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - zero rows", assay_name))
                log_warn("Assay '{assay_name}': Contains zero rows (features). Skipping.", .logr = TRUE)
                return(NULL)
            }

            # --- Identify Sample Columns ---
            design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) {
                character(0)
            })
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': design_samples count = %d", assay_name, length(design_samples)))
            if (length(design_samples) == 0) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - no design samples found", assay_name))
                log_warn("Assay '{assay_name}': Could not extract valid sample IDs from design matrix column '{sample_id}'. Skipping.", .logr = TRUE)
                return(NULL)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': sample_cols count = %d (matched)", assay_name, length(sample_cols)))
            if (length(sample_cols) < 2) { # Need at least 2 samples for ANOVA
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - fewer than 2 sample columns", assay_name))
                log_warn("Assay '{assay_name}': Fewer than 2 sample columns identified matching design matrix. Skipping ANOVA.", .logr = TRUE)
                return(NULL)
            }
            # Ensure sample columns are numeric
            non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
            if (length(non_numeric_samples) > 0) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Coercing %d non-numeric sample columns", assay_name, length(non_numeric_samples)))
                log_warn("Assay '{assay_name}': Non-numeric sample columns found: {paste(non_numeric_samples, collapse=', ')}. Attempting coercion.", .logr = TRUE)
                assay_tibble <- assay_tibble |>
                    dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # --- Prepare Matrix for Helper ---
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Converting to matrix...", assay_name))
            assay_matrix <- tryCatch(
                {
                    assay_tibble |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |>
                        tibble::column_to_rownames(var = metabolite_id_col_name) |>
                        as.matrix()
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Matrix conversion FAILED - %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error converting tibble to matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(NULL)
            }

            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Matrix dims = %d x %d", assay_name, nrow(assay_matrix), ncol(assay_matrix)))
            assay_matrix[!is.finite(assay_matrix)] <- NA

            # Check for sufficient valid data
            valid_rows <- rowSums(!is.na(assay_matrix)) > 1
            valid_cols <- colSums(!is.na(assay_matrix)) > 1
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': valid_rows = %d, valid_cols = %d", assay_name, sum(valid_rows), sum(valid_cols)))
            if (sum(valid_rows) < 2 || sum(valid_cols) < 2) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - insufficient non-NA data", assay_name))
                log_warn("Assay '{assay_name}': Insufficient non-NA data points (<2 features or <2 samples with data) for ANOVA. Skipping.", .logr = TRUE)
                return(NULL)
            }
            assay_matrix_filt <- assay_matrix[valid_rows, valid_cols, drop = FALSE]
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Filtered matrix dims = %d x %d", assay_name, nrow(assay_matrix_filt), ncol(assay_matrix_filt)))

            # Filter design matrix to match valid columns in assay_matrix_filt
            design_matrix_filtered <- design_matrix |>
                dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix_filt)) |>
                as.data.frame() # Helper might expect data.frame
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Filtered design_matrix rows = %d", assay_name, nrow(design_matrix_filtered)))

            # Check if grouping variable has enough levels/samples after filtering
            if (!ruv_grouping_variable_final %in% colnames(design_matrix_filtered)) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - grouping variable '%s' not found in design matrix", assay_name, ruv_grouping_variable_final))
                log_warn("Assay '{assay_name}': Grouping variable '{ruv_grouping_variable_final}' not found in filtered design matrix. Skipping ANOVA.", .logr = TRUE)
                return(NULL)
            }
            group_counts <- table(design_matrix_filtered[[ruv_grouping_variable_final]])
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Group counts: %s", assay_name, paste(names(group_counts), "=", group_counts, collapse = ", ")))
            if (length(group_counts) < 2 || any(group_counts < 2)) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - insufficient groups (%d) or samples per group", assay_name, length(group_counts)))
                log_warn("Assay '{assay_name}': Insufficient groups ({length(group_counts)}) or samples per group (<2) for ANOVA based on '{ruv_grouping_variable_final}' after filtering. Skipping.", .logr = TRUE)
                return(NULL)
            }


            # --- Resolve Assay-Specific Parameters (Revised for flexible percentage) ---

            # Determine the percentage for *this specific assay*
            percentage_to_use_for_assay <- NULL

            # Check 1: Explicit function argument provided?
            if (!is.null(percentage_as_neg_ctrl)) {
                if ((is.list(percentage_as_neg_ctrl) || is.vector(percentage_as_neg_ctrl)) && !is.null(names(percentage_as_neg_ctrl))) {
                    # Check 1a: Named list/vector provided - try to match name
                    if (assay_name %in% names(percentage_as_neg_ctrl)) {
                        percentage_to_use_for_assay <- percentage_as_neg_ctrl[[assay_name]]
                        log_info("   Assay '{assay_name}': Using percentage from named argument: {percentage_to_use_for_assay}", .logr = TRUE)
                    }
                } else if (is.vector(percentage_as_neg_ctrl) && is.null(names(percentage_as_neg_ctrl)) && length(percentage_as_neg_ctrl) == length(assay_list)) {
                    # Check 1b: Unnamed vector of correct length provided - use position
                    percentage_to_use_for_assay <- percentage_as_neg_ctrl[[i]]
                    log_info("   Assay '{assay_name}': Using percentage from positional argument: {percentage_to_use_for_assay}", .logr = TRUE)
                } else if (is.numeric(percentage_as_neg_ctrl) && length(percentage_as_neg_ctrl) == 1) {
                    # Check 1c: Single numeric value provided
                    percentage_to_use_for_assay <- percentage_as_neg_ctrl
                    log_info("   Assay '{assay_name}': Using single percentage value from argument: {percentage_to_use_for_assay}", .logr = TRUE)
                }
            }

            # Check 2: Fallback to config/default if not found in explicit args
            if (is.null(percentage_to_use_for_assay)) {
                percentage_to_use_for_assay <- checkParamsObjectFunctionSimplify(
                    theObject, "percentage_as_neg_ctrl",
                    default_value = 10
                ) # Generic key, hardcoded default 10
                log_info("   Assay '{assay_name}': Using percentage from config/default: {percentage_to_use_for_assay}", .logr = TRUE)
                # Update object args only if we resolved from config/default
                # Avoid overwriting if a specific value was passed via function arg
                if (is.null(percentage_as_neg_ctrl)) { # Only update args if function call arg was NULL
                    theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
                }
            }

            # Validate the resolved percentage
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': percentage_to_use_for_assay = %s", assay_name, ifelse(is.null(percentage_to_use_for_assay), "NULL", percentage_to_use_for_assay)))
            if (!is.numeric(percentage_to_use_for_assay) || length(percentage_to_use_for_assay) != 1 || is.na(percentage_to_use_for_assay) || percentage_to_use_for_assay < 0 || percentage_to_use_for_assay > 100) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - Invalid percentage", assay_name))
                log_warn("   Assay '{assay_name}': Invalid percentage resolved ({percentage_to_use_for_assay}). Must be numeric between 0 and 100. Skipping assay.", .logr = TRUE)
                return(NULL)
            }

            # Calculate default num_neg_ctrl based on resolved percentage and *filtered* matrix
            # We use the *resolved* percentage for this assay now
            default_num_neg_ctrl <- round(nrow(assay_matrix_filt) * percentage_to_use_for_assay / 100, 0)
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': default_num_neg_ctrl = %d (from %d features * %.1f%%)", assay_name, default_num_neg_ctrl, nrow(assay_matrix_filt), percentage_to_use_for_assay))

            # Resolve num_neg_ctrl (prioritize function arg, then config, then calculated default)
            num_neg_ctrl_assay <- NULL
            if (!is.null(num_neg_ctrl)) { # Check explicit function arg first
                # Add similar logic here if you want num_neg_ctrl to also be per-assay via list/vector
                # For now, assume num_neg_ctrl function arg is single value if provided
                if (is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0) {
                    num_neg_ctrl_assay <- num_neg_ctrl
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Using num_neg_ctrl from argument = %d", assay_name, num_neg_ctrl_assay))
                    log_info("   Assay '{assay_name}': Using num_neg_ctrl from argument: {num_neg_ctrl_assay}", .logr = TRUE)
                } else {
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Invalid num_neg_ctrl argument, ignoring", assay_name))
                    log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl argument provided. Ignoring.", .logr = TRUE)
                }
            }
            if (is.null(num_neg_ctrl_assay)) { # If not provided or invalid in args, check config/default
                num_neg_ctrl_assay <- checkParamsObjectFunctionSimplify(
                    theObject, "num_neg_ctrl",
                    default_value = default_num_neg_ctrl
                ) # Generic key
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Using num_neg_ctrl from config/default = %d", assay_name, num_neg_ctrl_assay))
                log_info("   Assay '{assay_name}': Using num_neg_ctrl from config/calculated default: {num_neg_ctrl_assay}", .logr = TRUE)
                # Update object args only if resolved from config/default and function arg was NULL/invalid
                if (is.null(num_neg_ctrl) || !(is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0)) {
                    theObject <- updateParamInObject(theObject, "num_neg_ctrl")
                }
            }
            # Validate the resolved num_neg_ctrl
            if (!is.numeric(num_neg_ctrl_assay) || length(num_neg_ctrl_assay) != 1 || is.na(num_neg_ctrl_assay) || num_neg_ctrl_assay < 0) {
                message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - Invalid num_neg_ctrl = %s", assay_name, num_neg_ctrl_assay))
                log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl resolved ({num_neg_ctrl_assay}). Must be non-negative integer. Skipping assay.", .logr = TRUE)
                return(NULL)
            }
            # Ensure integer
            num_neg_ctrl_assay <- as.integer(num_neg_ctrl_assay)

            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Final num_neg_ctrl = %d", assay_name, num_neg_ctrl_assay))
            log_info("  Assay '{assay_name}': Final Neg Ctrl Count: {num_neg_ctrl_assay} (based on percentage: {percentage_to_use_for_assay}%)", .logr = TRUE)


            # --- Prepare Design Matrix for Helper ---
            # Helper expects rownames = sample IDs, and group_id column removed
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Preparing design_matrix_for_helper...", assay_name))
            design_matrix_for_helper <- tryCatch(
                {
                    design_matrix_filtered |>
                        tibble::column_to_rownames(var = sample_id) |>
                        dplyr::select(-dplyr::any_of(group_id)) # Remove group_id if it exists
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - design_matrix_for_helper error: %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error preparing design matrix for helper: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(design_matrix_for_helper)) {
                return(NULL)
            }
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': design_matrix_for_helper ready, dims = %d x %d", assay_name, nrow(design_matrix_for_helper), ncol(design_matrix_for_helper)))


            # --- Call Helper ---
            # **ASSUMPTION**: getNegCtrlProtAnovaHelper can handle the data
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': Calling getNegCtrlProtAnovaHelper...", assay_name))
            control_indices_assay <- tryCatch(
                {
                    getNegCtrlProtAnovaHelper(
                        assay_matrix_filt, # Matrix with features as rows, samples as cols
                        design_matrix = design_matrix_for_helper,
                        grouping_variable = ruv_grouping_variable_final,
                        # Pass the specifically resolved percentage and number for *this* assay
                        percentage_as_neg_ctrl = percentage_to_use_for_assay,
                        num_neg_ctrl = num_neg_ctrl_assay,
                        ruv_qval_cutoff = ruv_qval_cutoff_final,
                        ruv_fdr_method = ruv_fdr_method_final
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': FAIL - getNegCtrlProtAnovaHelper error: %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error calling getNegCtrlProtAnovaHelper: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL) # Return NULL for this assay on error
                }
            )

            num_ctrl_selected <- sum(control_indices_assay, na.rm = TRUE)
            message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Assay '%s': SUCCESS - Selected %d control features", assay_name, num_ctrl_selected))
            log_info("  Assay '{assay_name}': Selected {sum(control_indices_assay, na.rm = TRUE)} control features.", .logr = TRUE)
            return(control_indices_assay)
        })

        # Set names for the list of results
        names(control_features_list) <- assay_names

        # Remove NULL elements (skipped assays)
        final_control_list <- control_features_list[!sapply(control_features_list, is.null)]

        message(sprintf("   DEBUG66 [getNegCtrlMetabAnova] Finished. Returning %d assay results.", length(final_control_list)))
        message("+===========================================================================+")
        message("|  DEBUG66: Exiting getNegCtrlMetabAnova                                    |")
        message("+===========================================================================+")
        log_info("Finished Negative Control selection for {length(final_control_list)} assay(s).")
        return(final_control_list)
    }
)

#' @title RUV Canonical Correlation for MetaboliteAssayData
#' @name ruvCancor,MetaboliteAssayData-method
#' @importFrom methods slot
#' @importFrom purrr map set_names map_lgl
#' @importFrom tibble column_to_rownames is_tibble as_tibble
#' @importFrom dplyr pull select filter all_of any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom logger log_info log_warn log_error
#' @export
#' @export
setMethod(
    f = "ruvCancor",
    signature = "MetaboliteAssayData",
    definition = function(theObject, ctrl = NULL, num_components_to_impute = NULL, ruv_grouping_variable = NULL) {
        message("+===========================================================================+")
        message("|  DEBUG66: Entering ruvCancor (MetaboliteAssayData)                        |")
        message("+===========================================================================+")

        assay_list <- methods::slot(theObject, "metabolite_data")
        metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id <- methods::slot(theObject, "sample_id")
        # group_id is not directly used here but good practice to extract if needed later

        message(sprintf("   DEBUG66 [ruvCancor] Function args: ctrl is.null = %s, ruv_grouping_variable = %s", is.null(ctrl), ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable)))
        message(sprintf("   DEBUG66 [ruvCancor] Number of assays: %d", length(assay_list)))

        # --- Resolve Global Parameters ---
        # Use generic keys as per handover doc
        # Default ctrl=NULL means it MUST be provided in args or function call
        ctrl_final <- checkParamsObjectFunctionSimplify(theObject, "ctrl", default_value = NULL)
        num_components_to_impute_final <- checkParamsObjectFunctionSimplify(theObject, "num_components_to_impute", default_value = 2)
        # Default ruv_grouping_variable = NULL, MUST be provided
        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = NULL)

        message(sprintf(
            "   DEBUG66 [ruvCancor] Resolved ctrl_final is.null = %s, is.list = %s, class = '%s'",
            is.null(ctrl_final), is.list(ctrl_final), class(ctrl_final)[1]
        ))
        if (is.list(ctrl_final)) {
            message(sprintf("   DEBUG66 [ruvCancor] ctrl_final names: %s", paste(names(ctrl_final), collapse = ", ")))
        }
        message(sprintf("   DEBUG66 [ruvCancor] Resolved ruv_grouping_variable_final = '%s'", ifelse(is.null(ruv_grouping_variable_final), "NULL", ruv_grouping_variable_final)))
        message(sprintf("   DEBUG66 [ruvCancor] Resolved num_components_to_impute_final = %s", num_components_to_impute_final))

        # Update object args (using generic keys)
        theObject <- updateParamInObject(theObject, "ctrl")
        theObject <- updateParamInObject(theObject, "num_components_to_impute")
        theObject <- updateParamInObject(theObject, "ruv_grouping_variable")

        log_info("Starting RUV Canonical Correlation plot generation for metabolites.")
        log_info("Parameters (Resolved):")
        log_info("  - Control Features Key: 'ctrl' (Value type depends on input/config)")
        log_info("  - Num Imputation Components: {num_components_to_impute_final}")
        log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")

        # --- Input Validation ---
        if (is.null(ctrl_final)) {
            message("   DEBUG66 [ruvCancor] FAIL - ctrl_final is NULL, stopping!")
            log_error("Negative control features ('ctrl') must be provided either via function argument or object configuration ('args$ctrl').")
            stop("Missing required 'ctrl' parameter for ruvCancor.")
        }
        if (is.null(ruv_grouping_variable_final)) {
            message("   DEBUG66 [ruvCancor] FAIL - ruv_grouping_variable_final is NULL, stopping!")
            log_error("RUV grouping variable ('ruv_grouping_variable') must be provided either via function argument or object configuration.")
            stop("Missing required 'ruv_grouping_variable' parameter for ruvCancor.")
        }
        if (!ruv_grouping_variable_final %in% colnames(design_matrix)) {
            message(sprintf("   DEBUG66 [ruvCancor] FAIL - ruv_grouping_variable '%s' not in design_matrix columns!", ruv_grouping_variable_final))
            log_error("The 'ruv_grouping_variable' ('{ruv_grouping_variable_final}') is not a column in the design matrix.")
            stop(paste0("The 'ruv_grouping_variable = ", ruv_grouping_variable_final, "' is not a column in the design matrix."))
        }
        if (!is.numeric(num_components_to_impute_final) || is.na(num_components_to_impute_final) || num_components_to_impute_final < 1) {
            message(sprintf("   DEBUG66 [ruvCancor] FAIL - invalid num_components_to_impute = %s", num_components_to_impute_final))
            log_error("Invalid 'num_components_to_impute': {num_components_to_impute_final}. Must be a positive integer.")
            stop(paste0("The num_components_to_impute = ", num_components_to_impute_final, " value is invalid."))
        }


        if (length(assay_list) == 0) {
            message("   DEBUG66 [ruvCancor] WARNING - no assays found, returning empty list")
            log_warn("No assays found in `metabolite_data` slot. Returning empty list.")
            return(list())
        }
        # Ensure list is named
        assay_names <- names(assay_list)
        if (is.null(assay_names)) {
            assay_names <- paste0("Assay_", seq_along(assay_list))
            message("   DEBUG66 [ruvCancor] Assay list was unnamed. Using default names.")
            log_warn("Assay list was unnamed. Using default names.")
        }
        message(sprintf("   DEBUG66 [ruvCancor] Assay names: %s", paste(assay_names, collapse = ", ")))

        # --- Process Each Assay ---
        cancor_plots_list <- lapply(seq_along(assay_list), function(i) {
            assay_name <- assay_names[i]
            assay_tibble <- assay_list[[i]]
            message(sprintf("   DEBUG66 [ruvCancor] === Processing assay: %s ===", assay_name))
            message(sprintf("-- Processing assay for RUV Cancor Plot: %s", assay_name))

            # --- Basic Checks ---
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': tibble rows = %d, cols = %d", assay_name, nrow(assay_tibble), ncol(assay_tibble)))
            if (!tibble::is_tibble(assay_tibble)) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Not a tibble, attempting coercion", assay_name))
                log_warn("Assay '{assay_name}' is not a tibble. Attempting coercion.", .logr = TRUE)
                assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                    message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Coercion FAILED - %s", assay_name, e$message))
                    log_warn("Failed to coerce assay '{assay_name}' to tibble: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                })
                if (is.null(assay_tibble)) {
                    return(NULL)
                } # Skip assay
            }
            if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - metabolite ID column not found", assay_name))
                log_warn("Assay '{assay_name}': Metabolite ID column '{metabolite_id_col_name}' not found. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_tibble) == 0) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - zero rows", assay_name))
                log_warn("Assay '{assay_name}': Contains zero rows (features). Skipping.", .logr = TRUE)
                return(NULL)
            }

            # --- Identify Sample Columns ---
            design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) {
                character(0)
            })
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': design_samples count = %d", assay_name, length(design_samples)))
            if (length(design_samples) == 0) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - no design samples found", assay_name))
                log_warn("Assay '{assay_name}': Could not extract valid sample IDs from design matrix column '{sample_id}'. Skipping.", .logr = TRUE)
                return(NULL)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': sample_cols count = %d", assay_name, length(sample_cols)))
            if (length(sample_cols) < 2) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - fewer than 2 sample columns", assay_name))
                log_warn("Assay '{assay_name}': Fewer than 2 sample columns identified matching design matrix. Skipping RUV cancor plot.", .logr = TRUE)
                return(NULL)
            }
            # Ensure sample columns are numeric
            non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
            if (length(non_numeric_samples) > 0) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Coercing %d non-numeric sample columns", assay_name, length(non_numeric_samples)))
                log_warn("Assay '{assay_name}': Non-numeric sample columns found: {paste(non_numeric_samples, collapse=', ')}. Attempting coercion.", .logr = TRUE)
                assay_tibble <- assay_tibble |>
                    dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # --- Prepare Matrix (Features x Samples) ---
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Converting to matrix...", assay_name))
            assay_matrix <- tryCatch(
                {
                    assay_tibble |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |>
                        tibble::column_to_rownames(var = metabolite_id_col_name) |>
                        as.matrix()
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Matrix conversion FAILED - %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error converting tibble to matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(NULL)
            }

            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Matrix dims = %d x %d", assay_name, nrow(assay_matrix), ncol(assay_matrix)))
            assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf first

            # --- Filter Design Matrix ---
            # Ensure design matrix matches the actual columns used in the assay_matrix
            design_matrix_filtered <- design_matrix |>
                dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix)) |>
                as.data.frame() # Ensure it's a data.frame if needed
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Filtered design_matrix rows = %d", assay_name, nrow(design_matrix_filtered)))

            # --- Prepare Y (Samples x Features) ---
            Y_matrix <- t(assay_matrix[, as.character(design_matrix_filtered[[sample_id]]), drop = FALSE]) # Ensure column order matches filtered design
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Y_matrix dims = %d x %d (samples x features)", assay_name, nrow(Y_matrix), ncol(Y_matrix)))

            # --- Imputation (using mixOmics::impute.nipals) ---
            if (anyNA(Y_matrix)) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': NAs detected, performing NIPALS imputation", assay_name))
                log_info("   Assay '{assay_name}': Missing values detected. Performing NIPALS imputation with {num_components_to_impute_final} components.", .logr = TRUE)
                Y_imputed <- tryCatch(
                    {
                        mixOmics::impute.nipals(Y_matrix, ncomp = num_components_to_impute_final)
                    },
                    error = function(e) {
                        message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': NIPALS imputation FAILED - %s", assay_name, e$message))
                        log_warn("Assay '{assay_name}': Error during NIPALS imputation: {e$message}. Skipping.", .logr = TRUE)
                        return(NULL)
                    }
                )
                if (is.null(Y_imputed)) {
                    return(NULL)
                }
                Y_final <- Y_imputed
            } else {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': No NAs in Y_matrix, skipping imputation", assay_name))
                Y_final <- Y_matrix
            }

            # --- Prepare X (Grouping Variable) ---
            if (!ruv_grouping_variable_final %in% colnames(design_matrix_filtered)) {
                # This check should be redundant due to the initial check, but safe
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - grouping variable lost after filtering", assay_name))
                log_error("Assay '{assay_name}': Grouping variable '{ruv_grouping_variable_final}' lost after filtering design matrix. This shouldn't happen.", .logr = TRUE)
                return(NULL)
            }
            X_vector <- design_matrix_filtered[[ruv_grouping_variable_final]]
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': X_vector length = %d, unique values = %d", assay_name, length(X_vector), length(unique(X_vector))))

            # --- Prepare ctl (Control Features Indices/Logical) ---
            # Resolve control features specific to this assay
            ctrl_assay <- NULL
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Resolving ctrl_assay, ctrl_final is.list = %s", assay_name, is.list(ctrl_final)))
            if (is.list(ctrl_final)) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Checking if '%s' in names(ctrl_final) = %s", assay_name, assay_name, assay_name %in% names(ctrl_final)))
                if (assay_name %in% names(ctrl_final)) {
                    ctrl_assay <- ctrl_final[[assay_name]]
                    message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Found assay-specific controls, class = '%s'", assay_name, class(ctrl_assay)[1]))
                    log_info("   Assay '{assay_name}': Found assay-specific controls in 'ctrl' list.", .logr = TRUE)
                } else {
                    message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - assay not in ctrl_final names", assay_name))
                    log_warn("Assay '{assay_name}': 'ctrl' is a list, but does not contain an element named '{assay_name}'. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            } else {
                # If ctrl_final is not a list, assume it's a global vector (numeric, logical, character)
                # This maintains backwards compatibility if a global ctrl vector is provided
                ctrl_assay <- ctrl_final
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Using globally provided ctrl, class = '%s'", assay_name, class(ctrl_assay)[1]))
                log_info("   Assay '{assay_name}': Using the globally provided 'ctrl' vector.", .logr = TRUE)
            }

            if (is.null(ctrl_assay)) {
                # This case should ideally be caught by the list check above, but as a safeguard:
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - ctrl_assay is NULL", assay_name))
                log_warn("Assay '{assay_name}': Failed to resolve control features for this assay. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # `ruv_cancorplot` expects `ctl` relative to the *columns* of Y_final (features)
            feature_names_in_assay <- colnames(Y_final)
            control_indices_assay <- NULL # Initialize

            if (is.numeric(ctrl_assay)) {
                # If numeric indices are provided, check bounds
                if (any(ctrl_assay < 1) || any(ctrl_assay > length(feature_names_in_assay))) {
                    log_warn("Assay '{assay_name}': Numeric 'ctrl' indices are out of bounds for the features in this assay ({length(feature_names_in_assay)}). Skipping.", .logr = TRUE)
                    return(NULL)
                }
                control_indices_assay <- ctrl_assay
            } else if (is.logical(ctrl_assay)) {
                # If logical, check length relative to the features *in this specific assay*
                if (length(ctrl_assay) != length(feature_names_in_assay)) {
                    # We need to align the logical vector with the current assay's features if names are present
                    if (!is.null(names(ctrl_assay))) {
                        feature_match <- match(feature_names_in_assay, names(ctrl_assay))
                        if (anyNA(feature_match)) {
                            log_warn("Assay '{assay_name}': Some assay features not found in named logical 'ctrl' vector. Skipping.", .logr = TRUE)
                            return(NULL)
                        }
                        control_indices_assay <- ctrl_assay[feature_match]
                        if (length(control_indices_assay) != length(feature_names_in_assay)) {
                            log_warn("Assay '{assay_name}': Length mismatch after aligning named logical 'ctrl' vector ({length(control_indices_assay)}) with assay features ({length(feature_names_in_assay)}). Skipping.", .logr = TRUE)
                            return(NULL)
                        }
                        log_info("   Assay '{assay_name}': Aligned named logical 'ctrl' vector to assay features.", .logr = TRUE)
                    } else {
                        log_warn("Assay '{assay_name}': Unnamed logical 'ctrl' vector length ({length(ctrl_assay)}) does not match number of features ({length(feature_names_in_assay)}). Skipping.", .logr = TRUE)
                        return(NULL)
                    }
                } else {
                    # Length matches, assume order is correct
                    control_indices_assay <- ctrl_assay
                }
            } else if (is.character(ctrl_assay)) {
                # If character IDs, find which ones are in the current assay
                control_indices_assay <- feature_names_in_assay %in% ctrl_assay
                if (sum(control_indices_assay) == 0) {
                    log_warn("Assay '{assay_name}': None of the provided character 'ctrl' IDs were found in the features of this assay. Skipping.", .logr = TRUE)
                    return(NULL)
                }
                # ruv_cancorplot expects logical or numeric indices, convert the logical vector derived from character IDs
                # control_indices_assay remains logical here, which is valid for ruv_cancorplot ctl
            } else {
                log_warn("Assay '{assay_name}': Invalid type for resolved 'ctrl_assay' parameter. Expected numeric, logical, or character. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Final check on number of controls (use sum for logical, length for numeric)
            # Need to ensure control_indices_assay is not NULL before checking
            if (is.null(control_indices_assay)) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - control_indices_assay is NULL", assay_name))
                log_warn("Assay '{assay_name}': Control indices could not be determined. Skipping.", .logr = TRUE)
                return(NULL)
            }
            num_controls_found <- if (is.logical(control_indices_assay)) sum(control_indices_assay, na.rm = TRUE) else length(control_indices_assay)
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': num_controls_found = %d", assay_name, num_controls_found))
            if (num_controls_found < 5) {
                message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': FAIL - fewer than 5 controls (%d)", assay_name, num_controls_found))
                log_warn("Assay '{assay_name}': Fewer than 5 negative control features found/specified for this assay ({num_controls_found}). RUV results may be unreliable. Skipping cancor plot.", .logr = TRUE)
                # Proceeding might still work but is discouraged by the original protein code's check
                return(NULL) # Skip plot generation as per original check
            }
            log_info("   Assay '{assay_name}': Using {num_controls_found} control features for cancor plot.", .logr = TRUE)

            # --- Call ruv_cancorplot ---
            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': Calling ruv_cancorplot...", assay_name))
            cancor_plot_assay <- tryCatch(
                {
                    # Ensure ruv_cancorplot is loaded/available in the environment
                    # Requires Y = samples x features matrix
                    # Requires X = grouping vector (same length as nrow(Y))
                    # Requires ctl = logical/numeric vector identifying control columns in Y
                    ruv_cancorplot(
                        Y = Y_final,
                        X = X_vector,
                        ctl = control_indices_assay
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': ruv_cancorplot FAILED - %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error calling ruv_cancorplot: {e$message}. Check if the function exists and is loaded correctly. Skipping.", .logr = TRUE)
                    return(NULL) # Return NULL for this assay on error
                }
            )

            message(sprintf("   DEBUG66 [ruvCancor] Assay '%s': SUCCESS - cancor_plot created", assay_name))
            return(cancor_plot_assay)
        })

        # Set names for the list of plots
        names(cancor_plots_list) <- assay_names

        # Remove NULL elements (skipped assays)
        final_plots_list <- cancor_plots_list[!sapply(cancor_plots_list, is.null)]

        message(sprintf("   DEBUG66 [ruvCancor] Finished. Returning %d assay plots.", length(final_plots_list)))
        message("+===========================================================================+")
        message("|  DEBUG66: Exiting ruvCancor                                               |")
        message("+===========================================================================+")
        log_info("Finished RUV Canonical Correlation plot generation for {length(final_plots_list)} assay(s).")
        return(final_plots_list)
    }
)

#' Apply RUV-III Correction with Varying K
#'
#' Applies the RUV-III correction method to each assay within a
#' `MetaboliteAssayData` object. This method accounts for unwanted variation
#' using control features and a replicate structure matrix. It allows for
#' potentially different numbers of factors (`k`) to be removed for each assay.
#'
#' @param theObject A `MetaboliteAssayData` object.
#' @param ruv_grouping_variable Character string. The column name in the
#'   `design_matrix` that defines the replicate structure for RUV-III (e.g.,
#'   biological groups where variation *within* the group is considered noise).
#'   Defaults are looked up via `checkParamsObjectFunctionSimplify` using the
#'   key `"ruv_grouping_variable"`. Must be provided.
#' @param ruv_number_k An integer or a named list/vector.
#'   - If an integer, this number of factors (`k`) is removed from all assays.
#'   - If a named list/vector, the names must correspond to the assay names
#'     in `metabolite_data`. The value associated with each name specifies the
#'     `k` for that assay. Assays not named will use a default `k`.
#'   Defaults are looked up via `checkParamsObjectFunctionSimplify` using the key
#'   `"ruv_number_k"`.
#' @param ctrl A logical vector, numeric vector, character vector, or a named list.
#'   - If a vector, it specifies the control features used for all assays.
#'     Can be logical (matching features), numeric indices, or character IDs.
#'   - If a named list, the names must correspond to assay names. Each element
#'     should be a vector specifying controls for that specific assay.
#'   Defaults are looked up via `checkParamsObjectFunctionSimplify` using the key
#'   `"ctrl"`. Must be provided.
#' @param num_components_to_impute Integer. The number of principal components
#'   to use for NIPALS imputation if missing values are present before RUV.
#'   Defaults are looked up via `checkParamsObjectFunctionSimplify` using the key
#'   `"num_components_to_impute"`.
#'
#' @return A modified `MetaboliteAssayData` object where the `metabolite_data`
#'   slot contains the RUV-corrected assay data. Features or samples with only
#'   NA/NaN values after correction are removed. The `design_matrix` is updated
#'   via `cleanDesignMatrix`.
#'
#' @importFrom methods slot slot<-
#' @importFrom purrr map map_lgl map_chr set_names
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble is_tibble
#' @importFrom dplyr pull select filter all_of any_of mutate across left_join relocate distinct
#' @importFrom rlang sym !! :=
#' @importFrom logger log_info log_warn log_error
#' @importFrom stringr str_split
#' @describeIn ruvIII_C_Varying Method for MetaboliteAssayData
#' @export
setMethod(
    f = "ruvIII_C_Varying",
    signature = "MetaboliteAssayData",
    definition = function(theObject,
                          ruv_grouping_variable = NULL,
                          ruv_number_k = NULL,
                          ctrl = NULL) {
        message("+===========================================================================+")
        message("|  DEBUG66: Entering ruvIII_C_Varying (MetaboliteAssayData)                 |")
        message("+===========================================================================+")

        assay_list <- methods::slot(theObject, "metabolite_data")
        metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id <- methods::slot(theObject, "sample_id")
        group_id <- methods::slot(theObject, "group_id") # Extract for context, though not directly used below
        technical_replicate_id <- methods::slot(theObject, "technical_replicate_id") # Extract for context

        message(sprintf(
            "   DEBUG66 [ruvIII_C_Varying] Function args: ruv_grouping_variable = %s, ruv_number_k is.null = %s, ctrl is.null = %s",
            ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable), is.null(ruv_number_k), is.null(ctrl)
        ))
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Number of assays: %d", length(assay_list)))

        # --- Resolve Parameters (Prioritize function args, then object args) ---

        # 1. RUV Grouping Variable
        if (!is.null(ruv_grouping_variable)) {
            ruv_grouping_variable_final <- ruv_grouping_variable
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_grouping_variable from function arg: '%s'", ruv_grouping_variable_final))
            log_info("Using 'ruv_grouping_variable' from function argument: {ruv_grouping_variable_final}")
        } else {
            ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = NULL)
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_grouping_variable from object/default: '%s'", ifelse(is.null(ruv_grouping_variable_final), "NULL", ruv_grouping_variable_final)))
            log_info("Using 'ruv_grouping_variable' from object args or default: {ruv_grouping_variable_final}")
        }

        # 2. RUV Number K (k)
        if (!is.null(ruv_number_k)) {
            ruv_number_k_resolved <- ruv_number_k
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_number_k from function arg, class = '%s'", class(ruv_number_k_resolved)[1]))
            log_info("Using 'ruv_number_k' from function argument.")
        } else {
            ruv_number_k_resolved <- checkParamsObjectFunctionSimplify(theObject, "ruv_number_k", default_value = NULL)
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Using ruv_number_k from object/default: is.null = %s, class = '%s'",
                is.null(ruv_number_k_resolved), class(ruv_number_k_resolved)[1]
            ))
            log_info("Using 'ruv_number_k' from object args or default.")
        }

        # 3. Control Features (ctrl)
        if (!is.null(ctrl)) {
            ctrl_resolved <- ctrl
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ctrl from function arg, class = '%s', is.list = %s", class(ctrl_resolved)[1], is.list(ctrl_resolved)))
            log_info("Using 'ctrl' from function argument.")
        } else {
            ctrl_resolved <- checkParamsObjectFunctionSimplify(theObject, "ctrl", default_value = NULL)
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Using ctrl from object/default: is.null = %s, class = '%s', is.list = %s",
                is.null(ctrl_resolved), class(ctrl_resolved)[1], is.list(ctrl_resolved)
            ))
            log_info("Using 'ctrl' from object args or default.")
        }

        # Debug detailed ctrl info
        if (is.list(ctrl_resolved) && !is.null(ctrl_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ctrl_resolved is a list with names: %s", paste(names(ctrl_resolved), collapse = ", ")))
            for (nm in names(ctrl_resolved)) {
                ctrl_item <- ctrl_resolved[[nm]]
                num_ctrl <- if (is.logical(ctrl_item)) sum(ctrl_item, na.rm = TRUE) else length(ctrl_item)
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] ctrl_resolved[['%s']]: class = '%s', num_ctrl = %d", nm, class(ctrl_item)[1], num_ctrl))
            }
        }
        # Debug detailed k info
        if (is.list(ruv_number_k_resolved) && !is.null(ruv_number_k_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved is a list with names: %s", paste(names(ruv_number_k_resolved), collapse = ", ")))
            for (nm in names(ruv_number_k_resolved)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved[['%s']] = %s", nm, ruv_number_k_resolved[[nm]]))
            }
        } else if (!is.null(ruv_number_k_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved = %s (single value)", ruv_number_k_resolved))
        }

        # --- Update Object Args (Store the resolved values) ---
        # We store the *final* value used, regardless of source
        # Using direct slot assignment as updateParamInObject seemed problematic
        theObject@args$ruv_grouping_variable <- ruv_grouping_variable_final
        theObject@args$ruv_number_k <- ruv_number_k_resolved
        theObject@args$ctrl <- ctrl_resolved

        # --- Validation (Using the final resolved values) ---
        if (is.null(ruv_grouping_variable_final)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_grouping_variable_final is NULL!")
            log_error("Missing required parameter 'ruv_grouping_variable'. Must be provided in function call or object args.")
            stop("Missing required 'ruv_grouping_variable'")
        }
        if (!ruv_grouping_variable_final %in% colnames(design_matrix)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_grouping_variable '%s' not in design_matrix!", ruv_grouping_variable_final))
            log_error("Resolved 'ruv_grouping_variable' ('{ruv_grouping_variable_final}') not found as a column in the design matrix.", .logr = TRUE)
            stop("'ruv_grouping_variable' not in design matrix")
        }
        if (is.null(ruv_number_k_resolved)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_number_k_resolved is NULL!")
            log_error("Missing required parameter 'ruv_number_k' (K value). Must be provided in function call or object args.")
            stop("Missing required 'ruv_number_k'")
        }
        if (is.null(ctrl_resolved)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl_resolved is NULL!")
            log_error("Missing required parameter 'ctrl' (control features). Must be provided in function call or object args.")
            stop("Missing required 'ctrl'")
        }

        log_info("Starting RUV-III C Varying correction for metabolites.")
        log_info("Parameters (Resolved):")
        log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
        log_info("  - RUV Number K (k): Type '{class(ruv_number_k_resolved)}' (Value(s) resolved per assay)")
        log_info("  - Control Features (ctrl): Type '{class(ctrl_resolved)}' (Value(s) resolved per assay)")

        if (!is.list(assay_list)) {
            assay_list <- list(assay_list)
        }

        if (length(assay_list) == 0) {
            message("   DEBUG66 [ruvIII_C_Varying] WARNING - no assays found, returning object unchanged")
            log_warn("No assays found in `metabolite_data` slot. Returning object unchanged.")
            return(theObject)
        }
        # Ensure list is named
        assay_names <- names(assay_list)
        if (is.null(assay_names) || any(assay_names == "")) {
            needs_name <- which(is.null(assay_names) | assay_names == "")
            new_names <- paste0("Assay_", seq_along(assay_list))
            assay_names[needs_name] <- new_names[needs_name]
            names(assay_list) <- assay_names
            message("   DEBUG66 [ruvIII_C_Varying] Assay list had unnamed elements, using defaults")
            log_warn("Assay list contained unnamed or empty elements. Using default names (Assay_...).", immediate. = TRUE)
        }
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay names: %s", paste(assay_names, collapse = ", ")))

        # --- Validate structure of k and ctrl if they are lists ---
        is_k_list <- is.list(ruv_number_k_resolved)
        # Check if ctrl is a list, but NOT a data.frame (which could be passed accidentally)
        is_ctrl_list <- is.list(ctrl_resolved) && !is.data.frame(ctrl_resolved)
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] is_k_list = %s, is_ctrl_list = %s", is_k_list, is_ctrl_list))

        # Validate names if k is a list
        if (is_k_list && !all(assay_names %in% names(ruv_number_k_resolved))) {
            missing_k <- setdiff(assay_names, names(ruv_number_k_resolved))
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - k list missing assays: %s", paste(missing_k, collapse = ", ")))
            log_error("If 'ruv_number_k' is a list, its names must match assay names. Missing K for: {paste(missing_k, collapse=', ')}", .logr = TRUE)
            stop("Names in 'ruv_number_k' list do not match assay names.")
        }
        # Validate names if ctrl is a list
        if (is_ctrl_list && !all(assay_names %in% names(ctrl_resolved))) {
            missing_ctrl <- setdiff(assay_names, names(ctrl_resolved))
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl list missing assays: %s", paste(missing_ctrl, collapse = ", ")))
            log_error("If 'ctrl' is a list, its names must match assay names. Missing ctrl for: {paste(missing_ctrl, collapse=', ')}", .logr = TRUE)
            stop("Names in 'ctrl' list do not match assay names.")
        }
        # Validate type if k is NOT a list (must be single numeric for multiple assays)
        if (!is_k_list && length(assay_list) > 1 && !(is.numeric(ruv_number_k_resolved) && length(ruv_number_k_resolved) == 1 && !is.na(ruv_number_k_resolved))) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - k not a list but multiple assays, invalid format"))
            log_error("If multiple assays exist, 'ruv_number_k' must be a single non-NA numeric value or a named list.", .logr = TRUE)
            stop("Invalid format for 'ruv_number_k' for multiple assays.")
        }
        # Validate type if ctrl is NOT a list (must be vector for multiple assays)
        if (!is_ctrl_list && length(assay_list) > 1 && !(is.logical(ctrl_resolved) || is.numeric(ctrl_resolved) || is.character(ctrl_resolved))) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl not a list but multiple assays, invalid format"))
            log_error("If multiple assays exist and 'ctrl' is not a list, it must be a logical, numeric, or character vector (applied globally).", .logr = TRUE)
            stop("Invalid format for 'ctrl' for multiple assays.")
        }


        # --- Process Each Assay ---
        corrected_assay_list <- lapply(seq_along(assay_list), function(i) {
            assay_name <- assay_names[i]
            assay_tibble <- assay_list[[i]]
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] === Processing assay: %s ===", assay_name))
            message(sprintf("-- Processing assay for RUVIII: %s", assay_name))

            # --- Get Assay-Specific k and ctrl ---
            k_assay <- if (is_k_list) ruv_number_k_resolved[[assay_name]] else ruv_number_k_resolved
            ctrl_assay_input <- if (is_ctrl_list) ctrl_resolved[[assay_name]] else ctrl_resolved
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': k_assay = %s, ctrl_assay_input class = '%s'", assay_name, k_assay, class(ctrl_assay_input)[1]))
            if (is.logical(ctrl_assay_input)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': ctrl_assay_input sum(TRUE) = %d", assay_name, sum(ctrl_assay_input, na.rm = TRUE)))
            }

            # Validate k_assay
            if (!is.numeric(k_assay) || length(k_assay) != 1 || is.na(k_assay) || k_assay < 0) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - Invalid k_assay = %s", assay_name, k_assay))
                log_warn("Assay '{assay_name}': Invalid K value resolved ({k_assay}). Must be a non-negative integer. Skipping.", .logr = TRUE)
                return(NULL)
            }
            k_assay <- as.integer(k_assay) # Ensure integer

            # --- Basic Checks & Data Prep ---
            if (!tibble::is_tibble(assay_tibble)) {
                log_warn("Assay '{assay_name}' is not a tibble. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                log_warn("Assay '{assay_name}': ID column '{metabolite_id_col_name}' not found. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_tibble) < 1) {
                log_warn("Assay '{assay_name}' has no features. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Identify sample columns based on design matrix
            design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) {
                character(0)
            })
            if (length(design_samples) == 0) {
                log_warn("Assay '{assay_name}': No valid sample IDs in design matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) < 2) {
                log_warn("Assay '{assay_name}': Fewer than 2 sample columns found. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Ensure sample columns are numeric
            non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
            if (length(non_numeric_samples) > 0) {
                log_warn("Assay '{assay_name}': Coercing non-numeric sample columns to numeric: {paste(non_numeric_samples, collapse=', ')}", .logr = TRUE)
                assay_tibble <- assay_tibble |> dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # Convert to matrix (features x samples)
            # Handle potential duplicate feature IDs before converting to rownames
            assay_matrix <- tryCatch(
                {
                    n_initial <- nrow(assay_tibble)
                    assay_tibble_unique <- assay_tibble |>
                        dplyr::group_by(!!rlang::sym(metabolite_id_col_name)) |>
                        dplyr::filter(dplyr::row_number() == 1) |> # Keep only first instance
                        dplyr::ungroup()
                    n_final <- nrow(assay_tibble_unique)
                    if (n_final < n_initial) {
                        log_warn("Assay '{assay_name}': Duplicate feature IDs detected in '{metabolite_id_col_name}'. Keeping first instance only ({n_final}/{n_initial} features).", .logr = TRUE)
                    }
                    assay_tibble_unique |>
                        tibble::column_to_rownames(var = metabolite_id_col_name) |>
                        dplyr::select(dplyr::all_of(sample_cols)) |> # Select only sample columns
                        as.matrix()
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error converting to matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(NULL)
            }
            assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf AFTER conversion

            # Filter design matrix to match actual samples in matrix
            design_matrix_filtered <- design_matrix |>
                dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix)) |>
                as.data.frame() # Ensure data.frame

            if (nrow(design_matrix_filtered) < 2) {
                log_warn("Assay '{assay_name}': Fewer than 2 samples remain after filtering design matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_matrix) < 1) {
                log_warn("Assay '{assay_name}': Fewer than 1 feature remains. Skipping.", .logr = TRUE)
                return(NULL)
            }


            # --- Prepare Y (Samples x Features) --- NO IMPUTATION ---
            # Ensure column order matches filtered design matrix sample order
            Y_final <- t(assay_matrix[, as.character(design_matrix_filtered[[sample_id]]), drop = FALSE])

            # Check for NAs *before* RUV-III, as the helper might not handle them
            if (anyNA(Y_final)) {
                log_warn("   Assay '{assay_name}': Missing values (NA) detected in data matrix Y *before* RUV. RUVIII_C_Varying might fail or produce unexpected results if it doesn't handle NAs internally. Consider imputation *before* calling ruvIII_C_Varying if needed.", .logr = TRUE)
            }

            # Check dimensions after transpose
            if (nrow(Y_final) < 2 || ncol(Y_final) < 1) {
                log_warn("Assay '{assay_name}': Insufficient dimensions after preparing Y matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }


            # --- Prepare M Matrix (using filtered design matrix) ---
            M <- tryCatch(
                {
                    getRuvIIIReplicateMatrixHelper(
                        design_matrix_filtered,
                        !!rlang::sym(sample_id),
                        !!rlang::sym(ruv_grouping_variable_final)
                    )
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error getting RUV III Replicate Matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(M)) {
                return(NULL)
            } # Skip if M fails

            # Ensure M matrix dimensions match Y_final rows (samples)
            if (nrow(M) != nrow(Y_final) || !identical(rownames(M), rownames(Y_final))) {
                log_warn("Assay '{assay_name}': M matrix rownames do not match Y matrix rownames after filtering. Attempting to reorder.", .logr = TRUE)
                # Attempt to reorder M based on Y_final rownames if possible
                matched_m_rows <- match(rownames(Y_final), rownames(M))
                if (anyNA(matched_m_rows)) {
                    log_error("   Cannot reorder M matrix - rownames mismatch. Skipping.")
                    return(NULL)
                }
                M <- M[matched_m_rows, , drop = FALSE]
                log_info("   Reordered M matrix rows to match Y matrix.")
                if (nrow(M) != nrow(Y_final)) { # Double check after reorder
                    log_error("   M matrix row count still mismatch after reorder. Skipping.")
                    return(NULL)
                }
            }


            # --- Prepare potentialControls for RUVIII_C_Varying ---
            feature_names_in_assay <- colnames(Y_final) # Features present in Y_final

            # Resolve ctrl_assay_input into a logical vector aligned with feature_names_in_assay
            ctrl_logical_assay <- NULL
            if (is.null(ctrl_assay_input)) {
                log_warn("Assay '{assay_name}': Resolved control features ('ctrl') is NULL. Skipping.", .logr = TRUE)
                return(NULL)
            } else if (is.numeric(ctrl_assay_input)) {
                if (any(ctrl_assay_input < 1) || any(ctrl_assay_input > length(feature_names_in_assay))) {
                    log_warn("Assay '{assay_name}': Numeric 'ctrl' indices are out of bounds ({length(feature_names_in_assay)} features). Skipping.", .logr = TRUE)
                    return(NULL)
                }
                ctrl_logical_assay <- seq_along(feature_names_in_assay) %in% ctrl_assay_input
            } else if (is.logical(ctrl_assay_input)) {
                if (length(ctrl_assay_input) != length(feature_names_in_assay)) {
                    if (!is.null(names(ctrl_assay_input))) {
                        # Try to align based on names
                        feature_match <- match(feature_names_in_assay, names(ctrl_assay_input))
                        if (anyNA(feature_match)) {
                            log_warn("Assay '{assay_name}': Some assay features not found in named logical 'ctrl' vector. Skipping.", .logr = TRUE)
                            return(NULL)
                        }
                        ctrl_logical_assay <- ctrl_assay_input[feature_match]
                        if (length(ctrl_logical_assay) != length(feature_names_in_assay)) {
                            log_warn("Assay '{assay_name}': Length mismatch after aligning named logical 'ctrl' vector. Skipping.", .logr = TRUE)
                            return(NULL)
                        }
                        log_info("   Assay '{assay_name}': Aligned named logical 'ctrl' vector to assay features.", .logr = TRUE)
                    } else {
                        log_warn("Assay '{assay_name}': Unnamed logical 'ctrl' vector length ({length(ctrl_assay_input)}) does not match features ({length(feature_names_in_assay)}). Skipping.", .logr = TRUE)
                        return(NULL)
                    }
                } else {
                    ctrl_logical_assay <- ctrl_assay_input # Assume correct order
                }
            } else if (is.character(ctrl_assay_input)) {
                ctrl_logical_assay <- feature_names_in_assay %in% ctrl_assay_input
            } else {
                log_warn("Assay '{assay_name}': Invalid type for resolved 'ctrl' parameter. Expected numeric, logical, or character. Skipping.", .logr = TRUE)
                return(NULL)
            }

            if (is.null(ctrl_logical_assay)) {
                log_warn("Assay '{assay_name}': Failed to resolve control features to logical vector. Skipping.", .logr = TRUE)
                return(NULL)
            }
            num_controls_found <- sum(ctrl_logical_assay, na.rm = TRUE)
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': num_controls_found = %d", assay_name, num_controls_found))
            if (num_controls_found < 1) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - no control features", assay_name))
                log_warn("Assay '{assay_name}': No control features identified after resolution. Skipping.", .logr = TRUE)
                return(NULL)
            }
            log_info("   Assay '{assay_name}': Using {num_controls_found} control features.", .logr = TRUE)

            # Get the names of the control features
            potential_controls_names <- feature_names_in_assay[ctrl_logical_assay]
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': potential_controls_names count = %d", assay_name, length(potential_controls_names)))


            # --- Call RUVIII_C_Varying ---
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Assay '%s': Calling RUVIII_C_Varying with k = %d, Y dims = %dx%d, M dims = %dx%d",
                assay_name, k_assay, nrow(Y_final), ncol(Y_final), nrow(M), ncol(M)
            ))
            cln_mat <- tryCatch(
                {
                    # Check if RUVIII_C_Varying exists
                    if (!exists("RUVIII_C_Varying", mode = "function")) {
                        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - RUVIII_C_Varying function not found!", assay_name))
                        stop("Function 'RUVIII_C_Varying' not found. Ensure it is loaded from its package or source file.")
                    }
                    RUVIII_C_Varying(
                        k = k_assay,
                        Y = Y_final, # Use data potentially containing NAs
                        M = M,
                        toCorrect = colnames(Y_final), # Correct all features
                        potentialControls = potential_controls_names
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUVIII_C_Varying FAILED - %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error calling RUVIII_C_Varying: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(cln_mat) || !is.matrix(cln_mat)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - RUVIII_C_Varying returned NULL or non-matrix", assay_name))
                log_warn("Assay '{assay_name}': RUVIII_C_Varying did not return a valid matrix. Skipping.", .logr = TRUE)
                return(NULL) # Skip assay if RUV fails
            }
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUVIII_C_Varying SUCCESS - returned matrix dims = %dx%d", assay_name, nrow(cln_mat), ncol(cln_mat)))

            # --- Clean Corrected Matrix ---
            # Transpose result back to Features x Samples for cleaning
            corrected_matrix <- t(cln_mat)

            # Remove features (rows) with no finite values
            valid_features <- rowSums(is.finite(corrected_matrix), na.rm = TRUE) > 0
            corrected_matrix_filt_f <- corrected_matrix[valid_features, , drop = FALSE]
            if (nrow(corrected_matrix_filt_f) == 0) {
                log_warn("Assay '{assay_name}': No features remained after removing non-finite rows post-RUV. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Remove samples (columns) with no finite values
            valid_samples <- colSums(is.finite(corrected_matrix_filt_f), na.rm = TRUE) > 0
            corrected_matrix_filt_fs <- corrected_matrix_filt_f[, valid_samples, drop = FALSE]
            if (ncol(corrected_matrix_filt_fs) == 0) {
                log_warn("Assay '{assay_name}': No samples remained after removing non-finite columns post-RUV. Skipping.", .logr = TRUE)
                return(NULL)
            }

            log_info("   Assay '{assay_name}': RUV correction applied. Dimensions before cleaning: {nrow(corrected_matrix)}x{ncol(corrected_matrix)}, After: {nrow(corrected_matrix_filt_fs)}x{ncol(corrected_matrix_filt_fs)}", .logr = TRUE)

            # --- Reconstruct Tibble ---
            # Get original metadata columns relevant to the remaining features
            metadata_cols <- setdiff(colnames(assay_tibble), sample_cols) # All original non-sample columns
            # Filter the *original* tibble to get metadata for rows that remain
            original_metadata_tibble <- assay_tibble |>
                dplyr::filter(!!rlang::sym(metabolite_id_col_name) %in% rownames(corrected_matrix_filt_fs)) |>
                dplyr::select(dplyr::all_of(c(metabolite_id_col_name, metadata_cols))) # Ensure ID column is selected

            # Ensure metadata IDs are unique before join (should be due to earlier handling, but safe)
            original_metadata_tibble <- original_metadata_tibble |>
                dplyr::distinct(!!rlang::sym(metabolite_id_col_name), .keep_all = TRUE)


            reconstructed_tibble <- tryCatch(
                {
                    corrected_data_tibble <- corrected_matrix_filt_fs |>
                        as.data.frame() |>
                        tibble::rownames_to_column(var = metabolite_id_col_name) |>
                        tibble::as_tibble()

                    # Join corrected data with filtered original metadata
                    # Ensure join column types match (rownames_to_column is character)
                    original_metadata_tibble_char <- original_metadata_tibble |>
                        dplyr::mutate(!!rlang::sym(metabolite_id_col_name) := as.character(!!rlang::sym(metabolite_id_col_name)))
                    corrected_data_tibble_char <- corrected_data_tibble |>
                        dplyr::mutate(!!rlang::sym(metabolite_id_col_name) := as.character(!!rlang::sym(metabolite_id_col_name)))

                    final_tibble <- dplyr::left_join(original_metadata_tibble_char, corrected_data_tibble_char, by = metabolite_id_col_name) |>
                        # Ensure original column order (metadata first, then remaining samples)
                        dplyr::relocate(
                            dplyr::all_of(colnames(original_metadata_tibble)), # All metadata cols
                            dplyr::all_of(colnames(corrected_matrix_filt_fs))
                        ) # Remaining sample cols

                    # Check if join resulted in expected columns
                    if (!identical(sort(colnames(final_tibble)), sort(c(colnames(original_metadata_tibble), colnames(corrected_matrix_filt_fs))))) {
                        log_warn("Assay '{assay_name}': Column mismatch after joining corrected data and metadata.", .logr = TRUE)
                        # Potentially return NULL or the corrected_data_tibble only
                    }
                    final_tibble
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error reconstructing tibble after RUV correction: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL) # Return NULL on error
                }
            )

            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUV-III correction complete.", assay_name))
            message(sprintf("   Assay '%s' RUV-III correction complete.", assay_name))
            return(reconstructed_tibble)
        })

        # Set names and remove NULLs
        names(corrected_assay_list) <- assay_names
        final_corrected_list <- corrected_assay_list[!sapply(corrected_assay_list, is.null)]
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] final_corrected_list length = %d", length(final_corrected_list)))

        if (length(final_corrected_list) == 0) {
            message("   DEBUG66 [ruvIII_C_Varying] WARNING - no assays successfully processed, returning original object")
            log_warn("No assays were successfully processed by RUV-III. Returning original object.")
            return(theObject)
        }

        # Update the slot in the object
        methods::slot(theObject, "metabolite_data") <- final_corrected_list

        # --- Clean Design Matrix ---
        theObject <- tryCatch(
            {
                log_info("Cleaning design matrix to match remaining samples after RUV...")
                cleanDesignMatrix(theObject)
            },
            error = function(e) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] cleanDesignMatrix FAILED - %s", e$message))
                log_warn("Error running cleanDesignMatrix after RUV correction: {e$message}. Design matrix might not be fully synchronized.", .logr = TRUE)
                return(theObject)
            }
        )

        message(sprintf("   DEBUG66 [ruvIII_C_Varying] RUV-III correction finished for %d assay(s).", length(final_corrected_list)))
        message("+===========================================================================+")
        message("|  DEBUG66: Exiting ruvIII_C_Varying                                        |")
        message("+===========================================================================+")
        log_info("RUV-III correction process finished for {length(final_corrected_list)} assay(s).")
        return(theObject)
    }
)

