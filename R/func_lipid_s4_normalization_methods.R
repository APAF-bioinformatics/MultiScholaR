#' @title Normalize by Internal Standard (ITSD) for LipidomicsAssayData
#'
#' @description
#' Normalizes the untransformed lipid intensity data in each assay of a
#' `LipidomicsAssayData` object using Internal Standards (ITSDs).
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
#' log-transformed) data in the `lipid_data` slot. It aims to preserve the
#' overall scale of the data while adjusting relative differences between samples,
#' centering the normalized data around the average ITSD response observed.
#'
#' Assays where no ITSDs are found will issue a warning and remain unchanged. Samples
#' where the aggregate ITSD signal is zero or NA will also result in warnings,
#' and the normalized values for those samples will likely become NA.
#'
#' @describeIn normaliseUntransformedData Method for LipidomicsAssayData using ITSD
#'
#' @param theObject A `LipidomicsAssayData` object.
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
#' @return An updated `LipidomicsAssayData` object with ITSD-normalized values in the
#'   `lipid_data` slot. Normalization status is recorded in the `args` slot.
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
    signature = signature(theObject = "LipidomicsAssayData", method = "character"),
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
        assay_list <- methods::slot(theObject, "lipid_data")
        lipid_id_col <- methods::slot(theObject, "lipid_id_column")
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
            warning("LipidomicsAssayData object has no assays in 'lipid_data' slot. No normalization performed.")
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
        # KIT-242 (VAL-242 OD-2 (b), OD-3 (d)): the assays not normalised to standards, with the reason, and the samples
        # without standard signal, per assay.
        itsd_collector$not_applied <- list()
        # KIT-242 (OD-4 (c), the coordinator's ruling on N05): a record made before replays with the earlier ITSD record
        # (applied always TRUE, no per-assay fields), so that it replays exactly (`workflowReplayRuleMarkers$itsd_outcome`).
        itsd_outcome_marker <- workflowReplayRuleMarkers$itsd_outcome
        itsd_current <- !identical(getOption(itsd_outcome_marker$option), itsd_outcome_marker$earlier)
        itsd_collector$samples_without <- list()

        # --- Process Each Assay ---
        normalized_assay_list <- lapply(seq_along(assay_list), function(i) {
            assay_index_name <- assay_names[i] # Use the name from the corrected list
            assay_tibble <- assay_list[[i]]
            message(sprintf("-- Processing assay: %s", assay_index_name))
            # KIT-242 (VAL-242 OD-2 (b)): no standard in this assay: not normalised to standards, recorded, with a notice.
            notApplied <- function(reason, text) {
                if (!itsd_current) return(warning(text, immediate. = TRUE))
                itsd_collector$not_applied[[assay_index_name]] <- reason
                itsd_collector$counts_per_assay[[assay_index_name]] <- 0L
                rlang::warn(text, class = "omics_itsd_not_applied", assay = assay_index_name, reason = reason)
            }
            # KIT-242 (VAL-242: a method failure or a malformed assay stops).
            itsdFailed <- function(reason, text, parent = NULL) {
                if (!itsd_current) return(warning(text, immediate. = TRUE))
                rlang::abort(text, class = "omics_itsd_failed", assay = assay_index_name, reason = reason, parent = parent,
                    call = NULL)
            }

            # --- Basic Checks ---
            if (!tibble::is_tibble(assay_tibble)) {
                warning(sprintf("Assay '%s' is not a tibble. Attempting to coerce.", assay_index_name), immediate. = TRUE)
                assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                    itsdFailed("malformed_assay", sprintf("Failed to coerce assay '%s' to tibble: %s. Skipping normalization.", assay_index_name, e$message))
                    return(NULL) # Signal to skip this assay
                })
                if (is.null(assay_tibble)) {
                    return(assay_list[[i]])
                } # Return original if coercion failed
            }
            if (!lipid_id_col %in% colnames(assay_tibble)) {
                itsdFailed("malformed_assay", sprintf("Assay '%s': Lipid ID column '%s' not found. Skipping normalization.", assay_index_name, lipid_id_col))
                return(assay_tibble)
            }
            # --- Validate ITSD Identification Prerequisites ---
            # Only check pattern columns if using regex mode
            actual_itsd_cols <- NULL
            if (!use_manual_itsd) {
                # Check if *any* of the itsd_pattern_columns exist
                if (!any(itsd_pattern_columns %in% colnames(assay_tibble))) {
                    notApplied("no_pattern_column", sprintf(
                        "Assay '%s': None of the specified ITSD pattern columns ('%s') found. Cannot identify ITSDs. Skipping normalization.",
                        assay_index_name, paste(itsd_pattern_columns, collapse = "', '")
                    ))
                    return(assay_tibble)
                }
                # Filter itsd_pattern_columns to only those present in the current assay
                actual_itsd_cols <- intersect(itsd_pattern_columns, colnames(assay_tibble))
                if (length(actual_itsd_cols) == 0) { # Should be caught above, but double check
                    notApplied("no_pattern_column", sprintf("Assay '%s': No ITSD identification columns found after checking existence. Skipping normalization.", assay_index_name))
                    return(assay_tibble)
                }
            }


            # --- Identify Sample Columns in this Assay ---
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) == 0) {
                itsdFailed("malformed_assay", sprintf("Assay '%s': No sample columns identified matching design matrix sample IDs. Skipping normalization.", assay_index_name))
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
                    notApplied("no_standard_selected", sprintf("Assay '%s': No manual ITSD features provided in itsd_feature_ids. Skipping normalization.", assay_index_name))
                    return(assay_tibble)
                }
                # Match feature IDs against lipid_id_col
                feature_ids <- as.character(assay_tibble[[lipid_id_col]])
                itsd_rows_logical <- feature_ids %in% as.character(assay_itsd_ids)

                if (!any(itsd_rows_logical, na.rm = TRUE)) {
                    notApplied("no_standard_matched", sprintf(
                        "Assay '%s': None of the %d manual ITSD feature IDs matched rows in the data. Skipping normalization.",
                        assay_index_name, length(assay_itsd_ids)
                    ))
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
                    notApplied("no_standard_matched", sprintf(
                        "Assay '%s': No ITSD features identified using regex '%s' in columns '%s'. Skipping normalization.",
                        assay_index_name, itsd_regex, paste(actual_itsd_cols, collapse = "', '")
                    ))
                    return(assay_tibble)
                }
                message(sprintf("   Identified %d ITSD features via regex.", sum(itsd_rows_logical)))
            }

            itsd_data <- assay_tibble |> dplyr::filter(itsd_rows_logical)
            non_itsd_data <- assay_tibble |> dplyr::filter(!itsd_rows_logical)
            n_itsd <- nrow(itsd_data)
            message(sprintf("   Using %d ITSD features for normalization.", n_itsd))

            # --- Capture ITSD feature names for reporting ---
            itsd_feature_names <- as.character(itsd_data[[lipid_id_col]])
            itsd_collector$features_per_assay[[assay_index_name]] <- itsd_feature_names
            itsd_collector$counts_per_assay[[assay_index_name]] <- n_itsd

            # --- Calculate Normalization Factors ---
            # Step 1: Calculate aggregate ITSD per sample
            norm_factors_long <- tryCatch(
                {
                    itsd_data |>
                        dplyr::select(dplyr::all_of(c(lipid_id_col, sample_cols))) |>
                        tidyr::pivot_longer(cols = dplyr::all_of(sample_cols), names_to = "Sample", values_to = "Intensity") |>
                        dplyr::group_by(.data$Sample) |>
                        # Aggregate ITSD intensities per sample
                        # Use SampleNormFactor to be clear about the value's meaning
                        dplyr::summarise(SampleNormFactor = itsd_aggregation_func(.data$Intensity, na.rm = TRUE), .groups = "drop")
                },
                error = function(e) {
                    itsdFailed("factor_error", sprintf("Assay '%s': Error calculating normalization factors: %s. Skipping normalization.", assay_index_name, e$message))
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
                itsd_collector$samples_without[[assay_index_name]] <- as.character(problematic_factors$Sample)
                # KIT-242 (VAL-242 OD-3 (d)): the sample is set missing in this assay only (the between-sample step then
                # leaves it out of this assay), recorded per assay and sample, with a notice.
                itsd_text <- sprintf(
                    "Assay '%s': Aggregate ITSD signal (SampleNormFactor) is NA or zero for samples: %s. Normalized values will be NA for these samples.",
                    assay_index_name, paste(problematic_factors$Sample, collapse = ", ")
                )
                if (itsd_current) {
                    rlang::warn(itsd_text, class = "omics_itsd_sample_without_standard", assay = assay_index_name,
                        sample = as.character(problematic_factors$Sample))
                } else {
                    warning(itsd_text, immediate. = TRUE)
                }
                # Set factor to NA to ensure division results in NA
                norm_factors_long <- norm_factors_long |>
                    dplyr::mutate(SampleNormFactor = dplyr::if_else(is.na(.data$SampleNormFactor) | .data$SampleNormFactor == 0, NA_real_, .data$SampleNormFactor))
            }

            # Step 2: Calculate the average aggregate ITSD signal across all samples
            average_norm_factor <- mean(norm_factors_long$SampleNormFactor, na.rm = TRUE)
            if (is.na(average_norm_factor) || average_norm_factor == 0) {
                itsdFailed("no_standard_signal", sprintf("Assay '%s': Average aggregate ITSD signal is NA or zero. Cannot perform average-centered normalization. Skipping.", assay_index_name))
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
                    itsdFailed("apply_error", sprintf("Assay '%s': Error applying normalization factors: %s. Returning unnormalized data for non-ITSD features.", assay_index_name, e$message))
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
                    dplyr::arrange(!!rlang::sym(lipid_id_col))
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
        methods::slot(theObject, "lipid_data") <- normalized_assay_list

        # --- Update Args Slot ---
        if (!is.list(theObject@args)) {
            warning("Slot 'args' is not a list. Cannot record ITSD normalization status.", immediate. = TRUE)
        } else {
            # Ensure the specific list exists
            if (!"ITSDNormalization" %in% names(theObject@args)) {
                theObject@args$ITSDNormalization <- list()
            }
            # KIT-242 (VAL-242 OD-2 (b)): applied only when at least one assay was normalised to its standards; per assay
            # the outcome, the reason where not applied, and the samples without standard signal (current rule only).
            if (itsd_current) {
                theObject@args$ITSDNormalization$applied <- length(itsd_collector$features_per_assay) > 0
                theObject@args$ITSDNormalization$applied_per_assay <- as.list(stats::setNames(
                    assay_names %in% names(itsd_collector$features_per_assay), assay_names))
                theObject@args$ITSDNormalization$not_applied_reason_per_assay <- as.list(itsd_collector$not_applied)
                theObject@args$ITSDNormalization$samples_without_standard_per_assay <- as.list(itsd_collector$samples_without)
            } else {
                theObject@args$ITSDNormalization$applied <- TRUE
            }
            theObject@args$ITSDNormalization$method_type <- "average_centered"
            theObject@args$ITSDNormalization$itsd_aggregation <- tolower(itsd_aggregation)
            theObject@args$ITSDNormalization$itsd_pattern_columns <- itsd_pattern_columns # Record actual columns used (potentially default)
            theObject@args$ITSDNormalization$removed_itsd <- remove_itsd_after_norm
            theObject@args$ITSDNormalization$timestamp <- Sys.time()

            # Store per-assay ITSD feature information for reporting
            if (length(itsd_collector$features_per_assay) > 0) {
                theObject@args$ITSDNormalization$itsd_features_per_assay <- as.list(itsd_collector$features_per_assay)
                if (!itsd_current) theObject@args$ITSDNormalization$itsd_counts_per_assay <- as.list(itsd_collector$counts_per_assay)
                message(sprintf("   Stored ITSD feature names for %d assays", length(itsd_collector$features_per_assay)))
            }
            # KIT-242: the count of standards used per assay, 0 where none matched (the count of standards matched, not of
            # standards observed in each sample).
            if (itsd_current && length(itsd_collector$counts_per_assay) > 0) {
                theObject@args$ITSDNormalization$itsd_counts_per_assay <- as.list(itsd_collector$counts_per_assay)[
                    intersect(assay_names, names(itsd_collector$counts_per_assay))]
            }
        }

        message("ITSD normalization process complete for all applicable assays.")
        return(theObject)
    }
)


#' @title Clean Design Matrix for LipidomicsAssayData
#' @name cleanDesignMatrix,LipidomicsAssayData-method
#' @importFrom dplyr inner_join select rename filter all_of any_of
#' @importFrom rlang sym !!
#' @importFrom methods slot
#' @importFrom tibble tibble
#' @export
#' @param theObject Runtime inputs used by this function; see the usage section for accepted values.
setMethod(
    f = "cleanDesignMatrix",
    signature = "LipidomicsAssayData",
    definition = function(theObject) {
        assay_list <- methods::slot(theObject, "lipid_data")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col_name <- methods::slot(theObject, "sample_id")
        lipid_id_col_name <- methods::slot(theObject, "lipid_id_column") # Needed to exclude from sample cols

        if (length(assay_list) == 0) {
            warning("cleanDesignMatrix: No assays found in `lipid_data`. Returning object unchanged.")
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


#' @title Log2 Transform Assay Data for LipidomicsAssayData
#'
#' @description
#' Applies a log2 transformation to the numeric sample columns in all assay tibbles stored within the `lipid_data`
#' slot of a `LipidomicsAssayData` object: plain `log2(x)` by default, `log2(x + offset)` for an offset above 0.
#'
#' @details
#' At the default offset of 0 a missing value stays missing, and a value of 0 or below or a value that is not finite
#' (`NaN`, `Inf`, `-Inf`) becomes missing; these are counted per assay in `args$log_transform_nonpositive`
#' (non-positive and non-finite) and reported in one warning of class `omics_log2_nonpositive`. The import already
#' reads a zero as not detected and refuses negative values, so in the application every value reaching this step
#' is positive. An offset above 0 runs the transformation of earlier releases unchanged (negative values set to 0,
#' then `log2(x + offset)`), so a record made with an offset replays exactly; the result then depends on the unit of
#' the table (a value v moves by log2(1 + offset / v)). An offset below 0, missing or not a single number is refused
#' with class `omics_log2_offset_invalid`; an error of the transformation stops with class `omics_log2_failed`.
#' Non-numeric columns and the lipid ID column are preserved.
#'
#' @title Log Transform Assays for LipidomicsAssayData
#' @name logTransformAssays,LipidomicsAssayData-method
#' @param theObject A `LipidomicsAssayData` object.
#' @param offset A number of 0 or above added before the log transformation (default: 0, plain log2).
#' @param ... Currently unused.
#'
#' @return An updated `LipidomicsAssayData` object with log2 transformed values in the `lipid_data` slot.
#'
#' @importFrom dplyr mutate across all_of select filter intersect union setdiff
#' @importFrom rlang sym !!
#' @importFrom purrr map set_names
#' @importFrom methods slot slot<- is
#' @importFrom tibble as_tibble is_tibble
#' @export
setMethod(
    f = "logTransformAssays",
    signature = "LipidomicsAssayData",
    definition = function(theObject, offset = 0, ...) {
        # --- Input Validation ---
        omicsLog2CheckOffset(offset)
        if (offset > 0) {
            message(sprintf("Applying log2(x + %s) transformation to assays.", as.character(offset)))
        } else {
            message("Applying log2(x) transformation to assays.")
        }
        tally <- new.env(parent = emptyenv())
        tally$counts <- list()

        # --- Get Object Slots ---
        assay_list <- methods::slot(theObject, "lipid_data")
        lipid_id_col_name <- methods::slot(theObject, "lipid_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col_name <- methods::slot(theObject, "sample_id")

        if (length(assay_list) == 0) {
            warning("LipidomicsAssayData object has no assays in 'lipid_data' slot. No transformation performed.")
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

            # Check for lipid ID column
            if (!lipid_id_col_name %in% colnames(assay_tibble)) {
                warning(sprintf("Assay '%s': Lipid ID column '%s' not found. Skipping transformation.", assay_index_name, lipid_id_col_name), immediate. = TRUE)
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
            # Ensure lipid ID is not treated as a sample column
            metadata_cols <- union(metadata_cols, lipid_id_col_name)
            sample_cols <- setdiff(all_assay_cols, metadata_cols)
            # ----------------------------- #

            if (length(sample_cols) == 0) {
                warning(sprintf("Assay '%s': No numeric sample columns identified matching design matrix sample IDs. Skipping transformation.", assay_index_name), immediate. = TRUE)
                return(assay_tibble) # Return original
            }

            # --- Apply Log2 Transformation (KIT-172: the shared step; an error stops) ---
            step <- omicsLog2SampleColumns(assay_tibble, sample_cols, offset, assay_index_name)
            transformed_tibble <- step$assay
            if (!is.null(step$counts)) tally$counts[[assay_index_name]] <- step$counts

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
        methods::slot(theObject, "lipid_data") <- transformed_assay_list

        # Add log transformation status to args (KIT-172: at offset 0 also the counts, and the warning)
        theObject <- omicsLog2RecordArgs(theObject, offset, tally$counts)


        message("Log2 transformation complete for all assays.")
        return(theObject)
    }
)


#' @title Normalize Between Samples for LipidomicsAssayData
#' @name normaliseBetweenSamples,LipidomicsAssayData-method
#' @param theObject Object of class LipidomicsAssayData
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
    signature = "LipidomicsAssayData",
    definition = function(theObject, normalisation_method = NULL) {
        assay_list <- methods::slot(theObject, "lipid_data")
        lipid_id_col_name <- methods::slot(theObject, "lipid_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id_col_name <- methods::slot(theObject, "sample_id")

        # --- Get Normalization Method ---
        # Use the general parameter name as in protein version for consistency
        normalisation_method_final <- checkParamsObjectFunctionSimplify(
            theObject,
            "normalisation_method",
            default_value = "cyclicloess" # Default if not found in args or user override
        )
        # Store the *actually used* method back into args
        theObject@args$normalisation_method <- normalisation_method_final
        log_info("Applying between-sample normalization method: {normalisation_method_final}")
        # KIT-242 (VAL-242): the rule this run is made under (the earlier one only when replaying an older record).
        current_rule <- omicsNormalisationCurrentRule()
        norm_collector <- new.env(parent = emptyenv())
        norm_collector$excluded <- list()


        if (length(assay_list) == 0) {
            if (current_rule) abortOmicsNormalisation(NA_character_, "malformed_assay", normalisation_method_final,
                "Between-sample normalisation was not run: the object holds no assay.")
            warning("No assays found in `lipid_data` slot. Skipping normalization.")
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
            # KIT-242 (VAL-242, UVR-082): a malformed assay stops under the current rule; earlier, it was skipped.
            malformed <- function(text) {
                if (current_rule) abortOmicsNormalisation(assay_name, "malformed_assay", normalisation_method_final, text)
                warning(text, immediate. = TRUE)
            }

            # --- Basic Checks ---
            if (!tibble::is_tibble(assay_tibble)) {
                warning(sprintf("Assay '%s' is not a tibble. Attempting coercion.", assay_name), immediate. = TRUE)
                assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                    malformed(sprintf("Failed to coerce assay '%s' to tibble: %s. Skipping normalization.", assay_name, e$message))
                    return(NULL) # Signal to skip
                })
                if (is.null(assay_tibble)) {
                    return(assay_list[[i]])
                } # Return original if coercion failed
            }
            if (!lipid_id_col_name %in% colnames(assay_tibble)) {
                malformed(sprintf("Assay '%s': Lipid ID column '%s' not found. Skipping normalization.", assay_name, lipid_id_col_name))
                return(assay_tibble)
            }
            # KIT-242 (VAL-242: a malformed assay stops): duplicate feature identifiers (the QC tab's duplicate resolution
            # not run) stop with the cause and what to do; earlier, the matrix conversion failed and the assay was skipped.
            if (current_rule) {
                feature_ids <- as.character(assay_tibble[[lipid_id_col_name]])
                duplicated_ids <- unique(feature_ids[duplicated(feature_ids)])
                if (length(duplicated_ids)) {
                    abortOmicsNormalisation(assay_name, "malformed_assay", normalisation_method_final, sprintf(paste(
                        "Between-sample normalisation was not run: assay '%s' has %d duplicated feature identifiers (%s).",
                        "Run 'Resolve Duplicates' in the Quality Control tab, then run the normalisation again."), assay_name,
                        length(duplicated_ids), paste(utils::head(duplicated_ids, 3L), collapse = ", ")))
                }
            }

            # --- Identify Sample Columns ---
            design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) {
                character(0)
            })
            if (length(design_samples) == 0) {
                malformed(sprintf("Assay '%s': Could not extract valid sample IDs from design matrix column '%s'. Skipping normalization.", assay_name, sample_id_col_name))
                return(assay_tibble)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) == 0) {
                malformed(sprintf("Assay '%s': No sample columns identified matching design matrix sample IDs. Skipping normalization.", assay_name))
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
                        dplyr::select(dplyr::all_of(c(lipid_id_col_name, sample_cols))) |> # Select ID + Samples
                        tibble::column_to_rownames(var = lipid_id_col_name) |>
                        as.matrix()
                },
                error = function(e) {
                    malformed(sprintf("Assay '%s': Error converting tibble to matrix: %s. Skipping normalization.", assay_name, e$message))
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(assay_tibble)
            } # Return original if matrix conversion failed

            assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf

            if (current_rule) {
                # KIT-242 (VAL-242 OD-1 (c)): a sample with no value is left out of the method and kept empty in this
                # assay only; KIT-145's guard counts the other samples (item viii); a method failure stops.
                if (nrow(assay_matrix) < 1 || ncol(assay_matrix) < 1) {
                    malformed(sprintf("Assay '%s': Matrix is empty after preparation. Skipping normalization.", assay_name))
                }
                normalised <- normaliseOmicsAssayMatrix(assay_matrix, assay_name, normalisation_method_final, "lipids")
                normalized_matrix <- normalised$matrix
                if (length(normalised$excluded)) norm_collector$excluded[[assay_name]] <- normalised$excluded
            } else {
                # KIT-145 (M17): cyclic loess is refused below its minimum number of features. The refusal is raised here,
                # before the skips below and outside the tryCatch around the fit, so that it is never turned into
                # unnormalised data with a warning.
                if (identical(normalisation_method_final, "cyclicloess")) {
                    checkOmicsCyclicLoessMinimum(assay_matrix, assay_name, "lipids")
                }

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
            }

            normalized_matrix[!is.finite(normalized_matrix)] <- NA # Ensure NAs remain NAs

            # --- Reconstruct Tibble ---
            # Get original metadata columns
            metadata_cols <- setdiff(colnames(assay_tibble), sample_cols)

            reconstructed_tibble <- tryCatch(
                {
                    normalized_data_tibble <- normalized_matrix |>
                        as.data.frame() |> # Convert matrix to data frame
                        tibble::rownames_to_column(var = lipid_id_col_name) |>
                        tibble::as_tibble()

                    original_metadata_tibble <- assay_tibble |>
                        dplyr::select(dplyr::any_of(metadata_cols)) # Use any_of in case some metadata cols were dynamic

                    # Ensure join column types match (rownames_to_column creates character)
                    original_metadata_tibble_char <- original_metadata_tibble |>
                        dplyr::mutate(!!rlang::sym(lipid_id_col_name) := as.character(!!rlang::sym(lipid_id_col_name)))
                    normalized_data_tibble_char <- normalized_data_tibble |>
                        dplyr::mutate(!!rlang::sym(lipid_id_col_name) := as.character(!!rlang::sym(lipid_id_col_name)))

                    # Join normalized data with original metadata
                    dplyr::left_join(original_metadata_tibble_char, normalized_data_tibble_char, by = lipid_id_col_name) |>
                        # Ensure original column order (metadata first, then samples in original order)
                        dplyr::relocate(dplyr::all_of(metadata_cols), dplyr::all_of(sample_cols))
                },
                error = function(e) {
                    malformed(sprintf("Assay '%s': Error reconstructing tibble after normalization: %s. Returning original data.", assay_name, e$message))
                    return(assay_tibble) # Return original on error
                }
            )

            message(sprintf("   Assay '%s' normalization complete.", assay_name))
            return(reconstructed_tibble)
        })

        # Restore original names
        names(normalized_assay_list) <- assay_names

        # KIT-242 (item 1a): a sample with no value in any assay leaves the design, with a notice.
        if (current_rule && length(norm_collector$excluded)) {
            pruned <- omicsRemoveSamplesEmptyEverywhere(normalized_assay_list, theObject@design_matrix, sample_id_col_name,
                norm_collector$excluded)
            normalized_assay_list <- pruned$assays
            theObject@design_matrix <- pruned$design
        }

        # Update the slot in the object
        methods::slot(theObject, "lipid_data") <- normalized_assay_list

        # --- Clean Design Matrix (as done in protein version) ---
        # Ensure the cleanDesignMatrix method exists for LipidomicsAssayData
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

#' @title Lipid Intensity Filtering Method for LipidomicsAssayData
#'
#' @description
#' Filters lipids in *all* assays of a LipidomicsAssayData object.
#' It removes lipids that have intensities below a certain percentile threshold
#' in a proportion of the design's samples exceeding a defined cutoff. The threshold is calculated
#' independently for each assay.
#'
#' @describeIn lipidIntensityFiltering Method for LipidomicsAssayData
#'
#' @param theObject A LipidomicsAssayData object.
#' @param lipids_intensity_cutoff_percentile See generic definition.
#' @param lipids_proportion_of_samples_below_cutoff See generic definition.
#'
#' @importFrom dplyr pull select all_of across
#' @importFrom rlang sym
#' @importFrom stats quantile
#'
#' @return An updated LipidomicsAssayData object.
#' @export
setMethod(
    f = "lipidIntensityFiltering",
    signature = "LipidomicsAssayData",
    definition = function(theObject, lipids_intensity_cutoff_percentile = NULL, lipids_proportion_of_samples_below_cutoff = NULL) {
        # --- Parameter Resolution (Done once) ---
        config_intensity_percentile <- "lipids_intensity_cutoff_percentile"
        raw_intensity_percentile <- checkParamsObjectFunctionSimplify(
            theObject,
            config_intensity_percentile,
            lipids_intensity_cutoff_percentile
        )
        message("Raw intensity percentile from config/param: ", raw_intensity_percentile)
        cleaned_intensity_percentile <- trimws(sub("#.*$", "", raw_intensity_percentile))
        intensity_cutoff_percentile_final <- as.numeric(cleaned_intensity_percentile)

        config_proportion_cutoff <- "lipids_proportion_of_samples_below_cutoff"
        raw_proportion_cutoff <- checkParamsObjectFunctionSimplify(
            theObject,
            config_proportion_cutoff,
            lipids_proportion_of_samples_below_cutoff
        )
        message("Raw proportion cutoff from config/param: ", raw_proportion_cutoff)
        cleaned_proportion_cutoff <- trimws(sub("#.*$", "", raw_proportion_cutoff))
        proportion_of_samples_below_cutoff_final <- as.numeric(cleaned_proportion_cutoff)

        if (is.na(intensity_cutoff_percentile_final)) {
            stop("Failed to convert cleaned lipids_intensity_cutoff_percentile ('", cleaned_intensity_percentile, "' from raw '", raw_intensity_percentile, "') to numeric. Check config.ini or parameter value.")
        }
        if (is.na(proportion_of_samples_below_cutoff_final)) {
            stop("Failed to convert cleaned lipids_proportion_of_samples_below_cutoff ('", cleaned_proportion_cutoff, "' from raw '", raw_proportion_cutoff, "') to numeric. Check config.ini or parameter value.")
        }

        # --- Update Object Parameters (Done once) ---
        theObject <- updateParamInObject(theObject, config_intensity_percentile)
        theObject <- updateParamInObject(theObject, config_proportion_cutoff)

        # --- Process Each Assay in the List ---
        lipid_id_col <- theObject@lipid_id_column
        original_assay_list <- theObject@lipid_data
        original_assay_names <- names(original_assay_list)

        if (length(original_assay_list) == 0) {
            warning("LipidomicsAssayData object has no assays in 'lipid_data' slot. No filtering performed.")
            return(theObject)
        }

        design_samples <- as.character(theObject@design_matrix[[theObject@sample_id]])
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

            if (!(lipid_id_col %in% names(assay_table))) {
                warning("Lipid ID column '", lipid_id_col, "' not found in assay '", assay_name_for_msg, "'. Skipping this assay.")
                return(assay_table) # Return the original table if ID is missing
            }

            # KIT-202 (VAL-118 (1)): the design's samples of this assay, less the id column (every numeric column under
            # the earlier rule of a record made before); an absent or non-numeric design sample is refused.
            sample_cols <- omicsIntensityFilterSampleColumns(names(assay_table),
                names(assay_table)[sapply(assay_table, is.numeric)], design_samples, lipid_id_col, assay_name_for_msg)

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

            # Calculate threshold specifically for this assay. KIT-031 (user decision 2026-10-28): the threshold is the
            # quantile itself; the ceiling() that rounded it up to the next integer intensity is removed.
            min_lipid_intensity_threshold <- quantile(all_intensity_values,
                na.rm = TRUE,
                probs = c(intensity_cutoff_percentile_final / 100)
            )[1]

            message("Calculated minimum intensity threshold for assay '", assay_name_for_msg, "': ", min_lipid_intensity_threshold)

            # Filter using Helper
            filtered_assay <- lipidIntensityFilteringHelper(
                assay_table = assay_table,
                min_lipid_intensity_threshold = min_lipid_intensity_threshold,
                lipids_proportion_of_samples_below_cutoff = proportion_of_samples_below_cutoff_final,
                lipid_id_column = lipid_id_col,
                sample_columns = sample_cols
            )

            message("Filtered assay '", assay_name_for_msg, "'. Original rows: ", nrow(assay_table), ", Filtered rows: ", nrow(filtered_assay))
            return(filtered_assay)
        })

        # Restore original names if they existed
        if (!is.null(original_assay_names)) {
            names(filtered_assay_list) <- original_assay_names
        }

        # Assign the list of filtered assays back to the object
        theObject@lipid_data <- filtered_assay_list

        # Optional: Call a generic cleanup/design matrix function if applicable
        # theObject <- cleanDesignMatrix(theObject) # If a generic method exists

        return(theObject)
    }
)
