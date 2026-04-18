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

