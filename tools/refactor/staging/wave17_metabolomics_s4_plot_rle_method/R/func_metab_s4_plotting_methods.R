#' @title Plot RLE for MetaboliteAssayData
#' @name plotRle,MetaboliteAssayData-method
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames
#' @export
setMethod(
    f = "plotRle",
    signature = "MetaboliteAssayData",
    definition = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
        # --- Input Validation ---
        if (!is.character(grouping_variable) || length(grouping_variable) != 1 || is.na(grouping_variable)) {
            stop("`grouping_variable` must be a single non-NA character string.")
        }
        if (!grouping_variable %in% colnames(theObject@design_matrix)) {
            stop(sprintf("`grouping_variable` '%s' not found in design_matrix.", grouping_variable))
        }
        if (!is.null(sample_label) && (!is.character(sample_label) || length(sample_label) != 1 || sample_label == "")) {
            stop("`sample_label` must be NULL or a single non-empty character string.")
        }
        if (!is.null(sample_label) && !sample_label %in% colnames(theObject@design_matrix)) {
            stop(sprintf("`sample_label` '%s' not found in design_matrix.", sample_label))
        }

        design_matrix_df <- as.data.frame(theObject@design_matrix) # Ensure it's a data.frame for rownames
        sample_id_col_name <- theObject@sample_id
        metabolite_id_col_name <- theObject@metabolite_id_column
        assay_list <- theObject@metabolite_data

        if (!is.list(assay_list)) {
            assay_list <- list(assay_list)
        }

        if (length(assay_list) == 0) {
            warning("No assays found in `metabolite_data` slot. Returning empty list.")
            return(list())
        }

        # Ensure list is named
        if (is.null(names(assay_list))) {
            names(assay_list) <- paste0("Assay_", seq_along(assay_list))
            warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
        }

        message("--- Entering plotRle for MetaboliteAssayData ---")
        message(sprintf("   plotRle: grouping_variable = %s", grouping_variable))
        message(sprintf("   plotRle: yaxis_limit = %s", paste(yaxis_limit, collapse = ", ")))
        message(sprintf("   plotRle: sample_label = %s", ifelse(is.null(sample_label), "NULL", sample_label)))


        # --- Plotting Logic per Assay ---
        rle_plots_list <- purrr::map(seq_along(assay_list), function(i) {
            assay_name <- names(assay_list)[i]
            current_assay_data <- assay_list[[i]]

            # --- Correctly identify sample columns based on design matrix ---
            design_samples <- as.character(design_matrix_df[[sample_id_col_name]]) # Get sample IDs from design matrix
            all_assay_cols <- colnames(current_assay_data)
            sample_cols <- intersect(all_assay_cols, design_samples) # Find which design samples are columns in the assay
            metadata_cols <- setdiff(all_assay_cols, sample_cols) # All other columns are metadata/ID

            # Ensure the primary metabolite ID column is considered metadata
            if (metabolite_id_col_name %in% sample_cols) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' is also listed as a sample ID. Check configuration.", assay_name, metabolite_id_col_name))
            }
            metadata_cols <- union(metadata_cols, metabolite_id_col_name) # Ensure metabolite ID is not treated as a sample column
            sample_cols <- setdiff(all_assay_cols, metadata_cols) # Final list of sample columns


            if (length(sample_cols) == 0) {
                warning(sprintf("Assay '%s': No sample columns found in assay matching sample IDs in '%s' column of design matrix. Skipping RLE.", assay_name, sample_id_col_name))
                return(NULL)
            }
            # --- End Correction ---


            # Check sample consistency (now based on correctly identified sample_cols)
            design_samples_check <- design_matrix_df[[sample_id_col_name]] # Use original type for check
            missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
            if (length(missing_samples_in_design) > 0) {
                # This condition should theoretically not be met if sample_cols were derived from design_samples,
                # but keeping as a safeguard against type issues or unexpected data.
                warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping RLE.", assay_name, paste(missing_samples_in_design, collapse = ", ")))
                return(NULL)
            }


            # Filter design matrix to match assay samples
            design_matrix_filtered <- design_matrix_df[design_matrix_df[[sample_id_col_name]] %in% sample_cols, ]

            # Ensure metabolite ID column exists
            if (!metabolite_id_col_name %in% colnames(current_assay_data)) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping RLE.", assay_name, metabolite_id_col_name))
                return(NULL)
            }

            # Convert to matrix, ensuring correct sample columns
            frozen_metabolite_matrix <- current_assay_data |>
                tibble::column_to_rownames(metabolite_id_col_name) |>
                dplyr::select(all_of(sample_cols)) |> # Select only relevant sample columns
                as.matrix()

            # Check for sufficient data
            if (nrow(frozen_metabolite_matrix) < 1 || ncol(frozen_metabolite_matrix) < 1) {
                warning(sprintf("Assay '%s': Matrix has zero rows or columns after preparation. Skipping RLE.", assay_name))
                return(NULL)
            }

            # Handle sample labels
            rownames(design_matrix_filtered) <- design_matrix_filtered[[sample_id_col_name]] # Set rownames for indexing
            temp_matrix_colnames <- colnames(frozen_metabolite_matrix) # Store original colnames

            if (!is.null(sample_label)) {
                if (sample_label %in% colnames(design_matrix_filtered)) {
                    # Create a mapping from original sample ID to new label
                    label_map <- setNames(design_matrix_filtered[[sample_label]], design_matrix_filtered[[sample_id_col_name]])
                    # Apply the mapping to the matrix column names
                    colnames(frozen_metabolite_matrix) <- label_map[temp_matrix_colnames]
                    # Update the rownames of the filtered design matrix to match the new labels for lookup
                    rownames(design_matrix_filtered) <- design_matrix_filtered[[sample_label]]
                } # else: sample_label not found, already checked, but defensive
            }


            # Prepare rowinfo vector using the potentially updated colnames/rownames
            rowinfo_vector <- NA
            current_colnames <- colnames(frozen_metabolite_matrix)
            if (all(current_colnames %in% rownames(design_matrix_filtered))) {
                rowinfo_vector <- design_matrix_filtered[current_colnames, grouping_variable]
            } else {
                warning(sprintf("Assay '%s': Not all matrix column names ('%s') found in design matrix rownames after label application. Check sample_label consistency. Proceeding without fill.", assay_name, paste(head(current_colnames), collapse = ", ")))
                # Proceeding without rowinfo fill color if lookup fails
            }


            # Call the helper function
            tryCatch(
                {
                    plotRleHelper(
                        Y = t(frozen_metabolite_matrix), # Helper expects samples in rows
                        rowinfo = rowinfo_vector,
                        yaxis_limit = yaxis_limit
                    )
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error during RLE plotting: %s. Skipping.", assay_name, e$message))
                    return(NULL) # Skip on error
                }
            )
        })

        # Set names for the list of plots
        names(rle_plots_list) <- names(assay_list)

        # Remove NULL elements (skipped assays)
        rle_plots_list <- rle_plots_list[!sapply(rle_plots_list, is.null)]

        return(rle_plots_list)
    }
)

