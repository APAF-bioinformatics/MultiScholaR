#' @title Plot PCA for MetaboliteAssayData
#' @name plotPca,MetaboliteAssayData-method
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames
#' @export
setMethod(
    f = "plotPca",
    signature = "MetaboliteAssayData",
    definition = function(theObject, grouping_variable, shape_variable = NULL, label_column = NULL, title = NULL, font_size = 8) {
        # --- Input Validation ---
        if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
            stop("`grouping_variable` must be a single character string.")
        }
        if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
            stop("`shape_variable` must be NULL or a single character string.")
        }
        if (!grouping_variable %in% colnames(theObject@design_matrix)) {
            stop(sprintf("`grouping_variable` '%s' not found in design_matrix.", grouping_variable))
        }
        if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
            stop(sprintf("`shape_variable` '%s' not found in design_matrix.", shape_variable))
        }
        if (!is.null(label_column) && label_column != "" && !label_column %in% colnames(theObject@design_matrix)) {
            # Allow label_column to be empty/NULL, but if specified, it must exist
            stop(sprintf("`label_column` '%s' not found in design_matrix.", label_column))
        }

        design_matrix <- theObject@design_matrix
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

        # Ensure list is named, provide default names if not
        if (is.null(names(assay_list))) {
            names(assay_list) <- paste0("Assay_", seq_along(assay_list))
            warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
        }


        # --- Plotting Logic per Assay ---
        pca_plots_list <- purrr::map(seq_along(assay_list), function(i) {
            assay_name <- names(assay_list)[i]
            current_assay_data <- assay_list[[i]]

            # --- Correctly identify sample columns based on design matrix ---
            design_samples <- as.character(design_matrix[[sample_id_col_name]]) # Get sample IDs from design matrix
            all_assay_cols <- colnames(current_assay_data)
            sample_cols <- intersect(all_assay_cols, design_samples) # Find which design samples are columns in the assay
            metadata_cols <- setdiff(all_assay_cols, sample_cols) # All other columns are metadata/ID

            # Ensure the primary metabolite ID column is considered metadata if it's not a sample ID itself
            if (metabolite_id_col_name %in% sample_cols) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' is also listed as a sample ID. Check configuration.", assay_name, metabolite_id_col_name))
            }
            metadata_cols <- union(metadata_cols, metabolite_id_col_name) # Ensure metabolite ID is not treated as a sample column
            sample_cols <- setdiff(all_assay_cols, metadata_cols) # Final list of sample columns

            if (length(sample_cols) == 0) {
                warning(sprintf("Assay '%s': No sample columns found in assay matching sample IDs in '%s' column of design matrix. Skipping PCA.", assay_name, sample_id_col_name))
                return(NULL) # Skip this assay
            }
            # --- End Correction ---


            # Check if all identified sample columns exist in the design matrix (redundant check now, but safe)
            design_samples_check <- design_matrix[[sample_id_col_name]] # Use original type for check
            missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
            if (length(missing_samples_in_design) > 0) {
                warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping PCA.", assay_name, paste(missing_samples_in_design, collapse = ", ")))
                return(NULL)
            }

            # Filter design matrix to match assay samples
            design_matrix_filtered <- design_matrix[design_matrix[[sample_id_col_name]] %in% sample_cols, ]

            # Ensure metabolite ID column exists
            if (!metabolite_id_col_name %in% colnames(current_assay_data)) {
                warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping PCA.", assay_name, metabolite_id_col_name))
                return(NULL)
            }

            # Check for sufficient features after removing non-finite values
            frozen_metabolite_matrix_pca <- current_assay_data |>
                tibble::column_to_rownames(metabolite_id_col_name) |>
                dplyr::select(all_of(sample_cols)) |> # Ensure correct columns
                as.matrix()

            # Replace Inf/-Inf with NA
            frozen_metabolite_matrix_pca[!is.finite(frozen_metabolite_matrix_pca)] <- NA

            # Check for sufficient features and samples after NA handling
            valid_rows <- rowSums(is.finite(frozen_metabolite_matrix_pca)) > 1 # Need at least 2 points per feature for variance
            valid_cols <- colSums(is.finite(frozen_metabolite_matrix_pca)) > 1 # Need at least 2 points per sample for variance

            if (sum(valid_rows) < 2 || sum(valid_cols) < 2) {
                warning(sprintf("Assay '%s': Insufficient finite data points (< 2 features or < 2 samples with data) for PCA. Skipping.", assay_name))
                return(NULL)
            }

            frozen_metabolite_matrix_pca_final <- frozen_metabolite_matrix_pca[valid_rows, valid_cols, drop = FALSE]
            design_matrix_filtered_final <- design_matrix_filtered[design_matrix_filtered[[sample_id_col_name]] %in% colnames(frozen_metabolite_matrix_pca_final), ]

            # Generate title for this specific assay
            assay_title <- if (!is.null(title) && title != "") paste(title, "-", assay_name) else ""

            # --- Ensure consistent type for sample ID column before join ---
            design_matrix_filtered_final[[sample_id_col_name]] <- as.character(design_matrix_filtered_final[[sample_id_col_name]])
            # --- End type consistency fix ---

            # Call the helper function
            tryCatch(
                {
                    plotPcaHelper(
                        data = frozen_metabolite_matrix_pca_final,
                        design_matrix = design_matrix_filtered_final,
                        sample_id_column = sample_id_col_name,
                        grouping_variable = grouping_variable,
                        shape_variable = shape_variable,
                        label_column = label_column,
                        title = assay_title,
                        geom.text.size = font_size
                    )
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error during PCA plotting: %s. Skipping.", assay_name, e$message))
                    return(NULL) # Skip on error
                }
            )
        })

        # Set names for the list of plots
        names(pca_plots_list) <- names(assay_list)

        # Remove NULL elements (skipped assays)
        pca_plots_list <- pca_plots_list[!sapply(pca_plots_list, is.null)]

        return(pca_plots_list)
    }
)

