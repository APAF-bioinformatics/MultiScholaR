#' Metabolite Assay Data S4 Class
#'
#' An S4 class to store and manage multiple metabolomics quantitative datasets
#' derived from different platforms or ionization modes, along with associated
#' experimental design and metadata.
#'
#' @slot metabolite_data A named list of data frames. Each data frame contains
#'   quantitative data for one assay (e.g., LCMS-Pos). Rows should represent
#'   metabolites/features, and columns should represent samples. Each data frame
#'   MUST contain the column specified by `metabolite_id_column`. The names of
#'   the list elements should describe the assay source (e.g., "LCMS_Pos").
#' @slot metabolite_id_column Character string. The name of the column within
#'   each assay data frame that contains the **primary feature identifier** (used
#'   for uniquely identifying rows in quantitative analysis).
#' @slot annotation_id_column Character string. The name of the column (often
#'   in feature metadata, but could be in assay tables) that contains the
#'   **biological or database annotation identifier** (e.g., HMDB ID, KEGG ID,
#'   chemical name, CHEBI ID). This is the ID typically used for downstream
#'   biological mapping.
#' @slot database_identifier_type Character string. Describes the type or format
#'   of identifiers found in the column specified by `annotation_id_column`
#'   (e.g., "HMDB", "KEGG", "CHEBI", "Mixed_CHEBI_Unknown", "InternalName").
#' @slot internal_standard_regex Character string. A regular expression used to
#'   identify features that are internal standards based on their identifier in
#'   the `metabolite_id_column` or `annotation_id_column`. Set to `NA_character_` or `""`
#'   if not applicable or no internal standards are used.
#' @slot design_matrix A data frame containing the experimental design. Must
#'   include the column specified by `sample_id`.
#' @slot sample_id Character string. The name of the column in `design_matrix`
#'   that contains unique sample identifiers. These identifiers must correspond
#'   to the sample column names in the assay data frames.
#' @slot group_id Character string. The name of the column in `design_matrix`
#'   that defines the primary experimental groups or conditions.
#' @slot technical_replicate_id Character string. The name of the column in
#'   `design_matrix` that identifies technical replicates, if applicable. Use
#'   `NA_character_` if there are no technical replicates.
#' @slot args A list, typically populated from a configuration file, holding
#'   parameters used during processing.
#'
#' @importFrom methods setClass slotNames slot new
#' @importFrom dplyr pull distinct
#' @importFrom rlang sym !!
#' @exportClass MetaboliteAssayData
MetaboliteAssayData <- setClass("MetaboliteAssayData",
    slots = c(
        metabolite_data = "list",
        metabolite_id_column = "character",
        annotation_id_column = "character",
        database_identifier_type = "character",
        internal_standard_regex = "character",
        design_matrix = "data.frame",
        sample_id = "character",
        group_id = "character",
        technical_replicate_id = "character",
        args = "list"
    ),
    prototype = list(
        metabolite_data = list(),
        metabolite_id_column = "database_identifier",
        annotation_id_column = "metabolite_identification",
        database_identifier_type = "Unknown",
        internal_standard_regex = NA_character_,
        design_matrix = data.frame(),
        sample_id = "Sample_ID",
        group_id = "group",
        technical_replicate_id = NA_character_,
        args = list()
    ),
    validity = function(object) {
        errors <- character()

        # --- Basic Type Checks ---
        if (!is.list(object@metabolite_data)) {
            errors <- c(errors, "`metabolite_data` must be a list.")
        }
        if (length(object@metabolite_data) > 0 && !all(sapply(object@metabolite_data, is.data.frame))) {
            errors <- c(errors, "All elements in `metabolite_data` must be data frames.")
        }
        if (!is.character(object@metabolite_id_column) || length(object@metabolite_id_column) != 1) {
            errors <- c(errors, "`metabolite_id_column` must be a single character string.")
        }
        if (!is.character(object@annotation_id_column) || length(object@annotation_id_column) != 1) {
            errors <- c(errors, "`annotation_id_column` must be a single character string.")
        }
        if (!is.character(object@database_identifier_type) || length(object@database_identifier_type) != 1) {
            errors <- c(errors, "`database_identifier_type` must be a single character string.")
        }
        if (!is.character(object@internal_standard_regex) || length(object@internal_standard_regex) != 1) {
            errors <- c(errors, "`internal_standard_regex` must be a single character string (can be NA_character_).")
        }
        if (!is.data.frame(object@design_matrix)) {
            errors <- c(errors, "`design_matrix` must be a data frame.")
        }
        if (!is.character(object@sample_id) || length(object@sample_id) != 1) {
            errors <- c(errors, "`sample_id` must be a single character string.")
        }
        if (!is.character(object@group_id) || length(object@group_id) != 1) {
            errors <- c(errors, "`group_id` must be a single character string.")
        }
        if (!is.character(object@technical_replicate_id) || length(object@technical_replicate_id) != 1) {
            errors <- c(errors, "`technical_replicate_id` must be a single character string (can be NA_character_).")
        }
        if (!is.list(object@args)) {
            errors <- c(errors, "`args` must be a list.")
        }

        # --- Content Checks ---
        if (length(object@metabolite_data) > 0) {
            # Check primary metabolite ID column exists in all assays
            metabolite_col_exists <- sapply(object@metabolite_data, function(df) {
                object@metabolite_id_column %in% colnames(df)
            })
            if (!all(metabolite_col_exists)) {
                missing_in <- names(object@metabolite_data)[!metabolite_col_exists]
                errors <- c(errors, paste0("Primary `metabolite_id_column` ('", object@metabolite_id_column, "') not found in assay(s): ", paste(missing_in, collapse = ", ")))
            }

            # Check sample ID column exists in design matrix
            if (!(object@sample_id %in% colnames(object@design_matrix))) {
                errors <- c(errors, paste0("`sample_id` column ('", object@sample_id, "') not found in `design_matrix`."))
            } else {
                # Check consistency of samples between assays and design matrix
                samples_in_design <- tryCatch(
                    object@design_matrix |> dplyr::pull(!!rlang::sym(object@sample_id)) |> as.character() |> unique() |> sort(),
                    error = function(e) character(0) # Handle case where column doesn't exist
                )

                if (length(samples_in_design) == 0 && !(object@sample_id %in% colnames(object@design_matrix))) {
                    # Error already captured above
                } else {
                    first_assay_samples <- NULL
                    for (assay_name in names(object@metabolite_data)) {
                        assay_df <- object@metabolite_data[[assay_name]]
                        assay_samples <- setdiff(colnames(assay_df), object@metabolite_id_column) |> sort()

                        if (is.null(first_assay_samples)) {
                            first_assay_samples <- assay_samples
                        }

                        # Check all assays have the same sample columns
                        if (!identical(assay_samples, first_assay_samples)) {
                            errors <- c(errors, paste0("Sample columns differ between assays. Check assay '", assay_name, "'."))
                            # Stop checking further assays for sample consistency after first mismatch
                            break
                        }
                        # Check assay samples match design matrix samples (only need to check once if all assays are identical)
                        if (assay_name == names(object@metabolite_data)[1]) {
                            if (!identical(assay_samples, samples_in_design)) {
                                errors <- c(errors, paste0("Sample columns in assays do not exactly match unique sample IDs ('", object@sample_id, "') in `design_matrix`."))
                                # Stop checking further assays for sample consistency after first mismatch
                                break
                            }
                        }
                    }
                }
            }
        }

        # --- Final Check ---
        if (length(errors) == 0) TRUE else errors
    }
)
#' @export MetaboliteAssayData


#' Create MetaboliteAssayData Object
#'
#' Constructor function for the MetaboliteAssayData class.
#'
#' @param metabolite_data Named list of data frames (assays).
#' @param design_matrix Experimental design data frame.
#' @param metabolite_id_column Name of the **primary feature ID** column within assays (e.g., `"database_identifier"`).
#' @param annotation_id_column Name of the **annotation ID** column (e.g., `"metabolite_identification"`).
#' @param sample_id Name of the sample ID column in design_matrix and assays.
#' @param group_id Name of the group column in design_matrix.
#' @param technical_replicate_id Name of the technical replicate column in design_matrix (use NA_character_ if none).
#' @param database_identifier_type Type of identifier in the `annotation_id_column` (e.g., `"Mixed_CHEBI_Unknown"`).
#' @param internal_standard_regex Regex to identify internal standards. Use `NA_character_` or `""` if none.
#' @param args List of arguments (e.g., from config).
#'
#' @return A MetaboliteAssayData object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming lcms_pos_df, lcms_neg_df, gcms_df are data frames
#' # with 'Metabolite' as ID column and samples as other columns
#' # Assuming design_df has 'SampleID', 'Group', 'Replicate' columns
#' assays_list <- list(
#'     LCMS_Pos = lcms_pos_df,
#'     LCMS_Neg = lcms_neg_df,
#'     GCMS = gcms_df
#' )
#' config <- list(...) # Your config list
#'
#' met_assay_obj <- createMetaboliteAssayData(
#'     metabolite_data = assays_list,
#'     design_matrix = design_df,
#'     metabolite_id_column = "Metabolite",
#'     sample_id = "SampleID",
#'     group_id = "Group",
#'     technical_replicate_id = "Replicate",
#'     database_identifier_type = "InternalName",
#'     internal_standard_regex = "^IS_",
#'     args = config
#' )
#' }
createMetaboliteAssayData <- function(
    metabolite_data,
    design_matrix,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "metabolite_identification",
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    database_identifier_type = "Unknown",
    internal_standard_regex = NA_character_,
    args = list()) {
    # Perform basic checks before creating the object
    stopifnot(is.list(metabolite_data))
    stopifnot(all(sapply(metabolite_data, is.data.frame)))
    stopifnot(is.data.frame(design_matrix))
    # Add more checks as needed...

    obj <- new("MetaboliteAssayData",
        metabolite_data = metabolite_data,
        metabolite_id_column = metabolite_id_column,
        annotation_id_column = annotation_id_column,
        database_identifier_type = database_identifier_type,
        internal_standard_regex = internal_standard_regex,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
    )
    # Validity check is automatically called by 'new'
    return(obj)
}

##-----------------------------------------------------------------------------
## Plotting Methods for MetaboliteAssayData
##-----------------------------------------------------------------------------

#' @describeIn plotPca Method for MetaboliteAssayData
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames
#' @export
setMethod(f = "plotPca",
          signature = "MetaboliteAssayData",
          definition = function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size = 8) {
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
                if(length(missing_samples_in_design) > 0) {
                     warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping PCA.", assay_name, paste(missing_samples_in_design, collapse=", ")))
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
                assay_title <- paste(title, "-", assay_name)

                # --- Ensure consistent type for sample ID column before join ---
                design_matrix_filtered_final[[sample_id_col_name]] <- as.character(design_matrix_filtered_final[[sample_id_col_name]])
                # --- End type consistency fix ---

                # Call the helper function
                tryCatch({
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
                }, error = function(e) {
                    warning(sprintf("Assay '%s': Error during PCA plotting: %s. Skipping.", assay_name, e$message))
                    return(NULL) # Skip on error
                })
            })

            # Set names for the list of plots
            names(pca_plots_list) <- names(assay_list)

            # Remove NULL elements (skipped assays)
            pca_plots_list <- pca_plots_list[!sapply(pca_plots_list, is.null)]

            return(pca_plots_list)
          })

#' @describeIn plotRle Method for MetaboliteAssayData
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames
#' @export
setMethod(f = "plotRle",
          signature = "MetaboliteAssayData",
          definition = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
            # --- Input Validation ---
            if (!is.character(grouping_variable) || length(grouping_variable) != 1 || is.na(grouping_variable)) {
              stop("`grouping_variable` must be a single non-NA character string.")
            }
            if (!grouping_variable %in% colnames(theObject@design_matrix)) {
              stop(sprintf("`grouping_variable` '%s' not found in design_matrix.", grouping_variable))
            }
             if (!is.null(sample_label) && (!is.character(sample_label) || length(sample_label) != 1 || sample_label == "" )) {
                stop("`sample_label` must be NULL or a single non-empty character string.")
             }
             if (!is.null(sample_label) && !sample_label %in% colnames(theObject@design_matrix)) {
                stop(sprintf("`sample_label` '%s' not found in design_matrix.", sample_label))
             }

            design_matrix_df <- as.data.frame(theObject@design_matrix) # Ensure it's a data.frame for rownames
            sample_id_col_name <- theObject@sample_id
            metabolite_id_col_name <- theObject@metabolite_id_column
            assay_list <- theObject@metabolite_data

            if (length(assay_list) == 0) {
                warning("No assays found in `metabolite_data` slot. Returning empty list.")
                return(list())
            }

            # Ensure list is named
            if (is.null(names(assay_list))) {
                names(assay_list) <- paste0("Assay_", seq_along(assay_list))
                warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
            }

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
                if(length(missing_samples_in_design) > 0) {
                     # This condition should theoretically not be met if sample_cols were derived from design_samples,
                     # but keeping as a safeguard against type issues or unexpected data.
                     warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping RLE.", assay_name, paste(missing_samples_in_design, collapse=", ")))
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

                 if(!is.null(sample_label)) {
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
                   warning(sprintf("Assay '%s': Not all matrix column names ('%s') found in design matrix rownames after label application. Check sample_label consistency. Proceeding without fill.", assay_name, paste(head(current_colnames), collapse=", ")))
                   # Proceeding without rowinfo fill color if lookup fails
                }


                # Call the helper function
                tryCatch({
                    plotRleHelper(
                        Y = t(frozen_metabolite_matrix), # Helper expects samples in rows
                        rowinfo = rowinfo_vector,
                        yaxis_limit = yaxis_limit
                    )
                }, error = function(e) {
                    warning(sprintf("Assay '%s': Error during RLE plotting: %s. Skipping.", assay_name, e$message))
                    return(NULL) # Skip on error
                })
            })

            # Set names for the list of plots
            names(rle_plots_list) <- names(assay_list)

            # Remove NULL elements (skipped assays)
            rle_plots_list <- rle_plots_list[!sapply(rle_plots_list, is.null)]

            return(rle_plots_list)
          })


#' @describeIn plotDensity Method for MetaboliteAssayData
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw labs theme element_blank element_text margin
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom mixOmics pca
#' @importFrom dplyr left_join select all_of
#' @importFrom rlang sym !!
#' @export
setMethod(f = "plotDensity",
          signature = "MetaboliteAssayData",
          definition = function(theObject, grouping_variable, title = "", font_size = 8) {
              # --- Input is MetaboliteAssayData: Original logic (Calculate PCA) ---

              # --- Input Validation ---
              if (!is.character(grouping_variable) || length(grouping_variable) != 1 || is.na(grouping_variable)) {
                stop("`grouping_variable` must be a single non-NA character string.")
              }
              # Ensure design matrix exists and has the grouping variable
              if (is.null(slot(theObject, "design_matrix")) || nrow(slot(theObject, "design_matrix")) == 0) {
                  stop("`design_matrix` slot is missing or empty in the MetaboliteAssayData object.")
              }
              if (!grouping_variable %in% colnames(theObject@design_matrix)) {
                stop(sprintf("`grouping_variable` '%s' not found in design_matrix.", grouping_variable))
              }

              design_matrix <- theObject@design_matrix
              sample_id_col_name <- theObject@sample_id
              metabolite_id_col_name <- theObject@metabolite_id_column
              assay_list <- theObject@metabolite_data

               if (length(assay_list) == 0) {
                  warning("No assays found in `metabolite_data` slot. Returning empty list.")
                  return(list())
              }

              # Ensure list is named
              if (is.null(names(assay_list))) {
                  names(assay_list) <- paste0("Assay_", seq_along(assay_list))
                  warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
              }

              # --- Plotting Logic per Assay ---
              density_plots_list <- purrr::map(seq_along(assay_list), function(i) {
                   assay_name <- names(assay_list)[i]
                   current_assay_data <- assay_list[[i]]

                   # --- Correctly identify sample columns based on design matrix ---
                   design_samples <- as.character(design_matrix[[sample_id_col_name]]) # Get sample IDs from design matrix
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
                       warning(sprintf("Assay '%s': No sample columns found in assay matching sample IDs in '%s' column of design matrix. Skipping Density plot.", assay_name, sample_id_col_name))
                       return(NULL)
                   }
                   # --- End Correction ---


                   # Check sample consistency
                   design_samples_check <- design_matrix[[sample_id_col_name]] # Use original type for check
                   missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
                   if(length(missing_samples_in_design) > 0) {
                        warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping Density plot.", assay_name, paste(missing_samples_in_design, collapse=", ")))
                        return(NULL)
                   }

                   # Filter design matrix
                   design_matrix_filtered <- design_matrix[as.character(design_matrix[[sample_id_col_name]]) %in% sample_cols, ]

                   # Ensure metabolite ID column exists
                   if (!metabolite_id_col_name %in% colnames(current_assay_data)) {
                       warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping Density plot.", assay_name, metabolite_id_col_name))
                       return(NULL)
                   }

                   # Convert to matrix, handle non-finite, check dimensions (similar to plotPca)
                   frozen_metabolite_matrix_pca <- current_assay_data |>
                       tibble::column_to_rownames(metabolite_id_col_name) |>
                       dplyr::select(all_of(sample_cols)) |>
                       as.matrix()
                   frozen_metabolite_matrix_pca[!is.finite(frozen_metabolite_matrix_pca)] <- NA

                   valid_rows <- rowSums(is.finite(frozen_metabolite_matrix_pca)) > 1
                   valid_cols <- colSums(is.finite(frozen_metabolite_matrix_pca)) > 1
                   if (sum(valid_rows) < 2 || sum(valid_cols) < 2) {
                       warning(sprintf("Assay '%s': Insufficient finite data points for PCA. Skipping Density plot.", assay_name))
                       return(NULL)
                   }
                   frozen_metabolite_matrix_pca_final <- frozen_metabolite_matrix_pca[valid_rows, valid_cols, drop = FALSE]
                   design_matrix_filtered_final <- design_matrix_filtered[as.character(design_matrix_filtered[[sample_id_col_name]]) %in% colnames(frozen_metabolite_matrix_pca_final), ]
                   # Ensure consistent type for sample ID column before join
                    design_matrix_filtered_final[[sample_id_col_name]] <- as.character(design_matrix_filtered_final[[sample_id_col_name]])


                   # --- Perform PCA ---
                   pca_data_for_plot <- tryCatch({
                       pca.res <- mixOmics::pca(t(as.matrix(frozen_metabolite_matrix_pca_final)), ncomp = 2) # Ensure ncomp=2 for PC1/PC2
                       pca.res$variates$X |>
                           as.data.frame() |>
                           tibble::rownames_to_column(var = sample_id_col_name) |>
                           dplyr::left_join(design_matrix_filtered_final, by = sample_id_col_name) |>
                           tibble::as_tibble() # Ensure it's a tibble
                   }, error = function(e) {
                       warning(sprintf("Assay '%s': Error during PCA calculation for density plot: %s. Skipping.", assay_name, e$message))
                       return(NULL) # Skip this assay if PCA fails
                   })

                   if (is.null(pca_data_for_plot) || !("PC1" %in% colnames(pca_data_for_plot)) || !("PC2" %in% colnames(pca_data_for_plot))) {
                      warning(sprintf("Assay '%s': PCA result is invalid or missing PC1/PC2. Skipping Density plot.", assay_name))
                      return(NULL)
                   }
                   if (!grouping_variable %in% colnames(pca_data_for_plot)) {
                        warning(sprintf("Assay '%s': Grouping variable '%s' not found in PCA data after join. Skipping Density plot.", assay_name, grouping_variable))
                        return(NULL)
                   }


                   # --- Create Density/Box Plots ---
                   assay_title <- paste(title, "-", assay_name)

                   tryCatch({
                        # Create PC1 boxplot
                       pc1_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC1, fill = !!rlang::sym(grouping_variable))) +
                           geom_boxplot(notch = FALSE) +
                           theme_bw() +
                           labs(title = assay_title,
                                x = "",
                                y = "PC1") +
                           theme(
                               legend.position = "none",
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               text = element_text(size = font_size),
                               plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank()
                           )

                       # Create PC2 boxplot
                       pc2_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC2, fill = !!rlang::sym(grouping_variable))) +
                           geom_boxplot(notch = FALSE) +
                           theme_bw() +
                           labs(x = "",
                                y = "PC2") +
                           theme(
                               legend.position = "none",
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               text = element_text(size = font_size),
                               plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank()
                           )

                       # Combine plots with minimal spacing
                       combined_plot <- pc1_box / pc2_box +
                           patchwork::plot_layout(heights = c(1, 1)) +
                           patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

                       return(combined_plot)

                   }, error = function(e) {
                       warning(sprintf("Assay '%s': Error creating density boxplots: %s. Skipping.", assay_name, e$message))
                       return(NULL)
                   })
              })

              # Set names for the list of plots
              names(density_plots_list) <- names(assay_list)

              # Remove NULL elements (skipped assays)
              density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

              return(density_plots_list)
          })


#' @describeIn plotDensity Method for list of ggplot objects
#' @importFrom purrr map set_names
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw labs theme element_blank element_text margin
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom rlang sym !!
#' @export
setMethod(f = "plotDensity",
          signature = c(theObject = "list"), # Explicitly define signature argument
          definition = function(theObject, grouping_variable, title = "", font_size = 8) {
               # --- Input is a list: Assume list of ggplot objects (Use existing PCA data) ---

               # Basic validation
               if (!all(sapply(theObject, function(x) inherits(x, "ggplot")))) {
                  stop("If 'theObject' is a list, all its elements must be ggplot objects.")
               }
               if (!is.character(grouping_variable) || length(grouping_variable) != 1 || is.na(grouping_variable)) {
                 stop("`grouping_variable` must be a single non-NA character string.")
               }

               pca_plots_list <- theObject # Rename for clarity

               if (length(pca_plots_list) == 0) {
                   warning("Input list of ggplot objects is empty. Returning empty list.")
                   return(list())
               }

               # Ensure list is named, provide default names if not
               if (is.null(names(pca_plots_list))) {
                   names(pca_plots_list) <- paste0("Plot_", seq_along(pca_plots_list))
                   warning("Input ggplot list was unnamed. Using default names (Plot_1, Plot_2, ...).")
               }

               # --- Plotting Logic per Input ggplot ---
               density_plots_list <- purrr::map(seq_along(pca_plots_list), function(i) {
                   pca_plot <- pca_plots_list[[i]]
                   plot_name <- names(pca_plots_list)[i]
                   # Use title override if provided, otherwise use original plot title or name
                   plot_title_final <- if (!is.null(title) && title != "") paste(title, "-", plot_name) else tryCatch(pca_plot$labels$title, error = function(e) plot_name)

                   # --- Extract PCA data from the ggplot object ---
                   pca_data_for_plot <- NULL
                   if (!is.null(pca_plot$data) && is.data.frame(pca_plot$data) && all(c("PC1", "PC2", grouping_variable) %in% colnames(pca_plot$data))) {
                       pca_data_for_plot <- tibble::as_tibble(pca_plot$data)
                   } else {
                       warning(sprintf("Plot '%s': Could not reliably extract required data (PC1, PC2, %s) from the ggplot object's internal structure. Skipping density plot generation.", plot_name, grouping_variable))
                       return(NULL)
                   }

                   if (!grouping_variable %in% colnames(pca_data_for_plot)) {
                       warning(sprintf("Plot '%s': Grouping variable '%s' not found in extracted data. Skipping.", plot_name, grouping_variable))
                       return(NULL)
                   }

                   # --- Create Density/Box Plots (using extracted data) ---
                   tryCatch({
                      # Create PC1 boxplot
                     pc1_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC1, fill = !!rlang::sym(grouping_variable))) +
                         geom_boxplot(notch = FALSE) +
                         theme_bw() +
                         labs(title = plot_title_final, # Use final title
                              x = "",
                              y = "PC1") +
                         theme(
                             legend.position = "none",
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             text = element_text(size = font_size),
                             plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank()
                         )

                     # Create PC2 boxplot
                     pc2_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC2, fill = !!rlang::sym(grouping_variable))) +
                         geom_boxplot(notch = FALSE) +
                         theme_bw() +
                         labs(x = "",
                              y = "PC2") +
                         theme(
                             legend.position = "none",
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             text = element_text(size = font_size),
                             plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank()
                         )

                     # Combine plots with minimal spacing
                     combined_plot <- pc1_box / pc2_box +
                         patchwork::plot_layout(heights = c(1, 1)) +
                         patchwork::plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

                     return(combined_plot)

                   }, error = function(e) {
                       warning(sprintf("Plot '%s': Error creating density boxplots from ggplot input: %s. Skipping.", plot_name, e$message))
                       return(NULL)
                   })
               })

               # Set names for the list of plots
               names(density_plots_list) <- names(pca_plots_list)

               # Remove NULL elements (skipped plots)
               density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

               return(density_plots_list)
          })


#' @describeIn pearsonCorForSamplePairs Method for MetaboliteAssayData
#' @importFrom purrr map set_names map_df
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join select mutate filter distinct arrange group_by summarise ungroup pull if_else case_when row_number n rename add_column relocate
#' @importFrom stringr str_detect
#' @importFrom rlang sym !! :=
#' @export
setMethod(f = "pearsonCorForSamplePairs",
          signature = "MetaboliteAssayData",
          definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
            # --- Input Validation ---
            # tech_rep_remove_regex can be NULL, checked inside helper/later use
            # correlation_group can be NA, checked below

            design_matrix <- theObject@design_matrix
            sample_id_col_name <- theObject@sample_id
            metabolite_id_col_name <- theObject@metabolite_id_column
            tech_rep_col_name <- theObject@technical_replicate_id # Default grouping if correlation_group is NA
            assay_list <- theObject@metabolite_data

            if (length(assay_list) == 0) {
                warning("No assays found in `metabolite_data` slot. Returning empty list.")
                return(list())
            }

             # Ensure list is named
             if (is.null(names(assay_list))) {
                 names(assay_list) <- paste0("Assay_", seq_along(assay_list))
                 warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
             }

             # Determine the actual grouping column to use for pairing samples
             replicate_group_column_name <- correlation_group
             if (is.na(correlation_group)) {
                 if (is.na(tech_rep_col_name) || !tech_rep_col_name %in% colnames(design_matrix)) {
                     stop("`correlation_group` is NA and `technical_replicate_id` ('", tech_rep_col_name, "') is NA or not found in design_matrix. Cannot determine sample pairing.")
                 }
                 replicate_group_column_name <- tech_rep_col_name
             } else {
                 if (!correlation_group %in% colnames(design_matrix)) {
                     stop(sprintf("Specified `correlation_group` ('%s') not found in design_matrix.", correlation_group))
                 }
             }

            # Resolve tech_rep_remove_regex from config if needed (or use default)
            # Assuming the helper function or subsequent filtering handles NULL regex gracefully (meaning no filtering)
            tech_rep_remove_regex_final <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", tech_rep_remove_regex) # Allow override
            # theObject <- updateParamInObject(theObject, "tech_rep_remove_regex") # Update object if needed


            # --- Correlation Logic per Assay ---
            correlation_results_list <- purrr::map(seq_along(assay_list), function(i) {
                 assay_name <- names(assay_list)[i]
                 current_assay_data <- assay_list[[i]]

                 # --- Correctly identify sample columns based on design matrix ---
                 design_samples <- as.character(design_matrix[[sample_id_col_name]]) # Get sample IDs from design matrix
                 all_assay_cols <- colnames(current_assay_data)
                 sample_cols <- intersect(all_assay_cols, design_samples) # Find which design samples are columns in the assay
                 metadata_cols <- setdiff(all_assay_cols, sample_cols) # All other columns are metadata/ID

                 # Ensure the primary metabolite ID column is considered metadata
                 if (metabolite_id_col_name %in% sample_cols) {
                      warning(sprintf("Assay '%s': Metabolite ID column '%s' is also listed as a sample ID. Check configuration.", assay_name, metabolite_id_col_name))
                 }
                 metadata_cols <- union(metadata_cols, metabolite_id_col_name) # Ensure metabolite ID is not treated as a sample column
                 sample_cols <- setdiff(all_assay_cols, metadata_cols) # Final list of sample columns

                 if (length(sample_cols) < 2) { # Need at least 2 samples for correlation
                     warning(sprintf("Assay '%s': Fewer than 2 sample columns found matching design matrix sample IDs. Skipping Pearson correlation.", assay_name))
                     return(NULL)
                 }
                 # --- End Correction ---


                 # Check sample consistency (now based on correctly identified sample_cols)
                 design_samples_check <- design_matrix[[sample_id_col_name]] # Use original type for check
                 missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
                 if(length(missing_samples_in_design) > 0) {
                      # This condition should theoretically not be met if sample_cols were derived from design_samples,
                      # but keeping as a safeguard against type issues or unexpected data.
                      warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping Pearson correlation.", assay_name, paste(missing_samples_in_design, collapse=", ")))
                      return(NULL)
                 }


                 # Filter design matrix to match assay samples
                 design_matrix_filtered <- design_matrix[as.character(design_matrix[[sample_id_col_name]]) %in% sample_cols, ]

                 # Ensure metabolite ID column exists
                  if (!metabolite_id_col_name %in% colnames(current_assay_data)) {
                      warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping Pearson correlation.", assay_name, metabolite_id_col_name))
                      return(NULL)
                  }


                 # Prepare long data for helper
                 assay_long <- current_assay_data |>
                     tidyr::pivot_longer(
                         cols = all_of(sample_cols),
                         names_to = sample_id_col_name,
                         values_to = "abundance"
                     )


                 # Prepare the design matrix subset for the helper
                 design_subset <- design_matrix_filtered |>
                     dplyr::select(!!rlang::sym(sample_id_col_name), !!rlang::sym(replicate_group_column_name))


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
                  pairs_for_comparison <- tryCatch({
                      getPairsOfSamplesTable(design_subset, # Contains sample_id and replicate_group_column
                                            run_id_column = sample_id_col_name,
                                            replicate_group_column = replicate_group_column_name)
                  }, error = function(e) {
                       warning(sprintf("Assay '%s': Error getting sample pairs: %s. Skipping correlation.", assay_name, e$message))
                       return(NULL)
                  })

                  if(is.null(pairs_for_comparison) || nrow(pairs_for_comparison) == 0) {
                      warning(sprintf("Assay '%s': No valid sample pairs found for correlation. Skipping.", assay_name))
                      return(NULL)
                  }

                  # Get the names of the columns containing paired sample IDs (e.g., "Run.x", "Run.y")
                  run_id_col_x <- paste0(sample_id_col_name, ".x")
                  run_id_col_y <- paste0(sample_id_col_name, ".y")

                  # Check if these columns exist in the pairs table
                  if (!all(c(run_id_col_x, run_id_col_y) %in% colnames(pairs_for_comparison))) {
                      warning(sprintf("Assay '%s': Expected paired sample columns ('%s', '%s') not found in pairs table. Skipping correlation.", assay_name, run_id_col_x, run_id_col_y))
                      return(NULL)
                  }

                  # Calculate correlations as a separate vector first
                  calculated_correlations <- tryCatch({
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
                  }, error = function(e) {
                      warning(sprintf("Assay '%s': Error during map2_dbl correlation calculation: %s. Returning NULL results.", assay_name, e$message))
                      return(NULL) # Return NULL if map2_dbl fails
                  })

                  # Check if calculation succeeded and add the column
                  if (is.null(calculated_correlations)) {
                      correlation_results_raw <- NULL # Propagate failure
                  } else if (length(calculated_correlations) != nrow(pairs_for_comparison)) {
                       warning(sprintf("Assay '%s': Number of calculated correlations (%d) does not match number of pairs (%d). Skipping.", assay_name, length(calculated_correlations), nrow(pairs_for_comparison)))
                       correlation_results_raw <- NULL
                  } else {
                       correlation_results_raw <- pairs_for_comparison |>
                           dplyr::mutate(pearson_correlation = calculated_correlations)
                  }
                  # ---------------------------------------------------------- #

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
                        }
                   } else {
                       # If no regex, correlation_results_filtered remains unchanged (holds original data)
                   }

                   return(correlation_results_filtered)
             })

            # Set names for the list of results
            names(correlation_results_list) <- names(assay_list)

            # Remove NULL elements (skipped assays)
            correlation_results_list <- correlation_results_list[!sapply(correlation_results_list, is.null)]

            return(correlation_results_list)
           })


#' @describeIn plotPearson Method for MetaboliteAssayData
#' @importFrom purrr map set_names
#' @importFrom ggplot2 ggplot aes geom_histogram scale_y_continuous xlab ylab theme element_blank
#' @export
setMethod(f = "plotPearson",
          signature = "MetaboliteAssayData",
          definition = function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {

            # Get the list of correlation tibbles (one per assay)
            # tech_rep_remove_regex and correlation_group are passed down
            correlation_list <- pearsonCorForSamplePairs(theObject,
                                                        tech_rep_remove_regex = tech_rep_remove_regex,
                                                        correlation_group = correlation_group)

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
                tryCatch({
                    pearson_plot <- correlation_vec |>
                        ggplot(aes(pearson_correlation)) +
                        geom_histogram(breaks = hist_breaks, na.rm = TRUE) +
                        # Set x-axis limits and breaks
                        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
                        # Set fixed y-axis scale, similar to protein version
                        scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
                        xlab("Pearson Correlation") +
                        ylab("Counts") +
                        ggtitle(paste("Pearson Correlation -", assay_name)) +
                        theme_bw() + 
                        theme(
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()
                           )

                    return(pearson_plot)
                }, error = function(e) {
                     warning(sprintf("Assay '%s': Error creating Pearson histogram: %s. Skipping.", assay_name, e$message))
                     return(NULL)
                })
            })

            # Set names for the list of plots
            names(pearson_plots_list) <- names(correlation_list)

             # Remove NULL elements (skipped assays)
            pearson_plots_list <- pearson_plots_list[!sapply(pearson_plots_list, is.null)]

            return(pearson_plots_list)
          })


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
    wide_pair_table <- tryCatch({
        input_pair_table |>
            dplyr::select(!!rlang::sym(feature_id_column), !!rlang::sym(sample_id_column), !!rlang::sym(value_column)) |>
            tidyr::pivot_wider(
                names_from = !!rlang::sym(sample_id_column),
                values_from = !!rlang::sym(value_column)
            )
    }, error = function(e) {
        warning(sprintf("Error pivoting data wider for correlation between %s and %s: %s", sample_x_id, sample_y_id, e$message))
        return(NULL)
    })

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
    cor_result <- tryCatch({
        stats::cor(values_x, values_y, use = "pairwise.complete.obs")
    }, error = function(e) {
        warning(sprintf("calculateMetabolitePairCorrelation: Error in stats::cor for samples %s and %s: %s", sample_x_id, sample_y_id, e$message))
        return(NA_real_) # Returns NA_real_ on cor error
    })

    # --- Modified Check ---
    # Ensure the result is a single, finite numeric value
    if (length(cor_result) != 1 || !is.numeric(cor_result) || !is.finite(cor_result)) {
        return(NA_real_)
    }
    # --- End Modified Check ---

    return(cor_result)
}

# ------------------------------------------------------- #

#' @describeIn normaliseBetweenSamples Method for MetaboliteAssayData
#' @param theObject Object of class MetaboliteAssayData
#' @param normalisation_method Method to use for normalization. Options are
#'   "cyclicloess", "quantile", "scale", "none". If NULL, the value is retrieved
#'   from the object's configuration arguments (looking for "normalisation_method").
#'
#' @importFrom limma normalizeCyclicLoess normalizeQuantiles normalizeMedianAbsValues
#' @importFrom purrr map set_names
#' @importFrom tibble column_to_rownames rownames_to_column as_tibble
#' @importFrom dplyr select all_of left_join relocate any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom methods slot slot<- is
#' @importFrom logger log_info
#' @export
setMethod(f = "normaliseBetweenSamples",
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
                    if (is.null(assay_tibble)) return(assay_list[[i]]) # Return original if coercion failed
                }
                 if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                     warning(sprintf("Assay '%s': Metabolite ID column '%s' not found. Skipping normalization.", assay_name, metabolite_id_col_name), immediate. = TRUE)
                     return(assay_tibble)
                 }

                 # --- Identify Sample Columns ---
                 design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) { character(0) })
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
                    warning(sprintf("Assay '%s': Non-numeric sample columns found: %s. Attempting coercion, but this may indicate upstream issues.", assay_name, paste(non_numeric_samples, collapse=", ")), immediate. = TRUE)
                    assay_tibble <- assay_tibble |>
                        dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
                 }

                 # --- Prepare Matrix for Normalization ---
                 assay_matrix <- tryCatch({
                    assay_tibble |>
                       dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |> # Select ID + Samples
                       tibble::column_to_rownames(var = metabolite_id_col_name) |>
                       as.matrix()
                 }, error = function(e) {
                     warning(sprintf("Assay '%s': Error converting tibble to matrix: %s. Skipping normalization.", assay_name, e$message), immediate. = TRUE)
                     return(NULL)
                 })
                 if (is.null(assay_matrix)) return(assay_tibble) # Return original if matrix conversion failed

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
                    normalized_matrix <- tryCatch({
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
                    }, error = function(e) {
                        warning(sprintf("Assay '%s': Error during '%s' normalization: %s. Returning unnormalized data for this assay.", assay_name, normalisation_method_final, e$message), immediate. = TRUE)
                        return(assay_matrix) # Return original matrix on error
                    })
                } else {
                     message("   Normalization method is 'none'. Skipping application.")
                }

                normalized_matrix[!is.finite(normalized_matrix)] <- NA # Ensure NAs remain NAs

                # --- Reconstruct Tibble ---
                 # Get original metadata columns
                 metadata_cols <- setdiff(colnames(assay_tibble), sample_cols)

                 reconstructed_tibble <- tryCatch({
                     normalized_data_tibble <- normalized_matrix |>
                        as.data.frame() |> # Convert matrix to data frame
                        tibble::rownames_to_column(var = metabolite_id_col_name) |>
                        tibble::as_tibble()

                     original_metadata_tibble <- assay_tibble |>
                        dplyr::select(dplyr::any_of(metadata_cols)) # Use any_of in case some metadata cols were dynamic

                     # Join normalized data with original metadata
                     dplyr::left_join(original_metadata_tibble, normalized_data_tibble, by = metabolite_id_col_name) |>
                     # Ensure original column order (metadata first, then samples in original order)
                     dplyr::relocate(dplyr::all_of(metadata_cols), dplyr::all_of(sample_cols))

                 }, error = function(e) {
                      warning(sprintf("Assay '%s': Error reconstructing tibble after normalization: %s. Returning original data.", assay_name, e$message), immediate. = TRUE)
                      return(assay_tibble) # Return original on error
                 })

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
            theObject <- tryCatch({
                cleanDesignMatrix(theObject)
            }, error = function(e) {
                 warning(sprintf("Error running cleanDesignMatrix after normalization: %s. Design matrix might not be fully synchronized.", e$message))
                 return(theObject) # Return object even if cleaning fails
            })

            log_info("Between-sample normalization process finished for all assays.")
            return(theObject)
          })

#' @describeIn cleanDesignMatrix Method for MetaboliteAssayData
#' @importFrom dplyr inner_join select rename filter all_of any_of
#' @importFrom rlang sym !!
#' @importFrom methods slot
#' @importFrom tibble tibble
#' @export
setMethod(f = "cleanDesignMatrix",
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
            design_samples <- tryCatch(as.character(design_matrix[[sample_id_col_name]]), error = function(e) { character(0) })
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

            cleaned_design_matrix <- tryCatch({
                 # Create a tibble with just the sample IDs in the order they appear in the data
                 sample_order_tibble <- tibble::tibble(temp_sample_id = sample_cols_vector)

                 # Join with the design matrix to filter and reorder
                 sample_order_tibble |>
                     dplyr::inner_join(design_matrix_char_id,
                                      by = c("temp_sample_id" = sample_id_col_name))
             }, error = function(e) {
                 warning(sprintf("cleanDesignMatrix: Error during inner_join: %s. Returning object unchanged.", e$message))
                 return(NULL) # Signal error
             })

             if(is.null(cleaned_design_matrix)) {
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
          })

