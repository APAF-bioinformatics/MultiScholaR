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
        # --- Get required info ---
        sample_id_col <- object@sample_id
        metabolite_id_col <- object@metabolite_id_column
        design_matrix <- object@design_matrix
        metabolite_data <- object@metabolite_data

        # --- Basic slot type checks (as before) ---
        if (!is.list(metabolite_data)) {
            errors <- c(errors, "`metabolite_data` must be a list.")
        } else if (length(metabolite_data) > 0 && !all(sapply(metabolite_data, is.data.frame))) {
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
        # Check design matrix first
        if (!is.data.frame(design_matrix)) {
             # Error already added by basic checks, but prevent further processing
        } else if (!sample_id_col %in% colnames(design_matrix)) {
             errors <- c(errors, paste0("`sample_id` column ('", sample_id_col, "') not found in `design_matrix`."))
        } else {
             # Get unique, sorted sample IDs from design matrix (ensure character)
             samples_in_design <- tryCatch(
                 design_matrix[[sample_id_col]] |> as.character() |> unique() |> sort(),
                 error = function(e) { errors <- c(errors, "Error extracting sample IDs from design matrix."); character(0) }
             )

             if (length(samples_in_design) == 0 && length(errors) == 0) {
                 errors <- c(errors, "No valid sample IDs found in design matrix.")
             }

             # Proceed with assay checks only if design matrix looks okay so far
             if(length(metabolite_data) > 0 && length(errors) == 0) {
                 assay_names_vec <- names(metabolite_data)
                 if (is.null(assay_names_vec)) assay_names_vec <- paste0("Assay_", seq_along(metabolite_data))
                 names(metabolite_data) <- assay_names_vec # Ensure the list is named for lapply output

                 # Use lapply to check each assay and collect results/errors
                 assay_check_results <- lapply(assay_names_vec, function(assay_name) {
                      assay_df <- metabolite_data[[assay_name]]
                      assay_errors <- character()

                      # Check metabolite ID column exists
                      if (!metabolite_id_col %in% colnames(assay_df)) {
                           assay_errors <- c(assay_errors, paste0("Assay '", assay_name,"': `metabolite_id_column` ('", metabolite_id_col, "') not found."))
                      }

                      # Identify actual sample columns in the assay
                      assay_colnames <- colnames(assay_df)
                      actual_sample_cols_in_assay <- intersect(assay_colnames, samples_in_design) |> sort()

                      # Store results for later checks
                      list(
                          errors = assay_errors,
                          sample_cols = actual_sample_cols_in_assay
                      )
                 })

                 # Aggregate errors from individual assay checks
                 all_assay_errors <- unlist(lapply(assay_check_results, `[[`, "errors"))
                 errors <- c(errors, all_assay_errors)

                 # Perform cross-assay consistency checks if no individual errors found yet
                 if (length(errors) == 0 && length(assay_check_results) > 1) {
                      first_assay_samples <- assay_check_results[[1]]$sample_cols
                      consistency_check <- sapply(assay_check_results[-1], function(res) {
                           identical(res$sample_cols, first_assay_samples)
                      })
                      if (!all(consistency_check)) {
                           mismatched_assays <- assay_names_vec[c(FALSE, !consistency_check)] # Get names of inconsistent assays
                           errors <- c(errors, paste0("Actual sample columns differ between assays. First mismatch found in: ", mismatched_assays[1]))
                      }
                 }

                 # Perform comparison with design matrix if no errors found yet
                 if (length(errors) == 0 && length(assay_check_results) >= 1) {
                     first_assay_samples <- assay_check_results[[1]]$sample_cols # Get samples from first (or only) assay
                     if (!identical(first_assay_samples, samples_in_design)) {
                          errors <- c(errors, paste0("Sample columns in assays do not exactly match unique sample IDs ('", sample_id_col, "') in `design_matrix`."))
                          # Add more detail:
                          missing_in_assay <- setdiff(samples_in_design, first_assay_samples)
                          extra_in_assay <- setdiff(first_assay_samples, samples_in_design)
                          if(length(missing_in_assay) > 0) errors <- c(errors, paste0("   Samples in design_matrix missing from assay columns: ", paste(utils::head(missing_in_assay, 10), collapse=", "), ifelse(length(missing_in_assay)>10,"...","")))
                          if(length(extra_in_assay) > 0) errors <- c(errors, paste0("   Sample columns in assay not found in design_matrix: ", paste(utils::head(extra_in_assay, 10), collapse=", "), ifelse(length(extra_in_assay)>10,"...","")))
                     }
                 }
             }
        }

        # --- Final Check ---
        if (length(errors) == 0) TRUE else errors
    }
)


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
          signature = c("MetaboliteAssayData", "ANY", "ANY", "ANY", "ANY", "ANY"),
          definition = function(theObject, grouping_variable, shape_variable = NULL, label_column, title = NULL, font_size = 8) {
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


            if(!is.list( assay_list)) {
              assay_list <- list( assay_list)
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
                assay_title <- if (!is.null(title) && title != "") paste(title, "-", assay_name) else ""

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
#'
setMethod(f = "plotRle",
          signature = c("MetaboliteAssayData", "ANY", "ANY", "ANY"),
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

            if(!is.list( assay_list)) {
              assay_list <- list( assay_list)
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
            message(sprintf("   plotRle: yaxis_limit = %s", paste(yaxis_limit, collapse=", ")))
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

              if(!is.list( assay_list)) {
                assay_list <- list( assay_list)
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
                   # assay_title <- paste0( ifelse( is.null(title) | is.na(title) | title == "", "", paste0(title, " - ")), assay_name)

                   tryCatch({
                        # Create PC1 boxplot
                       pc1_box <- ggplot(pca_data_for_plot, aes(x = !!rlang::sym(grouping_variable), y = PC1, fill = !!rlang::sym(grouping_variable))) +
                           geom_boxplot(notch = TRUE) +
                           theme_bw() +
                           labs(title = paste(title, "-", assay_name),
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
                           geom_boxplot(notch = TRUE) +
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
                 design_matrix_filtered <- design_matrix[design_matrix[[sample_id_col_name]] %in% sample_cols, ]

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
                        # # Set fixed y-axis scale, similar to protein version
                        # scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
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


##-----------------------------------------------------------------------------
## Negative Control Selection Methods for MetaboliteAssayData
##-----------------------------------------------------------------------------

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
setMethod(f = "getNegCtrlMetabAnova",
          signature = "MetaboliteAssayData",
          definition = function(theObject,
                                ruv_grouping_variable = NULL, # These args are now effectively ignored for default resolution
                                percentage_as_neg_ctrl = NULL,
                                num_neg_ctrl = NULL,
                                ruv_qval_cutoff = NULL,
                                ruv_fdr_method = NULL) {

            assay_list <- methods::slot(theObject, "metabolite_data")
            metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
            design_matrix <- methods::slot(theObject, "design_matrix")
            group_id <- methods::slot(theObject, "group_id") # Needed for helper
            sample_id <- methods::slot(theObject, "sample_id")

            # --- Resolve Global Parameters (Mimicking Protein version exactly) ---
            # Get values from object args slot, falling back to hardcoded defaults.
            ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = "group") # Hardcoded default
            ruv_qval_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_qval_cutoff", default_value = 0.05)       # Hardcoded default
            ruv_fdr_method_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_fdr_method", default_value = "BH")         # Hardcoded default

            # Update object args with resolved values (Mimicking Protein version)
            theObject <- updateParamInObject(theObject, "ruv_grouping_variable") # Correct: Only 2 args
            theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")       # Correct: Only 2 args
            theObject <- updateParamInObject(theObject, "ruv_fdr_method")        # Correct: Only 2 args

            log_info("Starting Negative Control selection using ANOVA for metabolites.")
            log_info("Parameters (Resolved):")
            log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
            log_info("  - RUV Q-value Cutoff: {ruv_qval_cutoff_final}")
            log_info("  - RUV FDR Method: {ruv_fdr_method_final}")
            # Percentage/Num are resolved per assay

            if (length(assay_list) == 0) {
                log_warn("No assays found in `metabolite_data` slot. Returning empty list.")
                return(list())
            }
            # Ensure list is named
            assay_names <- names(assay_list)
            if (is.null(assay_names)) {
                 assay_names <- paste0("Assay_", seq_along(assay_list))
                 log_warn("Assay list was unnamed. Using default names.")
            }

            # --- Process Each Assay ---
            control_features_list <- lapply(seq_along(assay_list), function(i) {
                assay_name <- assay_names[i]
                assay_tibble <- assay_list[[i]]
                message(sprintf("-- Processing assay for NegCtrl ANOVA: %s", assay_name))

                # --- Basic Checks ---
                if (!tibble::is_tibble(assay_tibble)) {
                    log_warn("Assay '{assay_name}' is not a tibble. Attempting coercion.", .logr = TRUE)
                    assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                        log_warn("Failed to coerce assay '{assay_name}' to tibble: {e$message}. Skipping.", .logr = TRUE)
                        return(NULL)
                    })
                    if (is.null(assay_tibble)) return(NULL) # Skip assay
                }
                 if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                     log_warn("Assay '{assay_name}': Metabolite ID column '{metabolite_id_col_name}' not found. Skipping.", .logr = TRUE)
                     return(NULL)
                 }
                 if (nrow(assay_tibble) == 0) {
                     log_warn("Assay '{assay_name}': Contains zero rows (features). Skipping.", .logr = TRUE)
                     return(NULL)
                 }

                 # --- Identify Sample Columns ---
                 design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) { character(0) })
                 if (length(design_samples) == 0) {
                      log_warn("Assay '{assay_name}': Could not extract valid sample IDs from design matrix column '{sample_id}'. Skipping.", .logr = TRUE)
                      return(NULL)
                 }
                 all_assay_cols <- colnames(assay_tibble)
                 sample_cols <- intersect(all_assay_cols, design_samples)
                 if (length(sample_cols) < 2) { # Need at least 2 samples for ANOVA
                     log_warn("Assay '{assay_name}': Fewer than 2 sample columns identified matching design matrix. Skipping ANOVA.", .logr = TRUE)
                     return(NULL)
                 }
                 # Ensure sample columns are numeric
                 non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
                 if (length(non_numeric_samples) > 0) {
                    log_warn("Assay '{assay_name}': Non-numeric sample columns found: {paste(non_numeric_samples, collapse=', ')}. Attempting coercion.", .logr = TRUE)
                    assay_tibble <- assay_tibble |>
                        dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
                 }

                # --- Prepare Matrix for Helper ---
                assay_matrix <- tryCatch({
                    assay_tibble |>
                       dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |>
                       tibble::column_to_rownames(var = metabolite_id_col_name) |>
                       as.matrix()
                 }, error = function(e) {
                     log_warn("Assay '{assay_name}': Error converting tibble to matrix: {e$message}. Skipping.", .logr = TRUE)
                     return(NULL)
                 })
                 if (is.null(assay_matrix)) return(NULL)

                 assay_matrix[!is.finite(assay_matrix)] <- NA

                 # Check for sufficient valid data
                 valid_rows <- rowSums(!is.na(assay_matrix)) > 1
                 valid_cols <- colSums(!is.na(assay_matrix)) > 1
                 if (sum(valid_rows) < 2 || sum(valid_cols) < 2) {
                      log_warn("Assay '{assay_name}': Insufficient non-NA data points (<2 features or <2 samples with data) for ANOVA. Skipping.", .logr = TRUE)
                      return(NULL)
                 }
                 assay_matrix_filt <- assay_matrix[valid_rows, valid_cols, drop = FALSE]

                 # Filter design matrix to match valid columns in assay_matrix_filt
                 design_matrix_filtered <- design_matrix |>
                    dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix_filt)) |>
                    as.data.frame() # Helper might expect data.frame

                 # Check if grouping variable has enough levels/samples after filtering
                 if (!ruv_grouping_variable_final %in% colnames(design_matrix_filtered)) {
                     log_warn("Assay '{assay_name}': Grouping variable '{ruv_grouping_variable_final}' not found in filtered design matrix. Skipping ANOVA.", .logr = TRUE)
                     return(NULL)
                 }
                  group_counts <- table(design_matrix_filtered[[ruv_grouping_variable_final]])
                  if (length(group_counts) < 2 || any(group_counts < 2)) {
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
                        theObject, "percentage_as_neg_ctrl", default_value = 10) # Generic key, hardcoded default 10
                    log_info("   Assay '{assay_name}': Using percentage from config/default: {percentage_to_use_for_assay}", .logr = TRUE)
                    # Update object args only if we resolved from config/default
                    # Avoid overwriting if a specific value was passed via function arg
                     if(is.null(percentage_as_neg_ctrl)){ # Only update args if function call arg was NULL
                        theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
                     }

                }

                # Validate the resolved percentage
                if (!is.numeric(percentage_to_use_for_assay) || length(percentage_to_use_for_assay) != 1 || is.na(percentage_to_use_for_assay) || percentage_to_use_for_assay < 0 || percentage_to_use_for_assay > 100) {
                    log_warn("   Assay '{assay_name}': Invalid percentage resolved ({percentage_to_use_for_assay}). Must be numeric between 0 and 100. Skipping assay.", .logr=TRUE)
                    return(NULL)
                }

                # Calculate default num_neg_ctrl based on resolved percentage and *filtered* matrix
                # We use the *resolved* percentage for this assay now
                default_num_neg_ctrl <- round(nrow(assay_matrix_filt) * percentage_to_use_for_assay / 100, 0)

                # Resolve num_neg_ctrl (prioritize function arg, then config, then calculated default)
                num_neg_ctrl_assay <- NULL
                if(!is.null(num_neg_ctrl)){ # Check explicit function arg first
                    # Add similar logic here if you want num_neg_ctrl to also be per-assay via list/vector
                    # For now, assume num_neg_ctrl function arg is single value if provided
                    if(is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0){
                        num_neg_ctrl_assay <- num_neg_ctrl
                         log_info("   Assay '{assay_name}': Using num_neg_ctrl from argument: {num_neg_ctrl_assay}", .logr = TRUE)
                    } else {
                        log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl argument provided. Ignoring.", .logr=TRUE)
                    }
                }
                if(is.null(num_neg_ctrl_assay)) { # If not provided or invalid in args, check config/default
                    num_neg_ctrl_assay <- checkParamsObjectFunctionSimplify(
                        theObject, "num_neg_ctrl", default_value = default_num_neg_ctrl) # Generic key
                     log_info("   Assay '{assay_name}': Using num_neg_ctrl from config/calculated default: {num_neg_ctrl_assay}", .logr = TRUE)
                     # Update object args only if resolved from config/default and function arg was NULL/invalid
                      if(is.null(num_neg_ctrl) || !(is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0)){
                         theObject <- updateParamInObject(theObject, "num_neg_ctrl")
                      }
                }
                 # Validate the resolved num_neg_ctrl
                 if (!is.numeric(num_neg_ctrl_assay) || length(num_neg_ctrl_assay) != 1 || is.na(num_neg_ctrl_assay) || num_neg_ctrl_assay < 0) {
                     log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl resolved ({num_neg_ctrl_assay}). Must be non-negative integer. Skipping assay.", .logr=TRUE)
                     return(NULL)
                 }
                 # Ensure integer
                 num_neg_ctrl_assay <- as.integer(num_neg_ctrl_assay)


                log_info("  Assay '{assay_name}': Final Neg Ctrl Count: {num_neg_ctrl_assay} (based on percentage: {percentage_to_use_for_assay}%)", .logr = TRUE)


                # --- Prepare Design Matrix for Helper ---
                # Helper expects rownames = sample IDs, and group_id column removed
                design_matrix_for_helper <- tryCatch({
                    design_matrix_filtered |>
                        tibble::column_to_rownames(var = sample_id) |>
                        dplyr::select(-dplyr::any_of(group_id)) # Remove group_id if it exists
                }, error = function(e) {
                    log_warn("Assay '{assay_name}': Error preparing design matrix for helper: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                })
                if (is.null(design_matrix_for_helper)) return(NULL)


                # --- Call Helper ---
                # **ASSUMPTION**: getNegCtrlProtAnovaHelper can handle the data
                control_indices_assay <- tryCatch({
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
                }, error = function(e) {
                     log_warn("Assay '{assay_name}': Error calling getNegCtrlProtAnovaHelper: {e$message}. Skipping.", .logr = TRUE)
                     return(NULL) # Return NULL for this assay on error
                })

                log_info("  Assay '{assay_name}': Selected {sum(control_indices_assay, na.rm = TRUE)} control features.", .logr = TRUE)
                return(control_indices_assay)
            })

            # Set names for the list of results
            names(control_features_list) <- assay_names

            # Remove NULL elements (skipped assays)
            final_control_list <- control_features_list[!sapply(control_features_list, is.null)]

            log_info("Finished Negative Control selection for {length(final_control_list)} assay(s).")
            return(final_control_list)
          })

#' @describeIn ruvCancor Method for MetaboliteAssayData
#' @importFrom methods slot
#' @importFrom purrr map set_names map_lgl
#' @importFrom tibble column_to_rownames is_tibble as_tibble
#' @importFrom dplyr pull select filter all_of any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom logger log_info log_warn log_error
#' @importFrom mixOmics impute.nipals
#' @export
setMethod(f = "ruvCancor",
          signature = "MetaboliteAssayData",
          definition = function(theObject, ctrl = NULL, num_components_to_impute = NULL, ruv_grouping_variable = NULL) {

            assay_list <- methods::slot(theObject, "metabolite_data")
            metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
            design_matrix <- methods::slot(theObject, "design_matrix")
            sample_id <- methods::slot(theObject, "sample_id")
            # group_id is not directly used here but good practice to extract if needed later

            # --- Resolve Global Parameters ---
            # Use generic keys as per handover doc
            # Default ctrl=NULL means it MUST be provided in args or function call
            ctrl_final <- checkParamsObjectFunctionSimplify(theObject, "ctrl", default_value = NULL)
            num_components_to_impute_final <- checkParamsObjectFunctionSimplify(theObject, "num_components_to_impute", default_value = 2)
            # Default ruv_grouping_variable = NULL, MUST be provided
            ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = NULL)

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
                log_error("Negative control features ('ctrl') must be provided either via function argument or object configuration ('args$ctrl').")
                stop("Missing required 'ctrl' parameter for ruvCancor.")
            }
            if (is.null(ruv_grouping_variable_final)) {
                log_error("RUV grouping variable ('ruv_grouping_variable') must be provided either via function argument or object configuration.")
                stop("Missing required 'ruv_grouping_variable' parameter for ruvCancor.")
            }
            if (!ruv_grouping_variable_final %in% colnames(design_matrix)) {
                log_error("The 'ruv_grouping_variable' ('{ruv_grouping_variable_final}') is not a column in the design matrix.")
                stop(paste0("The 'ruv_grouping_variable = ", ruv_grouping_variable_final, "' is not a column in the design matrix."))
            }
            if (!is.numeric(num_components_to_impute_final) || is.na(num_components_to_impute_final) || num_components_to_impute_final < 1) {
                log_error("Invalid 'num_components_to_impute': {num_components_to_impute_final}. Must be a positive integer.")
                stop(paste0("The num_components_to_impute = ", num_components_to_impute_final, " value is invalid."))
            }


            if (length(assay_list) == 0) {
                log_warn("No assays found in `metabolite_data` slot. Returning empty list.")
                return(list())
            }
            # Ensure list is named
            assay_names <- names(assay_list)
            if (is.null(assay_names)) {
                 assay_names <- paste0("Assay_", seq_along(assay_list))
                 log_warn("Assay list was unnamed. Using default names.")
            }

            # --- Process Each Assay ---
            cancor_plots_list <- lapply(seq_along(assay_list), function(i) {
                 assay_name <- assay_names[i]
                 assay_tibble <- assay_list[[i]]
                 message(sprintf("-- Processing assay for RUV Cancor Plot: %s", assay_name))

                 # --- Basic Checks ---
                 if (!tibble::is_tibble(assay_tibble)) {
                     log_warn("Assay '{assay_name}' is not a tibble. Attempting coercion.", .logr = TRUE)
                     assay_tibble <- tryCatch(tibble::as_tibble(assay_tibble), error = function(e) {
                         log_warn("Failed to coerce assay '{assay_name}' to tibble: {e$message}. Skipping.", .logr = TRUE)
                         return(NULL)
                     })
                     if (is.null(assay_tibble)) return(NULL) # Skip assay
                 }
                  if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                      log_warn("Assay '{assay_name}': Metabolite ID column '{metabolite_id_col_name}' not found. Skipping.", .logr = TRUE)
                      return(NULL)
                  }
                  if (nrow(assay_tibble) == 0) {
                      log_warn("Assay '{assay_name}': Contains zero rows (features). Skipping.", .logr = TRUE)
                      return(NULL)
                  }

                 # --- Identify Sample Columns ---
                 design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) { character(0) })
                 if (length(design_samples) == 0) {
                      log_warn("Assay '{assay_name}': Could not extract valid sample IDs from design matrix column '{sample_id}'. Skipping.", .logr = TRUE)
                      return(NULL)
                 }
                 all_assay_cols <- colnames(assay_tibble)
                 sample_cols <- intersect(all_assay_cols, design_samples)
                 if (length(sample_cols) < 2) {
                     log_warn("Assay '{assay_name}': Fewer than 2 sample columns identified matching design matrix. Skipping RUV cancor plot.", .logr = TRUE)
                     return(NULL)
                 }
                 # Ensure sample columns are numeric
                 non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
                 if (length(non_numeric_samples) > 0) {
                    log_warn("Assay '{assay_name}': Non-numeric sample columns found: {paste(non_numeric_samples, collapse=', ')}. Attempting coercion.", .logr = TRUE)
                    assay_tibble <- assay_tibble |>
                        dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
                 }

                 # --- Prepare Matrix (Features x Samples) ---
                 assay_matrix <- tryCatch({
                     assay_tibble |>
                        dplyr::select(dplyr::all_of(c(metabolite_id_col_name, sample_cols))) |>
                        tibble::column_to_rownames(var = metabolite_id_col_name) |>
                        as.matrix()
                 }, error = function(e) {
                      log_warn("Assay '{assay_name}': Error converting tibble to matrix: {e$message}. Skipping.", .logr = TRUE)
                      return(NULL)
                 })
                 if (is.null(assay_matrix)) return(NULL)

                 assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf first

                 # --- Filter Design Matrix ---
                 # Ensure design matrix matches the actual columns used in the assay_matrix
                 design_matrix_filtered <- design_matrix |>
                    dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix)) |>
                    as.data.frame() # Ensure it's a data.frame if needed

                 # --- Prepare Y (Samples x Features) ---
                 Y_matrix <- t(assay_matrix[, as.character(design_matrix_filtered[[sample_id]]), drop = FALSE]) # Ensure column order matches filtered design

                 # --- Imputation (using mixOmics::impute.nipals) ---
                 if (anyNA(Y_matrix)) {
                      log_info("   Assay '{assay_name}': Missing values detected. Performing NIPALS imputation with {num_components_to_impute_final} components.", .logr = TRUE)
                      Y_imputed <- tryCatch({
                          mixOmics::impute.nipals(Y_matrix, ncomp = num_components_to_impute_final)
                      }, error = function(e) {
                           log_warn("Assay '{assay_name}': Error during NIPALS imputation: {e$message}. Skipping.", .logr = TRUE)
                           return(NULL)
                      })
                      if (is.null(Y_imputed)) return(NULL)
                      Y_final <- Y_imputed
                 } else {
                     Y_final <- Y_matrix
                 }

                 # --- Prepare X (Grouping Variable) ---
                 if (!ruv_grouping_variable_final %in% colnames(design_matrix_filtered)) {
                      # This check should be redundant due to the initial check, but safe
                      log_error("Assay '{assay_name}': Grouping variable '{ruv_grouping_variable_final}' lost after filtering design matrix. This shouldn't happen.", .logr = TRUE)
                      return(NULL)
                 }
                 X_vector <- design_matrix_filtered[[ruv_grouping_variable_final]]

                 # --- Prepare ctl (Control Features Indices/Logical) ---
                 # Resolve control features specific to this assay
                 ctrl_assay <- NULL
                 if (is.list(ctrl_final)) {
                    if (assay_name %in% names(ctrl_final)) {
                        ctrl_assay <- ctrl_final[[assay_name]]
                        log_info("   Assay '{assay_name}': Found assay-specific controls in 'ctrl' list.", .logr = TRUE)
                    } else {
                        log_warn("Assay '{assay_name}': 'ctrl' is a list, but does not contain an element named '{assay_name}'. Skipping.", .logr = TRUE)
                        return(NULL)
                    }
                 } else {
                    # If ctrl_final is not a list, assume it's a global vector (numeric, logical, character)
                    # This maintains backwards compatibility if a global ctrl vector is provided
                    ctrl_assay <- ctrl_final
                    log_info("   Assay '{assay_name}': Using the globally provided 'ctrl' vector.", .logr = TRUE)
                 }

                 if (is.null(ctrl_assay)) {
                    # This case should ideally be caught by the list check above, but as a safeguard:
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
                                    log_warn("Assay '{assay_name}': Some assay features not found in named logical 'ctrl' vector. Skipping.", .logr=TRUE)
                                    return(NULL)
                               }
                               control_indices_assay <- ctrl_assay[feature_match]
                               if (length(control_indices_assay) != length(feature_names_in_assay)){
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
                      log_warn("Assay '{assay_name}': Control indices could not be determined. Skipping.", .logr=TRUE)
                      return(NULL)
                 }
                 num_controls_found <- if (is.logical(control_indices_assay)) sum(control_indices_assay, na.rm=TRUE) else length(control_indices_assay)
                 if (num_controls_found < 5) {
                     log_warn("Assay '{assay_name}': Fewer than 5 negative control features found/specified for this assay ({num_controls_found}). RUV results may be unreliable. Skipping cancor plot.", .logr = TRUE)
                     # Proceeding might still work but is discouraged by the original protein code's check
                     return(NULL) # Skip plot generation as per original check
                 }
                  log_info("   Assay '{assay_name}': Using {num_controls_found} control features for cancor plot.", .logr = TRUE)

                 # --- Call ruv_cancorplot ---
                 cancor_plot_assay <- tryCatch({
                      # Ensure ruv_cancorplot is loaded/available in the environment
                      # Requires Y = samples x features matrix
                      # Requires X = grouping vector (same length as nrow(Y))
                      # Requires ctl = logical/numeric vector identifying control columns in Y
                      ruv_cancorplot(Y = Y_final,
                                     X = X_vector,
                                     ctl = control_indices_assay)
                 }, error = function(e) {
                      log_warn("Assay '{assay_name}': Error calling ruv_cancorplot: {e$message}. Check if the function exists and is loaded correctly. Skipping.", .logr = TRUE)
                      return(NULL) # Return NULL for this assay on error
                 })

                 return(cancor_plot_assay)
            })

            # Set names for the list of plots
            names(cancor_plots_list) <- assay_names

            # Remove NULL elements (skipped assays)
            final_plots_list <- cancor_plots_list[!sapply(cancor_plots_list, is.null)]

            log_info("Finished RUV Canonical Correlation plot generation for {length(final_plots_list)} assay(s).")
            return(final_plots_list)
          })

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
#' @importFrom mixOmics impute.nipals
#' @importFrom stringr str_split
#' @describeIn ruvIII_C_Varying Method for MetaboliteAssayData
#' @export
setMethod(f = "ruvIII_C_Varying",
          signature = "MetaboliteAssayData",
          definition = function(theObject,
                                ruv_grouping_variable = NULL,
                                ruv_number_k = NULL,
                                ctrl = NULL) {

            assay_list <- methods::slot(theObject, "metabolite_data")
            metabolite_id_col_name <- methods::slot(theObject, "metabolite_id_column")
            design_matrix <- methods::slot(theObject, "design_matrix")
            sample_id <- methods::slot(theObject, "sample_id")
            group_id <- methods::slot(theObject, "group_id") # Extract for context, though not directly used below
            technical_replicate_id <- methods::slot(theObject, "technical_replicate_id") # Extract for context

            # --- Resolve Parameters (Prioritize function args, then object args) ---

            # 1. RUV Grouping Variable
            if (!is.null(ruv_grouping_variable)) {
                ruv_grouping_variable_final <- ruv_grouping_variable
                log_info("Using 'ruv_grouping_variable' from function argument: {ruv_grouping_variable_final}")
            } else {
                ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = NULL)
                log_info("Using 'ruv_grouping_variable' from object args or default: {ruv_grouping_variable_final}")
            }

            # 2. RUV Number K (k)
            if (!is.null(ruv_number_k)) {
                ruv_number_k_resolved <- ruv_number_k
                log_info("Using 'ruv_number_k' from function argument.")
            } else {
                ruv_number_k_resolved <- checkParamsObjectFunctionSimplify(theObject, "ruv_number_k", default_value = NULL)
                log_info("Using 'ruv_number_k' from object args or default.")
            }

            # 3. Control Features (ctrl)
            if (!is.null(ctrl)) {
                ctrl_resolved <- ctrl
                log_info("Using 'ctrl' from function argument.")
            } else {
                ctrl_resolved <- checkParamsObjectFunctionSimplify(theObject, "ctrl", default_value = NULL)
                 log_info("Using 'ctrl' from object args or default.")
            }

            # --- Update Object Args (Store the resolved values) ---
            # We store the *final* value used, regardless of source
            # Using direct slot assignment as updateParamInObject seemed problematic
            theObject@args$ruv_grouping_variable <- ruv_grouping_variable_final
            theObject@args$ruv_number_k <- ruv_number_k_resolved
            theObject@args$ctrl <- ctrl_resolved

            # --- Validation (Using the final resolved values) ---
            if (is.null(ruv_grouping_variable_final)) {
                log_error("Missing required parameter 'ruv_grouping_variable'. Must be provided in function call or object args.")
                stop("Missing required 'ruv_grouping_variable'")
            }
            if (!ruv_grouping_variable_final %in% colnames(design_matrix)) {
                 log_error("Resolved 'ruv_grouping_variable' ('{ruv_grouping_variable_final}') not found as a column in the design matrix.", .logr = TRUE)
                 stop("'ruv_grouping_variable' not in design matrix")
            }
            if (is.null(ruv_number_k_resolved)) {
                 log_error("Missing required parameter 'ruv_number_k' (K value). Must be provided in function call or object args.")
                 stop("Missing required 'ruv_number_k'")
            }
            if (is.null(ctrl_resolved)) {
                 log_error("Missing required parameter 'ctrl' (control features). Must be provided in function call or object args.")
                 stop("Missing required 'ctrl'")
            }

            log_info("Starting RUV-III C Varying correction for metabolites.")
            log_info("Parameters (Resolved):")
            log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
            log_info("  - RUV Number K (k): Type '{class(ruv_number_k_resolved)}' (Value(s) resolved per assay)")
            log_info("  - Control Features (ctrl): Type '{class(ctrl_resolved)}' (Value(s) resolved per assay)")

            if(!is.list( assay_list)) {
              assay_list <- list( assay_list)
            }

            if (length(assay_list) == 0) {
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
                log_warn("Assay list contained unnamed or empty elements. Using default names (Assay_...).", immediate. = TRUE)
            }

            # --- Validate structure of k and ctrl if they are lists ---
             is_k_list <- is.list(ruv_number_k_resolved)
             # Check if ctrl is a list, but NOT a data.frame (which could be passed accidentally)
             is_ctrl_list <- is.list(ctrl_resolved) && !is.data.frame(ctrl_resolved)

             # Validate names if k is a list
             if (is_k_list && !all(assay_names %in% names(ruv_number_k_resolved))) {
                 missing_k <- setdiff(assay_names, names(ruv_number_k_resolved))
                 log_error("If 'ruv_number_k' is a list, its names must match assay names. Missing K for: {paste(missing_k, collapse=', ')}", .logr = TRUE)
                 stop("Names in 'ruv_number_k' list do not match assay names.")
             }
             # Validate names if ctrl is a list
             if (is_ctrl_list && !all(assay_names %in% names(ctrl_resolved))) {
                  missing_ctrl <- setdiff(assay_names, names(ctrl_resolved))
                  log_error("If 'ctrl' is a list, its names must match assay names. Missing ctrl for: {paste(missing_ctrl, collapse=', ')}", .logr = TRUE)
                  stop("Names in 'ctrl' list do not match assay names.")
             }
             # Validate type if k is NOT a list (must be single numeric for multiple assays)
             if (!is_k_list && length(assay_list) > 1 && !(is.numeric(ruv_number_k_resolved) && length(ruv_number_k_resolved) == 1 && !is.na(ruv_number_k_resolved))) {
                  log_error("If multiple assays exist, 'ruv_number_k' must be a single non-NA numeric value or a named list.", .logr = TRUE)
                  stop("Invalid format for 'ruv_number_k' for multiple assays.")
             }
             # Validate type if ctrl is NOT a list (must be vector for multiple assays)
             if (!is_ctrl_list && length(assay_list) > 1 && !(is.logical(ctrl_resolved) || is.numeric(ctrl_resolved) || is.character(ctrl_resolved))) {
                  log_error("If multiple assays exist and 'ctrl' is not a list, it must be a logical, numeric, or character vector (applied globally).", .logr = TRUE)
                  stop("Invalid format for 'ctrl' for multiple assays.")
             }


            # --- Process Each Assay ---
            corrected_assay_list <- lapply(seq_along(assay_list), function(i) {
                assay_name <- assay_names[i]
                assay_tibble <- assay_list[[i]]
                message(sprintf("-- Processing assay for RUVIII: %s", assay_name))

                # --- Get Assay-Specific k and ctrl ---
                k_assay <- if (is_k_list) ruv_number_k_resolved[[assay_name]] else ruv_number_k_resolved
                ctrl_assay_input <- if (is_ctrl_list) ctrl_resolved[[assay_name]] else ctrl_resolved

                # Validate k_assay
                if (!is.numeric(k_assay) || length(k_assay) != 1 || is.na(k_assay) || k_assay < 0) {
                    log_warn("Assay '{assay_name}': Invalid K value resolved ({k_assay}). Must be a non-negative integer. Skipping.", .logr = TRUE)
                    return(NULL)
                }
                k_assay <- as.integer(k_assay) # Ensure integer

                # --- Basic Checks & Data Prep ---
                if (!tibble::is_tibble(assay_tibble)) {
                    log_warn("Assay '{assay_name}' is not a tibble. Skipping.", .logr=TRUE); return(NULL)
                }
                if (!metabolite_id_col_name %in% colnames(assay_tibble)) {
                     log_warn("Assay '{assay_name}': ID column '{metabolite_id_col_name}' not found. Skipping.", .logr=TRUE); return(NULL)
                }
                if (nrow(assay_tibble) < 1) {
                    log_warn("Assay '{assay_name}' has no features. Skipping.", .logr=TRUE); return(NULL)
                }

                # Identify sample columns based on design matrix
                design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) { character(0) })
                if (length(design_samples) == 0) { log_warn("Assay '{assay_name}': No valid sample IDs in design matrix. Skipping.", .logr=TRUE); return(NULL) }
                all_assay_cols <- colnames(assay_tibble)
                sample_cols <- intersect(all_assay_cols, design_samples)
                if (length(sample_cols) < 2) { log_warn("Assay '{assay_name}': Fewer than 2 sample columns found. Skipping.", .logr=TRUE); return(NULL) }

                # Ensure sample columns are numeric
                non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
                if (length(non_numeric_samples) > 0) {
                   log_warn("Assay '{assay_name}': Coercing non-numeric sample columns to numeric: {paste(non_numeric_samples, collapse=', ')}", .logr = TRUE)
                   assay_tibble <- assay_tibble |> dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
                }

                # Convert to matrix (features x samples)
                # Handle potential duplicate feature IDs before converting to rownames
                assay_matrix <- tryCatch({
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
                }, error = function(e) {
                    log_warn("Assay '{assay_name}': Error converting to matrix: {e$message}. Skipping.", .logr=TRUE); return(NULL)
                })
                 if (is.null(assay_matrix)) return(NULL)
                 assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf AFTER conversion

                 # Filter design matrix to match actual samples in matrix
                 design_matrix_filtered <- design_matrix |>
                     dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix)) |>
                     as.data.frame() # Ensure data.frame

                 if (nrow(design_matrix_filtered) < 2) { log_warn("Assay '{assay_name}': Fewer than 2 samples remain after filtering design matrix. Skipping.", .logr=TRUE); return(NULL) }
                 if (nrow(assay_matrix) < 1) { log_warn("Assay '{assay_name}': Fewer than 1 feature remains. Skipping.", .logr=TRUE); return(NULL) }


                # --- Prepare Y (Samples x Features) --- NO IMPUTATION ---
                # Ensure column order matches filtered design matrix sample order
                Y_final <- t(assay_matrix[, as.character(design_matrix_filtered[[sample_id]]), drop = FALSE])

                 # Check for NAs *before* RUV-III, as the helper might not handle them
                 if (anyNA(Y_final)) {
                    log_warn("   Assay '{assay_name}': Missing values (NA) detected in data matrix Y *before* RUV. RUVIII_C_Varying might fail or produce unexpected results if it doesn't handle NAs internally. Consider imputation *before* calling ruvIII_C_Varying if needed.", .logr = TRUE)
                 }

                # Check dimensions after transpose
                if(nrow(Y_final) < 2 || ncol(Y_final) < 1) {
                    log_warn("Assay '{assay_name}': Insufficient dimensions after preparing Y matrix. Skipping.", .logr=TRUE); return(NULL)
                }


                # --- Prepare M Matrix (using filtered design matrix) ---
                 M <- tryCatch({
                     getRuvIIIReplicateMatrixHelper(design_matrix_filtered,
                                                  !!rlang::sym(sample_id),
                                                  !!rlang::sym(ruv_grouping_variable_final))
                 }, error = function(e){
                      log_warn("Assay '{assay_name}': Error getting RUV III Replicate Matrix: {e$message}. Skipping.", .logr = TRUE)
                      return(NULL)
                 })
                 if(is.null(M)) return(NULL) # Skip if M fails

                 # Ensure M matrix dimensions match Y_final rows (samples)
                 if(nrow(M) != nrow(Y_final) || !identical(rownames(M), rownames(Y_final))) {
                      log_warn("Assay '{assay_name}': M matrix rownames do not match Y matrix rownames after filtering. Attempting to reorder.", .logr = TRUE)
                      # Attempt to reorder M based on Y_final rownames if possible
                      matched_m_rows <- match(rownames(Y_final), rownames(M))
                      if (anyNA(matched_m_rows)) {
                         log_error("   Cannot reorder M matrix - rownames mismatch. Skipping.")
                         return(NULL)
                      }
                      M <- M[matched_m_rows, , drop=FALSE]
                      log_info("   Reordered M matrix rows to match Y matrix.")
                      if(nrow(M) != nrow(Y_final)) { # Double check after reorder
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
                  if (num_controls_found < 1) {
                      log_warn("Assay '{assay_name}': No control features identified after resolution. Skipping.", .logr = TRUE)
                       return(NULL)
                  }
                  log_info("   Assay '{assay_name}': Using {num_controls_found} control features.", .logr = TRUE)

                 # Get the names of the control features
                 potential_controls_names <- feature_names_in_assay[ctrl_logical_assay]


                # --- Call RUVIII_C_Varying ---
                cln_mat <- tryCatch({
                    # Check if RUVIII_C_Varying exists
                    if (!exists("RUVIII_C_Varying", mode = "function")) {
                         stop("Function 'RUVIII_C_Varying' not found. Ensure it is loaded from its package or source file.")
                    }
                    RUVIII_C_Varying(k = k_assay,
                                   Y = Y_final, # Use data potentially containing NAs
                                   M = M,
                                   toCorrect = colnames(Y_final), # Correct all features
                                   potentialControls = potential_controls_names)
                }, error = function(e) {
                    log_warn("Assay '{assay_name}': Error calling RUVIII_C_Varying: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                })
                if (is.null(cln_mat) || !is.matrix(cln_mat)) {
                    log_warn("Assay '{assay_name}': RUVIII_C_Varying did not return a valid matrix. Skipping.", .logr = TRUE)
                    return(NULL) # Skip assay if RUV fails
                }

                # --- Clean Corrected Matrix ---
                # Transpose result back to Features x Samples for cleaning
                corrected_matrix <- t(cln_mat)

                # Remove features (rows) with no finite values
                valid_features <- rowSums(is.finite(corrected_matrix), na.rm = TRUE) > 0
                corrected_matrix_filt_f <- corrected_matrix[valid_features, , drop = FALSE]
                 if(nrow(corrected_matrix_filt_f) == 0) {
                     log_warn("Assay '{assay_name}': No features remained after removing non-finite rows post-RUV. Skipping.", .logr = TRUE)
                     return(NULL)
                 }

                # Remove samples (columns) with no finite values
                valid_samples <- colSums(is.finite(corrected_matrix_filt_f), na.rm = TRUE) > 0
                corrected_matrix_filt_fs <- corrected_matrix_filt_f[, valid_samples, drop = FALSE]
                if(ncol(corrected_matrix_filt_fs) == 0) {
                     log_warn("Assay '{assay_name}': No samples remained after removing non-finite columns post-RUV. Skipping.", .logr = TRUE)
                     return(NULL)
                 }

                log_info("   Assay '{assay_name}': RUV correction applied. Dimensions before cleaning: {nrow(corrected_matrix)}x{ncol(corrected_matrix)}, After: {nrow(corrected_matrix_filt_fs)}x{ncol(corrected_matrix_filt_fs)}", .logr=TRUE)

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


                reconstructed_tibble <- tryCatch({
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
                    dplyr::relocate(dplyr::all_of(colnames(original_metadata_tibble)), # All metadata cols
                                    dplyr::all_of(colnames(corrected_matrix_filt_fs))) # Remaining sample cols

                     # Check if join resulted in expected columns
                     if(!identical(sort(colnames(final_tibble)), sort(c(colnames(original_metadata_tibble), colnames(corrected_matrix_filt_fs))))) {
                          log_warn("Assay '{assay_name}': Column mismatch after joining corrected data and metadata.", .logr=TRUE)
                          # Potentially return NULL or the corrected_data_tibble only
                     }
                     final_tibble

                }, error = function(e) {
                     log_warn("Assay '{assay_name}': Error reconstructing tibble after RUV correction: {e$message}. Skipping.", .logr = TRUE)
                     return(NULL) # Return NULL on error
                })

                message(sprintf("   Assay '%s' RUV-III correction complete.", assay_name))
                return(reconstructed_tibble)
            })

            # Set names and remove NULLs
            names(corrected_assay_list) <- assay_names
            final_corrected_list <- corrected_assay_list[!sapply(corrected_assay_list, is.null)]

            if(length(final_corrected_list) == 0) {
                log_warn("No assays were successfully processed by RUV-III. Returning original object.")
                return(theObject)
            }

            # Update the slot in the object
            methods::slot(theObject, "metabolite_data") <- final_corrected_list

            # --- Clean Design Matrix ---
            theObject <- tryCatch({
                log_info("Cleaning design matrix to match remaining samples after RUV...")
                cleanDesignMatrix(theObject)
            }, error = function(e) {
                 log_warn("Error running cleanDesignMatrix after RUV correction: {e$message}. Design matrix might not be fully synchronized.", .logr=TRUE)
                 return(theObject)
            })

            log_info("RUV-III correction process finished for {length(final_corrected_list)} assay(s).")
            return(theObject)
          })

#' plot number of significant differentially expressed metabolites
#'@export
setGeneric(name="plotNumSigDiffExpBarPlot",
           def=function(objectsList ) {
             standardGeneric("plotNumSigDiffExpBarPlot")
           },
           signature=c("objectsList"  ))


# MetabolomicsDifferentialAbundanceResults
#'@export
setMethod(f = "plotNumSigDiffExpBarPlot",
          signature = "list",
          definition = function(objectsList) {

           return_object_list <-  purrr::map( objectsList
                        , function(object ) {

                          ## Count the number of up or down significant differentially expressed metabolites.
                          # The contrasts_results_table is already a list of data frames (one per contrast)
                          # So we don't need to wrap it in another list
                          num_sig_de_molecules_first_go <- printCountDeGenesTableHelper(
                            list_of_de_tables = object@contrasts_results_table,
                            list_of_descriptions = names(object@contrasts_results_table)
                          )


                          object@num_sig_diff_exp_bar_plot <- num_sig_de_molecules_first_go$plot

                          object@num_sig_diff_table <- num_sig_de_molecules_first_go$table

                          object
                        })

           return_object_list

 })






#' Plot static volcano plot (without gene names)

#' plot number of significant differentially expressed metabolites
#'@export
setGeneric(name="plotVolcano",
           def=function(objectsList,
                        de_q_val_thresh = 0.05,
                        qvalue_column = "q_value",
                        log2fc_column = "logFC") {
             standardGeneric("plotVolcano")
           },
           signature=c("objectsList"))

#'@export
setMethod(f = "plotVolcano",
          signature = "list",
          definition = function(objectsList
                                , de_q_val_thresh = 0.05
                                , qvalue_column = "fdr_qvalue"
                                , log2fc_column = "logFC") {


            return_object_list <-  purrr::map( objectsList
                                               , function(object ) {

                                                 ## Plot static volcano plot
                                                 static_volcano_plot_data <- object@contrasts_results_table  |>
                                                   bind_rows(.id="comparison") |>
                                                   mutate( lqm = -log10(!!sym(qvalue_column) ) ) |>
                                                   dplyr::mutate(label = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "Significant",
                                                                                    TRUE ~ "Not sig.")) |>
                                                   dplyr::mutate(colour = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "purple",
                                                                                     TRUE ~ "black")) |>
                                                   dplyr::mutate(colour = factor(colour, levels = c("black", "purple"))) |>
                                                   dplyr::mutate( title =  str_split_i( comparison, "=", 1))

                                                 list_of_volcano_plots_tbl <- static_volcano_plot_data |>
                                                   group_by( comparison, title) |>
                                                   nest() |>
                                                    mutate( plot = purrr::map2( data, title, \(x,y) { plotOneVolcanoNoVerticalLines(x, y
                                                                                                                                   , log_q_value_column = "lqm"
                                                                                                                                   , log_fc_column = log2fc_column) } ) )

                                                  # THE FIX: Extract the 'plot' column to get a LIST of plots,
                                                  # and correctly name the list elements.
                                                  plots_list <- list_of_volcano_plots_tbl$plot
                                                  names(plots_list) <- list_of_volcano_plots_tbl$comparison

                                                  # Assign the LIST of plots to the slot, not the whole table
                                                  object@volcano_plot <- plots_list
                                                  return(object)
                                               })

            return(return_object_list)

          })


# Get the differential expression results in wide format
#'@export
setGeneric(name="getDeResultsWideFormat"
           , def=function(objectsList
                        , qvalue_column = "fdr_qvalue"
                        , raw_pvalue_column = "raw_pvalue"
                        , log2fc_column = "logFC") {
             standardGeneric("getDeResultsWideFormat")
           },
           signature=c("objectsList"))

#'@export
setMethod(f = "getDeResultsWideFormat"
          , signature = "list"
          , definition = function(objectsList
                                , qvalue_column = "fdr_qvalue"
                                , raw_pvalue_column = "raw_pvalue"
                                , log2fc_column = "logFC")  {


            return_object_list <-  purrr::map( objectsList, function(object) {

              # Correctly access the metabolite data from the nested 'theObject' slot.
              # This defensively handles cases where the slot might hold a list of
              # data frames (correct) or a single data frame (incorrect but handled).
              counts_data_slot <- object@theObject@metabolite_data
              counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
              } else {
                counts_data_slot
              }

              id_col_name <- object@theObject@metabolite_id_column

              # Bind the list of data frames into a single tidy data frame
              tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id="comparison") |>
                dplyr::mutate( comparision_short = str_split_i( comparison, "=", 1))

              # Pivot the tidy data frame to a wide format using the correct ID column
              wide_results <- tidy_results |>
                tidyr::pivot_wider(id_cols = c(!!sym(id_col_name)),
                                   names_from = c(comparision_short),
                                   names_sep = ":",
                                   values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))) |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by( !!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
                dplyr::distinct()

              # Assign to the correct slot and return the object
              object@results_table_wide <- wide_results
              return(object)
            })

            return(return_object_list)
          })




# Get the differential expression results in wide format
#'@export
setGeneric(name="getDeResultsLongFormat"
           , def=function(objectsList) {
             standardGeneric("getDeResultsLongFormat")
           },
           signature=c("objectsList"))

#'@export
setMethod(f = "getDeResultsLongFormat"
          , signature = "list"
          , definition = function(objectsList)  {


            return_object_list <-  purrr::map( objectsList, function(object) {

              # Correctly access the metabolite data from the nested 'theObject' slot.
              # This defensively handles cases where the slot might hold a list of
              # data frames (correct) or a single data frame (incorrect but handled).
              counts_data_slot <- object@theObject@metabolite_data
              counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
              } else {
                counts_data_slot
              }

              id_col_name <- object@theObject@metabolite_id_column

              # Bind the list of data frames into a single tidy data frame
              tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id="comparison") |>
                dplyr::mutate( comparision_short = str_split_i( comparison, "=", 1))

              # Pivot the tidy data frame to a wide format using the correct ID column
              long_results <- tidy_results |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by( !!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::distinct()

              print(head( long_results))

              # Assign to the correct slot and return the object
              object@results_table_long <- long_results
              return(object)
            })

            return(return_object_list)
          })





## Create proteomics interactive volcano plot
#' @export
setGeneric(name="plotInteractiveVolcano"
           , def=function(objectsList, anno_tbl = NULL) {
             standardGeneric("plotInteractiveVolcano")
           },
           signature=c("objectsList", "anno_tbl"))

#'@export
setMethod(f = "plotInteractiveVolcano"
          , signature = "list"
          , definition =
            function( objectsList, anno_tbl = NULL ) {

              list_of_objects <- purrr::map( objectsList,
                                             \(de_output_object) {



                                               updateWithQvalue <-function(r_obj) {
                                                 r_obj_output <- r_obj

                                                 # Prepare the data for the interactive volcano plot
                                                 for(coef in seq_len( ncol( de_analysis_long_table$LCMS_Pos@fit.eb$p.value))) {

                                                   r_obj_output$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues
                                                 }

                                                 r_obj_output

                                               }

                                               my_fit_eb <-   updateWithQvalue( de_output_object@fit.eb)

                                               counts_matrix <- de_output_object@theObject@metabolite_data |>
                                                 column_to_rownames( de_output_object@theObject@metabolite_id_column) |>
                                                 as.matrix()
                                               #
                                               # # anno_tbl <- de_output_object@theObject@metabolite_data |>
                                               # #   dplyr::select( !!sym(de_output_object@theObject@metabolite_id_column) ) |>
                                               # #   mutate( random_stuff = NA_real_ )
                                               #
                                               # anno_tbl$random_stuff <- rnorm( n=nrow( anno_tbl))


                                               groups <- data.frame( Run = colnames(counts_matrix) ) |>
                                                 left_join( de_output_object@theObject@design_matrix
                                                           , by=join_by( !!sym(de_output_object@theObject@sample_id) == !!sym(de_output_object@theObject@sample_id) )) |>
                                                 dplyr::pull(genotype_group)


                                               list_of_glimma_objs <- purrr::map( seq_len( ncol( de_analysis_long_table$LCMS_Pos@fit.eb$p.value))

                                                                                  , \(idx) {

                                                                                    anno_tbl_filt <-  de_output_object@contrasts_results_table[[idx]] |>
                                                                                      dplyr::select( !!sym(de_output_object@theObject@metabolite_id_column )) |>
                                                                                      left_join( anno_tbl, by=join_by( !!sym(de_output_object@theObject@metabolite_id_column)  == !!sym(de_output_object@theObject@metabolite_id_column) ) )

                                                                                    Glimma::glimmaVolcano(my_fit_eb
                                                                                                          , coef=idx
                                                                                                          , anno=anno_tbl_filt
                                                                                                          , counts = counts_matrix
                                                                                                          , groups = groups
                                                                                                          , display.columns = colnames(anno_tbl )
                                                                                                          , status=decideTests(my_fit_eb, adjust.method="none")
                                                                                                          , p.adj.method = "none"
                                                                                                          , transform.counts='none')

                                                                                  })

                                               de_output_object@interactive_volcano_plot <- list_of_glimma_objs

                                               de_output_object })

              list_of_objects

            } )



    # htmlwidgets::saveWidget( widget = Glimma::glimmaVolcano(r_obj
    #                                                         , coef=coef
    #                                                         , anno=anno_tbl
    #                                                         , counts = counts_tbl
    #                                                         , groups = groups
    #                                                         , display.columns = colnames(anno_tbl )
    #                                                         , status=decideTests(r_obj, adjust.method="none")
    #                                                         , p.adj.method = "none"
    #                                                         , transform.counts='none'
    # ) #the plotly object
    # , file = file.path( output_dir
    #                     , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
    # , selfcontained = TRUE #creates a single html file
    # )



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create a QC composite figure

#' @export
#' @export
setGeneric(name = "createGridQCMetabolomics",
           def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQCMetabolomics")
           },
           signature = c("theObject", "pca_titles", "density_titles", "rle_titles", "pearson_titles", "save_path", "file_name"))

#' @export
setMethod(f = "createGridQCMetabolomics",
          signature = "GridPlotData",
          definition = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {

            # --- Identify all unique assay names from the plot lists ---
            all_plot_lists <- c(theObject@pca_plots, theObject@density_plots, theObject@rle_plots, theObject@pearson_plots)
            all_plot_lists <- all_plot_lists[sapply(all_plot_lists, function(x) is.list(x) && length(x) > 0)]

            assay_names <- if (length(all_plot_lists) > 0) unique(unlist(lapply(all_plot_lists, names))) else character(0)

            if (length(assay_names) == 0 || all(sapply(assay_names, is.null))) {
                warning("No assays with named plots found. Cannot generate composite QC plot.", immediate. = TRUE)
                return(list())
            }

            # --- Determine the grid layout ---
            # Assume the layout is defined by the number of PCA plot types passed
            num_cols <- length(theObject@pca_plots)

            # --- Loop over each assay to create a composite plot ---
            composite_plots_list <- purrr::map(assay_names, function(current_assay_name) {
                message(sprintf("--- Generating composite QC plot for assay: %s ---", current_assay_name))

                # Helper to extract and prepare plots for the current assay
                prepare_plot_row <- function(plot_groups_list) {
                    plots <- purrr::map(plot_groups_list, ~ .x[[current_assay_name]])
                    # Replace any NULLs with a blank plot to maintain grid alignment
                    lapply(plots, function(p) if(is.null(p)) ggplot() + theme_void() else p)
                }

                # Extract and prepare plots for the current assay
                pca_plots_assay <- prepare_plot_row(theObject@pca_plots)
                density_plots_assay <- prepare_plot_row(theObject@density_plots)
                rle_plots_assay <- prepare_plot_row(theObject@rle_plots)
                pearson_plots_assay <- prepare_plot_row(theObject@pearson_plots)

                # --- Plot creation helper functions ---
                createLabelPlot <- function(title) {
                  ggplot() +
                    annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
                    xlim(0, 1) +
                    theme_void() +
                    theme(plot.margin = margin(5, 5, 5, 5), panel.background = element_blank())
                }
                # These functions now only apply theme, as titles are handled by labels
                applyTheme <- function(plot) {
                    if (inherits(plot, "ggplot") && !is.null(plot$data)) { # Check if it's not an empty plot
                         plot <- plot + theme(text = element_text(size = 15),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank())
                         if("patchwork" %in% class(plot)) { # for density plots
                             plot <- plot & theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
                         }
                    }
                    plot
                }


                # --- Generate and combine plots for the current assay ---
                plots_to_combine <- list()
                add_plot_row <- function(plots, titles, create_fn) {
                    # Only add row if there are titles provided for it
                    if (!is.null(titles) && length(titles) > 0) {
                        list(wrap_plots(lapply(titles, createLabelPlot), ncol = num_cols),
                             wrap_plots(lapply(plots, create_fn), ncol = num_cols))
                    } else {
                        list()
                    }
                }

                plots_to_combine <- c(plots_to_combine, add_plot_row(pca_plots_assay, pca_titles, applyTheme))
                plots_to_combine <- c(plots_to_combine, add_plot_row(density_plots_assay, density_titles, applyTheme))
                plots_to_combine <- c(plots_to_combine, add_plot_row(rle_plots_assay, rle_titles, applyTheme))
                plots_to_combine <- c(plots_to_combine, add_plot_row(pearson_plots_assay, pearson_titles, applyTheme))

                if (length(plots_to_combine) == 0) {
                    warning(paste("No plots to combine for assay:", current_assay_name))
                    return(NULL)
                }

                num_rows <- length(plots_to_combine) / 2
                layout_heights <- rep(c(0.1, 1), num_rows)

                combined_plot <- wrap_plots(plots_to_combine, ncol = 1) + plot_layout(heights = layout_heights)

                if (!is.null(save_path)) {
                  assay_file_name <- paste0(file_name, "_", current_assay_name)
                  sapply(c("png", "pdf", "svg"), function(ext) {
                    ggsave(
                      plot = combined_plot,
                      filename = file.path(save_path, paste0(assay_file_name, ".", ext)),
                      width = 3.5 * num_cols,
                      height = 4 * num_rows
                    )
                  })
                  message(paste("Plots saved for assay '", current_assay_name, "' in", save_path))
                }

                return(combined_plot)
            })

            names(composite_plots_list) <- assay_names
            composite_plots_list[!sapply(composite_plots_list, is.null)]
          })


