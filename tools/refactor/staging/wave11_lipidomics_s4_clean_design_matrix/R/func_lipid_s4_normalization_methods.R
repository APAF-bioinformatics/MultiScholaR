#' @title Clean Design Matrix for LipidomicsAssayData
#' @name cleanDesignMatrix,LipidomicsAssayData-method
#' @importFrom dplyr inner_join select rename filter all_of any_of
#' @importFrom rlang sym !!
#' @importFrom methods slot
#' @importFrom tibble tibble
#' @export
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

