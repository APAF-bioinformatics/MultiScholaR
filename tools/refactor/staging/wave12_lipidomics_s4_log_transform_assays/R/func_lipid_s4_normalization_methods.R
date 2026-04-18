#"
#' @description
#' Applies a log2 transformation (log2(x + offset)) to the numeric sample columns
#' in all assay tibbles stored within the `lipid_data` slot of a
#' `LipidomicsAssayData` object.
#'
#' @details
#' This transformation is often applied to reduce skewness and stabilize variance
#' in quantitative omics data, especially before visualization or statistical modeling.
#' Zeros in the data are handled by adding a small `offset` before transformation.
#' Non-numeric columns and the lipid ID column are preserved.
#'
#' @title Log Transform Assays for LipidomicsAssayData
#' @name logTransformAssays,LipidomicsAssayData-method
#' @param theObject A `LipidomicsAssayData` object.
#' @param offset A small positive number to add before log transformation (default: 1).
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
    definition = function(theObject, offset = 1, ...) {
        # --- Input Validation ---
        if (!is.numeric(offset) || length(offset) != 1 || offset <= 0) {
            stop("`offset` must be a single positive numeric value.")
        }
        message(sprintf("Applying log2(x + %s) transformation to assays.", as.character(offset)))

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
        methods::slot(theObject, "lipid_data") <- transformed_assay_list

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

