# ============================================================================
# func_metab_norm.R
# ============================================================================
# Purpose: Metabolomics normalization functions
# 
# This file contains functions for metabolomics normalization, including
# log transformation and between-sample normalization. Functions in this
# file are used by metabolomics normalization modules and related workflows.
#
# Functions to extract here:
# - logTransformAssays(): S4 method for log transformation (metabolite)
# - normaliseUntransformedData(): S4 method for untransformed data normalization
# - normaliseBetweenSamples(): S4 method for between-sample normalization
# - getNegCtrlMetabAnova(): S4 method for negative control selection
# - Additional metabolomics normalization helper functions
#
# Dependencies:
# - limma
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: logTransformAssays() (metabolite method)
# Current location: R/metabolite_normalization.R
# Type: S4 method (exportMethods)
# Description: Applies log2 transformation to metabolite assays
# setMethod(f = "logTransformAssays", signature = "MetaboliteAssayData", ...) {
#   # Extract from R/metabolite_normalization.R
# }

# Function 2: normaliseUntransformedData() (metabolite method)
# Current location: R/metabolite_normalization.R
# Type: S4 method (exportMethods)
# Description: Normalizes untransformed metabolite data
# setMethod(f = "normaliseUntransformedData", signature = "MetaboliteAssayData", ...) {
#   # Extract from R/metabolite_normalization.R
# }

# Function 3: normaliseBetweenSamples() (metabolite method)
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Normalizes between samples (shared with proteins)
# setMethod(f = "normaliseBetweenSamples", signature = "MetaboliteAssayData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R (if method exists) or create
# }

# Function 4: getNegCtrlMetabAnova()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets negative control metabolites using ANOVA
# setMethod(f = "getNegCtrlMetabAnova", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }


# ----------------------------------------------------------------------------
# logTransformAssays
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------
# getNegCtrlMetabAnova
# ----------------------------------------------------------------------------
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

