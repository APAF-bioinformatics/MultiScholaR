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

