#' Apply RUV-III Correction with Varying K
#'
#' Applies the RUV-III correction method to each assay within a
#' `LipidomicsAssayData` object. This method accounts for unwanted variation
#' using control features and a replicate structure matrix. It allows for
#' potentially different numbers of factors (`k`) to be removed for each assay.
#'
#' @param theObject A `LipidomicsAssayData` object.
#' @param ruv_grouping_variable Character string. The column name in the
#'   `design_matrix` that defines the replicate structure for RUV-III (e.g.,
#'   biological groups where variation *within* the group is considered noise).
#'   Defaults are looked up via `checkParamsObjectFunctionSimplify` using the
#'   key `"ruv_grouping_variable"`. Must be provided.
#' @param ruv_number_k An integer or a named list/vector.
#'   - If an integer, this number of factors (`k`) is removed from all assays.
#'   - If a named list/vector, the names must correspond to the assay names
#'     in `lipid_data`. The value associated with each name specifies the
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
#' @return A modified `LipidomicsAssayData` object where the `lipid_data`
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
#' @importFrom stringr str_split
#' @describeIn ruvIII_C_Varying Method for LipidomicsAssayData
#' @export
setMethod(
    f = "ruvIII_C_Varying",
    signature = "LipidomicsAssayData",
    definition = function(theObject,
                          ruv_grouping_variable = NULL,
                          ruv_number_k = NULL,
                          ctrl = NULL) {
        message("+---------------------------------------------------------------------------+")
        message("|  DEBUG66: Entering ruvIII_C_Varying (LipidomicsAssayData)                 |")
        message("+---------------------------------------------------------------------------+")

        assay_list <- methods::slot(theObject, "lipid_data")
        lipid_id_col_name <- methods::slot(theObject, "lipid_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        sample_id <- methods::slot(theObject, "sample_id")
        group_id <- methods::slot(theObject, "group_id") # Extract for context, though not directly used below
        technical_replicate_id <- methods::slot(theObject, "technical_replicate_id") # Extract for context

        message(sprintf(
            "   DEBUG66 [ruvIII_C_Varying] Function args: ruv_grouping_variable = %s, ruv_number_k is.null = %s, ctrl is.null = %s",
            ifelse(is.null(ruv_grouping_variable), "NULL", ruv_grouping_variable), is.null(ruv_number_k), is.null(ctrl)
        ))
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Number of assays: %d", length(assay_list)))

        # --- Resolve Parameters (Prioritize function args, then object args) ---

        # 1. RUV Grouping Variable
        if (!is.null(ruv_grouping_variable)) {
            ruv_grouping_variable_final <- ruv_grouping_variable
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_grouping_variable from function arg: '%s'", ruv_grouping_variable_final))
            log_info("Using 'ruv_grouping_variable' from function argument: {ruv_grouping_variable_final}")
        } else {
            ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = NULL)
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_grouping_variable from object/default: '%s'", ifelse(is.null(ruv_grouping_variable_final), "NULL", ruv_grouping_variable_final)))
            log_info("Using 'ruv_grouping_variable' from object args or default: {ruv_grouping_variable_final}")
        }

        # 2. RUV Number K (k)
        if (!is.null(ruv_number_k)) {
            ruv_number_k_resolved <- ruv_number_k
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ruv_number_k from function arg, class = '%s'", class(ruv_number_k_resolved)[1]))
            log_info("Using 'ruv_number_k' from function argument.")
        } else {
            ruv_number_k_resolved <- checkParamsObjectFunctionSimplify(theObject, "ruv_number_k", default_value = NULL)
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Using ruv_number_k from object/default: is.null = %s, class = '%s'",
                is.null(ruv_number_k_resolved), class(ruv_number_k_resolved)[1]
            ))
            log_info("Using 'ruv_number_k' from object args or default.")
        }

        # 3. Control Features (ctrl)
        if (!is.null(ctrl)) {
            ctrl_resolved <- ctrl
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Using ctrl from function arg, class = '%s', is.list = %s", class(ctrl_resolved)[1], is.list(ctrl_resolved)))
            log_info("Using 'ctrl' from function argument.")
        } else {
            ctrl_resolved <- checkParamsObjectFunctionSimplify(theObject, "ctrl", default_value = NULL)
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Using ctrl from object/default: is.null = %s, class = '%s', is.list = %s",
                is.null(ctrl_resolved), class(ctrl_resolved)[1], is.list(ctrl_resolved)
            ))
            log_info("Using 'ctrl' from object args or default.")
        }

        # Debug detailed ctrl info
        if (is.list(ctrl_resolved) && !is.null(ctrl_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ctrl_resolved is a list with names: %s", paste(names(ctrl_resolved), collapse = ", ")))
            for (nm in names(ctrl_resolved)) {
                ctrl_item <- ctrl_resolved[[nm]]
                num_ctrl <- if (is.logical(ctrl_item)) sum(ctrl_item, na.rm = TRUE) else length(ctrl_item)
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] ctrl_resolved[['%s']]: class = '%s', num_ctrl = %d", nm, class(ctrl_item)[1], num_ctrl))
            }
        }
        # Debug detailed k info
        if (is.list(ruv_number_k_resolved) && !is.null(ruv_number_k_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved is a list with names: %s", paste(names(ruv_number_k_resolved), collapse = ", ")))
            for (nm in names(ruv_number_k_resolved)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved[['%s']] = %s", nm, ruv_number_k_resolved[[nm]]))
            }
        } else if (!is.null(ruv_number_k_resolved)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] ruv_number_k_resolved = %s (single value)", ruv_number_k_resolved))
        }

        # --- Update Object Args (Store the resolved values) ---
        # We store the *final* value used, regardless of source
        # Using direct slot assignment as updateParamInObject seemed problematic
        theObject@args$ruv_grouping_variable <- ruv_grouping_variable_final
        theObject@args$ruv_number_k <- ruv_number_k_resolved
        theObject@args$ctrl <- ctrl_resolved

        # --- Validation (Using the final resolved values) ---
        if (is.null(ruv_grouping_variable_final)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_grouping_variable_final is NULL!")
            log_error("Missing required parameter 'ruv_grouping_variable'. Must be provided in function call or object args.")
            stop("Missing required 'ruv_grouping_variable'")
        }
        if (!ruv_grouping_variable_final %in% colnames(design_matrix)) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_grouping_variable '%s' not in design_matrix!", ruv_grouping_variable_final))
            log_error("Resolved 'ruv_grouping_variable' ('{ruv_grouping_variable_final}') not found as a column in the design matrix.", .logr = TRUE)
            stop("'ruv_grouping_variable' not in design matrix")
        }
        if (is.null(ruv_number_k_resolved)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ruv_number_k_resolved is NULL!")
            log_error("Missing required parameter 'ruv_number_k' (K value). Must be provided in function call or object args.")
            stop("Missing required 'ruv_number_k'")
        }
        if (is.null(ctrl_resolved)) {
            message("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl_resolved is NULL!")
            log_error("Missing required parameter 'ctrl' (control features). Must be provided in function call or object args.")
            stop("Missing required 'ctrl'")
        }

        log_info("Starting RUV-III C Varying correction for lipids.")
        log_info("Parameters (Resolved):")
        log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
        log_info("  - RUV Number K (k): Type '{class(ruv_number_k_resolved)}' (Value(s) resolved per assay)")
        log_info("  - Control Features (ctrl): Type '{class(ctrl_resolved)}' (Value(s) resolved per assay)")

        if (!is.list(assay_list)) {
            assay_list <- list(assay_list)
        }

        if (length(assay_list) == 0) {
            message("   DEBUG66 [ruvIII_C_Varying] WARNING - no assays found, returning object unchanged")
            log_warn("No assays found in `lipid_data` slot. Returning object unchanged.")
            return(theObject)
        }
        # Ensure list is named
        assay_names <- names(assay_list)
        if (is.null(assay_names) || any(assay_names == "")) {
            needs_name <- which(is.null(assay_names) | assay_names == "")
            new_names <- paste0("Assay_", seq_along(assay_list))
            assay_names[needs_name] <- new_names[needs_name]
            names(assay_list) <- assay_names
            message("   DEBUG66 [ruvIII_C_Varying] Assay list had unnamed elements, using defaults")
            log_warn("Assay list contained unnamed or empty elements. Using default names (Assay_...).", immediate. = TRUE)
        }
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay names: %s", paste(assay_names, collapse = ", ")))

        # --- Validate structure of k and ctrl if they are lists ---
        is_k_list <- is.list(ruv_number_k_resolved)
        # Check if ctrl is a list, but NOT a data.frame (which could be passed accidentally)
        is_ctrl_list <- is.list(ctrl_resolved) && !is.data.frame(ctrl_resolved)
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] is_k_list = %s, is_ctrl_list = %s", is_k_list, is_ctrl_list))

        # Validate names if k is a list
        if (is_k_list && !all(assay_names %in% names(ruv_number_k_resolved))) {
            missing_k <- setdiff(assay_names, names(ruv_number_k_resolved))
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - k list missing assays: %s", paste(missing_k, collapse = ", ")))
            log_error("If 'ruv_number_k' is a list, its names must match assay names. Missing K for: {paste(missing_k, collapse=', ')}", .logr = TRUE)
            stop("Names in 'ruv_number_k' list do not match assay names.")
        }
        # Validate names if ctrl is a list
        if (is_ctrl_list && !all(assay_names %in% names(ctrl_resolved))) {
            missing_ctrl <- setdiff(assay_names, names(ctrl_resolved))
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl list missing assays: %s", paste(missing_ctrl, collapse = ", ")))
            log_error("If 'ctrl' is a list, its names must match assay names. Missing ctrl for: {paste(missing_ctrl, collapse=', ')}", .logr = TRUE)
            stop("Names in 'ctrl' list do not match assay names.")
        }
        # Validate type if k is NOT a list (must be single numeric for multiple assays)
        if (!is_k_list && length(assay_list) > 1 && !(is.numeric(ruv_number_k_resolved) && length(ruv_number_k_resolved) == 1 && !is.na(ruv_number_k_resolved))) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - k not a list but multiple assays, invalid format"))
            log_error("If multiple assays exist, 'ruv_number_k' must be a single non-NA numeric value or a named list.", .logr = TRUE)
            stop("Invalid format for 'ruv_number_k' for multiple assays.")
        }
        # Validate type if ctrl is NOT a list (must be vector for multiple assays)
        if (!is_ctrl_list && length(assay_list) > 1 && !(is.logical(ctrl_resolved) || is.numeric(ctrl_resolved) || is.character(ctrl_resolved))) {
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] FAIL - ctrl not a list but multiple assays, invalid format"))
            log_error("If multiple assays exist and 'ctrl' is not a list, it must be a logical, numeric, or character vector (applied globally).", .logr = TRUE)
            stop("Invalid format for 'ctrl' for multiple assays.")
        }


        # --- Process Each Assay ---
        corrected_assay_list <- lapply(seq_along(assay_list), function(i) {
            assay_name <- assay_names[i]
            assay_tibble <- assay_list[[i]]
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] === Processing assay: %s ===", assay_name))
            message(sprintf("-- Processing assay for RUVIII: %s", assay_name))

            # --- Get Assay-Specific k and ctrl ---
            k_assay <- if (is_k_list) ruv_number_k_resolved[[assay_name]] else ruv_number_k_resolved
            ctrl_assay_input <- if (is_ctrl_list) ctrl_resolved[[assay_name]] else ctrl_resolved
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': k_assay = %s, ctrl_assay_input class = '%s'", assay_name, k_assay, class(ctrl_assay_input)[1]))
            if (is.logical(ctrl_assay_input)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': ctrl_assay_input sum(TRUE) = %d", assay_name, sum(ctrl_assay_input, na.rm = TRUE)))
            }

            # Validate k_assay
            if (!is.numeric(k_assay) || length(k_assay) != 1 || is.na(k_assay) || k_assay < 0) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - Invalid k_assay = %s", assay_name, k_assay))
                log_warn("Assay '{assay_name}': Invalid K value resolved ({k_assay}). Must be a non-negative integer. Skipping.", .logr = TRUE)
                return(NULL)
            }
            k_assay <- as.integer(k_assay) # Ensure integer

            # --- Basic Checks & Data Prep ---
            if (!tibble::is_tibble(assay_tibble)) {
                log_warn("Assay '{assay_name}' is not a tibble. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (!lipid_id_col_name %in% colnames(assay_tibble)) {
                log_warn("Assay '{assay_name}': ID column '{lipid_id_col_name}' not found. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_tibble) < 1) {
                log_warn("Assay '{assay_name}' has no features. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Identify sample columns based on design matrix
            design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) {
                character(0)
            })
            if (length(design_samples) == 0) {
                log_warn("Assay '{assay_name}': No valid sample IDs in design matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }
            all_assay_cols <- colnames(assay_tibble)
            sample_cols <- intersect(all_assay_cols, design_samples)
            if (length(sample_cols) < 2) {
                log_warn("Assay '{assay_name}': Fewer than 2 sample columns found. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Ensure sample columns are numeric
            non_numeric_samples <- sample_cols[!purrr::map_lgl(assay_tibble[sample_cols], is.numeric)]
            if (length(non_numeric_samples) > 0) {
                log_warn("Assay '{assay_name}': Coercing non-numeric sample columns to numeric: {paste(non_numeric_samples, collapse=', ')}", .logr = TRUE)
                assay_tibble <- assay_tibble |> dplyr::mutate(dplyr::across(dplyr::all_of(non_numeric_samples), as.numeric))
            }

            # Convert to matrix (features x samples)
            # Handle potential duplicate feature IDs before converting to rownames
            assay_matrix <- tryCatch(
                {
                    n_initial <- nrow(assay_tibble)
                    assay_tibble_unique <- assay_tibble |>
                        dplyr::group_by(!!rlang::sym(lipid_id_col_name)) |>
                        dplyr::filter(dplyr::row_number() == 1) |> # Keep only first instance
                        dplyr::ungroup()
                    n_final <- nrow(assay_tibble_unique)
                    if (n_final < n_initial) {
                        log_warn("Assay '{assay_name}': Duplicate feature IDs detected in '{lipid_id_col_name}'. Keeping first instance only ({n_final}/{n_initial} features).", .logr = TRUE)
                    }
                    assay_tibble_unique |>
                        tibble::column_to_rownames(var = lipid_id_col_name) |>
                        dplyr::select(dplyr::all_of(sample_cols)) |> # Select only sample columns
                        as.matrix()
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error converting to matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(NULL)
            }
            assay_matrix[!is.finite(assay_matrix)] <- NA # Handle Inf/-Inf AFTER conversion

            # Filter design matrix to match actual samples in matrix
            design_matrix_filtered <- design_matrix |>
                dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix)) |>
                as.data.frame() # Ensure data.frame

            if (nrow(design_matrix_filtered) < 2) {
                log_warn("Assay '{assay_name}': Fewer than 2 samples remain after filtering design matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_matrix) < 1) {
                log_warn("Assay '{assay_name}': Fewer than 1 feature remains. Skipping.", .logr = TRUE)
                return(NULL)
            }


            # --- Prepare Y (Samples x Features) --- NO IMPUTATION ---
            # Ensure column order matches filtered design matrix sample order
            Y_final <- t(assay_matrix[, as.character(design_matrix_filtered[[sample_id]]), drop = FALSE])

            # Check for NAs *before* RUV-III, as the helper might not handle them
            if (anyNA(Y_final)) {
                log_warn("   Assay '{assay_name}': Missing values (NA) detected in data matrix Y *before* RUV. RUVIII_C_Varying might fail or produce unexpected results if it doesn't handle NAs internally. Consider imputation *before* calling ruvIII_C_Varying if needed.", .logr = TRUE)
            }

            # Check dimensions after transpose
            if (nrow(Y_final) < 2 || ncol(Y_final) < 1) {
                log_warn("Assay '{assay_name}': Insufficient dimensions after preparing Y matrix. Skipping.", .logr = TRUE)
                return(NULL)
            }


            # --- Prepare M Matrix (using filtered design matrix) ---
            M <- tryCatch(
                {
                    getRuvIIIReplicateMatrixHelper(
                        design_matrix_filtered,
                        !!rlang::sym(sample_id),
                        !!rlang::sym(ruv_grouping_variable_final)
                    )
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error getting RUV III Replicate Matrix: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(M)) {
                return(NULL)
            } # Skip if M fails

            # Ensure M matrix dimensions match Y_final rows (samples)
            if (nrow(M) != nrow(Y_final) || !identical(rownames(M), rownames(Y_final))) {
                log_warn("Assay '{assay_name}': M matrix rownames do not match Y matrix rownames after filtering. Attempting to reorder.", .logr = TRUE)
                # Attempt to reorder M based on Y_final rownames if possible
                matched_m_rows <- match(rownames(Y_final), rownames(M))
                if (anyNA(matched_m_rows)) {
                    log_error("   Cannot reorder M matrix - rownames mismatch. Skipping.")
                    return(NULL)
                }
                M <- M[matched_m_rows, , drop = FALSE]
                log_info("   Reordered M matrix rows to match Y matrix.")
                if (nrow(M) != nrow(Y_final)) { # Double check after reorder
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
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': num_controls_found = %d", assay_name, num_controls_found))
            if (num_controls_found < 1) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - no control features", assay_name))
                log_warn("Assay '{assay_name}': No control features identified after resolution. Skipping.", .logr = TRUE)
                return(NULL)
            }
            log_info("   Assay '{assay_name}': Using {num_controls_found} control features.", .logr = TRUE)

            # Get the names of the control features
            potential_controls_names <- feature_names_in_assay[ctrl_logical_assay]
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': potential_controls_names count = %d", assay_name, length(potential_controls_names)))


            # --- Call RUVIII_C_Varying ---
            message(sprintf(
                "   DEBUG66 [ruvIII_C_Varying] Assay '%s': Calling RUVIII_C_Varying with k = %d, Y dims = %dx%d, M dims = %dx%d",
                assay_name, k_assay, nrow(Y_final), ncol(Y_final), nrow(M), ncol(M)
            ))
            cln_mat <- tryCatch(
                {
                    # Check if RUVIII_C_Varying exists
                    if (!exists("RUVIII_C_Varying", mode = "function")) {
                        message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - RUVIII_C_Varying function not found!", assay_name))
                        stop("Function 'RUVIII_C_Varying' not found. Ensure it is loaded from its package or source file.")
                    }
                    RUVIII_C_Varying(
                        k = k_assay,
                        Y = Y_final, # Use data potentially containing NAs
                        M = M,
                        toCorrect = colnames(Y_final), # Correct all features
                        potentialControls = potential_controls_names
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUVIII_C_Varying FAILED - %s", assay_name, e$message))
                    log_warn("Assay '{assay_name}': Error calling RUVIII_C_Varying: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(cln_mat) || !is.matrix(cln_mat)) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': FAIL - RUVIII_C_Varying returned NULL or non-matrix", assay_name))
                log_warn("Assay '{assay_name}': RUVIII_C_Varying did not return a valid matrix. Skipping.", .logr = TRUE)
                return(NULL) # Skip assay if RUV fails
            }
            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUVIII_C_Varying SUCCESS - returned matrix dims = %dx%d", assay_name, nrow(cln_mat), ncol(cln_mat)))

            # --- Clean Corrected Matrix ---
            # Transpose result back to Features x Samples for cleaning
            corrected_matrix <- t(cln_mat)

            # Remove features (rows) with no finite values
            valid_features <- rowSums(is.finite(corrected_matrix), na.rm = TRUE) > 0
            corrected_matrix_filt_f <- corrected_matrix[valid_features, , drop = FALSE]
            if (nrow(corrected_matrix_filt_f) == 0) {
                log_warn("Assay '{assay_name}': No features remained after removing non-finite rows post-RUV. Skipping.", .logr = TRUE)
                return(NULL)
            }

            # Remove samples (columns) with no finite values
            valid_samples <- colSums(is.finite(corrected_matrix_filt_f), na.rm = TRUE) > 0
            corrected_matrix_filt_fs <- corrected_matrix_filt_f[, valid_samples, drop = FALSE]
            if (ncol(corrected_matrix_filt_fs) == 0) {
                log_warn("Assay '{assay_name}': No samples remained after removing non-finite columns post-RUV. Skipping.", .logr = TRUE)
                return(NULL)
            }

            log_info("   Assay '{assay_name}': RUV correction applied. Dimensions before cleaning: {nrow(corrected_matrix)}x{ncol(corrected_matrix)}, After: {nrow(corrected_matrix_filt_fs)}x{ncol(corrected_matrix_filt_fs)}", .logr = TRUE)

            # --- Reconstruct Tibble ---
            # Get original metadata columns relevant to the remaining features
            metadata_cols <- setdiff(colnames(assay_tibble), sample_cols) # All original non-sample columns
            # Filter the *original* tibble to get metadata for rows that remain
            original_metadata_tibble <- assay_tibble |>
                dplyr::filter(!!rlang::sym(lipid_id_col_name) %in% rownames(corrected_matrix_filt_fs)) |>
                dplyr::select(dplyr::all_of(c(lipid_id_col_name, metadata_cols))) # Ensure ID column is selected

            # Ensure metadata IDs are unique before join (should be due to earlier handling, but safe)
            original_metadata_tibble <- original_metadata_tibble |>
                dplyr::distinct(!!rlang::sym(lipid_id_col_name), .keep_all = TRUE)


            reconstructed_tibble <- tryCatch(
                {
                    corrected_data_tibble <- corrected_matrix_filt_fs |>
                        as.data.frame() |>
                        tibble::rownames_to_column(var = lipid_id_col_name) |>
                        tibble::as_tibble()

                    # Join corrected data with filtered original metadata
                    # Ensure join column types match (rownames_to_column is character)
                    original_metadata_tibble_char <- original_metadata_tibble |>
                        dplyr::mutate(!!rlang::sym(lipid_id_col_name) := as.character(!!rlang::sym(lipid_id_col_name)))
                    corrected_data_tibble_char <- corrected_data_tibble |>
                        dplyr::mutate(!!rlang::sym(lipid_id_col_name) := as.character(!!rlang::sym(lipid_id_col_name)))

                    final_tibble <- dplyr::left_join(original_metadata_tibble_char, corrected_data_tibble_char, by = lipid_id_col_name) |>
                        # Ensure original column order (metadata first, then remaining samples)
                        dplyr::relocate(
                            dplyr::all_of(colnames(original_metadata_tibble)), # All metadata cols
                            dplyr::all_of(colnames(corrected_matrix_filt_fs))
                        ) # Remaining sample cols

                    # Check if join resulted in expected columns
                    if (!identical(sort(colnames(final_tibble)), sort(c(colnames(original_metadata_tibble), colnames(corrected_matrix_filt_fs))))) {
                        log_warn("Assay '{assay_name}': Column mismatch after joining corrected data and metadata.", .logr = TRUE)
                        # Potentially return NULL or the corrected_data_tibble only
                    }
                    final_tibble
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error reconstructing tibble after RUV correction: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL) # Return NULL on error
                }
            )

            message(sprintf("   DEBUG66 [ruvIII_C_Varying] Assay '%s': RUV-III correction complete.", assay_name))
            message(sprintf("   Assay '%s' RUV-III correction complete.", assay_name))
            return(reconstructed_tibble)
        })

        # Set names and remove NULLs
        names(corrected_assay_list) <- assay_names
        final_corrected_list <- corrected_assay_list[!sapply(corrected_assay_list, is.null)]
        message(sprintf("   DEBUG66 [ruvIII_C_Varying] final_corrected_list length = %d", length(final_corrected_list)))

        if (length(final_corrected_list) == 0) {
            message("   DEBUG66 [ruvIII_C_Varying] WARNING - no assays successfully processed, returning original object")
            log_warn("No assays were successfully processed by RUV-III. Returning original object.")
            return(theObject)
        }

        # Update the slot in the object
        methods::slot(theObject, "lipid_data") <- final_corrected_list

        # --- Clean Design Matrix ---
        theObject <- tryCatch(
            {
                log_info("Cleaning design matrix to match remaining samples after RUV...")
                cleanDesignMatrix(theObject)
            },
            error = function(e) {
                message(sprintf("   DEBUG66 [ruvIII_C_Varying] cleanDesignMatrix FAILED - %s", e$message))
                log_warn("Error running cleanDesignMatrix after RUV correction: {e$message}. Design matrix might not be fully synchronized.", .logr = TRUE)
                return(theObject)
            }
        )

        message(sprintf("   DEBUG66 [ruvIII_C_Varying] RUV-III correction finished for %d assay(s).", length(final_corrected_list)))
        message("+---------------------------------------------------------------------------+")
        message("|  DEBUG66: Exiting ruvIII_C_Varying                                        |")
        message("+---------------------------------------------------------------------------+")
        log_info("RUV-III correction process finished for {length(final_corrected_list)} assay(s).")
        return(theObject)
    }
)

