#' Get Negative Control Features using ANOVA (Lipids)
#'
#' Identifies potential negative control features (lipids) for RUV correction
#' based on ANOVA across a specified grouping variable. Features are ranked by a
#' one-way ANOVA across the grouping variable: features with a Benjamini-Hochberg
#' q-value above the cutoff, observed in at least two groups, ranked by p-value,
#' highest first, ties broken by feature identifier, so the table's row order
#' does not decide (features observed in one group score p = 1 in the family and
#' are not chosen).
#'
#' This method iterates through each assay in the `LipidomicsAssayData` object. Samples of
#' groups whose name contains "pool" are left out of the ANOVA. A grouping column
#' absent from the design, or fewer than two groups left after leaving out the
#' pools, is refused with an error of class `multischolar_ruv_grouping_unusable`.
#' Records made before this rule (which took the first rows of the table) replay
#' with the selection they recorded; that earlier rule is set only in process,
#' by the replay of such a record, and no argument selects it.
#'
#' @param theObject A `LipidomicsAssayData` object.
#' @param ruv_grouping_variable Character string. The column name in the
#'   `design_matrix` to use for grouping in the ANOVA (e.g., biological replicate
#'   group, batch). Defaults are looked up via `checkParamsObjectFunctionSimplify`
#'   using the key `"ruv_grouping_variable"`.
#' @param percentage_as_neg_ctrl Numeric (0-100). The percentage of total features
#'   to select as negative controls based on ANOVA p-value ranking. Overridden by
#'   `num_neg_ctrl` if provided. Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"lipids_percentage_as_neg_ctrl"`.
#' @param num_neg_ctrl Integer. The absolute number of features to select as
#'   negative controls. Overrides `percentage_as_neg_ctrl`. Defaults are looked up via
#'   `checkParamsObjectFunctionSimplify` using the key `"lipids_num_neg_ctrl"`.
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
#'   (lipids) in that assay were selected as negative controls. The vector
#'   is named with the feature IDs.
#'
#' @importFrom methods slot
#' @importFrom purrr map set_names map_lgl
#' @importFrom tibble column_to_rownames as_tibble is_tibble
#' @importFrom dplyr pull select filter all_of any_of mutate across
#' @importFrom rlang sym !!
#' @importFrom logger log_info log_warn
#' @describeIn getNegCtrlMetabAnova Method for LipidomicsAssayData
#' @export
setMethod(
    f = "getNegCtrlMetabAnova",
    signature = "LipidomicsAssayData",
    definition = function(theObject,
                          ruv_grouping_variable = NULL, # These args are now effectively ignored for default resolution
                          percentage_as_neg_ctrl = NULL,
                          num_neg_ctrl = NULL,
                          ruv_qval_cutoff = NULL,
                          ruv_fdr_method = NULL) {

        assay_list <- methods::slot(theObject, "lipid_data")
        lipid_id_col_name <- methods::slot(theObject, "lipid_id_column")
        design_matrix <- methods::slot(theObject, "design_matrix")
        group_id <- methods::slot(theObject, "group_id") # Needed for helper
        sample_id <- methods::slot(theObject, "sample_id")


        # --- Resolve Global Parameters (Mimicking Protein version exactly) ---
        # Get values from object args slot, falling back to hardcoded defaults.
        ruv_grouping_variable_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_grouping_variable", default_value = "group") # Hardcoded default
        ruv_qval_cutoff_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_qval_cutoff", default_value = 0.05) # Hardcoded default
        ruv_fdr_method_final <- checkParamsObjectFunctionSimplify(theObject, "ruv_fdr_method", default_value = "BH") # Hardcoded default


        # Update object args with resolved values (Mimicking Protein version)
        theObject <- updateParamInObject(theObject, "ruv_grouping_variable") # Correct: Only 2 args
        theObject <- updateParamInObject(theObject, "ruv_qval_cutoff") # Correct: Only 2 args
        theObject <- updateParamInObject(theObject, "ruv_fdr_method") # Correct: Only 2 args

        log_info("Starting Negative Control selection using ANOVA for lipids.")
        log_info("Parameters (Resolved):")
        log_info("  - RUV Grouping Variable: {ruv_grouping_variable_final}")
        log_info("  - RUV Q-value Cutoff: {ruv_qval_cutoff_final}")
        log_info("  - RUV FDR Method: {ruv_fdr_method_final}")
        # Percentage/Num are resolved per assay

        if (length(assay_list) == 0) {
            log_warn("No assays found in `lipid_data` slot. Returning empty list.")
            return(list())
        }
        # Ensure list is named
        assay_names <- names(assay_list)
        if (is.null(assay_names)) {
            assay_names <- paste0("Assay_", seq_along(assay_list))
            log_warn("Assay list was unnamed. Using default names.")
        }

        # KIT-242 (VAL-373 Q-A): why an assay was skipped, kept for the RUV search under the current rule.
        skip_collector <- new.env(parent = emptyenv())
        skip_collector$reasons <- list()

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
                if (is.null(assay_tibble)) {
                    return(NULL)
                } # Skip assay
            }
            if (!lipid_id_col_name %in% colnames(assay_tibble)) {
                log_warn("Assay '{assay_name}': Lipid ID column '{lipid_id_col_name}' not found. Skipping.", .logr = TRUE)
                return(NULL)
            }
            if (nrow(assay_tibble) == 0) {
                log_warn("Assay '{assay_name}': Contains zero rows (features). Skipping.", .logr = TRUE)
                return(NULL)
            }

            # --- Identify Sample Columns ---
            design_samples <- tryCatch(as.character(design_matrix[[sample_id]]), error = function(e) {
                character(0)
            })
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
            assay_matrix <- tryCatch(
                {
                    assay_tibble |>
                        dplyr::select(dplyr::all_of(c(lipid_id_col_name, sample_cols))) |>
                        tibble::column_to_rownames(var = lipid_id_col_name) |>
                        as.matrix()
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error converting tibble to matrix: {e$message}. Skipping.", .logr = TRUE)
                    skip_collector$reasons[[assay_name]] <- list(reason = "not_read")
                    return(NULL)
                }
            )
            if (is.null(assay_matrix)) {
                return(NULL)
            }

            assay_matrix[!is.finite(assay_matrix)] <- NA

            # Check for sufficient valid data
            valid_rows <- rowSums(!is.na(assay_matrix)) > 1
            valid_cols <- colSums(!is.na(assay_matrix)) > 1
            if (sum(valid_rows) < 2 || sum(valid_cols) < 2) {
                log_warn("Assay '{assay_name}': Insufficient non-NA data points (<2 features or <2 samples with data) for ANOVA. Skipping.", .logr = TRUE)
                skip_collector$reasons[[assay_name]] <- list(reason = "too_little_data")
                return(NULL)
            }
            assay_matrix_filt <- assay_matrix[valid_rows, valid_cols, drop = FALSE]

            # Filter design matrix to match valid columns in assay_matrix_filt
            design_matrix_filtered <- design_matrix |>
                dplyr::filter(!!rlang::sym(sample_id) %in% colnames(assay_matrix_filt)) |>
                as.data.frame() # Helper might expect data.frame

            # Check if grouping variable has enough levels/samples after filtering
            if (!ruv_grouping_variable_final %in% colnames(design_matrix_filtered)) {
                cli::cli_abort(
                    "Control features cannot be chosen for {assay_name}: the design has no column '{ruv_grouping_variable_final}', so the features cannot be compared between groups. Choose a grouping column that is in the design.",
                    class = c("multischolar_ruv_grouping_unusable", "multischolar_refusal")
                )
            }
            group_counts <- table(design_matrix_filtered[[ruv_grouping_variable_final]])
            if (length(group_counts) < 2 || any(group_counts < 2)) {
                log_warn("Assay '{assay_name}': Insufficient groups ({length(group_counts)}) or samples per group (<2) for ANOVA based on '{ruv_grouping_variable_final}' after filtering. Skipping.", .logr = TRUE)
                all_groups <- unique(as.character(design_matrix[[ruv_grouping_variable_final]]))
                skip_collector$reasons[[assay_name]] <- list(reason = "group_counts", groups = unique(c(names(group_counts)[group_counts < 2],
                    setdiff(all_groups[!is.na(all_groups)], names(group_counts)))))
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
                    theObject, "percentage_as_neg_ctrl",
                    default_value = 10
                ) # Generic key, hardcoded default 10
                log_info("   Assay '{assay_name}': Using percentage from config/default: {percentage_to_use_for_assay}", .logr = TRUE)
                # Update object args only if we resolved from config/default
                # Avoid overwriting if a specific value was passed via function arg
                if (is.null(percentage_as_neg_ctrl)) { # Only update args if function call arg was NULL
                    theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
                }
            }

            # Validate the resolved percentage
            if (!is.numeric(percentage_to_use_for_assay) || length(percentage_to_use_for_assay) != 1 || is.na(percentage_to_use_for_assay) || percentage_to_use_for_assay < 0 || percentage_to_use_for_assay > 100) {
                log_warn("   Assay '{assay_name}': Invalid percentage resolved ({percentage_to_use_for_assay}). Must be numeric between 0 and 100. Skipping assay.", .logr = TRUE)
                return(NULL)
            }

            # Calculate default num_neg_ctrl based on resolved percentage and *filtered* matrix
            # We use the *resolved* percentage for this assay now
            default_num_neg_ctrl <- round(nrow(assay_matrix_filt) * percentage_to_use_for_assay / 100, 0)

            # Resolve num_neg_ctrl (prioritize function arg, then config, then calculated default)
            num_neg_ctrl_assay <- NULL
            if (!is.null(num_neg_ctrl)) { # Check explicit function arg first
                # Add similar logic here if you want num_neg_ctrl to also be per-assay via list/vector
                # For now, assume num_neg_ctrl function arg is single value if provided
                if (is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0) {
                    num_neg_ctrl_assay <- num_neg_ctrl
                    log_info("   Assay '{assay_name}': Using num_neg_ctrl from argument: {num_neg_ctrl_assay}", .logr = TRUE)
                } else {
                    log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl argument provided. Ignoring.", .logr = TRUE)
                }
            }
            if (is.null(num_neg_ctrl_assay)) { # If not provided or invalid in args, check config/default
                num_neg_ctrl_assay <- checkParamsObjectFunctionSimplify(
                    theObject, "num_neg_ctrl",
                    default_value = default_num_neg_ctrl
                ) # Generic key
                log_info("   Assay '{assay_name}': Using num_neg_ctrl from config/calculated default: {num_neg_ctrl_assay}", .logr = TRUE)
                # Update object args only if resolved from config/default and function arg was NULL/invalid
                if (is.null(num_neg_ctrl) || !(is.numeric(num_neg_ctrl) && length(num_neg_ctrl) == 1 && !is.na(num_neg_ctrl) && num_neg_ctrl >= 0)) {
                    theObject <- updateParamInObject(theObject, "num_neg_ctrl")
                }
            }
            # Validate the resolved num_neg_ctrl
            if (!is.numeric(num_neg_ctrl_assay) || length(num_neg_ctrl_assay) != 1 || is.na(num_neg_ctrl_assay) || num_neg_ctrl_assay < 0) {
                log_warn("   Assay '{assay_name}': Invalid num_neg_ctrl resolved ({num_neg_ctrl_assay}). Must be non-negative integer. Skipping assay.", .logr = TRUE)
                return(NULL)
            }
            # Ensure integer
            num_neg_ctrl_assay <- as.integer(num_neg_ctrl_assay)

            log_info("  Assay '{assay_name}': Final Neg Ctrl Count: {num_neg_ctrl_assay} (based on percentage: {percentage_to_use_for_assay}%)", .logr = TRUE)


            # --- Prepare Design Matrix for Helper ---
            # Helper expects rownames = sample IDs and reads the grouping column, so only that column is passed
            # (as the proteomics method does). The earlier rule removed group_id, which is the grouping column on
            # every route, so the helper saw no grouping and took the first rows; it is reached only through the
            # option set when replaying records made with it (KIT-190).
            design_matrix_for_helper <- tryCatch(
                {
                    if (identical(getOption("MultiScholaR.ruv_control_grouping", "kept"), "dropped")) {
                        design_matrix_filtered |>
                            tibble::column_to_rownames(var = sample_id) |>
                            dplyr::select(-dplyr::any_of(group_id)) # Remove group_id if it exists
                    } else {
                        design_matrix_filtered |>
                            tibble::column_to_rownames(var = sample_id) |>
                            dplyr::select(dplyr::all_of(ruv_grouping_variable_final))
                    }
                },
                error = function(e) {
                    log_warn("Assay '{assay_name}': Error preparing design matrix for helper: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL)
                }
            )
            if (is.null(design_matrix_for_helper)) {
                return(NULL)
            }


            # --- Call Helper ---
            # **ASSUMPTION**: getNegCtrlProtAnovaHelper can handle the data
            control_indices_assay <- tryCatch(
                {
                    getNegCtrlProtAnovaHelper(
                        assay_matrix_filt, # Matrix with features as rows, samples as cols
                        design_matrix = design_matrix_for_helper,
                        grouping_variable = ruv_grouping_variable_final,
                        # Pass the specifically resolved percentage and number for *this* assay
                        percentage_as_neg_ctrl = percentage_to_use_for_assay,
                        num_neg_ctrl = num_neg_ctrl_assay,
                        ruv_qval_cutoff = ruv_qval_cutoff_final,
                        ruv_fdr_method = ruv_fdr_method_final,
                        # KIT-190 rule (d) (VAL-110, VAL-126); the earlier rule keeps the earlier ranking.
                        candidate_rule = if (identical(getOption("MultiScholaR.ruv_control_grouping", "kept"), "dropped"))
                            "q_table_order" else "p_id_two_groups"
                    )
                },
                error = function(e) {
                    # Within this method only the helper's pool check raises this class: the method refuses an
                    # absent column itself and passes the column to the helper.
                    if (inherits(e, "multischolar_ruv_grouping_unusable")) {
                        cli::cli_abort(
                            c(
                                "Control features cannot be chosen for {assay_name}: after leaving out the pooled samples, fewer than two groups of '{ruv_grouping_variable_final}' remain.",
                                i = "Choose a grouping column with at least two groups that are not pools."
                            ),
                            parent = e,
                            class = c("multischolar_ruv_grouping_unusable", "multischolar_refusal")
                        )
                    }
                    log_warn("Assay '{assay_name}': Error calling getNegCtrlProtAnovaHelper: {e$message}. Skipping.", .logr = TRUE)
                    return(NULL) # Return NULL for this assay on error
                }
            )

            num_ctrl_selected <- sum(control_indices_assay, na.rm = TRUE)
            log_info("  Assay '{assay_name}': Selected {sum(control_indices_assay, na.rm = TRUE)} control features.", .logr = TRUE)
            return(control_indices_assay)
        })

        # Set names for the list of results
        names(control_features_list) <- assay_names

        # Remove NULL elements (skipped assays)
        final_control_list <- control_features_list[!sapply(control_features_list, is.null)]
        # KIT-242 (VAL-373 Q-A): under the current rule the skipped assays' reasons travel with the list (none: no attribute).
        if (length(skip_collector$reasons) && omicsRuvFallbackCurrent()) attr(final_control_list, "skipped") <- skip_collector$reasons

        log_info("Finished Negative Control selection for {length(final_control_list)} assay(s).")
        return(final_control_list)
    }
)

