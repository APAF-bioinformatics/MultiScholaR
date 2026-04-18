#' @title Calculate Pearson Correlation for Sample Pairs
#' @name pearsonCorForSamplePairs,LipidomicsAssayData-method
#' @importFrom purrr map set_names map_df
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join select mutate filter distinct arrange group_by summarise ungroup pull if_else case_when row_number n rename relocate
#' @importFrom tibble add_column
#' @importFrom stringr str_detect
#' @importFrom rlang sym !! :=
#' @export
#' @export
setMethod(
    f = "pearsonCorForSamplePairs",
    signature = "LipidomicsAssayData",
    definition = function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
        message("+---------------------------------------------------------------------------+")
        message("|  DEBUG66: Entering pearsonCorForSamplePairs (LipidomicsAssayData)        |")
        message("+---------------------------------------------------------------------------+")

        # --- Input Validation ---
        # tech_rep_remove_regex can be NULL, checked inside helper/later use
        # correlation_group can be NA, checked below

        design_matrix <- theObject@design_matrix
        sample_id_col_name <- theObject@sample_id
        lipid_id_col_name <- theObject@lipid_id_column
        tech_rep_col_name <- theObject@technical_replicate_id # Default grouping if correlation_group is NA
        assay_list <- theObject@lipid_data

        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] sample_id_col = '%s', lipid_id_col = '%s'", sample_id_col_name, lipid_id_col_name))
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] tech_rep_col_name = '%s', correlation_group = '%s'", tech_rep_col_name, correlation_group))
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Number of assays: %d", length(assay_list)))

        if (length(assay_list) == 0) {
            warning("No assays found in `lipid_data` slot. Returning empty list.")
            return(list())
        }

        # Ensure list is named
        if (is.null(names(assay_list))) {
            names(assay_list) <- paste0("Assay_", seq_along(assay_list))
            warning("Assay list was unnamed. Using default names (Assay_1, Assay_2, ...).")
        }
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay names: %s", paste(names(assay_list), collapse = ", ")))

        # Determine the actual grouping column to use for pairing samples
        replicate_group_column_name <- correlation_group
        if (is.na(correlation_group)) {
            message("   DEBUG66 [pearsonCorForSamplePairs] correlation_group is NA, falling back to tech_rep_col")
            if (is.na(tech_rep_col_name) || !tech_rep_col_name %in% colnames(design_matrix)) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] FAIL - tech_rep_col '%s' not valid", tech_rep_col_name))
                stop("`correlation_group` is NA and `technical_replicate_id` ('", tech_rep_col_name, "') is NA or not found in design_matrix. Cannot determine sample pairing.")
            }
            replicate_group_column_name <- tech_rep_col_name
        } else {
            if (!correlation_group %in% colnames(design_matrix)) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] FAIL - correlation_group '%s' not in design_matrix", correlation_group))
                stop(sprintf("Specified `correlation_group` ('%s') not found in design_matrix.", correlation_group))
            }
        }
        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Using replicate_group_column_name = '%s'", replicate_group_column_name))

        # Resolve tech_rep_remove_regex from config if needed (or use default)
        # Assuming the helper function or subsequent filtering handles NULL regex gracefully (meaning no filtering)
        tech_rep_remove_regex_final <- checkParamsObjectFunctionSimplifyAcceptNull(theObject, "tech_rep_remove_regex", tech_rep_remove_regex) # Allow override
        # theObject <- updateParamInObject(theObject, "tech_rep_remove_regex") # Update object if needed


        # --- Correlation Logic per Assay ---
        correlation_results_list <- purrr::map(seq_along(assay_list), function(i) {
            assay_name <- names(assay_list)[i]
            current_assay_data <- assay_list[[i]]
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] === Processing assay: %s ===", assay_name))

            # --- Correctly identify sample columns based on design matrix ---
            design_samples <- as.character(design_matrix[[sample_id_col_name]]) # Get sample IDs from design matrix
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_samples count = %d", assay_name, length(design_samples)))
            all_assay_cols <- colnames(current_assay_data)
            sample_cols <- intersect(all_assay_cols, design_samples) # Find which design samples are columns in the assay
            metadata_cols <- setdiff(all_assay_cols, sample_cols) # All other columns are metadata/ID
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': sample_cols count = %d, metadata_cols count = %d", assay_name, length(sample_cols), length(metadata_cols)))

            # Ensure the primary lipid ID column is considered metadata
            if (lipid_id_col_name %in% sample_cols) {
                warning(sprintf("Assay '%s': Lipid ID column '%s' is also listed as a sample ID. Check configuration.", assay_name, lipid_id_col_name))
            }
            metadata_cols <- union(metadata_cols, lipid_id_col_name) # Ensure lipid ID is not treated as a sample column
            sample_cols <- setdiff(all_assay_cols, metadata_cols) # Final list of sample columns
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': final sample_cols count = %d", assay_name, length(sample_cols)))

            if (length(sample_cols) < 2) { # Need at least 2 samples for correlation
                warning(sprintf("Assay '%s': Fewer than 2 sample columns found matching design matrix sample IDs. Skipping Pearson correlation.", assay_name))
                return(NULL)
            }
            # --- End Correction ---


            # Check sample consistency (now based on correctly identified sample_cols)
            design_samples_check <- design_matrix[[sample_id_col_name]] # Use original type for check
            missing_samples_in_design <- setdiff(sample_cols, as.character(design_samples_check)) # Compare character versions
            if (length(missing_samples_in_design) > 0) {
                # This condition should theoretically not be met if sample_cols were derived from design_samples,
                # but keeping as a safeguard against type issues or unexpected data.
                warning(sprintf("Assay '%s': Identified sample columns missing in design_matrix (check for type mismatches?): %s. Skipping Pearson correlation.", assay_name, paste(missing_samples_in_design, collapse = ", ")))
                return(NULL)
            }


            # Filter design matrix to match assay samples
            design_matrix_filtered <- design_matrix[design_matrix[[sample_id_col_name]] %in% sample_cols, ]
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_matrix_filtered rows = %d", assay_name, nrow(design_matrix_filtered)))

            # Ensure lipid ID column exists
            if (!lipid_id_col_name %in% colnames(current_assay_data)) {
                warning(sprintf("Assay '%s': Lipid ID column '%s' not found. Skipping Pearson correlation.", assay_name, lipid_id_col_name))
                return(NULL)
            }


            # Prepare long data for helper
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Pivoting to long format...", assay_name))
            assay_long <- current_assay_data |>
                tidyr::pivot_longer(
                    cols = all_of(sample_cols),
                    names_to = sample_id_col_name,
                    values_to = "abundance"
                )
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': assay_long rows = %d", assay_name, nrow(assay_long)))


            # Prepare the design matrix subset for the helper
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Creating design_subset with cols '%s' and '%s'", assay_name, sample_id_col_name, replicate_group_column_name))
            design_subset <- design_matrix_filtered |>
                dplyr::select(!!rlang::sym(sample_id_col_name), !!rlang::sym(replicate_group_column_name))
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_subset rows = %d, cols = %s", assay_name, nrow(design_subset), paste(colnames(design_subset), collapse = ", ")))


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
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Calling getPairsOfSamplesTable...", assay_name))
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': design_subset unique groups = %s", assay_name, paste(unique(design_subset[[replicate_group_column_name]]), collapse = ", ")))
            pairs_for_comparison <- tryCatch(
                {
                    getPairsOfSamplesTable(design_subset, # Contains sample_id and replicate_group_column
                        run_id_column = sample_id_col_name,
                        replicate_group_column = replicate_group_column_name
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': ERROR in getPairsOfSamplesTable: %s", assay_name, e$message))
                    warning(sprintf("Assay '%s': Error getting sample pairs: %s. Skipping correlation.", assay_name, e$message))
                    return(NULL)
                }
            )

            if (is.null(pairs_for_comparison) || nrow(pairs_for_comparison) == 0) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': No valid sample pairs found. Skipping.", assay_name))
                warning(sprintf("Assay '%s': No valid sample pairs found for correlation. Skipping.", assay_name))
                return(NULL)
            }
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': pairs_for_comparison rows = %d", assay_name, nrow(pairs_for_comparison)))

            # Get the names of the columns containing paired sample IDs (e.g., "Run.x", "Run.y")
            run_id_col_x <- paste0(sample_id_col_name, ".x")
            run_id_col_y <- paste0(sample_id_col_name, ".y")

            # Check if these columns exist in the pairs table
            if (!all(c(run_id_col_x, run_id_col_y) %in% colnames(pairs_for_comparison))) {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Missing expected columns. pairs_for_comparison cols = %s", assay_name, paste(colnames(pairs_for_comparison), collapse = ", ")))
                warning(sprintf("Assay '%s': Expected paired sample columns ('%s', '%s') not found in pairs table. Skipping correlation.", assay_name, run_id_col_x, run_id_col_y))
                return(NULL)
            }

            # Calculate correlations as a separate vector first
            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Calculating correlations for %d pairs...", assay_name, nrow(pairs_for_comparison)))
            calculated_correlations <- tryCatch(
                {
                    purrr::map2_dbl(
                        .x = pairs_for_comparison[[run_id_col_x]], # Directly access columns
                        .y = pairs_for_comparison[[run_id_col_y]], # Directly access columns
                        .f = ~ {
                            # Filter the long assay data for the current pair
                            assay_pair_filtered <- assay_long |>
                                dplyr::filter(!!rlang::sym(sample_id_col_name) %in% c(.x, .y))

                            # Call the new lipid-specific helper
                            correlation_val <- calculateLipidPairCorrelation(
                                input_pair_table = assay_pair_filtered,
                                feature_id_column = lipid_id_col_name,
                                sample_id_column = sample_id_col_name,
                                value_column = "abundance"
                            )
                            return(correlation_val) # Explicit return
                        }
                    )
                },
                error = function(e) {
                    message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': ERROR in map2_dbl: %s", assay_name, e$message))
                    warning(sprintf("Assay '%s': Error during map2_dbl correlation calculation: %s. Returning NULL results.", assay_name, e$message))
                    return(NULL) # Return NULL if map2_dbl fails
                }
            )

            # Check if calculation succeeded and add the column
            if (is.null(calculated_correlations)) {
                correlation_results_raw <- NULL # Propagate failure
            } else if (length(calculated_correlations) != nrow(pairs_for_comparison)) {
                warning(sprintf("Assay '%s': Number of calculated correlations (%d) does not match number of pairs (%d). Skipping.", assay_name, length(calculated_correlations), nrow(pairs_for_comparison)))
                correlation_results_raw <- NULL
            } else {
                message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Successfully calculated %d correlations", assay_name, length(calculated_correlations)))
                correlation_results_raw <- pairs_for_comparison |>
                    dplyr::mutate(pearson_correlation = calculated_correlations)
            }
            # ---------------------------------------------------------- #

            # Initialize correlation_results_filtered with raw results (will be filtered if regex provided)
            correlation_results_filtered <- correlation_results_raw

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
                    correlation_results_filtered <- correlation_results_raw
                }
            }
            # If no regex, correlation_results_filtered already holds correlation_results_raw

            message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Assay '%s': Returning %d correlation results", assay_name, nrow(correlation_results_filtered)))
            return(correlation_results_filtered)
        })

        # Set names for the list of results
        names(correlation_results_list) <- names(assay_list)

        # Remove NULL elements (skipped assays)
        non_null_count_before <- sum(!sapply(correlation_results_list, is.null))
        correlation_results_list <- correlation_results_list[!sapply(correlation_results_list, is.null)]

        message(sprintf("   DEBUG66 [pearsonCorForSamplePairs] Finished. Returning %d assay results (removed %d NULL)", length(correlation_results_list), length(assay_list) - non_null_count_before))
        message("+---------------------------------------------------------------------------+")
        message("|  DEBUG66: Exiting pearsonCorForSamplePairs                               |")
        message("+---------------------------------------------------------------------------+")

        return(correlation_results_list)
    }
)

#' @title Plot Pearson Correlation
#' @name plotPearson,LipidomicsAssayData-method
#' @importFrom purrr map set_names
#' @importFrom ggplot2 ggplot aes geom_histogram scale_y_continuous xlab ylab theme element_blank
#' @export
setMethod(
    f = "plotPearson",
    signature = "LipidomicsAssayData",
    definition = function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {
        # Get the list of correlation tibbles (one per assay)
        # tech_rep_remove_regex and correlation_group are passed down
        correlation_list <- pearsonCorForSamplePairs(theObject,
            tech_rep_remove_regex = tech_rep_remove_regex,
            correlation_group = correlation_group
        )

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
            tryCatch(
                {
                    pearson_plot <- correlation_vec |>
                        ggplot(aes(pearson_correlation)) +
                        geom_histogram(breaks = hist_breaks, na.rm = TRUE) +
                        # Set x-axis limits and breaks
                        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
                        # # Set fixed y-axis scale, similar to protein version
                        # scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4), expand = c(0, 0)) +
                        xlab("Pearson Correlation") +
                        ylab("Counts") +
                        # ggtitle(paste(assay_name)) +
                        theme_bw() +
                        theme(
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()
                        )

                    return(pearson_plot)
                },
                error = function(e) {
                    warning(sprintf("Assay '%s': Error creating Pearson histogram: %s. Skipping.", assay_name, e$message))
                    return(NULL)
                }
            )
        })

        # Set names for the list of plots
        names(pearson_plots_list) <- names(correlation_list)

        # Remove NULL elements (skipped assays)
        pearson_plots_list <- pearson_plots_list[!sapply(pearson_plots_list, is.null)]

        return(pearson_plots_list)
    }
)

#' @title Filter Samples by Lipid Correlation Threshold
#' @name filterSamplesByLipidCorrelationThreshold,LipidomicsAssayData-method
#' @description Removes samples from a LipidomicsAssayData object based on
#'   Pearson correlation thresholds. Samples that do not have at least one
#'   replicate pair with correlation above the threshold are removed.
#' @param theObject A LipidomicsAssayData object
#' @param pearson_correlation_per_pair A list of data frames (one per assay)
#'   containing pair-wise correlation results from \code{pearsonCorForSamplePairs}.
#' @param min_pearson_correlation_threshold A numeric value (0-1). Samples with
#'   correlation below this threshold in their replicate group are removed.
#' @return An updated LipidomicsAssayData object with poorly correlated samples removed.
#' @importFrom dplyr filter select distinct pull
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @export
setMethod(
    f = "filterSamplesByLipidCorrelationThreshold",
    signature = "LipidomicsAssayData",
    definition = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = 0.5) {
        message("+---------------------------------------------------------------------------+")
        message("|        Lipid Sample Filtering by Correlation Threshold (S4)          |")
        message("+---------------------------------------------------------------------------+")

        if (is.null(pearson_correlation_per_pair) || !is.list(pearson_correlation_per_pair)) {
            stop("`pearson_correlation_per_pair` must be a list of correlation data frames (one per assay).")
        }

        if (is.null(min_pearson_correlation_threshold) || !is.numeric(min_pearson_correlation_threshold)) {
            stop("`min_pearson_correlation_threshold` must be a numeric value.")
        }

        design_matrix <- theObject@design_matrix
        sample_id_col_name <- theObject@sample_id
        assay_list <- theObject@lipid_data

        if (length(assay_list) == 0) {
            warning("No assays found in LipidomicsAssayData object.")
            return(theObject)
        }

        samples_to_remove_total <- character()

        filtered_assay_list <- purrr::map2(assay_list, pearson_correlation_per_pair, function(current_assay_data, correlation_results) {
            if (is.null(correlation_results) || nrow(correlation_results) == 0) {
                warning("No correlation results provided for assay. Skipping filtering.")
                return(current_assay_data)
            }

            run_id_col_x <- paste0(sample_id_col_name, ".x")
            run_id_col_y <- paste0(sample_id_col_name, ".y")

            if (!all(c(run_id_col_x, run_id_col_y, "pearson_correlation") %in% colnames(correlation_results))) {
                warning("Correlation results table missing expected columns. Skipping filtering.")
                return(current_assay_data)
            }

            all_samples_in_analysis <- correlation_results |>
                tidyr::pivot_longer(cols = c(!!rlang::sym(run_id_col_x), !!rlang::sym(run_id_col_y)), values_to = "sample_id") |>
                dplyr::distinct(sample_id) |>
                dplyr::pull(sample_id)

            passing_pairs <- correlation_results |>
                dplyr::filter(pearson_correlation >= min_pearson_correlation_threshold)

            samples_to_keep <- passing_pairs |>
                tidyr::pivot_longer(cols = c(!!rlang::sym(run_id_col_x), !!rlang::sym(run_id_col_y)), values_to = "sample_id") |>
                dplyr::distinct(sample_id) |>
                dplyr::pull(sample_id)

            samples_to_remove <- setdiff(all_samples_in_analysis, samples_to_keep)
            samples_to_remove_total <<- c(samples_to_remove_total, samples_to_remove)

            if (length(samples_to_remove) > 0) {
                message(sprintf(
                    "  Removing %d samples below correlation threshold: %s",
                    length(samples_to_remove), paste(samples_to_remove, collapse = ", ")
                ))
                cols_to_keep <- setdiff(colnames(current_assay_data), samples_to_remove)
                current_assay_data <- current_assay_data[, cols_to_keep, drop = FALSE]
            } else {
                message("  No samples below correlation threshold.")
            }

            return(current_assay_data)
        })

        names(filtered_assay_list) <- names(assay_list)
        theObject@lipid_data <- filtered_assay_list

        if (length(samples_to_remove_total) > 0) {
            samples_to_remove_unique <- unique(samples_to_remove_total)
            theObject@design_matrix <- design_matrix |>
                dplyr::filter(!(!!rlang::sym(sample_id_col_name) %in% samples_to_remove_unique))
            message(sprintf("Total samples removed across all assays: %d", length(samples_to_remove_unique)))
        }

        message("+---------------------------------------------------------------------------+")
        return(theObject)
    }
)
