#' @export
setMethod(
    f = "getDaResultsWideFormat",
    signature = "list",
    definition = function(objectsList,
                          qvalue_column = "fdr_qvalue",
                          raw_pvalue_column = "raw_pvalue",
                          log2fc_column = "logFC") {
        return_object_list <- purrr::map(objectsList, function(object) {
            # Correctly access the lipid data from the nested 'theObject' slot.
            # This defensively handles cases where the slot might hold a list of
            # data frames (correct) or a single data frame (incorrect but handled).
            counts_data_slot <- object@theObject@lipid_data
            counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
            } else {
                counts_data_slot
            }

            id_col_name <- object@theObject@lipid_id_column

            # Bind the list of data frames into a single tidy data frame
            tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id = "comparison") |>
                dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

            # Pivot the tidy data frame to a wide format using the correct ID column
            wide_results <- tidy_results |>
                tidyr::pivot_wider(
                    id_cols = c(!!sym(id_col_name)),
                    names_from = c(comparision_short),
                    names_sep = ":",
                    values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))
                ) |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
                dplyr::distinct()

            # Assign to the correct slot and return the object
            object@results_table_wide <- wide_results
            return(object)
        })

        return(return_object_list)
    }
)

# Get the differential expression results in wide format
#' @export
setMethod(
    f = "getDaResultsLongFormat",
    signature = "list",
    definition = function(objectsList) {
        return_object_list <- purrr::map(objectsList, function(object) {
            # Correctly access the lipid data from the nested 'theObject' slot.
            # This defensively handles cases where the slot might hold a list of
            # data frames (correct) or a single data frame (incorrect but handled).
            counts_data_slot <- object@theObject@lipid_data
            counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
                counts_data_slot[[1]]
            } else {
                counts_data_slot
            }

            id_col_name <- object@theObject@lipid_id_column

            # Bind the list of data frames into a single tidy data frame
            tidy_results <- object@contrasts_results_table |>
                dplyr::bind_rows(.id = "comparison") |>
                dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

            # Pivot the tidy data frame to a wide format using the correct ID column
            long_results <- tidy_results |>
                # Join with original counts using the correct ID column.
                dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
                dplyr::distinct()

            print(head(long_results))

            # Assign to the correct slot and return the object
            object@results_table_long <- long_results
            return(object)
        })

        return(return_object_list)
    }
)

