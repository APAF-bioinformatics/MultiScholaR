#' @title Filter Samples by Metabolite Correlation Threshold
#' @name filterSamplesByMetaboliteCorrelationThreshold,MetaboliteAssayData-method
#' @description Removes samples from a MetaboliteAssayData object based on
#'   Pearson correlation thresholds. Samples that do not have at least one
#'   replicate pair with correlation above the threshold are removed.
#' @param theObject A MetaboliteAssayData object
#' @param pearson_correlation_per_pair A list of data frames (one per assay)
#'   containing pair-wise correlation results from \code{pearsonCorForSamplePairs}.
#' @param min_pearson_correlation_threshold A numeric value (0-1). Samples with
#'   correlation below this threshold in their replicate group are removed.
#' @return An updated MetaboliteAssayData object with poorly correlated samples removed.
#' @importFrom dplyr filter select distinct pull
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @export
setMethod(
    f = "filterSamplesByMetaboliteCorrelationThreshold",
    signature = "MetaboliteAssayData",
    definition = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = 0.5) {
        message("+===========================================================================+")
        message("|        Metabolite Sample Filtering by Correlation Threshold (S4)          |")
        message("+===========================================================================+")

        if (is.null(pearson_correlation_per_pair) || !is.list(pearson_correlation_per_pair)) {
            stop("`pearson_correlation_per_pair` must be a list of correlation data frames (one per assay).")
        }

        if (is.null(min_pearson_correlation_threshold) || !is.numeric(min_pearson_correlation_threshold)) {
            stop("`min_pearson_correlation_threshold` must be a numeric value.")
        }

        design_matrix <- theObject@design_matrix
        sample_id_col_name <- theObject@sample_id
        assay_list <- theObject@metabolite_data

        if (length(assay_list) == 0) {
            warning("No assays found in MetaboliteAssayData object.")
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
        theObject@metabolite_data <- filtered_assay_list

        if (length(samples_to_remove_total) > 0) {
            samples_to_remove_unique <- unique(samples_to_remove_total)
            theObject@design_matrix <- design_matrix |>
                dplyr::filter(!(!!rlang::sym(sample_id_col_name) %in% samples_to_remove_unique))
            message(sprintf("Total samples removed across all assays: %d", length(samples_to_remove_unique)))
        }

        message("+===========================================================================+")

        return(theObject)
    }
)

