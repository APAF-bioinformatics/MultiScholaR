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
    wide_pair_table <- tryCatch(
        {
            input_pair_table |>
                dplyr::select(!!rlang::sym(feature_id_column), !!rlang::sym(sample_id_column), !!rlang::sym(value_column)) |>
                tidyr::pivot_wider(
                    names_from = !!rlang::sym(sample_id_column),
                    values_from = !!rlang::sym(value_column)
                )
        },
        error = function(e) {
            warning(sprintf("Error pivoting data wider for correlation between %s and %s: %s", sample_x_id, sample_y_id, e$message))
            return(NULL)
        }
    )

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
    cor_result <- tryCatch(
        {
            stats::cor(values_x, values_y, use = "pairwise.complete.obs")
        },
        error = function(e) {
            warning(sprintf("calculateMetabolitePairCorrelation: Error in stats::cor for samples %s and %s: %s", sample_x_id, sample_y_id, e$message))
            return(NA_real_) # Returns NA_real_ on cor error
        }
    )

    # --- Modified Check ---
    # Ensure the result is a single, finite numeric value
    if (length(cor_result) != 1 || !is.numeric(cor_result) || !is.finite(cor_result)) {
        return(NA_real_)
    }
    # --- End Modified Check ---

    return(cor_result)
}

