# ----------------------------------------------------------------------------
# resolveLipidDuplicateFeaturesByIntensity
# ----------------------------------------------------------------------------
#' Resolve Duplicate Features by Keeping Highest Average Intensity
#'
#' Within an assay tibble, identifies features with duplicate IDs and keeps only
#' the one with the highest average intensity across sample columns.
#'
#' @param assay_tibble A data frame or tibble representing one assay.
#' @param id_col Character string. The name of the column containing the feature IDs.
#' @param sample_cols Character vector. The names of the columns containing quantitative sample data.
#'
#' @return A tibble with duplicate features resolved based on highest average intensity.
#' @keywords internal
#' @importFrom dplyr group_by summarise ungroup filter slice_max select rowwise mutate c_across any_of
#' @importFrom rlang sym !!
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
resolveLipidDuplicateFeaturesByIntensity <- function(assay_tibble, id_col, sample_cols) {
    if (!id_col %in% colnames(assay_tibble)) {
        warning(sprintf("ID column '%s' not found in assay tibble. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    if (length(sample_cols) == 0) {
        warning("No sample columns provided. Returning original tibble.")
        return(assay_tibble)
    }

    # Check for duplicates first
    id_counts <- assay_tibble %>% dplyr::count(!!rlang::sym(id_col), name = "feature_count")
    duplicates_exist <- any(id_counts$feature_count > 1)

    if (!duplicates_exist) {
        # message(sprintf("No duplicates found in ID column '%s'. Returning original tibble.", id_col))
        return(assay_tibble)
    }

    message(sprintf("Resolving duplicates in ID column '%s' by keeping highest average intensity feature...", id_col))

    # Ensure sample columns are numeric for mean calculation
    assay_tibble_numeric <- assay_tibble %>%
        dplyr::mutate(dplyr::across(dplyr::any_of(sample_cols), as.numeric))

    # Calculate average intensity (handle NAs)
    # Using rowwise is more robust to non-numeric columns than converting to matrix first
    resolved_tibble <- assay_tibble_numeric %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            avg_intensity = mean(dplyr::c_across(dplyr::any_of(sample_cols)), na.rm = TRUE)
        ) %>%
        dplyr::ungroup() %>%
        # Handle cases where avg_intensity might be NaN (if all samples are NA)
        dplyr::mutate(avg_intensity = ifelse(is.nan(avg_intensity), -Inf, avg_intensity)) %>%
        # Group by the ID and keep the one with the highest average intensity
        dplyr::group_by(!!rlang::sym(id_col)) %>%
        # slice_max keeps ties by default; with_ties = FALSE ensures only one row per ID
        dplyr::slice_max(order_by = avg_intensity, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        # Remove the temporary average intensity column
        dplyr::select(-avg_intensity)

    # Report how many rows were removed
    rows_removed <- nrow(assay_tibble) - nrow(resolved_tibble)
    if (rows_removed > 0) {
        message(sprintf("Removed %d lower-intensity duplicate feature row(s).", rows_removed))
    }

    return(resolved_tibble)
}

# ----------------------------------------------------------------------------
# calculateLipidPairCorrelation
# ----------------------------------------------------------------------------
# --- Internal Helper for Lipid Pair Correlation --- #
#' Calculate Pearson correlation for a pair of samples from long-format data
#'
#' Internal helper function specifically for lipidomics data structure.
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
calculateLipidPairCorrelation <- function(input_pair_table, feature_id_column, sample_id_column, value_column) {
    # Get the two unique sample IDs from the input table
    sample_ids <- unique(input_pair_table[[sample_id_column]])
    if (length(sample_ids) != 2) {
        warning("calculateLipidPairCorrelation: Input table does not contain exactly two samples.")
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
            warning(sprintf("calculateLipidPairCorrelation: Error in stats::cor for samples %s and %s: %s", sample_x_id, sample_y_id, e$message))
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

