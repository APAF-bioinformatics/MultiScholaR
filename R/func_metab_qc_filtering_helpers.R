# ----------------------------------------------------------------------------
# metaboliteIntensityFilteringHelper
# ----------------------------------------------------------------------------
#' @title Helper function for metabolite intensity filtering
#' @name metaboliteIntensityFilteringHelper
#' @description Filter metabolites based on an intensity threshold and the proportion of samples below that threshold in a wide-format table.
#' @param assay_table A wide data frame where rows are metabolites and columns include a metabolite identifier and numeric sample intensities.
#' @param min_metabolite_intensity_threshold The calculated minimum intensity value. Metabolites in samples below this threshold are considered 'below threshold'.
#' @param metabolites_proportion_of_samples_below_cutoff The maximum allowed proportion (0 to 1) of samples where a metabolite can be below the threshold. If a metabolite exceeds this proportion, it's removed.
#' @param metabolite_id_column A string specifying the name of the column containing the unique metabolite identifiers.
#' @return A filtered wide data frame containing only the metabolites that pass the filter.
#' @export
metaboliteIntensityFilteringHelper <- function(
  assay_table,
  min_metabolite_intensity_threshold,
  metabolites_proportion_of_samples_below_cutoff,
  metabolite_id_column
) {
    # Identify numeric columns representing sample intensities
    sample_cols <- names(assay_table)[sapply(assay_table, is.numeric)]
    num_samples <- length(sample_cols)

    if (num_samples == 0) {
        warning("No numeric sample columns found in the assay table. Returning original table.")
        return(assay_table)
    }

    # Calculate the number of samples below threshold for each metabolite
    metabolites_below_threshold <- assay_table |>
        # Ensure id column is character for safe rowwise operations if needed
        # mutate({{metabolite_id_column}} := as.character({{metabolite_id_column}})) |>
        rowwise() |>
        mutate(
            num_below_threshold = sum(c_across(all_of(sample_cols)) < min_metabolite_intensity_threshold, na.rm = TRUE),
            proportion_below_threshold = num_below_threshold / num_samples
        ) |>
        ungroup()

    # Filter metabolites based on the proportion cutoff
    filtered_assay_table <- metabolites_below_threshold |>
        dplyr::filter(proportion_below_threshold < metabolites_proportion_of_samples_below_cutoff) |>
        # Remove the temporary calculation columns
        dplyr::select(-num_below_threshold, -proportion_below_threshold)

    return(filtered_assay_table)
}

# ----------------------------------------------------------------------------
# resolveDuplicateFeaturesByIntensity
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
resolveDuplicateFeaturesByIntensity <- function(assay_tibble, id_col, sample_cols) {
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

