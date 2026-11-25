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
#' @param sample_columns Optional. The sample columns to count. Without it, every numeric column other than the identifier column is taken as a sample: pass the design's samples when the table holds other numeric columns (retention time, m/z, scores).
#' @return A filtered wide data frame containing only the metabolites that pass the filter.
#' @export
metaboliteIntensityFilteringHelper <- function(
  assay_table,
  min_metabolite_intensity_threshold,
  metabolites_proportion_of_samples_below_cutoff,
  metabolite_id_column,
  sample_columns = NULL
) {
    # KIT-202 (VAL-118 (1), VAL-130): the columns counted as samples are the ones given; without them, every numeric
    # column except the identifier, with a classed warning naming them (the package always passes them).
    sample_cols <- omicsIntensityFilterHelperColumns(assay_table, metabolite_id_column, sample_columns)
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
