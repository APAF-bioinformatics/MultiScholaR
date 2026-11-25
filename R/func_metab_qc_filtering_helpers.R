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
#' @param group_labels Optional. The group of each sample column, in the same order as `sample_columns`. A metabolite is kept only if at least `min_groups` groups each have a value (not missing) in at least `min_values_per_group` of their samples (in every sample, if the group is smaller). Groups named like "pool" and samples without a group are not counted as groups. Without it, only metabolites with no value in any sample are removed, with a warning.
#' @param min_values_per_group Values required per group (a whole number of at least 1; default 2). Used with `group_labels`.
#' @param min_groups Groups required (a whole number from 1 to the number of groups; default 1). Used with `group_labels`.
#' @return A filtered wide data frame containing only the metabolites that pass the filter.
#' @export
metaboliteIntensityFilteringHelper <- function(
  assay_table,
  min_metabolite_intensity_threshold,
  metabolites_proportion_of_samples_below_cutoff,
  metabolite_id_column,
  sample_columns = NULL,
  group_labels = NULL,
  min_values_per_group = NULL,
  min_groups = NULL
) {
    # KIT-202 (VAL-118 (1), VAL-130): the columns counted as samples are the ones given; without them, every numeric
    # column except the identifier, with a classed warning naming them (the package always passes them).
    sample_cols <- omicsIntensityFilterHelperColumns(assay_table, metabolite_id_column, sample_columns)
    num_samples <- length(sample_cols)

    if (num_samples == 0) {
        warning("No numeric sample columns found in the assay table. Returning original table.")
        return(assay_table)
    }

    # KIT-217 (VAL-173 (2)): the per-group count of present values, then today's share (omicsIntensityFilterKeep).
    settings <- omicsIntensityFilterHelperSettings(sample_cols, group_labels, min_values_per_group, min_groups, "metabolites")
    keep_rows <- omicsIntensityFilterKeep(as.matrix(assay_table[, sample_cols, drop = FALSE]), settings$groups,
        min_metabolite_intensity_threshold, metabolites_proportion_of_samples_below_cutoff, settings$min_values_per_group,
        settings$min_groups)

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

    # Filter metabolites by the rule's decision (the share below the cutoff, and the per-group count)
    filtered_assay_table <- metabolites_below_threshold |>
        dplyr::filter(keep_rows) |>
        # Remove the temporary calculation columns
        dplyr::select(-num_below_threshold, -proportion_below_threshold)

    return(filtered_assay_table)
}
