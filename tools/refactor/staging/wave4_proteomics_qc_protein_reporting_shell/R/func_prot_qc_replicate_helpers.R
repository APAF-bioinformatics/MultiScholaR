# ----------------------------------------------------------------------------
# removeProteinWithOnlyOneReplicate
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @title Remove Proteins with Single Replicate
#' @description Remove proteins that only have data for one technical replicate for all sample.
#' This can be repurposed for removing proteins that only have one biological replicates for all experimental groups.
#' @export
removeProteinWithOnlyOneReplicate <- function(
  input_table,
  samples_id_tbl,
  input_table_sample_id_column = Run,
  sample_id_tbl_sample_id_column = ms_filename,
  replicate_group_column = general_sample_info,
  protein_id_column = Protein.Ids,
  core_utilisation
) {
  # Count the number of technical replicates per sample and protein combination
  num_tech_reps_per_sample_and_protein <- NA
  if (length(which(is.na(core_utilisation))) == 0) {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join(samples_id_tbl, by = join_by({{ input_table_sample_id_column }} == {{ sample_id_tbl_sample_id_column }})) |>
      group_by({{ replicate_group_column }}, {{ protein_id_column }}) |>
      # partition(core_utilisation) |>
      summarise(counts = n()) |>
      # collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join(samples_id_tbl, by = join_by({{ input_table_sample_id_column }} == {{ sample_id_tbl_sample_id_column }})) |>
      group_by({{ replicate_group_column }}, {{ protein_id_column }}) |>
      partition(core_utilisation) |>
      summarise(counts = n()) |>
      collect() |>
      ungroup()
  }

  # Any proteins found in more than one replicates in any patient will be kept for analysis
  removed_proteins_with_only_one_replicate <- input_table |>
    inner_join(
      num_tech_reps_per_sample_and_protein |>
        dplyr::filter(counts > 1) |>
        dplyr::select(-counts, -{{ replicate_group_column }}) |>
        distinct(),
      by = join_by({{ protein_id_column }})
    ) |>
    distinct()

  removed_proteins_with_only_one_replicate
}

