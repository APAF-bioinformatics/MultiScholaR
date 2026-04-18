# ----------------------------------------------------------------------------
# removeEmptyRows
# ----------------------------------------------------------------------------
#' Remove rows in the table where the columns specified by the column regular expression pattern are all zero or NA value.
#' @param input_table Input table with columns recording protein abundances for each sample. The name of these columns matches a regular expression pattern, defined by 'col_pattern'. Remove rows with all samples having no protein abundance.
#' @param col_pattern String representing regular expression pattern that matches the name of columns containing the protein abundance values.
#' @param row_id The column name with the row_id, tidyverse style name.
#' @return A data frame with the rows without abundance values removed.
#' @export
removeEmptyRows <- function(input_table, col_pattern, row_id) {
  temp_col_name <- paste0("temp_", as_string(as_name(enquo(row_id))))

  temp_input_table <- input_table |>
    dplyr::mutate(!!rlang::sym(temp_col_name) := row_number())

  sites_to_accept <- temp_input_table |>
    mutate(across(matches(col_pattern, perl = TRUE), \(x){
      (is.na(x) | x == 0)
    })) |>
    dplyr::filter(!if_all(matches(col_pattern, perl = TRUE), \(x){
      x == TRUE
    })) |>
    dplyr::select({{ temp_col_name }})

  ## Removing entries where all the "Reporter intensity corrected" rows are zero
  filtered_table <- temp_input_table |>
    inner_join(sites_to_accept, by = temp_col_name) |>
    dplyr::select(-temp_col_name)

  return(filtered_table)
}

# ----------------------------------------------------------------------------
# removeProteinsWithOnlyOneReplicateHelper
# ----------------------------------------------------------------------------
#' @title Remove Proteins with Single Replicate
#' @description
#' Remove proteins that have been identified in only one replicate across all patients
#' (e.g. identified no more than one relpicate in any patient)
#' @export
removeProteinsWithOnlyOneReplicateHelper <- function(
  input_table,
  samples_id_tbl,
  input_table_sample_id_column = Run,
  sample_id_tbl_sample_id_column = ms_filename,
  replicate_group_column = general_sample_info,
  protein_id_column = Protein.Ids,
  quantity_column = Protein.Normalised,
  core_utilisation
) {
  # Count the number of technical replicates per sample and peptide combination
  num_tech_reps_per_sample_and_protein <- NA
  if (length(which(is.na(core_utilisation))) == 0) {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join(samples_id_tbl, by = join_by({{ input_table_sample_id_column }} == {{ sample_id_tbl_sample_id_column }})) |>
      dplyr::filter(!is.na({{ quantity_column }})) |>
      group_by({{ replicate_group_column }}, {{ protein_id_column }}) |>
      # partition(core_utilisation) |>
      summarise(counts = n()) |>
      # collect() |>
      ungroup()
  } else {
    num_tech_reps_per_sample_and_protein <- input_table |>
      left_join(samples_id_tbl, by = join_by({{ input_table_sample_id_column }} == {{ sample_id_tbl_sample_id_column }})) |>
      dplyr::filter(!is.na({{ quantity_column }})) |>
      group_by({{ replicate_group_column }}, {{ protein_id_column }}) |>
      partition(core_utilisation) |>
      summarise(counts = n()) |>
      collect() |>
      ungroup()
  }

  ## Need to have two or more replicates in at least two groups to be included
  proteins_in_two_or_more_groups_with_two_or_more_replicates <- num_tech_reps_per_sample_and_protein |>
    dplyr::filter(counts > 1) |>
    group_by({{ protein_id_column }}) |>
    summarise(num_groups = n()) |>
    ungroup() |>
    dplyr::filter(num_groups > 1) |>
    dplyr::select(-num_groups) |>
    distinct()


  removed_proteins_with_only_one_replicate <- input_table |>
    inner_join(proteins_in_two_or_more_groups_with_two_or_more_replicates,
      by = join_by({{ protein_id_column }})
    ) |>
    distinct()

  removed_proteins_with_only_one_replicate
}

# ----------------------------------------------------------------------------
# removeRowsWithMissingValues
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that have more than accepted number of missing values per group.
#' @param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#' @param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param sample_id The name of the column in design_matrix table that has the sample ID.
#' @param row_id A unique ID for each row of the 'input_table' variable.
#' @param grouping_variable The name of the column in design_matrix table that has the experimental group.
#' @param max_num_samples_miss_per_group An integer representing the maximum number of samples with missing values per group.
#' @param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#' @param temporary_abundance_column The name of a temporary column, as a string, to keep the abundance value you want to filter upon
#' @return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#' @export
removeRowsWithMissingValues <- function(
  input_table, cols, design_matrix, sample_id, row_id, grouping_variable, max_num_samples_miss_per_group, abundance_threshold,
  temporary_abundance_column = "Abundance"
) {
  abundance_long <- input_table |>
    pivot_longer(
      cols = {{ cols }},
      names_to = as_string(as_name(enquo(sample_id))),
      values_to = temporary_abundance_column
    ) |>
    mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)) |>
    left_join(
      design_matrix |>
        mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)),
      by = as_string(as_name(enquo(sample_id)))
    )

  count_missing_values_per_group <- abundance_long |>
    mutate(is_missing = ifelse(!is.na(!!sym(temporary_abundance_column)) & !!sym(temporary_abundance_column) > abundance_threshold, 0, 1)) |>
    group_by({{ row_id }}, {{ grouping_variable }}) |>
    summarise(num_missing_values = sum(is_missing)) |>
    ungroup()

  remove_rows_temp <- count_missing_values_per_group |>
    dplyr::filter(max_num_samples_miss_per_group < num_missing_values) |>
    dplyr::select(-num_missing_values, -{{ grouping_variable }}) |>
    distinct({{ row_id }})

  filtered_tbl <- input_table |>
    dplyr::anti_join(remove_rows_temp, by = as_string(as_name(enquo(row_id))))

  return(filtered_tbl)
}

# ----------------------------------------------------------------------------
# removeRowsWithMissingValuesPercentHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Remove rows with missing values
#' @param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#' @param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param sample_id The name of the column in design_matrix table that has the sample ID.
#' @param row_id A unique ID for each row of the 'input_table' variable.
#' @param grouping_variable The name of the column in design_matrix table that has the experimental group.
#' @param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a protein .
#' @param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with protein abundance values below threshold.
#' @param temporary_abundance_column The name of a temporary column to keep the abundance value you want to filter upon
#' @param proteins_intensity_cutoff_percentile The percentile of the protein intensity values to be used as the minimum threshold for protein intensity.
#' @return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#' @export
removeRowsWithMissingValuesPercentHelper <- function(
  input_table,
  cols,
  design_matrix,
  sample_id # symbol, e.g. Run
  , row_id # symbol
  , grouping_variable # symbol
  , groupwise_percentage_cutoff = 1,
  max_groups_percentage_cutoff = 50,
  proteins_intensity_cutoff_percentile = 1,
  temporary_abundance_column = "Abundance"
) {
  message("+===========================================================================+")
  message("|  DEBUG66: Entering removeRowsWithMissingValuesPercentHelper (OPTIMIZED)   |")
  message("+===========================================================================+")

  # 1. Setup strings and symbols
  sample_id_str <- rlang::as_string(rlang::ensym(sample_id))
  row_id_str <- rlang::as_string(rlang::ensym(row_id))
  group_var_str <- rlang::as_string(rlang::ensym(grouping_variable))

  message(sprintf(
    "   DEBUG66: Pivoting %d rows x %d cols. Memory: %.1f MB",
    nrow(input_table), ncol(input_table), sum(gc()[, 2])
  ))

  # 2. Prepare minimal design matrix (Select only needed columns to save memory in join)
  design_matrix_minimal <- design_matrix |>
    dplyr::select(!!rlang::ensym(sample_id), !!rlang::ensym(grouping_variable)) |>
    dplyr::mutate(!!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)))

  # 3. Pivot and Join in one pipeline to reduce intermediate copies
  # Use .data[[]] for robust column access within dplyr
  abundance_long <- input_table |>
    tidyr::pivot_longer(
      cols = !all_of(cols),
      names_to = sample_id_str,
      values_to = temporary_abundance_column
    ) |>
    dplyr::mutate(
      !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
      !!rlang::sym(temporary_abundance_column) := dplyr::if_else(
        is.nan(!!rlang::sym(temporary_abundance_column)),
        NA_real_,
        !!rlang::sym(temporary_abundance_column)
      )
    ) |>
    dplyr::left_join(design_matrix_minimal, by = sample_id_str)

  message(sprintf(
    "   DEBUG66: Pivot/Join complete. Long table: %d rows. Memory: %.1f MB",
    nrow(abundance_long), sum(gc()[, 2])
  ))

  # Force garbage collection after big join
  gc()

  # 4. Calculate Threshold
  # Extract vector directly to avoid data.frame overhead for simple stats
  valid_values <- abundance_long[[temporary_abundance_column]]
  valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]

  min_protein_intensity_threshold <- if (length(valid_values) > 0) {
    ceiling(quantile(valid_values, probs = proteins_intensity_cutoff_percentile / 100, na.rm = TRUE))[1]
  } else {
    0
  }
  message(sprintf("   DEBUG66: Threshold calculated: %g", min_protein_intensity_threshold))

  # 5. Count samples per group
  count_values_per_group <- abundance_long |>
    dplyr::distinct(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::count(!!rlang::sym(group_var_str), name = "num_per_group")

  # 6. Calculate missing stats
  # Perform aggregation directly
  count_percent_missing_per_group <- abundance_long |>
    dplyr::mutate(is_missing = dplyr::if_else(
      !is.na(!!rlang::sym(temporary_abundance_column)) &
        !!rlang::sym(temporary_abundance_column) > min_protein_intensity_threshold,
      0, 1
    )) |>
    dplyr::group_by(!!rlang::sym(row_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::summarise(num_missing_per_group = sum(is_missing), .groups = "drop") |>
    dplyr::left_join(count_values_per_group, by = group_var_str) |>
    dplyr::mutate(perc_missing_per_group = num_missing_per_group / num_per_group * 100)

  # 7. Identify proteins to remove
  total_num_of_groups <- nrow(count_values_per_group)

  proteins_to_remove <- count_percent_missing_per_group |>
    dplyr::filter(perc_missing_per_group > groupwise_percentage_cutoff) |>
    dplyr::group_by(!!rlang::sym(row_id_str)) |>
    dplyr::summarise(percent_groups_failed = n() / total_num_of_groups * 100, .groups = "drop") |>
    dplyr::filter(percent_groups_failed > max_groups_percentage_cutoff)

  num_removed <- nrow(proteins_to_remove)
  message(sprintf("   DEBUG66: Proteins to remove: %d", num_removed))

  # 8. Filter original table
  if (num_removed > 0) {
    filtered_tbl <- input_table |>
      dplyr::anti_join(proteins_to_remove, by = row_id_str)
  } else {
    filtered_tbl <- input_table
  }

  message(sprintf(
    "   DEBUG66: Final table: %d rows. Memory: %.1f MB",
    nrow(filtered_tbl), sum(gc()[, 2])
  ))

  # Final GC
  gc()

  message("+===========================================================================+")
  message("|  DEBUG66: Exiting removeRowsWithMissingValuesPercentHelper                |")
  message("+===========================================================================+")

  return(filtered_tbl)
}

