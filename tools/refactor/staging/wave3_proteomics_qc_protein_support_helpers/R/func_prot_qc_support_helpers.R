# ----------------------------------------------------------------------------
# avgReplicateProteinIntensity
# ----------------------------------------------------------------------------
#' @title Average Replicate Protein Intensity
#' @description
#' Protein average values from replicate samples
#' @export
avgReplicateProteinIntensity <- function(
  input_table,
  metadata_table,
  protein_id_column = protein_id_column,
  input_table_sample_id_column = Run,
  sample_id_tbl_sample_id_column = Run,
  replicate_group_column = collaborator_patient_id,
  quantity_column = Log2.Protein.Imputed,
  avg_quantity_column = Avg.Log2.Protein.Imputed
) {
  avg_log2_protein_intensity_imputed <- input_table |>
    inner_join(metadata_table,
      by = join_by({{ input_table_sample_id_column }} == {{ sample_id_tbl_sample_id_column }})
    ) |>
    group_by({{ protein_id_column }}, {{ replicate_group_column }}) |>
    summarise({{ avg_quantity_column }} := mean({{ quantity_column }}, na.rm = TRUE)) |>
    ungroup()

  avg_log2_protein_intensity_imputed
}

# ----------------------------------------------------------------------------
# calculatePercentMissingPeptidePerReplicate
# ----------------------------------------------------------------------------
#' calculatePercentMissingPeptidePerReplicate
#' @description Calculate percentage of peptides from each sample that is missing and merge with metadata
#' @export
calculatePercentMissingPeptidePerReplicate <- function(
  input_table,
  metadata_table,
  protein_id_column = Protein.Ids,
  intensity_column = Peptide.Normalised,
  replicate_id_column = Run,
  peptide_sequence_column = Stripped.Sequence
) {
  # Total number of peptides with values per run
  total_num_of_peptides_with_values_per_run <- input_table |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }})) |>
    dplyr::filter(!is.na({{ intensity_column }})) |>
    group_by({{ replicate_id_column }}) |>
    summarise(counts = n()) |>
    ungroup()

  # Total number of peptides
  total_num_of_peptides <- input_table |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }})) |>
    distinct({{ protein_id_column }}, {{ peptide_sequence_column }}) |>
    nrow()

  percent_missing_per_run <- total_num_of_peptides_with_values_per_run |>
    mutate(percent_missing = (1 - (counts / total_num_of_peptides)) * 100) |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }}))

  return(percent_missing_per_run)
}

# ----------------------------------------------------------------------------
# calculatePercentMissingProteinPerReplicate
# ----------------------------------------------------------------------------
#' calculatePercentMissingProteinPerReplicate
#' @description Calculate percentage of proteins from each sample that is missing and merge with metadata
#' @export
calculatePercentMissingProteinPerReplicate <- function(
  input_table,
  metadata_table,
  protein_id_column = Protein.Ids,
  intensity_column = Log2.Protein.Imputed,
  replicate_id_column = Run
) {
  # Total number of peptides with values per run
  total_num_of_proteins_with_values_per_run <- input_table |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }})) |>
    dplyr::filter(!is.na({{ intensity_column }})) |>
    group_by({{ replicate_id_column }}) |>
    summarise(num_proteins_with_values = n()) |>
    ungroup()

  # Total number of peptides
  total_num_of_proteins <- input_table |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }})) |>
    distinct({{ protein_id_column }}) |>
    nrow()

  percent_missing_per_run <- total_num_of_proteins_with_values_per_run |>
    mutate(percent_missing = (1 - (num_proteins_with_values / total_num_of_proteins)) * 100) |>
    left_join(metadata_table, by = join_by({{ replicate_id_column }}))

  return(percent_missing_per_run)
}

# ----------------------------------------------------------------------------
# calculatePercentMissingPerProtein
# ----------------------------------------------------------------------------
#' @export
calculatePercentMissingPerProtein <- function(
  intensity_wide_table,
  protein_id = "uniprot_acc",
  pattern = !tidyselect::matches(protein_id),
  experimental_design_table,
  names_to = "sample_collaborator_sample_id",
  values_to = "Avg.Log2.Protein.Imputed",
  is_missing_column = is_missing
) {
  # print(deparse1(substitute(!!sym({{protein_id}}) )) )

  intensity_long_table <- intensity_wide_table |>
    pivot_longer(
      cols = {{ pattern }},
      names_to = names_to,
      values_to = values_to
    )


  intensity_vs_design_matrix <- intensity_long_table |>
    mutate({{ names_to }} := purrr::map_chr(!!rlang::sym(names_to), as.character)) |>
    left_join(experimental_design_table,
      by = join_by({{ names_to }})
    )


  list_of_columns_to_pivot <- setdiff(
    colnames(experimental_design_table),
    c(
      protein_id,
      names_to,
      values_to,
      as_string(as_name(enquo(is_missing)))
    )
  )

  intensity_vs_design_matrix_cln <- intensity_vs_design_matrix |>
    mutate({{ is_missing_column }} := case_when(
      is.nan(!!sym(values_to)) |
        is.na(Avg.Log2.Protein.Imputed) ~ TRUE,
      TRUE ~ FALSE
    )) |>
    relocate({{ is_missing_column }}, .after = !!sym(values_to)) |>
    pivot_longer(
      cols = all_of(list_of_columns_to_pivot),
      names_to = "parameter_name",
      values_to = "values"
    )

  missing_value_per_category <- intensity_vs_design_matrix_cln |>
    group_by(
      !!sym(protein_id),
      parameter_name,
      values
    ) |>
    summarise(
      num_values = n(),
      num_missing = sum(is_missing)
    ) |>
    ungroup() |>
    mutate(perc_missing = num_missing / num_values * 100) |>
    mutate(num_present = num_values - num_missing) |>
    mutate(perc_present = 100 - perc_missing) |>
    dplyr::select(
      uniprot_acc,
      parameter_name,
      values,
      num_missing,
      num_present,
      num_values,
      perc_missing,
      perc_present
    ) |>
    mutate(compare_column = paste0(parameter_name, as.character(values)))

  missing_value_per_category
}

# ----------------------------------------------------------------------------
# calculateMissingValuesPerProteinFishersTest
# ----------------------------------------------------------------------------
#' @export
calculateMissingValuesPerProteinFishersTest <- function(contrasts_table, missing_value_per_category) {
  contrasts_table_separated <- contrasts_table |>
    separate(col = contrasts, sep = "[=-]", into = c("contrast_name", "left", "right"))

  runFisherTest <- function(a1, b1, a2, b2) {
    fisher.test(matrix(c(a1, b1, a2, b2), 2, 2, byrow = TRUE))$p.value
  }

  plan(multisession, workers = 8)


  contasts_missing_counts_tbl <- contrasts_table_separated |>
    left_join(missing_value_per_category,
      by = join_by(left == compare_column)
    ) |>
    left_join(missing_value_per_category,
      by = join_by(
        right == compare_column,
        uniprot_acc == uniprot_acc
      ),
      suffix = c(".left", ".right")
    ) |>
    dplyr::filter(!(is.na(num_missing.left) &
      is.na(num_present.left) &
      is.na(num_missing.right) &
      is.na(num_present.right))) |>
    mutate(fisher_test = furrr::future_pmap_dbl(
      list(
        a1 = num_missing.left,
        b1 = num_present.left,
        a2 = num_missing.right,
        b2 = num_present.right
      ),
      \(a1, a2, b1, b2){
        runFisherTest(a1 = a1, b1 = b1, a2 = a2, b2 = b2)
      }
    ))

  # fisher.test(matrix( c(8, 19, 27, 73), 2,2, byrow=TRUE))
  # fisher.test(matrix( c(8, 19, 27, 73), 2,2, byrow=FALSE))

  contasts_missing_fdr_tbl <- contasts_missing_counts_tbl |>
    nest(.by = contrast_name, .key = "tables") |>
    dplyr::mutate(updated_tables = purrr::map(
      tables,
      \(x){
        x |> bind_cols(data.frame(fdr = p.adjust(x$fisher_test, method = "fdr")))
      }
    )) |>
    dplyr::select(-tables) |>
    unnest(updated_tables)

  contasts_missing_fdr_tbl
}

# ----------------------------------------------------------------------------
# getRowsToKeepList
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' For each experimental group, identify proteins that does have enough number of samples with abundance values.
#' @param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#' @param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#' @param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#' @param sample_id The name of the column in design_matrix table that has the sample ID.
#' @param row_id A unique ID for each row of the 'input_table' variable.
#' @param grouping_variable The name of the column in design_matrix table that has the experimental group.
#' @param min_num_samples_per_group An integer representing the minimum number of samples per group.
#' @param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#' @return A list, the name of each element is the sample ID and each element is a vector containing the protein accessions (e.g. row_id) with enough number of values.
#' @export
getRowsToKeepList <- function(input_table, cols, design_matrix, sample_id, row_id, grouping_variable, min_num_samples_per_group, abundance_threshold) {
  abundance_long <- input_table |>
    pivot_longer(
      cols = {{ cols }},
      names_to = as_string(as_name(enquo(sample_id))),
      values_to = "Abundance"
    ) |>
    left_join(design_matrix, by = as_string(as_name(enquo(sample_id))))


  count_values_per_group <- abundance_long |>
    mutate(has_value = ifelse(!is.na(Abundance) & Abundance > abundance_threshold, 1, 0)) |>
    group_by({{ row_id }}, {{ grouping_variable }}) |>
    summarise(num_values = sum(has_value)) |>
    ungroup()


  kept_rows_temp <- count_values_per_group |>
    dplyr::filter(num_values >= min_num_samples_per_group) |>
    dplyr::select(-num_values) |>
    group_by({{ grouping_variable }}) |>
    nest(data = c({{ row_id }})) |>
    ungroup() |>
    mutate(data = purrr::map(data, \(x){
      x[, as_name(enquo(row_id))][[1]]
    }))


  sample_rows_lists <- kept_rows_temp$data
  names(sample_rows_lists) <- kept_rows_temp[, as_name(enquo(grouping_variable))][[1]]

  return(sample_rows_lists)
}

# ----------------------------------------------------------------------------
# averageValuesFromReplicates
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Average values from replicates
#' @param design_matrix Contains the sample_id column and the average_replicates_id column
#' @export
averageValuesFromReplicates <- function(input_table, design_matrix, group_pattern, row_id, sample_id, average_replicates_id) {
  output_table <- input_table |>
    as.data.frame() |>
    rownames_to_column(row_id) |>
    pivot_longer(
      cols = matches(group_pattern),
      names_to = sample_id,
      values_to = "value"
    ) |>
    left_join(design_matrix, by = sample_id) |>
    group_by(!!rlang::sym(average_replicates_id), !!rlang::sym(row_id)) |>
    summarise(value = mean(value, na.rm = TRUE)) |>
    ungroup() |>
    pivot_wider(
      names_from = !!rlang::sym(average_replicates_id),
      values_from = "value"
    ) |>
    column_to_rownames(row_id) |>
    as.matrix()

  return(output_table)
}

# ----------------------------------------------------------------------------
# proteinTechRepCorrelationHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Calculate protein technical replicate correlation
#' @param design_matrix_tech_rep: design matrix with the technical replicates
#' @param data_matrix: input data matrix
#' @param sample_id_column: column name of the sample ID. This is the unique identifier for each sample.
#' @param tech_rep_column: column name of the technical replicates. Technical replicates of the same sample will have the same value.
#' @param tech_rep_num_column: column name of the technical replicate number. This is a unique number for each technical replicate for each sample.
#' @export
proteinTechRepCorrelationHelper <- function(
  design_matrix_tech_rep, data_matrix,
  protein_id_column = "Protein.Ids",
  sample_id_column = "Sample_ID", tech_rep_column = "replicates", tech_rep_num_column = "tech_rep_num", tech_rep_remove_regex = "pool"
) {
  tech_reps_list <- design_matrix_tech_rep |>
    dplyr::pull(!!sym(tech_rep_num_column)) |>
    unique()

  frozen_protein_matrix_tech_rep <- data_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(
      cols = !matches(protein_id_column),
      values_to = "log2_intensity",
      names_to = sample_id_column
    ) |>
    left_join(design_matrix_tech_rep,
      by = join_by(!!sym(sample_id_column) == !!sym(sample_id_column))
    ) |>
    dplyr::filter(!str_detect(!!sym(tech_rep_column), tech_rep_remove_regex)) |>
    dplyr::select(!!sym(protein_id_column), !!sym(tech_rep_column), log2_intensity, !!sym(tech_rep_num_column)) |>
    dplyr::filter(!!sym(tech_rep_num_column) %in% tech_reps_list) |>
    pivot_wider(
      id_cols = c(!!sym(protein_id_column), !!sym(tech_rep_column)),
      names_from = !!sym(tech_rep_num_column),
      values_from = log2_intensity
    ) |>
    nest(data = !matches(protein_id_column)) |>
    mutate(data = purrr::map(data, \(x){
      x |> column_to_rownames(tech_rep_column)
    })) |>
    mutate(pearson = purrr::map_dbl(data, \(x){
      if (length(which(!is.na(x[, 1]))) > 0 & length(which(!is.na(x[, 2]))) > 0) {
        cor(x, use = "pairwise.complete.obs")[1, 2]
      } else {
        NA_real_
      }
    })) |>
    mutate(spearman = purrr::map_dbl(data, \(x){
      if (length(which(!is.na(x[, 1]))) > 0 & length(which(!is.na(x[, 2]))) > 0) {
        cor(x, use = "pairwise.complete.obs", method = "spearman")[1, 2]
      } else {
        NA_real_
      }
    }))

  frozen_protein_matrix_tech_rep
}

