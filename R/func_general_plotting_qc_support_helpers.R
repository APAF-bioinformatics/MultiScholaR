# ----------------------------------------------------------------------------
# plotDensityOfProteinIntensityPerSample
# ----------------------------------------------------------------------------
#' @export
plotDensityOfProteinIntensityPerSample <- function(
  protein_intensity_long_tbl,
  number_of_peptides_per_protein_per_sample,
  protein_id_column = Protein.Ids,
  sample_id_column = Run,
  num_peptides_column = num_peptides_after_impute,
  protein_intensity_column = Log2.Protein.Imputed
) {
  protein_intensity_vs_num_peptides_for_replicates <- protein_intensity_long_tbl |>
    left_join(number_of_peptides_per_protein_per_sample,
      by = join_by({{ protein_id_column }}, {{ sample_id_column }})
    ) |>
    mutate(peptides_status = ifelse({{ num_peptides_column }} == 1, "Multiple Peptides",
      "Single Peptide"
    )) |>
    ggplot(aes({{ protein_intensity_column }}, group = peptides_status, fill = peptides_status, alpha = 0.5)) +
    geom_density() +
    scale_alpha(guide = "none") +
    apafTheme() +
    xlab("log2 Protein Intensity") +
    ylab("Density") +
    labs(fill = "Peptide") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}

# ----------------------------------------------------------------------------
# plotPercentSamplesVsProteinQuantified
# ----------------------------------------------------------------------------
#' @export
plotPercentSamplesVsProteinQuantified <- function(
  protein_intensity_long_tbl = frozen_protein_table,
  number_of_peptides_per_protein_per_sample = number_of_peptides_per_protein_per_sample,
  protein_id_column = Protein.Ids,
  sample_id_column = Run,
  num_peptides_column = num_peptides_after_impute,
  protein_intensity_column = Log2.Protein.Imputed
) {
  samples_vs_intensity <- protein_intensity_long_tbl |>
    left_join(number_of_peptides_per_protein_per_sample,
      by = join_by({{ protein_id_column }}, {{ sample_id_column }})
    ) |>
    mutate(peptides_status = ifelse(num_peptides_after_impute == 1, "Multiple Peptides",
      "Single Peptide"
    ))

  total_num_samples <- samples_vs_intensity |>
    distinct({{ sample_id_column }}) |>
    nrow()

  summarise_peptide_status <- function(input_vector) {
    if ("Multiple Peptides" %in% input_vector) {
      return("Multiple Peptides")
    } else {
      return("Single Peptide")
    }
  }

  num_samples_per_protein <- samples_vs_intensity |>
    dplyr::filter(!is.na({{ protein_intensity_column }})) |>
    group_by({{ protein_id_column }}) |>
    summarise(
      num_values = n(),
      peptides_status = summarise_peptide_status(peptides_status)
    ) |>
    ungroup() |>
    mutate(percentage = num_values / total_num_samples * 100)

  num_samples_per_protein |>
    mutate(percentage_bin = cut(percentage, breaks = c(0, 20, 40, 60, 80, 100))) |>
    ggplot(aes(percentage_bin, fill = peptides_status, group = peptides_status)) +
    geom_bar(position = "dodge") +
    apafTheme() +
    xlab("Percentage of Samples") +
    ylab("Num. Quantified Proteins ") +
    labs(fill = "Peptide") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}

# ----------------------------------------------------------------------------
# plotNumMissingValues
# ----------------------------------------------------------------------------
#' Plot the number of missing values in each sample
#' @param input_table  Data matrix with each row as a protein and each column a sample.
#' @return A ggplot2 bar plot showing the number of missing values per column.
#' @export
plotNumMissingValues <- function(input_table) {
  plot_num_missing_values <- apply(
    data.matrix(log2(input_table)), 2,
    function(x) {
      length(which(!is.finite(x)))
    }
  ) |>
    t() |>
    t() |>
    set_colnames("No. of Missing Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  # Strip captured environment to prevent memory bloat
  plot_num_missing_values$plot_env <- rlang::base_env()
  plot_num_missing_values
}

# ----------------------------------------------------------------------------
# plotNumOfValues
# ----------------------------------------------------------------------------
#' Plot the number of values in each sample
#' @param input_table  Data matrix with each row as a protein and each column a sample.
#' @return A ggplot2 bar plot showing the number of missing values per column.
#' @export
plotNumOfValues <- function(input_table) {
  plot_num_missing_values <- apply(
    data.matrix(log2(input_table)), 2,
    function(x) {
      length(which(!is.na(x)))
    }
  ) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  # Strip captured environment to prevent memory bloat
  plot_num_missing_values$plot_env <- rlang::base_env()
  plot_num_missing_values
}

# ----------------------------------------------------------------------------
# plotNumOfValuesNoLog
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Plot the number of values in each sample
#' @param input_table  Data matrix with each row as a protein and each column a sample.
#' @return A ggplot2 bar plot showing the number of missing values per column.
#' @export
plotNumOfValuesNoLog <- function(input_table) {
  plot_num_missing_values <- apply(
    data.matrix(input_table), 2,
    function(x) {
      length(which(!is.na(x)))
    }
  ) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  # Strip captured environment to prevent memory bloat
  plot_num_missing_values$plot_env <- rlang::base_env()
  plot_num_missing_values
}

