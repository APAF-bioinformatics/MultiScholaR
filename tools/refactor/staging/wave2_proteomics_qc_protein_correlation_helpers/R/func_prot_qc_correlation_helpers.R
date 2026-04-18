# ----------------------------------------------------------------------------
# getPairsOfSamplesTable
# ----------------------------------------------------------------------------
#' @title Get Pairs of Samples for Comparison
#' @param input_table A table with two columns, the Run ID column and the technical replicate group column
#' @param run_id_column A string representing the name of the column with the Run ID (or sample ID)
#' @param replicate_group_column A string representing the name of the column with the technical replicate group
#' @return A table with three columns, the technical replicate group column, run ID X column, and run ID Y column
#' @export
getPairsOfSamplesTable <- function(
  input_table,
  run_id_column,
  replicate_group_column
) {
  pairs_for_comparison <- input_table |>
    inner_join(input_table, by = join_by(!!rlang::sym(replicate_group_column))) |>
    dplyr::filter(!!rlang::sym(paste0(run_id_column, ".x")) > !!rlang::sym(paste0(run_id_column, ".y"))) |>
    arrange(!!rlang::sym(replicate_group_column)) |>
    relocate(!!rlang::sym(replicate_group_column),
      .before = paste0(run_id_column, ".x")
    )

  pairs_for_comparison
}

# ----------------------------------------------------------------------------
# calulatePearsonCorrelation
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation
#' @description Calculate the Pearson correlation of the abundances of peptides between two samples X and Y.
#' @param ms_filename_x A string representing the sample file name X (for a pair of sample in the same technical replicate group) for correlation score calculation.
#' @param ms_filename_y A string representing the sample file name Y (for a pair of sample in the same technical replicate group) for correlation score calculation.
#' @param input_table A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#' @param sample_id_column Sample ID column name as string (default = "Run").
#' @param protein_id_column Protein accession column name as string (default = "Protein.Ids").
#' @param peptide_sequence_column Peptide sequence column name as string (default = "Stripped.Sequence").
#' @param peptide_normalised_column Normalised peptide abundance column name as string (default = "Peptide.Normalised").
#' @return The pearson correlation value of the abundances of peptides between two samples X and Y.
#' @export
calulatePearsonCorrelation <- function(
  ms_filename_x, ms_filename_y, input_table,
  sample_id_column = "Run",
  protein_id_column = "Protein.Ids",
  peptide_sequence_column = "Stripped.Sequence",
  peptide_normalised_column = "Peptide.Normalised"
) {
  tab_x <- input_table |>
    dplyr::filter(!!rlang::sym(sample_id_column) == ms_filename_x)

  tab_y <- input_table |>
    dplyr::filter(!!rlang::sym(sample_id_column) == ms_filename_y)

  merged_tbl <- tab_x |>
    inner_join(tab_y, by = join_by(!!rlang::sym(protein_id_column), !!rlang::sym(peptide_sequence_column)))

  # merged_tbl |>
  #   dplyr::filter(!is.na( !!sym(paste0(peptide_normalised_column, ".x")) ) & !is.na( !!sym(paste0(peptide_normalised_column, ".x")))) |>
  #   head() |> print()

  # print( paste(ms_filename_x, ms_filename_y))
  input_x <- merged_tbl[[paste0(peptide_normalised_column, ".x")]]
  input_y <- merged_tbl[[paste0(peptide_normalised_column, ".y")]]
  if (length(input_x) > 0 & length(input_y) > 0) {
    cor_result <- cor(input_x,
      input_y,
      use = "pairwise.complete.obs"
    )

    cor_result
  } else {
    return(NA)
  }
}

# ----------------------------------------------------------------------------
# calculatePearsonCorrelationMatrix
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation Matrix (Optimized)
#' @description Computes all pairwise sample correlations using a single matrix operation.
#' Much faster than iterative approach for large datasets. Pivots data to wide format
#' (proteins x samples) and computes the full correlation matrix in one vectorized call.
#' @param input_table Long-format data table with sample, protein, and value columns
#' @param sample_id_column String name of sample ID column
#' @param protein_id_column String name of protein ID column
#' @param peptide_normalised_column String name of normalized value column
#' @return Correlation matrix (samples x samples)
#' @importFrom tidyr pivot_wider
#' @export
calculatePearsonCorrelationMatrix <- function(input_table,
                                              sample_id_column,
                                              protein_id_column,
                                              peptide_normalised_column) {
  message("*** PEARSON MATRIX: Pivoting data to wide format... ***")
  pivot_start <- Sys.time()


  # Pivot to wide: rows = proteins, columns = samples
  wide_matrix <- input_table |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(protein_id_column),
      names_from = dplyr::all_of(sample_id_column),
      values_from = dplyr::all_of(peptide_normalised_column),
      values_fn = mean # Handle any duplicates
    ) |>
    dplyr::select(-dplyr::all_of(protein_id_column)) |>
    as.matrix()

  pivot_elapsed <- as.numeric(difftime(Sys.time(), pivot_start, units = "secs"))
  message(sprintf(
    "*** PEARSON MATRIX: Pivot completed in %.2f seconds (%d proteins x %d samples) ***",
    pivot_elapsed, nrow(wide_matrix), ncol(wide_matrix)
  ))

  message("*** PEARSON MATRIX: Computing correlation matrix... ***")
  cor_start <- Sys.time()

  # Single cor() call computes ALL pairwise correlations
  cor_matrix <- cor(wide_matrix, use = "pairwise.complete.obs")

  cor_elapsed <- as.numeric(difftime(Sys.time(), cor_start, units = "secs"))
  message(sprintf("*** PEARSON MATRIX: Correlation matrix computed in %.2f seconds ***", cor_elapsed))

  cor_matrix
}

# ----------------------------------------------------------------------------
# calculatePearsonCorrelationOptimized
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation (Optimized for pre-filtered data)
#' @description Fast correlation calculation for pre-filtered sample data.
#' This function is optimized to work with pre-partitioned data and avoids
#' redundant filtering operations.
#' @param data_x Pre-filtered data for sample X
#' @param data_y Pre-filtered data for sample Y
#' @param protein_id_column String name of protein ID column
#' @param peptide_sequence_column String name of peptide sequence column
#' @param peptide_normalised_column String name of normalized peptide column
#' @return Numeric Pearson correlation value, or NA_real_ if insufficient data
#' @export
calculatePearsonCorrelationOptimized <- function(data_x, data_y,
                                                 protein_id_column,
                                                 peptide_sequence_column,
                                                 peptide_normalised_column) {
  # Quick validation
  if (is.null(data_x) || is.null(data_y) || nrow(data_x) == 0 || nrow(data_y) == 0) {
    return(NA_real_)
  }

  # Inner join on protein + peptide (already filtered, just need matching)
  merged_tbl <- data_x |>
    dplyr::inner_join(
      data_y,
      by = c(protein_id_column, peptide_sequence_column),
      suffix = c(".x", ".y")
    )

  if (nrow(merged_tbl) == 0) {
    return(NA_real_)
  }

  # Extract values
  values_x <- merged_tbl[[paste0(peptide_normalised_column, ".x")]]
  values_y <- merged_tbl[[paste0(peptide_normalised_column, ".y")]]

  # Calculate correlation
  if (length(values_x) > 0 && length(values_y) > 0) {
    cor(values_x, values_y, use = "pairwise.complete.obs")
  } else {
    NA_real_
  }
}

# ----------------------------------------------------------------------------
# calulatePearsonCorrelationForSamplePairsHelper
# ----------------------------------------------------------------------------
#' @title Calculate Pearson Correlation for Sample Pairs
#' @description Helper function to calculate pairwise Pearson correlations between samples
#' using optimized sequential processing with pre-partitioned data for fast performance.
#' Data is pre-partitioned by sample ID to eliminate redundant filtering operations.
#' @param samples_id_tbl Table of sample IDs
#' @param run_id_column String name of run ID column (default: "ms_filename")
#' @param replicate_group_column String name of replicate group column (default: "general_sample_info")
#' @param input_table Input data table
#' @param num_of_cores Number of CPU cores (parameter retained for compatibility but not used)
#' @param sample_id_column String name of sample ID column (default: "Run")
#' @param protein_id_column String name of protein ID column (default: "Protein.Ids")
#' @param peptide_sequence_column String name of peptide sequence column (default: "Stripped.Sequence")
#' @param peptide_normalised_column String name of normalized peptide column (default: "Peptide.Normalised")
#' @return Data frame with pairwise correlations
#' @importFrom purrr map2_dbl map_chr
#' @export
calulatePearsonCorrelationForSamplePairsHelper <- function(
  samples_id_tbl,
  run_id_column = "ms_filename",
  replicate_group_column = "general_sample_info",
  input_table,
  num_of_cores = 1,
  sample_id_column = "Run",
  protein_id_column = "Protein.Ids",
  peptide_sequence_column = "Stripped.Sequence",
  peptide_normalised_column = "Peptide.Normalised"
) {
  pairs_for_comparison <- getPairsOfSamplesTable(samples_id_tbl,
    run_id_column = run_id_column,
    replicate_group_column = replicate_group_column
  )

  # Log pair generation
  num_pairs <- nrow(pairs_for_comparison)
  message(sprintf("*** PEARSON HELPER: Generated %d sample pairs for correlation analysis ***", num_pairs))

  message(sprintf("--- DEBUG66 [calulatePearsonCorrelationForSamplePairsHelper]: Processing %d pairs ---", num_pairs))

  # Track total calculation time
  total_start_time <- Sys.time()

  # MATRIX-BASED OPTIMIZATION: Compute all correlations in one vectorized operation
  # This replaces the slow group_split + iterative approach
  message("*** PEARSON HELPER: Using matrix-based correlation (fast vectorized approach) ***")

  # Compute full correlation matrix in one call
  message(sprintf("   [calulatePearsonCorrelationForSamplePairsHelper] Step: Computing correlation matrix..."))
  mat_start <- Sys.time()
  cor_matrix <- calculatePearsonCorrelationMatrix(
    input_table = input_table,
    sample_id_column = sample_id_column,
    protein_id_column = protein_id_column,
    peptide_normalised_column = peptide_normalised_column
  )
  mat_end <- Sys.time()
  message(sprintf(
    "   [calulatePearsonCorrelationForSamplePairsHelper] Matrix computed. Dim: %d x %d. Duration: %.2f secs",
    nrow(cor_matrix), ncol(cor_matrix), as.numeric(difftime(mat_end, mat_start, units = "secs"))
  ))

  # Extract correlations for the specific pairs we need
  message("*** PEARSON HELPER: Extracting correlations for specified pairs... ***")
  extract_start <- Sys.time()

  sample_pairs_x <- pairs_for_comparison[[paste0(run_id_column, ".x")]]
  sample_pairs_y <- pairs_for_comparison[[paste0(run_id_column, ".y")]]

  message(sprintf("   [calulatePearsonCorrelationForSamplePairsHelper] Extracting %d values from matrix...", length(sample_pairs_x)))

  # Vectorized extraction from correlation matrix
  correlations <- purrr::map2_dbl(
    sample_pairs_x,
    sample_pairs_y,
    \(x, y) {
      # Direct lookup from correlation matrix - O(1) per pair
      if (x %in% colnames(cor_matrix) && y %in% colnames(cor_matrix)) {
        cor_matrix[x, y]
      } else {
        NA_real_
      }
    }
  )

  extract_elapsed <- as.numeric(difftime(Sys.time(), extract_start, units = "secs"))
  message(sprintf("*** PEARSON HELPER: Pair extraction completed in %.2f seconds ***", extract_elapsed))

  # Add results back to dataframe
  pearson_correlation_per_pair <- pairs_for_comparison |>
    mutate(pearson_correlation = correlations)

  # Log completion with timing statistics
  total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  message(sprintf(
    "*** PEARSON HELPER: Completed %d correlations in %.1f seconds (matrix approach) ***",
    num_pairs, total_elapsed
  ))

  pearson_correlation_per_pair
}

# ----------------------------------------------------------------------------
# filterSamplesByPeptideCorrelationThreshold
# ----------------------------------------------------------------------------
#' @title Filter Samples by Peptide Correlation Threshold
#' @description Remove samples which is correlated with any technical replicate samples
#' @param pearson_correlation_per_pair A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#' @param peptide_keep_samples_with_min_num_peptides A data frame with the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#' @param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#' @param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#' @param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#' @param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#' @param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#' @return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#' @export
filterSamplesByPeptideCorrelationThreshold <- function(
  pearson_correlation_per_pair,
  peptide_keep_samples_with_min_num_peptides,
  min_pearson_correlation_threshold = 0.75,
  filename_column_x = ms_filename.x,
  filename_column_y = ms_filename.y,
  correlation_column = pearson_correlation,
  filename_id_column = "Run"
) {
  # Samples to keep include all those pairs of samples with correlation score passing threshold
  samples_to_keep <- pearson_correlation_per_pair |>
    dplyr::filter({{ correlation_column }} >= min_pearson_correlation_threshold) |>
    pivot_longer(
      cols = c({{ filename_column_x }}, {{ filename_column_y }}),
      values_to = filename_id_column
    ) |>
    dplyr::distinct(!!rlang::sym(filename_id_column))

  samples_above_correlation_theshold <- peptide_keep_samples_with_min_num_peptides |>
    inner_join(samples_to_keep,
      by = join_by(!!rlang::sym(filename_id_column) == !!rlang::sym(filename_id_column))
    ) |>
    distinct()

  samples_above_correlation_theshold
}

# ----------------------------------------------------------------------------
# findSamplesPairBelowPeptideCorrelationThreshold
# ----------------------------------------------------------------------------
#' @export
findSamplesPairBelowPeptideCorrelationThreshold <- function(
  pearson_correlation_per_pair,
  peptide_keep_samples_with_min_num_peptides,
  min_pearson_correlation_threshold = 0.75,
  filename_column_x = ms_filename.x,
  filename_column_y = ms_filename.y,
  correlation_column = pearson_correlation,
  filename_id_column = "Run"
) {
  samples_to_keep <- pearson_correlation_per_pair |>
    dplyr::filter({{ correlation_column }} >= min_pearson_correlation_threshold) |>
    pivot_longer(
      cols = c({{ filename_column_x }}, {{ filename_column_y }}),
      values_to = filename_id_column
    ) |>
    dplyr::distinct(!!rlang::sym(filename_id_column))

  samples_below_correlation_theshold <- pearson_correlation_per_pair |>
    pivot_longer(
      cols = c({{ filename_column_x }}, {{ filename_column_y }}),
      values_to = filename_id_column
    ) |>
    dplyr::distinct(!!rlang::sym(filename_id_column)) |>
    dplyr::inner_join(samples_to_keep,
      by = join_by(!!rlang::sym(filename_id_column) == !!rlang::sym(filename_id_column))
    )

  samples_below_correlation_theshold
}

# ----------------------------------------------------------------------------
# filterSamplesByProteinCorrelationThresholdHelper
# ----------------------------------------------------------------------------
#' @title Filter Samples by Protein Correlation
#' @description Remove samples which is correlated with any technical replicate samples
#' @param protein_intensity_table A data frame with the following columns: 1. ID of technical replicate group, 2. sample file name X, 3. sample file name Y, 4. Pearson correlation of the abundances of peptides between sample X and Y.
#' @param peptide_keep_samples_with_min_num_peptides A data frame with the proteins as rows and samples ID as columns.
#' @param min_pearson_correlation_threshold Minimum pearson correlation for a pair of files to be considered to be consistent and kept for further analysis
#' @param filename_column_x Name of column containing the sample file name X (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#' @param filename_column_y Name of column containing the sample file name Y (for a pair of sample in the same technical replicate group). Tidyverse column header format, not a string.
#' @param protein_id_column Name of column containing the protein ID. Tidyverse column header format, not a string.
#' @param correlation_column Name of column containing the Pearson's correlation score between Sample X and Y. Tidyverse column header format, not a string.
#' @param filename_id_column A string indicating the name of column that contains the sample ID or Run ID in the data frame `peptide_keep_samples_with_min_num_peptides`.
#' @return A table without samples that are poorly correlated with the rest of the samples in the technical replicate group. Contains the following columns: 1. Sample file name or Run name, 2. Protein IDs, 3. Stripped peptide sequence, 4. Normalised peptide abundances
#' @export
filterSamplesByProteinCorrelationThresholdHelper <- function(
  pearson_correlation_per_pair,
  protein_intensity_table,
  min_pearson_correlation_threshold = 0.75,
  filename_column_x = ms_filename.x,
  filename_column_y = ms_filename.y,
  protein_id_column = Protein.Ids,
  correlation_column = pearson_correlation
) {
  message("--- DEBUG66 [filterSamplesByProteinCorrelationThresholdHelper]: Entry ---")
  message(sprintf(
    "   [filterSamplesByProteinCorrelationThresholdHelper] protein_intensity_table dims: %d x %d",
    nrow(protein_intensity_table), ncol(protein_intensity_table)
  ))
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Threshold: %.2f", min_pearson_correlation_threshold))

  # All Samples
  all_samples <- pearson_correlation_per_pair |>
    pivot_longer(
      cols = c({{ filename_column_x }}, {{ filename_column_y }}),
      values_to = "temp_column"
    ) |>
    dplyr::distinct(temp_column)

  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Total unique samples in pairs: %d", nrow(all_samples)))

  # Samples to keep include all those pairs of samples with correlation score passing threshold
  samples_to_keep <- pearson_correlation_per_pair |>
    dplyr::filter({{ correlation_column }} >= min_pearson_correlation_threshold) |>
    pivot_longer(
      cols = c({{ filename_column_x }}, {{ filename_column_y }}),
      values_to = "temp_column"
    ) |>
    dplyr::distinct(temp_column)

  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Samples passing correlation threshold: %d", nrow(samples_to_keep)))

  # Samples to keep anyway
  # Use robust method to get protein ID column name as string to prevent bad_alloc on setdiff
  pid_col_name <- tryCatch(
    {
      if (is.character(protein_id_column)) protein_id_column else rlang::as_string(rlang::ensym(protein_id_column))
    },
    error = function(e) rlang::as_string(rlang::ensym(protein_id_column))
  )

  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Protein ID Column Name resolved to: '%s'", pid_col_name))

  # DEBUG: Trace setdiff inputs
  all_cols <- colnames(protein_intensity_table)
  pair_samples <- all_samples |> dplyr::pull(temp_column)

  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Calculating 'samples_to_keep_anyway'..."))
  message(sprintf("      - Total columns: %d", length(all_cols)))
  message(sprintf("      - Pair samples: %d", length(pair_samples)))

  cols_not_in_pairs <- setdiff(all_cols, pair_samples)
  message(sprintf("      - Cols NOT in pairs: %d", length(cols_not_in_pairs)))

  samples_to_keep_anyway <- setdiff(cols_not_in_pairs, pid_col_name)
  message(sprintf("      - Final 'samples_to_keep_anyway' (excluding ID col): %d", length(samples_to_keep_anyway)))

  print(samples_to_keep_anyway)

  # Samples in the table to keep
  samples_to_keep_subset <- colnames(protein_intensity_table)[colnames(protein_intensity_table) %in% (samples_to_keep |> dplyr::pull(temp_column))]
  message(sprintf("   [filterSamplesByProteinCorrelationThresholdHelper] Samples from pairs to keep: %d", length(samples_to_keep_subset)))

  samples_above_correlation_threshold <- protein_intensity_table |>
    dplyr::select({{ protein_id_column }}, all_of(c(samples_to_keep_anyway, samples_to_keep_subset)))

  message(sprintf(
    "   [filterSamplesByProteinCorrelationThresholdHelper] Result dims: %d x %d",
    nrow(samples_above_correlation_threshold), ncol(samples_above_correlation_threshold)
  ))

  samples_above_correlation_threshold
}

