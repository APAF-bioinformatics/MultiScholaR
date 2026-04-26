# fidelity-coverage-compare: shared
library(testthat)

getPairsOfSamplesTable <- get(
  "getPairsOfSamplesTable",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calulatePearsonCorrelation <- get(
  "calulatePearsonCorrelation",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calculatePearsonCorrelationMatrix <- get(
  "calculatePearsonCorrelationMatrix",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calculatePearsonCorrelationOptimized <- get(
  "calculatePearsonCorrelationOptimized",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
calulatePearsonCorrelationForSamplePairsHelper <- get(
  "calulatePearsonCorrelationForSamplePairsHelper",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
filterSamplesByPeptideCorrelationThreshold <- get(
  "filterSamplesByPeptideCorrelationThreshold",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
findSamplesPairBelowPeptideCorrelationThreshold <- get(
  "findSamplesPairBelowPeptideCorrelationThreshold",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
filterSamplesByProteinCorrelationThresholdHelper <- get(
  "filterSamplesByProteinCorrelationThresholdHelper",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("protein QC correlation helpers preserve pair generation, matrix correlation, and filtering behavior", {
  samples_id_tbl <- data.frame(
    general_sample_info = c("G1", "G1", "G2"),
    ms_filename = c("S1", "S2", "S3"),
    stringsAsFactors = FALSE
  )
  peptide_table <- data.frame(
    Run = c("S1", "S1", "S2", "S2", "S3", "S3"),
    Protein.Ids = c("P1", "P2", "P1", "P2", "P1", "P2"),
    Stripped.Sequence = c("pep1", "pep2", "pep1", "pep2", "pep1", "pep2"),
    Peptide.Normalised = c(10, 20, 11, 19, 30, 60),
    stringsAsFactors = FALSE
  )

  pairs_tbl <- getPairsOfSamplesTable(
    input_table = samples_id_tbl,
    run_id_column = "ms_filename",
    replicate_group_column = "general_sample_info"
  )
  expect_identical(pairs_tbl$general_sample_info, "G1")
  expect_identical(pairs_tbl$ms_filename.x, "S2")
  expect_identical(pairs_tbl$ms_filename.y, "S1")

  pair_cor <- calulatePearsonCorrelation(
    ms_filename_x = "S1",
    ms_filename_y = "S2",
    input_table = peptide_table
  )
  expect_true(pair_cor > 0.9)

  missing_cor <- calulatePearsonCorrelation(
    ms_filename_x = "S1",
    ms_filename_y = "missing",
    input_table = peptide_table
  )
  expect_true(is.na(missing_cor))

  cor_matrix <- calculatePearsonCorrelationMatrix(
    input_table = peptide_table,
    sample_id_column = "Run",
    protein_id_column = "Protein.Ids",
    peptide_normalised_column = "Peptide.Normalised"
  )
  expect_identical(colnames(cor_matrix), c("S1", "S2", "S3"))
  expect_equal(unname(diag(cor_matrix)), c(1, 1, 1))

  optimized_cor <- calculatePearsonCorrelationOptimized(
    data_x = subset(peptide_table, Run == "S1"),
    data_y = subset(peptide_table, Run == "S2"),
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_normalised_column = "Peptide.Normalised"
  )
  expect_equal(optimized_cor, pair_cor)
  expect_true(is.na(calculatePearsonCorrelationOptimized(
    data_x = subset(peptide_table, Run == "S1"),
    data_y = peptide_table[0, ],
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_normalised_column = "Peptide.Normalised"
  )))

  helper_result <- calulatePearsonCorrelationForSamplePairsHelper(
    samples_id_tbl = samples_id_tbl,
    run_id_column = "ms_filename",
    replicate_group_column = "general_sample_info",
    input_table = peptide_table,
    sample_id_column = "Run",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_normalised_column = "Peptide.Normalised"
  )
  expect_true("pearson_correlation" %in% colnames(helper_result))
  expect_equal(helper_result$pearson_correlation[[1]], pair_cor)

  pearson_pairs <- data.frame(
    general_sample_info = c("G1", "G1"),
    ms_filename.x = c("S1", "S1"),
    ms_filename.y = c("S2", "S3"),
    pearson_correlation = c(0.99, 0.2),
    stringsAsFactors = FALSE
  )

  kept_peptides <- filterSamplesByPeptideCorrelationThreshold(
    pearson_correlation_per_pair = pearson_pairs,
    peptide_keep_samples_with_min_num_peptides = peptide_table,
    min_pearson_correlation_threshold = 0.95,
    filename_id_column = "Run"
  )
  expect_true(all(unique(kept_peptides$Run) %in% c("S1", "S2")))
  expect_false("S3" %in% kept_peptides$Run)

  below_pairs <- findSamplesPairBelowPeptideCorrelationThreshold(
    pearson_correlation_per_pair = pearson_pairs,
    peptide_keep_samples_with_min_num_peptides = peptide_table,
    min_pearson_correlation_threshold = 0.95,
    filename_id_column = "Run"
  )
  expect_true(all(below_pairs$Run %in% c("S1", "S2")))

  protein_intensity_table <- data.frame(
    Protein.Ids = c("P1", "P2"),
    S1 = c(10, 20),
    S2 = c(11, 19),
    S3 = c(30, 60),
    S4 = c(5, 6),
    stringsAsFactors = FALSE
  )

  protein_filtered <- filterSamplesByProteinCorrelationThresholdHelper(
    pearson_correlation_per_pair = pearson_pairs,
    protein_intensity_table = protein_intensity_table,
    min_pearson_correlation_threshold = 0.95,
    protein_id_column = "Protein.Ids"
  )
  expect_true(all(c("Protein.Ids", "S1", "S2", "S4") %in% colnames(protein_filtered)))
  expect_false("S3" %in% colnames(protein_filtered))
})
