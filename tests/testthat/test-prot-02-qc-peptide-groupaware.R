library(testthat)
library(dplyr)

# fidelity-coverage-compare: shared

build_peptide_filter_design <- function() {
  data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    Group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
}

build_peptide_filter_wide <- function() {
  data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("AAA", "BBB", "CCC", "DDD"),
    S1 = c(20, 20, 20, 10),
    S2 = c(20, 20, 20, 10),
    S3 = c(20, 10, 10, 10),
    S4 = c(20, 10, 10, 10),
    stringsAsFactors = FALSE
  )
}

build_peptide_filter_long <- function(quantity_column = "Peptide.Normalised") {
  values <- c(20, 20, 20, 20, 20, 20, 10, NaN, 20, 20, 10, 10, 10, 10, 10, 10)
  long_input <- tidyr::expand_grid(
    peptide_idx = seq_len(4),
    Run = c("S1", "S2", "S3", "S4")
  ) |>
    dplyr::mutate(
      Protein.Ids = c("P1", "P1", "P2", "P2")[peptide_idx],
      Stripped.Sequence = c("AAA", "BBB", "CCC", "DDD")[peptide_idx],
      "{quantity_column}" := values
    ) |>
    dplyr::select(Protein.Ids, Stripped.Sequence, Run, dplyr::all_of(quantity_column)) |>
    as.data.frame()

  long_input
}

test_that("peptideIntensityFilteringHelper filters by group correctly", {
  # Mock input table (WIDE format)
  input_table <- build_peptide_filter_wide()
  
  # Mock design matrix
  design_matrix <- build_peptide_filter_design()
  
  # Case 1: Threshold = 15
  # P1-AAA: G1 (20,20 - 0% below), G2 (20,20 - 0% below) -> Keep
  # P1-BBB: G1 (20,20 - 0% below), G2 (10,10 - 100% below) -> Failed in 1 group (50%). 
  #         If max_groups_percentage_cutoff = 40, it should be removed.
  # P2-CCC: G1 (20,20 - 0% below), G2 (10,10 - 100% below) -> Same as BBB.
  # P2-DDD: G1 (10,10 - 100% below), G2 (10,10 - 100% below) -> Failed in 2 groups (100%).
  
  result <- peptideIntensityFilteringHelper(
    input_table = input_table,
    design_matrix = design_matrix,
    min_peptide_intensity_threshold = 15,
    sample_id_column = "Run",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 50, # More than 50% below threshold means group fails
    max_groups_percentage_cutoff = 40, # More than 40% of groups failing means peptide removed
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_quantity_column = "Peptide.Normalised",
    core_utilisation = NA
  )
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("Protein.Ids", "Stripped.Sequence") %in% names(result)))
  expect_lte(nrow(result), nrow(input_table))
  
  # Case 2: max_groups_percentage_cutoff = 60 (allow 1 group to fail)
  result2 <- peptideIntensityFilteringHelper(
    input_table = input_table,
    design_matrix = design_matrix,
    min_peptide_intensity_threshold = 15,
    sample_id_column = "Run",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 60,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_quantity_column = "Peptide.Normalised",
    core_utilisation = NA
  )
  
  expect_s3_class(result2, "data.frame")
  expect_true(all(c("Protein.Ids", "Stripped.Sequence") %in% names(result2)))
  expect_lte(nrow(result2), nrow(input_table))
})

test_that("peptideIntensityFilteringHelper preserves aliases, long input, and validation behavior", {
  design_matrix <- build_peptide_filter_design()
  long_input <- build_peptide_filter_long("Intensity")

  alias_result <- peptideIntensityFilteringHelper(
    input_data = long_input,
    design_matrix = design_matrix,
    min_peptide_intensity_threshold = 15,
    sample_id_column = "run",
    grouping_variable = "group",
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 60,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    raw_quantity_column = "Intensity",
    core_utilisation = NA
  )

  expect_s3_class(alias_result, "data.frame")
  expect_true("AAA" %in% alias_result$Stripped.Sequence)
  expect_lte(nrow(alias_result), nrow(long_input))

  no_removal_result <- peptideIntensityFilteringHelper(
    input_table = build_peptide_filter_wide(),
    design_matrix = design_matrix,
    min_peptide_intensity_threshold = 0,
    sample_id_column = "Run",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 100,
    max_groups_percentage_cutoff = 100,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    peptide_quantity_column = "Peptide.Normalised",
    core_utilisation = NA
  )
  expect_equal(nrow(no_removal_result), nrow(build_peptide_filter_wide()))

  expect_error(
    peptideIntensityFilteringHelper(
      input_table = NULL,
      design_matrix = design_matrix,
      protein_id_column = "Protein.Ids",
      peptide_sequence_column = "Stripped.Sequence"
    ),
    "input_table or input_data"
  )

  expect_error(
    peptideIntensityFilteringHelper(
      input_table = build_peptide_filter_wide(),
      design_matrix = design_matrix,
      sample_id_column = "Run",
      grouping_variable = "MissingGroup",
      protein_id_column = "Protein.Ids",
      peptide_sequence_column = "Stripped.Sequence",
      peptide_quantity_column = "Peptide.Normalised"
    ),
    "not found"
  )
})

test_that("removePeptidesWithMissingValuesPercentHelper filters by group correctly", {
  # Mock input table (WIDE format)
  input_table <- data.frame(
    Protein.Ids = c("P1", "P1", "P2"),
    Stripped.Sequence = c("AAA", "BBB", "CCC"),
    S1 = c(20, 20, 20), # Group G1
    S2 = c(20, 20, NA), # Group G1
    S3 = c(20, NA, NA), # Group G2
    S4 = c(20, NA, NA), # Group G2
    stringsAsFactors = FALSE
  )
  
  design_matrix <- build_peptide_filter_design()
  
  # Case 1: Threshold = 0 (anything > 0 is NOT missing)
  # AAA: G1 (20,20 - 0% missing), G2 (20,20 - 0% missing) -> Keep
  # BBB: G1 (20,20 - 0% missing), G2 (NA,NA - 100% missing) -> Failed 1 group (50%)
  # CCC: G1 (20,NA - 50% missing), G2 (NA,NA - 100% missing) -> Failed 2 groups if cutoff < 50
  
  result <- removePeptidesWithMissingValuesPercentHelper(
    input_table = input_table,
    design_matrix = design_matrix,
    sample_id = "Run",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 25, # More than 25% missing means group fails
    max_groups_percentage_cutoff = 40,
    abundance_threshold = 0,
    abundance_column = "Abundance"
  )
  
  expect_equal(nrow(result), 1)
  expect_equal(result$Stripped.Sequence, "AAA")
})

test_that("removePeptidesWithMissingValuesPercentHelper preserves long input and no-removal branches", {
  design_matrix <- build_peptide_filter_design()
  long_input <- build_peptide_filter_long("Abundance")

  long_result <- removePeptidesWithMissingValuesPercentHelper(
    input_table = long_input,
    design_matrix = design_matrix,
    sample_id = "run",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    grouping_variable = "group",
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 60,
    abundance_threshold = 15,
    abundance_column = "Abundance"
  )
  expect_s3_class(long_result, "data.frame")
  expect_true("AAA" %in% long_result$Stripped.Sequence)
  expect_lte(nrow(long_result), nrow(long_input))

  no_removal_result <- removePeptidesWithMissingValuesPercentHelper(
    input_table = build_peptide_filter_wide(),
    design_matrix = design_matrix,
    sample_id = "Run",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 100,
    max_groups_percentage_cutoff = 100,
    abundance_threshold = 0,
    abundance_column = "Abundance"
  )
  expect_equal(nrow(no_removal_result), nrow(build_peptide_filter_wide()))

  expect_error(
    removePeptidesWithMissingValuesPercentHelper(
      input_table = build_peptide_filter_wide(),
      design_matrix = design_matrix,
      sample_id = "Run",
      protein_id_column = "Protein.Ids",
      peptide_sequence_column = "Stripped.Sequence",
      grouping_variable = "MissingGroup",
      abundance_column = "Abundance"
    ),
    "not found"
  )
})

# <!-- APAF Bioinformatics | test-prot-02-qc-peptide-groupaware.R | Approved | 2026-03-14 -->
