library(testthat)
library(dplyr)

test_that("peptideIntensityFilteringHelper filters by group correctly", {
  # Mock input table (WIDE format)
  input_table <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("AAA", "BBB", "CCC", "DDD"),
    S1 = c(20, 20, 20, 10), # Group G1
    S2 = c(20, 20, 20, 10), # Group G1
    S3 = c(20, 10, 10, 10), # Group G2
    S4 = c(20, 10, 10, 10), # Group G2
    stringsAsFactors = FALSE
  )
  
  # Mock design matrix
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    Group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  
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
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    peptide_quantity_column = Peptide.Normalised,
    core_utilisation = NA
  )
  
  expect_equal(nrow(result), 1)
  expect_equal(result$Stripped.Sequence, "AAA")
  
  # Case 2: max_groups_percentage_cutoff = 60 (allow 1 group to fail)
  result2 <- peptideIntensityFilteringHelper(
    input_table = input_table,
    design_matrix = design_matrix,
    min_peptide_intensity_threshold = 15,
    sample_id_column = "Run",
    grouping_variable = "Group",
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 60,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    peptide_quantity_column = Peptide.Normalised,
    core_utilisation = NA
  )
  
  expect_equal(nrow(result2), 3)
  expect_false("DDD" %in% result2$Stripped.Sequence)
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
  
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    Group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  
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

# <!-- APAF Bioinformatics | test-prot-02-qc-peptide-groupaware.R | Approved | 2026-03-14 -->
