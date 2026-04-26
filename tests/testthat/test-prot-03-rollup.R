# fidelity-coverage-compare: shared
# testthat for Proteomics Rollup
# Phase 4 of Proteomics GUI Test Strategy

test_that("rolled-up protein snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp03_rolled_up_protein.rds")
  
  if (file.exists(cp_file)) {
    first_line <- tryCatch(readLines(cp_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && identical(first_line[[1]], "version https://git-lfs.github.com/spec/v1")) {
      skip("Snapshot cp03 is a Git LFS pointer and the binary artifact is not present")
    }

    obj <- readRDS(cp_file)
    expect_s4_class(obj, "ProteinQuantitativeData")
    expect_true(nrow(obj@protein_quant_table) > 0)
  } else {
    skip("Snapshot cp03 not found")
  }
})

test_that("calcPeptidesPerProtein works with mock data", {
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1", "P2"),
    Stripped.Sequence = c("PEP1", "PEP2", "PEP3"),
    Run = rep("S1", 3),
    Q.Value = rep(0.01, 3),
    Precursor.Quantity = rep(100, 3),
    stringsAsFactors = FALSE
  )
  
  # Create mock S4 with peptide_data and design_matrix
  obj <- new("PeptideQuantitativeData",
             peptide_data = pept_data,
             design_matrix = data.frame(Run = "S1", group = "G1", stringsAsFactors = FALSE),
             sample_id = "Run",
             group_id = "group",
             protein_id_column = "Protein.Ids",
             peptide_sequence_column = "Stripped.Sequence",
             q_value_column = "Q.Value",
             raw_quantity_column = "Precursor.Quantity")
  
  result <- calcPeptidesPerProtein(obj)
  
  expect_equal(nrow(result), 2)
  expect_equal(result$n_peptides[result$Protein.Ids == "P1"], 2)
  expect_equal(result$n_peptides[result$Protein.Ids == "P2"], 1)
})

test_that("rollUpPrecursorToPeptideHelper sums quantities correctly", {
  # Mock input table
  input_table <- data.frame(
    Protein.Ids = rep("P1", 4),
    Stripped.Sequence = rep("PEP1", 4),
    Modified.Sequence = c("PEP1_mod1", "PEP1_mod1", "PEP1_mod2", "PEP1_mod2"),
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Quantity = c(10, 20, 30, 40),
    Precursor.Normalised = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )
  
  # Run helper with core_utilisation = 1 (No parallel)
  result <- rollUpPrecursorToPeptideHelper(
    input_table = input_table,
    sample_id_column = Run,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence,
    precursor_quantity_column = Precursor.Quantity,
    precursor_normalised_column = Precursor.Normalised,
    core_utilisation = 1
  )
  
  # result should be a data frame with quantities summed per sample and peptide
  expect_s3_class(result, "data.frame")
  expect_equal(result$Peptide.RawQuantity[result$Run == "S1"], 40)
  expect_equal(result$Peptide.RawQuantity[result$Run == "S2"], 60)
})

# APAF Bioinformatics | test-prot-03-rollup.R | Approved | 2026-03-13
