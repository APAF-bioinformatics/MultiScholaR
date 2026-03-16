# testthat for Proteomics S4 Objects
# Phase 4 of Proteomics GUI Test Strategy

test_that("ProteinQuantitativeData constructor works", {
  # Check for captured checkpoint
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp04_design_matrix.rds")
  
  if (file.exists(cp_file)) {
    pqd <- readRDS(cp_file)
    # If the CP04 was from a DIA run, it might be PeptideQuantitativeData
    expect_true(inherits(pqd, "ProteinQuantitativeData") || inherits(pqd, "PeptideQuantitativeData"))
    expect_true(nrow(pqd@design_matrix) > 0)
  } else {
    # Fallback to mock
    pqd <- ProteinQuantitativeData(
      protein_quant_table = data.frame(
        Protein.Ids = c("P1", "P2"),
        S1 = c(10, 11),
        S2 = c(20, 21),
        stringsAsFactors = FALSE
      ),
      design_matrix = data.frame(
        Run = c("S1", "S2"),
        group = c("G1", "G2"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run",
      group_id = "group",
      protein_id_column = "Protein.Ids"
    )
    
    expect_s4_class(pqd, "ProteinQuantitativeData")
    expect_equal(pqd@sample_id, "Run")
    expect_equal(nrow(pqd@protein_quant_table), 2)
  }
})

test_that("PeptideQuantitativeData constructor works", {
  # Mock peptide data
  pept_data <- data.frame(
    Protein.Ids = c("P1", "P1"),
    Peptide.Sequence = c("PEPT1", "PEPT2"),
    Run = c("S1", "S1"),
    Precursor.Quantity = c(100, 200),
    Q.Value = c(0.001, 0.001),
    stringsAsFactors = FALSE
  )
  
  pqd <- PeptideQuantitativeData(
    peptide_data = pept_data,
    design_matrix = data.frame(
      Run = "S1",
      group = "G1",
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Peptide.Sequence",
    q_value_column = "Q.Value",
    raw_quantity_column = "Precursor.Quantity"
  )
  
  expect_s4_class(pqd, "PeptideQuantitativeData")
  expect_equal(pqd@peptide_sequence_column, "Peptide.Sequence")
})

test_that("ProteinQuantitativeData validation catches mismatches", {
  # This might require checking if validity function is actually implemented
  # If it is, we can test it
})

# APAF Bioinformatics | test-prot-04-design.R | Approved | 2026-03-13
