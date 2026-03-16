# testthat for Proteomics QC Filtering
# Phase 4 of Proteomics GUI Test Strategy

test_that("QC filtered peptide snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp02_qc_filtered_peptide.rds")
  
  if (file.exists(cp_file)) {
    obj <- readRDS(cp_file)
    expect_s4_class(obj, "PeptideQuantitativeData")
    expect_true(nrow(obj@peptide_data) > 0)
  } else {
    skip("Snapshot cp02 not found")
  }
})

test_that("removeEmptyRows removes all-zero or all-NA rows", {
  df <- data.frame(
    ID = c("P1", "P2", "P3"),
    S1 = c(10, 0, NA),
    S2 = c(20, 0, 0),
    Other = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  
  # Columns matching "S\\d" are S1 and S2
  # P1: 10, 20 (Keep)
  # P2: 0, 0 (Remove)
  # P3: NA, 0 (Remove)
  
  result <- removeEmptyRows(df, "S\\d", ID)
  
  expect_equal(nrow(result), 1)
  expect_equal(result$ID, "P1")
})

test_that("removeRowsWithMissingValuesPercentHelper filters correctly", {
  # Mock input table (WIDE format)
  input_table <- data.frame(
    Protein.Ids = c("P1", "P2", "P3"),
    S1 = c(10, NA, NA),
    S2 = c(11, 11, NA),
    S3 = c(12, NA, NA),
    S4 = c(13, 13, 13),
    stringsAsFactors = FALSE
  )
  
  # Mock design matrix
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  
  # cols argument: column names NOT to pivot (the ID columns)
  cols <- "Protein.Ids"
  
  # Run helper
  result <- removeRowsWithMissingValuesPercentHelper(
    input_table = input_table,
    cols = cols,
    design_matrix = design_matrix,
    sample_id = Run,
    row_id = Protein.Ids,
    grouping_variable = group,
    groupwise_percentage_cutoff = 50,
    max_groups_percentage_cutoff = 50,
    temporary_abundance_column = "Abundance"
  )
  
  # According to the source, it returns a data.frame if it succeeded in Step 9
  # Wait, let's check what it actually returns.
  expect_s3_class(result, "data.frame")
  expect_true("P1" %in% result$Protein.Ids)
})

# APAF Bioinformatics | test-prot-02-qc-filtering.R | Approved | 2026-03-13
