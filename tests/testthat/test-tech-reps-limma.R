library(testthat)
library(dplyr)
library(limma)

# Load the package
devtools::load_all(".")

test_that("runTestsContrasts handles technical replicates correctly", {
  # 1. Create mock data with technical replicates
  # 4 biological samples, each with 2 technical replicates = 8 total runs, two biological samples per group.
  # (KIT-035: with one biological sample per group the block equals the group, the design spans it and limma
  # cannot estimate a within-block correlation; such a design is standard, with the reason recorded.)
  data_matrix <- matrix(
    c(10, 10.2, 10.5, 10.4, 12, 12.1, 11.8, 11.9,  # Protein 1
      15, 14.8, 15.1, 15.3, 13, 13.2, 12.9, 13.1), # Protein 2
    nrow = 2, byrow = TRUE
  )
  colnames(data_matrix) <- paste0("S", 1:8)
  rownames(data_matrix) <- c("P1", "P2")
  
  # Design matrix with replicates column
  design_matrix <- data.frame(
    Run = paste0("S", 1:8),
    group = rep(c("A", "B"), each = 4),
    replicates = c(1, 1, 2, 2, 3, 3, 4, 4), # Biological sample ids; each run twice (technical replicates)
    stringsAsFactors = FALSE
  )
  rownames(design_matrix) <- design_matrix$Run
  
  contrast_strings <- "groupB-groupA"
  formula_string <- "~ 0 + group"
  
  # 2. Run runTestsContrasts and capture messages
  # We expect "Detected technical replicates. Calculating duplicateCorrelation..."
  
  expect_message(
    results <- runTestsContrasts(
      data = data_matrix,
      contrast_strings = contrast_strings,
      design_matrix = design_matrix,
      formula_string = formula_string
    ),
    "Detected technical replicates"
  )
  
  # 3. Verify results structure
  expect_type(results, "list")
  expect_named(results, c("results", "fit.eb", "qvalue_warnings"))
  expect_s3_class(results$results[[1]], "data.frame")
  expect_true(nrow(results$results[[1]]) == 2)
})

test_that("runTestsContrasts handles NO technical replicates correctly", {
  # 1. Create mock data WITHOUT technical replicates
  data_matrix <- matrix(
    c(10, 11, 12, 13,
      15, 14, 13, 12),
    nrow = 2, byrow = TRUE
  )
  colnames(data_matrix) <- c("S1", "S2", "S3", "S4")
  rownames(data_matrix) <- c("P1", "P2")
  
  # Design matrix with unique replicates
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    replicates = c(1, 2, 3, 4), # Unique = No tech reps
    stringsAsFactors = FALSE
  )
  rownames(design_matrix) <- design_matrix$Run
  
  contrast_strings <- "groupB-groupA"
  formula_string <- "~ 0 + group"
  
  # 2. Run runTestsContrasts and capture messages
  # We expect "No technical replicates detected"
  
  expect_message(
    results <- runTestsContrasts(
      data = data_matrix,
      contrast_strings = contrast_strings,
      design_matrix = design_matrix,
      formula_string = formula_string
    ),
    "No technical replicates detected"
  )
})
# APAF Bioinformatics | tests/testthat/test-tech-reps-limma.R | Approved | 2026-03-14
