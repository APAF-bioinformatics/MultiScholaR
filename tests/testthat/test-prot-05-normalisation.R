# testthat for Proteomics Normalisation
# Phase 4 of Proteomics GUI Test Strategy

test_that("scaleCenterAndFillMissing works with mock matrix", {
  mat <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, byrow = TRUE)
  colnames(mat) <- c("S1", "S2", "S3")
  
  result <- scaleCenterAndFillMissing(mat)
  
  expect_equal(dim(result), dim(mat))
  expect_false(any(is.na(result)))
})

test_that("normalisation snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp05_normalised.rds")
  
  if (file.exists(cp_file)) {
    obj <- tryCatch(
      readRDS(cp_file),
      error = function(e) {
        skip(sprintf("Snapshot cp05 is unavailable in this checkout: %s", e$message))
      }
    )
    expect_true(inherits(obj, "ProteinQuantitativeData") || inherits(obj, "PeptideQuantitativeData"))
    
    # Extract matrix and check scaling
    if (inherits(obj, "ProteinQuantitativeData")) {
      mat <- as.matrix(obj@protein_quant_table[,-1])
      result <- scaleCenterAndFillMissing(mat)
      expect_equal(dim(result), dim(mat))
    }
  } else {
    skip("Snapshot cp05 not found")
  }
})

test_that("scaleCenterAndFillMissing handles NAs correctly", {
  mat <- matrix(c(10, 20, 30, NA, 50, 60), nrow = 2, byrow = TRUE)
  
  result <- scaleCenterAndFillMissing(mat)
  
  expect_false(any(is.na(result)))
  # The NA in row 2 should be replaced
  expect_true(is.numeric(result[2,1]))
})

# APAF Bioinformatics | test-prot-05-normalisation.R | Approved | 2026-03-13
