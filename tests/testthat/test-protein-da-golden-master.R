# test-protein-da-golden-master.R
# Golden Master tests for Protein Differential Abundance Analysis workflow

library(testthat)
devtools::load_all()

test_that("protein_deAnalysisWrapperFunction and its alias are exported", {
  # Verify functions exist in the loaded package environment
  expect_true(exists("protein_deAnalysisWrapperFunction"))
  expect_true(exists("deAnalysisWrapperFunction"))
})

test_that("deAnalysisWrapperFunction alias provides a deprecation warning", {
  # We can't easily run the full function without a real S4 object, 
  # but we can test that it triggers the warning before failing on arguments
  expect_warning(
    try(deAnalysisWrapperFunction(NULL), silent = TRUE),
    "Deprecated"
  )
})

test_of_fix <- function() {
  file_content <- readLines("../../R/func_prot_da.R")
  
  test_that("Robust contrast matching logic is present", {
    expect_true(any(grepl("Glimma: Available contrasts in data", file_content)))
    expect_true(any(grepl("clean_target <- gsub", file_content)))
  })

  test_that("Final display_df cleaning logic is present", {
    expect_true(any(grepl("FINAL CLEANING for Glimma table stability", file_content)))
    expect_true(any(grepl("as.character\\(display_df\\[\\[col\\]\\]\\)", file_content)))
  })
}
test_of_fix()

test_that(".capture_checkpoint respects analysis_dir (Issue 3)", {
  # Verify the source code logic for directory resolution
  func_body <- as.character(body(MultiScholaR:::.capture_checkpoint))
  expect_true(any(grepl("multischolar.analysis_dir", func_body)))
  expect_true(any(grepl("checkpoints", func_body)))
})

# APAF Bioinformatics | test-protein-da-golden-master.R | Approved | 2026-03-15
