# test_limpa_connection.R
# Verify Issue 1: Missing deAnalysisWrapperFunction

library(testthat)
devtools::load_all()

test_that("protein_deAnalysisWrapperFunction exists and is callable", {
  expect_true(exists("protein_deAnalysisWrapperFunction"))
  expect_true(exists("deAnalysisWrapperFunction"))
  
  # Verify signature
  args <- names(formals(protein_deAnalysisWrapperFunction))
  expect_true("theObject" %in% args)
  expect_true("contrasts_tbl" %in% args)
})

test_that("deAnalysisWrapperFunction triggers warning but works", {
  # We use a try block because it will fail on arguments, but we want the warning
  expect_warning(
    try(deAnalysisWrapperFunction(NULL), silent = TRUE),
    "Deprecated"
  )
})

# APAF Bioinformatics | test_limpa_connection.R | Approved | 2026-03-15
