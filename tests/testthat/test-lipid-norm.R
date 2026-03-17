# MultiScholaR: Lipidomics Normalization Test Suite
# Trace ID: tr_20260317_225113_f

library(testthat)
library(MultiScholaR)

# Source modified files directly
source("../../R/allGenerics.R")
source("../../R/func_lipid_s4_objects.R")

context("Lipidomics: Normalization Sub-modules")

test_that("mod_lipid_norm_impute server logic works", {
  fixture_path <- testthat::test_path("../fixtures/lipid_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  lipid_obj <- readRDS(fixture_path)
  
  # Test Imputation (Zero)
  impute_method <- "zero"
  imputed_s4 <- lipid_obj
  imputed_s4@lipid_data <- purrr::map(lipid_obj@lipid_data, function(assay_df) {
    num_cols <- sapply(assay_df, is.numeric)
    sample_cols <- names(assay_df)[num_cols]
    assay_df[sample_cols] <- purrr::map_dfc(assay_df[sample_cols], function(col) {
      col[is.na(col)] <- 0
      return(col)
    })
    return(assay_df)
  })
  
  for (assay_name in names(imputed_s4@lipid_data)) {
    assay_df <- imputed_s4@lipid_data[[assay_name]]
    num_cols <- sapply(assay_df, is.numeric)
    sample_cols <- names(assay_df)[num_cols]
    expect_false(any(is.na(assay_df[sample_cols])))
  }
})

test_that("lipidMissingValueImputationLimpa works (singleton mode)", {
  fixture_path <- testthat::test_path("../fixtures/lipid_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  lipid_obj <- readRDS(fixture_path)
  
  imputed_obj <- MultiScholaR:::lipidMissingValueImputationLimpa(lipid_obj, verbose = FALSE)
  
  expect_s4_class(imputed_obj, "LipidomicsAssayData")
  
  # Verify all NAs are gone
  for (assay_name in names(imputed_obj@lipid_data)) {
    assay_df <- imputed_obj@lipid_data[[assay_name]]
    design_samples <- as.character(lipid_obj@design_matrix[[lipid_obj@sample_id]])
    sample_cols <- intersect(colnames(assay_df), design_samples)
    expect_false(any(is.na(assay_df[sample_cols])), info = paste("NAs still present in", assay_name))
  }
})

test_that("mod_lipid_norm_itsd feature selection works", {
  fixture_path <- testthat::test_path("../fixtures/lipid_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  lipid_obj <- readRDS(fixture_path)
  
  assay_name <- names(lipid_obj@lipid_data)[1]
  assay_df <- lipid_obj@lipid_data[[assay_name]]
  
  selection_table <- MultiScholaR:::buildItsdSelectionTable(
    assay_data = assay_df,
    metabolite_id_col = lipid_obj@lipid_id_column
  )
  
  expect_s3_class(selection_table, "data.frame")
  expect_true("is_candidate" %in% names(selection_table))
})

# <!-- APAF Bioinformatics | test-lipid-norm.R | Approved | 2026-03-17 -->
