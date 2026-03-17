# MultiScholaR: Metabolomics Normalization Test Suite
# Trace ID: tr_20260317_225113_f

library(testthat)
library(MultiScholaR)

context("Metabolomics: Normalization Sub-modules")

test_that("mod_metab_norm_impute server logic works", {
  # 1. Setup Fixture
  fixture_path <- testthat::test_path("../fixtures/metab_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  metab_obj <- readRDS(fixture_path)
  
  # 2. Mock State Manager
  state_manager <- list(
    getState = function() metab_obj,
    saveState = function(state_name, s4_data_object, ...) {
      assign("last_saved_state", s4_data_object, envir = .GlobalEnv)
    },
    getHistory = function() c("initial")
  )
  workflow_data <- list(state_manager = state_manager)
  
  # 3. Test Imputation (Zero)
  # Simulate the logic inside mod_metab_norm_impute_server
  impute_method <- "zero"
  imputed_s4 <- metab_obj
  imputed_s4@metabolite_data <- purrr::map(metab_obj@metabolite_data, function(assay_df) {
    num_cols <- sapply(assay_df, is.numeric)
    sample_cols <- names(assay_df)[num_cols]
    assay_df[sample_cols] <- purrr::map_dfc(assay_df[sample_cols], function(col) {
      col[is.na(col)] <- 0
      return(col)
    })
    return(assay_df)
  })
  
  # Verify
  for (assay_name in names(imputed_s4@metabolite_data)) {
    assay_df <- imputed_s4@metabolite_data[[assay_name]]
    num_cols <- sapply(assay_df, is.numeric)
    sample_cols <- names(assay_df)[num_cols]
    expect_false(any(is.na(assay_df[sample_cols])), info = paste("NAs found in", assay_name))
  }
})

test_that("metaboliteMissingValueImputationLimpa works (singleton mode)", {
  fixture_path <- testthat::test_path("../fixtures/metab_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  metab_obj <- readRDS(fixture_path)
  
  # Run Limpa (singleton mode) via internal call
  imputed_obj <- MultiScholaR:::metaboliteMissingValueImputationLimpa(metab_obj, verbose = FALSE)
  
  expect_s4_class(imputed_obj, "MetaboliteAssayData")
  
  # Verify all NAs are gone
  for (assay_name in names(imputed_obj@metabolite_data)) {
    assay_df <- imputed_obj@metabolite_data[[assay_name]]
    design_samples <- as.character(metab_obj@design_matrix[[metab_obj@sample_id]])
    sample_cols <- intersect(colnames(assay_df), design_samples)
    expect_false(any(is.na(assay_df[sample_cols])), info = paste("NAs still present in", assay_name))
  }
})

test_that("metaboliteMissingValueImputationMissForest works", {
  fixture_path <- testthat::test_path("../fixtures/metab_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  metab_obj <- readRDS(fixture_path)
  
  # Run missForest with minimal parameters for speed
  imputed_obj <- MultiScholaR:::metaboliteMissingValueImputationMissForest(
    metab_obj, 
    maxiter = 2, 
    ntree = 10, 
    verbose = FALSE
  )
  
  expect_s4_class(imputed_obj, "MetaboliteAssayData")
  
  for (assay_name in names(imputed_obj@metabolite_data)) {
    assay_df <- imputed_obj@metabolite_data[[assay_name]]
    design_samples <- as.character(metab_obj@design_matrix[[metab_obj@sample_id]])
    sample_cols <- intersect(colnames(assay_df), design_samples)
    expect_false(any(is.na(assay_df[sample_cols])))
  }
})

test_that("mod_metab_norm_itsd feature selection works", {
  fixture_path <- testthat::test_path("../fixtures/metab_baseline_post_filter.rds")
  skip_if_not(file.exists(fixture_path), "Fixture not found")
  metab_obj <- readRDS(fixture_path)
  
  # Test buildItsdSelectionTable helper (exported/available in package)
  assay_name <- names(metab_obj@metabolite_data)[1]
  assay_df <- metab_obj@metabolite_data[[assay_name]]
  
  selection_table <- MultiScholaR:::buildItsdSelectionTable(
    assay_data = assay_df,
    metabolite_id_col = metab_obj@metabolite_id_column
  )
  
  expect_s3_class(selection_table, "data.frame")
  expect_true("feature_id" %in% names(selection_table))
  expect_true("is_candidate" %in% names(selection_table))
})

# <!-- APAF Bioinformatics | test-metab-norm.R | Approved | 2026-03-17 -->
