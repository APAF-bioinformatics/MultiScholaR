# MultiScholaR: Metabolomics Filtering Test Suite
# Trace ID: tr_20260317_225113_f

library(testthat)
library(MultiScholaR)

# Source updated code directly since we haven't reinstalled the package
source("../../R/allGenerics.R")
source("../../R/func_metab_s4_objects.R")

context("Metabolomics: Group-Aware Filtering")

test_that("metaboliteIntensityFiltering group-aware logic works", {
  # 1. Setup Mock Data
  # We create a simple MetaboliteAssayData object with 2 groups, 4 samples each
  # Group A: S1, S2, S3, S4
  # Group B: S5, S6, S7, S8
  
  design <- data.frame(
    Run = paste0("S", 1:8),
    group = c(rep("GroupA", 4), rep("GroupB", 4)),
    stringsAsFactors = FALSE
  )
  
  # Create a feature that should PASS (above threshold in 3 samples of GroupA and 3 samples of GroupB)
  # Create a feature that should FAIL (above threshold in only 2 samples of GroupA)
  # Threshold will be based on 1st percentile of all data.
  
  # Base values: high (1000), low (10)
  data_vals <- list(
    positive = data.frame(
      Metabolite = c("PassFeature", "FailFeature"),
      S1 = c(1000, 1000), # A
      S2 = c(1000, 1000), # A
      S3 = c(1000, 0),    # A - FailFeature low here
      S4 = c(1000, 0),    # A - FailFeature low here
      S5 = c(1000, 0),    # B
      S6 = c(1000, 0),    # B
      S7 = c(1000, 0),    # B
      S8 = c(1000, 0),    # B
      stringsAsFactors = FALSE
    )
  )
  
  metab_obj <- createMetaboliteAssayData(
    metabolite_data = data_vals,
    design_matrix = design,
    metabolite_id_column = "Metabolite",
    sample_id = "Run",
    group_id = "group"
  )
  
  # 2. Run Filtering
  # Requirement: 3 samples per group, at least 2 groups
  # Percentile = 1. Min value is 10, so threshold will likely be around 10.
  filtered_obj <- metaboliteIntensityFiltering(
    metab_obj,
    grouping_variable = "group",
    min_samples_per_group = 3,
    min_groups = 2,
    metabolites_intensity_cutoff_percentile = 1
  )
  
  # 3. Verify
  # PassFeature: GroupA (4/4 pass), GroupB (4/4 pass) -> 2 groups pass -> KEEP
  # FailFeature: GroupA (2/4 pass), GroupB (0/4 pass) -> 0 groups pass -> REMOVE
  
  res_data <- filtered_obj@metabolite_data$positive
  expect_true("PassFeature" %in% res_data$Metabolite)
  expect_false("FailFeature" %in% res_data$Metabolite)
})

# <!-- APAF Bioinformatics | test-metab-filtering.R | Approved | 2026-03-18 -->
