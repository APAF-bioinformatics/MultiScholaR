# MultiScholaR: Lipidomics Filtering Test Suite
# Trace ID: tr_20260317_225113_f

library(testthat)
library(MultiScholaR)

# Source updated code directly since we haven't reinstalled the package
source("../../R/allGenerics.R")
source("../../R/func_lipid_s4_objects.R")

context("Lipidomics: Group-Aware Filtering")

test_that("lipidIntensityFiltering group-aware logic works", {
  # 1. Setup Mock Data
  design <- data.frame(
    Run = paste0("S", 1:8),
    group = c(rep("GroupA", 4), rep("GroupB", 4)),
    stringsAsFactors = FALSE
  )
  
  data_vals <- list(
    positive = data.frame(
      Lipid = c("PassLipid", "FailLipid"),
      S1 = c(1000, 1000), # A
      S2 = c(1000, 1000), # A
      S3 = c(1000, 0),    # A
      S4 = c(1000, 0),    # A
      S5 = c(1000, 0),    # B
      S6 = c(1000, 0),    # B
      S7 = c(1000, 0),    # B
      S8 = c(1000, 0),    # B
      stringsAsFactors = FALSE
    )
  )
  
  lipid_obj <- createLipidomicsAssayData(
    lipid_data = data_vals,
    design_matrix = design,
    lipid_id_column = "Lipid",
    sample_id = "Run",
    group_id = "group"
  )
  
  # 2. Run Filtering
  filtered_obj <- lipidIntensityFiltering(
    lipid_obj,
    grouping_variable = "group",
    min_samples_per_group = 3,
    min_groups = 2,
    lipids_intensity_cutoff_percentile = 1
  )
  
  # 3. Verify
  res_data <- filtered_obj@lipid_data$positive
  expect_true("PassLipid" %in% res_data$Lipid)
  expect_false("FailLipid" %in% res_data$Lipid)
})

# <!-- APAF Bioinformatics | test-lipid-filtering.R | Approved | 2026-03-18 -->
