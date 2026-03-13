# testthat for Proteomics Volcano Plot
# Phase 4 of Proteomics GUI Test Strategy

test_that("generateProtDAVolcanoPlotGlimma works with captured snapshot", {
  cp_file <- test_path("..", "testdata", "prot_checkpoints", "cp08_volcano_input.rds")
  
  if (file.exists(cp_file)) {
    data <- readRDS(cp_file)
    
    # Run the function
    # Note: we need to handle potential NULLs in captured data if any
    p <- do.call(generateProtDAVolcanoPlotGlimma, data)
    
    # Assertions
    expect_s3_class(p, "htmlwidget")
    expect_s3_class(p, "glimmaXY")
  } else {
    skip("Snapshot cp08 not found")
  }
})

test_that("generateProtDAVolcanoPlotGlimma handles mock data", {
  # Mock CP08 input
  mock_cp08 <- list(
    da_results_list = list(
      da_proteins_long = data.frame(
        Protein.Ids = c("P1", "P2", "P3", "P4"),
        comparison = rep("T_vs_C", 4),
        log2FC = c(2.5, -3.0, 0.5, -0.2),
        fdr_qvalue = c(0.001, 0.005, 0.5, 0.8),
        raw_pvalue = c(0.0001, 0.0005, 0.1, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    selected_contrast = "T_vs_C",
    da_q_val_thresh = 0.05,
    args_row_id = "Protein.Ids"
  )
  
  # Run function
  result <- do.call(generateProtDAVolcanoPlotGlimma, mock_cp08)
  
  # Assertions
  expect_s3_class(result, "htmlwidget")
  expect_s3_class(result, "glimmaXY")
})

test_that("generateProtDAVolcanoStatic handles mock data", {
  # Mock inputs
  da_results_list <- list(
    da_proteins_long = data.frame(
      Protein.Ids = paste0("P", 1:100),
      comparison = rep("T_vs_C", 100),
      log2FC = rnorm(100, 0, 2),
      fdr_qvalue = runif(100, 0, 1),
      raw_pvalue = runif(100, 0, 1),
      stringsAsFactors = FALSE
    )
  )
  
  # Run function
  p <- generateProtDAVolcanoStatic(
    da_results_list = da_results_list,
    selected_contrast = "T_vs_C",
    da_q_val_thresh = 0.05,
    lfc_threshold = 1.0
  )
  
  # Assertions
  expect_s3_class(p, "ggplot")
})

# APAF Bioinformatics | test-prot-08-volcano.R | Approved | 2026-03-13
