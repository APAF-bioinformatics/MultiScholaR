# testthat for Proteomics Heatmap
# Phase 4 of Proteomics GUI Test Strategy

test_that("generateProtDAHeatmap works with captured snapshot", {
  cp_file <- test_path("..", "testdata", "prot_checkpoints", "cp09_heatmap_input.rds")
  
  if (file.exists(cp_file)) {
    data <- readRDS(cp_file)
    
    # Run function
    result <- do.call(generateProtDAHeatmap, data)
    
    expect_type(result, "list")
    expect_true("plot" %in% names(result))
    expect_s4_class(result$plot, "Heatmap")
  } else {
    skip("Snapshot cp09 not found")
  }
})

test_that("generateProtDAHeatmap works with mock data", {
  # Mock S4 object
  pqd <- new("ProteinQuantitativeData",
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2", "P3"),
      S1 = c(10, 11, 12),
      S2 = c(10.5, 11.5, 12.5),
      S3 = c(15, 16, 17),
      S4 = c(15.5, 16.5, 17.5),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      Group = c("G1", "G1", "G2", "G2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "Group",
    protein_id_column = "Protein.Ids",
    args = list()
  )
  
  # Mock DA results
  da_results <- list(
    da_proteins_long = data.frame(
      Protein.Ids = c("P1", "P2", "P3"),
      comparison = rep("G2_vs_G1", 3),
      log2FC = c(5, -5, 0.1),
      fdr_qvalue = c(0.001, 0.001, 0.8),
      stringsAsFactors = FALSE
    ),
    theObject = pqd
  )
  
  # Run function
  result <- generateProtDAHeatmap(
    da_results_list = da_results,
    selected_contrast = "G2_vs_G1",
    top_n_genes = 10,
    da_q_val_thresh = 0.05
  )
  
  expect_type(result, "list")
  expect_true("plot" %in% names(result))
  expect_s4_class(result$plot, "Heatmap")
})

# APAF Bioinformatics | test-prot-09-heatmap.R | Approved | 2026-03-13
