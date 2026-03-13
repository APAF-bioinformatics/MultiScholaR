# testthat for Proteomics DA Analysis
# Phase 4 of Proteomics GUI Test Strategy

test_that("DA analysis snapshot is valid", {
  cp_file <- test_path("..", "testdata", "prot_checkpoints", "cp07_da_results.rds")
  
  if (file.exists(cp_file)) {
    results <- readRDS(cp_file)
    expect_type(results, "list")
    expect_true("da_proteins_long" %in% names(results))
    expect_true(nrow(results$da_proteins_long) > 0)
  } else {
    skip("Snapshot cp07 not found")
  }
})

test_that("differentialAbundanceAnalysis works with mock data", {
  # Mock ProteinQuantitativeData
  # Increase number of proteins to 20 to avoid PCA filtering issues
  n_prot <- 20
  pqd <- new("ProteinQuantitativeData",
    protein_quant_table = data.frame(
      Protein.Ids = paste0("P", 1:n_prot),
      # Random abundances with some difference between groups
      S1 = rnorm(n_prot, 10, 1),
      S2 = rnorm(n_prot, 10, 1),
      S3 = rnorm(n_prot, 15, 1),
      S4 = rnorm(n_prot, 15, 1),
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
    protein_id_table = data.frame(
      Protein.Ids = paste0("P", 1:n_prot),
      Gene.Names = paste0("GENE", 1:n_prot),
      stringsAsFactors = FALSE
    ),
    args = list()
  )
  
  # Mock contrasts table
  contrasts <- data.frame(
    full_format = "G2_vs_G1=GroupG2-GroupG1",
    contrasts = "GroupG2-GroupG1",
    friendly_names = "G2_vs_G1",
    stringsAsFactors = FALSE
  )
  
  # Run analysis
  # Note: internal functions might expect specific column names in design matrix
  result <- differentialAbundanceAnalysis(
    theObject = pqd,
    contrasts_tbl = contrasts,
    formula_string = "~ 0 + Group",
    args_row_id = "Protein.Ids",
    group_id = "Group"
  )
  
  expect_type(result, "list")
  expect_true("da_proteins_long" %in% names(result))
  expect_true("G2_vs_G1" %in% result$da_proteins_long$comparison)
})

# APAF Bioinformatics | test-prot-07-da-analysis.R | Approved | 2026-03-13
