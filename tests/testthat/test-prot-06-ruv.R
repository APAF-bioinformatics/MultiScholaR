# testthat for Proteomics RUV
# Phase 4 of Proteomics GUI Test Strategy

test_that("RUV snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp06_ruv_corrected.rds")
  
  if (file.exists(cp_file)) {
    obj <- readRDS(cp_file)
    expect_true(inherits(obj, "ProteinQuantitativeData") || inherits(obj, "PeptideQuantitativeData"))
    expect_true(!is.null(obj@args$ruvIII_C_Varying))
  } else {
    skip("Snapshot cp06 not found")
  }
})

test_that("getRuvIIIReplicateMatrixHelper creates correct matrix", {
  design <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  
  result <- getRuvIIIReplicateMatrixHelper(design, Sample, Group)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 2))
  expect_equal(result["S1", "G1"], 1)
  expect_equal(result["S1", "G2"], 0)
})

test_that("findBestK works with mock cancorplot data", {
  # Mock data frame matching the structure used in findBestK
  mock_data <- data.frame(
    featureset = c(rep("Control", 3), rep("All", 3)),
    cc = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
    K = c(1, 2, 3, 1, 2, 3),
    stringsAsFactors = FALSE
  )
  
  # Difference (All - Control):
  # K=1: 0.5 - 0.1 = 0.4
  # K=2: 0.7 - 0.2 = 0.5
  # K=3: 0.9 - 0.3 = 0.6 (Max)
  
  mock_plot <- list(data = mock_data)
  
  result <- findBestK(mock_plot)
  
  expect_equal(result, 3)
})

test_that("getNegCtrlProtAnovaHelper identifies controls", {
  # Mock data matrix (log2 abundances)
  mat <- matrix(rnorm(100*4), nrow=100, ncol=4)
  rownames(mat) <- paste0("P", 1:100)
  colnames(mat) <- c("S1", "S2", "S3", "S4")
  
  # Make some proteins NOT significant (good controls)
  # In ANOVA, large p-values/q-values mean not significant
  # getNegCtrlProtAnovaHelper selects genes with HIGHEST q-values
  
  design <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  
  # We need to mock the internal ANOVA/qvalue calls?
  # Actually, let's just run it and see if it returns a logical vector
  
  # Note: This function calls qvalue() which might fail on small mock data
  # Use BH method to be safer
  result <- getNegCtrlProtAnovaHelper(
    data_matrix = mat,
    design_matrix = design,
    grouping_variable = "group",
    percentage_as_neg_ctrl = 20,
    ruv_fdr_method = "BH"
  )
  
  expect_type(result, "logical")
  expect_equal(length(result), 100)
  expect_equal(sum(result), 20) # 20% of 100
})

# APAF Bioinformatics | test-prot-06-ruv.R | Approved | 2026-03-13
