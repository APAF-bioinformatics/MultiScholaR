# testthat for Proteomics RUV
# Phase 4 of Proteomics GUI Test Strategy

test_that("RUV snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp06_ruv_corrected.rds")
  
  if (file.exists(cp_file)) {
    obj <- tryCatch(
      readRDS(cp_file),
      error = function(e) {
        skip(sprintf("Snapshot cp06 is unavailable in this checkout: %s", e$message))
      }
    )
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

test_that("findBestK works with mock cancorplot data (deprecated wrapper)", {
  mock_data <- data.frame(
    featureset = c(rep("Control", 3), rep("All", 3)),
    cc = c(0.1, 0.2, 0.3, 0.5, 0.7, 0.9),
    K = c(1, 2, 3, 1, 2, 3),
    stringsAsFactors = FALSE
  )
  # delta: K=1:0.4, K=2:0.5, K=3:0.6  max=0.6, target=0.57; only K=3 qualifies
  mock_plot <- list(data = mock_data)
  result <- suppressWarnings(findBestK(mock_plot))
  expect_equal(result, 3L)
})

# ---------------------------------------------------------------------------
# Phase 1 shared math contract tests (RUV-001)
# ---------------------------------------------------------------------------

makeMockCancorPlot <- function(featuresets, k_values, cc_values) {
  list(data = data.frame(
    featureset = featuresets
    , K = k_values
    , cc = cc_values
    , stringsAsFactors = FALSE
  ))
}

makeDeltaPlot <- function(k_all, delta_all, ctrl_base = 0.2) {
  all_cc <- ctrl_base + delta_all
  makeMockCancorPlot(
    featuresets = c(rep("Control", length(k_all)), rep("All", length(k_all)))
    , k_values = c(k_all, k_all)
    , cc_values = c(rep(ctrl_base, length(k_all)), all_cc)
  )
}

makeGgplotCancorPlot <- function(k_all, delta_all, ctrl_base = 0.2) {
  all_cc <- ctrl_base + delta_all
  pd <- data.frame(
    featureset = c(rep("Control", length(k_all)), rep("All", length(k_all)))
    , K = c(k_all, k_all)
    , cc = c(rep(ctrl_base, length(k_all)), all_cc)
    , stringsAsFactors = FALSE
  )
  ggplot2::ggplot(pd, ggplot2::aes(K, cc, colour = featureset)) + ggplot2::geom_point()
}

test_that("findBestKElbow returns smallest K at a plateau", {
  # delta: K=1:0.5, K=2:0.91, K=3:0.95 — both K=2 and K=3 satisfy target=0.9025
  plt <- makeDeltaPlot(1:3, c(0.5, 0.91, 0.95))
  expect_identical(findBestKElbow(plt), 2L)
})

test_that("findBestKElbow returns smallest K on exact ties", {
  # delta: K=1:0.5, K=2:0.9, K=3:0.9 — exact tie at max
  plt <- makeDeltaPlot(1:3, c(0.5, 0.9, 0.9))
  expect_identical(findBestKElbow(plt), 2L)
})

test_that("findBestKElbow returns largest K on strictly increasing gap curves", {
  # delta: K=1:0.1, K=2:0.5, K=3:1.0 — strictly increasing; only K=3 >= target 0.95
  plt <- makeDeltaPlot(1:3, c(0.1, 0.5, 1.0))
  expect_identical(findBestKElbow(plt), 3L)
})

test_that("findBestKElbow returns 1L for valid weak-effect curves", {
  plt <- makeDeltaPlot(1:3, c(0.01, 0.02, 0.03))
  expect_identical(findBestKElbow(plt), 1L)
})

test_that("findBestKElbow returns NA_integer_ for malformed input", {
  expect_identical(findBestKElbow(list()), NA_integer_)
  expect_identical(findBestKElbow(NULL), NA_integer_)
  expect_identical(findBestKElbow(list(data = data.frame(x = 1))), NA_integer_)
})

test_that("findBestKElbow returns NA_integer_ when All is missing", {
  plt <- makeMockCancorPlot(
    featuresets = c("Control", "Control", "Control")
    , k_values = 1:3
    , cc_values = c(0.2, 0.3, 0.4)
  )
  expect_identical(findBestKElbow(plt), NA_integer_)
})

test_that("findBestKElbow returns NA_integer_ when Control is missing", {
  plt <- makeMockCancorPlot(
    featuresets = c("All", "All", "All")
    , k_values = 1:3
    , cc_values = c(0.5, 0.7, 0.9)
  )
  expect_identical(findBestKElbow(plt), NA_integer_)
})

test_that("findBestKElbow is invariant to shuffled row order", {
  plt_ordered <- makeDeltaPlot(1:3, c(0.5, 0.91, 0.95))
  shuffled_data <- plt_ordered$data[c(4, 1, 5, 3, 6, 2), ]
  plt_shuffled <- list(data = shuffled_data)
  expect_identical(findBestKElbow(plt_ordered), findBestKElbow(plt_shuffled))
})

test_that("calculateSeparationScore aligns by K not by row position", {
  # Control at K=2,1 order; All at K=1,2 order
  # K-aligned: K=1: 0.8-0.2=0.6, K=2: 0.6-0.5=0.1 → max=0.6
  # Positional (old): (0.8-0.5, 0.6-0.2) = (0.3, 0.4) → max=0.4
  misaligned <- list(data = data.frame(
    featureset = c("Control", "Control", "All", "All")
    , cc = c(0.5, 0.2, 0.8, 0.6)
    , K = c(2, 1, 1, 2)
    , stringsAsFactors = FALSE
  ))
  expect_equal(calculateSeparationScore(misaligned, "max_difference"), 0.6)
})

test_that("calculateSeparationScore auc matches ordered and shuffled input", {
  ordered_data <- data.frame(
    featureset = c("Control", "Control", "Control", "All", "All", "All")
    , cc = c(0.2, 0.3, 0.4, 0.3, 0.7, 0.4)
    , K = c(1, 2, 3, 1, 2, 3)
    , stringsAsFactors = FALSE
  )
  shuffled_data <- ordered_data[c(4, 1, 5, 3, 6, 2), ]
  ordered_auc <- calculateSeparationScore(list(data = ordered_data), "auc")
  shuffled_auc <- calculateSeparationScore(list(data = shuffled_data), "auc")
  expect_equal(ordered_auc, shuffled_auc)
})

test_that("calculateCompositeScore handles max_acceptable_k == 1", {
  # k=1 at boundary: no penalty → composite equals separation_score
  expect_equal(calculateCompositeScore(0.5, 1, 0.5, 1), 0.5)
  # k=2 above max: heavy penalty → composite = 0
  expect_equal(calculateCompositeScore(0.5, 2, 0.5, 1), 0)
})

test_that("calculateCompositeScore returns NA_real_ on non-finite inputs", {
  expect_identical(calculateCompositeScore(Inf, 2, 0.5, 3), NA_real_)
  expect_identical(calculateCompositeScore(NaN, 2, 0.5, 3), NA_real_)
  expect_identical(calculateCompositeScore(0.5, Inf, 0.5, 3), NA_real_)
  expect_identical(calculateCompositeScore(0.5, NA_real_, 0.5, 3), NA_real_)
  expect_identical(calculateCompositeScore(NA_real_, 2, 0.5, 3), NA_real_)
})

test_that("calculateCompositeScore clamps penalties into [0, 1]", {
  # k=100 far above max → penalty clamped to 1, composite = 0
  result <- calculateCompositeScore(1.0, 100, 0.5, 3)
  expect_gte(result, 0)
  expect_lte(result, 1.0)
})

test_that("findBestK emits a deprecation warning and delegates correctly", {
  plt <- makeDeltaPlot(1:3, c(0.5, 0.91, 0.95))
  expect_warning(result <- findBestK(plt), "findBestKElbow")
  expect_identical(result, findBestKElbow(plt))
})

test_that("findBestKForAssayList propagates scalar K and preserves names", {
  plt <- makeGgplotCancorPlot(1:3, c(0.5, 0.91, 0.95))
  result <- findBestKForAssayList(list(Assay1 = plt, Assay2 = plt))
  expect_identical(names(result), c("Assay1", "Assay2"))
  expect_identical(result$Assay1, 2L)
  expect_identical(result$Assay2, 2L)
})

test_that("duplicate K rows within a featureset are treated as invalid input", {
  dup_data <- data.frame(
    featureset = c("Control", "Control", "Control", "All", "All", "All")
    , cc = c(0.2, 0.3, 0.25, 0.5, 0.7, 0.9)
    , K = c(1, 2, 2, 1, 2, 3)  # Control has K=2 twice
    , stringsAsFactors = FALSE
  )
  plt_dup <- list(data = dup_data)
  expect_identical(findBestKElbow(plt_dup), NA_integer_)
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
