# test-protein-da-golden-master.R
# Golden Master tests for Protein Differential Abundance Analysis workflow

library(testthat)

test_that("protein_deAnalysisWrapperFunction and its alias are exported", {
  # Verify functions exist in the loaded package environment
  expect_true(exists("protein_deAnalysisWrapperFunction"))
  expect_true(exists("deAnalysisWrapperFunction"))
})

test_that("deAnalysisWrapperFunction alias provides a deprecation warning", {
  # We can't easily run the full function without a real S4 object, 
  # but we can test that it triggers the warning before failing on arguments
  expect_warning(
    try(deAnalysisWrapperFunction(NULL), silent = TRUE),
    "Deprecated"
  )
})

test_that("Robust contrast matching logic is present", {
  func_body <- as.character(body(MultiScholaR::generateProtDAVolcanoPlotGlimma))
  expect_true(any(grepl("Glimma: Available contrasts in data", func_body, fixed = TRUE)))
  expect_true(any(grepl("clean_target <- gsub", func_body, fixed = TRUE)))
})

test_that("Final consumed annotation cleaning remains and discarded preparation is absent", {
  captured <- NULL
  testthat::local_mocked_bindings(glimmaXY = function(...) {
    captured <<- list(...)
    htmlwidgets::createWidget("glimmaXY", list(), package = "Glimma")
  }, .package = "Glimma")
  rows <- data.frame(Protein.Ids = c("P1", "P2"), comparison = "T_vs_C",
    log2FC = c(2, -1), fdr_qvalue = c(0.01, 0.5), raw_pvalue = c(0.001, 0.2),
    custom.factor = factor(c("label", NA), levels = c("label", "unused")))
  result <- generateProtDAVolcanoPlotGlimma(list(da_proteins_long = rows),
    selected_contrast = "T_vs_C", args_row_id = "Protein.Ids")
  expect_s3_class(result, "htmlwidget")
  expect_identical(captured$anno$Protein_Ids, rows$Protein.Ids)
  expect_identical(captured$anno$custom_factor, c("label", ""))
  expect_true(all(vapply(captured$anno, is.character, logical(1))))
  expect_false(anyNA(captured$anno))
  expect_false(any(c("logFC", "FDR", "Status") %in% names(captured$anno)))
  expect_identical(captured$x, rows$log2FC)
  expect_false("display_df" %in% all.names(body(generateProtDAVolcanoPlotGlimma)))
})

test_that(".capture_checkpoint respects analysis_dir (Issue 3)", {
  # Verify the source code logic for directory resolution
  func_body <- as.character(body(MultiScholaR:::.capture_checkpoint))
  expect_true(any(grepl("multischolar.analysis_dir", func_body)))
  expect_true(any(grepl("checkpoints", func_body)))
})

# APAF Bioinformatics | test-protein-da-golden-master.R | Approved | 2026-03-15
