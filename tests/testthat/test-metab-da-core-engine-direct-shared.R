# fidelity-coverage-compare: shared
library(testthat)

runTestsContrastsMetabDA <- get(
  "runTestsContrastsMetabDA",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("metabolomics DA core engine preserves validation guards", {
  skip_if_not_installed("limma")

  assay_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("Sample_1", "Sample_2"))
  )

  no_match_design <- data.frame(
    sample_id = c("Other_1", "Other_2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    runTestsContrastsMetabDA(
      data = assay_matrix,
      contrast_strings = "groupB-groupA",
      design_matrix = no_match_design,
      formula_string = "~ 0 + group",
      sample_id_col = "sample_id"
    ),
    "No common samples between data matrix and design matrix"
  )

  valid_design <- data.frame(
    sample_id = c("Sample_1", "Sample_2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    runTestsContrastsMetabDA(
      data = assay_matrix,
      contrast_strings = "Treatment-Control",
      design_matrix = valid_design,
      formula_string = "~ 0 + group",
      sample_id_col = "sample_id"
    ),
    "references undefined levels"
  )
})

test_that("metabolomics DA core engine preserves standard limma result assembly", {
  skip_if_not_installed("limma")
  skip_if_not_installed("qvalue")

  assay_matrix <- matrix(
    c(
      10, 11, 20, 21,
      12, 13, 24, 25,
      30, 29, 15, 14,
      28, 27, 13, 12
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("M1", "M2", "M3", "M4"), c("Sample_1", "Sample_2", "Sample_3", "Sample_4"))
  )

  design_matrix <- data.frame(
    sample_id = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- runTestsContrastsMetabDA(
    data = assay_matrix,
    contrast_strings = "groupB-groupA",
    design_matrix = design_matrix,
    formula_string = "~ 0 + group",
    sample_id_col = "sample_id",
    treat_lfc_cutoff = NA
  )

  expect_named(result, c("results", "fit.eb", "qvalue_warnings"))
  expect_named(result$results, "groupB-groupA")
  expect_true(all(c("logFC", "P.Value", "raw_pvalue", "fdr_qvalue", "fdr_value_bh") %in% names(result$results[[1L]])))
  expect_identical(sort(rownames(result$results[[1L]])), c("M1", "M2", "M3", "M4"))
  expect_true(inherits(result$fit.eb, "MArrayLM"))
})

test_that("metabolomics DA core engine preserves TREAT branch assembly", {
  skip_if_not_installed("limma")
  skip_if_not_installed("qvalue")

  assay_matrix <- matrix(
    c(
      10, 10.5, 20, 20.5,
      11, 11.5, 21, 21.5,
      30, 30.5, 15, 15.5,
      28, 28.5, 13, 13.5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("M1", "M2", "M3", "M4"), c("Sample_1", "Sample_2", "Sample_3", "Sample_4"))
  )

  design_matrix <- data.frame(
    sample_id = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- runTestsContrastsMetabDA(
    data = assay_matrix,
    contrast_strings = "groupB-groupA",
    design_matrix = design_matrix,
    formula_string = "~ 0 + group",
    sample_id_col = "sample_id",
    treat_lfc_cutoff = 1
  )

  expect_named(result$results, "groupB-groupA")
  expect_true(all(c("logFC", "P.Value", "raw_pvalue", "fdr_qvalue", "fdr_value_bh") %in% names(result$results[[1L]])))
  expect_true(inherits(result$fit.eb, "MArrayLM"))
})
