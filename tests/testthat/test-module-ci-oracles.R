library(testthat)

test_that("module CI identity oracles detect unexplained sample and feature drift", {
  expect_no_error(module_ci_assert_sample_identity(
    before = c("WT_1", "WT_2", "KO_1"),
    after = c("WT_1", "KO_1"),
    expected_dropped = "WT_2"
  ))
  expect_failure(module_ci_assert_sample_identity(
    before = c("WT_1", "WT_2"),
    after = "WT_1"
  ))

  expect_no_error(module_ci_assert_feature_identity(
    before = c("P1", "P2", "P3"),
    after = c("P1", "P3"),
    expected_dropped = "P2"
  ))
})

test_that("module CI assay provenance oracle supports lists, state payloads, and result tables", {
  expect_no_error(module_ci_assert_assay_provenance(
    list(LCMS_Pos = data.frame(x = 1), LCMS_Neg = data.frame(x = 2)),
    c("LCMS_Pos", "LCMS_Neg")
  ))
  expect_no_error(module_ci_assert_assay_provenance(
    list(assay_names = c("LCMS_Pos", "GCMS")),
    c("LCMS_Pos", "GCMS")
  ))
  expect_no_error(module_ci_assert_assay_provenance(
    data.frame(assay = c("A", "A", "B")),
    c("A", "B")
  ))
})

test_that("module CI table and numeric oracles validate common module outputs", {
  fixture <- file.path(module_ci_fixture_root(), "fx", "metabolomics", "metabolite_assay.tsv")
  table <- module_ci_assert_table_schema(
    fixture,
    required_columns = c("Feature.Name", "WT_1", "KO_1"),
    min_rows = 3L
  )

  module_ci_assert_numeric_finite(table, columns = c("WT_1", "WT_2", "KO_1", "KO_2"))
  module_ci_assert_no_duplicate_keys(table, "Feature.Name")
  expect_failure(module_ci_assert_no_duplicate_keys(
    data.frame(id = c("A", "A")),
    "id"
  ))
})

test_that("module CI design alignment oracle validates sample contracts", {
  design <- data.frame(sample = c("WT_1", "WT_2", "KO_1"), group = c("WT", "WT", "KO"))
  expect_no_error(module_ci_assert_design_alignment(design, c("KO_1", "WT_1", "WT_2")))
  expect_error(
    module_ci_assert_design_alignment(data.frame(id = "WT_1"), "WT_1"),
    "sample_col not found"
  )
})
