# fidelity-coverage-compare: shared
library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingMetabImportSplitFiles <- function() {
  required_paths <- c(
    "R/mod_metab_import_column_helpers.R",
    "R/mod_metab_import_display_helpers.R",
    "R/mod_metab_import_processing_helpers.R",
    "R/mod_metab_import_server_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only metab import split file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingMetabImportSplitFiles()

test_that("metabolomics import column resolver preserves case-insensitive lookup and fallback", {
  headers <- c("Peak ID", "Metabolite_Name", "Sample_A")

  expect_identical(resolveMetabImportColumnName(headers, "peak id"), "Peak ID")
  expect_identical(resolveMetabImportColumnName(headers, "METABOLITE_name"), "Metabolite_Name")
  expect_identical(resolveMetabImportColumnName(headers, "missing_column"), "missing_column")
  expect_identical(resolveMetabImportColumnName(headers, ""), "")
  expect_null(resolveMetabImportColumnName(headers, NULL))
})

test_that("metabolomics import sample-column resolver preserves custom pattern and fallback order", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2"),
    Sample_A = c(10, 11),
    sample_b = c(20, 21),
    Annotation = c("A", "B"),
    check.names = FALSE
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "custom",
      sampleColsPattern = "^sample",
      importResult = list(sample_columns = c("ignored"))
    ),
    c("Sample_A", "sample_b")
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "custom",
      sampleColsPattern = "^missing",
      importResult = list(sample_columns = c("Imported_A", "Imported_B"))
    ),
    c("Imported_A", "Imported_B")
  )
})

test_that("metabolomics import workflow payload builder preserves assay assembly and sanitization contracts", {
  assay_one <- data.frame(
    metabolite_id = c("M1", "M2", "M2"),
    "Sample A" = c(10, 11, 12),
    "Sample-B" = c(20, 21, 22),
    Annotation = c("A", "B", "B"),
    check.names = FALSE
  )
  assay_two <- data.frame(
    feature_id = c("N1", "N2"),
    "Sample A" = c(5, 6),
    "Sample-B" = c(7, 8),
    check.names = FALSE
  )
  fixed_time <- as.POSIXct("2026-04-16 10:15:00", tz = "UTC")

  payload <- buildMetabImportWorkflowPayload(
    assay1Name = "LCMS_Pos",
    assay1Data = assay_one,
    assay2File = "assay2.tsv",
    assay2Name = "LCMS_Neg",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    metaboliteCol = "metabolite_id",
    annotationCol = "",
    sampleCols = c("Sample A", "Sample-B"),
    sanitizeNames = TRUE,
    isPattern = "",
    assay2Importer = function(path) {
      expect_identical(path, "assay2.tsv")
      list(data = assay_two)
    },
    cleanNamesFn = function(x) {
      expect_identical(x, c("Sample A", "Sample-B"))
      c("sample_a", "sample_b")
    },
    mapAssaysFn = lapply,
    timestampFn = function() fixed_time
  )

  expect_true(payload$sampleNamesSanitized)
  expect_identical(payload$dataFormat, "custom")
  expect_identical(payload$sampleCols, c("sample_a", "sample_b"))
  expect_identical(payload$columnMapping$metabolite_id_col, "metabolite_id")
  expect_null(payload$columnMapping$annotation_col)
  expect_identical(payload$columnMapping$sample_columns, c("sample_a", "sample_b"))
  expect_true(is.na(payload$columnMapping$is_pattern))
  expect_identical(names(payload$assayList), c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(
    names(payload$assayList$LCMS_Pos),
    c("metabolite_id", "sample_a", "sample_b", "Annotation")
  )
  expect_identical(
    names(payload$assayList$LCMS_Neg),
    c("feature_id", "sample_a", "sample_b")
  )
  expect_identical(payload$processingLog$timestamp, fixed_time)
  expect_identical(payload$processingLog$n_assays, 2L)
  expect_identical(payload$processingLog$assay_names, c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(payload$processingLog$detected_format, "custom")
  expect_identical(unname(payload$processingLog$n_metabolites), c(2L, 2L))
  expect_identical(payload$processingLog$n_samples, 2L)
})

test_that("metabolomics import status builders preserve key display branches", {
  detected <- htmltools::renderTags(
    buildMetabImportFormatDetectionStatus("msdial", 0.88)
  )$html
  sample_cols <- buildMetabImportSampleColumnsDisplay(
    list(sample_columns = c("Sample_A", "Sample_B"))
  )
  available_cols <- buildMetabImportAvailableColumnsDisplay(c("Peak ID", "Annotation"))
  summary_html <- htmltools::renderTags(
    buildMetabImportValidationSummary(
      assayData = data.frame(
        metabolite_id = c("M1", "M2"),
        Sample_A = c(10, 11),
        Sample_B = c(20, 21),
        check.names = FALSE
      ),
      getMetaboliteIdColFn = function() "metabolite_id",
      getSampleColumnsFn = function() c("Sample_A", "Sample_B")
    )
  )$html

  failure_html <- htmltools::renderTags(
    buildMetabImportValidationSummary(
      assayData = data.frame(
        metabolite_id = c("M1", "M2"),
        check.names = FALSE
      ),
      getMetaboliteIdColFn = function() "metabolite_id",
      getSampleColumnsFn = function() character(0)
    )
  )$html

  expect_match(detected, "MS-DIAL", fixed = TRUE)
  expect_match(detected, "Confidence: 88%", fixed = TRUE)
  expect_identical(sample_cols, "Sample_A, Sample_B")
  expect_identical(available_cols, "Peak ID, Annotation")
  expect_match(summary_html, "Validation Passed", fixed = TRUE)
  expect_match(summary_html, "Metabolites: 2", fixed = TRUE)
  expect_match(summary_html, "Samples: 2", fixed = TRUE)
  expect_match(failure_html, "Validation Failed", fixed = TRUE)
})
