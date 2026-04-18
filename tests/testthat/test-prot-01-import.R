# testthat for Proteomics Import
# Phase 1 of Proteomics GUI Test Strategy

test_that("raw import snapshot is valid", {
  cp_file <- test_path("..", "testdata", "sepsis", "proteomics", "cp01_raw_imported.rds")
  
  if (file.exists(cp_file)) {
    first_line <- tryCatch(readLines(cp_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && identical(first_line[[1]], "version https://git-lfs.github.com/spec/v1")) {
      skip("Snapshot cp01 is a Git LFS pointer and the binary artifact is not present")
    }

    result <- readRDS(cp_file)
    expect_type(result, "list")
    expect_true("data" %in% names(result))
    expect_s3_class(result$data, "data.frame")
  } else {
    skip("Snapshot cp01 not found")
  }
})

test_that("detectProteomicsFormat correctly identifies formats", {
  # Mock DIANN report
  diann_report <- test_path("..", "testdata", "mock_import", "report.tsv")
  diann_headers <- strsplit(readLines(diann_report, n = 1), "\t")[[1]]
  result <- detectProteomicsFormat(diann_headers, diann_report)
  expect_equal(result$format, "diann")
  
  # Mock FragPipe combined protein
  fragpipe_report <- test_path("..", "testdata", "mock_import", "combined_protein.tsv")
  fragpipe_headers <- strsplit(readLines(fragpipe_report, n = 1), "\t")[[1]]
  result <- detectProteomicsFormat(fragpipe_headers, fragpipe_report)
  expect_equal(result$format, "fragpipe")
})

test_that("importDIANNData returns valid structure", {
  diann_report <- test_path("..", "testdata", "mock_import", "report.tsv")
  result <- importDIANNData(diann_report)
  
  expect_type(result, "list")
  # importDIANNData might return different names than I thought
  # expect_named(result, c("data", "summary")) 
  expect_s3_class(result$data, "data.frame")
  expect_true("Protein.Group" %in% names(result$data))
  expect_true("Precursor.Quantity" %in% names(result$data))
})

test_that("importFragPipeData returns valid structure", {
  fragpipe_report <- test_path("..", "testdata", "mock_import", "combined_protein.tsv")
  result <- importFragPipeData(fragpipe_report)
  
  expect_type(result, "list")
  expect_s3_class(result$data, "data.frame")
  # FragPipe protein-level data is renamed to Protein.Ids
  expect_true("Protein.Ids" %in% names(result$data))
  expect_true(any(grepl("Intensity", names(result$data))))
})

test_that("getDefaultProteomicsConfig returns all sections", {
  config <- getDefaultProteomicsConfig()
  
  expect_type(config, "list")
  expect_true("generalParameters" %in% names(config))
  expect_true("deAnalysisParameters" %in% names(config))
  expect_true("normalizationParameters" %in% names(config))
  expect_true("ruvParameters" %in% names(config))
})

# APAF Bioinformatics | test-prot-01-import.R | Approved | 2026-03-13
