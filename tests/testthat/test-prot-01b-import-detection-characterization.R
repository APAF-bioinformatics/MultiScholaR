test_that("detectProteomicsFormat boosts DIA-NN confidence for parquet inputs", {
  diann_headers <- c(
    "Protein.Group", "Protein.Ids", "Protein.Names",
    "Precursor.Id", "Modified.Sequence", "Stripped.Sequence",
    "Precursor.Charge", "Q.Value", "PG.Q.Value", "Run"
  )

  tsv_result <- detectProteomicsFormat(diann_headers, "report.tsv")
  parquet_result <- detectProteomicsFormat(diann_headers, "report.parquet")

  expect_equal(tsv_result$format, "diann")
  expect_equal(parquet_result$format, "diann")
  expect_gt(parquet_result$confidence, tsv_result$confidence)
})

test_that("detectProteomicsFormat returns unknown for weak unmatched headers", {
  weak_headers <- c("sample", "value", "description")

  result <- detectProteomicsFormat(weak_headers, "random.csv")

  expect_equal(result$format, "unknown")
  expect_lt(result$confidence, 0.3)
})

test_that("detectProteomicsFormat uses FragPipe filename bonus without misclassifying weak input", {
  fragpipe_like_headers <- c(
    "Protein ID", "Protein", "Gene", "Description",
    "Spectral Count", "Sample1 Intensity", "Sample2 Intensity"
  )

  neutral_result <- detectProteomicsFormat(fragpipe_like_headers, "protein.tsv")
  fragpipe_named_result <- detectProteomicsFormat(fragpipe_like_headers, "msfragger_protein.tsv")

  expect_equal(neutral_result$format, "fragpipe")
  expect_equal(fragpipe_named_result$format, "fragpipe")
  expect_gt(fragpipe_named_result$confidence, neutral_result$confidence)
})

# APAF Bioinformatics | test-prot-01b-import-detection-characterization.R | Approved | 2026-04-11
