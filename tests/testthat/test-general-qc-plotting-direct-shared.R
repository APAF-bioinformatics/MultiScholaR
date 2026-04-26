# fidelity-coverage-compare: shared
library(testthat)

test_that("general QC plotting helpers preserve current plot assembly", {
  protein_long <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2", "P3", "P3"),
    Stripped.Sequence = c("pep1", "pep1", "pep2", "pep2", "pep3", "pep3"),
    Run = c("S1", "S2", "S1", "S2", "S1", "S2"),
    Log2.Protein.Imputed = c(10, 12, NA, 9, 8, 11),
    Peptide.RawQuantity = c(100, 120, 80, NA, 60, 65),
    percent_missing = c(0, 0, 50, 50, 0, 0),
    stringsAsFactors = FALSE
  )

  peptide_counts <- data.frame(
    Protein.Ids = c("P1", "P1", "P2", "P2", "P3", "P3"),
    Run = c("S1", "S2", "S1", "S2", "S1", "S2"),
    num_peptides_after_impute = c(1, 2, 1, 2, 2, 2),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    plotDensityOfProteinIntensityPerSample(
      protein_intensity_long_tbl = protein_long,
      number_of_peptides_per_protein_per_sample = peptide_counts
    ),
    "ggplot"
  )

  expect_s3_class(
    plotPercentSamplesVsProteinQuantified(
      protein_intensity_long_tbl = protein_long,
      number_of_peptides_per_protein_per_sample = peptide_counts
    ),
    "ggplot"
  )

  expect_s3_class(
    plotPeptidesProteinsCountsPerSampleHelper(protein_long),
    "ggplot"
  )

  expect_s3_class(
    plotHistogramOfPercentMissingPerIndvidual(
      unique(protein_long[c("Run", "percent_missing")])
    ),
    "ggplot"
  )
})

test_that("general QC utility plots preserve current missingness and RLE helpers", {
  input_matrix <- data.frame(
    S1 = c(1, 2, NA),
    S2 = c(2, NaN, 4),
    S3 = c(3, 4, Inf),
    check.names = FALSE
  )

  expect_s3_class(plotNumMissingValues(input_matrix), "ggplot")
  expect_s3_class(plotNumOfValues(input_matrix), "ggplot")
  expect_s3_class(plotNumOfValuesNoLog(input_matrix), "ggplot")

  rle_data <- getOneRlePlotData(as.matrix(input_matrix))

  expect_true(all(c("Quantiles", "name", "value") %in% names(rle_data)))
  expect_s3_class(plotRleQc(rle_data), "ggplot")
})
