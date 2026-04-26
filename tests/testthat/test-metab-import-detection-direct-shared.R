# fidelity-coverage-compare: shared
library(testthat)

test_that("metabolomics import detection helpers preserve vendor scoring and defaults", {
  msdial_headers <- c(
    "Peak ID", "Name", "RT (min)", "Precursor m/z", "Adduct",
    "Formula", "Ontology", "Total Score", "Sample_1"
  )
  progenesis_headers <- c(
    "Compound", "Neutral mass (Da)", "m/z", "Charge",
    "Retention time (min)", "Raw abundance"
  )
  xcms_headers <- c(
    "featureid", "mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "into"
  )
  cd_headers <- c(
    "Compound Name", "Molecular Formula", "RT [min]", "Area", "Height", "MZ"
  )

  msdial_detected <- detectMetabolomicsFormat(msdial_headers, filename = "study_msdial_height.tsv")
  progenesis_detected <- detectMetabolomicsFormat(progenesis_headers, filename = "progenesis_qi_export.csv")
  xcms_detected <- detectMetabolomicsFormat(xcms_headers, filename = "xcms_peak_list.tsv")
  cd_detected <- detectMetabolomicsFormat(cd_headers, filename = "compound_discoverer_cd_export.csv")
  unknown_detected <- detectMetabolomicsFormat(
    headers = c("feature_label", "sample_a", "sample_b"),
    filename = "custom_export.csv"
  )

  expect_identical(msdial_detected$format, "msdial")
  expect_identical(progenesis_detected$format, "progenesis")
  expect_identical(xcms_detected$format, "xcms")
  expect_identical(cd_detected$format, "compound_discoverer")
  expect_identical(unknown_detected$format, "unknown")

  msdial_defaults <- getMetabolomicsColumnDefaults("msdial")
  progenesis_defaults <- getMetabolomicsColumnDefaults("progenesis")
  xcms_defaults <- getMetabolomicsColumnDefaults("xcms")
  cd_defaults <- getMetabolomicsColumnDefaults("compound_discoverer")
  custom_defaults <- getMetabolomicsColumnDefaults("custom")
  unknown_defaults <- getMetabolomicsColumnDefaults("unknown")

  expect_identical(findMetabMatchingColumn(msdial_headers, msdial_defaults$metabolite_id), "Peak ID")
  expect_identical(findMetabMatchingColumn(msdial_headers, msdial_defaults$annotation), "Name")
  expect_identical(findMetabMatchingColumn(progenesis_headers, progenesis_defaults$rt), "Retention time (min)")
  expect_identical(findMetabMatchingColumn(xcms_headers, xcms_defaults$mz), "mz")
  expect_identical(findMetabMatchingColumn(cd_headers, cd_defaults$annotation), "Compound Name")
  expect_identical(msdial_defaults$is_pattern, "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:")
  expect_true(is.null(custom_defaults$metabolite_id))
  expect_true(is.na(custom_defaults$is_pattern))
  expect_true(is.null(unknown_defaults$annotation))
})

test_that("metabolomics import validation helpers preserve valid, warning, and error branches", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2", "M2", NA),
    Sample1 = c(100, 0, 50, NA),
    Sample2 = c(200, NA, 0, 25),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validated <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = "metabolite_id",
    sample_columns = c("Sample1", "Sample2")
  )

  high_missing <- validateMetabColumnMapping(
    data = data.frame(
      metabolite_id = c("M1", "M2"),
      Sample1 = c(NA_real_, 0),
      Sample2 = c(0, NA_real_),
      check.names = FALSE
    ),
    metabolite_id_column = "metabolite_id",
    sample_columns = c("Sample1", "Sample2")
  )

  missing_id <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = NULL,
    sample_columns = "Sample1"
  )
  no_samples <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = "metabolite_id",
    sample_columns = character()
  )
  missing_cols <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = "missing_id",
    sample_columns = c("MissingSample")
  )

  expect_true(validated$valid)
  expect_identical(validated$summary$n_metabolites, 3L)
  expect_identical(validated$summary$n_samples, 2L)
  expect_equal(validated$summary$pct_missing, 50)
  expect_match(validated$warnings[[1]], "duplicate metabolite IDs", fixed = TRUE)
  expect_match(high_missing$warnings[[1]], "High proportion of missing values", fixed = TRUE)

  expect_false(missing_id$valid)
  expect_match(missing_id$errors[[1]], "Metabolite ID column is not specified", fixed = TRUE)

  expect_false(no_samples$valid)
  expect_match(no_samples$errors[[1]], "No sample columns specified", fixed = TRUE)

  expect_false(missing_cols$valid)
  expect_match(missing_cols$errors[[1]], "Metabolite ID column not found", fixed = TRUE)
  expect_match(missing_cols$errors[[2]], "Sample columns not found", fixed = TRUE)
})
