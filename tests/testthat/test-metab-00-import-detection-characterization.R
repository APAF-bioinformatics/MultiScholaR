library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_general_helpers.R"),
    file.path(repo_root, "R", "func_metab_import_detection.R"),
    file.path(repo_root, "R", "func_metab_import.R")
  ),
  symbols = c(
    "findMatchingColumn",
    "detectMetabolomicsFormat",
    "findMetabMatchingColumn",
    "getMetabolomicsColumnDefaults",
    "validateColumnMapping",
    "validateMetabColumnMapping"
  ),
  env = environment()
)

test_that("metabolomics import detection helpers preserve format scoring and column defaults", {
  headers <- c(
    "Peak ID",
    "Name",
    "RT (min)",
    "Precursor m/z",
    "Adduct",
    "Formula",
    "Ontology",
    "Total Score",
    "Sample_1"
  )

  detected <- detectMetabolomicsFormat(headers, filename = "study_msdial_height.tsv")
  defaults <- getMetabolomicsColumnDefaults("msdial")

  expect_identical(detected$format, "msdial")
  expect_equal(detected$confidence, 0.95, tolerance = 1e-9)
  expect_gt(detected$all_scores[["msdial"]], detected$all_scores[["xcms"]])
  expect_identical(findMetabMatchingColumn(headers, defaults$metabolite_id), "Peak ID")
  expect_identical(findMetabMatchingColumn(headers, defaults$annotation), "Name")
  expect_identical(findMetabMatchingColumn(headers, defaults$rt), "RT (min)")
  expect_identical(defaults$is_pattern, "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:")

  custom_defaults <- getMetabolomicsColumnDefaults("custom")
  expect_null(custom_defaults$metabolite_id)
  expect_true(is.na(custom_defaults$is_pattern))

  unknown <- detectMetabolomicsFormat(
    headers = c("feature_label", "sample_a", "sample_b"),
    filename = "custom_export.csv"
  )
  expect_identical(unknown$format, "unknown")
})

test_that("metabolomics import validation helpers preserve current summary and alias contracts", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2", "M2", NA),
    Sample1 = c(100, 0, 50, NA),
    Sample2 = c(200, NA, 0, 25),
    check.names = FALSE
  )

  validated <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = "metabolite_id",
    sample_columns = c("Sample1", "Sample2")
  )

  aliased <- validateColumnMapping(
    data = assay_data,
    metabolite_id_column = "metabolite_id",
    sample_columns = c("Sample1", "Sample2")
  )

  expect_true(validated$valid)
  expect_identical(aliased, validated)
  expect_identical(validated$summary$n_metabolites, 3L)
  expect_identical(validated$summary$n_samples, 2L)
  expect_equal(validated$summary$pct_missing, 50)
  expect_match(validated$warnings[[1]], "duplicate metabolite IDs")

  missing_cols <- validateMetabColumnMapping(
    data = assay_data,
    metabolite_id_column = "missing_id",
    sample_columns = c("MissingSample")
  )

  expect_false(missing_cols$valid)
  expect_match(missing_cols$errors[[1]], "Metabolite ID column not found")
  expect_match(missing_cols$errors[[2]], "Sample columns not found")
})
