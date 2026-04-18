library(testthat)
library(dplyr)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) && length(expr) >= 3 && as.character(expr[[1]]) %in% c("<-", "=")
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
    file.path(repo_root, "R", "func_metab_qc_filtering_helpers.R"),
    file.path(repo_root, "R", "func_metab_qc.R")
  ),
  symbols = c(
    "metaboliteIntensityFilteringHelper",
    "resolveDuplicateFeaturesByIntensity"
  ),
  env = environment()
)

test_that("metabolomics QC filtering helpers preserve current row filtering behavior", {
  assay_table <- data.frame(
    metabolite_id = c("M1", "M2", "M3"),
    annotation_score = c(10, 20, 30),
    Sample1 = c(100, 20, 55),
    Sample2 = c(200, 30, NA),
    Sample3 = c(10, 40, 60),
    check.names = FALSE
  )

  filtered <- metaboliteIntensityFilteringHelper(
    assay_table = assay_table,
    min_metabolite_intensity_threshold = 50,
    metabolites_proportion_of_samples_below_cutoff = 0.5,
    metabolite_id_column = "metabolite_id"
  )

  expect_equal(filtered$metabolite_id, "M3")
  expect_equal(filtered$annotation_score, 30)
})

test_that("metabolomics QC filtering helpers return stable fallbacks for degenerate inputs", {
  no_numeric <- data.frame(
    metabolite_id = c("M1", "M2"),
    label = c("A", "B"),
    check.names = FALSE
  )

  expect_warning(
    passthrough <- metaboliteIntensityFilteringHelper(
      assay_table = no_numeric,
      min_metabolite_intensity_threshold = 50,
      metabolites_proportion_of_samples_below_cutoff = 0.5,
      metabolite_id_column = "metabolite_id"
    ),
    "No numeric sample columns found"
  )
  expect_identical(passthrough, no_numeric)

  duplicate_free <- data.frame(
    metabolite_id = c("M1", "M2"),
    Sample1 = c(10, 20),
    Sample2 = c(30, 40),
    check.names = FALSE
  )

  expect_identical(
    resolveDuplicateFeaturesByIntensity(duplicate_free, "metabolite_id", c("Sample1", "Sample2")),
    duplicate_free
  )
})

test_that("metabolomics QC duplicate resolution keeps the highest-intensity feature per metabolite", {
  assay_tibble <- data.frame(
    metabolite_id = c("M1", "M1", "M2"),
    Sample1 = c(10, 30, 5),
    Sample2 = c(20, 40, 15),
    annotation_score = c(1, 2, 3),
    check.names = FALSE
  )

  resolved <- resolveDuplicateFeaturesByIntensity(
    assay_tibble = assay_tibble,
    id_col = "metabolite_id",
    sample_cols = c("Sample1", "Sample2")
  )
  resolved <- resolved[order(resolved$metabolite_id), , drop = FALSE]

  expect_equal(resolved$metabolite_id, c("M1", "M2"))
  expect_equal(resolved$annotation_score, c(2, 3))
  expect_equal(resolved$Sample1, c(30, 5))
  expect_equal(resolved$Sample2, c(40, 15))
})
