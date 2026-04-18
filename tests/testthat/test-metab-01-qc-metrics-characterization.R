library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "func_metab_import.R"), local = environment())

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
    file.path(repo_root, "R", "func_metab_qc_metrics_helpers.R"),
    file.path(repo_root, "R", "func_metab_qc.R")
  ),
  symbols = c(
    "countUniqueMetabolites",
    "countMetabolitesPerSample",
    "calculateMissingness",
    "calculateSumIntensityPerSample",
    "calculateTotalUniqueMetabolitesAcrossAssays"
  ),
  env = environment()
)

test_that("metabolomics QC metric helpers preserve the current count and summary contracts", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2", "M2", NA),
    annotation_score = c(10, 20, 30, 40),
    Sample1 = c(100, 0, 50, NA),
    Sample2 = c(200, NA, 0, 25),
    Sample3 = c(0, 10, 60, 0),
    check.names = FALSE
  )

  counts <- countMetabolitesPerSample(
    assay_data = assay_data,
    sample_id_col = "Run",
    metabolite_id_col = "metabolite_id",
    sample_columns = c("Sample1", "Sample2", "Sample3")
  )
  counts <- counts[order(counts$Run), , drop = FALSE]
  expect_equal(counts$Run, c("Sample1", "Sample2", "Sample3"))
  expect_equal(counts$n_detected, c(2L, 2L, 2L))

  missingness <- calculateMissingness(
    assay_data = assay_data,
    sample_id_col = "Run",
    sample_columns = c("Sample1", "Sample2", "Sample3")
  )
  expect_equal(missingness, 6 / 12 * 100)

  summed <- calculateSumIntensityPerSample(
    assay_data = assay_data,
    sample_id_col = "Run",
    sample_columns = c("Sample1", "Sample2", "Sample3")
  )
  summed <- summed[order(summed$Run), , drop = FALSE]
  expect_equal(summed$Run, c("Sample1", "Sample2", "Sample3"))
  expect_equal(summed$sum_intensity, c(150, 225, 70))

  unique_count <- countUniqueMetabolites(
    assay_data = assay_data,
    metabolite_id_col = "metabolite_id"
  )
  expect_equal(unique_count, 2L)
})

test_that("metabolomics QC metric helpers return stable fallback values for missing inputs", {
  assay_data <- data.frame(Sample1 = c(1, 2), Sample2 = c(3, 4))

  expect_warning(
    result <- countUniqueMetabolites(assay_data, "metabolite_id"),
    "Metabolite ID column"
  )
  expect_equal(result, 0)

  expect_warning(
    missingness <- calculateMissingness(
      assay_data = assay_data,
      sample_id_col = "Run",
      sample_columns = c("MissingSample")
    ),
    "None of the provided sample_columns exist"
  )
  expect_equal(missingness, 0)

  summed <- calculateSumIntensityPerSample(
    assay_data = assay_data,
    sample_id_col = "Run",
    sample_columns = character()
  )
  expect_equal(summed$Run, c("Sample1", "Sample2"))
  expect_equal(summed$sum_intensity, c(3, 7))
})

test_that("metabolomics QC metric helpers deduplicate metabolite IDs across assays", {
  assay_list <- list(
    Assay1 = data.frame(metabolite_id = c("M1", "M2", "M2", NA)),
    Assay2 = data.frame(metabolite_id = c("M2", "M3", "M4", NA)),
    Assay3 = data.frame(other_id = c("X1", "X2"))
  )

  total_unique <- calculateTotalUniqueMetabolitesAcrossAssays(
    assay_list = assay_list,
    metabolite_id_col = "metabolite_id"
  )

  expect_equal(total_unique, 4L)
})
