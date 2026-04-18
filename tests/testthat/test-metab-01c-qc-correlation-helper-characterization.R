library(testthat)

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
    file.path(repo_root, "R", "func_metab_qc_correlation_helpers.R"),
    file.path(repo_root, "R", "func_metab_qc.R")
  ),
  symbols = c("calculateMetabolitePairCorrelation"),
  env = environment()
)

test_that("metabolomics QC correlation helper preserves pairwise Pearson behavior", {
  input_pair_table <- data.frame(
    metabolite_id = c("M1", "M1", "M2", "M2", "M3", "M3"),
    Run = c("Sample1", "Sample2", "Sample1", "Sample2", "Sample1", "Sample2"),
    abundance = c(1, 2, 2, 4, 3, 6),
    check.names = FALSE
  )

  correlation <- calculateMetabolitePairCorrelation(
    input_pair_table = input_pair_table,
    feature_id_column = "metabolite_id",
    sample_id_column = "Run",
    value_column = "abundance"
  )

  expect_equal(correlation, 1)
})

test_that("metabolomics QC correlation helper keeps current fallback behavior for invalid pairs", {
  single_sample_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    Run = c("Sample1", "Sample1"),
    abundance = c(10, 20),
    check.names = FALSE
  )

  expect_warning(
    invalid_pair_result <- calculateMetabolitePairCorrelation(
      input_pair_table = single_sample_table,
      feature_id_column = "metabolite_id",
      sample_id_column = "Run",
      value_column = "abundance"
    ),
    "does not contain exactly two samples"
  )
  expect_true(is.na(invalid_pair_result))

  one_feature_pair <- data.frame(
    metabolite_id = c("M1", "M1"),
    Run = c("Sample1", "Sample2"),
    abundance = c(10, 15),
    check.names = FALSE
  )

  insufficient_features <- calculateMetabolitePairCorrelation(
    input_pair_table = one_feature_pair,
    feature_id_column = "metabolite_id",
    sample_id_column = "Run",
    value_column = "abundance"
  )
  expect_true(is.na(insufficient_features))
})
