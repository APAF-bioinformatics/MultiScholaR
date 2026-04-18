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
    file.path(repo_root, "R", "func_metab_da_core_engine.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "runTestsContrastsMetabDA",
  env = environment()
)

test_that("metabolomics DA core engine preserves the no-common-samples guard", {
  assay_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("Sample_1", "Sample_2"))
  )

  design_matrix <- data.frame(
    sample_id = c("Other_1", "Other_2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    runTestsContrastsMetabDA(
      data = assay_matrix,
      contrast_strings = "groupB-groupA",
      design_matrix = design_matrix,
      formula_string = "~ 0 + group",
      sample_id_col = "sample_id"
    ),
    "No common samples between data matrix and design matrix"
  )
})

test_that("metabolomics DA core engine preserves contrast-level validation", {
  assay_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("Sample_1", "Sample_2"))
  )

  design_matrix <- data.frame(
    sample_id = c("Sample_1", "Sample_2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    runTestsContrastsMetabDA(
      data = assay_matrix,
      contrast_strings = "Treatment-Control",
      design_matrix = design_matrix,
      formula_string = "~ 0 + group",
      sample_id_col = "sample_id"
    ),
    "references undefined levels"
  )
})
