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
    file.path(repo_root, "R", "func_metab_da_results_long_format.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "createMetabDaResultsLongFormat",
  env = environment()
)

test_that("metabolomics DA results long-format helper preserves intensity widening and ordering", {
  lfc_qval_tbl <- data.frame(
    metabolite_id = c("M1", "M2"),
    comparison = c("groupB-groupA", "groupB-groupA"),
    logFC = c(1.25, -0.5),
    fdr_qvalue = c(0.01, 0.20),
    stringsAsFactors = FALSE
  )

  expr_matrix <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("S1", "S2"))
  )

  design_matrix <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("groupA", "groupB"),
    stringsAsFactors = FALSE
  )

  result <- createMetabDaResultsLongFormat(
    lfc_qval_tbl = lfc_qval_tbl,
    expr_matrix = expr_matrix,
    design_matrix = design_matrix,
    sample_id_col = "sample_id",
    group_id_col = "group",
    metabolite_id_col = "metabolite_id"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(result$metabolite_id, c("M1", "M2"))
  expect_equal(result$numerator, rep("groupB", 2))
  expect_equal(result$denominator, rep("groupA", 2))
  expect_equal(result$intensity.S1.groupA, c(10, 20))
  expect_equal(result$intensity.S2.groupB, c(30, 40))
})

test_that("metabolomics DA results long-format helper keeps fallback numerator and denominator columns", {
  lfc_qval_tbl <- data.frame(
    metabolite_id = c("M1", "M2"),
    comparison = c("no_delimiter", "no_delimiter"),
    logFC = c(1.25, -0.5),
    fdr_qvalue = c(0.01, 0.20),
    stringsAsFactors = FALSE
  )

  expr_matrix <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("M1", "M2"), c("S1", "S2"))
  )

  design_matrix <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("groupA", "groupB"),
    stringsAsFactors = FALSE
  )

  expect_output(
    result <- createMetabDaResultsLongFormat(
      lfc_qval_tbl = lfc_qval_tbl,
      expr_matrix = expr_matrix,
      design_matrix = design_matrix,
      sample_id_col = "sample_id",
      group_id_col = "group",
      metabolite_id_col = "metabolite_id"
    ),
    "Could not parse contrast groups"
  )

  expect_true(all(is.na(result$numerator)))
  expect_true(all(is.na(result$denominator)))
  expect_equal(result$metabolite_id, c("M1", "M2"))
  expect_equal(result$intensity.S1.groupA, c(10, 20))
  expect_equal(result$intensity.S2.groupB, c(30, 40))
})
