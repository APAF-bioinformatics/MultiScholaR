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
    file.path(repo_root, "R", "func_metab_da_quant_data.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "getMetaboliteQuantData",
  env = environment()
)

test_that("metabolomics DA quant-data helper preserves default metadata filtering", {
  assay_df <- data.frame(
    `Alignment ID` = c("A1", "A2"),
    `Metabolite name` = c("Met One", "Met Two"),
    Sample_1 = c(10, 11),
    Sample_2 = c(20, 21),
    batch = c("B1", "B2"),
    check.names = FALSE
  )

  quant_info <- getMetaboliteQuantData(assay_df)

  expect_identical(quant_info$sample_cols, c("Sample_1", "Sample_2"))
  expect_identical(
    colnames(quant_info$quant_data),
    c("Sample_1", "Sample_2")
  )
  expect_identical(
    colnames(quant_info$meta_data),
    c("Alignment ID", "Metabolite name")
  )
})

test_that("metabolomics DA quant-data helper preserves explicit metadata exclusions", {
  assay_df <- data.frame(
    feature_id = c("A1", "A2"),
    retention_time = c(1.2, 1.5),
    Sample_1 = c(10, 11),
    Sample_2 = c(20, 21),
    label = c("Met One", "Met Two"),
    qc_flag = c(TRUE, FALSE),
    check.names = FALSE
  )

  quant_info <- getMetaboliteQuantData(
    assay_df = assay_df,
    metabolite_id_col = "feature_id",
    annotation_col = "label",
    additional_meta_cols = c("retention_time", "missing_column")
  )

  expect_identical(quant_info$sample_cols, c("Sample_1", "Sample_2"))
  expect_identical(
    colnames(quant_info$meta_data),
    c("feature_id", "label", "retention_time")
  )
  expect_false("retention_time" %in% colnames(quant_info$quant_data))
})
