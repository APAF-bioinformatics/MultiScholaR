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

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(metabolite_data = "list")
  )
}

if (!methods::isClass("ProteinQuantitativeData")) {
  methods::setClass(
    "ProteinQuantitativeData",
    slots = c(protein_quant_table = "data.frame")
  )
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "getCountsTable",
  env = environment()
)

test_that("metabolomics DA counts-table helper preserves metabolite assay access", {
  assay_list <- list(
    LCMS_Pos = data.frame(
      metabolite_id = c("M1", "M2"),
      Sample_1 = c(10, 20),
      stringsAsFactors = FALSE
    )
  )

  theObject <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list
  )

  expect_identical(getCountsTable(theObject), assay_list)
})

test_that("metabolomics DA counts-table helper preserves protein fallback access", {
  protein_quant_table <- data.frame(
    protein_id = c("P1", "P2"),
    Sample_1 = c(5, 6),
    stringsAsFactors = FALSE
  )

  theObject <- methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = protein_quant_table
  )

  expect_identical(getCountsTable(theObject), protein_quant_table)
})

test_that("metabolomics DA counts-table helper preserves unsupported-object guard", {
  expect_error(getCountsTable(list()), "Unsupported object type")
})
