library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedExpressions <- function(paths, matcher, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        eval(expr, envir = env)
      }
    }
  }
}

normalizeSelectorValue <- function(x) {
  if (is.character(x)) {
    return(x[[1]])
  }

  if (is.symbol(x)) {
    return(as.character(x))
  }

  as.character(x)
}

isTargetSetClass <- function(expr, class_name) {
  is.call(expr) &&
    identical(as.character(expr[[1]]), "setClass") &&
    length(expr) >= 2 &&
    identical(normalizeSelectorValue(expr[[2]]), class_name)
}

isTargetSetMethod <- function(expr, method_name) {
  if (!is.call(expr) || !identical(as.character(expr[[1]]), "setMethod")) {
    return(FALSE)
  }

  expr_parts <- as.list(expr)
  part_names <- names(expr_parts)
  method_arg <- expr_parts[[which(part_names == "f")[1]]]

  !is.null(method_arg) && identical(normalizeSelectorValue(method_arg), method_name)
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character"
    )
  )
}

if (!methods::isGeneric("getDaResultsLongFormat")) {
  methods::setGeneric(
    "getDaResultsLongFormat",
    function(objectsList) {
      standardGeneric("getDaResultsLongFormat")
    }
  )
}

join_by <- dplyr::join_by
sym <- rlang::sym
str_split_i <- stringr::str_split_i

loadSelectedExpressions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R")
  ),
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults") ||
      isTargetSetMethod(expr, "getDaResultsLongFormat")
  },
  env = environment()
)

newDaLongResultObject <- function(metabolite_data) {
  methods::new(
    "MetabolomicsDifferentialAbundanceResults",
    theObject = methods::new(
      "MetaboliteAssayData",
      metabolite_data = metabolite_data,
      metabolite_id_column = "metabolite_id"
    ),
    contrasts_results_table = list(
      "groupB-groupA = B vs A" = data.frame(
        metabolite_id = c("M1", "M2"),
        fdr_qvalue = c(0.20, 0.01),
        raw_pvalue = c(0.40, 0.02),
        logFC = c(1.5, -0.5),
        stringsAsFactors = FALSE
      ),
      "groupC-groupA = C vs A" = data.frame(
        metabolite_id = c("M1", "M2"),
        fdr_qvalue = c(0.05, 0.50),
        raw_pvalue = c(0.06, 0.70),
        logFC = c(0.7, 0.1),
        stringsAsFactors = FALSE
      )
    )
  )
}

test_that("metabolomics DA long-format method preserves joined counts for each contrast row", {
  counts_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    intensity = c(100, 200),
    annotation = c("alpha", "beta"),
    stringsAsFactors = FALSE
  )

  output <- getDaResultsLongFormat(list(primary = newDaLongResultObject(list(counts_table))))

  expect_identical(names(output), "primary")
  expect_s3_class(output$primary@results_table_long, "data.frame")
  expect_identical(
    names(output$primary@results_table_long),
    c(
      "comparison",
      "metabolite_id",
      "fdr_qvalue",
      "raw_pvalue",
      "logFC",
      "comparision_short",
      "intensity",
      "annotation"
    )
  )
  expect_equal(
    output$primary@results_table_long$comparison,
    c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A",
      "groupC-groupA = C vs A"
    )
  )
  expect_equal(output$primary@results_table_long$metabolite_id, c("M1", "M2", "M1", "M2"))
  expect_equal(output$primary@results_table_long$comparision_short, c(
    "groupB-groupA ",
    "groupB-groupA ",
    "groupC-groupA ",
    "groupC-groupA "
  ))
  expect_equal(output$primary@results_table_long$intensity, c(100, 200, 100, 200))
  expect_equal(output$primary@results_table_long$annotation, c("alpha", "beta", "alpha", "beta"))
})

test_that("metabolomics DA long-format method currently errors on a bare-data-frame counts slot", {
  counts_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    intensity = c(100, 200),
    stringsAsFactors = FALSE
  )

  expect_error(
    getDaResultsLongFormat(list(primary = newDaLongResultObject(counts_table))),
    "must share the same src"
  )
})
