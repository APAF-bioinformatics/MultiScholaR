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

if (!methods::isGeneric("getDaResultsWideFormat")) {
  methods::setGeneric(
    "getDaResultsWideFormat",
    function(
      objectsList,
      qvalue_column = "fdr_qvalue",
      raw_pvalue_column = "raw_pvalue",
      log2fc_column = "logFC"
    ) {
      standardGeneric("getDaResultsWideFormat")
    }
  )
}

sym <- rlang::sym
str_split_i <- stringr::str_split_i
join_by <- dplyr::join_by
matches <- tidyselect::matches

loadSelectedExpressions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R")
  ),
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults") ||
      isTargetSetMethod(expr, "getDaResultsWideFormat")
  },
  env = environment()
)

newDaWideResultObject <- function(metabolite_data) {
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

test_that("metabolomics DA wide-format method preserves joined counts and q-value ordering", {
  counts_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    intensity = c(100, 200),
    annotation = c("alpha", "beta"),
    stringsAsFactors = FALSE
  )

  output <- getDaResultsWideFormat(list(primary = newDaWideResultObject(list(counts_table))))

  expect_identical(names(output), "primary")
  expect_s3_class(output$primary@results_table_wide, "data.frame")
  expect_equal(output$primary@results_table_wide$metabolite_id, c("M2", "M1"))
  expect_equal(output$primary@results_table_wide$intensity, c(200, 100))
  expect_equal(output$primary@results_table_wide$annotation, c("beta", "alpha"))
  expect_identical(
    names(output$primary@results_table_wide),
    c(
      "metabolite_id",
      "logFC:groupB-groupA ",
      "logFC:groupC-groupA ",
      "fdr_qvalue:groupB-groupA ",
      "fdr_qvalue:groupC-groupA ",
      "raw_pvalue:groupB-groupA ",
      "raw_pvalue:groupC-groupA ",
      "intensity",
      "annotation"
    )
  )
  expect_equal(
    output$primary@results_table_wide[["fdr_qvalue:groupB-groupA "]],
    c(0.01, 0.20)
  )
})

test_that("metabolomics DA wide-format method currently errors on a bare-data-frame counts slot", {
  counts_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    intensity = c(100, 200),
    stringsAsFactors = FALSE
  )

  expect_error(
    getDaResultsWideFormat(list(primary = newDaWideResultObject(counts_table))),
    "must share the same src"
  )
})
