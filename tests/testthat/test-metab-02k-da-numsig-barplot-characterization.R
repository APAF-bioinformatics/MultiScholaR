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

if (!methods::isGeneric("plotNumSigDiffExpBarPlot")) {
  methods::setGeneric(
    "plotNumSigDiffExpBarPlot",
    function(objectsList) {
      standardGeneric("plotNumSigDiffExpBarPlot")
    }
  )
}

printCountDaGenesTable <- function(list_of_da_tables, list_of_descriptions, formula_string) {
  expect_true(is.na(formula_string))

  list(
    plot = list(
      descriptions = list_of_descriptions,
      counts = vapply(list_of_da_tables, nrow, integer(1))
    ),
    table = data.frame(
      comparison = list_of_descriptions,
      significant = vapply(list_of_da_tables, nrow, integer(1)),
      stringsAsFactors = FALSE
    )
  )
}

loadSelectedExpressions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R")
  ),
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults") ||
      isTargetSetMethod(expr, "plotNumSigDiffExpBarPlot")
  },
  env = environment()
)

test_that("metabolomics DA numsig barplot method preserves slot updates and list names", {
  contrasts_results_table <- list(
    "groupB-groupA = B vs A" = data.frame(logFC = c(1.2, -0.8)),
    "groupC-groupA = C vs A" = data.frame(logFC = c(0.5))
  )

  result_object <- methods::new(
    "MetabolomicsDifferentialAbundanceResults",
    theObject = methods::new(
      "MetaboliteAssayData",
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id"
    ),
    contrasts_results_table = contrasts_results_table
  )

  output <- plotNumSigDiffExpBarPlot(list(primary = result_object))

  expect_equal(names(output), "primary")
  expect_identical(
    output$primary@num_sig_diff_exp_bar_plot,
    list(
      descriptions = names(contrasts_results_table),
      counts = stats::setNames(c(2L, 1L), names(contrasts_results_table))
    )
  )
  expect_identical(
    output$primary@num_sig_diff_table,
    structure(
      data.frame(
        comparison = names(contrasts_results_table),
        significant = c(2L, 1L),
        stringsAsFactors = FALSE
      ),
      row.names = names(contrasts_results_table)
    )
  )
})
