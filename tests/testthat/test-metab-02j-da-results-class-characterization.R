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

isTargetSetClass <- function(expr, class_name) {
  is.call(expr) &&
    identical(as.character(expr[[1]]), "setClass") &&
    length(expr) >= 2 &&
    identical(as.character(expr[[2]]), class_name)
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(metabolite_data = "list")
  )
}

loadSelectedExpressions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R")
  ),
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults")
  },
  env = environment()
)

test_that("metabolomics DA results S4 class preserves slot layout and defaults", {
  class_def <- methods::getClassDef("MetabolomicsDifferentialAbundanceResults")

  expect_false(is.null(class_def))
  expect_identical(
    methods::slotNames(class_def),
    c(
      "theObject",
      "fit.eb",
      "contrasts_results_table",
      "num_sig_diff_exp_bar_plot",
      "num_sig_diff_table",
      "volcano_plot",
      "interactive_volcano_plot",
      "p_value_dist_plot",
      "results_table_long",
      "results_table_wide"
    )
  )

  result_object <- methods::new(
    "MetabolomicsDifferentialAbundanceResults",
    theObject = methods::new("MetaboliteAssayData", metabolite_data = list())
  )

  expect_s4_class(result_object, "MetabolomicsDifferentialAbundanceResults")
  expect_type(result_object@fit.eb, "NULL")
  expect_identical(result_object@contrasts_results_table, list())
  expect_identical(result_object@num_sig_diff_exp_bar_plot, list())
  expect_s3_class(result_object@num_sig_diff_table, "data.frame")
  expect_identical(nrow(result_object@num_sig_diff_table), 0L)
  expect_identical(result_object@volcano_plot, list())
  expect_identical(result_object@interactive_volcano_plot, list())
  expect_identical(result_object@p_value_dist_plot, list())
  expect_s3_class(result_object@results_table_long, "data.frame")
  expect_identical(nrow(result_object@results_table_long), 0L)
  expect_s3_class(result_object@results_table_wide, "data.frame")
  expect_identical(nrow(result_object@results_table_wide), 0L)
})
