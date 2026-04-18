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

findSelectedExpression <- function(paths, matcher) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        return(expr)
      }
    }
  }

  NULL
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
    slots = c(label = "character"),
    prototype = list(label = character())
  )
}

if (!methods::isGeneric("differentialAbundanceAnalysis")) {
  methods::setGeneric(
    "differentialAbundanceAnalysis",
    function(
      theObject,
      contrasts_tbl = NULL,
      formula_string = NULL,
      group_id = NULL,
      da_q_val_thresh = NULL,
      treat_lfc_cutoff = NULL,
      eBayes_trend = NULL,
      eBayes_robust = NULL,
      args_group_pattern = NULL,
      args_row_id = NULL,
      qvalue_column = NULL,
      raw_pvalue_column = NULL
    ) {
      standardGeneric("differentialAbundanceAnalysis")
    }
  )
}

helper_calls <- list()
differentialAbundanceAnalysisHelper <- function(
  obj,
  contrasts_tbl = NULL,
  formula_string = NULL,
  group_id = NULL,
  da_q_val_thresh = NULL,
  treat_lfc_cutoff = NULL,
  eBayes_trend = NULL,
  eBayes_robust = NULL,
  args_group_pattern = NULL
) {
  helper_calls[[length(helper_calls) + 1]] <<- list(
    label = obj@label,
    contrasts_tbl = contrasts_tbl,
    formula_string = formula_string,
    group_id = group_id,
    da_q_val_thresh = da_q_val_thresh,
    treat_lfc_cutoff = treat_lfc_cutoff,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust,
    args_group_pattern = args_group_pattern
  )

  list(label = obj@label, forwarded_group = group_id)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_da_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "differentialAbundanceAnalysis")
  },
  env = environment()
)

test_that("metabolomics S4 DA analysis list method validates inputs and preserves names", {
  object_list <- list(
    primary = methods::new("MetaboliteAssayData", label = "alpha"),
    secondary = methods::new("MetaboliteAssayData", label = "beta")
  )
  helper_calls <<- list()

  output <- differentialAbundanceAnalysis(
    object_list,
    contrasts_tbl = data.frame(contrasts = "B-A", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group",
    group_id = "group",
    da_q_val_thresh = 0.1,
    treat_lfc_cutoff = 0.5,
    eBayes_trend = FALSE,
    eBayes_robust = TRUE,
    args_group_pattern = "(A|B)",
    args_row_id = "ignored_here",
    qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue"
  )

  expect_identical(names(output), c("primary", "secondary"))
  expect_length(helper_calls, 2)
  expect_identical(vapply(helper_calls, `[[`, character(1), "label"), c("alpha", "beta"))
  expect_equal(helper_calls[[1]]$contrasts_tbl$contrasts, "B-A")
  expect_identical(helper_calls[[1]]$formula_string, "~ 0 + group")
  expect_identical(helper_calls[[1]]$group_id, "group")
  expect_identical(helper_calls[[1]]$da_q_val_thresh, 0.1)
  expect_identical(helper_calls[[1]]$treat_lfc_cutoff, 0.5)
  expect_identical(helper_calls[[1]]$eBayes_trend, FALSE)
  expect_identical(helper_calls[[1]]$eBayes_robust, TRUE)
  expect_identical(helper_calls[[1]]$args_group_pattern, "(A|B)")
  expect_identical(output$primary$forwarded_group, "group")
  expect_identical(output$secondary$label, "beta")
})

test_that("metabolomics S4 DA analysis list method rejects non-MetaboliteAssayData entries", {
  expect_error(
    differentialAbundanceAnalysis(
      list(methods::new("MetaboliteAssayData", label = "alpha"), "bad-entry")
    ),
    "All objects in objectsList must be of class MetaboliteAssayData"
  )
})

test_that("metabolomics S4 DA analysis definition still delegates through the helper", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "differentialAbundanceAnalysis")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "differentialAbundanceAnalysisHelper\\(obj,")
  expect_match(target_text, "names\\(results_list\\) <- names\\(objectsList\\)")
  expect_match(target_text, "All objects in objectsList must be of class MetaboliteAssayData", fixed = TRUE)
})
