library(methods)
library(testthat)
library(ggplot2)

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

  if (is.null(method_arg) && length(expr_parts) >= 2) {
    method_arg <- expr_parts[[2]]
  }

  !is.null(method_arg) && identical(normalizeSelectorValue(method_arg), method_name)
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(args = "list"),
    prototype = list(args = list())
  )
}

if (!methods::isGeneric("plotPearson")) {
  methods::setGeneric(
    "plotPearson",
    function(theObject, tech_rep_remove_regex = "pool", correlation_group = NA) {
      standardGeneric("plotPearson")
    }
  )
}

if (!methods::isGeneric("pearsonCorForSamplePairs")) {
  methods::setGeneric(
    "pearsonCorForSamplePairs",
    function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
      standardGeneric("pearsonCorForSamplePairs")
    }
  )
}

pearsonStubCalls <- list()
pearsonStubResults <- NULL

methods::setMethod(
  "pearsonCorForSamplePairs",
  "MetaboliteAssayData",
  function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
    pearsonStubCalls[[length(pearsonStubCalls) + 1]] <<- list(
      tech_rep_remove_regex = tech_rep_remove_regex,
      correlation_group = correlation_group
    )

    pearsonStubResults
  }
)

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "plotPearson")
  },
  env = environment()
)

newMetabPlotPearsonObject <- function() {
  methods::new("MetaboliteAssayData", args = list())
}

test_that("metabolomics S4 plotPearson delegates correlation args and preserves plot labels", {
  pearsonStubCalls <<- list()
  pearsonStubResults <<- list(
    LCMS_Pos = tibble::tibble(pearson_correlation = c(0.25, 0.5, NA_real_))
  )

  plots <- suppressWarnings(
    plotPearson(
      newMetabPlotPearsonObject(),
      tech_rep_remove_regex = "pool",
      correlation_group = "Batch"
    )
  )

  expect_length(pearsonStubCalls, 1)
  expect_identical(pearsonStubCalls[[1]]$tech_rep_remove_regex, "pool")
  expect_identical(pearsonStubCalls[[1]]$correlation_group, "Batch")
  expect_named(plots, "LCMS_Pos")
  expect_s3_class(plots$LCMS_Pos, "ggplot")
  expect_identical(plots$LCMS_Pos$labels$x, "Pearson Correlation")
  expect_identical(plots$LCMS_Pos$labels$y, "Counts")
  expect_equal(plots$LCMS_Pos$scales$get_scales("x")$limits, c(0, 1))
})

test_that("metabolomics S4 plotPearson keeps unnamed fallback names and drops invalid assays", {
  pearsonStubCalls <<- list()
  pearsonStubResults <<- list(
    tibble::tibble(pearson_correlation = c(0.9, 0.8)),
    tibble::tibble(other_column = 1),
    tibble::tibble(pearson_correlation = c(NA_real_, NA_real_))
  )

  plots <- suppressWarnings(
    plotPearson(newMetabPlotPearsonObject())
  )

  expect_named(plots, "Assay_1")
  expect_length(plots, 1)
})

test_that("metabolomics S4 plotPearson keeps empty-result fallback", {
  pearsonStubCalls <<- list()
  pearsonStubResults <<- list()

  expect_warning(
    plots <- plotPearson(newMetabPlotPearsonObject()),
    "No correlation results generated"
  )

  expect_identical(plots, list())
})

test_that("metabolomics S4 plotPearson source retains histogram guards and correlation delegation", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "plotPearson")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "pearsonCorForSamplePairs\\(")
  expect_match(target_text, "geom_histogram\\(breaks = hist_breaks, na\\.rm = TRUE\\)")
  expect_match(target_text, "scale_x_continuous\\(limits = c\\(0, 1\\)")
  expect_match(
    target_text,
    "pearson_plots_list <- pearson_plots_list\\[!sapply\\(pearson_plots_list,\\s+is\\.null\\)\\]"
  )
})
