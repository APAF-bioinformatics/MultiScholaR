library(methods)
library(testthat)
library(dplyr)
library(tibble)

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
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      metabolite_id_column = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      design_matrix = data.frame(),
      sample_id = "Run",
      metabolite_id_column = "Name"
    )
  )
}

if (!methods::isGeneric("plotRle")) {
  methods::setGeneric(
    "plotRle",
    function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
      standardGeneric("plotRle")
    }
  )
}

plotRleCalls <- list()
plotRleHelper <- function(Y, rowinfo = NULL, yaxis_limit = c()) {
  plotRleCalls[[length(plotRleCalls) + 1]] <<- list(
    Y = Y,
    rowinfo = rowinfo,
    yaxis_limit = yaxis_limit
  )

  structure(
    list(
      samples = rownames(Y),
      metabolites = colnames(Y),
      rowinfo = rowinfo,
      yaxis_limit = yaxis_limit
    ),
    class = "mock_plot_rle"
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_plotting_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "plotRle")
  },
  env = environment()
)

newMetabPlotObject <- function(named = TRUE) {
  assay_tbl <- tibble::tibble(
    Name = c("M1", "M2", "M3"),
    annotation = c("alpha", "beta", "gamma"),
    S1 = c(10, 20, 30),
    S2 = c(15, NA, 35),
    S3 = c(40, 50, 60)
  )

  assay_list <- list(assay_tbl)
  if (named) {
    names(assay_list) <- "LCMS_Pos"
  }

  methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list,
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3"),
      Group = c("Control", "Control", "Treatment"),
      Label = c("Sample 1", "Sample 2", "Sample 3"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    metabolite_id_column = "Name"
  )
}

test_that("metabolomics S4 plotRle preserves helper routing and sample label remapping", {
  plotRleCalls <<- list()

  expect_message(
    plots <- plotRle(
      newMetabPlotObject(),
      grouping_variable = "Group",
      yaxis_limit = c(-3, 3),
      sample_label = "Label"
    ),
    "--- Entering plotRle for MetaboliteAssayData ---"
  )

  expect_named(plots, "LCMS_Pos")
  expect_length(plotRleCalls, 1)

  helper_call <- plotRleCalls[[1]]
  expect_equal(unname(rownames(helper_call$Y)), c("Sample 1", "Sample 2", "Sample 3"))
  expect_equal(colnames(helper_call$Y), c("M1", "M2", "M3"))
  expect_equal(unname(helper_call$rowinfo), c("Control", "Control", "Treatment"))
  expect_equal(helper_call$yaxis_limit, c(-3, 3))
  expect_equal(unname(plots$LCMS_Pos$samples), c("Sample 1", "Sample 2", "Sample 3"))
})

test_that("metabolomics S4 plotRle preserves unnamed assay fallback naming", {
  plotRleCalls <<- list()

  expect_warning(
    plots <- plotRle(
      newMetabPlotObject(named = FALSE),
      grouping_variable = "Group"
    ),
    "Assay list was unnamed. Using default names"
  )

  expect_named(plots, "Assay_1")
  expect_equal(unname(plotRleCalls[[1]]$rowinfo), c("Control", "Control", "Treatment"))
})

test_that("metabolomics S4 plotRle rejects missing sample labels", {
  expect_error(
    plotRle(
      newMetabPlotObject(),
      grouping_variable = "Group",
      sample_label = "MissingLabel"
    ),
    "`sample_label` 'MissingLabel' not found in design_matrix.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 plotRle source retains label remapping and helper delegation", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "plotRle")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "label_map <- setNames\\(")
  expect_match(
    target_text,
    "rowinfo_vector <- design_matrix_filtered\\[current_colnames,\\s+grouping_variable\\]"
  )
  expect_match(target_text, "plotRleHelper\\(")
})
