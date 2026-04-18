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

if (!methods::isGeneric("plotPca")) {
  methods::setGeneric(
    "plotPca",
    function(
      theObject,
      grouping_variable,
      shape_variable = NULL,
      label_column = NULL,
      title = NULL,
      font_size = 8
    ) {
      standardGeneric("plotPca")
    }
  )
}

plotPcaCalls <- list()
plotPcaHelper <- function(
  data,
  design_matrix,
  sample_id_column,
  grouping_variable,
  shape_variable = NULL,
  label_column = NULL,
  title = NULL,
  geom.text.size = 8
) {
  plotPcaCalls[[length(plotPcaCalls) + 1]] <<- list(
    data = data,
    design_matrix = design_matrix,
    sample_id_column = sample_id_column,
    grouping_variable = grouping_variable,
    shape_variable = shape_variable,
    label_column = label_column,
    title = title,
    geom.text.size = geom.text.size
  )

  structure(
    list(
      title = title,
      samples = colnames(data),
      groups = design_matrix[[grouping_variable]]
    ),
    class = "mock_plot_pca"
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_plotting_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "plotPca")
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
      Batch = c("B1", "B1", "B2"),
      Label = c("Sample 1", "Sample 2", "Sample 3"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    metabolite_id_column = "Name"
  )
}

test_that("metabolomics S4 plotPca preserves helper routing and assay title composition", {
  plotPcaCalls <<- list()

  plots <- plotPca(
    newMetabPlotObject(),
    grouping_variable = "Group",
    shape_variable = "Batch",
    label_column = "Label",
    title = "QC PCA",
    font_size = 11
  )

  expect_named(plots, "LCMS_Pos")
  expect_length(plotPcaCalls, 1)

  helper_call <- plotPcaCalls[[1]]
  expect_equal(helper_call$sample_id_column, "Run")
  expect_equal(helper_call$grouping_variable, "Group")
  expect_equal(helper_call$shape_variable, "Batch")
  expect_equal(helper_call$label_column, "Label")
  expect_equal(helper_call$title, "QC PCA - LCMS_Pos")
  expect_equal(helper_call$geom.text.size, 11)
  expect_equal(rownames(helper_call$data), c("M1", "M2", "M3"))
  expect_equal(colnames(helper_call$data), c("S1", "S2", "S3"))
  expect_equal(helper_call$design_matrix$Run, c("S1", "S2", "S3"))
  expect_equal(plots$LCMS_Pos$title, "QC PCA - LCMS_Pos")
})

test_that("metabolomics S4 plotPca preserves unnamed assay fallback naming", {
  plotPcaCalls <<- list()

  expect_warning(
    plots <- plotPca(
      newMetabPlotObject(named = FALSE),
      grouping_variable = "Group"
    ),
    "Assay list was unnamed. Using default names"
  )

  expect_named(plots, "Assay_1")
  expect_equal(plotPcaCalls[[1]]$title, "")
})

test_that("metabolomics S4 plotPca rejects missing design variables", {
  expect_error(
    plotPca(newMetabPlotObject(), grouping_variable = "MissingGroup"),
    "`grouping_variable` 'MissingGroup' not found in design_matrix.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 plotPca source retains finite-data guards and helper delegation", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "plotPca")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(
    target_text,
    "valid_rows <- rowSums\\(is\\.finite\\(frozen_metabolite_matrix_pca\\)\\) >\\s+1"
  )
  expect_match(target_text, "plotPcaHelper\\(")
  expect_match(
    target_text,
    "assay_title <- if \\(!is\\.null\\(title\\) && title != \"\"\\)\\s+paste\\(title, \"-\", assay_name\\)\\s+else \"\""
  )
})
