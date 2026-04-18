library(methods)
library(testthat)
library(ggplot2)
library(patchwork)
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

findSelectedExpressions <- function(paths, matcher) {
  matches <- list()

  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        matches[[length(matches) + 1]] <- expr
      }
    }
  }

  matches
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

if (!methods::isGeneric("plotDensity")) {
  methods::setGeneric(
    "plotDensity",
    function(theObject, grouping_variable, title = "", font_size = 8) {
      standardGeneric("plotDensity")
    }
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_plotting_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "plotDensity")
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
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    metabolite_id_column = "Name"
  )
}

newPcaPlot <- function() {
  ggplot(
    data.frame(
      PC1 = c(1, 2, 3),
      PC2 = c(3, 4, 5),
      Group = c("A", "A", "B")
    ),
    aes(PC1, PC2)
  ) +
    geom_point()
}

test_that("metabolomics S4 plotDensity preserves current PCA failure fallback without mixOmics", {
  skip_if(requireNamespace("mixOmics", quietly = TRUE), "Current characterization covers the no-mixOmics fallback path.")

  expect_warning(
    expect_warning(
      plots <- plotDensity(
        newMetabPlotObject(),
        grouping_variable = "Group",
        title = "QC Density",
        font_size = 10
      ),
      "Error during PCA calculation for density plot",
      fixed = TRUE
    ),
    "PCA result is invalid or missing PC1/PC2. Skipping Density plot.",
    fixed = TRUE
  )

  expect_length(plots, 0)
})

test_that("metabolomics S4 plotDensity list method preserves plot title composition and extracted PCA data", {
  plots <- plotDensity(
    list(LCMS_Pos = newPcaPlot()),
    grouping_variable = "Group",
    title = "QC Density",
    font_size = 11
  )

  expect_named(plots, "LCMS_Pos")
  expect_s3_class(plots[[1]], "patchwork")
  expect_equal(plots[[1]][[1]]$labels$title, "QC Density - LCMS_Pos")
  expect_equal(plots[[1]][[1]]$labels$y, "PC1")
  expect_equal(plots[[1]][[2]]$labels$y, "PC2")
  expect_equal(plots[[1]][[1]]$data$Group, c("A", "A", "B"))
  expect_equal(plots[[1]][[1]]$theme$text$size, 11)
})

test_that("metabolomics S4 plotDensity list method preserves unnamed plot fallback naming", {
  expect_warning(
    plots <- plotDensity(
      list(newPcaPlot()),
      grouping_variable = "Group",
      title = "QC Density"
    ),
    "Input ggplot list was unnamed. Using default names",
    fixed = TRUE
  )

  expect_named(plots, "Plot_1")
  expect_equal(plots[[1]][[1]]$labels$title, "QC Density - Plot_1")
})

test_that("metabolomics S4 plotDensity source retains both PCA and list-method plotting paths", {
  target_exprs <- findSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "plotDensity")
    }
  )

  expect_length(target_exprs, 2)

  target_text <- vapply(target_exprs, function(expr) {
    paste(deparse(expr), collapse = "\n")
  }, character(1))

  metab_text <- target_text[grepl('signature = "MetaboliteAssayData"', target_text, fixed = TRUE)]
  list_text <- target_text[grepl('signature = c(theObject = "list")', target_text, fixed = TRUE)]

  expect_length(metab_text, 1)
  expect_length(list_text, 1)
  expect_match(metab_text, "mixOmics::pca\\(")
  expect_match(metab_text, "patchwork::plot_layout\\(")
  expect_match(list_text, "pca_plot\\$data")
  expect_match(list_text, "plot_title_final <- if \\(!is\\.null\\(title\\) && title !=\\s+\"\"\\)")
  expect_match(list_text, "paste\\(title, \"-\", plot_name\\)")
})
