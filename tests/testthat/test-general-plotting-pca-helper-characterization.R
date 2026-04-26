# fidelity-coverage-compare: shared
library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

plotPcaHelperSourcePath <- function() {
  candidates <- file.path(
    repo_root,
    "R",
    c("func_general_plotting_pca_rle_helpers.R", "func_general_plotting.R")
  )
  candidates[file.exists(candidates)][[1]]
}

findPlotPcaHelperExpression <- function(path) {
  exprs <- parse(file = path, keep.source = TRUE)

  for (expr in exprs) {
    is_assignment <- is.call(expr) &&
      length(expr) >= 3 &&
      as.character(expr[[1]]) %in% c("<-", "=")

    if (!is_assignment || !is.symbol(expr[[2]])) {
      next
    }

    if (identical(as.character(expr[[2]]), "plotPcaHelper")) {
      return(expr)
    }
  }

  NULL
}

plotPcaHelperHasNoMixOmicsFallback <- function(path = plotPcaHelperSourcePath()) {
  source_text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  grepl("requireNamespace\\(\"mixOmics\", quietly = TRUE\\)", source_text) &&
    grepl("stats::prcomp\\(", source_text)
}

plotPcaHelperInputs <- function() {
  data_matrix <- matrix(
    c(
      10, 11, 30, 31,
      20, 18, 25, 26,
      5, 9, 12, 13,
      30, 10, 21, 22,
      15, 14, 18, 19,
      7, 8, 20, 21
    ),
    nrow = 6,
    byrow = TRUE
  )
  rownames(data_matrix) <- paste0("feature_", seq_len(nrow(data_matrix)))
  colnames(data_matrix) <- c("S1", "S2", "S3", "S4")

  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    batch = c("B1", "B1", "B2", "B2"),
    label = c("Sample 1", "Sample 2", "Sample 3", "Sample 4"),
    stringsAsFactors = FALSE
  )

  list(data = data_matrix, design_matrix = design_matrix)
}

test_that("plotPcaHelper source retains the mixOmics fast path with an explicit no-mixOmics fallback", {
  helper_source_path <- plotPcaHelperSourcePath()
  target_expr <- findPlotPcaHelperExpression(helper_source_path)

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "mixOmics::pca\\(")
  if (basename(helper_source_path) == "func_general_plotting_pca_rle_helpers.R") {
    expect_match(target_text, "requireNamespace\\(\"mixOmics\", quietly = TRUE\\)")
    expect_match(target_text, "stats::prcomp\\(")
  } else {
    expect_false(grepl("stats::prcomp\\(", target_text))
  }
})

test_that("plotPcaHelper preserves a usable PCA plot when mixOmics is unavailable", {
  skip_if_not(
    plotPcaHelperHasNoMixOmicsFallback(),
    "Baseline helper predates the no-mixOmics fallback."
  )

  plot_env <- new.env(parent = globalenv())
  plot_env$checkMemoryBoth <- function(...) invisible(NULL)
  plot_env$message <- function(...) invisible(NULL)
  plot_env$requireNamespace <- function(pkg, quietly = TRUE) FALSE
  plot_env$ggplot <- ggplot2::ggplot
  plot_env$aes <- ggplot2::aes
  plot_env$geom_point <- ggplot2::geom_point
  plot_env$xlab <- ggplot2::xlab
  plot_env$ylab <- ggplot2::ylab
  plot_env$labs <- ggplot2::labs
  plot_env$theme <- ggplot2::theme
  plot_env$element_blank <- ggplot2::element_blank
  plot_env$scale_x_continuous <- ggplot2::scale_x_continuous
  plot_env$scale_y_continuous <- ggplot2::scale_y_continuous
  plot_env$scale_color_manual <- ggplot2::scale_color_manual
  plot_env$scale_shape_manual <- ggplot2::scale_shape_manual
  plot_env$coord_cartesian <- ggplot2::coord_cartesian
  plot_env$rownames_to_column <- tibble::rownames_to_column
  plot_env$left_join <- dplyr::left_join
  plot_env$sym <- rlang::sym
  plot_env$getCategoricalColourPalette <- function() {
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")
  }

  target_expr <- findPlotPcaHelperExpression(
    plotPcaHelperSourcePath()
  )
  eval(target_expr, envir = plot_env)

  data_matrix <- matrix(
    c(
      10, 11, 30,
      20, 18, 25,
      5, 9, 12,
      30, 10, 21
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(data_matrix) <- paste0("feature_", seq_len(nrow(data_matrix)))
  colnames(data_matrix) <- c("S1", "S2", "S3")

  design_matrix <- data.frame(
    Sample_ID = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    batch = c("B1", "B1", "B2"),
    stringsAsFactors = FALSE
  )

  plot_obj <- plot_env$plotPcaHelper(
    data = data_matrix,
    design_matrix = design_matrix,
    sample_id_column = "Sample_ID",
    grouping_variable = "group",
    shape_variable = "batch",
    title = "Fallback PCA",
    ncomp = 2,
    cv_percentile = 0
  )

  expect_s3_class(plot_obj, "ggplot")
  expect_identical(plot_obj$labels$title, "Fallback PCA")
  expect_true(all(c("PC1", "PC2", "Sample_ID", "group", "batch") %in% names(plot_obj$data)))
  expect_identical(sort(unique(as.character(plot_obj$data$group))), c("A", "B"))
})

test_that("plotPcaHelper package fallback returns grouped PCA plots across aesthetic branches", {
  skip_if_not(
    plotPcaHelperHasNoMixOmicsFallback(),
    "Baseline helper predates the no-mixOmics fallback."
  )

  inputs <- plotPcaHelperInputs()

  labeled_plot <- plotPcaHelper(
    data = inputs$data,
    design_matrix = inputs$design_matrix,
    sample_id_column = "Sample_ID",
    grouping_variable = "group",
    shape_variable = "batch",
    label_column = "label",
    title = "Fallback PCA",
    geom.text.size = 3,
    ncomp = 2,
    cv_percentile = 0
  )

  expect_s3_class(labeled_plot, "ggplot")
  expect_identical(labeled_plot$labels$title, "Fallback PCA")
  expect_true(all(c("PC1", "PC2", "Sample_ID", "group", "batch", "label") %in% names(labeled_plot$data)))
  expect_identical(sort(unique(as.character(labeled_plot$data$group))), c("A", "B"))

  shape_plot <- plotPcaHelper(
    data = inputs$data,
    design_matrix = inputs$design_matrix,
    sample_id_column = "Sample_ID",
    grouping_variable = "group",
    shape_variable = "batch",
    title = "Shape PCA",
    ncomp = 2,
    cv_percentile = 0
  )

  expect_s3_class(shape_plot, "ggplot")
  expect_identical(shape_plot$labels$title, "Shape PCA")

  plain_plot <- plotPcaHelper(
    data = inputs$data,
    design_matrix = inputs$design_matrix,
    sample_id_column = "Sample_ID",
    grouping_variable = "group",
    title = "Plain PCA",
    ncomp = 2,
    cv_percentile = 0
  )

  expect_s3_class(plain_plot, "ggplot")
  expect_identical(plain_plot$labels$title, "Plain PCA")
})

test_that("plotPcaHelper package fallback validates joined aesthetic columns", {
  skip_if_not(
    plotPcaHelperHasNoMixOmicsFallback(),
    "Baseline helper predates the no-mixOmics fallback."
  )

  inputs <- plotPcaHelperInputs()

  expect_error(
    plotPcaHelper(
      data = inputs$data,
      design_matrix = inputs$design_matrix,
      sample_id_column = "Sample_ID",
      grouping_variable = "missing_group",
      title = "Invalid PCA",
      ncomp = 2,
      cv_percentile = 0
    ),
    "Grouping variable 'missing_group' not found"
  )

  expect_error(
    plotPcaHelper(
      data = inputs$data,
      design_matrix = inputs$design_matrix,
      sample_id_column = "Sample_ID",
      grouping_variable = "group",
      shape_variable = "missing_shape",
      title = "Invalid PCA",
      ncomp = 2,
      cv_percentile = 0
    ),
    "Shape variable 'missing_shape' not found"
  )

  expect_error(
    plotPcaHelper(
      data = inputs$data,
      design_matrix = inputs$design_matrix,
      sample_id_column = "Sample_ID",
      grouping_variable = "group",
      label_column = "missing_label",
      title = "Invalid PCA",
      ncomp = 2,
      cv_percentile = 0
    ),
    "Label column 'missing_label' not found"
  )
})
