library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

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

test_that("plotPcaHelper source retains the mixOmics fast path with an explicit no-mixOmics fallback", {
  target_expr <- findPlotPcaHelperExpression(
    file.path(repo_root, "R", "func_general_plotting_pca_rle_helpers.R")
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "requireNamespace\\(\"mixOmics\", quietly = TRUE\\)")
  expect_match(target_text, "mixOmics::pca\\(")
  expect_match(target_text, "stats::prcomp\\(")
})

test_that("plotPcaHelper preserves a usable PCA plot when mixOmics is unavailable", {
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
    file.path(repo_root, "R", "func_general_plotting_pca_rle_helpers.R")
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
