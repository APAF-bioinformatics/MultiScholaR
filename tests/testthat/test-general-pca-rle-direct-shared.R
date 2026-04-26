# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localFakeMixOmics <- function(env = parent.frame()) {
  old_lib_paths <- .libPaths()
  fake_lib <- file.path(tempdir(), "multischolar-fake-mixomics-lib")
  fake_pkg <- file.path(tempdir(), "multischolar-fake-mixomics-src")
  installed_pkg <- file.path(fake_lib, "mixOmics")

  if (!dir.exists(installed_pkg)) {
    unlink(fake_pkg, recursive = TRUE)
    dir.create(file.path(fake_pkg, "R"), recursive = TRUE, showWarnings = FALSE)
    dir.create(fake_lib, recursive = TRUE, showWarnings = FALSE)

    writeLines(
      c(
        "Package: mixOmics",
        "Version: 0.0.1",
        "Title: Fake mixOmics test double",
        "Description: Minimal namespace used by MultiScholaR PCA coverage tests.",
        "License: MIT",
        "Encoding: UTF-8"
      ),
      file.path(fake_pkg, "DESCRIPTION")
    )
    writeLines("export(pca)", file.path(fake_pkg, "NAMESPACE"))
    writeLines(
      c(
        "pca <- function(x, ncomp = 2) {",
        "  x <- as.matrix(x)",
        "  scores <- data.frame(",
        "    PC1 = seq_len(nrow(x)),",
        "    PC2 = seq_len(nrow(x)) * -1,",
        "    row.names = rownames(x)",
        "  )",
        "  list(",
        "    variates = list(X = scores),",
        "    prop_expl_var = list(X = c(PC1 = 0.6, PC2 = 0.3))",
        "  )",
        "}"
      ),
      file.path(fake_pkg, "R", "mixOmics.R")
    )
    utils::install.packages(fake_pkg, lib = fake_lib, repos = NULL, type = "source", quiet = TRUE)
  }

  .libPaths(c(fake_lib, old_lib_paths))
  withr::defer(.libPaths(old_lib_paths), envir = env)
  invisible(fake_lib)
}

buildPcaRleInputs <- function() {
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

  list(data = data_matrix, design = design_matrix)
}

test_that("general PCA/RLE helpers preserve fallback, list, ggpairs, and layout branches", {
  inputs <- buildPcaRleInputs()

  pca_fallback <- makeFunctionWithOverrides(
    plotPcaHelper,
    list(
      checkMemoryBoth = function(...) invisible(NULL),
      message = function(...) invisible(NULL),
      requireNamespace = function(pkg, quietly = TRUE) FALSE,
      getCategoricalColourPalette = function() c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")
    )
  )

  labeled_plot <- pca_fallback(
    data = inputs$data,
    design_matrix = inputs$design,
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
  expect_true(all(c("PC1", "PC2", "Sample_ID", "group", "batch", "label") %in% names(labeled_plot$data)))

  expect_error(
    pca_fallback(
      data = inputs$data,
      design_matrix = inputs$design,
      sample_id_column = "Sample_ID",
      grouping_variable = "missing_group",
      title = "Invalid PCA",
      ncomp = 2,
      cv_percentile = 0
    ),
    "Grouping variable 'missing_group' not found",
    fixed = TRUE
  )
  expect_error(
    pca_fallback(
      data = inputs$data,
      design_matrix = inputs$design,
      sample_id_column = "Sample_ID",
      grouping_variable = "group",
      shape_variable = "missing_shape",
      title = "Invalid PCA",
      ncomp = 2,
      cv_percentile = 0
    ),
    "Shape variable 'missing_shape' not found",
    fixed = TRUE
  )

  localFakeMixOmics()

  pca_list_plot <- plotPcaListHelper(
    data = inputs$data,
    design_matrix = inputs$design,
    sample_id_column = "Sample_ID",
    grouping_variables_list = c("group", "batch"),
    label_column = "label",
    title = "PCA List",
    geom.text.size = 2,
    ncomp = 2,
    cv_percentile = 0
  )
  expect_length(pca_list_plot, 2L)
  expect_true(all(vapply(pca_list_plot, inherits, logical(1), what = "ggplot")))

  ggpairs_helper <- makeFunctionWithOverrides(
    plotPcaGgpairs,
    list(
      ggpairs = function(data, columns, aes, legend) {
        list(data = data, columns = columns, legend = legend)
      }
    )
  )
  ggpairs_plot <- ggpairs_helper(
    data_matrix = inputs$data,
    design_matrix = inputs$design,
    grouping_variable = "group",
    sample_id_column = "Sample_ID",
    ncomp = 2
  )
  expect_true(is.list(ggpairs_plot))
  expect_identical(ggpairs_plot$columns, c("PC1 (60%)", "PC2 (30%)"))

  rle_plot <- plotRleHelper(
    Y = t(inputs$data),
    rowinfo = inputs$design$group,
    yaxis_limit = c(-1, 1)
  )
  expect_s3_class(rle_plot, "ggplot")

  min_max <- getMaxMinBoxplot(
    list(data = data.frame(min = c(-2, -1), max = c(3, 4))),
    adjust_factor = 0.10
  )
  expect_equal(min_max, c(-2.2, 4.4))

  arranged <- makeFunctionWithOverrides(
    rlePcaPlotList,
    list(
      plotPcaHelper = function(data_matrix, design_matrix, sample_id_column, grouping_variable, title, cex = 7) {
        ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point()
      },
      ggarrange = function(plotlist, nrow, ncol, common.legend, legend, widths, heights) {
        list(
          plot_count = length(plotlist),
          nrow = nrow,
          ncol = ncol,
          legend = legend
        )
      }
    )
  )(
    list_of_data_matrix = list(inputs$data, inputs$data),
    list_of_design_matrix = list(inputs$design, inputs$design),
    sample_id_column = Sample_ID,
    grouping_variable = group,
    list_of_descriptions = c("First", "Second")
  )

  expect_identical(arranged$plot_count, 4L)
  expect_identical(arranged$ncol, 2L)
  expect_identical(arranged$legend, "bottom")
})
