# fidelity-coverage-compare: shared
library(testthat)
suppressPackageStartupMessages(library(ggplot2))

makeSharedPeptidePcaData <- function() {
  data.frame(
    PC1 = c(-2.0, -1.4, -0.8, -0.3, 0.4, 0.9, 1.5, 2.1),
    PC2 = c(1.5, 1.1, 0.6, 0.2, -0.2, -0.7, -1.0, -1.4),
    group = rep(c("A", "B"), each = 4),
    batch = rep(c("Batch1", "Batch2"), times = 4),
    stringsAsFactors = FALSE
  )
}

expectSharedPeptideDensityPatchwork <- function(density_plot,
                                                expected_title,
                                                expected_font_size) {
  expect_true(any(inherits(density_plot, c("patchwork", "ggplot"))))

  density_panels <- c(list(density_plot), density_plot$patches$plots)
  panel_x_labels <- vapply(density_panels, function(panel) panel$labels$x, character(1))
  pc1_panel <- density_panels[[which(panel_x_labels == "PC1")[[1]]]]
  pc2_panel <- density_panels[[which(panel_x_labels == "PC2")[[1]]]]

  expect_identical(pc1_panel$labels$title, expected_title)
  expect_identical(pc1_panel$labels$y, "Density")
  expect_identical(pc2_panel$labels$y, "Density")
  expect_true(all(c("PC1", "PC2", "group") %in% colnames(pc1_panel$data)))
  expect_true(all(c("PC1", "PC2", "group") %in% colnames(pc2_panel$data)))
  expect_equal(sort(unique(as.character(pc1_panel$data$group))), c("A", "B"))
  expect_equal(sort(unique(as.character(pc2_panel$data$group))), c("A", "B"))
  expect_identical(pc1_panel$theme$text$size, expected_font_size)
  expect_identical(pc2_panel$theme$text$size, expected_font_size)
}

test_that("peptide ggplot plotDensity converts embedded PCA data into density patchworks", {
  pca_data <- makeSharedPeptidePcaData()
  pca_plot <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(PC1, PC2, colour = group)
  ) +
    ggplot2::geom_point()

  density_plot <- plotDensity(
    pca_plot,
    grouping_variable = "group",
    title = "Peptide PCA density",
    font_size = 6
  )

  expectSharedPeptideDensityPatchwork(
    density_plot,
    expected_title = "Peptide PCA density",
    expected_font_size = 6
  )
})

test_that("peptide ggplot plotDensity falls back to PCA data from the plot environment", {
  pca_plot <- local({
    data <- makeSharedPeptidePcaData()
    ggplot2::ggplot(mapping = ggplot2::aes(PC1, PC2, colour = group)) +
      ggplot2::geom_point(data = data)
  })

  density_plot <- plotDensity(
    pca_plot,
    grouping_variable = "group",
    title = "Fallback density",
    font_size = 7
  )

  expectSharedPeptideDensityPatchwork(
    density_plot,
    expected_title = "Fallback density",
    expected_font_size = 7
  )
})

test_that("peptide ggplot plotDensity reports missing grouping data", {
  pca_data <- makeSharedPeptidePcaData()

  expect_error(
    plotDensity(
      ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(),
      grouping_variable = "missing_group"
    ),
    "grouping_variable 'missing_group' not found in the data",
    fixed = TRUE
  )
})
