# fidelity-coverage-compare: shared
library(testthat)

buildMetabVolcanoResultsDirect <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "", "Met Three", "Met Four"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA", "groupC-groupA"),
    friendly_name = c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A",
      "groupC-groupA = C vs A"
    ),
    assay = c("LCMS_Pos", "LCMS_Neg", "LCMS_Pos", "LCMS_Neg"),
    logFC = c(1.8, -1.2, 0.4, -2.2),
    fdr_qvalue = c(0.01, 0.02, 0.20, 0.03),
    significant = c("Up", "Down", "NS", "Down"),
    stringsAsFactors = FALSE
  )
}

test_that("metabolomics static volcano helper preserves null exits and filtered ggplot output", {
  expect_null(generateMetabDAVolcanoStatic(NULL, selected_contrast = "groupB-groupA"))
  expect_null(generateMetabDAVolcanoStatic(
    list(da_metabolites_long = buildMetabVolcanoResultsDirect()),
    selected_contrast = NULL
  ))
  expect_null(generateMetabDAVolcanoStatic(
    list(da_metabolites_long = buildMetabVolcanoResultsDirect()),
    selected_contrast = "missing-contrast"
  ))

  combined_plot <- generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabVolcanoResultsDirect()),
    selected_contrast = "groupB-groupA = B vs A",
    selected_assay = "Combined",
    show_labels = TRUE,
    n_labels = 2
  )

  expect_s3_class(combined_plot, "ggplot")
  expect_true("GeomTextRepel" %in% vapply(combined_plot$layers, function(x) class(x$geom)[1], character(1)))
  expect_true(inherits(combined_plot$facet, "FacetNull"))

  faceted_plot <- generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabVolcanoResultsDirect()),
    selected_contrast = "groupC-groupA",
    selected_assay = NULL,
    show_labels = FALSE,
    n_labels = 1
  )

  expect_s3_class(faceted_plot, "ggplot")
  expect_false("GeomTextRepel" %in% vapply(faceted_plot$layers, function(x) class(x$geom)[1], character(1)))
  expect_false(inherits(faceted_plot$facet, "FacetNull"))

  assay_plot <- generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabVolcanoResultsDirect()),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Pos",
    show_labels = TRUE,
    n_labels = 1
  )

  expect_s3_class(assay_plot, "ggplot")
  expect_identical(unique(ggplot2::ggplot_build(assay_plot)$data[[1]]$colour), "#E64B35")
})
