# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("MockMetabDaHeatmapObjectShared")) {
  methods::setClass(
    "MockMetabDaHeatmapObjectShared",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      metabolite_id_column = "character",
      sample_id = "character",
      group_id = "character"
    )
  )
}

buildMetabDaHeatmapResultsShared <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "Met Two", "Met Three", "Met Four"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupB-groupA", "groupC-groupA"),
    friendly_name = c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A"
    ),
    assay = c("LCMS_Pos", "LCMS_Pos", "LCMS_Neg", "LCMS_Neg"),
    logFC = c(1.5, -0.8, 0.6, -0.4),
    fdr_qvalue = c(0.01, 0.02, 0.03, 0.20),
    stringsAsFactors = FALSE
  )
}

buildMetabDaHeatmapObjectShared <- function() {
  methods::new(
    "MockMetabDaHeatmapObjectShared",
    metabolite_data = list(
      LCMS_Pos = data.frame(
        `Alignment ID` = c("M1", "M2"),
        S1 = c(10, 12),
        S2 = c(20, 18),
        S3 = c(25, 15),
        check.names = FALSE
      ),
      LCMS_Neg = data.frame(
        `Alignment ID` = c("M3", "M4"),
        S1 = c(8, 9),
        S2 = c(14, 11),
        S3 = c(17, 13),
        check.names = FALSE
      )
    ),
    design_matrix = data.frame(
      sample_id = c("S1", "S2", "S3"),
      group = c("groupA", "groupB", "groupB"),
      stringsAsFactors = FALSE
    ),
    metabolite_id_column = "Alignment ID",
    sample_id = "sample_id",
    group_id = "group"
  )
}

test_that("metabolomics DA heatmap helper preserves early null exits", {
  expect_null(generateMetabDAHeatmap(
    da_results_list = NULL,
    selected_contrast = "groupB-groupA"
  ))

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(da_metabolites_long = buildMetabDaHeatmapResultsShared()),
    selected_contrast = NULL
  ))

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(da_metabolites_long = buildMetabDaHeatmapResultsShared()),
    selected_contrast = "missing-contrast"
  ))
})

test_that("metabolomics DA heatmap helper preserves heatmap construction and clustering output", {
  heatmap_output <- generateMetabDAHeatmap(
    da_results_list = list(
      da_metabolites_long = buildMetabDaHeatmapResultsShared(),
      theObject = buildMetabDaHeatmapObjectShared()
    ),
    selected_contrast = "groupB-groupA",
    selected_assay = "Combined",
    top_n = 3,
    tree_cut_method = "k_clusters",
    n_clusters = 2,
    show_metabolite_names = TRUE
  )

  expect_true(is.list(heatmap_output))
  expect_s4_class(heatmap_output$plot, "Heatmap")
  expect_length(heatmap_output$row_clusters, 3L)
  expect_s3_class(heatmap_output$col_clusters, "hclust")
})

test_that("metabolomics DA heatmap helper preserves assay filtering before matrix construction", {
  missing_assay_object <- methods::new(
    "MockMetabDaHeatmapObjectShared",
    metabolite_data = list(
      LCMS_Neg = data.frame(
        `Alignment ID` = "M3",
        S1 = 10,
        S2 = 11,
        check.names = FALSE
      )
    ),
    design_matrix = data.frame(
      sample_id = c("S1", "S2"),
      group = c("groupA", "groupB"),
      stringsAsFactors = FALSE
    ),
    metabolite_id_column = "Alignment ID",
    sample_id = "sample_id",
    group_id = "group"
  )

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(
      da_metabolites_long = buildMetabDaHeatmapResultsShared(),
      theObject = missing_assay_object
    ),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Pos",
    top_n = NA_integer_,
    cluster_rows = NA,
    cluster_cols = NA,
    tree_cut_method = NA_character_
  ))

  assay_heatmap <- generateMetabDAHeatmap(
    da_results_list = list(
      da_metabolites_long = buildMetabDaHeatmapResultsShared(),
      theObject = missing_assay_object
    ),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Neg"
  )
  expect_true(is.list(assay_heatmap))
  expect_s4_class(assay_heatmap$plot, "Heatmap")
})
