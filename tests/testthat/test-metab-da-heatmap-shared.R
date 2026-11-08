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

# KIT-101 F-03: the metab DA UI's tree-cut defaults, read from the rendered UI.
metabDaTreeCutUiDefaults <- function(ui_html, prefix) {
  select_html <- regmatches(ui_html, regexpr(
    sprintf('(?s)<select id="%s-heatmap_tree_cut_method"[^>]*>.*?</select>', prefix),
    ui_html, perl = TRUE
  ))
  numeric_default <- function(id) {
    match <- regmatches(ui_html, regexec(
      sprintf('<input id="%s-%s" type="number"[^>]*value="([^"]+)"', prefix, id),
      ui_html, perl = TRUE
    ))[[1]]
    if (length(match) < 2L) NULL else as.numeric(match[[2L]])
  }
  method <- if (length(select_html)) {
    regmatches(select_html, regexec('<option value="([^"]+)" selected>', select_html))[[1]][2]
  }
  list(
    heatmap_tree_cut_method = method,
    heatmap_n_clusters = numeric_default("heatmap_n_clusters"),
    heatmap_cut_height = numeric_default("heatmap_cut_height"),
    heatmap_min_cluster_size = numeric_default("heatmap_min_cluster_size")
  )
}

test_that("KIT-101 F-03: the metab DA UI defines lipid's four tree-cut inputs and the heatmap server receives them", {
  ui_defaults <- metabDaTreeCutUiDefaults(as.character(mod_metab_da_ui("f03")), "f03")
  expect_identical(ui_defaults, list(
    heatmap_tree_cut_method = "none",
    heatmap_n_clusters = 4,
    heatmap_cut_height = 10,
    heatmap_min_cluster_size = 5
  ))

  # KIT-106 W4: the tab is the shared DA frame (mod_omic_da.R) under metab's configuration, mounted through
  # `mod_metab_da_server`. Its heatmap render hands the inputs to `omicDaHeatmapState` (captured here); its cluster
  # summary reads the tree-cut method itself (captured as the server reads it, with the summary's value).
  captured <- new.env(parent = emptyenv())
  heatmap_arguments <- names(formals(omicDaHeatmapState))
  testthat::local_mocked_bindings(omicDaHeatmapState = function(...) {
    captured$heatmap <- stats::setNames(list(...), heatmap_arguments)
    NULL
  })
  workflow <- list2env(list(state_manager = NULL, contrasts_tbl = NULL, tab_status = list()), parent = emptyenv())
  suppressMessages(shiny::testServer(mod_metab_da_server, args = list(workflow_data = workflow, experiment_paths = list(),
    omic_type = "metabolomics", experiment_label = "f03"), {
    do.call(session$setInputs, Filter(Negate(is.null), ui_defaults))
    output$heatmap_plot
    captured$cluster_summary <- list(treeCutMethod = input$heatmap_tree_cut_method,
      value = tryCatch(output$cluster_summary, shiny.silent.error = function(e) "<req>"))
  }))
  expect_identical(captured$heatmap$treeCutMethod, "none")
  expect_identical(captured$heatmap$nClusters, 4)
  expect_identical(captured$heatmap$cutHeight, 10)
  expect_identical(captured$heatmap$minClusterSize, 5)
  expect_identical(captured$cluster_summary$treeCutMethod, "none")
})

test_that("KIT-101 F-03: with tree cutting 'none' the metab heatmap is identical to the one drawn without the inputs", {
  da_results <- list(
    da_metabolites_long = buildMetabDaHeatmapResultsShared(),
    theObject = buildMetabDaHeatmapObjectShared()
  )
  # KIT-106 W4: the shared frame's heatmap view (`omicDaHeatmapState`, then `omicDaStoreHeatmap` on the tab) under metab's
  # configuration, with the start commit's builder defaults for the other inputs.
  heatmap_view <- omicGuiConfig("metabolomics")$da$views$heatmap
  render <- function(treeCutMethod, nClusters, cutHeight, minClusterSize) {
    daData <- new.env(parent = emptyenv())
    state <- omicDaHeatmapState(
      daResultsList = da_results,
      selectedContrast = "groupB-groupA",
      selectedAssay = "Combined",
      topN = 3,
      heatmapClustering = "both",
      scaleData = "row",
      clusteringMethod = "complete",
      distanceMethod = "euclidean",
      colorScheme = "RdBu",
      showLabels = TRUE,
      daQValThresh = 0.05,
      treeCutMethod = treeCutMethod,
      nClusters = nClusters,
      cutHeight = cutHeight,
      minClusterSize = minClusterSize,
      heatmap = heatmap_view
    )
    plot <- omicDaStoreHeatmap(state, daData)
    png_path <- tempfile(fileext = ".png")
    grDevices::png(png_path, width = 600, height = 500)
    tryCatch(ComplexHeatmap::draw(plot), finally = grDevices::dev.off())
    list(
      matrix = plot@matrix,
      row_clusters = daData$current_row_clusters,
      col_clusters = daData$current_col_clusters,
      png = unname(tools::md5sum(png_path))
    )
  }
  # Today the four inputs are absent, so the server passes NULL for each.
  without_inputs <- render(NULL, NULL, NULL, NULL)
  with_defaults <- render("none", 4, 10, 5)
  expect_identical(with_defaults, without_inputs)
  expect_null(with_defaults$row_clusters)
})
