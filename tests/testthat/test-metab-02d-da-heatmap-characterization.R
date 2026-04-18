library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_metab_da_heatmap.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "generateMetabDAHeatmap",
  env = environment()
)

if (!methods::isClass("MockMetabDaHeatmapObject")) {
  methods::setClass(
    "MockMetabDaHeatmapObject",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      metabolite_id_column = "character",
      sample_id = "character",
      group_id = "character"
    )
  )
}

buildMetabDaHeatmapResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3"),
    metabolite_name = c("Met One", "Met Two", "Met Three"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA"),
    friendly_name = c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A"
    ),
    assay = c("LCMS_Pos", "LCMS_Pos", "LCMS_Neg"),
    logFC = c(1.5, -0.8, 0.6),
    fdr_qvalue = c(0.01, 0.20, 0.03),
    stringsAsFactors = FALSE
  )
}

buildMetabDaHeatmapObject <- function(metabolite_data) {
  methods::new(
    "MockMetabDaHeatmapObject",
    metabolite_data = metabolite_data,
    design_matrix = data.frame(
      sample_id = c("S1", "S2"),
      group = c("groupA", "groupB"),
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
    da_results_list = list(da_metabolites_long = buildMetabDaHeatmapResults()),
    selected_contrast = NULL
  ))

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(da_metabolites_long = buildMetabDaHeatmapResults()),
    selected_contrast = "missing-contrast"
  ))
})

test_that("metabolomics DA heatmap helper preserves assay filtering before optional package calls", {
  missing_assay_object <- buildMetabDaHeatmapObject(
    metabolite_data = list(
      LCMS_Neg = data.frame(
        `Alignment ID` = "M3",
        S1 = 10,
        S2 = 11,
        check.names = FALSE
      )
    )
  )

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(
      da_metabolites_long = buildMetabDaHeatmapResults(),
      theObject = missing_assay_object
    ),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Pos",
    top_n = NA_integer_,
    cluster_rows = NA,
    cluster_cols = NA,
    tree_cut_method = NA_character_
  ))

  expect_null(generateMetabDAHeatmap(
    da_results_list = list(
      da_metabolites_long = buildMetabDaHeatmapResults(),
      theObject = missing_assay_object
    ),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Neg"
  ))
})
