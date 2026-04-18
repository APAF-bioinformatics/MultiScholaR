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
    file.path(repo_root, "R", "func_metab_da_volcano_static.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "generateMetabDAVolcanoStatic",
  env = environment()
)

buildMetabDaLongResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M1", "M3", "M4"),
    metabolite_name = c("Met One", "", "Met One", "Met Three", "Met Four"),
    comparison = c(
      "groupB-groupA",
      "groupB-groupA",
      "groupB-groupA",
      "groupB-groupA",
      "groupC-groupA"
    ),
    friendly_name = c(
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupB-groupA = B vs A",
      "groupC-groupA = C vs A"
    ),
    assay = c("LCMS_Pos", "LCMS_Pos", "LCMS_Neg", "LCMS_Neg", "LCMS_Pos"),
    logFC = c(1.8, -1.1, 0.9, -1.5, 0.5),
    fdr_qvalue = c(0.001, 0.02, 0.04, 0.03, 0.2),
    significant = c("Up", "Down", "Up", "Down", "NS"),
    stringsAsFactors = FALSE
  )
}

test_that("metabolomics DA static volcano helper preserves combined-view faceting and labels", {
  plot <- generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabDaLongResults()),
    selected_contrast = "groupB-groupA",
    selected_assay = "Combined",
    da_q_val_thresh = 0.05,
    show_labels = TRUE,
    n_labels = 2
  )

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$facet, "FacetWrap")
  expect_equal(sort(unique(plot$data$assay)), c("LCMS_Neg", "LCMS_Pos"))
  expect_equal(plot$labels$title, "Volcano Plot: groupB-groupA")
  expect_equal(
    sort(unique(plot$data$display_name)),
    c("M2", "Met One", "Met Three")
  )
  expect_equal(
    sum(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomTextRepel"), logical(1))),
    1
  )
})

test_that("metabolomics DA static volcano helper preserves assay filtering and null exits", {
  plot <- generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabDaLongResults()),
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Pos",
    show_labels = FALSE
  )

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$facet, "FacetNull")
  expect_true(all(plot$data$assay == "LCMS_Pos"))
  expect_equal(
    sum(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomTextRepel"), logical(1))),
    0
  )

  expect_null(generateMetabDAVolcanoStatic(
    da_results_list = list(da_metabolites_long = buildMetabDaLongResults()),
    selected_contrast = "missing-contrast"
  ))
  expect_null(generateMetabDAVolcanoStatic(
    da_results_list = NULL,
    selected_contrast = "groupB-groupA"
  ))
})
