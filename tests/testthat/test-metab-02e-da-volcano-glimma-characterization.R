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
    file.path(repo_root, "R", "func_metab_da_volcano_glimma.R"),
    file.path(repo_root, "R", "func_metab_da.R")
  ),
  symbols = "generateMetabDAVolcanoPlotGlimma",
  env = environment()
)

if (!methods::isClass("MockMetabDaGlimmaObject")) {
  methods::setClass(
    "MockMetabDaGlimmaObject",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      metabolite_id_column = "character",
      sample_id = "character",
      group_id = "character"
    )
  )
}

buildMetabDaGlimmaLongResults <- function() {
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
    raw_pvalue = c(0.01, 0.20, 0.03),
    fdr_qvalue = c(0.01, 0.20, 0.03),
    significant = c("Up", "NS", "Up"),
    stringsAsFactors = FALSE
  )
}

buildMetabDaGlimmaObject <- function() {
  methods::new(
    "MockMetabDaGlimmaObject",
    metabolite_data = list(
      LCMS_Pos = data.frame(
        metabolite_id = c("M1", "M2"),
        S1 = c(10, 12),
        S2 = c(11, 13),
        check.names = FALSE
      )
    ),
    design_matrix = data.frame(
      sample_id = c("S1", "S2"),
      group = c("groupA", "groupB"),
      stringsAsFactors = FALSE
    ),
    metabolite_id_column = "metabolite_id",
    sample_id = "sample_id",
    group_id = "group"
  )
}

buildMetabDaFitObject <- function(colnames_vec) {
  coefficients <- matrix(
    c(1.2, -0.7),
    nrow = 2,
    dimnames = list(c("M1", "M2"), colnames_vec)
  )

  list(coefficients = coefficients)
}

test_that("metabolomics DA Glimma helper preserves early null exits", {
  expect_null(generateMetabDAVolcanoPlotGlimma(
    da_results_list = NULL,
    selected_contrast = "groupB-groupA"
  ))

  expect_null(generateMetabDAVolcanoPlotGlimma(
    da_results_list = list(da_metabolites_long = buildMetabDaGlimmaLongResults()),
    selected_contrast = NULL
  ))

  expect_null(generateMetabDAVolcanoPlotGlimma(
    da_results_list = list(da_metabolites_long = buildMetabDaGlimmaLongResults()),
    selected_contrast = "groupB-groupA",
    selected_assay = "Combined"
  ))
})

test_that("metabolomics DA Glimma helper preserves missing-data and missing-coefficient guards", {
  base_results <- list(
    da_metabolites_long = buildMetabDaGlimmaLongResults(),
    theObject = buildMetabDaGlimmaObject(),
    contrasts_results = list(
      LCMS_Pos = list(fit.eb = buildMetabDaFitObject("otherContrast"))
    )
  )

  expect_null(generateMetabDAVolcanoPlotGlimma(
    da_results_list = base_results,
    selected_contrast = "missing-contrast",
    selected_assay = "LCMS_Pos"
  ))

  expect_null(generateMetabDAVolcanoPlotGlimma(
    da_results_list = base_results,
    selected_contrast = "groupB-groupA",
    selected_assay = "LCMS_Pos"
  ))
})
