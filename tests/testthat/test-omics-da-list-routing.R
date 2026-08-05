library(testthat)

test_that("shared DA list routing identifies homogeneous omics domains", {
  metabolite <- makeMetabCharacterizationObject()
  lipid <- createLipidomicsAssayData(
    lipid_data = list(),
    design_matrix = data.frame(Run = "S1", group = "A"),
    sample_id = "Run",
    group_id = "group"
  )

  expect_identical(
    .resolveOmicsDaListDomain(list(metabolite), "MetaboliteAssayData", "LipidomicsAssayData"),
    "metabolomics"
  )
  expect_identical(
    .resolveOmicsDaListDomain(list(lipid), "MetaboliteAssayData", "LipidomicsAssayData"),
    "lipidomics"
  )
  expect_identical(
    .resolveOmicsDaListDomain(list(), "MetaboliteAssayData", "LipidomicsAssayData"),
    "lipidomics"
  )
  expect_error(
    .resolveOmicsDaListDomain(list(metabolite, lipid), "MetaboliteAssayData", "LipidomicsAssayData"),
    "must all inherit"
  )
})

test_that("shared DA list methods are registered exactly once", {
  methods <- c(
    "differentialAbundanceAnalysis",
    "plotNumSigDiffExpBarPlot",
    "plotVolcanoS4",
    "plotInteractiveVolcano",
    "getDaResultsLongFormat",
    "getDaResultsWideFormat"
  )

  expect_true(all(vapply(methods, methods::hasMethod, logical(1L), signature = "list")))
})
