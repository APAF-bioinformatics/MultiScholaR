# fidelity-coverage-compare: shared
library(testthat)

makeSharedMetabQcPlottingFunction <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

buildSharedMetabQcProgress <- function() {
  methods::new(
    "FilteringProgressMetabolomics",
    steps = c("raw", "filtered"),
    assay_names = list(c("AssayA", "AssayB"), "AssayA"),
    n_metabolites_per_assay = list(
      list(AssayA = 3, AssayB = 2),
      list(AssayA = 2)
    ),
    n_metabolites_total = c(5, 2),
    detected_per_sample = list(
      list(
        AssayA = data.frame(Run = c("S1", "S2"), n_detected = c(3, 2)),
        AssayB = data.frame(Run = c("S1", "S2"), n_detected = c(2, 1))
      ),
      list(
        AssayA = data.frame(Run = c("S1", "S2"), n_detected = c(2, 1))
      )
    ),
    missingness_per_assay = list(
      list(AssayA = 10, AssayB = 25),
      list(AssayA = 40)
    ),
    sum_intensity_per_sample = list(
      list(
        AssayA = data.frame(Run = c("S1", "S2"), sum_intensity = c(8, 4)),
        AssayB = data.frame(Run = c("S1", "S2"), sum_intensity = c(16, 8))
      ),
      list(
        AssayA = data.frame(Run = c("S1", "S2"), sum_intensity = c(4, 2))
      )
    ),
    cv_distribution_per_assay = list(
      list(
        AssayA = data.frame(group = c("G1", "G2"), cv = c(5, 7)),
        AssayB = data.frame(group = "G1", cv = 10)
      ),
      list(
        AssayA = data.frame(group = "G1", cv = 4)
      )
    ),
    is_metrics_per_assay = list(
      list(
        AssayA = data.frame(is_id = "IS1", mean_intensity = 8, cv = 5),
        AssayB = data.frame(is_id = "IS2", mean_intensity = 16, cv = 8)
      ),
      list(
        AssayA = data.frame(is_id = "IS1", mean_intensity = 4, cv = 6)
      )
    )
  )
}

test_that("metabolomics QC plotting helper preserves getter fallback and empty-progress exit", {
  empty_progress <- methods::new("FilteringProgressMetabolomics")
  plotting_fun <- makeSharedMetabQcPlottingFunction(
    generateMetaboliteFilteringPlots,
    list(getFilteringProgressMetabolomics = function() empty_progress)
  )

  expect_message(
    plots <- plotting_fun(),
    "No metabolomics filtering steps have been tracked yet."
  )
  expect_identical(plots, list())
})

test_that("metabolomics QC plotting helper preserves current plot surface contracts", {
  plots <- generateMetaboliteFilteringPlots(buildSharedMetabQcProgress())

  expect_equal(
    names(plots),
    c(
      "total_metabolites",
      "metabolites_per_assay",
      "detected_per_sample",
      "missingness",
      "sum_intensity",
      "cv_distribution",
      "is_cv",
      "is_intensity"
    )
  )

  expect_s3_class(plots$total_metabolites, "ggplot")
  expect_s3_class(plots$metabolites_per_assay, "ggplot")
  expect_s3_class(plots$detected_per_sample, "ggplot")
  expect_s3_class(plots$missingness, "ggplot")
  expect_s3_class(plots$sum_intensity, "ggplot")
  expect_s3_class(plots$cv_distribution, "ggplot")
  expect_s3_class(plots$is_cv, "ggplot")
  expect_s3_class(plots$is_intensity, "ggplot")

  expect_equal(levels(plots$total_metabolites$data$step), c("raw", "filtered"))
  expect_equal(plots$total_metabolites$data$total_metabolites, c(5, 2))
  expect_equal(plots$metabolites_per_assay$data$n_metabolites, c(3, 2, 2))
  expect_true(all(c("Run", "n_detected", "step", "assay", "avg_detected") %in% names(plots$detected_per_sample$data)))
  expect_equal(plots$missingness$data$missingness, c(10, 25, 40))
  expect_equal(plots$sum_intensity$data$log2_intensity, c(3, 2, 4, 3, 2, 1))
  expect_equal(plots$cv_distribution$labels$title, "CV Distribution by Group")
  expect_equal(sort(unique(plots$cv_distribution$data$group)), c("G1", "G2"))
  expect_equal(plots$is_cv$labels$title, "Internal Standards CV")
  expect_equal(plots$is_intensity$data$log2_intensity, c(3, 4, 2))
})
