library(methods)
library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedDefinitions <- function(paths, symbols = character(), classes = character(), env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) && length(expr) >= 3 && as.character(expr[[1]]) %in% c("<-", "=")
      if (is_assignment && is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% symbols) {
        eval(expr, envir = env)
        next
      }

      is_class_definition <- is.call(expr) &&
        is.symbol(expr[[1]]) &&
        identical(as.character(expr[[1]]), "setClass") &&
        length(expr) >= 2 &&
        is.character(expr[[2]]) &&
        expr[[2]] %in% classes

      if (is_class_definition && !methods::isClass(expr[[2]])) {
        eval(expr, envir = env)
      }
    }
  }
}

clearFilteringProgress <- function() {
  if (exists("filtering_progress_metabolomics", envir = .GlobalEnv, inherits = FALSE)) {
    rm("filtering_progress_metabolomics", envir = .GlobalEnv)
  }
}

loadSelectedDefinitions(
  paths = c(
    file.path(repo_root, "R", "func_metab_qc_progress_helpers.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R"),
    file.path(repo_root, "R", "func_metab_qc.R")
  ),
  symbols = c(
    "getFilteringProgressMetabolomics",
    "updateFilteringProgressMetabolomics"
  ),
  classes = "FilteringProgressMetabolomics",
  env = environment()
)

test_that("metabolomics QC progress helpers preserve global initialization behavior", {
  clearFilteringProgress()
  on.exit(clearFilteringProgress(), add = TRUE)

  expect_message(
    prog_met <- getFilteringProgressMetabolomics(),
    "Initialized global 'filtering_progress_metabolomics' object"
  )
  expect_s4_class(prog_met, "FilteringProgressMetabolomics")
  expect_identical(
    get("filtering_progress_metabolomics", envir = .GlobalEnv),
    prog_met
  )

  expect_no_message(
    second_prog_met <- getFilteringProgressMetabolomics()
  )
  expect_identical(second_prog_met, prog_met)
})

test_that("metabolomics QC progress helpers keep append and overwrite slot updates stable", {
  clearFilteringProgress()
  on.exit(clearFilteringProgress(), add = TRUE)

  prog_met <- getFilteringProgressMetabolomics()
  metrics_list <- list(
    AssayA = list(
      n_metabolites = 3,
      detected_per_sample = data.frame(Run = "S1", n_detected = 3),
      missingness = 12.5,
      sum_intensity_per_sample = data.frame(Run = "S1", sum_intensity = 300),
      cv_distribution = data.frame(group = "G1", cv = 0.1),
      is_metrics = list(found = TRUE)
    ),
    AssayB = list(
      n_metabolites = 5,
      detected_per_sample = data.frame(Run = "S2", n_detected = 5),
      missingness = 7.5,
      sum_intensity_per_sample = data.frame(Run = "S2", sum_intensity = 500),
      cv_distribution = data.frame(group = "G1", cv = 0.2),
      is_metrics = list(found = FALSE)
    )
  )

  expect_invisible(
    updateFilteringProgressMetabolomics(
      prog_met = prog_met,
      step_name = "raw",
      current_assay_names = c("AssayA", "AssayB"),
      metrics_list = metrics_list,
      total_metabolites = 8
    )
  )

  stored_progress <- get("filtering_progress_metabolomics", envir = .GlobalEnv)
  expect_equal(stored_progress@steps, "raw")
  expect_equal(stored_progress@assay_names[[1]], c("AssayA", "AssayB"))
  expect_equal(stored_progress@n_metabolites_per_assay[[1]], list(AssayA = 3, AssayB = 5))
  expect_equal(stored_progress@n_metabolites_total[[1]], 8)
  expect_equal(stored_progress@missingness_per_assay[[1]], list(AssayA = 12.5, AssayB = 7.5))

  expect_error(
    updateFilteringProgressMetabolomics(
      prog_met = stored_progress,
      step_name = "raw",
      current_assay_names = "AssayA",
      metrics_list = metrics_list[1],
      total_metabolites = 3
    ),
    "already exists"
  )

  overwrite_metrics <- list(
    AssayA = list(
      n_metabolites = 4,
      detected_per_sample = data.frame(Run = "S1", n_detected = 4),
      missingness = 2.5,
      sum_intensity_per_sample = data.frame(Run = "S1", sum_intensity = 400),
      cv_distribution = data.frame(group = "G1", cv = 0.05),
      is_metrics = list(found = TRUE)
    )
  )

  expect_invisible(
    updateFilteringProgressMetabolomics(
      prog_met = stored_progress,
      step_name = "raw",
      current_assay_names = "AssayA",
      metrics_list = overwrite_metrics,
      total_metabolites = 4,
      overwrite = TRUE
    )
  )

  overwritten_progress <- get("filtering_progress_metabolomics", envir = .GlobalEnv)
  expect_equal(overwritten_progress@steps, "raw")
  expect_equal(overwritten_progress@assay_names[[1]], "AssayA")
  expect_equal(overwritten_progress@n_metabolites_per_assay[[1]], list(AssayA = 4))
  expect_equal(overwritten_progress@n_metabolites_total[[1]], 4)
  expect_equal(overwritten_progress@missingness_per_assay[[1]], list(AssayA = 2.5))
})
