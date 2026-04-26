# fidelity-coverage-compare: shared
library(testthat)

skipIfMissingLipidFilteringPublicBindings <- function() {
  missing <- character()
  package_ns <- asNamespace("MultiScholaR")
  if (!methods::isClass("FilteringProgressLipidomics")) {
    missing <- c(missing, "FilteringProgressLipidomics")
  }
  for (name in c("createLipidomicsAssayData", "updateLipidFiltering")) {
    if (!exists(name, envir = package_ns, inherits = FALSE)) {
      missing <- c(missing, name)
    }
  }
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Lipid filtering public contract binding(s) unavailable in this ref: %s",
        paste(missing, collapse = ", ")
      )
    )
  }
}

skipIfMissingLipidFilteringPublicBindings()

test_that("updateLipidFiltering preserves the public progress-and-plot contract", {
  tracked_globals <- c(
    "project_dirs",
    "omic_type",
    "experiment_label",
    "filtering_progress_lipidomics"
  )
  had_globals <- vapply(tracked_globals, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
  original_globals <- lapply(tracked_globals[had_globals], get, envir = .GlobalEnv, inherits = FALSE)

  on.exit({
    for (name in tracked_globals[had_globals]) {
      assign(name, original_globals[[name]], envir = .GlobalEnv)
    }
    for (name in tracked_globals[!had_globals]) {
      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  for (name in c("project_dirs", "omic_type", "experiment_label")) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }
  assign("filtering_progress_lipidomics", new("FilteringProgressLipidomics"), envir = .GlobalEnv)

  assay_data <- data.frame(
    lipid_id = c("L1", "L2", "L3"),
    S1 = c(10, 0, 5),
    S2 = c(12, NA, 15),
    S3 = c(20, 0, 1),
    check.names = FALSE
  )
  design_matrix <- data.frame(
    Run = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  lipid_object <- createLipidomicsAssayData(
    lipid_data = list(Assay1 = assay_data),
    design_matrix = design_matrix,
    lipid_id_column = "lipid_id",
    sample_id = "Run",
    group_id = "group",
    internal_standard_regex = "^IS_"
  )

  result <- NULL
  expect_warning(
    result <- suppressMessages(
      updateLipidFiltering(
        theObject = lipid_object,
        step_name = "raw",
        overwrite = TRUE,
        return_grid = TRUE
      )
    ),
    "global variables"
  )

  progress_object <- get("filtering_progress_lipidomics", envir = .GlobalEnv)

  expect_s3_class(result, "gtable")
  expect_identical(progress_object@steps, "raw")
  expect_identical(progress_object@assay_names[[1]], "Assay1")
  expect_equal(progress_object@n_lipids_total, 3)
})
