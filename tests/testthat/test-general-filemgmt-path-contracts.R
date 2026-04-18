library(testthat)

runSetupDirectoriesQuietly <- function(...) {
  result <- NULL
  suppressMessages(suppressWarnings(capture.output(
    result <- setupDirectories(...)
  )))
  result
}

runSetupAndShowDirectoriesQuietly <- function(...) {
  result <- NULL
  suppressMessages(suppressWarnings(capture.output(
    result <- setupAndShowDirectories(...)
  )))
  result
}

test_that("setupDirectories creates labeled proteomics paths and copies non-Rmd scripts", {
  base_dir <- tempfile("general-filemgmt-paths-")
  dir.create(base_dir, recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

  dir.create(file.path(base_dir, "scripts", "proteomics", "nested"), recursive = TRUE)
  writeLines("cat('run')", file.path(base_dir, "scripts", "proteomics", "run_me.R"))
  writeLines("cat('nested')", file.path(base_dir, "scripts", "proteomics", "nested", "nested.R"))
  writeLines("---", file.path(base_dir, "scripts", "proteomics", "draft.Rmd"))

  paths <- runSetupDirectoriesQuietly(
    base_dir = base_dir,
    omic_types = "proteomics",
    label = "Alpha",
    force = TRUE
  )

  expect_named(paths, "proteomics_Alpha")

  proteomics_paths <- paths$proteomics_Alpha

  expect_equal(proteomics_paths$omic_type, "proteomics")
  expect_equal(proteomics_paths$omic_label, "proteomics_Alpha")
  expect_equal(proteomics_paths$data_dir, file.path(base_dir, "data", "proteomics"))
  expect_true(dir.exists(proteomics_paths$results_dir))
  expect_true(dir.exists(proteomics_paths$results_summary_dir))
  expect_true(dir.exists(proteomics_paths$source_dir))
  expect_true(dir.exists(proteomics_paths$protein_qc_dir))
  expect_true(dir.exists(proteomics_paths$peptide_qc_dir))
  expect_true(dir.exists(proteomics_paths$publication_graphs_dir))
  expect_true(dir.exists(proteomics_paths$time_dir))
  expect_true(dir.exists(proteomics_paths$uniprot_annotation_dir))
  expect_true(file.exists(file.path(proteomics_paths$source_dir, "run_me.R")))
  expect_true(file.exists(file.path(proteomics_paths$source_dir, "nested", "nested.R")))
  expect_false(file.exists(file.path(proteomics_paths$source_dir, "draft.Rmd")))
})

test_that("setupDirectories reuses an existing directory structure without recopying scripts", {
  base_dir <- tempfile("general-filemgmt-reuse-")
  dir.create(base_dir, recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

  dir.create(file.path(base_dir, "results", "proteomics_Reuse"), recursive = TRUE)
  dir.create(file.path(base_dir, "results_summary", "proteomics_Reuse"), recursive = TRUE)
  dir.create(file.path(base_dir, "scripts", "proteomics_Reuse"), recursive = TRUE)
  writeLines("cat('keep')", file.path(base_dir, "scripts", "proteomics_Reuse", "existing.R"))

  dir.create(file.path(base_dir, "scripts", "proteomics"), recursive = TRUE)
  writeLines("cat('new')", file.path(base_dir, "scripts", "proteomics", "new.R"))

  paths <- runSetupDirectoriesQuietly(
    base_dir = base_dir,
    omic_types = "proteomics",
    label = "Reuse",
    reuse_existing = TRUE
  )

  expect_named(paths, "proteomics_Reuse")
  expect_equal(paths$proteomics_Reuse$source_dir, file.path(base_dir, "scripts", "proteomics_Reuse"))
  expect_true(file.exists(file.path(paths$proteomics_Reuse$source_dir, "existing.R")))
  expect_false(file.exists(file.path(paths$proteomics_Reuse$source_dir, "new.R")))
  expect_true(dir.exists(paths$proteomics_Reuse$time_dir))
})

test_that("setupAndShowDirectories creates proteomics directories and assigns globals", {
  base_dir <- tempfile("general-filemgmt-show-")
  dir.create(base_dir, recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

  dir.create(file.path(base_dir, "scripts", "proteomics", "nested"), recursive = TRUE)
  writeLines("cat('run')", file.path(base_dir, "scripts", "proteomics", "run_me.R"))
  writeLines("cat('nested')", file.path(base_dir, "scripts", "proteomics", "nested", "nested.R"))
  writeLines("---", file.path(base_dir, "scripts", "proteomics", "draft.Rmd"))

  binding_names <- c(
    "base_dir", "integration_dir", "results_dir", "data_dir", "source_dir",
    "da_output_dir", "publication_graphs_dir", "timestamp", "qc_dir",
    "time_dir", "results_summary_dir", "pathway_dir", "protein_qc_dir",
    "peptide_qc_dir", "clean_proteins_dir", "qc_figures_dir",
    "publication_figures_dir", "publication_tables_dir", "study_report_dir"
  )
  had_existing <- vapply(binding_names, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
  previous_values <- if (any(had_existing)) {
    mget(binding_names[had_existing], envir = .GlobalEnv, inherits = FALSE)
  } else {
    list()
  }
  on.exit({
    new_bindings <- binding_names[!had_existing]
    if (length(new_bindings)) {
      rm(list = new_bindings, envir = .GlobalEnv)
    }
    if (length(previous_values)) {
      list2env(previous_values, envir = .GlobalEnv)
    }
  }, add = TRUE)

  paths <- runSetupAndShowDirectoriesQuietly(
    base_dir = base_dir,
    label = "Alpha",
    force = TRUE
  )

  expect_equal(paths$results_dir, file.path(base_dir, "results", "proteomics_Alpha"))
  expect_equal(paths$results_summary_dir, file.path(base_dir, "results_summary", "proteomics_Alpha"))
  expect_equal(paths$data_dir, file.path(base_dir, "data"))
  expect_equal(paths$source_dir, file.path(base_dir, "scripts", "proteomics_Alpha"))
  expect_true(dir.exists(paths$integration_dir))
  expect_true(dir.exists(paths$results_dir))
  expect_true(dir.exists(paths$results_summary_dir))
  expect_true(dir.exists(paths$source_dir))
  expect_true(dir.exists(paths$protein_qc_dir))
  expect_true(dir.exists(paths$peptide_qc_dir))
  expect_true(dir.exists(paths$publication_graphs_dir))
  expect_true(dir.exists(paths$qc_dir))
  expect_true(dir.exists(paths$time_dir))
  expect_true(startsWith(paths$time_dir, paths$qc_dir))
  expect_true(file.exists(file.path(paths$source_dir, "run_me.R")))
  expect_true(file.exists(file.path(paths$source_dir, "nested", "nested.R")))
  expect_false(file.exists(file.path(paths$source_dir, "draft.Rmd")))
  expect_identical(get("results_dir", envir = .GlobalEnv), paths$results_dir)
  expect_identical(get("time_dir", envir = .GlobalEnv), paths$time_dir)
})

test_that("getProjectPaths resolves direct and fallback keys from a supplied environment", {
  env <- new.env(parent = emptyenv())
  env$project_dirs <- list(
    proteomics_Project42 = list(results_dir = "/tmp/proteomics-project42"),
    metabolomics = list(results_dir = "/tmp/metabolomics")
  )

  expect_identical(
    getProjectPaths(
      omic_type = "proteomics",
      experiment_label = "Project42",
      env = env
    ),
    env$project_dirs$proteomics_Project42
  )

  expect_identical(
    getProjectPaths(
      omic_type = "proteomics",
      experiment_label = "MissingLabel",
      env = env
    ),
    env$project_dirs$proteomics_Project42
  )

  expect_identical(
    getProjectPaths(
      omic_type = "metabolomics",
      env = env
    ),
    env$project_dirs$metabolomics
  )

  expect_error(
    getProjectPaths(
      omic_type = "lipidomics",
      env = env
    ),
    "Could not find project paths for omic_type='lipidomics'"
  )
})

test_that("directory helper wrappers create nested directories idempotently", {
  root_dir <- tempfile("general-filemgmt-helper-")
  on.exit(unlink(root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  nested_one <- file.path(root_dir, "one", "two")
  nested_two <- file.path(root_dir, "alpha", "beta")

  expect_false(dir.exists(nested_one))
  expect_false(dir.exists(nested_two))

  expect_invisible(createDirectoryIfNotExists(nested_one))
  expect_true(dir.exists(nested_one))

  expect_invisible(createDirectoryIfNotExists(nested_one))
  expect_true(dir.exists(nested_one))

  expect_invisible(createDirIfNotExists(nested_two))
  expect_true(dir.exists(nested_two))
})
