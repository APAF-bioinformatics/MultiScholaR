# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

skipIfMissingMultiScholaRBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    skip(sprintf("Target-only extracted helper(s) not present: %s", paste(missing, collapse = ", ")))
  }
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

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

test_that("setup directory extracted helpers preserve parsing and omic config behavior", {
  skipIfMissingMultiScholaRBindings(
    "parseSetupDirectoriesOmicTypes",
    "getSetupDirectoriesOmicConfig"
  )

  parse_omic_types <- getMultiScholaRBinding("parseSetupDirectoriesOmicTypes")
  get_omic_config <- getMultiScholaRBinding("getSetupDirectoriesOmicConfig")

  expect_identical(
    parse_omic_types("proteomics, metabolomics  lipidomics"),
    c("proteomics", "metabolomics", "lipidomics")
  )
  expect_identical(
    parse_omic_types(c(" transcriptomics ", "integration", "")),
    c("transcriptomics", "integration")
  )
  expect_error(parse_omic_types(42), "`omic_types` must be a character vector")
  expect_error(parse_omic_types(" , "), "No valid `omic_types` provided")
  expect_error(parse_omic_types("proteomics unknownomics"), "Invalid omic_type")

  proteomics_config <- get_omic_config("proteomics")
  metabolomics_config <- get_omic_config("metabolomics")
  transcriptomics_config <- get_omic_config("transcriptomics")
  lipidomics_config <- get_omic_config("lipidomics")
  integration_config <- get_omic_config("integration")

  expect_true("peptide_qc" %in% proteomics_config$results_subdirs)
  expect_identical(metabolomics_config$global_vars$subfeature_qc_leaf, NULL)
  expect_identical(transcriptomics_config$global_vars$raw_counts_leaf, "count_data")
  expect_identical(lipidomics_config$global_vars$da_output_leaf, "da_lipids")
  expect_identical(integration_config$global_vars$mofa_inputs_leaf, file.path("mofa", "inputs"))
  expect_error(get_omic_config("unknownomics"), "Configuration not found")
})

test_that("setup directory path-list helper preserves omic-specific path aliases", {
  skipIfMissingMultiScholaRBindings(
    "buildSetupDirectoriesPathList",
    "getSetupDirectoriesOmicConfig"
  )

  build_path_list <- getMultiScholaRBinding("buildSetupDirectoriesPathList")
  get_omic_config <- getMultiScholaRBinding("getSetupDirectoriesOmicConfig")

  base_dir <- tempfile("general-filemgmt-path-list-")
  dir.create(base_dir, recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE, force = TRUE), add = TRUE)

  build_for_type <- function(omic_type) {
    omic_config <- get_omic_config(omic_type)
    paths_def <- list(
      results_base = file.path(base_dir, "results", omic_type),
      data_dir = file.path(base_dir, "data", omic_type),
      scripts_dest_dir = file.path(base_dir, "scripts", omic_type),
      publication_graphs_dir = file.path(base_dir, "results", omic_type, "publication_graphs"),
      qc_dir = file.path(base_dir, "results", omic_type, "publication_graphs", "filtering_qc"),
      time_dir = file.path(base_dir, "results", omic_type, "publication_graphs", "filtering_qc", "20260422_000000"),
      results_summary_base = file.path(base_dir, "results_summary", omic_type)
    )
    build_path_list(
      base_dir = base_dir,
      current_omic_type = omic_type,
      omic_label_dirname = paste0(omic_type, "_Study"),
      timestamp = "20260422_000000",
      current_omic_paths_def = paths_def,
      omic_config = omic_config
    )
  }

  proteomics_paths <- build_for_type("proteomics")
  metabolomics_paths <- build_for_type("metabolomics")
  transcriptomics_paths <- build_for_type("transcriptomics")
  lipidomics_paths <- build_for_type("lipidomics")
  integration_paths <- build_for_type("integration")

  expect_identical(proteomics_paths$protein_qc_dir, proteomics_paths$feature_qc_dir)
  expect_identical(proteomics_paths$peptide_qc_dir, proteomics_paths$subfeature_qc_dir)
  expect_true(dir.exists(proteomics_paths$uniprot_annotation_dir))
  expect_identical(metabolomics_paths$metabolite_qc_dir, metabolomics_paths$feature_qc_dir)
  expect_identical(metabolomics_paths$clean_metabolites_dir, metabolomics_paths$clean_features_dir)
  expect_identical(transcriptomics_paths$gene_qc_dir, transcriptomics_paths$feature_qc_dir)
  expect_identical(transcriptomics_paths$normalized_counts_dir, transcriptomics_paths$clean_features_dir)
  expect_identical(lipidomics_paths$lipid_qc_dir, lipidomics_paths$feature_qc_dir)
  expect_identical(lipidomics_paths$clean_lipids_dir, lipidomics_paths$clean_features_dir)
  expect_true(endsWith(integration_paths$mofa_inputs_dir, file.path("mofa", "inputs")))
  expect_true(dir.exists(integration_paths$mofa_inputs_dir))
})

test_that("setup directory existing-dir helper preserves overwrite, reuse, cancel, and force decisions", {
  skipIfMissingMultiScholaRBindings("handleSetupDirectoriesExistingDirs")

  handle_existing <- getMultiScholaRBinding("handleSetupDirectoriesExistingDirs")
  results_path <- "/tmp/existing-results"
  summary_path <- "/tmp/existing-summary"
  scripts_path <- "/tmp/existing-scripts"
  existing_paths <- c(results_path, summary_path, scripts_path)

  unlink_calls <- character()
  overwrite_helper <- makeFunctionWithOverrides(handle_existing, list(
    dir.exists = function(path) path %in% existing_paths,
    readline = function(prompt) "y",
    unlink = function(path, ...) {
      unlink_calls <<- c(unlink_calls, path)
      0L
    }
  ))
  overwrite_result <- overwrite_helper(
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_Study",
    results_path = results_path,
    results_summary_path = summary_path,
    scripts_path = scripts_path,
    force = FALSE,
    reuse_existing = FALSE
  )
  expect_true(overwrite_result$process_current_omic)
  expect_false(overwrite_result$reuse_current_omic_dirs)
  expect_identical(unlink_calls, existing_paths)

  prompts <- c("n", "y")
  reuse_helper <- makeFunctionWithOverrides(handle_existing, list(
    dir.exists = function(path) path %in% existing_paths,
    readline = function(prompt) {
      value <- prompts[[1L]]
      prompts <<- prompts[-1L]
      value
    },
    unlink = function(...) stop("reuse path should not unlink")
  ))
  reuse_result <- reuse_helper(
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_Study",
    results_path = results_path,
    results_summary_path = summary_path,
    scripts_path = scripts_path,
    force = FALSE,
    reuse_existing = FALSE
  )
  expect_true(reuse_result$process_current_omic)
  expect_true(reuse_result$reuse_current_omic_dirs)

  prompts <- c("n", "n")
  cancel_helper <- makeFunctionWithOverrides(handle_existing, list(
    dir.exists = function(path) path %in% existing_paths,
    readline = function(prompt) {
      value <- prompts[[1L]]
      prompts <<- prompts[-1L]
      value
    },
    unlink = function(...) stop("cancel path should not unlink")
  ))
  cancel_result <- cancel_helper(
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_Study",
    results_path = results_path,
    results_summary_path = summary_path,
    scripts_path = scripts_path,
    force = FALSE,
    reuse_existing = FALSE
  )
  expect_false(cancel_result$process_current_omic)
  expect_false(cancel_result$reuse_current_omic_dirs)

  force_unlink_calls <- character()
  force_helper <- makeFunctionWithOverrides(handle_existing, list(
    dir.exists = function(path) path %in% existing_paths,
    unlink = function(path, ...) {
      force_unlink_calls <<- c(force_unlink_calls, path)
      0L
    }
  ))
  force_result <- force_helper(
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_Study",
    results_path = results_path,
    results_summary_path = summary_path,
    scripts_path = scripts_path,
    force = TRUE,
    reuse_existing = FALSE
  )
  expect_true(force_result$process_current_omic)
  expect_false(force_result$reuse_current_omic_dirs)
  expect_identical(force_unlink_calls, existing_paths)
})

test_that("setup directory materializer preserves reuse and script-copy edge branches", {
  skipIfMissingMultiScholaRBindings(
    "materializeSetupDirectoriesStructure",
    "getSetupDirectoriesOmicConfig"
  )

  materialize_structure <- getMultiScholaRBinding("materializeSetupDirectoriesStructure")
  get_omic_config <- getMultiScholaRBinding("getSetupDirectoriesOmicConfig")
  omic_config <- get_omic_config("proteomics")

  make_paths_def <- function(root, suffix, scripts_source_exists = TRUE) {
    scripts_source_dir <- file.path(root, "scripts-source", suffix)
    if (scripts_source_exists) {
      dir.create(scripts_source_dir, recursive = TRUE)
    }
    list(
      results_base = file.path(root, "results", suffix),
      results_summary_base = file.path(root, "results_summary", suffix),
      scripts_dest_dir = file.path(root, "scripts-dest", suffix),
      scripts_source_dir = scripts_source_dir
    )
  }

  root_dir <- tempfile("general-filemgmt-materialize-")
  dir.create(root_dir, recursive = TRUE)
  on.exit(unlink(root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  reuse_missing <- make_paths_def(root_dir, "reuse-missing")
  materialize_structure(
    current_omic_paths_def = reuse_missing,
    omic_config = omic_config,
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_ReuseMissing",
    reuse_current_omic_dirs = TRUE
  )
  expect_true(dir.exists(reuse_missing$results_base))

  no_source <- make_paths_def(root_dir, "no-source", scripts_source_exists = FALSE)
  materialize_structure(
    current_omic_paths_def = no_source,
    omic_config = omic_config,
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_NoSource",
    reuse_current_omic_dirs = FALSE
  )
  expect_true(dir.exists(no_source$scripts_dest_dir))

  no_copyable_scripts <- make_paths_def(root_dir, "rmd-only")
  writeLines("---", file.path(no_copyable_scripts$scripts_source_dir, "draft.Rmd"))
  dir.create(file.path(no_copyable_scripts$scripts_source_dir, "renv"), recursive = TRUE)
  writeLines("lock", file.path(no_copyable_scripts$scripts_source_dir, "renv", "renv.lock"))
  materialize_structure(
    current_omic_paths_def = no_copyable_scripts,
    omic_config = omic_config,
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_RmdOnly",
    reuse_current_omic_dirs = FALSE
  )
  expect_true(dir.exists(no_copyable_scripts$scripts_dest_dir))
  expect_length(list.files(no_copyable_scripts$scripts_dest_dir, recursive = TRUE), 0)

  copy_failure <- make_paths_def(root_dir, "copy-failure")
  writeLines("message('copy me')", file.path(copy_failure$scripts_source_dir, "copy_me.R"))
  failed_copy_calls <- character()
  copy_failure_materializer <- makeFunctionWithOverrides(materialize_structure, list(
    file.copy = function(from, to, overwrite = TRUE) {
      failed_copy_calls <<- c(failed_copy_calls, basename(from))
      FALSE
    }
  ))
  copy_failure_materializer(
    current_omic_paths_def = copy_failure,
    omic_config = omic_config,
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_CopyFailure",
    reuse_current_omic_dirs = FALSE
  )
  expect_identical(failed_copy_calls, "copy_me.R")

  reuse_existing <- make_paths_def(root_dir, "reuse-existing")
  dir.create(reuse_existing$results_base, recursive = TRUE)
  materialize_structure(
    current_omic_paths_def = reuse_existing,
    omic_config = omic_config,
    current_omic_type = "proteomics",
    omic_label_dirname = "proteomics_ReuseExisting",
    reuse_current_omic_dirs = TRUE
  )
  expect_true(dir.exists(file.path(reuse_existing$results_base, "protein_qc")))
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
