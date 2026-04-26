# fidelity-coverage-compare: shared
library(testthat)

getProtProteinQcServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_qc_protein_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_qc_protein_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_server
}

getProtProteinQcUi <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_qc_protein_ui", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_qc_protein_ui", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_ui
}

localProtProteinQcNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

newSharedProtProteinQcCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$server_calls <- list()
  captured$ui_ids <- character()
  captured$log_info <- character()
  captured
}

makeSharedProtProteinWorkflow <- function(workflow_type) {
  list(state_manager = list(workflow_type = workflow_type))
}

withSharedProtProteinSubmoduleMocks <- function(server_env, captured) {
  mock_frame <- parent.frame()

  record_server <- function(name, needs_paths = FALSE) {
    force(name)
    force(needs_paths)
    if (isTRUE(needs_paths)) {
      return(function(id, workflowData, experimentPaths, omicType, experimentLabel) {
        captured$server_calls[[length(captured$server_calls) + 1L]] <- list(
          name = name,
          id = id,
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      })
    }

    function(id, workflowData, omicType, experimentLabel) {
      captured$server_calls[[length(captured$server_calls) + 1L]] <- list(
        name = name,
        id = id,
        workflowData = workflowData,
        omicType = omicType,
        experimentLabel = experimentLabel
      )
      invisible(NULL)
    }
  }

  record_ui <- function(name, title) {
    force(name)
    force(title)
    function(id) {
      captured$ui_ids <- c(captured$ui_ids, paste(name, id, sep = ":"))
      shiny::tabPanel(title, shiny::div(`data-module-id` = id, name))
    }
  }

  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_rollup_server",
    record_server("rollup", needs_paths = TRUE), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_cleanup_server",
    record_server("cleanup"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_intensity_server",
    record_server("intensity_filter"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_dedup_server",
    record_server("duplicate_removal"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_replicate_server",
    record_server("replicate_filter", needs_paths = TRUE), mock_frame
  )

  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_rollup_ui",
    record_ui("rollup", "IQ Protein Rollup"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_cleanup_ui",
    record_ui("cleanup", "Accession Cleanup"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_intensity_ui",
    record_ui("intensity_filter", "Protein Intensity Filter"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_dedup_ui",
    record_ui("duplicate_removal", "Duplicate Removal"), mock_frame
  )
  localProtProteinQcNamespaceBinding(
    server_env, "mod_prot_qc_protein_replicate_ui",
    record_ui("replicate_filter", "Protein Replicate Filter"), mock_frame
  )
}

withSharedProtProteinShinyMocks <- function(captured) {
  mock_frame <- parent.frame()

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      captured$module_id <- id
      captured$output <- output
      invisible(NULL)
    },
    isolate = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    .package = "shiny",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    .package = "logger",
    .env = mock_frame
  )
}

test_that("mod_prot_qc_protein_server preserves DIA server fan-out", {
  captured <- newSharedProtProteinQcCapture()
  server_fn <- getProtProteinQcServer()
  server_env <- environment(server_fn)

  withSharedProtProteinSubmoduleMocks(server_env, captured)
  withSharedProtProteinShinyMocks(captured)

  workflow_data <- makeSharedProtProteinWorkflow("DIA")
  experiment_paths <- list(protein_qc_dir = tempdir())
  server_fn(
    id = "protein-qc",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$module_id, "protein-qc")
  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "name"),
    c("rollup", "cleanup", "intensity_filter", "duplicate_removal", "replicate_filter")
  )
  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "id"),
    c("rollup", "cleanup", "intensity_filter", "duplicate_removal", "replicate_filter")
  )
  expect_true(all(vapply(
    captured$server_calls,
    function(call) identical(call$workflowData, workflow_data),
    logical(1)
  )))
  expect_identical(captured$server_calls[[1L]]$experimentPaths, experiment_paths)
  expect_identical(captured$server_calls[[5L]]$experimentPaths, experiment_paths)
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "omicType") == "proteomics"))
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "experimentLabel") == "DIA Experiment"))
  expect_true("Protein QC orchestrator: workflow type is DIA" %in% captured$log_info)
  expect_true("Running DIA protein rollup sub-module." %in% captured$log_info)
  expect_true("Running common protein processing sub-modules." %in% captured$log_info)
})

test_that("mod_prot_qc_protein_server preserves non-DIA common fan-out", {
  captured <- newSharedProtProteinQcCapture()
  server_fn <- getProtProteinQcServer()
  server_env <- environment(server_fn)

  withSharedProtProteinSubmoduleMocks(server_env, captured)
  withSharedProtProteinShinyMocks(captured)

  workflow_data <- makeSharedProtProteinWorkflow("LFQ")
  server_fn(
    id = "protein-qc",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "LFQ Experiment"
  )

  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "name"),
    c("cleanup", "intensity_filter", "duplicate_removal", "replicate_filter")
  )
  expect_false("Running DIA protein rollup sub-module." %in% captured$log_info)
  expect_true("Protein QC orchestrator: workflow type is LFQ" %in% captured$log_info)
})

test_that("mod_prot_qc_protein_ui preserves DIA tab wiring", {
  captured <- newSharedProtProteinQcCapture()
  ui_fn <- getProtProteinQcUi()
  server_env <- environment(ui_fn)

  withSharedProtProteinSubmoduleMocks(server_env, captured)
  html <- htmltools::renderTags(ui_fn("protein-qc", workflow_type = "DIA"))$html

  expect_match(html, "id=\"protein-qc-protein_filter_tabs\"", fixed = TRUE)
  expect_false(grepl("Module not loaded", html, fixed = TRUE))
  expect_identical(
    captured$ui_ids,
    c(
      "rollup:protein-qc-rollup",
      "cleanup:protein-qc-cleanup",
      "intensity_filter:protein-qc-intensity_filter",
      "duplicate_removal:protein-qc-duplicate_removal",
      "replicate_filter:protein-qc-replicate_filter"
    )
  )
})

test_that("mod_prot_qc_protein_ui preserves non-DIA tab wiring", {
  captured <- newSharedProtProteinQcCapture()
  ui_fn <- getProtProteinQcUi()
  server_env <- environment(ui_fn)

  withSharedProtProteinSubmoduleMocks(server_env, captured)
  html <- htmltools::renderTags(ui_fn("protein-qc", workflow_type = "LFQ"))$html

  expect_false(grepl("protein-qc-rollup", html, fixed = TRUE))
  expect_identical(
    captured$ui_ids,
    c(
      "cleanup:protein-qc-cleanup",
      "intensity_filter:protein-qc-intensity_filter",
      "duplicate_removal:protein-qc-duplicate_removal",
      "replicate_filter:protein-qc-replicate_filter"
    )
  )
})
