# fidelity-coverage-compare: shared
library(testthat)

getProtPeptideQcServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_qc_peptide_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_qc_peptide_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_server
}

getProtPeptideQcUi <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_qc_peptide_ui", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_qc_peptide_ui", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_ui
}

localProtPeptideQcNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

newSharedProtPeptideQcCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$server_calls <- list()
  captured$ui_ids <- character()
  captured$log_info <- character()
  captured
}

makeSharedProtPeptideWorkflow <- function(workflow_type) {
  list(state_manager = list(workflow_type = workflow_type))
}

withSharedProtPeptideSubmoduleMocks <- function(server_env, captured) {
  mock_frame <- parent.frame()

  record_server <- function(name) {
    force(name)
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

  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_qvalue_server",
    record_server("qvalue_filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_rollup_server",
    record_server("rollup"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_intensity_server",
    record_server("intensity_filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_protein_server",
    record_server("protein_peptide_filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_sample_server",
    record_server("sample_filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_replicate_server",
    record_server("replicate_filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_impute_server",
    record_server("imputation"), mock_frame
  )

  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_qvalue_ui",
    record_ui("qvalue_filter", "Q-Value Filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_rollup_ui",
    record_ui("rollup", "Precursor Rollup"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_intensity_ui",
    record_ui("intensity_filter", "Intensity Filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_protein_ui",
    record_ui("protein_peptide_filter", "Protein Peptides"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_sample_ui",
    record_ui("sample_filter", "Sample Quality"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_replicate_ui",
    record_ui("replicate_filter", "Replicate Filter"), mock_frame
  )
  localProtPeptideQcNamespaceBinding(
    server_env, "mod_prot_qc_peptide_impute_ui",
    record_ui("imputation", "Imputation"), mock_frame
  )
}

withSharedProtPeptideShinyMocks <- function(captured) {
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

test_that("proteomics peptide QC orchestrator preserves DIA server fan-out", {
  captured <- newSharedProtPeptideQcCapture()
  server_fn <- getProtPeptideQcServer()
  server_env <- environment(server_fn)

  withSharedProtPeptideSubmoduleMocks(server_env, captured)
  withSharedProtPeptideShinyMocks(captured)

  workflow_data <- makeSharedProtPeptideWorkflow("DIA")
  server_fn(
    id = "peptide-qc",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$module_id, "peptide-qc")
  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "name"),
    c(
      "qvalue_filter",
      "rollup",
      "intensity_filter",
      "protein_peptide_filter",
      "sample_filter",
      "replicate_filter",
      "imputation"
    )
  )
  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "id"),
    c(
      "qvalue_filter",
      "rollup",
      "intensity_filter",
      "protein_peptide_filter",
      "sample_filter",
      "replicate_filter",
      "imputation"
    )
  )
  expect_true(all(vapply(
    captured$server_calls,
    function(call) identical(call$workflowData, workflow_data),
    logical(1)
  )))
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "omicType") == "proteomics"))
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "experimentLabel") == "DIA Experiment"))
  expect_true("Peptide QC orchestrator: workflow type is DIA" %in% captured$log_info)
  expect_true("Running DIA peptide processing sub-modules." %in% captured$log_info)
})

test_that("proteomics peptide QC orchestrator skips non-DIA server fan-out", {
  captured <- newSharedProtPeptideQcCapture()
  server_fn <- getProtPeptideQcServer()
  server_env <- environment(server_fn)

  withSharedProtPeptideSubmoduleMocks(server_env, captured)
  withSharedProtPeptideShinyMocks(captured)

  workflow_data <- makeSharedProtPeptideWorkflow("LFQ")
  server_fn(
    id = "peptide-qc",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "LFQ Experiment"
  )

  expect_length(captured$server_calls, 0L)
  expect_true("Peptide QC orchestrator: workflow type is LFQ" %in% captured$log_info)
  expect_true("Skipping peptide processing for workflow type: LFQ" %in% captured$log_info)
})

test_that("proteomics peptide QC UI preserves tab wiring", {
  captured <- newSharedProtPeptideQcCapture()
  ui_fn <- getProtPeptideQcUi()
  server_env <- environment(ui_fn)

  withSharedProtPeptideSubmoduleMocks(server_env, captured)
  html <- htmltools::renderTags(ui_fn("peptide-qc"))$html

  expect_match(html, "id=\"peptide-qc-peptide_filter_tabs\"", fixed = TRUE)
  expect_false(grepl("Module not loaded", html, fixed = TRUE))
  expect_identical(
    captured$ui_ids,
    c(
      "qvalue_filter:peptide-qc-qvalue_filter",
      "rollup:peptide-qc-rollup",
      "intensity_filter:peptide-qc-intensity_filter",
      "protein_peptide_filter:peptide-qc-protein_peptide_filter",
      "sample_filter:peptide-qc-sample_filter",
      "replicate_filter:peptide-qc-replicate_filter",
      "imputation:peptide-qc-imputation"
    )
  )
})
