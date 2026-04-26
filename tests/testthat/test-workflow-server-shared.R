# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeWorkflowStateStub <- function(default_type = "DIA") {
  list(
    new = function() {
      list(workflow_type = default_type)
    }
  )
}

makeRecordingServerStub <- function(captured, key) {
  force(captured)
  force(key)

  function(id, ...) {
    captured$calls[[key]] <- list(id = id, args = list(...))
    invisible(NULL)
  }
}

makeWorkflowStepperStub <- function(steps, tab_status) {
  shiny::tags$div(
    class = "workflow-stepper-stub",
    paste(vapply(steps, `[[`, character(1), "name"), collapse = "|"),
    paste(names(tab_status), collapse = "|")
  )
}

renderSharedOutput <- function(value) {
  paste(capture.output(print(value)), collapse = "\n")
}

test_that("proteomics workflow server wires modules and advances workflow state", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  fixture_dir <- tempfile("prot-workflow-server-")
  source_dir <- file.path(fixture_dir, "source")
  results_dir <- file.path(fixture_dir, "results")
  dir.create(source_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE)
  saveRDS(
    data.frame(sequence = c("PEPTIDE_1", "PEPTIDE_2"), stringsAsFactors = FALSE),
    file.path(source_dir, "aa_seq_tbl_final.RDS")
  )

  had_global_aa <- exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE)
  old_global_aa <- if (had_global_aa) get("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE) else NULL
  withr::defer({
    if (had_global_aa) {
      assign("aa_seq_tbl_final", old_global_aa, envir = .GlobalEnv)
    } else if (exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "aa_seq_tbl_final", envir = .GlobalEnv)
    }
    unlink(fixture_dir, recursive = TRUE, force = TRUE)
  })

  localNamespaceBindings(
    package_ns,
    list(
      WorkflowState = makeWorkflowStateStub("DIA"),
      mod_prot_import_server = makeRecordingServerStub(captured, "import"),
      mod_prot_design_server = makeRecordingServerStub(captured, "design"),
      mod_prot_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_prot_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_prot_da_server = makeRecordingServerStub(captured, "da"),
      mod_prot_enrich_server = makeRecordingServerStub(captured, "enrich"),
      mod_prot_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_proteomics_server,
    args = list(
      project_dirs = list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = source_dir,
          results_dir = results_dir
        )
      ),
      omic_type = "proteomics",
      experiment_label = "workflow-demo"
    ),
    {
      session$setInputs(workflow_tabs = "normalization")
      session$flushReact()

      expect_equal(captured$calls$import$id, "setup_import")
      expect_equal(captured$calls$design$id, "design_matrix")
      expect_equal(captured$calls$qc$id, "quality_control")
      expect_equal(captured$calls$norm$id, "normalization")
      expect_equal(captured$calls$da$id, "differential_expression")
      expect_equal(captured$calls$enrich$id, "enrichment_analysis")
      expect_equal(captured$calls$summary$id, "session_summary")
      expect_equal(captured$calls$summary$args[[1]], list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = source_dir,
          results_dir = results_dir
        )
      ))
      expect_true(is.function(captured$calls$norm$args[[5]]))
      expect_equal(captured$calls$norm$args[[5]](), "normalization")
      expect_true(exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE))

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Import", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_expression, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_expression <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$enrichment_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$enrichment_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")
    }
  )
})

test_that("proteomics workflow server bypasses QC for TMT workflows", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  fixture_dir <- tempfile("prot-workflow-tmt-")
  dir.create(file.path(fixture_dir, "source"), recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))

  localNamespaceBindings(
    package_ns,
    list(
      WorkflowState = makeWorkflowStateStub("TMT"),
      mod_prot_import_server = makeRecordingServerStub(captured, "import"),
      mod_prot_design_server = makeRecordingServerStub(captured, "design"),
      mod_prot_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_prot_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_prot_da_server = makeRecordingServerStub(captured, "da"),
      mod_prot_enrich_server = makeRecordingServerStub(captured, "enrich"),
      mod_prot_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_proteomics_server,
    args = list(
      project_dirs = list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = file.path(fixture_dir, "source"),
          results_dir = file.path(fixture_dir, "results")
        )
      ),
      omic_type = "proteomics",
      experiment_label = "workflow-demo"
    ),
    {
      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()

      expect_equal(workflow_data$tab_status$quality_control, "complete")
      expect_equal(workflow_data$tab_status$normalization, "pending")
    }
  )
})

test_that("lipidomics workflow server wires modules, status progression, and launcher shell", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      WorkflowState = makeWorkflowStateStub("DIA"),
      mod_lipid_import_server = makeRecordingServerStub(captured, "import"),
      mod_lipid_design_server = makeRecordingServerStub(captured, "design"),
      mod_lipid_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_lipid_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_lipid_da_server = makeRecordingServerStub(captured, "da"),
      mod_lipid_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_lipidomics_server,
    args = list(
      project_dirs = list(
        lipidomics = list(
          base_dir = tempfile("lipid-workflow-base-"),
          source_dir = tempfile("lipid-workflow-source-"),
          results_dir = tempfile("lipid-workflow-results-")
        )
      ),
      omic_type = "lipidomics",
      experiment_label = "lipid-demo"
    ),
    {
      session$setInputs(lipidomics_tabs = "norm")
      session$flushReact()

      expect_equal(captured$calls$import$id, "import")
      expect_equal(captured$calls$design$id, "design")
      expect_equal(captured$calls$qc$id, "qc")
      expect_equal(captured$calls$norm$id, "norm")
      expect_equal(captured$calls$da$id, "de")
      expect_equal(captured$calls$summary$id, "summary")
      expect_equal(captured$calls$norm$args$selected_tab(), "norm")

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Normalize", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")
    }
  )

  launcher_capture <- new.env(parent = emptyenv())
  launcher_capture$server_calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      mod_lipidomics_ui = function(id) {
        launcher_capture$ui_id <- id
        shiny::div("lipidomics launcher ui")
      },
      mod_lipidomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        launcher_capture$server_calls[[length(launcher_capture$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        invisible(NULL)
      }
    )
  )

  local_mocked_bindings(
    shinyApp = function(ui, server, ...) {
      list(ui = ui, server = server, dots = list(...))
    },
    .package = "shiny"
  )

  launcher_dir <- tempfile("lipidomics-launcher-")
  withr::defer(unlink(launcher_dir, recursive = TRUE, force = TRUE))
  app <- run_lipidomics_app(base_dir = launcher_dir)

  expect_equal(launcher_capture$ui_id, "lipidomics_app")
  expect_true(dir.exists(file.path(launcher_dir, "data")))
  expect_true(dir.exists(file.path(launcher_dir, "results")))
  expect_true(dir.exists(file.path(launcher_dir, "scripts")))

  app$server(list(), new.env(parent = emptyenv()), list())

  expect_equal(launcher_capture$server_calls[[1]]$id, "lipidomics_app")
  expect_equal(launcher_capture$server_calls[[1]]$omic_type, "lipidomics")
  expect_equal(launcher_capture$server_calls[[1]]$experiment_label, "Lipidomics Test")
})

test_that("metabolomics workflow server wires modules, status progression, and launcher shell", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      WorkflowState = makeWorkflowStateStub("DIA"),
      mod_metab_import_server = makeRecordingServerStub(captured, "import"),
      mod_metab_design_server = makeRecordingServerStub(captured, "design"),
      mod_metab_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_metab_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_metab_da_server = makeRecordingServerStub(captured, "da"),
      mod_metab_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_metabolomics_server,
    args = list(
      project_dirs = list(
        metabolomics = list(
          base_dir = tempfile("metab-workflow-base-"),
          source_dir = tempfile("metab-workflow-source-"),
          results_dir = tempfile("metab-workflow-results-")
        )
      ),
      omic_type = "metabolomics",
      experiment_label = "metab-demo"
    ),
    {
      session$setInputs(metabolomics_tabs = "norm")
      session$flushReact()

      expect_equal(captured$calls$import$id, "import")
      expect_equal(captured$calls$design$id, "design")
      expect_equal(captured$calls$qc$id, "qc")
      expect_equal(captured$calls$norm$id, "norm")
      expect_equal(captured$calls$da$id, "de")
      expect_equal(captured$calls$summary$id, "summary")
      expect_equal(captured$calls$norm$args$selected_tab(), "norm")

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Normalize", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")
    }
  )

  launcher_capture <- new.env(parent = emptyenv())
  launcher_capture$server_calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      mod_metabolomics_ui = function(id) {
        launcher_capture$ui_id <- id
        shiny::div("metabolomics launcher ui")
      },
      mod_metabolomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        launcher_capture$server_calls[[length(launcher_capture$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        invisible(NULL)
      }
    )
  )

  local_mocked_bindings(
    shinyApp = function(ui, server, ...) {
      list(ui = ui, server = server, dots = list(...))
    },
    .package = "shiny"
  )

  launcher_dir <- tempfile("metabolomics-launcher-")
  withr::defer(unlink(launcher_dir, recursive = TRUE, force = TRUE))
  app <- run_metabolomics_app(base_dir = launcher_dir)

  expect_equal(launcher_capture$ui_id, "metabolomics_app")
  expect_true(dir.exists(file.path(launcher_dir, "data")))
  expect_true(dir.exists(file.path(launcher_dir, "results")))
  expect_true(dir.exists(file.path(launcher_dir, "scripts")))

  app$server(list(), new.env(parent = emptyenv()), list())

  expect_equal(launcher_capture$server_calls[[1]]$id, "metabolomics_app")
  expect_equal(launcher_capture$server_calls[[1]]$omic_type, "metabolomics")
  expect_equal(launcher_capture$server_calls[[1]]$experiment_label, "Metabolomics Test")
})
