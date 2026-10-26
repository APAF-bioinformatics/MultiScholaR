# Extracted from test-workflow-server-shared.R:308

# prequel ----------------------------------------------------------------------
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
    new = function(audit_enabled = TRUE) {
      list(workflow_type = default_type, audit_enabled = audit_enabled)
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
expectWorkflowGatesAreMonotonic <- function(workflow_data, session, differential_key) {
  reset_status <- workflow_data$tab_status
  reset_status$design_matrix <- "disabled"
  reset_status$quality_control <- "complete"
  reset_status$normalization <- "pending"
  reset_status[[differential_key]] <- "disabled"
  reset_status$session_summary <- "disabled"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$design_matrix <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$quality_control, "complete")

  reset_status <- workflow_data$tab_status
  reset_status$quality_control <- "pending"
  reset_status$normalization <- "complete"
  reset_status[[differential_key]] <- "pending"
  reset_status$session_summary <- "disabled"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$quality_control <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$normalization, "complete")

  reset_status <- workflow_data$tab_status
  reset_status$quality_control <- "complete"
  reset_status$normalization <- "pending"
  reset_status[[differential_key]] <- "complete"
  reset_status$session_summary <- "pending"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$normalization <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status[[differential_key]], "complete")

  reset_status <- workflow_data$tab_status
  reset_status$normalization <- "complete"
  reset_status[[differential_key]] <- "pending"
  reset_status$session_summary <- "complete"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs[[differential_key]] <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$session_summary, "complete")
}

# test -------------------------------------------------------------------------
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
