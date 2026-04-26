# fidelity-coverage-compare: shared
library(testthat)

getMetabQcServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_metab_qc_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_metab_qc_server", envir = package_ns, inherits = FALSE))
  }

  mod_metab_qc_server
}

newSharedMetabQcCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$server_calls <- list()
  captured$ui_ids <- character()
  captured$log_info <- character()
  captured
}

makeSharedMetabQcWorkflow <- function(state = NULL, state_error = NULL) {
  state_manager <- list(
    getState = function() {
      if (!is.null(state_error)) {
        stop(state_error)
      }
      state
    }
  )

  list(state_manager = state_manager)
}

localSharedMetabQcNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedMetabQcSubmoduleMocks <- function(server_env, captured) {
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
  record_ui <- function(name) {
    force(name)
    function(id) {
      captured$ui_ids <- c(captured$ui_ids, paste(name, id, sep = ":"))
      paste(name, id, sep = ":")
    }
  }

  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_intensity_server",
    record_server("intensity"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_duplicates_server",
    record_server("duplicates"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_itsd_server",
    record_server("itsd"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_s4_server",
    record_server("s4_finalize"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_intensity_ui",
    record_ui("intensity"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_duplicates_ui",
    record_ui("duplicates"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_itsd_ui",
    record_ui("itsd"), mock_frame
  )
  localSharedMetabQcNamespaceBinding(
    server_env, "mod_metab_qc_s4_ui",
    record_ui("s4_finalize"), mock_frame
  )
}

withSharedMetabQcShinyMocks <- function(captured,
                                        evaluate_observe = FALSE,
                                        input_values = list()) {
  mock_frame <- parent.frame()

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      module(
        input_values,
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      captured$output <- output
      captured$module_id <- id
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      stored <- value
      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    observe = function(expr, ...) {
      if (isTRUE(evaluate_observe)) {
        eval(substitute(expr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value) || identical(value, FALSE)) {
          stop("required value missing", call. = FALSE)
        }
      }
      values[[1]]
    },
    tabsetPanel = function(..., id = NULL) {
      list(id = id, tabs = list(...))
    },
    .package = "shiny",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    log_error = function(...) invisible(NULL),
    .package = "logger",
    .env = mock_frame
  )
}

test_that("metabolomics QC orchestrator preserves the unavailable-state info panel", {
  captured <- newSharedMetabQcCapture()
  server_fn <- getMetabQcServer()
  server_env <- environment(server_fn)

  withSharedMetabQcSubmoduleMocks(server_env, captured)
  withSharedMetabQcShinyMocks(captured, evaluate_observe = TRUE)

  workflow_data <- makeSharedMetabQcWorkflow(state_error = "state unavailable")
  server_fn(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Metab Experiment",
    qc_trigger = function() NULL
  )

  expect_identical(captured$module_id, "qc")
  expect_match(
    as.character(htmltools::renderTags(captured$output$dynamic_qc_tabs)$html),
    "Please complete the 'Design Matrix' step first",
    fixed = TRUE
  )
  expect_length(captured$server_calls, 0L)
})

test_that("metabolomics QC orchestrator preserves populated tab rendering", {
  captured <- newSharedMetabQcCapture()
  server_fn <- getMetabQcServer()
  server_env <- environment(server_fn)

  withSharedMetabQcSubmoduleMocks(server_env, captured)
  withSharedMetabQcShinyMocks(captured, evaluate_observe = FALSE)

  workflow_data <- makeSharedMetabQcWorkflow(state = list(status = "ready"))
  server_fn(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Metab Experiment",
    qc_trigger = function() NULL
  )

  expect_identical(captured$output$dynamic_qc_tabs$id, "qc-metab_qc_tabs")
  expect_identical(
    captured$ui_ids,
    c(
      "intensity:qc-intensity",
      "duplicates:qc-duplicates",
      "itsd:qc-itsd",
      "s4_finalize:qc-s4_finalize"
    )
  )
})

test_that("metabolomics QC orchestrator preserves qc-trigger submodule initialization", {
  captured <- newSharedMetabQcCapture()
  server_fn <- getMetabQcServer()
  server_env <- environment(server_fn)

  withSharedMetabQcSubmoduleMocks(server_env, captured)
  withSharedMetabQcShinyMocks(captured, evaluate_observe = FALSE)

  workflow_data <- makeSharedMetabQcWorkflow(state = NULL)
  server_fn(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Metab Experiment",
    qc_trigger = function() TRUE
  )

  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "name"),
    c("intensity", "duplicates", "itsd", "s4_finalize")
  )
  expect_true(all(vapply(
    captured$server_calls,
    function(call) identical(call$workflowData, workflow_data),
    logical(1)
  )))
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "omicType") == "metabolomics"))
  expect_true(all(vapply(captured$server_calls, `[[`, character(1), "experimentLabel") == "Metab Experiment"))
  expect_true("Metabolomics QC: Initializing sub-module servers" %in% captured$log_info)
  expect_true("Metabolomics QC: Sub-modules initialized successfully" %in% captured$log_info)
})

test_that("metabolomics QC orchestrator preserves state-detected auto-initialization", {
  captured <- newSharedMetabQcCapture()
  server_fn <- getMetabQcServer()
  server_env <- environment(server_fn)

  withSharedMetabQcSubmoduleMocks(server_env, captured)
  withSharedMetabQcShinyMocks(captured, evaluate_observe = TRUE)

  workflow_data <- makeSharedMetabQcWorkflow(
    state = structure(list(status = "ready"), class = "MetaboliteAssayData")
  )
  server_fn(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Metab Experiment",
    qc_trigger = function() NULL
  )

  expect_identical(
    vapply(captured$server_calls, `[[`, character(1), "id"),
    c("intensity", "duplicates", "itsd", "s4_finalize")
  )
  expect_true("Metabolomics QC: Auto-initializing sub-modules (state detected)" %in% captured$log_info)
})
