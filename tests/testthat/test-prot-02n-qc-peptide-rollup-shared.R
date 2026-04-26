# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideRollupState")) {
  methods::setClass(
    "FakeSharedPeptideRollupState",
    slots = c(peptide_data = "data.frame")
  )
}

getProtPeptideRollupServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_rollup_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_rollup_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_rollup_server
}

makeSharedPeptideRollupState <- function(proteins = c("P0", "P0"),
                                         peptides = c("PEP_A", "PEP_B")) {
  methods::new(
    "FakeSharedPeptideRollupState",
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Stripped.Sequence = peptides,
      Intensity = seq_along(proteins) * 10,
      stringsAsFactors = FALSE
    )
  )
}

makeSharedPeptideRollupWorkflow <- function(current_state,
                                            captured,
                                            history = c("qvalue_filtered", "precursor_rollup")) {
  state_manager <- list(
    getState = function(...) current_state,
    saveState = function(...) {
      captured$save_state <- list(...)
      invisible(NULL)
    },
    getHistory = function() history,
    revertToState = function(state_name) {
      captured$reverted_state <- state_name
      current_state
    }
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data
}

newSharedPeptideRollupCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedPeptideRollupBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localSharedPeptideRollupBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localSharedPeptideRollupBinding(
      env,
      name,
      bindings[[name]],
      .local_envir = .local_envir
    )
  }

  invisible(NULL)
}

withSharedPeptideRollupPackageMocks <- function(server_env,
                                                captured,
                                                rolled_up_state,
                                                rollup_error = NULL) {
  localSharedPeptideRollupBindings(
    server_env,
    list(
      rollUpPrecursorToPeptide = function(theObject) {
      captured$rollup_input <- theObject
      if (!is.null(rollup_error)) {
        stop(rollup_error)
      }
      rolled_up_state
      },
      updateProteinFiltering = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$plot_update <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot-token"
      }
    ),
    .local_envir = parent.frame()
  )
}

withSharedPeptideRollupUiMocks <- function(captured,
                                           input_values,
                                           evaluate_plot = TRUE) {
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) list(...)[[1]],
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      if (evaluate_plot) {
        eval(substitute(expr), parent.frame())
      }
      "rendered-plot"
    },
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notifications <- c(captured$removed_notifications, id)
      invisible(NULL)
    },
    .package = "shiny",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    log_info = function(message) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    log_error = function(message) {
      captured$log_error <- c(captured$log_error, message)
      invisible(NULL)
    },
    .package = "logger",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    grid.draw = function(value) {
      captured$drawn_plot <- value
      invisible(NULL)
    },
    .package = "grid",
    .env = mock_frame
  )
}

test_that("mod_prot_qc_peptide_rollup_server preserves successful apply behavior", {
  captured <- newSharedPeptideRollupCapture()
  server_fn <- getProtPeptideRollupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideRollupState()
  rolled_up_state <- makeSharedPeptideRollupState(
    proteins = c("P1", "P1", "P2"),
    peptides = c("PEP_A", "PEP_B", "PEP_C")
  )

  withSharedPeptideRollupPackageMocks(server_env, captured, rolled_up_state)
  withSharedPeptideRollupUiMocks(
    captured,
    input_values = list(
      apply_rollup = TRUE,
      revert_rollup = FALSE
    )
  )

  workflow_data <- makeSharedPeptideRollupWorkflow(current_state, captured)
  server_fn(
    id = "peptide_rollup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$rollup_input, current_state)
  expect_true(workflow_data$qc_params$peptide_qc$precursor_rollup$applied)
  expect_s3_class(workflow_data$qc_params$peptide_qc$precursor_rollup$timestamp, "POSIXct")
  expect_identical(captured$save_state$state_name, "precursor_rollup")
  expect_identical(captured$save_state$s4_data_object, rolled_up_state)
  expect_identical(captured$save_state$config_object, list())
  expect_identical(
    captured$save_state$description,
    "Applied precursor to peptide rollup"
  )
  expect_identical(captured$plot_update$step_name, "3_precursor_rollup")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$rollup_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$rollup_results, "State saved as: 'precursor_rollup'", fixed = TRUE)
  expect_identical(captured$output$rollup_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "rollup_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("mod_prot_qc_peptide_rollup_server preserves apply error behavior", {
  captured <- newSharedPeptideRollupCapture()
  server_fn <- getProtPeptideRollupServer()
  server_env <- environment(server_fn)

  withSharedPeptideRollupPackageMocks(
    server_env,
    captured,
    rolled_up_state = makeSharedPeptideRollupState(),
    rollup_error = "mock peptide rollup failure"
  )
  withSharedPeptideRollupUiMocks(
    captured,
    input_values = list(
      apply_rollup = TRUE,
      revert_rollup = FALSE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideRollupWorkflow(
    makeSharedPeptideRollupState(),
    captured
  )
  server_fn(
    id = "peptide_rollup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock peptide rollup failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock peptide rollup failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "rollup_working")
})

test_that("mod_prot_qc_peptide_rollup_server preserves revert behavior", {
  captured <- newSharedPeptideRollupCapture()
  server_fn <- getProtPeptideRollupServer()
  server_env <- environment(server_fn)

  withSharedPeptideRollupPackageMocks(
    server_env,
    captured,
    rolled_up_state = makeSharedPeptideRollupState()
  )
  withSharedPeptideRollupUiMocks(
    captured,
    input_values = list(
      apply_rollup = FALSE,
      revert_rollup = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideRollupWorkflow(
    makeSharedPeptideRollupState(),
    captured,
    history = c("qvalue_filtered", "precursor_rollup")
  )
  server_fn(
    id = "peptide_rollup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "qvalue_filtered")
  expect_identical(
    captured$output$rollup_results,
    "Reverted to previous state: qvalue_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("mod_prot_qc_peptide_rollup_server preserves revert error behavior", {
  captured <- newSharedPeptideRollupCapture()
  server_fn <- getProtPeptideRollupServer()
  server_env <- environment(server_fn)

  withSharedPeptideRollupPackageMocks(
    server_env,
    captured,
    rolled_up_state = makeSharedPeptideRollupState()
  )
  withSharedPeptideRollupUiMocks(
    captured,
    input_values = list(
      apply_rollup = FALSE,
      revert_rollup = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideRollupWorkflow(
    makeSharedPeptideRollupState(),
    captured,
    history = "precursor_rollup"
  )
  server_fn(
    id = "peptide_rollup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("reverted_state", envir = captured, inherits = FALSE))
  expect_true(any(grepl("No previous state to revert to.", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("No previous state to revert to.", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
})
