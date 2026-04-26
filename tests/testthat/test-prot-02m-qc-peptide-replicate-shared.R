# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideReplicateState")) {
  methods::setClass(
    "FakeSharedPeptideReplicateState",
    slots = c(peptide_data = "data.frame")
  )
}

getProtPeptideReplicateServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_replicate_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_replicate_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_replicate_server
}

makeSharedPeptideReplicateState <- function(proteins = c("P0", "P0"),
                                            peptides = c("PEP_A", "PEP_B")) {
  methods::new(
    "FakeSharedPeptideReplicateState",
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Stripped.Sequence = peptides,
      replicates = seq_along(proteins),
      Intensity = seq_along(proteins) * 10,
      stringsAsFactors = FALSE
    )
  )
}

makeSharedPeptideReplicateWorkflow <- function(current_state,
                                               captured,
                                               history = c("sample_filtered", "replicate_filtered")) {
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

newSharedPeptideReplicateCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedPeptideReplicateBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localSharedPeptideReplicateBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localSharedPeptideReplicateBinding(
      env,
      name,
      bindings[[name]],
      .local_envir = .local_envir
    )
  }

  invisible(NULL)
}

withSharedPeptideReplicatePackageMocks <- function(server_env,
                                                   captured,
                                                   filtered_state,
                                                   filter_error = NULL) {
  localSharedPeptideReplicateBindings(
    server_env,
    list(
      removePeptidesWithOnlyOneReplicate = function(theObject, replicate_group_column) {
      captured$filter_input <- list(
        theObject = theObject,
        replicate_group_column = replicate_group_column
      )
      if (!is.null(filter_error)) {
        stop(filter_error)
      }
      filtered_state
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

withSharedPeptideReplicateUiMocks <- function(captured,
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

test_that("proteomics peptide replicate module preserves successful apply behavior", {
  captured <- newSharedPeptideReplicateCapture()
  server_fn <- getProtPeptideReplicateServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideReplicateState()
  filtered_state <- makeSharedPeptideReplicateState(
    proteins = c("P1", "P1", "P2"),
    peptides = c("PEP_A", "PEP_B", "PEP_C")
  )

  withSharedPeptideReplicatePackageMocks(server_env, captured, filtered_state)
  withSharedPeptideReplicateUiMocks(
    captured,
    input_values = list(
      apply_replicate_filter = TRUE,
      revert_replicate = FALSE,
      replicate_group_column = "replicates"
    )
  )

  workflow_data <- makeSharedPeptideReplicateWorkflow(current_state, captured)
  server_fn(
    id = "peptide_replicate",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$filter_input$theObject, current_state)
  expect_identical(captured$filter_input$replicate_group_column, "replicates")
  expect_identical(workflow_data$qc_params$peptide_qc$replicate_filter$replicate_group_column, "replicates")
  expect_s3_class(workflow_data$qc_params$peptide_qc$replicate_filter$timestamp, "POSIXct")
  expect_identical(captured$save_state$state_name, "replicate_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_identical(captured$save_state$config_object$replicate_group_column, "replicates")
  expect_identical(
    captured$save_state$description,
    "Applied replicate filter (removed single-replicate peptides)"
  )
  expect_identical(captured$plot_update$step_name, "7_replicate_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$replicate_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$replicate_results, "Replicate group column: replicates", fixed = TRUE)
  expect_identical(captured$output$replicate_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "replicate_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide replicate module preserves apply error behavior", {
  captured <- newSharedPeptideReplicateCapture()
  server_fn <- getProtPeptideReplicateServer()
  server_env <- environment(server_fn)

  withSharedPeptideReplicatePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideReplicateState(),
    filter_error = "mock peptide replicate failure"
  )
  withSharedPeptideReplicateUiMocks(
    captured,
    input_values = list(
      apply_replicate_filter = TRUE,
      revert_replicate = FALSE,
      replicate_group_column = "replicates"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideReplicateWorkflow(
    makeSharedPeptideReplicateState(),
    captured
  )
  server_fn(
    id = "peptide_replicate",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock peptide replicate failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock peptide replicate failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "replicate_working")
})

test_that("proteomics peptide replicate module preserves revert behavior", {
  captured <- newSharedPeptideReplicateCapture()
  server_fn <- getProtPeptideReplicateServer()
  server_env <- environment(server_fn)

  withSharedPeptideReplicatePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideReplicateState()
  )
  withSharedPeptideReplicateUiMocks(
    captured,
    input_values = list(
      apply_replicate_filter = FALSE,
      revert_replicate = TRUE,
      replicate_group_column = "replicates"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideReplicateWorkflow(
    makeSharedPeptideReplicateState(),
    captured,
    history = c("sample_filtered", "replicate_filtered")
  )
  server_fn(
    id = "peptide_replicate",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "sample_filtered")
  expect_identical(
    captured$output$replicate_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide replicate module preserves revert error behavior", {
  captured <- newSharedPeptideReplicateCapture()
  server_fn <- getProtPeptideReplicateServer()
  server_env <- environment(server_fn)

  withSharedPeptideReplicatePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideReplicateState()
  )
  withSharedPeptideReplicateUiMocks(
    captured,
    input_values = list(
      apply_replicate_filter = FALSE,
      revert_replicate = TRUE,
      replicate_group_column = "replicates"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideReplicateWorkflow(
    makeSharedPeptideReplicateState(),
    captured,
    history = "replicate_filtered"
  )
  server_fn(
    id = "peptide_replicate",
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
