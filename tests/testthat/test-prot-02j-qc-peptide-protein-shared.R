# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinPeptideState")) {
  methods::setClass(
    "FakeSharedProteinPeptideState",
    slots = c(args = "list", peptide_data = "data.frame")
  )
}

getProtPeptideProteinServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_protein_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_protein_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_protein_server
}

makeSharedProteinPeptideState <- function(proteins = c("P0"),
                                          args = list(
                                            filterMinNumPeptidesPerProtein = list(
                                              peptides_per_protein_cutoff = 1,
                                              peptidoforms_per_protein_cutoff = 1
                                            )
                                          )) {
  methods::new(
    "FakeSharedProteinPeptideState",
    args = args,
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Intensity = seq_along(proteins),
      stringsAsFactors = FALSE
    )
  )
}

copySharedProteinPeptideState <- function(theObject,
                                          args = theObject@args,
                                          peptide_data = theObject@peptide_data) {
  methods::new(
    "FakeSharedProteinPeptideState",
    args = args,
    peptide_data = peptide_data
  )
}

makeSharedProteinPeptideWorkflow <- function(current_state,
                                             captured,
                                             history = c("raw_data_s4", "protein_peptide_filtered")) {
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

newSharedProteinPeptideCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$config_updates <- list()
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedProteinPeptideBinding <- function(env, name, value, .local_envir = parent.frame()) {
  target_env <- env
  while (!identical(target_env, emptyenv()) &&
         !exists(name, envir = target_env, inherits = FALSE)) {
    target_env <- parent.env(target_env)
  }
  if (identical(target_env, emptyenv())) {
    target_env <- env
  }

  had_binding <- exists(name, envir = target_env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = target_env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, target_env)

  if (was_locked) {
    unlockBinding(name, target_env)
  }
  assign(name, value, envir = target_env)
  if (was_locked) {
    lockBinding(name, target_env)
  }

  withr::defer({
    if (exists(name, envir = target_env, inherits = FALSE) &&
        bindingIsLocked(name, target_env)) {
      unlockBinding(name, target_env)
    }
    if (had_binding) {
      assign(name, old_value, envir = target_env)
    } else if (exists(name, envir = target_env, inherits = FALSE)) {
      rm(list = name, envir = target_env)
    }
    if (was_locked && exists(name, envir = target_env, inherits = FALSE)) {
      lockBinding(name, target_env)
    }
  }, envir = .local_envir)
}

localSharedProteinPeptideBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localSharedProteinPeptideBinding(
      env,
      name,
      bindings[[name]],
      .local_envir = .local_envir
    )
  }

  invisible(NULL)
}

withSharedProteinPeptidePackageMocks <- function(server_env,
                                                 captured,
                                                 filtered_state,
                                                 filter_error = NULL) {
  localSharedProteinPeptideBindings(
    server_env,
    list(
      updateConfigParameter = function(theObject, function_name, parameter_name, new_value) {
      captured$config_updates[[length(captured$config_updates) + 1L]] <- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )

      updated_args <- theObject@args
      updated_args[[function_name]][[parameter_name]] <- new_value
      copySharedProteinPeptideState(theObject, args = updated_args)
      },
      filterMinNumPeptidesPerProtein = function(theObject) {
      captured$filter_input <- theObject
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

withSharedProteinPeptideUiMocks <- function(server_env,
                                            captured,
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

  localSharedProteinPeptideBinding(
    server_env,
    "renderText",
    function(expr) {
      eval(substitute(expr), parent.frame())
    },
    .local_envir = parent.frame()
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

test_that("proteomics peptide protein module preserves successful apply behavior", {
  captured <- newSharedProteinPeptideCapture()
  server_fn <- getProtPeptideProteinServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinPeptideState()
  filtered_state <- makeSharedProteinPeptideState(
    proteins = c("P1", "P1", "P2")
  )

  withSharedProteinPeptidePackageMocks(server_env, captured, filtered_state)
  withSharedProteinPeptideUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_peptide_filter = TRUE,
      revert_protein_peptide = FALSE,
      min_peptides_per_protein = 3,
      min_peptidoforms_per_protein = 2
    )
  )

  workflow_data <- makeSharedProteinPeptideWorkflow(current_state, captured)
  server_fn(
    id = "protein_peptide_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(
    vapply(captured$config_updates, `[[`, character(1), "parameter_name"),
    c("peptides_per_protein_cutoff", "peptidoforms_per_protein_cutoff")
  )
  expect_identical(
    vapply(captured$config_updates, `[[`, numeric(1), "new_value"),
    c(3, 2)
  )
  expect_identical(
    captured$filter_input@args$filterMinNumPeptidesPerProtein$peptides_per_protein_cutoff,
    3
  )
  expect_identical(
    captured$filter_input@args$filterMinNumPeptidesPerProtein$peptidoforms_per_protein_cutoff,
    2
  )
  expect_identical(captured$save_state$state_name, "protein_peptide_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_identical(captured$save_state$config_object$min_peptides_per_protein, 3)
  expect_identical(captured$save_state$config_object$min_peptidoforms_per_protein, 2)
  expect_identical(
    captured$save_state$description,
    "Applied minimum peptides per protein filter"
  )
  expect_identical(workflow_data$qc_params$peptide_qc$protein_peptide_filter$min_peptides_per_protein, 3)
  expect_identical(workflow_data$qc_params$peptide_qc$protein_peptide_filter$min_peptidoforms_per_protein, 2)
  expect_s3_class(workflow_data$qc_params$peptide_qc$protein_peptide_filter$timestamp, "POSIXct")
  expect_identical(captured$plot_update$step_name, "5_protein_peptide_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$protein_peptida_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$protein_peptida_results, "Min peptides per protein: 3", fixed = TRUE)
  expect_match(captured$output$protein_peptida_results, "Min peptidoforms per protein: 2", fixed = TRUE)
  expect_identical(captured$output$protein_peptide_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "protein_peptide_working")
  expect_true(any(grepl("Applying protein peptide count filter", captured$log_info, fixed = TRUE)))
  expect_true(any(grepl("applied successfully", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide protein module preserves apply error behavior", {
  captured <- newSharedProteinPeptideCapture()
  server_fn <- getProtPeptideProteinServer()
  server_env <- environment(server_fn)

  withSharedProteinPeptidePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedProteinPeptideState(),
    filter_error = "mock protein peptide failure"
  )
  withSharedProteinPeptideUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_peptide_filter = TRUE,
      revert_protein_peptide = FALSE,
      min_peptides_per_protein = 3,
      min_peptidoforms_per_protein = 2
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinPeptideWorkflow(
    makeSharedProteinPeptideState(),
    captured
  )
  server_fn(
    id = "protein_peptide_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock protein peptide failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock protein peptide failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "protein_peptide_working")
})

test_that("proteomics peptide protein module preserves revert behavior", {
  captured <- newSharedProteinPeptideCapture()
  server_fn <- getProtPeptideProteinServer()
  server_env <- environment(server_fn)

  withSharedProteinPeptideUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_peptide_filter = FALSE,
      revert_protein_peptide = TRUE,
      min_peptides_per_protein = 3,
      min_peptidoforms_per_protein = 2
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinPeptideWorkflow(
    makeSharedProteinPeptideState(),
    captured,
    history = c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered")
  )
  server_fn(
    id = "protein_peptide_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "qvalue_filtered")
  expect_identical(
    captured$output$protein_peptida_results,
    "Reverted to previous state: qvalue_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide protein module preserves revert error behavior", {
  captured <- newSharedProteinPeptideCapture()
  server_fn <- getProtPeptideProteinServer()
  server_env <- environment(server_fn)

  withSharedProteinPeptideUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_peptide_filter = FALSE,
      revert_protein_peptide = TRUE,
      min_peptides_per_protein = 3,
      min_peptidoforms_per_protein = 2
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinPeptideWorkflow(
    makeSharedProteinPeptideState(),
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "protein_peptide_filter",
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
