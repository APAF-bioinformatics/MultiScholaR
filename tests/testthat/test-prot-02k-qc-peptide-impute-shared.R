# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideImputationState")) {
  methods::setClass(
    "FakeSharedPeptideImputationState",
    slots = c(args = "list", peptide_data = "data.frame")
  )
}

getProtPeptideImputationServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_impute_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_impute_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_impute_server
}

makeSharedPeptideImputationState <- function(proteins = c("P0"),
                                             args = list(
                                               peptideMissingValueImputation = list(
                                                 proportion_missing_values = 0.5
                                               )
                                             )) {
  methods::new(
    "FakeSharedPeptideImputationState",
    args = args,
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Intensity = seq_along(proteins),
      stringsAsFactors = FALSE
    )
  )
}

copySharedPeptideImputationState <- function(theObject,
                                             args = theObject@args,
                                             peptide_data = theObject@peptide_data) {
  methods::new(
    "FakeSharedPeptideImputationState",
    args = args,
    peptide_data = peptide_data
  )
}

makeSharedPeptideImputationWorkflow <- function(current_state,
                                                captured,
                                                history = c("sample_filtered", "imputed")) {
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

newSharedPeptideImputationCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$log_error <- character()
  captured
}

localSharedPeptideImputationBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localSharedPeptideImputationBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localSharedPeptideImputationBinding(
      env,
      name,
      bindings[[name]],
      .local_envir = .local_envir
    )
  }

  invisible(NULL)
}

withSharedPeptideImputationPackageMocks <- function(server_env,
                                                    captured,
                                                    imputed_state,
                                                    imputation_error = NULL) {
  localSharedPeptideImputationBindings(
    server_env,
    list(
      updateConfigParameter = function(theObject, function_name, parameter_name, new_value) {
      captured$config_update <- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )

      updated_args <- theObject@args
      updated_args[[function_name]][[parameter_name]] <- new_value
      copySharedPeptideImputationState(theObject, args = updated_args)
      },
      peptideMissingValueImputation = function(theObject) {
      captured$imputation_input <- theObject
      if (!is.null(imputation_error)) {
        stop(imputation_error)
      }
      imputed_state
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
      },
      .capture_checkpoint = function(object, checkpoint_id, checkpoint_name) {
      captured$checkpoint <- list(
        object = object,
        checkpoint_id = checkpoint_id,
        checkpoint_name = checkpoint_name
      )
      invisible(NULL)
      },
      renderText = function(expr) {
      eval(substitute(expr), parent.frame())
      },
      renderPlot = function(expr) {
      eval(substitute(expr), parent.frame())
      "rendered-plot"
      },
      req = function(...) list(...)[[1]],
      grid.draw = function(value) {
      captured$drawn_plot <- value
      invisible(NULL)
      }
    ),
    .local_envir = parent.frame()
  )
}

withSharedPeptideImputationUiMocks <- function(captured,
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
    log_warn = function(message) {
      captured$log_warn <- c(captured$log_warn, message)
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

test_that("proteomics peptide imputation module preserves successful apply behavior", {
  captured <- newSharedPeptideImputationCapture()
  server_fn <- getProtPeptideImputationServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideImputationState()
  imputed_state <- makeSharedPeptideImputationState(
    proteins = c("P1", "P1", "P2")
  )

  withSharedPeptideImputationPackageMocks(server_env, captured, imputed_state)
  withSharedPeptideImputationUiMocks(
    captured,
    input_values = list(
      apply_imputation = TRUE,
      revert_imputation = FALSE,
      proportion_missing_values = 0.4
    )
  )

  workflow_data <- makeSharedPeptideImputationWorkflow(current_state, captured)
  server_fn(
    id = "imputation",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$config_update$function_name, "peptideMissingValueImputation")
  expect_identical(captured$config_update$parameter_name, "proportion_missing_values")
  expect_identical(captured$config_update$new_value, 0.4)
  expect_identical(
    captured$imputation_input@args$peptideMissingValueImputation$proportion_missing_values,
    0.4
  )
  expect_identical(workflow_data$qc_params$peptide_qc$imputation$proportion_missing_values, 0.4)
  expect_s3_class(workflow_data$qc_params$peptide_qc$imputation$timestamp, "POSIXct")
  expect_identical(captured$save_state$state_name, "imputed")
  expect_identical(captured$save_state$s4_data_object, imputed_state)
  expect_identical(captured$save_state$config_object$proportion_missing_values, 0.4)
  expect_identical(
    captured$save_state$description,
    "Applied missing value imputation using technical replicates"
  )
  expect_identical(captured$checkpoint$object, imputed_state)
  expect_identical(captured$checkpoint$checkpoint_id, "cp02")
  expect_identical(captured$checkpoint$checkpoint_name, "qc_filtered_peptide")
  expect_identical(captured$plot_update$step_name, "8_imputed")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$imputation_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$imputation_results, "Max proportion missing: 0.4", fixed = TRUE)
  expect_match(captured$output$imputation_results, "State saved as: 'imputed'", fixed = TRUE)
  expect_identical(captured$output$imputation_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "imputation_working")
  expect_true(any(grepl("Applying missing value imputation", captured$log_info, fixed = TRUE)))
  expect_true(any(grepl("applied successfully", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide imputation module preserves apply error behavior", {
  captured <- newSharedPeptideImputationCapture()
  server_fn <- getProtPeptideImputationServer()
  server_env <- environment(server_fn)

  withSharedPeptideImputationPackageMocks(
    server_env,
    captured,
    imputed_state = makeSharedPeptideImputationState(),
    imputation_error = "mock imputation failure"
  )
  withSharedPeptideImputationUiMocks(
    captured,
    input_values = list(
      apply_imputation = TRUE,
      revert_imputation = FALSE,
      proportion_missing_values = 0.4
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideImputationWorkflow(
    makeSharedPeptideImputationState(),
    captured
  )
  server_fn(
    id = "imputation",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock imputation failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock imputation failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "imputation_working")
})

test_that("proteomics peptide imputation module preserves revert behavior", {
  captured <- newSharedPeptideImputationCapture()
  server_fn <- getProtPeptideImputationServer()
  server_env <- environment(server_fn)

  withSharedPeptideImputationPackageMocks(
    server_env,
    captured,
    imputed_state = makeSharedPeptideImputationState()
  )
  withSharedPeptideImputationUiMocks(
    captured,
    input_values = list(
      apply_imputation = FALSE,
      revert_imputation = TRUE,
      proportion_missing_values = 0.4
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideImputationWorkflow(
    makeSharedPeptideImputationState(),
    captured,
    history = c("raw_data_s4", "sample_filtered", "imputed")
  )
  server_fn(
    id = "imputation",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "sample_filtered")
  expect_identical(
    captured$output$imputation_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide imputation module preserves revert error behavior", {
  captured <- newSharedPeptideImputationCapture()
  server_fn <- getProtPeptideImputationServer()
  server_env <- environment(server_fn)

  withSharedPeptideImputationPackageMocks(
    server_env,
    captured,
    imputed_state = makeSharedPeptideImputationState()
  )
  withSharedPeptideImputationUiMocks(
    captured,
    input_values = list(
      apply_imputation = FALSE,
      revert_imputation = TRUE,
      proportion_missing_values = 0.4
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideImputationWorkflow(
    makeSharedPeptideImputationState(),
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "imputation",
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
