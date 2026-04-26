# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinIntensityState")) {
  methods::setClass(
    "FakeSharedProteinIntensityState",
    slots = c(args = "list", protein_quant_table = "data.frame")
  )
}

getProtProteinIntensityServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_protein_intensity_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_protein_intensity_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_intensity_server
}

makeSharedProteinIntensityState <- function(groupwise = 1,
                                            max_groups = 2,
                                            cutoff = 0.5,
                                            proteins = "P0") {
  methods::new(
    "FakeSharedProteinIntensityState",
    args = list(
      removeRowsWithMissingValuesPercent = list(
        groupwise_percentage_cutoff = groupwise,
        max_groups_percentage_cutoff = max_groups,
        proteins_intensity_cutoff_percentile = cutoff
      )
    ),
    protein_quant_table = data.frame(
      Protein.Ids = proteins,
      stringsAsFactors = FALSE
    )
  )
}

copySharedProteinIntensityState <- function(theObject,
                                            args = theObject@args,
                                            protein_quant_table = theObject@protein_quant_table) {
  methods::new(
    "FakeSharedProteinIntensityState",
    args = args,
    protein_quant_table = protein_quant_table
  )
}

makeSharedProteinIntensityWorkflow <- function(current_state,
                                               captured,
                                               history = c("raw_data_s4", "protein_intensity_filtered")) {
  state_manager <- list(
    getState = function() current_state,
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

localSharedProteinIntensityBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedProteinIntensityPackageMocks <- function(server_env, captured, filtered_state) {
  localSharedProteinIntensityBinding(
    server_env,
    "updateProteinFiltering",
    function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
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
    .local_envir = parent.frame()
  )

  testthat::local_mocked_bindings(
    updateMissingValueParameters = function(theObject,
                                            min_reps_per_group,
                                            min_groups) {
      captured$missing_updates[[length(captured$missing_updates) + 1L]] <- list(
        theObject = theObject,
        min_reps_per_group = min_reps_per_group,
        min_groups = min_groups
      )

      updated_args <- theObject@args
      updated_args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- 12.345
      updated_args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- 67.89
      copySharedProteinIntensityState(theObject, args = updated_args)
    },
    updateConfigParameter = function(theObject, function_name, parameter_name, new_value) {
      captured$config_updates[[length(captured$config_updates) + 1L]] <- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )

      updated_args <- theObject@args
      updated_args[[function_name]][[parameter_name]] <- new_value
      copySharedProteinIntensityState(theObject, args = updated_args)
    },
    removeRowsWithMissingValuesPercent = function(theObject) {
      captured$filter_input <- theObject
      filtered_state
    },
    .env = server_env
  )
}

withSharedProteinIntensityUiMocks <- function(server_env,
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
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      if (evaluate_plot) {
        eval(substitute(expr), parent.frame())
      }
      "rendered-plot"
    },
    req = function(...) list(...)[[1]],
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
    renderPlot = function(expr) {
      if (evaluate_plot) {
        eval(substitute(expr), parent.frame())
      }
      "rendered-plot"
    },
    req = function(...) list(...)[[1]],
    grid.draw = function(value) {
      captured$drawn_plot <- value
      invisible(NULL)
    },
    .env = server_env
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

newSharedProteinIntensityCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$missing_updates <- list()
  captured$config_updates <- list()
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

test_that("proteomics protein intensity module preserves flexible apply behavior", {
  captured <- newSharedProteinIntensityCapture()
  server_fn <- getProtProteinIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinIntensityState()
  filtered_state <- makeSharedProteinIntensityState(
    groupwise = 12.345,
    max_groups = 67.89,
    cutoff = 1.5,
    proteins = c("P1", "P1", "P2")
  )

  withSharedProteinIntensityPackageMocks(server_env, captured, filtered_state)
  withSharedProteinIntensityUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_intensity_filter = TRUE,
      revert_protein_intensity_filter = FALSE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      proteins_intensity_cutoff_percentile = 1.5
    )
  )

  workflow_data <- makeSharedProteinIntensityWorkflow(current_state, captured)
  server_fn(
    id = "protein_intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_length(captured$missing_updates, 1L)
  expect_identical(captured$missing_updates[[1L]]$min_reps_per_group, 2)
  expect_identical(captured$missing_updates[[1L]]$min_groups, 3)
  expect_length(captured$config_updates, 1L)
  expect_identical(
    captured$config_updates[[1L]]$function_name,
    "removeRowsWithMissingValuesPercent"
  )
  expect_identical(
    captured$config_updates[[1L]]$parameter_name,
    "proteins_intensity_cutoff_percentile"
  )
  expect_equal(captured$config_updates[[1L]]$new_value, 1.5)
  expect_identical(
    captured$filter_input@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    12.345
  )
  expect_identical(
    captured$filter_input@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
    67.89
  )
  expect_equal(
    captured$filter_input@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile,
    1.5
  )
  expect_false(workflow_data$qc_params$protein_qc$intensity_filter$strict_mode)
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$min_reps_per_group, 2)
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$min_groups, 3)
  expect_equal(
    workflow_data$qc_params$protein_qc$intensity_filter$proteins_intensity_cutoff_percentile,
    1.5
  )
  expect_identical(captured$save_state$state_name, "protein_intensity_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_false(captured$save_state$config_object$strict_mode)
  expect_equal(captured$save_state$config_object$proteins_intensity_cutoff_percentile, 1.5)
  expect_identical(
    captured$save_state$description,
    "Applied FLEXIBLE protein intensity filter (adaptive thresholds)"
  )
  expect_identical(captured$plot_update$step_name, "11_protein_intensity_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(
    captured$output$protein_intensity_filter_results,
    "Mode: FLEXIBLE",
    fixed = TRUE
  )
  expect_match(
    captured$output$protein_intensity_filter_results,
    "Proteins remaining: 2",
    fixed = TRUE
  )
  expect_identical(captured$output$protein_intensity_filter_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "protein_intensity_filter_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein intensity module preserves strict apply behavior", {
  captured <- newSharedProteinIntensityCapture()
  server_fn <- getProtProteinIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinIntensityState(groupwise = 55, max_groups = 66)
  filtered_state <- makeSharedProteinIntensityState(
    groupwise = 0,
    max_groups = 0,
    cutoff = 2.5,
    proteins = c("P1", "P2", "P3")
  )

  withSharedProteinIntensityPackageMocks(server_env, captured, filtered_state)
  withSharedProteinIntensityUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_intensity_filter = TRUE,
      revert_protein_intensity_filter = FALSE,
      use_strict_mode = TRUE,
      min_reps_per_group = 4,
      min_groups = 5,
      proteins_intensity_cutoff_percentile = 2.5
    )
  )

  workflow_data <- makeSharedProteinIntensityWorkflow(current_state, captured)
  server_fn(
    id = "protein_intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(
    captured$log_info[1],
    "Protein Processing: Using STRICT MODE (no missing values allowed)"
  )
  expect_identical(
    vapply(captured$config_updates, `[[`, character(1), "parameter_name"),
    c(
      "groupwise_percentage_cutoff",
      "max_groups_percentage_cutoff",
      "proteins_intensity_cutoff_percentile"
    )
  )
  expect_equal(
    vapply(captured$config_updates, `[[`, numeric(1), "new_value"),
    c(0, 0, 2.5)
  )
  expect_identical(
    captured$filter_input@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    0
  )
  expect_identical(
    captured$filter_input@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
    0
  )
  expect_true(workflow_data$qc_params$protein_qc$intensity_filter$strict_mode)
  expect_true(is.na(workflow_data$qc_params$protein_qc$intensity_filter$min_reps_per_group))
  expect_true(is.na(workflow_data$qc_params$protein_qc$intensity_filter$min_groups))
  expect_identical(
    captured$save_state$description,
    "Applied STRICT protein intensity filter (no missing values)"
  )
  expect_match(
    captured$output$protein_intensity_filter_results,
    "Mode: STRICT",
    fixed = TRUE
  )
  expect_match(
    captured$output$protein_intensity_filter_results,
    "Proteins remaining: 3",
    fixed = TRUE
  )
  expect_identical(captured$output$protein_intensity_filter_plot, "rendered-plot")
})

test_that("proteomics protein intensity module preserves apply error behavior", {
  captured <- newSharedProteinIntensityCapture()
  server_fn <- getProtProteinIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinIntensityState()

  withSharedProteinIntensityPackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedProteinIntensityState()
  )
  testthat::local_mocked_bindings(
    removeRowsWithMissingValuesPercent = function(theObject) {
      captured$filter_input <- theObject
      stop("mock protein intensity failure")
    },
    .env = server_env
  )
  withSharedProteinIntensityUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_intensity_filter = TRUE,
      revert_protein_intensity_filter = FALSE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      proteins_intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinIntensityWorkflow(current_state, captured)
  server_fn(
    id = "protein_intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock protein intensity failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock protein intensity failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "protein_intensity_filter_working")
})

test_that("proteomics protein intensity module preserves revert behavior", {
  captured <- newSharedProteinIntensityCapture()
  server_fn <- getProtProteinIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinIntensityState()

  withSharedProteinIntensityUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_intensity_filter = FALSE,
      revert_protein_intensity_filter = TRUE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      proteins_intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinIntensityWorkflow(
    current_state,
    captured,
    history = c("raw_data_s4", "qvalue_filtered", "protein_intensity_filtered")
  )
  server_fn(
    id = "protein_intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "qvalue_filtered")
  expect_match(
    captured$output$protein_intensity_filter_results,
    "Reverted to previous state: qvalue_filtered",
    fixed = TRUE
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein intensity module preserves revert error behavior", {
  captured <- newSharedProteinIntensityCapture()
  server_fn <- getProtProteinIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinIntensityState()

  withSharedProteinIntensityUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_intensity_filter = FALSE,
      revert_protein_intensity_filter = TRUE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      proteins_intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinIntensityWorkflow(
    current_state,
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "protein_intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("reverted_state", envir = captured, inherits = FALSE))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("No previous state to revert to.", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
})
