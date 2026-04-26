# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideIntensityState")) {
  methods::setClass(
    "FakeSharedPeptideIntensityState",
    slots = c(args = "list", peptide_data = "data.frame")
  )
}

getProtPeptideIntensityServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_intensity_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_intensity_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_intensity_server
}

makeSharedPeptideIntensityState <- function(groupwise = 1,
                                            max_groups = 2,
                                            cutoff = 0.5,
                                            proteins = "P0") {
  methods::new(
    "FakeSharedPeptideIntensityState",
    args = list(
      peptideIntensityFiltering = list(
        groupwise_percentage_cutoff = groupwise,
        max_groups_percentage_cutoff = max_groups,
        peptides_intensity_cutoff_percentile = cutoff
      )
    ),
    peptide_data = data.frame(
      Protein.Ids = proteins,
      stringsAsFactors = FALSE
    )
  )
}

copySharedPeptideIntensityState <- function(theObject,
                                            args = theObject@args,
                                            peptide_data = theObject@peptide_data) {
  methods::new(
    "FakeSharedPeptideIntensityState",
    args = args,
    peptide_data = peptide_data
  )
}

makeSharedPeptideIntensityWorkflow <- function(current_state,
                                               captured,
                                               history = c("raw_data_s4", "intensity_filtered")) {
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

localSharedPeptideIntensityBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedPeptideIntensityPackageMocks <- function(server_env, captured, filtered_state) {
  localSharedPeptideIntensityBinding(
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
                                            min_groups,
                                            function_name,
                                            grouping_variable) {
      captured$missing_updates[[length(captured$missing_updates) + 1L]] <- list(
        theObject = theObject,
        min_reps_per_group = min_reps_per_group,
        min_groups = min_groups,
        function_name = function_name,
        grouping_variable = grouping_variable
      )

      updated_args <- theObject@args
      updated_args$peptideIntensityFiltering$groupwise_percentage_cutoff <- 12.345
      updated_args$peptideIntensityFiltering$max_groups_percentage_cutoff <- 67.89
      copySharedPeptideIntensityState(theObject, args = updated_args)
    },
    updateConfigParameter = function(theObject, function_name, parameter_name, new_value) {
      captured$config_updates[[length(captured$config_updates) + 1L]] <- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )

      updated_args <- theObject@args
      updated_args[[function_name]][[parameter_name]] <- new_value
      copySharedPeptideIntensityState(theObject, args = updated_args)
    },
    peptideIntensityFiltering = function(theObject) {
      captured$filter_input <- theObject
      filtered_state
    },
    .env = server_env
  )
}

withSharedPeptideIntensityUiMocks <- function(captured, input_values, evaluate_plot = TRUE) {
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
    reactive = function(expr) {
      expr_quo <- substitute(expr)
      expr_env <- parent.frame()
      function() eval(expr_quo, expr_env)
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

newSharedPeptideIntensityCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$missing_updates <- list()
  captured$config_updates <- list()
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

test_that("proteomics peptide intensity module preserves flexible apply behavior", {
  captured <- newSharedPeptideIntensityCapture()
  server_fn <- getProtPeptideIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideIntensityState()
  filtered_state <- makeSharedPeptideIntensityState(
    groupwise = 12.345,
    max_groups = 67.89,
    cutoff = 1.5,
    proteins = c("P1", "P1", "P2")
  )

  withSharedPeptideIntensityPackageMocks(server_env, captured, filtered_state)
  withSharedPeptideIntensityUiMocks(
    captured,
    input_values = list(
      apply_intensity_filter = TRUE,
      revert_intensity = FALSE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      intensity_cutoff_percentile = 1.5
    )
  )

  workflow_data <- makeSharedPeptideIntensityWorkflow(current_state, captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_length(captured$missing_updates, 3L)
  expect_identical(captured$missing_updates[[1L]]$function_name, "peptideIntensityFiltering")
  expect_identical(captured$missing_updates[[1L]]$grouping_variable, "group")
  expect_identical(captured$missing_updates[[1L]]$min_reps_per_group, 2)
  expect_identical(captured$missing_updates[[1L]]$min_groups, 3)
  expect_length(captured$config_updates, 1L)
  expect_identical(
    captured$config_updates[[1L]]$parameter_name,
    "peptides_intensity_cutoff_percentile"
  )
  expect_equal(captured$config_updates[[1L]]$new_value, 1.5)
  expect_identical(
    captured$filter_input@args$peptideIntensityFiltering$groupwise_percentage_cutoff,
    12.345
  )
  expect_identical(
    captured$filter_input@args$peptideIntensityFiltering$max_groups_percentage_cutoff,
    67.89
  )
  expect_equal(
    captured$filter_input@args$peptideIntensityFiltering$peptides_intensity_cutoff_percentile,
    1.5
  )
  expect_false(workflow_data$qc_params$peptide_qc$intensity_filter$strict_mode)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$min_reps_per_group, 2)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$min_groups, 3)
  expect_equal(workflow_data$qc_params$peptide_qc$intensity_filter$intensity_cutoff_percentile, 1.5)
  expect_identical(captured$save_state$state_name, "intensity_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_false(captured$save_state$config_object$strict_mode)
  expect_equal(captured$save_state$config_object$intensity_cutoff_percentile, 1.5)
  expect_identical(captured$save_state$description, "Applied FLEXIBLE peptide intensity filter")
  expect_identical(captured$plot_update$step_name, "4_intensity_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_identical(captured$output$calculated_groupwise_percent, "Groupwise % cutoff: 12.345%")
  expect_identical(captured$output$calculated_max_groups_percent, "Max groups % cutoff: 67.890%")
  expect_match(captured$output$intensity_results, "Mode: FLEXIBLE", fixed = TRUE)
  expect_match(captured$output$intensity_results, "Proteins remaining: 2", fixed = TRUE)
  expect_identical(captured$output$intensity_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "intensity_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide intensity module preserves strict apply behavior", {
  captured <- newSharedPeptideIntensityCapture()
  server_fn <- getProtPeptideIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideIntensityState(groupwise = 55, max_groups = 66)
  filtered_state <- makeSharedPeptideIntensityState(
    groupwise = 0,
    max_groups = 0,
    cutoff = 2.5,
    proteins = c("P1", "P2", "P3")
  )

  withSharedPeptideIntensityPackageMocks(server_env, captured, filtered_state)
  withSharedPeptideIntensityUiMocks(
    captured,
    input_values = list(
      apply_intensity_filter = TRUE,
      revert_intensity = FALSE,
      use_strict_mode = TRUE,
      min_reps_per_group = 4,
      min_groups = 5,
      intensity_cutoff_percentile = 2.5
    )
  )

  workflow_data <- makeSharedPeptideIntensityWorkflow(current_state, captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$log_info, "Peptide Processing: Using STRICT MODE")
  expect_identical(
    vapply(captured$config_updates, `[[`, character(1), "parameter_name"),
    c(
      "groupwise_percentage_cutoff",
      "max_groups_percentage_cutoff",
      "peptides_intensity_cutoff_percentile"
    )
  )
  expect_equal(
    vapply(captured$config_updates, `[[`, numeric(1), "new_value"),
    c(0, 0, 2.5)
  )
  expect_identical(
    captured$filter_input@args$peptideIntensityFiltering$groupwise_percentage_cutoff,
    0
  )
  expect_identical(
    captured$filter_input@args$peptideIntensityFiltering$max_groups_percentage_cutoff,
    0
  )
  expect_true(workflow_data$qc_params$peptide_qc$intensity_filter$strict_mode)
  expect_true(is.na(workflow_data$qc_params$peptide_qc$intensity_filter$min_reps_per_group))
  expect_true(is.na(workflow_data$qc_params$peptide_qc$intensity_filter$min_groups))
  expect_identical(captured$save_state$description, "Applied STRICT peptide intensity filter")
  expect_match(captured$output$intensity_results, "Mode: STRICT", fixed = TRUE)
  expect_match(captured$output$intensity_results, "Proteins remaining: 3", fixed = TRUE)
  expect_identical(captured$output$intensity_plot, "rendered-plot")
})

test_that("proteomics peptide intensity module preserves apply error behavior", {
  captured <- newSharedPeptideIntensityCapture()
  server_fn <- getProtPeptideIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideIntensityState()

  withSharedPeptideIntensityPackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideIntensityState()
  )
  testthat::local_mocked_bindings(
    peptideIntensityFiltering = function(theObject) {
      captured$filter_input <- theObject
      stop("mock intensity failure")
    },
    .env = server_env
  )
  withSharedPeptideIntensityUiMocks(
    captured,
    input_values = list(
      apply_intensity_filter = TRUE,
      revert_intensity = FALSE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideIntensityWorkflow(current_state, captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock intensity failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock intensity failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "intensity_working")
})

test_that("proteomics peptide intensity module preserves revert behavior", {
  captured <- newSharedPeptideIntensityCapture()
  server_fn <- getProtPeptideIntensityServer()
  current_state <- makeSharedPeptideIntensityState()

  withSharedPeptideIntensityUiMocks(
    captured,
    input_values = list(
      apply_intensity_filter = FALSE,
      revert_intensity = TRUE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideIntensityWorkflow(
    current_state,
    captured,
    history = c("raw_data_s4", "qvalue_filtered", "intensity_filtered")
  )
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "qvalue_filtered")
  expect_match(
    captured$output$intensity_results,
    "Reverted to previous state: qvalue_filtered",
    fixed = TRUE
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide intensity module preserves revert error behavior", {
  captured <- newSharedPeptideIntensityCapture()
  server_fn <- getProtPeptideIntensityServer()
  current_state <- makeSharedPeptideIntensityState()

  withSharedPeptideIntensityUiMocks(
    captured,
    input_values = list(
      apply_intensity_filter = FALSE,
      revert_intensity = TRUE,
      use_strict_mode = FALSE,
      min_reps_per_group = 2,
      min_groups = 3,
      intensity_cutoff_percentile = 1.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideIntensityWorkflow(
    current_state,
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "intensity",
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
