# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideSampleState")) {
  methods::setClass(
    "FakeSharedPeptideSampleState",
    slots = c(args = "list", peptide_data = "data.frame", sample_id = "character")
  )
}

getProtPeptideSampleServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_sample_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_sample_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_sample_server
}

makeSharedPeptideSampleState <- function(runs = c("S1", "S1", "S2", "S3"),
                                         proteins = c("P1", "P1", "P2", "P3"),
                                         args = list()) {
  methods::new(
    "FakeSharedPeptideSampleState",
    args = args,
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Run = runs,
      stringsAsFactors = FALSE
    ),
    sample_id = "Run"
  )
}

makeSharedPeptideSampleWorkflow <- function(current_state,
                                            captured,
                                            history = c("raw_data_s4", "sample_filtered")) {
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

newSharedPeptideSampleCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedPeptideSampleBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedPeptideSamplePackageMocks <- function(server_env,
                                                captured,
                                                filtered_state,
                                                filter_error = NULL) {
  localSharedPeptideSampleBinding(
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
    updateConfigParameter = function(theObject, function_name, parameter_name, new_value) {
      captured$config_update <- list(
        theObject = theObject,
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      theObject
    },
    filterMinNumPeptidesPerSample = function(theObject) {
      captured$filter_input <- theObject
      if (!is.null(filter_error)) {
        stop(filter_error)
      }
      filtered_state
    },
    .env = server_env
  )
}

withSharedPeptideSampleUiMocks <- function(server_env,
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

test_that("proteomics peptide sample module preserves successful apply behavior", {
  captured <- newSharedPeptideSampleCapture()
  server_fn <- getProtPeptideSampleServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideSampleState()
  filtered_state <- makeSharedPeptideSampleState(
    runs = c("S1", "S1", "S3"),
    proteins = c("P1", "P1", "P3")
  )

  withSharedPeptideSamplePackageMocks(server_env, captured, filtered_state)
  withSharedPeptideSampleUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_sample_filter = TRUE,
      revert_sample = FALSE,
      min_peptides_per_sample = 500
    )
  )

  workflow_data <- makeSharedPeptideSampleWorkflow(current_state, captured)
  server_fn(
    id = "sample_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$config_update$theObject, current_state)
  expect_identical(captured$config_update$function_name, "filterMinNumPeptidesPerSample")
  expect_identical(captured$config_update$parameter_name, "peptides_per_sample_cutoff")
  expect_identical(captured$config_update$new_value, 500)
  expect_identical(captured$filter_input, current_state)
  expect_identical(
    captured$save_state$state_name,
    "sample_filtered"
  )
  expect_identical(
    captured$save_state$s4_data_object@args$filterMinNumPeptidesPerSample$samples_removed,
    "S2"
  )
  expect_identical(
    captured$save_state$s4_data_object@args$filterMinNumPeptidesPerSample$samples_removed_count,
    1L
  )
  expect_identical(
    captured$save_state$s4_data_object@args$filterMinNumPeptidesPerSample$samples_before_count,
    3L
  )
  expect_identical(
    captured$save_state$s4_data_object@args$filterMinNumPeptidesPerSample$samples_after_count,
    2L
  )
  expect_identical(captured$save_state$config_object$min_peptides_per_sample, 500)
  expect_identical(captured$save_state$config_object$samples_removed, "S2")
  expect_identical(captured$save_state$config_object$samples_removed_count, 1L)
  expect_identical(
    captured$save_state$description,
    "Applied minimum peptides per sample filter"
  )
  expect_identical(workflow_data$qc_params$peptide_qc$sample_filter$min_peptides_per_sample, 500)
  expect_identical(workflow_data$qc_params$peptide_qc$sample_filter$samples_removed, "S2")
  expect_identical(workflow_data$qc_params$peptide_qc$sample_filter$samples_removed_count, 1L)
  expect_identical(workflow_data$qc_params$peptide_qc$sample_filter$samples_before_count, 3L)
  expect_identical(workflow_data$qc_params$peptide_qc$sample_filter$samples_after_count, 2L)
  expect_s3_class(workflow_data$qc_params$peptide_qc$sample_filter$timestamp, "POSIXct")
  expect_identical(captured$plot_update$step_name, "6_sample_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$sample_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$sample_results, "Samples remaining: 2", fixed = TRUE)
  expect_match(captured$output$sample_results, "Samples removed: 1", fixed = TRUE)
  expect_match(captured$output$sample_results, "Min peptides per sample: 500", fixed = TRUE)
  expect_match(captured$output$sample_results, "Removed samples:\nS2", fixed = TRUE)
  expect_identical(captured$output$sample_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "sample_working")
  expect_true(any(grepl("Applying sample quality filter", captured$log_info, fixed = TRUE)))
  expect_true(any(grepl("applied successfully", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide sample module preserves apply error behavior", {
  captured <- newSharedPeptideSampleCapture()
  server_fn <- getProtPeptideSampleServer()
  server_env <- environment(server_fn)

  withSharedPeptideSamplePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideSampleState(),
    filter_error = "mock sample failure"
  )
  withSharedPeptideSampleUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_sample_filter = TRUE,
      revert_sample = FALSE,
      min_peptides_per_sample = 500
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideSampleWorkflow(
    makeSharedPeptideSampleState(),
    captured
  )
  server_fn(
    id = "sample_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock sample failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock sample failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "sample_working")
})

test_that("proteomics peptide sample module preserves revert behavior", {
  captured <- newSharedPeptideSampleCapture()
  server_fn <- getProtPeptideSampleServer()
  server_env <- environment(server_fn)

  withSharedPeptideSampleUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_sample_filter = FALSE,
      revert_sample = TRUE,
      min_peptides_per_sample = 500
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideSampleWorkflow(
    makeSharedPeptideSampleState(),
    captured,
    history = c("raw_data_s4", "sample_filtered", "replicate_filtered")
  )
  server_fn(
    id = "sample_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "sample_filtered")
  expect_identical(
    captured$output$sample_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide sample module preserves revert error behavior", {
  captured <- newSharedPeptideSampleCapture()
  server_fn <- getProtPeptideSampleServer()
  server_env <- environment(server_fn)

  withSharedPeptideSampleUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_sample_filter = FALSE,
      revert_sample = TRUE,
      min_peptides_per_sample = 500
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideSampleWorkflow(
    makeSharedPeptideSampleState(),
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "sample_filter",
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
