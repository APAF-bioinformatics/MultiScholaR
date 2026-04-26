# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedPeptideQvalueState")) {
  methods::setClass(
    "FakeSharedPeptideQvalueState",
    slots = c(args = "list", peptide_data = "data.frame")
  )
}

getProtPeptideQvalueServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_peptide_qvalue_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_peptide_qvalue_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_peptide_qvalue_server
}

makeSharedPeptideQvalueState <- function(ids = c(" sample1", "", "sample3 "),
                                         qvalue = 0.05,
                                         global_qvalue = 0.05,
                                         proteotypic = 0,
                                         proteins = c("P0", "P0")) {
  methods::new(
    "FakeSharedPeptideQvalueState",
    args = list(
      srlQvalueProteotypicPeptideClean = list(
        input_matrix_column_ids = ids,
        qvalue_threshold = qvalue,
        global_qvalue_threshold = global_qvalue,
        choose_only_proteotypic_peptide = proteotypic
      )
    ),
    peptide_data = data.frame(
      Protein.Ids = proteins,
      Intensity = seq_along(proteins),
      stringsAsFactors = FALSE
    )
  )
}

copySharedPeptideQvalueState <- function(theObject,
                                         args = theObject@args,
                                         peptide_data = theObject@peptide_data) {
  methods::new(
    "FakeSharedPeptideQvalueState",
    args = args,
    peptide_data = peptide_data
  )
}

makeSharedPeptideQvalueWorkflow <- function(current_state,
                                            captured,
                                            history = c("raw_data_s4", "qvalue_filtered")) {
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

newSharedPeptideQvalueCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$config_updates <- list()
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$log_error <- character()
  captured
}

localSharedPeptideQvalueBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedPeptideQvaluePackageMocks <- function(server_env,
                                                captured,
                                                filtered_state,
                                                filter_error = NULL) {
  localSharedPeptideQvalueBinding(
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
      captured$config_updates[[length(captured$config_updates) + 1L]] <- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )

      updated_args <- theObject@args
      updated_args[[function_name]][[parameter_name]] <- new_value
      copySharedPeptideQvalueState(theObject, args = updated_args)
    },
    srlQvalueProteotypicPeptideClean = function(theObject) {
      captured$filter_input <- theObject
      if (!is.null(filter_error)) {
        stop(filter_error)
      }
      filtered_state
    },
    .env = server_env
  )
}

withSharedPeptideQvalueUiMocks <- function(server_env,
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

test_that("proteomics peptide qvalue module preserves successful apply behavior", {
  captured <- newSharedPeptideQvalueCapture()
  server_fn <- getProtPeptideQvalueServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedPeptideQvalueState()
  filtered_state <- makeSharedPeptideQvalueState(
    ids = c("Run", "Precursor.Id", "Intensity"),
    qvalue = 0.01,
    global_qvalue = 0.02,
    proteotypic = 1,
    proteins = c("P1", "P1", "P2")
  )

  withSharedPeptideQvaluePackageMocks(server_env, captured, filtered_state)
  withSharedPeptideQvalueUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_qvalue_filter = TRUE,
      revert_qvalue = FALSE,
      qvalue_threshold = 0.01,
      global_qvalue_threshold = 0.02,
      proteotypic_only = TRUE
    )
  )

  workflow_data <- makeSharedPeptideQvalueWorkflow(current_state, captured)
  server_fn(
    id = "qvalue_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_true(any(grepl("S4 class = FakeSharedPeptideQvalueState", captured$log_info, fixed = TRUE)))
  expect_true(any(grepl("input_matrix_column_ids length = 3", captured$log_info, fixed = TRUE)))
  expect_identical(
    captured$log_warn,
    c(
      "Q-value filter: input_matrix_column_ids contains values with leading/trailing whitespace!",
      "Q-value filter: input_matrix_column_ids contains empty strings!"
    )
  )
  expect_identical(
    vapply(captured$config_updates, `[[`, character(1), "parameter_name"),
    c(
      "qvalue_threshold",
      "global_qvalue_threshold",
      "choose_only_proteotypic_peptide"
    )
  )
  expect_equal(
    vapply(captured$config_updates, `[[`, numeric(1), "new_value"),
    c(0.01, 0.02, 1)
  )
  expect_identical(
    captured$filter_input@args$srlQvalueProteotypicPeptideClean$qvalue_threshold,
    0.01
  )
  expect_identical(
    captured$filter_input@args$srlQvalueProteotypicPeptideClean$global_qvalue_threshold,
    0.02
  )
  expect_identical(
    captured$filter_input@args$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide,
    1
  )
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$qvalue_threshold, 0.01)
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$global_qvalue_threshold, 0.02)
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$proteotypic_only, TRUE)
  expect_s3_class(workflow_data$qc_params$peptide_qc$qvalue_filter$timestamp, "POSIXct")
  expect_identical(captured$save_state$state_name, "qvalue_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_identical(captured$save_state$config_object$qvalue_threshold, 0.01)
  expect_identical(captured$save_state$config_object$global_qvalue_threshold, 0.02)
  expect_identical(captured$save_state$config_object$proteotypic_only, TRUE)
  expect_identical(
    captured$save_state$description,
    "Applied Q-value and proteotypic peptide filter"
  )
  expect_identical(captured$plot_update$step_name, "2_qval_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$qvalue_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$qvalue_results, "Q-value threshold: 0.01", fixed = TRUE)
  expect_match(captured$output$qvalue_results, "Global Q-value threshold: 0.02", fixed = TRUE)
  expect_match(captured$output$qvalue_results, "Proteotypic only: TRUE", fixed = TRUE)
  expect_identical(captured$output$qvalue_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "qvalue_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide qvalue module preserves apply error behavior", {
  captured <- newSharedPeptideQvalueCapture()
  server_fn <- getProtPeptideQvalueServer()
  server_env <- environment(server_fn)

  withSharedPeptideQvaluePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedPeptideQvalueState(),
    filter_error = "mock qvalue failure"
  )
  withSharedPeptideQvalueUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_qvalue_filter = TRUE,
      revert_qvalue = FALSE,
      qvalue_threshold = 0.01,
      global_qvalue_threshold = 0.02,
      proteotypic_only = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideQvalueWorkflow(
    makeSharedPeptideQvalueState(),
    captured
  )
  server_fn(
    id = "qvalue_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock qvalue failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock qvalue failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "qvalue_working")
})

test_that("proteomics peptide qvalue module preserves revert behavior", {
  captured <- newSharedPeptideQvalueCapture()
  server_fn <- getProtPeptideQvalueServer()
  server_env <- environment(server_fn)

  withSharedPeptideQvalueUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_qvalue_filter = FALSE,
      revert_qvalue = TRUE,
      qvalue_threshold = 0.01,
      global_qvalue_threshold = 0.02,
      proteotypic_only = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideQvalueWorkflow(
    makeSharedPeptideQvalueState(),
    captured,
    history = c("raw_data_s4", "qvalue_filtered")
  )
  server_fn(
    id = "qvalue_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "raw_data_s4")
  expect_identical(captured$output$qvalue_results, "Reverted to raw data state")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics peptide qvalue module preserves revert error behavior", {
  captured <- newSharedPeptideQvalueCapture()
  server_fn <- getProtPeptideQvalueServer()
  server_env <- environment(server_fn)

  withSharedPeptideQvalueUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_qvalue_filter = FALSE,
      revert_qvalue = TRUE,
      qvalue_threshold = 0.01,
      global_qvalue_threshold = 0.02,
      proteotypic_only = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedPeptideQvalueWorkflow(
    makeSharedPeptideQvalueState(),
    captured,
    history = "qvalue_filtered"
  )
  server_fn(
    id = "qvalue_filter",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("reverted_state", envir = captured, inherits = FALSE))
  expect_true(any(grepl("Cannot revert", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("Cannot revert", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
})
