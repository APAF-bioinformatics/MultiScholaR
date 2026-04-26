# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinReplicateState")) {
  methods::setClass(
    "FakeSharedProteinReplicateState",
    slots = c(protein_quant_table = "data.frame")
  )
}

getProtProteinReplicateServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_protein_replicate_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_protein_replicate_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_replicate_server
}

makeSharedProteinReplicateState <- function(proteins = c("P0", "P0")) {
  methods::new(
    "FakeSharedProteinReplicateState",
    protein_quant_table = data.frame(
      Protein.Ids = proteins,
      sample_a = seq_along(proteins),
      stringsAsFactors = FALSE
    )
  )
}

makeSharedProteinReplicateWorkflow <- function(current_state,
                                               captured,
                                               history = c("sample_filtered", "protein_replicate_filtered")) {
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

newSharedProteinReplicateCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$log_error <- character()
  captured
}

localSharedProteinReplicateBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedProteinReplicatePackageMocks <- function(server_env,
                                                   captured,
                                                   filtered_state,
                                                   filter_error = NULL) {
  localSharedProteinReplicateBinding(
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
    removeProteinsWithOnlyOneReplicate = function(theObject,
                                                  core_utilisation,
                                                  grouping_variable) {
      captured$filter_input <- list(
        theObject = theObject,
        core_utilisation = core_utilisation,
        grouping_variable = grouping_variable
      )
      if (!is.null(filter_error)) {
        stop(filter_error)
      }
      filtered_state
    },
    .env = server_env
  )
}

withSharedProteinReplicateUiMocks <- function(server_env,
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

test_that("proteomics protein replicate module preserves successful apply behavior", {
  captured <- newSharedProteinReplicateCapture()
  server_fn <- getProtProteinReplicateServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinReplicateState()
  filtered_state <- makeSharedProteinReplicateState(proteins = c("P1", "P1", "P2"))
  protein_qc_dir <- tempfile("protein-qc-")
  source_dir <- tempfile("source-")
  dir.create(protein_qc_dir)
  dir.create(source_dir)

  withr::defer(unlink(c(protein_qc_dir, source_dir), recursive = TRUE))
  withSharedProteinReplicatePackageMocks(server_env, captured, filtered_state)
  withSharedProteinReplicateUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_replicate_filter = TRUE,
      revert_protein_replicate_filter = FALSE,
      protein_grouping_variable = "group",
      parallel_cores = 3
    )
  )

  workflow_data <- makeSharedProteinReplicateWorkflow(current_state, captured)
  server_fn(
    id = "protein_replicate_filter",
    workflow_data = workflow_data,
    experiment_paths = list(
      protein_qc_dir = protein_qc_dir,
      source_dir = source_dir
    ),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$filter_input$theObject, current_state)
  expect_identical(captured$filter_input$grouping_variable, "group")
  expect_identical(workflow_data$qc_params$protein_qc$replicate_filter$grouping_variable, "group")
  expect_identical(workflow_data$qc_params$protein_qc$replicate_filter$parallel_cores, 3)
  expect_s3_class(workflow_data$qc_params$protein_qc$replicate_filter$timestamp, "POSIXct")
  expect_identical(captured$save_state$state_name, "protein_replicate_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_identical(captured$save_state$config_object$grouping_variable, "group")
  expect_identical(captured$save_state$config_object$parallel_cores, 3)
  expect_identical(
    captured$save_state$config_object$output_file,
    file.path(protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
  )
  expect_true(file.exists(captured$save_state$config_object$output_file))
  expect_true(file.exists(file.path(source_dir, "qc_params.RDS")))
  expect_identical(
    captured$save_state$description,
    "Applied protein replicate filter (removed single-replicate proteins)"
  )
  expect_identical(workflow_data$protein_counts$after_qc_filtering, 2L)
  expect_identical(captured$plot_update$step_name, "13_protein_replicate_filtered")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$protein_replicate_filter_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$protein_replicate_filter_results, "Grouping variable: group", fixed = TRUE)
  expect_match(captured$output$protein_replicate_filter_results, "Parallel cores used: 3", fixed = TRUE)
  expect_identical(captured$output$protein_replicate_filter_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "protein_replicate_filter_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein replicate module preserves fallback output and qc-params warning behavior", {
  captured <- newSharedProteinReplicateCapture()
  server_fn <- getProtProteinReplicateServer()
  server_env <- environment(server_fn)
  filtered_state <- makeSharedProteinReplicateState(proteins = c("P1", "P2"))
  bad_source <- tempfile("source-file-")
  writeLines("not a directory", bad_source)
  old_dir <- getwd()

  withr::defer(setwd(old_dir))
  withr::defer(unlink(bad_source))
  setwd(tempdir())
  withSharedProteinReplicatePackageMocks(server_env, captured, filtered_state)
  withSharedProteinReplicateUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_replicate_filter = TRUE,
      revert_protein_replicate_filter = FALSE,
      protein_grouping_variable = "condition",
      parallel_cores = 2
    )
  )

  workflow_data <- makeSharedProteinReplicateWorkflow(
    makeSharedProteinReplicateState(),
    captured
  )
  suppressWarnings(
    server_fn(
      id = "protein_replicate_filter",
      workflow_data = workflow_data,
      experiment_paths = list(
        protein_qc_dir = NULL,
        source_dir = bad_source
      ),
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    )
  )

  expect_identical(
    captured$save_state$config_object$output_file,
    "remove_proteins_with_only_one_rep.tsv"
  )
  expect_true(any(grepl("Could not save QC parameters file", captured$log_warn, fixed = TRUE)))
  expect_match(captured$output$protein_replicate_filter_results, "Output file: remove_proteins_with_only_one_rep.tsv", fixed = TRUE)
})

test_that("proteomics protein replicate module preserves apply error behavior", {
  captured <- newSharedProteinReplicateCapture()
  server_fn <- getProtProteinReplicateServer()
  server_env <- environment(server_fn)

  withSharedProteinReplicatePackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedProteinReplicateState(),
    filter_error = "mock replicate failure"
  )
  withSharedProteinReplicateUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_replicate_filter = TRUE,
      revert_protein_replicate_filter = FALSE,
      protein_grouping_variable = "group",
      parallel_cores = 3
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinReplicateWorkflow(
    makeSharedProteinReplicateState(),
    captured
  )
  server_fn(
    id = "protein_replicate_filter",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir(), source_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock replicate failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock replicate failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "protein_replicate_filter_working")
})

test_that("proteomics protein replicate module preserves revert behavior", {
  captured <- newSharedProteinReplicateCapture()
  server_fn <- getProtProteinReplicateServer()
  server_env <- environment(server_fn)

  withSharedProteinReplicateUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_replicate_filter = FALSE,
      revert_protein_replicate_filter = TRUE,
      protein_grouping_variable = "group",
      parallel_cores = 3
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinReplicateWorkflow(
    makeSharedProteinReplicateState(),
    captured,
    history = c("sample_filtered", "protein_replicate_filtered")
  )
  server_fn(
    id = "protein_replicate_filter",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir(), source_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "sample_filtered")
  expect_identical(
    captured$output$protein_replicate_filter_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein replicate module preserves revert error behavior", {
  captured <- newSharedProteinReplicateCapture()
  server_fn <- getProtProteinReplicateServer()
  server_env <- environment(server_fn)

  withSharedProteinReplicateUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_protein_replicate_filter = FALSE,
      revert_protein_replicate_filter = TRUE,
      protein_grouping_variable = "group",
      parallel_cores = 3
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinReplicateWorkflow(
    makeSharedProteinReplicateState(),
    captured,
    history = "protein_replicate_filtered"
  )
  server_fn(
    id = "protein_replicate_filter",
    workflow_data = workflow_data,
    experiment_paths = list(protein_qc_dir = tempdir(), source_dir = tempdir()),
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
