# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinDedupState")) {
  methods::setClass(
    "FakeSharedProteinDedupState",
    slots = c(protein_quant_table = "data.frame")
  )
}

getProtProteinDedupServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_protein_dedup_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_protein_dedup_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_dedup_server
}

makeSharedProteinDedupState <- function(proteins = c("P1", "P1", "P2"),
                                        sample1 = c(2, 4, 8),
                                        sample2 = c(10, 14, 18)) {
  methods::new(
    "FakeSharedProteinDedupState",
    protein_quant_table = data.frame(
      Protein.Ids = proteins,
      sample1 = sample1,
      sample2 = sample2,
      stringsAsFactors = FALSE
    )
  )
}

makeSharedProteinDedupWorkflow <- function(current_state,
                                           captured,
                                           history = c("raw_data_s4", "duplicates_removed")) {
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

newSharedProteinDedupCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedProteinDedupBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

localSharedProteinDedupBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localSharedProteinDedupBinding(
      env,
      name,
      bindings[[name]],
      .local_envir = .local_envir
    )
  }

  invisible(NULL)
}

withSharedProteinDedupPackageMocks <- function(server_env, captured) {
  localSharedProteinDedupBindings(
    server_env,
    list(
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

withSharedProteinDedupUiMocks <- function(captured,
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

test_that("proteomics protein dedup module preserves successful apply behavior", {
  captured <- newSharedProteinDedupCapture()
  server_fn <- getProtProteinDedupServer()
  server_env <- environment(server_fn)

  withSharedProteinDedupPackageMocks(server_env, captured)
  withSharedProteinDedupUiMocks(
    captured,
    input_values = list(
      apply_duplicate_removal = TRUE,
      revert_duplicate_removal = FALSE,
      duplicate_aggregation_method = "max"
    )
  )

  workflow_data <- makeSharedProteinDedupWorkflow(
    makeSharedProteinDedupState(),
    captured
  )
  server_fn(
    id = "protein_dedup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  deduped_table <- as.data.frame(captured$save_state$s4_data_object@protein_quant_table)
  expect_identical(captured$save_state$state_name, "duplicates_removed")
  expect_identical(captured$save_state$config_object$aggregation_method, "max")
  expect_identical(captured$save_state$config_object$duplicates_found, "P1")
  expect_identical(captured$save_state$config_object$num_duplicates, 1L)
  expect_identical(
    captured$save_state$description,
    "Removed duplicate proteins by aggregation"
  )
  expect_identical(deduped_table$Protein.Ids, c("P1", "P2"))
  expect_identical(deduped_table$sample1, c(4, 8))
  expect_identical(deduped_table$sample2, c(14, 18))
  expect_identical(captured$plot_update$step_name, "12_duplicates_removed")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(captured$output$duplicate_removal_results, "Proteins remaining: 2", fixed = TRUE)
  expect_match(captured$output$duplicate_removal_results, "Duplicates found: 1", fixed = TRUE)
  expect_match(captured$output$duplicate_removal_results, "Aggregation method: max", fixed = TRUE)
  expect_identical(captured$output$duplicate_removal_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "duplicate_removal_working")
  expect_true(any(grepl("Duplicate protein removal completed successfully", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein dedup module preserves apply error behavior", {
  captured <- newSharedProteinDedupCapture()
  server_fn <- getProtProteinDedupServer()
  server_env <- environment(server_fn)

  withSharedProteinDedupPackageMocks(server_env, captured)
  withSharedProteinDedupUiMocks(
    captured,
    input_values = list(
      apply_duplicate_removal = TRUE,
      revert_duplicate_removal = FALSE,
      duplicate_aggregation_method = "missing_aggregation_function"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinDedupWorkflow(
    makeSharedProteinDedupState(),
    captured
  )
  server_fn(
    id = "protein_dedup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("Error removing duplicate proteins:", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("Error removing duplicate proteins:", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "duplicate_removal_working")
})

test_that("proteomics protein dedup module preserves revert behavior", {
  captured <- newSharedProteinDedupCapture()
  server_fn <- getProtProteinDedupServer()
  server_env <- environment(server_fn)

  withSharedProteinDedupPackageMocks(server_env, captured)
  withSharedProteinDedupUiMocks(
    captured,
    input_values = list(
      apply_duplicate_removal = FALSE,
      revert_duplicate_removal = TRUE,
      duplicate_aggregation_method = "mean"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinDedupWorkflow(
    makeSharedProteinDedupState(),
    captured,
    history = c("raw_data_s4", "duplicates_removed")
  )
  server_fn(
    id = "protein_dedup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "raw_data_s4")
  expect_identical(
    captured$output$duplicate_removal_results,
    "Reverted to previous state: raw_data_s4"
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein dedup module preserves revert error behavior", {
  captured <- newSharedProteinDedupCapture()
  server_fn <- getProtProteinDedupServer()
  server_env <- environment(server_fn)

  withSharedProteinDedupPackageMocks(server_env, captured)
  withSharedProteinDedupUiMocks(
    captured,
    input_values = list(
      apply_duplicate_removal = FALSE,
      revert_duplicate_removal = TRUE,
      duplicate_aggregation_method = "mean"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinDedupWorkflow(
    makeSharedProteinDedupState(),
    captured,
    history = "duplicates_removed"
  )
  server_fn(
    id = "protein_dedup",
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
