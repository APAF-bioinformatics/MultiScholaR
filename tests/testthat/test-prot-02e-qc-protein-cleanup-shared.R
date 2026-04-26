# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinCleanupState")) {
  methods::setClass(
    "FakeSharedProteinCleanupState",
    slots = c(protein_quant_table = "data.frame")
  )
}

getProtProteinCleanupServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_protein_cleanup_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_protein_cleanup_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_cleanup_server
}

makeSharedProteinCleanupState <- function(proteins = "P0") {
  methods::new(
    "FakeSharedProteinCleanupState",
    protein_quant_table = data.frame(
      Protein.Ids = proteins,
      stringsAsFactors = FALSE
    )
  )
}

makeSharedProteinCleanupWorkflow <- function(current_state,
                                             captured,
                                             history = c("raw_data_s4", "protein_accession_cleaned"),
                                             fasta_metadata = list(
                                               has_protein_evidence = TRUE,
                                               has_gene_names = TRUE
                                             )) {
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
  workflow_data$fasta_metadata <- fasta_metadata
  workflow_data
}

setSharedProteinCleanupGlobal <- function(name, value) {
  old_exists <- exists(name, envir = .GlobalEnv, inherits = FALSE)
  old_value <- if (old_exists) get(name, envir = .GlobalEnv, inherits = FALSE) else NULL

  assign(name, value, envir = .GlobalEnv)
  withr::defer({
    if (old_exists) {
      assign(name, old_value, envir = .GlobalEnv)
    } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }, envir = testthat::teardown_env())
}

withoutSharedProteinCleanupGlobal <- function(name) {
  old_exists <- exists(name, envir = .GlobalEnv, inherits = FALSE)
  old_value <- if (old_exists) get(name, envir = .GlobalEnv, inherits = FALSE) else NULL

  if (old_exists) {
    rm(list = name, envir = .GlobalEnv)
  }
  withr::defer({
    if (old_exists) {
      assign(name, old_value, envir = .GlobalEnv)
    } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = name, envir = .GlobalEnv)
    }
  }, envir = testthat::teardown_env())
}

localSharedProteinCleanupBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedProteinCleanupPackageMocks <- function(server_env, captured, cleaned_state) {
  localSharedProteinCleanupBinding(
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
    chooseBestProteinAccession = function(theObject,
                                          delim,
                                          seqinr_obj,
                                          seqinr_accession_column,
                                          replace_zero_with_na,
                                          aggregation_method) {
      captured$cleanup_call <- list(
        theObject = theObject,
        delim = delim,
        seqinr_obj = seqinr_obj,
        seqinr_accession_column = seqinr_accession_column,
        replace_zero_with_na = replace_zero_with_na,
        aggregation_method = aggregation_method
      )
      cleaned_state
    },
    .env = server_env
  )
}

withSharedProteinCleanupUiMocks <- function(server_env,
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
    log_warn = function(message) {
      captured$log_warn <- c(captured$log_warn, message)
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

newSharedProteinCleanupCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured$log_warn <- character()
  captured
}

test_that("proteomics protein cleanup module preserves accession-cleanup apply behavior", {
  captured <- newSharedProteinCleanupCapture()
  server_fn <- getProtProteinCleanupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinCleanupState(c("P1;P2", "P1;P2", "P3"))
  cleaned_state <- makeSharedProteinCleanupState(c("P1", "P1"))
  source_dir <- tempfile("protein-cleanup-shared-")
  dir.create(source_dir)
  withr::defer(unlink(source_dir, recursive = TRUE, force = TRUE), envir = teardown_env())

  setSharedProteinCleanupGlobal(
    "aa_seq_tbl_final",
    data.frame(database_id = c("P1", "P2"), gene = c("G1", "G2"))
  )
  setSharedProteinCleanupGlobal("experiment_paths", list(source_dir = source_dir))
  withSharedProteinCleanupPackageMocks(server_env, captured, cleaned_state)
  withSharedProteinCleanupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_accession_cleanup = TRUE,
      revert_accession_cleanup = FALSE,
      delimiter = ";",
      aggregation_method = "median"
    )
  )

  workflow_data <- makeSharedProteinCleanupWorkflow(current_state, captured)
  server_fn(
    id = "cleanup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$cleanup_call$theObject, current_state)
  expect_identical(captured$cleanup_call$delim, ";")
  expect_identical(captured$cleanup_call$seqinr_accession_column, "uniprot_acc")
  expect_true(captured$cleanup_call$replace_zero_with_na)
  expect_identical(captured$cleanup_call$aggregation_method, "median")
  expect_true("uniprot_acc" %in% colnames(captured$cleanup_call$seqinr_obj))
  expect_identical(workflow_data$accession_cleanup_results$cleanup_applied, TRUE)
  expect_identical(workflow_data$accession_cleanup_results$delimiter_used, ";")
  expect_identical(workflow_data$accession_cleanup_results$aggregation_method, "median")
  expect_identical(workflow_data$accession_cleanup_results$proteins_before, 2L)
  expect_identical(workflow_data$accession_cleanup_results$proteins_after, 1L)
  expect_identical(workflow_data$accession_cleanup_results$had_full_metadata, TRUE)
  expect_identical(
    workflow_data$qc_params$protein_qc$accession_cleanup,
    workflow_data$accession_cleanup_results
  )
  expect_true(file.exists(file.path(source_dir, "accession_cleanup_results.RDS")))
  expect_identical(captured$save_state$state_name, "protein_accession_cleaned")
  expect_identical(captured$save_state$s4_data_object, cleaned_state)
  expect_identical(captured$save_state$config_object$delimiter, ";")
  expect_identical(captured$save_state$config_object$aggregation_method, "median")
  expect_true(captured$save_state$config_object$cleanup_applied)
  expect_identical(captured$save_state$description, "Applied protein accession cleanup")
  expect_match(
    captured$output$accession_cleanup_results,
    "Protein Accession Cleanup Applied Successfully",
    fixed = TRUE
  )
  expect_match(captured$output$accession_cleanup_results, "Proteins remaining: 1", fixed = TRUE)
  expect_identical(captured$plot_update$step_name, "10_protein_accession_cleaned")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_identical(captured$output$accession_cleanup_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "accession_cleanup_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein cleanup module preserves no-FASTA skip behavior", {
  captured <- newSharedProteinCleanupCapture()
  server_fn <- getProtProteinCleanupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinCleanupState(c("P1;P2", "P3"))

  withoutSharedProteinCleanupGlobal("aa_seq_tbl_final")
  withSharedProteinCleanupPackageMocks(server_env, captured, cleaned_state = current_state)
  withSharedProteinCleanupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_accession_cleanup = TRUE,
      revert_accession_cleanup = FALSE,
      delimiter = ";",
      aggregation_method = "mean"
    )
  )

  workflow_data <- makeSharedProteinCleanupWorkflow(current_state, captured, fasta_metadata = NULL)
  server_fn(
    id = "cleanup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("cleanup_call", envir = captured, inherits = FALSE))
  expect_identical(workflow_data$accession_cleanup_results$cleanup_applied, FALSE)
  expect_identical(workflow_data$accession_cleanup_results$had_full_metadata, FALSE)
  expect_identical(captured$save_state$s4_data_object, current_state)
  expect_false(captured$save_state$config_object$cleanup_applied)
  expect_identical(captured$save_state$description, "Skipped accession cleanup (no FASTA data)")
  expect_match(
    captured$output$accession_cleanup_results,
    "Protein Accession Cleanup Skipped",
    fixed = TRUE
  )
})

test_that("proteomics protein cleanup module preserves apply error behavior", {
  captured <- newSharedProteinCleanupCapture()
  server_fn <- getProtProteinCleanupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinCleanupState(c("P1;P2", "P3"))

  setSharedProteinCleanupGlobal(
    "aa_seq_tbl_final",
    data.frame(database_id = c("P1", "P2"), gene = c("G1", "G2"))
  )
  withSharedProteinCleanupPackageMocks(
    server_env,
    captured,
    cleaned_state = makeSharedProteinCleanupState()
  )
  testthat::local_mocked_bindings(
    chooseBestProteinAccession = function(...) {
      stop("mock cleanup failure")
    },
    .env = server_env
  )
  withSharedProteinCleanupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_accession_cleanup = TRUE,
      revert_accession_cleanup = FALSE,
      delimiter = ";",
      aggregation_method = "median"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinCleanupWorkflow(current_state, captured)
  server_fn(
    id = "cleanup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("mock cleanup failure", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("mock cleanup failure", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "accession_cleanup_working")
})

test_that("proteomics protein cleanup module preserves revert behavior", {
  captured <- newSharedProteinCleanupCapture()
  server_fn <- getProtProteinCleanupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinCleanupState()

  withSharedProteinCleanupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_accession_cleanup = FALSE,
      revert_accession_cleanup = TRUE,
      delimiter = ";",
      aggregation_method = "mean"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinCleanupWorkflow(
    current_state,
    captured,
    history = c("raw_data_s4", "rollup_complete", "protein_accession_cleaned")
  )
  server_fn(
    id = "cleanup",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "rollup_complete")
  expect_match(
    captured$output$accession_cleanup_results,
    "Reverted to previous state: rollup_complete",
    fixed = TRUE
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein cleanup module preserves revert error behavior", {
  captured <- newSharedProteinCleanupCapture()
  server_fn <- getProtProteinCleanupServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedProteinCleanupState()

  withSharedProteinCleanupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_accession_cleanup = FALSE,
      revert_accession_cleanup = TRUE,
      delimiter = ";",
      aggregation_method = "mean"
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinCleanupWorkflow(
    current_state,
    captured,
    history = "raw_data_s4"
  )
  server_fn(
    id = "cleanup",
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
