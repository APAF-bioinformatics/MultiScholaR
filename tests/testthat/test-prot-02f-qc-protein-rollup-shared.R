# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinRollupPeptideState")) {
  methods::setClass(
    "FakeSharedProteinRollupPeptideState",
    slots = c(
      peptide_data = "data.frame",
      design_matrix = "data.frame",
      args = "list"
    )
  )
}

if (!methods::isClass("FakeSharedProteinRollupProteinState")) {
  methods::setClass(
    "FakeSharedProteinRollupProteinState",
    slots = c(
      protein_quant_table = "data.frame",
      design_matrix = "data.frame",
      args = "list"
    )
  )
}

getProtProteinRollupServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists(
    "mod_prot_qc_protein_rollup_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    return(get("mod_prot_qc_protein_rollup_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_rollup_server
}

makeSharedProteinRollupFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

makeSharedProteinRollupPeptideState <- function() {
  methods::new(
    "FakeSharedProteinRollupPeptideState",
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1", "P2"),
      Stripped.Sequence = c("AAA", "BBB", "CCC"),
      Run = c("sample_a", "sample_b", "sample_a"),
      Peptide.Imputed = c(10, NA, 20),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("sample_a", "sample_b"),
      group = c("case", "control"),
      replicates = c(1, 1),
      stringsAsFactors = FALSE
    ),
    args = list(source = "shared-rollup-test")
  )
}

makeSharedProteinRollupProteinState <- function(protein_quant_table,
                                                design_matrix,
                                                args) {
  methods::new(
    "FakeSharedProteinRollupProteinState",
    protein_quant_table = protein_quant_table,
    design_matrix = design_matrix,
    args = args
  )
}

makeSharedProteinRollupWorkflow <- function(current_state,
                                            captured,
                                            history = c("raw_data_s4", "qvalue_filtered", "protein_s4_created")) {
  state_manager <- list(
    current_state = "imputed",
    getState = function(state_name) {
      captured$get_state <- state_name
      current_state
    },
    saveState = function(...) {
      captured$save_state <- list(...)
      invisible(NULL)
    },
    getHistory = function() history,
    revertToState = function(state_name) {
      captured$reverted_state <- state_name
      if (is.na(state_name)) {
        stop("no previous peptide state")
      }
      current_state
    }
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data
}

newSharedProteinRollupCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$log_error <- character()
  captured
}

localSharedProteinRollupBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

withSharedProteinRollupPackageMocks <- function(server_env,
                                                captured,
                                                create_output = TRUE) {
  localSharedProteinRollupBinding(
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
    ProteinQuantitativeData = function(protein_quant_table,
                                       protein_id_column,
                                       protein_id_table,
                                       design_matrix,
                                       sample_id,
                                       group_id,
                                       technical_replicate_id,
                                       args) {
      captured$create <- list(
        protein_quant_table = protein_quant_table,
        protein_id_column = protein_id_column,
        protein_id_table = protein_id_table,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
      )
      makeSharedProteinRollupProteinState(
        protein_quant_table = protein_quant_table,
        design_matrix = design_matrix,
        args = args
      )
    },
    .capture_checkpoint = function(object, checkpointId, checkpointLabel) {
      captured$checkpoint <- list(
        object = object,
        checkpointId = checkpointId,
        checkpointLabel = checkpointLabel
      )
      invisible(NULL)
    },
    .env = server_env
  )
}

`%||%` <- function(left, right) {
  if (is.null(left)) {
    right
  } else {
    left
  }
}

withSharedProteinRollupUiMocks <- function(server_env,
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

test_that("proteomics protein rollup module preserves successful apply behavior", {
  captured <- newSharedProteinRollupCapture()
  server_fn <- getProtProteinRollupServer()
  server_env <- environment(server_fn)

  withSharedProteinRollupPackageMocks(server_env, captured)
  withSharedProteinRollupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_iq_rollup = TRUE,
      revert_iq_rollup = FALSE
    )
  )

  workflow_data <- makeSharedProteinRollupWorkflow(
    makeSharedProteinRollupPeptideState(),
    captured
  )
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())

  server_fn(
    id = "protein_rollup",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$get_state, "imputed")
  expect_identical(colnames(captured$create$protein_quant_table), c("Protein.Ids", "sample_a"))
  expect_identical(captured$create$protein_id_column, "Protein.Ids")
  expect_identical(captured$create$design_matrix$Run, "sample_a")
  expect_identical(captured$create$sample_id, "Run")
  expect_identical(captured$create$group_id, "group")
  expect_identical(captured$create$technical_replicate_id, "replicates")
  expect_identical(captured$create$args, list(source = "shared-rollup-test"))
  expect_identical(captured$save_state$state_name, "protein_s4_created")
  expect_identical(captured$save_state$s4_data_object, captured$checkpoint$object)
  expect_identical(captured$save_state$config_object$s4_class, "ProteinQuantitativeData")
  expect_identical(captured$checkpoint$checkpointId, "cp03")
  expect_identical(captured$checkpoint$checkpointLabel, "rolled_up_protein")
  expect_identical(captured$plot_update$step_name, "9_protein_s4_created")
  expect_identical(captured$plot_update$omic_type, "proteomics")
  expect_identical(captured$plot_update$experiment_label, "DIA Experiment")
  expect_true(captured$plot_update$return_grid)
  expect_true(captured$plot_update$overwrite)
  expect_match(
    captured$output$iq_rollup_results,
    "IQ Protein Rollup & S4 Object Creation Completed Successfully",
    fixed = TRUE
  )
  expect_match(captured$output$iq_rollup_results, "Proteins quantified: 2", fixed = TRUE)
  expect_match(captured$output$iq_rollup_results, "Samples: 1", fixed = TRUE)
  expect_identical(captured$output$iq_rollup_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, "plot-token")
  expect_identical(captured$removed_notifications, "iq_rollup_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "warning"),
    logical(1)
  )))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein rollup module preserves IQ timeout behavior", {
  captured <- newSharedProteinRollupCapture()
  server_fn <- getProtProteinRollupServer()
  server_env <- environment(server_fn)
  if (exists("runProteinIqRollupApplyStep", envir = server_env, inherits = FALSE)) {
    skip("Timeout branch lives in the baseline inline wrapper, not the extracted target wrapper.")
  }
  server_fn <- makeSharedProteinRollupFunctionWithOverrides(
    server_fn,
    list(
      file.exists = function(...) FALSE,
      Sys.sleep = function(...) {
        captured$sleep_count <- (captured$sleep_count %||% 0L) + 1L
        invisible(NULL)
      }
    )
  )

  withSharedProteinRollupPackageMocks(server_env, captured, create_output = FALSE)
  withSharedProteinRollupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_iq_rollup = TRUE,
      revert_iq_rollup = FALSE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinRollupWorkflow(
    makeSharedProteinRollupPeptideState(),
    captured
  )

  server_fn(
    id = "protein_rollup",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$sleep_count, 30L)
  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_true(any(grepl(
    "IQ output file not created within timeout period",
    captured$log_error,
    fixed = TRUE
  )))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("IQ output file not created within timeout period", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "iq_rollup_working")
})

test_that("proteomics protein rollup module preserves generic apply error behavior", {
  captured <- newSharedProteinRollupCapture()
  server_fn <- getProtProteinRollupServer()
  server_env <- environment(server_fn)

  withSharedProteinRollupPackageMocks(server_env, captured, create_output = FALSE)
  withSharedProteinRollupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_iq_rollup = TRUE,
      revert_iq_rollup = FALSE
    ),
    evaluate_plot = FALSE
  )

  malformed_state <- makeSharedProteinRollupPeptideState()
  malformed_state@peptide_data <- data.frame(
    Protein.Ids = "P1",
    Stripped.Sequence = "AAA",
    Run = "sample_a",
    stringsAsFactors = FALSE
  )
  workflow_data <- makeSharedProteinRollupWorkflow(malformed_state, captured)

  server_fn(
    id = "protein_rollup",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_update", envir = captured, inherits = FALSE))
  expect_true(any(grepl("Peptide.Imputed", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("Peptide.Imputed", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "iq_rollup_working")
})

test_that("proteomics protein rollup module preserves revert behavior", {
  captured <- newSharedProteinRollupCapture()
  server_fn <- getProtProteinRollupServer()
  server_env <- environment(server_fn)

  withSharedProteinRollupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_iq_rollup = FALSE,
      revert_iq_rollup = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinRollupWorkflow(
    makeSharedProteinRollupPeptideState(),
    captured,
    history = c("raw_data_s4", "qvalue_filtered", "protein_s4_created")
  )

  server_fn(
    id = "protein_rollup",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_identical(captured$reverted_state, "qvalue_filtered")
  expect_identical(captured$output$iq_rollup_results, "Reverted to qvalue_filtered state")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("proteomics protein rollup module preserves revert error behavior", {
  captured <- newSharedProteinRollupCapture()
  server_fn <- getProtProteinRollupServer()
  server_env <- environment(server_fn)

  withSharedProteinRollupUiMocks(
    server_env,
    captured,
    input_values = list(
      apply_iq_rollup = FALSE,
      revert_iq_rollup = TRUE
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedProteinRollupWorkflow(
    makeSharedProteinRollupPeptideState(),
    captured,
    history = "protein_s4_created"
  )

  server_fn(
    id = "protein_rollup",
    workflow_data = workflow_data,
    experiment_paths = list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir()),
    omic_type = "proteomics",
    experiment_label = "DIA Experiment"
  )

  expect_true(any(grepl("Error reverting:", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("Error reverting:", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
})
