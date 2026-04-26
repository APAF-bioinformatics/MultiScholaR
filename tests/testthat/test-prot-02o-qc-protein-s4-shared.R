# fidelity-coverage-compare: shared
library(testthat)

if (!methods::isClass("FakeSharedProteinS4State")) {
  methods::setClass(
    "FakeSharedProteinS4State",
    slots = c(protein_quant_table = "data.frame")
  )
}

getProtProteinS4Server <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_qc_protein_s4_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_qc_protein_s4_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_qc_protein_s4_server
}

makeSharedProteinS4Object <- function(proteins = c("P1", "P2", "P2")) {
  methods::new(
    "FakeSharedProteinS4State",
    protein_quant_table = data.frame(
      Protein.Ids = proteins,
      S1 = seq_along(proteins) * 10,
      stringsAsFactors = FALSE
    )
  )
}

makeSharedProteinS4Workflow <- function(captured,
                                        data_cln = data.frame(
                                          Protein.Ids = c("P1", "P2"),
                                          S1 = c(10, 20),
                                          stringsAsFactors = FALSE
                                        ),
                                        revert_error = NULL) {
  state_manager <- list(
    saveState = function(...) {
      captured$save_state <- list(...)
      invisible(NULL)
    },
    revertToState = function(state_name) {
      captured$reverted_state <- state_name
      if (!is.null(revert_error)) {
        stop(revert_error)
      }
      "initial-state"
    }
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_cln <- data_cln
  workflow_data$design_matrix <- data.frame(
    Run = "S1",
    group = "Control",
    replicates = "R1",
    stringsAsFactors = FALSE
  )
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "TMT"))
  workflow_data$state_manager <- state_manager
  workflow_data
}

newSharedProteinS4Capture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_error <- character()
  captured
}

localSharedProteinS4Binding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

withSharedProteinS4PackageMocks <- function(server_env,
                                            captured,
                                            protein_s4,
                                            constructor_error = NULL) {
  localSharedProteinS4Binding(
    server_env,
    "ProteinQuantitativeData",
    function(...) {
      captured$constructor_args <- list(...)
      if (!is.null(constructor_error)) {
        stop(constructor_error)
      }
      protein_s4
    },
    .local_envir = parent.frame()
  )
}

withSharedProteinS4UiMocks <- function(captured, input_values) {
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value)) {
          stop("required value missing", call. = FALSE)
        }
      }
      values[[1]]
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
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
    log_info = function(message, ...) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    log_error = function(message, ...) {
      captured$log_error <- c(captured$log_error, message)
      invisible(NULL)
    },
    .package = "logger",
    .env = mock_frame
  )
}

test_that("proteomics protein S4 module preserves successful creation behavior", {
  captured <- newSharedProteinS4Capture()
  server_fn <- getProtProteinS4Server()
  server_env <- environment(server_fn)
  protein_s4 <- makeSharedProteinS4Object()

  withSharedProteinS4PackageMocks(server_env, captured, protein_s4)
  withSharedProteinS4UiMocks(
    captured,
    input_values = list(create_protein_s4 = TRUE, revert_s4_creation = FALSE)
  )

  workflow_data <- makeSharedProteinS4Workflow(captured)
  server_fn(
    id = "protein_s4",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "TMT Experiment"
  )

  expect_identical(captured$constructor_args$protein_quant_table, workflow_data$data_cln)
  expect_identical(captured$constructor_args$protein_id_column, "Protein.Ids")
  expect_identical(captured$constructor_args$sample_id, "Run")
  expect_identical(captured$constructor_args$group_id, "group")
  expect_identical(captured$constructor_args$technical_replicate_id, "replicates")
  expect_identical(captured$constructor_args$args, workflow_data$config_list)
  expect_identical(captured$save_state$state_name, "protein_s4_initial")
  expect_identical(captured$save_state$s4_data_object, protein_s4)
  expect_identical(
    captured$save_state$config_object,
    list(s4_class = "ProteinQuantitativeData", protein_id_column = "Protein.Ids")
  )
  expect_identical(
    captured$save_state$description,
    "Created initial ProteinQuantitativeData S4 object from protein-level data."
  )
  expect_match(captured$output$s4_creation_results, "Proteins loaded: 2", fixed = TRUE)
  expect_match(captured$output$s4_creation_results, "S4 Class: FakeSharedProteinS4State", fixed = TRUE)
  expect_match(captured$output$s4_creation_results, "State saved as: 'protein_s4_initial'", fixed = TRUE)
  expect_identical(captured$removed_notifications, "s4_creation_working")
  expect_identical(captured$notifications[[1L]]$message, "Creating Protein S4 object from imported data...")
  expect_identical(captured$notifications[[2L]]$message, "Protein S4 object created successfully")
  expect_identical(captured$notifications[[2L]]$type, "message")
})

test_that("proteomics protein S4 module preserves creation error behavior", {
  captured <- newSharedProteinS4Capture()
  server_fn <- getProtProteinS4Server()
  server_env <- environment(server_fn)

  withSharedProteinS4PackageMocks(
    server_env,
    captured,
    protein_s4 = makeSharedProteinS4Object(),
    constructor_error = "mock creation failure"
  )
  withSharedProteinS4UiMocks(
    captured,
    input_values = list(create_protein_s4 = TRUE, revert_s4_creation = FALSE)
  )

  workflow_data <- makeSharedProteinS4Workflow(captured)
  server_fn(
    id = "protein_s4",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "TMT Experiment"
  )

  expect_match(captured$log_error, "Error creating Protein S4 object: mock creation failure", fixed = TRUE)
  expect_identical(captured$notifications[[2L]]$type, "error")
  expect_identical(captured$notifications[[2L]]$duration, 15)
  expect_match(captured$notifications[[2L]]$message, "mock creation failure", fixed = TRUE)
  expect_identical(captured$removed_notifications, "s4_creation_working")
})

test_that("proteomics protein S4 module preserves successful revert behavior", {
  captured <- newSharedProteinS4Capture()
  server_fn <- getProtProteinS4Server()
  server_env <- environment(server_fn)

  withSharedProteinS4PackageMocks(server_env, captured, makeSharedProteinS4Object())
  withSharedProteinS4UiMocks(
    captured,
    input_values = list(create_protein_s4 = FALSE, revert_s4_creation = TRUE)
  )

  workflow_data <- makeSharedProteinS4Workflow(captured)
  server_fn(
    id = "protein_s4",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "TMT Experiment"
  )

  expect_identical(captured$reverted_state, "initial")
  expect_identical(
    captured$output$s4_creation_results,
    "Reverted to initial empty state. You may need to re-run previous steps."
  )
  expect_identical(captured$notifications[[1L]]$message, "Reverted successfully")
  expect_identical(captured$notifications[[1L]]$type, "message")
})

test_that("proteomics protein S4 module preserves revert error behavior", {
  captured <- newSharedProteinS4Capture()
  server_fn <- getProtProteinS4Server()
  server_env <- environment(server_fn)

  withSharedProteinS4PackageMocks(server_env, captured, makeSharedProteinS4Object())
  withSharedProteinS4UiMocks(
    captured,
    input_values = list(create_protein_s4 = FALSE, revert_s4_creation = TRUE)
  )

  workflow_data <- makeSharedProteinS4Workflow(captured, revert_error = "mock revert failure")
  server_fn(
    id = "protein_s4",
    workflow_data = workflow_data,
    omic_type = "proteomics",
    experiment_label = "TMT Experiment"
  )

  expect_identical(captured$reverted_state, "initial")
  expect_match(captured$log_error, "Error reverting: mock revert failure", fixed = TRUE)
  expect_identical(captured$notifications[[1L]]$type, "error")
  expect_match(captured$notifications[[1L]]$message, "mock revert failure", fixed = TRUE)
})
