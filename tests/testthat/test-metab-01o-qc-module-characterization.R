library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "mod_metab_qc.R")
  ),
  symbols = c(
    "initializeMetabQcSubmodules",
    "mod_metab_qc_server"
  ),
  env = environment()
)

test_that("metabolomics QC initialization seam preserves sub-module registration order", {
  helper_env <- environment(initializeMetabQcSubmodules)
  binding_names <- c(
    "mod_metab_qc_intensity_server",
    "mod_metab_qc_duplicates_server",
    "mod_metab_qc_itsd_server",
    "mod_metab_qc_s4_server"
  )
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = helper_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = helper_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = helper_env)
      } else if (exists(name, envir = helper_env, inherits = FALSE)) {
        rm(list = name, envir = helper_env)
      }
    }
  }, add = TRUE)

  capture <- new.env(parent = emptyenv())
  capture$calls <- list()

  assign(
    "mod_metab_qc_intensity_server",
    function(id, workflowData, omicType, experimentLabel) {
      capture$calls[[length(capture$calls) + 1L]] <- list(
        id = id,
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
    },
    envir = helper_env
  )
  assign(
    "mod_metab_qc_duplicates_server",
    function(id, workflowData, omicType, experimentLabel) {
      capture$calls[[length(capture$calls) + 1L]] <- list(
        id = id,
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
    },
    envir = helper_env
  )
  assign(
    "mod_metab_qc_itsd_server",
    function(id, workflowData, omicType, experimentLabel) {
      capture$calls[[length(capture$calls) + 1L]] <- list(
        id = id,
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
    },
    envir = helper_env
  )
  assign(
    "mod_metab_qc_s4_server",
    function(id, workflowData, omicType, experimentLabel) {
      capture$calls[[length(capture$calls) + 1L]] <- list(
        id = id,
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
    },
    envir = helper_env
  )

  visible <- withVisible(
    initializeMetabQcSubmodules(
      workflowData = "workflow-state",
      omicType = "metabolomics",
      experimentLabel = "Experiment A"
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(
    vapply(capture$calls, `[[`, character(1), "id"),
    c("intensity", "duplicates", "itsd", "s4_finalize")
  )
  expect_true(all(vapply(capture$calls, `[[`, character(1), "workflow_data") == "workflow-state"))
  expect_true(all(vapply(capture$calls, `[[`, character(1), "omic_type") == "metabolomics"))
  expect_true(all(vapply(capture$calls, `[[`, character(1), "experiment_label") == "Experiment A"))
})

test_that("metabolomics QC server delegates qc-trigger initialization through the seam", {
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_server)
  had_initializer <- exists("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)

  if (had_initializer) {
    original_initializer <- get("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_initializer) {
      assign("initializeMetabQcSubmodules", original_initializer, envir = server_env)
    } else if (exists("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)) {
      rm(list = "initializeMetabQcSubmodules", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "initializeMetabQcSubmodules",
    function(workflowData, omicType, experimentLabel) {
      captured$initializer_call <- list(
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
      invisible(NULL)
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(list(), new.env(parent = emptyenv()), list(ns = function(value) value))
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
    renderUI = function(expr) substitute(expr),
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    observe = function(expr, ...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_server(
    id = "qc",
    workflow_data = list(state_manager = NULL),
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Experiment A",
    qc_trigger = function() TRUE
  )

  expect_identical(captured$module_id, "qc")
  expect_identical(
    captured$initializer_call,
    list(
      workflow_data = list(state_manager = NULL),
      omic_type = "metabolomics",
      experiment_label = "Experiment A"
    )
  )
})

test_that("metabolomics QC server delegates auto-initialization through the seam when state appears", {
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_server)
  had_initializer <- exists("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)

  if (had_initializer) {
    original_initializer <- get("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_initializer) {
      assign("initializeMetabQcSubmodules", original_initializer, envir = server_env)
    } else if (exists("initializeMetabQcSubmodules", envir = server_env, inherits = FALSE)) {
      rm(list = "initializeMetabQcSubmodules", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "initializeMetabQcSubmodules",
    function(workflowData, omicType, experimentLabel) {
      captured$initializer_call <- list(
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
      invisible(NULL)
    },
    envir = server_env
  )

  state_manager <- list(
    getState = function() structure(list(status = "ready"), class = "MetaboliteAssayData")
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(list(), new.env(parent = emptyenv()), list(ns = function(value) value))
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
    renderUI = function(expr) substitute(expr),
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    observe = function(expr, ...) {
      eval(substitute(expr), parent.frame())
      invisible(NULL)
    },
    req = function(...) invisible(list(...)[[1]]),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(state_manager = state_manager)

  mod_metab_qc_server(
    id = "qc",
    workflow_data = workflow_data,
    experiment_paths = list(),
    omic_type = "metabolomics",
    experiment_label = "Experiment A",
    qc_trigger = function() NULL
  )

  expect_identical(captured$module_id, "qc")
  expect_identical(
    captured$initializer_call,
    list(
      workflow_data = workflow_data,
      omic_type = "metabolomics",
      experiment_label = "Experiment A"
    )
  )
})
