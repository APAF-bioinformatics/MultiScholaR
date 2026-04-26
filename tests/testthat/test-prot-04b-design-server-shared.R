# fidelity-coverage-compare: shared
library(testthat)

getSharedProtDesignServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_prot_design_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_prot_design_server", envir = package_ns, inherits = FALSE))
  }

  mod_prot_design_server
}

test_that("proteomics design module server preserves public shell wiring", {
  server_fn <- getSharedProtDesignServer()
  server_env <- environment(server_fn)
  captured <- new.env(parent = emptyenv())
  captured$events <- character()

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$design_matrix <- data.frame(Run = "S1", group = "A", stringsAsFactors = FALSE)
  workflow_data$contrasts_tbl <- data.frame(contrasts = "groupA-groupB", stringsAsFactors = FALSE)
  workflow_data$tab_status <- list(design_matrix = "pending")

  experiment_paths <- list(
    base_dir = tempdir(),
    source_dir = tempdir(),
    results_dir = tempdir()
  )
  volumes <- c(Home = tempdir())
  qc_value <- FALSE
  qc_trigger <- function(value) {
    if (!missing(value)) {
      qc_value <<- value
    }
    qc_value
  }

  input <- list()
  output <- new.env(parent = emptyenv())
  session <- list(ns = function(id) paste("design", id, sep = "-"))

  has_target_shell <- exists(
    "registerProtDesignServerShells",
    envir = server_env,
    inherits = FALSE
  )
  if (has_target_shell) {
    local_mocked_bindings(
      registerProtDesignServerShells = function(input,
                                                output,
                                                session,
                                                workflowData,
                                                experimentPaths,
                                                volumes = NULL,
                                                qcTrigger = NULL,
                                                ...) {
        captured$shell_args <- list(
          input = input,
          output = output,
          session = session,
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          volumes = volumes,
          qcTrigger = qcTrigger
        )
        "builder-results-rv"
      },
      .env = server_env
    )
  } else if (exists("mod_prot_design_builder_server", envir = server_env, inherits = FALSE)) {
    local_mocked_bindings(
      mod_prot_design_builder_server = function(...) {
        captured$builder_registered <- TRUE
        function() NULL
      },
      .env = server_env
    )
  }

  local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(input, output, session)
      invisible(NULL)
    },
    isolate = function(expr) eval(substitute(expr), parent.frame()),
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
      captured$events <- c(
        captured$events,
        paste(deparse(substitute(eventExpr)), collapse = "")
      )
      invisible(NULL)
    },
    reactive = function(expr) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
    },
    renderText = function(expr) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
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
    outputOptions = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    showModal = function(...) invisible(NULL),
    removeModal = function(...) invisible(NULL),
    .package = "shiny"
  )
  local_mocked_bindings(
    shinyDirChoose = function(...) invisible(NULL),
    shinyFileChoose = function(...) invisible(NULL),
    .package = "shinyFiles"
  )
  local_mocked_bindings(
    renderDT = function(expr, ...) {
      quoted <- substitute(expr)
      env <- parent.frame()
      function() eval(quoted, env)
    },
    .package = "DT"
  )

  server_fn(
    id = "design",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = volumes,
    qc_trigger = qc_trigger
  )

  expect_identical(captured$module_id, "design")

  if (has_target_shell) {
    expect_identical(captured$shell_args$input, input)
    expect_identical(captured$shell_args$output, output)
    expect_identical(captured$shell_args$session, session)
    expect_identical(captured$shell_args$workflowData, workflow_data)
    expect_identical(captured$shell_args$experimentPaths, experiment_paths)
    expect_identical(captured$shell_args$volumes, volumes)
    expect_identical(captured$shell_args$qcTrigger, qc_trigger)
  } else {
    expect_true("input$confirm_import" %in% captured$events)
    expect_true(is.function(output$data_available))
    expect_true(is.function(output$design_matrix_exists))
    expect_true(is.function(output$design_matrix_preview))
    expect_true(is.function(output$contrasts_preview))
  }
})
