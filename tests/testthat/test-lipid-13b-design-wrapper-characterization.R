# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

test_that("lipid design server preserves module registration handoff", {
  namespace <- asNamespace("MultiScholaR")
  skip_if_not(
    exists("mod_lipid_design_server", envir = namespace, inherits = FALSE),
    "mod_lipid_design_server is unavailable in this ref."
  )

  server <- get("mod_lipid_design_server", envir = namespace, inherits = FALSE)
  server_env <- environment(server)
  helper_names <- c(
    "initializeLipidDesignImportBootstrap",
    "registerLipidDesignImportModalShell",
    "registerLipidDesignImportConfirmationObserver",
    "registerLipidDesignBuilderModule",
    "registerLipidDesignBuilderResultsObserver",
    "registerLipidDesignPreviewOutputs"
  )
  has_body_helpers <- all(vapply(
    helper_names,
    exists,
    logical(1),
    envir = server_env,
    inherits = FALSE
  ))

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- list(LCMS_Pos = data.frame(Sample_1 = 1))
  workflow_data$config_list <- list(config = TRUE)
  workflow_data$design_matrix <- data.frame(Run = "Sample_1", group = "A")
  experiment_paths <- list(base_dir = tempfile("lipid-design-base-"))
  qc_trigger <- shiny::reactiveVal(0)
  captured <- new.env(parent = emptyenv())
  captured$helper_calls <- character(0)

  if (has_body_helpers) {
    testthat::local_mocked_bindings(
      initializeLipidDesignImportBootstrap = function(input, session, experimentPaths, volumes, ...) {
        captured$helper_calls <- c(captured$helper_calls, "bootstrap")
        captured$bootstrap <- list(
          input = input,
          session = session,
          experimentPaths = experimentPaths,
          volumes = volumes
        )
        list(resolvedVolumes = c(Project = experimentPaths$base_dir))
      },
      registerLipidDesignImportModalShell = function(input, output, session, resolvedVolumes, ...) {
        captured$helper_calls <- c(captured$helper_calls, "modal")
        captured$modal <- list(
          input = input,
          output = output,
          session = session,
          resolvedVolumes = resolvedVolumes
        )
        output
      },
      registerLipidDesignImportConfirmationObserver = function(input, resolvedVolumes, workflowData, experimentPaths, qcTrigger, ...) {
        captured$helper_calls <- c(captured$helper_calls, "import")
        captured$import <- list(
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          qcTrigger = qcTrigger,
          resolvedVolumes = resolvedVolumes
        )
        invisible(NULL)
      },
      registerLipidDesignBuilderModule = function(workflowData, ...) {
        captured$helper_calls <- c(captured$helper_calls, "builder")
        captured$builder <- list(workflowData = workflowData)
        "builder-results-rv"
      },
      registerLipidDesignBuilderResultsObserver = function(builderResultsRv, workflowData, experimentPaths, qcTrigger, ...) {
        captured$helper_calls <- c(captured$helper_calls, "builder-observer")
        captured$builder_observer <- list(
          builderResultsRv = builderResultsRv,
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          qcTrigger = qcTrigger
        )
        invisible(NULL)
      },
      registerLipidDesignPreviewOutputs = function(output, workflowData, ...) {
        captured$helper_calls <- c(captured$helper_calls, "preview")
        captured$preview <- list(output = output, workflowData = workflowData)
        output
      },
      .env = server_env
    )
  }

  testthat::local_mocked_bindings(
    moduleServer = function(id, module, ...) {
      captured$id <- id
      captured$module <- module

      if (isTRUE(has_body_helpers)) {
        output <- new.env(parent = emptyenv())
        module(
          input = list(show_import_modal = 0, confirm_import = 0),
          output = output,
          session = list(ns = function(value) paste(id, value, sep = "-"))
        )
        captured$output <- output
      }

      invisible("registered")
    },
    reactive = function(expr, ...) {
      force(expr)
      function() expr
    },
    outputOptions = function(output, name, suspendWhenHidden = TRUE, ...) {
      captured$output_options <- c(captured$output_options, name)
      invisible(output)
    },
    .package = "shiny"
  )

  result <- server(
    id = "design",
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    volumes = c(Project = experiment_paths$base_dir),
    qc_trigger = qc_trigger
  )

  expect_identical(result, "registered")
  expect_identical(captured$id, "design")
  expect_true(is.function(captured$module))

  if (has_body_helpers) {
    expect_identical(
      captured$helper_calls,
      c("bootstrap", "modal", "import", "builder", "builder-observer", "preview")
    )
    expect_identical(captured$bootstrap$experimentPaths, experiment_paths)
    expect_identical(captured$import$workflowData, workflow_data)
    expect_identical(captured$builder$workflowData, workflow_data)
    expect_identical(captured$builder_observer$builderResultsRv, "builder-results-rv")
    expect_identical(captured$builder_observer$qcTrigger, qc_trigger)
    expect_identical(captured$preview$workflowData, workflow_data)
    expect_identical(captured$output$data_available(), TRUE)
    expect_identical(captured$output$design_matrix_exists(), TRUE)
    expect_true(all(c("data_available", "design_matrix_exists") %in% captured$output_options))
  }
})
