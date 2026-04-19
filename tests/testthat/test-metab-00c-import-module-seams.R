library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadMetabImportModuleEnv <- function() {
  module_env <- new.env(parent = globalenv())
  sys.source(file.path(repo_root, "R", "mod_metab_import_server_helpers.R"), envir = module_env)
  sys.source(file.path(repo_root, "R", "mod_metab_import_server.R"), envir = module_env)
  module_env
}

test_that("runMetabImportModuleServerShell preserves helper wiring and shared local state", {
  module_env <- loadMetabImportModuleEnv()
  local_data <- new.env(parent = emptyenv())
  local_data$all_headers <- NULL
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())

  column_accessors <- list(
    getMetaboliteIdCol = function() "metabolite_id",
    getAnnotationCol = function() "annotation",
    getSampleColumns = function() c("sample_1", "sample_2")
  )

  result <- module_env$runMetabImportModuleServerShell(
    input = list(vendor_format = "custom"),
    output = output,
    session = "import-session",
    id = "import",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    volumes = c(home = "/tmp/home"),
    requireNamespaceFn = function(pkg, quietly = TRUE) TRUE,
    createReactiveValuesFn = function(...) {
      captured$reactive_values <- list(...)
      local_data
    },
    setupAssaySelectionCallbackFn = function(...) {
      captured$selection_callback <- list(...)
      function() "import-data-callback"
    },
    setupShinyFilesFn = function(...) {
      captured$shiny_files <- list(...)
      c(data = "/tmp/data")
    },
    setupColumnAccessorsFn = function(...) {
      captured$column_accessors <- list(...)
      column_accessors
    },
    setupFileLoadedOutputFn = function(...) {
      captured$file_loaded <- list(...)
      invisible(NULL)
    },
    setupFormatDetectionStatusOutputFn = function(...) {
      captured$format_detection <- list(...)
      invisible(NULL)
    },
    setupMetaboliteIdStatusOutputFn = function(...) {
      captured$metabolite_id_status <- list(...)
      invisible(NULL)
    },
    setupAnnotationStatusOutputFn = function(...) {
      captured$annotation_status <- list(...)
      invisible(NULL)
    },
    setupSampleColumnsDisplayOutputFn = function(...) {
      captured$sample_columns_display <- list(...)
      invisible(NULL)
    },
    setupAvailableColumnsDisplayOutputFn = function(...) {
      captured$available_columns_display <- list(...)
      invisible(NULL)
    },
    setupCustomMetaboliteIdStatusOutputFn = function(...) {
      captured$custom_metabolite_id_status <- list(...)
      invisible(NULL)
    },
    setupCustomAnnotationStatusOutputFn = function(...) {
      captured$custom_annotation_status <- list(...)
      invisible(NULL)
    },
    setupValidationSummaryOutputFn = function(...) {
      captured$validation_summary <- list(...)
      invisible(NULL)
    },
    setupProcessingObserverFn = function(...) {
      captured$processing_observer <- list(...)
      invisible(NULL)
    },
    setupStatusOutputFn = function(...) {
      captured$status_output <- list(...)
      invisible(NULL)
    },
    logMessageFn = function(...) invisible(NULL)
  )

  expect_identical(result, local_data)
  expect_equal(captured$reactive_values$detected_format, NULL)
  expect_identical(captured$selection_callback$localData, local_data)
  expect_identical(captured$selection_callback$session, "import-session")
  expect_identical(captured$shiny_files$volumes, c(home = "/tmp/home"))
  expect_identical(captured$shiny_files$localData, local_data)
  expect_true(is.function(captured$shiny_files$importDataFn))
  expect_identical(captured$column_accessors$input$vendor_format, "custom")
  expect_identical(captured$column_accessors$localData, local_data)
  expect_identical(captured$file_loaded$localData, local_data)
  expect_identical(captured$format_detection$localData, local_data)
  expect_identical(captured$metabolite_id_status$input$vendor_format, "custom")
  expect_identical(captured$annotation_status$input$vendor_format, "custom")
  expect_identical(captured$validation_summary$columnAccessors, column_accessors)
  expect_identical(captured$processing_observer$workflowData, "workflow-state")
  expect_identical(captured$processing_observer$columnAccessors, column_accessors)
  expect_identical(captured$status_output$workflowData, "workflow-state")
})

test_that("runMetabImportModuleServerEntryShell preserves moduleServer handoff", {
  module_env <- loadMetabImportModuleEnv()
  captured <- new.env(parent = emptyenv())
  logged_messages <- character()

  result <- module_env$runMetabImportModuleServerEntryShell(
    id = "import",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    volumes = "volumes-state",
    moduleServerFn = function(id, module) {
      captured$id <- id
      captured$module <- module
      "module-server-state"
    },
    runModuleServerShellFn = function(...) {
      captured$shell_args <- list(...)
      "module-shell-state"
    },
    logMessageFn = function(...) {
      logged_messages <<- c(logged_messages, paste(..., collapse = " "))
      invisible(NULL)
    }
  )

  expect_identical(result, "module-server-state")
  expect_identical(captured$id, "import")
  expect_true(is.function(captured$module))
  expect_length(logged_messages, 2L)

  captured$module(
    input = list(),
    output = new.env(parent = emptyenv()),
    session = "import-session"
  )

  expect_identical(captured$shell_args$id, "import")
  expect_identical(captured$shell_args$workflowData, "workflow-state")
  expect_identical(captured$shell_args$experimentPaths, "experiment-paths")
  expect_identical(captured$shell_args$volumes, "volumes-state")
  expect_identical(captured$shell_args$session, "import-session")
})

test_that("runMetabImportModuleServerPublicWrapper preserves breadcrumb public-wrapper handoff", {
  module_env <- loadMetabImportModuleEnv()
  captured <- new.env(parent = emptyenv())

  result <- module_env$runMetabImportModuleServerPublicWrapper(
    id = "import",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    volumes = "volumes-state",
    runModuleServerEntryShellFn = function(...) {
      captured$args <- list(...)
      "entry-shell-state"
    }
  )

  expect_identical(result, "entry-shell-state")
  expect_identical(captured$args$id, "import")
  expect_identical(captured$args$workflowData, "workflow-state")
  expect_identical(captured$args$experimentPaths, "experiment-paths")
  expect_identical(captured$args$volumes, "volumes-state")
})

test_that("mod_metab_import_server routes through the top-level public seam", {
  module_env <- loadMetabImportModuleEnv()
  captured <- new.env(parent = emptyenv())

  module_env$runMetabImportModuleServerPublicWrapper <- function(...) {
    captured$args <- list(...)
    "public-wrapper-state"
  }

  result <- module_env$mod_metab_import_server(
    id = "import",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    volumes = "volumes-state"
  )

  expect_identical(result, "public-wrapper-state")
  expect_identical(captured$args$id, "import")
  expect_identical(captured$args$workflow_data, "workflow-state")
  expect_identical(captured$args$experiment_paths, "experiment-paths")
  expect_identical(captured$args$volumes, "volumes-state")
})
