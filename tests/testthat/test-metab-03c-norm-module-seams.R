library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadMetabNormModuleEnv <- function() {
  module_env <- new.env(parent = globalenv())
  sys.source(file.path(repo_root, "R", "mod_metab_norm_server_helpers.R"), envir = module_env)
  sys.source(file.path(repo_root, "R", "mod_metab_norm_server.R"), envir = module_env)
  module_env
}

test_that("runMetabNormModuleServerShell preserves wrapper orchestration through helper seams", {
  module_env <- loadMetabNormModuleEnv()
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  norm_state <- new.env(parent = emptyenv())
  norm_state$assay_names <- NULL

  input <- list(
    color_variable = "group",
    shape_variable = "batch",
    run_normalization = 1,
    reset_normalization = 1,
    apply_correlation_filter = 1,
    skip_correlation_filter = 1,
    export_session = 1
  )
  session <- list(ns = identity)
  workflow_data <- list(
    state_manager = "state-manager",
    design_matrix = "design-matrix"
  )
  experiment_paths <- list(metabolite_qc_dir = tempdir())
  selected_tab <- function() "normalization"

  observe_fn <- function(expr) {
    force(expr)
    invisible(NULL)
  }
  observe_event_fn <- function(eventExpr, handlerExpr, ignoreInit = FALSE) {
    force(eventExpr)
    force(handlerExpr)
    invisible(NULL)
  }

  result <- module_env$runMetabNormModuleServerShell(
    input = input,
    output = output,
    session = session,
    id = "norm",
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    omicType = "metabolomics",
    experimentLabel = "Metabolomics",
    selectedTab = selected_tab,
    logInfoFn = function(...) invisible(NULL),
    createReactiveStateFn = function() {
      captured$create_reactive_state <- TRUE
      norm_state
    },
    appendNormalizationLogFn = function(...) {
      captured$append_log <- list(...)
      invisible(NULL)
    },
    observeFn = observe_fn,
    observeEventFn = observe_event_fn,
    initializeAssayNamesFn = function(...) {
      captured$initialize_assay_names <- list(...)
      invisible(NULL)
    },
    runAutoPreNormalizationQcObserverShellFn = function(...) {
      captured$auto_pre_qc <- list(...)
      invisible(NULL)
    },
    updateDesignDrivenChoicesFn = function(...) {
      captured$update_design_choices <- list(...)
      invisible(NULL)
    },
    renderNormalizationLogFn = function(...) {
      captured$render_log <- list(...)
      "render-log-state"
    },
    renderItsdSelectionUiFn = function(...) {
      captured$render_itsd_ui <- list(...)
      "render-itsd-ui-state"
    },
    renderRuvQcUiFn = function(...) {
      captured$render_ruv_ui <- list(...)
      "render-ruv-ui-state"
    },
    runQcImageBindingShellFn = function(...) {
      captured$qc_image_binding <- list(...)
      invisible(NULL)
    },
    runAssayLabelBindingShellFn = function(...) {
      captured$assay_label_binding <- list(...)
      invisible(NULL)
    },
    runItsdSelectionTableObserverShellFn = function(...) {
      captured$itsd_table_observer <- list(...)
      invisible(NULL)
    },
    runItsdSelectionTrackingObserverShellFn = function(...) {
      captured$itsd_tracking_observer <- list(...)
      invisible(NULL)
    },
    runNormalizationObserverWrapperFn = function(...) {
      captured$run_normalization <- list(...)
      invisible(NULL)
    },
    runResetNormalizationObserverWrapperFn = function(...) {
      captured$reset_normalization <- list(...)
      invisible(NULL)
    },
    runRuvBindingObserverShellFn = function(...) {
      captured$ruv_binding <- list(...)
      invisible(NULL)
    },
    runApplyCorrelationObserverWrapperFn = function(...) {
      captured$apply_correlation <- list(...)
      invisible(NULL)
    },
    runSkipCorrelationObserverWrapperFn = function(...) {
      captured$skip_correlation <- list(...)
      invisible(NULL)
    },
    renderCorrelationFilterSummaryFn = function(...) {
      captured$render_correlation_summary <- list(...)
      "render-correlation-summary"
    },
    renderFinalQcPlotFn = function(...) {
      captured$render_final_qc <- list(...)
      "render-final-qc"
    },
    runExportSessionObserverWrapperFn = function(...) {
      captured$export_session <- list(...)
      invisible(NULL)
    }
  )

  expect_identical(result, norm_state)
  expect_true(isTRUE(captured$create_reactive_state))
  expect_identical(captured$initialize_assay_names$stateManager, "state-manager")
  expect_identical(captured$initialize_assay_names$normData, norm_state)
  expect_identical(captured$auto_pre_qc$selectedTab, "normalization")
  expect_identical(captured$auto_pre_qc$experimentPaths, experiment_paths)
  expect_identical(captured$update_design_choices$session, session)
  expect_identical(captured$update_design_choices$designMatrix, "design-matrix")
  expect_identical(captured$render_log$normData, norm_state)
  expect_identical(output$norm_log, "render-log-state")
  expect_identical(captured$render_itsd_ui$ns, identity)
  expect_identical(output$itsd_selection_ui, "render-itsd-ui-state")
  expect_identical(output$ruv_qc_ui, "render-ruv-ui-state")
  expect_identical(captured$qc_image_binding$qcDir, experiment_paths$metabolite_qc_dir)
  expect_true(is.function(captured$assay_label_binding$getAssayNamesFn))
  expect_identical(captured$itsd_table_observer$workflowData, workflow_data)
  expect_identical(captured$itsd_tracking_observer$input$shape_variable, "batch")
  expect_identical(captured$run_normalization$omicType, "metabolomics")
  expect_identical(captured$reset_normalization$workflowData, workflow_data)
  expect_identical(captured$ruv_binding$normData, norm_state)
  expect_identical(captured$apply_correlation$workflowData, workflow_data)
  expect_identical(captured$skip_correlation$workflowData, workflow_data)
  expect_identical(output$correlation_filter_summary, "render-correlation-summary")
  expect_identical(captured$render_final_qc$colorVariableFn(), "group")
  expect_identical(captured$render_final_qc$shapeVariableFn(), "batch")
  expect_identical(output$final_qc_plot, "render-final-qc")
  expect_identical(captured$export_session$experimentLabel, "Metabolomics")
  expect_true(is.function(captured$run_normalization$addLogFn))
})

test_that("runMetabNormModuleServerEntryShell preserves moduleServer handoff", {
  module_env <- loadMetabNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  result <- module_env$runMetabNormModuleServerEntryShell(
    id = "norm",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    omicType = "metabolomics",
    experimentLabel = "Metabolomics",
    selectedTab = function() "norm",
    moduleServerFn = function(id, module) {
      captured$id <- id
      captured$module <- module
      "module-server-state"
    },
    runModuleServerShellFn = function(...) {
      captured$shell_args <- list(...)
      "module-shell-state"
    }
  )

  expect_identical(result, "module-server-state")
  expect_identical(captured$id, "norm")
  expect_true(is.function(captured$module))

  captured$module(
    input = list(),
    output = new.env(parent = emptyenv()),
    session = list(ns = identity)
  )

  expect_identical(captured$shell_args$id, "norm")
  expect_identical(captured$shell_args$workflowData, "workflow-state")
  expect_identical(captured$shell_args$experimentPaths, "experiment-paths")
  expect_identical(captured$shell_args$omicType, "metabolomics")
  expect_identical(captured$shell_args$experimentLabel, "Metabolomics")
  expect_identical(captured$shell_args$selectedTab(), "norm")
})

test_that("runMetabNormModuleServerPublicWrapper preserves breadcrumb public-wrapper handoff", {
  module_env <- loadMetabNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  result <- module_env$runMetabNormModuleServerPublicWrapper(
    id = "norm",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    omic_type = "metabolomics",
    experiment_label = "Metabolomics",
    selected_tab = function() "norm",
    runModuleServerEntryShellFn = function(...) {
      captured$args <- list(...)
      "entry-shell-state"
    }
  )

  expect_identical(result, "entry-shell-state")
  expect_identical(captured$args$id, "norm")
  expect_identical(captured$args$workflowData, "workflow-state")
  expect_identical(captured$args$experimentPaths, "experiment-paths")
  expect_identical(captured$args$omicType, "metabolomics")
  expect_identical(captured$args$experimentLabel, "Metabolomics")
  expect_identical(captured$args$selectedTab(), "norm")
})

test_that("mod_metab_norm_server routes through the top-level public seam", {
  module_env <- loadMetabNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  module_env$runMetabNormModuleServerPublicWrapper <- function(...) {
    captured$args <- list(...)
    "public-wrapper-state"
  }

  result <- module_env$mod_metab_norm_server(
    id = "norm",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    omic_type = "metabolomics",
    experiment_label = "Metabolomics",
    selected_tab = function() "norm"
  )

  expect_identical(result, "public-wrapper-state")
  expect_identical(captured$args$id, "norm")
  expect_identical(captured$args$workflow_data, "workflow-state")
  expect_identical(captured$args$experiment_paths, "experiment-paths")
  expect_identical(captured$args$omic_type, "metabolomics")
  expect_identical(captured$args$experiment_label, "Metabolomics")
  expect_identical(captured$args$selected_tab(), "norm")
})
