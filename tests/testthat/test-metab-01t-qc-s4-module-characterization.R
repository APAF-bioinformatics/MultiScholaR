# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

source("helpers-scoped-mocked-bindings.R")

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

register_binding_teardown(
  asNamespace("shiny"),
  c("showNotification"),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("logger"),
  c("log_info", "log_error"),
  .local_envir = environment()
)
register_binding_teardown(
  multiScholaRNamespace(),
  c("updateMetaboliteFiltering"),
  .local_envir = environment()
)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id",
      args = list()
    )
  )
}

skip_if_missing_metab_qc_s4_helpers <- function() {
  skip_if_not(
    hasMultiScholaRBinding("runMetabQcS4ServerBody"),
    "mod_metab_qc_s4 helper seams are only available on the refactored branch"
  )
}

makeMetabQcS4CharacterizationData <- function(label = "qc_s4_fixture") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c(paste0(label, "_m1"), paste0(label, "_m2")),
        annotation_id = c("ann_a", "ann_b"),
        S1 = c(10, NA_real_),
        S2 = c(20, 30),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id",
    args = list(label = label)
  )
}

makeMetabQcS4WorkflowHarness <- function(
  current_s4 = makeMetabQcS4CharacterizationData(),
  history = c("import_complete", "duplicates_resolved", "itsd_applied")
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_s4
  state_manager$current_history <- history
  state_manager$getState <- function() state_manager$current_state
  state_manager$getHistory <- function() state_manager$current_history
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1L]] <<- list(...)
    invisible(NULL)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(qc = "complete")
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "locked"
  )

  list(
    workflow_data = workflow_data,
    state_manager = state_manager,
    capture = capture,
    current_s4 = current_s4
  )
}

test_that("mod_metab_qc_s4_server preserves public finalize behavior", {
  harness <- makeMetabQcS4WorkflowHarness()
  mod_metab_qc_s4_server <- getMultiScholaRBinding("mod_metab_qc_s4_server")

  scoped_mocked_bindings(
    showNotification = function(...) invisible(NULL),
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .env = asNamespace("logger")
  )
  scoped_mocked_bindings(
    updateMetaboliteFiltering = function(...) grid::nullGrob(),
    .env = multiScholaRNamespace()
  )

  testServer(
    mod_metab_qc_s4_server,
    args = list(
      workflow_data = harness$workflow_data,
      omic_type = "metabolomics",
      experiment_label = "Metabolomics"
    ),
    {
      session$setInputs(finalize_qc = 1)
    }
  )

  expect_identical(harness$workflow_data$tab_status$quality_control, "complete")
  expect_length(harness$capture$saved_states, 1L)
  expect_identical(harness$capture$saved_states[[1]]$state_name, "metab_qc_complete")
})

test_that("metabolomics QC S4 helper builders preserve summary and assay statistics behavior", {
  skip_if_missing_metab_qc_s4_helpers()

  current_s4 <- makeMetabQcS4CharacterizationData()
  build_data_summary_ui <- getMultiScholaRBinding("buildMetabQcS4DataSummaryUi")
  build_assay_stats_datatable <- getMultiScholaRBinding("buildMetabQcS4AssayStatsDatatable")

  summary_html <- htmltools::renderTags(build_data_summary_ui(current_s4))$html
  expect_match(summary_html, "Number of Assays", fixed = TRUE)
  expect_match(summary_html, "Total Metabolites", fixed = TRUE)
  expect_match(summary_html, "metabolite_id", fixed = TRUE)

  stats_df <- build_assay_stats_datatable(
    current_s4,
    datatableFn = function(data, ...) data
  )

  expect_identical(stats_df$Assay, "Plasma")
  expect_identical(stats_df$Metabolites, 2)
  expect_identical(stats_df$Samples, 2)
  expect_identical(stats_df$Missingness, 25)
})

test_that("metabolomics QC S4 history and state helpers preserve fallback behavior", {
  skip_if_missing_metab_qc_s4_helpers()

  build_state_history_ui <- getMultiScholaRBinding("buildMetabQcS4StateHistoryUi")
  get_state_history <- getMultiScholaRBinding("getMetabQcS4StateHistory")
  build_state_history_render_output <- getMultiScholaRBinding("buildMetabQcS4StateHistoryRenderOutput")
  get_finalize_state <- getMultiScholaRBinding("getMetabQcS4FinalizeState")
  validate_finalize_state <- getMultiScholaRBinding("validateMetabQcS4FinalizeState")
  get_finalize_history <- getMultiScholaRBinding("getMetabQcS4FinalizeHistory")

  state_manager <- new.env(parent = emptyenv())
  state_manager$getHistory <- function() c("import_complete", "duplicates_resolved")
  state_manager$getState <- function() makeMetabQcS4CharacterizationData()

  history_html <- htmltools::renderTags(
    build_state_history_ui(state_manager$getHistory())
  )$html
  expect_match(history_html, "duplicates_resolved", fixed = TRUE)
  expect_match(history_html, "(current)", fixed = TRUE)

  fallback_history <- get_state_history(
    stateManager = state_manager,
    historyGetterFn = function(manager) stop("history unavailable", call. = FALSE),
    errorHistory = "fallback_state"
  )
  expect_identical(fallback_history, "fallback_state")

  rendered_history <- build_state_history_render_output(
    stateManager = state_manager,
    getStateHistoryFn = function(stateManager) c("a", "b"),
    buildStateHistoryUiFn = function(history) paste(history, collapse = " -> ")
  )
  expect_identical(rendered_history, "a -> b")

  expect_s4_class(
    get_finalize_state(stateManager = state_manager),
    "MetaboliteAssayData"
  )
  expect_s4_class(
    validate_finalize_state(
      currentS4 = state_manager$getState(),
      reqFn = function(value) value
    ),
    "MetaboliteAssayData"
  )
  expect_identical(
    get_finalize_history(stateManager = state_manager),
    c("import_complete", "duplicates_resolved")
  )
})

test_that("metabolomics QC S4 finalize helpers preserve persistence and completion behavior", {
  skip_if_missing_metab_qc_s4_helpers()

  save_completed_state <- getMultiScholaRBinding("saveMetabQcS4CompletedState")
  complete_tab_status <- getMultiScholaRBinding("completeMetabQcS4TabStatus")
  build_finalize_results_text <- getMultiScholaRBinding("buildMetabQcS4FinalizeResultsText")

  harness <- makeMetabQcS4WorkflowHarness()
  current_s4 <- harness$current_s4

  saved_state <- save_completed_state(
    stateManager = harness$state_manager,
    currentS4 = current_s4,
    configObject = harness$workflow_data$config_list
  )
  expect_identical(saved_state, "metab_qc_complete")
  expect_length(harness$capture$saved_states, 1L)

  completed_status <- complete_tab_status(harness$workflow_data)
  expect_identical(completed_status$quality_control, "complete")
  expect_identical(harness$workflow_data$tab_status$quality_control, "complete")

  result_text <- build_finalize_results_text(
    currentS4 = current_s4,
    history = c("import_complete", "duplicates_resolved")
  )
  expect_match(result_text, "QC Finalization Complete", fixed = TRUE)
  expect_match(result_text, "State saved as: 'metab_qc_complete'", fixed = TRUE)
})

test_that("metabolomics QC S4 render helpers preserve render-path behavior", {
  skip_if_missing_metab_qc_s4_helpers()

  get_data_summary_state <- getMultiScholaRBinding("getMetabQcS4DataSummaryState")
  get_assay_stats_state <- getMultiScholaRBinding("getMetabQcS4AssayStatsState")
  build_data_summary_render_output <- getMultiScholaRBinding("buildMetabQcS4DataSummaryRenderOutput")
  build_assay_stats_render_output <- getMultiScholaRBinding("buildMetabQcS4AssayStatsRenderOutput")
  build_filter_plot_render_output <- getMultiScholaRBinding("buildMetabQcS4FilterPlotRenderOutput")
  render_filter_plot <- getMultiScholaRBinding("renderMetabQcS4FilterPlot")

  state_manager <- new.env(parent = emptyenv())
  state_manager$getState <- function() makeMetabQcS4CharacterizationData()

  expect_s4_class(
    get_data_summary_state(state_manager),
    "MetaboliteAssayData"
  )
  expect_null(
    get_data_summary_state(
      state_manager,
      stateGetterFn = function(manager) stop("summary state missing", call. = FALSE)
    )
  )

  expect_s4_class(
    get_assay_stats_state(state_manager),
    "MetaboliteAssayData"
  )
  expect_null(
    get_assay_stats_state(
      state_manager,
      stateGetterFn = function(manager) stop("stats state missing", call. = FALSE)
    )
  )

  rendered_summary <- build_data_summary_render_output(
    stateManager = state_manager,
    buildDataSummaryUiFn = function(currentS4) currentS4@sample_id
  )
  expect_identical(rendered_summary, "Run")

  rendered_stats <- build_assay_stats_render_output(
    stateManager = state_manager,
    buildAssayStatsDatatableFn = function(currentS4) names(currentS4@metabolite_data)
  )
  expect_identical(rendered_stats, "Plasma")

  capture <- new.env(parent = emptyenv())
  rendered_plot <- build_filter_plot_render_output(
    filterPlot = function() grid::nullGrob(),
    renderFilterPlotFn = function(filterPlot) {
      capture$plot <- filterPlot()
      "filter-rendered"
    }
  )
  expect_identical(rendered_plot, "filter-rendered")
  expect_true(inherits(capture$plot, "grob"))

  draw_capture <- new.env(parent = emptyenv())
  visible <- withVisible(
    render_filter_plot(
      filterPlot = function() grid::nullGrob(),
      reqFn = function(value) value,
      gridDrawFn = function(plotObject) {
        draw_capture$plot <- plotObject
        invisible(NULL)
      }
    )
  )
  expect_false(visible$visible)
  expect_null(visible$value)
  expect_true(inherits(draw_capture$plot, "grob"))
})

test_that("metabolomics QC S4 finalize workflow preserves success and error orchestration", {
  skip_if_missing_metab_qc_s4_helpers()

  run_finalize_workflow <- getMultiScholaRBinding("runMetabQcS4FinalizeWorkflow")

  harness <- makeMetabQcS4WorkflowHarness()
  capture <- new.env(parent = emptyenv())
  capture$steps <- character()

  filter_plot <- local({
    stored <- NULL
    function(value) {
      if (missing(value)) {
        stored
      } else {
        stored <<- value
        invisible(NULL)
      }
    }
  })

  run_finalize_workflow(
    workflowData = harness$workflow_data,
    omicType = "metabolomics",
    filterPlot = filter_plot,
    getFinalizeStateFn = function(stateManager) {
      capture$steps <- c(capture$steps, "get")
      harness$current_s4
    },
    validateFinalizeStateFn = function(currentS4, ...) {
      capture$steps <- c(capture$steps, "validate")
      currentS4
    },
    saveCompletedStateFn = function(...) {
      capture$steps <- c(capture$steps, "save")
      invisible(NULL)
    },
    updateTrackingPlotFn = function(...) {
      capture$steps <- c(capture$steps, "plot")
      invisible(NULL)
    },
    completeTabStatusFn = function(...) {
      capture$steps <- c(capture$steps, "tab")
      invisible(NULL)
    },
    getFinalizeHistoryFn = function(...) {
      capture$steps <- c(capture$steps, "history")
      c("import_complete", "duplicates_resolved")
    },
    reportFinalizeSuccessFn = function(...) {
      capture$steps <- c(capture$steps, "success")
      invisible(NULL)
    },
    reportFinalizeErrorFn = function(error) {
      capture$error <- error$message
      capture$steps <- c(capture$steps, "error")
      invisible(NULL)
    }
  )

  expect_identical(
    capture$steps,
    c("get", "validate", "save", "plot", "tab", "history", "success")
  )

  capture$steps <- character()
  capture$error <- NULL

  run_finalize_workflow(
    workflowData = harness$workflow_data,
    omicType = "metabolomics",
    filterPlot = filter_plot,
    getFinalizeStateFn = function(...) stop("finalize exploded", call. = FALSE),
    reportFinalizeErrorFn = function(error) {
      capture$error <- error$message
      capture$steps <- c(capture$steps, "error")
      invisible(NULL)
    }
  )

  expect_identical(capture$error, "finalize exploded")
  expect_identical(capture$steps, "error")
})

test_that("metabolomics QC S4 reporting and plot helpers preserve success and error tails", {
  skip_if_missing_metab_qc_s4_helpers()

  report_finalize_success <- getMultiScholaRBinding("reportMetabQcS4FinalizeSuccess")
  report_finalize_error <- getMultiScholaRBinding("reportMetabQcS4FinalizeError")
  update_tracking_plot <- getMultiScholaRBinding("updateMetabQcS4TrackingPlot")

  capture <- new.env(parent = emptyenv())
  current_s4 <- makeMetabQcS4CharacterizationData()
  output <- new.env(parent = emptyenv())

  result_text <- report_finalize_success(
    currentS4 = current_s4,
    history = c("import_complete", "duplicates_resolved"),
    output = output,
    buildResultsTextFn = function(...) "finalization complete",
    renderTextFn = function(value) value,
    logInfoFn = function(message) {
      capture$info <- message
      invisible(NULL)
    },
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      capture$notification <- list(message = message, type = type, duration = duration)
      invisible(NULL)
    }
  )

  expect_identical(result_text, "finalization complete")
  expect_identical(output$finalize_results, "finalization complete")
  expect_identical(capture$info, "Metabolomics QC finalized successfully")
  expect_identical(
    capture$notification,
    list(
      message = "QC complete! Proceed to Normalization.",
      type = "message",
      duration = 5
    )
  )

  error_text <- report_finalize_error(
    error = simpleError("save failed"),
    logErrorFn = function(message) {
      capture$error_log <- message
      invisible(NULL)
    },
    showNotificationFn = function(message, type = NULL) {
      capture$error_notification <- list(message = message, type = type)
      invisible(NULL)
    }
  )

  expect_identical(error_text, "Error finalizing QC: save failed")
  expect_identical(capture$error_log, "Error finalizing QC: save failed")
  expect_identical(
    capture$error_notification,
    list(message = "Error finalizing QC: save failed", type = "error")
  )

  updated_plot <- update_tracking_plot(
    currentS4 = current_s4,
    omicType = "metabolomics",
    setFilterPlotFn = function(plotObject) {
      capture$plot <- plotObject
      invisible(NULL)
    },
    updateMetaboliteFilteringFn = function(...) grid::nullGrob()
  )
  expect_true(inherits(updated_plot, "grob"))
  expect_true(inherits(capture$plot, "grob"))

  failed_plot <- update_tracking_plot(
    currentS4 = current_s4,
    omicType = "metabolomics",
    setFilterPlotFn = function(plotObject) {
      capture$failed_plot <- plotObject
      invisible(NULL)
    },
    updateMetaboliteFilteringFn = function(...) stop("plot failed", call. = FALSE)
  )
  expect_null(failed_plot)
  expect_null(capture$failed_plot)
})

test_that("metabolomics QC S4 server body preserves render and observer wiring", {
  skip_if_missing_metab_qc_s4_helpers()

  run_server_body <- getMultiScholaRBinding("runMetabQcS4ServerBody")

  harness <- makeMetabQcS4WorkflowHarness()
  capture <- new.env(parent = emptyenv())
  input <- list(finalize_qc = 1)
  output <- new.env(parent = emptyenv())

  reactive_val_fn <- function(value = NULL) {
    stored <- value

    function(new_value) {
      if (missing(new_value)) {
        stored
      } else {
        stored <<- new_value
        invisible(NULL)
      }
    }
  }

  result <- run_server_body(
    input = input,
    output = output,
    session = list(ns = function(id) id),
    workflowData = harness$workflow_data,
    omicType = "metabolomics",
    experimentLabel = "Metabolomics",
    reactiveValFn = reactive_val_fn,
    renderUiFn = function(expr) eval(substitute(expr), parent.frame()),
    renderDtFn = function(expr) eval(substitute(expr), parent.frame()),
    observeEventFn = function(eventExpr, handlerExpr, ...) {
      capture$event <- eval(substitute(eventExpr), parent.frame())
      eval(substitute(handlerExpr), parent.frame())
      invisible(NULL)
    },
    renderPlotFn = function(expr) eval(substitute(expr), parent.frame()),
    reqFn = function(value) value,
    buildStateHistoryRenderOutputFn = function(stateManager) "history-ui",
    buildDataSummaryRenderOutputFn = function(stateManager) "summary-ui",
    buildAssayStatsRenderOutputFn = function(stateManager) "assay-stats",
    runFinalizeWorkflowFn = function(workflowData, omicType, filterPlot, output) {
      capture$finalize <- list(
        workflowData = workflowData,
        omicType = omicType,
        filterPlot = filterPlot,
        output = output
      )
      filterPlot(grid::nullGrob())
      invisible(NULL)
    },
    buildFilterPlotRenderOutputFn = function(filterPlot) {
      capture$filter_plot_value <- filterPlot()
      "filter-plot"
    }
  )

  expect_null(result)
  expect_identical(output$state_history, "history-ui")
  expect_identical(output$data_summary, "summary-ui")
  expect_identical(output$assay_stats_table, "assay-stats")
  expect_identical(output$filter_plot, "filter-plot")
  expect_identical(capture$event, 1)
  expect_identical(capture$finalize$workflowData, harness$workflow_data)
  expect_identical(capture$finalize$omicType, "metabolomics")
  expect_true(is.function(capture$finalize$filterPlot))
  expect_true(inherits(capture$filter_plot_value, "grob"))
})
