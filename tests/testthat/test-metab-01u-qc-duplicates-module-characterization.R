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
  c("showNotification", "removeNotification"),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("logger"),
  c("log_info", "log_warn", "log_error"),
  .local_envir = environment()
)
register_binding_teardown(
  multiScholaRNamespace(),
  c(
    "findMetabDuplicateFeatureIDs",
    "resolveDuplicateFeaturesByIntensity",
    "updateMetaboliteFiltering"
  ),
  .local_envir = environment()
)

skip_if_missing_metab_duplicate_helpers <- function() {
  skip_if_not(
    hasMultiScholaRBinding("runMetabDuplicateResolutionWorkflow"),
    "mod_metab_qc_duplicates helper seams are only available on the refactored branch"
  )
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      database_identifier_type = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id",
      database_identifier_type = "Unknown",
      internal_standard_regex = NA_character_,
      design_matrix = data.frame(),
      sample_id = "Sample_ID",
      group_id = "group",
      technical_replicate_id = NA_character_,
      args = list()
    )
  )
}

makeMetabDuplicateCharacterizationData <- function(label = "metab_dup_fixture") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c(
          paste0(label, "_m1"),
          paste0(label, "_m1"),
          paste0(label, "_m2")
        ),
        annotation_id = c("ann_1", "ann_1", "ann_2"),
        Sample1 = c(10, 4, 7),
        Sample2 = c(9, 3, 6),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id",
    database_identifier_type = "Mock",
    internal_standard_regex = "",
    design_matrix = data.frame(
      Sample_ID = c("Sample1", "Sample2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list(label = label)
  )
}

makeMetabDuplicateHarness <- function(
  current_s4 = makeMetabDuplicateCharacterizationData(),
  history = c("import_complete", "qc_ready", "metab_duplicates_resolved")
) {
  capture <- new.env(parent = emptyenv())
  capture$notifications <- list()
  capture$removed_notifications <- character()
  capture$saved_states <- list()
  capture$reverted_to <- NULL

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_s4
  state_manager$current_history <- history
  state_manager$getState <- function() state_manager$current_state
  state_manager$getHistory <- function() state_manager$current_history
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1L]] <<- list(...)
    invisible(NULL)
  }
  state_manager$revertToState <- function(state_name) {
    capture$reverted_to <<- state_name
    invisible(state_name)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$config_list <- list(qc = "ready")

  list(
    capture = capture,
    workflow_data = workflow_data,
    state_manager = state_manager,
    current_s4 = current_s4
  )
}

notificationText <- function(notification) {
  paste(as.character(notification$ui), collapse = "")
}

recordNotification <- function(capture, ui, id = NULL, type = "default", duration = 5) {
  capture$notifications[[length(capture$notifications) + 1L]] <- list(
    ui = ui,
    id = id,
    type = type,
    duration = duration
  )
  invisible(if (is.null(id)) paste0("notification-", length(capture$notifications)) else id)
}

test_that("mod_metab_qc_duplicates_server preserves detect-resolve-revert public behavior", {
  harness <- makeMetabDuplicateHarness()
  mod_metab_qc_duplicates_server <- getMultiScholaRBinding("mod_metab_qc_duplicates_server")

  scoped_mocked_bindings(
    showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL, type = "default", session = getDefaultReactiveDomain()) {
      recordNotification(harness$capture, ui = ui, id = id, type = type, duration = duration)
    },
    removeNotification = function(id, session = getDefaultReactiveDomain()) {
      harness$capture$removed_notifications <<- c(harness$capture$removed_notifications, id)
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .env = asNamespace("logger")
  )
  scoped_mocked_bindings(
    findMetabDuplicateFeatureIDs = function(theObject) {
      expect_s4_class(theObject, "MetaboliteAssayData")
      list(
        Plasma = data.frame(
          metabolite_id = c("dup_1", "dup_2"),
          duplicate_count = c(2L, 3L),
          stringsAsFactors = FALSE
        )
      )
    },
    resolveDuplicateFeaturesByIntensity = function(assay_tibble, id_col, sample_cols) {
      expect_identical(id_col, "metabolite_id")
      expect_identical(sample_cols, c("Sample1", "Sample2"))
      assay_tibble[c(1, 3), , drop = FALSE]
    },
    updateMetaboliteFiltering = function(...) grid::nullGrob(),
    .env = multiScholaRNamespace()
  )

  testServer(
    mod_metab_qc_duplicates_server,
    args = list(
      workflow_data = harness$workflow_data,
      omic_type = "metabolomics",
      experiment_label = "Metabolomics"
    ),
    {
      session$setInputs(detect_duplicates = 1)
      session$setInputs(resolve_duplicates = 1)
      session$setInputs(revert_duplicates = 1)
    }
  )

  notice_text <- vapply(harness$capture$notifications, notificationText, character(1))
  expect_length(harness$capture$saved_states, 1L)
  expect_identical(harness$capture$saved_states[[1]]$state_name, "metab_duplicates_resolved")
  expect_identical(harness$capture$reverted_to, "qc_ready")
  expect_true(any(grepl("Detection complete: 2 duplicate IDs found", notice_text, fixed = TRUE)))
  expect_true(any(grepl("Duplicates resolved: 1 rows removed", notice_text, fixed = TRUE)))
  expect_true(any(grepl("Reverted successfully", notice_text, fixed = TRUE)))
  expect_true("metab_dup_resolve_working" %in% harness$capture$removed_notifications)
})

test_that("metabolomics duplicate helper builders preserve assay summaries and renderers", {
  skip_if_missing_metab_duplicate_helpers()

  build_summary_ui <- getMultiScholaRBinding("buildMetabDuplicateSummaryUi")
  build_tables_ui <- getMultiScholaRBinding("buildMetabDuplicateTablesUi")
  register_table_renderers <- getMultiScholaRBinding("registerMetabDuplicateTableRenderers")
  render_filter_plot <- getMultiScholaRBinding("renderMetabDuplicateFilterPlot")

  summary_html <- htmltools::renderTags(
    build_summary_ui(list(
      Plasma = data.frame(metabolite_id = c("dup_1", "dup_2")),
      Serum = NULL
    ))
  )$html
  expect_match(summary_html, "Plasma: 2 duplicate IDs", fixed = TRUE)
  expect_match(summary_html, "Serum: 0 duplicate IDs", fixed = TRUE)

  tables_html <- htmltools::renderTags(
    build_tables_ui(
      dupList = list(
        Plasma = data.frame(metabolite_id = c("dup_1", "dup_2")),
        Serum = NULL
      ),
      nsFn = function(id) paste0("metab-", id)
    )
  )$html
  expect_match(tables_html, "Plasma", fixed = TRUE)
  expect_match(tables_html, "metab-dup_table_Plasma", fixed = TRUE)
  expect_false(grepl("Serum", tables_html, fixed = TRUE))

  output_env <- new.env(parent = emptyenv())
  register_table_renderers(
    dupList = list(
      Plasma = data.frame(metabolite_id = c("dup_1", "dup_2")),
      Serum = NULL
    ),
    output = output_env,
    renderDtFn = function(expr) {
      force(expr)
      "datatable-rendered"
    },
    datatableFn = function(data, ...) data
  )
  expect_identical(output_env$dup_table_Plasma, "datatable-rendered")
  expect_false(exists("dup_table_Serum", envir = output_env, inherits = FALSE))

  drawn_grob <- NULL
  render_filter_plot(
    filterPlot = function() grid::nullGrob(),
    reqFn = function(value) value,
    inheritsFn = inherits,
    gridDrawFn = function(value) {
      drawn_grob <<- value
      invisible(NULL)
    },
    printFn = function(...) stop("printFn should not be called for grobs")
  )
  expect_true(inherits(drawn_grob, "grob"))

  printed_plot <- NULL
  render_filter_plot(
    filterPlot = function() structure(list(token = "ggplot"), class = "ggplot"),
    reqFn = function(value) value,
    inheritsFn = inherits,
    gridDrawFn = function(...) stop("gridDrawFn should not be called for ggplot"),
    printFn = function(value) {
      printed_plot <<- value
      invisible(NULL)
    }
  )
  expect_identical(printed_plot$token, "ggplot")
})

test_that("metabolomics duplicate misc helpers preserve resolution workflow behavior", {
  skip_if_missing_metab_duplicate_helpers()

  resolve_assay_data <- getMultiScholaRBinding("resolveMetabDuplicateAssayData")
  build_resolution_summary <- getMultiScholaRBinding("buildMetabDuplicateResolutionSummary")
  detect_features <- getMultiScholaRBinding("detectMetabDuplicateFeatures")
  prepare_resolution_state <- getMultiScholaRBinding("prepareMetabDuplicateResolutionState")
  apply_resolution_state <- getMultiScholaRBinding("applyMetabDuplicateResolutionState")
  run_resolution_workflow <- getMultiScholaRBinding("runMetabDuplicateResolutionWorkflow")
  revert_resolution <- getMultiScholaRBinding("revertMetabDuplicateResolution")

  assay_result <- resolve_assay_data(
    assayList = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m1", "m2"),
        Sample1 = c(1, 5, 3),
        Sample2 = c(2, 4, 7),
        stringsAsFactors = FALSE
      ),
      MetadataOnly = data.frame(
        metabolite_id = c("m3", "m4"),
        note = c("x", "y"),
        stringsAsFactors = FALSE
      )
    ),
    metaboliteIdCol = "metabolite_id",
    resolveDuplicateFeaturesByIntensityFn = function(assay_tibble, id_col, sample_cols) {
      expect_identical(id_col, "metabolite_id")
      expect_identical(sample_cols, c("Sample1", "Sample2"))
      assay_tibble[!duplicated(assay_tibble[[id_col]]), , drop = FALSE]
    },
    logWarnFn = function(...) invisible(NULL)
  )
  expect_identical(assay_result$statsList$Plasma$removed, 1L)
  expect_identical(assay_result$statsList$MetadataOnly$removed, 0)

  summary <- build_resolution_summary(
    statsList = assay_result$statsList,
    stateName = "custom_state"
  )
  expect_identical(summary$totalRemoved, 1)
  expect_match(summary$resultText, "State saved as: 'custom_state'", fixed = TRUE)

  current_s4 <- makeMetabDuplicateCharacterizationData()
  detection <- detect_features(
    stateManager = list(getState = function() current_s4),
    duplicateFinderFn = function(theObject) {
      expect_identical(theObject, current_s4)
      list(Plasma = data.frame(metabolite_id = c("dup_1", "dup_2")))
    },
    reqFn = function(value) value,
    inheritsFn = inherits,
    expectedClass = "MetaboliteAssayData"
  )
  expect_identical(detection$totalDuplicates, 2L)

  harness <- makeMetabDuplicateHarness()
  resolution_stats_value <- NULL
  filter_plot_value <- "old"

  preflight <- prepare_resolution_state(
    stateManager = harness$state_manager,
    resolveDuplicateAssayDataFn = function(assayList, metaboliteIdCol) {
      expect_identical(metaboliteIdCol, "metabolite_id")
      list(
        resolvedAssayList = list(
          Plasma = assayList$Plasma[c(1, 3), , drop = FALSE]
        ),
        statsList = list(Plasma = list(original = 3, resolved = 2, removed = 1))
      )
    },
    reqFn = function(value) value,
    inheritsFn = inherits,
    expectedClass = "MetaboliteAssayData"
  )
  expect_identical(preflight$statsList$Plasma$removed, 1)

  apply_resolution_state(
    currentS4 = preflight$currentS4,
    statsList = preflight$statsList,
    workflowData = harness$workflow_data,
    omicType = "metabolomics",
    setResolutionStatsFn = function(value) {
      resolution_stats_value <<- value
    },
    setFilterPlotFn = function(value) {
      filter_plot_value <<- value
    },
    updateMetaboliteFilteringFn = function(...) grid::nullGrob(),
    logWarnFn = function(...) invisible(NULL)
  )
  expect_identical(harness$capture$saved_states[[1]]$state_name, "metab_duplicates_resolved")
  expect_identical(resolution_stats_value$Plasma$removed, 1)
  expect_true(inherits(filter_plot_value, "grob"))

  workflow_result <- run_resolution_workflow(
    workflowData = harness$workflow_data,
    omicType = "metabolomics",
    setResolutionStatsFn = function(value) {
      resolution_stats_value <<- value
    },
    setFilterPlotFn = function(value) {
      filter_plot_value <<- value
    },
    prepareResolutionStateFn = function(...) preflight,
    applyResolutionStateFn = function(...) list(stateName = "workflow_state", qcPlot = NULL),
    buildResolutionSummaryFn = function(statsList, stateName) {
      list(totalRemoved = 1, resultText = paste("saved", stateName))
    }
  )
  expect_identical(workflow_result$totalRemoved, 1)
  expect_identical(workflow_result$resultText, "saved workflow_state")

  revert_result <- revert_resolution(
    stateManager = harness$state_manager,
    reqFn = function(value) value,
    historyGetterFn = function(manager) manager$getHistory(),
    revertStateFn = function(manager, stateName) manager$revertToState(stateName)
  )
  expect_identical(revert_result$previousStateName, "qc_ready")
  expect_identical(harness$capture$reverted_to, "qc_ready")
})

test_that("metabolomics duplicate observer helpers preserve notification and output behavior", {
  skip_if_missing_metab_duplicate_helpers()

  report_detection <- getMultiScholaRBinding("reportMetabDuplicateDetection")
  report_detection_error <- getMultiScholaRBinding("reportMetabDuplicateDetectionError")
  run_detection_observer <- getMultiScholaRBinding("runMetabDuplicateDetectionObserverShell")
  run_resolution_observer_shell <- getMultiScholaRBinding("runMetabDuplicateResolutionObserverShell")
  run_revert_observer_shell <- getMultiScholaRBinding("runMetabDuplicateRevertObserverShell")
  run_resolution_observer <- getMultiScholaRBinding("runMetabDuplicateResolutionObserver")

  notifications <- list()
  report_detection(
    totalDuplicates = 2L,
    logInfoFn = function(...) invisible(NULL),
    showNotificationFn = function(message, type) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    },
    sprintfFn = sprintf
  )
  expect_identical(notifications[[1]]$message, "Detection complete: 2 duplicate IDs found")
  expect_identical(notifications[[1]]$type, "warning")

  error_notifications <- list()
  report_detection_error(
    errorCondition = simpleError("boom"),
    logErrorFn = function(...) invisible(NULL),
    showNotificationFn = function(message, type) {
      error_notifications <<- c(error_notifications, list(list(message = message, type = type)))
    }
  )
  expect_identical(error_notifications[[1]]$message, "Error detecting duplicates: boom")
  expect_identical(error_notifications[[1]]$type, "error")

  duplicate_value <- NULL
  detection_dispatch <- run_detection_observer(
    stateManager = "state-manager",
    setDuplicateInfoFn = function(value) {
      duplicate_value <<- value
    },
    detectDuplicatesFn = function(...) list(
      duplicatesList = list(Plasma = data.frame(metabolite_id = c("dup_1", "dup_2"))),
      totalDuplicates = 2L
    ),
    reportDetectionFn = function(...) invisible(NULL),
    reportDetectionErrorFn = function(...) invisible(NULL)
  )
  expect_identical(detection_dispatch$status, "success")
  expect_identical(nrow(duplicate_value$Plasma), 2L)

  output_env <- new.env(parent = emptyenv())
  removed_notifications <- character()
  resolution_dispatch <- run_resolution_observer_shell(
    runResolutionFn = function() list(
      resultText = "Duplicate Resolution Complete",
      totalRemoved = 1L
    ),
    output = output_env,
    setDuplicateInfoFn = function(value) {
      duplicate_value <<- value
    },
    renderTextFn = function(text) text,
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL),
    showNotificationFn = function(message, type, duration = NULL) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    },
    removeNotificationFn = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    }
  )
  expect_identical(resolution_dispatch$status, "success")
  expect_identical(output_env$resolution_results, "Duplicate Resolution Complete")
  expect_identical(duplicate_value, NULL)
  expect_true("metab_dup_resolve_working" %in% removed_notifications)

  resolution_stats_value <- "stats"
  filter_plot_value <- "plot"
  revert_dispatch <- run_revert_observer_shell(
    runRevertFn = function() list(
      previousStateName = "qc_ready",
      resultText = "Reverted to previous state: qc_ready"
    ),
    output = output_env,
    setResolutionStatsFn = function(value) {
      resolution_stats_value <<- value
    },
    setDuplicateInfoFn = function(value) {
      duplicate_value <<- value
    },
    setFilterPlotFn = function(value) {
      filter_plot_value <<- value
    },
    renderTextFn = function(text) text,
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL),
    showNotificationFn = function(message, type) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    }
  )
  expect_identical(revert_dispatch$status, "success")
  expect_identical(resolution_stats_value, NULL)
  expect_identical(duplicate_value, NULL)
  expect_identical(filter_plot_value, NULL)

  observer_capture <- new.env(parent = emptyenv())
  run_resolution_observer(
    workflowData = list(state_manager = "state-manager"),
    omicType = "metabolomics",
    output = output_env,
    setDuplicateInfoFn = function(...) invisible(NULL),
    setResolutionStatsFn = function(...) invisible(NULL),
    setFilterPlotFn = function(...) invisible(NULL),
    reqFn = function(value) value,
    showNotificationFn = function(message, id = NULL, duration = NULL) {
      observer_capture$message <- message
      observer_capture$id <- id
      observer_capture$duration <- duration
      invisible(NULL)
    },
    runResolutionObserverShellFn = function(runResolutionFn, output, setDuplicateInfoFn, ...) {
      observer_capture$shell <- runResolutionFn()
      invisible(NULL)
    },
    runResolutionWorkflowFn = function(...) list(
      resultText = "resolution ok",
      totalRemoved = 2L
    )
  )
  expect_identical(observer_capture$message, "Resolving duplicate features...")
  expect_identical(observer_capture$id, "metab_dup_resolve_working")
  expect_identical(observer_capture$shell$totalRemoved, 2L)
})
