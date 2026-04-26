# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

skip_if_missing_lipid_duplicate_helpers <- function() {
  skip_if_not(
    hasMultiScholaRBinding("runLipidDuplicateModuleServerShell"),
    "mod_lipid_qc_duplicates helper seams are only available on the refactored branch"
  )
}

if (!methods::isClass("LipidomicsAssayData")) {
  methods::setClass(
    "LipidomicsAssayData",
    slots = c(
      lipid_data = "list",
      lipid_id_column = "character",
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
      lipid_data = list(),
      lipid_id_column = "database_identifier",
      annotation_id_column = "lipid_identification",
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

makeLipidDuplicateCharacterizationData <- function(label = "lipid_dup_fixture") {
  methods::new(
    "LipidomicsAssayData",
    lipid_data = list(
      `Positive Mode` = data.frame(
        database_identifier = c(
          paste0(label, "_L1"),
          paste0(label, "_L1"),
          paste0(label, "_L2")
        ),
        lipid_identification = c("Annot 1", "Annot 1", "Annot 2"),
        Sample1 = c(10, 4, 7),
        Sample2 = c(9, 3, 6),
        stringsAsFactors = FALSE
      )
    ),
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
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

makeLipidDuplicateHarness <- function(
  current_s4 = makeLipidDuplicateCharacterizationData(),
  history = c("import_complete", "qc_ready", "lipid_duplicates_resolved")
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

test_that("mod_lipid_qc_duplicates_server preserves detect-resolve-revert public behavior", {
  harness <- makeLipidDuplicateHarness()
  mod_lipid_qc_duplicates_server <- getMultiScholaRBinding("mod_lipid_qc_duplicates_server")

  local_mocked_bindings(
    showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL, type = "default", session = getDefaultReactiveDomain()) {
      recordNotification(harness$capture, ui = ui, id = id, type = type, duration = duration)
    },
    removeNotification = function(id, session = getDefaultReactiveDomain()) {
      harness$capture$removed_notifications <<- c(harness$capture$removed_notifications, id)
      invisible(NULL)
    },
    .package = "shiny"
  )
  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )
  local_mocked_bindings(
    findLipidDuplicateFeatureIDs = function(theObject) {
      expect_s4_class(theObject, "LipidomicsAssayData")
      list(
        `Positive Mode` = data.frame(
          database_identifier = paste0("dup_", seq_len(2)),
          duplicate_count = c(2L, 3L),
          stringsAsFactors = FALSE
        )
      )
    },
    resolveLipidDuplicateFeaturesByIntensity = function(assay_tibble, id_col, sample_cols) {
      expect_identical(id_col, "database_identifier")
      expect_identical(sample_cols, c("Sample1", "Sample2"))
      assay_tibble[c(1, 3), , drop = FALSE]
    },
    updateLipidFiltering = function(...) grid::nullGrob(),
    .env = multiScholaRNamespace()
  )

  testServer(
    mod_lipid_qc_duplicates_server,
    args = list(
      workflow_data = harness$workflow_data,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics"
    ),
    {
      session$setInputs(detect_duplicates = 1)
      session$setInputs(resolve_duplicates = 1)
      session$setInputs(revert_duplicates = 1)
    }
  )

  notice_text <- vapply(harness$capture$notifications, notificationText, character(1))
  expect_length(harness$capture$saved_states, 1L)
  expect_identical(harness$capture$saved_states[[1]]$state_name, "lipid_duplicates_resolved")
  expect_identical(harness$capture$reverted_to, "qc_ready")
  expect_true(any(grepl("Detection complete: 2 duplicate IDs found", notice_text, fixed = TRUE)))
  expect_true(any(grepl("Duplicates resolved: 1 rows removed", notice_text, fixed = TRUE)))
  expect_true(any(grepl("Reverted successfully", notice_text, fixed = TRUE)))
  expect_true("lipid_dup_resolve_working" %in% harness$capture$removed_notifications)
})

test_that("lipid duplicate helper builders preserve initial-state and UI behavior", {
  skip_if_missing_lipid_duplicate_helpers()

  initialize_state <- getMultiScholaRBinding("initializeLipidDuplicateServerState")
  build_summary_ui <- getMultiScholaRBinding("buildLipidDuplicateSummaryUi")
  build_tables_ui <- getMultiScholaRBinding("buildLipidDuplicateTablesUi")

  state <- initialize_state(
    reactiveValFn = function(value) structure(list(value = value), class = "mock_reactive")
  )
  expect_identical(names(state), c("duplicateInfo", "resolutionStats", "filterPlot"))
  expect_s3_class(state$duplicateInfo, "mock_reactive")

  summary_html <- htmltools::renderTags(
    build_summary_ui(list(
      `Positive Mode` = data.frame(
        database_identifier = c("dup_1", "dup_2"),
        duplicate_count = c(2L, 3L),
        stringsAsFactors = FALSE
      ),
      `Negative Mode` = NULL
    ))
  )$html
  expect_match(summary_html, "Positive Mode: 2 duplicate IDs", fixed = TRUE)
  expect_match(summary_html, "Negative Mode: 0 duplicate IDs", fixed = TRUE)

  tables_html <- htmltools::renderTags(
    build_tables_ui(
      dupList = list(
        `Positive Mode` = data.frame(database_identifier = c("dup_1", "dup_2")),
        `Negative Mode` = NULL
      ),
      ns = function(id) paste0("dup-", id)
    )
  )$html
  expect_match(tables_html, "Positive Mode", fixed = TRUE)
  expect_match(tables_html, "dup-dup_table_Positive_Mode", fixed = TRUE)
  expect_false(grepl("Negative Mode", tables_html, fixed = TRUE))
})

test_that("lipid duplicate output helpers preserve render registration behavior", {
  skip_if_missing_lipid_duplicate_helpers()

  register_summary_output <- getMultiScholaRBinding("registerLipidDuplicateSummaryOutput")
  register_tables_output <- getMultiScholaRBinding("registerLipidDuplicateTablesOutput")
  register_filter_plot_output <- getMultiScholaRBinding("registerLipidDuplicateFilterPlotOutput")

  output_env <- new.env(parent = emptyenv())
  duplicate_info <- list(
    `Positive Mode` = data.frame(database_identifier = c("dup_1", "dup_2"))
  )

  register_summary_output(
    output = output_env,
    duplicateInfoVal = function() duplicate_info,
    renderUiFn = function(expr) htmltools::renderTags(expr)$html,
    buildSummaryUiFn = function(dupList) shiny::tags$span(paste(names(dupList), collapse = ","))
  )
  expect_match(output_env$duplicate_summary, "Positive Mode", fixed = TRUE)

  register_tables_output(
    output = output_env,
    duplicateInfoVal = function() duplicate_info,
    ns = function(id) paste0("dup-", id),
    renderUiFn = function(expr) htmltools::renderTags(expr)$html,
    buildTablesUiFn = function(dupList, ns) shiny::tags$span(ns(paste0("tables-", names(dupList)[[1]])))
  )
  expect_match(output_env$duplicate_tables, "dup-tables-Positive Mode", fixed = TRUE)

  drawn_grob <- NULL
  register_filter_plot_output(
    output = output_env,
    filterPlotVal = function() grid::nullGrob(),
    renderPlotFn = function(expr) {
      force(expr)
      "plot-grob-ok"
    },
    reqFn = function(value) value,
    gridDrawFn = function(value) {
      drawn_grob <<- value
      invisible(NULL)
    },
    printFn = function(...) stop("printFn should not be called for grobs")
  )
  expect_identical(output_env$filter_plot, "plot-grob-ok")
  expect_true(inherits(drawn_grob, "grob"))

  printed_plot <- NULL
  register_filter_plot_output(
    output = output_env,
    filterPlotVal = function() structure(list(token = "ggplot"), class = "ggplot"),
    renderPlotFn = function(expr) {
      force(expr)
      "plot-ggplot-ok"
    },
    reqFn = function(value) value,
    gridDrawFn = function(...) stop("gridDrawFn should not be called for ggplot"),
    printFn = function(value) {
      printed_plot <<- value
      invisible(NULL)
    }
  )
  expect_identical(output_env$filter_plot, "plot-ggplot-ok")
  expect_identical(printed_plot$token, "ggplot")
})

test_that("lipid duplicate detection helpers preserve duplicate summary behavior", {
  skip_if_missing_lipid_duplicate_helpers()

  handle_detection <- getMultiScholaRBinding("handleLipidDuplicateDetection")
  apply_detection <- getMultiScholaRBinding("applyLipidDuplicateDetectionResult")
  current_s4 <- makeLipidDuplicateCharacterizationData()
  duplicate_value <- NULL
  notifications <- list()

  detection_result <- handle_detection(
    workflowData = list(
      state_manager = list(getState = function() current_s4)
    ),
    reqFn = function(value) value,
    findDuplicatesFn = function(theObject) {
      expect_identical(theObject, current_s4)
      list(
        `Positive Mode` = data.frame(database_identifier = c("dup_1", "dup_2")),
        `Negative Mode` = NULL
      )
    },
    logInfoFn = function(...) invisible(NULL)
  )

  expect_identical(detection_result$totalDuplicates, 2)
  expect_identical(detection_result$notificationType, "warning")

  apply_detection(
    detectionResult = detection_result,
    duplicateInfoVal = function(value) {
      duplicate_value <<- value
    },
    showNotificationFn = function(message, type) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    }
  )

  expect_identical(duplicate_value, detection_result$duplicatesList)
  expect_identical(notifications[[1]]$message, "Detection complete: 2 duplicate IDs found")
  expect_identical(notifications[[1]]$type, "warning")
})

test_that("lipid duplicate resolution and revert helpers preserve save and fallback behavior", {
  skip_if_missing_lipid_duplicate_helpers()

  handle_resolution <- getMultiScholaRBinding("handleLipidDuplicateResolution")
  apply_resolution <- getMultiScholaRBinding("applyLipidDuplicateResolutionResult")
  handle_revert <- getMultiScholaRBinding("handleLipidDuplicateRevert")
  apply_revert <- getMultiScholaRBinding("applyLipidDuplicateRevertResult")

  harness <- makeLipidDuplicateHarness()
  current_s4 <- harness$current_s4
  resolution_stats_value <- NULL
  duplicate_value <- list(existing = TRUE)
  filter_plot_value <- "old-plot"
  output_env <- new.env(parent = emptyenv())
  notifications <- list()
  removed_notifications <- character()

  resolution_result <- handle_resolution(
    workflowData = harness$workflow_data,
    omicType = "lipidomics",
    reqFn = function(value) value,
    resolveDuplicatesFn = function(assay_tibble, id_col, sample_cols) {
      expect_identical(id_col, "database_identifier")
      expect_identical(sample_cols, c("Sample1", "Sample2"))
      assay_tibble[c(1, 3), , drop = FALSE]
    },
    updateFilteringFn = function(...) grid::nullGrob(),
    logWarnFn = function(...) invisible(NULL)
  )

  apply_resolution(
    resolutionResult = resolution_result,
    resolutionStatsVal = function(value) {
      resolution_stats_value <<- value
    },
    duplicateInfoVal = function(value) {
      duplicate_value <<- value
    },
    filterPlotVal = function(value) {
      filter_plot_value <<- value
    },
    output = output_env,
    renderTextFn = function(text) text,
    logInfoFn = function(...) invisible(NULL),
    removeNotificationFn = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    },
    showNotificationFn = function(message, type) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    }
  )

  expect_identical(resolution_result$totalRemoved, 1)
  expect_identical(output_env$resolution_results, resolution_result$resultText)
  expect_identical(duplicate_value, NULL)
  expect_true(inherits(filter_plot_value, "grob"))
  expect_identical(resolution_stats_value[["Positive Mode"]]$removed, 1L)
  expect_true("lipid_dup_resolve_working" %in% removed_notifications)
  expect_identical(notifications[[1]]$message, "Duplicates resolved: 1 rows removed")

  no_numeric_s4 <- makeLipidDuplicateCharacterizationData("no_numeric")
  no_numeric_s4@lipid_data[[1]] <- data.frame(
    database_identifier = c("L1", "L1"),
    lipid_identification = c("Annot 1", "Annot 1"),
    label = c("keep", "drop"),
    stringsAsFactors = FALSE
  )
  warn_messages <- character()

  fallback_result <- handle_resolution(
    workflowData = list(
      state_manager = list(
        getState = function() no_numeric_s4,
        saveState = function(...) invisible(NULL)
      ),
      config_list = list()
    ),
    omicType = "lipidomics",
    reqFn = function(value) value,
    resolveDuplicatesFn = function(...) stop("resolveDuplicatesFn should not be called"),
    updateFilteringFn = function(...) stop("plot failed"),
    logWarnFn = function(message) {
      warn_messages <<- c(warn_messages, message)
    }
  )

  expect_identical(fallback_result$totalRemoved, 0)
  expect_true(any(grepl("No numeric columns found in assay: Positive Mode", warn_messages, fixed = TRUE)))
  expect_true(any(grepl("Could not generate QC plot: plot failed", warn_messages, fixed = TRUE)))

  revert_result <- handle_revert(
    workflowData = harness$workflow_data,
    reqFn = function(value) value,
    logInfoFn = function(...) invisible(NULL)
  )
  apply_revert(
    revertResult = revert_result,
    resolutionStatsVal = function(value) {
      resolution_stats_value <<- value
    },
    duplicateInfoVal = function(value) {
      duplicate_value <<- value
    },
    filterPlotVal = function(value) {
      filter_plot_value <<- value
    },
    output = output_env,
    renderTextFn = function(text) text,
    showNotificationFn = function(message, type) {
      notifications <<- c(notifications, list(list(message = message, type = type)))
    }
  )

  expect_identical(harness$capture$reverted_to, "qc_ready")
  expect_identical(resolution_stats_value, NULL)
  expect_identical(duplicate_value, NULL)
  expect_identical(filter_plot_value, NULL)
  expect_match(output_env$resolution_results, "Reverted to previous state: qc_ready", fixed = TRUE)
  expect_identical(notifications[[2]]$message, "Reverted successfully")
})

test_that("lipid duplicate observer and shell helpers preserve registration wiring", {
  skip_if_missing_lipid_duplicate_helpers()

  register_table_observer <- getMultiScholaRBinding("registerLipidDuplicateTableObserver")
  register_detect_observer <- getMultiScholaRBinding("registerLipidDuplicateDetectObserver")
  register_resolve_observer <- getMultiScholaRBinding("registerLipidDuplicateResolveObserver")
  register_revert_observer <- getMultiScholaRBinding("registerLipidDuplicateRevertObserver")
  register_server_bindings <- getMultiScholaRBinding("registerLipidDuplicateServerBindings")
  run_server_shell <- getMultiScholaRBinding("runLipidDuplicateModuleServerShell")

  output_env <- new.env(parent = emptyenv())
  duplicate_value <- list(`Positive Mode` = data.frame(database_identifier = "dup_1"))
  register_args <- NULL

  register_table_observer(
    output = output_env,
    duplicateInfoVal = function() duplicate_value,
    observeFn = function(handler) {
      force(handler)
      invisible(NULL)
    },
    reqFn = function(value) value,
    registerTableOutputsFn = function(output, dupList) {
      register_args <<- list(output = output, dupList = dupList)
      "dup_table_Positive_Mode"
    }
  )
  expect_identical(register_args$output, output_env)
  expect_identical(register_args$dupList, duplicate_value)

  event_trace <- character()
  register_detect_observer(
    input = list(detect_duplicates = 1L),
    workflowData = list(state_manager = "state-manager"),
    duplicateInfoVal = function(...) invisible(NULL),
    observeEventFn = function(event, handler) {
      event_trace <<- c(event_trace, paste0("detect:", event))
      force(handler)
      invisible(NULL)
    },
    reqFn = function(value) value,
    handleDetectionFn = function(...) list(
      duplicatesList = duplicate_value,
      totalDuplicates = 1L,
      notificationMessage = "Detection complete: 1 duplicate IDs found",
      notificationType = "warning"
    ),
    applyDetectionResultFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL)
  )

  register_resolve_observer(
    input = list(resolve_duplicates = 2L),
    workflowData = list(state_manager = "state-manager"),
    omicType = "lipidomics",
    resolutionStatsVal = function(...) invisible(NULL),
    duplicateInfoVal = function(...) invisible(NULL),
    filterPlotVal = function(...) invisible(NULL),
    output = output_env,
    observeEventFn = function(event, handler) {
      event_trace <<- c(event_trace, paste0("resolve:", event))
      force(handler)
      invisible(NULL)
    },
    reqFn = function(value) value,
    showNotificationFn = function(...) invisible(NULL),
    handleResolutionFn = function(...) list(
      statsList = list(),
      qcPlot = NULL,
      totalRemoved = 0L,
      resultText = "Duplicate Resolution Complete"
    ),
    applyResolutionResultFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL),
    removeNotificationFn = function(...) invisible(NULL)
  )

  register_revert_observer(
    input = list(revert_duplicates = 3L),
    workflowData = list(state_manager = "state-manager"),
    resolutionStatsVal = function(...) invisible(NULL),
    duplicateInfoVal = function(...) invisible(NULL),
    filterPlotVal = function(...) invisible(NULL),
    output = output_env,
    observeEventFn = function(event, handler) {
      event_trace <<- c(event_trace, paste0("revert:", event))
      force(handler)
      invisible(NULL)
    },
    handleRevertFn = function(...) list(
      resultText = "Reverted to previous state: qc_ready",
      notificationMessage = "Reverted successfully",
      notificationType = "message"
    ),
    applyRevertResultFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL)
  )

  expect_identical(event_trace, c("detect:1", "resolve:2", "revert:3"))

  binding_trace <- character()
  register_server_bindings(
    input = list(detect_duplicates = 1L, resolve_duplicates = 2L, revert_duplicates = 3L),
    output = output_env,
    workflowData = list(state_manager = "state-manager"),
    omicType = "lipidomics",
    ns = function(id) paste0("dup-", id),
    duplicateInfoVal = function(...) invisible(NULL),
    resolutionStatsVal = function(...) invisible(NULL),
    filterPlotVal = function(...) invisible(NULL),
    registerDetectObserverFn = function(...) binding_trace <<- c(binding_trace, "detect"),
    registerSummaryOutputFn = function(...) binding_trace <<- c(binding_trace, "summary"),
    registerTablesOutputFn = function(...) binding_trace <<- c(binding_trace, "tables"),
    registerTableObserverFn = function(...) binding_trace <<- c(binding_trace, "table-observer"),
    registerResolveObserverFn = function(...) binding_trace <<- c(binding_trace, "resolve"),
    registerRevertObserverFn = function(...) binding_trace <<- c(binding_trace, "revert"),
    registerFilterPlotOutputFn = function(...) binding_trace <<- c(binding_trace, "plot")
  )
  expect_identical(
    binding_trace,
    c("detect", "summary", "tables", "table-observer", "resolve", "revert", "plot")
  )

  shell_trace <- list()
  shell_state <- run_server_shell(
    input = list(),
    output = output_env,
    session = list(ns = function(id) paste0("dup-", id)),
    workflowData = list(state_manager = "state-manager"),
    omicType = "lipidomics",
    initializeServerStateFn = function() list(
      duplicateInfo = "dup-state",
      resolutionStats = "stats-state",
      filterPlot = "plot-state"
    ),
    registerServerBindingsFn = function(input, output, workflowData, omicType, ns, duplicateInfoVal, resolutionStatsVal, filterPlotVal) {
      shell_trace <<- list(
        input = input,
        output = output,
        workflowData = workflowData,
        omicType = omicType,
        ns_value = ns("tables"),
        duplicateInfoVal = duplicateInfoVal,
        resolutionStatsVal = resolutionStatsVal,
        filterPlotVal = filterPlotVal
      )
      invisible(output)
    }
  )

  expect_identical(shell_state$duplicateInfo, "dup-state")
  expect_identical(shell_trace$ns_value, "dup-tables")
  expect_identical(shell_trace$omicType, "lipidomics")
  expect_identical(shell_trace$duplicateInfoVal, "dup-state")
})
