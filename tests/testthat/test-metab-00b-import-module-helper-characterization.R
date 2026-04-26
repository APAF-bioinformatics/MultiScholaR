library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingMetabImportTargetFiles <- function() {
  required_paths <- c(
    "R/mod_metab_import_column_helpers.R",
    "R/mod_metab_import_display_helpers.R",
    "R/mod_metab_import_registration_helpers.R",
    "R/mod_metab_import_processing_helpers.R",
    "R/mod_metab_import_server_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only metab import helper file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingMetabImportTargetFiles()

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
    file.path(repo_root, "R", "mod_metab_import_column_helpers.R"),
    file.path(repo_root, "R", "mod_metab_import_display_helpers.R"),
    file.path(repo_root, "R", "mod_metab_import_registration_helpers.R"),
    file.path(repo_root, "R", "mod_metab_import_processing_helpers.R"),
    file.path(repo_root, "R", "mod_metab_import_server_helpers.R"),
    file.path(repo_root, "R", "mod_metab_import.R")
  ),
  symbols = c(
    "resolveMetabImportColumnName",
    "resolveMetabImportSampleColumns",
    "buildMetabImportWorkflowPayload",
    "applyMetabImportWorkflowPayload",
    "prepareMetabImportAssaySelectionState",
    "applyMetabImportAssaySelectionState",
    "handleMetabImportAssayFileSelection",
    "finalizeMetabImportProcessingFeedback",
    "finalizeMetabImportAssaySelectionError",
    "runMetabImportAssaySelection",
    "setupMetabImportShinyFiles",
    "setupMetabImportColumnAccessors",
    "setupMetabImportAssaySelectionCallback",
    "setupMetabImportProcessingObserver",
    "setupMetabImportFileLoadedOutput",
    "setupMetabImportStatusOutput",
    "setupMetabImportValidationSummaryOutput",
    "setupMetabImportFormatDetectionStatusOutput",
    "setupMetabImportMetaboliteIdStatusOutput",
    "setupMetabImportAnnotationStatusOutput",
    "setupMetabImportSampleColumnsDisplayOutput",
    "setupMetabImportAvailableColumnsDisplayOutput",
    "setupMetabImportCustomMetaboliteIdStatusOutput",
    "setupMetabImportCustomAnnotationStatusOutput",
    "buildMetabImportFormatDetectionStatus",
    "buildMetabImportMetaboliteIdStatus",
    "buildMetabImportCustomMetaboliteIdStatus",
    "buildMetabImportCustomAnnotationStatus",
    "buildMetabImportAnnotationStatus",
    "buildMetabImportSampleColumnsDisplay",
    "buildMetabImportAvailableColumnsDisplay",
    "buildMetabImportValidationSummary",
    "buildMetabImportStatus",
    "runMetabImportProcessing"
  ),
  env = environment()
)

test_that("metabolomics import column resolver preserves case-insensitive lookup and fallback", {
  headers <- c("Peak ID", "Metabolite_Name", "Sample_A")

  expect_identical(resolveMetabImportColumnName(headers, "peak id"), "Peak ID")
  expect_identical(resolveMetabImportColumnName(headers, "METABOLITE_name"), "Metabolite_Name")
  expect_identical(resolveMetabImportColumnName(headers, "missing_column"), "missing_column")
  expect_identical(resolveMetabImportColumnName(headers, ""), "")
  expect_null(resolveMetabImportColumnName(headers, NULL))
})

test_that("metabolomics import sample-column resolver preserves custom pattern and fallback order", {
  assay_data <- data.frame(
    metabolite_id = c("M1", "M2"),
    Sample_A = c(10, 11),
    sample_b = c(20, 21),
    Annotation = c("A", "B"),
    check.names = FALSE
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "custom",
      sampleColsPattern = "^sample",
      importResult = list(sample_columns = c("ignored"))
    ),
    c("Sample_A", "sample_b")
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "custom",
      sampleColsPattern = "^missing",
      importResult = list(sample_columns = c("Imported_A", "Imported_B"))
    ),
    c("Imported_A", "Imported_B")
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "msdial",
      sampleColsPattern = "^sample",
      importResult = list(sample_columns = c("Imported_A", "Imported_B"))
    ),
    c("Imported_A", "Imported_B")
  )

  expect_identical(
    resolveMetabImportSampleColumns(
      assayData = assay_data,
      vendorFormat = "custom",
      sampleColsPattern = "^missing",
      importResult = NULL
    ),
    c("Sample_A", "sample_b")
  )
})

test_that("metabolomics import workflow payload builder preserves assay assembly and sanitization contracts", {
  assay_one <- data.frame(
    metabolite_id = c("M1", "M2", "M2"),
    "Sample A" = c(10, 11, 12),
    "Sample-B" = c(20, 21, 22),
    Annotation = c("A", "B", "B"),
    check.names = FALSE
  )
  assay_two <- data.frame(
    feature_id = c("N1", "N2"),
    "Sample A" = c(5, 6),
    "Sample-B" = c(7, 8),
    check.names = FALSE
  )
  fixed_time <- as.POSIXct("2026-04-16 10:15:00", tz = "UTC")

  payload <- buildMetabImportWorkflowPayload(
    assay1Name = "LCMS_Pos",
    assay1Data = assay_one,
    assay2File = "assay2.tsv",
    assay2Name = "LCMS_Neg",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    metaboliteCol = "metabolite_id",
    annotationCol = "",
    sampleCols = c("Sample A", "Sample-B"),
    sanitizeNames = TRUE,
    isPattern = "",
    assay2Importer = function(path) {
      expect_identical(path, "assay2.tsv")
      list(data = assay_two)
    },
    cleanNamesFn = function(x) {
      expect_identical(x, c("Sample A", "Sample-B"))
      c("sample_a", "sample_b")
    },
    mapAssaysFn = lapply,
    timestampFn = function() fixed_time
  )

  expect_true(payload$sampleNamesSanitized)
  expect_identical(payload$dataFormat, "custom")
  expect_identical(payload$sampleCols, c("sample_a", "sample_b"))
  expect_identical(payload$columnMapping$metabolite_id_col, "metabolite_id")
  expect_null(payload$columnMapping$annotation_col)
  expect_identical(payload$columnMapping$sample_columns, c("sample_a", "sample_b"))
  expect_true(is.na(payload$columnMapping$is_pattern))
  expect_identical(names(payload$assayList), c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(
    names(payload$assayList$LCMS_Pos),
    c("metabolite_id", "sample_a", "sample_b", "Annotation")
  )
  expect_identical(
    names(payload$assayList$LCMS_Neg),
    c("feature_id", "sample_a", "sample_b")
  )
  expect_identical(payload$processingLog$timestamp, fixed_time)
  expect_identical(payload$processingLog$n_assays, 2L)
  expect_identical(payload$processingLog$assay_names, c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(payload$processingLog$detected_format, "custom")
  expect_identical(unname(payload$processingLog$n_metabolites), c(2L, 2L))
  expect_identical(payload$processingLog$n_samples, 2L)
})

test_that("metabolomics import workflow payload application preserves workflow commit and finalization contracts", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- NULL
  workflow_data$data_format <- NULL
  workflow_data$data_type <- NULL
  workflow_data$column_mapping <- NULL
  workflow_data$processing_log <- list(setup_import = NULL)
  workflow_data$tab_status <- list(setup_import = "incomplete", downstream = "pending")

  workflow_type_calls <- character()
  workflow_data$state_manager <- list(
    setWorkflowType = function(type) {
      workflow_type_calls <<- c(workflow_type_calls, type)
      invisible(type)
    }
  )

  assay_one <- data.frame(metabolite_id = c("M1", "M2"), sample_a = c(1, 2))
  assay_two <- data.frame(feature_id = c("N1"), sample_a = 3)
  workflow_payload <- list(
    assayList = list(LCMS_Pos = assay_one, LCMS_Neg = assay_two),
    dataFormat = "custom",
    columnMapping = list(
      metabolite_id_col = "metabolite_id",
      annotation_col = NULL,
      sample_columns = "sample_a",
      is_pattern = NA_character_
    ),
    processingLog = list(
      timestamp = as.POSIXct("2026-04-16 11:00:00", tz = "UTC"),
      n_assays = 2L,
      assay_names = c("LCMS_Pos", "LCMS_Neg"),
      detected_format = "custom",
      n_metabolites = c(2L, 1L),
      n_samples = 1L
    )
  )
  log_messages <- character()

  result <- applyMetabImportWorkflowPayload(
    workflowData = workflow_data,
    workflowPayload = workflow_payload,
    logInfo = function(message) {
      log_messages <<- c(log_messages, message)
    }
  )

  expect_identical(workflow_data$data_tbl, workflow_payload$assayList)
  expect_identical(workflow_data$data_format, "custom")
  expect_identical(workflow_data$data_type, "metabolite")
  expect_identical(workflow_data$column_mapping, workflow_payload$columnMapping)
  expect_identical(workflow_type_calls, "metabolomics_standard")
  expect_identical(workflow_data$processing_log$setup_import, workflow_payload$processingLog)
  expect_identical(workflow_data$tab_status$setup_import, "complete")
  expect_identical(workflow_data$tab_status$downstream, "pending")
  expect_identical(result$assayCount, 2L)
  expect_identical(result$totalMetabolites, 3L)
  expect_identical(result$updatedStatus, workflow_data$tab_status)
  expect_identical(
    log_messages,
    c(
      "Workflow type set to: metabolomics_standard",
      "Metabolomics import complete: 2 assays, 3 total metabolites"
    )
  )
})

test_that("metabolomics import assay-selection helper preserves header parsing and input-default contracts", {
  assay_file <- tempfile(fileext = ".csv")
  writeLines('"Peak ID","Name","Sample_1"', assay_file)
  on.exit(unlink(assay_file), add = TRUE)

  detect_calls <- list()
  import_calls <- character()

  selection_state <- prepareMetabImportAssaySelectionState(
    assay1File = assay_file,
    detectFormatFn = function(headers, filename) {
      detect_calls <<- list(headers = headers, filename = filename)
      list(format = "msdial", confidence = 0.95)
    },
    importers = list(
      msdial = function(path) {
        import_calls <<- c(import_calls, path)
        list(
          data = data.frame(
            "Peak ID" = "P1",
            Name = "Met1",
            Sample_1 = 123,
            check.names = FALSE
          ),
          detected_columns = list(
            metabolite_id = "Peak ID",
            annotation = "Name"
          ),
          sample_columns = "Sample_1",
          is_pattern = "^IS_"
        )
      }
    ),
    defaultImporter = function(path) {
      stop(sprintf("unexpected default importer call for %s", path))
    }
  )

  expect_identical(selection_state$headers, c("Peak ID", "Name", "Sample_1"))
  expect_identical(detect_calls$headers, c("Peak ID", "Name", "Sample_1"))
  expect_identical(detect_calls$filename, basename(assay_file))
  expect_identical(import_calls, assay_file)
  expect_identical(selection_state$formatInfo$format, "msdial")
  expect_equal(selection_state$formatInfo$confidence, 0.95, tolerance = 1e-9)
  expect_identical(selection_state$metaboliteIdChoices, c("Peak ID", "Name", "Sample_1"))
  expect_identical(selection_state$selectedMetaboliteId, "Peak ID")
  expect_identical(unname(selection_state$annotationChoices), c("", "Peak ID", "Name", "Sample_1"))
  expect_identical(names(selection_state$annotationChoices)[1], "(None)")
  expect_identical(selection_state$selectedAnnotation, "Name")
  expect_identical(selection_state$isPattern, "^IS_")
})

test_that("metabolomics import assay-selection helper preserves header-read failure and default-import fallback contracts", {
  expect_error(
    prepareMetabImportAssaySelectionState(
      assay1File = "broken.tsv",
      readHeadersFn = function(path) character(0),
      detectFormatFn = function(headers, filename) {
        stop("detectFormatFn should not be called when headers are missing")
      }
    ),
    "Could not read headers from file"
  )

  default_import_calls <- character()

  fallback_state <- prepareMetabImportAssaySelectionState(
    assay1File = "fallback.tsv",
    readHeadersFn = function(path) c("Feature", "Sample_1"),
    detectFormatFn = function(headers, filename) {
      expect_identical(headers, c("Feature", "Sample_1"))
      expect_identical(filename, "fallback.tsv")
      list(format = "unknown", confidence = 0.2)
    },
    importers = list(
      msdial = function(path) {
        stop(sprintf("unexpected named importer call for %s", path))
      }
    ),
    defaultImporter = function(path) {
      default_import_calls <<- c(default_import_calls, path)
      list(
        data = data.frame(Feature = "F1", Sample_1 = 1, check.names = FALSE),
        detected_columns = list(
          metabolite_id = "Feature",
          annotation = ""
        ),
        sample_columns = "Sample_1",
        is_pattern = NA_character_
      )
    }
  )

  expect_identical(default_import_calls, "fallback.tsv")
  expect_identical(fallback_state$selectedMetaboliteId, "Feature")
  expect_identical(fallback_state$selectedAnnotation, "")
  expect_null(fallback_state$isPattern)
})

test_that("metabolomics import assay-selection application helper preserves local-data commit and update contracts", {
  local_data <- new.env(parent = emptyenv())
  local_data$all_headers <- NULL
  local_data$detected_format <- NULL
  local_data$format_confidence <- NULL
  local_data$assay1_import_result <- NULL
  local_data$assay1_data <- NULL

  assay_data <- data.frame(
    "Peak ID" = c("P1", "P2"),
    Sample_1 = c(10, 12),
    Sample_2 = c(11, 13),
    check.names = FALSE
  )
  import_state <- list(
    headers = c("Peak ID", "Sample_1", "Sample_2"),
    formatInfo = list(format = "msdial", confidence = 0.93),
    importResult = list(
      data = assay_data,
      detected_columns = list(
        metabolite_id = "Peak ID",
        annotation = ""
      ),
      sample_columns = c("Sample_1", "Sample_2"),
      is_pattern = "^IS_"
    ),
    metaboliteIdChoices = c("Peak ID", "Sample_1", "Sample_2"),
    selectedMetaboliteId = "Peak ID",
    annotationChoices = c("(None)" = "", "Peak ID", "Sample_1", "Sample_2"),
    selectedAnnotation = "",
    isPattern = "^IS_"
  )
  select_updates <- list()
  text_updates <- list()
  log_messages <- character()

  result <- applyMetabImportAssaySelectionState(
    localData = local_data,
    importState = import_state,
    session = "fake-session",
    updateSelectInputFn = function(session, inputId, choices, selected) {
      select_updates[[length(select_updates) + 1]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
    },
    updateTextInputFn = function(session, inputId, value) {
      text_updates[[length(text_updates) + 1]] <<- list(
        session = session,
        inputId = inputId,
        value = value
      )
    },
    logInfoFn = function(message) {
      log_messages <<- c(log_messages, message)
    }
  )

  expect_identical(local_data$all_headers, import_state$headers)
  expect_identical(local_data$detected_format, "msdial")
  expect_equal(local_data$format_confidence, 0.93, tolerance = 1e-9)
  expect_identical(local_data$assay1_import_result, import_state$importResult)
  expect_identical(local_data$assay1_data, assay_data)
  expect_identical(
    select_updates,
    list(
      list(
        session = "fake-session",
        inputId = "metabolite_id_col",
        choices = c("Peak ID", "Sample_1", "Sample_2"),
        selected = "Peak ID"
      ),
      list(
        session = "fake-session",
        inputId = "annotation_col",
        choices = c("(None)" = "", "Peak ID", "Sample_1", "Sample_2"),
        selected = ""
      )
    )
  )
  expect_identical(
    text_updates,
    list(
      list(
        session = "fake-session",
        inputId = "is_pattern",
        value = "^IS_"
      )
    )
  )
  expect_identical(
    log_messages,
    "Imported assay 1: 2 rows, 3 columns, format: msdial"
  )
  expect_identical(result$importedRows, 2L)
  expect_identical(result$importedCols, 3L)
  expect_identical(result$detectedFormat, "msdial")
  expect_equal(result$confidence, 0.93, tolerance = 1e-9)
})

test_that("metabolomics import assay-file selection helper preserves parse, path update, and import trigger contracts", {
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_file <- NULL
  local_data$assay2_file <- NULL
  output <- new.env(parent = emptyenv())

  parse_calls <- list()
  selected_paths <- character()

  result <- handleMetabImportAssayFileSelection(
    selectedInput = "selected-assay",
    volumes = c(home = "/tmp"),
    localData = local_data,
    localField = "assay1_file",
    output = output,
    outputId = "assay1_path",
    onSelected = function(path) {
      selected_paths <<- c(selected_paths, path)
    },
    parseFilePathsFn = function(volumes, selectedInput) {
      parse_calls <<- list(volumes = volumes, selectedInput = selectedInput)
      data.frame(datapath = "/tmp/assay1.csv", stringsAsFactors = FALSE)
    },
    renderTextFn = function(path) {
      paste("rendered:", path)
    }
  )

  expect_identical(parse_calls$volumes, c(home = "/tmp"))
  expect_identical(parse_calls$selectedInput, "selected-assay")
  expect_identical(result, "/tmp/assay1.csv")
  expect_identical(local_data$assay1_file, "/tmp/assay1.csv")
  expect_identical(output$assay1_path, "rendered: /tmp/assay1.csv")
  expect_identical(selected_paths, "/tmp/assay1.csv")

  assay2_result <- handleMetabImportAssayFileSelection(
    selectedInput = "selected-assay-2",
    volumes = c(home = "/tmp"),
    localData = local_data,
    localField = "assay2_file",
    output = output,
    outputId = "assay2_path",
    parseFilePathsFn = function(volumes, selectedInput) {
      expect_identical(volumes, c(home = "/tmp"))
      expect_identical(selectedInput, "selected-assay-2")
      data.frame(datapath = "/tmp/assay2.csv", stringsAsFactors = FALSE)
    },
    renderTextFn = function(path) {
      paste("rendered:", path)
    }
  )

  expect_identical(assay2_result, "/tmp/assay2.csv")
  expect_identical(local_data$assay2_file, "/tmp/assay2.csv")
  expect_identical(output$assay2_path, "rendered: /tmp/assay2.csv")
  expect_identical(selected_paths, "/tmp/assay1.csv")

  expect_null(
    handleMetabImportAssayFileSelection(
      selectedInput = NULL,
      volumes = c(home = "/tmp"),
      localData = local_data,
      localField = "assay1_file",
      output = output,
      outputId = "assay1_path",
      parseFilePathsFn = function(...) {
        stop("parseFilePathsFn should not be called for NULL input")
      }
    )
  )
})

test_that("metabolomics shinyFiles setup helper preserves volume fallback, chooser wiring, and observer delegation contracts", {
  input <- new.env(parent = emptyenv())
  input$assay1_file <- NULL
  input$assay2_file <- NULL
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())

  choose_calls <- list()
  handle_calls <- list()
  observer_handlers <- list()
  import_calls <- integer()
  log_messages <- character()
  error_logs <- character()

  setup_result <- setupMetabImportShinyFiles(
    input = input,
    output = output,
    session = "fake-session",
    volumes = NULL,
    localData = local_data,
    importDataFn = function() {
      import_calls <<- c(import_calls, 1L)
    },
    getVolumesFn = function() {
      c(home = "/tmp", data = "/data")
    },
    shinyFileChooseFn = function(input, id, roots, session, filetypes) {
      choose_calls[[length(choose_calls) + 1]] <<- list(
        input = input,
        id = id,
        roots = roots,
        session = session,
        filetypes = filetypes
      )
    },
    observeEventFn = function(eventExpr, handlerExpr) {
      handler_env <- parent.frame()
      handler_expr <- substitute(handlerExpr)
      observer_handlers[[length(observer_handlers) + 1]] <<- function() {
        eval(handler_expr, envir = handler_env)
      }
      invisible(NULL)
    },
    handleAssayFileSelectionFn = function(
        selectedInput,
        volumes,
        localData,
        localField,
        output,
        outputId,
        onSelected = NULL,
        ...
    ) {
      handle_calls[[length(handle_calls) + 1]] <<- list(
        selectedInput = selectedInput,
        volumes = volumes,
        localData = localData,
        localField = localField,
        output = output,
        outputId = outputId,
        hasOnSelected = is.function(onSelected)
      )

      if (identical(selectedInput, "assay2-selected")) {
        stop("bad assay 2 path")
      }

      if (is.function(onSelected)) {
        onSelected("/tmp/assay1.csv")
      }

      invisible(selectedInput)
    },
    logMessageFn = function(message) {
      log_messages <<- c(log_messages, message)
    },
    logErrorFn = function(message) {
      error_logs <<- c(error_logs, message)
    }
  )

  expect_identical(setup_result, c(home = "/tmp", data = "/data"))
  expect_identical(
    log_messages,
    c(
      "   mod_metab_import_server: volumes is NULL, creating from getVolumes()",
      "   mod_metab_import_server: volumes type = character, length = 2",
      "   mod_metab_import_server: volumes names = home, data"
    )
  )
  expect_length(choose_calls, 2L)
  expect_identical(choose_calls[[1]]$id, "assay1_file")
  expect_identical(choose_calls[[2]]$id, "assay2_file")
  expect_identical(choose_calls[[1]]$roots, c(home = "/tmp", data = "/data"))
  expect_identical(choose_calls[[2]]$roots, c(home = "/tmp", data = "/data"))
  expect_identical(choose_calls[[1]]$session, "fake-session")
  expect_identical(
    choose_calls[[1]]$filetypes,
    c("tsv", "tab", "txt", "csv", "xlsx", "parquet")
  )
  expect_length(observer_handlers, 2L)

  input$assay1_file <- "assay1-selected"
  observer_handlers[[1]]()

  input$assay2_file <- "assay2-selected"
  observer_handlers[[2]]()

  expect_identical(import_calls, 1L)
  expect_length(handle_calls, 2L)
  expect_identical(handle_calls[[1]]$selectedInput, "assay1-selected")
  expect_identical(handle_calls[[1]]$volumes, c(home = "/tmp", data = "/data"))
  expect_identical(handle_calls[[1]]$localData, local_data)
  expect_identical(handle_calls[[1]]$output, output)
  expect_identical(handle_calls[[1]]$localField, "assay1_file")
  expect_identical(handle_calls[[1]]$outputId, "assay1_path")
  expect_true(handle_calls[[1]]$hasOnSelected)
  expect_identical(handle_calls[[2]]$selectedInput, "assay2-selected")
  expect_identical(handle_calls[[2]]$localField, "assay2_file")
  expect_identical(handle_calls[[2]]$outputId, "assay2_path")
  expect_false(handle_calls[[2]]$hasOnSelected)
  expect_identical(error_logs, "Error parsing file path: bad assay 2 path")
})

test_that("metabolomics import processing feedback helper preserves success and error finalization contracts", {
  removed_ids <- character()
  shown_notifications <- list()
  error_logs <- character()

  record_notification <- function(message, type = NULL, duration = NULL) {
    shown_notifications[[length(shown_notifications) + 1]] <<- list(
      message = message,
      type = type,
      duration = duration
    )
  }

  success_result <- finalizeMetabImportProcessingFeedback(
    status = "success",
    removeNotificationFn = function(id) {
      removed_ids <<- c(removed_ids, id)
    },
    showNotificationFn = record_notification,
    logErrorFn = function(message) {
      error_logs <<- c(error_logs, message)
    }
  )

  expect_identical(success_result$status, "success")
  expect_identical(success_result$notificationId, "metab_import_working")
  expect_identical(success_result$message, "Data imported successfully!")
  expect_identical(removed_ids, "metab_import_working")
  expect_length(shown_notifications, 1L)
  expect_identical(
    shown_notifications[[1]],
    list(
      message = "Data imported successfully!",
      type = "message",
      duration = NULL
    )
  )
  expect_identical(error_logs, character())

  removed_ids <- character()
  shown_notifications <- list()

  error_result <- finalizeMetabImportProcessingFeedback(
    status = "error",
    error = simpleError("bad import state"),
    removeNotificationFn = function(id) {
      removed_ids <<- c(removed_ids, id)
    },
    showNotificationFn = record_notification,
    logErrorFn = function(message) {
      error_logs <<- c(error_logs, message)
    }
  )

  expect_identical(error_result$status, "error")
  expect_identical(error_result$notificationId, "metab_import_working")
  expect_identical(error_result$message, "Error: bad import state")
  expect_identical(removed_ids, "metab_import_working")
  expect_length(shown_notifications, 1L)
  expect_identical(
    shown_notifications[[1]],
    list(
      message = "Error: bad import state",
      type = "error",
      duration = 10
    )
  )
  expect_identical(
    error_logs,
    "Error processing import: bad import state"
  )
})

test_that("metabolomics import assay-selection error helper preserves logging and notification contracts", {
  shown_notifications <- list()
  error_logs <- character()

  result <- finalizeMetabImportAssaySelectionError(
    error = simpleError("bad assay headers"),
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      shown_notifications[[length(shown_notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    logErrorFn = function(message) {
      error_logs <<- c(error_logs, message)
    }
  )

  expect_identical(result$status, "error")
  expect_identical(result$message, "Error importing data: bad assay headers")
  expect_identical(error_logs, "Error importing data: bad assay headers")
  expect_identical(
    shown_notifications,
    list(
      list(
        message = "Error importing data: bad assay headers",
        type = "error",
        duration = NULL
      )
    )
  )
})

test_that("metabolomics import import-data orchestration helper preserves req, delegation, and error-finalization contracts", {
  req_calls <- list()
  prepare_calls <- list()
  apply_calls <- list()
  finalize_calls <- list()
  local_data <- new.env(parent = emptyenv())
  import_state <- list(
    formatInfo = list(format = "msdial", confidence = 0.9)
  )

  result <- runMetabImportAssaySelection(
    assay1File = "assay1.csv",
    localData = local_data,
    session = "fake-session",
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    prepareImportStateFn = function(assay1File) {
      prepare_calls[[length(prepare_calls) + 1]] <<- assay1File
      import_state
    },
    applyImportStateFn = function(localData, importState, session) {
      apply_calls[[length(apply_calls) + 1]] <<- list(
        localData = localData,
        importState = importState,
        session = session
      )
      invisible(list(status = "applied", session = session))
    },
    finalizeErrorFn = function(error) {
      finalize_calls[[length(finalize_calls) + 1]] <<- error
      invisible(list(status = "error"))
    }
  )

  expect_identical(req_calls, list("assay1.csv"))
  expect_identical(prepare_calls, list("assay1.csv"))
  expect_length(apply_calls, 1L)
  expect_identical(apply_calls[[1]]$localData, local_data)
  expect_identical(apply_calls[[1]]$importState, import_state)
  expect_identical(apply_calls[[1]]$session, "fake-session")
  expect_identical(result$status, "applied")
  expect_identical(finalize_calls, list())

  error_result <- runMetabImportAssaySelection(
    assay1File = "broken.csv",
    localData = local_data,
    session = "fake-session",
    reqFn = function(value) invisible(value),
    prepareImportStateFn = function(assay1File) {
      stop(sprintf("cannot import %s", assay1File))
    },
    applyImportStateFn = function(...) {
      stop("applyImportStateFn should not run on prepare failure")
    },
    finalizeErrorFn = function(error) {
      finalize_calls[[length(finalize_calls) + 1]] <<- conditionMessage(error)
      invisible(list(status = "error", message = conditionMessage(error)))
    }
  )

  expect_identical(error_result$status, "error")
  expect_identical(error_result$message, "cannot import broken.csv")
  expect_identical(finalize_calls, list("cannot import broken.csv"))
})

test_that("metabolomics import import-data callback seam preserves current file lookup and delegation", {
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_file <- "initial.csv"
  run_calls <- list()

  import_data <- setupMetabImportAssaySelectionCallback(
    localData = local_data,
    session = "fake-session",
    runAssaySelectionFn = function(assay1File, localData, session) {
      run_calls[[length(run_calls) + 1]] <<- list(
        assay1File = assay1File,
        localData = localData,
        session = session
      )

      invisible(list(status = "delegated", assay1File = assay1File))
    }
  )

  local_data$assay1_file <- "updated.csv"
  result <- import_data()

  expect_length(run_calls, 1L)
  expect_identical(run_calls[[1]]$assay1File, "updated.csv")
  expect_identical(run_calls[[1]]$localData, local_data)
  expect_identical(run_calls[[1]]$session, "fake-session")
  expect_identical(result$status, "delegated")
  expect_identical(result$assay1File, "updated.csv")
})

test_that("metabolomics process-import observer setup seam preserves current input lookup and delegation", {
  input <- new.env(parent = emptyenv())
  input$process_import <- 0L
  input$assay1_name <- "Initial_Assay"
  input$assay2_name <- "Initial_Assay_2"
  input$vendor_format <- "msdial"
  input$sanitize_names <- FALSE
  input$is_pattern <- "^OLD_"

  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay-data"
  local_data$assay2_file <- "initial-assay2.csv"
  local_data$detected_format <- "initial-format"

  column_accessors <- list(
    getMetaboliteIdCol = function() "Peak ID",
    getAnnotationCol = function() "Annotation",
    getSampleColumns = function() c("Sample_A", "Sample_B")
  )
  workflow_data <- new.env(parent = emptyenv())
  run_calls <- list()
  observer_handler <- NULL

  setup_result <- setupMetabImportProcessingObserver(
    input = input,
    localData = local_data,
    columnAccessors = column_accessors,
    workflowData = workflow_data,
    observeEventFn = function(eventExpr, handlerExpr) {
      handler_env <- parent.frame()
      handler_expr <- substitute(handlerExpr)
      observer_handler <<- function() {
        eval(handler_expr, envir = handler_env)
      }

      invisible(NULL)
    },
    runProcessingFn = function(...) {
      run_calls[[length(run_calls) + 1]] <<- list(...)
      invisible(list(status = "delegated"))
    }
  )

  expect_null(setup_result)
  expect_true(is.function(observer_handler))

  input$assay1_name <- "Updated_Assay"
  input$assay2_name <- "Updated_Assay_2"
  input$vendor_format <- "custom"
  input$sanitize_names <- TRUE
  input$is_pattern <- "^UPDATED_"
  local_data$assay1_data <- "updated-assay-data"
  local_data$assay2_file <- "updated-assay2.csv"
  local_data$detected_format <- "updated-format"

  observer_handler()

  expect_length(run_calls, 1L)
  expect_identical(run_calls[[1]]$assay1Data, "updated-assay-data")
  expect_identical(run_calls[[1]]$assay1Name, "Updated_Assay")
  expect_identical(run_calls[[1]]$assay2File, "updated-assay2.csv")
  expect_identical(run_calls[[1]]$assay2Name, "Updated_Assay_2")
  expect_identical(run_calls[[1]]$vendorFormat, "custom")
  expect_identical(run_calls[[1]]$detectedFormat, "updated-format")
  expect_true(run_calls[[1]]$sanitizeNames)
  expect_identical(run_calls[[1]]$isPattern, "^UPDATED_")
  expect_identical(run_calls[[1]]$getMetaboliteIdColFn, column_accessors$getMetaboliteIdCol)
  expect_identical(run_calls[[1]]$getAnnotationColFn, column_accessors$getAnnotationCol)
  expect_identical(run_calls[[1]]$getSampleColumnsFn, column_accessors$getSampleColumns)
  expect_identical(run_calls[[1]]$workflowData, workflow_data)
})

test_that("metabolomics file-loaded output setup seam preserves current data lookup and unsuspended registration", {
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- NULL
  output <- new.env(parent = emptyenv())
  output_options_calls <- list()

  setup_result <- setupMetabImportFileLoadedOutput(
    output = output,
    localData = local_data,
    reactiveFn = function(expr) {
      reactive_env <- parent.frame()
      reactive_expr <- substitute(expr)

      function() {
        eval(reactive_expr, envir = reactive_env)
      }
    },
    outputOptionsFn = function(output, name, suspendWhenHidden) {
      output_options_calls[[length(output_options_calls) + 1]] <<- list(
        output = output,
        name = name,
        suspendWhenHidden = suspendWhenHidden
      )

      invisible(output)
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$file_loaded))

  local_data$assay1_data <- "updated-assay"

  expect_true(output$file_loaded())
  expect_length(output_options_calls, 1L)
  expect_identical(output_options_calls[[1]]$output, output)
  expect_identical(output_options_calls[[1]]$name, "file_loaded")
  expect_false(output_options_calls[[1]]$suspendWhenHidden)

  local_data$assay1_data <- NULL

  expect_false(output$file_loaded())
})

test_that("metabolomics import-status output setup seam preserves current workflow lookup and delegation", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "incomplete")
  workflow_data$processing_log <- list(
    setup_import = list(
      detected_format = "initial",
      n_assays = 1L,
      n_samples = 2L
    )
  )

  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportStatusOutput(
    output = output,
    workflowData = workflow_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildStatusFn = function(setupImportStatus, setupImportLog) {
      build_calls[[length(build_calls) + 1]] <<- list(
        setupImportStatus = setupImportStatus,
        setupImportLog = setupImportLog
      )

      list(
        status = setupImportStatus,
        log = setupImportLog
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$import_status))

  workflow_data$tab_status$setup_import <- "complete"
  workflow_data$processing_log$setup_import <- list(
    detected_format = "updated",
    n_assays = 2L,
    n_samples = 3L
  )

  render_result <- output$import_status()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$setupImportStatus, "complete")
  expect_identical(
    build_calls[[1]]$setupImportLog,
    workflow_data$processing_log$setup_import
  )
  expect_identical(render_result$status, "complete")
  expect_identical(render_result$log, workflow_data$processing_log$setup_import)
})

test_that("metabolomics validation-summary output setup seam preserves current data lookup and delegation", {
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay"

  column_accessors <- list(
    getMetaboliteIdCol = function() "initial-metabolite",
    getSampleColumns = function() c("initial-sample")
  )
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportValidationSummaryOutput(
    output = output,
    localData = local_data,
    columnAccessors = column_accessors,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildValidationSummaryFn = function(assayData, getMetaboliteIdColFn, getSampleColumnsFn) {
      build_calls[[length(build_calls) + 1]] <<- list(
        assayData = assayData,
        getMetaboliteIdColFn = getMetaboliteIdColFn,
        getSampleColumnsFn = getSampleColumnsFn
      )

      list(
        assayData = assayData,
        metabolite = getMetaboliteIdColFn(),
        samples = getSampleColumnsFn()
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$validation_summary))

  local_data$assay1_data <- "updated-assay"
  column_accessors$getMetaboliteIdCol <- function() "updated-metabolite"
  column_accessors$getSampleColumns <- function() c("updated-sample-a", "updated-sample-b")

  render_result <- output$validation_summary()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assayData, "updated-assay")
  expect_identical(render_result$assayData, "updated-assay")
  expect_identical(render_result$metabolite, "updated-metabolite")
  expect_identical(render_result$samples, c("updated-sample-a", "updated-sample-b"))
})

test_that("metabolomics format-detection output setup seam preserves current format lookup and delegation", {
  local_data <- new.env(parent = emptyenv())
  local_data$detected_format <- "initial-format"
  local_data$format_confidence <- 0.12
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportFormatDetectionStatusOutput(
    output = output,
    localData = local_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildFormatDetectionStatusFn = function(detectedFormat, formatConfidence) {
      build_calls[[length(build_calls) + 1]] <<- list(
        detectedFormat = detectedFormat,
        formatConfidence = formatConfidence
      )

      list(
        detectedFormat = detectedFormat,
        formatConfidence = formatConfidence
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$format_detection_status))

  local_data$detected_format <- "updated-format"
  local_data$format_confidence <- 0.87

  render_result <- output$format_detection_status()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$detectedFormat, "updated-format")
  expect_equal(build_calls[[1]]$formatConfidence, 0.87, tolerance = 1e-9)
  expect_identical(render_result$detectedFormat, "updated-format")
  expect_equal(render_result$formatConfidence, 0.87, tolerance = 1e-9)
})

test_that("metabolomics metabolite-id output setup seam preserves current data lookup and delegation", {
  input <- list(metabolite_id_col = "initial-metabolite")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay"
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportMetaboliteIdStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildMetaboliteIdStatusFn = function(assayData, metaboliteIdCol) {
      build_calls[[length(build_calls) + 1]] <<- list(
        assayData = assayData,
        metaboliteIdCol = metaboliteIdCol
      )

      list(
        assayData = assayData,
        metaboliteIdCol = metaboliteIdCol
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$metabolite_id_status))

  local_data$assay1_data <- "updated-assay"
  input$metabolite_id_col <- "updated-metabolite"

  render_result <- output$metabolite_id_status()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assayData, "updated-assay")
  expect_identical(build_calls[[1]]$metaboliteIdCol, "updated-metabolite")
  expect_identical(render_result$assayData, "updated-assay")
  expect_identical(render_result$metaboliteIdCol, "updated-metabolite")
})

test_that("metabolomics annotation-status output setup seam preserves current data lookup and delegation", {
  input <- list(annotation_col = "initial-annotation")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay"
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportAnnotationStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildAnnotationStatusFn = function(assayData, annotationCol) {
      build_calls[[length(build_calls) + 1]] <<- list(
        assayData = assayData,
        annotationCol = annotationCol
      )

      list(
        assayData = assayData,
        annotationCol = annotationCol
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$annotation_status))

  local_data$assay1_data <- "updated-assay"
  input$annotation_col <- "updated-annotation"

  render_result <- output$annotation_status()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assayData, "updated-assay")
  expect_identical(build_calls[[1]]$annotationCol, "updated-annotation")
  expect_identical(render_result$assayData, "updated-assay")
  expect_identical(render_result$annotationCol, "updated-annotation")
})

test_that("metabolomics custom-annotation output setup seam preserves current data lookup and delegation", {
  input <- list(annotation_col_custom = "initial-annotation")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay"
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportCustomAnnotationStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildCustomAnnotationStatusFn = function(assayData, columnName) {
      build_calls[[length(build_calls) + 1]] <<- list(
        assayData = assayData,
        columnName = columnName
      )

      list(
        assayData = assayData,
        columnName = columnName
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$annotation_status_custom))

  local_data$assay1_data <- "updated-assay"
  input$annotation_col_custom <- "updated-annotation"

  render_result <- output$annotation_status_custom()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assayData, "updated-assay")
  expect_identical(build_calls[[1]]$columnName, "updated-annotation")
  expect_identical(render_result$assayData, "updated-assay")
  expect_identical(render_result$columnName, "updated-annotation")
})

test_that("metabolomics custom-metabolite output setup seam preserves current data lookup and delegation", {
  input <- list(metabolite_id_col_custom = "initial-metabolite")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- "initial-assay"
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportCustomMetaboliteIdStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUIFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildCustomMetaboliteIdStatusFn = function(assayData, columnName) {
      build_calls[[length(build_calls) + 1]] <<- list(
        assayData = assayData,
        columnName = columnName
      )

      list(
        assayData = assayData,
        columnName = columnName
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$metabolite_id_status_custom))

  local_data$assay1_data <- "updated-assay"
  input$metabolite_id_col_custom <- "updated-metabolite"

  render_result <- output$metabolite_id_status_custom()

  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assayData, "updated-assay")
  expect_identical(build_calls[[1]]$columnName, "updated-metabolite")
  expect_identical(render_result$assayData, "updated-assay")
  expect_identical(render_result$columnName, "updated-metabolite")
})

test_that("metabolomics available-columns output setup seam preserves current header lookup and delegation", {
  local_data <- new.env(parent = emptyenv())
  local_data$all_headers <- c("initial_a", "initial_b")
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportAvailableColumnsDisplayOutput(
    output = output,
    localData = local_data,
    renderTextFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildAvailableColumnsDisplayFn = function(allHeaders) {
      build_calls[[length(build_calls) + 1]] <<- list(
        allHeaders = allHeaders
      )

      list(
        allHeaders = allHeaders
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$available_columns_display))

  local_data$all_headers <- c("updated_a", "updated_b", "updated_c")

  render_result <- output$available_columns_display()

  expect_length(build_calls, 1L)
  expect_identical(
    build_calls[[1]]$allHeaders,
    c("updated_a", "updated_b", "updated_c")
  )
  expect_identical(
    render_result$allHeaders,
    c("updated_a", "updated_b", "updated_c")
  )
})

test_that("metabolomics sample-columns output setup seam preserves current import-result lookup and delegation", {
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_import_result <- list(sample_columns = c("initial_a", "initial_b"))
  output <- new.env(parent = emptyenv())
  build_calls <- list()

  setup_result <- setupMetabImportSampleColumnsDisplayOutput(
    output = output,
    localData = local_data,
    renderTextFn = function(expr) {
      render_env <- parent.frame()
      render_expr <- substitute(expr)

      function() {
        eval(render_expr, envir = render_env)
      }
    },
    buildSampleColumnsDisplayFn = function(importResult) {
      build_calls[[length(build_calls) + 1]] <<- list(
        importResult = importResult
      )

      list(
        importResult = importResult
      )
    }
  )

  expect_null(setup_result)
  expect_true(is.function(output$sample_columns_display))

  local_data$assay1_import_result <- list(
    sample_columns = c("updated_a", "updated_b", "updated_c")
  )

  render_result <- output$sample_columns_display()

  expect_length(build_calls, 1L)
  expect_identical(
    build_calls[[1]]$importResult,
    list(sample_columns = c("updated_a", "updated_b", "updated_c"))
  )
  expect_identical(
    render_result$importResult,
    list(sample_columns = c("updated_a", "updated_b", "updated_c"))
  )
})

test_that("metabolomics import column-accessor setup seam preserves custom resolution and sample-column delegation", {
  input <- list(
    vendor_format = "custom",
    metabolite_id_col = "Peak ID",
    metabolite_id_col_custom = "peak id",
    annotation_col = "Annotation",
    annotation_col_custom = "annotation",
    sample_cols_pattern = "^sample"
  )
  local_data <- list(
    assay1_data = data.frame(
      "Peak ID" = c("M1", "M2"),
      Annotation = c("A", "B"),
      Sample_A = c(10, 11),
      sample_b = c(20, 21),
      check.names = FALSE
    ),
    assay1_import_result = list(sample_columns = c("Imported_A", "Imported_B"))
  )
  column_name_calls <- list()
  sample_column_calls <- list()
  req_calls <- list()

  fakeReactive <- function(expr) {
    expr_sub <- substitute(expr)
    env <- parent.frame()

    function() eval(expr_sub, env)
  }

  accessors <- setupMetabImportColumnAccessors(
    input = input,
    localData = local_data,
    reactiveFn = fakeReactive,
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    resolveColumnNameFn = function(headers, columnName) {
      column_name_calls[[length(column_name_calls) + 1]] <<- list(
        headers = headers,
        columnName = columnName
      )
      paste0("resolved:", columnName)
    },
    resolveSampleColumnsFn = function(assayData, vendorFormat, sampleColsPattern, importResult) {
      sample_column_calls[[length(sample_column_calls) + 1]] <<- list(
        assayData = assayData,
        vendorFormat = vendorFormat,
        sampleColsPattern = sampleColsPattern,
        importResult = importResult
      )
      c("Sample_A", "sample_b")
    }
  )

  expect_identical(accessors$getMetaboliteIdCol(), "resolved:peak id")
  expect_identical(accessors$getAnnotationCol(), "resolved:annotation")
  expect_identical(accessors$getSampleColumns(), c("Sample_A", "sample_b"))
  expect_identical(req_calls, list(local_data$assay1_data))
  expect_identical(
    column_name_calls,
    list(
      list(headers = names(local_data$assay1_data), columnName = "peak id"),
      list(headers = names(local_data$assay1_data), columnName = "annotation")
    )
  )
  expect_length(sample_column_calls, 1L)
  expect_identical(sample_column_calls[[1]]$assayData, local_data$assay1_data)
  expect_identical(sample_column_calls[[1]]$vendorFormat, "custom")
  expect_identical(sample_column_calls[[1]]$sampleColsPattern, "^sample")
  expect_identical(
    sample_column_calls[[1]]$importResult,
    local_data$assay1_import_result
  )
})

test_that("metabolomics import column-accessor setup seam preserves dropdown bypass and req gating", {
  input <- list(
    vendor_format = "msdial",
    metabolite_id_col = "Peak ID",
    metabolite_id_col_custom = "ignored",
    annotation_col = "Annotation",
    annotation_col_custom = "ignored",
    sample_cols_pattern = "^sample"
  )
  local_data <- list(
    assay1_data = NULL,
    assay1_import_result = list(sample_columns = c("Imported_A", "Imported_B"))
  )
  column_name_calls <- list()

  fakeReactive <- function(expr) {
    expr_sub <- substitute(expr)
    env <- parent.frame()

    function() eval(expr_sub, env)
  }

  accessors <- setupMetabImportColumnAccessors(
    input = input,
    localData = local_data,
    reactiveFn = fakeReactive,
    reqFn = function(value) {
      if (is.null(value)) {
        stop("missing assay data")
      }

      invisible(value)
    },
    resolveColumnNameFn = function(headers, columnName) {
      column_name_calls[[length(column_name_calls) + 1]] <<- list(
        headers = headers,
        columnName = columnName
      )
      columnName
    }
  )

  expect_identical(accessors$getMetaboliteIdCol(), "Peak ID")
  expect_identical(accessors$getAnnotationCol(), "Annotation")
  expect_identical(column_name_calls, list())
  expect_error(accessors$getSampleColumns(), "missing assay data")
})

test_that("metabolomics format-detection status helper preserves confidence threshold and label mapping contracts", {
  req_calls <- list()

  status_ui <- buildMetabImportFormatDetectionStatus(
    detectedFormat = "compound_discoverer",
    formatConfidence = 0.66,
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    divFn = function(...) list(...),
    strongFn = function(text) list(tag = "strong", text = text),
    brFn = function() list(tag = "br"),
    smallFn = function(text) list(tag = "small", text = text)
  )

  expect_identical(req_calls, list("compound_discoverer"))
  expect_identical(status_ui$class, "alert alert-warning")
  expect_identical(status_ui[[2]], list(tag = "strong", text = "Detected format: "))
  expect_identical(status_ui[[3]], "Compound Discoverer")
  expect_identical(status_ui[[4]], list(tag = "br"))
  expect_identical(status_ui[[5]], list(tag = "small", text = "Confidence: 66%"))

  unknown_ui <- buildMetabImportFormatDetectionStatus(
    detectedFormat = "unexpected_format",
    formatConfidence = 0.05,
    reqFn = function(value) invisible(value),
    divFn = function(...) list(...),
    strongFn = function(text) list(tag = "strong", text = text),
    brFn = function() list(tag = "br"),
    smallFn = function(text) list(tag = "small", text = text)
  )

  expect_identical(unknown_ui$class, "alert alert-danger")
  expect_identical(unknown_ui[[3]], "Unknown")
  expect_identical(unknown_ui[[5]], list(tag = "small", text = "Confidence: 5%"))
})

test_that("metabolomics format-detection status helper preserves req gating", {
  expect_error(
    buildMetabImportFormatDetectionStatus(
      detectedFormat = NULL,
      formatConfidence = 0.9,
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing detected format")
        }
        invisible(value)
      }
    ),
    "missing detected format"
  )
})

test_that("metabolomics metabolite-id status helper preserves success and missing-column contracts", {
  req_calls <- list()
  assay_data <- data.frame(
    feature_id = c("M1", "M2", "M2", NA_character_),
    intensity = c(10, 20, 30, 40),
    check.names = FALSE
  )

  success_ui <- buildMetabImportMetaboliteIdStatus(
    assayData = assay_data,
    metaboliteIdCol = "feature_id",
    reqFn = function(...) {
      req_calls[[length(req_calls) + 1]] <<- list(...)
      invisible(TRUE)
    },
    namesFn = names,
    uniqueFn = unique,
    lengthFn = length,
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style),
    sprintfFn = sprintf
  )

  expect_identical(req_calls, list(list(assay_data, "feature_id")))
  expect_identical(success_ui[[1]], list(name = "check-circle", style = "color: green;"))
  expect_identical(success_ui[[2]], " 3 unique IDs")

  missing_ui <- buildMetabImportMetaboliteIdStatus(
    assayData = assay_data,
    metaboliteIdCol = "missing_col",
    reqFn = function(...) invisible(TRUE),
    namesFn = names,
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(missing_ui[[1]], list(name = "times-circle", style = "color: red;"))
  expect_identical(missing_ui[[2]], " Column not found")
})

test_that("metabolomics metabolite-id status helper preserves req gating", {
  assay_data <- data.frame(feature_id = "M1", check.names = FALSE)

  expect_error(
    buildMetabImportMetaboliteIdStatus(
      assayData = assay_data,
      metaboliteIdCol = NULL,
      reqFn = function(...) {
        values <- list(...)
        if (any(vapply(values, is.null, logical(1)))) {
          stop("missing metabolite id state")
        }
        invisible(TRUE)
      }
    ),
    "missing metabolite id state"
  )
})

test_that("metabolomics custom metabolite-id status helper preserves prompt, resolution, and missing-column contracts", {
  req_calls <- list()
  resolve_calls <- list()
  assay_data <- data.frame(
    "Peak ID" = c("M1", "M2", "M2", NA_character_),
    intensity = c(10, 20, 30, 40),
    check.names = FALSE
  )

  prompt_ui <- buildMetabImportCustomMetaboliteIdStatus(
    assayData = assay_data,
    columnName = "",
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(req_calls, list(assay_data))
  expect_identical(prompt_ui[[1]], list(name = "question-circle", style = "color: gray;"))
  expect_identical(prompt_ui[[2]], " Enter column name")

  success_ui <- buildMetabImportCustomMetaboliteIdStatus(
    assayData = assay_data,
    columnName = "peak id",
    reqFn = function(value) invisible(value),
    namesFn = names,
    resolveColumnNameFn = function(headers, columnName) {
      resolve_calls[[length(resolve_calls) + 1]] <<- list(headers = headers, columnName = columnName)
      "Peak ID"
    },
    uniqueFn = unique,
    lengthFn = length,
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style),
    sprintfFn = sprintf
  )

  expect_identical(
    resolve_calls,
    list(list(headers = c("Peak ID", "intensity"), columnName = "peak id"))
  )
  expect_identical(success_ui[[1]], list(name = "check-circle", style = "color: green;"))
  expect_identical(success_ui[[2]], " Found: 3 unique IDs")

  missing_ui <- buildMetabImportCustomMetaboliteIdStatus(
    assayData = assay_data,
    columnName = "missing_col",
    reqFn = function(value) invisible(value),
    resolveColumnNameFn = function(headers, columnName) columnName,
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(missing_ui[[1]], list(name = "times-circle", style = "color: red;"))
  expect_identical(missing_ui[[2]], " Column not found")
})

test_that("metabolomics custom metabolite-id status helper preserves req gating", {
  expect_error(
    buildMetabImportCustomMetaboliteIdStatus(
      assayData = NULL,
      columnName = "Peak ID",
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing custom metabolite id state")
        }
        invisible(value)
      }
    ),
    "missing custom metabolite id state"
  )
})

test_that("metabolomics custom annotation status helper preserves optional, resolution, and missing-column contracts", {
  req_calls <- list()
  resolve_calls <- list()
  assay_data <- data.frame(
    Annotation = c("A", "B"),
    intensity = c(10, 20),
    check.names = FALSE
  )

  optional_ui <- buildMetabImportCustomAnnotationStatus(
    assayData = assay_data,
    columnName = "",
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(req_calls, list(assay_data))
  expect_identical(optional_ui[[1]], list(name = "minus-circle", style = "color: gray;"))
  expect_identical(optional_ui[[2]], " Optional")

  found_ui <- buildMetabImportCustomAnnotationStatus(
    assayData = assay_data,
    columnName = "annotation",
    reqFn = function(value) invisible(value),
    namesFn = names,
    resolveColumnNameFn = function(headers, columnName) {
      resolve_calls[[length(resolve_calls) + 1]] <<- list(headers = headers, columnName = columnName)
      "Annotation"
    },
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(
    resolve_calls,
    list(list(headers = c("Annotation", "intensity"), columnName = "annotation"))
  )
  expect_identical(found_ui[[1]], list(name = "check-circle", style = "color: green;"))
  expect_identical(found_ui[[2]], " Found")

  missing_ui <- buildMetabImportCustomAnnotationStatus(
    assayData = assay_data,
    columnName = "missing_col",
    reqFn = function(value) invisible(value),
    resolveColumnNameFn = function(headers, columnName) columnName,
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(missing_ui[[1]], list(name = "times-circle", style = "color: red;"))
  expect_identical(missing_ui[[2]], " Column not found")
})

test_that("metabolomics custom annotation status helper preserves req gating", {
  expect_error(
    buildMetabImportCustomAnnotationStatus(
      assayData = NULL,
      columnName = "Annotation",
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing custom annotation state")
        }
        invisible(value)
      }
    ),
    "missing custom annotation state"
  )
})

test_that("metabolomics annotation status helper preserves found, missing, and optional contracts", {
  req_calls <- list()
  assay_data <- data.frame(
    annotation = c("A", "B"),
    intensity = c(10, 20),
    check.names = FALSE
  )

  found_ui <- buildMetabImportAnnotationStatus(
    assayData = assay_data,
    annotationCol = "annotation",
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(req_calls, list(assay_data))
  expect_identical(found_ui[[1]], list(name = "check-circle", style = "color: green;"))
  expect_identical(found_ui[[2]], " OK")

  missing_ui <- buildMetabImportAnnotationStatus(
    assayData = assay_data,
    annotationCol = "missing_col",
    reqFn = function(value) invisible(value),
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(missing_ui[[1]], list(name = "times-circle", style = "color: red;"))
  expect_identical(missing_ui[[2]], " Column not found")

  optional_ui <- buildMetabImportAnnotationStatus(
    assayData = assay_data,
    annotationCol = "",
    reqFn = function(value) invisible(value),
    spanFn = function(...) list(...),
    iconFn = function(name, style = NULL) list(name = name, style = style)
  )

  expect_identical(optional_ui[[1]], list(name = "minus-circle", style = "color: gray;"))
  expect_identical(optional_ui[[2]], " Optional")
})

test_that("metabolomics annotation status helper preserves req gating", {
  expect_error(
    buildMetabImportAnnotationStatus(
      assayData = NULL,
      annotationCol = "annotation",
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing annotation status state")
        }
        invisible(value)
      }
    ),
    "missing annotation status state"
  )
})

test_that("metabolomics sample-columns display helper preserves truncation and full-list contracts", {
  long_display <- buildMetabImportSampleColumnsDisplay(
    importResult = list(sample_columns = paste0("Sample_", seq_len(12))),
    reqFn = function(value) invisible(value)
  )

  expect_identical(
    long_display,
    "Sample_1, Sample_2, Sample_3, Sample_4, Sample_5, Sample_6, Sample_7, Sample_8, Sample_9, Sample_10 ... and 2 more"
  )

  short_display <- buildMetabImportSampleColumnsDisplay(
    importResult = list(sample_columns = c("Sample_A", "Sample_B")),
    reqFn = function(value) invisible(value)
  )

  expect_identical(short_display, "Sample_A, Sample_B")
})

test_that("metabolomics sample-columns display helper preserves req gating", {
  expect_error(
    buildMetabImportSampleColumnsDisplay(
      importResult = NULL,
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing sample-columns display state")
        }
        invisible(value)
      }
    ),
    "missing sample-columns display state"
  )
})

test_that("metabolomics available-columns display helper preserves collapse contract", {
  display <- buildMetabImportAvailableColumnsDisplay(
    allHeaders = c("Peak ID", "Retention_Time", "Sample_A"),
    reqFn = function(value) invisible(value)
  )

  expect_identical(display, "Peak ID, Retention_Time, Sample_A")
})

test_that("metabolomics available-columns display helper preserves req gating", {
  expect_error(
    buildMetabImportAvailableColumnsDisplay(
      allHeaders = NULL,
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing available-columns display state")
        }
        invisible(value)
      }
    ),
    "missing available-columns display state"
  )
})

test_that("metabolomics validation summary helper preserves delegation and report assembly contracts", {
  make_div <- function(...) {
    args <- list(...)
    arg_names <- names(args)
    children <- unname(args[is.null(arg_names) | arg_names != "class"])

    list(
      tag = "div",
      class = args$class,
      children = children
    )
  }
  make_ul <- function(...) {
    children <- list(...)
    if (length(children) == 1 && is.list(children[[1]]) && length(children[[1]]) > 0) {
      children <- children[[1]]
    }

    list(tag = "ul", children = children)
  }

  assay_data <- data.frame(
    "Peak ID" = c("M1", "M2"),
    Sample_A = c(10, 11),
    Sample_B = c(20, 21),
    check.names = FALSE
  )
  req_calls <- list()
  validate_calls <- list()

  success_ui <- buildMetabImportValidationSummary(
    assayData = assay_data,
    getMetaboliteIdColFn = function() "Peak ID",
    getSampleColumnsFn = function() c("Sample_A", "Sample_B"),
    validateColumnMappingFn = function(data, metabolite_id_column, sample_columns) {
      validate_calls[[length(validate_calls) + 1]] <<- list(
        data = data,
        metabolite_id_column = metabolite_id_column,
        sample_columns = sample_columns
      )

      list(
        valid = TRUE,
        errors = character(0),
        warnings = "Duplicate metabolite IDs detected",
        summary = list(
          n_metabolites = 2L,
          n_samples = 2L,
          pct_missing = 12.5
        )
      )
    },
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    tagListFn = function(...) list(tag = "tagList", children = list(...)),
    divFn = make_div,
    iconFn = function(name, style = NULL) list(tag = "icon", name = name, style = style),
    strongFn = function(text) list(tag = "strong", text = text),
    ulFn = make_ul,
    liFn = function(text) list(tag = "li", text = text),
    lapplyFn = lapply,
    sprintfFn = sprintf
  )

  expect_identical(req_calls, list(assay_data, "Peak ID"))
  expect_identical(validate_calls[[1]]$data, assay_data)
  expect_identical(validate_calls[[1]]$metabolite_id_column, "Peak ID")
  expect_identical(validate_calls[[1]]$sample_columns, c("Sample_A", "Sample_B"))
  expect_identical(success_ui$tag, "tagList")
  expect_identical(success_ui$children[[1]]$class, "alert alert-success")
  expect_identical(success_ui$children[[1]]$children[[1]], list(tag = "icon", name = "check-circle", style = NULL))
  expect_identical(success_ui$children[[1]]$children[[2]], list(tag = "strong", text = " Validation Passed"))
  expect_identical(success_ui$children[[2]]$children[[1]], list(tag = "li", text = "Metabolites: 2"))
  expect_identical(success_ui$children[[2]]$children[[2]], list(tag = "li", text = "Samples: 2"))
  expect_identical(success_ui$children[[2]]$children[[3]], list(tag = "li", text = "Missing values: 12.5%"))
  expect_identical(success_ui$children[[3]]$class, "alert alert-warning")
  expect_identical(success_ui$children[[3]]$children[[1]], list(tag = "icon", name = "exclamation-triangle", style = NULL))
  expect_identical(success_ui$children[[3]]$children[[2]], " Warnings:")
  expect_identical(success_ui$children[[3]]$children[[3]]$children[[1]], list(tag = "li", text = "Duplicate metabolite IDs detected"))

  failure_ui <- buildMetabImportValidationSummary(
    assayData = assay_data,
    getMetaboliteIdColFn = function() "Peak ID",
    getSampleColumnsFn = function() character(0),
    validateColumnMappingFn = function(data, metabolite_id_column, sample_columns) {
      list(
        valid = FALSE,
        errors = c("No sample columns specified", "Metabolite ID column not found"),
        warnings = character(0),
        summary = list()
      )
    },
    reqFn = function(value) invisible(value),
    divFn = make_div,
    iconFn = function(name, style = NULL) list(tag = "icon", name = name, style = style),
    strongFn = function(text) list(tag = "strong", text = text),
    ulFn = make_ul,
    liFn = function(text) list(tag = "li", text = text),
    lapplyFn = lapply
  )

  expect_identical(failure_ui$class, "alert alert-danger")
  expect_identical(failure_ui$children[[1]], list(tag = "icon", name = "times-circle", style = NULL))
  expect_identical(failure_ui$children[[2]], list(tag = "strong", text = " Validation Failed"))
  expect_identical(failure_ui$children[[3]]$children[[1]], list(tag = "li", text = "No sample columns specified"))
  expect_identical(failure_ui$children[[3]]$children[[2]], list(tag = "li", text = "Metabolite ID column not found"))
})

test_that("metabolomics validation summary helper preserves req gating", {
  assay_data <- data.frame("Peak ID" = "M1", check.names = FALSE)

  expect_error(
    buildMetabImportValidationSummary(
      assayData = NULL,
      getMetaboliteIdColFn = function() "Peak ID",
      getSampleColumnsFn = function() "Sample_A",
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing validation summary state")
        }
        invisible(value)
      }
    ),
    "missing validation summary state"
  )

  expect_error(
    buildMetabImportValidationSummary(
      assayData = assay_data,
      getMetaboliteIdColFn = function() NULL,
      getSampleColumnsFn = function() "Sample_A",
      validateColumnMappingFn = function(...) {
        stop("validation should not run")
      },
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing validation summary metabolite column")
        }
        invisible(value)
      }
    ),
    "missing validation summary metabolite column"
  )
})

test_that("metabolomics import-status helper preserves completion-summary assembly and incomplete fallback", {
  make_div <- function(...) {
    args <- list(...)
    arg_names <- names(args)
    children <- unname(args[is.null(arg_names) | arg_names != "class"])

    list(
      tag = "div",
      class = args$class,
      children = children
    )
  }

  toupper_calls <- character()

  status_ui <- buildMetabImportStatus(
    setupImportStatus = "complete",
    setupImportLog = list(
      detected_format = "msdial",
      n_assays = 2L,
      n_samples = 12L
    ),
    divFn = make_div,
    iconFn = function(name, style = NULL) list(tag = "icon", name = name, style = style),
    strongFn = function(text) list(tag = "strong", text = text),
    brFn = function() list(tag = "br"),
    sprintfFn = sprintf,
    toupperFn = function(value) {
      toupper_calls <<- c(toupper_calls, value)
      toupper(value)
    }
  )

  expect_identical(toupper_calls, "msdial")
  expect_identical(status_ui$class, "alert alert-success")
  expect_identical(status_ui$children[[1]], list(tag = "icon", name = "check-circle", style = NULL))
  expect_identical(status_ui$children[[2]], list(tag = "strong", text = " Import Complete"))
  expect_identical(status_ui$children[[3]], list(tag = "br"))
  expect_identical(status_ui$children[[4]], "Format: MSDIAL | Assays: 2 | Samples: 12")

  expect_null(buildMetabImportStatus(
    setupImportStatus = NULL,
    setupImportLog = list()
  ))

  expect_null(buildMetabImportStatus(
    setupImportStatus = "pending",
    setupImportLog = list(
      detected_format = "msdial",
      n_assays = 1L,
      n_samples = 1L
    )
  ))
})

test_that("metabolomics process-import helper preserves accessor delegation and processing finalization contracts", {
  req_calls <- list()
  accessor_calls <- character()
  build_calls <- list()
  apply_calls <- list()
  finalize_calls <- list()
  shown_notifications <- list()
  info_logs <- character()
  assay_data <- data.frame(
    "Peak ID" = c("M1", "M2"),
    Sample_A = c(10, 11),
    check.names = FALSE
  )
  workflow_data <- new.env(parent = emptyenv())

  result <- runMetabImportProcessing(
    assay1Data = assay_data,
    assay1Name = "LCMS_Pos",
    assay2File = "assay2.csv",
    assay2Name = "LCMS_Neg",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    sanitizeNames = TRUE,
    isPattern = "^IS_",
    getMetaboliteIdColFn = function() {
      accessor_calls <<- c(accessor_calls, "metabolite")
      "Peak ID"
    },
    getAnnotationColFn = function() {
      accessor_calls <<- c(accessor_calls, "annotation")
      "Annotation"
    },
    getSampleColumnsFn = function() {
      accessor_calls <<- c(accessor_calls, "samples")
      c("Sample_A", "Sample_B")
    },
    workflowData = workflow_data,
    reqFn = function(value) {
      req_calls[[length(req_calls) + 1]] <<- value
      invisible(value)
    },
    showNotificationFn = function(message, type = NULL, duration = NULL, id = NULL) {
      shown_notifications[[length(shown_notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration,
        id = id
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      info_logs <<- c(info_logs, message)
      invisible(NULL)
    },
    buildWorkflowPayloadFn = function(...) {
      build_calls[[length(build_calls) + 1]] <<- list(...)
      list(
        sampleNamesSanitized = TRUE,
        sampleCols = c("Sample_A", "Sample_B"),
        payloadId = "workflow-payload"
      )
    },
    applyWorkflowPayloadFn = function(workflowData, workflowPayload) {
      apply_calls[[length(apply_calls) + 1]] <<- list(
        workflowData = workflowData,
        workflowPayload = workflowPayload
      )
      invisible(list(status = "applied"))
    },
    finalizeFeedbackFn = function(status, error = NULL) {
      finalize_calls[[length(finalize_calls) + 1]] <<- list(
        status = status,
        error = if (is.null(error)) NULL else conditionMessage(error)
      )
      invisible(list(status = status))
    }
  )

  expect_identical(req_calls, list(assay_data, "Peak ID"))
  expect_identical(accessor_calls, c("metabolite", "annotation", "samples"))
  expect_length(build_calls, 1L)
  expect_identical(build_calls[[1]]$assay1Name, "LCMS_Pos")
  expect_identical(build_calls[[1]]$assay2File, "assay2.csv")
  expect_identical(build_calls[[1]]$annotationCol, "Annotation")
  expect_identical(build_calls[[1]]$sampleCols, c("Sample_A", "Sample_B"))
  expect_true(build_calls[[1]]$sanitizeNames)
  expect_length(apply_calls, 1L)
  expect_identical(apply_calls[[1]]$workflowData, workflow_data)
  expect_identical(apply_calls[[1]]$workflowPayload$payloadId, "workflow-payload")
  expect_identical(
    shown_notifications,
    list(
      list(
        message = "Processing imported data...",
        type = NULL,
        duration = NULL,
        id = "metab_import_working"
      ),
      list(
        message = "Sample names sanitized for R compatibility.",
        type = "message",
        duration = NULL,
        id = NULL
      )
    )
  )
  expect_identical(
    info_logs,
    c(
      "Sanitizing sample names in metabolomics data...",
      "Sanitized 2 sample column names."
    )
  )
  expect_identical(finalize_calls, list(list(status = "success", error = NULL)))
  expect_identical(result$status, "success")
  expect_identical(result$finalizeResult$status, "success")

  error_result <- runMetabImportProcessing(
    assay1Data = assay_data,
    assay1Name = "LCMS_Pos",
    assay2File = NULL,
    assay2Name = "",
    vendorFormat = "msdial",
    detectedFormat = "msdial",
    sanitizeNames = FALSE,
    isPattern = NULL,
    getMetaboliteIdColFn = function() "Peak ID",
    getAnnotationColFn = function() NULL,
    getSampleColumnsFn = function() "Sample_A",
    workflowData = workflow_data,
    reqFn = function(value) invisible(value),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    buildWorkflowPayloadFn = function(...) {
      stop("cannot build workflow payload")
    },
    applyWorkflowPayloadFn = function(...) {
      stop("applyWorkflowPayloadFn should not run on payload failure")
    },
    finalizeFeedbackFn = function(status, error = NULL) {
      invisible(list(
        status = status,
        error = if (is.null(error)) NULL else conditionMessage(error)
      ))
    }
  )

  expect_identical(error_result$status, "error")
  expect_identical(error_result$error$message, "cannot build workflow payload")
  expect_identical(error_result$finalizeResult$status, "error")
  expect_identical(error_result$finalizeResult$error, "cannot build workflow payload")
})

test_that("metabolomics process-import helper preserves req gating", {
  expect_error(
    runMetabImportProcessing(
      assay1Data = NULL,
      assay1Name = "LCMS_Pos",
      assay2File = NULL,
      assay2Name = "",
      vendorFormat = "msdial",
      detectedFormat = "msdial",
      sanitizeNames = FALSE,
      isPattern = NULL,
      getMetaboliteIdColFn = function() "Peak ID",
      getAnnotationColFn = function() NULL,
      getSampleColumnsFn = function() "Sample_A",
      workflowData = list(),
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing process import assay")
        }
        invisible(value)
      }
    ),
    "missing process import assay"
  )

  assay_data <- data.frame("Peak ID" = "M1", check.names = FALSE)

  expect_error(
    runMetabImportProcessing(
      assay1Data = assay_data,
      assay1Name = "LCMS_Pos",
      assay2File = NULL,
      assay2Name = "",
      vendorFormat = "msdial",
      detectedFormat = "msdial",
      sanitizeNames = FALSE,
      isPattern = NULL,
      getMetaboliteIdColFn = function() NULL,
      getAnnotationColFn = function() NULL,
      getSampleColumnsFn = function() "Sample_A",
      workflowData = list(),
      reqFn = function(value) {
        if (is.null(value)) {
          stop("missing process import metabolite column")
        }
        invisible(value)
      }
    ),
    "missing process import metabolite column"
  )
})
