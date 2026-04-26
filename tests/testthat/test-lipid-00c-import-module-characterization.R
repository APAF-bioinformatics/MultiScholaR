# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

makeLipidImportCharacterizationOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("mod_lipid_import_ui preserves the wrapper shell contract", {
  rendered <- htmltools::renderTags(mod_lipid_import_ui("import"))$html

  expect_match(rendered, "Lipidomics Data Import", fixed = TRUE)
  expect_match(rendered, "import-process_import", fixed = TRUE)
  expect_match(rendered, "import-validation_summary", fixed = TRUE)
})

test_that("mod_lipid_import_server preserves the wrapper wiring contract", {
  workflow_data <- shiny::reactiveValues()

  testServer(
    mod_lipid_import_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = list(source_dir = tempdir()),
      volumes = c(Home = tempdir())
    ),
    {
      expect_true(exists("output"))
      expect_true(exists("session"))
      expect_true(exists("input"))
      expect_s3_class(session, "session_proxy")
    }
  )
})

test_that("buildLipidImportAssayInputPanel preserves shinyFiles and standard-upload wiring", {
  skip_if_not(
    exists("buildLipidImportAssayInputPanel", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_assay_input_panel <- get("buildLipidImportAssayInputPanel", envir = asNamespace("MultiScholaR"))
  shiny_files_html <- htmltools::renderTags(
    build_assay_input_panel(
      ns = function(id) paste0("lipid-", id),
      assayNameInputId = "assay1_name",
      assayNameLabel = "Assay Name",
      assayNameValue = "LCMS_Pos",
      shinyFilesButtonId = "assay1_file",
      standardInputId = "assay1_file_std",
      pathOutputId = "assay1_path",
      useShinyFiles = TRUE,
      buildShinyFilesButton = function(...) shiny::tags$span("shiny-files-button"),
      buildPathOutput = function(...) shiny::tags$span("path-output"),
      buildFileInput = function(...) shiny::tags$span("standard-file-input")
    )
  )$html

  expect_match(shiny_files_html, "shiny-files-button", fixed = TRUE)
  expect_match(shiny_files_html, "path-output", fixed = TRUE)

  standard_upload_html <- htmltools::renderTags(
    build_assay_input_panel(
      ns = function(id) paste0("lipid-", id),
      assayNameInputId = "assay2_name",
      assayNameLabel = "Assay Name (Optional)",
      assayNameValue = "",
      assayNamePlaceholder = "e.g., LCMS_Neg",
      shinyFilesButtonId = "assay2_file",
      standardInputId = "assay2_file_std",
      pathOutputId = "assay2_path",
      useShinyFiles = FALSE,
      buildShinyFilesButton = function(...) shiny::tags$span("shiny-files-button"),
      buildPathOutput = function(...) shiny::tags$span("path-output"),
      buildFileInput = function(...) shiny::tags$span("standard-file-input")
    )
  )$html

  expect_match(standard_upload_html, "standard-file-input", fixed = TRUE)
})

test_that("buildLipidImportFileImportSection and footer preserve UI wiring", {
  skip_if_not(
    exists("buildLipidImportFileImportSection", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportProcessFooterSection", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_file_import_section <- get("buildLipidImportFileImportSection", envir = asNamespace("MultiScholaR"))
  build_footer_section <- get("buildLipidImportProcessFooterSection", envir = asNamespace("MultiScholaR"))
  import_html <- htmltools::renderTags(
    build_file_import_section(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = TRUE,
      buildAssayInputPanel = function(...) shiny::tags$span("assay-panel"),
      buildRadioButtons = function(...) shiny::tags$span("vendor-radio")
    )
  )$html

  expect_match(import_html, "vendor-radio", fixed = TRUE)
  expect_match(import_html, "assay-panel", fixed = TRUE)
  expect_match(import_html, "Step 1: Select Vendor Format", fixed = TRUE)
  expect_match(import_html, "Step 2: Import Data Files", fixed = TRUE)

  footer_html <- htmltools::renderTags(
    build_footer_section(
      ns = function(id) paste0("lipid-", id),
      buildActionButton = function(...) shiny::tags$span("process-button"),
      buildUiOutput = function(...) shiny::tags$span("status-output")
    )
  )$html

  expect_match(footer_html, "process-button", fixed = TRUE)
  expect_match(footer_html, "status-output", fixed = TRUE)
})

test_that("buildLipidImportColumnMappingSection and module shell preserve layout wiring", {
  skip_if_not(
    exists("buildLipidImportColumnMappingSection", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportModulePanel", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportUiShell", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_column_mapping_section <- get("buildLipidImportColumnMappingSection", envir = asNamespace("MultiScholaR"))
  build_module_panel <- get("buildLipidImportModulePanel", envir = asNamespace("MultiScholaR"))
  build_ui_shell <- get("buildLipidImportUiShell", envir = asNamespace("MultiScholaR"))
  captured <- list()

  mapping_html <- htmltools::renderTags(
    build_column_mapping_section(function(id) paste0("lipid-", id))
  )$html

  expect_match(mapping_html, "Step 3: Column Mapping", fixed = TRUE)
  expect_match(mapping_html, "lipid-sample_columns_display", fixed = TRUE)
  expect_match(mapping_html, "lipid-validation_summary", fixed = TRUE)

  panel_html <- htmltools::renderTags(
    build_module_panel(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = FALSE,
      buildFileImportSection = function(...) {
        captured$file_section <- list(...)
        shiny::tags$span("file-section")
      },
      buildColumnMappingSection = function(...) {
        captured$mapping_section <- list(...)
        shiny::tags$span("mapping-section")
      },
      buildProcessFooterSection = function(...) {
        captured$footer <- list(...)
        shiny::tags$span("footer-section")
      }
    )
  )$html

  expect_match(panel_html, "Lipidomics Data Import", fixed = TRUE)
  expect_match(panel_html, "file-section", fixed = TRUE)
  expect_match(panel_html, "mapping-section", fixed = TRUE)
  expect_match(panel_html, "footer-section", fixed = TRUE)

  shell_html <- htmltools::renderTags(
    build_ui_shell(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = FALSE,
      buildModulePanel = function(...) {
        captured$module <- list(...)
        shiny::tags$span("module-panel")
      },
      useShinyjs = function() shiny::tags$span("shinyjs-stub")
    )
  )$html

  expect_match(shell_html, "shinyjs-stub", fixed = TRUE)
  expect_match(shell_html, "module-panel", fixed = TRUE)
})

test_that("lipid import format and column helpers preserve detection and fallback behavior", {
  skip_if_not(
    exists("buildLipidImportFormatDetectionStatus", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportColumnValidationStatus", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("resolveLipidImportEffectiveColumn", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("resolveLipidImportSampleColumns", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("formatLipidImportColumnPreviewText", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  build_format_status <- get("buildLipidImportFormatDetectionStatus", envir = asNamespace("MultiScholaR"))
  build_validation_status <- get("buildLipidImportColumnValidationStatus", envir = asNamespace("MultiScholaR"))
  resolve_effective_column <- get("resolveLipidImportEffectiveColumn", envir = asNamespace("MultiScholaR"))
  resolve_sample_columns <- get("resolveLipidImportSampleColumns", envir = asNamespace("MultiScholaR"))
  format_preview <- get("formatLipidImportColumnPreviewText", envir = asNamespace("MultiScholaR"))

  success_html <- htmltools::renderTags(build_format_status("msdial", 0.72))$html
  fallback_html <- htmltools::renderTags(build_format_status("unexpected_vendor", 0.12))$html
  expect_match(success_html, "MS-DIAL", fixed = TRUE)
  expect_match(fallback_html, "Unknown", fixed = TRUE)

  assay_data <- data.frame(
    `Lipid Name` = c("PC 34:1", "TG 52:2", "TG 52:2"),
    Annotation = c("class_a", "class_b", "class_b"),
    Sample_A_Norm = c(10, 11, 12),
    Sample_B = c(20, 21, 22),
    check.names = FALSE
  )

  dropdown_html <- htmltools::renderTags(
    build_validation_status(
      assayData = assay_data,
      columnName = "Lipid Name",
      successMode = "unique_ids"
    )
  )$html
  optional_html <- htmltools::renderTags(
    build_validation_status(
      assayData = assay_data,
      columnName = "",
      successMode = "ok",
      emptyMode = "optional"
    )
  )$html

  expect_match(dropdown_html, "2 unique IDs", fixed = TRUE)
  expect_match(optional_html, "Optional", fixed = TRUE)

  expect_identical(
    resolve_effective_column(
      assayData = assay_data,
      vendorFormat = "custom",
      selectedColumn = "unused",
      customColumn = "lipid name"
    ),
    "Lipid Name"
  )

  expect_identical(
    resolve_sample_columns(
      assayData = assay_data,
      assayImportResult = list(sample_columns = c("Sample_A_Norm", "Sample_B")),
      vendorFormat = "msdial",
      sampleColsPattern = "",
      excludeNormalized = TRUE
    ),
    "Sample_B"
  )

  expect_identical(
    format_preview(paste0("Sample_", seq_len(12)), maxColumns = 10),
    "Sample_1, Sample_2, Sample_3, Sample_4, Sample_5, Sample_6, Sample_7, Sample_8, Sample_9, Sample_10 ... and 2 more"
  )
})

test_that("lipid import workflow helpers preserve workflow commit and setup finalization", {
  skip_if_not(
    exists("applyLipidImportResultToWorkflow", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("finalizeLipidImportSetupState", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  apply_import_result <- get("applyLipidImportResultToWorkflow", envir = asNamespace("MultiScholaR"))
  finalize_setup <- get("finalizeLipidImportSetupState", envir = asNamespace("MultiScholaR"))
  workflow_data <- new.env(parent = emptyenv())
  workflow_type_calls <- character()
  log_messages <- character()

  workflow_data$state_manager <- list(
    setWorkflowType = function(workflow_type) {
      workflow_type_calls <<- c(workflow_type_calls, workflow_type)
    }
  )
  assay_list <- list(
    LCMS_Pos = data.frame(
      lipid_id = c("L1", "L2", "L2"),
      sample_a = c(10, 11, 12)
    ),
    LCMS_Neg = data.frame(
      lipid_id = c("L3", "L4"),
      sample_b = c(20, 21)
    )
  )

  result <- apply_import_result(
    workflowData = workflow_data,
    assayList = assay_list,
    dataFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = c("sample_a", "sample_b"),
    isPattern = "^ISTD_",
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(result, "lipidomics_standard")
  expect_identical(workflow_type_calls, "lipidomics_standard")
  expect_identical(workflow_data$data_format, "msdial")
  expect_identical(workflow_data$data_type, "lipid")
  expect_identical(workflow_data$column_mapping$annotation_col, "annotation")

  workflow_data$processing_log <- list()
  workflow_data$tab_status <- list(other_tab = "incomplete")
  log_messages <- character()
  timestamp <- as.POSIXct("2026-04-15 11:20:00", tz = "UTC")

  finalize_setup(
    workflowData = workflow_data,
    assayList = assay_list,
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    sampleColumns = c("sample_a", "sample_b"),
    now = function() timestamp,
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(workflow_data$processing_log$setup_import$timestamp, timestamp)
  expect_identical(workflow_data$processing_log$setup_import$n_assays, 2L)
  expect_identical(as.integer(workflow_data$processing_log$setup_import$n_lipids), c(2L, 2L))
  expect_identical(workflow_data$tab_status$setup_import, "complete")
  expect_identical(log_messages, "Lipidomics import complete: 2 assays, 5 total lipids")
})

test_that("lipid import sample-name sanitation and assay assembly preserve helper behavior", {
  skip_if_not(
    exists("sanitizeLipidImportSampleNames", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("assembleLipidImportAssayList", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  sanitize_sample_names <- get("sanitizeLipidImportSampleNames", envir = asNamespace("MultiScholaR"))
  assemble_assay_list <- get("assembleLipidImportAssayList", envir = asNamespace("MultiScholaR"))

  assay_list <- list(
    LCMS_Pos = data.frame(
      lipid_id = c("L1", "L2"),
      `123-Sample!` = c(10, 11),
      `Control (A)` = c(20, 21),
      check.names = FALSE
    )
  )
  log_messages <- character()
  notifications <- list()

  sanitized <- sanitize_sample_names(
    assayList = assay_list,
    sampleColumns = c("123-Sample!", "Control (A)"),
    sanitizeNames = TRUE,
    makeCleanNames = function(x) c("x123_sample", "control_a"),
    logInfo = function(message) log_messages <<- c(log_messages, message),
    notify = function(message, type) {
      notifications <<- append(notifications, list(list(message = message, type = type)))
      invisible(NULL)
    }
  )

  expect_identical(sanitized$sampleColumns, c("x123_sample", "control_a"))
  expect_identical(names(sanitized$assayList$LCMS_Pos), c("lipid_id", "x123_sample", "control_a"))
  expect_identical(length(notifications), 1L)

  assay2_data <- data.frame(
    lipid_id = c("L3", "L4"),
    sample_b = c(20, 21)
  )
  assembled <- assemble_assay_list(
    assay1Name = "LCMS_Pos",
    assay1Data = assay_list$LCMS_Pos,
    assay2File = "/tmp/lcms-neg.tsv",
    assay2Name = "LCMS_Neg",
    importSecondAssay = function(path) {
      expect_identical(path, "/tmp/lcms-neg.tsv")
      list(data = assay2_data)
    }
  )

  expect_identical(names(assembled), c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(assembled$LCMS_Neg, assay2_data)
})

test_that("runLipidImportProcessing preserves success and failure cleanup behavior", {
  skip_if_not(
    exists("runLipidImportProcessing", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_import_processing <- get("runLipidImportProcessing", envir = asNamespace("MultiScholaR"))
  workflow_data <- new.env(parent = emptyenv())
  notifications <- list()
  removed_notifications <- character()
  apply_args <- NULL
  finalize_args <- NULL

  success_result <- run_import_processing(
    workflowData = workflow_data,
    assay1Name = "LCMS_Pos",
    assay1Data = data.frame(lipid_id = c("L1", "L2"), sample_a = c(10, 11)),
    assay2File = NULL,
    assay2Name = "",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = "sample_a",
    isPattern = "^ISTD_",
    sanitizeNames = TRUE,
    assembleAssayList = function(...) list(LCMS_Pos = data.frame(lipid_id = c("L1", "L2"), sample_a = c(10, 11))),
    sanitizeSampleNames = function(...) {
      list(
        assayList = list(LCMS_Pos = data.frame(lipid_id = c("L1", "L2"), clean_sample_a = c(10, 11))),
        sampleColumns = "clean_sample_a"
      )
    },
    applyResultToWorkflow = function(...) {
      apply_args <<- list(...)
      invisible("lipidomics_standard")
    },
    finalizeSetupState = function(...) {
      finalize_args <<- list(...)
      invisible(NULL)
    },
    notify = function(message, type = NULL, id = NULL, duration = NULL) {
      notifications <<- append(notifications, list(list(message = message, type = type, id = id, duration = duration)))
      invisible(NULL)
    },
    removeNotify = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    },
    logError = function(...) stop("unexpected error path")
  )

  expect_identical(success_result$sampleColumns, "clean_sample_a")
  expect_identical(removed_notifications, "lipid_import_working")
  expect_length(notifications, 2L)
  expect_true(!is.null(apply_args))
  expect_true(!is.null(finalize_args))

  notifications <- list()
  removed_notifications <- character()
  error_messages <- character()
  failure_result <- run_import_processing(
    workflowData = workflow_data,
    assay1Name = "LCMS_Pos",
    assay1Data = data.frame(lipid_id = "L1"),
    assay2File = NULL,
    assay2Name = "",
    vendorFormat = "msdial",
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "",
    sampleColumns = "sample_a",
    isPattern = "",
    sanitizeNames = FALSE,
    assembleAssayList = function(...) stop("assembly exploded"),
    sanitizeSampleNames = function(...) stop("should not sanitize"),
    applyResultToWorkflow = function(...) stop("should not apply"),
    finalizeSetupState = function(...) stop("should not finalize"),
    notify = function(message, type = NULL, id = NULL, duration = NULL) {
      notifications <<- append(notifications, list(list(message = message, type = type, id = id, duration = duration)))
      invisible(NULL)
    },
    removeNotify = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    },
    logError = function(message) {
      error_messages <<- c(error_messages, message)
      invisible(NULL)
    }
  )

  expect_null(failure_result)
  expect_identical(removed_notifications, "lipid_import_working")
  expect_match(error_messages[[1]], "assembly exploded", fixed = TRUE)
  expect_identical(notifications[[2]]$type, "error")
})

test_that("lipid import preview and status helpers preserve format and display branches", {
  skip_if_not(
    exists("loadLipidImportAssayPreview", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportValidationSummary", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildLipidImportStatusDisplay", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  load_assay_preview <- get("loadLipidImportAssayPreview", envir = asNamespace("MultiScholaR"))
  build_validation_summary <- get("buildLipidImportValidationSummary", envir = asNamespace("MultiScholaR"))
  build_status_display <- get("buildLipidImportStatusDisplay", envir = asNamespace("MultiScholaR"))
  log_messages <- character()

  msdial_preview <- load_assay_preview(
    assay1File = "assay1.tsv",
    readHeaders = function(path) {
      expect_identical(path, "assay1.tsv")
      c("LipidName", "Annotation", "Sample_A")
    },
    detectFormat = function(headers, filename) {
      expect_identical(headers, c("LipidName", "Annotation", "Sample_A"))
      expect_identical(filename, "assay1.tsv")
      list(format = "msdial", confidence = 0.72)
    },
    importMsdial = function(path) {
      expect_identical(path, "assay1.tsv")
      list(
        data = data.frame(LipidName = "L1", Annotation = "A", Sample_A = 10),
        detected_columns = list(lipid_id = "LipidName", annotation = "Annotation"),
        is_pattern = "^ISTD_"
      )
    },
    importLipidSearch = function(...) stop("unexpected lipidsearch import"),
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  lipidsearch_preview <- load_assay_preview(
    assay1File = "assay2.tsv",
    readHeaders = function(path) {
      expect_identical(path, "assay2.tsv")
      c("Lipid", "Sample_B")
    },
    detectFormat = function(headers, filename) {
      expect_identical(headers, c("Lipid", "Sample_B"))
      expect_identical(filename, "assay2.tsv")
      list(format = "lipidsearch", confidence = 0.91)
    },
    importMsdial = function(...) stop("unexpected msdial import"),
    importLipidSearch = function(path) {
      expect_identical(path, "assay2.tsv")
      list(
        data = data.frame(Lipid = "L2", Sample_B = 20),
        detected_columns = list(lipid_id = "Lipid", annotation = ""),
        is_pattern = NA_character_
      )
    },
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(msdial_preview$detectedFormat, "msdial")
  expect_identical(msdial_preview$updates$lipidId$selected, "LipidName")
  expect_identical(msdial_preview$updates$isPattern, "^ISTD_")
  expect_identical(lipidsearch_preview$detectedFormat, "lipidsearch")
  expect_null(lipidsearch_preview$updates$isPattern)
  expect_length(log_messages, 2L)

  expect_error(
    load_assay_preview(
      assay1File = "broken.tsv",
      readHeaders = function(path) character(0),
      detectFormat = function(...) stop("should not detect format"),
      importMsdial = function(...) stop("should not import"),
      importLipidSearch = function(...) stop("should not import"),
      logInfo = function(...) invisible(NULL)
    ),
    "Could not read headers from file"
  )

  success_html <- htmltools::renderTags(
    build_validation_summary(
      list(
        valid = TRUE,
        summary = list(n_lipids = 2L, n_samples = 3L, pct_missing = 12.5),
        warnings = c("Check duplicate IDs"),
        errors = character()
      )
    )
  )$html
  failure_html <- htmltools::renderTags(
    build_validation_summary(
      list(
        valid = FALSE,
        summary = list(),
        warnings = character(),
        errors = c("Missing Lipid ID column")
      )
    )
  )$html

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = NULL)
  expect_null(build_status_display(workflow_data))

  workflow_data$tab_status$setup_import <- "complete"
  workflow_data$processing_log <- list(
    setup_import = list(detected_format = "msdial", n_assays = 2L, n_samples = 4L)
  )
  status_html <- htmltools::renderTags(build_status_display(workflow_data))$html

  expect_match(success_html, "Validation Passed", fixed = TRUE)
  expect_match(success_html, "Warnings:", fixed = TRUE)
  expect_match(failure_html, "Validation Failed", fixed = TRUE)
  expect_match(failure_html, "Missing Lipid ID column", fixed = TRUE)
  expect_match(status_html, "Import Complete", fixed = TRUE)
  expect_match(status_html, "MSDIAL | Assays: 2 | Samples: 4", fixed = TRUE)
})

test_that("lipid import selection and process helpers preserve event handling behavior", {
  skip_if_not(
    exists("handleLipidImportFileSelection", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportProcessRequest", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  handle_file_selection <- get("handleLipidImportFileSelection", envir = asNamespace("MultiScholaR"))
  handle_process_request <- get("handleLipidImportProcessRequest", envir = asNamespace("MultiScholaR"))
  selected_paths <- character()
  error_messages <- character()

  expect_false(
    handle_file_selection(
      fileInput = NULL,
      volumes = c(Home = tempdir()),
      onPathSelected = function(path) selected_paths <<- c(selected_paths, path)
    )
  )
  expect_false(
    handle_file_selection(
      fileInput = 0L,
      volumes = c(Home = tempdir()),
      onPathSelected = function(path) selected_paths <<- c(selected_paths, path)
    )
  )
  expect_true(
    handle_file_selection(
      fileInput = list(path = "/tmp/selected.tsv"),
      volumes = c(Home = tempdir()),
      onPathSelected = function(path) selected_paths <<- c(selected_paths, path),
      parseFilePaths = function(volumes, fileInput) data.frame(datapath = "/tmp/selected.tsv"),
      logError = function(message) error_messages <<- c(error_messages, message)
    )
  )
  expect_identical(selected_paths, "/tmp/selected.tsv")
  expect_false(
    handle_file_selection(
      fileInput = list(path = "/tmp/empty.tsv"),
      volumes = c(Home = tempdir()),
      onPathSelected = function(path) selected_paths <<- c(selected_paths, path),
      parseFilePaths = function(volumes, fileInput) data.frame(datapath = character()),
      logError = function(message) error_messages <<- c(error_messages, message)
    )
  )
  expect_false(
    handle_file_selection(
      fileInput = list(path = "/tmp/error.tsv"),
      volumes = c(Home = tempdir()),
      onPathSelected = function(path) selected_paths <<- c(selected_paths, path),
      parseFilePaths = function(...) stop("parse exploded"),
      logError = function(message) error_messages <<- c(error_messages, message)
    )
  )
  expect_match(error_messages[[1]], "Error parsing file path: parse exploded", fixed = TRUE)

  process_args <- NULL
  handle_process_request(
    workflowData = list(stage = "import"),
    assay1Name = "LCMS_Pos",
    assay1Data = data.frame(lipid_id = "L1", sample_a = 10),
    assay2File = NULL,
    assay2Name = "",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = "sample_a",
    isPattern = "^ISTD_",
    sanitizeNames = TRUE,
    processImport = function(...) {
      process_args <<- list(...)
      invisible("processed")
    }
  )

  expect_identical(process_args$assay1Name, "LCMS_Pos")
  expect_identical(process_args$vendorFormat, "custom")
  expect_identical(process_args$lipidIdCol, "lipid_id")
  expect_identical(process_args$sanitizeNames, TRUE)
})

test_that("lipid import output render helpers preserve forwarded payloads", {
  skip_if_not(
    exists("handleLipidImportValidationSummaryRender", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportSampleColumnsDisplayRender", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportAvailableColumnsDisplayRender", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportCustomLipidIdStatusRender", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportCustomAnnotationStatusRender", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleLipidImportStatusRender", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  handle_validation_summary <- get("handleLipidImportValidationSummaryRender", envir = asNamespace("MultiScholaR"))
  handle_sample_columns <- get("handleLipidImportSampleColumnsDisplayRender", envir = asNamespace("MultiScholaR"))
  handle_available_columns <- get("handleLipidImportAvailableColumnsDisplayRender", envir = asNamespace("MultiScholaR"))
  handle_custom_lipid_id <- get("handleLipidImportCustomLipidIdStatusRender", envir = asNamespace("MultiScholaR"))
  handle_custom_annotation <- get("handleLipidImportCustomAnnotationStatusRender", envir = asNamespace("MultiScholaR"))
  handle_status <- get("handleLipidImportStatusRender", envir = asNamespace("MultiScholaR"))

  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    sample_a = c(10, 11),
    sample_b = c(20, 21)
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "complete")
  workflow_data$processing_log <- list(setup_import = list(detected_format = "msdial", n_assays = 1L, n_samples = 2L))

  expect_identical(
    handle_validation_summary(
      assay1Data = assay_data,
      lipidIdCol = "lipid_id",
      sampleColumns = c("sample_a", "sample_b"),
      validateMapping = function(data, lipid_id_column, sample_columns) {
        list(data = data, lipid_id_column = lipid_id_column, sample_columns = sample_columns)
      },
      buildSummary = function(validation) list(kind = "summary", validation = validation)
    ),
    list(kind = "summary", validation = list(data = assay_data, lipid_id_column = "lipid_id", sample_columns = c("sample_a", "sample_b")))
  )
  expect_identical(
    handle_sample_columns(
      assay1ImportResult = list(sample_columns = c("sample_a", "sample_b")),
      formatPreviewText = function(columnNames, maxColumns = Inf) paste(columnNames, collapse = ", ")
    ),
    "sample_a, sample_b"
  )
  expect_identical(
    handle_available_columns(
      allHeaders = c("lipid_id", "annotation", "sample_a", "sample_b"),
      formatPreviewText = function(columnNames, maxColumns = Inf) paste(columnNames, collapse = ", ")
    ),
    "lipid_id, annotation, sample_a, sample_b"
  )
  expect_identical(
    handle_custom_lipid_id(
      assay1Data = assay_data,
      lipidIdColCustom = "lipid_id",
      buildValidationStatus = function(...) "custom-lipid-id-status"
    ),
    "custom-lipid-id-status"
  )
  expect_identical(
    handle_custom_annotation(
      assay1Data = assay_data,
      annotationColCustom = "annotation",
      buildValidationStatus = function(...) "custom-annotation-status"
    ),
    "custom-annotation-status"
  )
  expect_identical(
    handle_status(
      workflowData = workflow_data,
      buildStatusDisplay = function(workflowData) list(kind = "status", status = workflowData$tab_status$setup_import)
    ),
    list(kind = "status", status = "complete")
  )
})
