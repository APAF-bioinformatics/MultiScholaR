library(testthat)

source(test_path("..", "..", "R", "mod_lipid_import_ui_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_import_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_import_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_import_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_import.R"), local = environment())

test_that("buildLipidImportFormatDetectionStatus renders the vendor label and success class", {
  rendered <- htmltools::renderTags(
    buildLipidImportFormatDetectionStatus("msdial", 0.72)
  )$html

  expect_match(rendered, "alert-success", fixed = TRUE)
  expect_match(rendered, "MS-DIAL", fixed = TRUE)
  expect_match(rendered, "Confidence: 72%", fixed = TRUE)
})

test_that("buildLipidImportFormatDetectionStatus keeps warning and fallback mappings stable", {
  warning_html <- htmltools::renderTags(
    buildLipidImportFormatDetectionStatus("lipidsearch", 0.41)
  )$html
  fallback_html <- htmltools::renderTags(
    buildLipidImportFormatDetectionStatus("unexpected_vendor", 0.12)
  )$html

  expect_match(warning_html, "alert-warning", fixed = TRUE)
  expect_match(warning_html, "LipidSearch", fixed = TRUE)
  expect_match(warning_html, "Confidence: 41%", fixed = TRUE)

  expect_match(fallback_html, "alert-danger", fixed = TRUE)
  expect_match(fallback_html, "Unknown", fixed = TRUE)
  expect_match(fallback_html, "Confidence: 12%", fixed = TRUE)
})

test_that("buildLipidImportColumnValidationStatus preserves dropdown and custom ID messaging", {
  assay_data <- data.frame(
    `Lipid Name` = c("PC 34:1", "TG 52:2", "TG 52:2"),
    Annotation = c("class_a", "class_b", "class_b"),
    check.names = FALSE
  )

  dropdown_html <- htmltools::renderTags(
    buildLipidImportColumnValidationStatus(
      assayData = assay_data,
      columnName = "Lipid Name",
      successMode = "unique_ids"
    )
  )$html
  custom_html <- htmltools::renderTags(
    buildLipidImportColumnValidationStatus(
      assayData = assay_data,
      columnName = "lipid name",
      successMode = "found_unique_ids",
      emptyMode = "prompt",
      allowCaseInsensitive = TRUE
    )
  )$html

  expect_match(dropdown_html, "circle-check icon", fixed = TRUE)
  expect_match(dropdown_html, "2 unique IDs", fixed = TRUE)
  expect_match(custom_html, "circle-check icon", fixed = TRUE)
  expect_match(custom_html, "Found: 2 unique IDs", fixed = TRUE)
})

test_that("buildLipidImportColumnValidationStatus keeps optional and missing-column states stable", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )

  optional_html <- htmltools::renderTags(
    buildLipidImportColumnValidationStatus(
      assayData = assay_data,
      columnName = "",
      successMode = "ok",
      emptyMode = "optional"
    )
  )$html
  missing_html <- htmltools::renderTags(
    buildLipidImportColumnValidationStatus(
      assayData = assay_data,
      columnName = "missing_annotation",
      successMode = "found",
      emptyMode = "optional",
      allowCaseInsensitive = TRUE
    )
  )$html

  expect_match(optional_html, "circle-minus icon", fixed = TRUE)
  expect_match(optional_html, "Optional", fixed = TRUE)
  expect_match(missing_html, "circle-xmark icon", fixed = TRUE)
  expect_match(missing_html, "Column not found", fixed = TRUE)
})

test_that("resolveLipidImportEffectiveColumn preserves dropdown passthrough and custom lookup", {
  assay_data <- data.frame(
    `Lipid Name` = c("PC 34:1", "TG 52:2"),
    Annotation = c("class_a", "class_b"),
    check.names = FALSE
  )

  expect_identical(
    resolveLipidImportEffectiveColumn(
      assayData = assay_data,
      vendorFormat = "msdial",
      selectedColumn = "Peak ID",
      customColumn = "ignored"
    ),
    "Peak ID"
  )

  expect_identical(
    resolveLipidImportEffectiveColumn(
      assayData = assay_data,
      vendorFormat = "custom",
      selectedColumn = "unused",
      customColumn = "lipid name"
    ),
    "Lipid Name"
  )
})

test_that("resolveLipidImportEffectiveColumn keeps unmatched and empty custom values stable", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )

  expect_identical(
    resolveLipidImportEffectiveColumn(
      assayData = assay_data,
      vendorFormat = "custom",
      selectedColumn = "unused",
      customColumn = "missing_annotation"
    ),
    "missing_annotation"
  )

  expect_identical(
    resolveLipidImportEffectiveColumn(
      assayData = assay_data,
      vendorFormat = "custom",
      selectedColumn = "unused",
      customColumn = ""
    ),
    ""
  )
})

test_that("resolveLipidImportSampleColumns preserves custom pattern matches without exclusion filtering", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    Sample_A_Norm = c(10, 11),
    Sample_B = c(20, 21),
    Meta = c("m1", "m2"),
    check.names = FALSE
  )

  expect_identical(
    resolveLipidImportSampleColumns(
      assayData = assay_data,
      assayImportResult = list(sample_columns = c("detected_norm", "detected_raw")),
      vendorFormat = "custom",
      sampleColsPattern = "^sample_",
      excludeNormalized = TRUE
    ),
    c("Sample_A_Norm", "Sample_B")
  )
})

test_that("resolveLipidImportSampleColumns preserves detected and fallback exclusion behavior", {
  detected_assay <- data.frame(
    lipid_id = c("L1", "L2"),
    Sample_A = c(10, 11),
    Sample_B_Normalized = c(20, 21),
    check.names = FALSE
  )
  fallback_assay <- data.frame(
    lipid_id = c("L3", "L4"),
    Sample_C = c(30, 31),
    Sample_D_normallized = c(40, 41),
    Batch = c("b1", "b1"),
    check.names = FALSE
  )

  expect_identical(
    resolveLipidImportSampleColumns(
      assayData = detected_assay,
      assayImportResult = list(sample_columns = c("Sample_A", "Sample_B_Normalized")),
      vendorFormat = "msdial",
      sampleColsPattern = "",
      excludeNormalized = TRUE
    ),
    "Sample_A"
  )

  expect_identical(
    resolveLipidImportSampleColumns(
      assayData = fallback_assay,
      assayImportResult = NULL,
      vendorFormat = "custom",
      sampleColsPattern = "^missing_",
      excludeNormalized = TRUE
    ),
    "Sample_C"
  )
})

test_that("formatLipidImportColumnPreviewText preserves sample preview truncation after ten columns", {
  columns <- paste0("Sample_", seq_len(12))

  expect_identical(
    formatLipidImportColumnPreviewText(columns, maxColumns = 10),
    paste(
      paste(columns[seq_len(10)], collapse = ", "),
      "... and 2 more"
    )
  )
})

test_that("formatLipidImportColumnPreviewText preserves full and empty preview text", {
  expect_identical(
    formatLipidImportColumnPreviewText(c("lipid_id", "annotation", "sample_a")),
    "lipid_id, annotation, sample_a"
  )

  expect_identical(
    formatLipidImportColumnPreviewText(character(0)),
    ""
  )
})

test_that("buildLipidImportAssayInputPanel preserves shinyFiles assay-panel wiring", {
  captured <- list()

  rendered <- htmltools::renderTags(
    buildLipidImportAssayInputPanel(
      ns = function(id) paste0("lipid-", id),
      assayNameInputId = "assay1_name",
      assayNameLabel = "Assay Name",
      assayNameValue = "LCMS_Pos",
      shinyFilesButtonId = "assay1_file",
      standardInputId = "assay1_file_std",
      pathOutputId = "assay1_path",
      useShinyFiles = TRUE,
      buildShinyFilesButton = function(...) {
        captured$button <<- list(...)
        shiny::tags$span("shiny-files-button")
      },
      buildPathOutput = function(...) {
        captured$path <<- list(...)
        shiny::tags$span("path-output")
      },
      buildFileInput = function(...) {
        captured$fileInput <<- list(...)
        shiny::tags$span("standard-file-input")
      }
    )
  )$html

  expect_match(rendered, "Assay Name", fixed = TRUE)
  expect_match(rendered, "LCMS_Pos", fixed = TRUE)
  expect_match(rendered, "shiny-files-button", fixed = TRUE)
  expect_match(rendered, "path-output", fixed = TRUE)
  expect_null(captured$fileInput)
  expect_identical(captured$button$id, "lipid-assay1_file")
  expect_identical(captured$button$label, "Select File")
  expect_identical(captured$button$title, "Choose Data File")
  expect_false(captured$button$multiple)
  expect_true(inherits(captured$button$icon, "shiny.tag"))
  expect_identical(captured$path$outputId, "lipid-assay1_path")
  expect_true(captured$path$placeholder)
})

test_that("buildLipidImportAssayInputPanel preserves standard-upload assay-panel wiring", {
  captured <- list()

  rendered <- htmltools::renderTags(
    buildLipidImportAssayInputPanel(
      ns = function(id) paste0("lipid-", id),
      assayNameInputId = "assay2_name",
      assayNameLabel = "Assay Name (Optional)",
      assayNameValue = "",
      assayNamePlaceholder = "e.g., LCMS_Neg",
      shinyFilesButtonId = "assay2_file",
      standardInputId = "assay2_file_std",
      pathOutputId = "assay2_path",
      useShinyFiles = FALSE,
      buildShinyFilesButton = function(...) {
        captured$button <<- list(...)
        shiny::tags$span("shiny-files-button")
      },
      buildPathOutput = function(...) {
        captured$path <<- list(...)
        shiny::tags$span("path-output")
      },
      buildFileInput = function(...) {
        captured$fileInput <<- list(...)
        shiny::tags$span("standard-file-input")
      }
    )
  )$html

  expect_match(rendered, "Assay Name (Optional)", fixed = TRUE)
  expect_match(rendered, "e.g., LCMS_Neg", fixed = TRUE)
  expect_match(rendered, "standard-file-input", fixed = TRUE)
  expect_null(captured$button)
  expect_null(captured$path)
  expect_identical(captured$fileInput$inputId, "lipid-assay2_file_std")
  expect_null(captured$fileInput$label)
  expect_identical(
    captured$fileInput$accept,
    c(".tsv", ".tab", ".txt", ".csv", ".xlsx", ".parquet")
  )
})

test_that("buildLipidImportFileImportSection preserves vendor-selector and assay-panel wiring", {
  captured <- list(assays = list())

  rendered <- htmltools::renderTags(
    buildLipidImportFileImportSection(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = FALSE,
      buildAssayInputPanel = function(...) {
        captured$assays[[length(captured$assays) + 1]] <<- list(...)
        shiny::tags$div(sprintf("assay-panel-%d", length(captured$assays)))
      },
      buildRadioButtons = function(...) {
        captured$radio <<- list(...)
        shiny::tags$div("vendor-radios")
      }
    )
  )$html

  expect_match(rendered, "Step 1: Select Vendor Format", fixed = TRUE)
  expect_match(rendered, "Step 2: Import Data Files", fixed = TRUE)
  expect_match(rendered, "vendor-radios", fixed = TRUE)
  expect_length(captured$assays, 2)
  expect_identical(captured$radio[[1]], "lipid-vendor_format")
  expect_null(captured$radio[[2]])
  expect_identical(
    captured$radio$choices,
    c(
      "MS-DIAL" = "msdial",
      "Progenesis QI" = "progenesis",
      "XCMS" = "xcms",
      "Compound Discoverer" = "compound_discoverer",
      "LipidSearch" = "lipidsearch",
      "Other/Custom" = "custom"
    )
  )
  expect_identical(captured$radio$selected, "msdial")
  expect_true(captured$radio$inline)
  expect_identical(captured$assays[[1]]$assayNameInputId, "assay1_name")
  expect_identical(captured$assays[[1]]$assayNameLabel, "Assay Name")
  expect_identical(captured$assays[[1]]$assayNameValue, "LCMS_Pos")
  expect_identical(captured$assays[[1]]$shinyFilesButtonId, "assay1_file")
  expect_identical(captured$assays[[1]]$standardInputId, "assay1_file_std")
  expect_identical(captured$assays[[1]]$pathOutputId, "assay1_path")
  expect_false(captured$assays[[1]]$useShinyFiles)
  expect_identical(captured$assays[[2]]$assayNameInputId, "assay2_name")
  expect_identical(captured$assays[[2]]$assayNameLabel, "Assay Name (Optional)")
  expect_identical(captured$assays[[2]]$assayNameValue, "")
  expect_identical(captured$assays[[2]]$assayNamePlaceholder, "e.g., LCMS_Neg")
  expect_identical(captured$assays[[2]]$shinyFilesButtonId, "assay2_file")
  expect_identical(captured$assays[[2]]$standardInputId, "assay2_file_std")
  expect_identical(captured$assays[[2]]$pathOutputId, "assay2_path")
  expect_false(captured$assays[[2]]$useShinyFiles)
})

test_that("buildLipidImportColumnMappingSection preserves mapping and validation wiring", {
  rendered <- htmltools::renderTags(
    buildLipidImportColumnMappingSection(
      ns = function(id) paste0("lipid-", id)
    )
  )$html

  expect_match(rendered, "Step 3: Column Mapping", fixed = TRUE)
  expect_match(rendered, "Step 4: Validation", fixed = TRUE)
  expect_match(rendered, "lipid-format_detection_status", fixed = TRUE)
  expect_match(rendered, "lipid-lipid_id_col", fixed = TRUE)
  expect_match(rendered, "lipid-lipid_id_status", fixed = TRUE)
  expect_match(rendered, "lipid-lipid_id_col_custom", fixed = TRUE)
  expect_match(rendered, "lipid-lipid_id_status_custom", fixed = TRUE)
  expect_match(rendered, "lipid-annotation_col", fixed = TRUE)
  expect_match(rendered, "lipid-annotation_status", fixed = TRUE)
  expect_match(rendered, "lipid-annotation_col_custom", fixed = TRUE)
  expect_match(rendered, "lipid-annotation_status_custom", fixed = TRUE)
  expect_match(rendered, "lipid-sample_cols_pattern", fixed = TRUE)
  expect_match(rendered, "Exclude Vendor-Normalized Columns", fixed = TRUE)
  expect_match(rendered, "Sanitize Sample Names", fixed = TRUE)
  expect_match(rendered, "lipid-is_pattern", fixed = TRUE)
  expect_match(rendered, "lipid-sample_columns_display", fixed = TRUE)
  expect_match(rendered, "lipid-available_columns_display", fixed = TRUE)
  expect_match(rendered, "lipid-validation_summary", fixed = TRUE)
  expect_match(rendered, "Import a data file to configure column mappings.", fixed = TRUE)
  expect_match(rendered, "output[&#39;lipid-file_loaded&#39;]", fixed = TRUE)
  expect_match(rendered, "input[&#39;lipid-vendor_format&#39;] != &#39;custom&#39;", fixed = TRUE)
  expect_match(rendered, "input[&#39;lipid-vendor_format&#39;] == &#39;custom&#39;", fixed = TRUE)
  expect_match(rendered, "!output[&#39;lipid-file_loaded&#39;]", fixed = TRUE)
})

test_that("buildLipidImportProcessFooterSection preserves process-button and status wiring", {
  captured <- list()

  rendered <- htmltools::renderTags(
    buildLipidImportProcessFooterSection(
      ns = function(id) paste0("lipid-", id),
      buildActionButton = function(...) {
        captured$button <<- list(...)
        shiny::tags$div("process-button")
      },
      buildUiOutput = function(...) {
        captured$status <<- list(...)
        shiny::tags$div("import-status-output")
      }
    )
  )$html

  expect_match(rendered, "process-button", fixed = TRUE)
  expect_match(rendered, "import-status-output", fixed = TRUE)
  expect_identical(captured$button[[1]], "lipid-process_import")
  expect_identical(captured$button[[2]], "Process Imported Data")
  expect_identical(captured$button$class, "btn-success")
  expect_identical(captured$button$width, "100%")
  expect_s3_class(captured$button$icon, "shiny.tag")
  expect_identical(captured$status[[1]], "lipid-import_status")
})

test_that("buildLipidImportModulePanel preserves section wiring and footer layout", {
  captured <- list(columns = integer())

  rendered <- htmltools::renderTags(
    buildLipidImportModulePanel(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = TRUE,
      buildFileImportSection = function(...) {
        captured$file_import <<- list(...)
        shiny::tags$div("file-import-section")
      },
      buildColumnMappingSection = function(...) {
        captured$column_mapping <<- list(...)
        shiny::tags$div("column-mapping-section")
      },
      buildProcessFooterSection = function(...) {
        captured$footer <<- list(...)
        shiny::tags$div("process-footer-section")
      },
      buildWellPanel = function(...) shiny::tags$section(...),
      buildHeader = function(...) shiny::tags$h3(...),
      buildFluidRow = function(...) shiny::tags$div(class = "fluid-row", ...),
      buildColumn = function(width, ...) {
        captured$columns <<- c(captured$columns, width)
        shiny::tags$div(class = sprintf("column-%s", width), ...)
      },
      buildHr = function(...) shiny::tags$hr(...)
    )
  )$html

  expect_match(rendered, "Lipidomics Data Import", fixed = TRUE)
  expect_match(rendered, "file-import-section", fixed = TRUE)
  expect_match(rendered, "column-mapping-section", fixed = TRUE)
  expect_match(rendered, "process-footer-section", fixed = TRUE)
  expect_identical(captured$file_import$ns("vendor_format"), "lipid-vendor_format")
  expect_true(captured$file_import$useShinyFiles)
  expect_identical(captured$column_mapping$ns("validation_summary"), "lipid-validation_summary")
  expect_identical(captured$footer$ns("import_status"), "lipid-import_status")
  expect_identical(captured$columns, c(6, 6))
})

test_that("buildLipidImportUiShell preserves shinyjs wrapping and module-panel wiring", {
  captured <- list(columns = integer(), shinyjs_calls = 0L)

  rendered <- htmltools::renderTags(
    buildLipidImportUiShell(
      ns = function(id) paste0("lipid-", id),
      useShinyFiles = FALSE,
      buildModulePanel = function(...) {
        captured$panel <<- list(...)
        shiny::tags$div("module-panel")
      },
      buildTagList = function(...) shiny::tags$div(class = "tag-list", ...),
      useShinyjs = function(...) {
        captured$shinyjs_calls <<- captured$shinyjs_calls + 1L
        shiny::tags$div("use-shinyjs")
      },
      buildFluidRow = function(...) shiny::tags$div(class = "fluid-row", ...),
      buildColumn = function(width, ...) {
        captured$columns <<- c(captured$columns, width)
        shiny::tags$div(class = sprintf("column-%s", width), ...)
      }
    )
  )$html

  expect_match(rendered, "use-shinyjs", fixed = TRUE)
  expect_match(rendered, "module-panel", fixed = TRUE)
  expect_identical(captured$shinyjs_calls, 1L)
  expect_identical(captured$panel$ns("process_import"), "lipid-process_import")
  expect_false(captured$panel$useShinyFiles)
  expect_identical(captured$columns, 12)
})

test_that("loadLipidImportAssayPreview preserves the detected lipidsearch import payload", {
  captured <- list()
  log_messages <- character()
  assay_data <- data.frame(
    `Lipid Name` = c("PC 34:1", "TG 52:2"),
    Annotation = c("class_a", "class_b"),
    `Sample A` = c(10, 20),
    check.names = FALSE
  )

  preview <- loadLipidImportAssayPreview(
    assay1File = "/tmp/lipidsearch.csv",
    readHeaders = function(path) {
      captured$header_path <<- path
      c("Lipid Name", "Annotation", "Sample A")
    },
    detectFormat = function(headers, filename) {
      captured$headers <<- headers
      captured$filename <<- filename
      list(format = "lipidsearch", confidence = 0.83)
    },
    importMsdial = function(path) {
      captured$msdial_path <<- path
      stop("MS-DIAL importer should not be used")
    },
    importLipidSearch = function(path) {
      captured$lipidsearch_path <<- path
      list(
        data = assay_data,
        detected_columns = list(
          lipid_id = "Lipid Name",
          annotation = "Annotation"
        ),
        is_pattern = "^ISTD_",
        sample_columns = "Sample A"
      )
    },
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(captured$header_path, "/tmp/lipidsearch.csv")
  expect_identical(captured$headers, c("Lipid Name", "Annotation", "Sample A"))
  expect_identical(captured$filename, "lipidsearch.csv")
  expect_false("msdial_path" %in% names(captured))
  expect_identical(captured$lipidsearch_path, "/tmp/lipidsearch.csv")
  expect_identical(preview$headers, c("Lipid Name", "Annotation", "Sample A"))
  expect_identical(preview$detectedFormat, "lipidsearch")
  expect_identical(preview$formatConfidence, 0.83)
  expect_identical(preview$importResult$data, assay_data)
  expect_identical(preview$assayData, assay_data)
  expect_identical(preview$updates$lipidId$selected, "Lipid Name")
  expect_identical(preview$updates$annotation$selected, "Annotation")
  expect_identical(
    unname(preview$updates$annotation$choices),
    c("", "Lipid Name", "Annotation", "Sample A")
  )
  expect_identical(preview$updates$isPattern, "^ISTD_")
  expect_identical(
    log_messages,
    "Imported assay 1: 2 rows, 3 columns, format: lipidsearch"
  )
})

test_that("loadLipidImportAssayPreview preserves fallback import and omitted IS-pattern updates", {
  captured <- list()
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 11)
  )

  preview <- loadLipidImportAssayPreview(
    assay1File = "/tmp/custom.tsv",
    readHeaders = function(path) {
      captured$header_path <<- path
      c("lipid_id", "sample_a")
    },
    detectFormat = function(headers, filename) {
      captured$headers <<- headers
      captured$filename <<- filename
      list(format = "custom_vendor", confidence = 0.27)
    },
    importMsdial = function(path) {
      captured$msdial_path <<- path
      list(
        data = assay_data,
        detected_columns = list(
          lipid_id = "lipid_id",
          annotation = ""
        ),
        is_pattern = NA_character_,
        sample_columns = "sample_a"
      )
    },
    importLipidSearch = function(path) {
      captured$lipidsearch_path <<- path
      stop("LipidSearch importer should not be used")
    },
    logInfo = function(...) NULL
  )

  expect_identical(captured$header_path, "/tmp/custom.tsv")
  expect_identical(captured$headers, c("lipid_id", "sample_a"))
  expect_identical(captured$filename, "custom.tsv")
  expect_identical(captured$msdial_path, "/tmp/custom.tsv")
  expect_false("lipidsearch_path" %in% names(captured))
  expect_identical(preview$detectedFormat, "custom_vendor")
  expect_identical(preview$assayData, assay_data)
  expect_identical(preview$updates$lipidId$selected, "lipid_id")
  expect_identical(preview$updates$annotation$selected, "")
  expect_null(preview$updates$isPattern)
})

test_that("loadLipidImportAssayPreview preserves the unreadable-header failure", {
  expect_error(
    loadLipidImportAssayPreview(
      assay1File = "/tmp/missing.tsv",
      readHeaders = function(path) character(0),
      logInfo = function(...) NULL
    ),
    "Could not read headers from file"
  )
})

test_that("handleLipidImportDataPreviewLoad preserves preview-load orchestration and state application", {
  local_data <- list(assay1_file = "/tmp/assay1.tsv")
  captured <- list()
  import_preview <- list(headers = c("lipid_id", "sample_a"))

  result <- handleLipidImportDataPreviewLoad(
    session = "mock-session",
    localData = local_data,
    loadPreview = function(assay1File) {
      captured$assay1_file <<- assay1File
      import_preview
    },
    applyPreview = function(session, localData, importPreview) {
      captured$session <<- session
      captured$local_data <<- localData
      captured$import_preview <<- importPreview
      "applied"
    },
    handleImportError = function(error) {
      stop(error)
    }
  )

  expect_identical(result, "applied")
  expect_identical(captured$assay1_file, "/tmp/assay1.tsv")
  expect_identical(captured$session, "mock-session")
  expect_identical(captured$local_data, local_data)
  expect_identical(captured$import_preview, import_preview)
})

test_that("handleLipidImportDataPreviewLoad preserves error forwarding and req guarding", {
  local_data <- list(assay1_file = "/tmp/assay1.tsv")
  apply_called <- FALSE
  forwarded_error <- NULL

  result <- handleLipidImportDataPreviewLoad(
    session = "mock-session",
    localData = local_data,
    loadPreview = function(assay1File) {
      stop("preview failed")
    },
    applyPreview = function(...) {
      apply_called <<- TRUE
      invisible(NULL)
    },
    handleImportError = function(error) {
      forwarded_error <<- error$message
      FALSE
    }
  )

  expect_false(result)
  expect_false(apply_called)
  expect_identical(forwarded_error, "preview failed")

  local_data_missing <- list(assay1_file = NULL)
  load_called <- FALSE

  expect_error(
    handleLipidImportDataPreviewLoad(
      session = "mock-session",
      localData = local_data_missing,
      loadPreview = function(...) {
        load_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(load_called)
})

test_that("buildLipidImportAssay1SelectedCallback preserves assay-1 selection dispatch", {
  local_data <- shiny::reactiveValues(assay1_file = "/tmp/assay1.csv")
  captured <- list()

  callback <- buildLipidImportAssay1SelectedCallback(
    session = "mock-session",
    localData = local_data,
    runPreviewLoad = function(session, localData) {
      captured$session <<- session
      captured$local_data <<- localData
      "preview-loaded"
    }
  )

  expect_true(is.function(callback))
  expect_identical(callback(), "preview-loaded")
  expect_identical(captured$session, "mock-session")
  expect_identical(captured$local_data, local_data)
})

test_that("applyLipidImportPreviewToModuleState preserves preview application and input updates", {
  local_data <- shiny::reactiveValues()
  select_updates <- list()
  text_updates <- list()
  import_preview <- list(
    headers = c("Lipid Name", "Annotation", "Sample A"),
    detectedFormat = "lipidsearch",
    formatConfidence = 0.83,
    importResult = list(
      data = data.frame(
        `Lipid Name` = c("PC 34:1", "TG 52:2"),
        check.names = FALSE
      ),
      detected_columns = list(
        lipid_id = "Lipid Name",
        annotation = "Annotation"
      )
    ),
    assayData = data.frame(
      `Lipid Name` = c("PC 34:1", "TG 52:2"),
      check.names = FALSE
    ),
    updates = list(
      lipidId = list(
        choices = c("Lipid Name", "Annotation", "Sample A"),
        selected = "Lipid Name"
      ),
      annotation = list(
        choices = c("(None)" = "", "Lipid Name", "Annotation", "Sample A"),
        selected = "Annotation"
      ),
      isPattern = "^ISTD_"
    )
  )

  result <- applyLipidImportPreviewToModuleState(
    session = "mock-session",
    localData = local_data,
    importPreview = import_preview,
    updateSelectInput = function(session, inputId, choices, selected) {
      select_updates[[length(select_updates) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextInput = function(session, inputId, value) {
      text_updates[[length(text_updates) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_identical(result, import_preview)
  expect_identical(shiny::isolate(local_data$all_headers), import_preview$headers)
  expect_identical(shiny::isolate(local_data$detected_format), "lipidsearch")
  expect_identical(shiny::isolate(local_data$format_confidence), 0.83)
  expect_identical(shiny::isolate(local_data$assay1_import_result), import_preview$importResult)
  expect_identical(shiny::isolate(local_data$assay1_data), import_preview$assayData)
  expect_length(select_updates, 2L)
  expect_identical(select_updates[[1]]$session, "mock-session")
  expect_identical(select_updates[[1]]$inputId, "lipid_id_col")
  expect_identical(select_updates[[1]]$choices, import_preview$updates$lipidId$choices)
  expect_identical(select_updates[[1]]$selected, "Lipid Name")
  expect_identical(select_updates[[2]]$inputId, "annotation_col")
  expect_identical(select_updates[[2]]$choices, import_preview$updates$annotation$choices)
  expect_identical(select_updates[[2]]$selected, "Annotation")
  expect_identical(
    text_updates,
    list(list(
      session = "mock-session",
      inputId = "is_pattern",
      value = "^ISTD_"
    ))
  )
})

test_that("applyLipidImportPreviewToModuleState keeps optional IS-pattern updates omitted", {
  local_data <- shiny::reactiveValues()
  text_updates <- list()
  import_preview <- list(
    headers = c("lipid_id", "sample_a"),
    detectedFormat = "msdial",
    formatConfidence = 0.27,
    importResult = list(
      data = data.frame(lipid_id = c("L1", "L2"), sample_a = c(10, 11)),
      detected_columns = list(lipid_id = "lipid_id", annotation = "")
    ),
    assayData = data.frame(lipid_id = c("L1", "L2"), sample_a = c(10, 11)),
    updates = list(
      lipidId = list(
        choices = c("lipid_id", "sample_a"),
        selected = "lipid_id"
      ),
      annotation = list(
        choices = c("(None)" = "", "lipid_id", "sample_a"),
        selected = ""
      ),
      isPattern = NULL
    )
  )

  applyLipidImportPreviewToModuleState(
    session = "mock-session",
    localData = local_data,
    importPreview = import_preview,
    updateSelectInput = function(...) invisible(NULL),
    updateTextInput = function(...) {
      text_updates[[length(text_updates) + 1L]] <<- list(...)
      invisible(NULL)
    }
  )

  expect_identical(shiny::isolate(local_data$detected_format), "msdial")
  expect_identical(shiny::isolate(local_data$assay1_import_result), import_preview$importResult)
  expect_identical(shiny::isolate(local_data$assay1_data), import_preview$assayData)
  expect_length(text_updates, 0L)
})

test_that("handleLipidImportDataImportError preserves the import error logging and notification payload", {
  log_messages <- character()
  notifications <- list()

  result <- handleLipidImportDataImportError(
    error = simpleError("preview failed"),
    logError = function(message) {
      log_messages <<- c(log_messages, message)
      invisible(NULL)
    },
    notify = function(message, type) {
      notifications[[length(notifications) + 1L]] <<- list(
        message = message,
        type = type
      )
      invisible(NULL)
    }
  )

  expect_false(result)
  expect_identical(log_messages, "Error importing data: preview failed")
  expect_identical(
    notifications,
    list(list(
      message = "Error importing data: preview failed",
      type = "error"
    ))
  )
})

test_that("handleLipidImportFileSelection preserves the chooser observer success path", {
  captured <- list()
  selected_paths <- character()

  result <- handleLipidImportFileSelection(
    fileInput = list(path = "assay1"),
    volumes = c(home = "/tmp"),
    onPathSelected = function(path) {
      selected_paths <<- c(selected_paths, path)
      invisible(NULL)
    },
    parseFilePaths = function(volumes, file_input) {
      captured$volumes <<- volumes
      captured$file_input <<- file_input
      data.frame(datapath = "/tmp/lipid-assay-1.tsv", stringsAsFactors = FALSE)
    },
    logError = function(...) stop("unexpected error path")
  )

  expect_true(result)
  expect_identical(captured$volumes, c(home = "/tmp"))
  expect_identical(captured$file_input, list(path = "assay1"))
  expect_identical(selected_paths, "/tmp/lipid-assay-1.tsv")
})

test_that("handleLipidImportFileSelection keeps ignored and parse-error chooser paths stable", {
  parse_calls <- 0L
  selected_paths <- character()
  log_messages <- character()

  ignored_result <- handleLipidImportFileSelection(
    fileInput = 0L,
    volumes = c(home = "/tmp"),
    onPathSelected = function(path) {
      selected_paths <<- c(selected_paths, path)
      invisible(NULL)
    },
    parseFilePaths = function(...) {
      parse_calls <<- parse_calls + 1L
      data.frame(datapath = "/tmp/ignored.tsv", stringsAsFactors = FALSE)
    },
    logError = function(message) {
      log_messages <<- c(log_messages, message)
      invisible(NULL)
    }
  )

  error_result <- handleLipidImportFileSelection(
    fileInput = list(path = "assay2"),
    volumes = c(home = "/tmp"),
    onPathSelected = function(path) {
      selected_paths <<- c(selected_paths, path)
      invisible(NULL)
    },
    parseFilePaths = function(...) {
      parse_calls <<- parse_calls + 1L
      stop("boom")
    },
    logError = function(message) {
      log_messages <<- c(log_messages, message)
      invisible(NULL)
    }
  )

  expect_false(ignored_result)
  expect_false(error_result)
  expect_identical(parse_calls, 1L)
  expect_identical(selected_paths, character())
  expect_identical(log_messages, "Error parsing file path: boom")
})

test_that("handleLipidImportProcessRequest preserves the process observer payload", {
  workflow_data <- new.env(parent = emptyenv())
  assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 11)
  )
  process_args <- NULL

  result <- handleLipidImportProcessRequest(
    workflowData = workflow_data,
    assay1Name = "LCMS_Pos",
    assay1Data = assay1_data,
    assay2File = "/tmp/lcms-neg.tsv",
    assay2Name = "LCMS_Neg",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = "sample_a",
    isPattern = "^ISTD_",
    sanitizeNames = TRUE,
    processImport = function(...) {
      process_args <<- list(...)
      "processed"
    }
  )

  expect_identical(result, "processed")
  expect_identical(process_args$workflowData, workflow_data)
  expect_identical(process_args$assay1Name, "LCMS_Pos")
  expect_identical(process_args$assay1Data, assay1_data)
  expect_identical(process_args$assay2File, "/tmp/lcms-neg.tsv")
  expect_identical(process_args$assay2Name, "LCMS_Neg")
  expect_identical(process_args$vendorFormat, "custom")
  expect_identical(process_args$detectedFormat, "msdial")
  expect_identical(process_args$lipidIdCol, "lipid_id")
  expect_identical(process_args$annotationCol, "annotation")
  expect_identical(process_args$sampleColumns, "sample_a")
  expect_identical(process_args$isPattern, "^ISTD_")
  expect_true(process_args$sanitizeNames)
})

test_that("handleLipidImportProcessRequest keeps missing prerequisites guarded by req", {
  process_called <- FALSE

  expect_error(
    handleLipidImportProcessRequest(
      workflowData = new.env(parent = emptyenv()),
      assay1Name = "LCMS_Pos",
      assay1Data = NULL,
      assay2File = NULL,
      assay2Name = "",
      vendorFormat = "msdial",
      detectedFormat = "msdial",
      lipidIdCol = "lipid_id",
      annotationCol = "",
      sampleColumns = "sample_a",
      isPattern = "",
      sanitizeNames = FALSE,
      processImport = function(...) {
        process_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_error(
    handleLipidImportProcessRequest(
      workflowData = new.env(parent = emptyenv()),
      assay1Name = "LCMS_Pos",
      assay1Data = data.frame(lipid_id = "L1"),
      assay2File = NULL,
      assay2Name = "",
      vendorFormat = "msdial",
      detectedFormat = "msdial",
      lipidIdCol = NULL,
      annotationCol = "",
      sampleColumns = "sample_a",
      isPattern = "",
      sanitizeNames = FALSE,
      processImport = function(...) {
        process_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(process_called)
})

test_that("registerLipidImportProcessObserver preserves the process observer shell payload", {
  input <- list(
    process_import = 3L,
    assay1_name = "LCMS_Pos",
    assay2_name = "LCMS_Neg",
    vendor_format = "custom",
    is_pattern = "^ISTD_",
    sanitize_names = TRUE
  )
  local_data <- list(
    assay1_data = data.frame(
      lipid_id = c("L1", "L2"),
      sample_a = c(10, 11)
    ),
    assay2_file = "/tmp/lcms-neg.tsv",
    detected_format = "msdial"
  )
  workflow_data <- new.env(parent = emptyenv())
  observer_event <- NULL
  process_args <- NULL
  getter_calls <- character()

  registerLipidImportProcessObserver(
    input = input,
    workflowData = workflow_data,
    localData = local_data,
    getLipidIdCol = function() {
      getter_calls <<- c(getter_calls, "lipid")
      "lipid_id"
    },
    getAnnotationCol = function() {
      getter_calls <<- c(getter_calls, "annotation")
      "annotation"
    },
    getSampleColumns = function() {
      getter_calls <<- c(getter_calls, "sample")
      "sample_a"
    },
    handleProcessRequest = function(...) {
      process_args <<- list(...)
      "processed"
    },
    observeEvent = function(event, handler) {
      observer_event <<- event
      force(handler)
      "registered"
    }
  )

  expect_identical(observer_event, 3L)
  expect_identical(getter_calls, c("lipid", "annotation", "sample"))
  expect_identical(process_args$workflowData, workflow_data)
  expect_identical(process_args$assay1Name, "LCMS_Pos")
  expect_identical(process_args$assay1Data, local_data$assay1_data)
  expect_identical(process_args$assay2File, "/tmp/lcms-neg.tsv")
  expect_identical(process_args$assay2Name, "LCMS_Neg")
  expect_identical(process_args$vendorFormat, "custom")
  expect_identical(process_args$detectedFormat, "msdial")
  expect_identical(process_args$lipidIdCol, "lipid_id")
  expect_identical(process_args$annotationCol, "annotation")
  expect_identical(process_args$sampleColumns, "sample_a")
  expect_identical(process_args$isPattern, "^ISTD_")
  expect_true(process_args$sanitizeNames)
})

test_that("handleLipidImportValidationSummaryRender preserves the validation summary payload", {
  assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 20)
  )
  validation_args <- list()
  built_validation <- NULL
  validation_result <- list(
    valid = TRUE,
    summary = list(n_lipids = 2L, n_samples = 1L, pct_missing = 0),
    warnings = "warn",
    errors = character(0)
  )

  result <- handleLipidImportValidationSummaryRender(
    assay1Data = assay1_data,
    lipidIdCol = "lipid_id",
    sampleColumns = "sample_a",
    validateMapping = function(...) {
      validation_args <<- list(...)
      validation_result
    },
    buildSummary = function(validation) {
      built_validation <<- validation
      "rendered-summary"
    }
  )

  expect_identical(result, "rendered-summary")
  expect_identical(validation_args$data, assay1_data)
  expect_identical(validation_args$lipid_id_column, "lipid_id")
  expect_identical(validation_args$sample_columns, "sample_a")
  expect_identical(built_validation, validation_result)
})

test_that("handleLipidImportValidationSummaryRender keeps missing prerequisites guarded by req", {
  validate_called <- FALSE

  expect_error(
    handleLipidImportValidationSummaryRender(
      assay1Data = NULL,
      lipidIdCol = "lipid_id",
      sampleColumns = "sample_a",
      validateMapping = function(...) {
        validate_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_error(
    handleLipidImportValidationSummaryRender(
      assay1Data = data.frame(lipid_id = "L1"),
      lipidIdCol = NULL,
      sampleColumns = "sample_a",
      validateMapping = function(...) {
        validate_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(validate_called)
})

test_that("handleLipidImportSampleColumnsDisplayRender preserves the sample-column preview payload", {
  assay1_import_result <- list(sample_columns = c("Sample_A", "Sample_B", "Sample_C"))
  preview_args <- NULL

  result <- handleLipidImportSampleColumnsDisplayRender(
    assay1ImportResult = assay1_import_result,
    formatPreviewText = function(...) {
      preview_args <<- list(...)
      "sample-preview"
    }
  )

  expect_identical(result, "sample-preview")
  expect_identical(preview_args$columnNames, c("Sample_A", "Sample_B", "Sample_C"))
  expect_identical(preview_args$maxColumns, 10)
})

test_that("handleLipidImportSampleColumnsDisplayRender keeps missing prerequisites guarded by req", {
  preview_called <- FALSE

  expect_error(
    handleLipidImportSampleColumnsDisplayRender(
      assay1ImportResult = NULL,
      formatPreviewText = function(...) {
        preview_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(preview_called)
})

test_that("handleLipidImportAssayPathRender preserves the selected assay path text", {
  expect_identical(
    handleLipidImportAssayPathRender("/tmp/assay1.tsv"),
    "/tmp/assay1.tsv"
  )
})

test_that("handleLipidImportAssayPathRender keeps empty and NULL passthrough stable", {
  expect_identical(handleLipidImportAssayPathRender(""), "")
  expect_null(handleLipidImportAssayPathRender(NULL))
})

test_that("handleLipidImportSelectedAssayPathAssignment preserves selected-path assignment and follow-up payloads", {
  assigned_path <- NULL
  rendered_path <- NULL
  rendered_inputs <- character()
  on_selected_calls <- 0L

  result <- handleLipidImportSelectedAssayPathAssignment(
    selectedPath = "/tmp/assay1.tsv",
    assignSelectedPath = function(path) {
      assigned_path <<- path
      invisible(NULL)
    },
    setRenderedPath = function(pathText) {
      rendered_path <<- pathText
      invisible(NULL)
    },
    onSelected = function() {
      on_selected_calls <<- on_selected_calls + 1L
      invisible(NULL)
    },
    renderText = function(text) {
      paste0("rendered:", text)
    },
    renderAssayPath = function(path) {
      rendered_inputs <<- c(rendered_inputs, path)
      paste0("path=", path)
    }
  )

  expect_identical(result, "/tmp/assay1.tsv")
  expect_identical(assigned_path, "/tmp/assay1.tsv")
  expect_identical(rendered_inputs, "/tmp/assay1.tsv")
  expect_identical(rendered_path, "rendered:path=/tmp/assay1.tsv")
  expect_identical(on_selected_calls, 1L)
})

test_that("handleLipidImportSelectedAssayPathAssignment keeps optional follow-up and empty path passthrough stable", {
  assigned_path <- "unset"
  rendered_path <- "unset"
  on_selected_calls <- 0L

  result <- handleLipidImportSelectedAssayPathAssignment(
    selectedPath = "",
    assignSelectedPath = function(path) {
      assigned_path <<- path
      invisible(NULL)
    },
    setRenderedPath = function(pathText) {
      rendered_path <<- pathText
      invisible(NULL)
    },
    onSelected = NULL,
    renderText = function(text) {
      text
    },
    renderAssayPath = function(path) {
      path
    }
  )

  expect_identical(result, "")
  expect_identical(assigned_path, "")
  expect_identical(rendered_path, "")
  expect_identical(on_selected_calls, 0L)
})

test_that("handleLipidImportAssayFileSelectionEvent preserves paired chooser observer payloads", {
  captured <- list()
  assigned_path <- NULL
  rendered_path <- NULL
  on_selected_calls <- 0L

  result <- handleLipidImportAssayFileSelectionEvent(
    fileInput = list(path = "assay1"),
    volumes = c(home = "/tmp"),
    assignSelectedPath = function(path) {
      assigned_path <<- path
      invisible(NULL)
    },
    setRenderedPath = function(pathText) {
      rendered_path <<- pathText
      invisible(NULL)
    },
    onSelected = function() {
      on_selected_calls <<- on_selected_calls + 1L
      invisible(NULL)
    },
    handleFileSelection = function(fileInput, volumes, onPathSelected) {
      captured$file_input <<- fileInput
      captured$volumes <<- volumes
      onPathSelected("/tmp/lipid-assay-1.tsv")
      TRUE
    },
    assignAssayPath = function(selectedPath, assignSelectedPath, setRenderedPath, onSelected) {
      captured$selected_path <<- selectedPath
      assignSelectedPath(selectedPath)
      setRenderedPath(paste0("rendered:", selectedPath))
      onSelected()
      "assigned"
    }
  )

  expect_true(result)
  expect_identical(captured$file_input, list(path = "assay1"))
  expect_identical(captured$volumes, c(home = "/tmp"))
  expect_identical(captured$selected_path, "/tmp/lipid-assay-1.tsv")
  expect_identical(assigned_path, "/tmp/lipid-assay-1.tsv")
  expect_identical(rendered_path, "rendered:/tmp/lipid-assay-1.tsv")
  expect_identical(on_selected_calls, 1L)
})

test_that("handleLipidImportAssayFileSelectionEvent keeps false and optional follow-up paths stable", {
  assign_calls <- 0L

  result <- handleLipidImportAssayFileSelectionEvent(
    fileInput = list(path = "assay2"),
    volumes = c(home = "/tmp"),
    assignSelectedPath = function(path) {
      assign_calls <<- assign_calls + 1L
      invisible(NULL)
    },
    setRenderedPath = function(pathText) {
      assign_calls <<- assign_calls + 1L
      invisible(NULL)
    },
    onSelected = NULL,
    handleFileSelection = function(fileInput, volumes, onPathSelected) {
      FALSE
    },
    assignAssayPath = function(...) {
      assign_calls <<- assign_calls + 1L
      invisible(NULL)
    }
  )

  expect_false(result)
  expect_identical(assign_calls, 0L)
})

test_that("registerLipidImportAssayFileSelectionObserver preserves observer trigger and chooser payload forwarding", {
  input <- list(assay1_file = list(path = "assay1"))
  observer_event <- NULL
  handled_args <- NULL

  registerLipidImportAssayFileSelectionObserver(
    input = input,
    fileInputId = "assay1_file",
    volumes = c(home = "/tmp"),
    assignSelectedPath = function(path) {
      invisible(path)
    },
    setRenderedPath = function(pathText) {
      invisible(pathText)
    },
    onSelected = function() {
      invisible(NULL)
    },
    handleSelectionEvent = function(...) {
      handled_args <<- list(...)
      "handled"
    },
    observeEvent = function(event, handler) {
      observer_event <<- event
      force(handler)
      "registered"
    }
  )

  expect_identical(observer_event, input$assay1_file)
  expect_identical(handled_args$fileInput, input$assay1_file)
  expect_identical(handled_args$volumes, c(home = "/tmp"))
  expect_true(is.function(handled_args$assignSelectedPath))
  expect_true(is.function(handled_args$setRenderedPath))
  expect_true(is.function(handled_args$onSelected))
})

test_that("registerLipidImportAssayFileSelectionObserver keeps optional follow-up forwarding stable", {
  input <- list(assay2_file = list(path = "assay2"))
  handled_args <- NULL

  registerLipidImportAssayFileSelectionObserver(
    input = input,
    fileInputId = "assay2_file",
    volumes = c(home = "/tmp"),
    assignSelectedPath = function(path) {
      invisible(path)
    },
    setRenderedPath = function(pathText) {
      invisible(pathText)
    },
    onSelected = NULL,
    handleSelectionEvent = function(...) {
      handled_args <<- list(...)
      "handled"
    },
    observeEvent = function(event, handler) {
      force(handler)
      "registered"
    }
  )

  expect_identical(handled_args$fileInput, input$assay2_file)
  expect_null(handled_args$onSelected)
})

test_that("registerLipidImportShinyFileInputs preserves paired chooser setup and observer routing", {
  local_data <- shiny::reactiveValues()
  output <- new.env(parent = emptyenv())
  setup_calls <- list()
  observer_calls <- list()
  assay1_selected_calls <- 0L

  result <- registerLipidImportShinyFileInputs(
    input = structure(list(), class = "mock_input"),
    output = output,
    session = "mock-session",
    localData = local_data,
    volumes = c(home = "/tmp"),
    onAssay1Selected = function() {
      assay1_selected_calls <<- assay1_selected_calls + 1L
      invisible(NULL)
    },
    setupAssayFileChooser = function(...) {
      setup_calls[[length(setup_calls) + 1L]] <<- list(...)
      invisible(NULL)
    },
    registerSelectionObserver = function(...) {
      observer_calls[[length(observer_calls) + 1L]] <<- list(...)
      invisible(NULL)
    }
  )

  expect_identical(result, c(home = "/tmp"))
  expect_length(setup_calls, 2L)
  expect_identical(setup_calls[[1]]$fileInputId, "assay1_file")
  expect_identical(setup_calls[[2]]$fileInputId, "assay2_file")
  expect_length(observer_calls, 2L)
  expect_identical(observer_calls[[1]]$fileInputId, "assay1_file")
  expect_identical(observer_calls[[2]]$fileInputId, "assay2_file")
  expect_true(is.function(observer_calls[[1]]$assignSelectedPath))
  expect_true(is.function(observer_calls[[1]]$setRenderedPath))
  expect_true(is.function(observer_calls[[1]]$onSelected))
  expect_true(is.function(observer_calls[[2]]$assignSelectedPath))
  expect_true(is.function(observer_calls[[2]]$setRenderedPath))
  expect_null(observer_calls[[2]]$onSelected)

  observer_calls[[1]]$assignSelectedPath("/tmp/lcms-pos.tsv")
  observer_calls[[1]]$setRenderedPath("LCMS Pos")
  observer_calls[[1]]$onSelected()
  observer_calls[[2]]$assignSelectedPath("/tmp/lcms-neg.tsv")
  observer_calls[[2]]$setRenderedPath("LCMS Neg")

  expect_identical(shiny::isolate(local_data$assay1_file), "/tmp/lcms-pos.tsv")
  expect_identical(shiny::isolate(local_data$assay2_file), "/tmp/lcms-neg.tsv")
  expect_identical(output$assay1_path, "LCMS Pos")
  expect_identical(output$assay2_path, "LCMS Neg")
  expect_identical(assay1_selected_calls, 1L)
})

test_that("setupLipidImportShinyFileInputs preserves the shinyFiles setup shell when enabled", {
  local_data <- shiny::reactiveValues()
  prepared_volumes <- c(cache = "/tmp/cache")
  prepare_calls <- list()
  callback_args <- NULL
  register_args <- NULL
  assay1_selected <- function() "assay1-selected"

  result <- setupLipidImportShinyFileInputs(
    useShinyFiles = TRUE,
    input = structure(list(), class = "mock_input"),
    output = new.env(parent = emptyenv()),
    session = "mock-session",
    localData = local_data,
    volumes = c(home = "/tmp/home"),
    prepareVolumes = function(volumes) {
      prepare_calls[[length(prepare_calls) + 1L]] <<- volumes
      prepared_volumes
    },
    buildAssay1Selected = function(session, localData) {
      callback_args <<- list(session = session, localData = localData)
      assay1_selected
    },
    registerInputs = function(...) {
      register_args <<- list(...)
      invisible(NULL)
    }
  )

  expect_identical(result, prepared_volumes)
  expect_identical(prepare_calls, list(c(home = "/tmp/home")))
  expect_identical(
    callback_args,
    list(session = "mock-session", localData = local_data)
  )
  expect_identical(register_args$volumes, prepared_volumes)
  expect_identical(register_args$onAssay1Selected, assay1_selected)
})

test_that("setupLipidImportShinyFileInputs skips setup work when shinyFiles is unavailable", {
  prepare_called <- FALSE
  build_called <- FALSE
  register_called <- FALSE
  existing_volumes <- c(home = "/tmp/home")

  result <- setupLipidImportShinyFileInputs(
    useShinyFiles = FALSE,
    input = structure(list(), class = "mock_input"),
    output = new.env(parent = emptyenv()),
    session = "mock-session",
    localData = shiny::reactiveValues(),
    volumes = existing_volumes,
    prepareVolumes = function(volumes) {
      prepare_called <<- TRUE
      volumes
    },
    buildAssay1Selected = function(session, localData) {
      build_called <<- TRUE
      function() invisible(list(session = session, localData = localData))
    },
    registerInputs = function(...) {
      register_called <<- TRUE
      invisible(NULL)
    }
  )

  expect_identical(result, existing_volumes)
  expect_false(prepare_called)
  expect_false(build_called)
  expect_false(register_called)
})

test_that("setupLipidImportAssayFileChooser preserves chooser setup payloads", {
  captured <- list()

  result <- setupLipidImportAssayFileChooser(
    input = structure(list(), class = "mock_input"),
    fileInputId = "assay1_file",
    volumes = c(home = "/tmp"),
    session = "mock-session",
    shinyFileChoose = function(input, id, roots, session, filetypes) {
      captured$input <<- input
      captured$id <<- id
      captured$roots <<- roots
      captured$session <<- session
      captured$filetypes <<- filetypes
      "chooser-setup"
    }
  )

  expect_identical(result, "chooser-setup")
  expect_s3_class(captured$input, "mock_input")
  expect_identical(captured$id, "assay1_file")
  expect_identical(captured$roots, c(home = "/tmp"))
  expect_identical(captured$session, "mock-session")
  expect_identical(
    captured$filetypes,
    c("tsv", "tab", "txt", "csv", "xlsx", "parquet")
  )
})

test_that("setupLipidImportAssayFileChooser keeps custom filetype overrides stable", {
  captured_filetypes <- NULL

  result <- setupLipidImportAssayFileChooser(
    input = NULL,
    fileInputId = "assay2_file",
    volumes = c(data = "/srv/data"),
    session = "mock-session",
    filetypes = c("csv", "parquet"),
    shinyFileChoose = function(input, id, roots, session, filetypes) {
      captured_filetypes <<- filetypes
      list(id = id, roots = roots)
    }
  )

  expect_identical(result, list(id = "assay2_file", roots = c(data = "/srv/data")))
  expect_identical(captured_filetypes, c("csv", "parquet"))
})

test_that("prepareLipidImportShinyFileVolumes initializes missing volumes through getVolumes and logs diagnostics", {
  emitted_messages <- character()
  get_volumes_calls <- 0L

  result <- prepareLipidImportShinyFileVolumes(
    volumes = NULL,
    getVolumes = function() {
      get_volumes_calls <<- get_volumes_calls + 1L
      function() {
        c(home = "/tmp", data = "/srv/data")
      }
    },
    emitMessage = function(messageText) {
      emitted_messages <<- c(emitted_messages, messageText)
      invisible(NULL)
    }
  )

  expect_identical(result, c(home = "/tmp", data = "/srv/data"))
  expect_identical(get_volumes_calls, 1L)
  expect_identical(
    emitted_messages,
    c(
      "   mod_lipid_import_server: volumes is NULL, creating from getVolumes()",
      "   mod_lipid_import_server: volumes type = character, length = 2",
      "   mod_lipid_import_server: volumes names = home, data"
    )
  )
})

test_that("prepareLipidImportShinyFileVolumes keeps existing and empty volumes logging stable", {
  emitted_messages <- character()

  existing_result <- prepareLipidImportShinyFileVolumes(
    volumes = c(cache = "/var/cache"),
    getVolumes = function() {
      stop("getVolumes should not be used when volumes are already provided")
    },
    emitMessage = function(messageText) {
      emitted_messages <<- c(emitted_messages, messageText)
      invisible(NULL)
    }
  )

  empty_result <- prepareLipidImportShinyFileVolumes(
    volumes = structure(character(), names = character()),
    emitMessage = function(messageText) {
      emitted_messages <<- c(emitted_messages, messageText)
      invisible(NULL)
    }
  )

  expect_identical(existing_result, c(cache = "/var/cache"))
  expect_identical(empty_result, structure(character(), names = character()))
  expect_identical(
    emitted_messages,
    c(
      "   mod_lipid_import_server: volumes type = character, length = 1",
      "   mod_lipid_import_server: volumes names = cache",
      "   mod_lipid_import_server: volumes type = character, length = 0",
      "   mod_lipid_import_server: WARNING - volumes is empty!"
    )
  )
})

test_that("registerLipidImportFileLoadedOutput preserves the file-loaded reactive and outputOptions contract", {
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(lipid_id = "L1")
  reactive_values <- list()
  output_options_calls <- list()

  result <- registerLipidImportFileLoadedOutput(
    output = output,
    localData = local_data,
    reactiveFn = function(expr) {
      reactive_values <<- c(reactive_values, list(eval.parent(substitute(expr))))
      reactive_values[[length(reactive_values)]]
    },
    outputOptionsFn = function(output, name, suspendWhenHidden) {
      output_options_calls <<- c(output_options_calls, list(list(
        output = output,
        name = name,
        suspendWhenHidden = suspendWhenHidden
      )))
      invisible(NULL)
    }
  )

  expect_identical(result, output)
  expect_identical(output$file_loaded, TRUE)
  expect_identical(reactive_values, list(TRUE))
  expect_length(output_options_calls, 1L)
  expect_identical(output_options_calls[[1]]$output, output)
  expect_identical(output_options_calls[[1]]$name, "file_loaded")
  expect_false(output_options_calls[[1]]$suspendWhenHidden)
})

test_that("registerLipidImportFileLoadedOutput keeps the unloaded state stable", {
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- NULL

  registerLipidImportFileLoadedOutput(
    output = output,
    localData = local_data,
    reactiveFn = function(expr) eval.parent(substitute(expr)),
    outputOptionsFn = function(...) invisible(NULL)
  )

  expect_false(output$file_loaded)
})

test_that("registerLipidImportValidationSummaryOutput preserves forwarded reactive payload", {
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 20)
  )
  captured <- NULL
  lipid_id_calls <- 0L
  sample_column_calls <- 0L

  result <- registerLipidImportValidationSummaryOutput(
    output = output,
    localData = local_data,
    getLipidIdCol = function() {
      lipid_id_calls <<- lipid_id_calls + 1L
      "resolved-lipid-id"
    },
    getSampleColumns = function() {
      sample_column_calls <<- sample_column_calls + 1L
      c("sample_a", "sample_b")
    },
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "validation-summary-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$validation_summary, "validation-summary-ui")
  expect_identical(captured$assay1Data, local_data$assay1_data)
  expect_identical(captured$lipidIdCol, "resolved-lipid-id")
  expect_identical(captured$sampleColumns, c("sample_a", "sample_b"))
  expect_identical(lipid_id_calls, 1L)
  expect_identical(sample_column_calls, 1L)
})

test_that("registerLipidImportStatusOutput preserves forwarded workflow payload", {
  output <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "complete")
  captured <- NULL

  result <- registerLipidImportStatusOutput(
    output = output,
    workflowData = workflow_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "import-status-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$import_status, "import-status-ui")
  expect_identical(captured$workflowData, workflow_data)
})

test_that("registerLipidImportFormatDetectionStatusOutput preserves forwarded format payload", {
  output <- new.env(parent = emptyenv())
  local_data <- list(
    detected_format = "lipidsearch",
    format_confidence = 0.83
  )
  captured <- NULL

  result <- registerLipidImportFormatDetectionStatusOutput(
    output = output,
    localData = local_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "format-detection-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$format_detection_status, "format-detection-ui")
  expect_identical(captured$detectedFormat, "lipidsearch")
  expect_identical(captured$formatConfidence, 0.83)
})

test_that("registerLipidImportLipidIdStatusOutput preserves forwarded dropdown lipid-ID payload", {
  output <- new.env(parent = emptyenv())
  input <- list(lipid_id_col = "lipid_id")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 20)
  )
  captured <- NULL

  result <- registerLipidImportLipidIdStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "lipid-id-status-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$lipid_id_status, "lipid-id-status-ui")
  expect_identical(captured$assay1Data, local_data$assay1_data)
  expect_identical(captured$lipidIdCol, "lipid_id")
})

test_that("registerLipidImportAnnotationStatusOutput preserves forwarded dropdown annotation payload", {
  output <- new.env(parent = emptyenv())
  input <- list(annotation_col = "annotation")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B")
  )
  captured <- NULL

  result <- registerLipidImportAnnotationStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "annotation-status-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$annotation_status, "annotation-status-ui")
  expect_identical(captured$assay1Data, local_data$assay1_data)
  expect_identical(captured$annotationCol, "annotation")
})

test_that("registerLipidImportSampleColumnsDisplayOutput preserves forwarded sample-column payload", {
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_import_result <- list(
    sample_columns = c("Sample_A", "Sample_B", "Sample_C")
  )
  captured <- NULL

  result <- registerLipidImportSampleColumnsDisplayOutput(
    output = output,
    localData = local_data,
    renderText = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "sample-columns-preview"
    }
  )

  expect_identical(result, output)
  expect_identical(output$sample_columns_display, "sample-columns-preview")
  expect_identical(captured$assay1ImportResult, local_data$assay1_import_result)
})

test_that("registerLipidImportAvailableColumnsDisplayOutput preserves forwarded header payload", {
  output <- new.env(parent = emptyenv())
  local_data <- new.env(parent = emptyenv())
  local_data$all_headers <- c("LipidIon", "Class", "Sample_A", "Sample_B")
  captured <- NULL

  result <- registerLipidImportAvailableColumnsDisplayOutput(
    output = output,
    localData = local_data,
    renderText = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "available-columns-preview"
    }
  )

  expect_identical(result, output)
  expect_identical(output$available_columns_display, "available-columns-preview")
  expect_identical(captured$allHeaders, local_data$all_headers)
})

test_that("registerLipidImportCustomLipidIdStatusOutput preserves forwarded custom lipid-ID payload", {
  output <- new.env(parent = emptyenv())
  input <- list(lipid_id_col_custom = "lipid_name")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(
    lipid_name = c("L1", "L2"),
    sample_a = c(10, 20)
  )
  captured <- NULL

  result <- registerLipidImportCustomLipidIdStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "custom-lipid-id-status-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$lipid_id_status_custom, "custom-lipid-id-status-ui")
  expect_identical(captured$assay1Data, local_data$assay1_data)
  expect_identical(captured$lipidIdColCustom, "lipid_name")
})

test_that("registerLipidImportCustomAnnotationStatusOutput preserves forwarded custom annotation payload", {
  output <- new.env(parent = emptyenv())
  input <- list(annotation_col_custom = "lipid_class")
  local_data <- new.env(parent = emptyenv())
  local_data$assay1_data <- data.frame(
    lipid_name = c("L1", "L2"),
    lipid_class = c("PC", "TG")
  )
  captured <- NULL

  result <- registerLipidImportCustomAnnotationStatusOutput(
    output = output,
    input = input,
    localData = local_data,
    renderUi = function(expr) eval.parent(substitute(expr)),
    handleRender = function(...) {
      captured <<- list(...)
      "custom-annotation-status-ui"
    }
  )

  expect_identical(result, output)
  expect_identical(output$annotation_status_custom, "custom-annotation-status-ui")
  expect_identical(captured$assay1Data, local_data$assay1_data)
  expect_identical(captured$annotationColCustom, "lipid_class")
})

test_that("registerLipidImportModuleOutputs preserves bundled output registration payloads", {
  output <- new.env(parent = emptyenv())
  input <- list(vendor_format = "custom")
  local_data <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  get_lipid_id_col <- function() "resolved-lipid-id"
  get_annotation_col <- function() "resolved-annotation"
  get_sample_columns <- function() c("sample_a", "sample_b")
  captured <- list()
  call_order <- character()

  result <- registerLipidImportModuleOutputs(
    output = output,
    input = input,
    localData = local_data,
    workflowData = workflow_data,
    getLipidIdCol = get_lipid_id_col,
    getSampleColumns = get_sample_columns,
    registerFileLoadedOutput = function(...) {
      call_order <<- c(call_order, "file_loaded")
      captured$file_loaded <<- list(...)
      output
    },
    registerFormatDetectionStatusOutput = function(...) {
      call_order <<- c(call_order, "format_detection")
      captured$format_detection <<- list(...)
      output
    },
    registerLipidIdStatusOutput = function(...) {
      call_order <<- c(call_order, "lipid_id_status")
      captured$lipid_id_status <<- list(...)
      output
    },
    registerAnnotationStatusOutput = function(...) {
      call_order <<- c(call_order, "annotation_status")
      captured$annotation_status <<- list(...)
      output
    },
    registerSampleColumnsDisplayOutput = function(...) {
      call_order <<- c(call_order, "sample_columns_display")
      captured$sample_columns_display <<- list(...)
      output
    },
    registerAvailableColumnsDisplayOutput = function(...) {
      call_order <<- c(call_order, "available_columns_display")
      captured$available_columns_display <<- list(...)
      output
    },
    registerCustomLipidIdStatusOutput = function(...) {
      call_order <<- c(call_order, "custom_lipid_id_status")
      captured$custom_lipid_id_status <<- list(...)
      output
    },
    registerCustomAnnotationStatusOutput = function(...) {
      call_order <<- c(call_order, "custom_annotation_status")
      captured$custom_annotation_status <<- list(...)
      output
    },
    registerValidationSummaryOutput = function(...) {
      call_order <<- c(call_order, "validation_summary")
      captured$validation_summary <<- list(...)
      output
    },
    registerStatusOutput = function(...) {
      call_order <<- c(call_order, "status")
      captured$status <<- list(...)
      output
    }
  )

  expect_identical(result, output)
  expect_identical(
    call_order,
    c(
      "file_loaded",
      "format_detection",
      "lipid_id_status",
      "annotation_status",
      "sample_columns_display",
      "available_columns_display",
      "custom_lipid_id_status",
      "custom_annotation_status",
      "validation_summary",
      "status"
    )
  )
  expect_identical(captured$file_loaded$output, output)
  expect_identical(captured$file_loaded$localData, local_data)
  expect_identical(captured$format_detection$output, output)
  expect_identical(captured$format_detection$localData, local_data)
  expect_identical(captured$lipid_id_status$input, input)
  expect_identical(captured$lipid_id_status$localData, local_data)
  expect_identical(captured$annotation_status$input, input)
  expect_identical(captured$annotation_status$localData, local_data)
  expect_identical(captured$sample_columns_display$output, output)
  expect_identical(captured$sample_columns_display$localData, local_data)
  expect_identical(captured$available_columns_display$output, output)
  expect_identical(captured$available_columns_display$localData, local_data)
  expect_identical(captured$custom_lipid_id_status$input, input)
  expect_identical(captured$custom_lipid_id_status$localData, local_data)
  expect_identical(captured$custom_annotation_status$input, input)
  expect_identical(captured$custom_annotation_status$localData, local_data)
  expect_identical(captured$validation_summary$output, output)
  expect_identical(captured$validation_summary$localData, local_data)
  expect_identical(captured$validation_summary$getLipidIdCol, get_lipid_id_col)
  expect_identical(captured$validation_summary$getSampleColumns, get_sample_columns)
  expect_identical(captured$status$output, output)
  expect_identical(captured$status$workflowData, workflow_data)
})

test_that("initializeLipidImportLocalData preserves the seeded reactive-values payload", {
  captured <- NULL
  local_data <- new.env(parent = emptyenv())

  result <- initializeLipidImportLocalData(
    reactiveValues = function(...) {
      captured <<- list(...)
      local_data
    }
  )

  expect_identical(result, local_data)
  expect_named(
    captured,
    c(
      "assay1_file",
      "assay1_data",
      "assay1_import_result",
      "assay2_file",
      "assay2_data",
      "assay2_import_result",
      "detected_format",
      "format_confidence",
      "all_headers"
    ),
    ignore.order = FALSE
  )
  expect_true(all(vapply(captured, is.null, logical(1))))
})

test_that("emitLipidImportModuleServerEntryDiagnostics preserves entry-phase and module-phase logging", {
  messages <- character()

  entry_result <- emitLipidImportModuleServerEntryDiagnostics(
    volumesIsNull = TRUE,
    emitMessage = function(text) {
      messages <<- c(messages, text)
      invisible(NULL)
    }
  )
  module_result <- emitLipidImportModuleServerEntryDiagnostics(
    insideModuleServer = TRUE,
    emitMessage = function(text) {
      messages <<- c(messages, text)
      invisible(NULL)
    }
  )

  expect_null(entry_result)
  expect_null(module_result)
  expect_identical(
    messages,
    c(
      "--- Entering mod_lipid_import_server ---",
      "   mod_lipid_import_server: volumes param is NULL = TRUE",
      "   mod_lipid_import_server: Inside moduleServer function"
    )
  )
})

test_that("probeLipidImportShinyFilesAvailability preserves the availability probe and diagnostic log", {
  captured <- NULL
  messages <- character()

  result <- probeLipidImportShinyFilesAvailability(
    requireNamespaceFn = function(package, quietly = FALSE) {
      captured <<- list(package = package, quietly = quietly)
      FALSE
    },
    emitMessage = function(text) {
      messages <<- c(messages, text)
      invisible(NULL)
    }
  )

  expect_false(result)
  expect_identical(
    captured,
    list(package = "shinyFiles", quietly = TRUE)
  )
  expect_identical(
    messages,
    "   mod_lipid_import_server: shinyFiles available = FALSE"
  )
})

test_that("handleLipidImportFormatDetectionStatusRender preserves the detected-format payload", {
  status_args <- NULL

  result <- handleLipidImportFormatDetectionStatusRender(
    detectedFormat = "msdial",
    formatConfidence = 0.72,
    buildStatus = function(...) {
      status_args <<- list(...)
      "format-detection-status"
    }
  )

  expect_identical(result, "format-detection-status")
  expect_identical(status_args$detectedFormat, "msdial")
  expect_identical(status_args$formatConfidence, 0.72)
})

test_that("handleLipidImportFormatDetectionStatusRender keeps missing format guarded by req", {
  status_called <- FALSE

  expect_error(
    handleLipidImportFormatDetectionStatusRender(
      detectedFormat = NULL,
      formatConfidence = 0.72,
      buildStatus = function(...) {
        status_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(status_called)
})

test_that("handleLipidImportAvailableColumnsDisplayRender preserves the available-header preview payload", {
  preview_args <- NULL

  result <- handleLipidImportAvailableColumnsDisplayRender(
    allHeaders = c("lipid_id", "annotation", "sample_a"),
    formatPreviewText = function(...) {
      preview_args <<- list(...)
      "available-columns-preview"
    }
  )

  expect_identical(result, "available-columns-preview")
  expect_identical(preview_args$columnNames, c("lipid_id", "annotation", "sample_a"))
})

test_that("handleLipidImportAvailableColumnsDisplayRender keeps missing prerequisites guarded by req", {
  preview_called <- FALSE

  expect_error(
    handleLipidImportAvailableColumnsDisplayRender(
      allHeaders = NULL,
      formatPreviewText = function(...) {
        preview_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(preview_called)
})

test_that("resolveLipidImportLipidIdColumn preserves the dropdown passthrough without assay data", {
  expect_identical(
    resolveLipidImportLipidIdColumn(
      assay1Data = NULL,
      vendorFormat = "msdial",
      lipidIdCol = "lipid_id",
      lipidIdColCustom = "ignored"
    ),
    "lipid_id"
  )
})

test_that("resolveLipidImportLipidIdColumn preserves the forwarded lipid-ID payload", {
  assay_data <- data.frame(
    `Lipid Name` = c("PC 34:1", "TG 52:2"),
    check.names = FALSE
  )
  resolver_args <- NULL

  result <- resolveLipidImportLipidIdColumn(
    assay1Data = assay_data,
    vendorFormat = "custom",
    lipidIdCol = "unused_dropdown",
    lipidIdColCustom = "lipid name",
    resolveEffectiveColumn = function(...) {
      resolver_args <<- list(...)
      "resolved-lipid-id-column"
    }
  )

  expect_identical(result, "resolved-lipid-id-column")
  expect_identical(resolver_args$assayData, assay_data)
  expect_identical(resolver_args$vendorFormat, "custom")
  expect_identical(resolver_args$selectedColumn, "unused_dropdown")
  expect_identical(resolver_args$customColumn, "lipid name")
})

test_that("resolveLipidImportAnnotationColumn preserves the dropdown passthrough without assay data", {
  expect_identical(
    resolveLipidImportAnnotationColumn(
      assay1Data = NULL,
      vendorFormat = "msdial",
      annotationCol = "annotation",
      annotationColCustom = "ignored"
    ),
    "annotation"
  )
})

test_that("resolveLipidImportAnnotationColumn preserves the forwarded annotation payload", {
  assay_data <- data.frame(
    `Lipid Class` = c("PC", "TG"),
    check.names = FALSE
  )
  resolver_args <- NULL

  result <- resolveLipidImportAnnotationColumn(
    assay1Data = assay_data,
    vendorFormat = "custom",
    annotationCol = "unused_dropdown",
    annotationColCustom = "lipid class",
    resolveEffectiveColumn = function(...) {
      resolver_args <<- list(...)
      "resolved-annotation-column"
    }
  )

  expect_identical(result, "resolved-annotation-column")
  expect_identical(resolver_args$assayData, assay_data)
  expect_identical(resolver_args$vendorFormat, "custom")
  expect_identical(resolver_args$selectedColumn, "unused_dropdown")
  expect_identical(resolver_args$customColumn, "lipid class")
})

test_that("resolveLipidImportSelectedSampleColumns keeps missing assay data guarded by req", {
  resolver_called <- FALSE

  expect_error(
    resolveLipidImportSelectedSampleColumns(
      assay1Data = NULL,
      assay1ImportResult = list(sample_columns = "Sample_A"),
      vendorFormat = "msdial",
      sampleColsPattern = "",
      excludeNormalized = FALSE,
      resolveSampleColumns = function(...) {
        resolver_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(resolver_called)
})

test_that("resolveLipidImportSelectedSampleColumns preserves the forwarded sample-column payload", {
  assay_data <- data.frame(
    Sample_A = c(10, 11),
    Sample_B_norm = c(20, 21),
    check.names = FALSE
  )
  assay_import_result <- list(sample_columns = c("Sample_A", "Sample_B_norm"))
  resolver_args <- NULL

  result <- resolveLipidImportSelectedSampleColumns(
    assay1Data = assay_data,
    assay1ImportResult = assay_import_result,
    vendorFormat = "custom",
    sampleColsPattern = "^sample_",
    excludeNormalized = TRUE,
    resolveSampleColumns = function(...) {
      resolver_args <<- list(...)
      "resolved-sample-columns"
    }
  )

  expect_identical(result, "resolved-sample-columns")
  expect_identical(resolver_args$assayData, assay_data)
  expect_identical(resolver_args$assayImportResult, assay_import_result)
  expect_identical(resolver_args$vendorFormat, "custom")
  expect_identical(resolver_args$sampleColsPattern, "^sample_")
  expect_true(resolver_args$excludeNormalized)
})

test_that("buildLipidImportColumnSelectionReactives preserves reactive resolver wiring", {
  local_data <- shiny::reactiveValues(
    assay1_data = data.frame(
      lipid_id = c("L1", "L2"),
      annotation = c("A", "B"),
      Sample_A = c(10, 11),
      check.names = FALSE
    ),
    assay1_import_result = list(sample_columns = "Sample_A")
  )
  input <- list(
    vendor_format = "custom",
    lipid_id_col = "lipid_id",
    lipid_id_col_custom = "lipid id",
    annotation_col = "annotation",
    annotation_col_custom = "annotation custom",
    sample_cols_pattern = "^sample_",
    exclude_norm = TRUE
  )
  captured <- list()

  reactives <- buildLipidImportColumnSelectionReactives(
    input = input,
    localData = local_data,
    resolveLipidIdColumn = function(...) {
      captured$lipid_id <<- list(...)
      "resolved-lipid-id"
    },
    resolveAnnotationColumn = function(...) {
      captured$annotation <<- list(...)
      "resolved-annotation"
    },
    resolveSampleColumns = function(...) {
      captured$sample_columns <<- list(...)
      "resolved-sample-columns"
    }
  )

  expect_named(
    reactives,
    c("lipidIdCol", "annotationCol", "sampleColumns"),
    ignore.order = FALSE
  )
  expect_true(is.function(reactives$lipidIdCol))
  expect_true(is.function(reactives$annotationCol))
  expect_true(is.function(reactives$sampleColumns))
  expect_identical(shiny::isolate(reactives$lipidIdCol()), "resolved-lipid-id")
  expect_identical(shiny::isolate(reactives$annotationCol()), "resolved-annotation")
  expect_identical(shiny::isolate(reactives$sampleColumns()), "resolved-sample-columns")
  expect_identical(captured$lipid_id$assay1Data, shiny::isolate(local_data$assay1_data))
  expect_identical(captured$lipid_id$vendorFormat, "custom")
  expect_identical(captured$lipid_id$lipidIdCol, "lipid_id")
  expect_identical(captured$lipid_id$lipidIdColCustom, "lipid id")
  expect_identical(captured$annotation$assay1Data, shiny::isolate(local_data$assay1_data))
  expect_identical(captured$annotation$vendorFormat, "custom")
  expect_identical(captured$annotation$annotationCol, "annotation")
  expect_identical(captured$annotation$annotationColCustom, "annotation custom")
  expect_identical(captured$sample_columns$assay1Data, shiny::isolate(local_data$assay1_data))
  expect_identical(
    captured$sample_columns$assay1ImportResult,
    shiny::isolate(local_data$assay1_import_result)
  )
  expect_identical(captured$sample_columns$vendorFormat, "custom")
  expect_identical(captured$sample_columns$sampleColsPattern, "^sample_")
  expect_true(captured$sample_columns$excludeNormalized)
})

test_that("handleLipidImportLipidIdStatusRender preserves the dropdown lipid-ID validation payload", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )
  validation_args <- NULL

  result <- handleLipidImportLipidIdStatusRender(
    assay1Data = assay_data,
    lipidIdCol = "lipid_id",
    buildValidationStatus = function(...) {
      validation_args <<- list(...)
      "dropdown-lipid-id-status"
    }
  )

  expect_identical(result, "dropdown-lipid-id-status")
  expect_identical(validation_args$assayData, assay_data)
  expect_identical(validation_args$columnName, "lipid_id")
  expect_identical(validation_args$successMode, "unique_ids")
})

test_that("handleLipidImportLipidIdStatusRender keeps missing prerequisites guarded by req", {
  validation_called <- FALSE

  expect_error(
    handleLipidImportLipidIdStatusRender(
      assay1Data = NULL,
      lipidIdCol = "lipid_id",
      buildValidationStatus = function(...) {
        validation_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_error(
    handleLipidImportLipidIdStatusRender(
      assay1Data = data.frame(lipid_id = "L1"),
      lipidIdCol = NULL,
      buildValidationStatus = function(...) {
        validation_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(validation_called)
})

test_that("handleLipidImportAnnotationStatusRender preserves the dropdown annotation validation payload", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )
  validation_args <- NULL

  result <- handleLipidImportAnnotationStatusRender(
    assay1Data = assay_data,
    annotationCol = "annotation",
    buildValidationStatus = function(...) {
      validation_args <<- list(...)
      "dropdown-annotation-status"
    }
  )

  expect_identical(result, "dropdown-annotation-status")
  expect_identical(validation_args$assayData, assay_data)
  expect_identical(validation_args$columnName, "annotation")
  expect_identical(validation_args$successMode, "ok")
  expect_identical(validation_args$emptyMode, "optional")
})

test_that("handleLipidImportAnnotationStatusRender leaves optional empty state with the validation builder", {
  validation_args <- NULL

  result <- handleLipidImportAnnotationStatusRender(
    assay1Data = data.frame(lipid_id = "L1", check.names = FALSE),
    annotationCol = "",
    buildValidationStatus = function(...) {
      validation_args <<- list(...)
      "optional-annotation-status"
    }
  )

  expect_identical(result, "optional-annotation-status")
  expect_identical(validation_args$columnName, "")
  expect_identical(validation_args$emptyMode, "optional")
})

test_that("handleLipidImportAnnotationStatusRender keeps missing prerequisites guarded by req", {
  validation_called <- FALSE

  expect_error(
    handleLipidImportAnnotationStatusRender(
      assay1Data = NULL,
      annotationCol = "annotation",
      buildValidationStatus = function(...) {
        validation_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(validation_called)
})

test_that("handleLipidImportCustomLipidIdStatusRender preserves the custom lipid-ID validation payload", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )
  validation_args <- NULL

  result <- handleLipidImportCustomLipidIdStatusRender(
    assay1Data = assay_data,
    lipidIdColCustom = "lipid_id",
    buildValidationStatus = function(...) {
      validation_args <<- list(...)
      "custom-lipid-id-status"
    }
  )

  expect_identical(result, "custom-lipid-id-status")
  expect_identical(validation_args$assayData, assay_data)
  expect_identical(validation_args$columnName, "lipid_id")
  expect_identical(validation_args$successMode, "found_unique_ids")
  expect_identical(validation_args$emptyMode, "prompt")
  expect_true(validation_args$allowCaseInsensitive)
})

test_that("handleLipidImportCustomLipidIdStatusRender keeps missing prerequisites guarded by req", {
  validation_called <- FALSE

  expect_error(
    handleLipidImportCustomLipidIdStatusRender(
      assay1Data = NULL,
      lipidIdColCustom = "lipid_id",
      buildValidationStatus = function(...) {
        validation_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(validation_called)
})

test_that("handleLipidImportCustomAnnotationStatusRender preserves the custom annotation validation payload", {
  assay_data <- data.frame(
    lipid_id = c("L1", "L2"),
    annotation = c("A", "B"),
    check.names = FALSE
  )
  validation_args <- NULL

  result <- handleLipidImportCustomAnnotationStatusRender(
    assay1Data = assay_data,
    annotationColCustom = "annotation",
    buildValidationStatus = function(...) {
      validation_args <<- list(...)
      "custom-annotation-status"
    }
  )

  expect_identical(result, "custom-annotation-status")
  expect_identical(validation_args$assayData, assay_data)
  expect_identical(validation_args$columnName, "annotation")
  expect_identical(validation_args$successMode, "found")
  expect_identical(validation_args$emptyMode, "optional")
  expect_true(validation_args$allowCaseInsensitive)
})

test_that("handleLipidImportCustomAnnotationStatusRender keeps missing prerequisites guarded by req", {
  validation_called <- FALSE

  expect_error(
    handleLipidImportCustomAnnotationStatusRender(
      assay1Data = NULL,
      annotationColCustom = "annotation",
      buildValidationStatus = function(...) {
        validation_called <<- TRUE
        invisible(NULL)
      }
    ),
    class = "shiny.silent.error"
  )

  expect_false(validation_called)
})

test_that("handleLipidImportStatusRender preserves the import-status workflow payload", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "complete")
  workflow_data$processing_log <- list(
    setup_import = list(
      detected_format = "msdial",
      n_assays = 2L,
      n_samples = 8L
    )
  )
  captured_workflow <- NULL

  result <- handleLipidImportStatusRender(
    workflowData = workflow_data,
    buildStatusDisplay = function(workflowData) {
      captured_workflow <<- workflowData
      "rendered-import-status"
    }
  )

  expect_identical(result, "rendered-import-status")
  expect_identical(captured_workflow, workflow_data)
})

test_that("buildLipidImportValidationSummary preserves passing summary and warning rendering", {
  rendered <- htmltools::renderTags(
    buildLipidImportValidationSummary(list(
      valid = TRUE,
      summary = list(
        n_lipids = 12L,
        n_samples = 4L,
        pct_missing = 8.5
      ),
      warnings = c("Some sample columns contain zeros"),
      errors = character(0)
    ))
  )$html

  expect_match(rendered, "alert-success", fixed = TRUE)
  expect_match(rendered, "Validation Passed", fixed = TRUE)
  expect_match(rendered, "Lipids: 12", fixed = TRUE)
  expect_match(rendered, "Samples: 4", fixed = TRUE)
  expect_match(rendered, "Missing values: 8.5%", fixed = TRUE)
  expect_match(rendered, "alert-warning", fixed = TRUE)
  expect_match(rendered, "Warnings:", fixed = TRUE)
  expect_match(rendered, "Some sample columns contain zeros", fixed = TRUE)
})

test_that("buildLipidImportValidationSummary preserves failing error rendering", {
  rendered <- htmltools::renderTags(
    buildLipidImportValidationSummary(list(
      valid = FALSE,
      summary = NULL,
      warnings = character(0),
      errors = c("Lipid ID column not found", "At least one sample column is required")
    ))
  )$html

  expect_match(rendered, "alert-danger", fixed = TRUE)
  expect_match(rendered, "Validation Failed", fixed = TRUE)
  expect_match(rendered, "Lipid ID column not found", fixed = TRUE)
  expect_match(rendered, "At least one sample column is required", fixed = TRUE)
})

test_that("buildLipidImportStatusDisplay preserves the completed import banner", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "complete")
  workflow_data$processing_log <- list(setup_import = list(
    detected_format = "msdial",
    n_assays = 2L,
    n_samples = 5L
  ))

  rendered <- htmltools::renderTags(
    buildLipidImportStatusDisplay(workflow_data)
  )$html

  expect_match(rendered, "alert-success", fixed = TRUE)
  expect_match(rendered, "Import Complete", fixed = TRUE)
  expect_match(rendered, "Format: MSDIAL | Assays: 2 | Samples: 5", fixed = TRUE)
})

test_that("buildLipidImportStatusDisplay keeps incomplete import state hidden", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(setup_import = "pending")
  workflow_data$processing_log <- list(setup_import = list(
    detected_format = "lipidsearch",
    n_assays = 1L,
    n_samples = 3L
  ))

  expect_null(buildLipidImportStatusDisplay(workflow_data))
})

test_that("applyLipidImportResultToWorkflow stores lipid import state and initializes workflow type", {
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
      lipid_id = c("L1", "L2"),
      sample_a = c(10, 11)
    )
  )

  result <- applyLipidImportResultToWorkflow(
    workflowData = workflow_data,
    assayList = assay_list,
    dataFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = "sample_a",
    isPattern = "^ISTD_",
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(result, "lipidomics_standard")
  expect_identical(workflow_data$data_tbl, assay_list)
  expect_identical(workflow_data$data_format, "msdial")
  expect_identical(workflow_data$data_type, "lipid")
  expect_identical(
    workflow_data$column_mapping,
    list(
      lipid_id_col = "lipid_id",
      annotation_col = "annotation",
      sample_columns = "sample_a",
      is_pattern = "^ISTD_"
    )
  )
  expect_identical(workflow_type_calls, "lipidomics_standard")
  expect_identical(log_messages, "Workflow type set to: lipidomics_standard")
})

test_that("applyLipidImportResultToWorkflow preserves blank option fallbacks without a state manager", {
  workflow_data <- new.env(parent = emptyenv())
  assay_list <- list(
    LCMS_Neg = data.frame(
      lipid_id = c("L3", "L4"),
      sample_b = c(20, 21)
    )
  )
  log_messages <- character()

  result <- applyLipidImportResultToWorkflow(
    workflowData = workflow_data,
    assayList = assay_list,
    dataFormat = "custom",
    lipidIdCol = "lipid_id",
    annotationCol = "",
    sampleColumns = "sample_b",
    isPattern = "",
    logInfo = function(message) log_messages <<- c(log_messages, message)
  )

  expect_identical(result, "lipidomics_standard")
  expect_identical(workflow_data$data_tbl, assay_list)
  expect_identical(workflow_data$data_format, "custom")
  expect_identical(workflow_data$data_type, "lipid")
  expect_identical(
    workflow_data$column_mapping,
    list(
      lipid_id_col = "lipid_id",
      annotation_col = NULL,
      sample_columns = "sample_b",
      is_pattern = NA_character_
    )
  )
  expect_false(exists("state_manager", envir = workflow_data, inherits = FALSE))
  expect_identical(log_messages, character())
})

test_that("finalizeLipidImportSetupState records completion metrics and log output", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$processing_log <- list()
  workflow_data$tab_status <- list(other_tab = "incomplete")

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

  log_messages <- character()
  timestamp <- as.POSIXct("2026-04-15 11:20:00", tz = "UTC")

  finalizeLipidImportSetupState(
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
  expect_identical(workflow_data$processing_log$setup_import$assay_names, c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(workflow_data$processing_log$setup_import$detected_format, "msdial")
  expect_identical(as.integer(workflow_data$processing_log$setup_import$n_lipids), c(2L, 2L))
  expect_identical(workflow_data$processing_log$setup_import$n_samples, 2L)
  expect_identical(workflow_data$tab_status$setup_import, "complete")
  expect_identical(workflow_data$tab_status$other_tab, "incomplete")
  expect_identical(log_messages, "Lipidomics import complete: 2 assays, 5 total lipids")
})

test_that("finalizeLipidImportSetupState falls back to row counts when the lipid column is absent", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$processing_log <- list()
  workflow_data$tab_status <- list()

  assay_list <- list(
    GCMS = data.frame(
      feature = c("F1", "F2", "F3"),
      sample_only = c(1, 2, 3)
    )
  )

  finalizeLipidImportSetupState(
    workflowData = workflow_data,
    assayList = assay_list,
    detectedFormat = "custom",
    lipidIdCol = "lipid_id",
    sampleColumns = "sample_only",
    now = function() as.POSIXct("2026-04-15 11:25:00", tz = "UTC"),
    logInfo = function(...) NULL
  )

  expect_identical(as.integer(workflow_data$processing_log$setup_import$n_lipids), 3L)
  expect_identical(workflow_data$processing_log$setup_import$assay_names, "GCMS")
  expect_identical(workflow_data$processing_log$setup_import$detected_format, "custom")
  expect_identical(workflow_data$tab_status$setup_import, "complete")
})

test_that("sanitizeLipidImportSampleNames renames matched sample columns across assays", {
  assay_list <- list(
    LCMS_Pos = data.frame(
      lipid_id = c("L1", "L2"),
      `123-Sample!` = c(10, 11),
      `Control (A)` = c(20, 21),
      check.names = FALSE
    ),
    LCMS_Neg = data.frame(
      lipid_id = c("L3", "L4"),
      `123-Sample!` = c(30, 31),
      batch = c("b1", "b1"),
      check.names = FALSE
    )
  )

  log_messages <- character()
  notifications <- list()
  cleaned_columns <- c("x123_sample", "control_a")

  result <- sanitizeLipidImportSampleNames(
    assayList = assay_list,
    sampleColumns = c("123-Sample!", "Control (A)"),
    sanitizeNames = TRUE,
    makeCleanNames = function(x) cleaned_columns,
    logInfo = function(message) log_messages <<- c(log_messages, message),
    notify = function(message, type) {
      notifications <<- append(notifications, list(list(message = message, type = type)))
      invisible(NULL)
    }
  )

  expect_identical(result$sampleColumns, cleaned_columns)
  expect_identical(
    names(result$assayList$LCMS_Pos),
    c("lipid_id", "x123_sample", "control_a")
  )
  expect_identical(
    names(result$assayList$LCMS_Neg),
    c("lipid_id", "x123_sample", "batch")
  )
  expect_identical(
    log_messages,
    c(
      "Sanitizing sample names in lipidomics data...",
      "Sanitized 2 sample column names."
    )
  )
  expect_identical(
    notifications,
    list(list(
      message = "Sample names sanitized for R compatibility.",
      type = "message"
    ))
  )
})

test_that("sanitizeLipidImportSampleNames keeps assay data unchanged when sanitization is disabled", {
  assay_list <- list(
    LCMS_Pos = data.frame(
      lipid_id = c("L1", "L2"),
      `123-Sample!` = c(10, 11),
      check.names = FALSE
    )
  )

  log_messages <- character()
  notifications <- list()

  result <- sanitizeLipidImportSampleNames(
    assayList = assay_list,
    sampleColumns = "123-Sample!",
    sanitizeNames = FALSE,
    makeCleanNames = function(x) paste0("clean_", x),
    logInfo = function(message) log_messages <<- c(log_messages, message),
    notify = function(message, type) {
      notifications <<- append(notifications, list(list(message = message, type = type)))
      invisible(NULL)
    }
  )

  expect_identical(result$assayList, assay_list)
  expect_identical(result$sampleColumns, "123-Sample!")
  expect_identical(log_messages, character())
  expect_identical(notifications, list())
})

test_that("assembleLipidImportAssayList preserves the primary assay-only import path", {
  assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 11)
  )

  result <- assembleLipidImportAssayList(
    assay1Name = "LCMS_Pos",
    assay1Data = assay1_data,
    assay2File = NULL,
    assay2Name = ""
  )

  expect_identical(result, list(LCMS_Pos = assay1_data))
})

test_that("assembleLipidImportAssayList appends the optional second assay through the importer seam", {
  assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 11)
  )
  importer_calls <- character()
  assay2_data <- data.frame(
    lipid_id = c("L3", "L4"),
    sample_b = c(20, 21)
  )

  result <- assembleLipidImportAssayList(
    assay1Name = "LCMS_Pos",
    assay1Data = assay1_data,
    assay2File = "/tmp/lcms-neg.tsv",
    assay2Name = "LCMS_Neg",
    importSecondAssay = function(path) {
      importer_calls <<- c(importer_calls, path)
      list(data = assay2_data)
    }
  )

  expect_identical(importer_calls, "/tmp/lcms-neg.tsv")
  expect_identical(
    result,
    list(
      LCMS_Pos = assay1_data,
      LCMS_Neg = assay2_data
    )
  )
})

test_that("runLipidImportProcessing preserves notification lifecycle and downstream seam inputs", {
  workflow_data <- new.env(parent = emptyenv())
  assay1_data <- data.frame(
    lipid_id = c("L1", "L2"),
    sample_a = c(10, 11)
  )
  assay_list <- list(LCMS_Pos = assay1_data)
  sanitized_assays <- list(
    LCMS_Pos = data.frame(
      lipid_id = c("L1", "L2"),
      clean_sample_a = c(10, 11)
    )
  )
  notifications <- list()
  removed_notifications <- character()
  apply_args <- NULL
  finalize_args <- NULL

  result <- runLipidImportProcessing(
    workflowData = workflow_data,
    assay1Name = "LCMS_Pos",
    assay1Data = assay1_data,
    assay2File = NULL,
    assay2Name = "",
    vendorFormat = "custom",
    detectedFormat = "msdial",
    lipidIdCol = "lipid_id",
    annotationCol = "annotation",
    sampleColumns = "sample_a",
    isPattern = "^ISTD_",
    sanitizeNames = TRUE,
    assembleAssayList = function(...) assay_list,
    sanitizeSampleNames = function(...) {
      list(
        assayList = sanitized_assays,
        sampleColumns = "clean_sample_a"
      )
    },
    applyResultToWorkflow = function(
      workflowData,
      assayList,
      dataFormat,
      lipidIdCol,
      annotationCol,
      sampleColumns,
      isPattern
    ) {
      workflowData$data_format <- dataFormat
      apply_args <<- list(
        assayList = assayList,
        dataFormat = dataFormat,
        lipidIdCol = lipidIdCol,
        annotationCol = annotationCol,
        sampleColumns = sampleColumns,
        isPattern = isPattern
      )
      invisible("lipidomics_standard")
    },
    finalizeSetupState = function(
      workflowData,
      assayList,
      detectedFormat,
      lipidIdCol,
      sampleColumns
    ) {
      finalize_args <<- list(
        assayList = assayList,
        detectedFormat = detectedFormat,
        lipidIdCol = lipidIdCol,
        sampleColumns = sampleColumns
      )
      invisible(NULL)
    },
    notify = function(message, type = NULL, id = NULL, duration = NULL) {
      notifications <<- append(notifications, list(list(
        message = message,
        type = type,
        id = id,
        duration = duration
      )))
      invisible(NULL)
    },
    removeNotify = function(id) {
      removed_notifications <<- c(removed_notifications, id)
      invisible(NULL)
    },
    logError = function(...) stop("unexpected error path")
  )

  expect_identical(
    result,
    list(
      assayList = sanitized_assays,
      sampleColumns = "clean_sample_a"
    )
  )
  expect_identical(apply_args$assayList, sanitized_assays)
  expect_identical(apply_args$dataFormat, "custom")
  expect_identical(apply_args$lipidIdCol, "lipid_id")
  expect_identical(apply_args$annotationCol, "annotation")
  expect_identical(apply_args$sampleColumns, "clean_sample_a")
  expect_identical(apply_args$isPattern, "^ISTD_")
  expect_identical(finalize_args$assayList, sanitized_assays)
  expect_identical(finalize_args$detectedFormat, "custom")
  expect_identical(finalize_args$lipidIdCol, "lipid_id")
  expect_identical(finalize_args$sampleColumns, "clean_sample_a")
  expect_identical(removed_notifications, "lipid_import_working")
  expect_identical(
    notifications,
    list(
      list(
        message = "Processing imported data...",
        type = NULL,
        id = "lipid_import_working",
        duration = NULL
      ),
      list(
        message = "Data imported successfully!",
        type = "message",
        id = NULL,
        duration = NULL
      )
    )
  )
})

test_that("runLipidImportProcessing keeps the working notification cleanup on failures", {
  workflow_data <- new.env(parent = emptyenv())
  notifications <- list()
  removed_notifications <- character()
  error_messages <- character()
  apply_called <- FALSE
  finalize_called <- FALSE

  result <- runLipidImportProcessing(
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
    assembleAssayList = function(...) {
      stop("assembly exploded")
    },
    sanitizeSampleNames = function(...) stop("should not sanitize"),
    applyResultToWorkflow = function(...) {
      apply_called <<- TRUE
      invisible(NULL)
    },
    finalizeSetupState = function(...) {
      finalize_called <<- TRUE
      invisible(NULL)
    },
    notify = function(message, type = NULL, id = NULL, duration = NULL) {
      notifications <<- append(notifications, list(list(
        message = message,
        type = type,
        id = id,
        duration = duration
      )))
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

  expect_null(result)
  expect_false(apply_called)
  expect_false(finalize_called)
  expect_identical(error_messages, "Error processing import: assembly exploded")
  expect_identical(removed_notifications, "lipid_import_working")
  expect_identical(
    notifications,
    list(
      list(
        message = "Processing imported data...",
        type = NULL,
        id = "lipid_import_working",
        duration = NULL
      ),
      list(
        message = "Error: assembly exploded",
        type = "error",
        id = NULL,
        duration = 10
      )
    )
  )
})
