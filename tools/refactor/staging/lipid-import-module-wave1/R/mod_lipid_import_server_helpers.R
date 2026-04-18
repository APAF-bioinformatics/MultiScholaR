buildLipidImportFormatDetectionStatus <- function(detectedFormat, formatConfidence) {
  confidence_pct <- round(formatConfidence * 100)
  color_class <- if (confidence_pct >= 70) "success" else if (confidence_pct >= 40) "warning" else "danger"

  format_display <- switch(detectedFormat,
    "msdial" = "MS-DIAL",
    "lipidsearch" = "LipidSearch",
    "progenesis" = "Progenesis QI",
    "xcms" = "XCMS",
    "compound_discoverer" = "Compound Discoverer",
    "Unknown"
  )

  shiny::tags$div(
    class = paste("alert", paste0("alert-", color_class)),
    shiny::tags$strong("Detected format: "),
    format_display,
    shiny::tags$br(),
    shiny::tags$small(sprintf("Confidence: %d%%", confidence_pct))
  )
}

buildLipidImportColumnValidationStatus <- function(
  assayData,
  columnName,
  successMode = c("unique_ids", "found_unique_ids", "ok", "found"),
  emptyMode = c("optional", "prompt"),
  allowCaseInsensitive = FALSE
) {
  successMode <- match.arg(successMode)
  emptyMode <- match.arg(emptyMode)

  if (is.null(columnName) || !nzchar(columnName)) {
    if (identical(emptyMode, "prompt")) {
      return(shiny::tags$span(
        shiny::icon("question-circle", style = "color: gray;"),
        " Enter column name"
      ))
    }

    return(shiny::tags$span(
      shiny::icon("minus-circle", style = "color: gray;"),
      " Optional"
    ))
  }

  actual_column <- columnName
  if (allowCaseInsensitive) {
    headers <- names(assayData)
    header_index <- match(tolower(columnName), tolower(headers))
    actual_column <- if (is.na(header_index)) NULL else headers[[header_index]]
  } else if (!columnName %in% names(assayData)) {
    actual_column <- NULL
  }

  if (is.null(actual_column)) {
    return(shiny::tags$span(
      shiny::icon("times-circle", style = "color: red;"),
      " Column not found"
    ))
  }

  success_text <- switch(successMode,
    "unique_ids" = sprintf(" %d unique IDs", length(unique(assayData[[actual_column]]))),
    "found_unique_ids" = sprintf(" Found: %d unique IDs", length(unique(assayData[[actual_column]]))),
    "ok" = " OK",
    "found" = " Found"
  )

  shiny::tags$span(
    shiny::icon("check-circle", style = "color: green;"),
    success_text
  )
}

resolveLipidImportEffectiveColumn <- function(
  assayData,
  vendorFormat,
  selectedColumn,
  customColumn
) {
  if (!identical(vendorFormat, "custom")) {
    return(selectedColumn)
  }

  if (is.null(customColumn) || !nzchar(customColumn)) {
    return(customColumn)
  }

  headers <- names(assayData)
  header_index <- match(tolower(customColumn), tolower(headers))

  if (is.na(header_index)) {
    return(customColumn)
  }

  headers[[header_index]]
}

resolveLipidImportSampleColumns <- function(
  assayData,
  assayImportResult = NULL,
  vendorFormat,
  sampleColsPattern = "",
  excludeNormalized = FALSE,
  exclusionPattern = "_norm(ali[sz]ed|allized)?$|normalized"
) {
  if (identical(vendorFormat, "custom") &&
      !is.null(sampleColsPattern) &&
      nzchar(sampleColsPattern)) {
    all_cols <- names(assayData)
    matched <- all_cols[grepl(sampleColsPattern, all_cols, ignore.case = TRUE)]

    if (length(matched) > 0) {
      return(matched)
    }
  }

  if (!is.null(assayImportResult)) {
    sample_cols <- assayImportResult$sample_columns
  } else {
    sample_cols <- names(assayData)[sapply(assayData, is.numeric)]
  }

  if (isTRUE(excludeNormalized)) {
    sample_cols <- sample_cols[!grepl(exclusionPattern, sample_cols, ignore.case = TRUE)]
  }

  sample_cols
}

formatLipidImportColumnPreviewText <- function(columnNames, maxColumns = Inf) {
  if (length(columnNames) > maxColumns) {
    return(paste(
      paste(head(columnNames, maxColumns), collapse = ", "),
      sprintf("... and %d more", length(columnNames) - maxColumns)
    ))
  }

  paste(columnNames, collapse = ", ")
}

buildLipidImportValidationSummary <- function(validation) {
  if (validation$valid) {
    return(shiny::tagList(
      shiny::tags$div(
        class = "alert alert-success",
        shiny::icon("check-circle"),
        shiny::tags$strong(" Validation Passed")
      ),
      shiny::tags$ul(
        shiny::tags$li(sprintf("Lipids: %d", validation$summary$n_lipids)),
        shiny::tags$li(sprintf("Samples: %d", validation$summary$n_samples)),
        shiny::tags$li(sprintf("Missing values: %.1f%%", validation$summary$pct_missing))
      ),
      if (length(validation$warnings) > 0) {
        shiny::tags$div(
          class = "alert alert-warning",
          shiny::icon("exclamation-triangle"),
          " Warnings:",
          shiny::tags$ul(
            lapply(validation$warnings, shiny::tags$li)
          )
        )
      }
    ))
  }

  shiny::tags$div(
    class = "alert alert-danger",
    shiny::icon("times-circle"),
    shiny::tags$strong(" Validation Failed"),
    shiny::tags$ul(
      lapply(validation$errors, shiny::tags$li)
    )
  )
}

buildLipidImportStatusDisplay <- function(workflowData) {
  if (is.null(workflowData$tab_status$setup_import) ||
    !identical(workflowData$tab_status$setup_import, "complete")) {
    return(NULL)
  }

  log_info <- workflowData$processing_log$setup_import

  shiny::tags$div(
    class = "alert alert-success",
    shiny::icon("check-circle"),
    shiny::tags$strong(" Import Complete"),
    shiny::tags$br(),
    sprintf(
      "Format: %s | Assays: %d | Samples: %d",
      toupper(log_info$detected_format),
      log_info$n_assays,
      log_info$n_samples
    )
  )
}

applyLipidImportResultToWorkflow <- function(
  workflowData,
  assayList,
  dataFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  workflowType = "lipidomics_standard",
  logInfo = logger::log_info
) {
  workflowData$data_tbl <- assayList
  workflowData$data_format <- dataFormat
  workflowData$data_type <- "lipid"
  workflowData$column_mapping <- list(
    lipid_id_col = lipidIdCol,
    annotation_col = if (!is.null(annotationCol) && nzchar(annotationCol)) annotationCol else NULL,
    sample_columns = sampleColumns,
    is_pattern = if (!is.null(isPattern) && nzchar(isPattern)) isPattern else NA_character_
  )

  if (!is.null(workflowData$state_manager)) {
    workflowData$state_manager$setWorkflowType(workflowType)
    logInfo(sprintf("Workflow type set to: %s", workflowType))
  }

  invisible(workflowType)
}

finalizeLipidImportSetupState <- function(
  workflowData,
  assayList,
  detectedFormat,
  lipidIdCol,
  sampleColumns,
  now = Sys.time,
  logInfo = logger::log_info
) {
  workflowData$processing_log$setup_import <- list(
    timestamp = now(),
    n_assays = length(assayList),
    assay_names = names(assayList),
    detected_format = detectedFormat,
    n_lipids = sapply(assayList, function(assayData) {
      if (lipidIdCol %in% names(assayData)) {
        length(unique(assayData[[lipidIdCol]]))
      } else {
        nrow(assayData)
      }
    }),
    n_samples = length(sampleColumns)
  )

  updated_status <- workflowData$tab_status
  updated_status$setup_import <- "complete"
  workflowData$tab_status <- updated_status

  logInfo(sprintf(
    "Lipidomics import complete: %d assays, %d total lipids",
    length(assayList),
    sum(sapply(assayList, nrow))
  ))

  invisible(workflowData$processing_log$setup_import)
}

sanitizeLipidImportSampleNames <- function(
  assayList,
  sampleColumns,
  sanitizeNames,
  makeCleanNames = janitor::make_clean_names,
  logInfo = logger::log_info,
  notify = shiny::showNotification
) {
  if (!isTRUE(sanitizeNames)) {
    return(list(
      assayList = assayList,
      sampleColumns = sampleColumns
    ))
  }

  logInfo("Sanitizing sample names in lipidomics data...")

  original_sample_cols <- sampleColumns
  cleaned_sample_cols <- makeCleanNames(original_sample_cols)

  sanitized_assays <- purrr::map(assayList, function(assayData) {
    all_cols <- names(assayData)
    match_indices <- match(original_sample_cols, all_cols)
    valid_matches <- !is.na(match_indices)

    if (any(valid_matches)) {
      names(assayData)[match_indices[valid_matches]] <- cleaned_sample_cols[valid_matches]
    }

    assayData
  })

  logInfo(sprintf("Sanitized %d sample column names.", length(cleaned_sample_cols)))
  notify("Sample names sanitized for R compatibility.", type = "message")

  list(
    assayList = sanitized_assays,
    sampleColumns = cleaned_sample_cols
  )
}

assembleLipidImportAssayList <- function(
  assay1Name,
  assay1Data,
  assay2File = NULL,
  assay2Name = NULL,
  importSecondAssay = importLipidMSDIALData
) {
  assay_list <- list()
  assay_list[[assay1Name]] <- assay1Data

  if (!is.null(assay2File) && nzchar(assay2Name)) {
    assay2_import <- importSecondAssay(assay2File)
    assay_list[[assay2Name]] <- assay2_import$data
  }

  assay_list
}

runLipidImportProcessing <- function(
  workflowData,
  assay1Name,
  assay1Data,
  assay2File = NULL,
  assay2Name = NULL,
  vendorFormat,
  detectedFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  sanitizeNames,
  assembleAssayList = assembleLipidImportAssayList,
  sanitizeSampleNames = sanitizeLipidImportSampleNames,
  applyResultToWorkflow = applyLipidImportResultToWorkflow,
  finalizeSetupState = finalizeLipidImportSetupState,
  notify = shiny::showNotification,
  removeNotify = shiny::removeNotification,
  logError = logger::log_error
) {
  notify(
    "Processing imported data...",
    id = "lipid_import_working",
    duration = NULL
  )

  tryCatch(
    {
      assay_list <- assembleAssayList(
        assay1Name = assay1Name,
        assay1Data = assay1Data,
        assay2File = assay2File,
        assay2Name = assay2Name
      )

      sanitized_import <- sanitizeSampleNames(
        assayList = assay_list,
        sampleColumns = sampleColumns,
        sanitizeNames = sanitizeNames
      )
      assay_list <- sanitized_import$assayList
      sampleColumns <- sanitized_import$sampleColumns

      applyResultToWorkflow(
        workflowData = workflowData,
        assayList = assay_list,
        dataFormat = if (identical(vendorFormat, "custom")) "custom" else detectedFormat,
        lipidIdCol = lipidIdCol,
        annotationCol = annotationCol,
        sampleColumns = sampleColumns,
        isPattern = isPattern
      )

      finalizeSetupState(
        workflowData = workflowData,
        assayList = assay_list,
        detectedFormat = workflowData$data_format,
        lipidIdCol = lipidIdCol,
        sampleColumns = sampleColumns
      )

      removeNotify("lipid_import_working")
      notify("Data imported successfully!", type = "message")

      invisible(list(
        assayList = assay_list,
        sampleColumns = sampleColumns
      ))
    },
    error = function(e) {
      logError(paste("Error processing import:", e$message))
      removeNotify("lipid_import_working")
      notify(paste("Error:", e$message), type = "error", duration = 10)
      invisible(NULL)
    }
  )
}

loadLipidImportAssayPreview <- function(
  assay1File,
  readHeaders = function(path) {
    tryCatch(
      {
        con <- file(path, "r")
        on.exit(close(con), add = TRUE)

        first_line <- readLines(con, n = 1)
        raw_headers <- if (grepl(",", first_line)) {
          strsplit(first_line, ",")[[1]]
        } else {
          strsplit(first_line, "\t")[[1]]
        }

        gsub('^"|"$', "", raw_headers)
      },
      error = function(e) {
        character(0)
      }
    )
  },
  detectFormat = detectLipidomicsFormat,
  importMsdial = importLipidMSDIALData,
  importLipidSearch = importLipidSearchData,
  logInfo = logger::log_info
) {
  headers <- readHeaders(assay1File)

  if (length(headers) == 0) {
    stop("Could not read headers from file")
  }

  format_info <- detectFormat(
    headers = headers,
    filename = basename(assay1File)
  )

  import_result <- switch(format_info$format,
    "msdial" = importMsdial(assay1File),
    "lipidsearch" = importLipidSearch(assay1File),
    importMsdial(assay1File)
  )

  logInfo(sprintf(
    "Imported assay 1: %d rows, %d columns, format: %s",
    nrow(import_result$data),
    ncol(import_result$data),
    format_info$format
  ))

  list(
    headers = headers,
    detectedFormat = format_info$format,
    formatConfidence = format_info$confidence,
    importResult = import_result,
    assayData = import_result$data,
    updates = list(
      lipidId = list(
        choices = headers,
        selected = import_result$detected_columns$lipid_id
      ),
      annotation = list(
        choices = c("(None)" = "", headers),
        selected = import_result$detected_columns$annotation
      ),
      isPattern = if (!is.null(import_result$is_pattern) && !is.na(import_result$is_pattern)) {
        import_result$is_pattern
      } else {
        NULL
      }
    )
  )
}

handleLipidImportFileSelection <- function(
  fileInput,
  volumes,
  onPathSelected,
  parseFilePaths = shinyFiles::parseFilePaths,
  logError = logger::log_error
) {
  if (is.null(fileInput) || is.integer(fileInput)) {
    return(FALSE)
  }

  tryCatch(
    {
      file_info <- parseFilePaths(volumes, fileInput)
      if (nrow(file_info) == 0) {
        return(FALSE)
      }

      onPathSelected(as.character(file_info$datapath[1]))
      TRUE
    },
    error = function(e) {
      logError(paste("Error parsing file path:", e$message))
      FALSE
    }
  )
}

handleLipidImportProcessRequest <- function(
  workflowData,
  assay1Name,
  assay1Data,
  assay2File = NULL,
  assay2Name = NULL,
  vendorFormat,
  detectedFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  sanitizeNames,
  processImport = runLipidImportProcessing
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  processImport(
    workflowData = workflowData,
    assay1Name = assay1Name,
    assay1Data = assay1Data,
    assay2File = assay2File,
    assay2Name = assay2Name,
    vendorFormat = vendorFormat,
    detectedFormat = detectedFormat,
    lipidIdCol = lipidIdCol,
    annotationCol = annotationCol,
    sampleColumns = sampleColumns,
    isPattern = isPattern,
    sanitizeNames = sanitizeNames
  )
}

handleLipidImportValidationSummaryRender <- function(
  assay1Data,
  lipidIdCol,
  sampleColumns,
  validateMapping = validateLipidColumnMapping,
  buildSummary = buildLipidImportValidationSummary
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  validation <- validateMapping(
    data = assay1Data,
    lipid_id_column = lipidIdCol,
    sample_columns = sampleColumns
  )

  buildSummary(validation)
}

handleLipidImportSampleColumnsDisplayRender <- function(
  assay1ImportResult,
  formatPreviewText = formatLipidImportColumnPreviewText
) {
  shiny::req(assay1ImportResult)

  formatPreviewText(
    columnNames = assay1ImportResult$sample_columns,
    maxColumns = 10
  )
}

handleLipidImportAvailableColumnsDisplayRender <- function(
  allHeaders,
  formatPreviewText = formatLipidImportColumnPreviewText
) {
  shiny::req(allHeaders)

  formatPreviewText(columnNames = allHeaders)
}

handleLipidImportCustomLipidIdStatusRender <- function(
  assay1Data,
  lipidIdColCustom,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = lipidIdColCustom,
    successMode = "found_unique_ids",
    emptyMode = "prompt",
    allowCaseInsensitive = TRUE
  )
}

handleLipidImportCustomAnnotationStatusRender <- function(
  assay1Data,
  annotationColCustom,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = annotationColCustom,
    successMode = "found",
    emptyMode = "optional",
    allowCaseInsensitive = TRUE
  )
}

handleLipidImportStatusRender <- function(
  workflowData,
  buildStatusDisplay = buildLipidImportStatusDisplay
) {
  buildStatusDisplay(workflowData)
}

