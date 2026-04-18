buildMetabImportWorkflowPayload <- function(
    assay1Name,
    assay1Data,
    assay2File = NULL,
    assay2Name = NULL,
    vendorFormat,
    detectedFormat,
    metaboliteCol,
    annotationCol,
    sampleCols,
    sanitizeNames = FALSE,
    isPattern = "",
    assay2Importer = importMetabMSDIALData,
    cleanNamesFn = janitor::make_clean_names,
    mapAssaysFn = purrr::map,
    timestampFn = Sys.time
) {
  assayList <- list()
  assayList[[assay1Name]] <- assay1Data

  if (!is.null(assay2File) && nzchar(assay2Name)) {
    assay2Import <- assay2Importer(assay2File)
    assayList[[assay2Name]] <- assay2Import$data
  }

  finalSampleCols <- sampleCols

  if (isTRUE(sanitizeNames)) {
    originalSampleCols <- sampleCols
    cleanedSampleCols <- cleanNamesFn(originalSampleCols)

    assayList <- mapAssaysFn(assayList, function(assayDf) {
      allCols <- names(assayDf)
      matchIndices <- match(originalSampleCols, allCols)
      validMatches <- !is.na(matchIndices)

      if (any(validMatches)) {
        names(assayDf)[matchIndices[validMatches]] <- cleanedSampleCols[validMatches]
      }

      assayDf
    })

    finalSampleCols <- cleanedSampleCols
  }

  dataFormat <- if (identical(vendorFormat, "custom")) "custom" else detectedFormat

  list(
    assayList = assayList,
    sampleCols = finalSampleCols,
    sampleNamesSanitized = isTRUE(sanitizeNames),
    dataFormat = dataFormat,
    columnMapping = list(
      metabolite_id_col = metaboliteCol,
      annotation_col = if (!is.null(annotationCol) && nzchar(annotationCol)) annotationCol else NULL,
      sample_columns = finalSampleCols,
      is_pattern = if (nzchar(isPattern)) isPattern else NA_character_
    ),
    processingLog = list(
      timestamp = timestampFn(),
      n_assays = length(assayList),
      assay_names = names(assayList),
      detected_format = dataFormat,
      n_metabolites = sapply(assayList, function(assayDf) {
        if (metaboliteCol %in% names(assayDf)) {
          length(unique(assayDf[[metaboliteCol]]))
        } else {
          nrow(assayDf)
        }
      }),
      n_samples = length(finalSampleCols)
    )
  )
}

applyMetabImportWorkflowPayload <- function(
    workflowData,
    workflowPayload,
    workflowType = "metabolomics_standard",
    logInfo = logger::log_info
) {
  workflowData$data_tbl <- workflowPayload$assayList
  workflowData$data_format <- workflowPayload$dataFormat
  workflowData$data_type <- "metabolite"
  workflowData$column_mapping <- workflowPayload$columnMapping

  if (!is.null(workflowData$state_manager)) {
    workflowData$state_manager$setWorkflowType(workflowType)
    logInfo(sprintf("Workflow type set to: %s", workflowType))
  }

  workflowData$processing_log$setup_import <- workflowPayload$processingLog

  updatedStatus <- workflowData$tab_status
  updatedStatus$setup_import <- "complete"
  workflowData$tab_status <- updatedStatus

  assayRowCounts <- vapply(workflowPayload$assayList, nrow, integer(1))
  totalMetabolites <- sum(assayRowCounts)

  logInfo(sprintf(
    "Metabolomics import complete: %d assays, %d total metabolites",
    length(workflowPayload$assayList),
    totalMetabolites
  ))

  invisible(list(
    assayCount = length(workflowPayload$assayList),
    totalMetabolites = totalMetabolites,
    updatedStatus = updatedStatus
  ))
}

prepareMetabImportAssaySelectionState <- function(
    assay1File,
    detectFormatFn = detectMetabolomicsFormat,
    defaultImporter = importMetabMSDIALData,
    importers = list(msdial = importMetabMSDIALData),
    readHeadersFn = function(path) {
      tryCatch(
        {
          con <- file(path, "r")
          on.exit(close(con), add = TRUE)
          firstLine <- readLines(con, n = 1)

          rawHeaders <- if (grepl(",", firstLine)) {
            strsplit(firstLine, ",")[[1]]
          } else {
            strsplit(firstLine, "\t")[[1]]
          }

          gsub('^"|"$', "", rawHeaders)
        },
        error = function(e) {
          character(0)
        }
      )
    }
) {
  headers <- readHeadersFn(assay1File)

  if (length(headers) == 0) {
    stop("Could not read headers from file")
  }

  formatInfo <- detectFormatFn(
    headers = headers,
    filename = basename(assay1File)
  )

  importFn <- importers[[formatInfo$format]]
  if (is.null(importFn)) {
    importFn <- defaultImporter
  }

  importResult <- importFn(assay1File)

  list(
    headers = headers,
    formatInfo = formatInfo,
    importResult = importResult,
    metaboliteIdChoices = headers,
    selectedMetaboliteId = importResult$detected_columns$metabolite_id,
    annotationChoices = c("(None)" = "", headers),
    selectedAnnotation = importResult$detected_columns$annotation,
    isPattern = if (!is.null(importResult$is_pattern) && !is.na(importResult$is_pattern)) {
      importResult$is_pattern
    } else {
      NULL
    }
  )
}

applyMetabImportAssaySelectionState <- function(
    localData,
    importState,
    session,
    updateSelectInputFn = shiny::updateSelectInput,
    updateTextInputFn = shiny::updateTextInput,
    logInfoFn = logger::log_info
) {
  localData$all_headers <- importState$headers
  localData$detected_format <- importState$formatInfo$format
  localData$format_confidence <- importState$formatInfo$confidence
  localData$assay1_import_result <- importState$importResult
  localData$assay1_data <- importState$importResult$data

  updateSelectInputFn(
    session,
    "metabolite_id_col",
    choices = importState$metaboliteIdChoices,
    selected = importState$selectedMetaboliteId
  )

  updateSelectInputFn(
    session,
    "annotation_col",
    choices = importState$annotationChoices,
    selected = importState$selectedAnnotation
  )

  if (!is.null(importState$isPattern)) {
    updateTextInputFn(session, "is_pattern", value = importState$isPattern)
  }

  importedRows <- nrow(importState$importResult$data)
  importedCols <- ncol(importState$importResult$data)

  logInfoFn(sprintf(
    "Imported assay 1: %d rows, %d columns, format: %s",
    importedRows,
    importedCols,
    importState$formatInfo$format
  ))

  invisible(list(
    importedRows = importedRows,
    importedCols = importedCols,
    detectedFormat = importState$formatInfo$format,
    confidence = importState$formatInfo$confidence
  ))
}

finalizeMetabImportProcessingFeedback <- function(
    status = c("success", "error"),
    error = NULL,
    workingNotificationId = "metab_import_working",
    successMessage = "Data imported successfully!",
    removeNotificationFn = shiny::removeNotification,
    showNotificationFn = shiny::showNotification,
    logErrorFn = logger::log_error
) {
  status <- match.arg(status)
  removeNotificationFn(workingNotificationId)

  if (identical(status, "success")) {
    showNotificationFn(successMessage, type = "message")

    return(invisible(list(
      status = status,
      notificationId = workingNotificationId,
      message = successMessage
    )))
  }

  errorMessage <- if (inherits(error, "condition")) {
    conditionMessage(error)
  } else {
    as.character(error)[1]
  }

  logErrorFn(paste("Error processing import:", errorMessage))

  notificationMessage <- paste("Error:", errorMessage)
  showNotificationFn(notificationMessage, type = "error", duration = 10)

  invisible(list(
    status = status,
    notificationId = workingNotificationId,
    message = notificationMessage
  ))
}

finalizeMetabImportAssaySelectionError <- function(
    error,
    showNotificationFn = shiny::showNotification,
    logErrorFn = logger::log_error
) {
  errorMessage <- if (inherits(error, "condition")) {
    conditionMessage(error)
  } else {
    as.character(error)[1]
  }

  notificationMessage <- paste("Error importing data:", errorMessage)
  logErrorFn(notificationMessage)
  showNotificationFn(notificationMessage, type = "error")

  invisible(list(
    status = "error",
    message = notificationMessage
  ))
}

runMetabImportAssaySelection <- function(
    assay1File,
    localData,
    session,
    reqFn = shiny::req,
    prepareImportStateFn = prepareMetabImportAssaySelectionState,
    applyImportStateFn = applyMetabImportAssaySelectionState,
    finalizeErrorFn = finalizeMetabImportAssaySelectionError
) {
  reqFn(assay1File)

  tryCatch(
    {
      importState <- prepareImportStateFn(assay1File = assay1File)

      applyImportStateFn(
        localData = localData,
        importState = importState,
        session = session
      )
    },
    error = function(e) {
      finalizeErrorFn(error = e)
    }
  )
}

runMetabImportProcessing <- function(
    assay1Data,
    assay1Name,
    assay2File,
    assay2Name,
    vendorFormat,
    detectedFormat,
    sanitizeNames,
    isPattern,
    getMetaboliteIdColFn,
    getAnnotationColFn,
    getSampleColumnsFn,
    workflowData,
    reqFn = shiny::req,
    isTRUEFn = isTRUE,
    showNotificationFn = shiny::showNotification,
    logInfoFn = logger::log_info,
    buildWorkflowPayloadFn = buildMetabImportWorkflowPayload,
    applyWorkflowPayloadFn = applyMetabImportWorkflowPayload,
    finalizeFeedbackFn = finalizeMetabImportProcessingFeedback,
    sprintfFn = sprintf,
    lengthFn = length
) {
  reqFn(assay1Data)

  metaboliteCol <- getMetaboliteIdColFn()
  annotationCol <- getAnnotationColFn()
  sampleCols <- getSampleColumnsFn()

  reqFn(metaboliteCol)

  showNotificationFn(
    "Processing imported data...",
    id = "metab_import_working",
    duration = NULL
  )

  tryCatch(
    {
      if (isTRUEFn(sanitizeNames)) {
        logInfoFn("Sanitizing sample names in metabolomics data...")
      }

      workflowPayload <- buildWorkflowPayloadFn(
        assay1Name = assay1Name,
        assay1Data = assay1Data,
        assay2File = assay2File,
        assay2Name = assay2Name,
        vendorFormat = vendorFormat,
        detectedFormat = detectedFormat,
        metaboliteCol = metaboliteCol,
        annotationCol = annotationCol,
        sampleCols = sampleCols,
        sanitizeNames = isTRUEFn(sanitizeNames),
        isPattern = isPattern
      )

      if (isTRUEFn(workflowPayload$sampleNamesSanitized)) {
        logInfoFn(sprintfFn(
          "Sanitized %d sample column names.",
          lengthFn(workflowPayload$sampleCols)
        ))
        showNotificationFn(
          "Sample names sanitized for R compatibility.",
          type = "message"
        )
      }

      applyResult <- applyWorkflowPayloadFn(
        workflowData = workflowData,
        workflowPayload = workflowPayload
      )
      finalizeResult <- finalizeFeedbackFn(status = "success")

      invisible(list(
        status = "success",
        workflowPayload = workflowPayload,
        applyResult = applyResult,
        finalizeResult = finalizeResult
      ))
    },
    error = function(e) {
      finalizeResult <- finalizeFeedbackFn(status = "error", error = e)

      invisible(list(
        status = "error",
        error = e,
        finalizeResult = finalizeResult
      ))
    }
  )
}

