metabImportFunctionAcceptsArg <- function(fn, arg) {
  formalNames <- names(formals(fn))
  arg %in% formalNames || "..." %in% formalNames
}

buildMetabImportWorkflowPayload <- function(
    assay1Name,
    assay1Data,
    assay1File = NULL,
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
    timestampFn = Sys.time,
    assay2PreparationFn = prepareMetabImportAssaySelectionState
) {
  assayNames <- c(assay1Name, if (!is.null(assay2File) && nzchar(assay2Name)) assay2Name else character(0))
  if (length(assayNames) == 0 || any(!nzchar(assayNames))) {
    stop("Metabolomics import requires non-empty assay names.")
  }
  if (anyDuplicated(assayNames)) {
    stop("Duplicate metabolomics assay names are not allowed: ", paste(unique(assayNames[duplicated(assayNames)]), collapse = ", "))
  }

  assayList <- list()
  assayList[[assay1Name]] <- assay1Data

  if (!is.null(assay2File) && nzchar(assay2Name)) {
    assay2Import <- if (missing(assay2Importer)) {
      assay2PreparationFn(
        assay1File = assay2File,
        vendorFormat = vendorFormat
      )$importResult
    } else {
      assay2Importer(assay2File)
    }
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
    ),
    sourceFiles = list(
      assay1 = assay1File,
      assay2 = if (!is.null(assay2File) && nzchar(assay2Name)) assay2File else NULL
    )
  )
}

metabImportValuesAreNumericLike <- function(values) {
  if (is.numeric(values) || is.integer(values)) {
    return(all(is.na(values) | is.finite(values)))
  }

  if (is.factor(values)) {
    values <- as.character(values)
  }

  if (!is.character(values)) {
    return(FALSE)
  }

  trimmed <- trimws(values)
  missingLike <- is.na(trimmed) | !nzchar(trimmed)
  coerced <- suppressWarnings(as.numeric(trimmed))
  all(missingLike | (!is.na(coerced) & is.finite(coerced)))
}

validateMetabImportWorkflowPayload <- function(workflowPayload,
                                               validateColumnMappingFn = validateMetabColumnMapping) {
  if (is.null(workflowPayload$assayList)) {
    return(list(valid = TRUE, errors = character(), warnings = "No assayList found; payload validation skipped."))
  }

  errors <- character()
  warnings <- character()
  assayList <- workflowPayload$assayList
  assayNames <- names(assayList)
  columnMapping <- workflowPayload$columnMapping
  metaboliteCol <- columnMapping$metabolite_id_col
  sampleCols <- columnMapping$sample_columns

  if (length(assayList) == 0) {
    errors <- c(errors, "No assay data provided.")
  }
  if (is.null(assayNames) || any(!nzchar(assayNames))) {
    errors <- c(errors, "All assays must have non-empty names.")
  }
  if (anyDuplicated(assayNames)) {
    errors <- c(errors, sprintf("Duplicate assay names are not allowed: %s", paste(unique(assayNames[duplicated(assayNames)]), collapse = ", ")))
  }
  if (any(grepl("[/\\\\\r\n\t]", assayNames))) {
    errors <- c(errors, "Assay names must not contain path separators or control characters.")
  }
  if (is.null(metaboliteCol) || !nzchar(metaboliteCol)) {
    errors <- c(errors, "Metabolite ID column is not specified.")
  }
  if (is.null(sampleCols) || length(sampleCols) == 0) {
    errors <- c(errors, "No sample columns specified.")
  }

  assaySummaries <- lapply(assayNames, function(assayName) {
    assayData <- assayList[[assayName]]
    assayErrors <- character()
    assayWarnings <- character()

    if (!is.data.frame(assayData)) {
      assayErrors <- c(assayErrors, sprintf("Assay '%s' is not a data frame.", assayName))
    } else {
      validation <- validateColumnMappingFn(
        data = assayData,
        metabolite_id_column = metaboliteCol,
        sample_columns = sampleCols
      )
      assayErrors <- c(assayErrors, validation$errors)
      assayWarnings <- c(assayWarnings, validation$warnings)

      presentSamples <- intersect(sampleCols, names(assayData))
      nonNumericSamples <- presentSamples[!vapply(assayData[presentSamples], metabImportValuesAreNumericLike, logical(1))]
      if (length(nonNumericSamples) > 0) {
        assayErrors <- c(
          assayErrors,
          sprintf(
            "Assay '%s' has non-numeric sample column(s): %s",
            assayName,
            paste(nonNumericSamples, collapse = ", ")
          )
        )
      }
    }

    list(
      assay_name = assayName,
      errors = assayErrors,
      warnings = assayWarnings
    )
  })

  for (summary in assaySummaries) {
    errors <- c(errors, summary$errors)
    warnings <- c(warnings, summary$warnings)
  }

  list(
    valid = length(errors) == 0,
    errors = unique(errors),
    warnings = unique(warnings),
    assay_summaries = assaySummaries
  )
}

coerceMetabImportWorkflowPayloadSamples <- function(workflowPayload) {
  if (is.null(workflowPayload$assayList) || is.null(workflowPayload$columnMapping$sample_columns)) {
    return(workflowPayload)
  }

  sampleCols <- workflowPayload$columnMapping$sample_columns
  workflowPayload$assayList <- lapply(workflowPayload$assayList, function(assayData) {
    presentSamples <- intersect(sampleCols, names(assayData))
    for (sampleCol in presentSamples) {
      if (!is.numeric(assayData[[sampleCol]])) {
        assayData[[sampleCol]] <- suppressWarnings(as.numeric(trimws(as.character(assayData[[sampleCol]]))))
      }
    }
    assayData
  })
  workflowPayload
}

writeMetabImportSourceArtifacts <- function(workflowPayload,
                                            experimentPaths = NULL,
                                            writeTableFn = utils::write.table,
                                            writeLinesFn = writeLines,
                                            writeJsonFn = jsonlite::write_json,
                                            fileCopyFn = file.copy,
                                            dirCreateFn = dir.create,
                                            dirExistsFn = dir.exists,
                                            fileExistsFn = file.exists) {
  sourceDir <- experimentPaths$source_dir
  if (is.null(sourceDir) || !nzchar(sourceDir) || !dirExistsFn(sourceDir)) {
    return(list(written = FALSE, reason = "source_dir unavailable"))
  }

  assayList <- workflowPayload$assayList
  assayNames <- names(assayList)
  artifactPaths <- list()

  for (assayName in assayNames) {
    assayPath <- file.path(sourceDir, paste0("data_cln_", assayName, ".tab"))
    writeTableFn(assayList[[assayName]], file = assayPath, sep = "\t", row.names = FALSE, quote = FALSE)
    artifactPaths[[paste0("data_cln_", assayName)]] <- assayPath
  }

  manifestPath <- file.path(sourceDir, "assay_manifest.txt")
  writeLinesFn(assayNames, manifestPath)
  artifactPaths$assay_manifest <- manifestPath

  columnMappingPath <- file.path(sourceDir, "column_mapping.json")
  writeJsonFn(workflowPayload$columnMapping, columnMappingPath, auto_unbox = TRUE, pretty = TRUE, null = "null")
  artifactPaths$column_mapping <- columnMappingPath

  sourceDirCopies <- file.path(sourceDir, "import_sources")
  dirCreateFn(sourceDirCopies, recursive = TRUE, showWarnings = FALSE)
  sourceFiles <- workflowPayload$sourceFiles
  copiedSources <- character()
  if (!is.null(sourceFiles)) {
    sourceFiles <- sourceFiles[!vapply(sourceFiles, is.null, logical(1))]
    for (sourceName in names(sourceFiles)) {
      sourceFile <- sourceFiles[[sourceName]]
      if (!is.null(sourceFile) && nzchar(sourceFile) && fileExistsFn(sourceFile)) {
        destination <- file.path(sourceDirCopies, paste0(sourceName, "_", basename(sourceFile)))
        fileCopyFn(sourceFile, destination, overwrite = TRUE)
        copiedSources <- c(copiedSources, destination)
      }
    }
  }

  summaryRows <- data.frame(
    assay_name = assayNames,
    feature_count = vapply(assayList, nrow, integer(1)),
    sample_count = length(workflowPayload$columnMapping$sample_columns),
    sample_columns = paste(workflowPayload$columnMapping$sample_columns, collapse = ","),
    stringsAsFactors = FALSE
  )
  summaryPath <- file.path(sourceDir, "metabolomics_import_summary.tsv")
  writeTableFn(summaryRows, file = summaryPath, sep = "\t", row.names = FALSE, quote = FALSE)
  artifactPaths$import_summary <- summaryPath

  list(
    written = TRUE,
    paths = artifactPaths,
    source_copies = copiedSources
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
    },
    vendorFormat = NULL,
    resolveFormatSupportFn = resolveWorkflowFormatSupport
) {
  headers <- readHeadersFn(assay1File)

  if (length(headers) == 0) {
    stop("Could not read headers from file")
  }

  formatInfo <- detectFormatFn(
    headers = headers,
    filename = basename(assay1File)
  )

  requestedFormat <- vendorFormat %||% formatInfo$format
  decision <- resolveFormatSupportFn(
    omicType = "metabolomics",
    requestedFormat = requestedFormat,
    detectedFormat = formatInfo$format,
    detectionConfidence = formatInfo$confidence
  )
  observedFormat <- formatInfo$format
  formatInfo$observed_format <- observedFormat
  formatInfo$format <- decision$format
  formatInfo$support_status <- decision$support_status

  importFn <- importers[[decision$format]]
  if (is.null(importFn) && identical(decision$format, "custom")) {
    importFn <- defaultImporter
  }
  if (is.null(importFn)) {
    workflowFormatSupportAbort(
      sprintf("No reader is registered for metabolomics format '%s'", decision$format),
      "multischolar_format_unsupported",
      "metabolomics",
      requestedFormat,
      observedFormat,
      decision$support_status
    )
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
    finalizeErrorFn = finalizeMetabImportAssaySelectionError,
    vendorFormat = NULL
) {
  reqFn(assay1File)

  tryCatch(
    {
      prepareArgs <- list(assay1File = assay1File)
      if (!is.null(vendorFormat) &&
          metabImportFunctionAcceptsArg(prepareImportStateFn, "vendorFormat")) {
        prepareArgs$vendorFormat <- vendorFormat
      }
      importState <- do.call(prepareImportStateFn, prepareArgs)

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
    assay1File = NULL,
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
    experimentPaths = NULL,
    reqFn = shiny::req,
    isTRUEFn = isTRUE,
    showNotificationFn = shiny::showNotification,
    logInfoFn = logger::log_info,
    buildWorkflowPayloadFn = buildMetabImportWorkflowPayload,
    validateWorkflowPayloadFn = validateMetabImportWorkflowPayload,
    coerceWorkflowPayloadFn = coerceMetabImportWorkflowPayloadSamples,
    writeImportArtifactsFn = writeMetabImportSourceArtifacts,
    applyWorkflowPayloadFn = applyMetabImportWorkflowPayload,
    persistArtifactFn = persistMetabImportArtifacts,
    finalizeFeedbackFn = finalizeMetabImportProcessingFeedback,
    sprintfFn = sprintf,
    lengthFn = length,
    formatConfidence = NULL,
    resolveFormatSupportFn = resolveWorkflowFormatSupport
) {
  formatDecision <- resolveFormatSupportFn(
    omicType = "metabolomics",
    requestedFormat = vendorFormat,
    detectedFormat = detectedFormat,
    detectionConfidence = formatConfidence
  )
  detectedFormat <- formatDecision$format

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

      workflowPayloadArgs <- list(
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
      if (metabImportFunctionAcceptsArg(buildWorkflowPayloadFn, "assay1File")) {
        workflowPayloadArgs$assay1File <- assay1File
      }
      workflowPayload <- do.call(buildWorkflowPayloadFn, workflowPayloadArgs)

      validationResult <- validateWorkflowPayloadFn(workflowPayload)
      if (!isTRUE(validationResult$valid)) {
        stop(sprintf(
          "Invalid metabolomics import: %s",
          paste(validationResult$errors, collapse = "; ")
        ))
      }

      workflowPayload <- coerceWorkflowPayloadFn(workflowPayload)
      artifactResult <- writeImportArtifactsFn(
        workflowPayload = workflowPayload,
        experimentPaths = experimentPaths
      )
      workflowPayload$processingLog$artifacts <- artifactResult

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
      artifactStageResult <- persistArtifactFn(
        workflow_data = workflowData,
        workflow_payload = workflowPayload
      )
      workflowData$processing_log$setup_import$artifact_stage <-
        artifactStageResult
      finalizeResult <- finalizeFeedbackFn(status = "success")

      invisible(list(
        status = "success",
        workflowPayload = workflowPayload,
        validationResult = validationResult,
        artifactResult = artifactResult,
        artifactStageResult = artifactStageResult,
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
