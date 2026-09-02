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
  logInfo = logger::log_info,
  vendorFormat = NULL,
  resolveFormatSupport = resolveWorkflowFormatSupport,
  workflowData = NULL,
  preIngressResolver = resolveLipidPreIngress,
  previewRows = 1000L
) {
  headers <- readHeaders(assay1File)

  if (length(headers) == 0) {
    stop("Could not read headers from file")
  }

  format_info <- detectFormat(
    headers = headers,
    filename = basename(assay1File)
  )

  requested_format <- vendorFormat %||% format_info$format
  decision <- resolveFormatSupport(
    omicType = "lipidomics",
    requestedFormat = requested_format,
    detectedFormat = format_info$format,
    detectionConfidence = format_info$confidence
  )
  observed_format <- format_info$format
  format_info$observed_format <- observed_format
  format_info$format <- decision$format
  format_info$support_status <- decision$support_status
  pre_ingress <- if (!is.null(workflowData) &&
      inherits(workflowData$workflow_context, "WorkflowContext")) {
    preIngressResolver(
      workflowData$workflow_context,
      assay1File,
      decision$format,
      headers
    )
  } else {
    list(
      status = "not_requested",
      defer_full_import = FALSE,
      outcome = NULL
    )
  }

  import_result <- switch(decision$format,
    "msdial" = importMsdial(assay1File),
    "lipidsearch" = {
      args <- list(assay1File)
      if (isTRUE(pre_ingress$defer_full_import) &&
          lipidImportFunctionAcceptsArg(importLipidSearch, "n_max")) {
        args$n_max <- as.integer(previewRows)
      }
      do.call(importLipidSearch, args)
    },
    "custom" = importMsdial(assay1File),
    workflowFormatSupportAbort(
      sprintf("No reader is registered for lipidomics format '%s'", decision$format),
      "multischolar_format_unsupported",
      "lipidomics",
      requested_format,
      observed_format,
      decision$support_status
    )
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
    deferred = isTRUE(pre_ingress$defer_full_import),
    preIngress = pre_ingress,
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

handleLipidImportDataPreviewLoad <- function(
  session,
  localData,
  loadPreview = loadLipidImportAssayPreview,
  applyPreview = applyLipidImportPreviewToModuleState,
  handleImportError = handleLipidImportDataImportError,
  vendorFormat = NULL,
  workflowData = NULL
) {
  shiny::req(localData$assay1_file)

  tryCatch(
    {
      previewArgs <- list(assay1File = localData$assay1_file)
      if (!is.null(vendorFormat) &&
          lipidImportFunctionAcceptsArg(loadPreview, "vendorFormat")) {
        previewArgs$vendorFormat <- vendorFormat
      }
      if (!is.null(workflowData) &&
          lipidImportFunctionAcceptsArg(loadPreview, "workflowData")) {
        previewArgs$workflowData <- workflowData
      }
      importPreview <- do.call(loadPreview, previewArgs)

      applyPreview(
        session = session,
        localData = localData,
        importPreview = importPreview
      )
    },
    error = function(e) {
      handleImportError(error = e)
    }
  )
}

buildLipidImportAssay1SelectedCallback <- function(
  session,
  localData,
  runPreviewLoad = handleLipidImportDataPreviewLoad,
  vendorFormat = NULL,
  workflowData = NULL
) {
  force(session)
  force(localData)
  force(vendorFormat)
  force(runPreviewLoad)

  function() {
    previewArgs <- list(
      session = session,
      localData = localData
    )
    if (!is.null(vendorFormat) &&
        lipidImportFunctionAcceptsArg(runPreviewLoad, "vendorFormat")) {
      previewArgs$vendorFormat <- if (is.function(vendorFormat)) {
        vendorFormat()
      } else {
        vendorFormat
      }
    }
    if (!is.null(workflowData) &&
        lipidImportFunctionAcceptsArg(runPreviewLoad, "workflowData")) {
      previewArgs$workflowData <- workflowData
    }
    do.call(runPreviewLoad, previewArgs)
  }
}

applyLipidImportPreviewToModuleState <- function(
  session,
  localData,
  importPreview,
  updateSelectInput = shiny::updateSelectInput,
  updateTextInput = shiny::updateTextInput
) {
  localData$all_headers <- importPreview$headers
  localData$detected_format <- importPreview$detectedFormat
  localData$format_confidence <- importPreview$formatConfidence
  localData$assay1_import_result <- importPreview$importResult
  localData$assay1_data <- importPreview$assayData
  localData$assay1_deferred <- isTRUE(importPreview$deferred)
  localData$preingress <- importPreview$preIngress

  updateSelectInput(
    session,
    "lipid_id_col",
    choices = importPreview$updates$lipidId$choices,
    selected = importPreview$updates$lipidId$selected
  )

  updateSelectInput(
    session,
    "annotation_col",
    choices = importPreview$updates$annotation$choices,
    selected = importPreview$updates$annotation$selected
  )

  if (!is.null(importPreview$updates$isPattern)) {
    updateTextInput(
      session,
      "is_pattern",
      value = importPreview$updates$isPattern
    )
  }

  invisible(importPreview)
}

handleLipidImportDataImportError <- function(
  error,
  logError = logger::log_error,
  notify = shiny::showNotification
) {
  errorMessage <- paste("Error importing data:", error$message)

  logError(errorMessage)
  notify(errorMessage, type = "error")

  invisible(FALSE)
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

registerLipidImportAssayFileSelectionObserver <- function(
  input,
  fileInputId,
  volumes,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  handleSelectionEvent = handleLipidImportAssayFileSelectionEvent,
  observeEvent = shiny::observeEvent
) {
  observeEvent(input[[fileInputId]], {
    handleSelectionEvent(
      fileInput = input[[fileInputId]],
      volumes = volumes,
      assignSelectedPath = assignSelectedPath,
      setRenderedPath = setRenderedPath,
      onSelected = onSelected
    )
  })
}

registerLipidImportShinyFileInputs <- function(
  input,
  output,
  session,
  localData,
  volumes,
  onAssay1Selected = NULL,
  setupAssayFileChooser = setupLipidImportAssayFileChooser,
  registerSelectionObserver = registerLipidImportAssayFileSelectionObserver
) {
  setupAssayFileChooser(
    input = input,
    fileInputId = "assay1_file",
    volumes = volumes,
    session = session
  )

  registerSelectionObserver(
    input = input,
    fileInputId = "assay1_file",
    volumes = volumes,
    assignSelectedPath = function(value) {
      localData$assay1_file <- value
    },
    setRenderedPath = function(renderedPath) {
      output$assay1_path <- renderedPath
    },
    onSelected = onAssay1Selected
  )

  setupAssayFileChooser(
    input = input,
    fileInputId = "assay2_file",
    volumes = volumes,
    session = session
  )

  registerSelectionObserver(
    input = input,
    fileInputId = "assay2_file",
    volumes = volumes,
    assignSelectedPath = function(value) {
      localData$assay2_file <- value
    },
    setRenderedPath = function(renderedPath) {
      output$assay2_path <- renderedPath
    }
  )

  invisible(volumes)
}

setupLipidImportShinyFileInputs <- function(
  useShinyFiles,
  input,
  output,
  session,
  localData,
  volumes,
  workflowData = NULL,
  prepareVolumes = prepareLipidImportShinyFileVolumes,
  buildAssay1Selected = buildLipidImportAssay1SelectedCallback,
  registerInputs = registerLipidImportShinyFileInputs
) {
  if (!isTRUE(useShinyFiles)) {
    return(invisible(volumes))
  }

  preparedVolumes <- prepareVolumes(volumes = volumes)

  assay1SelectedArgs <- list(
    session = session,
    localData = localData
  )
  if (lipidImportFunctionAcceptsArg(buildAssay1Selected, "vendorFormat")) {
    assay1SelectedArgs$vendorFormat <- function() input$vendor_format
  }
  if (!is.null(workflowData) && lipidImportFunctionAcceptsArg(
      buildAssay1Selected,
      "workflowData"
  )) {
    assay1SelectedArgs$workflowData <- workflowData
  }

  registerInputs(
    input = input,
    output = output,
    session = session,
    localData = localData,
    volumes = preparedVolumes,
    onAssay1Selected = do.call(buildAssay1Selected, assay1SelectedArgs)
  )

  invisible(preparedVolumes)
}

handleLipidImportStandardFileSelection <- function(
  fileInput,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  assignAssayPath = handleLipidImportSelectedAssayPathAssignment
) {
  if (is.null(fileInput) || is.null(fileInput$datapath) || length(fileInput$datapath) == 0) {
    return(invisible(NULL))
  }

  selectedPath <- as.character(fileInput$datapath[1])
  if (!nzchar(selectedPath)) {
    return(invisible(NULL))
  }

  assignAssayPath(
    selectedPath = selectedPath,
    assignSelectedPath = assignSelectedPath,
    setRenderedPath = setRenderedPath,
    onSelected = onSelected
  )
}

registerLipidImportStandardFileInputs <- function(
  input,
  output,
  localData,
  onAssay1Selected = NULL,
  observeEvent = shiny::observeEvent,
  handleStandardSelection = handleLipidImportStandardFileSelection
) {
  observeEvent(input$assay1_file_std, {
    handleStandardSelection(
      fileInput = input$assay1_file_std,
      assignSelectedPath = function(value) {
        localData$assay1_file <- value
      },
      setRenderedPath = function(renderedPath) {
        output$assay1_path <- renderedPath
      },
      onSelected = onAssay1Selected
    )
  })

  observeEvent(input$assay2_file_std, {
    handleStandardSelection(
      fileInput = input$assay2_file_std,
      assignSelectedPath = function(value) {
        localData$assay2_file <- value
      },
      setRenderedPath = function(renderedPath) {
        output$assay2_path <- renderedPath
      }
    )
  })

  invisible(NULL)
}

setupLipidImportStandardFileInputs <- function(
  useShinyFiles,
  input,
  output,
  session,
  localData,
  workflowData = NULL,
  buildAssay1Selected = buildLipidImportAssay1SelectedCallback,
  registerInputs = registerLipidImportStandardFileInputs
) {
  if (isTRUE(useShinyFiles)) {
    return(invisible(NULL))
  }

  assay1SelectedArgs <- list(
    session = session,
    localData = localData
  )
  if (lipidImportFunctionAcceptsArg(buildAssay1Selected, "vendorFormat")) {
    assay1SelectedArgs$vendorFormat <- function() input$vendor_format
  }
  if (!is.null(workflowData) && lipidImportFunctionAcceptsArg(
      buildAssay1Selected,
      "workflowData"
  )) {
    assay1SelectedArgs$workflowData <- workflowData
  }

  registerInputs(
    input = input,
    output = output,
    localData = localData,
    onAssay1Selected = do.call(buildAssay1Selected, assay1SelectedArgs)
  )

  invisible(NULL)
}

registerLipidImportProcessObserver <- function(
  input,
  workflowData,
  localData,
  getLipidIdCol,
  getAnnotationCol,
  getSampleColumns,
  experimentPaths = NULL,
  handleProcessRequest = handleLipidImportProcessRequest,
  observeEvent = shiny::observeEvent,
  releaseMemory = artifactReleaseTransientMemory
) {
  observeEvent(input$process_import, {
    request_args <- list(
      workflowData = workflowData,
      assay1Name = input$assay1_name,
      assay1Data = localData$assay1_data,
      assay1File = localData$assay1_file,
      assay2File = localData$assay2_file,
      assay2Name = input$assay2_name,
      vendorFormat = input$vendor_format,
      detectedFormat = localData$detected_format,
      formatConfidence = localData$format_confidence,
      lipidIdCol = getLipidIdCol(),
      annotationCol = getAnnotationCol(),
      sampleColumns = getSampleColumns(),
      isPattern = input$is_pattern,
      sanitizeNames = input$sanitize_names,
      experimentPaths = experimentPaths
    )
    if (lipidImportFunctionAcceptsArg(
        handleProcessRequest,
        "assay1Deferred"
    )) {
      request_args$assay1Deferred <- isTRUE(localData$assay1_deferred)
    }
    result <- do.call(handleProcessRequest, request_args)
    if (is.list(result) && identical(result$status, "success")) {
      localData$assay1_data <- NULL
      localData$assay1_import_result <- NULL
      localData$assay2_data <- NULL
      localData$assay2_import_result <- NULL
      localData$assay1_deferred <- FALSE
      if (isTRUE(request_args$assay1Deferred)) {
        result$assayList <- NULL
        releaseMemory(full = TRUE)
      }
    }
    invisible(result)
  })
}

handleLipidImportProcessRequest <- function(
  workflowData,
  assay1Name,
  assay1Data,
  assay1File = NULL,
  assay2File = NULL,
  assay2Name = NULL,
  vendorFormat,
  detectedFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  sanitizeNames,
  experimentPaths = NULL,
  processImport = runLipidImportProcessing,
  formatConfidence = NULL,
  assay1Deferred = FALSE
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  process_args <- list(
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
  if (lipidImportFunctionAcceptsArg(processImport, "formatConfidence")) {
    process_args$formatConfidence <- formatConfidence
  }
  if (lipidImportFunctionAcceptsArg(processImport, "assay1File")) {
    process_args$assay1File <- assay1File
  }
  if (lipidImportFunctionAcceptsArg(processImport, "assay1Deferred")) {
    process_args$assay1Deferred <- isTRUE(assay1Deferred)
  }
  if (lipidImportFunctionAcceptsArg(processImport, "experimentPaths")) {
    process_args$experimentPaths <- experimentPaths
  }

  do.call(processImport, process_args)
}
