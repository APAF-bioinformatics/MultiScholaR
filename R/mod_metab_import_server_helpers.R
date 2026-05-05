metabImportServerFunctionAcceptsArg <- function(fn, arg) {
  formalNames <- names(formals(fn))
  arg %in% formalNames || "..." %in% formalNames
}

setupMetabImportShinyFiles <- function(
    input,
    output,
    session,
    volumes,
    localData,
    importDataFn,
    getVolumesFn = function() shinyFiles::getVolumes()(),
    shinyFileChooseFn = shinyFiles::shinyFileChoose,
    observeEventFn = shiny::observeEvent,
    handleAssayFileSelectionFn = handleMetabImportAssayFileSelection,
    logMessageFn = message,
    logErrorFn = logger::log_error,
    typeofFn = typeof,
    lengthFn = length,
    namesFn = names,
    pasteFn = paste
) {
  resolvedVolumes <- volumes

  if (is.null(resolvedVolumes)) {
    logMessageFn("   mod_metab_import_server: volumes is NULL, creating from getVolumes()")
    resolvedVolumes <- getVolumesFn()
  }

  logMessageFn(sprintf(
    "   mod_metab_import_server: volumes type = %s, length = %d",
    typeofFn(resolvedVolumes), lengthFn(resolvedVolumes)
  ))

  if (lengthFn(resolvedVolumes) > 0) {
    logMessageFn(sprintf(
      "   mod_metab_import_server: volumes names = %s",
      pasteFn(namesFn(resolvedVolumes), collapse = ", ")
    ))
  } else {
    logMessageFn("   mod_metab_import_server: WARNING - volumes is empty!")
  }

  shinyFileChooseFn(
    input,
    "assay1_file",
    roots = resolvedVolumes,
    session = session,
    filetypes = c("tsv", "tab", "txt", "csv", "xlsx", "parquet")
  )

  observeEventFn(input$assay1_file, {
    if (!is.null(input$assay1_file) && !is.integer(input$assay1_file)) {
      tryCatch(
        {
          handleAssayFileSelectionFn(
            selectedInput = input$assay1_file,
            volumes = resolvedVolumes,
            localData = localData,
            localField = "assay1_file",
            output = output,
            outputId = "assay1_path",
            onSelected = function(path) {
              importDataFn()
            }
          )
        },
        error = function(e) {
          logErrorFn(paste("Error parsing file path:", e$message))
        }
      )
    }
  })

  shinyFileChooseFn(
    input,
    "assay2_file",
    roots = resolvedVolumes,
    session = session,
    filetypes = c("tsv", "tab", "txt", "csv", "xlsx", "parquet")
  )

  observeEventFn(input$assay2_file, {
    if (!is.null(input$assay2_file) && !is.integer(input$assay2_file)) {
      tryCatch(
        {
          handleAssayFileSelectionFn(
            selectedInput = input$assay2_file,
            volumes = resolvedVolumes,
            localData = localData,
            localField = "assay2_file",
            output = output,
            outputId = "assay2_path"
          )
        },
        error = function(e) {
          logErrorFn(paste("Error parsing file path:", e$message))
        }
      )
    }
  })

  invisible(resolvedVolumes)
}

handleMetabImportStandardAssayFileSelection <- function(
    fileInput,
    localData,
    localField,
    output,
    outputId,
    onSelected = NULL,
    renderTextFn = shiny::renderText
) {
  if (is.null(fileInput) || is.null(fileInput$datapath) || length(fileInput$datapath) == 0) {
    return(invisible(NULL))
  }

  selectedPath <- as.character(fileInput$datapath[1])
  if (!nzchar(selectedPath)) {
    return(invisible(NULL))
  }

  localData[[localField]] <- selectedPath
  output[[outputId]] <- renderTextFn(selectedPath)

  if (is.function(onSelected)) {
    onSelected(selectedPath)
  }

  invisible(selectedPath)
}

setupMetabImportStandardFileInputs <- function(
    input,
    output,
    localData,
    importDataFn,
    observeEventFn = shiny::observeEvent,
    handleStandardSelectionFn = handleMetabImportStandardAssayFileSelection
) {
  observeEventFn(input$assay1_file_std, {
    handleStandardSelectionFn(
      fileInput = input$assay1_file_std,
      localData = localData,
      localField = "assay1_file",
      output = output,
      outputId = "assay1_path",
      onSelected = function(path) {
        importDataFn()
      }
    )
  })

  observeEventFn(input$assay2_file_std, {
    handleStandardSelectionFn(
      fileInput = input$assay2_file_std,
      localData = localData,
      localField = "assay2_file",
      output = output,
      outputId = "assay2_path"
    )
  })

  invisible(NULL)
}

setupMetabImportColumnAccessors <- function(
    input,
    localData,
    reactiveFn = shiny::reactive,
    reqFn = shiny::req,
    resolveColumnNameFn = resolveMetabImportColumnName,
    resolveSampleColumnsFn = resolveMetabImportSampleColumns
) {
  getMetaboliteIdCol <- reactiveFn({
    if (input$vendor_format == "custom") {
      return(resolveColumnNameFn(
        headers = names(localData$assay1_data),
        columnName = input$metabolite_id_col_custom
      ))
    }

    input$metabolite_id_col
  })

  getAnnotationCol <- reactiveFn({
    if (input$vendor_format == "custom") {
      return(resolveColumnNameFn(
        headers = names(localData$assay1_data),
        columnName = input$annotation_col_custom
      ))
    }

    input$annotation_col
  })

  getSampleColumns <- reactiveFn({
    reqFn(localData$assay1_data)

    resolveSampleColumnsFn(
      assayData = localData$assay1_data,
      vendorFormat = input$vendor_format,
      sampleColsPattern = input$sample_cols_pattern,
      importResult = localData$assay1_import_result
    )
  })

  list(
    getMetaboliteIdCol = getMetaboliteIdCol,
    getAnnotationCol = getAnnotationCol,
    getSampleColumns = getSampleColumns
  )
}

setupMetabImportAssaySelectionCallback <- function(
    localData,
    session,
    runAssaySelectionFn = runMetabImportAssaySelection
) {
  function() {
    runAssaySelectionFn(
      assay1File = localData$assay1_file,
      localData = localData,
      session = session
    )
  }
}

setupMetabImportProcessingObserver <- function(
    input,
    localData,
    columnAccessors,
    workflowData,
    experimentPaths = NULL,
    observeEventFn = shiny::observeEvent,
    runProcessingFn = runMetabImportProcessing
) {
  observeEventFn(input$process_import, {
    processingArgs <- list(
      assay1Data = localData$assay1_data,
      assay1Name = input$assay1_name,
      assay1File = localData$assay1_file,
      assay2File = localData$assay2_file,
      assay2Name = input$assay2_name,
      vendorFormat = input$vendor_format,
      detectedFormat = localData$detected_format,
      sanitizeNames = input$sanitize_names,
      isPattern = input$is_pattern,
      getMetaboliteIdColFn = columnAccessors$getMetaboliteIdCol,
      getAnnotationColFn = columnAccessors$getAnnotationCol,
      getSampleColumnsFn = columnAccessors$getSampleColumns,
      workflowData = workflowData
    )
    if (metabImportServerFunctionAcceptsArg(runProcessingFn, "experimentPaths")) {
      processingArgs$experimentPaths <- experimentPaths
    }
    do.call(runProcessingFn, processingArgs)
  })

  invisible(NULL)
}

runMetabImportModuleServerShell <- function(
    input,
    output,
    session,
    id,
    workflowData,
    experimentPaths,
    volumes = NULL,
    requireNamespaceFn = requireNamespace,
    createReactiveValuesFn = shiny::reactiveValues,
    setupAssaySelectionCallbackFn = setupMetabImportAssaySelectionCallback,
    setupShinyFilesFn = setupMetabImportShinyFiles,
    setupStandardFileInputsFn = setupMetabImportStandardFileInputs,
    setupColumnAccessorsFn = setupMetabImportColumnAccessors,
    setupFileLoadedOutputFn = setupMetabImportFileLoadedOutput,
    setupFormatDetectionStatusOutputFn = setupMetabImportFormatDetectionStatusOutput,
    setupMetaboliteIdStatusOutputFn = setupMetabImportMetaboliteIdStatusOutput,
    setupAnnotationStatusOutputFn = setupMetabImportAnnotationStatusOutput,
    setupSampleColumnsDisplayOutputFn = setupMetabImportSampleColumnsDisplayOutput,
    setupAvailableColumnsDisplayOutputFn = setupMetabImportAvailableColumnsDisplayOutput,
    setupCustomMetaboliteIdStatusOutputFn = setupMetabImportCustomMetaboliteIdStatusOutput,
    setupCustomAnnotationStatusOutputFn = setupMetabImportCustomAnnotationStatusOutput,
    setupValidationSummaryOutputFn = setupMetabImportValidationSummaryOutput,
    setupProcessingObserverFn = setupMetabImportProcessingObserver,
    setupStatusOutputFn = setupMetabImportStatusOutput,
    logMessageFn = message,
    isTestModeFn = is_test_mode
) {
  logMessageFn("   mod_metab_import_server: Inside moduleServer function")

  useShinyFiles <- requireNamespaceFn("shinyFiles", quietly = TRUE) && !isTRUE(isTestModeFn())
  logMessageFn(sprintf("   mod_metab_import_server: shinyFiles available = %s", useShinyFiles))

  localData <- createReactiveValuesFn(
    assay1_file = NULL,
    assay1_data = NULL,
    assay1_import_result = NULL,
    assay2_file = NULL,
    assay2_data = NULL,
    assay2_import_result = NULL,
    detected_format = NULL,
    format_confidence = NULL,
    all_headers = NULL
  )

  importData <- setupAssaySelectionCallbackFn(
    localData = localData,
    session = session
  )

  if (useShinyFiles) {
    volumes <- setupShinyFilesFn(
      input = input,
      output = output,
      session = session,
      volumes = volumes,
      localData = localData,
      importDataFn = importData
    )
  } else {
    setupStandardFileInputsFn(
      input = input,
      output = output,
      localData = localData,
      importDataFn = importData
    )
  }

  columnAccessors <- setupColumnAccessorsFn(
    input = input,
    localData = localData
  )

  setupFileLoadedOutputFn(
    output = output,
    localData = localData
  )

  setupFormatDetectionStatusOutputFn(
    output = output,
    localData = localData
  )

  setupMetaboliteIdStatusOutputFn(
    output = output,
    input = input,
    localData = localData
  )

  setupAnnotationStatusOutputFn(
    output = output,
    input = input,
    localData = localData
  )

  setupSampleColumnsDisplayOutputFn(
    output = output,
    localData = localData
  )

  setupAvailableColumnsDisplayOutputFn(
    output = output,
    localData = localData
  )

  setupCustomMetaboliteIdStatusOutputFn(
    output = output,
    input = input,
    localData = localData
  )

  setupCustomAnnotationStatusOutputFn(
    output = output,
    input = input,
    localData = localData
  )

  setupValidationSummaryOutputFn(
    output = output,
    localData = localData,
    columnAccessors = columnAccessors
  )

  processingObserverArgs <- list(
    input = input,
    localData = localData,
    columnAccessors = columnAccessors,
    workflowData = workflowData
  )
  if (metabImportServerFunctionAcceptsArg(setupProcessingObserverFn, "experimentPaths")) {
    processingObserverArgs$experimentPaths <- experimentPaths
  }
  do.call(setupProcessingObserverFn, processingObserverArgs)

  setupStatusOutputFn(
    output = output,
    workflowData = workflowData
  )

  invisible(localData)
}

runMetabImportModuleServerEntryShell <- function(
    id,
    workflowData,
    experimentPaths,
    volumes = NULL,
    moduleServerFn = shiny::moduleServer,
    runModuleServerShellFn = runMetabImportModuleServerShell,
    logMessageFn = message
) {
  logMessageFn("--- Entering mod_metab_import_server ---")
  logMessageFn(sprintf("   mod_metab_import_server: volumes param is NULL = %s", is.null(volumes)))

  moduleServerFn(id, function(input, output, session) {
    runModuleServerShellFn(
      input = input,
      output = output,
      session = session,
      id = id,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      volumes = volumes
    )
  })
}

runMetabImportModuleServerPublicWrapper <- function(
    id,
    workflow_data,
    experiment_paths,
    volumes = NULL,
    runModuleServerEntryShellFn = runMetabImportModuleServerEntryShell
) {
  runModuleServerEntryShellFn(
    id = id,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    volumes = volumes
  )
}
