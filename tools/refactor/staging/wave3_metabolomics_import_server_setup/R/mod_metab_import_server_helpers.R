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
    observeEventFn = shiny::observeEvent,
    runProcessingFn = runMetabImportProcessing
) {
  observeEventFn(input$process_import, {
    runProcessingFn(
      assay1Data = localData$assay1_data,
      assay1Name = input$assay1_name,
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
  })

  invisible(NULL)
}

