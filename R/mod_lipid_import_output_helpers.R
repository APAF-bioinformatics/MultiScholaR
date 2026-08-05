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

handleLipidImportAssayPathRender <- function(selectedPath) {
  selectedPath
}

handleLipidImportSelectedAssayPathAssignment <- function(
  selectedPath,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  renderText = shiny::renderText,
  renderAssayPath = handleLipidImportAssayPathRender
) {
  assignSelectedPath(selectedPath)
  setRenderedPath(renderText({
    renderAssayPath(selectedPath)
  }))

  if (is.function(onSelected)) {
    onSelected()
  }

  invisible(selectedPath)
}

handleLipidImportAssayFileSelectionEvent <- function(
  fileInput,
  volumes,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  handleFileSelection = handleLipidImportFileSelection,
  assignAssayPath = handleLipidImportSelectedAssayPathAssignment
) {
  handleFileSelection(
    fileInput = fileInput,
    volumes = volumes,
    onPathSelected = function(selectedPath) {
      assignAssayPath(
        selectedPath = selectedPath,
        assignSelectedPath = assignSelectedPath,
        setRenderedPath = setRenderedPath,
        onSelected = onSelected
      )
    }
  )
}

setupLipidImportAssayFileChooser <- function(
  input,
  fileInputId,
  volumes,
  session,
  filetypes = c("tsv", "tab", "txt", "csv", "xlsx", "parquet"),
  shinyFileChoose = shinyFiles::shinyFileChoose
) {
  shinyFileChoose(
    input = input,
    id = fileInputId,
    roots = volumes,
    session = session,
    filetypes = filetypes
  )
}

prepareLipidImportShinyFileVolumes <- function(
  volumes,
  getVolumes = shinyFiles::getVolumes,
  emitMessage = message
) {
  if (is.null(volumes)) {
    emitMessage("   mod_lipid_import_server: volumes is NULL, creating from getVolumes()")
    volumes <- getVolumes()()
  }

  emitMessage(sprintf(
    "   mod_lipid_import_server: volumes type = %s, length = %d",
    typeof(volumes), length(volumes)
  ))

  if (length(volumes) > 0) {
    emitMessage(sprintf(
      "   mod_lipid_import_server: volumes names = %s",
      paste(names(volumes), collapse = ", ")
    ))
  } else {
    emitMessage("   mod_lipid_import_server: WARNING - volumes is empty!")
  }

  volumes
}

registerLipidImportFileLoadedOutput <- function(
  output,
  localData,
  reactiveFn = shiny::reactive,
  outputOptionsFn = shiny::outputOptions
) {
  output$file_loaded <- reactiveFn({
    !is.null(localData$assay1_data)
  })
  outputOptionsFn(output, "file_loaded", suspendWhenHidden = FALSE)

  invisible(output)
}

registerLipidImportValidationSummaryOutput <- function(
  output,
  localData,
  getLipidIdCol,
  getSampleColumns,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportValidationSummaryRender
) {
  output$validation_summary <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdCol = getLipidIdCol(),
      sampleColumns = getSampleColumns()
    )
  })

  invisible(output)
}

registerLipidImportFormatDetectionStatusOutput <- function(
  output,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportFormatDetectionStatusRender
) {
  output$format_detection_status <- renderUi({
    handleRender(
      detectedFormat = localData$detected_format,
      formatConfidence = localData$format_confidence
    )
  })

  invisible(output)
}

registerLipidImportLipidIdStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportLipidIdStatusRender
) {
  output$lipid_id_status <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdCol = input$lipid_id_col
    )
  })

  invisible(output)
}

registerLipidImportAnnotationStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportAnnotationStatusRender
) {
  output$annotation_status <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      annotationCol = input$annotation_col
    )
  })

  invisible(output)
}

registerLipidImportSampleColumnsDisplayOutput <- function(
  output,
  localData,
  renderText = shiny::renderText,
  handleRender = handleLipidImportSampleColumnsDisplayRender
) {
  output$sample_columns_display <- renderText({
    handleRender(
      assay1ImportResult = localData$assay1_import_result
    )
  })

  invisible(output)
}

registerLipidImportAvailableColumnsDisplayOutput <- function(
  output,
  localData,
  renderText = shiny::renderText,
  handleRender = handleLipidImportAvailableColumnsDisplayRender
) {
  output$available_columns_display <- renderText({
    handleRender(
      allHeaders = localData$all_headers
    )
  })

  invisible(output)
}

registerLipidImportCustomLipidIdStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportCustomLipidIdStatusRender
) {
  output$lipid_id_status_custom <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdColCustom = input$lipid_id_col_custom
    )
  })

  invisible(output)
}

registerLipidImportCustomAnnotationStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportCustomAnnotationStatusRender
) {
  output$annotation_status_custom <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      annotationColCustom = input$annotation_col_custom
    )
  })

  invisible(output)
}

registerLipidImportStatusOutput <- function(
  output,
  workflowData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportStatusRender
) {
  output$import_status <- renderUi({
    handleRender(workflowData = workflowData)
  })

  invisible(output)
}

registerLipidImportModuleOutputs <- function(
  output,
  input,
  localData,
  workflowData,
  getLipidIdCol,
  getSampleColumns,
  registerFileLoadedOutput = registerLipidImportFileLoadedOutput,
  registerFormatDetectionStatusOutput = registerLipidImportFormatDetectionStatusOutput,
  registerLipidIdStatusOutput = registerLipidImportLipidIdStatusOutput,
  registerAnnotationStatusOutput = registerLipidImportAnnotationStatusOutput,
  registerSampleColumnsDisplayOutput = registerLipidImportSampleColumnsDisplayOutput,
  registerAvailableColumnsDisplayOutput = registerLipidImportAvailableColumnsDisplayOutput,
  registerCustomLipidIdStatusOutput = registerLipidImportCustomLipidIdStatusOutput,
  registerCustomAnnotationStatusOutput = registerLipidImportCustomAnnotationStatusOutput,
  registerValidationSummaryOutput = registerLipidImportValidationSummaryOutput,
  registerStatusOutput = registerLipidImportStatusOutput
) {
  registerFileLoadedOutput(
    output = output,
    localData = localData
  )
  registerFormatDetectionStatusOutput(
    output = output,
    localData = localData
  )
  registerLipidIdStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerAnnotationStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerSampleColumnsDisplayOutput(
    output = output,
    localData = localData
  )
  registerAvailableColumnsDisplayOutput(
    output = output,
    localData = localData
  )
  registerCustomLipidIdStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerCustomAnnotationStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerValidationSummaryOutput(
    output = output,
    localData = localData,
    getLipidIdCol = getLipidIdCol,
    getSampleColumns = getSampleColumns
  )
  registerStatusOutput(
    output = output,
    workflowData = workflowData
  )

  invisible(output)
}

handleLipidImportFormatDetectionStatusRender <- function(
  detectedFormat,
  formatConfidence,
  buildStatus = buildLipidImportFormatDetectionStatus
) {
  shiny::req(detectedFormat)

  buildStatus(
    detectedFormat = detectedFormat,
    formatConfidence = formatConfidence
  )
}

handleLipidImportAvailableColumnsDisplayRender <- function(
  allHeaders,
  formatPreviewText = formatLipidImportColumnPreviewText
) {
  shiny::req(allHeaders)

  formatPreviewText(columnNames = allHeaders)
}

handleLipidImportLipidIdStatusRender <- function(
  assay1Data,
  lipidIdCol,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = lipidIdCol,
    successMode = "unique_ids"
  )
}

handleLipidImportAnnotationStatusRender <- function(
  assay1Data,
  annotationCol,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = annotationCol,
    successMode = "ok",
    emptyMode = "optional"
  )
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

resolveLipidImportLipidIdColumn <- function(
  assay1Data,
  vendorFormat,
  lipidIdCol,
  lipidIdColCustom,
  resolveEffectiveColumn = resolveLipidImportEffectiveColumn
) {
  resolveEffectiveColumn(
    assayData = assay1Data,
    vendorFormat = vendorFormat,
    selectedColumn = lipidIdCol,
    customColumn = lipidIdColCustom
  )
}

resolveLipidImportAnnotationColumn <- function(
  assay1Data,
  vendorFormat,
  annotationCol,
  annotationColCustom,
  resolveEffectiveColumn = resolveLipidImportEffectiveColumn
) {
  resolveEffectiveColumn(
    assayData = assay1Data,
    vendorFormat = vendorFormat,
    selectedColumn = annotationCol,
    customColumn = annotationColCustom
  )
}

resolveLipidImportSelectedSampleColumns <- function(
  assay1Data,
  assay1ImportResult,
  vendorFormat,
  sampleColsPattern,
  excludeNormalized,
  resolveSampleColumns = resolveLipidImportSampleColumns
) {
  shiny::req(assay1Data)

  resolveSampleColumns(
    assayData = assay1Data,
    assayImportResult = assay1ImportResult,
    vendorFormat = vendorFormat,
    sampleColsPattern = sampleColsPattern,
    excludeNormalized = excludeNormalized
  )
}

buildLipidImportColumnSelectionReactives <- function(
  input,
  localData,
  reactiveFn = shiny::reactive,
  resolveLipidIdColumn = resolveLipidImportLipidIdColumn,
  resolveAnnotationColumn = resolveLipidImportAnnotationColumn,
  resolveSampleColumns = resolveLipidImportSelectedSampleColumns
) {
  list(
    lipidIdCol = reactiveFn({
      resolveLipidIdColumn(
        assay1Data = localData$assay1_data,
        vendorFormat = input$vendor_format,
        lipidIdCol = input$lipid_id_col,
        lipidIdColCustom = input$lipid_id_col_custom
      )
    }),
    annotationCol = reactiveFn({
      resolveAnnotationColumn(
        assay1Data = localData$assay1_data,
        vendorFormat = input$vendor_format,
        annotationCol = input$annotation_col,
        annotationColCustom = input$annotation_col_custom
      )
    }),
    sampleColumns = reactiveFn({
      resolveSampleColumns(
        assay1Data = localData$assay1_data,
        assay1ImportResult = localData$assay1_import_result,
        vendorFormat = input$vendor_format,
        sampleColsPattern = input$sample_cols_pattern,
        excludeNormalized = input$exclude_norm
      )
    })
  )
}

handleLipidImportStatusRender <- function(
  workflowData,
  buildStatusDisplay = buildLipidImportStatusDisplay
) {
  buildStatusDisplay(workflowData)
}

