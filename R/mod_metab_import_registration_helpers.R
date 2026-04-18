setupMetabImportFileLoadedOutput <- function(
    output,
    localData,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions
) {
  output$file_loaded <- reactiveFn({
    !is.null(localData$assay1_data)
  })
  outputOptionsFn(output, "file_loaded", suspendWhenHidden = FALSE)

  invisible(NULL)
}

setupMetabImportStatusOutput <- function(
    output,
    workflowData,
    renderUIFn = shiny::renderUI,
    buildStatusFn = buildMetabImportStatus
) {
  output$import_status <- renderUIFn({
    buildStatusFn(
      setupImportStatus = workflowData$tab_status$setup_import,
      setupImportLog = workflowData$processing_log$setup_import
    )
  })

  invisible(NULL)
}

setupMetabImportValidationSummaryOutput <- function(
    output,
    localData,
    columnAccessors,
    renderUIFn = shiny::renderUI,
    buildValidationSummaryFn = buildMetabImportValidationSummary
) {
  output$validation_summary <- renderUIFn({
    buildValidationSummaryFn(
      assayData = localData$assay1_data,
      getMetaboliteIdColFn = columnAccessors$getMetaboliteIdCol,
      getSampleColumnsFn = columnAccessors$getSampleColumns
    )
  })

  invisible(NULL)
}

setupMetabImportFormatDetectionStatusOutput <- function(
    output,
    localData,
    renderUIFn = shiny::renderUI,
    buildFormatDetectionStatusFn = buildMetabImportFormatDetectionStatus
) {
  output$format_detection_status <- renderUIFn({
    buildFormatDetectionStatusFn(
      detectedFormat = localData$detected_format,
      formatConfidence = localData$format_confidence
    )
  })

  invisible(NULL)
}

setupMetabImportMetaboliteIdStatusOutput <- function(
    output,
    input,
    localData,
    renderUIFn = shiny::renderUI,
    buildMetaboliteIdStatusFn = buildMetabImportMetaboliteIdStatus
) {
  output$metabolite_id_status <- renderUIFn({
    buildMetaboliteIdStatusFn(
      assayData = localData$assay1_data,
      metaboliteIdCol = input$metabolite_id_col
    )
  })

  invisible(NULL)
}

setupMetabImportAnnotationStatusOutput <- function(
    output,
    input,
    localData,
    renderUIFn = shiny::renderUI,
    buildAnnotationStatusFn = buildMetabImportAnnotationStatus
) {
  output$annotation_status <- renderUIFn({
    buildAnnotationStatusFn(
      assayData = localData$assay1_data,
      annotationCol = input$annotation_col
    )
  })

  invisible(NULL)
}

setupMetabImportSampleColumnsDisplayOutput <- function(
    output,
    localData,
    renderTextFn = shiny::renderText,
    buildSampleColumnsDisplayFn = buildMetabImportSampleColumnsDisplay
) {
  output$sample_columns_display <- renderTextFn({
    buildSampleColumnsDisplayFn(
      importResult = localData$assay1_import_result
    )
  })

  invisible(NULL)
}

setupMetabImportAvailableColumnsDisplayOutput <- function(
    output,
    localData,
    renderTextFn = shiny::renderText,
    buildAvailableColumnsDisplayFn = buildMetabImportAvailableColumnsDisplay
) {
  output$available_columns_display <- renderTextFn({
    buildAvailableColumnsDisplayFn(
      allHeaders = localData$all_headers
    )
  })

  invisible(NULL)
}

setupMetabImportCustomMetaboliteIdStatusOutput <- function(
    output,
    input,
    localData,
    renderUIFn = shiny::renderUI,
    buildCustomMetaboliteIdStatusFn = buildMetabImportCustomMetaboliteIdStatus
) {
  output$metabolite_id_status_custom <- renderUIFn({
    buildCustomMetaboliteIdStatusFn(
      assayData = localData$assay1_data,
      columnName = input$metabolite_id_col_custom
    )
  })

  invisible(NULL)
}

setupMetabImportCustomAnnotationStatusOutput <- function(
    output,
    input,
    localData,
    renderUIFn = shiny::renderUI,
    buildCustomAnnotationStatusFn = buildMetabImportCustomAnnotationStatus
) {
  output$annotation_status_custom <- renderUIFn({
    buildCustomAnnotationStatusFn(
      assayData = localData$assay1_data,
      columnName = input$annotation_col_custom
    )
  })

  invisible(NULL)
}

