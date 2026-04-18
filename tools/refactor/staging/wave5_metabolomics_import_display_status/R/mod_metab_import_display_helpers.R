handleMetabImportAssayFileSelection <- function(
    selectedInput,
    volumes,
    localData,
    localField,
    output,
    outputId,
    onSelected = NULL,
    parseFilePathsFn = shinyFiles::parseFilePaths,
    renderTextFn = shiny::renderText
) {
  if (is.null(selectedInput) || is.integer(selectedInput)) {
    return(invisible(NULL))
  }

  fileInfo <- parseFilePathsFn(volumes, selectedInput)
  if (nrow(fileInfo) == 0) {
    return(invisible(NULL))
  }

  selectedPath <- as.character(fileInfo$datapath[1])
  localData[[localField]] <- selectedPath
  output[[outputId]] <- renderTextFn(selectedPath)

  if (is.function(onSelected)) {
    onSelected(selectedPath)
  }

  invisible(selectedPath)
}

buildMetabImportFormatDetectionStatus <- function(
    detectedFormat,
    formatConfidence,
    reqFn = shiny::req,
    divFn = shiny::tags$div,
    strongFn = shiny::tags$strong,
    brFn = shiny::tags$br,
    smallFn = shiny::tags$small
) {
  reqFn(detectedFormat)

  confidencePct <- round(formatConfidence * 100)
  colorClass <- if (confidencePct >= 70) "success" else if (confidencePct >= 40) "warning" else "danger"

  formatDisplay <- switch(detectedFormat,
    "msdial" = "MS-DIAL",
    "progenesis" = "Progenesis QI",
    "xcms" = "XCMS",
    "compound_discoverer" = "Compound Discoverer",
    "Unknown"
  )

  divFn(
    class = paste("alert", paste0("alert-", colorClass)),
    strongFn("Detected format: "),
    formatDisplay,
    brFn(),
    smallFn(sprintf("Confidence: %d%%", confidencePct))
  )
}

buildMetabImportMetaboliteIdStatus <- function(
    assayData,
    metaboliteIdCol,
    reqFn = shiny::req,
    namesFn = names,
    uniqueFn = unique,
    lengthFn = length,
    spanFn = shiny::tags$span,
    iconFn = shiny::icon,
    sprintfFn = sprintf
) {
  reqFn(assayData, metaboliteIdCol)

  if (metaboliteIdCol %in% namesFn(assayData)) {
    uniqueCount <- lengthFn(uniqueFn(assayData[[metaboliteIdCol]]))

    return(spanFn(
      iconFn("check-circle", style = "color: green;"),
      sprintfFn(" %d unique IDs", uniqueCount)
    ))
  }

  spanFn(
    iconFn("times-circle", style = "color: red;"),
    " Column not found"
  )
}

buildMetabImportCustomMetaboliteIdStatus <- function(
    assayData,
    columnName,
    reqFn = shiny::req,
    nzcharFn = nzchar,
    namesFn = names,
    resolveColumnNameFn = resolveMetabImportColumnName,
    uniqueFn = unique,
    lengthFn = length,
    spanFn = shiny::tags$span,
    iconFn = shiny::icon,
    sprintfFn = sprintf
) {
  reqFn(assayData)

  if (is.null(columnName) || !nzcharFn(columnName)) {
    return(spanFn(
      iconFn("question-circle", style = "color: gray;"),
      " Enter column name"
    ))
  }

  headers <- namesFn(assayData)
  actualCol <- resolveColumnNameFn(headers = headers, columnName = columnName)

  if (actualCol %in% headers) {
    uniqueCount <- lengthFn(uniqueFn(assayData[[actualCol]]))

    return(spanFn(
      iconFn("check-circle", style = "color: green;"),
      sprintfFn(" Found: %d unique IDs", uniqueCount)
    ))
  }

  spanFn(
    iconFn("times-circle", style = "color: red;"),
    " Column not found"
  )
}

buildMetabImportCustomAnnotationStatus <- function(
    assayData,
    columnName,
    reqFn = shiny::req,
    nzcharFn = nzchar,
    namesFn = names,
    resolveColumnNameFn = resolveMetabImportColumnName,
    spanFn = shiny::tags$span,
    iconFn = shiny::icon
) {
  reqFn(assayData)

  if (is.null(columnName) || !nzcharFn(columnName)) {
    return(spanFn(
      iconFn("minus-circle", style = "color: gray;"),
      " Optional"
    ))
  }

  headers <- namesFn(assayData)
  actualCol <- resolveColumnNameFn(headers = headers, columnName = columnName)

  if (actualCol %in% headers) {
    return(spanFn(
      iconFn("check-circle", style = "color: green;"),
      " Found"
    ))
  }

  spanFn(
    iconFn("times-circle", style = "color: red;"),
    " Column not found"
  )
}

buildMetabImportAnnotationStatus <- function(
    assayData,
    annotationCol,
    reqFn = shiny::req,
    nzcharFn = nzchar,
    namesFn = names,
    spanFn = shiny::tags$span,
    iconFn = shiny::icon
) {
  reqFn(assayData)

  if (!is.null(annotationCol) && nzcharFn(annotationCol)) {
    if (annotationCol %in% namesFn(assayData)) {
      return(spanFn(
        iconFn("check-circle", style = "color: green;"),
        " OK"
      ))
    }

    return(spanFn(
      iconFn("times-circle", style = "color: red;"),
      " Column not found"
    ))
  }

  spanFn(
    iconFn("minus-circle", style = "color: gray;"),
    " Optional"
  )
}

buildMetabImportSampleColumnsDisplay <- function(
    importResult,
    reqFn = shiny::req,
    lengthFn = length,
    headFn = utils::head,
    pasteFn = paste,
    sprintfFn = sprintf
) {
  reqFn(importResult)

  sampleCols <- importResult$sample_columns
  if (lengthFn(sampleCols) > 10) {
    return(pasteFn(
      pasteFn(headFn(sampleCols, 10), collapse = ", "),
      sprintfFn("... and %d more", lengthFn(sampleCols) - 10)
    ))
  }

  pasteFn(sampleCols, collapse = ", ")
}

buildMetabImportAvailableColumnsDisplay <- function(
    allHeaders,
    reqFn = shiny::req,
    pasteFn = paste
) {
  reqFn(allHeaders)
  pasteFn(allHeaders, collapse = ", ")
}

buildMetabImportValidationSummary <- function(
    assayData,
    getMetaboliteIdColFn,
    getSampleColumnsFn,
    validateColumnMappingFn = validateMetabColumnMapping,
    reqFn = shiny::req,
    lengthFn = length,
    tagListFn = shiny::tagList,
    divFn = shiny::tags$div,
    iconFn = shiny::icon,
    strongFn = shiny::tags$strong,
    ulFn = shiny::tags$ul,
    liFn = shiny::tags$li,
    lapplyFn = lapply,
    sprintfFn = sprintf
) {
  reqFn(assayData)

  metaboliteCol <- getMetaboliteIdColFn()
  sampleCols <- getSampleColumnsFn()

  reqFn(metaboliteCol)

  validation <- validateColumnMappingFn(
    data = assayData,
    metabolite_id_column = metaboliteCol,
    sample_columns = sampleCols
  )

  if (validation$valid) {
    return(tagListFn(
      divFn(
        class = "alert alert-success",
        iconFn("check-circle"),
        strongFn(" Validation Passed")
      ),
      ulFn(
        liFn(sprintfFn("Metabolites: %d", validation$summary$n_metabolites)),
        liFn(sprintfFn("Samples: %d", validation$summary$n_samples)),
        liFn(sprintfFn("Missing values: %.1f%%", validation$summary$pct_missing))
      ),
      if (lengthFn(validation$warnings) > 0) {
        divFn(
          class = "alert alert-warning",
          iconFn("exclamation-triangle"),
          " Warnings:",
          ulFn(
            lapplyFn(validation$warnings, liFn)
          )
        )
      }
    ))
  }

  divFn(
    class = "alert alert-danger",
    iconFn("times-circle"),
    strongFn(" Validation Failed"),
    ulFn(
      lapplyFn(validation$errors, liFn)
    )
  )
}

buildMetabImportStatus <- function(
    setupImportStatus,
    setupImportLog,
    divFn = shiny::tags$div,
    iconFn = shiny::icon,
    strongFn = shiny::tags$strong,
    brFn = shiny::tags$br,
    sprintfFn = sprintf,
    toupperFn = toupper
) {
  if (is.null(setupImportStatus) || !identical(setupImportStatus, "complete")) {
    return(NULL)
  }

  divFn(
    class = "alert alert-success",
    iconFn("check-circle"),
    strongFn(" Import Complete"),
    brFn(),
    sprintfFn(
      "Format: %s | Assays: %d | Samples: %d",
      toupperFn(setupImportLog$detected_format),
      setupImportLog$n_assays,
      setupImportLog$n_samples
    )
  )
}

