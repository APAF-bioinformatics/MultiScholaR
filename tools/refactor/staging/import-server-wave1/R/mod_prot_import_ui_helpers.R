buildProtImportStatusUi <- function(isProcessing, dataTbl, dataFormat, dataType) {
  if (isProcessing) {
    return(shiny::tags$div(
      class = "alert alert-info",
      shiny::tags$strong("Processing..."),
      " Please wait while data is being imported and validated."
    ))
  }

  if (!is.null(dataTbl)) {
    return(shiny::tags$div(
      class = "alert alert-success",
      shiny::tags$strong("[OK] Data imported successfully!"),
      shiny::tags$br(),
      paste("Format:", dataFormat, "| Type:", dataType)
    ))
  }

  NULL
}

buildProtImportDataSummaryUi <- function(dataTbl, setupImportLog) {
  shiny::req(dataTbl)
  shiny::req(setupImportLog)

  summaryInfo <- setupImportLog

  listItems <- list(
    shiny::tags$li(paste("Data format:", toupper(summaryInfo$detected_format))),
    shiny::tags$li(paste("Data type:", summaryInfo$data_type)),
    shiny::tags$li(paste("Total rows:", format(summaryInfo$n_rows, big.mark = ","))),
    shiny::tags$li(paste("Number of runs:", summaryInfo$n_runs)),
    shiny::tags$li(paste("Number of protein groups:", format(summaryInfo$n_proteins, big.mark = ",")))
  )

  if (!is.null(summaryInfo$n_peptides) && !is.na(summaryInfo$n_peptides)) {
    listItems <- append(listItems, list(
      shiny::tags$li(paste("Number of peptides:", format(summaryInfo$n_peptides, big.mark = ",")))
    ))
  }

  listItems <- append(listItems, list(
    shiny::tags$li(paste("Organism:", summaryInfo$organism, "(taxon:", summaryInfo$taxon_id, ")"))
  ))

  shiny::tags$div(
    class = "well",
    shiny::h4("Data Summary"),
    do.call(shiny::tags$ul, listItems)
  )
}

buildProtImportFormatDetectionUi <- function(detectedFormat, formatConfidence) {
  shiny::req(detectedFormat)

  confidenceColor <- if (formatConfidence >= 0.8) {
    "success"
  } else if (formatConfidence >= 0.5) {
    "warning"
  } else {
    "danger"
  }

  formatDisplay <- switch(detectedFormat,
    "diann" = "DIA-NN",
    "spectronaut" = "Spectronaut DIA",
    "fragpipe" = "FragPipe LFQ",
    "maxquant" = "MaxQuant LFQ",
    "pd_tmt" = "Proteome Discoverer TMT",
    "unknown" = "Unknown format"
  )

  shiny::tags$div(
    class = paste("alert", paste0("alert-", confidenceColor)),
    shiny::tags$strong("Detected format: "),
    formatDisplay,
    shiny::tags$br(),
    shiny::tags$small(
      paste0("Confidence: ", round(formatConfidence * 100), "%")
    )
  )
}

buildProtImportFormatSpecificOptionsUi <- function(format, ns) {
  shiny::req(format)

  switch(format,
    "diann" = shiny::tagList(
      shiny::h5("DIA-NN Specific Options"),
      shiny::checkboxInput(ns("diann_use_precursor_norm"),
        "Use Precursor.Normalised values",
        value = TRUE
      )
    ),
    "spectronaut" = shiny::tagList(
      shiny::h5("Spectronaut Specific Options"),
      shiny::radioButtons(ns("spectronaut_quantity"),
        "Quantity type:",
        choices = c("PG.Quantity" = "pg",
          "PEP.Quantity" = "pep"),
        selected = "pg"
      )
    ),
    "fragpipe" = shiny::tagList(
      shiny::h5("FragPipe LFQ Options"),
      shiny::checkboxInput(ns("fragpipe_use_maxlfq"),
        "Use MaxLFQ intensities",
        value = FALSE
      )
    ),
    "maxquant" = shiny::tagList(
      shiny::h5("MaxQuant LFQ Options"),
      shiny::checkboxInput(ns("maxquant_use_lfq"),
        "Use LFQ intensities",
        value = TRUE
      ),
      shiny::checkboxInput(ns("maxquant_filter_contaminants"),
        "Filter contaminants",
        value = TRUE
      )
    ),
    shiny::tags$div(
      class = "alert alert-warning",
      "Format-specific options not available"
    )
  )
}

