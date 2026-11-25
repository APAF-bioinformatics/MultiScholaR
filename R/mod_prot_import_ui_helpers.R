# KIT-250: the organism dialog's texts (drafts; the owner sets the final wording, VAL-41).
PROT_IMPORT_ORGANISM_CANCELLED_TEXT <- paste(
  "Import not finished: no organism was chosen. Press 'Process Imported Data' again and choose the organism in the",
  "window that follows, or untick 'Mixed species FASTA' and enter the organism yourself."
)
PROT_IMPORT_ORGANISM_WAITING_TEXT <- "Choose the organism in the window to finish the import."
# The Cancel warning's fixed id: one copy at a time, removed when the next import starts (KIT-175's pattern).
PROT_IMPORT_ORGANISM_CANCELLED_ID <- "prot_import_organism_cancelled"

buildProtImportStatusUi <- function(isProcessing, dataTbl, dataFormat, dataType, waitingForOrganism = FALSE,
                                    organismUnconfirmed = FALSE) {
  if (isProcessing) {
    return(shiny::tags$div(
      class = "alert alert-info",
      shiny::tags$strong("Processing..."),
      " Please wait while data is being imported and validated."
    ))
  }

  # KIT-250: with 'Mixed species FASTA' the import is finished only when the organism is confirmed.
  if (isTRUE(waitingForOrganism)) {
    return(shiny::tags$div(class = "alert alert-info", PROT_IMPORT_ORGANISM_WAITING_TEXT))
  }
  if (isTRUE(organismUnconfirmed)) {
    return(shiny::tags$div(class = "alert alert-warning", PROT_IMPORT_ORGANISM_CANCELLED_TEXT))
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

buildProtImportDataSummaryUi <- function(dataTbl, setupImportLog, organismPending = FALSE) {
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

  # KIT-250: no organism is shown before one is confirmed in the organism dialog.
  organismLine <- if (isTRUE(organismPending)) {
    "Organism: not confirmed"
  } else {
    paste("Organism:", summaryInfo$organism, "(taxon:", summaryInfo$taxon_id, ")")
  }
  listItems <- append(listItems, list(shiny::tags$li(organismLine)))

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
      # KIT-118: the "Use Precursor.Normalised values" checkbox is gone: it changed nothing in the
      # analysis, which reads both columns whatever it said.
      shiny::helpText(paste("Both DIA-NN quantity columns are read. Protein quantification (IQ) and",
        "the technical-replicate imputation use Precursor.Quantity; the peptide filters and limpa use",
        "Precursor.Normalised."))
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
        value = TRUE
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
