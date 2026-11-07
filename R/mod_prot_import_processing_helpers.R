buildProtImportProcessingModal <- function(ns, initialStatus = "Initializing...") {
  shiny::modalDialog(
    title = shiny::tagList(
      shiny::icon("cog", class = "fa-spin")
      , " Processing Data"
    )
    , shiny::tags$div(
        style = "text-align: center; padding: 20px;"
        , shiny::tags$div(
            style = "font-size: 48px; color: #3c8dbc; margin-bottom: 20px;"
            , shiny::icon("spinner", class = "fa-spin fa-pulse")
          )
        , shiny::tags$div(
            id = ns("processing_status_text")
            , style = "font-size: 16px; color: #333;"
            , initialStatus
          )
        , shiny::tags$div(
            style = "margin-top: 15px; font-size: 12px; color: #999;"
            , "Please wait while your data is being processed."
          )
      )
    , footer = NULL
    , size = "s"
    , easyClose = FALSE
  )
}

updateProtImportProcessingStatus <- function(
    statusText
    , updateHtml = shinyjs::html
    , messageFn = message
) {
  messageFn(sprintf("[mod_prot_import] %s", statusText))
  updateHtml("processing_status_text", statusText)
}

importProtImportDataByFormat <- function(
    format
    , searchResultsPath
    , input
    , importDiann = importDIANNData
    , importSpectronaut = importSpectronautData
    , importFragpipe = importFragPipeData
    , importMaxquant = importMaxQuantData
    , importPdTmt = importProteomeDiscovererTMTData
    , logError = logger::log_error
) {
  # KIT-102 P2: the reading is the core's importProteomicsData; this GUI-side name maps the inputs to its options.
  importProteomicsData(
    format = format
    , path = searchResultsPath
    , options = list(
        use_precursor_norm = input$diann_use_precursor_norm
        , quantity_type = input$spectronaut_quantity
        , use_maxlfq = input$fragpipe_use_maxlfq
        , use_lfq = input$maxquant_use_lfq
        , filter_contaminants = input$maxquant_filter_contaminants
      )
    , readers = list(
        diann = importDiann
        , spectronaut = importSpectronaut
        , fragpipe = importFragpipe
        , maxquant = importMaxquant
        , pd_tmt = importPdTmt
      )
    , logError = logError
  )
}

sanitizeProtImportRunNames <- function(
    workflowData
    , runCol
    , makeCleanNames = janitor::make_clean_names
    , logInfo = logger::log_info
    , showNotification = shiny::showNotification
) {
  table <- getDIAOwnedPayload(workflowData, "data_tbl")
  if (is.null(runCol) || !(runCol %in% names(table))) {
    return(FALSE)
  }

  sanitized <- sanitizeProtImportTable(table, runCol, makeCleanNames)
  workflowData <- setDIAOwnedPayload(workflowData, "data_tbl", sanitized$data)

  logInfo(sprintf("Sanitized %d unique sample names.", sanitized$count))
  showNotification("Sample names sanitized for R compatibility.", type = "message")

  TRUE
}

applyProtImportResultToWorkflow <- function(
    workflowData
    , dataImportResult
    , format
    , fastaPath
    , sanitizeNames = FALSE
    , logInfo = logger::log_info
    , showNotification = shiny::showNotification
    , sanitizeRunNames = sanitizeProtImportRunNames
) {
  manager <- getDIAActiveManager()
  # KIT-020: every format the owner serves (DIA-NN and the protein-level readers) takes the owned path.
  candidate <- inherits(manager, "DIAStateAdapter") && !is.na(workflowTypeForImportFormat(format, dataImportResult$data_type))
  if (candidate) {
    table <- prepareDIAParserTable(dataImportResult$data)
    if (isTRUE(sanitizeNames)) {
      sanitized <- sanitizeProtImportTable(table, dataImportResult$column_mapping$run_col)
      table <- sanitized$data
      if (!is.null(sanitized$count)) {
        logInfo(sprintf("Sanitized %d unique sample names.", sanitized$count))
        showNotification("Sample names sanitized for R compatibility.", type = "message")
      }
      sanitized <- NULL
    }
    workflowData <- setDIAOwnedPayload(workflowData, "data_tbl", table)
    recordDIAImportMetadata(manager, table, dataImportResult, workflowData$data_tbl)
    table <- NULL
  } else workflowData <- setDIAOwnedPayload(workflowData, "data_tbl", dataImportResult$data)
  workflowData$data_format <- format
  workflowData$data_type <- dataImportResult$data_type
  workflowData$column_mapping <- dataImportResult$column_mapping

  if (isTRUE(sanitizeNames) && !candidate) {
    logInfo("Sanitizing sample names using janitor::make_clean_names()...")
    sanitizeRunNames(
      workflowData = workflowData
      , runCol = dataImportResult$column_mapping$run_col
      , logInfo = logInfo
      , showNotification = showNotification
    )
  }

  workflowData$data_cln <- workflowData$data_tbl

  # KIT-046: a DIA-NN import keeps the route type the entry point set on the manager (`DIA` or `DIA_LIMPA`;
  # on an owner it is the store's identity); every other format's type follows the format.
  manager_type <- tryCatch(workflowData$state_manager$workflow_type, error = function(e) NULL)
  workflowType <- if (identical(format, "diann") && identical(manager_type, "DIA_LIMPA")) "DIA_LIMPA" else workflowTypeForImportFormat(format)
  if (is.na(workflowType)) workflowType <- "LFQ"
  workflowData$state_manager$setWorkflowType(workflowType)
  workflowData$fasta_file_path <- fastaPath

  invisible(workflowType)
}

finalizeProtImportSetupState <- function(
    workflowData
    , dataImportResult
    , format
    , searchFilename
    , fastaFilename
    , taxonId
    , organismName
    , mixedSpeciesFasta = FALSE
    , now = Sys.time
) {
  summary <- dataImportResult$import_summary %||% describeProtImportData(dataImportResult)

  workflowData$taxon_id <- taxonId
  workflowData$organism_name <- organismName

  timestamp <- now()

  if (!isTRUE(mixedSpeciesFasta)) {
    workflowData$mixed_species_analysis <- list(
      enabled = FALSE,
      selected_taxon_id = taxonId,
      selected_organism = organismName,
      organism_distribution = NULL,
      organism_mapping = NULL,
      filter_applied_at_import = FALSE,
      timestamp = timestamp
    )
  }

  workflowData$processing_log$setup_import <- list(
    timestamp = timestamp,
    search_file = searchFilename,
    detected_format = format,
    data_type = dataImportResult$data_type,
    fasta_file = fastaFilename,
    taxon_id = taxonId,
    organism = organismName,
    n_rows = summary$n_rows,
    n_runs = summary$n_runs,
    n_proteins = summary$n_proteins,
    n_peptides = summary$n_peptides
  )

  invisible(workflowData$processing_log$setup_import)
}

completeProtImportSuccessState <- function(
    workflowData
    , localData
    , dataImportResult
    , format
    , complete = function(notify) completeOmicImport(workflowData, omicGuiConfig("proteomics"),
        omicImportSourcePanel("proteomics")$design_status, notify)
    , removeModal = shiny::removeModal
    , showNotification = shiny::showNotification
    , messageFn = message
) {
  # KIT-102 W6: the import frame's completion helper makes today's single write (setup_import complete, design_matrix
  # pending). Its notification is not used: the modal closes and the success text shows below, after the log lines.
  updated_status <- complete(notify = FALSE)

  messageFn("========================================")
  messageFn("[mod_prot_import] Data import completed successfully!")
  messageFn(sprintf(
    "[mod_prot_import] Rows: %d, Proteins: %s, Format: %s",
    dataImportResult$import_summary$n_rows %||% nrow(dataImportResult$data),
    workflowData$processing_log$setup_import$n_proteins,
    format
  ))
  messageFn("========================================")

  if (!isTRUE(localData$waiting_for_organism_selection)) {
    removeModal()
    showNotification("Data import successful!", type = "message")
  } else {
    messageFn("[mod_prot_import] Waiting for organism selection from user...")
  }

  localData$processing <- FALSE

  invisible(updated_status)
}

resetProtImportWorkflowStateOnError <- function(
    workflowData
    , localData
    , errorMessage
    , logError = logger::log_error
    , removeModal = shiny::removeModal
    , showNotification = shiny::showNotification
    , messageFn = message
) {
  messageFn("========================================")
  messageFn(sprintf("[mod_prot_import] ERROR: %s", errorMessage))
  messageFn("========================================")

  logError(paste("Error during data import:", errorMessage))
  removeModal()
  showNotification(paste("Error:", errorMessage), type = "error", duration = 10)
  localData$processing <- FALSE

  workflowData$data_tbl <- NULL
  workflowData$data_format <- NULL
  workflowData$data_type <- NULL
  workflowData$column_mapping <- NULL
  workflowData$data_cln <- NULL
  workflowData$fasta_file_path <- NULL
  workflowData$aa_seq_tbl_final <- NULL
  workflowData$config_list <- NULL
  workflowData$processing_log$setup_import <- NULL
  updated_status <- workflowData$tab_status
  updated_status$setup_import <- "incomplete"
  workflowData$tab_status <- updated_status

  invisible(updated_status)
}
