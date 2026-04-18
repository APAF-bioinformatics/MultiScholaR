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
  tryCatch({
    switch(format,
      "diann" = importDiann(
        filepath = searchResultsPath,
        use_precursor_norm = input$diann_use_precursor_norm %||% TRUE
      ),
      "spectronaut" = importSpectronaut(
        filepath = searchResultsPath,
        quantity_type = input$spectronaut_quantity %||% "pg"
      ),
      "fragpipe" = importFragpipe(
        filepath = searchResultsPath,
        use_maxlfq = input$fragpipe_use_maxlfq %||% TRUE
      ),
      "maxquant" = importMaxquant(
        filepath = searchResultsPath,
        use_lfq = input$maxquant_use_lfq %||% TRUE,
        filter_contaminants = input$maxquant_filter_contaminants %||% TRUE
      ),
      "pd_tmt" = importPdTmt(
        filepath = searchResultsPath
      ),
      stop("Unsupported format: ", format)
    )
  }, error = function(e) {
    logError(paste("Error in data import function:", e$message))
    stop("Failed to import data: ", e$message)
  })
}

sanitizeProtImportRunNames <- function(
    workflowData
    , runCol
    , makeCleanNames = janitor::make_clean_names
    , logInfo = logger::log_info
    , showNotification = shiny::showNotification
) {
  if (is.null(runCol) || !(runCol %in% names(workflowData$data_tbl))) {
    return(FALSE)
  }

  original_runs <- unique(as.character(workflowData$data_tbl[[runCol]]))
  clean_runs <- makeCleanNames(original_runs)
  name_mapping <- stats::setNames(clean_runs, original_runs)

  workflowData$data_tbl[[runCol]] <- name_mapping[as.character(workflowData$data_tbl[[runCol]])]

  logInfo(sprintf("Sanitized %d unique sample names.", length(original_runs)))
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
  workflowData$data_tbl <- dataImportResult$data
  workflowData$data_format <- format
  workflowData$data_type <- dataImportResult$data_type
  workflowData$column_mapping <- dataImportResult$column_mapping

  if (isTRUE(sanitizeNames)) {
    logInfo("Sanitizing sample names using janitor::make_clean_names()...")
    sanitizeRunNames(
      workflowData = workflowData
      , runCol = dataImportResult$column_mapping$run_col
      , logInfo = logInfo
      , showNotification = showNotification
    )
  }

  workflowData$data_cln <- workflowData$data_tbl

  workflowType <- switch(format,
    "pd_tmt" = "TMT",
    "diann" = "DIA",
    "LFQ"
  )
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
  run_col <- dataImportResult$column_mapping$run_col
  protein_col <- dataImportResult$column_mapping$protein_col
  peptide_col <- dataImportResult$column_mapping$peptide_col

  n_runs <- if (!is.null(run_col) && run_col %in% names(dataImportResult$data)) {
    length(unique(dataImportResult$data[[run_col]]))
  } else {
    NA
  }

  n_proteins <- if (!is.null(protein_col) && protein_col %in% names(dataImportResult$data)) {
    length(unique(dataImportResult$data[[protein_col]]))
  } else {
    NA
  }

  n_peptides <- if (!is.null(peptide_col) && peptide_col %in% names(dataImportResult$data)) {
    length(unique(dataImportResult$data[[peptide_col]]))
  } else {
    NA
  }

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
    n_rows = nrow(dataImportResult$data),
    n_runs = n_runs,
    n_proteins = n_proteins,
    n_peptides = n_peptides
  )

  invisible(workflowData$processing_log$setup_import)
}

processProtImportFastaData <- function(
    workflowData
    , fastaPath
    , organismName
    , experimentPaths = NULL
    , processFasta = processFastaFile
    , assignFn = assign
    , assignEnv = .GlobalEnv
    , saveRds = saveRDS
    , logInfo = logger::log_info
    , logWarn = logger::log_warn
    , dirExists = dir.exists
    , dirCreate = dir.create
    , tempdirFn = tempdir
) {
  logInfo("Processing FASTA file...")

  uniprot_mapping <- if (!is.null(workflowData$uniprot_mapping)) {
    workflowData$uniprot_mapping
  } else {
    NULL
  }

  uniparc_mapping <- if (!is.null(workflowData$uniparc_mapping)) {
    workflowData$uniparc_mapping
  } else {
    NULL
  }

  logInfo("Checking experiment_paths...")
  if (!is.null(experimentPaths) && !is.null(experimentPaths$results_dir)) {
    fasta_meta_file <- file.path(
      experimentPaths$results_dir,
      "cache",
      "aa_seq_tbl.RDS"
    )
    cache_dir <- file.path(experimentPaths$results_dir, "cache")
  } else {
    logWarn("experiment_paths not properly initialized, using temp directory for cache")
    cache_dir <- file.path(tempdirFn(), "proteomics_cache")
    fasta_meta_file <- file.path(cache_dir, "aa_seq_tbl.RDS")
  }

  if (!dirExists(cache_dir)) {
    dirCreate(cache_dir, recursive = TRUE)
  }

  tryCatch({
    fasta_result <- processFasta(
      fasta_file_path = fastaPath,
      uniprot_search_results = uniprot_mapping,
      uniparc_search_results = uniparc_mapping,
      fasta_meta_file = fasta_meta_file,
      organism_name = organismName
    )

    aa_seq_tbl_final <- fasta_result$aa_seq_tbl_final
    fasta_metadata <- fasta_result$fasta_metadata

    workflowData$aa_seq_tbl_final <- aa_seq_tbl_final
    workflowData$fasta_metadata <- fasta_metadata

    assignFn("aa_seq_tbl_final", aa_seq_tbl_final, envir = assignEnv)

    if (!is.null(experimentPaths) && !is.null(experimentPaths$source_dir)) {
      scripts_aa_seq_path <- file.path(experimentPaths$source_dir, "aa_seq_tbl_final.RDS")
      saveRds(aa_seq_tbl_final, scripts_aa_seq_path)
      logInfo(sprintf("Saved aa_seq_tbl_final to scripts directory: %s", scripts_aa_seq_path))

      fasta_metadata_path <- file.path(experimentPaths$source_dir, "fasta_metadata.RDS")
      saveRds(fasta_metadata, fasta_metadata_path)
      logInfo(sprintf("Saved FASTA metadata to scripts directory: %s", fasta_metadata_path))
      logInfo(sprintf("FASTA Format: %s, Sequences: %d", fasta_metadata$fasta_format, fasta_metadata$num_sequences))
    }

    logInfo(sprintf("FASTA file processed successfully. Found %d sequences", nrow(aa_seq_tbl_final)))

    invisible(list(
      success = TRUE,
      cacheDir = cache_dir,
      fastaMetaFile = fasta_meta_file
    ))
  }, error = function(e) {
    logWarn(paste("Error processing FASTA file:", e$message))
    logWarn("Continuing without FASTA processing - protein ID conversion will be skipped")
    workflowData$aa_seq_tbl_final <- NULL
    workflowData$fasta_metadata <- NULL

    invisible(list(
      success = FALSE,
      cacheDir = cache_dir,
      fastaMetaFile = fasta_meta_file,
      error = e$message
    ))
  })
}

completeProtImportSuccessState <- function(
    workflowData
    , localData
    , dataImportResult
    , format
    , removeModal = shiny::removeModal
    , showNotification = shiny::showNotification
    , messageFn = message
) {
  updated_status <- workflowData$tab_status
  updated_status$setup_import <- "complete"
  workflowData$tab_status <- updated_status

  messageFn("========================================")
  messageFn("[mod_prot_import] Data import completed successfully!")
  messageFn(sprintf(
    "[mod_prot_import] Rows: %d, Proteins: %s, Format: %s",
    nrow(dataImportResult$data),
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

