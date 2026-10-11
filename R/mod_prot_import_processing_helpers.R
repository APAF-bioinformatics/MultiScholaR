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

#' Apply the existing run-name mapping before a durable GUI handoff.
#' @param table Imported table; all non-Run values and attributes are retained.
#' @param run_col Configured run column.
#' @param make_clean_names Existing sanitizer, injectable for parity tests.
#' @return Table and number of unique names; NULL count means no run column.
#' @noRd
sanitizeProtImportTable <- function(table, run_col, make_clean_names = janitor::make_clean_names) {
    if (is.null(run_col) || !run_col %in% names(table)) return(list(data = table, count = NULL))
    original <- unique(as.character(table[[run_col]]))
    mapping <- stats::setNames(make_clean_names(original), original)
    table[[run_col]] <- mapping[as.character(table[[run_col]])]
    list(data = table, count = length(original))
}

#' Describe import results without retaining the source table in the summary.
#' @param result Existing public import result.
#' @return Small counts used by the existing setup log and completion message.
#' @noRd
describeProtImportData <- function(result) {
    count <- function(column) {
        if (is.null(column) || !column %in% names(result$data)) return(NA)
        length(unique(result$data[[column]]))
    }
    list(n_rows = nrow(result$data), n_runs = count(result$column_mapping$run_col),
        n_proteins = count(result$column_mapping$protein_col),
        n_peptides = count(result$column_mapping$peptide_col))
}

#' Store bounded import/design metadata beside its immutable source reference.
#' @param manager Active DIA manager.
#' @param table Sanitized imported table, already validated for persistence.
#' @param result Original import result; used only to preserve original counts.
#' @param reference Committed imported-table reference.
#' @return Bounded metadata, with no abundance columns or parser pointers.
#' @noRd
recordDIAImportMetadata <- function(manager, table, result, reference) {
    mapping <- result$column_mapping
    columns <- intersect(c(mapping$run_col, "Batch"), names(table))
    samples <- dplyr::distinct(table[, columns, drop = FALSE])
    samples <- structure(lapply(samples, identity), names = names(samples), class = "data.frame",
        row.names = .set_row_names(nrow(samples)))
    quantities <- unique(Filter(Negate(is.null), mapping[c("norm_quantity_col", "raw_quantity_col", "quantity_col")]))
    quantities <- unlist(quantities, use.names = FALSE)
    quantities <- quantities[quantities %in% names(table)]
    metadata <- list(data = NULL, reference_id = reference$entity_id, data_type = result$data_type,
        column_mapping = mapping, import_summary = describeProtImportData(result),
        samples = samples, numeric_columns = quantities,
        needs_numeric = quantities[vapply(table[quantities], function(column)
            !is.double(column) || !is.null(attributes(column)), logical(1))],
        protein_ids = if (!is.null(mapping$protein_col) && mapping$protein_col %in% names(table))
            unique(table[[mapping$protein_col]]) else NULL)
    if (as.numeric(object.size(metadata)) > manager$.context$store$max_metadata_bytes) {
        cli::cli_abort("DIA import metadata exceeds its bounded allowance", class = "saved_entity_limit")
    }
    snapshot <- manager$.context$snapshot
    snapshot$import_metadata <- saveDIAStateValue(manager$.context, metadata,
        "import-metadata", codec = "rds")
    publishDIAStateHead(manager$.context, snapshot)
    manager$.context$import_metadata_cache <- list(reference = snapshot$import_metadata, value = metadata)
    metadata
}

#' Read the bounded import summary without reading its source assay.
#' @param manager Current DIA owner.
#' @return Metadata or NULL; only this immutable bounded record is retained.
#' @noRd
getDIAImportMetadata <- function(manager) {
    context <- manager$.context
    reference <- context$snapshot$import_metadata
    if (is.null(reference)) return(NULL)
    cached <- context$import_metadata_cache
    if (identical(cached$reference, reference)) return(cached$value)
    metadata <- resolveDIAStateValue(context, reference)
    if (as.numeric(object.size(metadata)) > context$store$max_metadata_bytes) {
        cli::cli_abort("DIA import metadata exceeds its bounded allowance", class = "saved_entity_limit")
    }
    context$import_metadata_cache <- list(reference = reference, value = metadata)
    metadata
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
  candidate <- inherits(manager, "DIAStateAdapter") && identical(format, "diann")
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

processProtImportFastaData <- function(
    workflowData
    , fastaPath
    , organismName
    , experimentPaths = NULL
    , processFasta = processFastaFile
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
  updated_status$design_matrix <- "pending"
  workflowData$tab_status <- updated_status

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
