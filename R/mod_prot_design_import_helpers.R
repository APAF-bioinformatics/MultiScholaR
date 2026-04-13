readProtDesignImportedContrasts <- function(contrastFile) {
  if (!file.exists(contrastFile)) {
    return(NULL)
  }

  contrastStrings <- readLines(contrastFile)
  if (!length(contrastStrings)) {
    return(NULL)
  }

  cleanStrings <- gsub("^group", "", contrastStrings)
  cleanStrings <- gsub("-group", "-", cleanStrings)
  friendlyNames <- gsub("-", "_vs_", cleanStrings)

  data.frame(
    contrasts = contrastStrings,
    friendly_names = friendlyNames,
    full_format = paste0(friendlyNames, "=", contrastStrings),
    stringsAsFactors = FALSE
  )
}

hydrateProtDesignImportedUniprotSidecar <- function(workflowData, importPath, sourceDir) {
  uniprotFileImport <- file.path(importPath, "uniprot_dat_cln.RDS")
  uniprotFileScripts <- file.path(sourceDir, "uniprot_dat_cln.RDS")

  if (file.exists(uniprotFileImport)) {
    logger::log_info("Loading uniprot_dat_cln from import directory.")
    uniprotDatCln <- readRDS(uniprotFileImport)
    workflowData$uniprot_dat_cln <- uniprotDatCln
    assign("uniprot_dat_cln", uniprotDatCln, envir = .GlobalEnv)

    saveRDS(uniprotDatCln, uniprotFileScripts)
    logger::log_info(sprintf(
      "Copied uniprot_dat_cln to scripts directory (%d annotations).",
      nrow(uniprotDatCln)
    ))

    shiny::showNotification(
      sprintf("UniProt annotations loaded: %d protein annotations available", nrow(uniprotDatCln)),
      type = "message",
      duration = 5
    )
    return(invisible(uniprotDatCln))
  }

  if (file.exists(uniprotFileScripts)) {
    logger::log_info("Loading existing uniprot_dat_cln from scripts directory.")
    uniprotDatCln <- readRDS(uniprotFileScripts)
    workflowData$uniprot_dat_cln <- uniprotDatCln
    assign("uniprot_dat_cln", uniprotDatCln, envir = .GlobalEnv)

    shiny::showNotification(
      sprintf("UniProt annotations loaded: %d protein annotations available", nrow(uniprotDatCln)),
      type = "message",
      duration = 5
    )
    return(invisible(uniprotDatCln))
  }

  logger::log_warn("No uniprot_dat_cln found. DE analysis will use accession IDs as gene names.")
  workflowData$uniprot_dat_cln <- NULL

  shiny::showNotification(
    "Warning: No UniProt annotations found. Gene names in DE results will be limited.",
    type = "warning",
    duration = 8
  )

  invisible(NULL)
}

hydrateProtDesignImportedFastaSidecar <- function(
    workflowData,
    importPath,
    sourceDir,
    resultsDir = NULL,
    fastaPath = NULL,
    organismName = NULL
) {
  aaSeqFileImport <- file.path(importPath, "aa_seq_tbl_final.RDS")
  aaSeqFileScripts <- file.path(sourceDir, "aa_seq_tbl_final.RDS")
  fastaMetadataFileImport <- file.path(importPath, "fasta_metadata.RDS")
  fastaMetadataFileScripts <- file.path(sourceDir, "fasta_metadata.RDS")

  if (file.exists(aaSeqFileImport)) {
    logger::log_info("Loading aa_seq_tbl_final from import directory.")
    aaSeqTblFinal <- readRDS(aaSeqFileImport)
    workflowData$aa_seq_tbl_final <- aaSeqTblFinal
    assign("aa_seq_tbl_final", aaSeqTblFinal, envir = .GlobalEnv)

    saveRDS(aaSeqTblFinal, aaSeqFileScripts)
    logger::log_info("Copied aa_seq_tbl_final to scripts directory for persistence.")

    if (file.exists(fastaMetadataFileImport)) {
      workflowData$fasta_metadata <- readRDS(fastaMetadataFileImport)
      saveRDS(workflowData$fasta_metadata, fastaMetadataFileScripts)
      logger::log_info("Loaded and copied FASTA metadata from import directory.")
    }
  } else if (file.exists(aaSeqFileScripts)) {
    logger::log_info("Loading existing aa_seq_tbl_final from scripts directory.")
    aaSeqTblFinal <- readRDS(aaSeqFileScripts)
    workflowData$aa_seq_tbl_final <- aaSeqTblFinal
    assign("aa_seq_tbl_final", aaSeqTblFinal, envir = .GlobalEnv)

    if (file.exists(fastaMetadataFileScripts)) {
      workflowData$fasta_metadata <- readRDS(fastaMetadataFileScripts)
      logger::log_info("Loaded FASTA metadata from scripts directory.")
    }
  } else {
    logger::log_warn("No aa_seq_tbl_final found. Protein accession cleanup will be skipped.")
    workflowData$aa_seq_tbl_final <- NULL
    workflowData$fasta_metadata <- NULL
  }

  if (is.null(workflowData$aa_seq_tbl_final) && !is.null(fastaPath)) {
    logger::log_info("Processing FASTA file for accession cleanup...")

    tryCatch({
      cacheDir <- if (!is.null(resultsDir)) {
        file.path(resultsDir, "cache")
      } else {
        file.path(tempdir(), "proteomics_cache")
      }

      if (!dir.exists(cacheDir)) {
        dir.create(cacheDir, recursive = TRUE)
      }

      fastaMetaFile <- file.path(cacheDir, "aa_seq_tbl.RDS")
      fastaResult <- processFastaFile(
        fasta_file_path = fastaPath,
        uniprot_search_results = NULL,
        uniparc_search_results = NULL,
        fasta_meta_file = fastaMetaFile,
        organism_name = organismName
      )

      aaSeqTblFinal <- fastaResult$aa_seq_tbl_final
      fastaMetadata <- fastaResult$fasta_metadata

      workflowData$aa_seq_tbl_final <- aaSeqTblFinal
      workflowData$fasta_metadata <- fastaMetadata
      assign("aa_seq_tbl_final", aaSeqTblFinal, envir = .GlobalEnv)

      if (!is.null(sourceDir)) {
        scriptsAaSeqPath <- file.path(sourceDir, "aa_seq_tbl_final.RDS")
        saveRDS(aaSeqTblFinal, scriptsAaSeqPath)
        logger::log_info(sprintf("Saved aa_seq_tbl_final to scripts: %s", scriptsAaSeqPath))

        fastaMetadataPath <- file.path(sourceDir, "fasta_metadata.RDS")
        saveRDS(fastaMetadata, fastaMetadataPath)
        logger::log_info(sprintf("Saved FASTA metadata: %s", fastaMetadataPath))
        logger::log_info(sprintf(
          "FASTA Format: %s, Sequences: %d",
          fastaMetadata$fasta_format,
          fastaMetadata$num_sequences
        ))
      }

      logger::log_info(sprintf("FASTA processed successfully: %d sequences", nrow(aaSeqTblFinal)))
      shiny::showNotification(
        sprintf(
          "FASTA file processed: %d protein sequences available for accession cleanup (Format: %s)",
          nrow(aaSeqTblFinal),
          fastaMetadata$fasta_format
        ),
        type = "message",
        duration = 5
      )
    }, error = function(e) {
      logger::log_warn(paste("Error processing FASTA file:", e$message))
      logger::log_warn("Continuing without FASTA - accession cleanup will be skipped")
      workflowData$aa_seq_tbl_final <- NULL
      workflowData$fasta_metadata <- NULL
      shiny::showNotification(
        paste("Warning: Could not process FASTA file:", e$message),
        type = "warning",
        duration = 8
      )
    })
  }

  invisible(workflowData$aa_seq_tbl_final)
}

loadProtDesignImportedConfigAndTables <- function(
    workflowData,
    experimentPaths,
    designFile,
    dataClnFile,
    contrastFile,
    readConfig = readConfigFile,
    readTabular = function(path) vroom::vroom(path, show_col_types = FALSE),
    systemFileFn = system.file,
    fileExists = file.exists,
    fileCopy = file.copy,
    downloadFile = download.file,
    showNotification = shiny::showNotification,
    assignFn = assign
) {
  configPath <- file.path(experimentPaths$source_dir, "config.ini")

  if (!fileExists(configPath)) {
    logger::log_info("config.ini not found in project. Retrieving default config.")
    tryCatch({
      pkgConfig <- systemFileFn("config", "config.ini", package = "MultiScholaR")

      if (fileExists(pkgConfig) && pkgConfig != "") {
        logger::log_info(paste("Found default config.ini in package:", pkgConfig))
        fileCopy(pkgConfig, configPath)
        logger::log_info(paste("Default config.ini copied to:", configPath))
        showNotification("Using default config.ini (from package).", type = "message")
      } else {
        logger::log_info("Default config.ini not found in package, downloading from GitHub...")
        defaultConfigUrl <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/config/config.ini"
        downloadFile(defaultConfigUrl, destfile = configPath, quiet = TRUE)
        logger::log_info(paste("Default config.ini downloaded to:", configPath))
        showNotification("Using default config.ini (downloaded).", type = "message")
      }
    }, error = function(e_download) {
      msg <- paste("Failed to retrieve default config.ini:", e_download$message)
      logger::log_error(msg)
      showNotification(msg, type = "error", duration = 10)
      stop("Could not obtain a configuration file.")
    })
  }

  logger::log_info("Loading config.ini from project.")
  workflowData$config_list <- readConfig(file = configPath)
  assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
  logger::log_info("Created global config_list for updateConfigParameter compatibility")

  logger::log_info("Loading imported design/data/contrast artifacts from disk.")
  list(
    importedDesign = readTabular(designFile),
    importedDataCln = readTabular(dataClnFile),
    importedContrasts = readProtDesignImportedContrasts(contrastFile)
  )
}

resolveProtDesignImportArtifacts <- function(
    importPath,
    selectedFastaPath = NULL,
    listFiles = list.files,
    fileExists = file.exists,
    basenameFn = basename
) {
  fastaFiles <- listFiles(
    importPath,
    pattern = "\\.fasta$|\\.fa$|\\.faa$",
    ignore.case = TRUE,
    full.names = TRUE
  )

  fastaPath <- NULL
  if (!is.null(selectedFastaPath) && fileExists(selectedFastaPath)) {
    logger::log_info(sprintf("Using user-selected FASTA file: %s", basenameFn(selectedFastaPath)))
    fastaPath <- selectedFastaPath
  } else if (length(fastaFiles) > 0) {
    logger::log_info(sprintf("Auto-detected FASTA file: %s", basenameFn(fastaFiles[1])))
    fastaPath <- fastaFiles[1]
  } else {
    logger::log_info("No FASTA file provided - accession cleanup will be skipped if aa_seq_tbl_final.RDS not found")
  }

  designFile <- file.path(importPath, "design_matrix.tab")
  dataClnFile <- file.path(importPath, "data_cln.tab")
  contrastFile <- file.path(importPath, "contrast_strings.tab")

  if (!fileExists(designFile) || !fileExists(dataClnFile)) {
    return(list(
      ok = FALSE,
      errorMessage = "Import failed: 'design_matrix.tab' and/or 'data_cln.tab' not found in the selected directory."
    ))
  }

  list(
    ok = TRUE,
    importPath = importPath,
    fastaPath = fastaPath,
    designFile = designFile,
    dataClnFile = dataClnFile,
    contrastFile = contrastFile
  )
}

initializeProtDesignImportedWorkflowState <- function(
    workflowData,
    importedDesign,
    importedDataCln,
    importedContrasts,
    taxonId,
    organismName
) {
  workflowData$design_matrix <- importedDesign
  workflowData$data_cln <- importedDataCln
  workflowData$contrasts_tbl <- importedContrasts

  if (!is.null(importedContrasts)) {
    assign("contrasts_tbl", importedContrasts, envir = .GlobalEnv)
    logger::log_info("Saved contrasts_tbl to global environment for DE analysis.")
  }

  workflowData$taxon_id <- taxonId
  workflowData$organism_name <- organismName
  logger::log_info(sprintf(
    "Set organism info from import modal: %s (taxon: %d)",
    organismName,
    taxonId
  ))

  workflowType <- NULL
  if (!is.null(workflowData$config_list$globalParameters$workflow_type)) {
    workflowType <- workflowData$config_list$globalParameters$workflow_type
    logger::log_info(paste("Loaded workflow_type from config.ini:", workflowType))
  } else {
    logger::log_warn("workflow_type not found in config.ini. Attempting to detect from data structure.")

    hasPrecursorCols <- any(c("Precursor.Id", "Precursor.Charge", "Precursor.Quantity") %in% names(importedDataCln))
    hasPeptideCols <- any(c("Stripped.Sequence", "Modified.Sequence") %in% names(importedDataCln))

    if (hasPrecursorCols && hasPeptideCols) {
      workflowType <- "DIA"
      logger::log_info("Detected peptide-level data. Setting workflow_type to DIA.")
    } else if (workflowType == "LFQ" || workflowType == "TMT") {
      logger::log_info(sprintf("Keeping workflow_type as %s (protein-level workflow).", workflowType))
    } else {
      workflowType <- "TMT"
      logger::log_info("Detected protein-level data. Setting workflow_type to TMT.")
    }
  }

  workflowData$state_manager$setWorkflowType(workflowType)
  logger::log_info(paste("Workflow type set to:", workflowType))

  if (is.null(workflowData$column_mapping)) {
    logger::log_warn("column_mapping not found in workflow_data. Inferring from data structure.")

    if (workflowType == "TMT") {
      workflowData$column_mapping <- list(
        protein_col = "Protein.Ids",
        run_col = "Run",
        quantity_col = "Abundance"
      )
    } else if (workflowType == "LFQ") {
      workflowData$column_mapping <- list(
        protein_col = "Protein.Ids",
        run_col = "Run",
        quantity_col = "Intensity"
      )
    } else {
      workflowData$column_mapping <- list(
        protein_col = "Protein.Ids",
        run_col = "Run",
        quantity_col = "Precursor.Quantity"
      )
    }
    logger::log_info("Inferred column_mapping for workflow operations.")
  }

  workflowType
}

runProtDesignImportConfirmationFlow <- function(
    workflowData,
    experimentPaths,
    importArtifacts,
    importedArtifacts,
    taxonId,
    organismName,
    session,
    qcTrigger = NULL,
    successMessage = "Design imported successfully!",
    successNotificationId = "importing_design",
    hydrateFastaSidecar = hydrateProtDesignImportedFastaSidecar,
    hydrateUniprotSidecar = hydrateProtDesignImportedUniprotSidecar,
    initializeWorkflowState = initializeProtDesignImportedWorkflowState,
    buildStateCheckpointFn = buildProtDesignStateCheckpoint,
    completePostCheckpointFn = completeProtDesignPostCheckpoint
) {
  hydrateFastaSidecar(
    workflowData = workflowData,
    importPath = importArtifacts$importPath,
    sourceDir = experimentPaths$source_dir,
    resultsDir = experimentPaths$results_dir,
    fastaPath = importArtifacts$fastaPath,
    organismName = organismName
  )

  hydrateUniprotSidecar(
    workflowData = workflowData,
    importPath = importArtifacts$importPath,
    sourceDir = experimentPaths$source_dir
  )

  workflowType <- initializeWorkflowState(
    workflowData = workflowData,
    importedDesign = importedArtifacts$importedDesign,
    importedDataCln = importedArtifacts$importedDataCln,
    importedContrasts = importedArtifacts$importedContrasts,
    taxonId = taxonId,
    organismName = organismName
  )

  stateName <- buildStateCheckpointFn(
    workflowData = workflowData,
    workflowType = workflowType,
    actionLabel = "Import",
    validateColumnMapping = TRUE
  )

  logger::log_info(sprintf("Import: S4 object saved to R6 state manager as '%s'", stateName))
  logger::log_info("Import: This state is now ACTIVE - QC modules will read from it")
  logger::log_info("Import: User can proceed to QC -> Accession Cleanup")

  completePostCheckpointFn(
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    session = session,
    qcTrigger = qcTrigger,
    successMessage = successMessage,
    successNotificationId = successNotificationId
  )

  invisible(list(
    workflowType = workflowType,
    stateName = stateName
  ))
}

runProtDesignImportObserverShell <- function(
    workflowData,
    experimentPaths,
    importPath,
    selectedFastaPath = NULL,
    taxonId,
    organismName,
    session,
    qcTrigger = NULL,
    resolveImportArtifacts = resolveProtDesignImportArtifacts,
    loadImportedArtifacts = loadProtDesignImportedConfigAndTables,
    runImportConfirmationFlow = runProtDesignImportConfirmationFlow,
    showNotification = shiny::showNotification,
    removeNotification = shiny::removeNotification
) {
  importArtifacts <- resolveImportArtifacts(
    importPath = importPath,
    selectedFastaPath = selectedFastaPath
  )

  if (!isTRUE(importArtifacts$ok)) {
    logger::log_error(importArtifacts$errorMessage)
    showNotification(importArtifacts$errorMessage, type = "error", duration = 10)
    return(invisible(list(
      ok = FALSE,
      stage = "resolve",
      errorMessage = importArtifacts$errorMessage
    )))
  }

  showNotification("Importing design files...", id = "importing_design", duration = NULL)

  tryCatch({
    importedArtifacts <- loadImportedArtifacts(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      designFile = importArtifacts$designFile,
      dataClnFile = importArtifacts$dataClnFile,
      contrastFile = importArtifacts$contrastFile
    )

    flowResult <- runImportConfirmationFlow(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      importArtifacts = importArtifacts,
      importedArtifacts = importedArtifacts,
      taxonId = taxonId,
      organismName = organismName,
      session = session,
      qcTrigger = qcTrigger
    )

    invisible(list(
      ok = TRUE,
      importArtifacts = importArtifacts,
      importedArtifacts = importedArtifacts,
      flowResult = flowResult
    ))
  }, error = function(e) {
    msg <- paste("Error during import:", e$message)
    logger::log_error(msg)
    showNotification(msg, type = "error", duration = 15)
    removeNotification("importing_design")
    invisible(list(
      ok = FALSE,
      stage = "run",
      errorMessage = msg
    ))
  })
}

registerProtDesignImportModalShell <- function(
    input,
    output,
    session,
    resolvedVolumes,
    importFastaPath
) {
  shiny::observeEvent(input$import_fasta_file, {
    if (!is.null(input$import_fasta_file) && !is.integer(input$import_fasta_file)) {
      fileInfo <- shinyFiles::parseFilePaths(resolvedVolumes, input$import_fasta_file)
      if (nrow(fileInfo) > 0) {
        selectedPath <- as.character(fileInfo$datapath[1])
        importFastaPath(selectedPath)
        output$import_fasta_file_path <- shiny::renderText(selectedPath)
      }
    }
  })

  shiny::observeEvent(input$show_import_modal, {
    importFastaPath(NULL)

    ns <- session$ns
    shiny::showModal(shiny::modalDialog(
      title = "Import Existing Design Matrix",
      shiny::p("Select the folder containing 'design_matrix.tab' and 'data_cln.tab' files."),
      shinyFiles::shinyDirButton(ns("import_dir"), "Select Folder", "Choose a directory"),
      shiny::verbatimTextOutput(ns("import_dir_path"), placeholder = TRUE),
      shiny::hr(),
      shiny::h5("FASTA File (Optional)"),
      shiny::p("Required for protein accession cleanup. Will be auto-detected from import folder if available."),
      shinyFiles::shinyFilesButton(
        ns("import_fasta_file"),
        "Select FASTA file (if not in import folder)",
        "Choose File",
        multiple = FALSE,
        icon = shiny::icon("file")
      ),
      shiny::br(),
      shiny::verbatimTextOutput(ns("import_fasta_file_path"), placeholder = TRUE),
      shiny::textOutput(ns("fasta_detection_status")),
      shiny::hr(),
      shiny::h5("Organism Information"),
      shiny::p("Required for UniProt annotation lookups:"),
      shiny::numericInput(ns("import_taxon_id"), "Taxonomy ID:", value = 9606, min = 1),
      shiny::helpText("e.g., 9606 for Homo sapiens, 10090 for Mus musculus"),
      shiny::textInput(ns("import_organism_name"), "Organism Name:", value = "Homo sapiens"),
      footer = shiny::tagList(
        shiny::modalButton("Cancel"),
        shiny::actionButton(ns("confirm_import"), "Import", class = "btn-primary")
      )
    ))
  })

  output$import_dir_path <- shiny::renderText({
    shiny::req(input$import_dir)
    shinyFiles::parseDirPath(resolvedVolumes, input$import_dir)
  })

  output$fasta_detection_status <- shiny::renderText({
    shiny::req(input$import_dir)
    importPath <- shinyFiles::parseDirPath(resolvedVolumes, input$import_dir)

    if (length(importPath) > 0 && dir.exists(importPath)) {
      fastaFiles <- list.files(importPath, pattern = "\\.fasta$|\\.fa$|\\.faa$", ignore.case = TRUE)
      if (length(fastaFiles) > 0) {
        paste("[OK] Found FASTA file:", fastaFiles[1])
      } else {
        "[WARNING] No FASTA file detected in folder. Upload one above if needed for accession cleanup."
      }
    } else {
      ""
    }
  })

  invisible(importFastaPath)
}

registerProtDesignImportConfirmationObserver <- function(
    input,
    resolvedVolumes,
    importFastaPath,
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    parseDirPathFn = shinyFiles::parseDirPath,
    removeModalFn = shiny::removeModal,
    runImportObserverShell = runProtDesignImportObserverShell
) {
  observeEventFn(input$confirm_import, {
    reqFn(input$import_dir)

    importPath <- parseDirPathFn(resolvedVolumes, input$import_dir)
    reqFn(importPath)

    removeModalFn()

    runImportObserverShell(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      importPath = importPath,
      selectedFastaPath = importFastaPath(),
      taxonId = input$import_taxon_id,
      organismName = input$import_organism_name,
      session = session,
      qcTrigger = qcTrigger
    )
  })
}

initializeProtDesignImportBootstrap <- function(
    input,
    session,
    experimentPaths,
    volumes = NULL,
    dirChooseFn = shinyFiles::shinyDirChoose,
    fileChooseFn = shinyFiles::shinyFileChoose,
    reactiveValFn = shiny::reactiveVal,
    isolateFn = shiny::isolate,
    getVolumesFn = shinyFiles::getVolumes,
    dirExistsFn = dir.exists,
    logInfo = logger::log_info
) {
  resolvedVolumes <- isolateFn({
    baseVolumes <- if (is.function(volumes)) {
      volumes()
    } else if (is.null(volumes)) {
      getVolumesFn()()
    } else {
      volumes
    }

    if (!is.null(experimentPaths) &&
        !is.null(experimentPaths$base_dir) &&
        dirExistsFn(experimentPaths$base_dir)) {
      enhancedVolumes <- c("Project Base Dir" = experimentPaths$base_dir, baseVolumes)
      logInfo(paste("Added base_dir to volumes for easier navigation:", experimentPaths$base_dir))
      enhancedVolumes
    } else {
      baseVolumes
    }
  })

  dirChooseFn(input, "import_dir", roots = resolvedVolumes, session = session)
  fileChooseFn(
    input,
    "import_fasta_file",
    roots = resolvedVolumes,
    session = session,
    filetypes = c("fasta", "fa", "faa")
  )

  list(
    resolvedVolumes = resolvedVolumes,
    importFastaPath = reactiveValFn(NULL)
  )
}

