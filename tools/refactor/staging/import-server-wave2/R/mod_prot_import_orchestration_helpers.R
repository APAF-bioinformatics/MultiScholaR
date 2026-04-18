announceProtImportStart <- function(
    searchResultsPath
    , fastaPath
    , format
    , messageFn = message
) {
  messageFn("========================================")
  messageFn("[mod_prot_import] Starting data import process")
  messageFn("========================================")
  messageFn(sprintf("[mod_prot_import] Search results: %s", searchResultsPath))
  messageFn(sprintf("[mod_prot_import] FASTA file: %s", fastaPath))
  messageFn(sprintf("[mod_prot_import] Detected format: %s", format))

  invisible(list(
    searchResultsPath = searchResultsPath
    , fastaPath = fastaPath
    , format = format
  ))
}

recordProtImportImportedData <- function(
    dataImportResult
    , captureCheckpoint = .capture_checkpoint
    , logInfo = logger::log_info
    , messageFn = message
) {
  logInfo(sprintf("Data imported successfully. Rows: %d", nrow(dataImportResult$data)))
  captureCheckpoint(dataImportResult, "cp01", "raw_imported")
  messageFn(sprintf("[mod_prot_import] Data imported: %d rows", nrow(dataImportResult$data)))

  invisible(dataImportResult)
}

recordProtImportSetupLog <- function(
    workflowData
    , dataImportResult
    , format
    , useShinyFiles
    , searchResultsPath
    , fastaPath
    , searchResultsStandard = NULL
    , fastaFileStandard = NULL
    , taxonId
    , organismName
    , mixedSpeciesFasta = FALSE
    , resolveFilename = resolveProtImportInputFilename
    , finalizeSetupState = finalizeProtImportSetupState
    , logInfo = logger::log_info
) {
  logInfo("Creating processing log...")

  search_filename <- resolveFilename(
    useShinyFiles = useShinyFiles
    , resolvedPath = searchResultsPath
    , standardInput = searchResultsStandard
  )

  fasta_filename <- resolveFilename(
    useShinyFiles = useShinyFiles
    , resolvedPath = fastaPath
    , standardInput = fastaFileStandard
  )

  finalizeSetupState(
    workflowData = workflowData
    , dataImportResult = dataImportResult
    , format = format
    , searchFilename = search_filename
    , fastaFilename = fasta_filename
    , taxonId = taxonId
    , organismName = organismName
    , mixedSpeciesFasta = mixedSpeciesFasta
  )

  logInfo("Processing log created successfully")

  invisible(list(
    searchFilename = search_filename
    , fastaFilename = fasta_filename
  ))
}

startProtImportProcessing <- function(
    searchResultsPath
    , fastaPath
    , format
    , localData
    , ns
    , requireValue = shiny::req
    , announceStart = announceProtImportStart
    , buildProcessingModal = buildProtImportProcessingModal
    , showModal = shiny::showModal
) {
  requireValue(searchResultsPath)
  requireValue(fastaPath)

  localData$processing <- TRUE

  announceStart(
    searchResultsPath = searchResultsPath
    , fastaPath = fastaPath
    , format = format
  )

  showModal(buildProcessingModal(ns))

  invisible(list(
    searchResultsPath = searchResultsPath
    , fastaPath = fastaPath
    , format = format
  ))
}

loadProtImportConfigurationResources <- function(
    workflowData
    , experimentPaths
    , useShinyFiles
    , configFilePath = NULL
    , configFileStandard = NULL
    , uniprotMappingFile = NULL
    , uniprotMappingStandard = NULL
    , uniparcMappingFile = NULL
    , uniparcMappingStandard = NULL
    , updateStatus = updateProtImportProcessingStatus
    , messageFn = message
    , resolveUploadPath = resolveProtImportOptionalUploadPath
    , loadConfiguration = loadProtImportConfiguration
    , storeConfiguration = storeProtImportConfiguration
    , loadOptionalMappings = loadProtImportOptionalMappings
    , readConfig = readConfigFile
    , showNotification = shiny::showNotification
    , logInfo = logger::log_info
    , logError = logger::log_error
    , getDefaultConfig = getDefaultProteomicsConfig
    , fileExists = file.exists
    , fileCopy = file.copy
    , downloadFile = download.file
    , systemFileFn = system.file
    , assignFn = assign
    , assignEnv = .GlobalEnv
    , readOptionalMapping = readProtImportOptionalMapping
) {
  updateStatus("Loading configuration...")
  messageFn("[mod_prot_import] Loading configuration")

  config_path <- resolveUploadPath(
    useShinyFiles = useShinyFiles
    , shinyPath = configFilePath
    , standardInput = configFileStandard
  )

  config_list <- loadConfiguration(
    configPath = config_path
    , experimentPaths = experimentPaths
    , readConfig = readConfig
    , showNotification = showNotification
    , logInfo = logInfo
    , logError = logError
    , getDefaultConfig = getDefaultConfig
    , fileExists = fileExists
    , fileCopy = fileCopy
    , downloadFile = downloadFile
    , systemFileFn = systemFileFn
  )

  storeConfiguration(
    workflowData = workflowData
    , configList = config_list
    , assignFn = assignFn
    , assignEnv = assignEnv
    , logInfo = logInfo
  )

  uniprot_path <- resolveUploadPath(
    useShinyFiles = useShinyFiles
    , shinyPath = uniprotMappingFile
    , standardInput = uniprotMappingStandard
  )

  uniparc_path <- resolveUploadPath(
    useShinyFiles = useShinyFiles
    , shinyPath = uniparcMappingFile
    , standardInput = uniparcMappingStandard
  )

  loadOptionalMappings(
    workflowData = workflowData
    , uniprotPath = uniprot_path
    , uniparcPath = uniparc_path
    , readOptionalMapping = readOptionalMapping
    , logInfo = logInfo
    , logError = logError
  )

  invisible(list(
    configPath = config_path
    , configList = config_list
    , uniprotPath = uniprot_path
    , uniparcPath = uniparc_path
  ))
}

finalizeProtImportProcessing <- function(
    workflowData
    , localData
    , dataImportResult
    , format
    , updateStatus = updateProtImportProcessingStatus
    , messageFn = message
    , completeSuccessState = completeProtImportSuccessState
) {
  updateStatus("Finalizing import...")
  messageFn("[mod_prot_import] Finalizing import")

  updated_status <- completeSuccessState(
    workflowData = workflowData
    , localData = localData
    , dataImportResult = dataImportResult
    , format = format
    , messageFn = messageFn
  )

  invisible(list(
    updatedStatus = updated_status
    , format = format
  ))
}

readProtImportDataWithStatus <- function(
    format
    , searchResultsPath
    , input
    , updateStatus = updateProtImportProcessingStatus
    , logInfo = logger::log_info
    , importDataByFormat = importProtImportDataByFormat
    , recordImportedData = recordProtImportImportedData
    , importDiann = importDIANNData
    , importSpectronaut = importSpectronautData
    , importFragpipe = importFragPipeData
    , importMaxquant = importMaxQuantData
    , importPdTmt = importProteomeDiscovererTMTData
    , logError = logger::log_error
    , captureCheckpoint = .capture_checkpoint
    , messageFn = message
) {
  updateStatus(sprintf("Reading %s data...", toupper(format)))
  logInfo(sprintf("Reading %s data from %s", format, searchResultsPath))

  data_import_result <- importDataByFormat(
    format = format
    , searchResultsPath = searchResultsPath
    , input = input
    , importDiann = importDiann
    , importSpectronaut = importSpectronaut
    , importFragpipe = importFragpipe
    , importMaxquant = importMaxquant
    , importPdTmt = importPdTmt
    , logError = logError
  )

  recordImportedData(
    dataImportResult = data_import_result
    , captureCheckpoint = captureCheckpoint
    , logInfo = logInfo
    , messageFn = messageFn
  )

  invisible(data_import_result)
}

applyProtImportWorkflowWithStatus <- function(
    workflowData
    , dataImportResult
    , format
    , fastaPath
    , organismName
    , experimentPaths
    , sanitizeNames = FALSE
    , updateStatus = updateProtImportProcessingStatus
    , applyImportResult = applyProtImportResultToWorkflow
    , processFastaData = processProtImportFastaData
    , logInfo = logger::log_info
    , logWarn = logger::log_warn
    , showNotification = shiny::showNotification
    , sanitizeRunNames = sanitizeProtImportRunNames
    , processFasta = processFastaFile
) {
  updateStatus("Storing imported data...")

  if (isTRUE(sanitizeNames)) {
    updateStatus("Sanitizing sample names...")
  }

  applyImportResult(
    workflowData = workflowData
    , dataImportResult = dataImportResult
    , format = format
    , fastaPath = fastaPath
    , sanitizeNames = sanitizeNames
    , logInfo = logInfo
    , showNotification = showNotification
    , sanitizeRunNames = sanitizeRunNames
  )

  updateStatus("Processing FASTA file...")

  processFastaData(
    workflowData = workflowData
    , fastaPath = fastaPath
    , organismName = organismName
    , experimentPaths = experimentPaths
    , processFasta = processFasta
    , logInfo = logInfo
    , logWarn = logWarn
  )

  invisible(list(
    sanitizeNames = sanitizeNames
    , fastaPath = fastaPath
  ))
}

runProtImportMixedSpeciesAnalysisIfNeeded <- function(
    workflowData
    , localData
    , dataImportResult
    , fastaPath
    , session
    , mixedSpeciesFasta = FALSE
    , analyzeMixedSpeciesData = analyzeProtImportMixedSpeciesData
    , updateStatus = updateProtImportProcessingStatus
    , messageFn = message
    , logInfo = logger::log_info
    , logWarn = logger::log_warn
    , showNotification = shiny::showNotification
    , extractOrganisms = extractOrganismsFromFasta
    , analyzeDistribution = analyzeOrganismDistribution
    , buildSelectionModal = buildProtImportOrganismSelectionModal
    , showModal = shiny::showModal
) {
  if (!isTRUE(mixedSpeciesFasta) || is.null(workflowData$aa_seq_tbl_final)) {
    return(FALSE)
  }

  analyzeMixedSpeciesData(
    workflowData = workflowData
    , localData = localData
    , dataImportResult = dataImportResult
    , fastaPath = fastaPath
    , session = session
    , updateStatus = updateStatus
    , messageFn = messageFn
    , logInfo = logInfo
    , logWarn = logWarn
    , showNotification = showNotification
    , extractOrganisms = extractOrganisms
    , analyzeDistribution = analyzeDistribution
    , buildSelectionModal = buildSelectionModal
    , showModal = showModal
  )
}

runProtImportProcessingSafely <- function(
    runProcessing
    , workflowData
    , localData
    , logError = logger::log_error
    , resetWorkflowState = resetProtImportWorkflowStateOnError
) {
  tryCatch(
    runProcessing()
    , error = function(e) {
      resetWorkflowState(
        workflowData = workflowData
        , localData = localData
        , errorMessage = e$message
        , logError = logError
      )
    }
  )
}

