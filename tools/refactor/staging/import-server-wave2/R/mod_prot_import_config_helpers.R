loadProtImportConfiguration <- function(
    configPath = NULL
    , experimentPaths
    , readConfig = readConfigFile
    , showNotification = shiny::showNotification
    , logInfo = logger::log_info
    , logError = logger::log_error
    , getDefaultConfig = getDefaultProteomicsConfig
    , fileExists = file.exists
    , fileCopy = file.copy
    , downloadFile = download.file
    , systemFileFn = system.file
) {
  if (!is.null(configPath)) {
    logInfo(paste("Reading configuration from", configPath))
    return(readConfig(file = configPath))
  }

  logInfo("Using default configuration")
  default_config_path <- file.path(experimentPaths$source_dir, "config.ini")

  if (fileExists(default_config_path)) {
    return(readConfig(file = default_config_path))
  }

  logInfo("config.ini not found in project. Retrieving default config.")
  tryCatch({
    pkg_config <- systemFileFn("config", "config.ini", package = "MultiScholaR")

    if (fileExists(pkg_config) && pkg_config != "") {
      logInfo(paste("Found default config.ini in package:", pkg_config))
      fileCopy(pkg_config, default_config_path)
      logInfo(paste("Default config.ini copied to:", default_config_path))
      showNotification("Copied default config.ini to scripts directory.", type = "message")
      readConfig(file = default_config_path)
    } else {
      logInfo("Default config.ini not found in package, downloading from GitHub...")
      default_config_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/config/config.ini"
      downloadFile(default_config_url, destfile = default_config_path, quiet = TRUE)
      logInfo(paste("Default config.ini downloaded to:", default_config_path))
      showNotification("Downloaded default config.ini to scripts directory.", type = "message")
      readConfig(file = default_config_path)
    }
  }, error = function(e_download) {
    msg <- paste("Failed to retrieve default config.ini:", e_download$message)
    logError(msg)
    showNotification(msg, type = "warning", duration = 10)
    logInfo("Using minimal fallback configuration")
    getDefaultConfig()
  })
}

storeProtImportConfiguration <- function(
    workflowData
    , configList
    , assignFn = assign
    , assignEnv = .GlobalEnv
    , logInfo = logger::log_info
) {
  workflowData$config_list <- configList
  assignFn("config_list", configList, envir = assignEnv)
  logInfo("Created global config_list for compatibility")

  invisible(configList)
}

readProtImportOptionalMapping <- function(
    mappingPath = NULL
    , mappingLabel
    , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
    , logInfo = logger::log_info
    , logError = logger::log_error
) {
  if (is.null(mappingPath)) {
    return(NULL)
  }

  logInfo(paste("Reading", mappingLabel, "from", mappingPath))

  tryCatch({
    readMapping(mappingPath)
  }, error = function(e) {
    logError(paste("Failed to read", mappingLabel, "file:", e$message))
    stop("Failed to read ", mappingLabel, " file: ", e$message)
  })
}

resolveProtImportOptionalUploadPath <- function(
    useShinyFiles
    , shinyPath = NULL
    , standardInput = NULL
) {
  if (isTRUE(useShinyFiles)) {
    return(shinyPath)
  }

  if (is.null(standardInput)) {
    return(NULL)
  }

  standardInput$datapath
}

resolveProtImportInputFilename <- function(
    useShinyFiles
    , resolvedPath
    , standardInput = NULL
) {
  if (isTRUE(useShinyFiles)) {
    if (is.null(resolvedPath)) {
      return(NULL)
    }

    return(basename(resolvedPath))
  }

  if (is.null(standardInput)) {
    return(NULL)
  }

  standardInput$name
}

loadProtImportOptionalMappings <- function(
    workflowData
    , uniprotPath = NULL
    , uniparcPath = NULL
    , readOptionalMapping = readProtImportOptionalMapping
    , logInfo = logger::log_info
    , logError = logger::log_error
) {
  if (!is.null(uniprotPath)) {
    workflowData$uniprot_mapping <- readOptionalMapping(
      mappingPath = uniprotPath
      , mappingLabel = "UniProt mapping"
      , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
      , logInfo = logInfo
      , logError = logError
    )
  }

  if (!is.null(uniparcPath)) {
    workflowData$uniparc_mapping <- readOptionalMapping(
      mappingPath = uniparcPath
      , mappingLabel = "UniParc mapping"
      , readMapping = function(path) vroom::vroom(path, show_col_types = FALSE)
      , logInfo = logInfo
      , logError = logError
    )
  }

  invisible(list(
    uniprotMapping = workflowData$uniprot_mapping
    , uniparcMapping = workflowData$uniparc_mapping
  ))
}

