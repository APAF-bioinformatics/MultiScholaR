storeProtImportConfiguration <- function(
    workflowData
    , configList
    , logInfo = logger::log_info
) {
  workflowData$config_list <- configList
  logInfo("Stored config_list on workflow_data")

  invisible(configList)
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
