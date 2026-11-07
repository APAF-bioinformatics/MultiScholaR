resolveProtImportDetectionFilename <- function(
    useShinyFiles
    , filePath
    , searchResultsStandard = NULL
    , basenameFn = basename
) {
  if (isTRUE(useShinyFiles)) {
    return(basenameFn(filePath))
  }

  if (is.null(searchResultsStandard) || is.null(searchResultsStandard$name)) {
    return(NULL)
  }

  searchResultsStandard$name
}

applyProtImportDetectedFormat <- function(
    localData
    , formatInfo
    , logInfo = logger::log_info
) {
  localData$detected_format <- formatInfo$format
  localData$format_confidence <- formatInfo$confidence

  logInfo(sprintf("Detected format: %s (confidence: %s)", formatInfo$format, formatInfo$confidence))

  invisible(formatInfo)
}

resetProtImportFormatDetectionState <- function(
    localData
    , errorMessage
    , logError = logger::log_error
) {
  logError(paste("Error detecting file format:", errorMessage))
  localData$detected_format <- "unknown"
  localData$format_confidence <- 0

  invisible(list(
    format = "unknown"
    , confidence = 0
    , errorMessage = errorMessage
  ))
}

runProtImportFormatDetection <- function(
    filePath
    , useShinyFiles
    , localData
    , searchResultsStandard = NULL
    , readHeaders = readProtImportHeaders
    , resolveFilename = resolveProtImportDetectionFilename
    , detectFormat = detectProteomicsFormat
    , applyDetectedFormat = applyProtImportDetectedFormat
    , resetDetectionState = resetProtImportFormatDetectionState
    , logInfo = logger::log_info
    , logError = logger::log_error
) {
  tryCatch({
    headers <- readHeaders(filePath)
    filename <- resolveFilename(
      useShinyFiles = useShinyFiles
      , filePath = filePath
      , searchResultsStandard = searchResultsStandard
    )

    format_info <- detectFormat(
      headers = headers
      , filename = filename
    )

    applyDetectedFormat(
      localData = localData
      , formatInfo = format_info
      , logInfo = logInfo
    )
  }, error = function(e) {
    resetDetectionState(
      localData = localData
      , errorMessage = e$message
      , logError = logError
    )
  })
}

resolveActiveProtImportFormat <- function(formatOverride, detectedFormat) {
  if (identical(formatOverride, "auto")) {
    return(detectedFormat)
  }

  formatOverride
}
