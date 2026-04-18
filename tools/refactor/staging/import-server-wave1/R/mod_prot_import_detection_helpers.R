readProtImportHeaders <- function(
    filePath
    , unzipList = function(path) unzip(path, list = TRUE)
    , unzipFiles = unzip
    , fileExt = tools::file_ext
    , requireNamespaceFn = requireNamespace
    , readExcel = readxl::read_excel
    , readParquet = arrow::read_parquet
    , readLinesFn = readLines
    , tempfileFn = tempfile
    , dirCreate = dir.create
    , unlinkFn = unlink
    , basenameFn = basename
) {
  headers <- NULL

  if (tolower(fileExt(filePath)) == "zip") {
    zip_contents <- unzipList(filePath)
    first_data_file <- zip_contents$Name[grep("\\.(xlsx|tsv|csv)$", zip_contents$Name, ignore.case = TRUE)[1]]

    if (is.na(first_data_file)) {
      stop("No data files (.xlsx, .tsv, .csv) found inside the ZIP archive.")
    }

    if (tolower(fileExt(first_data_file)) == "xlsx") {
      if (!requireNamespaceFn("readxl", quietly = TRUE)) {
        stop("Package 'readxl' needed for .xlsx files.")
      }

      temp_dir <- tempfileFn()
      dirCreate(temp_dir)
      on.exit(unlinkFn(temp_dir, recursive = TRUE), add = TRUE)

      unzipFiles(filePath, files = first_data_file, exdir = temp_dir, junkpaths = TRUE)
      unzipped_file_path <- file.path(temp_dir, basenameFn(first_data_file))
      headers <- names(readExcel(unzipped_file_path, n_max = 0))
    } else {
      con <- unz(filePath, first_data_file)
      on.exit(close(con), add = TRUE)
      headers <- strsplit(readLinesFn(con, n = 1), "\t|,")[[1]]
    }
  } else if (grepl("\\.parquet$", filePath, ignore.case = TRUE)) {
    headers <- names(readParquet(filePath, col_select = NULL))
  } else {
    preview_lines <- readLinesFn(filePath, n = 10)
    headers <- strsplit(preview_lines[1], "\t|,")[[1]]
  }

  if (is.null(headers)) {
    stop("Could not read headers from the provided file.")
  }

  headers
}

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

