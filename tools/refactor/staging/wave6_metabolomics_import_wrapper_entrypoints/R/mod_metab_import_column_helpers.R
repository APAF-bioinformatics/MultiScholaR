resolveMetabImportColumnName <- function(headers, columnName) {
  if (is.null(columnName) || !nzchar(columnName) || is.null(headers)) {
    return(columnName)
  }

  headersLower <- tolower(headers)
  columnLower <- tolower(columnName)

  if (columnLower %in% headersLower) {
    return(headers[which(headersLower == columnLower)[1]])
  }

  columnName
}

resolveMetabImportSampleColumns <- function(
    assayData,
    vendorFormat,
    sampleColsPattern = NULL,
    importResult = NULL
) {
  if (identical(vendorFormat, "custom") && !is.null(sampleColsPattern) && nzchar(sampleColsPattern)) {
    allCols <- names(assayData)
    matched <- allCols[grepl(sampleColsPattern, allCols, ignore.case = TRUE)]

    if (length(matched) > 0) {
      return(matched)
    }
  }

  if (!is.null(importResult)) {
    return(importResult$sample_columns)
  }

  names(assayData)[sapply(assayData, is.numeric)]
}

