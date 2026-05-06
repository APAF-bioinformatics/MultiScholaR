buildLipidImportFormatDetectionStatus <- function(detectedFormat, formatConfidence) {
  confidence_pct <- round(formatConfidence * 100)
  color_class <- if (confidence_pct >= 70) "success" else if (confidence_pct >= 40) "warning" else "danger"

  format_display <- switch(detectedFormat,
    "msdial" = "MS-DIAL",
    "lipidsearch" = "LipidSearch",
    "progenesis" = "Progenesis QI",
    "xcms" = "XCMS",
    "compound_discoverer" = "Compound Discoverer",
    "Unknown"
  )

  shiny::tags$div(
    class = paste("alert", paste0("alert-", color_class)),
    shiny::tags$strong("Detected format: "),
    format_display,
    shiny::tags$br(),
    shiny::tags$small(sprintf("Confidence: %d%%", confidence_pct))
  )
}

buildLipidImportColumnValidationStatus <- function(
  assayData,
  columnName,
  successMode = c("unique_ids", "found_unique_ids", "ok", "found"),
  emptyMode = c("optional", "prompt"),
  allowCaseInsensitive = FALSE
) {
  successMode <- match.arg(successMode)
  emptyMode <- match.arg(emptyMode)

  if (is.null(columnName) || !nzchar(columnName)) {
    if (identical(emptyMode, "prompt")) {
      return(shiny::tags$span(
        shiny::icon("question-circle", style = "color: gray;"),
        " Enter column name"
      ))
    }

    return(shiny::tags$span(
      shiny::icon("minus-circle", style = "color: gray;"),
      " Optional"
    ))
  }

  actual_column <- columnName
  if (allowCaseInsensitive) {
    headers <- names(assayData)
    header_index <- match(tolower(columnName), tolower(headers))
    actual_column <- if (is.na(header_index)) NULL else headers[[header_index]]
  } else if (!columnName %in% names(assayData)) {
    actual_column <- NULL
  }

  if (is.null(actual_column)) {
    return(shiny::tags$span(
      shiny::icon("times-circle", style = "color: red;"),
      " Column not found"
    ))
  }

  success_text <- switch(successMode,
    "unique_ids" = sprintf(" %d unique IDs", length(unique(assayData[[actual_column]]))),
    "found_unique_ids" = sprintf(" Found: %d unique IDs", length(unique(assayData[[actual_column]]))),
    "ok" = " OK",
    "found" = " Found"
  )

  shiny::tags$span(
    shiny::icon("check-circle", style = "color: green;"),
    success_text
  )
}

resolveLipidImportEffectiveColumn <- function(
  assayData,
  vendorFormat,
  selectedColumn,
  customColumn
) {
  if (!identical(vendorFormat, "custom")) {
    return(selectedColumn)
  }

  if (is.null(customColumn) || !nzchar(customColumn)) {
    return(customColumn)
  }

  headers <- names(assayData)
  header_index <- match(tolower(customColumn), tolower(headers))

  if (is.na(header_index)) {
    return(customColumn)
  }

  headers[[header_index]]
}

resolveLipidImportSampleColumns <- function(
  assayData,
  assayImportResult = NULL,
  vendorFormat,
  sampleColsPattern = "",
  excludeNormalized = FALSE,
  exclusionPattern = "_norm(ali[sz]ed|allized)?$|normalized"
) {
  if (identical(vendorFormat, "custom") &&
      !is.null(sampleColsPattern) &&
      nzchar(sampleColsPattern)) {
    all_cols <- names(assayData)
    matched <- all_cols[grepl(sampleColsPattern, all_cols, ignore.case = TRUE)]

    if (length(matched) > 0) {
      return(matched)
    }
  }

  if (!is.null(assayImportResult)) {
    sample_cols <- assayImportResult$sample_columns
  } else {
    sample_cols <- names(assayData)[sapply(assayData, is.numeric)]
  }

  if (isTRUE(excludeNormalized)) {
    sample_cols <- sample_cols[!grepl(exclusionPattern, sample_cols, ignore.case = TRUE)]
  }

  sample_cols
}

formatLipidImportColumnPreviewText <- function(columnNames, maxColumns = Inf) {
  if (length(columnNames) > maxColumns) {
    return(paste(
      paste(head(columnNames, maxColumns), collapse = ", "),
      sprintf("... and %d more", length(columnNames) - maxColumns)
    ))
  }

  paste(columnNames, collapse = ", ")
}

buildLipidImportValidationSummary <- function(validation) {
  if (validation$valid) {
    return(shiny::tagList(
      shiny::tags$div(
        class = "alert alert-success",
        shiny::icon("check-circle"),
        shiny::tags$strong(" Validation Passed")
      ),
      shiny::tags$ul(
        shiny::tags$li(sprintf("Lipids: %d", validation$summary$n_lipids)),
        shiny::tags$li(sprintf("Samples: %d", validation$summary$n_samples)),
        shiny::tags$li(sprintf("Missing values: %.1f%%", validation$summary$pct_missing))
      ),
      if (length(validation$warnings) > 0) {
        shiny::tags$div(
          class = "alert alert-warning",
          shiny::icon("exclamation-triangle"),
          " Warnings:",
          shiny::tags$ul(
            lapply(validation$warnings, shiny::tags$li)
          )
        )
      }
    ))
  }

  shiny::tags$div(
    class = "alert alert-danger",
    shiny::icon("times-circle"),
    shiny::tags$strong(" Validation Failed"),
    shiny::tags$ul(
      lapply(validation$errors, shiny::tags$li)
    )
  )
}

buildLipidImportStatusDisplay <- function(workflowData) {
  if (is.null(workflowData$tab_status$setup_import) ||
    !identical(workflowData$tab_status$setup_import, "complete")) {
    return(NULL)
  }

  log_info <- workflowData$processing_log$setup_import

  shiny::tags$div(
    class = "alert alert-success",
    shiny::icon("check-circle"),
    shiny::tags$strong(" Import Complete"),
    shiny::tags$br(),
    sprintf(
      "Format: %s | Assays: %d | Samples: %d",
      toupper(log_info$detected_format),
      log_info$n_assays,
      log_info$n_samples
    )
  )
}

emitLipidImportModuleServerEntryDiagnostics <- function(
  volumesIsNull = NULL,
  insideModuleServer = FALSE,
  emitMessage = message
) {
  if (insideModuleServer) {
    emitMessage("   mod_lipid_import_server: Inside moduleServer function")
    return(invisible(NULL))
  }

  emitMessage("--- Entering mod_lipid_import_server ---")
  emitMessage(sprintf(
    "   mod_lipid_import_server: volumes param is NULL = %s",
    volumesIsNull
  ))

  invisible(NULL)
}

probeLipidImportShinyFilesAvailability <- function(
  requireNamespaceFn = requireNamespace,
  isTestModeFn = is_test_mode,
  emitMessage = message
) {
  useShinyFiles <- requireNamespaceFn("shinyFiles", quietly = TRUE) && !isTRUE(isTestModeFn())
  emitMessage(sprintf(
    "   mod_lipid_import_server: shinyFiles available = %s",
    useShinyFiles
  ))
  useShinyFiles
}

initializeLipidImportLocalData <- function(
  reactiveValues = shiny::reactiveValues
) {
  reactiveValues(
    assay1_file = NULL,
    assay1_data = NULL,
    assay1_import_result = NULL,
    assay2_file = NULL,
    assay2_data = NULL,
    assay2_import_result = NULL,
    detected_format = NULL,
    format_confidence = NULL,
    all_headers = NULL
  )
}

applyLipidImportResultToWorkflow <- function(
  workflowData,
  assayList,
  dataFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  workflowType = "lipidomics_standard",
  logInfo = logger::log_info
) {
  workflowData$data_tbl <- assayList
  workflowData$data_format <- dataFormat
  workflowData$data_type <- "lipid"
  workflowData$column_mapping <- buildLipidImportColumnMapping(
    lipidIdCol = lipidIdCol,
    annotationCol = annotationCol,
    sampleColumns = sampleColumns,
    isPattern = isPattern
  )

  if (!is.null(workflowData$state_manager)) {
    workflowData$state_manager$setWorkflowType(workflowType)
    logInfo(sprintf("Workflow type set to: %s", workflowType))
  }

  invisible(workflowType)
}

finalizeLipidImportSetupState <- function(
  workflowData,
  assayList,
  detectedFormat,
  lipidIdCol,
  sampleColumns,
  artifactResult = NULL,
  now = Sys.time,
  logInfo = logger::log_info
) {
  workflowData$processing_log$setup_import <- list(
    timestamp = now(),
    n_assays = length(assayList),
    assay_names = names(assayList),
    detected_format = detectedFormat,
    n_lipids = sapply(assayList, function(assayData) {
      if (lipidIdCol %in% names(assayData)) {
        length(unique(assayData[[lipidIdCol]]))
      } else {
        nrow(assayData)
      }
    }),
    n_samples = length(sampleColumns),
    artifacts = artifactResult
  )

  updated_status <- workflowData$tab_status
  updated_status$setup_import <- "complete"
  if (identical(updated_status$design_matrix, "disabled")) {
    updated_status$design_matrix <- "pending"
  }
  workflowData$tab_status <- updated_status

  logInfo(sprintf(
    "Lipidomics import complete: %d assays, %d total lipids",
    length(assayList),
    sum(sapply(assayList, nrow))
  ))

  invisible(workflowData$processing_log$setup_import)
}

sanitizeLipidImportSampleNames <- function(
  assayList,
  sampleColumns,
  sanitizeNames,
  makeCleanNames = janitor::make_clean_names,
  logInfo = logger::log_info,
  notify = shiny::showNotification
) {
  if (!isTRUE(sanitizeNames)) {
    return(list(
      assayList = assayList,
      sampleColumns = sampleColumns
    ))
  }

  logInfo("Sanitizing sample names in lipidomics data...")

  original_sample_cols <- sampleColumns
  cleaned_sample_cols <- makeCleanNames(original_sample_cols)

  sanitized_assays <- purrr::map(assayList, function(assayData) {
    all_cols <- names(assayData)
    match_indices <- match(original_sample_cols, all_cols)
    valid_matches <- !is.na(match_indices)

    if (any(valid_matches)) {
      names(assayData)[match_indices[valid_matches]] <- cleaned_sample_cols[valid_matches]
    }

    assayData
  })

  logInfo(sprintf("Sanitized %d sample column names.", length(cleaned_sample_cols)))
  notify("Sample names sanitized for R compatibility.", type = "message")

  list(
    assayList = sanitized_assays,
    sampleColumns = cleaned_sample_cols
  )
}

assembleLipidImportAssayList <- function(
  assay1Name,
  assay1Data,
  assay2File = NULL,
  assay2Name = NULL,
  dataFormat = "msdial",
  lipidIdCol = NULL,
  annotationCol = NULL,
  importSecondAssay = NULL,
  resolveSecondAssayReader = resolveLipidImportSecondAssayReader,
  callSecondAssayReader = callLipidImportSecondAssayReader
) {
  assay_list <- list()
  assay_list[[assay1Name]] <- assay1Data

  if (!is.null(assay2File) && nzchar(assay2Name)) {
    if (is.null(importSecondAssay)) {
      importSecondAssay <- resolveSecondAssayReader(dataFormat = dataFormat)
    }
    assay2_import <- callSecondAssayReader(
      importSecondAssay = importSecondAssay,
      assay2File = assay2File,
      lipidIdCol = lipidIdCol,
      annotationCol = annotationCol
    )
    assay_list[[assay2Name]] <- assay2_import$data
  }

  assay_list
}

resolveLipidImportSecondAssayReader <- function(
  dataFormat,
  importMsdial = importLipidMSDIALData,
  importLipidSearch = importLipidSearchData
) {
  switch(
    tolower(dataFormat %||% "msdial"),
    lipidsearch = importLipidSearch,
    msdial = importMsdial,
    importMsdial
  )
}

callLipidImportSecondAssayReader <- function(
  importSecondAssay,
  assay2File,
  lipidIdCol = NULL,
  annotationCol = NULL
) {
  importer_formals <- tryCatch(names(formals(importSecondAssay)), error = function(e) character(0))
  args <- list(assay2File)

  if ("lipid_id_column" %in% importer_formals) {
    args$lipid_id_column <- lipidIdCol
  }
  if ("annotation_column" %in% importer_formals) {
    args$annotation_column <- annotationCol
  }

  do.call(importSecondAssay, args)
}

lipidImportFunctionAcceptsArg <- function(fn, arg) {
  formal_names <- names(formals(fn))
  arg %in% formal_names || "..." %in% formal_names
}

lipidImportScalarString <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

lipidImportScalarValue <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return("")
  }
  as.character(value[[1]])
}

validateLipidImportAssayNames <- function(assay1Name, assay2Name = NULL, assay2File = NULL) {
  assay1Name <- lipidImportScalarValue(assay1Name)
  assay2Name <- lipidImportScalarValue(assay2Name)
  assay2File <- lipidImportScalarValue(assay2File)
  has_assay2_file <- lipidImportScalarString(assay2File)
  has_assay2_name <- lipidImportScalarString(assay2Name)

  if (!lipidImportScalarString(assay1Name)) {
    stop("Primary lipidomics assay name is required", call. = FALSE)
  }

  assay_names <- c(assay1Name, if (has_assay2_name) assay2Name else character())
  if (any(grepl("[/\\\\]", assay_names))) {
    stop("Assay names must not contain path separators", call. = FALSE)
  }

  if (has_assay2_file != has_assay2_name) {
    stop("Second lipidomics assay requires both a file and an assay name", call. = FALSE)
  }

  if (has_assay2_file && identical(assay1Name, assay2Name)) {
    stop("Duplicate lipidomics assay names are not allowed", call. = FALSE)
  }

  invisible(TRUE)
}

buildLipidImportColumnMapping <- function(lipidIdCol, annotationCol, sampleColumns, isPattern) {
  list(
    lipid_id_col = lipidIdCol,
    annotation_col = if (lipidImportScalarString(annotationCol)) annotationCol else NULL,
    sample_columns = sampleColumns,
    is_pattern = if (lipidImportScalarString(isPattern)) isPattern else NA_character_
  )
}

validateLipidImportProcessingState <- function(
  assayList,
  lipidIdCol,
  sampleColumns,
  validateMapping = validateLipidColumnMapping
) {
  errors <- character()
  warnings <- character()
  summaries <- list()

  if (!is.list(assayList) || length(assayList) == 0L) {
    errors <- c(errors, "No lipidomics assays were assembled")
  }
  if (is.null(names(assayList)) || any(!nzchar(names(assayList)))) {
    errors <- c(errors, "All lipidomics assays must have non-empty names")
  }
  if (anyDuplicated(names(assayList)) > 0L) {
    errors <- c(errors, "Duplicate lipidomics assay names are not allowed")
  }
  if (is.null(sampleColumns) || length(sampleColumns) == 0L) {
    errors <- c(errors, "No sample columns specified")
  }

  for (assay_name in names(assayList)) {
    assay_data <- assayList[[assay_name]]
    if (!is.data.frame(assay_data)) {
      errors <- c(errors, sprintf("Assay '%s' is not a data frame", assay_name))
      next
    }

    mapping_result <- validateMapping(
      data = assay_data,
      lipid_id_column = lipidIdCol,
      sample_columns = sampleColumns
    )
    if (!isTRUE(mapping_result$valid)) {
      errors <- c(errors, sprintf(
        "Assay '%s': %s",
        assay_name,
        paste(mapping_result$errors, collapse = "; ")
      ))
    }
    if (length(mapping_result$warnings) > 0L) {
      warnings <- c(warnings, sprintf(
        "Assay '%s': %s",
        assay_name,
        mapping_result$warnings
      ))
    }
    summaries[[assay_name]] <- mapping_result$summary

    present_samples <- intersect(sampleColumns, names(assay_data))
    non_numeric <- present_samples[!vapply(assay_data[present_samples], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      errors <- c(errors, sprintf(
        "Assay '%s' has non-numeric sample columns: %s",
        assay_name,
        paste(non_numeric, collapse = ", ")
      ))
    }
  }

  list(
    valid = length(errors) == 0L,
    errors = errors,
    warnings = warnings,
    summary = summaries
  )
}

writeLipidImportSourceArtifacts <- function(
  assayList,
  columnMapping,
  experimentPaths = NULL,
  sourceFiles = NULL,
  writeTableFn = utils::write.table,
  writeLinesFn = writeLines,
  writeJsonFn = jsonlite::write_json,
  fileCopyFn = file.copy,
  dirCreateFn = dir.create,
  dirExistsFn = dir.exists,
  fileExistsFn = file.exists
) {
  source_dir <- experimentPaths$source_dir
  if (is.null(source_dir) || !nzchar(source_dir) || !dirExistsFn(source_dir)) {
    return(list(written = FALSE, reason = "source_dir unavailable"))
  }

  assay_names <- names(assayList)
  artifact_paths <- list()

  for (assay_name in assay_names) {
    assay_path <- file.path(source_dir, paste0("data_cln_", assay_name, ".tab"))
    writeTableFn(assayList[[assay_name]], file = assay_path, sep = "\t", row.names = FALSE, quote = FALSE)
    artifact_paths[[paste0("data_cln_", assay_name)]] <- assay_path
  }

  manifest_path <- file.path(source_dir, "assay_manifest.txt")
  writeLinesFn(assay_names, manifest_path)
  artifact_paths$assay_manifest <- manifest_path

  column_mapping_path <- file.path(source_dir, "column_mapping.json")
  writeJsonFn(columnMapping, column_mapping_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  artifact_paths$column_mapping <- column_mapping_path

  source_copy_dir <- file.path(source_dir, "import_sources")
  dirCreateFn(source_copy_dir, recursive = TRUE, showWarnings = FALSE)
  copied_sources <- character()
  if (!is.null(sourceFiles)) {
    sourceFiles <- sourceFiles[!vapply(sourceFiles, is.null, logical(1))]
    for (source_name in names(sourceFiles)) {
      source_file <- sourceFiles[[source_name]]
      if (lipidImportScalarString(source_file) && fileExistsFn(source_file)) {
        destination <- file.path(source_copy_dir, paste0(source_name, "_", basename(source_file)))
        fileCopyFn(source_file, destination, overwrite = TRUE)
        copied_sources <- c(copied_sources, destination)
      }
    }
  }

  summary_rows <- data.frame(
    assay_name = assay_names,
    feature_count = vapply(assayList, nrow, integer(1)),
    sample_count = length(columnMapping$sample_columns),
    sample_columns = paste(columnMapping$sample_columns, collapse = ","),
    stringsAsFactors = FALSE
  )
  summary_path <- file.path(source_dir, "lipidomics_import_summary.tsv")
  writeTableFn(summary_rows, file = summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
  artifact_paths$import_summary <- summary_path

  list(
    written = TRUE,
    paths = artifact_paths,
    source_copies = copied_sources
  )
}

runLipidImportProcessing <- function(
  workflowData,
  assay1Name,
  assay1Data,
  assay1File = NULL,
  assay2File = NULL,
  assay2Name = NULL,
  vendorFormat,
  detectedFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  sanitizeNames,
  experimentPaths = NULL,
  assembleAssayList = assembleLipidImportAssayList,
  sanitizeSampleNames = sanitizeLipidImportSampleNames,
  validateAssayNames = validateLipidImportAssayNames,
  validateImportState = validateLipidImportProcessingState,
  writeImportArtifacts = writeLipidImportSourceArtifacts,
  applyResultToWorkflow = applyLipidImportResultToWorkflow,
  finalizeSetupState = finalizeLipidImportSetupState,
  notify = shiny::showNotification,
  removeNotify = shiny::removeNotification,
  logError = logger::log_error
) {
  notify(
    "Processing imported data...",
    id = "lipid_import_working",
    duration = NULL
  )

  tryCatch(
    {
      validateAssayNames(
        assay1Name = assay1Name,
        assay2Name = assay2Name,
        assay2File = assay2File
      )

      assay_list <- assembleAssayList(
        assay1Name = assay1Name,
        assay1Data = assay1Data,
        assay2File = assay2File,
        assay2Name = assay2Name,
        dataFormat = if (identical(vendorFormat, "custom")) detectedFormat else detectedFormat %||% vendorFormat,
        lipidIdCol = lipidIdCol,
        annotationCol = annotationCol
      )

      sanitized_import <- sanitizeSampleNames(
        assayList = assay_list,
        sampleColumns = sampleColumns,
        sanitizeNames = sanitizeNames
      )
      assay_list <- sanitized_import$assayList
      sampleColumns <- sanitized_import$sampleColumns

      validation_result <- validateImportState(
        assayList = assay_list,
        lipidIdCol = lipidIdCol,
        sampleColumns = sampleColumns
      )
      if (!isTRUE(validation_result$valid)) {
        stop(sprintf(
          "Invalid lipidomics import: %s",
          paste(validation_result$errors, collapse = "; ")
        ), call. = FALSE)
      }

      column_mapping <- buildLipidImportColumnMapping(
        lipidIdCol = lipidIdCol,
        annotationCol = annotationCol,
        sampleColumns = sampleColumns,
        isPattern = isPattern
      )
      artifact_result <- writeImportArtifacts(
        assayList = assay_list,
        columnMapping = column_mapping,
        experimentPaths = experimentPaths,
        sourceFiles = c(assay1 = assay1File, assay2 = assay2File)
      )

      applyResultToWorkflow(
        workflowData = workflowData,
        assayList = assay_list,
        dataFormat = if (identical(vendorFormat, "custom")) "custom" else detectedFormat,
        lipidIdCol = lipidIdCol,
        annotationCol = annotationCol,
        sampleColumns = sampleColumns,
        isPattern = isPattern
      )

      finalizeSetupState(
        workflowData = workflowData,
        assayList = assay_list,
        detectedFormat = workflowData$data_format,
        lipidIdCol = lipidIdCol,
        sampleColumns = sampleColumns,
        artifactResult = artifact_result
      )

      removeNotify("lipid_import_working")
      notify("Data imported successfully!", type = "message")

      invisible(list(
        status = "success",
        assayList = assay_list,
        sampleColumns = sampleColumns,
        validationResult = validation_result,
        artifactResult = artifact_result
      ))
    },
    error = function(e) {
      logError(paste("Error processing import:", e$message))
      removeNotify("lipid_import_working")
      notify(paste("Error:", e$message), type = "error", duration = 10)
      invisible(NULL)
    }
  )
}

loadLipidImportAssayPreview <- function(
  assay1File,
  readHeaders = function(path) {
    tryCatch(
      {
        con <- file(path, "r")
        on.exit(close(con), add = TRUE)

        first_line <- readLines(con, n = 1)
        raw_headers <- if (grepl(",", first_line)) {
          strsplit(first_line, ",")[[1]]
        } else {
          strsplit(first_line, "\t")[[1]]
        }

        gsub('^"|"$', "", raw_headers)
      },
      error = function(e) {
        character(0)
      }
    )
  },
  detectFormat = detectLipidomicsFormat,
  importMsdial = importLipidMSDIALData,
  importLipidSearch = importLipidSearchData,
  logInfo = logger::log_info
) {
  headers <- readHeaders(assay1File)

  if (length(headers) == 0) {
    stop("Could not read headers from file")
  }

  format_info <- detectFormat(
    headers = headers,
    filename = basename(assay1File)
  )

  import_result <- switch(format_info$format,
    "msdial" = importMsdial(assay1File),
    "lipidsearch" = importLipidSearch(assay1File),
    importMsdial(assay1File)
  )

  logInfo(sprintf(
    "Imported assay 1: %d rows, %d columns, format: %s",
    nrow(import_result$data),
    ncol(import_result$data),
    format_info$format
  ))

  list(
    headers = headers,
    detectedFormat = format_info$format,
    formatConfidence = format_info$confidence,
    importResult = import_result,
    assayData = import_result$data,
    updates = list(
      lipidId = list(
        choices = headers,
        selected = import_result$detected_columns$lipid_id
      ),
      annotation = list(
        choices = c("(None)" = "", headers),
        selected = import_result$detected_columns$annotation
      ),
      isPattern = if (!is.null(import_result$is_pattern) && !is.na(import_result$is_pattern)) {
        import_result$is_pattern
      } else {
        NULL
      }
    )
  )
}

handleLipidImportDataPreviewLoad <- function(
  session,
  localData,
  loadPreview = loadLipidImportAssayPreview,
  applyPreview = applyLipidImportPreviewToModuleState,
  handleImportError = handleLipidImportDataImportError
) {
  shiny::req(localData$assay1_file)

  tryCatch(
    {
      importPreview <- loadPreview(
        assay1File = localData$assay1_file
      )

      applyPreview(
        session = session,
        localData = localData,
        importPreview = importPreview
      )
    },
    error = function(e) {
      handleImportError(error = e)
    }
  )
}

buildLipidImportAssay1SelectedCallback <- function(
  session,
  localData,
  runPreviewLoad = handleLipidImportDataPreviewLoad
) {
  force(session)
  force(localData)
  force(runPreviewLoad)

  function() {
    runPreviewLoad(
      session = session,
      localData = localData
    )
  }
}

applyLipidImportPreviewToModuleState <- function(
  session,
  localData,
  importPreview,
  updateSelectInput = shiny::updateSelectInput,
  updateTextInput = shiny::updateTextInput
) {
  localData$all_headers <- importPreview$headers
  localData$detected_format <- importPreview$detectedFormat
  localData$format_confidence <- importPreview$formatConfidence
  localData$assay1_import_result <- importPreview$importResult
  localData$assay1_data <- importPreview$assayData

  updateSelectInput(
    session,
    "lipid_id_col",
    choices = importPreview$updates$lipidId$choices,
    selected = importPreview$updates$lipidId$selected
  )

  updateSelectInput(
    session,
    "annotation_col",
    choices = importPreview$updates$annotation$choices,
    selected = importPreview$updates$annotation$selected
  )

  if (!is.null(importPreview$updates$isPattern)) {
    updateTextInput(
      session,
      "is_pattern",
      value = importPreview$updates$isPattern
    )
  }

  invisible(importPreview)
}

handleLipidImportDataImportError <- function(
  error,
  logError = logger::log_error,
  notify = shiny::showNotification
) {
  errorMessage <- paste("Error importing data:", error$message)

  logError(errorMessage)
  notify(errorMessage, type = "error")

  invisible(FALSE)
}

handleLipidImportFileSelection <- function(
  fileInput,
  volumes,
  onPathSelected,
  parseFilePaths = shinyFiles::parseFilePaths,
  logError = logger::log_error
) {
  if (is.null(fileInput) || is.integer(fileInput)) {
    return(FALSE)
  }

  tryCatch(
    {
      file_info <- parseFilePaths(volumes, fileInput)
      if (nrow(file_info) == 0) {
        return(FALSE)
      }

      onPathSelected(as.character(file_info$datapath[1]))
      TRUE
    },
    error = function(e) {
      logError(paste("Error parsing file path:", e$message))
      FALSE
    }
  )
}

registerLipidImportAssayFileSelectionObserver <- function(
  input,
  fileInputId,
  volumes,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  handleSelectionEvent = handleLipidImportAssayFileSelectionEvent,
  observeEvent = shiny::observeEvent
) {
  observeEvent(input[[fileInputId]], {
    handleSelectionEvent(
      fileInput = input[[fileInputId]],
      volumes = volumes,
      assignSelectedPath = assignSelectedPath,
      setRenderedPath = setRenderedPath,
      onSelected = onSelected
    )
  })
}

registerLipidImportShinyFileInputs <- function(
  input,
  output,
  session,
  localData,
  volumes,
  onAssay1Selected = NULL,
  setupAssayFileChooser = setupLipidImportAssayFileChooser,
  registerSelectionObserver = registerLipidImportAssayFileSelectionObserver
) {
  setupAssayFileChooser(
    input = input,
    fileInputId = "assay1_file",
    volumes = volumes,
    session = session
  )

  registerSelectionObserver(
    input = input,
    fileInputId = "assay1_file",
    volumes = volumes,
    assignSelectedPath = function(value) {
      localData$assay1_file <- value
    },
    setRenderedPath = function(renderedPath) {
      output$assay1_path <- renderedPath
    },
    onSelected = onAssay1Selected
  )

  setupAssayFileChooser(
    input = input,
    fileInputId = "assay2_file",
    volumes = volumes,
    session = session
  )

  registerSelectionObserver(
    input = input,
    fileInputId = "assay2_file",
    volumes = volumes,
    assignSelectedPath = function(value) {
      localData$assay2_file <- value
    },
    setRenderedPath = function(renderedPath) {
      output$assay2_path <- renderedPath
    }
  )

  invisible(volumes)
}

setupLipidImportShinyFileInputs <- function(
  useShinyFiles,
  input,
  output,
  session,
  localData,
  volumes,
  prepareVolumes = prepareLipidImportShinyFileVolumes,
  buildAssay1Selected = buildLipidImportAssay1SelectedCallback,
  registerInputs = registerLipidImportShinyFileInputs
) {
  if (!isTRUE(useShinyFiles)) {
    return(invisible(volumes))
  }

  preparedVolumes <- prepareVolumes(volumes = volumes)

  registerInputs(
    input = input,
    output = output,
    session = session,
    localData = localData,
    volumes = preparedVolumes,
    onAssay1Selected = buildAssay1Selected(
      session = session,
      localData = localData
    )
  )

  invisible(preparedVolumes)
}

handleLipidImportStandardFileSelection <- function(
  fileInput,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  assignAssayPath = handleLipidImportSelectedAssayPathAssignment
) {
  if (is.null(fileInput) || is.null(fileInput$datapath) || length(fileInput$datapath) == 0) {
    return(invisible(NULL))
  }

  selectedPath <- as.character(fileInput$datapath[1])
  if (!nzchar(selectedPath)) {
    return(invisible(NULL))
  }

  assignAssayPath(
    selectedPath = selectedPath,
    assignSelectedPath = assignSelectedPath,
    setRenderedPath = setRenderedPath,
    onSelected = onSelected
  )
}

registerLipidImportStandardFileInputs <- function(
  input,
  output,
  localData,
  onAssay1Selected = NULL,
  observeEvent = shiny::observeEvent,
  handleStandardSelection = handleLipidImportStandardFileSelection
) {
  observeEvent(input$assay1_file_std, {
    handleStandardSelection(
      fileInput = input$assay1_file_std,
      assignSelectedPath = function(value) {
        localData$assay1_file <- value
      },
      setRenderedPath = function(renderedPath) {
        output$assay1_path <- renderedPath
      },
      onSelected = onAssay1Selected
    )
  })

  observeEvent(input$assay2_file_std, {
    handleStandardSelection(
      fileInput = input$assay2_file_std,
      assignSelectedPath = function(value) {
        localData$assay2_file <- value
      },
      setRenderedPath = function(renderedPath) {
        output$assay2_path <- renderedPath
      }
    )
  })

  invisible(NULL)
}

setupLipidImportStandardFileInputs <- function(
  useShinyFiles,
  input,
  output,
  session,
  localData,
  buildAssay1Selected = buildLipidImportAssay1SelectedCallback,
  registerInputs = registerLipidImportStandardFileInputs
) {
  if (isTRUE(useShinyFiles)) {
    return(invisible(NULL))
  }

  registerInputs(
    input = input,
    output = output,
    localData = localData,
    onAssay1Selected = buildAssay1Selected(
      session = session,
      localData = localData
    )
  )

  invisible(NULL)
}

registerLipidImportProcessObserver <- function(
  input,
  workflowData,
  localData,
  getLipidIdCol,
  getAnnotationCol,
  getSampleColumns,
  experimentPaths = NULL,
  handleProcessRequest = handleLipidImportProcessRequest,
  observeEvent = shiny::observeEvent
) {
  observeEvent(input$process_import, {
    handleProcessRequest(
      workflowData = workflowData,
      assay1Name = input$assay1_name,
      assay1Data = localData$assay1_data,
      assay1File = localData$assay1_file,
      assay2File = localData$assay2_file,
      assay2Name = input$assay2_name,
      vendorFormat = input$vendor_format,
      detectedFormat = localData$detected_format,
      lipidIdCol = getLipidIdCol(),
      annotationCol = getAnnotationCol(),
      sampleColumns = getSampleColumns(),
      isPattern = input$is_pattern,
      sanitizeNames = input$sanitize_names,
      experimentPaths = experimentPaths
    )
  })
}

handleLipidImportProcessRequest <- function(
  workflowData,
  assay1Name,
  assay1Data,
  assay1File = NULL,
  assay2File = NULL,
  assay2Name = NULL,
  vendorFormat,
  detectedFormat,
  lipidIdCol,
  annotationCol,
  sampleColumns,
  isPattern,
  sanitizeNames,
  experimentPaths = NULL,
  processImport = runLipidImportProcessing
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  process_args <- list(
    workflowData = workflowData,
    assay1Name = assay1Name,
    assay1Data = assay1Data,
    assay2File = assay2File,
    assay2Name = assay2Name,
    vendorFormat = vendorFormat,
    detectedFormat = detectedFormat,
    lipidIdCol = lipidIdCol,
    annotationCol = annotationCol,
    sampleColumns = sampleColumns,
    isPattern = isPattern,
    sanitizeNames = sanitizeNames
  )
  if (lipidImportFunctionAcceptsArg(processImport, "assay1File")) {
    process_args$assay1File <- assay1File
  }
  if (lipidImportFunctionAcceptsArg(processImport, "experimentPaths")) {
    process_args$experimentPaths <- experimentPaths
  }

  do.call(processImport, process_args)
}

handleLipidImportValidationSummaryRender <- function(
  assay1Data,
  lipidIdCol,
  sampleColumns,
  validateMapping = validateLipidColumnMapping,
  buildSummary = buildLipidImportValidationSummary
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  validation <- validateMapping(
    data = assay1Data,
    lipid_id_column = lipidIdCol,
    sample_columns = sampleColumns
  )

  buildSummary(validation)
}

handleLipidImportSampleColumnsDisplayRender <- function(
  assay1ImportResult,
  formatPreviewText = formatLipidImportColumnPreviewText
) {
  shiny::req(assay1ImportResult)

  formatPreviewText(
    columnNames = assay1ImportResult$sample_columns,
    maxColumns = 10
  )
}

handleLipidImportAssayPathRender <- function(selectedPath) {
  selectedPath
}

handleLipidImportSelectedAssayPathAssignment <- function(
  selectedPath,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  renderText = shiny::renderText,
  renderAssayPath = handleLipidImportAssayPathRender
) {
  assignSelectedPath(selectedPath)
  setRenderedPath(renderText({
    renderAssayPath(selectedPath)
  }))

  if (is.function(onSelected)) {
    onSelected()
  }

  invisible(selectedPath)
}

handleLipidImportAssayFileSelectionEvent <- function(
  fileInput,
  volumes,
  assignSelectedPath,
  setRenderedPath,
  onSelected = NULL,
  handleFileSelection = handleLipidImportFileSelection,
  assignAssayPath = handleLipidImportSelectedAssayPathAssignment
) {
  handleFileSelection(
    fileInput = fileInput,
    volumes = volumes,
    onPathSelected = function(selectedPath) {
      assignAssayPath(
        selectedPath = selectedPath,
        assignSelectedPath = assignSelectedPath,
        setRenderedPath = setRenderedPath,
        onSelected = onSelected
      )
    }
  )
}

setupLipidImportAssayFileChooser <- function(
  input,
  fileInputId,
  volumes,
  session,
  filetypes = c("tsv", "tab", "txt", "csv", "xlsx", "parquet"),
  shinyFileChoose = shinyFiles::shinyFileChoose
) {
  shinyFileChoose(
    input = input,
    id = fileInputId,
    roots = volumes,
    session = session,
    filetypes = filetypes
  )
}

prepareLipidImportShinyFileVolumes <- function(
  volumes,
  getVolumes = shinyFiles::getVolumes,
  emitMessage = message
) {
  if (is.null(volumes)) {
    emitMessage("   mod_lipid_import_server: volumes is NULL, creating from getVolumes()")
    volumes <- getVolumes()()
  }

  emitMessage(sprintf(
    "   mod_lipid_import_server: volumes type = %s, length = %d",
    typeof(volumes), length(volumes)
  ))

  if (length(volumes) > 0) {
    emitMessage(sprintf(
      "   mod_lipid_import_server: volumes names = %s",
      paste(names(volumes), collapse = ", ")
    ))
  } else {
    emitMessage("   mod_lipid_import_server: WARNING - volumes is empty!")
  }

  volumes
}

registerLipidImportFileLoadedOutput <- function(
  output,
  localData,
  reactiveFn = shiny::reactive,
  outputOptionsFn = shiny::outputOptions
) {
  output$file_loaded <- reactiveFn({
    !is.null(localData$assay1_data)
  })
  outputOptionsFn(output, "file_loaded", suspendWhenHidden = FALSE)

  invisible(output)
}

registerLipidImportValidationSummaryOutput <- function(
  output,
  localData,
  getLipidIdCol,
  getSampleColumns,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportValidationSummaryRender
) {
  output$validation_summary <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdCol = getLipidIdCol(),
      sampleColumns = getSampleColumns()
    )
  })

  invisible(output)
}

registerLipidImportFormatDetectionStatusOutput <- function(
  output,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportFormatDetectionStatusRender
) {
  output$format_detection_status <- renderUi({
    handleRender(
      detectedFormat = localData$detected_format,
      formatConfidence = localData$format_confidence
    )
  })

  invisible(output)
}

registerLipidImportLipidIdStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportLipidIdStatusRender
) {
  output$lipid_id_status <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdCol = input$lipid_id_col
    )
  })

  invisible(output)
}

registerLipidImportAnnotationStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportAnnotationStatusRender
) {
  output$annotation_status <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      annotationCol = input$annotation_col
    )
  })

  invisible(output)
}

registerLipidImportSampleColumnsDisplayOutput <- function(
  output,
  localData,
  renderText = shiny::renderText,
  handleRender = handleLipidImportSampleColumnsDisplayRender
) {
  output$sample_columns_display <- renderText({
    handleRender(
      assay1ImportResult = localData$assay1_import_result
    )
  })

  invisible(output)
}

registerLipidImportAvailableColumnsDisplayOutput <- function(
  output,
  localData,
  renderText = shiny::renderText,
  handleRender = handleLipidImportAvailableColumnsDisplayRender
) {
  output$available_columns_display <- renderText({
    handleRender(
      allHeaders = localData$all_headers
    )
  })

  invisible(output)
}

registerLipidImportCustomLipidIdStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportCustomLipidIdStatusRender
) {
  output$lipid_id_status_custom <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      lipidIdColCustom = input$lipid_id_col_custom
    )
  })

  invisible(output)
}

registerLipidImportCustomAnnotationStatusOutput <- function(
  output,
  input,
  localData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportCustomAnnotationStatusRender
) {
  output$annotation_status_custom <- renderUi({
    handleRender(
      assay1Data = localData$assay1_data,
      annotationColCustom = input$annotation_col_custom
    )
  })

  invisible(output)
}

registerLipidImportStatusOutput <- function(
  output,
  workflowData,
  renderUi = shiny::renderUI,
  handleRender = handleLipidImportStatusRender
) {
  output$import_status <- renderUi({
    handleRender(workflowData = workflowData)
  })

  invisible(output)
}

registerLipidImportModuleOutputs <- function(
  output,
  input,
  localData,
  workflowData,
  getLipidIdCol,
  getSampleColumns,
  registerFileLoadedOutput = registerLipidImportFileLoadedOutput,
  registerFormatDetectionStatusOutput = registerLipidImportFormatDetectionStatusOutput,
  registerLipidIdStatusOutput = registerLipidImportLipidIdStatusOutput,
  registerAnnotationStatusOutput = registerLipidImportAnnotationStatusOutput,
  registerSampleColumnsDisplayOutput = registerLipidImportSampleColumnsDisplayOutput,
  registerAvailableColumnsDisplayOutput = registerLipidImportAvailableColumnsDisplayOutput,
  registerCustomLipidIdStatusOutput = registerLipidImportCustomLipidIdStatusOutput,
  registerCustomAnnotationStatusOutput = registerLipidImportCustomAnnotationStatusOutput,
  registerValidationSummaryOutput = registerLipidImportValidationSummaryOutput,
  registerStatusOutput = registerLipidImportStatusOutput
) {
  registerFileLoadedOutput(
    output = output,
    localData = localData
  )
  registerFormatDetectionStatusOutput(
    output = output,
    localData = localData
  )
  registerLipidIdStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerAnnotationStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerSampleColumnsDisplayOutput(
    output = output,
    localData = localData
  )
  registerAvailableColumnsDisplayOutput(
    output = output,
    localData = localData
  )
  registerCustomLipidIdStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerCustomAnnotationStatusOutput(
    output = output,
    input = input,
    localData = localData
  )
  registerValidationSummaryOutput(
    output = output,
    localData = localData,
    getLipidIdCol = getLipidIdCol,
    getSampleColumns = getSampleColumns
  )
  registerStatusOutput(
    output = output,
    workflowData = workflowData
  )

  invisible(output)
}

handleLipidImportFormatDetectionStatusRender <- function(
  detectedFormat,
  formatConfidence,
  buildStatus = buildLipidImportFormatDetectionStatus
) {
  shiny::req(detectedFormat)

  buildStatus(
    detectedFormat = detectedFormat,
    formatConfidence = formatConfidence
  )
}

handleLipidImportAvailableColumnsDisplayRender <- function(
  allHeaders,
  formatPreviewText = formatLipidImportColumnPreviewText
) {
  shiny::req(allHeaders)

  formatPreviewText(columnNames = allHeaders)
}

handleLipidImportLipidIdStatusRender <- function(
  assay1Data,
  lipidIdCol,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)
  shiny::req(lipidIdCol)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = lipidIdCol,
    successMode = "unique_ids"
  )
}

handleLipidImportAnnotationStatusRender <- function(
  assay1Data,
  annotationCol,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = annotationCol,
    successMode = "ok",
    emptyMode = "optional"
  )
}

handleLipidImportCustomLipidIdStatusRender <- function(
  assay1Data,
  lipidIdColCustom,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = lipidIdColCustom,
    successMode = "found_unique_ids",
    emptyMode = "prompt",
    allowCaseInsensitive = TRUE
  )
}

handleLipidImportCustomAnnotationStatusRender <- function(
  assay1Data,
  annotationColCustom,
  buildValidationStatus = buildLipidImportColumnValidationStatus
) {
  shiny::req(assay1Data)

  buildValidationStatus(
    assayData = assay1Data,
    columnName = annotationColCustom,
    successMode = "found",
    emptyMode = "optional",
    allowCaseInsensitive = TRUE
  )
}

resolveLipidImportLipidIdColumn <- function(
  assay1Data,
  vendorFormat,
  lipidIdCol,
  lipidIdColCustom,
  resolveEffectiveColumn = resolveLipidImportEffectiveColumn
) {
  resolveEffectiveColumn(
    assayData = assay1Data,
    vendorFormat = vendorFormat,
    selectedColumn = lipidIdCol,
    customColumn = lipidIdColCustom
  )
}

resolveLipidImportAnnotationColumn <- function(
  assay1Data,
  vendorFormat,
  annotationCol,
  annotationColCustom,
  resolveEffectiveColumn = resolveLipidImportEffectiveColumn
) {
  resolveEffectiveColumn(
    assayData = assay1Data,
    vendorFormat = vendorFormat,
    selectedColumn = annotationCol,
    customColumn = annotationColCustom
  )
}

resolveLipidImportSelectedSampleColumns <- function(
  assay1Data,
  assay1ImportResult,
  vendorFormat,
  sampleColsPattern,
  excludeNormalized,
  resolveSampleColumns = resolveLipidImportSampleColumns
) {
  shiny::req(assay1Data)

  resolveSampleColumns(
    assayData = assay1Data,
    assayImportResult = assay1ImportResult,
    vendorFormat = vendorFormat,
    sampleColsPattern = sampleColsPattern,
    excludeNormalized = excludeNormalized
  )
}

buildLipidImportColumnSelectionReactives <- function(
  input,
  localData,
  reactiveFn = shiny::reactive,
  resolveLipidIdColumn = resolveLipidImportLipidIdColumn,
  resolveAnnotationColumn = resolveLipidImportAnnotationColumn,
  resolveSampleColumns = resolveLipidImportSelectedSampleColumns
) {
  list(
    lipidIdCol = reactiveFn({
      resolveLipidIdColumn(
        assay1Data = localData$assay1_data,
        vendorFormat = input$vendor_format,
        lipidIdCol = input$lipid_id_col,
        lipidIdColCustom = input$lipid_id_col_custom
      )
    }),
    annotationCol = reactiveFn({
      resolveAnnotationColumn(
        assay1Data = localData$assay1_data,
        vendorFormat = input$vendor_format,
        annotationCol = input$annotation_col,
        annotationColCustom = input$annotation_col_custom
      )
    }),
    sampleColumns = reactiveFn({
      resolveSampleColumns(
        assay1Data = localData$assay1_data,
        assay1ImportResult = localData$assay1_import_result,
        vendorFormat = input$vendor_format,
        sampleColsPattern = input$sample_cols_pattern,
        excludeNormalized = input$exclude_norm
      )
    })
  )
}

handleLipidImportStatusRender <- function(
  workflowData,
  buildStatusDisplay = buildLipidImportStatusDisplay
) {
  buildStatusDisplay(workflowData)
}
