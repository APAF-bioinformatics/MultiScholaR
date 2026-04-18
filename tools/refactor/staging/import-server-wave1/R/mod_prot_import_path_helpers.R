updateProtImportCheckpointCaptureOption <- function(
    captureCheckpoints
    , optionsFn = options
    , logInfo = logger::log_info
) {
  capture_enabled <- isTRUE(captureCheckpoints)
  optionsFn(multischolar.capture_test_checkpoints = capture_enabled)

  if (capture_enabled) {
    logInfo("Proteomics test checkpoint capture ENABLED")
  } else {
    logInfo("Proteomics test checkpoint capture DISABLED")
  }

  invisible(capture_enabled)
}

toggleProtImportMixedSpeciesInputs <- function(
    mixedSpeciesFasta
    , disableInput = shinyjs::disable
    , enableInput = shinyjs::enable
    , messageFn = message
) {
  if (isTRUE(mixedSpeciesFasta)) {
    disableInput("taxon_id")
    disableInput("organism_name")
    messageFn("[mod_prot_import] Mixed species FASTA enabled - organism inputs disabled")
    return(invisible(TRUE))
  }

  enableInput("taxon_id")
  enableInput("organism_name")
  messageFn("[mod_prot_import] Mixed species FASTA disabled - organism inputs enabled")

  invisible(FALSE)
}

resolveProtImportPrimaryUploadPath <- function(
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

resolveProtImportShinyFileVolumes <- function(
    volumes = NULL
    , getVolumes = shinyFiles::getVolumes
) {
  if (!is.null(volumes)) {
    return(volumes)
  }

  getVolumes()()
}

registerProtImportShinyFileChooser <- function(
    input
    , inputId
    , volumes
    , session
    , filetypes
    , chooseFile = shinyFiles::shinyFileChoose
) {
  chooseFile(
    input
    , inputId
    , roots = volumes
    , session = session
    , filetypes = filetypes
  )

  invisible(inputId)
}

handleProtImportShinyFileSelection <- function(
    selectedInput
    , volumes
    , localData
    , localField
    , output
    , outputId
    , parseFilePaths = shinyFiles::parseFilePaths
    , renderText = shiny::renderText
    , messageFn = message
    , catchErrors = FALSE
) {
  if (is.null(selectedInput) || is.integer(selectedInput)) {
    return(invisible(NULL))
  }

  bindSelection <- function() {
    file_info <- parseFilePaths(volumes, selectedInput)
    if (nrow(file_info) == 0) {
      return(invisible(NULL))
    }

    selected_path <- as.character(file_info$datapath[1])
    localData[[localField]] <- selected_path
    output[[outputId]] <- renderText(selected_path)

    invisible(selected_path)
  }

  if (!isTRUE(catchErrors)) {
    return(bindSelection())
  }

  tryCatch(
    bindSelection(),
    error = function(e) {
      messageFn(sprintf("   mod_prot_import_server ERROR in parseFilePaths: %s", e$message))
      invisible(NULL)
    }
  )
}

