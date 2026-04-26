# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

skipIfMissingMultiScholaRBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    testthat::skip(sprintf("Target-only extracted helper(s) not present: %s", paste(missing, collapse = ", ")))
  }
}

skipIfMissingMultiScholaRBindings(
  "buildProtImportStatusUi",
  "readProtImportHeaders",
  "mod_prot_import_server"
)

makeFunctionWithOverrides <- function(fun, replacements) {
  funOverride <- fun
  environment(funOverride) <- list2env(replacements, parent = environment(fun))
  funOverride
}

withCleanGlobalObjects <- function(objectNames, code) {
  hadExisting <- vapply(
    objectNames,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
  oldValues <- lapply(seq_along(objectNames), function(i) {
    if (hadExisting[[i]]) {
      get(objectNames[[i]], envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(oldValues) <- objectNames

  on.exit({
    for (name in rev(objectNames)) {
      if (hadExisting[[name]]) {
        assign(name, oldValues[[name]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

makeFileInput <- function(path, type = "text/plain") {
  data.frame(
    name = basename(path),
    size = unname(file.info(path)$size),
    type = type,
    datapath = path,
    stringsAsFactors = FALSE
  )
}

makeImportFixture <- function(withConfig = TRUE) {
  rootDir <- tempfile("prot-import-module-")
  dir.create(rootDir, recursive = TRUE)

  searchResultsPath <- file.path(rootDir, "report.tsv")
  writeLines(
    c(
      paste(
        c(
          "Protein.Group", "Protein.Ids", "Protein.Names", "Precursor.Id",
          "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge",
          "Q.Value", "PG.Q.Value", "Run"
        ),
        collapse = "\t"
      ),
      paste(
        c(
          "P1", "P1", "Protein 1", "PEPTIDE1", "PEPTIDE1", "PEPTIDE1",
          "2", "0.01", "0.01", "Sample 1"
        ),
        collapse = "\t"
      )
    ),
    searchResultsPath
  )

  fastaPath <- file.path(rootDir, "proteins.fasta")
  writeLines(
    c(
      ">sp|P1|PROT1 OS=Homo sapiens OX=9606",
      "MPEPTIDESEQ"
    ),
    fastaPath
  )

  configPath <- NULL
  if (withConfig) {
    configPath <- file.path(rootDir, "config.ini")
    writeLines("[generalParameters]\nmin_peptides_per_protein=2", configPath)
  }

  resultsDir <- file.path(rootDir, "results")
  sourceDir <- file.path(rootDir, "scripts")
  dir.create(resultsDir, recursive = TRUE)
  dir.create(sourceDir, recursive = TRUE)

  list(
    rootDir = rootDir,
    searchResultsPath = searchResultsPath,
    fastaPath = fastaPath,
    configPath = configPath,
    experimentPaths = list(
      results_dir = resultsDir,
      source_dir = sourceDir
    )
  )
}

makeWorkflowData <- function() {
  shiny::reactiveValues(
    data_tbl = NULL,
    fasta_file_path = NULL,
    config_list = NULL,
    taxon_id = NULL,
    organism_name = NULL,
    design_matrix = NULL,
    data_cln = NULL,
    contrasts_tbl = NULL,
    state_manager = WorkflowState$new(),
    peptide_data = NULL,
    protein_log2_quant = NULL,
    protein_data = NULL,
    ruv_normalised_for_da_analysis_obj = NULL,
    da_analysis_results_list = NULL,
    uniprot_dat_cln = NULL,
    enrichment_results = NULL,
    data_format = NULL,
    data_type = NULL,
    column_mapping = NULL,
    aa_seq_tbl_final = NULL,
    fasta_metadata = NULL,
    uniprot_mapping = NULL,
    uniparc_mapping = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(
      setup_import = "pending",
      design_matrix = "disabled",
      quality_control = "disabled",
      normalization = "disabled",
      differential_expression = "disabled",
      enrichment_analysis = "disabled",
      session_summary = "disabled"
    ),
    state_update_trigger = NULL,
    processing_log = list()
  )
}

mockImportResult <- function() {
  list(
    data = data.frame(
      Protein.Group = c("P1", "P2"),
      Protein.Ids = c("P1", "P2"),
      Precursor.Id = c("PEPTIDE1", "PEPTIDE2"),
      Run = c("Sample 1", "Sample 2"),
      Precursor.Quantity = c(1000, 2000),
      stringsAsFactors = FALSE
    ),
    data_type = "peptide",
    column_mapping = list(
      protein_col = "Protein.Group",
      peptide_col = "Precursor.Id",
      run_col = "Run",
      quantity_col = "Precursor.Quantity"
    )
  )
}

mockFastaResult <- function() {
  list(
    aa_seq_tbl_final = data.frame(
      uniprot_acc = c("P1", "P2"),
      sequence = c("MPEPTIDESEQ", "MPEPTIDESEQ2"),
      stringsAsFactors = FALSE
    ),
    fasta_metadata = list(
      fasta_format = "standard",
      num_sequences = 2
    )
  )
}

test_that("buildProtImportStatusUi reports processing and success states", {
  processingHtml <- htmltools::renderTags(
    buildProtImportStatusUi(
      isProcessing = TRUE,
      dataTbl = NULL,
      dataFormat = NULL,
      dataType = NULL
    )
  )$html

  expect_match(processingHtml, "Processing\\.\\.\\.")
  expect_match(processingHtml, "Please wait while data is being imported and validated\\.")

  successHtml <- htmltools::renderTags(
    buildProtImportStatusUi(
      isProcessing = FALSE,
      dataTbl = data.frame(x = 1),
      dataFormat = "diann",
      dataType = "peptide"
    )
  )$html

  expect_match(successHtml, "\\[OK\\] Data imported successfully!")
  expect_match(successHtml, "Format: diann \\| Type: peptide")
  expect_null(buildProtImportStatusUi(FALSE, NULL, NULL, NULL))
})

test_that("buildProtImportDataSummaryUi renders summary fields and peptide row conditionally", {
  summaryHtml <- htmltools::renderTags(
    buildProtImportDataSummaryUi(
      dataTbl = data.frame(x = 1),
      setupImportLog = list(
        detected_format = "diann",
        data_type = "peptide",
        n_rows = 1234,
        n_runs = 2,
        n_proteins = 345,
        n_peptides = 678,
        organism = "Homo sapiens",
        taxon_id = 9606
      )
    )
  )$html

  expect_match(summaryHtml, "Data Summary")
  expect_match(summaryHtml, "Data format: DIANN")
  expect_match(summaryHtml, "Total rows: 1,234")
  expect_match(summaryHtml, "Number of peptides: 678")
  expect_match(summaryHtml, "Organism: Homo sapiens \\(taxon: 9606 \\)")

  noPeptideHtml <- htmltools::renderTags(
    buildProtImportDataSummaryUi(
      dataTbl = data.frame(x = 1),
      setupImportLog = list(
        detected_format = "fragpipe",
        data_type = "protein",
        n_rows = 20,
        n_runs = 3,
        n_proteins = 10,
        n_peptides = NA,
        organism = "Mus musculus",
        taxon_id = 10090
      )
    )
  )$html

  expect_no_match(noPeptideHtml, "Number of peptides:")
})

test_that("buildProtImportFormatDetectionUi maps format labels and confidence states", {
  detectedHtml <- htmltools::renderTags(
    buildProtImportFormatDetectionUi(
      detectedFormat = "diann",
      formatConfidence = 0.91
    )
  )$html

  expect_match(detectedHtml, "alert-success")
  expect_match(detectedHtml, "Detected format:")
  expect_match(detectedHtml, "DIA-NN")
  expect_match(detectedHtml, "Confidence: 91%")

  unknownHtml <- htmltools::renderTags(
    buildProtImportFormatDetectionUi(
      detectedFormat = "unknown",
      formatConfidence = 0.2
    )
  )$html

  expect_match(unknownHtml, "alert-danger")
  expect_match(unknownHtml, "Unknown format")
})

test_that("buildProtImportFormatSpecificOptionsUi renders per-format controls", {
  diannHtml <- htmltools::renderTags(
    buildProtImportFormatSpecificOptionsUi(
      format = "diann",
      ns = shiny::NS("import")
    )
  )$html

  expect_match(diannHtml, "DIA-NN Specific Options")
  expect_match(diannHtml, "import-diann_use_precursor_norm")

  spectronautHtml <- htmltools::renderTags(
    buildProtImportFormatSpecificOptionsUi(
      format = "spectronaut",
      ns = shiny::NS("import")
    )
  )$html

  expect_match(spectronautHtml, "Spectronaut Specific Options")
  expect_match(spectronautHtml, "import-spectronaut_quantity")

  fallbackHtml <- htmltools::renderTags(
    buildProtImportFormatSpecificOptionsUi(
      format = "unsupported",
      ns = shiny::NS("import")
    )
  )$html

  expect_match(fallbackHtml, "Format-specific options not available")
})

test_that("readProtImportHeaders reads plain, parquet, and zip/xlsx headers", {
  plainFile <- tempfile(fileext = ".tsv")
  writeLines(c("Protein.Group\tRun\tIntensity", "P1\tS1\t10"), plainFile)

  plainHeaders <- readProtImportHeaders(
    filePath = plainFile,
    readParquet = function(...) stop("unexpected parquet read")
  )

  expect_identical(plainHeaders, c("Protein.Group", "Run", "Intensity"))

  parquetHeaders <- readProtImportHeaders(
    filePath = "report.parquet",
    readParquet = function(path, col_select = NULL) {
      expect_identical(path, "report.parquet")
      data.frame(Protein.Group = character(), Run = character(), check.names = FALSE)
    },
    readLinesFn = function(...) stop("unexpected text read")
  )

  expect_identical(parquetHeaders, c("Protein.Group", "Run"))

  zipCalls <- list()
  xlsxHeaders <- readProtImportHeaders(
    filePath = "report.zip",
    unzipList = function(path) {
      expect_identical(path, "report.zip")
      data.frame(Name = "nested/report.xlsx", stringsAsFactors = FALSE)
    },
    unzipFiles = function(zipfile, files, exdir, junkpaths) {
      zipCalls[[length(zipCalls) + 1]] <<- list(
        zipfile = zipfile,
        files = files,
        exdir = exdir,
        junkpaths = junkpaths
      )
      invisible(NULL)
    },
    requireNamespaceFn = function(pkg, quietly = TRUE) {
      expect_identical(pkg, "readxl")
      TRUE
    },
    readExcel = function(path, n_max = 0) {
      expect_match(path, "report\\.xlsx$")
      data.frame(Protein.Group = character(), Run = character(), check.names = FALSE)
    },
    tempfileFn = function(...) tempfile("prot-import-zip-"),
    dirCreate = function(path) dir.create(path, recursive = TRUE),
    unlinkFn = function(path, recursive = TRUE) unlink(path, recursive = recursive),
    readParquet = function(...) stop("unexpected parquet read"),
    readLinesFn = function(...) stop("unexpected text read")
  )

  expect_identical(xlsxHeaders, c("Protein.Group", "Run"))
  expect_length(zipCalls, 1)
  expect_identical(zipCalls[[1]]$zipfile, "report.zip")
  expect_identical(zipCalls[[1]]$files, "nested/report.xlsx")
  expect_true(isTRUE(zipCalls[[1]]$junkpaths))
})

test_that("resolveProtImportDetectionFilename and resolveActiveProtImportFormat respect inputs", {
  expect_identical(
    resolveProtImportDetectionFilename(
      useShinyFiles = TRUE,
      filePath = "/tmp/report.tsv"
    ),
    "report.tsv"
  )

  expect_identical(
    resolveProtImportDetectionFilename(
      useShinyFiles = FALSE,
      filePath = "/tmp/report.tsv",
      searchResultsStandard = list(name = "upload.tsv")
    ),
    "upload.tsv"
  )

  expect_identical(resolveActiveProtImportFormat("auto", "diann"), "diann")
  expect_identical(resolveActiveProtImportFormat("fragpipe", "diann"), "fragpipe")
})

test_that("applyProtImportDetectedFormat and resetProtImportFormatDetectionState update local state", {
  localData <- list2env(list(detected_format = NULL, format_confidence = NULL))
  logInfoMessages <- character()
  logErrorMessages <- character()

  applyResult <- applyProtImportDetectedFormat(
    localData = localData,
    formatInfo = list(format = "diann", confidence = 0.91),
    logInfo = function(text) {
      logInfoMessages <<- c(logInfoMessages, text)
    }
  )

  expect_identical(localData$detected_format, "diann")
  expect_identical(localData$format_confidence, 0.91)
  expect_identical(applyResult$format, "diann")
  expect_match(logInfoMessages, "Detected format: diann")

  resetResult <- resetProtImportFormatDetectionState(
    localData = localData,
    errorMessage = "broken header read",
    logError = function(text) {
      logErrorMessages <<- c(logErrorMessages, text)
    }
  )

  expect_identical(localData$detected_format, "unknown")
  expect_identical(localData$format_confidence, 0)
  expect_identical(resetResult$format, "unknown")
  expect_identical(resetResult$confidence, 0)
  expect_identical(resetResult$errorMessage, "broken header read")
  expect_match(logErrorMessages, "Error detecting file format: broken header read")
})

test_that("runProtImportFormatDetection delegates success and error paths", {
  localData <- list2env(list(detected_format = NULL, format_confidence = NULL))
  applyCalls <- list()
  resetCalls <- list()

  success <- runProtImportFormatDetection(
    filePath = "/tmp/report.tsv",
    useShinyFiles = FALSE,
    localData = localData,
    searchResultsStandard = list(name = "upload.tsv"),
    readHeaders = function(path) {
      expect_identical(path, "/tmp/report.tsv")
      c("Protein.Group", "Run")
    },
    resolveFilename = function(useShinyFiles, filePath, searchResultsStandard) {
      expect_false(useShinyFiles)
      expect_identical(filePath, "/tmp/report.tsv")
      searchResultsStandard$name
    },
    detectFormat = function(headers, filename) {
      expect_identical(headers, c("Protein.Group", "Run"))
      expect_identical(filename, "upload.tsv")
      list(format = "diann", confidence = 0.9)
    },
    applyDetectedFormat = function(localData, formatInfo, logInfo) {
      applyCalls[[length(applyCalls) + 1]] <<- list(
        format = formatInfo$format,
        confidence = formatInfo$confidence,
        logInfo = identical(logInfo, identity)
      )
      "applied"
    },
    resetDetectionState = function(...) stop("should not reset"),
    logInfo = identity,
    logError = identity
  )

  errorResult <- runProtImportFormatDetection(
    filePath = "/tmp/report.tsv",
    useShinyFiles = TRUE,
    localData = localData,
    readHeaders = function(...) stop("header failure"),
    resolveFilename = function(...) stop("should not resolve"),
    detectFormat = function(...) stop("should not detect"),
    applyDetectedFormat = function(...) stop("should not apply"),
    resetDetectionState = function(localData, errorMessage, logError) {
      resetCalls[[length(resetCalls) + 1]] <<- list(
        errorMessage = errorMessage,
        logError = identical(logError, identity)
      )
      "reset"
    },
    logInfo = identity,
    logError = identity
  )

  expect_identical(success, "applied")
  expect_length(applyCalls, 1)
  expect_identical(applyCalls[[1]]$format, "diann")
  expect_true(isTRUE(applyCalls[[1]]$logInfo))
  expect_identical(errorResult, "reset")
  expect_length(resetCalls, 1)
  expect_identical(resetCalls[[1]]$errorMessage, "header failure")
  expect_true(isTRUE(resetCalls[[1]]$logError))
})

test_that("updateProtImportCheckpointCaptureOption updates option state and logs", {
  logMessages <- character()

  enabled <- withr::local_options(list(multischolar.capture_test_checkpoints = NULL))
  enabled <- updateProtImportCheckpointCaptureOption(
    captureCheckpoints = TRUE,
    optionsFn = options,
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    }
  )

  expect_true(isTRUE(enabled))
  expect_true(isTRUE(getOption("multischolar.capture_test_checkpoints")))
  expect_identical(logMessages, "Proteomics test checkpoint capture ENABLED")

  logMessages <- character()

  disabled <- updateProtImportCheckpointCaptureOption(
    captureCheckpoints = FALSE,
    optionsFn = options,
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    }
  )

  expect_false(disabled)
  expect_false(isTRUE(getOption("multischolar.capture_test_checkpoints")))
  expect_identical(logMessages, "Proteomics test checkpoint capture DISABLED")
})

test_that("toggleProtImportMixedSpeciesInputs toggles organism controls and messages", {
  disabledInputs <- character()
  enabledInputs <- character()
  messages <- character()

  enabled <- toggleProtImportMixedSpeciesInputs(
    mixedSpeciesFasta = TRUE,
    disableInput = function(id) {
      disabledInputs <<- c(disabledInputs, id)
    },
    enableInput = function(id) {
      enabledInputs <<- c(enabledInputs, id)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_true(isTRUE(enabled))
  expect_identical(disabledInputs, c("taxon_id", "organism_name"))
  expect_identical(enabledInputs, character())
  expect_identical(messages, "[mod_prot_import] Mixed species FASTA enabled - organism inputs disabled")

  disabledInputs <- character()
  enabledInputs <- character()
  messages <- character()

  disabled <- toggleProtImportMixedSpeciesInputs(
    mixedSpeciesFasta = FALSE,
    disableInput = function(id) {
      disabledInputs <<- c(disabledInputs, id)
    },
    enableInput = function(id) {
      enabledInputs <<- c(enabledInputs, id)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_false(disabled)
  expect_identical(disabledInputs, character())
  expect_identical(enabledInputs, c("taxon_id", "organism_name"))
  expect_identical(messages, "[mod_prot_import] Mixed species FASTA disabled - organism inputs enabled")
})

test_that("resolveProtImportPrimaryUploadPath returns shiny and standard file paths", {
  expect_identical(
    resolveProtImportPrimaryUploadPath(
      useShinyFiles = TRUE,
      shinyPath = "/tmp/from-shiny.tsv",
      standardInput = list(datapath = "/tmp/from-standard.tsv")
    ),
    "/tmp/from-shiny.tsv"
  )

  expect_identical(
    resolveProtImportPrimaryUploadPath(
      useShinyFiles = FALSE,
      shinyPath = "/tmp/from-shiny.tsv",
      standardInput = list(datapath = "/tmp/from-standard.tsv")
    ),
    "/tmp/from-standard.tsv"
  )

  expect_null(
    resolveProtImportPrimaryUploadPath(
      useShinyFiles = FALSE,
      shinyPath = NULL,
      standardInput = NULL
    )
  )
})

test_that("resolveProtImportShinyFileVolumes preserves explicit volumes and loads defaults", {
  explicitVolumes <- list(home = "/tmp/home")

  expect_identical(
    resolveProtImportShinyFileVolumes(
      volumes = explicitVolumes,
      getVolumes = function() {
        stop("should not load fallback volumes")
      }
    ),
    explicitVolumes
  )

  expect_identical(
    resolveProtImportShinyFileVolumes(
      volumes = NULL,
      getVolumes = function() {
        function() list(project = "/tmp/project")
      }
    ),
    list(project = "/tmp/project")
  )
})

test_that("registerProtImportShinyFileChooser delegates shinyFiles registration", {
  chooserCalls <- list()

  registered <- registerProtImportShinyFileChooser(
    input = "input-proxy",
    inputId = "search_results",
    volumes = list(home = "/tmp/home"),
    session = "session-proxy",
    filetypes = c("tsv", "txt"),
    chooseFile = function(input, inputId, roots, session, filetypes) {
      chooserCalls[[length(chooserCalls) + 1]] <<- list(
        input = input,
        inputId = inputId,
        roots = roots,
        session = session,
        filetypes = filetypes
      )
    }
  )

  expect_identical(registered, "search_results")
  expect_length(chooserCalls, 1)
  expect_identical(chooserCalls[[1]]$input, "input-proxy")
  expect_identical(chooserCalls[[1]]$inputId, "search_results")
  expect_identical(chooserCalls[[1]]$roots, list(home = "/tmp/home"))
  expect_identical(chooserCalls[[1]]$session, "session-proxy")
  expect_identical(chooserCalls[[1]]$filetypes, c("tsv", "txt"))
})

test_that("handleProtImportShinyFileSelection binds selected paths and catches parse errors", {
  localData <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  rendered <- character()

  selectedPath <- handleProtImportShinyFileSelection(
    selectedInput = list(selection = "report.tsv"),
    volumes = list(home = "/tmp/home"),
    localData = localData,
    localField = "search_results_file",
    output = output,
    outputId = "search_results_path",
    parseFilePaths = function(volumes, selectedInput) {
      data.frame(datapath = "/tmp/report.tsv", stringsAsFactors = FALSE)
    },
    renderText = function(text) {
      rendered <<- c(rendered, text)
      paste0("rendered:", text)
    },
    messageFn = function(text) {
      stop(sprintf("unexpected message: %s", text))
    }
  )

  expect_identical(selectedPath, "/tmp/report.tsv")
  expect_identical(localData$search_results_file, "/tmp/report.tsv")
  expect_identical(output$search_results_path, "rendered:/tmp/report.tsv")
  expect_identical(rendered, "/tmp/report.tsv")

  errorMessages <- character()
  errorLocalData <- new.env(parent = emptyenv())
  errorOutput <- new.env(parent = emptyenv())

  errorResult <- handleProtImportShinyFileSelection(
    selectedInput = list(selection = "report.tsv"),
    volumes = list(home = "/tmp/home"),
    localData = errorLocalData,
    localField = "search_results_file",
    output = errorOutput,
    outputId = "search_results_path",
    parseFilePaths = function(volumes, selectedInput) {
      stop("parse failed")
    },
    renderText = identity,
    messageFn = function(text) {
      errorMessages <<- c(errorMessages, text)
    },
    catchErrors = TRUE
  )

  expect_null(errorResult)
  expect_false(exists("search_results_file", envir = errorLocalData, inherits = FALSE))
  expect_false(exists("search_results_path", envir = errorOutput, inherits = FALSE))
  expect_identical(
    errorMessages,
    "   mod_prot_import_server ERROR in parseFilePaths: parse failed"
  )
})

test_that("buildProtImportProcessingModal renders the processing modal shell", {
  modalHtml <- htmltools::renderTags(
    buildProtImportProcessingModal(
      ns = shiny::NS("import"),
      initialStatus = "Loading FASTA..."
    )
  )$html

  expect_match(modalHtml, "Processing Data")
  expect_match(modalHtml, "Loading FASTA\\.\\.\\.")
  expect_match(modalHtml, "import-processing_status_text")
  expect_match(modalHtml, "Please wait while your data is being processed\\.")
})

test_that("updateProtImportProcessingStatus logs and updates the processing label", {
  updateCalls <- list()
  messages <- character()

  updateProtImportProcessingStatus(
    "Finalizing import...",
    updateHtml = function(id, html) {
      updateCalls[[length(updateCalls) + 1]] <<- list(id = id, html = html)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_length(updateCalls, 1)
  expect_identical(updateCalls[[1]]$id, "processing_status_text")
  expect_identical(updateCalls[[1]]$html, "Finalizing import...")
  expect_identical(messages, "[mod_prot_import] Finalizing import...")
})

test_that("importProtImportDataByFormat dispatches DIANN options and wraps import errors", {
  input <- list(diann_use_precursor_norm = FALSE)

  diannResult <- importProtImportDataByFormat(
    format = "diann",
    searchResultsPath = "report.tsv",
    input = input,
    importDiann = function(filepath, use_precursor_norm) {
      expect_identical(filepath, "report.tsv")
      expect_false(use_precursor_norm)
      list(data = data.frame(x = 1), data_type = "peptide", column_mapping = list())
    },
    importSpectronaut = function(...) stop("unexpected spectronaut call"),
    importFragpipe = function(...) stop("unexpected fragpipe call"),
    importMaxquant = function(...) stop("unexpected maxquant call"),
    importPdTmt = function(...) stop("unexpected pd_tmt call"),
    logError = function(...) NULL
  )

  expect_s3_class(diannResult$data, "data.frame")

  expect_error(
    importProtImportDataByFormat(
      format = "unknown",
      searchResultsPath = "report.tsv",
      input = list(),
      importDiann = function(...) NULL,
      importSpectronaut = function(...) NULL,
      importFragpipe = function(...) NULL,
      importMaxquant = function(...) NULL,
      importPdTmt = function(...) NULL,
      logError = function(...) NULL
    ),
    "Failed to import data: Unsupported format: unknown"
  )
})

test_that("sanitizeProtImportRunNames rewrites run column values when present", {
  workflowData <- list2env(list(
    data_tbl = data.frame(
      Run = c("Sample 1", "Sample 1", "Sample-2!"),
      value = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  ))

  notifications <- list()
  logMessages <- character()

  result <- sanitizeProtImportRunNames(
    workflowData = workflowData,
    runCol = "Run",
    makeCleanNames = function(x) gsub("[^a-z0-9]+", "_", tolower(x)),
    logInfo = function(text) logMessages <<- c(logMessages, text),
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    }
  )

  expect_true(result)
  expect_equal(workflowData$data_tbl$Run, c("sample_1", "sample_1", "sample_2_"))
  expect_identical(logMessages, "Sanitized 2 unique sample names.")
  expect_identical(notifications[[1]]$message, "Sample names sanitized for R compatibility.")
  expect_identical(notifications[[1]]$type, "message")

  unchangedWorkflow <- list2env(list(
    data_tbl = data.frame(
      Other = c("a", "b"),
      stringsAsFactors = FALSE
    )
  ))

  expect_false(
    sanitizeProtImportRunNames(
      workflowData = unchangedWorkflow,
      runCol = "Run",
      makeCleanNames = function(x) x,
      logInfo = function(...) stop("should not log"),
      showNotification = function(...) stop("should not notify")
    )
  )
})

test_that("applyProtImportResultToWorkflow stores import state and delegates sanitization", {
  stateManager <- new.env(parent = emptyenv())
  stateManager$workflow_type <- NULL
  stateManager$setWorkflowType <- function(type) {
    stateManager$workflow_type <- type
  }

  workflowData <- list2env(list(
    data_tbl = NULL,
    data_format = NULL,
    data_type = NULL,
    column_mapping = NULL,
    data_cln = NULL,
    fasta_file_path = NULL,
    state_manager = stateManager
  ))

  sanitizeCalls <- list()
  dataImportResult <- list(
    data = data.frame(
      Run = c("Sample 1", "Sample 2"),
      value = c(1, 2),
      stringsAsFactors = FALSE
    ),
    data_type = "peptide",
    column_mapping = list(run_col = "Run", protein_col = "Protein.Group")
  )

  workflowType <- applyProtImportResultToWorkflow(
    workflowData = workflowData,
    dataImportResult = dataImportResult,
    format = "diann",
    fastaPath = "/tmp/proteins.fasta",
    sanitizeNames = TRUE,
    logInfo = function(...) NULL,
    showNotification = function(...) NULL,
    sanitizeRunNames = function(workflowData, runCol, logInfo, showNotification) {
      sanitizeCalls[[length(sanitizeCalls) + 1]] <<- list(
        runCol = runCol,
        rows = nrow(workflowData$data_tbl)
      )
      workflowData$data_tbl$Run <- c("sample_1", "sample_2")
      TRUE
    }
  )

  expect_identical(workflowType, "DIA")
  expect_equal(workflowData$data_format, "diann")
  expect_equal(workflowData$data_type, "peptide")
  expect_equal(workflowData$column_mapping$run_col, "Run")
  expect_equal(workflowData$data_tbl$Run, c("sample_1", "sample_2"))
  expect_identical(workflowData$data_cln, workflowData$data_tbl)
  expect_equal(workflowData$fasta_file_path, "/tmp/proteins.fasta")
  expect_identical(stateManager$workflow_type, "DIA")
  expect_length(sanitizeCalls, 1)
  expect_identical(sanitizeCalls[[1]]$runCol, "Run")
  expect_identical(sanitizeCalls[[1]]$rows, 2L)

  workflowTypeNoSanitize <- applyProtImportResultToWorkflow(
    workflowData = workflowData,
    dataImportResult = dataImportResult,
    format = "fragpipe",
    fastaPath = "/tmp/proteins.fasta",
    sanitizeNames = FALSE,
    logInfo = function(...) stop("should not log"),
    showNotification = function(...) stop("should not notify"),
    sanitizeRunNames = function(...) stop("should not sanitize")
  )

  expect_identical(workflowTypeNoSanitize, "LFQ")
  expect_identical(stateManager$workflow_type, "LFQ")
})

test_that("finalizeProtImportSetupState records import summary and single-species metadata", {
  workflowData <- list2env(list(
    taxon_id = NULL,
    organism_name = NULL,
    mixed_species_analysis = NULL,
    processing_log = list()
  ))

  fixedTime <- as.POSIXct("2026-04-11 14:00:00", tz = "UTC")
  dataImportResult <- list(
    data = data.frame(
      Protein.Group = c("P1", "P2", "P2"),
      Precursor.Id = c("pep1", "pep2", "pep3"),
      Run = c("sample_1", "sample_1", "sample_2"),
      stringsAsFactors = FALSE
    ),
    data_type = "peptide",
    column_mapping = list(
      run_col = "Run",
      protein_col = "Protein.Group",
      peptide_col = "Precursor.Id"
    )
  )

  setupLog <- finalizeProtImportSetupState(
    workflowData = workflowData,
    dataImportResult = dataImportResult,
    format = "diann",
    searchFilename = "report.tsv",
    fastaFilename = "proteins.fasta",
    taxonId = 9606,
    organismName = "Homo sapiens",
    mixedSpeciesFasta = FALSE,
    now = function() fixedTime
  )

  expect_identical(workflowData$taxon_id, 9606)
  expect_identical(workflowData$organism_name, "Homo sapiens")
  expect_false(workflowData$mixed_species_analysis$enabled)
  expect_identical(workflowData$mixed_species_analysis$timestamp, fixedTime)
  expect_identical(setupLog$timestamp, fixedTime)
  expect_identical(setupLog$search_file, "report.tsv")
  expect_identical(setupLog$fasta_file, "proteins.fasta")
  expect_identical(setupLog$detected_format, "diann")
  expect_identical(setupLog$data_type, "peptide")
  expect_identical(setupLog$n_rows, 3L)
  expect_identical(setupLog$n_runs, 2L)
  expect_identical(setupLog$n_proteins, 2L)
  expect_identical(setupLog$n_peptides, 3L)

  workflowDataMixed <- list2env(list(
    taxon_id = NULL,
    organism_name = NULL,
    mixed_species_analysis = "keep-existing",
    processing_log = list()
  ))

  setupLogMixed <- finalizeProtImportSetupState(
    workflowData = workflowDataMixed,
    dataImportResult = list(
      data = data.frame(x = 1),
      data_type = "protein",
      column_mapping = list(run_col = NULL, protein_col = NULL, peptide_col = NULL)
    ),
    format = "fragpipe",
    searchFilename = "combined.tsv",
    fastaFilename = "proteins.fasta",
    taxonId = 10090,
    organismName = "Mus musculus",
    mixedSpeciesFasta = TRUE,
    now = function() fixedTime
  )

  expect_identical(workflowDataMixed$mixed_species_analysis, "keep-existing")
  expect_true(is.na(setupLogMixed$n_runs))
  expect_true(is.na(setupLogMixed$n_proteins))
  expect_true(is.na(setupLogMixed$n_peptides))
})

test_that("buildProtImportOrganismDistributionTable formats the organism distribution datatable", {
  tableWidget <- buildProtImportOrganismDistributionTable(
    data.frame(
      organism_name = c("Homo sapiens", "Mus musculus"),
      taxon_id = c(9606L, 10090L),
      protein_count = c(12L, 3L),
      percentage = c(80, 20),
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(tableWidget, "datatables")
  expect_named(tableWidget$x$data, c("Organism", "Taxon ID", "Proteins", "%"))
  expect_equal(as.character(tableWidget$x$data$Organism), c("Homo sapiens", "Mus musculus"))
  expect_equal(tableWidget$x$options$pageLength, 10)
  expect_false(tableWidget$x$options$searching)
  expect_false(tableWidget$x$options$ordering)
})

test_that("prepareProtImportOrganismSelectionChoices filters invalid rows and selects the top taxon", {
  selectionChoices <- prepareProtImportOrganismSelectionChoices(
    data.frame(
      organism_name = c("Homo sapiens", "[Unmatched/Unknown]", "Mus musculus"),
      taxon_id = c(9606L, NA, 10090L),
      protein_count = c(12L, 5L, 3L),
      percentage = c(80, 10, 20),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(selectionChoices$validOrganisms$organism_name, c("Homo sapiens", "Mus musculus"))
  expect_identical(unname(selectionChoices$choices), c("9606", "10090"))
  expect_identical(selectionChoices$selectedTaxon, "9606")
  expect_true(any(grepl("Homo sapiens \\(Taxon: 9606\\) - 12 proteins \\(80%\\)", names(selectionChoices$choices))))

  expect_null(
    prepareProtImportOrganismSelectionChoices(
      data.frame(
        organism_name = "[Unmatched/Unknown]",
        taxon_id = NA_integer_,
        protein_count = 1L,
        percentage = 100,
        stringsAsFactors = FALSE
      )
    )
  )
})

test_that("buildProtImportOrganismSelectionModal renders the organism-selection modal shell", {
  modalHtml <- htmltools::renderTags(
    buildProtImportOrganismSelectionModal(
      ns = shiny::NS("import"),
      organismDistribution = data.frame(
        organism_name = c("Homo sapiens", "Mus musculus"),
        taxon_id = c(9606L, 10090L),
        protein_count = c(12L, 3L),
        percentage = c(80, 20),
        stringsAsFactors = FALSE
      )
    )
  )$html

  expect_match(modalHtml, "Select Primary Organism")
  expect_match(modalHtml, "Multiple organisms detected in your FASTA database\\.")
  expect_match(modalHtml, "import-organism_dist_table")
  expect_match(modalHtml, "import-selected_organism")
  expect_match(modalHtml, "Filter data to keep only proteins from selected organism")
})

test_that("analyzeProtImportMixedSpeciesData populates distribution state and opens the selection modal", {
  workflowData <- list2env(list(organism_distribution = NULL))
  localData <- list2env(list(
    organism_mapping = NULL,
    organism_distribution = NULL,
    waiting_for_organism_selection = FALSE
  ))
  statusUpdates <- character()
  messages <- character()
  logMessages <- character()
  modalShown <- NULL
  notifications <- list()

  result <- analyzeProtImportMixedSpeciesData(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    fastaPath = "proteins.fasta",
    session = list(ns = shiny::NS("import")),
    updateStatus = function(text) {
      statusUpdates <<- c(statusUpdates, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logWarn = function(...) stop("should not warn"),
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    extractOrganisms = function(path) {
      expect_identical(path, "proteins.fasta")
      data.frame(
        uniprot_acc = c("P1", "P2"),
        organism_name = c("Homo sapiens", "Mus musculus"),
        taxon_id = c(9606L, 10090L),
        stringsAsFactors = FALSE
      )
    },
    analyzeDistribution = function(proteinIds, organismMapping) {
      expect_equal(sort(proteinIds), c("P1", "P2"))
      expect_equal(nrow(organismMapping), 2)
      data.frame(
        organism_name = c("Homo sapiens", "Mus musculus"),
        taxon_id = c(9606L, 10090L),
        protein_count = c(1L, 1L),
        percentage = c(50, 50),
        stringsAsFactors = FALSE
      )
    },
    buildSelectionModal = function(ns, organismDistribution) {
      expect_identical(ns("selected_organism"), "import-selected_organism")
      expect_equal(nrow(organismDistribution), 2)
      structure(list(modal = TRUE), class = "stub_modal")
    },
    showModal = function(modal) {
      modalShown <<- modal
    }
  )

  expect_true(result)
  expect_identical(statusUpdates, "Analyzing organism distribution...")
  expect_true(any(grepl("Mixed species FASTA - analyzing organism distribution", messages, fixed = TRUE)))
  expect_true(any(grepl("Found 2 organisms in data", logMessages, fixed = TRUE)))
  expect_equal(nrow(localData$organism_mapping), 2)
  expect_equal(nrow(localData$organism_distribution), 2)
  expect_equal(nrow(workflowData$organism_distribution), 2)
  expect_true(localData$waiting_for_organism_selection)
  expect_identical(class(modalShown), "stub_modal")
  expect_length(notifications, 0)
})

test_that("analyzeProtImportMixedSpeciesData warns when no valid organism choices can be built", {
  workflowData <- list2env(list(organism_distribution = NULL))
  localData <- list2env(list(
    organism_mapping = NULL,
    organism_distribution = NULL,
    waiting_for_organism_selection = FALSE
  ))
  warnMessages <- character()
  notifications <- list()

  result <- analyzeProtImportMixedSpeciesData(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    fastaPath = "proteins.fasta",
    session = list(ns = shiny::NS("import")),
    updateStatus = function(...) NULL,
    messageFn = function(...) NULL,
    logInfo = function(...) NULL,
    logWarn = function(text) {
      warnMessages <<- c(warnMessages, text)
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    extractOrganisms = function(...) data.frame(x = 1),
    analyzeDistribution = function(...) {
      data.frame(
        organism_name = "[Unmatched/Unknown]",
        taxon_id = NA_integer_,
        protein_count = 1L,
        percentage = 100,
        stringsAsFactors = FALSE
      )
    },
    buildSelectionModal = function(...) NULL,
    showModal = function(...) stop("should not show modal")
  )

  expect_false(result)
  expect_false(localData$waiting_for_organism_selection)
  expect_true(any(grepl("No valid organisms with taxon IDs found in FASTA", warnMessages, fixed = TRUE)))
  expect_identical(
    notifications[[1]]$message,
    "Could not extract organism information from FASTA headers. Using default organism settings."
  )
  expect_identical(notifications[[1]]$type, "warning")
})

test_that("resolveProtImportOptionalUploadPath and resolveProtImportInputFilename handle shiny and standard inputs", {
  expect_identical(
    resolveProtImportOptionalUploadPath(
      useShinyFiles = TRUE,
      shinyPath = "/tmp/config.ini",
      standardInput = list(datapath = "/tmp/ignored.ini")
    ),
    "/tmp/config.ini"
  )

  expect_identical(
    resolveProtImportOptionalUploadPath(
      useShinyFiles = FALSE,
      shinyPath = "/tmp/ignored.ini",
      standardInput = list(datapath = "/tmp/config.ini")
    ),
    "/tmp/config.ini"
  )

  expect_null(
    resolveProtImportOptionalUploadPath(
      useShinyFiles = FALSE,
      shinyPath = NULL,
      standardInput = NULL
    )
  )

  expect_identical(
    resolveProtImportInputFilename(
      useShinyFiles = TRUE,
      resolvedPath = "/tmp/report.tsv",
      standardInput = list(name = "ignored.tsv")
    ),
    "report.tsv"
  )

  expect_identical(
    resolveProtImportInputFilename(
      useShinyFiles = FALSE,
      resolvedPath = "/tmp/ignored.tsv",
      standardInput = list(name = "report.tsv")
    ),
    "report.tsv"
  )

  expect_null(
    resolveProtImportInputFilename(
      useShinyFiles = TRUE,
      resolvedPath = NULL,
      standardInput = NULL
    )
  )
})

test_that("loadProtImportOptionalMappings stores requested mapping tables", {
  workflowData <- list2env(list(
    uniprot_mapping = NULL,
    uniparc_mapping = NULL
  ))
  calls <- list()

  result <- loadProtImportOptionalMappings(
    workflowData = workflowData,
    uniprotPath = "uniprot.tsv",
    uniparcPath = "uniparc.tsv",
    readOptionalMapping = function(mappingPath, mappingLabel, readMapping, logInfo, logError) {
      calls[[length(calls) + 1]] <<- list(
        mappingPath = mappingPath,
        mappingLabel = mappingLabel
      )
      data.frame(label = mappingLabel, stringsAsFactors = FALSE)
    },
    logInfo = function(...) NULL,
    logError = function(...) NULL
  )

  expect_length(calls, 2)
  expect_identical(calls[[1]]$mappingPath, "uniprot.tsv")
  expect_identical(calls[[1]]$mappingLabel, "UniProt mapping")
  expect_identical(calls[[2]]$mappingPath, "uniparc.tsv")
  expect_identical(calls[[2]]$mappingLabel, "UniParc mapping")
  expect_identical(workflowData$uniprot_mapping$label, "UniProt mapping")
  expect_identical(workflowData$uniparc_mapping$label, "UniParc mapping")
  expect_identical(result$uniprotMapping$label, "UniProt mapping")
  expect_identical(result$uniparcMapping$label, "UniParc mapping")

  emptyWorkflow <- list2env(list(
    uniprot_mapping = "keep-uniprot",
    uniparc_mapping = "keep-uniparc"
  ))

  emptyResult <- loadProtImportOptionalMappings(
    workflowData = emptyWorkflow,
    uniprotPath = NULL,
    uniparcPath = NULL,
    readOptionalMapping = function(...) stop("should not read"),
    logInfo = function(...) NULL,
    logError = function(...) NULL
  )

  expect_identical(emptyWorkflow$uniprot_mapping, "keep-uniprot")
  expect_identical(emptyWorkflow$uniparc_mapping, "keep-uniparc")
  expect_identical(emptyResult$uniprotMapping, "keep-uniprot")
  expect_identical(emptyResult$uniparcMapping, "keep-uniparc")
})

test_that("announceProtImportStart logs the import start banner and key paths", {
  messages <- character()

  result <- announceProtImportStart(
    searchResultsPath = "/tmp/report.tsv",
    fastaPath = "/tmp/proteins.fasta",
    format = "diann",
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_identical(result$searchResultsPath, "/tmp/report.tsv")
  expect_identical(result$fastaPath, "/tmp/proteins.fasta")
  expect_identical(result$format, "diann")
  expect_equal(
    messages[1:3],
    c(
      "========================================",
      "[mod_prot_import] Starting data import process",
      "========================================"
    )
  )
  expect_true(any(grepl("\\[mod_prot_import\\] Search results: /tmp/report.tsv", messages)))
  expect_true(any(grepl("\\[mod_prot_import\\] FASTA file: /tmp/proteins.fasta", messages)))
  expect_true(any(grepl("\\[mod_prot_import\\] Detected format: diann", messages)))
})

test_that("recordProtImportImportedData logs, checkpoints, and reports row count", {
  logMessages <- character()
  messages <- character()
  checkpointCalls <- list()

  result <- recordProtImportImportedData(
    dataImportResult = mockImportResult(),
    captureCheckpoint = function(result, checkpointId, checkpointName) {
      checkpointCalls[[length(checkpointCalls) + 1]] <<- list(
        rows = nrow(result$data),
        checkpointId = checkpointId,
        checkpointName = checkpointName
      )
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_equal(nrow(result$data), 2)
  expect_identical(logMessages, "Data imported successfully. Rows: 2")
  expect_length(checkpointCalls, 1)
  expect_identical(checkpointCalls[[1]]$rows, 2L)
  expect_identical(checkpointCalls[[1]]$checkpointId, "cp01")
  expect_identical(checkpointCalls[[1]]$checkpointName, "raw_imported")
  expect_identical(messages, "[mod_prot_import] Data imported: 2 rows")
})

test_that("recordProtImportSetupLog resolves filenames and finalizes setup state", {
  workflowData <- list2env(list())
  logMessages <- character()
  finalizeCalls <- list()

  result <- recordProtImportSetupLog(
    workflowData = workflowData,
    dataImportResult = mockImportResult(),
    format = "diann",
    useShinyFiles = FALSE,
    searchResultsPath = "/tmp/report.tsv",
    fastaPath = "/tmp/proteins.fasta",
    searchResultsStandard = list(name = "report.tsv"),
    fastaFileStandard = list(name = "proteins.fasta"),
    taxonId = 9606,
    organismName = "Homo sapiens",
    mixedSpeciesFasta = FALSE,
    resolveFilename = function(useShinyFiles, resolvedPath, standardInput) {
      filename <- if (!is.null(standardInput) && !is.null(standardInput$name)) {
        standardInput$name
      } else {
        basename(resolvedPath)
      }
      paste0(if (isTRUE(useShinyFiles)) "shiny:" else "standard:", filename)
    },
    finalizeSetupState = function(...) {
      finalizeCalls[[length(finalizeCalls) + 1]] <<- list(...)
      invisible(NULL)
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    }
  )

  expect_identical(result$searchFilename, "standard:report.tsv")
  expect_identical(result$fastaFilename, "standard:proteins.fasta")
  expect_identical(logMessages, c("Creating processing log...", "Processing log created successfully"))
  expect_length(finalizeCalls, 1)
  expect_identical(finalizeCalls[[1]]$format, "diann")
  expect_identical(finalizeCalls[[1]]$searchFilename, "standard:report.tsv")
  expect_identical(finalizeCalls[[1]]$fastaFilename, "standard:proteins.fasta")
  expect_identical(finalizeCalls[[1]]$taxonId, 9606)
  expect_identical(finalizeCalls[[1]]$organismName, "Homo sapiens")
  expect_false(finalizeCalls[[1]]$mixedSpeciesFasta)
})

test_that("startProtImportProcessing requires paths, marks state, and shows modal", {
  localData <- list2env(list(processing = FALSE))
  requiredValues <- list()
  announceCalls <- list()
  modalCalls <- list()

  result <- startProtImportProcessing(
    searchResultsPath = "/tmp/report.tsv",
    fastaPath = "/tmp/proteins.fasta",
    format = "diann",
    localData = localData,
    ns = function(id) paste0("ns-", id),
    requireValue = function(value) {
      requiredValues[[length(requiredValues) + 1]] <<- value
      invisible(value)
    },
    announceStart = function(searchResultsPath, fastaPath, format) {
      announceCalls[[length(announceCalls) + 1]] <<- list(
        searchResultsPath = searchResultsPath,
        fastaPath = fastaPath,
        format = format
      )
      invisible(NULL)
    },
    buildProcessingModal = function(ns) list(ns = ns("processing_status")),
    showModal = function(modal) {
      modalCalls[[length(modalCalls) + 1]] <<- modal
      invisible(NULL)
    }
  )

  expect_identical(requiredValues, list("/tmp/report.tsv", "/tmp/proteins.fasta"))
  expect_true(isTRUE(localData$processing))
  expect_length(announceCalls, 1)
  expect_identical(announceCalls[[1]]$format, "diann")
  expect_identical(modalCalls[[1]]$ns, "ns-processing_status")
  expect_identical(result$searchResultsPath, "/tmp/report.tsv")
  expect_identical(result$fastaPath, "/tmp/proteins.fasta")
  expect_identical(result$format, "diann")
})

test_that("loadProtImportConfigurationResources resolves paths and delegates config setup", {
  workflowData <- list2env(list())
  statusMessages <- character()
  messages <- character()
  resolvedCalls <- list()
  storedConfigs <- list()
  mappingCalls <- list()

  result <- loadProtImportConfigurationResources(
    workflowData = workflowData,
    experimentPaths = list(source_dir = "/tmp/project"),
    useShinyFiles = FALSE,
    configFilePath = "/ignored/config.ini",
    configFileStandard = list(datapath = "/tmp/config.ini"),
    uniprotMappingFile = "/ignored/uniprot.tsv",
    uniprotMappingStandard = list(datapath = "/tmp/uniprot.tsv"),
    uniparcMappingFile = "/ignored/uniparc.tsv",
    uniparcMappingStandard = list(datapath = "/tmp/uniparc.tsv"),
    updateStatus = function(text) {
      statusMessages <<- c(statusMessages, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    },
    resolveUploadPath = function(useShinyFiles, shinyPath, standardInput) {
      resolvedCalls[[length(resolvedCalls) + 1]] <<- list(
        useShinyFiles = useShinyFiles,
        shinyPath = shinyPath,
        standardInput = standardInput
      )
      standardInput$datapath
    },
    loadConfiguration = function(configPath, experimentPaths, ...) {
      list(path = configPath, sourceDir = experimentPaths$source_dir)
    },
    storeConfiguration = function(workflowData, configList, ...) {
      workflowData$config_list <- configList
      storedConfigs[[length(storedConfigs) + 1]] <<- configList
      invisible(NULL)
    },
    loadOptionalMappings = function(workflowData, uniprotPath, uniparcPath, ...) {
      mappingCalls[[length(mappingCalls) + 1]] <<- list(
        uniprotPath = uniprotPath,
        uniparcPath = uniparcPath
      )
      workflowData$uniprot_mapping <- data.frame(path = uniprotPath, stringsAsFactors = FALSE)
      workflowData$uniparc_mapping <- data.frame(path = uniparcPath, stringsAsFactors = FALSE)
      invisible(NULL)
    }
  )

  expect_identical(statusMessages, "Loading configuration...")
  expect_identical(messages, "[mod_prot_import] Loading configuration")
  expect_length(resolvedCalls, 3)
  expect_identical(storedConfigs[[1]]$path, "/tmp/config.ini")
  expect_identical(workflowData$config_list$sourceDir, "/tmp/project")
  expect_length(mappingCalls, 1)
  expect_identical(mappingCalls[[1]]$uniprotPath, "/tmp/uniprot.tsv")
  expect_identical(mappingCalls[[1]]$uniparcPath, "/tmp/uniparc.tsv")
  expect_identical(result$configPath, "/tmp/config.ini")
  expect_identical(result$uniprotPath, "/tmp/uniprot.tsv")
  expect_identical(result$uniparcPath, "/tmp/uniparc.tsv")
})

test_that("finalizeProtImportProcessing updates status and delegates success completion", {
  statusMessages <- character()
  messages <- character()
  completionCalls <- list()
  workflowData <- list2env(list())
  localData <- list2env(list(processing = TRUE))

  result <- finalizeProtImportProcessing(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    format = "diann",
    updateStatus = function(text) {
      statusMessages <<- c(statusMessages, text)
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    },
    completeSuccessState = function(workflowData, localData, dataImportResult, format, messageFn) {
      completionCalls[[length(completionCalls) + 1]] <<- list(
        rows = nrow(dataImportResult$data),
        format = format
      )
      messageFn("inner-success")
      invisible("complete")
    }
  )

  expect_identical(statusMessages, "Finalizing import...")
  expect_identical(messages, c("[mod_prot_import] Finalizing import", "inner-success"))
  expect_length(completionCalls, 1)
  expect_identical(completionCalls[[1]]$rows, 2L)
  expect_identical(completionCalls[[1]]$format, "diann")
  expect_identical(result$updatedStatus, "complete")
  expect_identical(result$format, "diann")
})

test_that("readProtImportDataWithStatus updates status and records imported data", {
  statusMessages <- character()
  logMessages <- character()
  importCalls <- list()
  recordCalls <- list()

  result <- readProtImportDataWithStatus(
    format = "diann",
    searchResultsPath = "/tmp/report.tsv",
    input = list(source = "mock"),
    updateStatus = function(text) {
      statusMessages <<- c(statusMessages, text)
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    importDataByFormat = function(format, searchResultsPath, input, ...) {
      importCalls[[length(importCalls) + 1]] <<- list(
        format = format,
        searchResultsPath = searchResultsPath,
        input = input
      )
      mockImportResult()
    },
    recordImportedData = function(dataImportResult, captureCheckpoint, logInfo, messageFn) {
      recordCalls[[length(recordCalls) + 1]] <<- list(
        rows = nrow(dataImportResult$data),
        checkpoint = identical(captureCheckpoint, identity),
        logInfo = identical(logInfo, environment()$logInfo),
        messageFn = identical(messageFn, identity)
      )
      invisible(dataImportResult)
    },
    captureCheckpoint = identity,
    messageFn = identity
  )

  expect_identical(statusMessages, "Reading DIANN data...")
  expect_identical(logMessages, "Reading diann data from /tmp/report.tsv")
  expect_length(importCalls, 1)
  expect_identical(importCalls[[1]]$format, "diann")
  expect_identical(importCalls[[1]]$searchResultsPath, "/tmp/report.tsv")
  expect_length(recordCalls, 1)
  expect_identical(recordCalls[[1]]$rows, 2L)
  expect_equal(nrow(result$data), 2)
})

test_that("applyProtImportWorkflowWithStatus updates status and delegates apply plus FASTA", {
  statusMessages <- character()
  applyCalls <- list()
  fastaCalls <- list()

  result <- applyProtImportWorkflowWithStatus(
    workflowData = list2env(list()),
    dataImportResult = mockImportResult(),
    format = "diann",
    fastaPath = "/tmp/proteins.fasta",
    organismName = "Homo sapiens",
    experimentPaths = list(source_dir = "/tmp/project"),
    sanitizeNames = TRUE,
    updateStatus = function(text) {
      statusMessages <<- c(statusMessages, text)
    },
    applyImportResult = function(...) {
      applyCalls[[length(applyCalls) + 1]] <<- list(...)
      invisible(NULL)
    },
    processFastaData = function(...) {
      fastaCalls[[length(fastaCalls) + 1]] <<- list(...)
      invisible(NULL)
    },
    logInfo = function(...) NULL,
    logWarn = function(...) NULL,
    showNotification = function(...) NULL,
    sanitizeRunNames = function(x) x,
    processFasta = function(...) NULL
  )

  expect_identical(
    statusMessages,
    c("Storing imported data...", "Sanitizing sample names...", "Processing FASTA file...")
  )
  expect_length(applyCalls, 1)
  expect_true(isTRUE(applyCalls[[1]]$sanitizeNames))
  expect_identical(applyCalls[[1]]$fastaPath, "/tmp/proteins.fasta")
  expect_length(fastaCalls, 1)
  expect_identical(fastaCalls[[1]]$organismName, "Homo sapiens")
  expect_identical(result$sanitizeNames, TRUE)
  expect_identical(result$fastaPath, "/tmp/proteins.fasta")
})

test_that("runProtImportMixedSpeciesAnalysisIfNeeded only runs when analysis is eligible", {
  workflowData <- list2env(list(aa_seq_tbl_final = NULL))
  localData <- list2env(list())
  analyzeCalls <- list()

  skipped <- runProtImportMixedSpeciesAnalysisIfNeeded(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    fastaPath = "/tmp/proteins.fasta",
    session = structure(list(), class = "MockShinySession"),
    mixedSpeciesFasta = TRUE,
    analyzeMixedSpeciesData = function(...) {
      analyzeCalls[[length(analyzeCalls) + 1]] <<- list(...)
      TRUE
    }
  )

  workflowData$aa_seq_tbl_final <- data.frame(id = "P1", stringsAsFactors = FALSE)

  analyzed <- runProtImportMixedSpeciesAnalysisIfNeeded(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    fastaPath = "/tmp/proteins.fasta",
    session = structure(list(), class = "MockShinySession"),
    mixedSpeciesFasta = TRUE,
    analyzeMixedSpeciesData = function(...) {
      analyzeCalls[[length(analyzeCalls) + 1]] <<- list(...)
      "analyzed"
    }
  )

  expect_false(skipped)
  expect_identical(analyzed, "analyzed")
  expect_length(analyzeCalls, 1)
  expect_identical(analyzeCalls[[1]]$fastaPath, "/tmp/proteins.fasta")
})

test_that("runProtImportProcessingSafely returns success and resets state on error", {
  workflowData <- list2env(list())
  localData <- list2env(list())
  resetCalls <- list()

  success <- runProtImportProcessingSafely(
    runProcessing = function() "ok",
    workflowData = workflowData,
    localData = localData,
    logError = function(...) NULL,
    resetWorkflowState = function(...) stop("should not reset")
  )

  errorResult <- runProtImportProcessingSafely(
    runProcessing = function() stop("boom"),
    workflowData = workflowData,
    localData = localData,
    logError = function(text) text,
    resetWorkflowState = function(...) {
      resetCalls[[length(resetCalls) + 1]] <<- list(...)
      "reset"
    }
  )

  expect_identical(success, "ok")
  expect_identical(errorResult, "reset")
  expect_length(resetCalls, 1)
  expect_identical(resetCalls[[1]]$errorMessage, "boom")
})

test_that("confirmProtImportOrganismSelection updates state and filters selected organism", {
  importTable <- data.frame(
    Protein.Group = c("P1;P3", "P3", "P2-2"),
    value = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  workflowData <- list2env(list(
    data_tbl = importTable,
    data_cln = importTable,
    column_mapping = list(protein_col = "Protein.Group"),
    processing_log = list(
      setup_import = list(
        taxon_id = 10090L,
        organism = "Mus musculus"
      )
    ),
    mixed_species_analysis = NULL,
    taxon_id = NULL,
    organism_name = NULL
  ))

  localData <- list2env(list(
    organism_distribution = data.frame(
      organism_name = c("Homo sapiens", "Mus musculus"),
      taxon_id = c(9606L, 10090L),
      protein_count = c(2L, 1L),
      percentage = c(66.67, 33.33),
      stringsAsFactors = FALSE
    ),
    organism_mapping = data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      organism_name = c("Homo sapiens", "Homo sapiens", "Mus musculus"),
      taxon_id = c(9606L, 9606L, 10090L),
      stringsAsFactors = FALSE
    ),
    waiting_for_organism_selection = TRUE
  ))

  notifications <- list()
  numericUpdates <- list()
  textUpdates <- list()
  modalRemoved <- FALSE
  fixedTime <- as.POSIXct("2026-04-11 12:00:00", tz = "UTC")

  result <- confirmProtImportOrganismSelection(
    selectedTaxon = 9606L,
    filterToOrganism = TRUE,
    workflowData = workflowData,
    localData = localData,
    session = structure(list(), class = "MockShinySession"),
    updateNumericInput = function(session, inputId, value) {
      numericUpdates[[length(numericUpdates) + 1]] <<- list(
        inputId = inputId,
        value = value
      )
    },
    updateTextInput = function(session, inputId, value) {
      textUpdates[[length(textUpdates) + 1]] <<- list(
        inputId = inputId,
        value = value
      )
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    removeModal = function() {
      modalRemoved <<- TRUE
    },
    now = function() fixedTime,
    logInfo = function(...) NULL,
    normalizeUniprot = function(acc, remove_isoform = TRUE) sub("-\\d+$", "", acc),
    cleanAcc = function(acc) sub("-\\d+$", "", acc)
  )

  expect_true(result)
  expect_identical(workflowData$taxon_id, 9606L)
  expect_identical(workflowData$organism_name, "Homo sapiens")
  expect_equal(nrow(workflowData$data_tbl), 2)
  expect_equal(nrow(workflowData$data_cln), 2)
  expect_equal(workflowData$data_tbl$Protein.Group, c("P1;P3", "P2-2"))
  expect_false(localData$waiting_for_organism_selection)
  expect_true(modalRemoved)

  expect_identical(workflowData$processing_log$setup_import$taxon_id, 9606L)
  expect_identical(workflowData$processing_log$setup_import$organism, "Homo sapiens")
  expect_true(workflowData$processing_log$setup_import$mixed_species_selection$filter_applied)
  expect_identical(workflowData$mixed_species_analysis$selected_taxon_id, 9606L)
  expect_identical(workflowData$mixed_species_analysis$selected_organism, "Homo sapiens")
  expect_identical(workflowData$mixed_species_analysis$timestamp, fixedTime)

  expect_equal(numericUpdates[[1]]$inputId, "taxon_id")
  expect_identical(numericUpdates[[1]]$value, 9606L)
  expect_equal(textUpdates[[1]]$inputId, "organism_name")
  expect_identical(textUpdates[[1]]$value, "Homo sapiens")

  notificationMessages <- vapply(notifications, `[[`, character(1), "message")
  expect_true(any(grepl("Filtered to Homo sapiens: kept 2 rows, removed 1 rows", notificationMessages, fixed = TRUE)))
  expect_true(any(grepl("Primary organism set to: Homo sapiens \\(Taxon: 9606\\)", notificationMessages)))
  expect_true(any(notificationMessages == "Data import successful!"))
})

test_that("processProtImportFastaData stores FASTA results and persists helper outputs", {
  fixture <- makeImportFixture(withConfig = FALSE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  workflowData <- list2env(list(
    uniprot_mapping = data.frame(uniprot_id = "P1", stringsAsFactors = FALSE),
    uniparc_mapping = data.frame(uniparc_id = "UPI1", stringsAsFactors = FALSE),
    aa_seq_tbl_final = NULL,
    fasta_metadata = NULL
  ))
  assignEnv <- new.env(parent = emptyenv())
  savePaths <- character()
  logMessages <- character()
  warnMessages <- character()

  result <- processProtImportFastaData(
    workflowData = workflowData,
    fastaPath = fixture$fastaPath,
    organismName = "Homo sapiens",
    experimentPaths = fixture$experimentPaths,
    processFasta = function(fasta_file_path,
                            uniprot_search_results,
                            uniparc_search_results,
                            fasta_meta_file,
                            organism_name) {
      expect_identical(fasta_file_path, fixture$fastaPath)
      expect_identical(organism_name, "Homo sapiens")
      expect_equal(uniprot_search_results$uniprot_id, "P1")
      expect_equal(uniparc_search_results$uniparc_id, "UPI1")
      expect_true(grepl("aa_seq_tbl\\.RDS$", fasta_meta_file))
      mockFastaResult()
    },
    assignFn = assign,
    assignEnv = assignEnv,
    saveRds = function(object, file) {
      savePaths <<- c(savePaths, file)
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logWarn = function(text) {
      warnMessages <<- c(warnMessages, text)
    }
  )

  expect_true(result$success)
  expect_true(dir.exists(result$cacheDir))
  expect_equal(basename(savePaths), c("aa_seq_tbl_final.RDS", "fasta_metadata.RDS"))
  expect_equal(workflowData$fasta_metadata$num_sequences, 2)
  expect_equal(workflowData$aa_seq_tbl_final$uniprot_acc, c("P1", "P2"))
  expect_identical(get("aa_seq_tbl_final", envir = assignEnv), workflowData$aa_seq_tbl_final)
  expect_true(any(grepl("Checking experiment_paths...", logMessages, fixed = TRUE)))
  expect_true(any(grepl("FASTA file processed successfully. Found 2 sequences", logMessages, fixed = TRUE)))
  expect_length(warnMessages, 0)
})

test_that("processProtImportFastaData clears FASTA state and reports failure when parsing fails", {
  workflowData <- list2env(list(
    uniprot_mapping = NULL,
    uniparc_mapping = NULL,
    aa_seq_tbl_final = data.frame(old = 1),
    fasta_metadata = list(old = TRUE)
  ))
  warnMessages <- character()

  result <- processProtImportFastaData(
    workflowData = workflowData,
    fastaPath = "broken.fasta",
    organismName = "Homo sapiens",
    experimentPaths = NULL,
    processFasta = function(...) {
      stop("parse failed")
    },
    assignFn = function(...) {
      stop("should not assign")
    },
    saveRds = function(...) {
      stop("should not save")
    },
    logInfo = function(...) NULL,
    logWarn = function(text) {
      warnMessages <<- c(warnMessages, text)
    },
    dirExists = function(...) TRUE
  )

  expect_false(result$success)
  expect_match(result$error, "parse failed")
  expect_true(grepl("proteomics_cache$", result$cacheDir))
  expect_null(workflowData$aa_seq_tbl_final)
  expect_null(workflowData$fasta_metadata)
  expect_true(any(grepl("Error processing FASTA file: parse failed", warnMessages, fixed = TRUE)))
  expect_true(any(grepl("Continuing without FASTA processing", warnMessages, fixed = TRUE)))
})

test_that("completeProtImportSuccessState marks setup complete and closes the modal when ready", {
  workflowData <- list2env(list(
    tab_status = list(setup_import = "pending"),
    processing_log = list(
      setup_import = list(n_proteins = 2L)
    )
  ))
  localData <- list2env(list(
    processing = TRUE,
    waiting_for_organism_selection = FALSE
  ))
  notifications <- list()
  modalRemoved <- FALSE
  messages <- character()

  result <- completeProtImportSuccessState(
    workflowData = workflowData,
    localData = localData,
    dataImportResult = mockImportResult(),
    format = "diann",
    removeModal = function() {
      modalRemoved <<- TRUE
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_equal(workflowData$tab_status$setup_import, "complete")
  expect_identical(result$setup_import, "complete")
  expect_false(localData$processing)
  expect_true(modalRemoved)
  expect_identical(notifications[[1]]$message, "Data import successful!")
  expect_identical(notifications[[1]]$type, "message")
  expect_true(any(grepl("\\[mod_prot_import\\] Data import completed successfully!", messages)))
  expect_true(any(grepl("Rows: 2, Proteins: 2, Format: diann", messages, fixed = TRUE)))
})

test_that("resetProtImportWorkflowStateOnError clears import state and marks setup incomplete", {
  workflowData <- list2env(list(
    data_tbl = data.frame(x = 1),
    data_format = "diann",
    data_type = "peptide",
    column_mapping = list(run_col = "Run"),
    data_cln = data.frame(x = 1),
    fasta_file_path = "proteins.fasta",
    aa_seq_tbl_final = data.frame(x = 1),
    config_list = list(a = 1),
    processing_log = list(setup_import = list(done = TRUE)),
    tab_status = list(setup_import = "complete")
  ))
  localData <- list2env(list(processing = TRUE))
  notifications <- list()
  loggedErrors <- character()
  messages <- character()
  modalRemoved <- FALSE

  result <- resetProtImportWorkflowStateOnError(
    workflowData = workflowData,
    localData = localData,
    errorMessage = "boom",
    logError = function(text) {
      loggedErrors <<- c(loggedErrors, text)
    },
    removeModal = function() {
      modalRemoved <<- TRUE
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    messageFn = function(text) {
      messages <<- c(messages, text)
    }
  )

  expect_null(workflowData$data_tbl)
  expect_null(workflowData$data_format)
  expect_null(workflowData$data_type)
  expect_null(workflowData$column_mapping)
  expect_null(workflowData$data_cln)
  expect_null(workflowData$fasta_file_path)
  expect_null(workflowData$aa_seq_tbl_final)
  expect_null(workflowData$config_list)
  expect_null(workflowData$processing_log$setup_import)
  expect_equal(workflowData$tab_status$setup_import, "incomplete")
  expect_identical(result$setup_import, "incomplete")
  expect_false(localData$processing)
  expect_true(modalRemoved)
  expect_identical(loggedErrors, "Error during data import: boom")
  expect_identical(notifications[[1]]$message, "Error: boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
  expect_true(any(grepl("\\[mod_prot_import\\] ERROR: boom", messages)))
})

test_that("loadProtImportConfiguration reads explicit config paths and package defaults", {
  fixture <- makeImportFixture(withConfig = TRUE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  logMessages <- character()
  notifications <- list()

  explicitConfig <- loadProtImportConfiguration(
    configPath = fixture$configPath,
    experimentPaths = fixture$experimentPaths,
    readConfig = function(file) {
      expect_identical(file, fixture$configPath)
      list(config_source = "explicit")
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logError = function(...) stop("should not log errors"),
    getDefaultConfig = function() stop("should not fall back")
  )

  expect_identical(explicitConfig$config_source, "explicit")
  expect_true(any(grepl(paste("Reading configuration from", fixture$configPath), logMessages, fixed = TRUE)))
  expect_length(notifications, 0)

  pkgConfigPath <- file.path(fixture$rootDir, "pkg-config.ini")
  writeLines("[generalParameters]\nmin_peptides_per_protein=4", pkgConfigPath)

  file.remove(fixture$configPath)
  logMessages <- character()

  packageDefaultConfig <- loadProtImportConfiguration(
    configPath = NULL,
    experimentPaths = fixture$experimentPaths,
    readConfig = function(file) {
      expect_identical(file, file.path(fixture$experimentPaths$source_dir, "config.ini"))
      list(config_source = "package-default")
    },
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logError = function(...) stop("should not log errors"),
    getDefaultConfig = function() stop("should not fall back"),
    fileExists = file.exists,
    fileCopy = function(from, to) {
      file.copy(from, to, overwrite = TRUE)
    },
    downloadFile = function(...) stop("should not download"),
    systemFileFn = function(...) pkgConfigPath
  )

  expect_identical(packageDefaultConfig$config_source, "package-default")
  expect_true(file.exists(file.path(fixture$experimentPaths$source_dir, "config.ini")))
  expect_true(any(grepl("Found default config.ini in package:", logMessages, fixed = TRUE)))
  expect_true(any(vapply(
    notifications,
    function(entry) identical(entry$message, "Copied default config.ini to scripts directory."),
    logical(1)
  )))
})

test_that("loadProtImportConfiguration falls back to minimal defaults when retrieval fails", {
  fixture <- makeImportFixture(withConfig = FALSE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  logMessages <- character()
  errorMessages <- character()
  notifications <- list()
  fallbackConfig <- list(config_source = "fallback")

  config <- loadProtImportConfiguration(
    configPath = NULL,
    experimentPaths = fixture$experimentPaths,
    readConfig = function(...) stop("should not read config"),
    showNotification = function(message, type = "default", duration = NULL) {
      notifications[[length(notifications) + 1]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logError = function(text) {
      errorMessages <<- c(errorMessages, text)
    },
    getDefaultConfig = function() fallbackConfig,
    fileExists = function(path) FALSE,
    fileCopy = function(...) stop("should not copy"),
    downloadFile = function(...) stop("offline"),
    systemFileFn = function(...) ""
  )

  expect_identical(config, fallbackConfig)
  expect_true(any(grepl("Using default configuration", logMessages, fixed = TRUE)))
  expect_true(any(grepl("Using minimal fallback configuration", logMessages, fixed = TRUE)))
  expect_identical(errorMessages, "Failed to retrieve default config.ini: offline")
  expect_true(any(vapply(
    notifications,
    function(entry) identical(entry$type, "warning") && identical(entry$duration, 10),
    logical(1)
  )))
})

test_that("storeProtImportConfiguration saves config in workflow and compatibility env", {
  workflowData <- list2env(list(config_list = NULL))
  assignEnv <- new.env(parent = emptyenv())
  logMessages <- character()
  configList <- list(generalParameters = list(min_peptides_per_protein = 2))

  result <- storeProtImportConfiguration(
    workflowData = workflowData,
    configList = configList,
    assignFn = assign,
    assignEnv = assignEnv,
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    }
  )

  expect_identical(workflowData$config_list, configList)
  expect_identical(get("config_list", envir = assignEnv), configList)
  expect_identical(result, configList)
  expect_identical(logMessages, "Created global config_list for compatibility")
})

test_that("readProtImportOptionalMapping returns NULL without a path and wraps reader errors", {
  expect_null(
    readProtImportOptionalMapping(
      mappingPath = NULL,
      mappingLabel = "UniProt mapping",
      readMapping = function(...) stop("should not read"),
      logInfo = function(...) stop("should not log"),
      logError = function(...) stop("should not log")
    )
  )

  logMessages <- character()
  mapping <- readProtImportOptionalMapping(
    mappingPath = "mapping.tsv",
    mappingLabel = "UniParc mapping",
    readMapping = function(path) {
      expect_identical(path, "mapping.tsv")
      data.frame(id = "UPI1", stringsAsFactors = FALSE)
    },
    logInfo = function(text) {
      logMessages <<- c(logMessages, text)
    },
    logError = function(...) stop("should not log errors")
  )

  expect_equal(mapping$id, "UPI1")
  expect_identical(logMessages, "Reading UniParc mapping from mapping.tsv")

  errorMessages <- character()
  expect_error(
    readProtImportOptionalMapping(
      mappingPath = "broken.tsv",
      mappingLabel = "UniProt mapping",
      readMapping = function(...) stop("bad file"),
      logInfo = function(...) NULL,
      logError = function(text) {
        errorMessages <<- c(errorMessages, text)
      }
    ),
    "Failed to read UniProt mapping file: bad file"
  )
  expect_identical(errorMessages, "Failed to read UniProt mapping file: bad file")
})

test_that("mod_prot_import_server stores successful import side effects", {
  fixture <- makeImportFixture(withConfig = TRUE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  workflowData <- makeWorkflowData()
  configStub <- getDefaultProteomicsConfig()

  withCleanGlobalObjects(c("aa_seq_tbl_final", "config_list"), {
    serverUnderTest <- makeFunctionWithOverrides(
      mod_prot_import_server,
      list(
        requireNamespace = function(package, quietly = FALSE) {
          if (identical(package, "shinyFiles")) {
            FALSE
          } else {
            base::requireNamespace(package, quietly = quietly)
          }
        },
        importDIANNData = function(filepath, use_precursor_norm = TRUE) {
          expect_equal(filepath, fixture$searchResultsPath)
          expect_true(isTRUE(use_precursor_norm))
          mockImportResult()
        },
        processFastaFile = function(...) {
          mockFastaResult()
        },
        readConfigFile = function(file) {
          expect_equal(file, fixture$configPath)
          configStub
        },
        .capture_checkpoint = function(...) invisible(NULL)
      )
    )

    testServer(
      serverUnderTest,
      args = list(
        workflow_data = workflowData,
        experiment_paths = fixture$experimentPaths,
        volumes = NULL
      ),
      {
        session$setInputs(
          search_results_standard = makeFileInput(fixture$searchResultsPath),
          fasta_file_standard = makeFileInput(fixture$fastaPath),
          config_file_standard = makeFileInput(fixture$configPath),
          format_override = "diann",
          diann_use_precursor_norm = TRUE,
          sanitize_names = TRUE,
          mixed_species_fasta = FALSE,
          taxon_id = 9606,
          organism_name = "Homo sapiens",
          capture_checkpoints = FALSE
        )
        session$flushReact()

        session$setInputs(process_data = 1)
        session$flushReact()

        expect_equal(workflowData$data_format, "diann")
        expect_equal(workflowData$data_type, "peptide")
        expect_equal(workflowData$state_manager$workflow_type, "DIA")
        expect_equal(workflowData$tab_status$setup_import, "complete")
        expect_equal(workflowData$fasta_file_path, fixture$fastaPath)
        expect_equal(sort(unique(workflowData$data_tbl$Run)), c("sample_1", "sample_2"))
        expect_identical(workflowData$data_cln, workflowData$data_tbl)
        expect_equal(workflowData$taxon_id, 9606)
        expect_equal(workflowData$organism_name, "Homo sapiens")
        expect_equal(workflowData$config_list, configStub)
        expect_equal(workflowData$fasta_metadata$num_sequences, 2)
        expect_false(workflowData$mixed_species_analysis$enabled)
        expect_equal(workflowData$processing_log$setup_import$detected_format, "diann")
        expect_equal(workflowData$processing_log$setup_import$n_runs, 2)
        expect_true(exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE))
        expect_true(exists("config_list", envir = .GlobalEnv, inherits = FALSE))
      }
    )
  })
})

test_that("mod_prot_import_server falls back to default config when optional files are absent", {
  fixture <- makeImportFixture(withConfig = FALSE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  workflowData <- makeWorkflowData()
  fallbackConfig <- list(
    generalParameters = list(min_peptides_per_protein = 3),
    deAnalysisParameters = list(formula_string = "~ 0 + group"),
    normalizationParameters = list(normalisation_method = "median"),
    ruvParameters = list(percentage_as_neg_ctrl = 25)
  )

  withCleanGlobalObjects(c("aa_seq_tbl_final", "config_list"), {
    serverUnderTest <- makeFunctionWithOverrides(
      mod_prot_import_server,
      list(
        requireNamespace = function(package, quietly = FALSE) {
          if (identical(package, "shinyFiles")) {
            FALSE
          } else {
            base::requireNamespace(package, quietly = quietly)
          }
        },
        importDIANNData = function(...) {
          mockImportResult()
        },
        processFastaFile = function(...) {
          mockFastaResult()
        },
        system.file = function(...) "",
        download.file = function(...) {
          stop("offline test fallback")
        },
        getDefaultProteomicsConfig = function() fallbackConfig,
        .capture_checkpoint = function(...) invisible(NULL)
      )
    )

    testServer(
      serverUnderTest,
      args = list(
        workflow_data = workflowData,
        experiment_paths = fixture$experimentPaths,
        volumes = NULL
      ),
      {
        session$setInputs(
          search_results_standard = makeFileInput(fixture$searchResultsPath),
          fasta_file_standard = makeFileInput(fixture$fastaPath),
          format_override = "diann",
          diann_use_precursor_norm = TRUE,
          sanitize_names = FALSE,
          mixed_species_fasta = FALSE,
          taxon_id = 9606,
          organism_name = "Homo sapiens",
          capture_checkpoints = FALSE
        )
        session$flushReact()

        session$setInputs(process_data = 1)
        session$flushReact()

        expect_equal(workflowData$tab_status$setup_import, "complete")
        expect_equal(workflowData$config_list, fallbackConfig)
        expect_null(workflowData$uniprot_mapping)
        expect_null(workflowData$uniparc_mapping)
        expect_equal(workflowData$processing_log$setup_import$organism, "Homo sapiens")
      }
    )
  })
})

test_that("mod_prot_import_server cleans workflow state after unsupported format failure", {
  fixture <- makeImportFixture(withConfig = FALSE)
  on.exit(unlink(fixture$rootDir, recursive = TRUE, force = TRUE), add = TRUE)

  workflowData <- makeWorkflowData()
  workflowData$data_tbl <- data.frame(old = 1)
  workflowData$data_format <- "stale"
  workflowData$data_type <- "stale"
  workflowData$column_mapping <- list(old = TRUE)
  workflowData$data_cln <- data.frame(old = 1)
  workflowData$fasta_file_path <- "stale.fasta"
  workflowData$aa_seq_tbl_final <- data.frame(old = 1)
  workflowData$config_list <- list(old = TRUE)
  workflowData$processing_log <- list(setup_import = list(old = TRUE))
  workflowData$tab_status <- list(
    setup_import = "complete",
    design_matrix = "disabled",
    quality_control = "disabled",
    normalization = "disabled",
    differential_expression = "disabled",
    enrichment_analysis = "disabled",
    session_summary = "disabled"
  )

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_import_server,
    list(
      requireNamespace = function(package, quietly = FALSE) {
        if (identical(package, "shinyFiles")) {
          FALSE
        } else {
          base::requireNamespace(package, quietly = quietly)
        }
      },
      .capture_checkpoint = function(...) invisible(NULL)
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      workflow_data = workflowData,
      experiment_paths = fixture$experimentPaths,
      volumes = NULL
    ),
    {
      session$setInputs(
        search_results_standard = makeFileInput(fixture$searchResultsPath),
        fasta_file_standard = makeFileInput(fixture$fastaPath),
        format_override = "unknown",
        mixed_species_fasta = FALSE,
        capture_checkpoints = FALSE
      )
      session$flushReact()

      session$setInputs(process_data = 1)
      session$flushReact()

      expect_null(workflowData$data_tbl)
      expect_null(workflowData$data_format)
      expect_null(workflowData$data_type)
      expect_null(workflowData$column_mapping)
      expect_null(workflowData$data_cln)
      expect_null(workflowData$fasta_file_path)
      expect_null(workflowData$aa_seq_tbl_final)
      expect_null(workflowData$config_list)
      expect_null(workflowData$processing_log$setup_import)
      expect_equal(workflowData$tab_status$setup_import, "incomplete")
    }
  )
})

# APAF Bioinformatics | test-prot-01c-import-module-contracts.R | Approved | 2026-04-11
