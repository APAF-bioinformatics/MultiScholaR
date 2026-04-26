# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

makeProtImportCharacterizationOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

withProtImportCharacterizationGlobals <- function(object_names, code) {
  had_existing <- vapply(
    object_names,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
  old_values <- lapply(seq_along(object_names), function(i) {
    if (had_existing[[i]]) {
      get(object_names[[i]], envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(old_values) <- object_names

  on.exit({
    for (name in rev(object_names)) {
      if (had_existing[[name]]) {
        assign(name, old_values[[name]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

makeProtImportCharacterizationFileInput <- function(path, type = "text/plain") {
  data.frame(
    name = basename(path),
    size = unname(file.info(path)$size),
    type = type,
    datapath = path,
    stringsAsFactors = FALSE
  )
}

makeProtImportCharacterizationFixture <- function(with_config = TRUE) {
  root_dir <- tempfile("prot-import-characterization-")
  dir.create(root_dir, recursive = TRUE)

  search_results_path <- file.path(root_dir, "report.tsv")
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
    search_results_path
  )

  fasta_path <- file.path(root_dir, "proteins.fasta")
  writeLines(
    c(
      ">sp|P1|PROT1 OS=Homo sapiens OX=9606",
      "MPEPTIDESEQ"
    ),
    fasta_path
  )

  config_path <- NULL
  if (with_config) {
    config_path <- file.path(root_dir, "config.ini")
    writeLines("[generalParameters]\nmin_peptides_per_protein=2", config_path)
  }

  results_dir <- file.path(root_dir, "results")
  source_dir <- file.path(root_dir, "scripts")
  dir.create(results_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  list(
    root_dir = root_dir,
    search_results_path = search_results_path,
    fasta_path = fasta_path,
    config_path = config_path,
    experiment_paths = list(
      results_dir = results_dir,
      source_dir = source_dir
    )
  )
}

makeProtImportCharacterizationWorkflow <- function() {
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

makeProtImportCharacterizationWorkflowEnv <- function() {
  state_manager <- new.env(parent = emptyenv())
  state_manager$workflow_type <- NULL
  state_manager$setWorkflowType <- function(workflow_type) {
    state_manager$workflow_type <- workflow_type
    invisible(workflow_type)
  }

  list2env(list(
    data_tbl = NULL,
    fasta_file_path = NULL,
    config_list = NULL,
    taxon_id = NULL,
    organism_name = NULL,
    design_matrix = NULL,
    data_cln = NULL,
    contrasts_tbl = NULL,
    state_manager = state_manager,
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
  ), parent = emptyenv())
}

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

skipIfMissingMultiScholaRBindings <- function(...) {
  names <- unlist(list(...), use.names = FALSE)
  missing <- names[!vapply(names, hasMultiScholaRBinding, logical(1))]
  if (length(missing) > 0) {
    skip(paste("requires extracted helper bindings:", paste(missing, collapse = ", ")))
  }
}

mockProtImportCharacterizationResult <- function() {
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

mockProtImportCharacterizationFasta <- function() {
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

test_that("prot import shared helpers preserve format detection and filename resolution behavior", {
  skipIfMissingMultiScholaRBindings(
    "readProtImportHeaders",
    "resolveProtImportDetectionFilename",
    "runProtImportFormatDetection",
    "applyProtImportDetectedFormat",
    "resetProtImportFormatDetectionState"
  )

  fixture <- makeProtImportCharacterizationFixture(with_config = FALSE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_identical(
    readProtImportHeaders(fixture$search_results_path, readExcel = function(...) stop("unexpected xlsx path")),
    c(
      "Protein.Group", "Protein.Ids", "Protein.Names", "Precursor.Id",
      "Modified.Sequence", "Stripped.Sequence", "Precursor.Charge",
      "Q.Value", "PG.Q.Value", "Run"
    )
  )

  local_state <- reactiveValues(detected_format = NULL, format_confidence = NULL)
  detected <- runProtImportFormatDetection(
    filePath = fixture$search_results_path,
    useShinyFiles = FALSE,
    localData = local_state,
    searchResultsStandard = list(name = "report.parquet"),
    readHeaders = function(path, ...) {
      expect_identical(path, fixture$search_results_path)
      c(
        "Protein.Group", "Protein.Ids", "Protein.Names",
        "Precursor.Id", "Modified.Sequence", "Stripped.Sequence",
        "Precursor.Charge", "Q.Value", "PG.Q.Value", "Run"
      )
    },
    detectFormat = function(headers, filename) {
      expect_identical(filename, "report.parquet")
      detectProteomicsFormat(headers, filename)
    },
    logInfo = function(...) invisible(NULL),
    logError = function(...) invisible(NULL)
  )

  applied <- applyProtImportDetectedFormat(
    localData = local_state,
    formatInfo = detected,
    logInfo = function(...) invisible(NULL)
  )

  expect_identical(applied$format, "diann")
  expect_identical(shiny::isolate(local_state$detected_format), "diann")
  expect_gt(shiny::isolate(local_state$format_confidence), 0.8)
  expect_identical(
    resolveProtImportDetectionFilename(
      useShinyFiles = FALSE,
      filePath = fixture$search_results_path,
      searchResultsStandard = list(name = "report.tsv")
    ),
    "report.tsv"
  )
  expect_identical(
    resolveProtImportDetectionFilename(
      useShinyFiles = TRUE,
      filePath = "/tmp/shiny/report.tsv",
      searchResultsStandard = NULL
    ),
    "report.tsv"
  )
  expect_null(
    resolveProtImportDetectionFilename(
      useShinyFiles = FALSE,
      filePath = fixture$search_results_path,
      searchResultsStandard = NULL
    )
  )

  reset <- resetProtImportFormatDetectionState(
    localData = local_state,
    errorMessage = "boom",
    logError = function(...) invisible(NULL)
  )

  expect_identical(reset$format, "unknown")
  expect_identical(shiny::isolate(local_state$detected_format), "unknown")
  expect_identical(shiny::isolate(local_state$format_confidence), 0)
})

test_that("prot import shared helpers preserve shiny path and optional upload resolution", {
  skipIfMissingMultiScholaRBindings(
    "handleProtImportShinyFileSelection",
    "resolveProtImportInputFilename",
    "resolveProtImportPrimaryUploadPath",
    "resolveProtImportOptionalUploadPath"
  )

  parse_calls <- list()
  selected_path <- handleProtImportShinyFileSelection(
    selectedInput = "token",
    parseFilePaths = function(volumes, selection, ...) {
      parse_calls[[length(parse_calls) + 1]] <<- list(volumes = volumes, selection = selection)
      data.frame(datapath = "/tmp/import/report.tsv", stringsAsFactors = FALSE)
    },
    volumes = c(Home = "/tmp"),
    localData = reactiveValues(),
    localField = "searchResultsPath",
    output = new.env(parent = emptyenv()),
    outputId = "selected_path",
    messageFn = function(...) invisible(NULL)
  )

  expect_identical(selected_path, "/tmp/import/report.tsv")
  expect_identical(parse_calls[[1]]$selection, "token")
  expect_identical(
    resolveProtImportPrimaryUploadPath(
      useShinyFiles = FALSE,
      shinyPath = "/ignored",
      standardInput = list(datapath = "/tmp/report.tsv")
    ),
    "/tmp/report.tsv"
  )
  expect_identical(
    resolveProtImportOptionalUploadPath(
      useShinyFiles = FALSE,
      shinyPath = NULL,
      standardInput = list(datapath = "/tmp/config.ini")
    ),
    "/tmp/config.ini"
  )
  expect_identical(
    resolveProtImportInputFilename(
      useShinyFiles = FALSE,
      resolvedPath = "/tmp/report.tsv",
      standardInput = list(name = "report.tsv")
    ),
    "report.tsv"
  )
  expect_identical(
    resolveProtImportInputFilename(
      useShinyFiles = TRUE,
      resolvedPath = "/tmp/import/report.tsv",
      standardInput = NULL
    ),
    "report.tsv"
  )
  expect_null(
    resolveProtImportInputFilename(
      useShinyFiles = FALSE,
      resolvedPath = "/tmp/report.tsv",
      standardInput = NULL
    )
  )
})

test_that("prot import shared helpers preserve configuration loading behavior", {
  skipIfMissingMultiScholaRBindings(
    "loadProtImportConfiguration",
    "readProtImportOptionalMapping",
    "loadProtImportOptionalMappings",
    "storeProtImportConfiguration",
    "loadProtImportConfigurationResources"
  )

  fixture <- makeProtImportCharacterizationFixture(with_config = TRUE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  uniprot_path <- file.path(fixture$root_dir, "uniprot.tsv")
  uniparc_path <- file.path(fixture$root_dir, "uniparc.tsv")
  writeLines("id\tvalue\nP1\tA", uniprot_path)
  writeLines("id\tvalue\nUPI1\tB", uniparc_path)

  config <- loadProtImportConfiguration(
    configPath = fixture$config_path,
    experimentPaths = fixture$experiment_paths,
    readConfig = function(file) {
      expect_identical(file, fixture$config_path)
      list(config_source = "explicit")
    },
    showNotification = function(...) invisible(NULL),
    logInfo = function(...) invisible(NULL),
    logError = function(...) stop("unexpected error"),
    getDefaultConfig = function() stop("unexpected fallback")
  )

  expect_identical(config$config_source, "explicit")

  workflow_data <- list2env(list())
  loaded <- loadProtImportConfigurationResources(
    workflowData = workflow_data,
    experimentPaths = fixture$experiment_paths,
    useShinyFiles = FALSE,
    configFilePath = fixture$config_path,
    configFileStandard = list(datapath = fixture$config_path),
    uniprotMappingFile = uniprot_path,
    uniprotMappingStandard = list(datapath = uniprot_path),
    uniparcMappingFile = uniparc_path,
    uniparcMappingStandard = list(datapath = uniparc_path),
    updateStatus = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL),
    resolveUploadPath = function(useShinyFiles, shinyPath, standardInput) standardInput$datapath,
    loadConfiguration = function(configPath, experimentPaths, ...) {
      list(path = configPath, sourceDir = experimentPaths$source_dir)
    },
    storeConfiguration = storeProtImportConfiguration,
    loadOptionalMappings = loadProtImportOptionalMappings
  )

  expect_identical(loaded$configPath, fixture$config_path)
  expect_true(is.data.frame(workflow_data$uniprot_mapping))
  expect_true(is.data.frame(workflow_data$uniparc_mapping))
  expect_identical(workflow_data$config_list$path, fixture$config_path)
})

test_that("prot import shared helpers preserve import application and setup logging", {
  skipIfMissingMultiScholaRBindings(
    "importProtImportDataByFormat",
    "applyProtImportResultToWorkflow",
    "finalizeProtImportSetupState",
    "recordProtImportSetupLog",
    "loadProtImportConfiguration"
  )

  workflow_data <- makeProtImportCharacterizationWorkflowEnv()
  result <- importProtImportDataByFormat(
    format = "diann",
    searchResultsPath = "/tmp/report.tsv",
    input = list(diann_use_precursor_norm = TRUE),
    importDiann = function(filepath, use_precursor_norm = TRUE) {
      expect_identical(filepath, "/tmp/report.tsv")
      expect_true(isTRUE(use_precursor_norm))
      mockProtImportCharacterizationResult()
    }
  )

  workflow_type <- applyProtImportResultToWorkflow(
    workflowData = workflow_data,
    dataImportResult = result,
    format = "diann",
    fastaPath = "/tmp/proteins.fasta",
    sanitizeNames = TRUE,
    sanitizeRunNames = function(workflowData, runCol, ...) {
      shiny::isolate({
        workflowData$data_tbl[[runCol]] <- c("sample_1", "sample_2")
      })
      TRUE
    }
  )

  expect_identical(workflow_type, "DIA")
  expect_identical(
    sort(unique(shiny::isolate(workflow_data$data_tbl$Run))),
    c("sample_1", "sample_2")
  )

  setup_log <- finalizeProtImportSetupState(
    workflowData = workflow_data,
    dataImportResult = result,
    format = "diann",
    searchFilename = "report.tsv",
    fastaFilename = "proteins.fasta",
    taxonId = 9606,
    organismName = "Homo sapiens",
    mixedSpeciesFasta = FALSE
  )

  expect_identical(setup_log$detected_format, "diann")
  expect_identical(setup_log$organism, "Homo sapiens")

  recorded <- recordProtImportSetupLog(
    workflowData = workflow_data,
    dataImportResult = result,
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
      if (!is.null(standardInput$name)) standardInput$name else basename(resolvedPath)
    },
    finalizeSetupState = finalizeProtImportSetupState,
    logInfo = function(...) invisible(NULL)
  )

  expect_identical(recorded$searchFilename, "report.tsv")
  expect_identical(recorded$fastaFilename, "proteins.fasta")
})

test_that("prot import shared helpers preserve mixed-species analysis and selection behavior", {
  skipIfMissingMultiScholaRBindings(
    "toggleProtImportMixedSpeciesInputs",
    "analyzeProtImportMixedSpeciesData",
    "runProtImportMixedSpeciesAnalysisIfNeeded",
    "confirmProtImportOrganismSelection",
    "buildProtImportOrganismSelectionModal"
  )

  selection <- toggleProtImportMixedSpeciesInputs(
    mixedSpeciesFasta = TRUE,
    disableInput = function(...) invisible(NULL),
    enableInput = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL)
  )
  expect_true(isTRUE(selection))

  workflow_data <- makeProtImportCharacterizationWorkflowEnv()
  workflow_data$aa_seq_tbl_final <- data.frame(
    Protein.Group = c("P1", "P2"),
    Organism = c("Human", "Mouse"),
    Taxon = c(9606L, 10090L),
    stringsAsFactors = FALSE
  )
  workflow_data$data_tbl <- data.frame(
    Protein.Group = c("P1", "P2"),
    value = c(1, 2),
    stringsAsFactors = FALSE
  )
  workflow_data$data_cln <- workflow_data$data_tbl
  workflow_data$column_mapping <- list(protein_col = "Protein.Group")

  local_data <- list2env(list(
    waiting_for_organism_selection = FALSE,
    organism_distribution = NULL,
    organism_mapping = NULL
  ), parent = emptyenv())
  modal_calls <- list()
  analysis <- analyzeProtImportMixedSpeciesData(
    workflowData = workflow_data,
    localData = local_data,
    dataImportResult = mockProtImportCharacterizationResult(),
    fastaPath = "/tmp/proteins.fasta",
    session = shiny::MockShinySession$new(),
    updateStatus = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL),
    logInfo = function(...) invisible(NULL),
    logWarn = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    extractOrganisms = function(...) {
      data.frame(
        uniprot_acc = c("P1", "P2"),
        taxon_id = c(9606L, 10090L),
        organism_name = c("Human", "Mouse"),
        stringsAsFactors = FALSE
      )
    },
    analyzeDistribution = function(...) {
      data.frame(
        taxon_id = c(9606L, 10090L),
        organism_name = c("Human", "Mouse"),
        protein_count = c(1L, 1L),
        percentage = c(50, 50),
        stringsAsFactors = FALSE
      )
    },
    showModal = function(modal) {
      modal_calls[[length(modal_calls) + 1]] <<- modal
      invisible(NULL)
    },
    buildSelectionModal = buildProtImportOrganismSelectionModal
  )

  expect_true(isTRUE(analysis))
  expect_true(isTRUE(local_data$waiting_for_organism_selection))
  expect_length(modal_calls, 1)

  rerun <- runProtImportMixedSpeciesAnalysisIfNeeded(
    workflowData = workflow_data,
    localData = local_data,
    dataImportResult = mockProtImportCharacterizationResult(),
    fastaPath = "/tmp/proteins.fasta",
    session = shiny::MockShinySession$new(),
    mixedSpeciesFasta = TRUE,
    analyzeMixedSpeciesData = function(...) "ran"
  )
  expect_identical(rerun, "ran")

  confirmed <- confirmProtImportOrganismSelection(
    selectedTaxon = 9606L,
    filterToOrganism = TRUE,
    workflowData = workflow_data,
    localData = local_data,
    session = shiny::MockShinySession$new(),
    updateNumericInput = function(...) invisible(NULL),
    updateTextInput = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    removeModal = function() invisible(NULL),
    now = function() as.POSIXct("2026-04-20 20:00:00", tz = "UTC"),
    normalizeUniprot = function(x, remove_isoform = TRUE) x,
    cleanAcc = function(x) x
  )

  expect_true(isTRUE(confirmed))
  expect_identical(workflow_data$taxon_id, 9606L)
  expect_true(all(grepl("P1", shiny::isolate(workflow_data$data_tbl$Protein.Group))))
})

test_that("prot import shared helpers preserve orchestration lifecycle behavior", {
  skipIfMissingMultiScholaRBindings(
    "startProtImportProcessing",
    "applyProtImportWorkflowWithStatus",
    "processProtImportFastaData",
    "completeProtImportSuccessState",
    "resetProtImportWorkflowStateOnError"
  )

  fixture <- makeProtImportCharacterizationFixture(with_config = FALSE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workflow_data <- makeProtImportCharacterizationWorkflowEnv()
  local_data <- list2env(
    list(processing = FALSE, waiting_for_organism_selection = FALSE),
    parent = emptyenv()
  )

  started <- startProtImportProcessing(
    searchResultsPath = fixture$search_results_path,
    fastaPath = fixture$fasta_path,
    format = "diann",
    localData = local_data,
    ns = function(id) paste0("ns-", id),
    requireValue = function(value) value,
    announceStart = function(...) invisible(NULL),
    buildProcessingModal = function(ns) list(ns = ns("processing_status")),
    showModal = function(...) invisible(NULL)
  )

  expect_true(isTRUE(local_data$processing))
  expect_identical(started$format, "diann")

  applied <- applyProtImportWorkflowWithStatus(
    workflowData = workflow_data,
    dataImportResult = mockProtImportCharacterizationResult(),
    format = "diann",
    fastaPath = fixture$fasta_path,
    organismName = "Homo sapiens",
    experimentPaths = fixture$experiment_paths,
    sanitizeNames = TRUE,
    updateStatus = function(...) invisible(NULL),
    applyImportResult = applyProtImportResultToWorkflow,
    processFastaData = processProtImportFastaData,
    logInfo = function(...) invisible(NULL),
    logWarn = function(...) invisible(NULL),
    showNotification = function(...) invisible(NULL),
    sanitizeRunNames = function(workflowData, runCol, ...) {
      shiny::isolate({
        workflowData$data_tbl[[runCol]] <- c("sample_1", "sample_2")
      })
      TRUE
    },
    processFasta = function(...) mockProtImportCharacterizationFasta()
  )

  expect_identical(applied$sanitizeNames, TRUE)
  expect_true(is.data.frame(workflow_data$aa_seq_tbl_final))

  completed <- completeProtImportSuccessState(
    workflowData = workflow_data,
    localData = local_data,
    dataImportResult = mockProtImportCharacterizationResult(),
    format = "diann",
    removeModal = function() invisible(NULL),
    showNotification = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL)
  )

  expect_identical(completed$setup_import, "complete")

  reset <- resetProtImportWorkflowStateOnError(
    workflowData = workflow_data,
    localData = local_data,
    errorMessage = "boom",
    logError = function(...) invisible(NULL),
    removeModal = function() invisible(NULL),
    showNotification = function(...) invisible(NULL),
    messageFn = function(...) invisible(NULL)
  )

  expect_identical(reset$setup_import, "incomplete")
  expect_null(workflow_data$data_tbl)
})

test_that("mod_prot_import_server preserves successful import side effects", {
  fixture <- makeProtImportCharacterizationFixture(with_config = TRUE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workflow_data <- makeProtImportCharacterizationWorkflow()
  config_stub <- getDefaultProteomicsConfig()

  withProtImportCharacterizationGlobals(c("aa_seq_tbl_final", "config_list"), {
    server_under_test <- makeProtImportCharacterizationOverrides(
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
          expect_equal(filepath, fixture$search_results_path)
          expect_true(isTRUE(use_precursor_norm))
          mockProtImportCharacterizationResult()
        },
        processFastaFile = function(...) {
          mockProtImportCharacterizationFasta()
        },
        readConfigFile = function(file) {
          expect_equal(file, fixture$config_path)
          config_stub
        },
        .capture_checkpoint = function(...) invisible(NULL)
      )
    )

    testServer(
      server_under_test,
      args = list(
        workflow_data = workflow_data,
        experiment_paths = fixture$experiment_paths,
        volumes = NULL
      ),
      {
        session$setInputs(
          search_results_standard = makeProtImportCharacterizationFileInput(fixture$search_results_path),
          fasta_file_standard = makeProtImportCharacterizationFileInput(fixture$fasta_path),
          config_file_standard = makeProtImportCharacterizationFileInput(fixture$config_path),
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

        expect_equal(workflow_data$data_format, "diann")
        expect_equal(workflow_data$data_type, "peptide")
        expect_equal(workflow_data$state_manager$workflow_type, "DIA")
        expect_equal(workflow_data$tab_status$setup_import, "complete")
        expect_equal(workflow_data$fasta_file_path, fixture$fasta_path)
        expect_equal(sort(unique(workflow_data$data_tbl$Run)), c("sample_1", "sample_2"))
        expect_identical(workflow_data$data_cln, workflow_data$data_tbl)
        expect_equal(workflow_data$taxon_id, 9606)
        expect_equal(workflow_data$organism_name, "Homo sapiens")
        expect_equal(workflow_data$config_list, config_stub)
        expect_equal(workflow_data$fasta_metadata$num_sequences, 2)
        expect_false(workflow_data$mixed_species_analysis$enabled)
        expect_equal(workflow_data$processing_log$setup_import$detected_format, "diann")
        expect_equal(workflow_data$processing_log$setup_import$n_runs, 2)
        expect_true(exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE))
        expect_true(exists("config_list", envir = .GlobalEnv, inherits = FALSE))
      }
    )
  })
})

test_that("mod_prot_import_server preserves unsupported-format cleanup behavior", {
  fixture <- makeProtImportCharacterizationFixture(with_config = FALSE)
  on.exit(unlink(fixture$root_dir, recursive = TRUE, force = TRUE), add = TRUE)

  workflow_data <- makeProtImportCharacterizationWorkflow()
  workflow_data$data_tbl <- data.frame(old = 1)
  workflow_data$data_format <- "stale"
  workflow_data$data_type <- "stale"
  workflow_data$column_mapping <- list(old = TRUE)
  workflow_data$data_cln <- data.frame(old = 1)
  workflow_data$fasta_file_path <- "stale.fasta"
  workflow_data$aa_seq_tbl_final <- data.frame(old = 1)
  workflow_data$config_list <- list(old = TRUE)
  workflow_data$processing_log <- list(setup_import = list(old = TRUE))
  workflow_data$tab_status <- list(
    setup_import = "complete",
    design_matrix = "disabled",
    quality_control = "disabled",
    normalization = "disabled",
    differential_expression = "disabled",
    enrichment_analysis = "disabled",
    session_summary = "disabled"
  )

  server_under_test <- makeProtImportCharacterizationOverrides(
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
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = fixture$experiment_paths,
      volumes = NULL
    ),
    {
      session$setInputs(
        search_results_standard = makeProtImportCharacterizationFileInput(fixture$search_results_path),
        fasta_file_standard = makeProtImportCharacterizationFileInput(fixture$fasta_path),
        format_override = "unknown",
        mixed_species_fasta = FALSE,
        capture_checkpoints = FALSE
      )
      session$flushReact()

      session$setInputs(process_data = 1)
      session$flushReact()

      expect_null(workflow_data$data_tbl)
      expect_null(workflow_data$data_format)
      expect_null(workflow_data$data_type)
      expect_null(workflow_data$column_mapping)
      expect_null(workflow_data$data_cln)
      expect_null(workflow_data$fasta_file_path)
      expect_null(workflow_data$aa_seq_tbl_final)
      expect_null(workflow_data$config_list)
      expect_null(workflow_data$processing_log$setup_import)
      expect_equal(workflow_data$tab_status$setup_import, "incomplete")
    }
  )
})
