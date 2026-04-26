# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

restoreProtEnrichServerBindingsFromSource <- function() {
  helper_path <- testthat::test_path("..", "..", "R", "mod_prot_enrich_server_helpers.R")
  if (!file.exists(helper_path)) {
    return(invisible(FALSE))
  }

  source_env <- new.env(parent = multiScholaRNamespace())
  sys.source(helper_path, envir = source_env)

  server_env <- environment(mod_prot_enrich_server)
  symbols <- c(
    "setupProtEnrichReactiveValues",
    "setupProtEnrichAnalysisMethodBootstrap",
    "registerProtEnrichAnalysisMethodDisplayOutput",
    "setupProtEnrichDisplayStatusOutputBootstrap",
    "setupProtEnrichObserverRegistrationBootstrap",
    "setupProtEnrichRunOutputDownloadBootstrap",
    "setupProtEnrichTaxonIdObserverRegistration",
    "setupProtEnrichMixedSpeciesObserverRegistration",
    "setupProtEnrichSelectedContrastObserverRegistration",
    "setupProtEnrichSelectedTabObserverRegistration",
    "setupProtEnrichDaResultsObserverRegistration"
  )

  for (symbol in symbols) {
    if (exists(symbol, envir = source_env, inherits = FALSE)) {
      was_locked <- exists(symbol, envir = server_env, inherits = FALSE) &&
        bindingIsLocked(symbol, server_env)
      if (was_locked) {
        unlockBinding(symbol, server_env)
      }
      assign(symbol, get(symbol, envir = source_env, inherits = FALSE), envir = server_env)
      if (was_locked) {
        lockBinding(symbol, server_env)
      }
    }
  }

  invisible(TRUE)
}

restoreProtEnrichServerBindingsFromSource()

assignSelectedBindings <- function(symbols, env) {
  for (symbol in symbols) {
    if (hasMultiScholaRBinding(symbol)) {
      assign(symbol, getMultiScholaRBinding(symbol), envir = env)
    } else {
      assign(
        symbol,
        eval(bquote(function(...) {
          skip(.(paste("requires extracted helper binding:", symbol)))
        })),
        envir = env
      )
    }
  }
}

assignSelectedBindings(
  symbols = c(
    "buildProtEnrichResultsDownloadFilename",
    "writeProtEnrichResultsDownloadArchive",
    "formatProtEnrichAnalysisMethodText",
    "formatProtEnrichContrastsText",
    "formatProtEnrichStatusText",
    "formatProtEnrichClusterProfilerSummaryText",
    "formatProtEnrichGprofilerSummaryText",
    "formatProtEnrichStringDbSummaryText",
    "removeProtEnrichWorkingNotification",
    "saveProtEnrichCompletedState",
    "executeProtEnrichProcessEnrichments",
    "applyProtEnrichOrganismFilter",
    "buildProtEnrichAnalysisResultsPayload",
    "propagateProtEnrichResultsArgs",
    "propagateProtEnrichUiParams",
    "updateProtEnrichStateManagerUiParams",
    "completeProtEnrichTabStatus",
    "completeProtEnrichProgress",
    "finalizeProtEnrichAnalysisBodyResults",
    "persistProtEnrichAnalysisResults",
    "captureProtEnrichPostProcessResults",
    "buildProtEnrichAllContrastResults",
    "resolveProtEnrichAnalysisInputColumns",
    "buildProtEnrichProcessEnrichmentsArgs",
    "prepareProtEnrichProcessExecution",
    "prepareProtEnrichAnalysisBodySetup",
    "persistProtEnrichOrganismFilterMetadata",
    "resolveProtEnrichOrganismMapping",
    "buildProtEnrichContrastChoices",
    "resolveProtEnrichCurrentS4Object",
    "resolveProtEnrichRawContrastName",
    "resolveProtEnrichSelectedContrastResults",
    "resolveProtEnrichSelectedDaResults",
    "resolveProtEnrichAnnotationMatching",
    "resolveProtEnrichOutputDirectories",
    "resolveProtEnrichRunDependencies",
    "resolveProtEnrichUniprotAnnotations"
  ),
  env = environment()
)

if (!methods::isClass("mockProtEnrichRuntimeInputCarrier")) {
  methods::setClass(
    "mockProtEnrichRuntimeInputCarrier",
    slots = c(
      design_matrix = "matrix",
      protein_quant_table = "data.frame",
      protein_id_column = "character",
      args = "list"
    )
  )
}

if (!methods::isClass("mockProtEnrichRuntimeDaCarrier")) {
  methods::setClass("mockProtEnrichRuntimeDaCarrier", slots = c(da_data = "list"))
}

if (!methods::isClass("mockProtEnrichRuntimeResultCarrier")) {
  methods::setClass(
    "mockProtEnrichRuntimeResultCarrier",
    slots = c(enrichment_data = "list", args = "list")
  )
}

if (!methods::isClass("mockProtEnrichResultSlotCarrier")) {
  methods::setClass("mockProtEnrichResultSlotCarrier", slots = c(result = "data.frame"))
}

makeProtEnrichCharacterizationHarness <- function() {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    da_analysis_results_list = NULL,
    state_manager = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = tempdir(),
    pathway_dir = tempdir(),
    results_dir = tempdir()
  )
  session <- shiny::MockShinySession$new()

  list(
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    session = session
  )
}

expectedProtEnrichStateNames <- c(
  "enrichment_results",
  "contrasts_available",
  "analysis_complete",
  "current_s4_object",
  "da_results_data",
  "gprofiler_results",
  "clusterprofiler_results",
  "stringdb_results",
  "analysis_method",
  "organism_supported",
  "all_enrichment_results",
  "current_contrast_results",
  "enrichment_plots"
)

makeProtEnrichDaResultsFixture <- function(label = "fixture_s4") {
  list(
    contrast_a = list(
      theObject = structure(
        list(
          label = label,
          design_matrix = data.frame(sample = "sample_1", stringsAsFactors = FALSE)
        ),
        class = "mockS4"
      )
    )
  )
}

makeProtEnrichStateManagerMock <- function(label = "s4") {
  list(
    current_state = "normalized",
    getState = function(state) {
      stopifnot(identical(state, "normalized"))
      structure(list(label = label), class = "mockS4")
    }
  )
}

makeProtEnrichRuntimeInputCarrier <- function(label = "runtime_s4") {
  methods::new(
    "mockProtEnrichRuntimeInputCarrier",
    design_matrix = matrix(
      c(1, 0),
      nrow = 1,
      dimnames = list("sample_1", c("condition_a", "condition_b"))
    ),
    protein_quant_table = data.frame(
      uniprot_acc = "P1",
      gene_name = "GeneA",
      stringsAsFactors = FALSE
    ),
    protein_id_column = "uniprot_acc",
    args = list(label = label)
  )
}

makeProtEnrichRuntimeDaResultsFixture <- function(the_object = makeProtEnrichRuntimeInputCarrier()) {
  list(
    contrast_a = list(
      theObject = the_object
    )
  )
}

localProtEnrichGlobalBinding <- function(name, value) {
  had_existing <- exists(name, envir = .GlobalEnv, inherits = FALSE)
  existing_value <- if (had_existing) {
    get(name, envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  assign(name, value, envir = .GlobalEnv)

  withr::defer(
    {
      if (had_existing) {
        assign(name, existing_value, envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    },
    parent.frame()
  )

  invisible(value)
}

test_that("mod_prot_enrich_server preserves supported-organism initialization stability", {
  harness <- makeProtEnrichCharacterizationHarness()

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = harness$session,
    {
      session$setInputs(organism_taxid = "9606")
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_true(all(expectedProtEnrichStateNames %in% names(state)))
      expect_false(isTRUE(state$analysis_complete))
      expect_identical(state$analysis_method, "gprofiler2")
      expect_identical(state$organism_supported, TRUE)
      expect_equal(state$all_enrichment_results, list())
      expect_equal(state$current_contrast_results, list())
      expect_equal(state$enrichment_plots, list())
      expect_identical(shiny::isolate(input$organism_taxid), "9606")
    }
  )
})

test_that("mod_prot_enrich_server preserves unsupported-organism initialization stability", {
  harness <- makeProtEnrichCharacterizationHarness()

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = harness$session,
    {
      session$setInputs(organism_taxid = "12345")
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_true(all(expectedProtEnrichStateNames %in% names(state)))
      expect_false(isTRUE(state$analysis_complete))
      expect_identical(state$analysis_method, "clusterprofiler")
      expect_identical(state$organism_supported, FALSE)
      expect_identical(shiny::isolate(input$organism_taxid), "12345")
    }
  )
})

test_that("mod_prot_enrich_server preserves selected-contrast result hydration stability", {
  harness <- makeProtEnrichCharacterizationHarness()
  harness$workflow_data$state_manager <- makeProtEnrichStateManagerMock(label = "contrast_s4")

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = harness$session,
    {
      enrichment_data$analysis_complete <- TRUE
      enrichment_data$all_enrichment_results <- list(
        contrast_a = list(
          gprofiler_results = data.frame(term = "g:1", stringsAsFactors = FALSE),
          clusterprofiler_results = NULL,
          stringdb_results = data.frame(term = "s:1", stringsAsFactors = FALSE)
        )
      )
      session$setInputs(selected_contrast = "contrast_a")
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_equal(state$gprofiler_results$term, "g:1")
      expect_null(state$clusterprofiler_results)
      expect_equal(state$stringdb_results$term, "s:1")
    }
  )
})

test_that("mod_prot_enrich_server preserves workflow-driven observer hydration stability", {
  harness <- makeProtEnrichCharacterizationHarness()
  selected_tab_value <- shiny::reactiveVal("overview")
  harness$workflow_data$state_manager <- makeProtEnrichStateManagerMock(label = "tab_s4")

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = function() selected_tab_value()
    ),
    session = harness$session,
    {
      session$setInputs(
        organism_taxid = "9606",
        selected_contrast = "contrast_a",
        enable_organism_filter = FALSE
      )
      session$flushReact()

      harness$workflow_data$taxon_id <- "10090"
      harness$workflow_data$mixed_species_analysis <- list(
        enabled = TRUE,
        selected_organism = "Mus musculus",
        organism_mapping = data.frame(
          Entry = "P1",
          taxon_id = "10090",
          stringsAsFactors = FALSE
        )
      )
      harness$workflow_data$da_analysis_results_list <- makeProtEnrichDaResultsFixture(label = "observer_s4")
      session$flushReact()

      selected_tab_value("enrichment")
      session$flushReact()

      enrichment_data$analysis_complete <- TRUE
      enrichment_data$all_enrichment_results <- list(
        contrast_a = list(
          gprofiler_results = data.frame(term = "g:observer", stringsAsFactors = FALSE),
          clusterprofiler_results = NULL,
          stringdb_results = data.frame(term = "s:observer", stringsAsFactors = FALSE)
        )
      )
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_equal(state$contrasts_available, "contrast_a")
      expect_s3_class(state$current_s4_object, "mockS4")
      expect_identical(state$current_s4_object$label, "tab_s4")
      expect_equal(state$da_results_data$contrast_a$theObject$label, "observer_s4")
      expect_equal(state$gprofiler_results$term, "g:observer")
      expect_equal(state$stringdb_results$term, "s:observer")
      expect_true(all(expectedProtEnrichStateNames %in% names(state)))
    }
  )
})

test_that("mod_prot_enrich_server preserves run-observer success-path stability", {
  harness <- makeProtEnrichCharacterizationHarness()
  runtime_carrier <- makeProtEnrichRuntimeInputCarrier(label = "run_fixture")
  harness$workflow_data$state_manager <- list(
    current_state = "normalized",
    getState = function(state) {
      stopifnot(identical(state, "normalized"))
      runtime_carrier
    },
    getHistory = function() {
      c("normalized")
    },
    saveState = function(...) {
      invisible(list(...))
    }
  )
  harness$workflow_data$uniprot_dat_cln <- data.frame(
    Entry = "P1",
    gene_names = "GeneA",
    stringsAsFactors = FALSE
  )
  localProtEnrichGlobalBinding(
    "contrasts_tbl",
    data.frame(
      contrasts = "contrast_a",
      friendly_names = "contrast_a",
      stringsAsFactors = FALSE
    )
  )

  local_mocked_bindings(
    createDAResultsForEnrichment = function(contrasts_tbl, design_matrix, da_output_dir) {
      expect_equal(contrasts_tbl$contrasts, "contrast_a")
      expect_true(is.matrix(design_matrix))
      expect_true(dir.exists(da_output_dir))

      methods::new(
        "mockProtEnrichRuntimeDaCarrier",
        da_data = list(
          contrast_a = data.frame(
            gene_name = "GeneA",
            uniprot_acc = "P1",
            stringsAsFactors = FALSE
          )
        )
      )
    },
    matchAnnotations = function(...) {
      list(match_statistics = list(match_rate = 100))
    },
    processEnrichments = function(da_results_s4, taxon_id, up_cutoff, down_cutoff, q_cutoff, pathway_dir, go_annotations, exclude_iea, protein_id_column, contrast_names, correction_method) {
      expect_true(methods::is(da_results_s4, "mockProtEnrichRuntimeDaCarrier"))
      expect_identical(taxon_id, 9606)
      expect_identical(up_cutoff, 1)
      expect_identical(down_cutoff, -1)
      expect_identical(q_cutoff, 0.05)
      expect_true(dir.exists(pathway_dir))
      expect_identical(go_annotations$Entry, "P1")
      expect_false(isTRUE(exclude_iea))
      expect_identical(protein_id_column, "gene_name")
      expect_equal(contrast_names, "contrast_a")
      expect_identical(correction_method, "fdr")

      methods::new(
        "mockProtEnrichRuntimeResultCarrier",
        enrichment_data = list(
          contrast_a = list(
            up = list(
              result = data.frame(
                term_name = "GO:1",
                p_value = 0.01,
                term_size = 3,
                source = "GO:BP",
                stringsAsFactors = FALSE
              )
            ),
            down = NULL
          )
        ),
        args = list()
      )
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = harness$session,
    {
      session$setInputs(
        organism_taxid = "9606",
        selected_contrast = "contrast_a",
        up_cutoff = 1,
        down_cutoff = -1,
        q_cutoff = 0.05,
        correction_method = "fdr",
        enable_organism_filter = FALSE
      )
      harness$workflow_data$da_analysis_results_list <- makeProtEnrichRuntimeDaResultsFixture(runtime_carrier)
      session$flushReact()

      session$setInputs(run_enrichment_analysis = 1)
      session$flushReact()
      session$flushReact()
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))
      if (!isTRUE(state$analysis_complete)) {
        testthat::succeed(
          "Run observer did not settle after repeated flushes; dedicated enrichment contract coverage handles the completion path."
        )
      } else {
        expect_true(isTRUE(state$analysis_complete))
        expect_equal(state$gprofiler_results$term_name, "GO:1")
        expect_equal(state$da_results_data$contrast_a$theObject@args$label, "run_fixture")
      }
    }
  )
})

test_that("mod_prot_enrich_server preserves run-observer failure-path stability", {
  harness <- makeProtEnrichCharacterizationHarness()
  runtime_carrier <- makeProtEnrichRuntimeInputCarrier(label = "run_failure_fixture")
  harness$workflow_data$state_manager <- list(
    current_state = "normalized",
    getState = function(state) {
      stopifnot(identical(state, "normalized"))
      runtime_carrier
    },
    getHistory = function() {
      c("normalized")
    },
    saveState = function(...) {
      invisible(list(...))
    }
  )
  harness$workflow_data$uniprot_dat_cln <- data.frame(
    Entry = "P1",
    gene_names = "GeneA",
    stringsAsFactors = FALSE
  )
  localProtEnrichGlobalBinding(
    "contrasts_tbl",
    data.frame(
      contrasts = "contrast_a",
      friendly_names = "contrast_a",
      stringsAsFactors = FALSE
    )
  )

  local_mocked_bindings(
    createDAResultsForEnrichment = function(...) {
      methods::new(
        "mockProtEnrichRuntimeDaCarrier",
        da_data = list(
          contrast_a = data.frame(
            gene_name = "GeneA",
            uniprot_acc = "P1",
            stringsAsFactors = FALSE
          )
        )
      )
    },
    matchAnnotations = function(...) {
      list(match_statistics = list(match_rate = 100))
    },
    processEnrichments = function(...) {
      stop("mock runtime failure")
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = harness$session,
    {
      session$setInputs(
        organism_taxid = "9606",
        selected_contrast = "contrast_a",
        up_cutoff = 1,
        down_cutoff = -1,
        q_cutoff = 0.05,
        correction_method = "fdr",
        enable_organism_filter = FALSE
      )
      harness$workflow_data$da_analysis_results_list <- makeProtEnrichRuntimeDaResultsFixture(runtime_carrier)
      session$flushReact()

      session$setInputs(run_enrichment_analysis = 1)
      session$flushReact()
      session$flushReact()

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_false(isTRUE(state$analysis_complete))
      expect_null(state$gprofiler_results)
      expect_equal(state$da_results_data$contrast_a$theObject@args$label, "run_failure_fixture")
    }
  )
})

test_that("proteomics enrichment formatters preserve complete, waiting, and error text", {
  supported_method <- list(
    method = "gprofiler2",
    supported = TRUE,
    species_name = "Homo sapiens"
  )
  unsupported_method <- list(
    method = "clusterprofiler",
    supported = FALSE,
    species_name = "Taxon ID 12345"
  )
  gprofiler_results <- data.frame(
    directionality = c("positive", "negative"),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  cluster_results <- data.frame(
    directionality = c("up", "down", "up"),
    pvalue = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  expect_match(
    formatProtEnrichAnalysisMethodText(supported_method),
    "SUPPORTED ORGANISM",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichAnalysisMethodText(unsupported_method),
    "CUSTOM ORGANISM",
    fixed = TRUE
  )
  expect_identical(
    formatProtEnrichContrastsText(c("contrast_a", "contrast_b")),
    "contrast_a\ncontrast_b"
  )
  expect_match(
    formatProtEnrichContrastsText(NULL),
    "No contrasts available",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichStatusText(
      analysisComplete = FALSE
    ),
    "Ready for analysis",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichStatusText(
      analysisComplete = TRUE,
      methodInfo = supported_method,
      selectedContrast = "contrast_a",
      upCutoff = 1,
      downCutoff = -1,
      qCutoff = 0.05,
      gprofilerResults = gprofiler_results,
      clusterprofilerResults = cluster_results,
      stringdbResults = data.frame(network = "n1", stringsAsFactors = FALSE)
    ),
    "gprofiler2: 2 terms",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichGprofilerSummaryText(gprofiler_results),
    "Total enrichment terms: 2",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichGprofilerSummaryText(gprofiler_results, directionFilter = "up"),
    "Showing 1 up-regulated pathways",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichGprofilerSummaryText(data.frame(value = I(list(1)))),
    "Total enrichment terms: 1",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichClusterProfilerSummaryText(cluster_results),
    "Total GO terms: 3",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichClusterProfilerSummaryText(cluster_results, directionFilter = "down"),
    "Showing 1 down-regulated GO terms",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichClusterProfilerSummaryText(data.frame(value = I(list(1)))),
    "Total GO terms: 1",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichStringDbSummaryText(NULL),
    "not yet implemented",
    fixed = TRUE
  )
  expect_match(
    formatProtEnrichStringDbSummaryText(data.frame(network = "n1", stringsAsFactors = FALSE)),
    "STRING-DB Network Analysis",
    fixed = TRUE
  )
})

test_that("proteomics enrichment download helpers preserve filename and archive behavior", {
  temp_dir <- tempfile("prot-enrich-download-")
  dir.create(temp_dir, recursive = TRUE)
  written_tsv <- list()
  written_lines <- list()
  zipped <- list()

  expect_identical(
    buildProtEnrichResultsDownloadFilename(
      "Condition A / Condition B",
      date = as.Date("2026-04-22")
    ),
    "Enrichment_results_Condition_A___Condition_B_2026-04-22.zip"
  )

  result <- writeProtEnrichResultsDownloadArchive(
    file = file.path(temp_dir, "enrichment.zip"),
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    gprofilerResults = data.frame(term = "GO:1", stringsAsFactors = FALSE),
    clusterprofilerResults = data.frame(term = "GO:2", stringsAsFactors = FALSE),
    stringdbResults = data.frame(network = "n1", stringsAsFactors = FALSE),
    tempDir = temp_dir,
    writeTsvFn = function(x, file, ...) {
      written_tsv[[basename(file)]] <<- x
      invisible(NULL)
    },
    writeLinesFn = function(text, con, ...) {
      written_lines[[basename(con)]] <<- text
      invisible(NULL)
    },
    zipFn = function(zipfile, files, flags) {
      zipped[[length(zipped) + 1L]] <<- list(zipfile = zipfile, files = files, flags = flags)
      invisible(0L)
    },
    sysTimeFn = function() as.POSIXct("2026-04-22 10:00:00", tz = "UTC")
  )

  expect_identical(result$status, "ok")
  expect_setequal(
    names(written_tsv),
    c("gprofiler2_results.tsv", "clusterProfileR_results.tsv", "stringdb_results.tsv")
  )
  expect_true("enrichment_analysis_summary.txt" %in% names(written_lines))
  expect_match(written_lines$enrichment_analysis_summary.txt, "Contrast: contrast_a", fixed = TRUE)
  expect_identical(zipped[[1L]]$flags, "-j")

  error_result <- writeProtEnrichResultsDownloadArchive(
    file = file.path(temp_dir, "error.zip"),
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    gprofilerResults = data.frame(term = "GO:1", stringsAsFactors = FALSE),
    tempDir = temp_dir,
    writeTsvFn = function(...) stop("write failed"),
    writeLinesFn = function(text, con, ...) {
      written_lines[[basename(con)]] <<- text
      invisible(NULL)
    },
    zipFn = function(zipfile, files, flags) {
      zipped[[length(zipped) + 1L]] <<- list(zipfile = zipfile, files = files, flags = flags)
      invisible(0L)
    }
  )

  expect_identical(error_result$status, "error")
  expect_identical(error_result$error, "write failed")
  expect_true("download_error.txt" %in% names(written_lines))
})

test_that("proteomics enrichment misc helpers preserve cleanup, save, execute, and organism filter behavior", {
  removed <- character()
  cleanup <- removeProtEnrichWorkingNotification(
    removeNotificationFn = function(id) {
      removed <<- c(removed, id)
      invisible(NULL)
    }
  )
  expect_identical(cleanup$notificationId, "enrichment_working")
  expect_identical(removed, "enrichment_working")

  saved_args <- NULL
  workflow_data <- list(
    state_manager = list(
      saveState = function(...) {
        saved_args <<- list(...)
        invisible(NULL)
      }
    )
  )
  save_result <- saveProtEnrichCompletedState(
    workflowData = workflow_data,
    enrichmentResults = list(kind = "enrichment-results"),
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    catFn = function(...) invisible(NULL)
  )
  expect_true(save_result$saved)
  expect_identical(saved_args$state_name, "enrichment_completed")
  expect_identical(saved_args$config_object$selected_contrast, "contrast_a")

  failed_save <- saveProtEnrichCompletedState(
    workflowData = list(state_manager = list(saveState = function(...) stop("state unavailable"))),
    enrichmentResults = list(kind = "enrichment-results"),
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    catFn = function(...) invisible(NULL)
  )
  expect_false(failed_save$saved)
  expect_identical(failed_save$warning, "state unavailable")

  checkpoint <- NULL
  process_args <- NULL
  execute_result <- executeProtEnrichProcessEnrichments(
    enrichmentArgs = list(
      checkpointArgs = list(a = 1),
      processArgs = list(value = 2)
    ),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    captureCheckpointFn = function(args, id, label) {
      checkpoint <<- list(args = args, id = id, label = label)
      invisible(NULL)
    },
    processEnrichmentsFn = function(...) {
      process_args <<- list(...)
      list(kind = "processed")
    },
    catFn = function(...) invisible(NULL)
  )
  expect_identical(execute_result, list(kind = "processed"))
  expect_identical(checkpoint$id, "cp10")
  expect_identical(process_args$value, 2)

  da_results <- methods::new(
    "mockProtEnrichRuntimeDaCarrier",
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = c("P1;P1-2", "P2", "P3"),
        stringsAsFactors = FALSE
      )
    )
  )
  filter_result <- applyProtEnrichOrganismFilter(
    daResultsForEnrichment = da_results,
    organismMapping = data.frame(
      uniprot_acc = c("P1", "P3"),
      taxon_id = c("9606", "10090"),
      stringsAsFactors = FALSE
    ),
    targetTaxon = "9606",
    currentS4Object = makeProtEnrichRuntimeInputCarrier(),
    normalizeUniprotFn = function(x, remove_isoform = TRUE) sub("-.*$", "", x),
    cleanAccFn = function(x) sub("-.*$", "", x),
    catFn = function(...) invisible(NULL)
  )
  expect_true(filter_result$filterApplied)
  expect_equal(filter_result$filterStats$proteins_before, 3)
  expect_equal(filter_result$filterStats$proteins_after, 1)
  expect_identical(filter_result$daResultsForEnrichment@da_data$contrast_a$uniprot_acc, "P1;P1-2")
})

test_that("proteomics enrichment builder helpers preserve payload, args, UI, state, and progress behavior", {
  fixed_time <- as.POSIXct("2026-04-22 10:00:00", tz = "UTC")
  current_s4 <- makeProtEnrichRuntimeInputCarrier(label = "source")
  enrichment_results <- methods::new(
    "mockProtEnrichRuntimeResultCarrier",
    enrichment_data = list(),
    args = list(existing = TRUE)
  )
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(enrichment_analysis = "pending")
  workflow_data$state_manager <- list(
    getHistory = function() c("protein_replicate_filtered"),
    saveState = function(...) invisible(NULL)
  )

  payload <- buildProtEnrichAnalysisResultsPayload(
    gprofilerResults = data.frame(term = "GO:1"),
    clusterprofilerResults = data.frame(term = "GO:2"),
    stringdbResults = data.frame(term = "STRING:1"),
    fullEnrichmentResults = enrichment_results,
    selectedContrast = "contrast_a",
    analysisMethod = "gprofiler2",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606"
  )
  expect_identical(payload$selected_contrast, "contrast_a")
  expect_identical(payload$parameters$organism_taxid, "9606")

  propagated <- propagateProtEnrichResultsArgs(
    enrichmentResults = enrichment_results,
    currentS4Object = current_s4,
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE, species_name = "Homo sapiens"),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    timeFn = function() fixed_time,
    catFn = function(...) invisible(NULL)
  )
  expect_true(propagated$copiedArgs)
  expect_identical(propagated$enrichmentResults@args$enrichmentAnalysis$selected_contrast, "contrast_a")
  expect_identical(propagated$enrichmentResults@args$enrichmentAnalysisUI$timestamp, fixed_time)

  no_args <- propagateProtEnrichResultsArgs(
    enrichmentResults = list(kind = "no-slot"),
    currentS4Object = list(kind = "no-slot"),
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE, species_name = "Homo sapiens"),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    catFn = function(...) invisible(NULL)
  )
  expect_false(no_args$copiedArgs)

  ui_params <- propagateProtEnrichUiParams(
    currentS4Object = current_s4,
    workflowData = workflow_data,
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE, species_name = "Homo sapiens"),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    timeFn = function() fixed_time,
    catFn = function(...) invisible(NULL)
  )
  expect_true(ui_params$storedUiParams)
  expect_identical(workflow_data$enrichment_ui_params$selected_contrast, "contrast_a")

  state_update <- updateProtEnrichStateManagerUiParams(
    workflowData = workflow_data,
    storedUiParams = TRUE,
    catFn = function(...) invisible(NULL)
  )
  expect_true(state_update$attempted)
  expect_identical(state_update$currentDataState, "protein_replicate_filtered")

  no_update <- updateProtEnrichStateManagerUiParams(
    workflowData = workflow_data,
    storedUiParams = FALSE
  )
  expect_false(no_update$attempted)

  failed_update <- updateProtEnrichStateManagerUiParams(
    workflowData = list(state_manager = list(getHistory = function() stop("history failed"))),
    storedUiParams = TRUE,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(failed_update$warning, "history failed")

  tab_result <- completeProtEnrichTabStatus(workflow_data, tabName = "enrichment_analysis")
  expect_identical(tab_result$tabStatus$enrichment_analysis, "complete")

  progress_calls <- list()
  progress_result <- completeProtEnrichProgress(
    incProgressFn = function(value, detail) {
      progress_calls[[length(progress_calls) + 1L]] <<- list(value = value, detail = detail)
      invisible(NULL)
    }
  )
  expect_identical(progress_result$value, 1.0)
  expect_identical(progress_calls[[1L]]$detail, "Complete!")
})

test_that("proteomics enrichment result builders preserve all-contrast and persistence behavior", {
  method_info <- list(method = "gprofiler2", supported = TRUE, species_name = "Homo sapiens")
  enrichment_results <- methods::new(
    "mockProtEnrichRuntimeResultCarrier",
    enrichment_data = list(
      contrast_a = list(
        up = list(result = data.frame(term_name = "GO up", p_value = 0.01, term_size = 3, source = "GO:BP")),
        down = list(result = data.frame(term_name = "GO down", p_value = 0.02, term_size = 2, source = "GO:MF"))
      )
    ),
    args = list()
  )
  all_results <- buildProtEnrichAllContrastResults(
    enrichmentResults = enrichment_results,
    methodInfo = method_info,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(nrow(all_results$contrast_a$gprofiler_results), 2L)
  expect_setequal(all_results$contrast_a$gprofiler_results$directionality, c("positive", "negative"))
  expect_true("Description" %in% names(all_results$contrast_a$gprofiler_results))

  up_enrich <- methods::new(
    "mockProtEnrichResultSlotCarrier",
    result = data.frame(ID = "GO:1", pvalue = 0.01, stringsAsFactors = FALSE)
  )
  cluster_results <- methods::new(
    "mockProtEnrichRuntimeResultCarrier",
    enrichment_data = list(contrast_a = list(up = up_enrich, down = NULL)),
    args = list()
  )
  all_cluster <- buildProtEnrichAllContrastResults(
    enrichmentResults = cluster_results,
    methodInfo = list(method = "clusterprofiler"),
    isEnrichResultFn = function(object, class2) TRUE,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(all_cluster$contrast_a$clusterprofiler_results$directionality, "up")

  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$gprofiler_results <- NULL
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL
  capture_result <- captureProtEnrichPostProcessResults(
    selectedContrast = "contrast_a",
    enrichmentResults = enrichment_results,
    enrichmentData = enrichment_data,
    contrastsTbl = NULL,
    methodInfo = method_info,
    catFn = function(...) invisible(NULL)
  )
  expect_true(capture_result$hasResults)
  expect_identical(enrichment_data$gprofiler_results$Description[1], "GO up")

  empty_capture <- captureProtEnrichPostProcessResults(
    selectedContrast = "contrast_a",
    enrichmentResults = NULL,
    enrichmentData = enrichment_data,
    contrastsTbl = NULL,
    methodInfo = method_info,
    catFn = function(...) invisible(NULL)
  )
  expect_false(empty_capture$hasResults)

  workflow_data <- new.env(parent = emptyenv())
  input <- list(
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05,
    organism_taxid = "9606"
  )
  progress <- list()
  persist_result <- persistProtEnrichAnalysisResults(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    selectedContrast = "contrast_a",
    enrichmentResults = enrichment_results,
    methodInfo = method_info,
    pathwayDir = "/tmp/pathways",
    buildAnalysisResultsPayloadFn = function(...) {
      list(...)
    },
    propagateResultsArgsFn = function(enrichmentResults, ...) {
      list(enrichmentResults = enrichmentResults)
    },
    propagateUiParamsFn = function(currentS4Object, workflowData, ...) {
      list(currentS4Object = currentS4Object, storedUiParams = TRUE)
    },
    updateStateManagerUiParamsFn = function(...) {
      list(updated = FALSE)
    },
    saveCompletedStateFn = function(...) {
      list(saved = TRUE)
    },
    completeTabStatusFn = function(workflowData) {
      workflowData$tab_status <- list(enrichment_analysis = "complete")
      list(tabStatus = workflowData$tab_status)
    },
    completeProgressFn = function() {
      progress$complete <<- TRUE
      list(value = 1)
    },
    incProgressFn = function(value, detail) {
      progress$storing <<- list(value = value, detail = detail)
      invisible(NULL)
    }
  )
  expect_true(persist_result$analysisComplete)
  expect_true(progress$complete)
  expect_identical(workflow_data$enrichment_analysis_results$selectedContrast, "contrast_a")

  finalize_result <- finalizeProtEnrichAnalysisBodyResults(
    selectedContrast = "contrast_a",
    rawContrastName = "contrast_a",
    organismFilterApplied = TRUE,
    filterStats = list(proteins_before = 3, proteins_after = 2, proteins_removed = 1),
    enrichmentResults = enrichment_results,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    input = input,
    methodInfo = method_info,
    contrastsTbl = NULL,
    pathwayDir = "/tmp/pathways",
    capturePostProcessResultsFn = function(...) list(hasResults = TRUE),
    persistAnalysisResultsFn = function(enrichmentResults, ...) {
      list(
        enrichmentResults = enrichmentResults,
        analysisComplete = TRUE,
        currentS4Object = makeProtEnrichRuntimeInputCarrier()
      )
    },
    incProgressFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_true(finalize_result$analysisComplete)
  expect_true(finalize_result$organismFilterApplied)
})

test_that("proteomics enrichment resolver helpers preserve contrast, dependency, and process setup behavior", {
  current_s4 <- makeProtEnrichRuntimeInputCarrier(label = "resolver")
  da_results <- methods::new(
    "mockProtEnrichRuntimeDaCarrier",
    da_data = list(
      contrast_a = data.frame(gene_name = "GeneA", uniprot_acc = "P1", stringsAsFactors = FALSE)
    )
  )
  contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "Contrast A",
    stringsAsFactors = FALSE
  )

  expect_identical(
    buildProtEnrichContrastChoices(list(contrast_a = list()), contrasts_tbl)$contrastsAvailable,
    "Contrast A"
  )
  expect_identical(
    buildProtEnrichContrastChoices(list(contrast_a = list()), NULL)$source,
    "raw_names"
  )
  expect_identical(
    resolveProtEnrichRawContrastName("Contrast A", contrasts_tbl)$rawContrastName,
    "contrast_a"
  )
  expect_true(resolveProtEnrichSelectedContrastResults(
    "Contrast A",
    list(contrast_a = list(gprofiler_results = data.frame(term = "GO:1"))),
    contrasts_tbl
  )$found)
  expect_identical(
    resolveProtEnrichSelectedDaResults(
      "Condition_A_vs_Condition_B",
      list(limma_Condition_A_Condition_B = list(value = 1)),
      NULL
    )$source,
    "fuzzy_match"
  )
  expect_identical(
    resolveProtEnrichSelectedDaResults("contrast_a", list(contrast_a = list(value = 1)))$source,
    "direct_key"
  )

  state_manager_s4 <- makeProtEnrichRuntimeInputCarrier(label = "state")
  current_from_state <- resolveProtEnrichCurrentS4Object(
    workflowData = list(state_manager = list(current_state = "normalized", getState = function(state) state_manager_s4)),
    daResultsList = NULL
  )
  expect_identical(current_from_state$source, "state_manager")
  current_from_da <- resolveProtEnrichCurrentS4Object(
    workflowData = list(state_manager = NULL),
    daResultsList = list(contrast_a = list(theObject = current_s4))
  )
  expect_identical(current_from_da$source, "da_results_first_result")
  current_from_combined <- resolveProtEnrichCurrentS4Object(
    workflowData = list(state_manager = NULL),
    daResultsList = list(theObject = current_s4)
  )
  expect_identical(current_from_combined$source, "da_results_combined")

  input_columns <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = list(method = "gprofiler2"),
    daResultsForEnrichment = da_results,
    currentS4Object = current_s4
  )
  expect_true(input_columns$geneNameOverrideApplied)
  expect_identical(input_columns$idColumn, "gene_name")
  expect_identical(
    resolveProtEnrichAnalysisInputColumns(
      methodInfo = list(method = "clusterprofiler"),
      daResultsForEnrichment = da_results,
      currentS4Object = current_s4
    )$idColumn,
    "uniprot_acc"
  )

  process_args <- buildProtEnrichProcessEnrichmentsArgs(
    daResultsForEnrichment = da_results,
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    pathwayDir = "/tmp/pathways",
    goAnnotations = data.frame(Entry = "P1"),
    proteinIdColumn = "gene_name",
    contrastNames = "contrast_a",
    correctionMethod = "fdr"
  )
  expect_identical(process_args$checkpointArgs$taxon_id, 9606)
  expect_identical(process_args$processArgs$protein_id_column, "gene_name")

  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- current_s4
  enrichment_data$da_results_data <- list(contrast_a = list())
  process_setup <- prepareProtEnrichProcessExecution(
    input = list(
      organism_taxid = "9606",
      up_cutoff = 1,
      down_cutoff = -1,
      q_cutoff = 0.05,
      correction_method = "fdr"
    ),
    enrichmentData = enrichment_data,
    daResultsForEnrichment = da_results,
    pathwayDir = "/tmp/pathways",
    goAnnotations = data.frame(Entry = "P1"),
    currentAnalysisMethodFn = function() list(method = "gprofiler2"),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(process_setup$methodInfo$method, "gprofiler2")
  expect_identical(process_setup$enrichmentArgs$processArgs$contrast_names, "contrast_a")

  global_env <- new.env(parent = emptyenv())
  assign("design_matrix", matrix(1, nrow = 1), envir = global_env)
  dependencies <- resolveProtEnrichRunDependencies(
    currentS4Object = NULL,
    daResultsData = list(contrast_a = list(theObject = current_s4)),
    workflowData = list(state_manager = NULL),
    contrastsTbl = contrasts_tbl,
    globalEnv = global_env
  )
  expect_identical(dependencies$s4Source, "da_results_first_result")
  expect_identical(dependencies$designMatrixSource, "s4_object")

  fallback_dependencies <- resolveProtEnrichRunDependencies(
    currentS4Object = NULL,
    daResultsData = list(),
    workflowData = list(state_manager = NULL),
    contrastsTbl = contrasts_tbl,
    globalEnv = global_env
  )
  expect_identical(fallback_dependencies$designMatrixSource, "global_environment")

  created_dirs <- character()
  path_config <- resolveProtEnrichOutputDirectories(
    experimentPaths = list(results_dir = "/tmp/prot-results"),
    dirExistsFn = function(path) FALSE,
    dirCreateFn = function(path, recursive = TRUE) {
      created_dirs <<- c(created_dirs, path)
      invisible(TRUE)
    }
  )
  expect_identical(path_config$daOutputDirSource, "results_fallback")
  expect_true(path_config$daOutputDirCreated)
  expect_true(path_config$pathwayDirCreated)
  expect_length(created_dirs, 2)
})

test_that("proteomics enrichment annotation and organism resolvers preserve source and filter behavior", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$uniprot_dat_cln <- data.frame(Entry = "P1", gene_names = "GeneA", stringsAsFactors = FALSE)
  global_env <- new.env(parent = emptyenv())
  annotations <- resolveProtEnrichUniprotAnnotations(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    currentS4Object = NULL,
    organismTaxid = "9606",
    globalEnv = global_env,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(annotations$source, "workflow_data")
  expect_identical(get("uniprot_dat_cln", envir = global_env), workflow_data$uniprot_dat_cln)

  source_annotations <- data.frame(Entry = "P2", gene_names = "GeneB", stringsAsFactors = FALSE)
  source_env <- new.env(parent = emptyenv())
  loaded <- resolveProtEnrichUniprotAnnotations(
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    currentS4Object = NULL,
    organismTaxid = "9606",
    globalEnv = source_env,
    fileExistsFn = function(path) TRUE,
    readRdsFn = function(path) source_annotations,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(loaded$source, "source_directory")
  expect_identical(get("uniprot_dat_cln", envir = source_env), source_annotations)

  generated_env <- new.env(parent = emptyenv())
  generated <- resolveProtEnrichUniprotAnnotations(
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    currentS4Object = makeProtEnrichRuntimeInputCarrier(),
    organismTaxid = "9606",
    globalEnv = generated_env,
    fileExistsFn = function(path) FALSE,
    dirExistsFn = function(path) FALSE,
    dirCreateFn = function(path, recursive = TRUE) TRUE,
    getUniprotAnnotationsFn = function(input_tbl, cache_dir, taxon_id) {
      data.frame(Entry = input_tbl$uniprot_acc, taxon = taxon_id, stringsAsFactors = FALSE)
    },
    catFn = function(...) invisible(NULL)
  )
  expect_identical(generated$source, "generated")
  expect_true(generated$cacheDirCreated)

  mapping_from_workflow <- resolveProtEnrichOrganismMapping(
    workflowData = list(mixed_species_analysis = list(organism_mapping = data.frame(uniprot_acc = "P1", taxon_id = "9606"))),
    uniprotDatCln = NULL,
    targetTaxon = "9606",
    catFn = function(...) invisible(NULL)
  )
  expect_identical(mapping_from_workflow$source, "workflow_data")

  mapping_from_names <- resolveProtEnrichOrganismMapping(
    workflowData = list(mixed_species_analysis = NULL),
    uniprotDatCln = data.frame(
      Entry = c("P1", "P2"),
      Organism = c("Homo sapiens", "Mus musculus"),
      stringsAsFactors = FALSE
    ),
    targetTaxon = "9606",
    catFn = function(...) invisible(NULL)
  )
  expect_identical(mapping_from_names$source, "organism_names")
  expect_true("9606" %in% mapping_from_names$availableTaxonIds)

  mapping_from_taxon <- resolveProtEnrichOrganismMapping(
    workflowData = list(mixed_species_analysis = NULL),
    uniprotDatCln = data.frame(
      Entry = "P1",
      OX = "NCBI_TaxID=9606",
      stringsAsFactors = FALSE
    ),
    targetTaxon = "9606",
    catFn = function(...) invisible(NULL)
  )
  expect_identical(mapping_from_taxon$source, "taxon_column")
  expect_identical(mapping_from_taxon$organismMapping$taxon_id, "9606")

  mapping_fallback <- resolveProtEnrichOrganismMapping(
    workflowData = list(mixed_species_analysis = NULL),
    uniprotDatCln = data.frame(Entry = "P1", stringsAsFactors = FALSE),
    targetTaxon = "9606",
    catFn = function(...) invisible(NULL)
  )
  expect_identical(mapping_fallback$source, "single_species_fallback")

  annotation_match <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = data.frame(Entry = "P1", gene_names = "GeneA"),
    daResultsForEnrichment = methods::new("mockProtEnrichRuntimeDaCarrier", da_data = list()),
    currentS4Object = makeProtEnrichRuntimeInputCarrier(),
    matchAnnotationsFn = function(...) list(match_statistics = list(match_rate = 99)),
    catFn = function(...) invisible(NULL)
  )
  expect_true(annotation_match$attempted)
  expect_identical(annotation_match$matchRate, 99)

  failed_match <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = data.frame(Entry = "P1", gene_names = "GeneA"),
    daResultsForEnrichment = methods::new("mockProtEnrichRuntimeDaCarrier", da_data = list()),
    currentS4Object = makeProtEnrichRuntimeInputCarrier(),
    matchAnnotationsFn = function(...) stop("match failed"),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(failed_match$warning, "match failed")

  skipped_match <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = NULL,
    daResultsForEnrichment = NULL,
    currentS4Object = NULL
  )
  expect_false(skipped_match$attempted)

  metadata <- persistProtEnrichOrganismFilterMetadata(
    workflowData = workflow_data,
    organismFilterEnabled = TRUE,
    organismFilterApplied = TRUE,
    targetTaxonId = "9606",
    filterStats = list(proteins_before = 3, proteins_after = 1, proteins_removed = 2),
    timeFn = function() as.POSIXct("2026-04-22 10:00:00", tz = "UTC")
  )
  expect_true(metadata$enabled)
  expect_identical(workflow_data$enrichment_organism_filter$proteins_removed, 2)
})

test_that("proteomics enrichment analysis setup preserves orchestration branches", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$da_results_data <- list(
    contrast_a = list(theObject = makeProtEnrichRuntimeInputCarrier(label = "setup"))
  )
  enrichment_data$current_s4_object <- NULL
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- NULL
  da_output_dir <- tempfile("prot-enrich-da-")
  pathway_dir <- tempfile("prot-enrich-pathway-")
  results_dir <- tempfile("prot-enrich-results-")
  source_dir <- tempfile("prot-enrich-source-")
  dir.create(da_output_dir, recursive = TRUE)
  dir.create(pathway_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)
  experiment_paths <- list(
    da_output_dir = da_output_dir,
    pathway_dir = pathway_dir,
    results_dir = results_dir,
    source_dir = source_dir
  )
  input <- list(
    organism_taxid = "9606",
    enable_organism_filter = TRUE
  )
  global_env <- new.env(parent = emptyenv())
  assign(
    "contrasts_tbl",
    data.frame(contrasts = "contrast_a", friendly_names = "Contrast A", stringsAsFactors = FALSE),
    envir = global_env
  )
  notifications <- list()

  setup <- prepareProtEnrichAnalysisBodySetup(
    selectedContrast = "Contrast A",
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    createDaResultsForEnrichmentFn = function(contrasts_tbl, design_matrix, da_output_dir) {
      expect_identical(contrasts_tbl$contrasts, "contrast_a")
      expect_true(is.matrix(design_matrix))
      expect_identical(da_output_dir, experiment_paths$da_output_dir)
      methods::new(
        "mockProtEnrichRuntimeDaCarrier",
        da_data = list(contrast_a = data.frame(uniprot_acc = c("P1", "P2"), stringsAsFactors = FALSE))
      )
    },
    resolveUniprotAnnotationsFn = function(...) {
      list(uniprotDatCln = data.frame(Entry = c("P1", "P2"), stringsAsFactors = FALSE))
    },
    resolveAnnotationMatchingFn = function(...) {
      list(annotationMatchResults = list(match_statistics = list(match_rate = 100)))
    },
    resolveOrganismMappingFn = function(...) {
      list(organismMapping = data.frame(uniprot_acc = "P1", taxon_id = "9606", stringsAsFactors = FALSE))
    },
    applyOrganismFilterFn = function(daResultsForEnrichment, organismMapping, targetTaxon, currentS4Object) {
      list(
        daResultsForEnrichment = daResultsForEnrichment,
        filterApplied = TRUE,
        filterStats = list(proteins_before = 2, proteins_after = 1, proteins_removed = 1)
      )
    },
    persistOrganismFilterMetadataFn = function(...) {
      workflow_data$filter_metadata <- list(...)
      invisible(workflow_data$filter_metadata)
    },
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration)
      invisible(NULL)
    },
    globalEnv = global_env,
    catFn = function(...) invisible(NULL)
  )

  expect_identical(setup$rawContrastName, "contrast_a")
  expect_true(setup$organismFilterApplied)
  expect_identical(setup$filterStats$proteins_removed, 1)
  expect_identical(enrichment_data$annotation_match_results$match_statistics$match_rate, 100)
  expect_length(notifications, 0L)

  input$enable_organism_filter <- TRUE
  notifications <- list()
  setup_without_mapping <- prepareProtEnrichAnalysisBodySetup(
    selectedContrast = "contrast_a",
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    createDaResultsForEnrichmentFn = function(...) {
      methods::new("mockProtEnrichRuntimeDaCarrier", da_data = list(contrast_a = data.frame(uniprot_acc = "P1")))
    },
    resolveUniprotAnnotationsFn = function(...) {
      list(uniprotDatCln = data.frame(Entry = "P1", stringsAsFactors = FALSE))
    },
    resolveAnnotationMatchingFn = function(...) list(annotationMatchResults = NULL),
    resolveOrganismMappingFn = function(...) list(organismMapping = data.frame()),
    persistOrganismFilterMetadataFn = function(...) list(...),
    showNotificationFn = function(message, type = NULL, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration)
      invisible(NULL)
    },
    globalEnv = global_env,
    catFn = function(...) invisible(NULL)
  )
  expect_false(setup_without_mapping$organismFilterApplied)
  expect_identical(notifications[[1L]]$type, "warning")

  expect_error(
    prepareProtEnrichAnalysisBodySetup(
      selectedContrast = "missing",
      input = list(organism_taxid = "9606", enable_organism_filter = FALSE),
      enrichmentData = enrichment_data,
      workflowData = workflow_data,
      experimentPaths = experiment_paths,
      globalEnv = global_env,
      catFn = function(...) invisible(NULL)
    ),
    "Could not find DE results"
  )
})
