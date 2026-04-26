# fidelity-coverage-compare: shared
# testthat for Proteomics Enrichment helper contracts

skipIfMissingMultiScholaRBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    testthat::skip(sprintf("Target-only extracted helper(s) not present: %s", paste(missing, collapse = ", ")))
  }
}

skipIfMissingMultiScholaRBindings(
  "resolveProtEnrichCurrentS4Object",
  "buildProtEnrichContrastChoices",
  "setupProtEnrichReactiveValues"
)

if (!methods::isClass("mockProtEnrichPlotContainer")) {
  methods::setClass("mockProtEnrichPlotContainer", slots = c(enrichment_plotly = "list"))
}

if (!methods::isClass("mockProtEnrichDesignCarrier")) {
  methods::setClass("mockProtEnrichDesignCarrier", slots = c(design_matrix = "matrix"))
}

if (!methods::isClass("mockProtEnrichProteinCarrier")) {
  methods::setClass("mockProtEnrichProteinCarrier", slots = c(protein_quant_table = "data.frame"))
}

if (!methods::isClass("mockProtEnrichAnnotationCarrier")) {
  methods::setClass("mockProtEnrichAnnotationCarrier", slots = c(protein_id_column = "character"))
}

if (!methods::isClass("mockProtEnrichAllContrastCarrier")) {
  methods::setClass("mockProtEnrichAllContrastCarrier", slots = c(enrichment_data = "list"))
}

if (!methods::isClass("mockProtEnrichClusterResultCarrier")) {
  methods::setClass("mockProtEnrichClusterResultCarrier", slots = c(result = "data.frame"))
}

if (!methods::isClass("mockProtEnrichArgsCarrier")) {
  methods::setClass("mockProtEnrichArgsCarrier", slots = c(args = "list"))
}

test_that("resolveProtEnrichCurrentS4Object prefers state manager state", {
  state_object <- structure(list(label = "state"), class = "mock_state_object")
  workflow_data <- list(
    state_manager = list(
      current_state = "normalized",
      getState = function(state_name) {
        expect_equal(state_name, "normalized")
        state_object
      }
    )
  )
  da_results_list <- list(
    contrast_a = list(theObject = structure(list(label = "fallback"), class = "mock_fallback_object"))
  )

  resolved <- resolveProtEnrichCurrentS4Object(workflow_data, da_results_list)

  expect_identical(resolved$currentS4, state_object)
  expect_equal(resolved$source, "state_manager")
})

test_that("resolveProtEnrichCurrentS4Object falls back to first DE result object", {
  fallback_object <- structure(list(label = "first_result"), class = "mock_first_result_object")
  workflow_data <- list(
    state_manager = list(
      current_state = "normalized",
      getState = function(state_name) {
        expect_equal(state_name, "normalized")
        NULL
      }
    )
  )
  da_results_list <- list(
    contrast_a = list(theObject = fallback_object),
    contrast_b = list(theObject = structure(list(label = "unused"), class = "mock_unused_object"))
  )

  resolved <- resolveProtEnrichCurrentS4Object(workflow_data, da_results_list)

  expect_identical(resolved$currentS4, fallback_object)
  expect_equal(resolved$source, "da_results_first_result")
})

test_that("resolveProtEnrichCurrentS4Object falls back to combined DE object", {
  combined_object <- structure(list(label = "combined"), class = "mock_combined_object")
  workflow_data <- list(state_manager = NULL)
  da_results_list <- list(theObject = combined_object)

  resolved <- resolveProtEnrichCurrentS4Object(workflow_data, da_results_list)

  expect_identical(resolved$currentS4, combined_object)
  expect_equal(resolved$source, "da_results_combined")
})

test_that("buildProtEnrichContrastChoices prefers friendly names when available", {
  da_results_list <- list(contrast_a = list(), contrast_b = list())
  contrasts_tbl <- data.frame(
    contrasts = c("contrast_a", "contrast_b"),
    friendly_names = c("A vs B", "C vs D"),
    stringsAsFactors = FALSE
  )

  contrast_config <- buildProtEnrichContrastChoices(da_results_list, contrasts_tbl)

  expect_equal(contrast_config$contrastsAvailable, c("A vs B", "C vs D"))
  expect_equal(unname(contrast_config$contrastChoices), c("A vs B", "C vs D"))
  expect_equal(names(contrast_config$contrastChoices), c("A vs B", "C vs D"))
  expect_equal(contrast_config$source, "friendly_names")
})

test_that("buildProtEnrichContrastChoices falls back to raw contrast names", {
  da_results_list <- list(contrast_a = list(), contrast_b = list())

  contrast_config <- buildProtEnrichContrastChoices(da_results_list)

  expect_equal(contrast_config$contrastsAvailable, c("contrast_a", "contrast_b"))
  expect_equal(unname(contrast_config$contrastChoices), c("contrast_a", "contrast_b"))
  expect_equal(names(contrast_config$contrastChoices), c("contrast_a", "contrast_b"))
  expect_equal(contrast_config$source, "raw_names")
})

test_that("resolveProtEnrichRawContrastName maps friendly names to raw keys", {
  contrasts_tbl <- data.frame(
    contrasts = c("contrast_a", "contrast_b"),
    friendly_names = c("A vs B", "C vs D"),
    stringsAsFactors = FALSE
  )

  resolved <- resolveProtEnrichRawContrastName("A vs B", contrasts_tbl)

  expect_equal(resolved$rawContrastName, "contrast_a")
  expect_equal(resolved$source, "friendly_name")
})

test_that("createProtEnrichRawContrastNameReactive reads contrasts_tbl and delegates raw-name resolution", {
  lookup_env <- new.env(parent = emptyenv())
  lookup_env$contrasts_tbl <- data.frame(
    contrasts = c("contrast_a", "contrast_b"),
    friendly_names = c("A vs B", "C vs D"),
    stringsAsFactors = FALSE
  )
  input <- shiny::reactiveValues(selected_contrast = "A vs B")
  recorded <- new.env(parent = emptyenv())
  recorded$selectedContrast <- NULL
  recorded$contrastsTbl <- NULL

  raw_contrast_name <- createProtEnrichRawContrastNameReactive(
    input = input,
    globalEnv = lookup_env,
    resolveRawContrastNameFn = function(selectedContrast, contrastsTbl) {
      recorded$selectedContrast <- selectedContrast
      recorded$contrastsTbl <- contrastsTbl
      list(rawContrastName = "contrast_a", source = "friendly_name")
    }
  )

  expect_equal(shiny::isolate(raw_contrast_name()), "contrast_a")
  expect_equal(recorded$selectedContrast, "A vs B")
  expect_identical(recorded$contrastsTbl, lookup_env$contrasts_tbl)
})

test_that("setupProtEnrichRawContrastNameReactive delegates creation through the wrapper-local seam", {
  input <- shiny::reactiveValues(selected_contrast = "A vs B")
  recorded <- new.env(parent = emptyenv())
  recorded$input <- NULL
  raw_contrast_name <- shiny::reactive("contrast_a")

  registration <- setupProtEnrichRawContrastNameReactive(
    input = input,
    createRawContrastNameReactiveFn = function(input) {
      recorded$input <- input
      raw_contrast_name
    }
  )

  expect_identical(registration$rawContrastName, raw_contrast_name)
  expect_equal(registration$reason, "created")
  expect_identical(recorded$input, input)
  expect_equal(shiny::isolate(registration$rawContrastName()), "contrast_a")
})

test_that("createProtEnrichReactiveValues seeds the wrapper state container", {
  enrichment_data <- createProtEnrichReactiveValues()

  expect_null(shiny::isolate(enrichment_data$enrichment_results))
  expect_null(shiny::isolate(enrichment_data$contrasts_available))
  expect_false(isTRUE(shiny::isolate(enrichment_data$analysis_complete)))
  expect_null(shiny::isolate(enrichment_data$current_s4_object))
  expect_null(shiny::isolate(enrichment_data$da_results_data))
  expect_null(shiny::isolate(enrichment_data$gprofiler_results))
  expect_null(shiny::isolate(enrichment_data$clusterprofiler_results))
  expect_null(shiny::isolate(enrichment_data$stringdb_results))
  expect_null(shiny::isolate(enrichment_data$analysis_method))
  expect_null(shiny::isolate(enrichment_data$organism_supported))
  expect_equal(shiny::isolate(enrichment_data$all_enrichment_results), list())
  expect_equal(shiny::isolate(enrichment_data$current_contrast_results), list())
  expect_equal(shiny::isolate(enrichment_data$enrichment_plots), list())
})

test_that("setupProtEnrichReactiveValues delegates creation through the wrapper-local seam", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- 0L
  recorded$messages <- character()
  enrichment_data <- shiny::reactiveValues(analysis_complete = FALSE)

  registration <- setupProtEnrichReactiveValues(
    createReactiveValuesFn = function() {
      recorded$calls <- recorded$calls + 1L
      enrichment_data
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$calls, 1L)
  expect_equal(registration$reason, "created")
  expect_identical(registration$enrichmentData, enrichment_data)
  expect_equal(
    recorded$messages,
    "   mod_prot_enrich_server Step: Reactive values initialized\n"
  )
})

test_that("createProtEnrichSupportedOrganismsReactive delegates the static lookup through a reactive wrapper", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- 0L

  supported_organisms <- createProtEnrichSupportedOrganismsReactive(
    buildSupportedOrganismsFn = function() {
      recorded$calls <- recorded$calls + 1L

      tibble::tribble(
        ~taxid, ~id, ~name,
        "9606", "hsapiens", "Homo sapiens"
      )
    }
  )

  resolved <- shiny::isolate(supported_organisms())

  expect_equal(recorded$calls, 1L)
  expect_s3_class(resolved, "tbl_df")
  expect_equal(resolved$taxid, "9606")
  expect_equal(resolved$id, "hsapiens")
})

test_that("setupProtEnrichSupportedOrganismsReactive delegates creation through the wrapper-local seam", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- 0L

  registration <- setupProtEnrichSupportedOrganismsReactive(
    createSupportedOrganismsReactiveFn = function() {
      recorded$calls <- recorded$calls + 1L
      shiny::reactive(tibble::tribble(
        ~taxid, ~id, ~name,
        "9606", "hsapiens", "Homo sapiens"
      ))
    }
  )

  resolved <- shiny::isolate(registration$supportedOrganisms())

  expect_equal(recorded$calls, 1L)
  expect_equal(registration$reason, "created")
  expect_s3_class(resolved, "tbl_df")
  expect_equal(resolved$taxid, "9606")
  expect_equal(resolved$id, "hsapiens")
})

test_that("createProtEnrichCurrentAnalysisMethodReactive resolves method metadata and updates state", {
  input <- shiny::reactiveValues(organism_taxid = "9606")
  enrichment_data <- shiny::reactiveValues(
    analysis_method = NULL,
    organism_supported = NULL
  )
  recorded <- new.env(parent = emptyenv())
  recorded$organismTaxid <- NULL
  recorded$supportedOrganisms <- NULL
  recorded$supportedOrganismsCalls <- 0L

  current_analysis_method <- createProtEnrichCurrentAnalysisMethodReactive(
    input = input,
    enrichmentData = enrichment_data,
    supportedOrganismsFn = function() {
      recorded$supportedOrganismsCalls <- recorded$supportedOrganismsCalls + 1L

      tibble::tribble(
        ~taxid, ~id, ~name,
        "9606", "hsapiens", "Homo sapiens"
      )
    },
    resolveAnalysisMethodFn = function(organismTaxid, supportedOrganisms) {
      recorded$organismTaxid <- organismTaxid
      recorded$supportedOrganisms <- supportedOrganisms

      list(
        methodInfo = list(
          method = "gprofiler2",
          supported = TRUE,
          species_name = "Homo sapiens"
        ),
        analysisMethod = "gprofiler2",
        organismSupported = TRUE
      )
    }
  )

  method_info <- shiny::isolate(current_analysis_method())

  expect_equal(recorded$supportedOrganismsCalls, 1L)
  expect_equal(recorded$organismTaxid, "9606")
  expect_s3_class(recorded$supportedOrganisms, "tbl_df")
  expect_identical(method_info$method, "gprofiler2")
  expect_true(isTRUE(method_info$supported))
  expect_identical(shiny::isolate(enrichment_data$analysis_method), "gprofiler2")
  expect_true(isTRUE(shiny::isolate(enrichment_data$organism_supported)))
})

test_that("setupProtEnrichCurrentAnalysisMethodReactive delegates creation through the wrapper-local seam", {
  input <- shiny::reactiveValues(organism_taxid = "9606")
  enrichment_data <- shiny::reactiveValues(
    analysis_method = NULL,
    organism_supported = NULL
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL
  current_analysis_method <- shiny::reactive(list(
    method = "gprofiler2",
    supported = TRUE,
    species_name = "Homo sapiens"
  ))

  registration <- setupProtEnrichCurrentAnalysisMethodReactive(
    input = input,
    enrichmentData = enrichment_data,
    supportedOrganismsFn = function() {
      tibble::tribble(
        ~taxid, ~id, ~name,
        "9606", "hsapiens", "Homo sapiens"
      )
    },
    createCurrentAnalysisMethodReactiveFn = function(input,
                                                     enrichmentData,
                                                     supportedOrganismsFn) {
      recorded$args <- list(
        input = input,
        enrichmentData = enrichmentData,
        supportedOrganismsFn = supportedOrganismsFn,
        supportedOrganisms = shiny::isolate(supportedOrganismsFn())
      )

      current_analysis_method
    }
  )

  expect_identical(registration$currentAnalysisMethod, current_analysis_method)
  expect_equal(registration$reason, "created")
  expect_identical(recorded$args$input, input)
  expect_identical(recorded$args$enrichmentData, enrichment_data)
  expect_true(is.function(recorded$args$supportedOrganismsFn))
  expect_equal(recorded$args$supportedOrganisms$taxid, "9606")
})

test_that("setupProtEnrichAnalysisMethodBootstrap delegates supported-organisms and analysis-method setup through the wrapper-local seam", {
  input <- shiny::reactiveValues(organism_taxid = "9606")
  enrichment_data <- shiny::reactiveValues(
    analysis_method = NULL,
    organism_supported = NULL
  )
  recorded <- new.env(parent = emptyenv())
  recorded$supportedOrganismsCalls <- 0L
  recorded$currentAnalysisMethodArgs <- NULL
  supported_organisms <- shiny::reactive(tibble::tribble(
    ~taxid, ~id, ~name,
    "9606", "hsapiens", "Homo sapiens"
  ))
  current_analysis_method <- shiny::reactive(list(
    method = "gprofiler2",
    supported = TRUE,
    species_name = "Homo sapiens"
  ))

  bootstrap <- setupProtEnrichAnalysisMethodBootstrap(
    input = input,
    enrichmentData = enrichment_data,
    setupSupportedOrganismsReactiveFn = function() {
      recorded$supportedOrganismsCalls <- recorded$supportedOrganismsCalls + 1L
      list(
        supportedOrganisms = supported_organisms,
        reason = "created"
      )
    },
    setupCurrentAnalysisMethodReactiveFn = function(input,
                                                    enrichmentData,
                                                    supportedOrganismsFn) {
      recorded$currentAnalysisMethodArgs <- list(
        input = input,
        enrichmentData = enrichmentData,
        supportedOrganismsFn = supportedOrganismsFn,
        supportedOrganisms = shiny::isolate(supportedOrganismsFn())
      )
      list(
        currentAnalysisMethod = current_analysis_method,
        reason = "created"
      )
    }
  )

  expect_equal(recorded$supportedOrganismsCalls, 1L)
  expect_identical(bootstrap$supportedOrganisms, supported_organisms)
  expect_identical(bootstrap$currentAnalysisMethod, current_analysis_method)
  expect_equal(bootstrap$reason, "created")
  expect_identical(recorded$currentAnalysisMethodArgs$input, input)
  expect_identical(recorded$currentAnalysisMethodArgs$enrichmentData, enrichment_data)
  expect_identical(
    recorded$currentAnalysisMethodArgs$supportedOrganismsFn,
    supported_organisms
  )
  expect_equal(recorded$currentAnalysisMethodArgs$supportedOrganisms$taxid, "9606")
})

test_that("resolveProtEnrichSelectedContrastResults hydrates mapped contrast payloads", {
  contrasts_tbl <- data.frame(
    contrasts = c("contrast_a", "contrast_b"),
    friendly_names = c("A vs B", "C vs D"),
    stringsAsFactors = FALSE
  )
  all_enrichment_results <- list(
    contrast_a = list(
      gprofiler_results = data.frame(term = "go:1", stringsAsFactors = FALSE),
      clusterprofiler_results = NULL,
      stringdb_results = data.frame(network = "ppi", stringsAsFactors = FALSE)
    )
  )

  resolved <- resolveProtEnrichSelectedContrastResults(
    "A vs B",
    all_enrichment_results,
    contrasts_tbl
  )

  expect_true(resolved$found)
  expect_equal(resolved$rawContrastName, "contrast_a")
  expect_equal(resolved$source, "friendly_name")
  expect_equal(resolved$availableContrasts, "contrast_a")
  expect_equal(nrow(resolved$gprofilerResults), 1)
  expect_null(resolved$clusterprofilerResults)
  expect_equal(nrow(resolved$stringdbResults), 1)
})

test_that("resolveProtEnrichSelectedContrastResults clears payload for missing contrast", {
  all_enrichment_results <- list(
    contrast_a = list(
      gprofiler_results = data.frame(term = "go:1", stringsAsFactors = FALSE),
      clusterprofiler_results = data.frame(term = "cp:1", stringsAsFactors = FALSE),
      stringdb_results = data.frame(network = "ppi", stringsAsFactors = FALSE)
    )
  )

  resolved <- resolveProtEnrichSelectedContrastResults(
    "missing_contrast",
    all_enrichment_results
  )

  expect_false(resolved$found)
  expect_equal(resolved$rawContrastName, "missing_contrast")
  expect_equal(resolved$source, "selected_contrast")
  expect_equal(resolved$availableContrasts, "contrast_a")
  expect_null(resolved$gprofilerResults)
  expect_null(resolved$clusterprofilerResults)
  expect_null(resolved$stringdbResults)
})

test_that("resolveProtEnrichSelectedDaResults uses friendly-name mapping before fallback strategies", {
  contrasts_tbl <- data.frame(
    contrasts = c("contrast_a", "contrast_b"),
    friendly_names = c("A vs B", "C vs D"),
    stringsAsFactors = FALSE
  )
  da_results_data <- list(
    contrast_a = list(result = "mapped"),
    contrast_b = list(result = "other")
  )

  resolved <- resolveProtEnrichSelectedDaResults(
    selectedContrast = "A vs B",
    daResultsData = da_results_data,
    contrastsTbl = contrasts_tbl
  )

  expect_identical(resolved$selectedDaResults, da_results_data$contrast_a)
  expect_equal(resolved$rawContrastName, "contrast_a")
  expect_equal(resolved$source, "friendly_name")
  expect_equal(resolved$mappedRawContrastName, "contrast_a")
  expect_equal(resolved$availableKeys, c("contrast_a", "contrast_b"))
})

test_that("resolveProtEnrichSelectedDaResults falls back to fuzzy key matching", {
  da_results_data <- list(
    analysis_GA_Elevated__GA_Control = list(result = "fuzzy"),
    other_key = list(result = "other")
  )

  resolved <- resolveProtEnrichSelectedDaResults(
    selectedContrast = "GA_Elevated_vs_GA_Control",
    daResultsData = da_results_data
  )

  expect_identical(resolved$selectedDaResults, da_results_data$analysis_GA_Elevated__GA_Control)
  expect_equal(resolved$rawContrastName, "analysis_GA_Elevated__GA_Control")
  expect_equal(resolved$source, "fuzzy_match")
  expect_null(resolved$mappedRawContrastName)
})

test_that("resolveProtEnrichSelectedDaResults falls back to direct key lookup", {
  da_results_data <- list(
    direct_key = list(result = "direct")
  )

  resolved <- resolveProtEnrichSelectedDaResults(
    selectedContrast = "direct_key",
    daResultsData = da_results_data
  )

  expect_identical(resolved$selectedDaResults, da_results_data$direct_key)
  expect_equal(resolved$rawContrastName, "direct_key")
  expect_equal(resolved$source, "direct_key")
  expect_null(resolved$mappedRawContrastName)
})

test_that("resolveProtEnrichSelectedDaResults reports missing selections with available keys", {
  contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )
  da_results_data <- list(
    other_contrast = list(result = "present")
  )

  resolved <- resolveProtEnrichSelectedDaResults(
    selectedContrast = "A vs B",
    daResultsData = da_results_data,
    contrastsTbl = contrasts_tbl
  )

  expect_null(resolved$selectedDaResults)
  expect_equal(resolved$rawContrastName, "contrast_a")
  expect_null(resolved$source)
  expect_equal(resolved$mappedRawContrastName, "contrast_a")
  expect_equal(resolved$availableKeys, "other_contrast")
})

test_that("resolveProtEnrichRunDependencies keeps the current S4 object and reads design_matrix from it", {
  design_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
  current_s4_object <- methods::new(
    "mockProtEnrichDesignCarrier",
    design_matrix = design_matrix
  )
  workflow_data <- list(state_manager = NULL)
  global_env <- new.env(parent = emptyenv())

  resolved <- resolveProtEnrichRunDependencies(
    currentS4Object = current_s4_object,
    daResultsData = list(),
    workflowData = workflow_data,
    contrastsTbl = data.frame(contrasts = "contrast_a", stringsAsFactors = FALSE),
    globalEnv = global_env
  )

  expect_identical(resolved$currentS4, current_s4_object)
  expect_equal(resolved$s4Source, "current_s4_object")
  expect_identical(resolved$designMatrix, design_matrix)
  expect_equal(resolved$designMatrixSource, "s4_object")
  expect_null(resolved$designMatrixError)
  expect_equal(resolved$contrastsTbl$contrasts, "contrast_a")
})

test_that("resolveProtEnrichRunDependencies prefers the first DE result object before state manager fallback", {
  first_result_object <- methods::new(
    "mockProtEnrichDesignCarrier",
    design_matrix = matrix(c(1, 1), nrow = 1)
  )
  state_manager_object <- methods::new(
    "mockProtEnrichDesignCarrier",
    design_matrix = matrix(c(2, 2), nrow = 1)
  )
  workflow_data <- list(
    state_manager = list(
      current_state = "normalized",
      getState = function(state_name) {
        expect_equal(state_name, "normalized")
        state_manager_object
      }
    )
  )
  da_results_data <- list(
    contrast_a = list(theObject = first_result_object)
  )

  resolved <- resolveProtEnrichRunDependencies(
    currentS4Object = NULL,
    daResultsData = da_results_data,
    workflowData = workflow_data,
    globalEnv = new.env(parent = emptyenv())
  )

  expect_identical(resolved$currentS4, first_result_object)
  expect_equal(resolved$s4Source, "da_results_first_result")
  expect_identical(resolved$designMatrix, first_result_object@design_matrix)
  expect_equal(resolved$designMatrixSource, "s4_object")
})

test_that("resolveProtEnrichRunDependencies falls back to global design_matrix when slot access fails", {
  global_env <- new.env(parent = emptyenv())
  global_env$design_matrix <- matrix(c(3, 4), nrow = 1)

  resolved <- resolveProtEnrichRunDependencies(
    currentS4Object = structure(list(), class = "not_s4"),
    daResultsData = list(),
    workflowData = list(state_manager = NULL),
    globalEnv = global_env
  )

  expect_equal(resolved$s4Source, "current_s4_object")
  expect_identical(resolved$designMatrix, global_env$design_matrix)
  expect_equal(resolved$designMatrixSource, "global_environment")
  expect_false(is.null(resolved$designMatrixError))
})

test_that("resolveProtEnrichOutputDirectories prefers experiment path directories when present", {
  created_dirs <- character()

  resolved <- resolveProtEnrichOutputDirectories(
    experimentPaths = list(
      da_output_dir = "/tmp/da-existing",
      pathway_dir = "/tmp/pathway-existing",
      results_dir = "/tmp/results"
    ),
    dirExistsFn = function(path) {
      path %in% c("/tmp/da-existing", "/tmp/pathway-existing")
    },
    dirCreateFn = function(path, recursive = FALSE) {
      created_dirs <<- c(created_dirs, paste(path, recursive, sep = "|"))
      TRUE
    }
  )

  expect_equal(resolved$daOutputDir, "/tmp/da-existing")
  expect_equal(resolved$daOutputDirSource, "experiment_paths")
  expect_false(resolved$daOutputDirCreated)
  expect_equal(resolved$pathwayDir, "/tmp/pathway-existing")
  expect_equal(resolved$pathwayDirSource, "experiment_paths")
  expect_false(resolved$pathwayDirCreated)
  expect_identical(created_dirs, character())
})

test_that("resolveProtEnrichOutputDirectories creates result fallbacks when experiment paths are unavailable", {
  created_dirs <- character()

  resolved <- resolveProtEnrichOutputDirectories(
    experimentPaths = list(
      da_output_dir = "/tmp/da-missing",
      pathway_dir = "/tmp/pathway-missing",
      results_dir = "/tmp/results-root"
    ),
    dirExistsFn = function(path) {
      FALSE
    },
    dirCreateFn = function(path, recursive = FALSE) {
      created_dirs <<- c(created_dirs, paste(path, recursive, sep = "|"))
      TRUE
    }
  )

  expect_equal(resolved$daOutputDir, file.path("/tmp/results-root", "da_proteins"))
  expect_equal(resolved$daOutputDirSource, "results_fallback")
  expect_true(resolved$daOutputDirCreated)
  expect_equal(resolved$pathwayDir, file.path("/tmp/results-root", "pathway_enrichment"))
  expect_equal(resolved$pathwayDirSource, "results_fallback")
  expect_true(resolved$pathwayDirCreated)
  expect_equal(
    created_dirs,
    c(
      paste(file.path("/tmp/results-root", "da_proteins"), "TRUE", sep = "|"),
      paste(file.path("/tmp/results-root", "pathway_enrichment"), "TRUE", sep = "|")
    )
  )
})

test_that("resolveProtEnrichUniprotAnnotations prefers the global cache when available", {
  global_env <- new.env(parent = emptyenv())
  global_env$uniprot_dat_cln <- data.frame(Entry = "P1", stringsAsFactors = FALSE)
  workflow_data <- new.env(parent = emptyenv())

  resolved <- resolveProtEnrichUniprotAnnotations(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    globalEnv = global_env,
    catFn = function(...) NULL
  )

  expect_identical(resolved$uniprotDatCln, global_env$uniprot_dat_cln)
  expect_equal(resolved$source, "global_environment")
  expect_null(resolved$sourcePath)
  expect_null(resolved$cacheDir)
  expect_false(resolved$cacheDirCreated)
  expect_null(resolved$loadError)
  expect_null(resolved$creationError)
  expect_false(exists("uniprot_dat_cln", envir = workflow_data, inherits = FALSE))
})

test_that("resolveProtEnrichUniprotAnnotations loads the source-directory sidecar and syncs caches", {
  global_env <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  source_tbl <- data.frame(Entry = c("P1", "P2"), stringsAsFactors = FALSE)
  read_paths <- character()

  resolved <- resolveProtEnrichUniprotAnnotations(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    globalEnv = global_env,
    fileExistsFn = function(path) {
      identical(path, file.path("/tmp/source", "uniprot_dat_cln.RDS"))
    },
    readRdsFn = function(path) {
      read_paths <<- c(read_paths, path)
      source_tbl
    },
    catFn = function(...) NULL
  )

  expect_identical(resolved$uniprotDatCln, source_tbl)
  expect_equal(resolved$source, "source_directory")
  expect_equal(resolved$sourcePath, file.path("/tmp/source", "uniprot_dat_cln.RDS"))
  expect_null(resolved$cacheDir)
  expect_false(resolved$cacheDirCreated)
  expect_null(resolved$loadError)
  expect_null(resolved$creationError)
  expect_equal(read_paths, file.path("/tmp/source", "uniprot_dat_cln.RDS"))
  expect_identical(workflow_data$uniprot_dat_cln, source_tbl)
  expect_identical(global_env$uniprot_dat_cln, source_tbl)
})

test_that("resolveProtEnrichUniprotAnnotations generates annotations into the results cache on fallback", {
  global_env <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  current_s4_object <- methods::new(
    "mockProtEnrichProteinCarrier",
    protein_quant_table = data.frame(accession = c("P1", "P2"), stringsAsFactors = FALSE)
  )
  created_dirs <- character()
  generated_tbl <- data.frame(Entry = c("P1", "P2"), stringsAsFactors = FALSE)
  annotation_calls <- list()

  resolved <- resolveProtEnrichUniprotAnnotations(
    workflowData = workflow_data,
    experimentPaths = list(source_dir = "/tmp/source", results_dir = "/tmp/results"),
    currentS4Object = current_s4_object,
    organismTaxid = "9606",
    globalEnv = global_env,
    fileExistsFn = function(path) FALSE,
    dirExistsFn = function(path) FALSE,
    dirCreateFn = function(path, recursive = FALSE) {
      created_dirs <<- c(created_dirs, paste(path, recursive, sep = "|"))
      TRUE
    },
    getUniprotAnnotationsFn = function(input_tbl, cache_dir, taxon_id) {
      annotation_calls <<- list(
        input_tbl = input_tbl,
        cache_dir = cache_dir,
        taxon_id = taxon_id
      )
      generated_tbl
    },
    catFn = function(...) NULL
  )

  expect_identical(resolved$uniprotDatCln, generated_tbl)
  expect_equal(resolved$source, "generated")
  expect_equal(resolved$sourcePath, file.path("/tmp/source", "uniprot_dat_cln.RDS"))
  expect_equal(resolved$cacheDir, file.path("/tmp/results", "cache", "uniprot_annotations"))
  expect_true(resolved$cacheDirCreated)
  expect_null(resolved$loadError)
  expect_null(resolved$creationError)
  expect_equal(
    created_dirs,
    paste(file.path("/tmp/results", "cache", "uniprot_annotations"), "TRUE", sep = "|")
  )
  expect_identical(annotation_calls$input_tbl, current_s4_object@protein_quant_table)
  expect_equal(annotation_calls$cache_dir, file.path("/tmp/results", "cache", "uniprot_annotations"))
  expect_equal(annotation_calls$taxon_id, 9606)
  expect_identical(workflow_data$uniprot_dat_cln, generated_tbl)
  expect_identical(global_env$uniprot_dat_cln, generated_tbl)
})

test_that("resolveProtEnrichAnnotationMatching runs matching with the S4 protein id column", {
  logs <- character()
  annotation_calls <- list()
  match_results <- list(match_statistics = list(match_rate = 87))
  current_s4_object <- methods::new(
    "mockProtEnrichAnnotationCarrier",
    protein_id_column = "primary_id"
  )
  da_results_for_enrichment <- structure(list(label = "da"), class = "mock_da_results_s4")
  uniprot_tbl <- data.frame(Entry = "P1", gene_names = "GENE1", stringsAsFactors = FALSE)

  resolved <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = uniprot_tbl,
    daResultsForEnrichment = da_results_for_enrichment,
    currentS4Object = current_s4_object,
    matchAnnotationsFn = function(da_results_s4,
                                  uniprot_annotations,
                                  protein_id_column,
                                  uniprot_id_column,
                                  gene_names_column) {
      annotation_calls <<- list(
        da_results_s4 = da_results_s4,
        uniprot_annotations = uniprot_annotations,
        protein_id_column = protein_id_column,
        uniprot_id_column = uniprot_id_column,
        gene_names_column = gene_names_column
      )
      match_results
    },
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_true(resolved$attempted)
  expect_equal(resolved$proteinIdColumn, "primary_id")
  expect_identical(resolved$annotationMatchResults, match_results)
  expect_equal(resolved$matchRate, 87)
  expect_null(resolved$warning)
  expect_identical(annotation_calls$da_results_s4, da_results_for_enrichment)
  expect_identical(annotation_calls$uniprot_annotations, uniprot_tbl)
  expect_equal(annotation_calls$protein_id_column, "primary_id")
  expect_equal(annotation_calls$uniprot_id_column, "Entry")
  expect_equal(annotation_calls$gene_names_column, "gene_names")
  expect_true(any(grepl("Attempting UniProt annotation matching", logs, fixed = TRUE)))
  expect_true(any(grepl("Annotation matching completed - 87% match rate", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichAnnotationMatching falls back to uniprot_acc and surfaces warnings", {
  logs <- character()
  current_s4_object <- methods::new(
    "mockProtEnrichProteinCarrier",
    protein_quant_table = data.frame(accession = "P1", stringsAsFactors = FALSE)
  )

  resolved <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = data.frame(Entry = "P1", stringsAsFactors = FALSE),
    daResultsForEnrichment = structure(list(label = "da"), class = "mock_da_results_s4"),
    currentS4Object = current_s4_object,
    matchAnnotationsFn = function(protein_id_column, ...) {
      expect_equal(protein_id_column, "uniprot_acc")
      stop("match failure")
    },
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_true(resolved$attempted)
  expect_equal(resolved$proteinIdColumn, "uniprot_acc")
  expect_null(resolved$annotationMatchResults)
  expect_null(resolved$matchRate)
  expect_equal(resolved$warning, "match failure")
  expect_true(any(grepl("Warning in annotation matching: match failure", logs, fixed = TRUE)))
  expect_true(any(grepl("Continuing with enrichment analysis", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichAnnotationMatching skips matching when required inputs are missing", {
  resolved <- resolveProtEnrichAnnotationMatching(
    uniprotDatCln = NULL,
    daResultsForEnrichment = structure(list(label = "da"), class = "mock_da_results_s4"),
    currentS4Object = methods::new(
      "mockProtEnrichAnnotationCarrier",
      protein_id_column = "primary_id"
    ),
    matchAnnotationsFn = function(...) {
      stop("should not be called")
    }
  )

  expect_false(resolved$attempted)
  expect_null(resolved$proteinIdColumn)
  expect_null(resolved$annotationMatchResults)
  expect_null(resolved$matchRate)
  expect_null(resolved$warning)
})

test_that("resolveProtEnrichOrganismMapping prefers workflow-data organism mappings", {
  logs <- character()
  organism_mapping <- data.frame(
    uniprot_acc = c("P1", "P2"),
    taxon_id = c("9606", "10090"),
    stringsAsFactors = FALSE
  )

  resolved <- resolveProtEnrichOrganismMapping(
    workflowData = list(
      mixed_species_analysis = list(
        organism_mapping = organism_mapping
      )
    ),
    uniprotDatCln = NULL,
    targetTaxon = "9606",
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_identical(resolved$organismMapping, organism_mapping)
  expect_equal(resolved$source, "workflow_data")
  expect_null(resolved$accessionColumn)
  expect_null(resolved$taxonColumn)
  expect_equal(resolved$availableTaxonIds, c("9606", "10090"))
  expect_null(resolved$warning)
  expect_true(any(grepl("Using organism_mapping from import module", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichOrganismMapping maps organism names from UniProt annotations", {
  logs <- character()
  uniprot_tbl <- data.frame(
    Entry = c("P1", "P2", "P3"),
    Organism = c("Homo sapiens (Human)", "Mus musculus (Mouse)", NA),
    stringsAsFactors = FALSE
  )

  resolved <- resolveProtEnrichOrganismMapping(
    workflowData = list(),
    uniprotDatCln = uniprot_tbl,
    targetTaxon = "9606",
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_equal(resolved$source, "organism_names")
  expect_equal(resolved$accessionColumn, "Entry")
  expect_null(resolved$taxonColumn)
  expect_equal(resolved$availableTaxonIds, c("9606", "10090"))
  expect_null(resolved$warning)
  expect_equal(
    resolved$organismMapping,
    data.frame(
      uniprot_acc = c("P1", "P2", "P3"),
      taxon_id = c("9606", "10090", NA),
      stringsAsFactors = FALSE
    )
  )
  expect_true(any(grepl("Checking uniprot_dat_cln columns: Entry, Organism", logs, fixed = TRUE)))
  expect_true(any(grepl("Found 'Organism' column with organism names", logs, fixed = TRUE)))
  expect_true(any(grepl("Found 1 proteins matching target taxon 9606", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichOrganismMapping parses numeric taxon columns from UniProt annotations", {
  logs <- character()
  uniprot_tbl <- data.frame(
    UniProt_Acc = c("P1", "P2"),
    OX = c("ID=9606", "10090"),
    stringsAsFactors = FALSE
  )

  resolved <- resolveProtEnrichOrganismMapping(
    workflowData = list(),
    uniprotDatCln = uniprot_tbl,
    targetTaxon = "9606",
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_equal(resolved$source, "taxon_column")
  expect_equal(resolved$accessionColumn, "UniProt_Acc")
  expect_equal(resolved$taxonColumn, "OX")
  expect_equal(resolved$availableTaxonIds, c("9606", "10090"))
  expect_null(resolved$warning)
  expect_equal(
    resolved$organismMapping,
    data.frame(
      uniprot_acc = c("P1", "P2"),
      taxon_id = c("9606", "10090"),
      stringsAsFactors = FALSE
    )
  )
  expect_true(any(grepl("Found taxon column: OX", logs, fixed = TRUE)))
  expect_true(any(grepl("Created organism_mapping from taxon column (2 entries)", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichOrganismMapping falls back to single-species mapping", {
  logs <- character()
  uniprot_tbl <- data.frame(
    accession = c("P1", "P2"),
    stringsAsFactors = FALSE
  )

  resolved <- resolveProtEnrichOrganismMapping(
    workflowData = list(),
    uniprotDatCln = uniprot_tbl,
    targetTaxon = "9606",
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_equal(resolved$source, "single_species_fallback")
  expect_equal(resolved$accessionColumn, "accession")
  expect_null(resolved$taxonColumn)
  expect_equal(resolved$availableTaxonIds, "9606")
  expect_null(resolved$warning)
  expect_equal(
    resolved$organismMapping,
    data.frame(
      uniprot_acc = c("P1", "P2"),
      taxon_id = c("9606", "9606"),
      stringsAsFactors = FALSE
    )
  )
  expect_true(any(grepl("No organism column found - creating single-species mapping", logs, fixed = TRUE)))
  expect_true(any(grepl("Created single-species organism_mapping (2 entries, all assigned to taxon 9606)", logs, fixed = TRUE)))
})

test_that("resolveProtEnrichOrganismMapping returns an empty config when no mapping source is available", {
  resolved <- resolveProtEnrichOrganismMapping(
    workflowData = list(),
    uniprotDatCln = NULL,
    targetTaxon = "9606",
    catFn = function(...) NULL
  )

  expect_null(resolved$organismMapping)
  expect_null(resolved$source)
  expect_null(resolved$accessionColumn)
  expect_null(resolved$taxonColumn)
  expect_identical(resolved$availableTaxonIds, character())
  expect_null(resolved$warning)
})

test_that("applyProtEnrichOrganismFilter filters each contrast to the target organism", {
  logs <- character()
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = c("contrast_a", "contrast_b")),
    da_data = list(
      contrast_a = data.frame(
        protein_group = c("P1;Q1", "P3;Q3", "P5"),
        stringsAsFactors = FALSE
      ),
      contrast_b = data.frame(
        protein_group = c("P2-2;Q2", "P4"),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )
  organism_mapping <- data.frame(
    uniprot_acc = c("P1", "P2-2", "P9"),
    taxon_id = c("9606", "9606", "10090"),
    stringsAsFactors = FALSE
  )
  current_s4_object <- methods::new(
    "mockProtEnrichAnnotationCarrier",
    protein_id_column = "protein_group"
  )

  resolved <- applyProtEnrichOrganismFilter(
    daResultsForEnrichment = da_results,
    organismMapping = organism_mapping,
    targetTaxon = "9606",
    currentS4Object = current_s4_object,
    normalizeUniprotFn = function(accessions, remove_isoform = TRUE) {
      expect_true(isTRUE(remove_isoform))
      sub("-.*$", "", accessions)
    },
    cleanAccFn = function(accessions) sub("-.*$", "", accessions),
    catFn = function(...) {
      logs <<- c(logs, paste0(..., collapse = ""))
    }
  )

  expect_true(isTRUE(resolved$filterApplied))
  expect_equal(resolved$proteinIdColumn, "protein_group")
  expect_equal(resolved$targetProteins, c("P1", "P2-2"))
  expect_equal(resolved$targetProteinsClean, c("P1", "P2"))
  expect_equal(resolved$filterStats$proteins_before, 5)
  expect_equal(resolved$filterStats$proteins_after, 2)
  expect_equal(resolved$filterStats$proteins_removed, 3)
  expect_equal(
    resolved$daResultsForEnrichment@da_data$contrast_a$protein_group,
    "P1;Q1"
  )
  expect_equal(
    resolved$daResultsForEnrichment@da_data$contrast_b$protein_group,
    "P2-2;Q2"
  )
  expect_true(any(grepl("Found 2 proteins for target taxon 9606", logs, fixed = TRUE)))
  expect_true(any(grepl("Organism filtering complete - kept 2/5 proteins (40.0%)", logs, fixed = TRUE)))
})

test_that("applyProtEnrichOrganismFilter falls back to uniprot_acc when the S4 id slot is unavailable", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = c("P1", "P2"),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )
  current_s4_object <- methods::new(
    "mockProtEnrichProteinCarrier",
    protein_quant_table = data.frame(accession = "P1", stringsAsFactors = FALSE)
  )

  resolved <- applyProtEnrichOrganismFilter(
    daResultsForEnrichment = da_results,
    organismMapping = data.frame(
      uniprot_acc = "P1",
      taxon_id = "9606",
      stringsAsFactors = FALSE
    ),
    targetTaxon = "9606",
    currentS4Object = current_s4_object,
    normalizeUniprotFn = function(accessions, remove_isoform = TRUE) accessions,
    cleanAccFn = function(accessions) accessions,
    catFn = function(...) NULL
  )

  expect_true(isTRUE(resolved$filterApplied))
  expect_equal(resolved$proteinIdColumn, "uniprot_acc")
  expect_equal(resolved$filterStats$proteins_before, 2)
  expect_equal(resolved$filterStats$proteins_after, 1)
  expect_equal(resolved$filterStats$proteins_removed, 1)
  expect_equal(resolved$daResultsForEnrichment@da_data$contrast_a$uniprot_acc, "P1")
})

test_that("persistProtEnrichOrganismFilterMetadata stores assembled reporting metadata", {
  workflow_data <- new.env(parent = emptyenv())
  timestamp <- as.POSIXct("2026-04-14 10:15:00", tz = "UTC")

  resolved <- persistProtEnrichOrganismFilterMetadata(
    workflowData = workflow_data,
    organismFilterEnabled = TRUE,
    organismFilterApplied = TRUE,
    targetTaxonId = "9606",
    filterStats = list(
      proteins_before = 42,
      proteins_after = 30,
      proteins_removed = 12
    ),
    timeFn = function() timestamp
  )

  expect_true(is.list(workflow_data$enrichment_organism_filter))
  expect_identical(resolved, workflow_data$enrichment_organism_filter)
  expect_true(isTRUE(resolved$enabled))
  expect_true(isTRUE(resolved$filter_applied))
  expect_equal(resolved$target_taxon_id, "9606")
  expect_equal(resolved$proteins_before, 42)
  expect_equal(resolved$proteins_after, 30)
  expect_equal(resolved$proteins_removed, 12)
  expect_identical(resolved$timestamp, timestamp)
})

test_that("prepareProtEnrichAnalysisBodySetup resolves the pre-process runner context", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$da_results_data <- list(contrast_a = list(result = "mapped"))
  enrichment_data$current_s4_object <- NULL
  enrichment_data$annotation_match_results <- NULL

  workflow_data <- new.env(parent = emptyenv())

  input <- list(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    enable_organism_filter = FALSE
  )
  global_env <- new.env(parent = emptyenv())
  global_env$contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$order <- character()
  recorded$logs <- character()

  prepared <- prepareProtEnrichAnalysisBodySetup(
    selectedContrast = "A vs B",
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = list(
      da_output_dir = "/tmp/da",
      pathway_dir = "/tmp/pathway",
      results_dir = "/tmp/results"
    ),
    resolveSelectedDaResultsFn = function(selectedContrast, daResultsData, contrastsTbl) {
      recorded$order <- c(recorded$order, "selected")
      expect_equal(selectedContrast, "A vs B")
      expect_identical(daResultsData, enrichment_data$da_results_data)
      expect_identical(contrastsTbl, global_env$contrasts_tbl)
      list(
        selectedDaResults = daResultsData$contrast_a,
        rawContrastName = "contrast_a",
        source = "friendly_name",
        mappedRawContrastName = "contrast_a",
        availableKeys = "contrast_a"
      )
    },
    resolveRunDependenciesFn = function(currentS4Object, daResultsData, workflowData, contrastsTbl) {
      recorded$order <- c(recorded$order, "dependencies")
      expect_null(currentS4Object)
      expect_identical(daResultsData, enrichment_data$da_results_data)
      expect_identical(workflowData, workflow_data)
      expect_identical(contrastsTbl, global_env$contrasts_tbl)
      list(
        contrastsTbl = contrastsTbl,
        designMatrix = matrix(1, nrow = 1),
        currentS4 = "current_s4",
        s4Source = "state_manager",
        designMatrixSource = "global_environment",
        designMatrixError = NULL
      )
    },
    resolveOutputDirectoriesFn = function(experimentPaths) {
      recorded$order <- c(recorded$order, "paths")
      expect_equal(experimentPaths$da_output_dir, "/tmp/da")
      expect_equal(experimentPaths$pathway_dir, "/tmp/pathway")
      list(
        daOutputDir = "/tmp/da",
        pathwayDir = "/tmp/pathway",
        daOutputDirSource = "experiment_paths",
        pathwayDirSource = "experiment_paths"
      )
    },
    createDaResultsForEnrichmentFn = function(contrasts_tbl, design_matrix, da_output_dir) {
      recorded$order <- c(recorded$order, "create")
      expect_identical(contrasts_tbl, global_env$contrasts_tbl)
      expect_equal(design_matrix[1, 1], 1)
      expect_equal(da_output_dir, "/tmp/da")
      structure(list(kind = "da_results_s4"), class = "mock_da_results_s4")
    },
    resolveUniprotAnnotationsFn = function(workflowData, experimentPaths, currentS4Object, organismTaxid) {
      recorded$order <- c(recorded$order, "annotations")
      expect_identical(workflowData, workflow_data)
      expect_equal(experimentPaths$results_dir, "/tmp/results")
      expect_equal(currentS4Object, "current_s4")
      expect_equal(organismTaxid, "9606")
      list(uniprotDatCln = data.frame(Entry = "P1", stringsAsFactors = FALSE))
    },
    resolveAnnotationMatchingFn = function(uniprotDatCln, daResultsForEnrichment, currentS4Object) {
      recorded$order <- c(recorded$order, "matching")
      expect_equal(uniprotDatCln$Entry, "P1")
      expect_equal(daResultsForEnrichment$kind, "da_results_s4")
      expect_equal(currentS4Object, "current_s4")
      list(annotationMatchResults = list(match_rate = 100))
    },
    resolveOrganismMappingFn = function(...) {
      stop("organism mapping should not be called when filtering is disabled")
    },
    applyOrganismFilterFn = function(...) {
      stop("organism filter should not be called when filtering is disabled")
    },
    persistOrganismFilterMetadataFn = function(workflowData,
                                               organismFilterEnabled,
                                               organismFilterApplied,
                                               targetTaxonId,
                                               filterStats) {
      recorded$order <- c(recorded$order, "persist")
      workflowData$enrichment_organism_filter <- list(
        enabled = organismFilterEnabled,
        applied = organismFilterApplied,
        target = targetTaxonId,
        stats = filterStats
      )
    },
    showNotificationFn = function(...) {
      stop("warning notification should not fire on the happy path")
    },
    globalEnv = global_env,
    catFn = function(message) {
      recorded$logs <- c(recorded$logs, message)
    }
  )

  expect_identical(
    recorded$order,
    c("selected", "dependencies", "paths", "create", "annotations", "matching", "persist")
  )
  expect_equal(prepared$rawContrastName, "contrast_a")
  expect_identical(prepared$contrastsTbl, global_env$contrasts_tbl)
  expect_equal(prepared$pathwayDir, "/tmp/pathway")
  expect_equal(prepared$daResultsForEnrichment$kind, "da_results_s4")
  expect_equal(prepared$goAnnotations$Entry, "P1")
  expect_false(prepared$organismFilterApplied)
  expect_equal(prepared$filterStats$proteins_before, 0)
  expect_equal(prepared$filterStats$proteins_after, 0)
  expect_equal(prepared$filterStats$proteins_removed, 0)
  expect_equal(enrichment_data$current_s4_object, "current_s4")
  expect_identical(enrichment_data$annotation_match_results, list(match_rate = 100))
  expect_false(workflow_data$enrichment_organism_filter$enabled)
  expect_false(workflow_data$enrichment_organism_filter$applied)
  expect_equal(workflow_data$enrichment_organism_filter$target, "9606")
})

test_that("resolveProtEnrichAnalysisInputColumns keeps the S4 protein id column by default", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = "P1",
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )
  current_s4_object <- methods::new(
    "mockProtEnrichAnnotationCarrier",
    protein_id_column = "primary_id"
  )

  resolved <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = list(method = "clusterprofiler"),
    daResultsForEnrichment = da_results,
    currentS4Object = current_s4_object
  )

  expect_equal(resolved$idColumn, "primary_id")
  expect_equal(resolved$source, "s4_object")
  expect_false(resolved$geneNameOverrideApplied)
})

test_that("resolveProtEnrichAnalysisInputColumns falls back to uniprot_acc without an S4 object", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = "P1",
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )

  resolved <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = list(method = "clusterprofiler"),
    daResultsForEnrichment = da_results,
    currentS4Object = NULL
  )

  expect_equal(resolved$idColumn, "uniprot_acc")
  expect_equal(resolved$source, "default")
  expect_false(resolved$geneNameOverrideApplied)
})

test_that("resolveProtEnrichAnalysisInputColumns prefers gene_name for gprofiler2 inputs", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        primary_id = "P1",
        gene_name = "GENE1",
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )
  current_s4_object <- methods::new(
    "mockProtEnrichAnnotationCarrier",
    protein_id_column = "primary_id"
  )

  resolved <- resolveProtEnrichAnalysisInputColumns(
    methodInfo = list(method = "gprofiler2"),
    daResultsForEnrichment = da_results,
    currentS4Object = current_s4_object
  )

  expect_equal(resolved$idColumn, "gene_name")
  expect_equal(resolved$source, "gprofiler_gene_name_override")
  expect_true(resolved$geneNameOverrideApplied)
})

test_that("buildProtEnrichProcessEnrichmentsArgs reuses one shared argument bundle for CP10 and processEnrichments", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = "P1",
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )
  go_annotations <- data.frame(
    Entry = "P1",
    stringsAsFactors = FALSE
  )

  bundled_args <- buildProtEnrichProcessEnrichmentsArgs(
    daResultsForEnrichment = da_results,
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    pathwayDir = "/tmp/pathways",
    goAnnotations = go_annotations,
    proteinIdColumn = "gene_name",
    contrastNames = c("contrast_a", "contrast_b"),
    correctionMethod = "BH"
  )

  expect_identical(bundled_args$checkpointArgs$da_results_s4, da_results)
  expect_identical(bundled_args$processArgs$da_results, da_results)
  expect_equal(bundled_args$checkpointArgs$taxon_id, 9606)
  expect_equal(bundled_args$processArgs$taxon_id, 9606)
  expect_identical(bundled_args$checkpointArgs$pathway_dir, "/tmp/pathways")
  expect_identical(bundled_args$processArgs$pathway_dir, "/tmp/pathways")
  expect_identical(bundled_args$checkpointArgs$go_annotations, go_annotations)
  expect_identical(bundled_args$processArgs$go_annotations, go_annotations)
  expect_identical(bundled_args$checkpointArgs$protein_id_column, "gene_name")
  expect_identical(bundled_args$processArgs$protein_id_column, "gene_name")
  expect_identical(bundled_args$checkpointArgs$contrast_names, c("contrast_a", "contrast_b"))
  expect_identical(bundled_args$processArgs$contrast_names, c("contrast_a", "contrast_b"))
  expect_false(bundled_args$checkpointArgs$exclude_iea)
  expect_false(bundled_args$processArgs$exclude_iea)
})

test_that("buildProtEnrichProcessEnrichmentsArgs allows explicit exclude_iea overrides", {
  da_results <- methods::new(
    "da_results_for_enrichment",
    contrasts = tibble::tibble(contrasts = "contrast_a"),
    da_data = list(
      contrast_a = data.frame(
        uniprot_acc = "P1",
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(sample = character(), stringsAsFactors = FALSE)
  )

  bundled_args <- buildProtEnrichProcessEnrichmentsArgs(
    daResultsForEnrichment = da_results,
    organismTaxid = 10090,
    upCutoff = 0,
    downCutoff = 0,
    qCutoff = 0.1,
    pathwayDir = "/tmp/pathways",
    goAnnotations = NULL,
    proteinIdColumn = "uniprot_acc",
    contrastNames = "contrast_a",
    correctionMethod = "gSCS",
    excludeIea = TRUE
  )

  expect_true(bundled_args$checkpointArgs$exclude_iea)
  expect_true(bundled_args$processArgs$exclude_iea)
  expect_equal(bundled_args$checkpointArgs$taxon_id, 10090)
  expect_equal(bundled_args$processArgs$taxon_id, 10090)
})

test_that("prepareProtEnrichProcessExecution resolves method, input columns, and shared process args", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- "current_s4"
  enrichment_data$da_results_data <- list(contrast_a = list(), contrast_b = list())

  input <- list(
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05,
    correction_method = "fdr"
  )
  da_results_for_enrichment <- structure(list(kind = "da_results_s4"), class = "mock_da_results_s4")
  go_annotations <- data.frame(Entry = "P1", stringsAsFactors = FALSE)

  recorded <- new.env(parent = emptyenv())
  recorded$order <- character()
  recorded$logs <- character()

  prepared <- prepareProtEnrichProcessExecution(
    input = input,
    enrichmentData = enrichment_data,
    daResultsForEnrichment = da_results_for_enrichment,
    pathwayDir = "/tmp/pathway",
    goAnnotations = go_annotations,
    currentAnalysisMethodFn = function() {
      recorded$order <- c(recorded$order, "method")
      list(method = "gprofiler2", supported = TRUE, species_name = "Human")
    },
    resolveAnalysisInputColumnsFn = function(methodInfo, daResultsForEnrichment, currentS4Object) {
      recorded$order <- c(recorded$order, "inputs")
      expect_equal(methodInfo$method, "gprofiler2")
      expect_identical(daResultsForEnrichment, da_results_for_enrichment)
      expect_identical(currentS4Object, "current_s4")
      list(idColumn = "gene_name", source = "override", geneNameOverrideApplied = TRUE)
    },
    buildProcessEnrichmentsArgsFn = function(daResultsForEnrichment,
                                             organismTaxid,
                                             upCutoff,
                                             downCutoff,
                                             qCutoff,
                                             pathwayDir,
                                             goAnnotations,
                                             proteinIdColumn,
                                             contrastNames,
                                             correctionMethod) {
      recorded$order <- c(recorded$order, "args")
      expect_identical(daResultsForEnrichment, da_results_for_enrichment)
      expect_equal(organismTaxid, "9606")
      expect_equal(upCutoff, 1)
      expect_equal(downCutoff, -1)
      expect_equal(qCutoff, 0.05)
      expect_equal(pathwayDir, "/tmp/pathway")
      expect_identical(goAnnotations, go_annotations)
      expect_equal(proteinIdColumn, "gene_name")
      expect_equal(contrastNames, c("contrast_a", "contrast_b"))
      expect_equal(correctionMethod, "fdr")
      list(
        checkpointArgs = list(kind = "checkpoint"),
        processArgs = list(kind = "process")
      )
    },
    catFn = function(message) {
      recorded$logs <- c(recorded$logs, trimws(message))
    }
  )

  expect_equal(recorded$order, c("method", "inputs", "args"))
  expect_equal(recorded$logs, "ENRICHMENT Step: Using analysis method: gprofiler2")
  expect_equal(prepared$methodInfo$method, "gprofiler2")
  expect_identical(
    prepared$inputColumnConfig,
    list(idColumn = "gene_name", source = "override", geneNameOverrideApplied = TRUE)
  )
  expect_identical(prepared$enrichmentArgs$checkpointArgs$kind, "checkpoint")
  expect_identical(prepared$enrichmentArgs$processArgs$kind, "process")
})

test_that("executeProtEnrichProcessEnrichments captures CP10, dispatches processEnrichments, and logs completion", {
  recorded <- new.env(parent = emptyenv())
  recorded$checkpoint <- NULL
  recorded$process <- NULL
  recorded$logs <- character()

  enrichment_results <- list(status = "ok")
  enrichment_args <- list(
    checkpointArgs = list(kind = "checkpoint"),
    processArgs = list(kind = "process", taxon_id = 9606)
  )

  result <- executeProtEnrichProcessEnrichments(
    enrichmentArgs = enrichment_args,
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    captureCheckpointFn = function(args, id, label) {
      recorded$checkpoint <- list(args = args, id = id, label = label)
    },
    processEnrichmentsFn = function(kind, taxon_id) {
      recorded$process <- list(kind = kind, taxon_id = taxon_id)
      enrichment_results
    },
    catFn = function(message) {
      recorded$logs <- c(recorded$logs, trimws(message))
    }
  )

  expect_identical(result, enrichment_results)
  expect_equal(recorded$checkpoint$id, "cp10")
  expect_equal(recorded$checkpoint$label, "enrichment_input")
  expect_equal(recorded$checkpoint$args$kind, "checkpoint")
  expect_identical(recorded$process, list(kind = "process", taxon_id = 9606))
  expect_equal(
    recorded$logs,
    "ENRICHMENT Step: processEnrichments completed with up_cutoff: 1.0, down_cutoff: -1.0, q_cutoff: 0.050"
  )
})

test_that("buildProtEnrichAllContrastResults collates and standardizes gprofiler2 payloads", {
  enrichment_results <- methods::new(
    "mockProtEnrichAllContrastCarrier",
    enrichment_data = list(
      contrast_a = list(
        up = list(result = data.frame(
          term_name = "up term",
          p_value = 0.001,
          term_size = 5,
          source = "GO:BP",
          stringsAsFactors = FALSE
        )),
        down = list(result = data.frame(
          term_name = "down term",
          p_value = 0.01,
          term_size = 3,
          source = "REAC",
          stringsAsFactors = FALSE
        ))
      )
    )
  )

  collated <- buildProtEnrichAllContrastResults(
    enrichmentResults = enrichment_results,
    methodInfo = list(method = "gprofiler2"),
    catFn = function(...) invisible(NULL)
  )

  contrast_results <- collated$contrast_a$gprofiler_results

  expect_equal(names(collated), "contrast_a")
  expect_equal(contrast_results$directionality, c("positive", "negative"))
  expect_true(all(contrast_results$analysis_method == "gprofiler2"))
  expect_equal(contrast_results$Description, c("up term", "down term"))
  expect_equal(contrast_results$pvalue, c(0.001, 0.01))
  expect_equal(contrast_results$p.adjust, c(0.001, 0.01))
  expect_equal(contrast_results$qvalue, c(0.001, 0.01))
  expect_equal(contrast_results$Count, c(5, 3))
  expect_equal(contrast_results$data_source, c("GO:BP", "REAC"))
  expect_null(collated$contrast_a$clusterprofiler_results)
  expect_null(collated$contrast_a$stringdb_results)
})

test_that("buildProtEnrichAllContrastResults collates clusterprofiler enrichResult payloads", {
  up_results <- methods::new(
    "mockProtEnrichClusterResultCarrier",
    result = data.frame(ID = "GO:0001", stringsAsFactors = FALSE)
  )
  down_results <- methods::new(
    "mockProtEnrichClusterResultCarrier",
    result = data.frame(ID = "GO:0002", stringsAsFactors = FALSE)
  )
  enrichment_results <- methods::new(
    "mockProtEnrichAllContrastCarrier",
    enrichment_data = list(
      contrast_a = list(
        up = up_results,
        down = down_results
      )
    )
  )

  collated <- buildProtEnrichAllContrastResults(
    enrichmentResults = enrichment_results,
    methodInfo = list(method = "clusterprofiler"),
    isEnrichResultFn = function(object, class_name) {
      expect_identical(class_name, "enrichResult")
      methods::is(object, "mockProtEnrichClusterResultCarrier")
    },
    catFn = function(...) invisible(NULL)
  )

  contrast_results <- collated$contrast_a$clusterprofiler_results

  expect_equal(names(collated), "contrast_a")
  expect_equal(contrast_results$ID, c("GO:0001", "GO:0002"))
  expect_equal(contrast_results$directionality, c("up", "down"))
  expect_true(all(contrast_results$analysis_method == "clusterprofiler"))
  expect_null(collated$contrast_a$gprofiler_results)
  expect_null(collated$contrast_a$stringdb_results)
})

test_that("buildProtEnrichAnalysisResultsPayload assembles the workflow payload from stored results", {
  gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  clusterprofiler_results <- data.frame(ID = "GO:0001", stringsAsFactors = FALSE)
  stringdb_results <- data.frame(term = "ppi", stringsAsFactors = FALSE)
  full_enrichment_results <- methods::new(
    "mockProtEnrichAllContrastCarrier",
    enrichment_data = list()
  )

  payload <- buildProtEnrichAnalysisResultsPayload(
    gprofilerResults = gprofiler_results,
    clusterprofilerResults = clusterprofiler_results,
    stringdbResults = stringdb_results,
    fullEnrichmentResults = full_enrichment_results,
    selectedContrast = "contrast_a",
    analysisMethod = "gprofiler2",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606"
  )

  expect_identical(payload$gprofiler_results, gprofiler_results)
  expect_identical(payload$clusterprofiler_results, clusterprofiler_results)
  expect_identical(payload$stringdb_results, stringdb_results)
  expect_identical(payload$full_enrichment_results, full_enrichment_results)
  expect_identical(payload$selected_contrast, "contrast_a")
  expect_identical(payload$analysis_method, "gprofiler2")
  expect_identical(payload$parameters$up_cutoff, 1)
  expect_identical(payload$parameters$down_cutoff, -1)
  expect_identical(payload$parameters$q_cutoff, 0.05)
  expect_identical(payload$parameters$organism_taxid, "9606")
})

test_that("propagateProtEnrichResultsArgs copies source args and enrichment UI metadata", {
  source_object <- methods::new(
    "mockProtEnrichArgsCarrier",
    args = list(globalParameters = list(workflow_type = "proteomics"))
  )
  enrichment_results <- methods::new(
    "mockProtEnrichArgsCarrier",
    args = list(existing = list(value = "keep"))
  )
  fixed_time <- as.POSIXct("2026-04-14 12:00:00", tz = "UTC")

  propagated <- propagateProtEnrichResultsArgs(
    enrichmentResults = enrichment_results,
    currentS4Object = source_object,
    selectedContrast = "contrast_a",
    methodInfo = list(
      method = "gprofiler2",
      supported = TRUE,
      species_name = "Homo sapiens"
    ),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathway_enrichment",
    timeFn = function() fixed_time,
    catFn = function(...) invisible(NULL)
  )

  expect_true(propagated$dataHasArgs)
  expect_true(propagated$resultsHasArgs)
  expect_true(propagated$copiedArgs)
  expect_identical(
    propagated$enrichmentResults@args$globalParameters$workflow_type,
    "proteomics"
  )
  expect_identical(
    propagated$enrichmentResults@args$enrichmentAnalysis,
    list(
      selected_contrast = "contrast_a",
      analysis_method = "gprofiler2",
      organism_supported = TRUE,
      up_cutoff = 1,
      down_cutoff = -1,
      q_cutoff = 0.05,
      organism_taxid = "9606",
      pathway_dir = "/tmp/pathway_enrichment"
    )
  )
  expect_identical(
    propagated$enrichmentResults@args$enrichmentAnalysisUI,
    list(
      up_log2fc_cutoff = 1,
      down_log2fc_cutoff = -1,
      q_value_cutoff = 0.05,
      organism_taxon_id = "9606",
      analysis_method = "gprofiler2",
      organism_name = "Homo sapiens",
      organism_supported = TRUE,
      selected_contrast = "contrast_a",
      timestamp = fixed_time
    )
  )
})

test_that("propagateProtEnrichResultsArgs leaves results untouched when source args are unavailable", {
  enrichment_results <- methods::new(
    "mockProtEnrichArgsCarrier",
    args = list(existing = list(value = "keep"))
  )

  propagated <- propagateProtEnrichResultsArgs(
    enrichmentResults = enrichment_results,
    currentS4Object = structure(list(label = "no_args"), class = "mockNoArgsCarrier"),
    selectedContrast = "contrast_a",
    methodInfo = list(
      method = "gprofiler2",
      supported = TRUE,
      species_name = "Homo sapiens"
    ),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathway_enrichment",
    catFn = function(...) invisible(NULL)
  )

  expect_false(propagated$dataHasArgs)
  expect_true(propagated$resultsHasArgs)
  expect_false(propagated$copiedArgs)
  expect_identical(
    propagated$enrichmentResults@args,
    list(existing = list(value = "keep"))
  )
})

test_that("propagateProtEnrichUiParams stores source-object and workflow UI metadata", {
  source_object <- methods::new(
    "mockProtEnrichArgsCarrier",
    args = list(globalParameters = list(workflow_type = "proteomics"))
  )
  workflow_data <- new.env(parent = emptyenv())
  fixed_times <- as.POSIXct(
    c("2026-04-14 12:00:00", "2026-04-14 12:00:01"),
    tz = "UTC"
  )
  time_index <- 0L

  propagated <- propagateProtEnrichUiParams(
    currentS4Object = source_object,
    workflowData = workflow_data,
    selectedContrast = "contrast_a",
    methodInfo = list(
      method = "gprofiler2",
      supported = TRUE,
      species_name = "Homo sapiens"
    ),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    timeFn = function() {
      time_index <<- time_index + 1L
      fixed_times[[time_index]]
    },
    catFn = function(...) invisible(NULL)
  )

  expect_true(propagated$dataHasArgs)
  expect_true(propagated$storedUiParams)
  expect_identical(
    propagated$currentS4Object@args$enrichmentAnalysisUI,
    list(
      up_log2fc_cutoff = 1,
      down_log2fc_cutoff = -1,
      q_value_cutoff = 0.05,
      organism_taxon_id = "9606",
      analysis_method = "gprofiler2",
      organism_name = "Homo sapiens",
      organism_supported = TRUE,
      selected_contrast = "contrast_a",
      timestamp = fixed_times[[1]]
    )
  )
  expect_identical(
    workflow_data$enrichment_ui_params,
    list(
      up_log2fc_cutoff = 1,
      down_log2fc_cutoff = -1,
      q_value_cutoff = 0.05,
      organism_selected = "9606",
      database_source = "gprofiler2",
      organism_name = "Homo sapiens",
      organism_supported = TRUE,
      selected_contrast = "contrast_a",
      timestamp = fixed_times[[2]]
    )
  )
  expect_identical(propagated$workflowUiParams, workflow_data$enrichment_ui_params)
})

test_that("propagateProtEnrichUiParams skips UI propagation when source args are unavailable", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$enrichment_ui_params <- list(existing = "keep")

  propagated <- propagateProtEnrichUiParams(
    currentS4Object = structure(list(label = "no_args"), class = "mockNoArgsCarrier"),
    workflowData = workflow_data,
    selectedContrast = "contrast_a",
    methodInfo = list(
      method = "gprofiler2",
      supported = TRUE,
      species_name = "Homo sapiens"
    ),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    catFn = function(...) invisible(NULL)
  )

  expect_false(propagated$dataHasArgs)
  expect_false(propagated$storedUiParams)
  expect_null(propagated$workflowUiParams)
  expect_identical(workflow_data$enrichment_ui_params, list(existing = "keep"))
})

test_that("updateProtEnrichStateManagerUiParams detects the current persisted data state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(
    getHistory = function() c("protein_replicate_filtered", "enrichment_completed")
  )

  updated <- updateProtEnrichStateManagerUiParams(
    workflowData = workflow_data,
    storedUiParams = TRUE,
    catFn = function(...) invisible(NULL)
  )

  expect_true(updated$attempted)
  expect_equal(updated$currentDataState, "protein_replicate_filtered")
  expect_equal(updated$availableStates, c("protein_replicate_filtered", "enrichment_completed"))
  expect_false(updated$updated)
  expect_null(updated$warning)
})

test_that("updateProtEnrichStateManagerUiParams skips state-manager lookup when UI params were not stored", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(
    getHistory = function() stop("getHistory should not be called")
  )

  updated <- updateProtEnrichStateManagerUiParams(
    workflowData = workflow_data,
    storedUiParams = FALSE,
    catFn = function(...) invisible(NULL)
  )

  expect_false(updated$attempted)
  expect_null(updated$currentDataState)
  expect_null(updated$availableStates)
  expect_false(updated$updated)
  expect_null(updated$warning)
})

test_that("updateProtEnrichStateManagerUiParams surfaces state-manager lookup warnings", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(
    getHistory = function() stop("history unavailable")
  )

  updated <- updateProtEnrichStateManagerUiParams(
    workflowData = workflow_data,
    storedUiParams = TRUE,
    catFn = function(...) invisible(NULL)
  )

  expect_true(updated$attempted)
  expect_null(updated$currentDataState)
  expect_null(updated$availableStates)
  expect_false(updated$updated)
  expect_equal(updated$warning, "history unavailable")
})

test_that("saveProtEnrichCompletedState stores enrichment results and config in the state manager", {
  workflow_data <- new.env(parent = emptyenv())
  saved_state <- NULL
  workflow_data$state_manager <- list(
    saveState = function(...) {
      saved_state <<- list(...)
    }
  )

  saved <- saveProtEnrichCompletedState(
    workflowData = workflow_data,
    enrichmentResults = "enrichment_results",
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    catFn = function(...) invisible(NULL)
  )

  expect_true(saved$attempted)
  expect_true(saved$saved)
  expect_null(saved$warning)
  expect_equal(saved_state$state_name, "enrichment_completed")
  expect_identical(saved_state$s4_data_object, "enrichment_results")
  expect_identical(
    saved_state$config_object,
    list(
      selected_contrast = "contrast_a",
      analysis_method = "gprofiler2",
      organism_supported = TRUE,
      up_cutoff = 1,
      down_cutoff = -1,
      q_cutoff = 0.05,
      organism_taxid = "9606",
      pathway_dir = "/tmp/pathways"
    )
  )
  expect_equal(
    saved_state$description,
    "Enrichment analysis completed using gprofiler2 for contrast: contrast_a"
  )
})

test_that("saveProtEnrichCompletedState surfaces state-manager save warnings", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(
    saveState = function(...) stop("save failed")
  )

  saved <- saveProtEnrichCompletedState(
    workflowData = workflow_data,
    enrichmentResults = "enrichment_results",
    selectedContrast = "contrast_a",
    methodInfo = list(method = "gprofiler2", supported = TRUE),
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    organismTaxid = "9606",
    pathwayDir = "/tmp/pathways",
    catFn = function(...) invisible(NULL)
  )

  expect_true(saved$attempted)
  expect_false(saved$saved)
  expect_equal(saved$warning, "save failed")
})

test_that("completeProtEnrichTabStatus marks enrichment complete and preserves peer tabs", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(
    setup_import = "complete",
    differential_expression = "pending",
    enrichment_analysis = "pending"
  )

  completed <- completeProtEnrichTabStatus(workflowData = workflow_data)

  expect_equal(completed$tabName, "enrichment_analysis")
  expect_equal(completed$status, "complete")
  expect_equal(workflow_data$tab_status$setup_import, "complete")
  expect_equal(workflow_data$tab_status$differential_expression, "pending")
  expect_equal(workflow_data$tab_status$enrichment_analysis, "complete")
  expect_identical(completed$tabStatus, workflow_data$tab_status)
})

test_that("completeProtEnrichProgress issues the final completion update", {
  recorded <- new.env(parent = emptyenv())
  recorded$value <- NULL
  recorded$detail <- NULL

  completed <- completeProtEnrichProgress(
    value = 1,
    detail = "Complete!",
    incProgressFn = function(value, detail) {
      recorded$value <- value
      recorded$detail <- detail
    }
  )

  expect_equal(recorded$value, 1)
  expect_equal(recorded$detail, "Complete!")
  expect_equal(completed$value, 1)
  expect_equal(completed$detail, "Complete!")
})

test_that("notifyProtEnrichCompletion emits the success notification payload", {
  recorded <- new.env(parent = emptyenv())
  recorded$message <- NULL
  recorded$type <- NULL
  recorded$duration <- NULL

  completed <- notifyProtEnrichCompletion(
    selectedContrast = "treated_vs_control",
    showNotificationFn = function(message, type, duration) {
      recorded$message <- message
      recorded$type <- type
      recorded$duration <- duration
    }
  )

  expect_equal(
    recorded$message,
    "Enrichment analysis completed successfully for contrast: treated_vs_control"
  )
  expect_equal(recorded$type, "message")
  expect_equal(recorded$duration, 5)
  expect_equal(
    completed$message,
    "Enrichment analysis completed successfully for contrast: treated_vs_control"
  )
  expect_equal(completed$type, "message")
  expect_equal(completed$duration, 5)
  expect_equal(completed$selectedContrast, "treated_vs_control")
})

test_that("notifyProtEnrichAnalysisError emits the error notification payload", {
  recorded <- new.env(parent = emptyenv())
  recorded$message <- NULL
  recorded$type <- NULL
  recorded$duration <- NULL

  completed <- notifyProtEnrichAnalysisError(
    errorMessage = "mock failure",
    showNotificationFn = function(message, type, duration) {
      recorded$message <- message
      recorded$type <- type
      recorded$duration <- duration
    }
  )

  expect_equal(recorded$message, "Error in enrichment analysis: mock failure")
  expect_equal(recorded$type, "error")
  expect_equal(recorded$duration, 10)
  expect_equal(completed$message, "Error in enrichment analysis: mock failure")
  expect_equal(completed$type, "error")
  expect_equal(completed$duration, 10)
  expect_equal(completed$errorMessage, "mock failure")
})

test_that("logProtEnrichAnalysisError emits the formatted error log marker", {
  recorded <- new.env(parent = emptyenv())
  recorded$messages <- character()

  completed <- logProtEnrichAnalysisError(
    errorMessage = "mock failure",
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$messages, "*** ERROR in enrichment analysis: mock failure ***\n")
  expect_equal(completed$message, "*** ERROR in enrichment analysis: mock failure ***\n")
  expect_equal(completed$errorMessage, "mock failure")
})

test_that("messageProtEnrichAnalysisError emits the formatted error message payload", {
  recorded <- new.env(parent = emptyenv())
  recorded$messages <- character()

  completed <- messageProtEnrichAnalysisError(
    errorMessage = "mock failure",
    messageFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$messages, "*** ERROR in enrichment analysis: mock failure")
  expect_equal(completed$message, "*** ERROR in enrichment analysis: mock failure")
  expect_equal(completed$errorMessage, "mock failure")
})

test_that("reportProtEnrichAnalysisError delegates the grouped error-reporting trio", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- reportProtEnrichAnalysisError(
    errorMessage = "mock failure",
    messageErrorFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("message:", errorMessage))
      list(kind = "message", errorMessage = errorMessage)
    },
    logErrorFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("log:", errorMessage))
      list(kind = "log", errorMessage = errorMessage)
    },
    notifyErrorFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("notify:", errorMessage))
      list(kind = "notify", errorMessage = errorMessage)
    }
  )

  expect_equal(
    recorded$calls,
    c(
      "message:mock failure",
      "log:mock failure",
      "notify:mock failure"
    )
  )
  expect_equal(completed$errorMessage, "mock failure")
  expect_equal(completed$messageResult$kind, "message")
  expect_equal(completed$logResult$kind, "log")
  expect_equal(completed$notificationResult$kind, "notify")
})

test_that("finalizeProtEnrichObserverFailure delegates error reporting and cleanup", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- finalizeProtEnrichObserverFailure(
    errorMessage = "mock failure",
    reportAnalysisErrorFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("report:", errorMessage))
      list(kind = "report", errorMessage = errorMessage)
    },
    removeWorkingNotificationFn = function() {
      recorded$calls <- c(recorded$calls, "cleanup")
      list(kind = "cleanup")
    }
  )

  expect_equal(recorded$calls, c("report:mock failure", "cleanup"))
  expect_equal(completed$errorMessage, "mock failure")
  expect_equal(completed$reportResult$kind, "report")
  expect_equal(completed$reportResult$errorMessage, "mock failure")
  expect_equal(completed$cleanupResult$kind, "cleanup")
})

test_that("reportProtEnrichCompletion delegates the grouped success-reporting pair", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- reportProtEnrichCompletion(
    selectedContrast = "treated_vs_control",
    notifyCompletionFn = function(selectedContrast) {
      recorded$calls <- c(recorded$calls, paste0("notify:", selectedContrast))
      list(kind = "notify", selectedContrast = selectedContrast)
    },
    logCompletionFn = function() {
      recorded$calls <- c(recorded$calls, "log")
      list(kind = "log")
    }
  )

  expect_equal(
    recorded$calls,
    c(
      "notify:treated_vs_control",
      "log"
    )
  )
  expect_equal(completed$selectedContrast, "treated_vs_control")
  expect_equal(completed$notificationResult$kind, "notify")
  expect_equal(completed$notificationResult$selectedContrast, "treated_vs_control")
  expect_equal(completed$logResult$kind, "log")
})

test_that("finalizeProtEnrichObserverRun delegates completion reporting and cleanup on success", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- finalizeProtEnrichObserverRun(
    completed = TRUE,
    selectedContrast = "treated_vs_control",
    reportCompletionFn = function(selectedContrast) {
      recorded$calls <- c(recorded$calls, paste0("report:", selectedContrast))
      list(kind = "report", selectedContrast = selectedContrast)
    },
    removeWorkingNotificationFn = function() {
      recorded$calls <- c(recorded$calls, "cleanup")
      list(kind = "cleanup")
    }
  )

  expect_equal(recorded$calls, c("report:treated_vs_control", "cleanup"))
  expect_true(completed$completed)
  expect_equal(completed$selectedContrast, "treated_vs_control")
  expect_equal(completed$reportResult$kind, "report")
  expect_equal(completed$reportResult$selectedContrast, "treated_vs_control")
  expect_equal(completed$cleanupResult$kind, "cleanup")
})

test_that("finalizeProtEnrichObserverRun still clears the working notification on failure", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- finalizeProtEnrichObserverRun(
    completed = FALSE,
    selectedContrast = "treated_vs_control",
    reportCompletionFn = function(selectedContrast) {
      recorded$calls <- c(recorded$calls, paste0("report:", selectedContrast))
      list(kind = "report", selectedContrast = selectedContrast)
    },
    removeWorkingNotificationFn = function() {
      recorded$calls <- c(recorded$calls, "cleanup")
      list(kind = "cleanup")
    }
  )

  expect_equal(recorded$calls, "cleanup")
  expect_false(completed$completed)
  expect_equal(completed$selectedContrast, "treated_vs_control")
  expect_null(completed$reportResult)
  expect_equal(completed$cleanupResult$kind, "cleanup")
})

test_that("runProtEnrichObserverShell delegates progress wrapping and success finalization", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- runProtEnrichObserverShell(
    selectedContrast = "treated_vs_control",
    runAnalysisBody = function() {
      recorded$calls <- c(recorded$calls, "run")
    },
    withProgressFn = function(message, value, expr) {
      recorded$calls <- c(recorded$calls, paste0("progress:", message, ":", value))
      force(expr)
    },
    finalizeFailureFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("failure:", errorMessage))
      list(kind = "failure", errorMessage = errorMessage)
    },
    finalizeObserverRunFn = function(completed, selectedContrast) {
      recorded$calls <- c(recorded$calls, paste0("success:", completed, ":", selectedContrast))
      list(kind = "success", completed = completed, selectedContrast = selectedContrast)
    }
  )

  expect_equal(
    recorded$calls,
    c(
      "progress:Running enrichment analysis...:0",
      "run",
      "success:TRUE:treated_vs_control"
    )
  )
  expect_true(completed$completed)
  expect_equal(completed$selectedContrast, "treated_vs_control")
  expect_null(completed$failureResult)
  expect_equal(completed$successResult$kind, "success")
  expect_true(completed$successResult$completed)
  expect_equal(completed$successResult$selectedContrast, "treated_vs_control")
})

test_that("runProtEnrichObserverShell delegates failure finalization when the body errors", {
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  completed <- runProtEnrichObserverShell(
    selectedContrast = "treated_vs_control",
    runAnalysisBody = function() {
      recorded$calls <- c(recorded$calls, "run")
      stop("mock failure")
    },
    withProgressFn = function(message, value, expr) {
      recorded$calls <- c(recorded$calls, paste0("progress:", message, ":", value))
      force(expr)
    },
    finalizeFailureFn = function(errorMessage) {
      recorded$calls <- c(recorded$calls, paste0("failure:", errorMessage))
      list(kind = "failure", errorMessage = errorMessage)
    },
    finalizeObserverRunFn = function(completed, selectedContrast) {
      recorded$calls <- c(recorded$calls, paste0("success:", completed, ":", selectedContrast))
      list(kind = "success", completed = completed, selectedContrast = selectedContrast)
    }
  )

  expect_equal(
    recorded$calls,
    c(
      "progress:Running enrichment analysis...:0",
      "run",
      "failure:mock failure"
    )
  )
  expect_false(completed$completed)
  expect_equal(completed$selectedContrast, "treated_vs_control")
  expect_equal(completed$failureResult$kind, "failure")
  expect_equal(completed$failureResult$errorMessage, "mock failure")
  expect_null(completed$successResult)
})

test_that("handoffProtEnrichObserverRun delegates notification, shell, and body handoff", {
  input <- list(selected_contrast = "A vs B")
  enrichment_data <- new.env(parent = emptyenv())
  workflow_data <- shiny::reactiveValues(tab_status = list(enrichment_analysis = "pending"))
  experiment_paths <- list(pathway_dir = "/tmp/pathway", results_dir = "/tmp/results")
  recorded <- new.env(parent = emptyenv())
  recorded$notification <- NULL
  recorded$shell <- NULL
  recorded$body <- NULL

  handoff <- handoffProtEnrichObserverRun(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    currentAnalysisMethodFn = function() list(method = "gprofiler2"),
    showNotificationFn = function(message, id, duration) {
      recorded$notification <- list(
        message = message,
        id = id,
        duration = duration
      )
    },
    runObserverShellFn = function(selectedContrast, runAnalysisBody) {
      recorded$shell <- list(
        selectedContrast = selectedContrast,
        runAnalysisBody = runAnalysisBody
      )
      bodyResult <- runAnalysisBody()
      list(
        kind = "shell",
        selectedContrast = selectedContrast,
        bodyResult = bodyResult
      )
    },
    runAnalysisBodyFn = function(input,
                                 enrichmentData,
                                 workflowData,
                                 experimentPaths,
                                 currentAnalysisMethodFn) {
      recorded$body <- list(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(kind = "body")
    }
  )

  expect_identical(
    recorded$notification,
    list(
      message = "Running enrichment analysis...",
      id = "enrichment_working",
      duration = NULL
    )
  )
  expect_equal(recorded$shell$selectedContrast, "A vs B")
  expect_true(is.function(recorded$shell$runAnalysisBody))
  expect_identical(recorded$body$input, input)
  expect_identical(recorded$body$enrichmentData, enrichment_data)
  expect_identical(recorded$body$workflowData, workflow_data)
  expect_identical(recorded$body$experimentPaths, experiment_paths)
  expect_true(is.function(recorded$body$currentAnalysisMethodFn))
  expect_equal(handoff$selectedContrast, "A vs B")
  expect_equal(handoff$notificationMessage, "Running enrichment analysis...")
  expect_equal(handoff$notificationId, "enrichment_working")
  expect_null(handoff$notificationDuration)
  expect_equal(handoff$shellResult$kind, "shell")
  expect_equal(handoff$shellResult$selectedContrast, "A vs B")
  expect_equal(handoff$shellResult$bodyResult$kind, "body")
})

test_that("runProtEnrichObserverPreflight delegates req validation before the handoff helper", {
  input <- list(selected_contrast = "A vs B")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$da_results_data <- list(contrast_a = list(result = "mapped"))
  workflow_data <- shiny::reactiveValues(tab_status = list(enrichment_analysis = "pending"))
  experiment_paths <- list(pathway_dir = "/tmp/pathway", results_dir = "/tmp/results")
  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()
  recorded$req <- NULL
  recorded$handoff <- NULL

  preflight <- runProtEnrichObserverPreflight(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    currentAnalysisMethodFn = function() list(method = "gprofiler2"),
    reqFn = function(selectedContrast, daResultsData) {
      recorded$calls <- c(recorded$calls, "req")
      recorded$req <- list(
        selectedContrast = selectedContrast,
        daResultsData = daResultsData
      )
      invisible(TRUE)
    },
    handoffObserverRunFn = function(input,
                                    enrichmentData,
                                    workflowData,
                                    experimentPaths,
                                    currentAnalysisMethodFn) {
      recorded$calls <- c(recorded$calls, "handoff")
      recorded$handoff <- list(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(kind = "handoff", selectedContrast = input$selected_contrast)
    }
  )

  expect_equal(recorded$calls, c("req", "handoff"))
  expect_equal(recorded$req$selectedContrast, "A vs B")
  expect_identical(recorded$req$daResultsData, enrichment_data$da_results_data)
  expect_identical(recorded$handoff$input, input)
  expect_identical(recorded$handoff$enrichmentData, enrichment_data)
  expect_identical(recorded$handoff$workflowData, workflow_data)
  expect_identical(recorded$handoff$experimentPaths, experiment_paths)
  expect_true(is.function(recorded$handoff$currentAnalysisMethodFn))
  expect_equal(preflight$selectedContrast, "A vs B")
  expect_identical(preflight$daResultsData, enrichment_data$da_results_data)
  expect_equal(preflight$handoffResult$kind, "handoff")
  expect_equal(preflight$handoffResult$selectedContrast, "A vs B")
})

test_that("registerProtEnrichRunObserver wires the start log and preflight helper through observeEvent", {
  input <- list(run_enrichment_analysis = 1, selected_contrast = "A vs B")
  enrichment_data <- new.env(parent = emptyenv())
  workflow_data <- shiny::reactiveValues(tab_status = list(enrichment_analysis = "pending"))
  experiment_paths <- list(pathway_dir = "/tmp/pathway", results_dir = "/tmp/results")
  current_analysis_method_fn <- function() list(method = "gprofiler2")
  recorded <- new.env(parent = emptyenv())
  recorded$eventExpr <- NULL
  recorded$messages <- character()
  recorded$preflight <- NULL

  observer <- registerProtEnrichRunObserver(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    currentAnalysisMethodFn = current_analysis_method_fn,
    observeEventFn = function(eventExpr, handlerExpr) {
      recorded$eventExpr <- paste(deparse(substitute(eventExpr)), collapse = "")
      eval(substitute(handlerExpr), parent.frame())
      "observer_registration"
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    },
    runObserverPreflightFn = function(input,
                                      enrichmentData,
                                      workflowData,
                                      experimentPaths,
                                      currentAnalysisMethodFn) {
      recorded$preflight <- list(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(kind = "preflight", selectedContrast = input$selected_contrast)
    }
  )

  expect_equal(recorded$eventExpr, "input$run_enrichment_analysis")
  expect_equal(recorded$messages, "=== STARTING ENRICHMENT ANALYSIS ===\n")
  expect_identical(recorded$preflight$input, input)
  expect_identical(recorded$preflight$enrichmentData, enrichment_data)
  expect_identical(recorded$preflight$workflowData, workflow_data)
  expect_identical(recorded$preflight$experimentPaths, experiment_paths)
  expect_identical(recorded$preflight$currentAnalysisMethodFn, current_analysis_method_fn)
  expect_equal(observer, "observer_registration")
})

test_that("registerProtEnrichSelectedTabObserver initializes enrichment state through the helper seam", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(current_state = "normalized")
  workflow_data$da_analysis_results_list <- list(
    contrast_a = data.frame(logFC = 1, stringsAsFactors = FALSE)
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- NULL
  enrichment_data$da_results_data <- NULL
  enrichment_data$contrasts_available <- NULL

  recorded <- new.env(parent = emptyenv())
  recorded$eventExpr <- NULL
  recorded$ignoreInit <- NULL
  recorded$messages <- character()
  recorded$update <- NULL
  recorded$notifications <- list()
  contrast_choices <- c("A vs B" = "contrast_a")

  observer <- registerProtEnrichSelectedTabObserver(
    selectedTabFn = function() "enrichment",
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = "mock_session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE) {
      recorded$eventExpr <- paste(deparse(substitute(eventExpr)), collapse = "")
      recorded$ignoreInit <- ignoreInit
      eval(substitute(handlerExpr), parent.frame())
    },
    resolveCurrentS4ObjectFn = function(workflowData, daResultsList) {
      expect_identical(workflowData$da_analysis_results_list, daResultsList)
      list(
        currentS4 = structure(list(label = "s4"), class = "mockS4"),
        source = "state_manager"
      )
    },
    buildContrastChoicesFn = function(daResultsList, contrastsTbl) {
      expect_identical(daResultsList, workflow_data$da_analysis_results_list)
      expect_equal(contrastsTbl$friendly_names, "A vs B")
      list(
        contrastsAvailable = "A vs B",
        contrastChoices = contrast_choices,
        source = "friendly_names"
      )
    },
    existsFn = function(name, envir) identical(name, "contrasts_tbl") && identical(envir, recorded$globalEnv),
    getFn = function(name, envir) {
      expect_identical(name, "contrasts_tbl")
      expect_identical(envir, recorded$globalEnv)
      envir$contrasts_tbl
    },
    globalEnv = {
      recorded$globalEnv <- new.env(parent = emptyenv())
      recorded$globalEnv$contrasts_tbl <- data.frame(
        contrasts = "contrast_a",
        friendly_names = "A vs B",
        stringsAsFactors = FALSE
      )
      recorded$globalEnv
    },
    updateSelectInputFn = function(session, inputId, choices) {
      recorded$update <- list(session = session, inputId = inputId, choices = choices)
    },
    showNotificationFn = function(...) {
      recorded$notifications <- c(recorded$notifications, list(list(...)))
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$eventExpr, "selectedTabFn()")
  expect_true(isTRUE(recorded$ignoreInit))
  expect_true(observer$initialized)
  expect_equal(observer$reason, "initialized")
  expect_equal(observer$currentState, "normalized")
  expect_equal(observer$source, "state_manager")
  expect_identical(enrichment_data$da_results_data, workflow_data$da_analysis_results_list)
  expect_s3_class(enrichment_data$current_s4_object, "mockS4")
  expect_equal(enrichment_data$contrasts_available, "A vs B")
  expect_identical(recorded$update$session, "mock_session")
  expect_equal(recorded$update$inputId, "selected_contrast")
  expect_identical(recorded$update$choices, contrast_choices)
  expect_length(recorded$notifications, 0)
  expect_true(any(grepl("ENRICHMENT TAB CLICKED", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("Using friendly names", recorded$messages, fixed = TRUE)))
})

test_that("setupProtEnrichSelectedTabObserverRegistration delegates setup through the wrapper-local seam", {
  workflow_data <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  selected_tab_fn <- function() "enrichment"
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL
  recorded$messages <- character()

  registration <- setupProtEnrichSelectedTabObserverRegistration(
    selectedTabFn = selected_tab_fn,
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = "mock_session",
    registerSelectedTabObserverFn = function(selectedTabFn,
                                             workflowData,
                                             enrichmentData,
                                             session,
                                             ...) {
      recorded$registration <- list(
        selectedTabFn = selectedTabFn,
        workflowData = workflowData,
        enrichmentData = enrichmentData,
        session = session
      )
      list(kind = "selected_tab_observer_registration")
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_true(isTRUE(registration$selectedTabProvided))
  expect_equal(registration$reason, "selected_tab_provided")
  expect_identical(
    registration$registration,
    list(kind = "selected_tab_observer_registration")
  )
  expect_identical(recorded$registration$selectedTabFn, selected_tab_fn)
  expect_identical(recorded$registration$workflowData, workflow_data)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$session, "mock_session")
  expect_true(any(grepl("Setting up tab selection observer", recorded$messages, fixed = TRUE)))
})

test_that("setupProtEnrichSelectedTabObserverRegistration skips setup when selected_tab is absent", {
  recorded <- new.env(parent = emptyenv())
  recorded$registrationCalls <- 0L
  recorded$messages <- character()

  registration <- setupProtEnrichSelectedTabObserverRegistration(
    selectedTabFn = NULL,
    workflowData = new.env(parent = emptyenv()),
    enrichmentData = new.env(parent = emptyenv()),
    session = "mock_session",
    registerSelectedTabObserverFn = function(...) {
      recorded$registrationCalls <- recorded$registrationCalls + 1L
      NULL
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_false(isTRUE(registration$selectedTabProvided))
  expect_equal(registration$reason, "selected_tab_missing")
  expect_null(registration$registration)
  expect_equal(recorded$registrationCalls, 0L)
  expect_true(any(grepl("No selected_tab parameter provided", recorded$messages, fixed = TRUE)))
})

test_that("setupProtEnrichDaResultsObserverRegistration delegates setup through the wrapper-local seam", {
  workflow_data <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichDaResultsObserverRegistration(
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = "mock_session",
    registerDaResultsObserverFn = function(workflowData,
                                           enrichmentData,
                                           session,
                                           ...) {
      recorded$registration <- list(
        workflowData = workflowData,
        enrichmentData = enrichmentData,
        session = session
      )
      list(kind = "da_results_observer_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "da_results_observer_registration")
  )
  expect_identical(recorded$registration$workflowData, workflow_data)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$session, "mock_session")
})

test_that("setupProtEnrichSelectedContrastObserverRegistration delegates setup through the wrapper-local seam", {
  input <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichSelectedContrastObserverRegistration(
    input = input,
    enrichmentData = enrichment_data,
    registerSelectedContrastObserverFn = function(input,
                                                  enrichmentData,
                                                  ...) {
      recorded$registration <- list(
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "selected_contrast_observer_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "selected_contrast_observer_registration")
  )
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichContrastsDisplayOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichContrastsDisplayOutputRegistration(
    output = output,
    enrichmentData = enrichment_data,
    registerContrastsDisplayOutputFn = function(output,
                                                enrichmentData,
                                                ...) {
      recorded$registration <- list(
        output = output,
        enrichmentData = enrichmentData
      )
      list(kind = "contrasts_display_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "contrasts_display_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichStatusOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL
  current_analysis_method <- function() {
    list(method = "gprofiler2", species_name = "Homo sapiens")
  }

  registration <- setupProtEnrichStatusOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    currentAnalysisMethodFn = current_analysis_method,
    registerStatusOutputFn = function(output,
                                      input,
                                      enrichmentData,
                                      currentAnalysisMethodFn,
                                      ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(kind = "status_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "status_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$currentAnalysisMethodFn, current_analysis_method)
})

test_that("setupProtEnrichDisplayStatusOutputBootstrap delegates display/status setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$analysisMethodDisplay <- NULL
  recorded$contrastsDisplay <- NULL
  recorded$statusOutput <- NULL
  current_analysis_method <- function() {
    list(method = "gprofiler2", species_name = "Homo sapiens")
  }

  bootstrap <- setupProtEnrichDisplayStatusOutputBootstrap(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    currentAnalysisMethodFn = current_analysis_method,
    setupAnalysisMethodDisplayOutputRegistrationFn = function(output,
                                                              currentAnalysisMethodFn,
                                                              ...) {
      recorded$analysisMethodDisplay <- list(
        output = output,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(
        registration = list(kind = "analysis_method_display_registration"),
        reason = "registered"
      )
    },
    setupContrastsDisplayOutputRegistrationFn = function(output,
                                                         enrichmentData,
                                                         ...) {
      recorded$contrastsDisplay <- list(
        output = output,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "contrasts_display_registration"),
        reason = "registered"
      )
    },
    setupStatusOutputRegistrationFn = function(output,
                                               input,
                                               enrichmentData,
                                               currentAnalysisMethodFn,
                                               ...) {
      recorded$statusOutput <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(
        registration = list(kind = "status_output_registration"),
        reason = "registered"
      )
    }
  )

  expect_equal(bootstrap$reason, "registered")
  expect_identical(
    bootstrap$analysisMethodDisplayOutputRegistration,
    list(kind = "analysis_method_display_registration")
  )
  expect_identical(
    bootstrap$contrastsDisplayOutputRegistration,
    list(kind = "contrasts_display_registration")
  )
  expect_identical(
    bootstrap$statusOutputRegistration,
    list(kind = "status_output_registration")
  )
  expect_identical(recorded$analysisMethodDisplay$output, output)
  expect_identical(recorded$analysisMethodDisplay$currentAnalysisMethodFn, current_analysis_method)
  expect_identical(recorded$contrastsDisplay$output, output)
  expect_identical(recorded$contrastsDisplay$enrichmentData, enrichment_data)
  expect_identical(recorded$statusOutput$output, output)
  expect_identical(recorded$statusOutput$input, input)
  expect_identical(recorded$statusOutput$enrichmentData, enrichment_data)
  expect_identical(recorded$statusOutput$currentAnalysisMethodFn, current_analysis_method)
})

test_that("setupProtEnrichRunObserverRegistration delegates setup through the wrapper-local seam", {
  input <- shiny::reactiveValues(selected_contrast = "A vs B")
  enrichment_data <- new.env(parent = emptyenv())
  workflow_data <- shiny::reactiveValues(tab_status = list(enrichment_analysis = "pending"))
  experiment_paths <- list(pathway_dir = "/tmp/pathway", results_dir = "/tmp/results")
  current_analysis_method <- function() {
    list(method = "gprofiler2", species_name = "Homo sapiens")
  }
  mock_preflight <- function(input,
                             enrichmentData,
                             workflowData,
                             experimentPaths,
                             currentAnalysisMethodFn,
                             ...) {
    list(kind = "preflight", selectedContrast = input$selected_contrast)
  }
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichRunObserverRegistration(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    currentAnalysisMethodFn = current_analysis_method,
    runObserverPreflightFn = mock_preflight,
    registerRunObserverFn = function(input,
                                     enrichmentData,
                                     workflowData,
                                     experimentPaths,
                                     currentAnalysisMethodFn,
                                     runObserverPreflightFn,
                                     ...) {
      recorded$registration <- list(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn,
        runObserverPreflightFn = runObserverPreflightFn
      )
      list(kind = "run_observer_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "run_observer_registration")
  )
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$workflowData, workflow_data)
  expect_identical(recorded$registration$experimentPaths, experiment_paths)
  expect_identical(recorded$registration$currentAnalysisMethodFn, current_analysis_method)
  expect_identical(recorded$registration$runObserverPreflightFn, mock_preflight)
})

test_that("setupProtEnrichGprofilerResultsTableOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(gprofiler_direction_filter = "up")
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichGprofilerResultsTableOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    registerGprofilerResultsTableOutputFn = function(output,
                                                     input,
                                                     enrichmentData,
                                                     ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "gprofiler_results_table_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "gprofiler_results_table_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichGprofilerSummaryOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(gprofiler_direction_filter = "up")
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichGprofilerSummaryOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    registerGprofilerSummaryOutputFn = function(output,
                                                input,
                                                enrichmentData,
                                                ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "gprofiler_summary_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "gprofiler_summary_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichClusterProfilerResultsTableOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(clusterprofiler_direction_filter = "down")
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichClusterProfilerResultsTableOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    registerClusterProfilerResultsTableOutputFn = function(output,
                                                           input,
                                                           enrichmentData,
                                                           ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "clusterprofiler_results_table_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "clusterprofiler_results_table_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichClusterProfilerSummaryOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(clusterprofiler_direction_filter = "down")
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichClusterProfilerSummaryOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    registerClusterProfilerSummaryOutputFn = function(output,
                                                      input,
                                                      enrichmentData,
                                                      ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "clusterprofiler_summary_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "clusterprofiler_summary_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichStringDbResultsTableOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    stringdb_filter_significant = TRUE,
    enrichment_p_val_thresh = 0.05,
    stringdb_max_results = 12
  )
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichStringDbResultsTableOutputRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    registerStringDbResultsTableOutputFn = function(output,
                                                    input,
                                                    enrichmentData,
                                                    ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(kind = "stringdb_results_table_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "stringdb_results_table_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichStringDbSummaryOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichStringDbSummaryOutputRegistration(
    output = output,
    enrichmentData = enrichment_data,
    registerStringDbSummaryOutputFn = function(output,
                                               enrichmentData,
                                               ...) {
      recorded$registration <- list(
        output = output,
        enrichmentData = enrichmentData
      )
      list(kind = "stringdb_summary_output_registration")
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(kind = "stringdb_summary_output_registration")
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichResultsSummaryOutputBootstrap delegates results and summary setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    gprofiler_direction_filter = "up",
    clusterprofiler_direction_filter = "down",
    stringdb_filter_significant = TRUE,
    enrichment_p_val_thresh = 0.05,
    stringdb_max_results = 12
  )
  enrichment_data <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$gprofilerResultsTable <- NULL
  recorded$gprofilerSummary <- NULL
  recorded$clusterProfilerResultsTable <- NULL
  recorded$clusterProfilerSummary <- NULL
  recorded$stringDbResultsTable <- NULL
  recorded$stringDbSummary <- NULL

  bootstrap <- setupProtEnrichResultsSummaryOutputBootstrap(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    setupGprofilerResultsTableOutputRegistrationFn = function(output,
                                                              input,
                                                              enrichmentData,
                                                              ...) {
      recorded$gprofilerResultsTable <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "gprofiler_results_table_registration"),
        reason = "registered"
      )
    },
    setupGprofilerSummaryOutputRegistrationFn = function(output,
                                                         input,
                                                         enrichmentData,
                                                         ...) {
      recorded$gprofilerSummary <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "gprofiler_summary_registration"),
        reason = "registered"
      )
    },
    setupClusterProfilerResultsTableOutputRegistrationFn = function(output,
                                                                    input,
                                                                    enrichmentData,
                                                                    ...) {
      recorded$clusterProfilerResultsTable <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "clusterprofiler_results_table_registration"),
        reason = "registered"
      )
    },
    setupClusterProfilerSummaryOutputRegistrationFn = function(output,
                                                               input,
                                                               enrichmentData,
                                                               ...) {
      recorded$clusterProfilerSummary <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "clusterprofiler_summary_registration"),
        reason = "registered"
      )
    },
    setupStringDbResultsTableOutputRegistrationFn = function(output,
                                                             input,
                                                             enrichmentData,
                                                             ...) {
      recorded$stringDbResultsTable <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "stringdb_results_table_registration"),
        reason = "registered"
      )
    },
    setupStringDbSummaryOutputRegistrationFn = function(output,
                                                        enrichmentData,
                                                        ...) {
      recorded$stringDbSummary <- list(
        output = output,
        enrichmentData = enrichmentData
      )
      list(
        registration = list(kind = "stringdb_summary_registration"),
        reason = "registered"
      )
    }
  )

  expect_equal(bootstrap$reason, "registered")
  expect_identical(
    bootstrap$gprofilerResultsTableOutputRegistration,
    list(kind = "gprofiler_results_table_registration")
  )
  expect_identical(
    bootstrap$gprofilerSummaryOutputRegistration,
    list(kind = "gprofiler_summary_registration")
  )
  expect_identical(
    bootstrap$clusterProfilerResultsTableOutputRegistration,
    list(kind = "clusterprofiler_results_table_registration")
  )
  expect_identical(
    bootstrap$clusterProfilerSummaryOutputRegistration,
    list(kind = "clusterprofiler_summary_registration")
  )
  expect_identical(
    bootstrap$stringDbResultsTableOutputRegistration,
    list(kind = "stringdb_results_table_registration")
  )
  expect_identical(
    bootstrap$stringDbSummaryOutputRegistration,
    list(kind = "stringdb_summary_registration")
  )
  expect_identical(recorded$gprofilerResultsTable$output, output)
  expect_identical(recorded$gprofilerResultsTable$input, input)
  expect_identical(recorded$gprofilerResultsTable$enrichmentData, enrichment_data)
  expect_identical(recorded$gprofilerSummary$output, output)
  expect_identical(recorded$gprofilerSummary$input, input)
  expect_identical(recorded$gprofilerSummary$enrichmentData, enrichment_data)
  expect_identical(recorded$clusterProfilerResultsTable$output, output)
  expect_identical(recorded$clusterProfilerResultsTable$input, input)
  expect_identical(recorded$clusterProfilerResultsTable$enrichmentData, enrichment_data)
  expect_identical(recorded$clusterProfilerSummary$output, output)
  expect_identical(recorded$clusterProfilerSummary$input, input)
  expect_identical(recorded$clusterProfilerSummary$enrichmentData, enrichment_data)
  expect_identical(recorded$stringDbResultsTable$output, output)
  expect_identical(recorded$stringDbResultsTable$input, input)
  expect_identical(recorded$stringDbResultsTable$enrichmentData, enrichment_data)
  expect_identical(recorded$stringDbSummary$output, output)
  expect_identical(recorded$stringDbSummary$enrichmentData, enrichment_data)
})

test_that("setupProtEnrichPlotOutputsRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    gprofiler_direction_filter = "all",
    clusterprofiler_direction_filter = "down"
  )
  enrichment_data <- new.env(parent = emptyenv())
  raw_contrast_name <- shiny::reactive("contrast_a")
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichPlotOutputsRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    rawContrastNameFn = raw_contrast_name,
    registerPlotOutputsFn = function(output,
                                     input,
                                     enrichmentData,
                                     rawContrastNameFn,
                                     ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        rawContrastNameFn = rawContrastNameFn
      )
      list(
        gprofilerPlot = "gprofiler-plot",
        clusterprofilerPlot = "clusterprofiler-plot"
      )
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    list(
      gprofilerPlot = "gprofiler-plot",
      clusterprofilerPlot = "clusterprofiler-plot"
    )
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$rawContrastNameFn, raw_contrast_name)
  expect_equal(shiny::isolate(recorded$registration$rawContrastNameFn()), "contrast_a")
})

test_that("setupProtEnrichStringDbPlotOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichStringDbPlotOutputRegistration(
    output = output,
    registerStringDbPlotOutputFn = function(output, ...) {
      recorded$registration <- list(output = output)
      "stringdb-plot-registration"
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(registration$registration, "stringdb-plot-registration")
  expect_identical(recorded$registration$output, output)
})

test_that("setupProtEnrichPlotOutputBootstrap delegates raw-contrast and plot-output setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    gprofiler_direction_filter = "all",
    clusterprofiler_direction_filter = "down"
  )
  enrichment_data <- new.env(parent = emptyenv())
  raw_contrast_name <- shiny::reactive("contrast_a")
  recorded <- new.env(parent = emptyenv())
  recorded$rawContrastInput <- NULL
  recorded$plotOutputsRegistration <- NULL
  recorded$stringDbPlotOutput <- NULL

  bootstrap <- setupProtEnrichPlotOutputBootstrap(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    setupRawContrastNameReactiveFn = function(input, ...) {
      recorded$rawContrastInput <- input
      list(
        rawContrastName = raw_contrast_name,
        reason = "created"
      )
    },
    setupPlotOutputsRegistrationFn = function(output,
                                              input,
                                              enrichmentData,
                                              rawContrastNameFn,
                                              ...) {
      recorded$plotOutputsRegistration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        rawContrastNameFn = rawContrastNameFn
      )
      list(
        registration = list(kind = "plot_outputs_registration"),
        reason = "registered"
      )
    },
    setupStringDbPlotOutputRegistrationFn = function(output, ...) {
      recorded$stringDbPlotOutput <- list(output = output)
      list(
        registration = "stringdb_plot_output_registration",
        reason = "registered"
      )
    }
  )

  expect_equal(bootstrap$reason, "registered")
  expect_identical(bootstrap$rawContrastName, raw_contrast_name)
  expect_identical(
    bootstrap$plotOutputsRegistration,
    list(kind = "plot_outputs_registration")
  )
  expect_identical(
    bootstrap$stringDbPlotOutputRegistration,
    "stringdb_plot_output_registration"
  )
  expect_identical(recorded$rawContrastInput, input)
  expect_identical(recorded$plotOutputsRegistration$output, output)
  expect_identical(recorded$plotOutputsRegistration$input, input)
  expect_identical(recorded$plotOutputsRegistration$enrichmentData, enrichment_data)
  expect_identical(recorded$plotOutputsRegistration$rawContrastNameFn, raw_contrast_name)
  expect_equal(shiny::isolate(recorded$plotOutputsRegistration$rawContrastNameFn()), "contrast_a")
  expect_identical(recorded$stringDbPlotOutput$output, output)
})

test_that("setupProtEnrichRunOutputDownloadBootstrap delegates run/results/plot/download setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    gprofiler_direction_filter = "all",
    clusterprofiler_direction_filter = "down",
    stringdb_filter_significant = TRUE,
    enrichment_p_val_thresh = 0.05,
    stringdb_max_results = 12
  )
  enrichment_data <- new.env(parent = emptyenv())
  workflow_data <- shiny::reactiveValues(tab_status = list(enrichment_analysis = "pending"))
  experiment_paths <- list(pathway_dir = "/tmp/pathway", results_dir = "/tmp/results")
  current_analysis_method <- shiny::reactive(list(method = "gprofiler2", supported = TRUE))
  raw_contrast_name <- shiny::reactive("contrast_a")
  mock_preflight <- function(input,
                             enrichmentData,
                             workflowData,
                             experimentPaths,
                             currentAnalysisMethodFn,
                             ...) {
    list(kind = "preflight", selectedContrast = input$selected_contrast)
  }
  recorded <- new.env(parent = emptyenv())
  recorded$runObserver <- NULL
  recorded$resultsSummary <- NULL
  recorded$plotOutput <- NULL
  recorded$resultsDownload <- NULL

  bootstrap <- setupProtEnrichRunOutputDownloadBootstrap(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    currentAnalysisMethodFn = current_analysis_method,
    runObserverPreflightFn = mock_preflight,
    setupRunObserverRegistrationFn = function(input,
                                              enrichmentData,
                                              workflowData,
                                              experimentPaths,
                                              currentAnalysisMethodFn,
                                              runObserverPreflightFn,
                                              ...) {
      recorded$runObserver <- list(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn,
        runObserverPreflightFn = runObserverPreflightFn
      )
      list(
        registration = list(kind = "run_observer_registration"),
        reason = "registered"
      )
    },
    setupResultsSummaryOutputBootstrapFn = function(output,
                                                    input,
                                                    enrichmentData,
                                                    ...) {
      recorded$resultsSummary <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        gprofilerResultsTableOutputRegistration = list(kind = "gprofiler_results_table_registration"),
        gprofilerSummaryOutputRegistration = list(kind = "gprofiler_summary_registration"),
        clusterProfilerResultsTableOutputRegistration = list(kind = "clusterprofiler_results_table_registration"),
        clusterProfilerSummaryOutputRegistration = list(kind = "clusterprofiler_summary_registration"),
        stringDbResultsTableOutputRegistration = list(kind = "stringdb_results_table_registration"),
        stringDbSummaryOutputRegistration = list(kind = "stringdb_summary_registration"),
        reason = "registered"
      )
    },
    setupPlotOutputBootstrapFn = function(output,
                                          input,
                                          enrichmentData,
                                          ...) {
      recorded$plotOutput <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        rawContrastName = raw_contrast_name,
        plotOutputsRegistration = list(kind = "plot_outputs_registration"),
        stringDbPlotOutputRegistration = "stringdb_plot_output_registration",
        reason = "registered"
      )
    },
    setupResultsDownloadHandlerRegistrationFn = function(output,
                                                         input,
                                                         enrichmentData,
                                                         currentAnalysisMethodFn,
                                                         ...) {
      recorded$resultsDownload <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(
        registration = "results_download_handler_registration",
        reason = "registered"
      )
    }
  )

  expect_equal(bootstrap$reason, "registered")
  expect_identical(bootstrap$runObserverRegistration, list(kind = "run_observer_registration"))
  expect_identical(
    bootstrap$gprofilerResultsTableOutputRegistration,
    list(kind = "gprofiler_results_table_registration")
  )
  expect_identical(
    bootstrap$gprofilerSummaryOutputRegistration,
    list(kind = "gprofiler_summary_registration")
  )
  expect_identical(
    bootstrap$clusterProfilerResultsTableOutputRegistration,
    list(kind = "clusterprofiler_results_table_registration")
  )
  expect_identical(
    bootstrap$clusterProfilerSummaryOutputRegistration,
    list(kind = "clusterprofiler_summary_registration")
  )
  expect_identical(
    bootstrap$stringDbResultsTableOutputRegistration,
    list(kind = "stringdb_results_table_registration")
  )
  expect_identical(
    bootstrap$stringDbSummaryOutputRegistration,
    list(kind = "stringdb_summary_registration")
  )
  expect_identical(bootstrap$rawContrastName, raw_contrast_name)
  expect_identical(bootstrap$plotOutputsRegistration, list(kind = "plot_outputs_registration"))
  expect_identical(bootstrap$stringDbPlotOutputRegistration, "stringdb_plot_output_registration")
  expect_identical(bootstrap$resultsDownloadHandlerRegistration, "results_download_handler_registration")
  expect_identical(recorded$runObserver$input, input)
  expect_identical(recorded$runObserver$enrichmentData, enrichment_data)
  expect_identical(recorded$runObserver$workflowData, workflow_data)
  expect_identical(recorded$runObserver$experimentPaths, experiment_paths)
  expect_identical(recorded$runObserver$currentAnalysisMethodFn, current_analysis_method)
  expect_identical(recorded$runObserver$runObserverPreflightFn, mock_preflight)
  expect_identical(recorded$resultsSummary$output, output)
  expect_identical(recorded$resultsSummary$input, input)
  expect_identical(recorded$resultsSummary$enrichmentData, enrichment_data)
  expect_identical(recorded$plotOutput$output, output)
  expect_identical(recorded$plotOutput$input, input)
  expect_identical(recorded$plotOutput$enrichmentData, enrichment_data)
  expect_identical(recorded$resultsDownload$output, output)
  expect_identical(recorded$resultsDownload$input, input)
  expect_identical(recorded$resultsDownload$enrichmentData, enrichment_data)
  expect_identical(recorded$resultsDownload$currentAnalysisMethodFn, current_analysis_method)
  expect_equal(shiny::isolate(bootstrap$rawContrastName()), "contrast_a")
})

test_that("setupProtEnrichAnalysisMethodDisplayOutputRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  current_analysis_method <- shiny::reactive(list(
    method = "gprofiler2",
    supported = TRUE
  ))
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichAnalysisMethodDisplayOutputRegistration(
    output = output,
    currentAnalysisMethodFn = current_analysis_method,
    registerAnalysisMethodDisplayOutputFn = function(output,
                                                     currentAnalysisMethodFn,
                                                     ...) {
      recorded$registration <- list(
        output = output,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      "analysis-method-display-registration"
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(
    registration$registration,
    "analysis-method-display-registration"
  )
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$currentAnalysisMethodFn, current_analysis_method)
})

test_that("setupProtEnrichTaxonIdObserverRegistration delegates setup through the wrapper-local seam", {
  workflow_data <- shiny::reactiveValues(taxon_id = "9606")
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichTaxonIdObserverRegistration(
    workflowData = workflow_data,
    session = mock_session,
    registerTaxonIdObserverFn = function(workflowData, session, ...) {
      recorded$registration <- list(
        workflowData = workflowData,
        session = session
      )
      "taxon-id-registration"
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(registration$registration, "taxon-id-registration")
  expect_identical(recorded$registration$workflowData, workflow_data)
  expect_identical(recorded$registration$session, mock_session)
})

test_that("setupProtEnrichMixedSpeciesObserverRegistration delegates setup through the wrapper-local seam", {
  workflow_data <- shiny::reactiveValues(mixed_species_analysis = list(enabled = TRUE))
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichMixedSpeciesObserverRegistration(
    workflowData = workflow_data,
    session = mock_session,
    registerMixedSpeciesObserverFn = function(workflowData, session, ...) {
      recorded$registration <- list(
        workflowData = workflowData,
        session = session
      )
      "mixed-species-registration"
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(registration$registration, "mixed-species-registration")
  expect_identical(recorded$registration$workflowData, workflow_data)
  expect_identical(recorded$registration$session, mock_session)
})

test_that("setupProtEnrichObserverRegistrationBootstrap delegates observer setup through the wrapper-local seam", {
  input <- shiny::reactiveValues(selected_contrast = "A vs B")
  workflow_data <- shiny::reactiveValues(
    taxon_id = "9606",
    mixed_species_analysis = list(enabled = TRUE)
  )
  enrichment_data <- new.env(parent = emptyenv())
  selected_tab_fn <- function() "enrichment"
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$taxon <- NULL
  recorded$mixedSpecies <- NULL
  recorded$selectedContrast <- NULL
  recorded$selectedTab <- NULL
  recorded$daResults <- NULL

  bootstrap <- setupProtEnrichObserverRegistrationBootstrap(
    selectedTabFn = selected_tab_fn,
    input = input,
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = mock_session,
    setupTaxonIdObserverRegistrationFn = function(workflowData, session, ...) {
      recorded$taxon <- list(
        workflowData = workflowData,
        session = session
      )
      list(
        registration = "taxon-id-registration",
        reason = "registered"
      )
    },
    setupMixedSpeciesObserverRegistrationFn = function(workflowData, session, ...) {
      recorded$mixedSpecies <- list(
        workflowData = workflowData,
        session = session
      )
      list(
        registration = "mixed-species-registration",
        reason = "registered"
      )
    },
    setupSelectedContrastObserverRegistrationFn = function(input, enrichmentData, ...) {
      recorded$selectedContrast <- list(
        input = input,
        enrichmentData = enrichmentData
      )
      list(
        registration = "selected-contrast-registration",
        reason = "registered"
      )
    },
    setupSelectedTabObserverRegistrationFn = function(selectedTabFn,
                                                      workflowData,
                                                      enrichmentData,
                                                      session,
                                                      ...) {
      recorded$selectedTab <- list(
        selectedTabFn = selectedTabFn,
        workflowData = workflowData,
        enrichmentData = enrichmentData,
        session = session
      )
      list(
        selectedTabProvided = TRUE,
        registration = "selected-tab-registration",
        reason = "selected_tab_provided"
      )
    },
    setupDaResultsObserverRegistrationFn = function(workflowData, enrichmentData, session, ...) {
      recorded$daResults <- list(
        workflowData = workflowData,
        enrichmentData = enrichmentData,
        session = session
      )
      list(
        registration = "da-results-registration",
        reason = "registered"
      )
    }
  )

  expect_equal(bootstrap$reason, "registered")
  expect_identical(bootstrap$taxonIdObserverRegistration, "taxon-id-registration")
  expect_identical(bootstrap$mixedSpeciesObserverRegistration, "mixed-species-registration")
  expect_identical(bootstrap$selectedContrastObserverRegistration, "selected-contrast-registration")
  expect_identical(bootstrap$selectedTabObserverRegistration, "selected-tab-registration")
  expect_true(isTRUE(bootstrap$selectedTabProvided))
  expect_identical(bootstrap$selectedTabReason, "selected_tab_provided")
  expect_identical(bootstrap$daResultsObserverRegistration, "da-results-registration")
  expect_identical(recorded$taxon$workflowData, workflow_data)
  expect_identical(recorded$taxon$session, mock_session)
  expect_identical(recorded$mixedSpecies$workflowData, workflow_data)
  expect_identical(recorded$mixedSpecies$session, mock_session)
  expect_identical(recorded$selectedContrast$input, input)
  expect_identical(recorded$selectedContrast$enrichmentData, enrichment_data)
  expect_identical(recorded$selectedTab$selectedTabFn, selected_tab_fn)
  expect_identical(recorded$selectedTab$workflowData, workflow_data)
  expect_identical(recorded$selectedTab$enrichmentData, enrichment_data)
  expect_identical(recorded$selectedTab$session, mock_session)
  expect_identical(recorded$daResults$workflowData, workflow_data)
  expect_identical(recorded$daResults$enrichmentData, enrichment_data)
  expect_identical(recorded$daResults$session, mock_session)
})

test_that("setupProtEnrichResultsDownloadHandlerRegistration delegates setup through the wrapper-local seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    organism_taxid = "9606"
  )
  enrichment_data <- new.env(parent = emptyenv())
  current_analysis_method <- shiny::reactive(list(
    method = "gprofiler2",
    supported = TRUE
  ))
  recorded <- new.env(parent = emptyenv())
  recorded$registration <- NULL

  registration <- setupProtEnrichResultsDownloadHandlerRegistration(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    currentAnalysisMethodFn = current_analysis_method,
    registerResultsDownloadHandlerFn = function(output,
                                                input,
                                                enrichmentData,
                                                currentAnalysisMethodFn,
                                                ...) {
      recorded$registration <- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      "results-download-registration"
    }
  )

  expect_equal(registration$reason, "registered")
  expect_identical(registration$registration, "results-download-registration")
  expect_identical(recorded$registration$output, output)
  expect_identical(recorded$registration$input, input)
  expect_identical(recorded$registration$enrichmentData, enrichment_data)
  expect_identical(recorded$registration$currentAnalysisMethodFn, current_analysis_method)
  expect_equal(shiny::isolate(recorded$registration$input$selected_contrast), "A vs B")
  expect_equal(shiny::isolate(recorded$registration$input$organism_taxid), "9606")
})

test_that("registerProtEnrichDaResultsObserver refreshes DE-backed contrast choices through the helper seam", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$da_analysis_results_list <- list(
    contrast_a = data.frame(logFC = 1, stringsAsFactors = FALSE)
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- NULL
  enrichment_data$da_results_data <- NULL
  enrichment_data$contrasts_available <- NULL

  recorded <- new.env(parent = emptyenv())
  recorded$eventExpr <- NULL
  recorded$ignoreInit <- NULL
  recorded$ignoreNULL <- NULL
  recorded$messages <- character()
  recorded$update <- NULL
  contrast_choices <- c("A vs B" = "contrast_a")

  observer <- registerProtEnrichDaResultsObserver(
    workflowData = workflow_data,
    enrichmentData = enrichment_data,
    session = "mock_session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE, ignoreNULL = FALSE) {
      recorded$eventExpr <- paste(deparse(substitute(eventExpr)), collapse = "")
      recorded$ignoreInit <- ignoreInit
      recorded$ignoreNULL <- ignoreNULL
      eval(substitute(handlerExpr), parent.frame())
    },
    resolveCurrentS4ObjectFn = function(workflowData, daResultsList) {
      expect_identical(workflowData$da_analysis_results_list, daResultsList)
      list(
        currentS4 = structure(list(label = "s4"), class = "mockS4"),
        source = "state_manager"
      )
    },
    buildContrastChoicesFn = function(daResultsList, contrastsTbl) {
      expect_identical(daResultsList, workflow_data$da_analysis_results_list)
      expect_equal(contrastsTbl$friendly_names, "A vs B")
      list(
        contrastsAvailable = "A vs B",
        contrastChoices = contrast_choices,
        source = "friendly_names"
      )
    },
    existsFn = function(name, envir) identical(name, "contrasts_tbl") && identical(envir, recorded$globalEnv),
    getFn = function(name, envir) {
      expect_identical(name, "contrasts_tbl")
      expect_identical(envir, recorded$globalEnv)
      envir$contrasts_tbl
    },
    globalEnv = {
      recorded$globalEnv <- new.env(parent = emptyenv())
      recorded$globalEnv$contrasts_tbl <- data.frame(
        contrasts = "contrast_a",
        friendly_names = "A vs B",
        stringsAsFactors = FALSE
      )
      recorded$globalEnv
    },
    updateSelectInputFn = function(session, inputId, choices) {
      recorded$update <- list(session = session, inputId = inputId, choices = choices)
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$eventExpr, "workflowData$da_analysis_results_list")
  expect_true(isTRUE(recorded$ignoreInit))
  expect_true(isTRUE(recorded$ignoreNULL))
  expect_true(observer$hasResults)
  expect_true(observer$currentS4Stored)
  expect_equal(observer$source, "state_manager")
  expect_equal(observer$reason, "updated")
  expect_identical(observer$contrastChoices, contrast_choices)
  expect_s3_class(enrichment_data$current_s4_object, "mockS4")
  expect_identical(enrichment_data$da_results_data, workflow_data$da_analysis_results_list)
  expect_equal(enrichment_data$contrasts_available, "A vs B")
  expect_identical(recorded$update$session, "mock_session")
  expect_equal(recorded$update$inputId, "selected_contrast")
  expect_identical(recorded$update$choices, contrast_choices)
  expect_true(any(grepl("DE results detected - updating contrasts", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("Using friendly names", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("Updated contrasts dropdown with 1 contrasts", recorded$messages, fixed = TRUE)))
})

test_that("registerProtEnrichSelectedContrastObserver hydrates mapped contrast results through the helper seam", {
  input <- list(selected_contrast = "A vs B")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$analysis_complete <- TRUE

  gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  clusterprofiler_results <- data.frame(term = "cp:1", stringsAsFactors = FALSE)
  stringdb_results <- data.frame(network = "ppi", stringsAsFactors = FALSE)
  enrichment_data$all_enrichment_results <- list(
    contrast_a = list(
      gprofiler_results = gprofiler_results,
      clusterprofiler_results = clusterprofiler_results,
      stringdb_results = stringdb_results
    )
  )
  enrichment_data$gprofiler_results <- NULL
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL

  global_env <- new.env(parent = emptyenv())
  global_env$contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )

  recorded <- new.env(parent = emptyenv())
  recorded$messages <- character()
  recorded$req <- list()

  observer <- registerProtEnrichSelectedContrastObserver(
    input = input,
    enrichmentData = enrichment_data,
    observeFn = function(handlerExpr) {
      eval(substitute(handlerExpr), parent.frame())
    },
    reqFn = function(value) {
      recorded$req <- c(recorded$req, list(value))
      value
    },
    globalEnv = global_env,
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_length(recorded$req, 3)
  expect_equal(observer$rawContrastName, "contrast_a")
  expect_true(observer$found)
  expect_identical(enrichment_data$gprofiler_results, gprofiler_results)
  expect_identical(enrichment_data$clusterprofiler_results, clusterprofiler_results)
  expect_identical(enrichment_data$stringdb_results, stringdb_results)
  expect_true(any(grepl("CONTRAST MAPPING", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("UPDATED: 1 gprofiler2 results", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("UPDATED: 1 clusterProfileR results", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("UPDATED: 1 stringDB results", recorded$messages, fixed = TRUE)))
})

test_that("registerProtEnrichSelectedContrastObserver clears stale result slots when a contrast is missing", {
  input <- list(selected_contrast = "missing_contrast")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$analysis_complete <- TRUE
  enrichment_data$all_enrichment_results <- list(
    contrast_a = list(
      gprofiler_results = data.frame(term = "go:1", stringsAsFactors = FALSE),
      clusterprofiler_results = data.frame(term = "cp:1", stringsAsFactors = FALSE),
      stringdb_results = data.frame(network = "ppi", stringsAsFactors = FALSE)
    )
  )
  enrichment_data$gprofiler_results <- data.frame(stale = "gp", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- data.frame(stale = "cp", stringsAsFactors = FALSE)
  enrichment_data$stringdb_results <- data.frame(stale = "sd", stringsAsFactors = FALSE)

  recorded <- new.env(parent = emptyenv())
  recorded$messages <- character()

  observer <- registerProtEnrichSelectedContrastObserver(
    input = input,
    enrichmentData = enrichment_data,
    observeFn = function(handlerExpr) {
      eval(substitute(handlerExpr), parent.frame())
    },
    reqFn = function(value) value,
    globalEnv = new.env(parent = emptyenv()),
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_false(observer$found)
  expect_equal(observer$availableContrasts, "contrast_a")
  expect_null(enrichment_data$gprofiler_results)
  expect_null(enrichment_data$clusterprofiler_results)
  expect_null(enrichment_data$stringdb_results)
  expect_true(any(grepl("WARNING: No results found", recorded$messages, fixed = TRUE)))
  expect_true(any(grepl("AVAILABLE CONTRASTS: contrast_a", recorded$messages, fixed = TRUE)))
})

test_that("runProtEnrichAnalysisBody delegates the extracted enrichment runner sequence", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$da_results_data <- list(contrast_a = list(result = "mapped"))
  enrichment_data$current_s4_object <- NULL
  enrichment_data$gprofiler_results <- NULL
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL
  enrichment_data$all_enrichment_results <- list()
  enrichment_data$annotation_match_results <- NULL
  enrichment_data$enrichment_results_full <- NULL
  enrichment_data$analysis_complete <- FALSE

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(
    setup_import = "complete",
    enrichment_analysis = "pending"
  )

  input <- list(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05,
    correction_method = "fdr",
    enable_organism_filter = FALSE
  )
  input_cfg <- input

  global_env <- new.env(parent = emptyenv())
  global_env$contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )

  recorded <- new.env(parent = emptyenv())
  recorded$progress <- character()
  recorded$setup <- NULL
  recorded$prepare <- NULL
  recorded$execute <- NULL
  recorded$finalize <- NULL

  enrichment_results <- methods::new(
    "mockProtEnrichAllContrastCarrier",
    enrichment_data = list(contrast_a = list())
  )
  all_contrast_results <- list(
    contrast_a = list(
      gprofiler_results = data.frame(term = "go:1", stringsAsFactors = FALSE),
      clusterprofiler_results = NULL,
      stringdb_results = NULL
    )
  )
  capture_post_process_stub <- function(...) {
    stop("analysis body should delegate post-process capture through finalizeAnalysisResultsFn")
  }
  persist_analysis_results_stub <- function(...) {
    stop("analysis body should delegate persistence through finalizeAnalysisResultsFn")
  }
  direct_selected_da_results_stub <- function(...) {
    stop("analysis body should delegate selected-contrast resolution through prepareAnalysisSetupFn")
  }
  direct_run_dependencies_stub <- function(...) {
    stop("analysis body should delegate dependency resolution through prepareAnalysisSetupFn")
  }
  direct_output_directories_stub <- function(...) {
    stop("analysis body should delegate output-directory resolution through prepareAnalysisSetupFn")
  }
  direct_create_da_results_stub <- function(...) {
    stop("analysis body should delegate DA-results S4 creation through prepareAnalysisSetupFn")
  }
  direct_uniprot_annotations_stub <- function(...) {
    stop("analysis body should delegate annotation resolution through prepareAnalysisSetupFn")
  }
  direct_annotation_matching_stub <- function(...) {
    stop("analysis body should delegate annotation matching through prepareAnalysisSetupFn")
  }
  direct_organism_mapping_stub <- function(...) {
    stop("analysis body should delegate organism mapping through prepareAnalysisSetupFn")
  }
  direct_organism_filter_stub <- function(...) {
    stop("analysis body should delegate organism filtering through prepareAnalysisSetupFn")
  }
  direct_persist_filter_metadata_stub <- function(...) {
    stop("analysis body should delegate filter metadata persistence through prepareAnalysisSetupFn")
  }
  show_notification_stub <- function(...) {
    stop("warning notification should not fire on the happy path")
  }

  completed <- runProtEnrichAnalysisBody(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    experimentPaths = list(
      da_output_dir = "/tmp/da",
      pathway_dir = "/tmp/pathway",
      results_dir = "/tmp/results"
    ),
    prepareAnalysisSetupFn = function(selectedContrast,
                                      input,
                                      enrichmentData,
                                      workflowData,
                                      experimentPaths,
                                      resolveSelectedDaResultsFn,
                                      resolveRunDependenciesFn,
                                      resolveOutputDirectoriesFn,
                                      createDaResultsForEnrichmentFn,
                                      resolveUniprotAnnotationsFn,
                                      resolveAnnotationMatchingFn,
                                      resolveOrganismMappingFn,
                                      applyOrganismFilterFn,
                                      persistOrganismFilterMetadataFn,
                                      showNotificationFn,
                                      globalEnv,
                                      catFn) {
      expect_equal(selectedContrast, "A vs B")
      expect_identical(input, input_cfg)
      expect_identical(enrichmentData, enrichment_data)
      expect_identical(workflowData, workflow_data)
      expect_equal(experimentPaths$results_dir, "/tmp/results")
      expect_identical(resolveSelectedDaResultsFn, direct_selected_da_results_stub)
      expect_identical(resolveRunDependenciesFn, direct_run_dependencies_stub)
      expect_identical(resolveOutputDirectoriesFn, direct_output_directories_stub)
      expect_identical(createDaResultsForEnrichmentFn, direct_create_da_results_stub)
      expect_identical(resolveUniprotAnnotationsFn, direct_uniprot_annotations_stub)
      expect_identical(resolveAnnotationMatchingFn, direct_annotation_matching_stub)
      expect_identical(resolveOrganismMappingFn, direct_organism_mapping_stub)
      expect_identical(applyOrganismFilterFn, direct_organism_filter_stub)
      expect_identical(persistOrganismFilterMetadataFn, direct_persist_filter_metadata_stub)
      expect_identical(showNotificationFn, show_notification_stub)
      expect_identical(globalEnv, global_env)
      expect_true(is.function(catFn))
      enrichmentData$current_s4_object <- "current_s4"
      enrichmentData$annotation_match_results <- list(match_rate = 100)
      workflowData$enrichment_organism_filter <- list(
        enabled = FALSE,
        applied = FALSE,
        target = "9606",
        stats = list(
          proteins_before = 0,
          proteins_after = 0,
          proteins_removed = 0
        )
      )
      recorded$setup <- list(
        rawContrastName = "contrast_a",
        pathwayDir = "/tmp/pathway",
        annotationEntry = "P1"
      )
      list(
        rawContrastName = "contrast_a",
        contrastsTbl = global_env$contrasts_tbl,
        pathwayDir = "/tmp/pathway",
        daResultsForEnrichment = structure(list(kind = "da_results_s4"), class = "mock_da_results_s4"),
        goAnnotations = data.frame(Entry = "P1", stringsAsFactors = FALSE),
        organismFilterApplied = FALSE,
        filterStats = list(
          proteins_before = 0,
          proteins_after = 0,
          proteins_removed = 0
        )
      )
    },
    resolveSelectedDaResultsFn = direct_selected_da_results_stub,
    resolveRunDependenciesFn = direct_run_dependencies_stub,
    resolveOutputDirectoriesFn = direct_output_directories_stub,
    createDaResultsForEnrichmentFn = direct_create_da_results_stub,
    resolveUniprotAnnotationsFn = direct_uniprot_annotations_stub,
    resolveAnnotationMatchingFn = direct_annotation_matching_stub,
    resolveOrganismMappingFn = direct_organism_mapping_stub,
    applyOrganismFilterFn = direct_organism_filter_stub,
    persistOrganismFilterMetadataFn = direct_persist_filter_metadata_stub,
    currentAnalysisMethodFn = function() {
      stop("analysis body should delegate method resolution through prepareProcessExecutionFn")
    },
    resolveAnalysisInputColumnsFn = function(...) {
      stop("analysis body should delegate input-column resolution through prepareProcessExecutionFn")
    },
    buildProcessEnrichmentsArgsFn = function(...) {
      stop("analysis body should delegate process-arg assembly through prepareProcessExecutionFn")
    },
    prepareProcessExecutionFn = function(input,
                                         enrichmentData,
                                         daResultsForEnrichment,
                                         pathwayDir,
                                         goAnnotations,
                                         currentAnalysisMethodFn,
                                         resolveAnalysisInputColumnsFn,
                                         buildProcessEnrichmentsArgsFn,
                                         catFn) {
      expect_equal(input$organism_taxid, "9606")
      expect_identical(enrichmentData, enrichment_data)
      expect_equal(daResultsForEnrichment$kind, "da_results_s4")
      expect_equal(pathwayDir, "/tmp/pathway")
      expect_equal(goAnnotations$Entry, "P1")
      expect_true(is.function(currentAnalysisMethodFn))
      expect_true(is.function(resolveAnalysisInputColumnsFn))
      expect_true(is.function(buildProcessEnrichmentsArgsFn))
      expect_true(is.function(catFn))
      recorded$prepare <- list(
        pathwayDir = pathwayDir,
        annotationEntry = goAnnotations$Entry,
        contrastNames = names(enrichmentData$da_results_data)
      )
      list(
        methodInfo = list(method = "gprofiler2", supported = TRUE, species_name = "Human"),
        inputColumnConfig = list(idColumn = "gene_name"),
        enrichmentArgs = list(
          checkpointArgs = list(kind = "checkpoint"),
          processArgs = list(kind = "process")
        )
      )
    },
    captureCheckpointFn = function(...) {
      stop("analysis body should delegate checkpoint capture through executeProcessEnrichmentsFn")
    },
    processEnrichmentsFn = function(...) {
      stop("analysis body should delegate processEnrichments through executeProcessEnrichmentsFn")
    },
    executeProcessEnrichmentsFn = function(enrichmentArgs,
                                           upCutoff,
                                           downCutoff,
                                           qCutoff,
                                           captureCheckpointFn,
                                           processEnrichmentsFn,
                                           catFn) {
      expect_equal(enrichmentArgs$checkpointArgs$kind, "checkpoint")
      expect_equal(enrichmentArgs$processArgs$kind, "process")
      expect_equal(upCutoff, 1)
      expect_equal(downCutoff, -1)
      expect_equal(qCutoff, 0.05)
      expect_true(is.function(captureCheckpointFn))
      expect_true(is.function(processEnrichmentsFn))
      expect_true(is.function(catFn))
      recorded$execute <- list(
        checkpointKind = enrichmentArgs$checkpointArgs$kind,
        processKind = enrichmentArgs$processArgs$kind,
        upCutoff = upCutoff,
        downCutoff = downCutoff,
        qCutoff = qCutoff
      )
      enrichment_results
    },
    buildAllContrastResultsFn = function(...) {
      stop("post-process helper should own all-contrast assembly")
    },
    resolveSelectedContrastResultsFn = function(...) {
      stop("post-process helper should own selected-contrast hydration")
    },
    capturePostProcessResultsFn = capture_post_process_stub,
    persistAnalysisResultsFn = persist_analysis_results_stub,
    finalizeAnalysisResultsFn = function(selectedContrast,
                                         rawContrastName,
                                         organismFilterApplied,
                                         filterStats,
                                         enrichmentResults,
                                         enrichmentData,
                                         workflowData,
                                         input,
                                         methodInfo,
                                         contrastsTbl,
                                         pathwayDir,
                                         buildAllContrastResultsFn,
                                         resolveSelectedContrastResultsFn,
                                         capturePostProcessResultsFn,
                                         persistAnalysisResultsFn,
                                         buildAnalysisResultsPayloadFn,
                                         propagateResultsArgsFn,
                                         propagateUiParamsFn,
                                         updateStateManagerUiParamsFn,
                                         saveCompletedStateFn,
                                         completeTabStatusFn,
                                         completeProgressFn,
                                         incProgressFn,
                                         catFn) {
      expect_equal(selectedContrast, "A vs B")
      expect_equal(rawContrastName, "contrast_a")
      expect_false(organismFilterApplied)
      expect_equal(filterStats$proteins_before, 0)
      expect_equal(filterStats$proteins_after, 0)
      expect_equal(filterStats$proteins_removed, 0)
      expect_identical(enrichmentResults, enrichment_results)
      expect_identical(enrichmentData, enrichment_data)
      expect_identical(workflowData, workflow_data)
      expect_identical(input, input_cfg)
      expect_identical(contrastsTbl, global_env$contrasts_tbl)
      expect_equal(methodInfo$method, "gprofiler2")
      expect_equal(pathwayDir, "/tmp/pathway")
      expect_true(is.function(buildAllContrastResultsFn))
      expect_true(is.function(resolveSelectedContrastResultsFn))
      expect_identical(capturePostProcessResultsFn, capture_post_process_stub)
      expect_identical(persistAnalysisResultsFn, persist_analysis_results_stub)
      expect_identical(buildAnalysisResultsPayloadFn, buildProtEnrichAnalysisResultsPayload)
      expect_identical(propagateResultsArgsFn, propagateProtEnrichResultsArgs)
      expect_identical(propagateUiParamsFn, propagateProtEnrichUiParams)
      expect_identical(updateStateManagerUiParamsFn, updateProtEnrichStateManagerUiParams)
      expect_identical(saveCompletedStateFn, saveProtEnrichCompletedState)
      expect_identical(completeTabStatusFn, completeProtEnrichTabStatus)
      expect_identical(completeProgressFn, completeProtEnrichProgress)
      expect_true(is.function(incProgressFn))
      expect_true(is.function(catFn))
      enrichmentData$all_enrichment_results <- all_contrast_results
      enrichmentData$gprofiler_results <- all_contrast_results$contrast_a$gprofiler_results
      enrichmentData$clusterprofiler_results <- NULL
      enrichmentData$stringdb_results <- NULL
      enrichmentData$enrichment_results_full <- enrichmentResults
      enrichmentData$analysis_complete <- TRUE
      enrichmentData$current_s4_object <- "updated_s4"
      workflowData$enrichment_analysis_results <- list(kind = "payload")
      workflowData$tab_status$enrichment_analysis <- "complete"
      incProgressFn(0.8, detail = "Storing results...")
      recorded$finalize <- list(
        selectedContrast = selectedContrast,
        rawContrastName = rawContrastName,
        organismFilterApplied = organismFilterApplied,
        filterStats = filterStats,
        methodInfo = methodInfo,
        pathwayDir = pathwayDir,
        gprofilerTerms = enrichmentData$gprofiler_results$term
      )
      list(
        selectedContrast = selectedContrast,
        rawContrastName = rawContrastName,
        analysisMethod = methodInfo$method,
        analysisComplete = TRUE,
        organismFilterApplied = organismFilterApplied,
        filterStats = filterStats,
        enrichmentResults = "results_with_args"
      )
    },
    showNotificationFn = show_notification_stub,
    incProgressFn = function(value, detail) {
      recorded$progress <- c(recorded$progress, paste(value, detail, sep = "|"))
    },
    globalEnv = global_env,
    catFn = function(...) invisible(NULL)
  )

  expect_equal(
    recorded$progress,
    c(
      "0.2|Transforming DE data...",
      "0.3|Creating DE results S4 object...",
      "0.5|Running enrichment analysis...",
      "0.8|Storing results..."
    )
  )
  expect_identical(
    recorded$setup,
    list(
      rawContrastName = "contrast_a",
      pathwayDir = "/tmp/pathway",
      annotationEntry = "P1"
    )
  )
  expect_identical(
    recorded$prepare,
    list(
      pathwayDir = "/tmp/pathway",
      annotationEntry = "P1",
      contrastNames = "contrast_a"
    )
  )
  expect_identical(
    recorded$execute,
    list(
      checkpointKind = "checkpoint",
      processKind = "process",
      upCutoff = 1,
      downCutoff = -1,
      qCutoff = 0.05
    )
  )
  expect_identical(enrichment_data$annotation_match_results, list(match_rate = 100))
  expect_identical(enrichment_data$all_enrichment_results, all_contrast_results)
  expect_equal(enrichment_data$gprofiler_results$term, "go:1")
  expect_null(enrichment_data$clusterprofiler_results)
  expect_null(enrichment_data$stringdb_results)
  expect_identical(enrichment_data$enrichment_results_full, enrichment_results)
  expect_true(enrichment_data$analysis_complete)
  expect_equal(enrichment_data$current_s4_object, "updated_s4")
  expect_identical(workflow_data$enrichment_analysis_results, list(kind = "payload"))
  expect_false(workflow_data$enrichment_organism_filter$enabled)
  expect_false(workflow_data$enrichment_organism_filter$applied)
  expect_equal(workflow_data$enrichment_organism_filter$target, "9606")
  expect_equal(workflow_data$tab_status$enrichment_analysis, "complete")
  expect_equal(recorded$finalize$selectedContrast, "A vs B")
  expect_equal(recorded$finalize$rawContrastName, "contrast_a")
  expect_false(recorded$finalize$organismFilterApplied)
  expect_equal(recorded$finalize$filterStats$proteins_before, 0)
  expect_equal(recorded$finalize$filterStats$proteins_after, 0)
  expect_equal(recorded$finalize$filterStats$proteins_removed, 0)
  expect_equal(recorded$finalize$methodInfo$method, "gprofiler2")
  expect_equal(recorded$finalize$pathwayDir, "/tmp/pathway")
  expect_equal(recorded$finalize$gprofilerTerms, "go:1")
  expect_equal(completed$selectedContrast, "A vs B")
  expect_equal(completed$rawContrastName, "contrast_a")
  expect_equal(completed$analysisMethod, "gprofiler2")
  expect_true(completed$analysisComplete)
  expect_false(completed$organismFilterApplied)
  expect_equal(completed$filterStats$proteins_before, 0)
  expect_equal(completed$filterStats$proteins_after, 0)
  expect_equal(completed$filterStats$proteins_removed, 0)
  expect_equal(completed$enrichmentResults, "results_with_args")
})

test_that("captureProtEnrichPostProcessResults stores all-contrast results and hydrates the selected display state", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$all_enrichment_results <- NULL
  enrichment_data$gprofiler_results <- NULL
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL

  enrichment_results <- methods::new(
    "mockProtEnrichAllContrastCarrier",
    enrichment_data = list(contrast_a = list())
  )
  all_contrast_results <- list(
    contrast_a = list(
      gprofiler_results = data.frame(term = "go:1", stringsAsFactors = FALSE),
      clusterprofiler_results = NULL,
      stringdb_results = data.frame(network = "ppi:1", stringsAsFactors = FALSE)
    )
  )

  recorded <- new.env(parent = emptyenv())
  recorded$calls <- character()

  result <- captureProtEnrichPostProcessResults(
    selectedContrast = "A vs B",
    enrichmentResults = enrichment_results,
    enrichmentData = enrichment_data,
    contrastsTbl = data.frame(
      contrasts = "contrast_a",
      friendly_names = "A vs B",
      stringsAsFactors = FALSE
    ),
    methodInfo = list(method = "gprofiler2"),
    buildAllContrastResultsFn = function(enrichmentResults, methodInfo) {
      recorded$calls <- c(recorded$calls, "build_all")
      expect_identical(enrichmentResults, enrichment_results)
      expect_equal(methodInfo$method, "gprofiler2")
      all_contrast_results
    },
    resolveSelectedContrastResultsFn = function(selectedContrast, allContrastResults, contrastsTbl) {
      recorded$calls <- c(recorded$calls, "resolve_selected")
      expect_equal(selectedContrast, "A vs B")
      expect_identical(allContrastResults, all_contrast_results)
      expect_equal(contrastsTbl$friendly_names, "A vs B")
      list(
        found = TRUE,
        rawContrastName = "contrast_a",
        gprofilerResults = allContrastResults$contrast_a$gprofiler_results,
        clusterprofilerResults = NULL,
        stringdbResults = allContrastResults$contrast_a$stringdb_results
      )
    },
    catFn = function(message) {
      recorded$calls <- c(recorded$calls, trimws(message))
    }
  )

  expect_equal(recorded$calls, c(
    "build_all",
    "resolve_selected",
    "ENRICHMENT Step: Set initial display to contrast 'contrast_a'"
  ))
  expect_identical(enrichment_data$all_enrichment_results, all_contrast_results)
  expect_equal(enrichment_data$gprofiler_results$term, "go:1")
  expect_null(enrichment_data$clusterprofiler_results)
  expect_equal(enrichment_data$stringdb_results$network, "ppi:1")
  expect_true(result$hasResults)
  expect_identical(result$allContrastResults, all_contrast_results)
  expect_true(result$initialContrastState$found)
  expect_equal(result$initialContrastState$rawContrastName, "contrast_a")
})

test_that("captureProtEnrichPostProcessResults logs empty results without mutating display state", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$all_enrichment_results <- list(existing = TRUE)
  enrichment_data$gprofiler_results <- data.frame(term = "existing", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- data.frame(term = "existing", stringsAsFactors = FALSE)
  enrichment_data$stringdb_results <- data.frame(network = "existing", stringsAsFactors = FALSE)

  recorded <- character()

  result <- captureProtEnrichPostProcessResults(
    selectedContrast = "A vs B",
    enrichmentResults = NULL,
    enrichmentData = enrichment_data,
    contrastsTbl = NULL,
    methodInfo = list(method = "gprofiler2"),
    buildAllContrastResultsFn = function(...) {
      stop("all-contrast collation should be skipped for empty results")
    },
    resolveSelectedContrastResultsFn = function(...) {
      stop("selected-contrast hydration should be skipped for empty results")
    },
    catFn = function(message) {
      recorded <<- c(recorded, trimws(message))
    }
  )

  expect_equal(recorded, "ENRICHMENT Step: processEnrichments returned NULL or empty results")
  expect_identical(enrichment_data$all_enrichment_results, list(existing = TRUE))
  expect_equal(enrichment_data$gprofiler_results$term, "existing")
  expect_equal(enrichment_data$clusterprofiler_results$term, "existing")
  expect_equal(enrichment_data$stringdb_results$network, "existing")
  expect_false(result$hasResults)
  expect_null(result$allContrastResults)
  expect_null(result$initialContrastState)
})

test_that("finalizeProtEnrichAnalysisBodyResults captures, persists, and returns the runner summary", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- "current_s4"
  enrichment_data$gprofiler_results <- data.frame(term = "existing", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL
  enrichment_data$enrichment_results_full <- NULL
  enrichment_data$analysis_complete <- FALSE

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(enrichment_analysis = "pending")

  input_cfg <- list(
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )
  contrasts_tbl <- data.frame(
    contrasts = "contrast_a",
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )
  filter_stats <- list(proteins_before = 10, proteins_after = 8, proteins_removed = 2)
  recorded <- new.env(parent = emptyenv())
  recorded$steps <- character()

  build_all_contrast_results_stub <- function(...) {
    list(kind = "build_all")
  }
  resolve_selected_contrast_results_stub <- function(...) {
    list(kind = "resolve_selected")
  }
  build_analysis_results_payload_stub <- function(...) {
    list(kind = "build_payload")
  }
  propagate_results_args_stub <- function(...) {
    list(kind = "propagate_results")
  }
  propagate_ui_params_stub <- function(...) {
    list(kind = "propagate_ui")
  }
  update_state_manager_ui_params_stub <- function(...) {
    list(kind = "update_state")
  }
  save_completed_state_stub <- function(...) {
    list(kind = "save_state")
  }
  complete_tab_status_stub <- function(...) {
    list(kind = "complete_tab")
  }
  complete_progress_stub <- function(...) {
    list(kind = "complete_progress")
  }

  result <- finalizeProtEnrichAnalysisBodyResults(
    selectedContrast = "A vs B",
    rawContrastName = "contrast_a",
    organismFilterApplied = TRUE,
    filterStats = filter_stats,
    enrichmentResults = "raw_results",
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    input = input_cfg,
    methodInfo = list(method = "gprofiler2"),
    contrastsTbl = contrasts_tbl,
    pathwayDir = "/tmp/pathway",
    buildAllContrastResultsFn = build_all_contrast_results_stub,
    resolveSelectedContrastResultsFn = resolve_selected_contrast_results_stub,
    capturePostProcessResultsFn = function(selectedContrast,
                                           enrichmentResults,
                                           enrichmentData,
                                           contrastsTbl,
                                           methodInfo,
                                           buildAllContrastResultsFn,
                                           resolveSelectedContrastResultsFn,
                                           catFn) {
      recorded$steps <- c(recorded$steps, "capture")
      expect_equal(selectedContrast, "A vs B")
      expect_equal(enrichmentResults, "raw_results")
      expect_identical(enrichmentData, enrichment_data)
      expect_identical(contrastsTbl, contrasts_tbl)
      expect_equal(methodInfo$method, "gprofiler2")
      expect_identical(buildAllContrastResultsFn, build_all_contrast_results_stub)
      expect_identical(resolveSelectedContrastResultsFn, resolve_selected_contrast_results_stub)
      expect_true(is.function(catFn))
      enrichmentData$gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
      list(kind = "capture")
    },
    persistAnalysisResultsFn = function(input,
                                        enrichmentData,
                                        workflowData,
                                        selectedContrast,
                                        enrichmentResults,
                                        methodInfo,
                                        pathwayDir,
                                        buildAnalysisResultsPayloadFn,
                                        propagateResultsArgsFn,
                                        propagateUiParamsFn,
                                        updateStateManagerUiParamsFn,
                                        saveCompletedStateFn,
                                        completeTabStatusFn,
                                        completeProgressFn,
                                        incProgressFn) {
      recorded$steps <- c(recorded$steps, "persist")
      expect_identical(input, input_cfg)
      expect_identical(enrichmentData, enrichment_data)
      expect_identical(workflowData, workflow_data)
      expect_equal(selectedContrast, "A vs B")
      expect_equal(enrichmentResults, "raw_results")
      expect_equal(methodInfo$method, "gprofiler2")
      expect_equal(pathwayDir, "/tmp/pathway")
      expect_identical(buildAnalysisResultsPayloadFn, build_analysis_results_payload_stub)
      expect_identical(propagateResultsArgsFn, propagate_results_args_stub)
      expect_identical(propagateUiParamsFn, propagate_ui_params_stub)
      expect_identical(updateStateManagerUiParamsFn, update_state_manager_ui_params_stub)
      expect_identical(saveCompletedStateFn, save_completed_state_stub)
      expect_identical(completeTabStatusFn, complete_tab_status_stub)
      expect_identical(completeProgressFn, complete_progress_stub)
      expect_true(is.function(incProgressFn))
      list(
        enrichmentResults = "results_with_args",
        analysisComplete = TRUE,
        currentS4Object = "updated_s4"
      )
    },
    buildAnalysisResultsPayloadFn = build_analysis_results_payload_stub,
    propagateResultsArgsFn = propagate_results_args_stub,
    propagateUiParamsFn = propagate_ui_params_stub,
    updateStateManagerUiParamsFn = update_state_manager_ui_params_stub,
    saveCompletedStateFn = save_completed_state_stub,
    completeTabStatusFn = complete_tab_status_stub,
    completeProgressFn = complete_progress_stub,
    incProgressFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )

  expect_equal(recorded$steps, c("capture", "persist"))
  expect_equal(result$selectedContrast, "A vs B")
  expect_equal(result$rawContrastName, "contrast_a")
  expect_equal(result$analysisMethod, "gprofiler2")
  expect_true(result$analysisComplete)
  expect_true(result$organismFilterApplied)
  expect_identical(result$filterStats, filter_stats)
  expect_equal(result$enrichmentResults, "results_with_args")
  expect_equal(enrichment_data$gprofiler_results$term, "go:1")
})

test_that("persistProtEnrichAnalysisResults stores payload, propagates args, and finalizes completion", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$current_s4_object <- "current_s4"
  enrichment_data$gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- NULL
  enrichment_data$stringdb_results <- NULL
  enrichment_data$enrichment_results_full <- NULL
  enrichment_data$analysis_complete <- FALSE

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$tab_status <- list(enrichment_analysis = "pending")

  input <- list(
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )

  recorded <- new.env(parent = emptyenv())
  recorded$steps <- character()
  recorded$save_state <- NULL
  recorded$update_state <- NULL

  result <- persistProtEnrichAnalysisResults(
    input = input,
    enrichmentData = enrichment_data,
    workflowData = workflow_data,
    selectedContrast = "A vs B",
    enrichmentResults = "raw_results",
    methodInfo = list(method = "gprofiler2"),
    pathwayDir = "/tmp/pathway",
    buildAnalysisResultsPayloadFn = function(gprofilerResults,
                                             clusterprofilerResults,
                                             stringdbResults,
                                             fullEnrichmentResults,
                                             selectedContrast,
                                             analysisMethod,
                                             upCutoff,
                                             downCutoff,
                                             qCutoff,
                                             organismTaxid) {
      recorded$steps <- c(recorded$steps, "build_payload")
      expect_equal(gprofilerResults$term, "go:1")
      expect_null(clusterprofilerResults)
      expect_null(stringdbResults)
      expect_equal(fullEnrichmentResults, "raw_results")
      expect_equal(selectedContrast, "A vs B")
      expect_equal(analysisMethod, "gprofiler2")
      expect_equal(upCutoff, 1)
      expect_equal(downCutoff, -1)
      expect_equal(qCutoff, 0.05)
      expect_equal(organismTaxid, "9606")
      list(kind = "payload")
    },
    propagateResultsArgsFn = function(enrichmentResults,
                                      currentS4Object,
                                      selectedContrast,
                                      methodInfo,
                                      upCutoff,
                                      downCutoff,
                                      qCutoff,
                                      organismTaxid,
                                      pathwayDir) {
      recorded$steps <- c(recorded$steps, "propagate_results")
      expect_equal(enrichmentResults, "raw_results")
      expect_equal(currentS4Object, "current_s4")
      expect_equal(selectedContrast, "A vs B")
      expect_equal(methodInfo$method, "gprofiler2")
      expect_equal(upCutoff, 1)
      expect_equal(downCutoff, -1)
      expect_equal(qCutoff, 0.05)
      expect_equal(organismTaxid, "9606")
      expect_equal(pathwayDir, "/tmp/pathway")
      list(enrichmentResults = "results_with_args")
    },
    propagateUiParamsFn = function(currentS4Object,
                                   workflowData,
                                   selectedContrast,
                                   methodInfo,
                                   upCutoff,
                                   downCutoff,
                                   qCutoff,
                                   organismTaxid) {
      recorded$steps <- c(recorded$steps, "propagate_ui")
      expect_equal(currentS4Object, "current_s4")
      expect_identical(workflowData, workflow_data)
      expect_equal(selectedContrast, "A vs B")
      expect_equal(methodInfo$method, "gprofiler2")
      expect_equal(upCutoff, 1)
      expect_equal(downCutoff, -1)
      expect_equal(qCutoff, 0.05)
      expect_equal(organismTaxid, "9606")
      list(currentS4Object = "updated_s4", storedUiParams = TRUE)
    },
    updateStateManagerUiParamsFn = function(workflowData, storedUiParams) {
      recorded$steps <- c(recorded$steps, "update_state")
      recorded$update_state <- list(workflowData = workflowData, storedUiParams = storedUiParams)
    },
    saveCompletedStateFn = function(workflowData,
                                    enrichmentResults,
                                    selectedContrast,
                                    methodInfo,
                                    upCutoff,
                                    downCutoff,
                                    qCutoff,
                                    organismTaxid,
                                    pathwayDir) {
      recorded$steps <- c(recorded$steps, "save_state")
      recorded$save_state <- list(
        workflowData = workflowData,
        enrichmentResults = enrichmentResults,
        selectedContrast = selectedContrast,
        methodInfo = methodInfo,
        upCutoff = upCutoff,
        downCutoff = downCutoff,
        qCutoff = qCutoff,
        organismTaxid = organismTaxid,
        pathwayDir = pathwayDir
      )
    },
    completeTabStatusFn = function(workflowData) {
      recorded$steps <- c(recorded$steps, "complete_tab")
      workflowData$tab_status$enrichment_analysis <- "complete"
      list(status = "complete")
    },
    completeProgressFn = function() {
      recorded$steps <- c(recorded$steps, "complete_progress")
      list(done = TRUE)
    },
    incProgressFn = function(value, detail) {
      recorded$steps <- c(recorded$steps, paste("progress", value, detail, sep = ":"))
    }
  )

  expect_equal(
    recorded$steps,
    c(
      "progress:0.8:Storing results...",
      "build_payload",
      "propagate_results",
      "propagate_ui",
      "update_state",
      "save_state",
      "complete_tab",
      "complete_progress"
    )
  )
  expect_equal(enrichment_data$enrichment_results_full, "raw_results")
  expect_true(enrichment_data$analysis_complete)
  expect_equal(enrichment_data$current_s4_object, "updated_s4")
  expect_identical(workflow_data$enrichment_analysis_results, list(kind = "payload"))
  expect_equal(workflow_data$tab_status$enrichment_analysis, "complete")
  expect_identical(recorded$update_state$workflowData, workflow_data)
  expect_true(recorded$update_state$storedUiParams)
  expect_identical(recorded$save_state$workflowData, workflow_data)
  expect_equal(recorded$save_state$enrichmentResults, "results_with_args")
  expect_equal(recorded$save_state$selectedContrast, "A vs B")
  expect_equal(recorded$save_state$methodInfo$method, "gprofiler2")
  expect_equal(recorded$save_state$upCutoff, 1)
  expect_equal(recorded$save_state$downCutoff, -1)
  expect_equal(recorded$save_state$qCutoff, 0.05)
  expect_equal(recorded$save_state$organismTaxid, "9606")
  expect_equal(recorded$save_state$pathwayDir, "/tmp/pathway")
  expect_equal(result$enrichmentResults, "results_with_args")
  expect_true(result$analysisComplete)
  expect_equal(result$currentS4Object, "updated_s4")
})

test_that("runProtEnrichAnalysisBody stops when the selected contrast cannot be resolved", {
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$da_results_data <- list(other_contrast = list(result = "present"))
  enrichment_data$current_s4_object <- NULL

  input <- list(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05,
    correction_method = "fdr",
    enable_organism_filter = FALSE
  )

  expect_error(
    runProtEnrichAnalysisBody(
      input = input,
      enrichmentData = enrichment_data,
      workflowData = new.env(parent = emptyenv()),
      experimentPaths = list(results_dir = "/tmp/results"),
      currentAnalysisMethodFn = function() {
        stop("analysis method should not be reached")
      },
      resolveSelectedDaResultsFn = function(...) {
        list(
          selectedDaResults = NULL,
          rawContrastName = "contrast_a",
          source = NULL,
          mappedRawContrastName = "contrast_a",
          availableKeys = "other_contrast"
        )
      },
      incProgressFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      globalEnv = new.env(parent = emptyenv())
    ),
    "Could not find DE results for contrast 'A vs B'. Available contrasts: other_contrast"
  )
})


test_that("mod_prot_enrich_server characterizes run/output/download bootstrap delegation into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$bootstrap <- NULL
  mock_preflight <- function(input,
                             enrichmentData,
                             workflowData,
                             experimentPaths,
                             currentAnalysisMethodFn,
                             ...) {
    list(kind = "preflight", selectedContrast = input$selected_contrast)
  }

  local_mocked_bindings(
    setupProtEnrichRunOutputDownloadBootstrap = function(output,
                                                         input,
                                                         enrichmentData,
                                                         workflowData,
                                                         experimentPaths,
                                                         currentAnalysisMethodFn,
                                                         runObserverPreflightFn,
                                                         ...) {
      recorded$bootstrap <<- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn,
        runObserverPreflightFn = runObserverPreflightFn
      )
      list(
        runObserverRegistration = list(kind = "run_observer_registration"),
        gprofilerResultsTableOutputRegistration = list(kind = "gprofiler_results_table_setup"),
        gprofilerSummaryOutputRegistration = list(kind = "gprofiler_summary_setup"),
        clusterProfilerResultsTableOutputRegistration = list(kind = "clusterprofiler_results_table_setup"),
        clusterProfilerSummaryOutputRegistration = list(kind = "clusterprofiler_summary_setup"),
        stringDbResultsTableOutputRegistration = list(kind = "stringdb_results_table_setup"),
        stringDbSummaryOutputRegistration = list(kind = "stringdb_summary_setup"),
        rawContrastName = shiny::reactive("contrast_a"),
        plotOutputsRegistration = list(kind = "plot_outputs_setup"),
        stringDbPlotOutputRegistration = list(kind = "stringdb_plot_setup"),
        resultsDownloadHandlerRegistration = list(kind = "results_download_setup"),
        reason = "registered"
      )
    },
    runProtEnrichObserverPreflight = mock_preflight,
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = mock_session,
    {
      session$setInputs(
        selected_contrast = "A vs B",
        organism_taxid = "9606",
        gprofiler_direction_filter = "all",
        clusterprofiler_direction_filter = "all",
        stringdb_filter_significant = TRUE,
        enrichment_p_val_thresh = 0.05,
        stringdb_max_results = 12
      )
      session$flushReact()

      expect_false(is.null(recorded$bootstrap))
      expect_identical(recorded$bootstrap$output, output)
      expect_identical(recorded$bootstrap$workflowData, workflow_data)
      expect_identical(recorded$bootstrap$experimentPaths, experiment_paths)
      expect_identical(recorded$bootstrap$enrichmentData, enrichment_data)
      expect_true(is.function(recorded$bootstrap$currentAnalysisMethodFn))
      expect_identical(recorded$bootstrap$runObserverPreflightFn, mock_preflight)
      expect_equal(
        shiny::isolate(recorded$bootstrap$input$selected_contrast),
        "A vs B"
      )
      expect_equal(
        shiny::isolate(recorded$bootstrap$input$organism_taxid),
        "9606"
      )
      expect_equal(
        shiny::isolate(recorded$bootstrap$input$gprofiler_direction_filter),
        "all"
      )
      expect_equal(
        shiny::isolate(recorded$bootstrap$input$clusterprofiler_direction_filter),
        "all"
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$stringdb_filter_significant),
        TRUE
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$enrichment_p_val_thresh),
        0.05
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$stringdb_max_results),
        12
      )
    }
  )
})

test_that("mod_prot_enrich_server characterizes analysis-method bootstrap delegation into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$bootstrapArgs <- NULL
  recorded$currentAnalysisMethodFn <- NULL
  recorded$registration <- NULL

  local_mocked_bindings(
    setupProtEnrichAnalysisMethodBootstrap = function(input,
                                                      enrichmentData,
                                                      ...) {
      recorded$bootstrapArgs <<- list(
        input = input,
        enrichmentData = enrichmentData
      )
      recorded$currentAnalysisMethodFn <<- shiny::reactive(list(
        method = "gprofiler2",
        supported = TRUE,
        species_name = "Homo sapiens"
      ))
      list(
        supportedOrganisms = shiny::reactive(tibble::tribble(
          ~taxid, ~id, ~name,
          "9606", "hsapiens", "Homo sapiens"
        )),
        currentAnalysisMethod = recorded$currentAnalysisMethodFn,
        reason = "created"
      )
    },
    registerProtEnrichAnalysisMethodDisplayOutput = function(output,
                                                             currentAnalysisMethodFn,
                                                             ...) {
      recorded$registration <<- list(
        output = output,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      output$analysis_method_display <- shiny::renderText("mock analysis method")
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = mock_session,
    {
      session$setInputs(organism_taxid = "9606")
      session$flushReact()

      expect_false(is.null(recorded$bootstrapArgs))
      expect_identical(recorded$bootstrapArgs$input, input)
      expect_identical(recorded$bootstrapArgs$enrichmentData, enrichment_data)
      expect_false(is.null(recorded$registration))
      expect_identical(recorded$registration$output, output)
      expect_identical(
        recorded$registration$currentAnalysisMethodFn,
        recorded$currentAnalysisMethodFn
      )
    }
  )
})

test_that("mod_prot_enrich_server characterizes reactive-values setup delegation into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$setupCalls <- 0L
  recorded$messages <- character()
  recorded$enrichmentData <- NULL
  custom_enrichment_data <- createProtEnrichReactiveValues()

  local_mocked_bindings(
    setupProtEnrichReactiveValues = function(...) {
      recorded$setupCalls <<- recorded$setupCalls + 1L
      recorded$messages <<- c(
        recorded$messages,
        "   mod_prot_enrich_server Step: Reactive values initialized\n"
      )
      list(
        enrichmentData = custom_enrichment_data,
        reason = "created"
      )
    },
    setupProtEnrichAnalysisMethodBootstrap = function(input,
                                                      enrichmentData,
                                                      ...) {
      recorded$enrichmentData <<- enrichmentData
      list(
        supportedOrganisms = shiny::reactive(tibble::tribble(
          ~taxid, ~id, ~name,
          "9606", "hsapiens", "Homo sapiens"
        )),
        currentAnalysisMethod = shiny::reactive(list(
          method = "gprofiler2",
          supported = TRUE,
          species_name = "Homo sapiens"
        )),
        reason = "created"
      )
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = mock_session,
    {
      session$setInputs(organism_taxid = "9606")
      session$flushReact()

      expect_equal(recorded$setupCalls, 1L)
      expect_identical(recorded$enrichmentData, custom_enrichment_data)
      expect_identical(enrichment_data, custom_enrichment_data)
      expect_equal(
        recorded$messages,
        "   mod_prot_enrich_server Step: Reactive values initialized\n"
      )
    }
  )
})

test_that("mod_prot_enrich_server characterizes current-analysis-method bootstrap handoff into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$bootstrap <- NULL
  recorded$currentAnalysisMethodFn <- NULL
  recorded$registration <- NULL

  local_mocked_bindings(
    setupProtEnrichAnalysisMethodBootstrap = function(input,
                                                      enrichmentData,
                                                      ...) {
      recorded$bootstrap <<- list(
        input = input,
        enrichmentData = enrichmentData
      )
      recorded$currentAnalysisMethodFn <<- shiny::reactive(list(
        method = "gprofiler2",
        supported = TRUE,
        species_name = "Homo sapiens"
      ))
      list(
        supportedOrganisms = shiny::reactive(tibble::tribble(
          ~taxid, ~id, ~name,
          "9606", "hsapiens", "Homo sapiens"
        )),
        currentAnalysisMethod = recorded$currentAnalysisMethodFn,
        reason = "created"
      )
    },
    registerProtEnrichAnalysisMethodDisplayOutput = function(output,
                                                             currentAnalysisMethodFn,
                                                             ...) {
      recorded$registration <<- list(
        output = output,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      output$analysis_method_display <- shiny::renderText("mock analysis method")
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = mock_session,
    {
      session$setInputs(organism_taxid = "9606")
      session$flushReact()

      expect_false(is.null(recorded$bootstrap))
      expect_identical(recorded$bootstrap$input, input)
      expect_identical(recorded$bootstrap$enrichmentData, enrichment_data)

      expect_false(is.null(recorded$registration))
      expect_identical(recorded$registration$output, output)
      expect_identical(
        recorded$registration$currentAnalysisMethodFn,
        recorded$currentAnalysisMethodFn
      )
    }
  )
})

test_that("mod_prot_enrich_server characterizes display/status output bootstrap delegation into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$bootstrap <- NULL

  local_mocked_bindings(
    setupProtEnrichDisplayStatusOutputBootstrap = function(output,
                                                           input,
                                                           enrichmentData,
                                                           currentAnalysisMethodFn,
                                                           ...) {
      recorded$bootstrap <<- list(
        output = output,
        input = input,
        enrichmentData = enrichmentData,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
      list(
        analysisMethodDisplayOutputRegistration = list(kind = "analysis_method_display_setup"),
        contrastsDisplayOutputRegistration = list(kind = "contrasts_display_setup"),
        statusOutputRegistration = list(kind = "status_output_setup"),
        reason = "registered"
      )
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = NULL
    ),
    session = mock_session,
    {
      session$setInputs(
        selected_contrast = "A vs B",
        organism_taxid = "9606",
        up_cutoff = 1,
        down_cutoff = -1,
        q_cutoff = 0.05
      )
      session$flushReact()

      expect_false(is.null(recorded$bootstrap))
      expect_identical(recorded$bootstrap$output, output)
      expect_identical(recorded$bootstrap$enrichmentData, enrichment_data)
      expect_true(is.function(recorded$bootstrap$currentAnalysisMethodFn))
      expect_equal(
        shiny::isolate(recorded$bootstrap$currentAnalysisMethodFn())$method,
        "gprofiler2"
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$selected_contrast),
        "A vs B"
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$organism_taxid),
        "9606"
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$up_cutoff),
        1
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$down_cutoff),
        -1
      )
      expect_identical(
        shiny::isolate(recorded$bootstrap$input$q_cutoff),
        0.05
      )
    }
  )
})

test_that("registerProtEnrichTaxonIdObserver updates organism_taxid through the helper seam", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$taxon_id <- "9606"
  recorded <- new.env(parent = emptyenv())
  recorded$eventValue <- NULL
  recorded$ignoreInit <- NULL
  recorded$ignoreNULL <- NULL
  recorded$update <- NULL

  observer <- registerProtEnrichTaxonIdObserver(
    workflowData = workflow_data,
    session = "mock_session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE, ignoreNULL = FALSE) {
      recorded$eventValue <- eventExpr
      recorded$ignoreInit <- ignoreInit
      recorded$ignoreNULL <- ignoreNULL
      eval(
        substitute(handlerExpr),
        envir = list2env(list(workflowData = workflow_data), parent = parent.frame())
      )
    },
    updateTextInputFn = function(session, inputId, value) {
      recorded$update <- list(session = session, inputId = inputId, value = value)
    }
  )

  expect_identical(recorded$eventValue, "9606")
  expect_true(isTRUE(recorded$ignoreInit))
  expect_true(isTRUE(recorded$ignoreNULL))
  expect_true(isTRUE(observer$updated))
  expect_identical(observer$taxonId, "9606")
  expect_equal(observer$reason, "updated")
  expect_identical(
    recorded$update,
    list(session = "mock_session", inputId = "organism_taxid", value = "9606")
  )
})

test_that("registerProtEnrichMixedSpeciesObserver auto-enables organism filtering through the helper seam", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$mixed_species_analysis <- list(
    enabled = TRUE,
    selected_organism = "Homo sapiens"
  )
  recorded <- new.env(parent = emptyenv())
  recorded$eventValue <- NULL
  recorded$ignoreInit <- NULL
  recorded$ignoreNULL <- NULL
  recorded$update <- NULL
  recorded$notification <- NULL
  recorded$messages <- character()

  observer <- registerProtEnrichMixedSpeciesObserver(
    workflowData = workflow_data,
    session = "mock_session",
    observeEventFn = function(eventExpr, handlerExpr, ignoreInit = FALSE, ignoreNULL = FALSE) {
      recorded$eventValue <- eventExpr
      recorded$ignoreInit <- ignoreInit
      recorded$ignoreNULL <- ignoreNULL
      eval(
        substitute(handlerExpr),
        envir = list2env(list(workflowData = workflow_data), parent = parent.frame())
      )
    },
    updateCheckboxInputFn = function(session, inputId, value) {
      recorded$update <- list(session = session, inputId = inputId, value = value)
    },
    showNotificationFn = function(message, type, duration) {
      recorded$notification <- list(message = message, type = type, duration = duration)
    },
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_identical(recorded$eventValue, workflow_data$mixed_species_analysis)
  expect_true(isTRUE(recorded$ignoreInit))
  expect_true(isTRUE(recorded$ignoreNULL))
  expect_true(isTRUE(observer$filterEnabled))
  expect_true(isTRUE(observer$notificationShown))
  expect_equal(observer$reason, "enabled_filter")
  expect_identical(observer$mixedSpeciesAnalysis, workflow_data$mixed_species_analysis)
  expect_identical(
    recorded$update,
    list(session = "mock_session", inputId = "enable_organism_filter", value = TRUE)
  )
  expect_identical(
    recorded$notification,
    list(
      message = "Multi-species FASTA detected. Filtering to Homo sapiens enabled.",
      type = "message",
      duration = 5
    )
  )
  expect_true(any(grepl("Auto-enabled organism filter", recorded$messages, fixed = TRUE)))
})

test_that("mod_prot_enrich_server characterizes observer bootstrap delegation into the helper seam", {
  workflow_data <- shiny::reactiveValues(
    taxon_id = NULL,
    mixed_species_analysis = NULL,
    tab_status = list(enrichment_analysis = "incomplete")
  )
  experiment_paths <- list(
    da_output_dir = "/tmp/da",
    pathway_dir = "/tmp/pathway",
    results_dir = "/tmp/results"
  )
  selected_tab_fn <- function() "enrichment"
  mock_session <- shiny::MockShinySession$new()
  recorded <- new.env(parent = emptyenv())
  recorded$bootstrap <- NULL

  local_mocked_bindings(
    setupProtEnrichObserverRegistrationBootstrap = function(selectedTabFn,
                                                            input,
                                                            workflowData,
                                                            enrichmentData,
                                                            session,
                                                            setupTaxonIdObserverRegistrationFn,
                                                            setupMixedSpeciesObserverRegistrationFn,
                                                            setupSelectedContrastObserverRegistrationFn,
                                                            setupSelectedTabObserverRegistrationFn,
                                                            setupDaResultsObserverRegistrationFn,
                                                            ...) {
      recorded$bootstrap <<- list(
        selectedTabFn = selectedTabFn,
        input = input,
        workflowData = workflowData,
        enrichmentData = enrichmentData,
        session = session,
        setupTaxonIdObserverRegistrationFn = setupTaxonIdObserverRegistrationFn,
        setupMixedSpeciesObserverRegistrationFn = setupMixedSpeciesObserverRegistrationFn,
        setupSelectedContrastObserverRegistrationFn = setupSelectedContrastObserverRegistrationFn,
        setupSelectedTabObserverRegistrationFn = setupSelectedTabObserverRegistrationFn,
        setupDaResultsObserverRegistrationFn = setupDaResultsObserverRegistrationFn
      )
      list(
        taxonIdObserverRegistration = list(kind = "taxon_id_setup"),
        mixedSpeciesObserverRegistration = list(kind = "mixed_species_setup"),
        selectedContrastObserverRegistration = list(kind = "selected_contrast_setup"),
        selectedTabObserverRegistration = list(kind = "selected_tab_setup"),
        selectedTabProvided = TRUE,
        selectedTabReason = "selected_tab_provided",
        daResultsObserverRegistration = list(kind = "da_results_setup"),
        reason = "registered"
      )
    },
    .env = environment(mod_prot_enrich_server)
  )

  testServer(
    mod_prot_enrich_server,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = selected_tab_fn
    ),
    session = mock_session,
    {
      session$setInputs(selected_contrast = "A vs B", organism_taxid = "9606")
      session$flushReact()

      expect_false(is.null(recorded$bootstrap))
      expect_identical(recorded$bootstrap$selectedTabFn, selected_tab_fn)
      expect_identical(recorded$bootstrap$input, input)
      expect_identical(recorded$bootstrap$workflowData, workflow_data)
      expect_identical(recorded$bootstrap$enrichmentData, enrichment_data)
      expect_identical(recorded$bootstrap$session, session)
      expect_identical(
        recorded$bootstrap$setupTaxonIdObserverRegistrationFn,
        setupProtEnrichTaxonIdObserverRegistration
      )
      expect_identical(
        recorded$bootstrap$setupMixedSpeciesObserverRegistrationFn,
        setupProtEnrichMixedSpeciesObserverRegistration
      )
      expect_identical(
        recorded$bootstrap$setupSelectedContrastObserverRegistrationFn,
        setupProtEnrichSelectedContrastObserverRegistration
      )
      expect_identical(
        recorded$bootstrap$setupSelectedTabObserverRegistrationFn,
        setupProtEnrichSelectedTabObserverRegistration
      )
      expect_identical(
        recorded$bootstrap$setupDaResultsObserverRegistrationFn,
        setupProtEnrichDaResultsObserverRegistration
      )
    }
  )
})

test_that("registerProtEnrichContrastsDisplayOutput wires renderText through the formatter seam", {
  output <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$contrasts_available <- c("A vs B", "C vs D")
  recorded <- new.env(parent = emptyenv())
  recorded$contrastsAvailable <- NULL

  contrasts_display <- registerProtEnrichContrastsDisplayOutput(
    output = output,
    enrichmentData = enrichment_data,
    renderTextFn = function(expr) expr,
    formatContrastsTextFn = function(contrastsAvailable) {
      recorded$contrastsAvailable <- contrastsAvailable
      "formatted contrasts"
    }
  )

  expect_identical(recorded$contrastsAvailable, c("A vs B", "C vs D"))
  expect_identical(contrasts_display, "formatted contrasts")
  expect_identical(output$contrasts_display, "formatted contrasts")
})

test_that("registerProtEnrichStatusOutput wires renderText through the formatter seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(
    selected_contrast = "A vs B",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$analysis_complete <- TRUE
  enrichment_data$gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- data.frame(term = "cp:1", stringsAsFactors = FALSE)
  enrichment_data$stringdb_results <- data.frame(network = "ppi:1", stringsAsFactors = FALSE)
  recorded <- new.env(parent = emptyenv())
  recorded$methodCalls <- 0L
  recorded$args <- NULL

  enrichment_status <- registerProtEnrichStatusOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    currentAnalysisMethodFn = function() {
      recorded$methodCalls <- recorded$methodCalls + 1L
      list(method = "gprofiler2", species_name = "Homo sapiens")
    },
    renderTextFn = function(expr) shiny::isolate(expr),
    formatStatusTextFn = function(...) {
      recorded$args <- list(...)
      "formatted status"
    }
  )

  expect_identical(recorded$methodCalls, 1L)
  expect_true(isTRUE(recorded$args$analysisComplete))
  expect_identical(recorded$args$methodInfo$method, "gprofiler2")
  expect_identical(recorded$args$selectedContrast, "A vs B")
  expect_identical(recorded$args$upCutoff, 1)
  expect_identical(recorded$args$downCutoff, -1)
  expect_identical(recorded$args$qCutoff, 0.05)
  expect_identical(recorded$args$gprofilerResults, enrichment_data$gprofiler_results)
  expect_identical(recorded$args$clusterprofilerResults, enrichment_data$clusterprofiler_results)
  expect_identical(recorded$args$stringdbResults, enrichment_data$stringdb_results)
  expect_identical(enrichment_status, "formatted status")
  expect_identical(output$enrichment_status, "formatted status")
})

test_that("registerProtEnrichGprofilerSummaryOutput wires renderText through the formatter seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(gprofiler_direction_filter = "up")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  gprofiler_summary <- registerProtEnrichGprofilerSummaryOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    renderTextFn = function(expr) shiny::isolate(expr),
    formatSummaryTextFn = function(gprofilerResults, directionFilter) {
      recorded$args <- list(
        gprofilerResults = gprofilerResults,
        directionFilter = directionFilter
      )
      "formatted gprofiler summary"
    }
  )

  expect_identical(recorded$args$gprofilerResults, enrichment_data$gprofiler_results)
  expect_identical(recorded$args$directionFilter, "up")
  expect_identical(gprofiler_summary, "formatted gprofiler summary")
  expect_identical(output$gprofiler_summary_stats, "formatted gprofiler summary")
})

test_that("registerProtEnrichGprofilerResultsTableOutput wires rendering through the table seam", {
  output <- new.env(parent = emptyenv())
  input <- list(gprofiler_direction_filter = "up")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$gprofiler_results <- data.frame(
    term = "go:1",
    directionality = "positive",
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  gprofiler_table <- registerProtEnrichGprofilerResultsTableOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    renderResultsTableFn = function(gprofilerResults, directionFilter) {
      recorded$args <- list(
        gprofilerResults = gprofilerResults,
        directionFilter = directionFilter
      )
      "formatted gprofiler table"
    }
  )

  expect_identical(recorded$args$gprofilerResults, enrichment_data$gprofiler_results)
  expect_identical(recorded$args$directionFilter, "up")
  expect_identical(gprofiler_table, "formatted gprofiler table")
  expect_identical(output$gprofiler_results_table, "formatted gprofiler table")
})

test_that("registerProtEnrichClusterProfilerSummaryOutput wires renderText through the formatter seam", {
  output <- new.env(parent = emptyenv())
  input <- shiny::reactiveValues(clusterprofiler_direction_filter = "down")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$clusterprofiler_results <- data.frame(
    Description = "Cell cycle",
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  clusterprofiler_summary <- registerProtEnrichClusterProfilerSummaryOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    renderTextFn = function(expr) shiny::isolate(expr),
    formatSummaryTextFn = function(clusterprofilerResults, directionFilter) {
      recorded$args <- list(
        clusterprofilerResults = clusterprofilerResults,
        directionFilter = directionFilter
      )
      "formatted clusterprofiler summary"
    }
  )

  expect_identical(recorded$args$clusterprofilerResults, enrichment_data$clusterprofiler_results)
  expect_identical(recorded$args$directionFilter, "down")
  expect_identical(clusterprofiler_summary, "formatted clusterprofiler summary")
  expect_identical(output$clusterprofiler_summary_stats, "formatted clusterprofiler summary")
})

test_that("registerProtEnrichClusterProfilerResultsTableOutput wires rendering through the table seam", {
  output <- new.env(parent = emptyenv())
  input <- list(clusterprofiler_direction_filter = "down")
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$clusterprofiler_results <- data.frame(
    Description = "Cell cycle",
    directionality = "down",
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  clusterprofiler_table <- registerProtEnrichClusterProfilerResultsTableOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    renderResultsTableFn = function(clusterprofilerResults, directionFilter) {
      recorded$args <- list(
        clusterprofilerResults = clusterprofilerResults,
        directionFilter = directionFilter
      )
      "formatted clusterprofiler table"
    }
  )

  expect_identical(recorded$args$clusterprofilerResults, enrichment_data$clusterprofiler_results)
  expect_identical(recorded$args$directionFilter, "down")
  expect_identical(clusterprofiler_table, "formatted clusterprofiler table")
  expect_identical(output$clusterprofiler_results_table, "formatted clusterprofiler table")
})

test_that("registerProtEnrichStringDbResultsTableOutput wires rendering through the table seam", {
  output <- new.env(parent = emptyenv())
  input <- list(
    stringdb_filter_significant = TRUE,
    enrichment_p_val_thresh = 0.05,
    stringdb_max_results = 10
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$stringdb_results <- data.frame(
    network = "ppi:1",
    p_value = 0.01,
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  stringdb_table <- registerProtEnrichStringDbResultsTableOutput(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    renderResultsTableFn = function(stringdbResults,
                                    filterSignificant,
                                    enrichmentPValThresh,
                                    maxResults) {
      recorded$args <- list(
        stringdbResults = stringdbResults,
        filterSignificant = filterSignificant,
        enrichmentPValThresh = enrichmentPValThresh,
        maxResults = maxResults
      )
      "formatted stringdb table"
    }
  )

  expect_identical(recorded$args$stringdbResults, enrichment_data$stringdb_results)
  expect_identical(recorded$args$filterSignificant, TRUE)
  expect_identical(recorded$args$enrichmentPValThresh, 0.05)
  expect_identical(recorded$args$maxResults, 10)
  expect_identical(stringdb_table, "formatted stringdb table")
  expect_identical(output$stringdb_results_table, "formatted stringdb table")
})

test_that("registerProtEnrichStringDbSummaryOutput wires renderText through the formatter seam", {
  output <- new.env(parent = emptyenv())
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$stringdb_results <- data.frame(
    term = "network cluster",
    stringsAsFactors = FALSE
  )
  recorded <- new.env(parent = emptyenv())
  recorded$args <- NULL

  stringdb_summary <- registerProtEnrichStringDbSummaryOutput(
    output = output,
    enrichmentData = enrichment_data,
    renderTextFn = function(expr) shiny::isolate(expr),
    formatSummaryTextFn = function(stringdbResults) {
      recorded$args <- list(
        stringdbResults = stringdbResults
      )
      "formatted stringdb summary"
    }
  )

  expect_identical(recorded$args$stringdbResults, enrichment_data$stringdb_results)
  expect_identical(stringdb_summary, "formatted stringdb summary")
  expect_identical(output$stringdb_summary_stats, "formatted stringdb summary")
})

test_that("registerProtEnrichPlotOutputs preserves the paired plot-output shape through the render seams", {
  output <- new.env(parent = emptyenv())
  input <- list(
    gprofiler_direction_filter = "up",
    clusterprofiler_direction_filter = "down"
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$analysis_complete <- TRUE
  enrichment_data$enrichment_results_full <- list(results = "full-results")
  recorded <- new.env(parent = emptyenv())
  recorded$gprofiler <- NULL
  recorded$clusterprofiler <- NULL

  plot_outputs <- registerProtEnrichPlotOutputs(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    rawContrastNameFn = function() "contrast_raw",
    renderGprofilerPlotFn = function(analysisComplete, enrichmentResultsFull, rawContrast, directionFilter) {
      recorded$gprofiler <- list(
        analysisComplete = analysisComplete,
        enrichmentResultsFull = enrichmentResultsFull,
        rawContrast = rawContrast,
        directionFilter = directionFilter
      )
      list(kind = "gprofiler_plot")
    },
    renderClusterProfilerPlotFn = function(analysisComplete, enrichmentResultsFull, rawContrast, directionFilter) {
      recorded$clusterprofiler <- list(
        analysisComplete = analysisComplete,
        enrichmentResultsFull = enrichmentResultsFull,
        rawContrast = rawContrast,
        directionFilter = directionFilter
      )
      list(kind = "clusterprofiler_plot")
    }
  )

  expect_identical(
    recorded$gprofiler,
    list(
      analysisComplete = TRUE,
      enrichmentResultsFull = enrichment_data$enrichment_results_full,
      rawContrast = "contrast_raw",
      directionFilter = "up"
    )
  )
  expect_identical(
    recorded$clusterprofiler,
    list(
      analysisComplete = TRUE,
      enrichmentResultsFull = enrichment_data$enrichment_results_full,
      rawContrast = "contrast_raw",
      directionFilter = "down"
    )
  )
  expect_named(plot_outputs, c("gprofilerPlot", "clusterprofilerPlot"))
  expect_identical(plot_outputs$gprofilerPlot, list(kind = "gprofiler_plot"))
  expect_identical(plot_outputs$clusterprofilerPlot, list(kind = "clusterprofiler_plot"))
  expect_identical(output$gprofiler_plot, plot_outputs$gprofilerPlot)
  expect_identical(output$clusterprofiler_plot, plot_outputs$clusterprofilerPlot)
})

test_that("registerProtEnrichStringDbPlotOutput preserves the rendered STRING-DB plot payload", {
  output <- new.env(parent = emptyenv())

  stringdb_plot <- registerProtEnrichStringDbPlotOutput(
    output = output,
    renderStringDbPlotFn = function() list(kind = "stringdb_plot", status = "placeholder")
  )

  expect_identical(stringdb_plot, list(kind = "stringdb_plot", status = "placeholder"))
  expect_identical(output$stringdb_plot, stringdb_plot)
})

test_that("registerProtEnrichResultsDownloadHandler preserves the download request payload", {
  output <- new.env(parent = emptyenv())
  input <- list(
    selected_contrast = "A vs B",
    organism_taxid = "9606",
    up_cutoff = 1,
    down_cutoff = -1,
    q_cutoff = 0.05
  )
  enrichment_data <- new.env(parent = emptyenv())
  enrichment_data$gprofiler_results <- data.frame(term = "go:1", stringsAsFactors = FALSE)
  enrichment_data$clusterprofiler_results <- data.frame(term = "cp:1", stringsAsFactors = FALSE)
  enrichment_data$stringdb_results <- data.frame(network = "ppi:1", stringsAsFactors = FALSE)
  recorded <- new.env(parent = emptyenv())
  recorded$filenameContrast <- NULL
  recorded$archiveRequest <- NULL
  recorded$methodCalls <- 0L

  download_handler <- registerProtEnrichResultsDownloadHandler(
    output = output,
    input = input,
    enrichmentData = enrichment_data,
    currentAnalysisMethodFn = function() {
      recorded$methodCalls <- recorded$methodCalls + 1L
      list(method = "gprofiler2", species_name = "Homo sapiens")
    },
    downloadHandlerFn = function(filename, content) {
      list(filename = filename, content = content)
    },
    buildFilenameFn = function(selectedContrast) {
      recorded$filenameContrast <- selectedContrast
      "enrichment_results.zip"
    },
    writeArchiveFn = function(...) {
      recorded$archiveRequest <- list(...)
      "archive-written"
    }
  )

  expect_true(is.function(download_handler$filename))
  expect_true(is.function(download_handler$content))
  expect_identical(output$download_enrichment_results, download_handler)
  expect_identical(download_handler$filename(), "enrichment_results.zip")
  expect_identical(recorded$filenameContrast, "A vs B")

  download_result <- download_handler$content("/tmp/mock-enrichment.zip")

  expect_identical(download_result, "archive-written")
  expect_identical(recorded$methodCalls, 1L)
  expect_identical(recorded$archiveRequest$file, "/tmp/mock-enrichment.zip")
  expect_identical(recorded$archiveRequest$selectedContrast, "A vs B")
  expect_identical(recorded$archiveRequest$methodInfo, list(method = "gprofiler2", species_name = "Homo sapiens"))
  expect_identical(recorded$archiveRequest$organismTaxid, "9606")
  expect_identical(recorded$archiveRequest$upCutoff, 1)
  expect_identical(recorded$archiveRequest$downCutoff, -1)
  expect_identical(recorded$archiveRequest$qCutoff, 0.05)
  expect_identical(recorded$archiveRequest$gprofilerResults, enrichment_data$gprofiler_results)
  expect_identical(recorded$archiveRequest$clusterprofilerResults, enrichment_data$clusterprofiler_results)
  expect_identical(recorded$archiveRequest$stringdbResults, enrichment_data$stringdb_results)
})

test_that("logProtEnrichCompletion emits the completion log marker", {
  recorded <- new.env(parent = emptyenv())
  recorded$messages <- character()

  completed <- logProtEnrichCompletion(
    catFn = function(message) {
      recorded$messages <- c(recorded$messages, message)
    }
  )

  expect_equal(recorded$messages, "=== ENRICHMENT ANALYSIS COMPLETED ===\n")
  expect_equal(completed$message, "=== ENRICHMENT ANALYSIS COMPLETED ===\n")
})

test_that("removeProtEnrichWorkingNotification clears the working notification id", {
  recorded <- new.env(parent = emptyenv())
  recorded$notificationId <- NULL

  completed <- removeProtEnrichWorkingNotification(
    removeNotificationFn = function(notificationId) {
      recorded$notificationId <- notificationId
    }
  )

  expect_equal(recorded$notificationId, "enrichment_working")
  expect_equal(completed$notificationId, "enrichment_working")
})

test_that("buildProtEnrichSupportedOrganisms returns the static supported-organism lookup", {
  supported_organisms <- buildProtEnrichSupportedOrganisms()

  expect_s3_class(supported_organisms, "tbl_df")
  expect_equal(names(supported_organisms), c("taxid", "id", "name"))
  expect_equal(nrow(supported_organisms), 13)
  expect_identical(anyDuplicated(supported_organisms$taxid), 0L)
  expect_true(all(c("9606", "10090", "9598") %in% supported_organisms$taxid))

  human_row <- supported_organisms[supported_organisms$taxid == "9606", , drop = FALSE]
  expect_equal(human_row$id[[1]], "hsapiens")
  expect_equal(human_row$name[[1]], "Homo sapiens")
})

test_that("resolveProtEnrichAnalysisMethod returns supported organism metadata", {
  supported_organisms <- tibble::tribble(
    ~taxid, ~id, ~name,
    "9606", "hsapiens", "Homo sapiens",
    "10090", "mmusculus", "Mus musculus"
  )

  resolved <- resolveProtEnrichAnalysisMethod("9606", supported_organisms)

  expect_equal(resolved$analysisMethod, "gprofiler2")
  expect_true(resolved$organismSupported)
  expect_equal(resolved$methodInfo$method, "gprofiler2")
  expect_equal(resolved$methodInfo$species_id, "hsapiens")
  expect_equal(resolved$methodInfo$species_name, "Homo sapiens")
  expect_match(resolved$methodInfo$description, "gprofiler2 analysis for Homo sapiens", perl = TRUE)
})

test_that("resolveProtEnrichAnalysisMethod returns custom organism fallback metadata", {
  supported_organisms <- tibble::tribble(
    ~taxid, ~id, ~name,
    "9606", "hsapiens", "Homo sapiens"
  )

  resolved <- resolveProtEnrichAnalysisMethod("7227", supported_organisms)

  expect_equal(resolved$analysisMethod, "clusterprofiler")
  expect_false(resolved$organismSupported)
  expect_equal(resolved$methodInfo$method, "clusterprofiler")
  expect_null(resolved$methodInfo$species_id)
  expect_equal(resolved$methodInfo$species_name, "Taxon ID 7227")
  expect_match(
    resolved$methodInfo$description,
    "clusterProfileR analysis with custom GO annotations for taxon 7227",
    perl = TRUE
  )
})

test_that("formatProtEnrichStatusText reports completed enrichment counts", {
  method_info <- list(
    method = "gprofiler2",
    species_name = "Homo sapiens"
  )

  status_text <- formatProtEnrichStatusText(
    analysisComplete = TRUE,
    methodInfo = method_info,
    selectedContrast = "A vs B",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    gprofilerResults = data.frame(term = c("go:1", "go:2"), stringsAsFactors = FALSE),
    clusterprofilerResults = data.frame(term = "cp:1", stringsAsFactors = FALSE),
    stringdbResults = data.frame(network = c("ppi:1", "ppi:2", "ppi:3"), stringsAsFactors = FALSE)
  )

  expect_match(status_text, "\\[OK\\] Analysis Complete", perl = TRUE)
  expect_match(status_text, "Method: gprofiler2", perl = TRUE)
  expect_match(status_text, "Contrast: A vs B", perl = TRUE)
  expect_match(status_text, "\\* gprofiler2: 2 terms", perl = TRUE)
  expect_match(status_text, "\\* clusterProfileR: 1 terms", perl = TRUE)
  expect_match(status_text, "\\* STRING-DB: 3 networks", perl = TRUE)
})

test_that("formatProtEnrichStatusText reports waiting instructions before analysis", {
  status_text <- formatProtEnrichStatusText(analysisComplete = FALSE)

  expect_match(status_text, "\\[WAITING\\] Ready for analysis", perl = TRUE)
  expect_match(status_text, "1\\. Select contrast from DE results", perl = TRUE)
  expect_match(status_text, "4\\. Click 'Run Enrichment Analysis'", perl = TRUE)
  expect_match(status_text, "Method automatically determined by organism\\.", perl = TRUE)
})

test_that("formatProtEnrichAnalysisMethodText reports supported organisms", {
  display_text <- formatProtEnrichAnalysisMethodText(list(
    method = "gprofiler2",
    supported = TRUE,
    species_name = "Homo sapiens"
  ))

  expect_match(display_text, "\\[OK\\] SUPPORTED ORGANISM", perl = TRUE)
  expect_match(display_text, "Method: gprofiler2", perl = TRUE)
  expect_match(display_text, "Species: Homo sapiens", perl = TRUE)
  expect_match(display_text, "All enrichment methods available", perl = TRUE)
})

test_that("formatProtEnrichAnalysisMethodText reports custom organisms", {
  display_text <- formatProtEnrichAnalysisMethodText(list(
    method = "clusterprofiler",
    supported = FALSE,
    species_name = "Taxon ID 7227"
  ))

  expect_match(display_text, "\\[WARNING\\] CUSTOM ORGANISM", perl = TRUE)
  expect_match(display_text, "Method: clusterprofiler", perl = TRUE)
  expect_match(display_text, "Organism: Taxon ID 7227", perl = TRUE)
  expect_match(display_text, "Using UniProt GO annotations", perl = TRUE)
})

test_that("formatProtEnrichContrastsText reports available contrasts line by line", {
  display_text <- formatProtEnrichContrastsText(c("A vs B", "C vs D"))

  expect_identical(display_text, "A vs B\nC vs D")
})

test_that("formatProtEnrichContrastsText reports the pre-analysis fallback", {
  display_text <- formatProtEnrichContrastsText(NULL)

  expect_identical(
    display_text,
    "No contrasts available.\nComplete differential expression\nanalysis first."
  )
})

test_that("formatProtEnrichContrastsText preserves empty character output", {
  display_text <- formatProtEnrichContrastsText(character())

  expect_identical(display_text, "")
})

test_that("formatProtEnrichGprofilerSummaryText reports total and directional counts", {
  summary_text <- formatProtEnrichGprofilerSummaryText(
    gprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("positive", "positive", "negative"),
      stringsAsFactors = FALSE
    )
  )

  expect_match(summary_text, "Total enrichment terms: 3", perl = TRUE)
  expect_match(summary_text, "Up-regulated pathways: 2", perl = TRUE)
  expect_match(summary_text, "Down-regulated pathways: 1", perl = TRUE)
  expect_match(summary_text, "Results displayed in table below\\.", perl = TRUE)
})

test_that("formatProtEnrichGprofilerSummaryText reports filtered up-regulated counts", {
  summary_text <- formatProtEnrichGprofilerSummaryText(
    gprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("positive", "negative", "positive"),
      stringsAsFactors = FALSE
    ),
    directionFilter = "up"
  )

  expect_match(summary_text, "Showing 2 up-regulated pathways", perl = TRUE)
  expect_match(summary_text, "Results displayed in table below\\.", perl = TRUE)
})

test_that("formatProtEnrichGprofilerSummaryText reports empty results before analysis", {
  summary_text <- formatProtEnrichGprofilerSummaryText(
    gprofilerResults = data.frame()
  )

  expect_identical(summary_text, "No gprofiler2 results available.")
})

test_that("renderProtEnrichGprofilerResultsTable returns placeholder table before analysis", {
  rendered <- renderProtEnrichGprofilerResultsTable(
    gprofilerResults = data.frame(),
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      c(x, list(formattedColumns = columns, digits = digits))
    }
  )

  expect_equal(rendered$data$Message, "No gprofiler2 results available. Run analysis first.")
  expect_null(rendered$options)
  expect_null(rendered$extensions)
})

test_that("renderProtEnrichGprofilerResultsTable filters directional rows and tracks rounded columns", {
  rendered <- renderProtEnrichGprofilerResultsTable(
    gprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("positive", "negative", "positive"),
      pvalue = c(0.00123, 0.02, 0.00045),
      qvalue = c(0.005, 0.04, 0.002),
      stringsAsFactors = FALSE
    ),
    directionFilter = "up",
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      c(x, list(formattedColumns = columns, digits = digits))
    }
  )

  expect_equal(rendered$data$term, c("go:1", "go:3"))
  expect_true(all(rendered$data$directionality == "positive"))
  expect_equal(rendered$options$pageLength, 25)
  expect_equal(rendered$options$buttons, c("copy", "csv", "excel"))
  expect_equal(rendered$extensions, "Buttons")
  expect_equal(rendered$formattedColumns, c("pvalue", "qvalue"))
  expect_equal(rendered$digits, 4)
})

test_that("renderProtEnrichGprofilerResultsTable surfaces table-building errors", {
  rendered <- renderProtEnrichGprofilerResultsTable(
    gprofilerResults = data.frame(
      term = "go:1",
      pvalue = 0.01,
      stringsAsFactors = FALSE
    ),
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      stop("format failure")
    },
    catFn = function(...) NULL
  )

  expect_match(rendered$data$Message, "^Error: format failure$", perl = TRUE)
})

test_that("renderProtEnrichGprofilerPlot returns the pre-analysis placeholder shell", {
  rendered <- renderProtEnrichGprofilerPlot(
    analysisComplete = FALSE,
    enrichmentResultsFull = NULL,
    rawContrast = "contrast_a",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_equal(rendered$plot$kind, "plot")
  expect_equal(rendered$text, "Run enrichment analysis first")
  expect_false(rendered$showlegend)
})

test_that("renderProtEnrichGprofilerPlot routes direction-filtered plots and all-mode fallback", {
  up_plot <- list(id = "up")
  down_plot <- list(id = "down")
  enrichment_results_full <- methods::new(
    "mockProtEnrichPlotContainer",
    enrichment_plotly = list(
      contrast_a = list(up = up_plot, down = down_plot),
      contrast_b = list(up = up_plot, down = NULL)
    )
  )

  rendered_down <- renderProtEnrichGprofilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_a",
    directionFilter = "down",
    renderPlotlyFn = function(expr) force(expr)
  )
  rendered_all <- renderProtEnrichGprofilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_a",
    directionFilter = "all",
    renderPlotlyFn = function(expr) force(expr)
  )
  rendered_missing_down <- renderProtEnrichGprofilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_b",
    directionFilter = "down",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_identical(rendered_down, down_plot)
  expect_identical(rendered_all, up_plot)
  expect_equal(rendered_missing_down$text, "No down-regulated enrichment data")
  expect_false(rendered_missing_down$showlegend)
})

test_that("renderProtEnrichGprofilerPlot surfaces plot-building errors", {
  rendered <- renderProtEnrichGprofilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = structure(list(), class = "not_s4"),
    rawContrast = "contrast_a",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_equal(rendered$plot$kind, "plot")
  expect_match(rendered$text, "^Plot error:", perl = TRUE)
  expect_false(rendered$showlegend)
})

test_that("formatProtEnrichClusterProfilerSummaryText reports total and directional counts", {
  summary_text <- formatProtEnrichClusterProfilerSummaryText(
    clusterprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("up", "up", "down"),
      stringsAsFactors = FALSE
    )
  )

  expect_match(summary_text, "Total GO terms: 3", perl = TRUE)
  expect_match(summary_text, "Up-regulated: 2", perl = TRUE)
  expect_match(summary_text, "Down-regulated: 1", perl = TRUE)
  expect_match(summary_text, "Results displayed in table below\\.", perl = TRUE)
})

test_that("formatProtEnrichClusterProfilerSummaryText reports filtered down-regulated counts", {
  summary_text <- formatProtEnrichClusterProfilerSummaryText(
    clusterprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("up", "down", "down"),
      stringsAsFactors = FALSE
    ),
    directionFilter = "down"
  )

  expect_match(summary_text, "Showing 2 down-regulated GO terms", perl = TRUE)
  expect_match(summary_text, "Results displayed in table below\\.", perl = TRUE)
})

test_that("formatProtEnrichClusterProfilerSummaryText reports empty results before analysis", {
  summary_text <- formatProtEnrichClusterProfilerSummaryText(
    clusterprofilerResults = data.frame()
  )

  expect_identical(summary_text, "No clusterProfileR results available.")
})

test_that("renderProtEnrichClusterProfilerResultsTable returns placeholder table before analysis", {
  rendered <- renderProtEnrichClusterProfilerResultsTable(
    clusterprofilerResults = data.frame(),
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      c(x, list(formattedColumns = columns, digits = digits))
    }
  )

  expect_equal(rendered$data$Message, "No clusterProfileR results available. Run analysis first.")
  expect_null(rendered$options)
  expect_null(rendered$extensions)
})

test_that("renderProtEnrichClusterProfilerResultsTable filters directional rows and tracks rounded columns", {
  rendered <- renderProtEnrichClusterProfilerResultsTable(
    clusterprofilerResults = data.frame(
      term = c("go:1", "go:2", "go:3"),
      directionality = c("up", "down", "up"),
      pvalue = c(0.00123, 0.02, 0.00045),
      qvalue = c(0.005, 0.04, 0.002),
      stringsAsFactors = FALSE
    ),
    directionFilter = "up",
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      c(x, list(formattedColumns = columns, digits = digits))
    }
  )

  expect_equal(rendered$data$term, c("go:1", "go:3"))
  expect_true(all(rendered$data$directionality == "up"))
  expect_equal(rendered$options$pageLength, 25)
  expect_equal(rendered$options$buttons, c("copy", "csv", "excel"))
  expect_equal(rendered$extensions, "Buttons")
  expect_equal(rendered$formattedColumns, c("pvalue", "qvalue"))
  expect_equal(rendered$digits, 4)
})

test_that("renderProtEnrichClusterProfilerResultsTable surfaces table-building errors", {
  rendered <- renderProtEnrichClusterProfilerResultsTable(
    clusterprofilerResults = data.frame(
      term = "go:1",
      pvalue = 0.01,
      stringsAsFactors = FALSE
    ),
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    },
    formatRoundFn = function(x, columns, digits) {
      stop("format failure")
    },
    catFn = function(...) NULL
  )

  expect_match(rendered$data$Message, "^Error: format failure$", perl = TRUE)
})

test_that("renderProtEnrichClusterProfilerPlot returns the pre-analysis placeholder shell", {
  rendered <- renderProtEnrichClusterProfilerPlot(
    analysisComplete = FALSE,
    enrichmentResultsFull = NULL,
    rawContrast = "contrast_a",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_equal(rendered$plot$kind, "plot")
  expect_equal(rendered$text, "Run enrichment analysis first")
  expect_false(rendered$showlegend)
})

test_that("renderProtEnrichClusterProfilerPlot routes direction-filtered plots and all-mode fallback", {
  up_plot <- list(id = "up")
  down_plot <- list(id = "down")
  enrichment_results_full <- methods::new(
    "mockProtEnrichPlotContainer",
    enrichment_plotly = list(
      contrast_a = list(up = up_plot, down = down_plot),
      contrast_b = list(up = up_plot, down = NULL)
    )
  )

  rendered_down <- renderProtEnrichClusterProfilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_a",
    directionFilter = "down",
    renderPlotlyFn = function(expr) force(expr)
  )
  rendered_all <- renderProtEnrichClusterProfilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_a",
    directionFilter = "all",
    renderPlotlyFn = function(expr) force(expr)
  )
  rendered_missing_down <- renderProtEnrichClusterProfilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = enrichment_results_full,
    rawContrast = "contrast_b",
    directionFilter = "down",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_identical(rendered_down, down_plot)
  expect_identical(rendered_all, up_plot)
  expect_equal(rendered_missing_down$text, "No down-regulated enrichment data")
  expect_false(rendered_missing_down$showlegend)
})

test_that("renderProtEnrichClusterProfilerPlot surfaces plot-building errors", {
  rendered <- renderProtEnrichClusterProfilerPlot(
    analysisComplete = TRUE,
    enrichmentResultsFull = structure(list(), class = "not_s4"),
    rawContrast = "contrast_a",
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_equal(rendered$plot$kind, "plot")
  expect_match(rendered$text, "^Plot error:", perl = TRUE)
  expect_false(rendered$showlegend)
})

test_that("formatProtEnrichStringDbSummaryText reports empty results before analysis", {
  summary_text <- formatProtEnrichStringDbSummaryText(
    stringdbResults = data.frame()
  )

  expect_identical(
    summary_text,
    "STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics."
  )
})

test_that("formatProtEnrichStringDbSummaryText reports placeholder feature summary", {
  summary_text <- formatProtEnrichStringDbSummaryText(
    stringdbResults = data.frame(network = "ppi:1", stringsAsFactors = FALSE)
  )

  expect_match(summary_text, "STRING-DB Network Analysis", perl = TRUE)
  expect_match(summary_text, "Status: Implementation pending", perl = TRUE)
  expect_match(summary_text, "\\* Protein-protein interaction networks", perl = TRUE)
  expect_match(summary_text, "\\* Interactive network visualization", perl = TRUE)
})

test_that("renderProtEnrichStringDbResultsTable returns placeholder table before analysis", {
  rendered <- renderProtEnrichStringDbResultsTable(
    stringdbResults = data.frame(),
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    }
  )

  expect_equal(rendered$data$Message, "STRING-DB analysis not yet implemented.")
  expect_equal(rendered$data$Status, "Coming soon...")
  expect_null(rendered$options)
  expect_null(rendered$extensions)
})

test_that("renderProtEnrichStringDbResultsTable filters significant rows and limits display count", {
  rendered <- renderProtEnrichStringDbResultsTable(
    stringdbResults = data.frame(
      network = c("ppi:1", "ppi:2", "ppi:3"),
      p_value = c(0.01, 0.20, 0.03),
      stringsAsFactors = FALSE
    ),
    filterSignificant = TRUE,
    enrichmentPValThresh = 0.05,
    maxResults = 1,
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    }
  )

  expect_equal(rendered$data$network, "ppi:1")
  expect_equal(rendered$data$p_value, 0.01)
  expect_equal(rendered$options$pageLength, 25)
  expect_equal(rendered$options$buttons, c("copy", "csv", "excel"))
  expect_equal(rendered$extensions, "Buttons")
})

test_that("renderProtEnrichStringDbResultsTable surfaces table-building errors", {
  rendered <- renderProtEnrichStringDbResultsTable(
    stringdbResults = data.frame(network = "ppi:1", stringsAsFactors = FALSE),
    filterSignificant = TRUE,
    enrichmentPValThresh = 0.05,
    renderDtFn = function(expr) force(expr),
    datatableFn = function(data, options = NULL, extensions = NULL) {
      list(data = data, options = options, extensions = extensions)
    }
  )

  expect_match(rendered$data$Message, "^STRING-DB Error:", perl = TRUE)
})

test_that("renderProtEnrichStringDbPlot returns the STRING-DB placeholder plot shell", {
  rendered <- renderProtEnrichStringDbPlot(
    renderPlotlyFn = function(expr) force(expr),
    plotLyFn = function() list(kind = "plot"),
    addTextFn = function(plot, x, y, text, showlegend) {
      list(
        plot = plot,
        x = x,
        y = y,
        text = text,
        showlegend = showlegend
      )
    }
  )

  expect_equal(rendered$plot$kind, "plot")
  expect_equal(rendered$x, 0.5)
  expect_equal(rendered$y, 0.5)
  expect_equal(rendered$text, "STRING-DB coming soon")
  expect_false(rendered$showlegend)
})

test_that("buildProtEnrichResultsDownloadFilename sanitizes contrast names", {
  filename <- buildProtEnrichResultsDownloadFilename(
    selectedContrast = "Tumor vs Normal / batch#1",
    date = as.Date("2026-04-14")
  )

  expect_identical(
    filename,
    "Enrichment_results_Tumor_vs_Normal___batch_1_2026-04-14.zip"
  )
})

test_that("buildProtEnrichResultsDownloadFilename keeps allowed filename characters", {
  filename <- buildProtEnrichResultsDownloadFilename(
    selectedContrast = "contrast.A-1_B",
    date = as.Date("2026-04-14"),
    prefix = "CustomPrefix"
  )

  expect_identical(
    filename,
    "CustomPrefix_contrast.A-1_B_2026-04-14.zip"
  )
})

test_that("writeProtEnrichResultsDownloadArchive exports result files and summary", {
  temp_dir <- file.path(tempdir(), "prot-enrich-download-success")
  unlink(temp_dir, recursive = TRUE, force = TRUE)
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  zip_call <- NULL
  result <- writeProtEnrichResultsDownloadArchive(
    file = file.path(temp_dir, "results.zip"),
    selectedContrast = "A vs B",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    gprofilerResults = data.frame(term = "go:1", stringsAsFactors = FALSE),
    clusterprofilerResults = data.frame(term = "cp:1", stringsAsFactors = FALSE),
    stringdbResults = data.frame(network = "ppi:1", stringsAsFactors = FALSE),
    tempDir = temp_dir,
    writeTsvFn = function(data, path) {
      utils::write.table(data, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
    },
    zipFn = function(zipfile, files, flags) {
      zip_call <<- list(zipfile = zipfile, files = files, flags = flags)
    },
    sysTimeFn = function() as.POSIXct("2026-04-14 12:34:56", tz = "UTC")
  )

  expect_identical(result$status, "ok")
  expect_identical(
    sort(basename(result$files)),
    sort(c(
      "gprofiler2_results.tsv",
      "clusterProfileR_results.tsv",
      "stringdb_results.tsv",
      "enrichment_analysis_summary.txt"
    ))
  )
  expect_identical(zip_call$zipfile, file.path(temp_dir, "results.zip"))
  expect_identical(sort(basename(zip_call$files)), sort(basename(result$files)))
  expect_identical(zip_call$flags, "-j")

  summary_text <- paste(
    readLines(file.path(temp_dir, "enrichment_analysis_summary.txt")),
    collapse = "\n"
  )
  expect_match(summary_text, "Date: 2026-04-14 12:34:56", perl = TRUE)
  expect_match(summary_text, "Contrast: A vs B", perl = TRUE)
  expect_match(summary_text, "Analysis Method: gprofiler2", perl = TRUE)
  expect_match(summary_text, "gprofiler2 terms: 1", perl = TRUE)
  expect_match(summary_text, "clusterProfileR terms: 1", perl = TRUE)
  expect_match(summary_text, "STRING-DB networks: 1", perl = TRUE)
})

test_that("writeProtEnrichResultsDownloadArchive falls back to an error note", {
  temp_dir <- file.path(tempdir(), "prot-enrich-download-error")
  unlink(temp_dir, recursive = TRUE, force = TRUE)
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  zip_call <- NULL
  result <- writeProtEnrichResultsDownloadArchive(
    file = file.path(temp_dir, "results.zip"),
    selectedContrast = "A vs B",
    methodInfo = list(method = "gprofiler2", species_name = "Homo sapiens"),
    organismTaxid = "9606",
    upCutoff = 1,
    downCutoff = -1,
    qCutoff = 0.05,
    gprofilerResults = data.frame(term = "go:1", stringsAsFactors = FALSE),
    tempDir = temp_dir,
    writeTsvFn = function(data, path) {
      stop("kaboom")
    },
    zipFn = function(zipfile, files, flags) {
      zip_call <<- list(zipfile = zipfile, files = files, flags = flags)
    }
  )

  expect_identical(result$status, "error")
  expect_match(result$error, "kaboom", perl = TRUE)
  expect_identical(basename(result$files), "download_error.txt")
  expect_identical(zip_call$zipfile, file.path(temp_dir, "results.zip"))
  expect_identical(zip_call$flags, "-j")
  expect_identical(basename(zip_call$files), "download_error.txt")
  expect_true(file.exists(file.path(temp_dir, "download_error.txt")))
  expect_identical(
    readLines(file.path(temp_dir, "download_error.txt")),
    "Error creating download: kaboom"
  )
})
