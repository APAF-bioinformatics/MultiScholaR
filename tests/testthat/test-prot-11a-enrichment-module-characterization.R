library(testthat)
library(shiny)

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

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_true(isTRUE(state$analysis_complete))
      expect_equal(state$gprofiler_results$term_name, "GO:1")
      expect_equal(state$da_results_data$contrast_a$theObject@args$label, "run_fixture")
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

      state <- shiny::isolate(reactiveValuesToList(enrichment_data, all.names = TRUE))

      expect_false(isTRUE(state$analysis_complete))
      expect_null(state$gprofiler_results)
      expect_equal(state$da_results_data$contrast_a$theObject@args$label, "run_failure_fixture")
    }
  )
})
