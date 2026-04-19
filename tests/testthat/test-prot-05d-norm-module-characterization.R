library(testthat)
library(shiny)

if (!methods::isClass("mockProtNormCharacterizationData")) {
  methods::setClass(
    "mockProtNormCharacterizationData",
    slots = c(
      protein_quant_table = "data.frame",
      protein_id_column = "character",
      args = "list"
    )
  )
}

makeProtNormCharacterizationData <- function(label = "prot_norm_fixture") {
  methods::new(
    "mockProtNormCharacterizationData",
    protein_quant_table = data.frame(
      Protein.Ids = c("p1", "p2"),
      S1 = c(10, 20),
      S2 = c(30, 40),
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    args = list(
      label = label,
      globalParameters = list(workflow_type = "DIA")
    )
  )
}

makeProtNormCharacterizationHarness <- function(
  current_state = "ruv_corrected",
  current_s4 = makeProtNormCharacterizationData()
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()
  capture$reverted_states <- character()

  root_dir <- tempfile("prot-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  protein_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(protein_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_state
  state_manager$getState <- function(state) current_s4
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1]] <<- list(...)
  }
  state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "protein_replicate_filtered", current_state)
  }
  state_manager$revertToState <- function(state) {
    capture$reverted_states <<- c(capture$reverted_states, state)
    current_s4
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$design_matrix <- NULL
  workflow_data$state_manager <- state_manager
  workflow_data$tab_status <- list(normalization = "pending", differential_expression = "disabled")
  workflow_data$state_update_trigger <- NULL
  workflow_data$protein_counts <- list()
  workflow_data$ruv_optimization_result <- list(ruv_skipped = FALSE)
  workflow_data$contrasts_tbl <- data.frame(friendly_names = "A vs B", stringsAsFactors = FALSE)
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "DIA"))
  workflow_data$fasta_metadata <- list(database = "ref")
  workflow_data$accession_cleanup_results <- list(clean = TRUE)
  workflow_data$qc_params <- list(minimum = 1)
  workflow_data$mixed_species_analysis <- NULL

  experiment_paths <- list(
    protein_qc_dir = protein_qc_dir,
    source_dir = source_dir
  )

  list(
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    session = shiny::MockShinySession$new(),
    current_s4 = current_s4,
    capture = capture
  )
}

launchProtNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_prot_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "proteomics",
      experiment_label = "Proteomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

setProtNormInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("norm-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("mod_prot_norm_server preserves apply-correlation public behavior", {
  harness <- makeProtNormCharacterizationHarness()
  final_s4 <- makeProtNormCharacterizationData(label = "final")

  local_mocked_bindings(
    pearsonCorForSamplePairs = function(...) {
      data.frame(sample1 = "S1", sample2 = "S2", pearson = 0.99, stringsAsFactors = FALSE)
    },
    filterSamplesByProteinCorrelationThreshold = function(...) final_s4,
    updateProteinFiltering = function(...) "filter_plot",
    plotPca = function(...) "pca_plot",
    .env = asNamespace("MultiScholaR")
  )

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- harness$current_s4
  state$ruv_optimization_result <- list(ruv_skipped = FALSE)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.8
  )
  setProtNormInputs(harness$session, apply_correlation_filter = 1)

  expect_true(shiny::isolate(state$correlation_filtering_complete))
  expect_identical(shiny::isolate(state$correlation_filtered_obj), final_s4)
  expect_identical(harness$workflow_data$ruv_normalised_for_da_analysis_obj, final_s4)
  expect_equal(harness$workflow_data$tab_status$normalization, "complete")
  expect_equal(harness$workflow_data$tab_status$differential_expression, "pending")
  expect_length(harness$capture$saved_states, 1)
  expect_equal(harness$capture$saved_states[[1]]$state_name, "correlation_filtered")
})

test_that("mod_prot_norm_server preserves skip-correlation public behavior", {
  harness <- makeProtNormCharacterizationHarness()

  local_mocked_bindings(
    updateProteinFiltering = function(...) "filter_plot",
    plotPca = function(...) "pca_plot",
    .env = asNamespace("MultiScholaR")
  )

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- harness$current_s4
  state$ruv_optimization_result <- list(ruv_skipped = FALSE)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group"
  )
  setProtNormInputs(harness$session, skip_correlation_filter = 1)

  expect_true(shiny::isolate(state$correlation_filtering_complete))
  expect_identical(shiny::isolate(state$correlation_filtered_obj), harness$current_s4)
  expect_identical(harness$workflow_data$ruv_normalised_for_da_analysis_obj, harness$current_s4)
  expect_equal(harness$workflow_data$tab_status$differential_expression, "pending")
  expect_length(harness$capture$saved_states, 1)
  expect_true(isTRUE(harness$capture$saved_states[[1]]$config_object$skipped))
})

test_that("mod_prot_norm_server preserves export prereq and success behavior", {
  harness <- makeProtNormCharacterizationHarness(current_state = "correlation_filtered")

  state <- launchProtNormModule(harness)

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    export_filtered_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))

  state$correlation_filtering_complete <- TRUE
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$correlation_threshold <- 0.75

  setProtNormInputs(harness$session, export_filtered_session = 2)

  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))
  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_summary.txt")))
})

test_that("mod_prot_norm_server preserves reset behavior", {
  harness <- makeProtNormCharacterizationHarness()

  state <- launchProtNormModule(harness)
  state$normalization_complete <- TRUE
  state$ruv_complete <- TRUE
  state$correlation_filtering_complete <- TRUE
  state$normalized_protein_obj <- harness$current_s4
  state$ruv_normalized_obj <- harness$current_s4
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$control_genes_index <- c(TRUE, FALSE)
  state$correlation_vector <- data.frame(value = 0.95)
  state$correlation_threshold <- 0.8
  state$final_qc_plot <- "pca_plot"
  state$final_filtering_plot <- "filter_plot"
  state$post_norm_filtering_plot <- "post_plot"
  state$filtering_summary_text <- "summary"
  state$ruv_optimization_result <- list(best_k = 2)

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    reset_normalization = 1
  )

  expect_false(shiny::isolate(state$normalization_complete))
  expect_false(shiny::isolate(state$ruv_complete))
  expect_false(shiny::isolate(state$correlation_filtering_complete))
  expect_null(shiny::isolate(state$normalized_protein_obj))
  expect_null(shiny::isolate(state$ruv_normalized_obj))
  expect_null(shiny::isolate(state$correlation_filtered_obj))
  expect_equal(harness$workflow_data$tab_status$normalization, "pending")
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
  expect_equal(harness$capture$reverted_states[[1]], "protein_replicate_filtered")
})

test_that("mod_prot_norm_server preserves invalid-tab warning behavior", {
  harness <- makeProtNormCharacterizationHarness(current_state = "intensity_filtered")
  selected_tab_value <- shiny::reactiveVal("overview")

  state <- launchProtNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("normalization")
  shiny:::flushReact()

  expect_false(shiny::isolate(state$pre_norm_qc_generated))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})

test_that("mod_prot_norm_server preserves pre-normalization error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "protein_replicate_filtered")
  selected_tab_value <- shiny::reactiveVal("overview")

  local_mocked_bindings(
    plotPca = function(...) stop("pre-qc boom"),
    .env = asNamespace("MultiScholaR")
  )

  state <- launchProtNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("normalization")
  shiny:::flushReact()

  expect_false(shiny::isolate(state$pre_norm_qc_generated))
  expect_false(file.exists(file.path(harness$experiment_paths$protein_qc_dir, "pre_norm_pca.png")))
})

test_that("mod_prot_norm_server preserves normalization error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "protein_replicate_filtered")

  local_mocked_bindings(
    normaliseBetweenSamples = function(...) stop("normalization boom"),
    .env = asNamespace("MultiScholaR")
  )

  state <- launchProtNormModule(harness)

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    run_normalization = 1
  )

  expect_false(shiny::isolate(state$normalization_complete))
  expect_false(shiny::isolate(state$ruv_complete))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})

test_that("mod_prot_norm_server preserves correlation error handling", {
  harness <- makeProtNormCharacterizationHarness()

  state <- launchProtNormModule(harness)
  state$ruv_normalized_obj <- NULL

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.8,
    apply_correlation_filter = 1
  )

  expect_false(shiny::isolate(state$correlation_filtering_complete))
  expect_length(harness$capture$saved_states, 0)
})

test_that("mod_prot_norm_server preserves export error handling", {
  harness <- makeProtNormCharacterizationHarness(current_state = "correlation_filtered")
  harness$experiment_paths$source_dir <- file.path(harness$experiment_paths$source_dir, "missing")

  state <- launchProtNormModule(harness)
  state$correlation_filtering_complete <- TRUE
  state$correlation_filtered_obj <- harness$current_s4
  state$best_k <- 2
  state$correlation_threshold <- 0.75

  setProtNormInputs(
    harness$session,
    norm_method = "cyclicloess",
    ruv_mode = "automatic",
    export_filtered_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "filtered_session_data_latest.rds")))
})

test_that("mod_prot_norm_server preserves reset error handling", {
  harness <- makeProtNormCharacterizationHarness()
  harness$workflow_data$state_manager$getHistory <- function() stop("reset boom")

  state <- launchProtNormModule(harness)
  state$normalization_complete <- TRUE
  state$ruv_complete <- TRUE

  setProtNormInputs(
    harness$session,
    ruv_mode = "automatic",
    reset_normalization = 1
  )

  expect_true(shiny::isolate(state$normalization_complete))
  expect_true(shiny::isolate(state$ruv_complete))
  expect_equal(harness$workflow_data$tab_status$differential_expression, "disabled")
})
