library(testthat)
library(shiny)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id",
      args = list()
    )
  )
}

makeMetabNormCharacterizationData <- function(label = "metab_norm_fixture") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c(paste0(label, "_m1"), paste0(label, "_m2")),
        annotation_id = c("a1", "a2"),
        S1 = c(10, 20),
        S2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    metabolite_id_column = "metabolite_id",
    annotation_id_column = "annotation_id",
    args = list(label = label)
  )
}

makeMetabNormCharacterizationHarness <- function(
  current_state = "metab_post_norm",
  current_s4 = makeMetabNormCharacterizationData()
) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()

  root_dir <- tempfile("metab-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  metabolite_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(metabolite_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- current_state
  state_manager$getState <- function(state = current_state) current_s4
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1L]] <<- list(...)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$design_matrix <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("A", "B"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )
  workflow_data$contrasts_tbl <- data.frame(
    friendly_names = "A vs B",
    stringsAsFactors = FALSE
  )
  workflow_data$config_list <- list(method = "TIC")
  workflow_data$metabolite_counts <- list(Plasma = list(features = 2, samples = 2))
  workflow_data$qc_params <- list(minimum = 1)
  workflow_data$tab_status <- list(
    quality_control = "pending",
    normalization = "pending",
    differential_expression = "locked"
  )

  experiment_paths <- list(
    metabolite_qc_dir = metabolite_qc_dir,
    source_dir = source_dir,
    export_dir = source_dir
  )

  list(
    workflow_data = workflow_data,
    experiment_paths = experiment_paths,
    session = shiny::MockShinySession$new(),
    current_s4 = current_s4,
    capture = capture
  )
}

launchMetabNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_metab_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "metabolomics",
      experiment_label = "Metabolomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

setMetabNormInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("norm-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

test_that("mod_metab_norm_server preserves selected-tab pre-QC hydration behavior", {
  harness <- makeMetabNormCharacterizationHarness(current_state = "metab_import_complete")
  capture <- new.env(parent = emptyenv())
  capture$pre_qc_calls <- 0L
  selected_tab_value <- shiny::reactiveVal("overview")

  local_mocked_bindings(
    generateMetabQcPlots = function(...) {
      capture$pre_qc_calls <- capture$pre_qc_calls + 1L
      invisible(NULL)
    },
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  expect_identical(capture$pre_qc_calls, 0L)
  selected_tab_value("norm")
  shiny:::flushReact()

  expect_identical(capture$pre_qc_calls, 1L)
})

test_that("mod_metab_norm_server preserves normalization observer success behavior", {
  harness <- makeMetabNormCharacterizationHarness()

  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    color_variable = "group",
    shape_variable = "batch",
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_normalized"), logical(1))))
})

test_that("mod_metab_norm_server preserves normalization observer error behavior", {
  harness <- makeMetabNormCharacterizationHarness()

  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(...) stop("normalization boom"),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  expect_false(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_normalized"), logical(1))))
})

test_that("mod_metab_norm_server preserves apply-correlation public behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  filtered_s4 <- makeMetabNormCharacterizationData(label = "filtered")

  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    pearsonCorForSamplePairs = function(...) {
      list(Plasma = data.frame(pearson_correlation = 0.96, stringsAsFactors = FALSE))
    },
    filterSamplesByMetaboliteCorrelationThreshold = function(...) filtered_s4,
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(
    harness$session,
    ruv_grouping_variable = "batch",
    min_pearson_correlation_threshold = 0.85,
    apply_correlation_filter = 1
  )

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_correlation_filtered"), logical(1))))
  expect_identical(
    harness$workflow_data$tab_status,
    list(
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "locked"
    )
  )
})

test_that("mod_metab_norm_server preserves skip-correlation public behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(harness$session, skip_correlation_filter = 1)

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_norm_complete"), logical(1))))
  expect_identical(
    harness$workflow_data$tab_status,
    list(
      quality_control = "complete",
      normalization = "complete",
      differential_expression = "locked"
    )
  )
})

test_that("mod_metab_norm_server preserves export prereq and success behavior", {
  harness <- makeMetabNormCharacterizationHarness(current_state = "metab_norm_complete")
  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )
  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    itsd_aggregation = "median",
    log_offset = 1,
    min_pearson_correlation_threshold = 0.8,
    ruv_grouping_variable = "group",
    export_session = 1
  )

  expect_false(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_data_latest.rds")))

  setMetabNormInputs(harness$session, run_normalization = 1)
  setMetabNormInputs(harness$session, skip_correlation_filter = 1)

  setMetabNormInputs(harness$session, export_session = 2)

  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_data_latest.rds")))
  expect_true(file.exists(file.path(harness$experiment_paths$source_dir, "metab_filtered_session_summary.txt")))
})

test_that("mod_metab_norm_server preserves reset behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )
  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )

  setMetabNormInputs(harness$session, reset_normalization = 1)

  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "metab_reset"), logical(1))))
})

test_that("mod_metab_norm_server preserves reset error behavior", {
  harness <- makeMetabNormCharacterizationHarness()
  local_mocked_bindings(
    generateMetabQcPlots = function(...) invisible(NULL),
    logTransformAssays = function(theObject, ...) theObject,
    savePlot = function(...) invisible(NULL),
    plotPca = function(...) list(Plasma = ggplot2::ggplot()),
    .env = asNamespace("MultiScholaR")
  )

  launchMetabNormModule(harness)

  setMetabNormInputs(
    harness$session,
    apply_itsd = FALSE,
    norm_method = "none",
    ruv_mode = "skip",
    run_normalization = 1
  )
  harness$workflow_data$state_manager$saveState <- function(...) stop("save failed")
  pre_reset_save_count <- length(harness$capture$saved_states)

  setMetabNormInputs(harness$session, reset_normalization = 1)

  expect_identical(length(harness$capture$saved_states), pre_reset_save_count)
})
