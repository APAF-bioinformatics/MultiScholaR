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

skipIfMissingMultiScholaRBindings <- function(...) {
  names <- unlist(list(...), use.names = FALSE)
  missing <- names[!vapply(names, hasMultiScholaRBinding, logical(1))]
  if (length(missing) > 0) {
    skip(paste("requires extracted helper bindings:", paste(missing, collapse = ", ")))
  }
}

if (!methods::isClass("LipidomicsAssayData")) {
  methods::setClass(
    "LipidomicsAssayData",
    slots = c(
      lipid_data = "list",
      lipid_id_column = "character",
      annotation_id_column = "character",
      database_identifier_type = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      lipid_data = list(),
      lipid_id_column = "database_identifier",
      annotation_id_column = "lipid_identification",
      database_identifier_type = "Unknown",
      internal_standard_regex = NA_character_,
      design_matrix = data.frame(),
      sample_id = "Sample_ID",
      group_id = "group",
      technical_replicate_id = NA_character_,
      args = list()
    )
  )
}

makeLipidNormCharacterizationData <- function() {
  methods::new(
    "LipidomicsAssayData",
    lipid_data = list(
      `Positive Mode` = data.frame(
        database_identifier = c("L1", "L2"),
        lipid_identification = c("A1", "A2"),
        Sample1 = c(10, 20),
        Sample2 = c(30, 40),
        stringsAsFactors = FALSE
      )
    ),
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
    database_identifier_type = "Mock",
    internal_standard_regex = "",
    design_matrix = data.frame(
      Sample_ID = c("Sample1", "Sample2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

makeLipidNormCharacterizationHarness <- function(current_s4 = makeLipidNormCharacterizationData()) {
  capture <- new.env(parent = emptyenv())
  capture$saved_states <- list()

  root_dir <- tempfile("lipid-norm-characterization-")
  dir.create(root_dir, recursive = TRUE)
  lipid_qc_dir <- file.path(root_dir, "qc")
  source_dir <- file.path(root_dir, "source")
  dir.create(lipid_qc_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "lipid_normalized"
  state_manager$getState <- function(state = state_manager$current_state) current_s4
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1]] <<- list(...)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$design_matrix <- data.frame(
    Sample_ID = c("Sample1", "Sample2"),
    group = c("A", "B"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )
  workflow_data$config_list <- list(globalParameters = list(workflow_type = "lipidomics"))
  workflow_data$contrasts_tbl <- data.frame(friendly_names = "A vs B", stringsAsFactors = FALSE)
  workflow_data$lipid_counts <- list(`Positive Mode` = list(features = 2, samples = 2))
  workflow_data$qc_params <- list(minimum_intensity = 1)
  workflow_data$tab_status <- list(quality_control = "pending", normalization = "pending")

  list(
    workflow_data = workflow_data,
    experiment_paths = list(lipid_qc_dir = lipid_qc_dir, source_dir = source_dir, export_dir = source_dir),
    session = shiny::MockShinySession$new(),
    current_s4 = current_s4,
    capture = capture
  )
}

launchLipidNormModule <- function(harness, selected_tab = NULL) {
  shiny::withReactiveDomain(harness$session, {
    mod_lipid_norm_server(
      id = "norm",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics",
      selected_tab = selected_tab
    )
  })

  harness$session$getReturned()
}

renderLipidNormOutput <- function(session, output_id) {
  htmltools::renderTags(session$getOutput(paste0("norm-", output_id)))$html
}

test_that("mod_lipid_norm_ui preserves the public normalization shell", {
  ui <- mod_lipid_norm_ui("norm")
  html <- htmltools::renderTags(ui)$html

  expect_match(html, "Normalization Options", fixed = TRUE)
  expect_match(html, "Run Normalization Pipeline", fixed = TRUE)
  expect_match(html, "Correlation Filtering", fixed = TRUE)
})

test_that("mod_lipid_norm_server preserves startup log placeholder behavior", {
  harness <- makeLipidNormCharacterizationHarness()

  local_mocked_bindings(
    generateLipidQcPlots = function(...) invisible(NULL),
    .env = asNamespace("MultiScholaR")
  )

  launchLipidNormModule(harness)

  norm_log_html <- renderLipidNormOutput(harness$session, "norm_log")

  expect_match(
    norm_log_html,
    "Normalization log will appear here as you apply steps...",
    fixed = TRUE
  )
})

test_that("mod_lipid_norm_server preserves selected-tab pre-QC error logging behavior", {
  harness <- makeLipidNormCharacterizationHarness()
  selected_tab_value <- shiny::reactiveVal("overview")

  local_mocked_bindings(
    generateLipidQcPlots = function(...) stop("pre qc boom"),
    .env = asNamespace("MultiScholaR")
  )

  launchLipidNormModule(
    harness,
    selected_tab = function() selected_tab_value()
  )

  selected_tab_value("norm")
  shiny:::flushReact()

  norm_log_html <- renderLipidNormOutput(harness$session, "norm_log")

  expect_match(norm_log_html, "Error generating Pre-QC: pre qc boom", fixed = TRUE)
})

test_that("lipid norm UI helpers preserve the generated control and QC shell", {
  skipIfMissingMultiScholaRBindings(
    "buildLipidNormOptionsControlPanel",
    "buildLipidNormQcTabsetPanel"
  )
  skip_if_not_installed("shinyjqui")

  build_options <- getMultiScholaRBinding("buildLipidNormOptionsControlPanel")
  build_tabs <- getMultiScholaRBinding("buildLipidNormQcTabsetPanel")

  html <- htmltools::renderTags(shiny::tagList(
    build_options(shiny::NS("norm")),
    build_tabs(shiny::NS("norm"))
  ))$html

  expect_match(html, "Normalization Options", fixed = TRUE)
  expect_match(html, "ITSD Selection", fixed = TRUE)
  expect_match(html, "Correlation Filtering", fixed = TRUE)
  expect_match(html, "Export Normalized Data", fixed = TRUE)
})

test_that("lipid norm runtime helpers wire startup and server orchestration layers", {
  skipIfMissingMultiScholaRBindings(
    "createLipidNormReactiveState",
    "createLipidNormStartupRuntime",
    "registerLipidNormPrimaryStartupOutputs",
    "registerLipidNormItsdSelectionRuntime",
    "registerLipidNormPostNormalizationOutputs",
    "registerLipidNormStartupObserverRuntime",
    "registerLipidNormServerRuntime",
    "runLipidNormModuleServerShell",
    "runLipidNormModuleServerEntryShell",
    "runLipidNormModuleServerPublicWrapper"
  )

  create_state <- getMultiScholaRBinding("createLipidNormReactiveState")
  create_startup <- getMultiScholaRBinding("createLipidNormStartupRuntime")
  register_primary <- getMultiScholaRBinding("registerLipidNormPrimaryStartupOutputs")
  register_itsd_runtime <- getMultiScholaRBinding("registerLipidNormItsdSelectionRuntime")
  register_post_outputs <- getMultiScholaRBinding("registerLipidNormPostNormalizationOutputs")
  register_startup_observers <- getMultiScholaRBinding("registerLipidNormStartupObserverRuntime")
  register_server_runtime <- getMultiScholaRBinding("registerLipidNormServerRuntime")
  run_shell <- getMultiScholaRBinding("runLipidNormModuleServerShell")
  run_entry_shell <- getMultiScholaRBinding("runLipidNormModuleServerEntryShell")
  run_public <- getMultiScholaRBinding("runLipidNormModuleServerPublicWrapper")

  norm_data <- create_state(reactiveValuesFn = function(...) {
    structure(list(...), class = "mock_reactive_values")
  })
  expect_identical(norm_data$normalization_log, character(0))
  expect_false(norm_data$normalization_complete)

  startup <- create_startup(
    input = list(color_variable = "batch", shape_variable = "group"),
    experimentPaths = list(lipid_qc_dir = tempdir()),
    normData = norm_data,
    ns = identity,
    buildAddLogFn = function(normData) function(message) {
      normData$normalization_log <- c(normData$normalization_log, message)
      invisible(normData$normalization_log)
    },
    buildPlotAestheticsGetterFn = function(input) function() {
      list(color_var = input$color_variable, shape_var = input$shape_variable)
    },
    buildCompositeFromFilesGeneratorFn = function() function(...) "composite",
    buildQcImageRendererFn = function(...) function(...) "image",
    buildLogRendererFn = function(...) function() "log",
    buildItsdSelectionUiRendererFn = function(...) function() "itsd",
    buildRuvQcUiRendererFn = function(...) function() "ruv",
    buildAssayLabelRendererFn = function(...) function(...) "assay",
    buildCorrelationFilterSummaryRendererFn = function(...) function() "summary",
    buildFinalQcPlotRendererFn = function(...) function() "plot"
  )
  expect_equal(sort(names(startup)), sort(c(
    "addLog", "getPlotAesthetics", "generateCompositeFromFiles",
    "renderQcImageForAssay", "renderNormLog", "renderItsdSelectionUi",
    "renderRuvQcUi", "renderAssayLabel", "renderCorrelationFilterSummary",
    "renderFinalQcPlot"
  )))
  expect_identical(startup$getPlotAesthetics(), list(color_var = "batch", shape_var = "group"))

  output <- new.env(parent = emptyenv())
  calls <- character()
  record_call <- function(name) {
    calls <<- c(calls, name)
  }

  register_primary(
    output = output,
    startupRuntime = startup,
    registerLogOutputFn = function(...) record_call("log"),
    registerItsdSelectionOutputFn = function(...) record_call("itsd_output"),
    registerRuvQcOutputFn = function(...) record_call("ruv_output"),
    registerStaticQcImageOutputsFn = function(...) record_call("static_qc"),
    registerAssayLabelOutputsFn = function(...) record_call("assay_labels")
  )
  register_itsd_runtime(
    input = list(),
    output = output,
    workflowData = new.env(parent = emptyenv()),
    normData = norm_data,
    registerItsdTableOutputsFn = function(...) record_call("itsd_tables"),
    registerItsdSelectionTrackingFn = function(...) record_call("itsd_tracking")
  )
  register_post_outputs(
    output = output,
    normData = norm_data,
    startupRuntime = startup,
    registerRuvCancorOutputsFn = function(...) record_call("ruv_cancor"),
    registerCorrelationFilterSummaryOutputFn = function(...) record_call("correlation_summary"),
    registerFinalQcPlotOutputFn = function(...) record_call("final_qc")
  )
  expect_setequal(
    calls,
    c(
      "log", "itsd_output", "ruv_output", "static_qc", "assay_labels",
      "itsd_tables", "itsd_tracking", "ruv_cancor", "correlation_summary",
      "final_qc"
    )
  )

  selected_tab <- function() "norm"
  startup_calls <- character()
  register_startup_observers(
    session = shiny::MockShinySession$new(),
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(lipid_qc_dir = tempdir()),
    normData = norm_data,
    addLog = startup$addLog,
    getPlotAestheticsFn = startup$getPlotAesthetics,
    selectedTab = selected_tab,
    registerAssayNameInitializationObserverFn = function(...) startup_calls <<- c(startup_calls, "assay"),
    registerSelectedTabPreNormalizationObserverFn = function(...) startup_calls <<- c(startup_calls, "selected_tab"),
    registerDesignDrivenChoiceObserverFn = function(...) startup_calls <<- c(startup_calls, "design")
  )
  expect_identical(startup_calls, c("assay", "selected_tab", "design"))

  server_calls <- character()
  register_server_runtime(
    input = list(),
    output = output,
    session = shiny::MockShinySession$new(),
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(lipid_qc_dir = tempdir()),
    omicType = "lipidomics",
    experimentLabel = "Lipidomics",
    normData = norm_data,
    startupRuntime = startup,
    addLog = startup$addLog,
    getPlotAestheticsFn = startup$getPlotAesthetics,
    generateCompositeFromFilesFn = startup$generateCompositeFromFiles,
    selectedTab = selected_tab,
    registerStartupObserverRuntimeFn = function(...) server_calls <<- c(server_calls, "startup"),
    registerPrimaryStartupOutputsFn = function(...) server_calls <<- c(server_calls, "primary"),
    registerItsdSelectionRuntimeFn = function(...) server_calls <<- c(server_calls, "itsd"),
    registerRunNormalizationObserverFn = function(...) server_calls <<- c(server_calls, "run"),
    registerResetNormalizationObserverFn = function(...) server_calls <<- c(server_calls, "reset"),
    registerPostNormalizationOutputsFn = function(...) server_calls <<- c(server_calls, "post"),
    registerApplyCorrelationFilterObserverFn = function(...) server_calls <<- c(server_calls, "apply"),
    registerSkipCorrelationFilterObserverFn = function(...) server_calls <<- c(server_calls, "skip"),
    registerExportSessionObserverFn = function(...) server_calls <<- c(server_calls, "export")
  )
  expect_identical(server_calls, c("startup", "primary", "itsd", "run", "reset", "post", "apply", "skip", "export"))

  shell_calls <- character()
  shell_result <- run_shell(
    input = list(),
    output = output,
    session = shiny::MockShinySession$new(),
    id = "norm",
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(lipid_qc_dir = tempdir()),
    omicType = "lipidomics",
    experimentLabel = "Lipidomics",
    logInfoFn = function(...) invisible(NULL),
    createReactiveStateFn = function() {
      shell_calls <<- c(shell_calls, "state")
      norm_data
    },
    createStartupRuntimeFn = function(...) {
      shell_calls <<- c(shell_calls, "startup")
      startup
    },
    registerServerRuntimeFn = function(...) shell_calls <<- c(shell_calls, "runtime")
  )
  expect_identical(shell_result, norm_data)
  expect_identical(shell_calls, c("state", "startup", "runtime"))

  entry_calls <- list()
  run_entry_shell(
    id = "norm",
    workflowData = new.env(parent = emptyenv()),
    experimentPaths = list(lipid_qc_dir = tempdir()),
    omicType = "lipidomics",
    experimentLabel = "Lipidomics",
    moduleServerFn = function(id, module) {
      entry_calls$id <<- id
      module(list(), new.env(parent = emptyenv()), shiny::MockShinySession$new())
    },
    runModuleServerShellFn = function(...) {
      entry_calls$args <<- list(...)
      "entry"
    }
  )
  expect_identical(entry_calls$id, "norm")
  expect_identical(entry_calls$args$id, "norm")

  public_args <- NULL
  run_public(
    id = "norm",
    workflow_data = new.env(parent = emptyenv()),
    experiment_paths = list(lipid_qc_dir = tempdir()),
    omic_type = "lipidomics",
    experiment_label = "Lipidomics",
    selected_tab = selected_tab,
    runModuleServerEntryShellFn = function(...) {
      public_args <<- list(...)
      "public"
    }
  )
  expect_identical(public_args$omicType, "lipidomics")
  expect_identical(public_args$selectedTab, selected_tab)
})

test_that("lipid norm observer helpers preserve assay and pre-QC behavior", {
  skipIfMissingMultiScholaRBindings(
    "registerLipidNormDesignDrivenChoiceObserver",
    "registerLipidNormAssayNameInitializationObserver",
    "registerLipidNormRuvCancorOutputs",
    "registerLipidNormItsdTableOutputs",
    "registerLipidNormItsdSelectionTracking",
    "handleLipidNormPreNormalizationQc",
    "handleLipidNormSelectedTabPreNormalizationTrigger"
  )
  skip_if_not_installed("ggplot2")

  observe_now <- function(expr, ...) force(expr)
  observe_event_now <- function(eventExpr, handlerExpr, ..., ignoreNULL = TRUE, ignoreInit = FALSE) {
    force(eventExpr)
    force(handlerExpr)
  }
  walk_now <- function(.x, .f) invisible(lapply(.x, .f))
  render_now <- function(expr, ...) eval.parent(substitute(expr))

  register_design <- getMultiScholaRBinding("registerLipidNormDesignDrivenChoiceObserver")
  register_assays <- getMultiScholaRBinding("registerLipidNormAssayNameInitializationObserver")
  register_ruv_outputs <- getMultiScholaRBinding("registerLipidNormRuvCancorOutputs")
  register_itsd_tables <- getMultiScholaRBinding("registerLipidNormItsdTableOutputs")
  register_itsd_tracking <- getMultiScholaRBinding("registerLipidNormItsdSelectionTracking")
  handle_pre_qc <- getMultiScholaRBinding("handleLipidNormPreNormalizationQc")
  handle_selected_tab <- getMultiScholaRBinding("handleLipidNormSelectedTabPreNormalizationTrigger")

  updates <- list()
  workflow_data <- makeLipidNormCharacterizationHarness()$workflow_data
  register_design(
    session = shiny::MockShinySession$new(),
    workflowData = workflow_data,
    observeFn = observe_now,
    updateSelectInputFn = function(session, inputId, choices, selected) {
      updates[[inputId]] <<- list(choices = choices, selected = selected)
    }
  )
  expect_identical(updates$color_variable$selected, "group")
  expect_true("batch" %in% updates$ruv_grouping_variable$choices)

  norm_data <- new.env(parent = emptyenv())
  norm_data$itsd_selections <- list()
  info_messages <- character()
  register_assays(
    workflowData = workflow_data,
    normData = norm_data,
    observeFn = observe_now,
    reqFn = force,
    logInfoFn = function(message) info_messages <<- c(info_messages, message),
    logWarnFn = function(message) info_messages <<- c(info_messages, message)
  )
  expect_identical(norm_data$assay_names, "Positive Mode")
  expect_identical(length(norm_data$itsd_selections), 0L)

  output <- new.env(parent = emptyenv())
  norm_data$ruv_optimization_results <- list(
    `Positive Mode` = list(
      success = TRUE,
      cancor_plot = ggplot2::ggplot(),
      best_k = 2,
      best_percentage = 5,
      separation_score = 0.81234,
      control_genes_index = c(TRUE, FALSE),
      optimization_results = data.frame(k = 2, score = 0.8)
    )
  )
  register_ruv_outputs(
    output = output,
    normData = norm_data,
    observeFn = observe_now,
    reqFn = force,
    walkFn = walk_now,
    renderPlotFn = render_now,
    renderTextFn = render_now,
    renderDataTableFn = render_now,
    datatableFn = function(data, ...) list(data = data, options = list(...))
  )
  expect_s3_class(output$cancor_plot_positive_mode, "ggplot")
  expect_match(output$ruv_summary_positive_mode, "Best k: 2", fixed = TRUE)
  expect_equal(output$ruv_table_positive_mode$data$score, 0.8)

  register_itsd_tables(
    output = output,
    workflowData = workflow_data,
    normData = norm_data,
    observeFn = observe_now,
    reqFn = force,
    walkFn = walk_now,
    renderDataTableFn = render_now,
    datatableFn = function(data, ...) structure(list(data = data), class = "mock_datatable"),
    formatStyleFn = function(x, ...) x,
    styleEqualFn = function(...) "style",
    formatRoundFn = function(x, ...) x,
    buildLipidItsdSelectionTableFn = function(...) data.frame(
      feature_id = c("L1", "L2"),
      is_candidate = c(TRUE, FALSE),
      mean_intensity = c(10, 20),
      cv_percent = c(1, 2)
    )
  )
  expect_s3_class(output$itsd_table_positive_mode, "mock_datatable")
  expect_identical(output$itsd_table_positive_mode$data$feature_id, c("L1", "L2"))

  input <- new.env(parent = emptyenv())
  input$itsd_table_positive_mode_rows_selected <- c(1L, 2L)
  register_itsd_tracking(
    input = input,
    normData = norm_data,
    observeFn = observe_now,
    reqFn = force,
    walkFn = walk_now,
    observeEventFn = observe_event_now,
    logInfoFn = function(...) invisible(NULL)
  )
  expect_identical(norm_data$itsd_selections[["Positive Mode"]], c(1L, 2L))

  log_messages <- character()
  norm_data$pre_norm_qc_generated <- FALSE
  norm_data$plot_refresh_trigger <- 0
  pre_qc_result <- handle_pre_qc(
    workflowData = workflow_data,
    experimentPaths = list(lipid_qc_dir = tempdir()),
    normData = norm_data,
    addLog = function(message) log_messages <<- c(log_messages, message),
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    reqFn = force,
    generateLipidQcPlotsFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  expect_true(pre_qc_result)
  expect_true(norm_data$pre_norm_qc_generated)
  expect_equal(norm_data$plot_refresh_trigger, 1)

  norm_data$pre_norm_qc_generated <- FALSE
  selected_result <- handle_selected_tab(
    selectedTabValue = "norm",
    workflowData = workflow_data,
    experimentPaths = list(lipid_qc_dir = tempdir()),
    normData = norm_data,
    addLog = function(message) log_messages <<- c(log_messages, message),
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    reqFn = force,
    withProgressFn = function(message, value, expr) force(expr),
    handlePreNormalizationQcFn = function(...) TRUE,
    logInfoFn = function(...) invisible(NULL)
  )
  expect_true(selected_result)
  expect_false(handle_selected_tab(
    selectedTabValue = "overview",
    workflowData = workflow_data,
    experimentPaths = list(lipid_qc_dir = tempdir()),
    normData = norm_data,
    addLog = function(...) invisible(NULL),
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    reqFn = force,
    withProgressFn = function(message, value, expr) force(expr),
    handlePreNormalizationQcFn = function(...) TRUE,
    logInfoFn = function(...) invisible(NULL)
  ))
})

test_that("lipid norm support helpers preserve output registration and composite behavior", {
  skipIfMissingMultiScholaRBindings(
    "registerLipidNormStaticQcImageOutputs",
    "registerLipidNormAssayLabelOutputs",
    "buildLipidNormCompositeFromFilesGenerator"
  )
  skip_if_not_installed("ggplot2")

  register_static <- getMultiScholaRBinding("registerLipidNormStaticQcImageOutputs")
  register_labels <- getMultiScholaRBinding("registerLipidNormAssayLabelOutputs")
  build_composite <- getMultiScholaRBinding("buildLipidNormCompositeFromFilesGenerator")

  output <- new.env(parent = emptyenv())
  register_static(output, function(assaySlot, plotType, stagePrefix) {
    paste(assaySlot, plotType, stagePrefix, sep = ":")
  })
  expect_identical(output$pca_post_filter_assay1, "1:pca:pre_norm")
  expect_identical(output$correlation_ruv_corrected_assay2, "2:correlation:ruv_corrected")

  register_labels(output, function(assaySlot) paste0("label-", assaySlot))
  expect_identical(output$assay1_label_pca, "label-1")
  expect_identical(output$assay2_label_correlation, "label-2")

  composite <- build_composite(
    requireNamespaceFn = function(package, quietly = TRUE) TRUE,
    warningFn = function(...) invisible(NULL),
    fileExistsFn = function(path) TRUE,
    readPngFn = function(path) array(1, dim = c(1, 1, 4)),
    rasterGrobFn = function(...) grid::nullGrob(),
    wrapPlotsFn = function(plots, ncol = 1) ggplot2::ggplot(),
    plotLayoutFn = function(...) NULL,
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  result <- composite(
    plot_files = c("a.png", "b.png", "c.png", "d.png"),
    ncol = 2,
    row_labels = list(row1 = c("a)", "b)"), row2 = c("c)", "d)")),
    column_labels = c("Pre", "Post")
  )
  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$width, 10)
  expect_true(result$height > 4)
})

test_that("lipid norm workflow helpers preserve normalization, export, reset, and correlation paths", {
  skipIfMissingMultiScholaRBindings(
    "handleLipidNormRunNormalization",
    "handleLipidNormExportSession",
    "handleLipidNormSkipCorrelationFilter",
    "handleLipidNormResetNormalization",
    "handleLipidNormApplyCorrelationFilter"
  )
  skip_if_not_installed("ggplot2")

  run_normalization <- getMultiScholaRBinding("handleLipidNormRunNormalization")
  export_session <- getMultiScholaRBinding("handleLipidNormExportSession")
  skip_correlation <- getMultiScholaRBinding("handleLipidNormSkipCorrelationFilter")
  reset_normalization <- getMultiScholaRBinding("handleLipidNormResetNormalization")
  apply_correlation <- getMultiScholaRBinding("handleLipidNormApplyCorrelationFilter")

  harness <- makeLipidNormCharacterizationHarness()
  norm_data <- new.env(parent = emptyenv())
  norm_data$assay_names <- "Positive Mode"
  norm_data$itsd_selections <- list(`Positive Mode` = 1L)
  norm_data$normalization_complete <- FALSE
  norm_data$ruv_complete <- FALSE
  norm_data$correlation_filtering_complete <- FALSE
  norm_data$ruv_optimization_results <- list()
  norm_data$correlation_results <- list()
  norm_data$normalization_log <- character()
  norm_data$plot_refresh_trigger <- 0

  input <- list(
    apply_itsd = TRUE,
    itsd_aggregation = "median",
    log_offset = 1,
    norm_method = "none",
    ruv_mode = "automatic",
    auto_percentage_min = 1,
    auto_percentage_max = 10,
    ruv_grouping_variable = "group",
    separation_metric = "max_difference",
    k_penalty_weight = 0.5,
    adaptive_k_penalty = TRUE,
    ruv_k = 2,
    ruv_percentage = 5,
    min_pearson_correlation_threshold = 0.5
  )

  log_messages <- character()
  add_log <- function(message) {
    log_messages <<- c(log_messages, message)
    invisible(log_messages)
  }
  imap_now <- function(.x, .f) {
    out <- lapply(names(.x), function(name) .f(.x[[name]], name))
    names(out) <- names(.x)
    out
  }
  compact_now <- function(.x) Filter(Negate(is.null), .x)

  result <- run_normalization(
    input = input,
    workflowData = harness$workflow_data,
    experimentPaths = harness$experiment_paths,
    omicType = "lipidomics",
    normData = norm_data,
    addLog = add_log,
    getPlotAestheticsFn = function() list(color_var = "group", shape_var = "batch"),
    reqFn = force,
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(...) invisible(NULL),
    generateLipidQcPlotsFn = function(...) invisible(NULL),
    buildLipidItsdSelectionTableFn = function(...) data.frame(
      feature_id = c("L1", "L2"),
      is_candidate = c(TRUE, FALSE)
    ),
    compactFn = compact_now,
    imapFn = imap_now,
    normaliseUntransformedDataFn = function(theObject, ...) theObject,
    logTransformAssaysFn = function(theObject, ...) theObject,
    normaliseBetweenSamplesFn = function(theObject, ...) theObject,
    runLipidPerAssayRuvOptimizationFn = function(...) list(
      `Positive Mode` = list(
        success = TRUE,
        best_k = 1,
        best_percentage = 5,
        control_genes_index = c(TRUE, FALSE),
        cancor_plot = ggplot2::ggplot(),
        optimization_results = data.frame(k = 1)
      )
    ),
    extractLipidBestKPerAssayFn = function(results) list(`Positive Mode` = 1),
    extractLipidCtrlPerAssayFn = function(results) list(`Positive Mode` = c(TRUE, FALSE)),
    ruvIII_C_VaryingFn = function(theObject, ...) theObject,
    generateCompositeFromFilesFn = function(...) list(plot = ggplot2::ggplot(), width = 8, height = 6),
    savePlotFn = function(...) invisible(NULL),
    dirExistsFn = function(path) TRUE,
    showNotificationFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  expect_true(result)
  expect_true(norm_data$normalization_complete)
  expect_true(norm_data$ruv_complete)
  expect_true(any(vapply(harness$capture$saved_states, \(x) identical(x$state_name, "lipid_ruv_corrected"), logical(1))))
  expect_true(any(grepl("Normalization pipeline complete", log_messages, fixed = TRUE)))

  exported <- export_session(
    input = input,
    workflowData = harness$workflow_data,
    experimentPaths = harness$experiment_paths,
    experimentLabel = "Lipidomics",
    normData = norm_data,
    addLog = add_log,
    reqFn = force,
    showNotificationFn = function(...) invisible(NULL),
    withProgressFn = function(message, value, expr) force(expr),
    incProgressFn = function(...) invisible(NULL),
    saveRdsFn = function(object, file) invisible(file),
    writeLinesFn = function(text, con) invisible(text),
    getTimeFn = function() as.POSIXct("2026-04-22 07:00:00", tz = "UTC"),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  expect_match(exported$session_filename, "lipid_filtered_session_data_20260422_070000.rds", fixed = TRUE)

  norm_data$ruv_corrected_obj <- harness$current_s4
  skip_result <- skip_correlation(
    workflowData = harness$workflow_data,
    normData = norm_data,
    addLog = add_log,
    reqFn = force,
    showNotificationFn = function(...) invisible(NULL)
  )
  expect_true(skip_result)
  expect_identical(harness$workflow_data$tab_status$normalization, "complete")

  reset_result <- reset_normalization(
    workflowData = harness$workflow_data,
    normData = norm_data,
    addLog = add_log,
    reqFn = force,
    showNotificationFn = function(...) invisible(NULL)
  )
  expect_true(reset_result)
  expect_false(norm_data$normalization_complete)
  expect_false(norm_data$ruv_complete)

  norm_data$normalization_complete <- TRUE
  norm_data$post_norm_obj <- harness$current_s4
  apply_result <- apply_correlation(
    input = input,
    workflowData = harness$workflow_data,
    normData = norm_data,
    addLog = add_log,
    reqFn = force,
    showNotificationFn = function(...) invisible(NULL),
    removeNotificationFn = function(...) invisible(NULL),
    pearsonCorForSamplePairsFn = function(...) list(
      `Positive Mode` = data.frame(pearson_correlation = c(0.8, 0.9))
    ),
    filterSamplesByLipidCorrelationThresholdFn = function(theObject, ...) theObject,
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  )
  expect_true(apply_result)
  expect_true(norm_data$correlation_filtering_complete)
  expect_identical(harness$workflow_data$tab_status$quality_control, "complete")

  incomplete_norm_data <- new.env(parent = emptyenv())
  incomplete_norm_data$normalization_complete <- FALSE
  expect_false(export_session(
    input = input,
    workflowData = harness$workflow_data,
    experimentPaths = harness$experiment_paths,
    experimentLabel = "Lipidomics",
    normData = incomplete_norm_data,
    addLog = add_log,
    reqFn = force,
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    logWarnFn = function(...) invisible(NULL),
    logErrorFn = function(...) invisible(NULL)
  ))
})
