library(testthat)
library(shiny)

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

makeLipidDaCurrentS4 <- function(label = "lipid_da_fixture") {
  methods::new(
    "LipidomicsAssayData",
    lipid_data = list(
      `Positive Mode` = data.frame(
        database_identifier = c(paste0(label, "_l1"), paste0(label, "_l2")),
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
    args = list(
      label = label,
      daAnalysisParameters = list(formula_string = "~ 0 + group")
    )
  )
}

makeLipidDaResults <- function() {
  data.frame(
    database_identifier = c("L1", "L2", "L3"),
    lipid_identification = c("Annot 1", "Annot 2", "Annot 3"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupB-groupA"),
    friendly_name = c("B vs A", "B vs A", "B vs A"),
    assay = c("Positive Mode", "Positive Mode", "Negative Mode"),
    logFC = c(1.2, -1.1, 0.5),
    raw_pvalue = c(0.001, 0.004, 0.02),
    fdr_qvalue = c(0.01, 0.02, 0.03),
    significant = c("Up", "Down", "Up"),
    stringsAsFactors = FALSE
  )
}

makeLipidDaSessionData <- function(current_s4 = makeLipidDaCurrentS4()) {
  list(
    current_s4_object = current_s4,
    r6_current_state_name = "filtered_ready",
    contrasts_tbl = data.frame(
      friendly_names = "B vs A",
      contrasts = "groupB-groupA",
      stringsAsFactors = FALSE
    ),
    assay_names = c("Positive Mode", "Negative Mode")
  )
}

makeLipidDaHarness <- function(session_file_mode = c("valid", "corrupt", "missing", "missing_source")) {
  session_file_mode <- match.arg(session_file_mode)

  capture <- new.env(parent = emptyenv())
  capture$notifications <- list()
  capture$removed_notifications <- character()
  capture$update_select_calls <- list()
  capture$update_textarea_calls <- list()
  capture$saved_states <- list()
  capture$analysis_calls <- list()
  capture$glimma_calls <- list()
  capture$heatmap_calls <- list()
  capture$heatmap_saves <- list()
  capture$output_results_calls <- list()

  root_dir <- tempfile("lipid-da-characterization-")
  dir.create(root_dir, recursive = TRUE)
  da_output_dir <- file.path(root_dir, "da_output")
  publication_graphs_dir <- file.path(root_dir, "publication_graphs")
  dir.create(da_output_dir, recursive = TRUE)
  dir.create(publication_graphs_dir, recursive = TRUE)

  session_file <- file.path(root_dir, "lipid_filtered_session_data_latest.rds")
  if (identical(session_file_mode, "valid")) {
    saveRDS(makeLipidDaSessionData(), session_file)
  } else if (identical(session_file_mode, "corrupt")) {
    writeLines("not an rds payload", session_file)
  }

  state_manager <- new.env(parent = emptyenv())
  state_manager$getState <- function(...) NULL
  state_manager$saveState <- function(...) {
    capture$saved_states[[length(capture$saved_states) + 1L]] <<- list(...)
    invisible(NULL)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager
  workflow_data$contrasts_tbl <- NULL
  workflow_data$tab_status <- list(differential_analysis = "pending", other = "keep")

  source_dir <- root_dir
  export_dir <- root_dir
  if (identical(session_file_mode, "missing_source")) {
    source_dir <- file.path(root_dir, "missing-source")
    export_dir <- file.path(root_dir, "missing-export")
  }

  list(
    capture = capture,
    workflow_data = workflow_data,
    experiment_paths = list(
      source_dir = source_dir,
      export_dir = export_dir,
      da_output_dir = da_output_dir,
      publication_graphs_dir = publication_graphs_dir
    ),
    session = shiny::MockShinySession$new(),
    cleanup = function() unlink(root_dir, recursive = TRUE, force = TRUE)
  )
}

makeLipidDaCharacterizationServer <- function() {
  server_under_test <- mod_lipid_da_server
  compatibility_env <- list2env(
    list(
      generateMetabVolcanoPlotGlimma = function(...) {
        get("generateLipidDAVolcanoPlotGlimma", envir = asNamespace("MultiScholaR"))(...)
      }
    ),
    parent = environment(mod_lipid_da_server)
  )
  environment(server_under_test) <- compatibility_env
  server_under_test
}

launchLipidDaModule <- function(harness, server_under_test = makeLipidDaCharacterizationServer()) {
  shiny::withReactiveDomain(harness$session, {
    server_under_test(
      id = "da",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "lipidomics",
      experiment_label = "Lipidomics"
    )
  })

  harness$session$getReturned()
}

setLipidDaInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("da-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

renderLipidDaOutput <- function(session, output_id) {
  htmltools::renderTags(session$getOutput(paste0("da-", output_id)))$html
}

notificationText <- function(notification) {
  paste(as.character(notification$ui), collapse = "")
}

installLipidDaModuleMocks <- function(
  harness,
  run_analysis_fn = function(...) {
    harness$capture$analysis_calls[[length(harness$capture$analysis_calls) + 1L]] <<- list(...)
    list(
      da_lipids_long = makeLipidDaResults(),
      significant_counts = list(`Positive Mode` = list(up = 1L, down = 1L, ns = 0L))
    )
  },
  glimma_fn = function(...) {
    harness$capture$glimma_calls[[length(harness$capture$glimma_calls) + 1L]] <<- list(...)
    htmltools::div("glimma widget")
  },
  heatmap_fn = function(...) {
    harness$capture$heatmap_calls[[length(harness$capture$heatmap_calls) + 1L]] <<- list(...)
    list(
      plot = ggplot2::ggplot(),
      row_clusters = c(L1 = 1, L2 = 2),
      col_clusters = NULL
    )
  },
  output_results_fn = function(...) {
    harness$capture$output_results_calls[[length(harness$capture$output_results_calls) + 1L]] <<- list(...)
    TRUE
  },
  save_heatmap_fn = function(...) {
    harness$capture$heatmap_saves[[length(harness$capture$heatmap_saves) + 1L]] <<- list(...)
    invisible(NULL)
  }
) {
  local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL) {
      harness$capture$update_select_calls[[length(harness$capture$update_select_calls) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextAreaInput = function(session, inputId, value = NULL, label = NULL, placeholder = NULL) {
      harness$capture$update_textarea_calls[[length(harness$capture$update_textarea_calls) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        value = value,
        label = label,
        placeholder = placeholder
      )
      invisible(NULL)
    },
    showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL, type = "default", session = getDefaultReactiveDomain()) {
      harness$capture$notifications[[length(harness$capture$notifications) + 1L]] <<- list(
        ui = ui,
        id = id,
        type = type,
        duration = duration
      )
      invisible(if (is.null(id)) paste0("notification-", length(harness$capture$notifications)) else id)
    },
    removeNotification = function(id, session = getDefaultReactiveDomain()) {
      harness$capture$removed_notifications <<- c(harness$capture$removed_notifications, id)
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )

  package_bindings <- list(
    runLipidsDA = run_analysis_fn,
    generateLipidDAVolcanoPlotGlimma = glimma_fn,
    generateLipidDAHeatmap = heatmap_fn,
    outputLipidDaResultsAllContrasts = output_results_fn,
    save_heatmap_products = save_heatmap_fn,
    generateLipidDAVolcanoStatic = function(...) ggplot2::ggplot()
  )

  do.call(
    local_mocked_bindings,
    c(package_bindings, list(.env = asNamespace("MultiScholaR")))
  )
}

test_that("mod_lipid_da_server preserves load-session public success behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "valid")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(harness)

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)

  expect_identical(harness$workflow_data$contrasts_tbl$friendly_names, "B vs A")
  expect_identical(harness$capture$saved_states[[1]]$state_name, "filtered_ready")
  expect_identical(
    harness$capture$saved_states[[1]]$description,
    "Loaded from filtered session for DA analysis"
  )
})

test_that("mod_lipid_da_server preserves missing-source-dir public error behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "missing_source")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(harness)

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)

  expect_length(harness$capture$saved_states, 0)
  expect_length(harness$capture$update_select_calls, 0)
  expect_identical(harness$workflow_data$contrasts_tbl, NULL)
})

test_that("mod_lipid_da_server preserves missing-session-file public error behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "missing")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(harness)

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)

  expect_length(harness$capture$saved_states, 0)
  expect_length(harness$capture$update_select_calls, 0)
  expect_identical(harness$workflow_data$contrasts_tbl, NULL)
})

test_that("mod_lipid_da_server preserves corrupt-session public error behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "corrupt")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(harness)

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)

  expect_length(harness$capture$saved_states, 0)
  expect_length(harness$capture$update_select_calls, 0)
  expect_identical(harness$workflow_data$contrasts_tbl, NULL)
})

test_that("mod_lipid_da_server preserves run-analysis public success behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "valid")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(harness)

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)
  setLipidDaInputs(
    harness$session,
    formula_string = "~ 0 + group",
    da_q_val_thresh = 0.01,
    treat_lfc_cutoff = 1.5,
    run_da_analysis = 1
  )

  expect_identical(harness$workflow_data$tab_status$differential_analysis, "complete")
  expect_identical(harness$capture$analysis_calls[[1]]$formula_string, "~ 0 + group")
  expect_identical(harness$capture$analysis_calls[[1]]$da_q_val_thresh, 0.01)
  expect_identical(harness$capture$analysis_calls[[1]]$treat_lfc_cutoff, 1.5)
  expect_length(harness$capture$output_results_calls, 1)
})

test_that("mod_lipid_da_server preserves run-analysis public error behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "valid")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(
    harness,
    run_analysis_fn = function(...) {
      harness$capture$analysis_calls[[length(harness$capture$analysis_calls) + 1L]] <<- list(...)
      stop("model matrix singular")
    }
  )

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)
  setLipidDaInputs(
    harness$session,
    formula_string = "~ 0 + group",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0,
    run_da_analysis = 1
  )

  expect_identical(harness$workflow_data$tab_status$differential_analysis, "pending")
  expect_length(harness$capture$output_results_calls, 0)
  expect_length(harness$capture$analysis_calls, 1)
})

test_that("mod_lipid_da_server preserves render and heatmap-save public behavior", {
  harness <- makeLipidDaHarness(session_file_mode = "valid")
  on.exit(harness$cleanup(), add = TRUE)
  installLipidDaModuleMocks(
    harness,
    glimma_fn = function(...) {
      harness$capture$glimma_calls[[length(harness$capture$glimma_calls) + 1L]] <<- list(...)
      stop("plot backend failed")
    }
  )

  launchLipidDaModule(harness)
  setLipidDaInputs(harness$session, load_filtered_session = 1)
  setLipidDaInputs(
    harness$session,
    formula_string = "~ 0 + group",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0,
    run_da_analysis = 1
  )
  setLipidDaInputs(
    harness$session,
    volcano_contrast = "B vs A",
    volcano_assay = "Positive Mode",
    heatmap_contrast = "B vs A",
    heatmap_assay = "Positive Mode",
    heatmap_top_n = 25,
    heatmap_cluster_method = "average",
    heatmap_distance_method = "manhattan",
    heatmap_clustering = "column",
    heatmap_scaling = "row",
    heatmap_color_scheme = "Viridis",
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 4,
    heatmap_cut_height = 0.8,
    heatmap_min_cluster_size = 3,
    table_contrast = "B vs A",
    table_assay = "All",
    table_significance = "all",
    table_max_rows = 50
  )

  glimma_html <- renderLipidDaOutput(harness$session, "volcano_glimma")
  expect_match(glimma_html, "Error generating plot: plot backend failed", fixed = TRUE)

  expect_no_error(harness$session$getOutput("da-heatmap_plot"))
  cluster_summary_html <- renderLipidDaOutput(harness$session, "cluster_summary")
  summary_stats_html <- renderLipidDaOutput(harness$session, "da_summary_stats")
  results_table_html <- renderLipidDaOutput(harness$session, "da_results_table")

  expect_true(nzchar(cluster_summary_html))
  expect_true(nzchar(summary_stats_html))
  expect_true(nzchar(results_table_html))

  setLipidDaInputs(harness$session, save_heatmap = 1)

  expect_length(harness$capture$glimma_calls, 1)
  expect_length(harness$capture$heatmap_calls, 1)
  expect_length(harness$capture$heatmap_saves, 1)
  expect_s3_class(harness$capture$heatmap_saves[[1]]$heatmap_obj, "ggplot")
  expect_identical(harness$capture$heatmap_saves[[1]]$file_prefix, "lipid_B_vs_A")
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$top_n, 25)
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$cluster_cols, TRUE)
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$cluster_rows, FALSE)
})
