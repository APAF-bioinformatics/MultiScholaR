library(testthat)
library(shiny)
library(htmltools)

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

makeMetabDaCurrentS4 <- function(label = "metab_da_fixture") {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = data.frame(
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
    args = list(
      label = label,
      daAnalysisParameters = list(formula_string = "~ 0 + group")
    )
  )
}

makeMetabDaResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3"),
    metabolite_name = c("Met One", "Met Two", "Met Three"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupB-groupA"),
    friendly_name = c("B vs A", "B vs A", "B vs A"),
    assay = c("LCMS_Pos", "LCMS_Pos", "LCMS_Neg"),
    logFC = c(1.2, -1.1, 0.5),
    raw_pvalue = c(0.001, 0.004, 0.02),
    fdr_qvalue = c(0.01, 0.02, 0.03),
    significant = c("Up", "Down", "Up"),
    stringsAsFactors = FALSE
  )
}

makeMetabDaSessionData <- function(current_s4 = makeMetabDaCurrentS4()) {
  list(
    current_s4_object = current_s4,
    r6_current_state_name = "filtered_ready",
    contrasts_tbl = data.frame(
      friendly_names = "B vs A",
      contrasts = "groupB-groupA",
      stringsAsFactors = FALSE
    ),
    assay_names = c("LCMS_Pos", "LCMS_Neg")
  )
}

makeMetabDaHarness <- function(session_file_mode = c("valid", "corrupt", "missing")) {
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

  root_dir <- tempfile("metab-da-characterization-")
  dir.create(root_dir, recursive = TRUE)
  da_output_dir <- file.path(root_dir, "da_output")
  publication_graphs_dir <- file.path(root_dir, "publication_graphs")
  dir.create(da_output_dir, recursive = TRUE)
  dir.create(publication_graphs_dir, recursive = TRUE)

  session_file <- file.path(root_dir, "metab_filtered_session_data_latest.rds")
  if (identical(session_file_mode, "valid")) {
    saveRDS(makeMetabDaSessionData(), session_file)
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

  list(
    capture = capture,
    workflow_data = workflow_data,
    experiment_paths = list(
      source_dir = root_dir,
      export_dir = root_dir,
      da_output_dir = da_output_dir,
      publication_graphs_dir = publication_graphs_dir
    ),
    session = shiny::MockShinySession$new()
  )
}

launchMetabDaModule <- function(harness) {
  shiny::withReactiveDomain(harness$session, {
    mod_metab_da_server(
      id = "da",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      omic_type = "metabolomics",
      experiment_label = "Metabolomics"
    )
  })

  harness$session$getReturned()
}

setMetabDaInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("da-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

renderMetabDaOutput <- function(session, output_id) {
  htmltools::renderTags(session$getOutput(paste0("da-", output_id)))$html
}

installMetabDaModuleMocks <- function(
  harness,
  run_analysis_fn = function(...) {
    harness$capture$analysis_calls[[length(harness$capture$analysis_calls) + 1L]] <<- list(...)
    list(
      da_metabolites_long = makeMetabDaResults(),
      significant_counts = list(LCMS_Pos = list(up = 1L, down = 1L, ns = 0L))
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
      row_clusters = c(M1 = 1, M2 = 2),
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

  local_mocked_bindings(
    runMetabolitesDA = run_analysis_fn,
    generateMetabDAVolcanoPlotGlimma = glimma_fn,
    generateMetabDAHeatmap = heatmap_fn,
    outputMetabDaResultsAllContrasts = output_results_fn,
    save_heatmap_products = save_heatmap_fn,
    generateMetabDAVolcanoStatic = function(...) ggplot2::ggplot(),
    .env = asNamespace("MultiScholaR")
  )
}

test_that("mod_metab_da_server preserves load-session public success behavior", {
  harness <- makeMetabDaHarness(session_file_mode = "valid")
  installMetabDaModuleMocks(harness)

  launchMetabDaModule(harness)
  setMetabDaInputs(harness$session, load_filtered_session = 1)

  expect_identical(harness$workflow_data$contrasts_tbl$friendly_names, "B vs A")
  expect_identical(harness$capture$saved_states[[1]]$state_name, "filtered_ready")
  expect_identical(
    harness$capture$saved_states[[1]]$description,
    "Loaded from filtered session for DA analysis"
  )
})

test_that("mod_metab_da_server preserves load-session public error behavior", {
  harness <- makeMetabDaHarness(session_file_mode = "corrupt")
  installMetabDaModuleMocks(harness)

  launchMetabDaModule(harness)
  setMetabDaInputs(harness$session, load_filtered_session = 1)

  expect_length(harness$capture$saved_states, 0)
  expect_length(harness$capture$update_select_calls, 0)
  expect_identical(harness$workflow_data$contrasts_tbl, NULL)
})

test_that("mod_metab_da_server preserves run-analysis public success behavior", {
  harness <- makeMetabDaHarness(session_file_mode = "valid")
  installMetabDaModuleMocks(harness)

  launchMetabDaModule(harness)
  setMetabDaInputs(harness$session, load_filtered_session = 1)
  setMetabDaInputs(
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

test_that("mod_metab_da_server preserves run-analysis public error behavior", {
  harness <- makeMetabDaHarness(session_file_mode = "valid")
  installMetabDaModuleMocks(
    harness,
    run_analysis_fn = function(...) {
      harness$capture$analysis_calls[[length(harness$capture$analysis_calls) + 1L]] <<- list(...)
      stop("model matrix singular")
    }
  )

  launchMetabDaModule(harness)
  setMetabDaInputs(harness$session, load_filtered_session = 1)
  setMetabDaInputs(
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

test_that("mod_metab_da_server preserves Glimma error and heatmap-save public behavior", {
  harness <- makeMetabDaHarness(session_file_mode = "valid")
  installMetabDaModuleMocks(
    harness,
    glimma_fn = function(...) {
      harness$capture$glimma_calls[[length(harness$capture$glimma_calls) + 1L]] <<- list(...)
      stop("plot backend failed")
    }
  )

  launchMetabDaModule(harness)
  setMetabDaInputs(harness$session, load_filtered_session = 1)
  setMetabDaInputs(
    harness$session,
    formula_string = "~ 0 + group",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0,
    run_da_analysis = 1
  )

  setMetabDaInputs(
    harness$session,
    volcano_contrast = "B vs A",
    volcano_assay = "LCMS_Pos"
  )
  glimma_html <- renderMetabDaOutput(harness$session, "volcano_glimma")
  expect_match(glimma_html, "Error generating plot: plot backend failed", fixed = TRUE)

  setMetabDaInputs(
    harness$session,
    heatmap_contrast = "B vs A",
    heatmap_assay = "LCMS_Pos",
    heatmap_top_n = 25,
    heatmap_cluster_method = "average",
    heatmap_distance_method = "manhattan",
    heatmap_clustering = "column",
    heatmap_scaling = "row",
    heatmap_color_scheme = "Viridis",
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 4,
    heatmap_cut_height = 0.8,
    heatmap_min_cluster_size = 3
  )
  expect_no_error(harness$session$getOutput("da-heatmap_plot"))

  setMetabDaInputs(harness$session, save_heatmap = 1)

  expect_length(harness$capture$glimma_calls, 1)
  expect_length(harness$capture$heatmap_calls, 1)
  expect_length(harness$capture$heatmap_saves, 1)
  expect_s3_class(harness$capture$heatmap_saves[[1]]$heatmap_obj, "ggplot")
  expect_identical(harness$capture$heatmap_saves[[1]]$file_prefix, "metab_B_vs_A")
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$top_n, 25)
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$cluster_cols, TRUE)
  expect_identical(harness$capture$heatmap_saves[[1]]$params_list$cluster_rows, FALSE)
})
