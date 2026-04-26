# fidelity-coverage-compare: shared
library(testthat)

msr <- function(name) {
  get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

da_server_heatmap_render_handler <- msr("da_server_heatmap_render_handler")
da_server_table_render_handler <- msr("da_server_table_render_handler")

if (!methods::isClass("FakeProtDaRenderObject")) {
  methods::setClass("FakeProtDaRenderObject", slots = c(protein_id_column = "character"))
}

newProtDaResultsEnv <- function() {
  da_data <- new.env(parent = emptyenv())
  da_data$da_results_list <- list(
    da_proteins_long = data.frame(
      comparison = c("A/B vs C", "A/B vs C", "A/B vs C", "D vs E"),
      Protein.Ids = c("P1", "P2", "P3", "P4"),
      log2FC = c(2.5, -2.2, 0.3, 1.1),
      raw_pvalue = c(0.001, 0.002, 0.5, 0.02),
      fdr_qvalue = c(0.01, 0.02, 0.6, 0.04),
      stringsAsFactors = FALSE
    ),
    theObject = methods::new("FakeProtDaRenderObject", protein_id_column = "Protein.Ids")
  )
  da_data$current_row_clusters <- NULL
  da_data$current_col_clusters <- NULL
  da_data$current_heatmap_plot <- NULL
  da_data
}

localShinyRenderMocks <- function(captured, .local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    observeEvent = function(eventExpr, handlerExpr, ..., ignoreInit = FALSE) {
      value <- eval(substitute(eventExpr), parent.frame())
      trigger <- TRUE

      if (is.list(value)) {
        trigger <- all(vapply(value, function(x) !is.null(x), logical(1)))
      } else if (is.null(value) || identical(value, FALSE)) {
        trigger <- FALSE
      }

      if (trigger) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value) || identical(value, FALSE)) {
          stop("required value missing", call. = FALSE)
        }
      }
      invisible(values[[1]])
    },
    renderPlot = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPrint = function(expr) {
      paste(capture.output(eval(substitute(expr), parent.frame())), collapse = "\n")
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    showNotification = function(message, type = NULL, duration = NULL, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    downloadHandler = function(filename, content, ...) {
      list(filename = filename, content = content)
    },
    .package = "shiny",
    .env = .local_envir
  )
}

localDtRenderMocks <- function(captured, .local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    datatable = function(data, options = NULL, extensions = NULL, ...) {
      captured$datatable_calls[[length(captured$datatable_calls) + 1L]] <- list(
        data = data,
        options = options,
        extensions = extensions
      )
      structure(
        list(data = data, options = options, extensions = extensions),
        class = "mock_datatable"
      )
    },
    formatRound = function(table, columns, digits, ...) {
      table$format_round <- list(columns = columns, digits = digits)
      table
    },
    .package = "DT",
    .env = .local_envir
  )
}

test_that("proteomics DA heatmap render handler preserves checkpoint capture, plot handoff, and save flow", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$checkpoints <- list()
  captured$heatmap_args <- NULL
  captured$save_args <- NULL
  captured$logs <- character()

  localShinyRenderMocks(captured)
  testthat::local_mocked_bindings(
    .capture_checkpoint = function(payload, checkpoint_id, label) {
      captured$checkpoints[[length(captured$checkpoints) + 1L]] <- list(
        payload = payload,
        checkpoint_id = checkpoint_id,
        label = label
      )
      invisible(NULL)
    },
    generateProtDAHeatmap = function(...) {
      captured$heatmap_args <- list(...)
      list(
        plot = "heatmap-plot",
        row_clusters = c(P1 = 1, P2 = 2),
        col_clusters = c(S1 = 1, S2 = 1)
      )
    },
    save_heatmap_products = function(...) {
      captured$save_args <- list(...)
      invisible(NULL)
    },
    .package = "MultiScholaR"
  )
  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$logs <- c(captured$logs, message)
      invisible(NULL)
    },
    .package = "logger"
  )

  input <- list(
    heatmap_contrast = "A/B vs C",
    heatmap_top_n = 25,
    heatmap_cluster_method = "average",
    heatmap_distance_method = "euclidean",
    heatmap_clustering = "both",
    heatmap_scaling = "row",
    heatmap_color_scheme = "Viridis",
    heatmap_show_labels = TRUE,
    da_q_val_thresh = 0.05,
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 3,
    heatmap_cut_height = 0.8,
    heatmap_min_cluster_size = 4,
    save_heatmap = 1
  )
  output <- new.env(parent = emptyenv())
  da_data <- newProtDaResultsEnv()
  experiment_paths <- list(publication_graphs_dir = tempdir())

  da_server_heatmap_render_handler(
    input = input,
    output = output,
    session = list(),
    ns = function(id) id,
    da_data = da_data,
    experiment_paths = experiment_paths
  )

  expect_identical(output$heatmap_plot, "heatmap-plot")
  expect_match(output$cluster_summary, "Total Clusters: 2", fixed = TRUE)
  expect_match(output$cluster_summary, "Cluster 1", fixed = TRUE)
  expect_identical(da_data$current_row_clusters, c(P1 = 1, P2 = 2))
  expect_identical(da_data$current_col_clusters, c(S1 = 1, S2 = 1))
  expect_identical(da_data$current_heatmap_plot, "heatmap-plot")
  expect_identical(captured$checkpoints[[1L]]$checkpoint_id, "cp09")
  expect_identical(captured$checkpoints[[1L]]$label, "heatmap_input")
  expect_true(isTRUE(captured$checkpoints[[1L]]$payload$cluster_rows))
  expect_true(isTRUE(captured$checkpoints[[1L]]$payload$cluster_cols))
  expect_identical(captured$heatmap_args$selected_contrast, "A/B vs C")
  expect_identical(captured$heatmap_args$top_n_genes, 25)
  expect_identical(captured$save_args$heatmap_obj, "heatmap-plot")
  expect_identical(captured$save_args$output_dir, tempdir())
  expect_identical(captured$save_args$file_prefix, "prot_A_B_vs_C")
  expect_identical(captured$notifications[[1L]]$message, "Heatmap and cluster info saved to publication_graphs/Heatmap")
  expect_true(any(grepl("Save Heatmap button clicked", captured$logs, fixed = TRUE)))
})

test_that("proteomics DA heatmap render handler preserves legacy passthrough results", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$checkpoints <- list()

  localShinyRenderMocks(captured)
  testthat::local_mocked_bindings(
    .capture_checkpoint = function(payload, checkpoint_id, label) {
      captured$checkpoints[[length(captured$checkpoints) + 1L]] <- list(
        payload = payload,
        checkpoint_id = checkpoint_id,
        label = label
      )
      invisible(NULL)
    },
    generateProtDAHeatmap = function(...) "legacy-heatmap-plot",
    .package = "MultiScholaR"
  )

  input <- list(
    heatmap_contrast = "A/B vs C",
    heatmap_top_n = 10,
    heatmap_cluster_method = "complete",
    heatmap_distance_method = "manhattan",
    heatmap_clustering = "column",
    heatmap_scaling = "none",
    heatmap_color_scheme = "Plasma",
    heatmap_show_labels = FALSE,
    da_q_val_thresh = 0.01,
    heatmap_tree_cut_method = "dynamic",
    heatmap_n_clusters = 2,
    heatmap_cut_height = 0.6,
    heatmap_min_cluster_size = 3,
    save_heatmap = NULL
  )
  output <- new.env(parent = emptyenv())
  da_data <- newProtDaResultsEnv()

  da_server_heatmap_render_handler(
    input = input,
    output = output,
    session = list(),
    ns = function(id) id,
    da_data = da_data,
    experiment_paths = list(publication_graphs_dir = tempdir())
  )

  expect_identical(output$heatmap_plot, "legacy-heatmap-plot")
  expect_match(output$cluster_summary, "No clusters defined", fixed = TRUE)
  expect_identical(da_data$current_row_clusters, NULL)
  expect_identical(captured$notifications, list())
})

test_that("proteomics DA table render handler preserves significance filtering, summaries, and downloads", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$datatable_calls <- list()

  localShinyRenderMocks(captured)
  localDtRenderMocks(captured)

  input <- list(
    table_contrast = "A/B vs C",
    table_significance = "up",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 1,
    table_max_rows = 1
  )
  output <- new.env(parent = emptyenv())
  da_data <- newProtDaResultsEnv()

  da_server_table_render_handler(
    input = input,
    output = output,
    session = list(),
    da_data = da_data
  )

  expect_identical(output$da_results_table$data$Protein.Ids, "P1")
  expect_identical(output$da_results_table$extensions, "Buttons")
  expect_identical(output$da_results_table$format_round$columns, c("log2FC", "raw_pvalue", "fdr_qvalue"))
  expect_identical(output$da_results_table$format_round$digits, 4)
  expect_match(output$da_summary_stats, "Total genes: 3", fixed = TRUE)
  expect_match(output$da_summary_stats, "Significant (q < 0.050): 2", fixed = TRUE)
  expect_match(output$da_summary_stats, "Up-regulated: 1", fixed = TRUE)
  expect_match(output$da_summary_stats, "Down-regulated: 1", fixed = TRUE)

  download_file <- tempfile("prot-da-results-", fileext = ".zip")
  output$download_da_results$content(download_file)
  expect_match(output$download_da_results$filename(), "^DA_results_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.zip$")
  expect_match(readLines(download_file), "packaged here", fixed = TRUE)
})

test_that("proteomics DA table render handler preserves empty and error fallbacks", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$datatable_calls <- list()

  localShinyRenderMocks(captured)
  localDtRenderMocks(captured)

  input <- list(
    table_contrast = "missing contrast",
    table_significance = "significant",
    da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 1,
    table_max_rows = 25
  )

  no_results_output <- new.env(parent = emptyenv())
  no_results_data <- newProtDaResultsEnv()
  da_server_table_render_handler(
    input = input,
    output = no_results_output,
    session = list(),
    da_data = no_results_data
  )
  expect_identical(
    no_results_output$da_results_table$data$Message,
    "No results available for selected contrast"
  )
  expect_identical(no_results_output$da_summary_stats, "No results available for selected contrast")

  missing_results_output <- new.env(parent = emptyenv())
  missing_results_data <- new.env(parent = emptyenv())
  missing_results_data$da_results_list <- list(theObject = methods::new("FakeProtDaRenderObject", protein_id_column = "Protein.Ids"))
  da_server_table_render_handler(
    input = input,
    output = missing_results_output,
    session = list(),
    da_data = missing_results_data
  )
  expect_identical(
    missing_results_output$da_results_table$data$Message,
    "No DE analysis results available"
  )
  expect_identical(missing_results_output$da_summary_stats, "No DE analysis results available")

  error_results_output <- new.env(parent = emptyenv())
  error_results_data <- new.env(parent = emptyenv())
  error_results_data$da_results_list <- list(
    da_proteins_long = data.frame(log2FC = 1, stringsAsFactors = FALSE),
    theObject = NULL
  )
  da_server_table_render_handler(
    input = input,
    output = error_results_output,
    session = list(),
    da_data = error_results_data
  )
  expect_match(error_results_output$da_results_table$data$Message, "Error:", fixed = TRUE)
  expect_match(error_results_output$da_summary_stats, "Error calculating statistics:", fixed = TRUE)
})
