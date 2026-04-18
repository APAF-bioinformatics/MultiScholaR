library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "mod_metab_qc_s4_server_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_s4.R")
  ),
  symbols = c(
    "buildMetabQcS4DataSummaryUi",
    "getMetabQcS4DataSummaryState",
    "buildMetabQcS4DataSummaryRenderOutput",
    "buildMetabQcS4StateHistoryUi",
    "getMetabQcS4StateHistory",
    "buildMetabQcS4StateHistoryRenderOutput",
    "buildMetabQcS4AssayStatsDatatable",
    "getMetabQcS4AssayStatsState",
    "buildMetabQcS4AssayStatsRenderOutput",
    "buildMetabQcS4FilterPlotRenderOutput",
    "renderMetabQcS4FilterPlot",
    "buildMetabQcS4FinalizeResultsText",
    "getMetabQcS4FinalizeState",
    "validateMetabQcS4FinalizeState",
    "saveMetabQcS4CompletedState",
    "completeMetabQcS4TabStatus",
    "getMetabQcS4FinalizeHistory",
    "updateMetabQcS4TrackingPlot",
    "reportMetabQcS4FinalizeSuccess",
    "reportMetabQcS4FinalizeError",
    "runMetabQcS4FinalizeWorkflow",
    "runMetabQcS4ServerBody",
    "mod_metab_qc_s4_server"
  ),
  env = environment()
)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      design_matrix = "data.frame",
      group_id = "character",
      sample_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      design_matrix = data.frame(),
      group_id = "group",
      sample_id = "sample_id"
    )
  )
}

test_that("metabolomics QC S4 data-summary seam preserves placeholder and summary rows", {
  placeholder <- buildMetabQcS4DataSummaryUi(
    currentS4 = NULL,
    paragraphFn = function(..., style = NULL) {
      list(type = "p", children = list(...), style = style)
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    }
  )

  expect_identical(placeholder$type, "p")
  expect_identical(placeholder$children[[1L]]$name, "exclamation-triangle")
  expect_match(
    placeholder$children[[2L]],
    "No MetaboliteAssayData object available.",
    fixed = TRUE
  )

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M1", "M2"),
        Sample1 = c(1, 2, 3),
        Sample2 = c(4, 5, 6),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        feature_id = c("M3", "M4"),
        Sample1 = c(7, 8),
        Sample2 = c(9, 10),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(
      group = c("A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    group_id = "group",
    sample_id = "sample_id"
  )

  summary_ui <- buildMetabQcS4DataSummaryUi(
    currentS4 = current_s4,
    paragraphFn = function(..., style = NULL) {
      list(type = "p", children = list(...), style = style)
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    },
    tagListFn = function(...) {
      list(type = "tagList", children = list(...))
    },
    tableTagFn = function(..., class = NULL, style = NULL) {
      list(type = "table", class = class, style = style, body = list(...))
    },
    tbodyTagFn = function(...) {
      list(type = "tbody", rows = list(...))
    },
    trTagFn = function(...) {
      list(type = "tr", cells = list(...))
    },
    tdTagFn = function(...) {
      list(type = "td", children = list(...))
    },
    strongFn = function(text) {
      list(type = "strong", text = text)
    }
  )

  expect_identical(summary_ui$type, "tagList")
  expect_identical(summary_ui$children[[1L]]$type, "table")
  expect_identical(summary_ui$children[[1L]]$class, "table table-condensed")
  expect_identical(summary_ui$children[[1L]]$style, "margin-bottom: 0;")
  expect_identical(
    vapply(
      summary_ui$children[[1L]]$body[[1L]]$rows,
      function(row) row$cells[[1L]]$children[[1L]]$text,
      character(1)
    ),
    c(
      "Number of Assays:",
      "Total Metabolites:",
      "Number of Samples:",
      "Experimental Groups:",
      "Metabolite ID Column:",
      "Sample ID Column:"
    )
  )
  expect_identical(
    vapply(
      summary_ui$children[[1L]]$body[[1L]]$rows,
      function(row) as.character(row$cells[[2L]]$children[[1L]]),
      character(1)
    ),
    c("2", "4", "2", "2", "feature_id", "sample_id")
  )
})

test_that("metabolomics QC S4 data-summary fetch seam returns current state", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    getMetabQcS4DataSummaryState(
      stateManager = state_manager,
      stateGetterFn = function(manager) {
        captured$state_manager <- manager
        current_s4
      }
    )
  )

  expect_identical(captured$state_manager, state_manager)
  expect_identical(visible$value, current_s4)
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 data-summary fetch seam falls back to NULL on getter errors", {
  visible <- withVisible(
    getMetabQcS4DataSummaryState(
      stateManager = list(token = "state-manager"),
      stateGetterFn = function(manager) {
        stop("state unavailable")
      }
    )
  )

  expect_null(visible$value)
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 data-summary render-output seam threads the current-state fetch and UI helpers", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    buildMetabQcS4DataSummaryRenderOutput(
      stateManager = state_manager,
      getDataSummaryStateFn = function(stateManager, ...) {
        captured$state_manager <- stateManager
        current_s4
      },
      buildDataSummaryUiFn = function(currentS4, ...) {
        captured$current_s4 <- currentS4
        "data-summary-token"
      }
    )
  )

  expect_identical(captured$state_manager, state_manager)
  expect_identical(captured$current_s4, current_s4)
  expect_identical(visible$value, "data-summary-token")
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 assay-stats fetch seam returns current state", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    getMetabQcS4AssayStatsState(
      stateManager = state_manager,
      stateGetterFn = function(manager) {
        captured$state_manager <- manager
        current_s4
      }
    )
  )

  expect_identical(captured$state_manager, state_manager)
  expect_identical(visible$value, current_s4)
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 assay-stats fetch seam falls back to NULL on getter errors", {
  visible <- withVisible(
    getMetabQcS4AssayStatsState(
      stateManager = list(token = "state-manager"),
      stateGetterFn = function(manager) {
        stop("state unavailable")
      }
    )
  )

  expect_null(visible$value)
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 assay-stats render-output seam threads the current-state fetch and datatable helpers", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    buildMetabQcS4AssayStatsRenderOutput(
      stateManager = state_manager,
      getAssayStatsStateFn = function(stateManager, ...) {
        captured$state_manager <- stateManager
        current_s4
      },
      buildAssayStatsDatatableFn = function(currentS4, ...) {
        captured$current_s4 <- currentS4
        "assay-stats-token"
      }
    )
  )

  expect_identical(captured$state_manager, state_manager)
  expect_identical(captured$current_s4, current_s4)
  expect_identical(visible$value, "assay-stats-token")
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 state-history seam preserves placeholder and current marker", {
  placeholder <- buildMetabQcS4StateHistoryUi(
    history = character(0),
    paragraphFn = function(..., style = NULL) {
      list(type = "p", children = list(...), style = style)
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    }
  )

  expect_identical(placeholder$type, "p")
  expect_identical(placeholder$children[[1L]]$name, "info-circle")
  expect_match(
    placeholder$children[[2L]],
    "No processing history available.",
    fixed = TRUE
  )
  expect_identical(placeholder$style, "color: #666;")

  state_history_ui <- buildMetabQcS4StateHistoryUi(
    history = c("raw_import", "qc_filtering", "metab_qc_complete"),
    paragraphFn = function(..., style = NULL) {
      list(type = "p", children = list(...), style = style)
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    },
    orderedListTagFn = function(..., style = NULL) {
      list(type = "ol", children = list(...), style = style)
    },
    listItemTagFn = function(...) {
      list(type = "li", children = list(...))
    },
    spanTagFn = function(text, style = NULL) {
      list(type = "span", text = text, style = style)
    }
  )

  expect_identical(state_history_ui$type, "ol")
  expect_identical(state_history_ui$style, "list-style: none; padding-left: 0;")
  expect_length(state_history_ui$children, 3)
  expect_identical(
    vapply(
      state_history_ui$children,
      function(item) item$children[[1L]]$name,
      character(1)
    ),
    c("check", "check", "arrow-right")
  )
  expect_identical(
    vapply(
      state_history_ui$children,
      function(item) item$children[[2L]]$text,
      character(1)
    ),
    c(
      " 1. raw_import",
      " 2. qc_filtering",
      " 3. metab_qc_complete"
    )
  )
  expect_identical(state_history_ui$children[[3L]]$children[[2L]]$style, "font-weight: bold;")
  expect_identical(state_history_ui$children[[3L]]$children[[3L]]$text, " (current)")
  expect_identical(state_history_ui$children[[3L]]$children[[3L]]$style, "color: blue;")
})

test_that("metabolomics QC S4 state-history fetch seam retrieves history through the state manager getter", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  history <- c("raw_import", "qc_filtering")

  visible <- withVisible(
    getMetabQcS4StateHistory(
      stateManager = state_manager,
      historyGetterFn = function(manager) {
        captured$state_manager <- manager
        history
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$state_manager, state_manager)
  expect_identical(visible$value, history)
})

test_that("metabolomics QC S4 state-history fetch seam falls back to empty history on getter errors", {
  visible <- withVisible(
    getMetabQcS4StateHistory(
      stateManager = list(token = "state-manager"),
      historyGetterFn = function(manager) {
        stop("history unavailable")
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(visible$value, character(0))
})

test_that("metabolomics QC S4 state-history render-output seam threads the history fetch and UI helpers", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  history <- c("raw_import", "qc_filtering")

  visible <- withVisible(
    buildMetabQcS4StateHistoryRenderOutput(
      stateManager = state_manager,
      getStateHistoryFn = function(stateManager, ...) {
        captured$state_manager <- stateManager
        history
      },
      buildStateHistoryUiFn = function(history, ...) {
        captured$history <- history
        "state-history-token"
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$state_manager, state_manager)
  expect_identical(captured$history, history)
  expect_identical(visible$value, "state-history-token")
})

test_that("metabolomics QC S4 assay-stats seam preserves null exits and assay metrics", {
  expect_null(buildMetabQcS4AssayStatsDatatable(NULL))

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M1", "M2"),
        Sample1 = c(1, 0, 3),
        Sample2 = c(NA, 5, 6),
        note = c("a", "b", "c"),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        feature_id = c("M3", "M4"),
        Sample1 = c(7, 8),
        Sample2 = c(9, 0),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = c("A", "B"), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  table_widget <- buildMetabQcS4AssayStatsDatatable(
    currentS4 = current_s4,
    datatableFn = function(data, options, rownames, class) {
      list(
        data = data,
        options = options,
        rownames = rownames,
        class = class
      )
    }
  )

  expect_identical(table_widget$options, list(dom = "t", paging = FALSE, ordering = FALSE))
  expect_false(table_widget$rownames)
  expect_identical(table_widget$class, "compact stripe")
  expect_identical(table_widget$data$Assay, c("Plasma", "Serum"))
  expect_identical(table_widget$data$Metabolites, c(2, 2))
  expect_identical(table_widget$data$Samples, c(2, 2))
  expect_equal(table_widget$data$Missingness, c(33.3, 25.0))
})

test_that("metabolomics QC S4 filter-plot seam dispatches grid and ggplot payloads", {
  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  grob_object <- structure(list(label = "grid"), class = c("gtable", "grob"))
  ggplot_object <- structure(list(label = "plot"), class = "ggplot")

  grob_visible <- withVisible(
    renderMetabQcS4FilterPlot(
      filterPlot = function() grob_object,
      reqFn = function(value) {
        captured$req_values[[length(captured$req_values) + 1L]] <- value
        value
      },
      inheritsFn = function(object, class_name) inherits(object, class_name),
      gridDrawFn = function(object) {
        captured$grid_draw_object <- object
        invisible(NULL)
      },
      printFn = function(object) {
        captured$unexpected_print_object <- object
        invisible(NULL)
      }
    )
  )

  expect_false(grob_visible$visible)
  expect_identical(captured$req_values[[1L]], grob_object)
  expect_identical(captured$grid_draw_object, grob_object)
  expect_false(exists("unexpected_print_object", envir = captured, inherits = FALSE))

  ggplot_visible <- withVisible(
    renderMetabQcS4FilterPlot(
      filterPlot = function() ggplot_object,
      reqFn = function(value) {
        captured$req_values[[length(captured$req_values) + 1L]] <- value
        value
      },
      inheritsFn = function(object, class_name) inherits(object, class_name),
      gridDrawFn = function(object) {
        captured$unexpected_grid_draw_object <- object
        invisible(NULL)
      },
      printFn = function(object) {
        captured$print_object <- object
        invisible(NULL)
      }
    )
  )

  expect_false(ggplot_visible$visible)
  expect_identical(captured$req_values[[2L]], ggplot_object)
  expect_identical(captured$print_object, ggplot_object)
  expect_false(exists("unexpected_grid_draw_object", envir = captured, inherits = FALSE))
})

test_that("metabolomics QC S4 filter-plot render-output seam threads the render helper", {
  captured <- new.env(parent = emptyenv())

  filter_plot <- function(new_value) {
    if (missing(new_value)) {
      "plot-token"
    } else {
      invisible(NULL)
    }
  }

  visible <- withVisible(
    buildMetabQcS4FilterPlotRenderOutput(
      filterPlot = filter_plot,
      renderFilterPlotFn = function(filterPlot, ...) {
        captured$filter_plot <- filterPlot
        "filter-plot-render-token"
      }
    )
  )

  expect_identical(captured$filter_plot, filter_plot)
  expect_identical(visible$value, "filter-plot-render-token")
  expect_true(visible$visible)
})

test_that("metabolomics QC S4 finalize-results seam preserves retained-count and history text", {
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M1", "M2"),
        Sample1 = c(1, 2, 3),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        feature_id = c("M3", "M4"),
        Sample1 = c(4, 5),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = c("A", "B"), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  result_text <- buildMetabQcS4FinalizeResultsText(
    currentS4 = current_s4,
    history = c("raw_import", "qc_filtering", "metab_qc_complete")
  )

  expect_match(result_text, "QC Finalization Complete", fixed = TRUE)
  expect_match(result_text, "Total metabolites retained: 4", fixed = TRUE)
  expect_match(result_text, "Number of assays: 2", fixed = TRUE)
  expect_match(result_text, "Processing steps completed: 3", fixed = TRUE)
  expect_match(
    result_text,
    "Processing History:\n  1\\. raw_import\n  2\\. qc_filtering\n  3\\. metab_qc_complete"
  )
  expect_match(result_text, "State saved as: 'metab_qc_complete'", fixed = TRUE)
  expect_match(
    result_text,
    "You can now proceed to the Normalization tab.",
    fixed = TRUE
  )
})

test_that("metabolomics QC S4 finalize-state seam fetches the current state through the state manager getter", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    getMetabQcS4FinalizeState(
      stateManager = state_manager,
      stateGetterFn = function(manager) {
        captured$state_manager <- manager
        current_s4
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$state_manager, state_manager)
  expect_identical(visible$value, current_s4)
})

test_that("metabolomics QC S4 finalize-state seam bubbles getter errors", {
  expect_error(
    getMetabQcS4FinalizeState(
      stateManager = list(token = "state-manager"),
      stateGetterFn = function(manager) {
        stop("state unavailable")
      }
    ),
    "state unavailable"
  )
})

test_that("metabolomics QC S4 finalize-state validation seam requires the current state and class", {
  captured <- new.env(parent = emptyenv())
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = character(), stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    validateMetabQcS4FinalizeState(
      currentS4 = current_s4,
      reqFn = function(value) {
        captured$req_value <- value
        value
      },
      inheritsFn = function(object, class_name) {
        captured$inherits_call <- list(
          object = object,
          class_name = class_name
        )
        TRUE
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$req_value, current_s4)
  expect_identical(captured$inherits_call$object, current_s4)
  expect_identical(captured$inherits_call$class_name, "MetaboliteAssayData")
  expect_identical(visible$value, current_s4)
})

test_that("metabolomics QC S4 finalize-state validation seam rejects invalid state classes", {
  expect_error(
    validateMetabQcS4FinalizeState(
      currentS4 = list(token = "not-s4"),
      reqFn = function(value) value,
      inheritsFn = function(object, class_name) FALSE
    ),
    "Current state is not a MetaboliteAssayData object"
  )
})

test_that("metabolomics QC S4 completed-state seam persists the finalized state and returns the saved state name", {
  captured <- new.env(parent = emptyenv())
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M2"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = "A", stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )
  state_manager <- list(
    saveState = function(...) {
      captured$save_state_call <- list(...)
      invisible(NULL)
    }
  )
  config_object <- list(config = "token")

  visible <- withVisible(
    saveMetabQcS4CompletedState(
      stateManager = state_manager,
      currentS4 = current_s4,
      configObject = config_object
    )
  )

  expect_false(visible$visible)
  expect_identical(visible$value, "metab_qc_complete")
  expect_identical(captured$save_state_call$state_name, "metab_qc_complete")
  expect_identical(captured$save_state_call$s4_data_object, current_s4)
  expect_identical(captured$save_state_call$config_object, config_object)
  expect_identical(
    captured$save_state_call$description,
    "QC processing complete - ready for normalization"
  )
})

test_that("metabolomics QC S4 tab-status seam completes quality control via list replacement", {
  workflow_data <- list2env(
    list(
      tab_status = list(
        setup_import = "complete",
        quality_control = "pending",
        normalization = "locked"
      )
    ),
    parent = emptyenv()
  )

  visible <- withVisible(
    completeMetabQcS4TabStatus(workflowData = workflow_data)
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    list(
      setup_import = "complete",
      quality_control = "complete",
      normalization = "locked"
    )
  )
  expect_identical(workflow_data$tab_status, visible$value)
})

test_that("metabolomics QC S4 finalize-history seam fetches history through the state manager getter", {
  captured <- new.env(parent = emptyenv())
  state_manager <- list(token = "state-manager")
  history <- c("raw_import", "qc_filtering", "metab_qc_complete")

  visible <- withVisible(
    getMetabQcS4FinalizeHistory(
      stateManager = state_manager,
      historyGetterFn = function(manager) {
        captured$state_manager <- manager
        history
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$state_manager, state_manager)
  expect_identical(visible$value, history)
})

test_that("metabolomics QC S4 tracking-plot seam refreshes and stores the QC plot", {
  captured <- new.env(parent = emptyenv())
  captured$plot_updates <- list()
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M2"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = "A", stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    updateMetabQcS4TrackingPlot(
      currentS4 = current_s4,
      omicType = "metabolomics",
      setFilterPlotFn = function(value) {
        captured$plot_updates <- c(captured$plot_updates, list(value))
        invisible(NULL)
      },
      updateMetaboliteFilteringFn = function(theObject, step_name, omics_type, return_grid, overwrite) {
        captured$qc_update <- list(
          the_object = theObject,
          step_name = step_name,
          omics_type = omics_type,
          return_grid = return_grid,
          overwrite = overwrite
        )
        "qc-plot-token"
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$qc_update$the_object, current_s4)
  expect_identical(captured$qc_update$step_name, "4_QC_Complete")
  expect_identical(captured$qc_update$omics_type, "metabolomics")
  expect_true(isTRUE(captured$qc_update$return_grid))
  expect_true(isTRUE(captured$qc_update$overwrite))
  expect_identical(captured$plot_updates, list("qc-plot-token"))
  expect_identical(visible$value, "qc-plot-token")
})

test_that("metabolomics QC S4 tracking-plot seam falls back to a null plot on refresh errors", {
  captured <- new.env(parent = emptyenv())
  captured$plot_updates <- list()
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M2"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = "A", stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )

  visible <- withVisible(
    updateMetabQcS4TrackingPlot(
      currentS4 = current_s4,
      omicType = "metabolomics",
      setFilterPlotFn = function(value) {
        captured$plot_updates <- c(captured$plot_updates, list(value))
        invisible(NULL)
      },
      updateMetaboliteFilteringFn = function(...) {
        stop("plot refresh failed")
      }
    )
  )

  expect_false(visible$visible)
  expect_length(captured$plot_updates, 1L)
  expect_null(captured$plot_updates[[1L]])
  expect_null(visible$value)
})

test_that("metabolomics QC S4 finalize-success seam renders results and reports completion", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M2"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = "A", stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )
  history <- c("raw_import", "qc_filtering", "metab_qc_complete")

  visible <- withVisible(
    reportMetabQcS4FinalizeSuccess(
      currentS4 = current_s4,
      history = history,
      output = output,
      buildResultsTextFn = function(currentS4, history, ...) {
        captured$current_s4 <- currentS4
        captured$history <- history
        "finalize-summary-token"
      },
      renderTextFn = function(text) {
        captured$render_text_value <- text
        paste0("rendered::", text)
      },
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, ...) {
        captured$notification <- c(list(message = message), list(...))
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$current_s4, current_s4)
  expect_identical(captured$history, history)
  expect_identical(captured$render_text_value, "finalize-summary-token")
  expect_identical(output$finalize_results, "rendered::finalize-summary-token")
  expect_identical(captured$log_message, "Metabolomics QC finalized successfully")
  expect_identical(
    captured$notification,
    list(
      message = "QC complete! Proceed to Normalization.",
      type = "message",
      duration = 5
    )
  )
  expect_identical(visible$value, "finalize-summary-token")
})

test_that("metabolomics QC S4 finalize-error seam formats, logs, and notifies failures", {
  captured <- new.env(parent = emptyenv())
  error_object <- simpleError("state unavailable")

  visible <- withVisible(
    reportMetabQcS4FinalizeError(
      error = error_object,
      logErrorFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, ...) {
        captured$notification <- c(list(message = message), list(...))
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$log_message, "Error finalizing QC: state unavailable")
  expect_identical(
    captured$notification,
    list(
      message = "Error finalizing QC: state unavailable",
      type = "error"
    )
  )
  expect_identical(visible$value, "Error finalizing QC: state unavailable")
})

test_that("metabolomics QC S4 server delegates the state-history render path through the render-output seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "buildMetabQcS4StateHistoryRenderOutput"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  assign(
    binding_name,
    function(stateManager, ...) {
      captured$render_output_call <- list(
        state_manager = stateManager
      )
      "state-history-token"
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      stored <- value

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    renderUI = function(expr) {
      rendered_expr <- substitute(expr)
      rendered_text <- paste(deparse(rendered_expr), collapse = "\n")

      if (grepl("buildStateHistoryRenderOutputFn", rendered_text, fixed = TRUE)) {
        eval(rendered_expr, parent.frame())
      } else {
        "other-ui-token"
      }
    },
    observeEvent = function(...) {
      invisible(NULL)
    },
    renderPlot = function(expr) {
      substitute(expr)
    },
    req = function(...) {
      list(...)[[1L]]
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      substitute(expr)
    },
    .package = "DT"
  )

  workflow_data <- list(
    state_manager = list(
      getHistory = function(...) stop("getHistory should not be called directly")
    )
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_identical(captured$render_output_call$state_manager, workflow_data$state_manager)
  expect_identical(output$state_history, "state-history-token")
})

test_that("metabolomics QC S4 server delegates the data-summary render path through the render-output seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "buildMetabQcS4DataSummaryRenderOutput"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  assign(
    binding_name,
    function(stateManager, ...) {
      captured$render_output_call <- list(
        state_manager = stateManager
      )
      "data-summary-token"
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      stored <- value

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    renderUI = function(expr) {
      rendered_expr <- substitute(expr)
      rendered_text <- paste(deparse(rendered_expr), collapse = "\n")

      if (grepl("buildDataSummaryRenderOutputFn", rendered_text, fixed = TRUE)) {
        eval(rendered_expr, parent.frame())
      } else {
        "other-ui-token"
      }
    },
    observeEvent = function(...) {
      invisible(NULL)
    },
    renderPlot = function(expr) {
      substitute(expr)
    },
    req = function(...) {
      list(...)[[1L]]
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      substitute(expr)
    },
    .package = "DT"
  )

  workflow_data <- list(
    state_manager = list(
      getState = function(...) stop("getState should not be called directly")
    )
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_identical(captured$render_output_call$state_manager, workflow_data$state_manager)
  expect_identical(output$data_summary, "data-summary-token")
})

test_that("metabolomics QC S4 server delegates the assay-stats render path through the render-output seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "buildMetabQcS4AssayStatsRenderOutput"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  assign(
    binding_name,
    function(stateManager, ...) {
      captured$render_output_call <- list(
        state_manager = stateManager
      )
      "assay-stats-token"
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      stored <- value

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    renderUI = function(expr) {
      "other-ui-token"
    },
    observeEvent = function(...) {
      invisible(NULL)
    },
    renderPlot = function(expr) {
      substitute(expr)
    },
    req = function(...) {
      list(...)[[1L]]
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      rendered_expr <- substitute(expr)
      rendered_text <- paste(deparse(rendered_expr), collapse = "\n")

      if (grepl("buildAssayStatsRenderOutputFn", rendered_text, fixed = TRUE)) {
        eval(rendered_expr, parent.frame())
      } else {
        "other-dt-token"
      }
    },
    .package = "DT"
  )

  workflow_data <- list(
    state_manager = list(
      getState = function(...) stop("getState should not be called directly")
    )
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_identical(captured$render_output_call$state_manager, workflow_data$state_manager)
  expect_identical(output$assay_stats_table, "assay-stats-token")
})

test_that("metabolomics QC S4 server delegates the filter-plot render path through the render-output seam", {
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "buildMetabQcS4FilterPlotRenderOutput"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  assign(
    binding_name,
    function(filterPlot, ...) {
      captured$render_output_call <- list(
        filter_plot = filterPlot
      )
      "filter-plot-render-token"
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      stored <- value

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    renderUI = function(expr) {
      substitute(expr)
    },
    observeEvent = function(...) {
      invisible(NULL)
    },
    renderPlot = function(expr) {
      rendered_expr <- substitute(expr)
      rendered_text <- paste(deparse(rendered_expr), collapse = "\n")

      if (grepl("buildFilterPlotRenderOutputFn", rendered_text, fixed = TRUE)) {
        eval(rendered_expr, parent.frame())
      } else {
        "other-plot-token"
      }
    },
    req = function(...) {
      list(...)[[1L]]
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      substitute(expr)
    },
    .package = "DT"
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = list(
      state_manager = list(
        getState = function() NULL
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_type(captured$render_output_call$filter_plot, "closure")
  expect_identical(captured$output$filter_plot, "filter-plot-render-token")
})

test_that("metabolomics QC S4 finalize workflow seam orchestrates the downstream helper chain", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("M1", "M2"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature_id",
    design_matrix = data.frame(group = "A", stringsAsFactors = FALSE),
    group_id = "group",
    sample_id = "sample_id"
  )
  history <- c("raw_import", "qc_filtering", "metab_qc_complete")

  workflow_data <- list2env(
    list(
      state_manager = list(
        token = "state-manager"
      ),
      config_list = list(config = "token"),
      tab_status = list(
        setup_import = "complete",
        quality_control = "pending",
        normalization = "locked"
      )
    ),
    parent = emptyenv()
  )

  visible <- withVisible(
    runMetabQcS4FinalizeWorkflow(
      workflowData = workflow_data,
      omicType = "metabolomics",
      filterPlot = function(value) {
        captured$steps <- c(captured$steps, "plot-set")
        captured$filter_plot_update <- value
        invisible(NULL)
      },
      output = output,
      getFinalizeStateFn = function(stateManager, ...) {
        captured$steps <- c(captured$steps, "state")
        captured$state_manager <- stateManager
        current_s4
      },
      validateFinalizeStateFn = function(currentS4, ...) {
        captured$steps <- c(captured$steps, "validate")
        captured$validated_s4 <- currentS4
        currentS4
      },
      saveCompletedStateFn = function(stateManager, currentS4, configObject, ...) {
        captured$steps <- c(captured$steps, "save")
        captured$save_call <- list(
          state_manager = stateManager,
          current_s4 = currentS4,
          config_object = configObject
        )
        invisible("metab_qc_complete")
      },
      updateTrackingPlotFn = function(currentS4, omicType, setFilterPlotFn, ...) {
        captured$steps <- c(captured$steps, "plot")
        captured$plot_call <- list(
          current_s4 = currentS4,
          omic_type = omicType,
          set_filter_plot_fn = setFilterPlotFn
        )
        setFilterPlotFn("qc-plot-token")
        invisible("qc-plot-token")
      },
      completeTabStatusFn = function(workflowData, ...) {
        captured$steps <- c(captured$steps, "status")
        updated_status <- workflowData$tab_status
        updated_status$quality_control <- "complete"
        workflowData$tab_status <- updated_status
        invisible(updated_status)
      },
      getFinalizeHistoryFn = function(stateManager, ...) {
        captured$steps <- c(captured$steps, "history")
        captured$history_state_manager <- stateManager
        history
      },
      reportFinalizeSuccessFn = function(currentS4, history, output, ...) {
        captured$steps <- c(captured$steps, "success")
        captured$success_call <- list(
          current_s4 = currentS4,
          history = history,
          output = output
        )
        output$finalize_results <- "finalize-success-token"
        invisible("finalize-success-token")
      },
      reportFinalizeErrorFn = function(error, ...) {
        captured$steps <- c(captured$steps, "error")
        captured$unexpected_error <- error
        invisible("handled-error-token")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(visible$value, "finalize-success-token")
  expect_identical(
    captured$steps,
    c("state", "validate", "save", "plot", "plot-set", "status", "history", "success")
  )
  expect_identical(captured$state_manager, workflow_data$state_manager)
  expect_identical(captured$validated_s4, current_s4)
  expect_identical(captured$save_call$state_manager, workflow_data$state_manager)
  expect_identical(captured$save_call$current_s4, current_s4)
  expect_identical(captured$save_call$config_object, list(config = "token"))
  expect_identical(captured$plot_call$current_s4, current_s4)
  expect_identical(captured$plot_call$omic_type, "metabolomics")
  expect_type(captured$plot_call$set_filter_plot_fn, "closure")
  expect_identical(captured$filter_plot_update, "qc-plot-token")
  expect_identical(captured$history_state_manager, workflow_data$state_manager)
  expect_identical(captured$success_call$current_s4, current_s4)
  expect_identical(captured$success_call$history, history)
  expect_identical(captured$success_call$output, output)
  expect_identical(
    workflow_data$tab_status,
    list(
      setup_import = "complete",
      quality_control = "complete",
      normalization = "locked"
    )
  )
  expect_identical(output$finalize_results, "finalize-success-token")
  expect_false(exists("unexpected_error", envir = captured, inherits = FALSE))
})

test_that("metabolomics QC S4 finalize workflow seam routes failures through the error reporter", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())

  workflow_data <- list2env(
    list(
      state_manager = list(
        token = "state-manager"
      ),
      config_list = list(config = "token"),
      tab_status = list(
        setup_import = "complete",
        quality_control = "pending",
        normalization = "locked"
      )
    ),
    parent = emptyenv()
  )

  visible <- withVisible(
    runMetabQcS4FinalizeWorkflow(
      workflowData = workflow_data,
      omicType = "metabolomics",
      filterPlot = function(value) {
        captured$steps <- c(captured$steps, "plot-set")
        invisible(NULL)
      },
      output = output,
      getFinalizeStateFn = function(stateManager, ...) {
        captured$steps <- c(captured$steps, "state")
        captured$state_manager <- stateManager
        stop("state unavailable")
      },
      validateFinalizeStateFn = function(...) {
        captured$steps <- c(captured$steps, "validate")
        invisible(NULL)
      },
      saveCompletedStateFn = function(...) {
        captured$steps <- c(captured$steps, "save")
        invisible(NULL)
      },
      updateTrackingPlotFn = function(...) {
        captured$steps <- c(captured$steps, "plot")
        invisible(NULL)
      },
      completeTabStatusFn = function(...) {
        captured$steps <- c(captured$steps, "status")
        invisible(NULL)
      },
      getFinalizeHistoryFn = function(...) {
        captured$steps <- c(captured$steps, "history")
        character(0)
      },
      reportFinalizeSuccessFn = function(...) {
        captured$steps <- c(captured$steps, "success")
        invisible("unexpected-success-token")
      },
      reportFinalizeErrorFn = function(error, ...) {
        captured$steps <- c(captured$steps, "error")
        captured$error <- error
        invisible("handled-error-token")
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(visible$value, "handled-error-token")
  expect_identical(captured$steps, c("state", "error"))
  expect_identical(captured$state_manager, workflow_data$state_manager)
  expect_identical(captured$error$message, "state unavailable")
  expect_identical(
    workflow_data$tab_status,
    list(
      setup_import = "complete",
      quality_control = "pending",
      normalization = "locked"
    )
  )
  expect_false(exists("finalize_results", envir = output, inherits = FALSE))
})

test_that("metabolomics QC S4 server-body seam preserves render handoffs and finalize wiring", {
  captured <- new.env(parent = emptyenv())
  captured$req_calls <- list()
  output <- new.env(parent = emptyenv())

  workflow_data <- list2env(
    list(
      state_manager = list(token = "state-manager"),
      config_list = list(config = "token"),
      tab_status = list(
        setup_import = "complete",
        quality_control = "pending",
        normalization = "locked"
      )
    ),
    parent = emptyenv()
  )

  visible <- withVisible(
    runMetabQcS4ServerBody(
      input = list(finalize_qc = TRUE),
      output = output,
      session = list(ns = function(value) paste("s4", value, sep = "-")),
      workflowData = workflow_data,
      omicType = "metabolomics",
      experimentLabel = "Experiment A",
      reactiveValFn = function(value = NULL) {
        stored <- value

        function(new_value) {
          if (missing(new_value)) {
            stored
          } else {
            stored <<- new_value
            captured$filter_plot_update <- new_value
            invisible(NULL)
          }
        }
      },
      renderUiFn = function(expr) {
        list(expr = substitute(expr), env = parent.frame())
      },
      renderDtFn = function(expr) {
        list(expr = substitute(expr), env = parent.frame())
      },
      observeEventFn = function(eventExpr, handlerExpr, ...) {
        if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
          eval(substitute(handlerExpr), parent.frame())
        }
        invisible(NULL)
      },
      renderPlotFn = function(expr) {
        list(expr = substitute(expr), env = parent.frame())
      },
      reqFn = function(...) {
        args <- list(...)
        captured$req_calls <- c(captured$req_calls, args)

        if (length(args) == 0L) {
          return(invisible(NULL))
        }

        args[[1L]]
      },
      buildStateHistoryRenderOutputFn = function(stateManager, ...) {
        captured$state_history_call <- list(state_manager = stateManager)
        "state-history-ui-token"
      },
      buildDataSummaryRenderOutputFn = function(stateManager, ...) {
        captured$data_summary_call <- list(state_manager = stateManager)
        "data-summary-ui-token"
      },
      buildAssayStatsRenderOutputFn = function(stateManager, ...) {
        captured$assay_stats_call <- list(state_manager = stateManager)
        "assay-stats-token"
      },
      runFinalizeWorkflowFn = function(workflowData, omicType, filterPlot, output, ...) {
        captured$finalize_call <- list(
          workflow_data = workflowData,
          omic_type = omicType,
          filter_plot = filterPlot,
          output = output
        )
        filterPlot("qc-plot-token")
        output$finalize_results <- "delegated-finalize-token"
        invisible("delegated-finalize-token")
      },
      buildFilterPlotRenderOutputFn = function(filterPlot, ...) {
        captured$filter_plot_call <- list(filter_plot = filterPlot)
        "filter-plot-token"
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)

  expect_identical(captured$finalize_call$workflow_data, workflow_data)
  expect_identical(captured$finalize_call$omic_type, "metabolomics")
  expect_type(captured$finalize_call$filter_plot, "closure")
  expect_identical(captured$finalize_call$output, output)
  expect_identical(captured$filter_plot_update, "qc-plot-token")
  expect_identical(output$finalize_results, "delegated-finalize-token")

  expect_identical(
    eval(output$state_history$expr, output$state_history$env),
    "state-history-ui-token"
  )
  expect_identical(captured$state_history_call$state_manager, workflow_data$state_manager)

  expect_identical(
    eval(output$data_summary$expr, output$data_summary$env),
    "data-summary-ui-token"
  )
  expect_identical(captured$data_summary_call$state_manager, workflow_data$state_manager)

  expect_identical(
    eval(output$assay_stats_table$expr, output$assay_stats_table$env),
    "assay-stats-token"
  )
  expect_identical(captured$assay_stats_call$state_manager, workflow_data$state_manager)

  expect_identical(
    eval(output$filter_plot$expr, output$filter_plot$env),
    "filter-plot-token"
  )
  expect_type(captured$filter_plot_call$filter_plot, "closure")
  expect_identical(captured$filter_plot_call$filter_plot, captured$finalize_call$filter_plot)

  expect_length(captured$req_calls, 4L)
  expect_true(all(vapply(captured$req_calls, function(arg) {
    identical(arg, workflow_data$state_manager)
  }, logical(1))))
})

test_that("metabolomics QC S4 module server delegates through the server body seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "runMetabQcS4ServerBody"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  expected_workflow_data <- list(
    state_manager = list(
      getState = function() "current-s4"
    )
  )

  assign(
    binding_name,
    function(input, output, session, workflowData, omicType, experimentLabel, ...) {
      captured$server_body_call <- list(
        input = input,
        output = output,
        session = session,
        workflow_data = workflowData,
        omic_type = omicType,
        experiment_label = experimentLabel
      )
      invisible(NULL)
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(finalize_qc = FALSE),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    .package = "shiny"
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = expected_workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_identical(captured$server_body_call$input$finalize_qc, FALSE)
  expect_identical(captured$server_body_call$output, output)
  expect_identical(captured$server_body_call$session$ns("filter_plot"), "s4-filter_plot")
  expect_identical(captured$server_body_call$workflow_data, expected_workflow_data)
  expect_identical(captured$server_body_call$omic_type, "metabolomics")
  expect_identical(captured$server_body_call$experiment_label, "Experiment A")
})

test_that("metabolomics QC S4 server delegates the finalize observer through the workflow seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_s4_server)
  binding_name <- "runMetabQcS4FinalizeWorkflow"
  had_binding <- exists(binding_name, envir = server_env, inherits = FALSE)

  if (had_binding) {
    original_binding <- get(binding_name, envir = server_env, inherits = FALSE)
  }

  on.exit({
    if (had_binding) {
      assign(binding_name, original_binding, envir = server_env)
    } else if (exists(binding_name, envir = server_env, inherits = FALSE)) {
      rm(list = binding_name, envir = server_env)
    }
  }, add = TRUE)

  assign(
    binding_name,
    function(workflowData, omicType, filterPlot, output, ...) {
      captured$workflow_call <- list(
        workflow_data = workflowData,
        omic_type = omicType,
        filter_plot = filterPlot,
        output = output
      )
      filterPlot("qc-plot-token")
      output$finalize_results <- "delegated-finalize-token"
      invisible("delegated-finalize-token")
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(finalize_qc = TRUE),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      stored <- value

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          captured$filter_plot_update <- new_value
          invisible(NULL)
        }
      }
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderUI = function(expr) {
      substitute(expr)
    },
    renderPlot = function(expr) {
      substitute(expr)
    },
    req = function(...) {
      args <- list(...)
      captured$req_calls <- c(captured$req_calls, args)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    renderDT = function(expr) {
      substitute(expr)
    },
    .package = "DT"
  )

  workflow_data <- list2env(
    list(
      state_manager = list(
        getState = function(...) stop("getState should not be called directly"),
        saveState = function(...) stop("saveState should not be called directly"),
        getHistory = function(...) stop("getHistory should not be called directly")
      ),
      config_list = list(config = "token"),
      tab_status = list(
        setup_import = "complete",
        quality_control = "pending",
        normalization = "locked"
      )
    ),
    parent = emptyenv()
  )

  mod_metab_qc_s4_server(
    id = "s4",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "s4")
  expect_identical(captured$workflow_call$workflow_data, workflow_data)
  expect_identical(captured$workflow_call$omic_type, "metabolomics")
  expect_type(captured$workflow_call$filter_plot, "closure")
  expect_identical(captured$workflow_call$output, output)
  expect_identical(captured$filter_plot_update, "qc-plot-token")
  expect_length(captured$req_calls, 1L)
  expect_identical(captured$req_calls[[1L]], workflow_data$state_manager)
  expect_identical(output$finalize_results, "delegated-finalize-token")
})
