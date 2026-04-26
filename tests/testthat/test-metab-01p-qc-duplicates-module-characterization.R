library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingMetabQcDuplicatesTargetFiles <- function() {
  required_paths <- c(
    "R/mod_metab_qc_duplicates_ui_helpers.R",
    "R/mod_metab_qc_duplicates_render_helpers.R",
    "R/mod_metab_qc_duplicates_server_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only metab QC duplicate helper file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingMetabQcDuplicatesTargetFiles()

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
    file.path(repo_root, "R", "mod_metab_qc_duplicates_ui_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_render_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_server_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_ui.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_server.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates.R")
  ),
  symbols = c(
    "reportMetabDuplicateDetection",
    "reportMetabDuplicateDetectionError",
    "runMetabDuplicateDetectionObserverShell",
    "runMetabDuplicateResolutionObserver",
    "runMetabDuplicateResolutionObserverShell",
    "runMetabDuplicateRevertObserverShell",
    "prepareMetabDuplicateResolutionState",
    "applyMetabDuplicateResolutionState",
    "runMetabDuplicateResolutionWorkflow",
    "resolveMetabDuplicateAssayData",
    "buildMetabDuplicateResolutionSummary",
    "buildMetabDuplicateSummaryUi",
    "buildMetabDuplicateTablesUi",
    "registerMetabDuplicateTableRenderers",
    "detectMetabDuplicateFeatures",
    "revertMetabDuplicateResolution",
    "renderMetabDuplicateFilterPlot",
    "mod_metab_qc_duplicates_server"
  ),
  env = environment()
)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id"
    )
  )
}

test_that("metabolomics duplicate assay seam preserves per-assay resolution stats and warnings", {
  capture <- new.env(parent = emptyenv())
  capture$calls <- list()
  capture$warnings <- character()

  assay_list <- list(
    Plasma = data.frame(
      metabolite_id = c("m1", "m1", "m2"),
      annotation_id = c("a1", "a1", "a2"),
      Sample1 = c(1, 5, 3),
      Sample2 = c(2, 4, 7),
      stringsAsFactors = FALSE
    ),
    MetadataOnly = data.frame(
      metabolite_id = c("m3", "m4"),
      annotation_id = c("a3", "a4"),
      note = c("x", "y"),
      stringsAsFactors = FALSE
    )
  )

  visible <- withVisible(
    resolveMetabDuplicateAssayData(
      assayList = assay_list,
      metaboliteIdCol = "metabolite_id",
      resolveDuplicateFeaturesByIntensityFn = function(assay_tibble, id_col, sample_cols) {
        capture$calls[[length(capture$calls) + 1L]] <- list(
          id_col = id_col,
          sample_cols = sample_cols,
          input_rows = nrow(assay_tibble)
        )
        assay_tibble[!duplicated(assay_tibble[[id_col]]), , drop = FALSE]
      },
      logWarnFn = function(message) {
        capture$warnings <- c(capture$warnings, message)
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_length(capture$calls, 1L)
  expect_identical(capture$calls[[1L]]$id_col, "metabolite_id")
  expect_identical(capture$calls[[1L]]$sample_cols, c("Sample1", "Sample2"))
  expect_identical(
    vapply(visible$value$statsList, `[[`, numeric(1), "removed"),
    c(Plasma = 1, MetadataOnly = 0)
  )
  expect_identical(names(visible$value$resolvedAssayList), c("Plasma", "MetadataOnly"))
  expect_identical(nrow(visible$value$resolvedAssayList$Plasma), 2L)
  expect_identical(visible$value$resolvedAssayList$MetadataOnly, assay_list$MetadataOnly)
  expect_identical(
    capture$warnings,
    "No numeric columns found in assay: MetadataOnly"
  )
})

test_that("metabolomics duplicate resolution summary seam assembles totals and result text", {
  summary <- buildMetabDuplicateResolutionSummary(
    statsList = list(
      Plasma = list(original = 4, resolved = 2, removed = 2),
      Serum = list(original = 3, resolved = 3, removed = 0)
    ),
    stateName = "custom_metab_state"
  )

  expect_identical(summary$totalRemoved, 2)
  expect_match(summary$resultText, "Duplicate Resolution Complete", fixed = TRUE)
  expect_match(
    summary$resultText,
    "  Plasma: 4 -> 2 rows \\(removed 2 duplicates\\)"
  )
  expect_match(
    summary$resultText,
    "  Serum: 3 -> 3 rows \\(removed 0 duplicates\\)"
  )
  expect_match(summary$resultText, "Total duplicate rows removed: 2", fixed = TRUE)
  expect_match(summary$resultText, "State saved as: 'custom_metab_state'", fixed = TRUE)
})

test_that("metabolomics duplicate resolution preflight seam validates state and dispatches assay resolution", {
  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  initial_assay_data <- list(
    Plasma = data.frame(
      metabolite_id = c("m1", "m1"),
      Sample1 = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = initial_assay_data,
    metabolite_id_column = "metabolite_id"
  )
  resolved_assay_data <- list(
    Plasma = data.frame(
      metabolite_id = "m1",
      Sample1 = 2,
      stringsAsFactors = FALSE
    )
  )
  expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )

  visible <- withVisible(
    prepareMetabDuplicateResolutionState(
      stateManager = list(
        getState = function() {
          captured$get_state_called <- TRUE
          current_s4
        }
      ),
      resolveDuplicateAssayDataFn = function(assayList, metaboliteIdCol) {
        captured$resolve_call <- list(
          assay_list = assayList,
          metabolite_id_col = metaboliteIdCol
        )
        list(
          resolvedAssayList = resolved_assay_data,
          statsList = expected_stats
        )
      },
      reqFn = function(value) {
        captured$req_values[[length(captured$req_values) + 1L]] <- value
        value
      },
      inheritsFn = function(object, class_name) {
        captured$inherits_call <- list(
          object = object,
          class_name = class_name
        )
        inherits(object, class_name)
      }
    )
  )

  expect_true(visible$visible)
  expect_true(isTRUE(captured$get_state_called))
  expect_identical(length(captured$req_values), 1L)
  expect_identical(captured$req_values[[1L]], current_s4)
  expect_identical(captured$inherits_call$class_name, "MetaboliteAssayData")
  expect_identical(captured$resolve_call$assay_list, initial_assay_data)
  expect_identical(captured$resolve_call$metabolite_id_col, "metabolite_id")
  expect_identical(visible$value$currentS4@metabolite_data, resolved_assay_data)
  expect_identical(visible$value$statsList, expected_stats)
})

test_that("metabolomics duplicate resolution apply seam updates reactive state, saves workflow state, and refreshes QC plot", {
  captured <- new.env(parent = emptyenv())
  captured$resolution_stats_updates <- list()
  captured$filter_plot_updates <- list()
  captured$warnings <- character()
  resolved_assay_data <- list(
    Plasma = data.frame(
      metabolite_id = "m1",
      Sample1 = 2,
      stringsAsFactors = FALSE
    )
  )
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = resolved_assay_data,
    metabolite_id_column = "metabolite_id"
  )
  expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )
  workflow_data <- list(
    state_manager = list(
      saveState = function(state_name, s4_data_object, config_object, description) {
        captured$save_state <- list(
          state_name = state_name,
          metabolite_data = s4_data_object@metabolite_data,
          config_object = config_object,
          description = description
        )
        invisible(NULL)
      }
    ),
    config_list = list(config = "value")
  )

  visible <- withVisible(
    applyMetabDuplicateResolutionState(
      currentS4 = current_s4,
      statsList = expected_stats,
      workflowData = workflow_data,
      omicType = "metabolomics",
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      updateMetaboliteFilteringFn = function(theObject, step_name, omics_type, return_grid, overwrite) {
        captured$qc_update <- list(
          metabolite_data = theObject@metabolite_data,
          step_name = step_name,
          omics_type = omics_type,
          return_grid = return_grid,
          overwrite = overwrite
        )
        "plot-token"
      },
      logWarnFn = function(message) {
        captured$warnings <- c(captured$warnings, message)
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$resolution_stats_updates, list(expected_stats))
  expect_identical(captured$save_state$state_name, "metab_duplicates_resolved")
  expect_identical(captured$save_state$metabolite_data, resolved_assay_data)
  expect_identical(captured$save_state$config_object, list(config = "value"))
  expect_identical(
    captured$save_state$description,
    "Resolved duplicate metabolite features by keeping highest intensity"
  )
  expect_identical(captured$qc_update$metabolite_data, resolved_assay_data)
  expect_identical(captured$qc_update$step_name, "3_Duplicates_Resolved")
  expect_identical(captured$qc_update$omics_type, "metabolomics")
  expect_true(isTRUE(captured$qc_update$return_grid))
  expect_true(isTRUE(captured$qc_update$overwrite))
  expect_identical(captured$filter_plot_updates, list("plot-token"))
  expect_identical(captured$warnings, character())
  expect_identical(visible$value$stateName, "metab_duplicates_resolved")
  expect_identical(visible$value$qcPlot, "plot-token")
})

test_that("metabolomics duplicate resolution apply seam falls back to a null QC plot when refresh fails", {
  captured <- new.env(parent = emptyenv())
  captured$resolution_stats_updates <- list()
  captured$filter_plot_updates <- list()
  captured$warnings <- character()
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m1",
        Sample1 = 2,
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "metabolite_id"
  )
  expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )
  workflow_data <- list(
    state_manager = list(
      saveState = function(state_name, s4_data_object, config_object, description) {
        captured$save_state <- list(
          state_name = state_name,
          metabolite_data = s4_data_object@metabolite_data,
          config_object = config_object,
          description = description
        )
        invisible(NULL)
      }
    ),
    config_list = list(config = "value")
  )

  visible <- withVisible(
    applyMetabDuplicateResolutionState(
      currentS4 = current_s4,
      statsList = expected_stats,
      workflowData = workflow_data,
      omicType = "metabolomics",
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      updateMetaboliteFilteringFn = function(...) {
        stop("plot refresh failed")
      },
      logWarnFn = function(message) {
        captured$warnings <- c(captured$warnings, message)
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$resolution_stats_updates, list(expected_stats))
  expect_identical(captured$save_state$state_name, "metab_duplicates_resolved")
  expect_identical(
    captured$warnings,
    "Could not generate QC plot: plot refresh failed"
  )
  expect_length(captured$filter_plot_updates, 1L)
  expect_null(captured$filter_plot_updates[[1L]])
  expect_identical(visible$value$stateName, "metab_duplicates_resolved")
  expect_null(visible$value$qcPlot)
})

test_that("metabolomics duplicate resolution workflow seam composes preflight, apply, and summary helpers", {
  captured <- new.env(parent = emptyenv())
  captured$resolution_stats_updates <- list()
  captured$filter_plot_updates <- list()

  resolved_assay_data <- list(
    Plasma = data.frame(
      metabolite_id = "m1",
      Sample1 = 2,
      stringsAsFactors = FALSE
    )
  )
  resolved_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = resolved_assay_data,
    metabolite_id_column = "metabolite_id"
  )
  expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )
  workflow_data <- list(
    state_manager = structure(list(id = "state-manager"), class = "mock_state_manager"),
    config_list = list(config = "value")
  )

  visible <- withVisible(
    runMetabDuplicateResolutionWorkflow(
      workflowData = workflow_data,
      omicType = "metabolomics",
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      prepareResolutionStateFn = function(stateManager) {
        captured$preflight_state_manager <- stateManager
        list(
          currentS4 = resolved_s4,
          statsList = expected_stats
        )
      },
      applyResolutionStateFn = function(currentS4, statsList, workflowData, omicType, setResolutionStatsFn, setFilterPlotFn) {
        captured$apply_call <- list(
          metabolite_data = currentS4@metabolite_data,
          stats_list = statsList,
          workflow_data = workflowData,
          omic_type = omicType
        )
        setResolutionStatsFn(statsList)
        setFilterPlotFn("plot-token")
        list(
          stateName = "metab_duplicates_resolved",
          qcPlot = "plot-token"
        )
      },
      buildResolutionSummaryFn = function(statsList, stateName) {
        captured$summary_call <- list(
          stats_list = statsList,
          state_name = stateName
        )
        list(
          totalRemoved = 1,
          resultText = "summary-token"
        )
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(captured$preflight_state_manager, workflow_data$state_manager)
  expect_identical(captured$apply_call$metabolite_data, resolved_assay_data)
  expect_identical(captured$apply_call$stats_list, expected_stats)
  expect_identical(captured$apply_call$workflow_data, workflow_data)
  expect_identical(captured$apply_call$omic_type, "metabolomics")
  expect_identical(captured$summary_call$stats_list, expected_stats)
  expect_identical(captured$summary_call$state_name, "metab_duplicates_resolved")
  expect_identical(captured$resolution_stats_updates, list(expected_stats))
  expect_identical(captured$filter_plot_updates, list("plot-token"))
  expect_identical(visible$value$resultText, "summary-token")
  expect_identical(visible$value$totalRemoved, 1)
})

test_that("metabolomics duplicate summary UI seam assembles placeholder and assay rows", {
  placeholder <- buildMetabDuplicateSummaryUi(
    dupList = NULL,
    paragraphFn = function(..., style = NULL) {
      list(type = "p", children = list(...), style = style)
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    }
  )

  expect_identical(placeholder$type, "p")
  expect_identical(placeholder$children[[1]]$name, "info-circle")
  expect_identical(placeholder$style, "color: #666;")
  expect_match(placeholder$children[[2]], "Detect Duplicates", fixed = TRUE)

  summary_ui <- buildMetabDuplicateSummaryUi(
    dupList = list(
      Plasma = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE),
      Serum = NULL
    ),
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    },
    listTagFn = function(items, style = NULL) {
      list(type = "ul", items = items, style = style)
    },
    itemTagFn = function(...) {
      list(type = "li", children = list(...))
    }
  )

  expect_identical(summary_ui$type, "ul")
  expect_identical(summary_ui$style, "list-style: none; padding-left: 0;")
  expect_length(summary_ui$items, 2L)
  expect_identical(summary_ui$items[[1]]$children[[1]]$name, "exclamation-triangle")
  expect_identical(summary_ui$items[[2]]$children[[1]]$name, "check-circle")
  expect_identical(
    summary_ui$items[[1]]$children[[2]],
    " Plasma: 2 duplicate IDs"
  )
  expect_identical(
    summary_ui$items[[2]]$children[[2]],
    " Serum: 0 duplicate IDs"
  )
})

test_that("metabolomics duplicate tables UI seam filters assay tabs and empty states", {
  empty_state <- buildMetabDuplicateTablesUi(
    dupList = list(
      Plasma = NULL,
      Serum = data.frame(feature = character(), stringsAsFactors = FALSE)
    ),
    wellPanelFn = function(...) {
      list(type = "well", children = list(...))
    },
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    },
    headerFn = function(text) {
      list(type = "h5", text = text)
    },
    paragraphFn = function(text) {
      list(type = "p", text = text)
    }
  )

  expect_identical(empty_state$type, "well")
  expect_identical(empty_state$children[[1]]$name, "check-circle")
  expect_identical(
    empty_state$children[[2]]$text,
    "No duplicates found in any assay!"
  )
  expect_identical(
    empty_state$children[[3]]$text,
    "All metabolite IDs are unique. No resolution needed."
  )

  tabbed_ui <- buildMetabDuplicateTablesUi(
    dupList = list(
      "Plasma/Serum" = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE),
      Control = NULL,
      Urine = data.frame(feature = "m2", stringsAsFactors = FALSE)
    ),
    nsFn = function(id) paste0("duplicates-", id),
    tabPanelFn = function(title, ...) {
      list(type = "tab", title = title, children = list(...))
    },
    breakFn = function() {
      list(type = "br")
    },
    tableOutputFn = function(outputId) {
      list(type = "dt", outputId = outputId)
    },
    tabsetPanelFn = function(id, ...) {
      list(type = "tabset", id = id, children = list(...))
    }
  )

  expect_identical(tabbed_ui$type, "tabset")
  expect_identical(tabbed_ui$id, "duplicates-dup_tables_tabs")
  expect_identical(
    vapply(tabbed_ui$children, `[[`, character(1), "title"),
    c("Plasma/Serum", "Urine")
  )
  expect_identical(
    tabbed_ui$children[[1]]$children[[2]]$outputId,
    "duplicates-dup_table_Plasma_Serum"
  )
  expect_identical(
    tabbed_ui$children[[2]]$children[[2]]$outputId,
    "duplicates-dup_table_Urine"
  )
})

test_that("metabolomics duplicate table renderer seam registers only populated assay tables", {
  output <- new.env(parent = emptyenv())
  detected_dup_list <- list(
    "Plasma/Serum" = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE),
    Control = NULL,
    Empty = data.frame(feature = character(), stringsAsFactors = FALSE),
    Urine = data.frame(feature = "m2", stringsAsFactors = FALSE)
  )

  visible <- withVisible(
    registerMetabDuplicateTableRenderers(
      dupList = detected_dup_list,
      output = output,
      renderDtFn = function(expr) {
        eval(substitute(expr), envir = parent.frame())
      },
      datatableFn = function(data, options, rownames, class) {
        list(
          data = data,
          options = options,
          rownames = rownames,
          class = class
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value,
    c("dup_table_Plasma_Serum", "dup_table_Urine")
  )
  expect_identical(
    ls(output),
    c("dup_table_Plasma_Serum", "dup_table_Urine")
  )
  expect_identical(
    output$dup_table_Plasma_Serum$data,
    detected_dup_list[["Plasma/Serum"]]
  )
  expect_identical(
    output$dup_table_Urine$data,
    detected_dup_list$Urine
  )
  expect_identical(
    output$dup_table_Plasma_Serum$options,
    list(pageLength = 10, scrollX = TRUE, dom = "frtip")
  )
  expect_false(output$dup_table_Plasma_Serum$rownames)
  expect_identical(output$dup_table_Plasma_Serum$class, "compact stripe")
})

test_that("metabolomics duplicate detection seam validates state and counts duplicates", {
  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  current_s4 <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m1"),
        Sample1 = c(1, 2),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "metabolite_id"
  )

  visible <- withVisible(
    detectMetabDuplicateFeatures(
      stateManager = list(
        getState = function() {
          captured$get_state_called <- TRUE
          current_s4
        }
      ),
      duplicateFinderFn = function(object) {
        captured$duplicate_finder_input <- object
        list(
          Plasma = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE),
          Control = NULL,
          Urine = data.frame(feature = "m2", stringsAsFactors = FALSE)
        )
      },
      reqFn = function(value) {
        captured$req_values[[length(captured$req_values) + 1L]] <- value
        value
      },
      inheritsFn = function(object, class_name) {
        captured$inherits_call <- list(
          object = object,
          class_name = class_name
        )
        inherits(object, class_name)
      }
    )
  )

  expect_true(visible$visible)
  expect_true(isTRUE(captured$get_state_called))
  expect_identical(length(captured$req_values), 2L)
  expect_identical(captured$duplicate_finder_input, current_s4)
  expect_identical(captured$inherits_call$class_name, "MetaboliteAssayData")
  expect_identical(names(visible$value$duplicatesList), c("Plasma", "Control", "Urine"))
  expect_identical(visible$value$totalDuplicates, 3L)
})

test_that("metabolomics duplicate detection reporting seam logs and notifies", {
  captured <- new.env(parent = emptyenv())
  captured$log_messages <- character()
  captured$notifications <- list()

  warning_visible <- withVisible(
    reportMetabDuplicateDetection(
      totalDuplicates = 2L,
      logInfoFn = function(message) {
        captured$log_messages <- c(captured$log_messages, message)
        invisible(NULL)
      },
      showNotificationFn = function(message, type) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type
        )
        invisible(NULL)
      }
    )
  )

  expect_false(warning_visible$visible)
  expect_identical(
    warning_visible$value$logMessage,
    "Detected 2 duplicate feature IDs across assays"
  )
  expect_identical(
    warning_visible$value$notificationMessage,
    "Detection complete: 2 duplicate IDs found"
  )
  expect_identical(warning_visible$value$notificationType, "warning")
  expect_identical(
    captured$log_messages,
    "Detected 2 duplicate feature IDs across assays"
  )
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Detection complete: 2 duplicate IDs found",
      type = "warning"
    )
  )

  message_visible <- withVisible(
    reportMetabDuplicateDetection(
      totalDuplicates = 0L,
      logInfoFn = function(message) {
        captured$log_messages <- c(captured$log_messages, message)
        invisible(NULL)
      },
      showNotificationFn = function(message, type) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type
        )
        invisible(NULL)
      }
    )
  )

  expect_false(message_visible$visible)
  expect_identical(message_visible$value$notificationType, "message")
  expect_identical(
    captured$notifications[[2L]],
    list(
      message = "Detection complete: 0 duplicate IDs found",
      type = "message"
    )
  )
})

test_that("metabolomics duplicate detection error seam logs and notifies", {
  captured <- new.env(parent = emptyenv())
  captured$log_messages <- character()
  captured$notifications <- list()

  visible <- withVisible(
    reportMetabDuplicateDetectionError(
      errorCondition = simpleError("detection exploded"),
      logErrorFn = function(message) {
        captured$log_messages <- c(captured$log_messages, message)
        invisible(NULL)
      },
      showNotificationFn = function(message, type) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    visible$value$errorMessage,
    "Error detecting duplicates: detection exploded"
  )
  expect_identical(visible$value$notificationType, "error")
  expect_identical(
    captured$log_messages,
    "Error detecting duplicates: detection exploded"
  )
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Error detecting duplicates: detection exploded",
      type = "error"
    )
  )
})

test_that("metabolomics duplicate detection observer shell updates duplicate state and reports success", {
  captured <- new.env(parent = emptyenv())
  captured$duplicate_info_updates <- list()
  detected_dup_list <- list(
    Plasma = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE)
  )
  state_manager <- list(label = "state-manager")

  visible <- withVisible(
    runMetabDuplicateDetectionObserverShell(
      stateManager = state_manager,
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates[[length(captured$duplicate_info_updates) + 1L]] <- value
        invisible(NULL)
      },
      detectDuplicatesFn = function(stateManager) {
        captured$detect_call <- list(
          state_manager = stateManager
        )
        list(
          duplicatesList = detected_dup_list,
          totalDuplicates = 2L
        )
      },
      reportDetectionFn = function(totalDuplicates) {
        captured$report_call <- list(
          total_duplicates = totalDuplicates
        )
        invisible(NULL)
      },
      reportDetectionErrorFn = function(errorCondition) {
        captured$error_call <- errorCondition
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$detect_call$state_manager, state_manager)
  expect_identical(captured$duplicate_info_updates, list(detected_dup_list))
  expect_identical(captured$report_call$total_duplicates, 2L)
  expect_false(exists("error_call", envir = captured, inherits = FALSE))
  expect_identical(visible$value$status, "success")
  expect_identical(visible$value$duplicatesList, detected_dup_list)
  expect_identical(visible$value$totalDuplicates, 2L)
})

test_that("metabolomics duplicate detection observer shell reports errors without mutating duplicate state", {
  captured <- new.env(parent = emptyenv())
  captured$duplicate_info_updates <- list()
  state_manager <- list(label = "state-manager")

  visible <- withVisible(
    runMetabDuplicateDetectionObserverShell(
      stateManager = state_manager,
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates[[length(captured$duplicate_info_updates) + 1L]] <- value
        invisible(NULL)
      },
      detectDuplicatesFn = function(stateManager) {
        captured$detect_call <- list(
          state_manager = stateManager
        )
        stop("detection exploded")
      },
      reportDetectionFn = function(totalDuplicates) {
        stop("reportDetectionFn should not be called on the error path")
      },
      reportDetectionErrorFn = function(errorCondition) {
        captured$error_call <- errorCondition
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$detect_call$state_manager, state_manager)
  expect_identical(captured$duplicate_info_updates, list())
  expect_s3_class(captured$error_call, "error")
  expect_identical(
    conditionMessage(captured$error_call),
    "detection exploded"
  )
  expect_identical(visible$value$status, "error")
  expect_s3_class(visible$value$errorCondition, "error")
  expect_identical(
    conditionMessage(visible$value$errorCondition),
    "detection exploded"
  )
})

test_that("metabolomics duplicate resolution observer shell renders results and clears duplicate state", {
  captured <- new.env(parent = emptyenv())
  captured$duplicate_info_updates <- list()
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabDuplicateResolutionObserverShell(
      runResolutionFn = function() {
        captured$run_resolution_called <- TRUE
        list(
          resultText = "resolution-summary-token",
          totalRemoved = 3L
        )
      },
      output = output,
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates <- c(
          captured$duplicate_info_updates,
          list(value)
        )
        invisible(NULL)
      },
      renderTextFn = function(text) text,
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error_log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, type, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- c(
          list(message = message, type = type),
          list(...)
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        captured$removed_notification <- id
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_true(isTRUE(captured$run_resolution_called))
  expect_identical(output$resolution_results, "resolution-summary-token")
  expect_identical(captured$duplicate_info_updates, list(NULL))
  expect_identical(captured$log_message, "Resolved duplicates: removed 3 rows")
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
  expect_identical(captured$removed_notification, "metab_dup_resolve_working")
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Duplicates resolved: 3 rows removed",
      type = "message"
    )
  )
  expect_identical(visible$value$status, "success")
  expect_identical(visible$value$totalRemoved, 3L)
  expect_identical(visible$value$resultText, "resolution-summary-token")
})

test_that("metabolomics duplicate resolution observer shell reports errors and tears down the working notification", {
  captured <- new.env(parent = emptyenv())
  captured$duplicate_info_updates <- list()
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabDuplicateResolutionObserverShell(
      runResolutionFn = function() {
        captured$run_resolution_called <- TRUE
        stop("resolution exploded")
      },
      output = output,
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates <- c(
          captured$duplicate_info_updates,
          list(value)
        )
        invisible(NULL)
      },
      renderTextFn = function(text) {
        stop("renderTextFn should not be called on the error path")
      },
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error_log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, type, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- c(
          list(message = message, type = type),
          list(...)
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        captured$removed_notification <- id
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_true(isTRUE(captured$run_resolution_called))
  expect_identical(captured$duplicate_info_updates, list())
  expect_false(exists("resolution_results", envir = output, inherits = FALSE))
  expect_false(exists("log_message", envir = captured, inherits = FALSE))
  expect_identical(
    captured$error_log_message,
    "Error resolving duplicates: resolution exploded"
  )
  expect_identical(captured$removed_notification, "metab_dup_resolve_working")
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Error resolving duplicates: resolution exploded",
      type = "error",
      duration = 15
    )
  )
  expect_identical(visible$value$status, "error")
  expect_s3_class(visible$value$errorCondition, "error")
  expect_identical(
    visible$value$errorMessage,
    "Error resolving duplicates: resolution exploded"
  )
})

test_that("metabolomics duplicate resolution observer seam reqs state, shows a working notification, and delegates through the shell workflow handoff", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$resolution_stats_updates <- list()
  captured$filter_plot_updates <- list()
  captured$duplicate_info_updates <- list()
  output <- new.env(parent = emptyenv())

  duplicate_info_setter <- function(value) {
    captured$duplicate_info_updates <- c(
      captured$duplicate_info_updates,
      list(value)
    )
    invisible(NULL)
  }
  workflow_data <- list(
    state_manager = structure(list(id = "state-manager"), class = "mock_state_manager"),
    config_list = list(config = "value")
  )

  visible <- withVisible(
    runMetabDuplicateResolutionObserver(
      workflowData = workflow_data,
      omicType = "metabolomics",
      output = output,
      setDuplicateInfoFn = duplicate_info_setter,
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        captured$req_value <- value
        value
      },
      showNotificationFn = function(ui, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- c(
          list(ui = ui),
          list(...)
        )
        invisible(NULL)
      },
      runResolutionObserverShellFn = function(runResolutionFn, output, setDuplicateInfoFn, ...) {
        captured$shell_call <- list(
          output = output,
          set_duplicate_info_fn = setDuplicateInfoFn
        )
        captured$resolution_dispatch <- runResolutionFn()
        invisible(list(status = "success"))
      },
      runResolutionWorkflowFn = function(workflowData, omicType, setResolutionStatsFn, setFilterPlotFn, ...) {
        captured$workflow_call <- list(
          workflow_data = workflowData,
          omic_type = omicType
        )
        setResolutionStatsFn("stats-token")
        setFilterPlotFn("plot-token")
        list(
          totalRemoved = 2L,
          resultText = "summary-token"
        )
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req_value, workflow_data$state_manager)
  expect_identical(captured$shell_call$output, output)
  expect_identical(captured$shell_call$set_duplicate_info_fn, duplicate_info_setter)
  expect_identical(captured$workflow_call$workflow_data, workflow_data)
  expect_identical(captured$workflow_call$omic_type, "metabolomics")
  expect_identical(captured$resolution_stats_updates, list("stats-token"))
  expect_identical(captured$filter_plot_updates, list("plot-token"))
  expect_identical(captured$duplicate_info_updates, list())
  expect_identical(
    captured$notifications[[1L]],
    list(
      ui = "Resolving duplicate features...",
      id = "metab_dup_resolve_working",
      duration = NULL
    )
  )
  expect_identical(captured$resolution_dispatch$totalRemoved, 2L)
  expect_identical(captured$resolution_dispatch$resultText, "summary-token")
})

test_that("metabolomics duplicate revert observer shell renders results and resets module state", {
  captured <- new.env(parent = emptyenv())
  captured$resolution_stats_updates <- list()
  captured$duplicate_info_updates <- list()
  captured$filter_plot_updates <- list()
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabDuplicateRevertObserverShell(
      runRevertFn = function() {
        captured$run_revert_called <- TRUE
        list(
          previousStateName = "after_import",
          resultText = "revert-result-token"
        )
      },
      output = output,
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates <- c(
          captured$duplicate_info_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      renderTextFn = function(text) text,
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error_log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, type, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- c(
          list(message = message, type = type),
          list(...)
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_true(isTRUE(captured$run_revert_called))
  expect_identical(output$resolution_results, "revert-result-token")
  expect_identical(captured$resolution_stats_updates, list(NULL))
  expect_identical(captured$duplicate_info_updates, list(NULL))
  expect_identical(captured$filter_plot_updates, list(NULL))
  expect_identical(
    captured$log_message,
    "Reverted duplicate resolution to after_import"
  )
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Reverted successfully",
      type = "message"
    )
  )
  expect_identical(visible$value$status, "success")
  expect_identical(visible$value$previousStateName, "after_import")
  expect_identical(visible$value$resultText, "revert-result-token")
})

test_that("metabolomics duplicate revert observer shell reports errors", {
  captured <- new.env(parent = emptyenv())
  captured$resolution_stats_updates <- list()
  captured$duplicate_info_updates <- list()
  captured$filter_plot_updates <- list()
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabDuplicateRevertObserverShell(
      runRevertFn = function() {
        captured$run_revert_called <- TRUE
        stop("revert exploded")
      },
      output = output,
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setDuplicateInfoFn = function(value) {
        captured$duplicate_info_updates <- c(
          captured$duplicate_info_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      renderTextFn = function(text) {
        stop("renderTextFn should not be called on the error path")
      },
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error_log_message <- message
        invisible(NULL)
      },
      showNotificationFn = function(message, type, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- c(
          list(message = message, type = type),
          list(...)
        )
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_true(isTRUE(captured$run_revert_called))
  expect_identical(captured$resolution_stats_updates, list())
  expect_identical(captured$duplicate_info_updates, list())
  expect_identical(captured$filter_plot_updates, list())
  expect_false(exists("resolution_results", envir = output, inherits = FALSE))
  expect_false(exists("log_message", envir = captured, inherits = FALSE))
  expect_identical(
    captured$error_log_message,
    "Error reverting: revert exploded"
  )
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Error reverting: revert exploded",
      type = "error"
    )
  )
  expect_identical(visible$value$status, "error")
  expect_s3_class(visible$value$errorCondition, "error")
  expect_identical(
    visible$value$errorMessage,
    "Error reverting: revert exploded"
  )
})

test_that("metabolomics duplicate revert seam selects the previous state and rejects missing history", {
  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  state_manager <- list(label = "state-manager")

  visible <- withVisible(
    revertMetabDuplicateResolution(
      stateManager = state_manager,
      reqFn = function(value) {
        captured$req_values[[length(captured$req_values) + 1L]] <- value
        value
      },
      historyGetterFn = function(manager) {
        captured$history_state_manager <- manager
        c("baseline", "after_import", "metab_duplicates_resolved")
      },
      revertStateFn = function(manager, stateName) {
        captured$revert_call <- list(
          state_manager = manager,
          state_name = stateName
        )
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(length(captured$req_values), 1L)
  expect_identical(captured$history_state_manager, state_manager)
  expect_identical(captured$revert_call$state_manager, state_manager)
  expect_identical(captured$revert_call$state_name, "after_import")
  expect_identical(visible$value$previousStateName, "after_import")
  expect_identical(
    visible$value$resultText,
    "Reverted to previous state: after_import"
  )

  expect_error(
    revertMetabDuplicateResolution(
      stateManager = state_manager,
      reqFn = function(value) value,
      historyGetterFn = function(manager) "baseline",
      revertStateFn = function(manager, stateName) {
        stop("revertStateFn should not be called")
      }
    ),
    "No previous state to revert to."
  )
})

test_that("metabolomics duplicate filter-plot seam dispatches grid and ggplot payloads", {
  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  grob_object <- structure(list(label = "grid"), class = c("gtable", "grob"))
  ggplot_object <- structure(list(label = "plot"), class = "ggplot")

  grob_visible <- withVisible(
    renderMetabDuplicateFilterPlot(
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
    renderMetabDuplicateFilterPlot(
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

test_that("metabolomics duplicate module resolve observer delegates through the resolution observer seam", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$reactive_values <- list()
  server_env <- environment(mod_metab_qc_duplicates_server)
  binding_names <- c(
    "runMetabDuplicateResolutionObserver",
    "runMetabDuplicateResolutionObserverShell",
    "runMetabDuplicateResolutionWorkflow",
    "prepareMetabDuplicateResolutionState",
    "applyMetabDuplicateResolutionState",
    "buildMetabDuplicateResolutionSummary",
    "buildMetabDuplicateSummaryUi",
    "updateMetaboliteFiltering"
  )
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)

  expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )

  assign(
    "runMetabDuplicateResolutionObserver",
    function(workflowData, omicType, output, setDuplicateInfoFn, setResolutionStatsFn, setFilterPlotFn, ...) {
      captured$observer_call <- list(
        workflow_data = workflowData,
        omic_type = omicType,
        output = output,
        set_duplicate_info_fn = setDuplicateInfoFn
      )
      setResolutionStatsFn(expected_stats)
      setFilterPlotFn("plot-token")
      output$resolution_results <- "summary-token"

      invisible(list(status = "success", resultText = "summary-token"))
    },
    envir = server_env
  )
  assign(
    "runMetabDuplicateResolutionWorkflow",
    function(...) {
      stop("runMetabDuplicateResolutionWorkflow should not be called directly")
    },
    envir = server_env
  )
  assign(
    "runMetabDuplicateResolutionObserverShell",
    function(...) {
      stop("runMetabDuplicateResolutionObserverShell should not be called directly")
    },
    envir = server_env
  )
  assign(
    "prepareMetabDuplicateResolutionState",
    function(...) {
      stop("prepareMetabDuplicateResolutionState should not be called directly")
    },
    envir = server_env
  )
  assign(
    "applyMetabDuplicateResolutionState",
    function(...) {
      stop("applyMetabDuplicateResolutionState should not be called directly")
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateResolutionSummary",
    function(...) {
      stop("buildMetabDuplicateResolutionSummary should not be called directly")
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateSummaryUi",
    function(dupList, ...) {
      captured$duplicate_summary_ui_call <- list(
        dup_list = dupList
      )
      "duplicate-summary-ui-token"
    },
    envir = server_env
  )
  assign(
    "updateMetaboliteFiltering",
    function(...) {
      stop("updateMetaboliteFiltering should not be called directly")
    },
    envir = server_env
  )

  state_manager <- list(
    getState = function() stop("getState should not be called directly"),
    saveState = function(...) stop("saveState should not be called directly"),
    getHistory = function() "baseline",
    revertToState = function(state_name) invisible(state_name)
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          detect_duplicates = FALSE,
          resolve_duplicates = TRUE,
          revert_duplicates = FALSE
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      store <- new.env(parent = emptyenv())
      store$set_count <- 0L
      store$last_value <- value
      captured$reactive_values[[length(captured$reactive_values) + 1L]] <- store

      function(new_value) {
        if (missing(new_value)) {
          store$last_value
        } else {
          store$set_count <- store$set_count + 1L
          store$last_value <- new_value
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
    observe = function(expr, ...) invisible(NULL),
    req = function(...) list(...)[[1]],
    showNotification = function(ui, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(ui = ui, ...)
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) substitute(expr),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(
    state_manager = state_manager,
    config_list = list(config = "value")
  )

  mod_metab_qc_duplicates_server(
    id = "duplicates",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "duplicates")
  expect_match(
    paste(deparse(captured$output$duplicate_summary), collapse = " "),
    "buildMetabDuplicateSummaryUi",
    fixed = TRUE
  )
  expect_match(
    paste(deparse(captured$output$duplicate_tables), collapse = " "),
    "buildMetabDuplicateTablesUi",
    fixed = TRUE
  )
  expect_identical(captured$observer_call$output, captured$output)
  expect_identical(captured$observer_call$workflow_data, workflow_data)
  expect_identical(captured$observer_call$omic_type, "metabolomics")
  expect_identical(length(captured$reactive_values), 3L)
  expect_identical(captured$reactive_values[[1L]]$set_count, 0L)
  expect_null(captured$reactive_values[[1L]]$last_value)
  expect_identical(captured$reactive_values[[2L]]$set_count, 1L)
  expect_identical(captured$reactive_values[[2L]]$last_value, expected_stats)
  expect_identical(captured$reactive_values[[3L]]$set_count, 1L)
  expect_identical(captured$reactive_values[[3L]]$last_value, "plot-token")
  expect_identical(captured$output$resolution_results, "summary-token")
  expect_identical(captured$notifications, list())
})

test_that("metabolomics duplicate module detect observer delegates through the detection observer shell seam", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  server_env <- environment(mod_metab_qc_duplicates_server)
  binding_names <- c(
    "runMetabDuplicateDetectionObserverShell",
    "buildMetabDuplicateSummaryUi",
    "buildMetabDuplicateTablesUi",
    "registerMetabDuplicateTableRenderers"
  )
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)

  detected_dup_list <- list(
    "Plasma/Serum" = data.frame(feature = c("m1", "m1"), stringsAsFactors = FALSE),
    Control = NULL
  )

  assign(
    "runMetabDuplicateDetectionObserverShell",
    function(stateManager, setDuplicateInfoFn, ...) {
      captured$shell_call <- list(
        state_manager = stateManager
      )
      setDuplicateInfoFn(detected_dup_list)
      invisible(list(status = "success"))
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateSummaryUi",
    function(dupList, ...) {
      captured$summary_ui_call <- dupList
      "duplicate-summary-ui-token"
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateTablesUi",
    function(dupList, nsFn, ...) {
      captured$tables_ui_call <- list(
        dup_list = dupList,
        sample_output_id = nsFn("dup_table_test")
      )
      "duplicate-tables-ui-token"
    },
    envir = server_env
  )
  assign(
    "registerMetabDuplicateTableRenderers",
    function(dupList, output, ...) {
      captured$renderer_call <- list(
        dup_list = dupList,
        output = output
      )
      invisible("registered")
    },
    envir = server_env
  )

  state_manager <- list(
    getState = function() stop("getState should not be called directly"),
    saveState = function(...) invisible(NULL),
    getHistory = function() "baseline",
    revertToState = function(state_name) invisible(state_name)
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          detect_duplicates = TRUE,
          resolve_duplicates = FALSE,
          revert_duplicates = FALSE
        ),
        output,
        list(ns = function(value) paste0("duplicates-", value))
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    observe = function(expr, ...) {
      eval(substitute(expr), parent.frame())
      invisible(NULL)
    },
    req = function(...) list(...)[[1]],
    showNotification = function(ui, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(ui = ui, ...)
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) eval(substitute(expr), parent.frame()),
    renderPlot = function(expr) substitute(expr),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(
    state_manager = state_manager,
    config_list = list(config = "value")
  )

  mod_metab_qc_duplicates_server(
    id = "duplicates",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "duplicates")
  expect_identical(captured$shell_call$state_manager, state_manager)
  expect_identical(captured$summary_ui_call, detected_dup_list)
  expect_identical(captured$tables_ui_call$dup_list, detected_dup_list)
  expect_identical(captured$renderer_call$dup_list, detected_dup_list)
  expect_identical(captured$renderer_call$output, captured$output)
  expect_identical(captured$output$duplicate_summary, "duplicate-summary-ui-token")
  expect_identical(captured$output$duplicate_tables, "duplicate-tables-ui-token")
})

test_that("metabolomics duplicate module detect observer error path delegates through the detection observer shell seam", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  server_env <- environment(mod_metab_qc_duplicates_server)
  binding_names <- c(
    "runMetabDuplicateDetectionObserverShell",
    "buildMetabDuplicateSummaryUi",
    "buildMetabDuplicateTablesUi",
    "registerMetabDuplicateTableRenderers"
  )
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)

  assign(
    "runMetabDuplicateDetectionObserverShell",
    function(stateManager, setDuplicateInfoFn, ...) {
      captured$shell_call <- list(
        state_manager = stateManager
      )
      invisible(list(status = "error"))
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateSummaryUi",
    function(dupList, ...) {
      captured$summary_ui_call <- dupList
      "duplicate-summary-ui-token"
    },
    envir = server_env
  )
  assign(
    "buildMetabDuplicateTablesUi",
    function(dupList, nsFn, ...) {
      captured$tables_ui_call <- list(
        dup_list = dupList,
        sample_output_id = nsFn("dup_table_test")
      )
      "duplicate-tables-ui-token"
    },
    envir = server_env
  )
  assign(
    "registerMetabDuplicateTableRenderers",
    function(dupList, output, ...) {
      captured$renderer_call <- list(
        dup_list = dupList,
        output = output
      )
      invisible("registered")
    },
    envir = server_env
  )

  state_manager <- list(
    getState = function() stop("getState should not be called directly"),
    saveState = function(...) invisible(NULL),
    getHistory = function() "baseline",
    revertToState = function(state_name) invisible(state_name)
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          detect_duplicates = TRUE,
          resolve_duplicates = FALSE,
          revert_duplicates = FALSE
        ),
        output,
        list(ns = function(value) paste0("duplicates-", value))
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    observe = function(expr, ...) {
      eval(substitute(expr), parent.frame())
      invisible(NULL)
    },
    req = function(...) list(...)[[1]],
    showNotification = function(ui, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(ui = ui, ...)
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) eval(substitute(expr), parent.frame()),
    renderPlot = function(expr) substitute(expr),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(
    state_manager = state_manager,
    config_list = list(config = "value")
  )

  mod_metab_qc_duplicates_server(
    id = "duplicates",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "duplicates")
  expect_identical(captured$shell_call$state_manager, state_manager)
  expect_identical(captured$summary_ui_call, NULL)
  expect_identical(captured$tables_ui_call$dup_list, NULL)
  expect_identical(captured$renderer_call$dup_list, NULL)
  expect_identical(captured$renderer_call$output, captured$output)
  expect_identical(captured$output$duplicate_summary, "duplicate-summary-ui-token")
  expect_identical(captured$output$duplicate_tables, "duplicate-tables-ui-token")
})

test_that("metabolomics duplicate module revert observer delegates through the revert observer shell seam", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$reactive_values <- list()
  server_env <- environment(mod_metab_qc_duplicates_server)
  binding_names <- c(
    "runMetabDuplicateRevertObserverShell",
    "revertMetabDuplicateResolution"
  )
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)

  state_manager <- list(
    getHistory = function() stop("getHistory should not be called directly"),
    revertToState = function(state_name) stop("revertToState should not be called directly")
  )

  assign(
    "runMetabDuplicateRevertObserverShell",
    function(runRevertFn, output, setResolutionStatsFn, setDuplicateInfoFn, setFilterPlotFn, ...) {
      captured$shell_call <- list(
        output = output,
        set_resolution_stats_fn = setResolutionStatsFn,
        set_duplicate_info_fn = setDuplicateInfoFn,
        set_filter_plot_fn = setFilterPlotFn
      )
      captured$revert_dispatch <- runRevertFn()
      output$resolution_results <- captured$revert_dispatch$resultText
      setResolutionStatsFn(NULL)
      setDuplicateInfoFn(NULL)
      setFilterPlotFn(NULL)
      invisible(list(status = "success"))
    },
    envir = server_env
  )
  assign(
    "revertMetabDuplicateResolution",
    function(stateManager, ...) {
      captured$revert_call <- list(
        state_manager = stateManager
      )
      list(
        previousStateName = "after_import",
        resultText = "revert-result-token"
      )
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          detect_duplicates = FALSE,
          resolve_duplicates = FALSE,
          revert_duplicates = TRUE
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      store <- new.env(parent = emptyenv())
      store$set_count <- 0L
      store$last_value <- value
      captured$reactive_values[[length(captured$reactive_values) + 1L]] <- store

      function(new_value) {
        if (missing(new_value)) {
          store$last_value
        } else {
          store$set_count <- store$set_count + 1L
          store$last_value <- new_value
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
    observe = function(expr, ...) invisible(NULL),
    req = function(...) list(...)[[1]],
    showNotification = function(ui, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(ui = ui, ...)
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) substitute(expr),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(
    state_manager = state_manager,
    config_list = list(config = "value")
  )

  mod_metab_qc_duplicates_server(
    id = "duplicates",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "duplicates")
  expect_identical(captured$shell_call$output, captured$output)
  expect_identical(captured$revert_call$state_manager, state_manager)
  expect_identical(captured$revert_dispatch$previousStateName, "after_import")
  expect_identical(captured$revert_dispatch$resultText, "revert-result-token")
  expect_identical(captured$output$resolution_results, "revert-result-token")
  expect_identical(length(captured$reactive_values), 3L)
  expect_true(all(vapply(captured$reactive_values, function(store) {
    identical(store$set_count, 1L) && is.null(store$last_value)
  }, logical(1))))
})

test_that("metabolomics duplicate module filter-plot render delegates through the filter-plot seam", {
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_duplicates_server)
  binding_names <- c("renderMetabDuplicateFilterPlot")
  had_bindings <- stats::setNames(
    vapply(
      binding_names,
      function(name) exists(name, envir = server_env, inherits = FALSE),
      logical(1)
    ),
    binding_names
  )
  old_bindings <- vector("list", length(binding_names))
  names(old_bindings) <- binding_names

  for (name in binding_names) {
    if (had_bindings[[name]]) {
      old_bindings[[name]] <- get(name, envir = server_env, inherits = FALSE)
    }
  }

  on.exit({
    for (name in binding_names) {
      if (had_bindings[[name]]) {
        assign(name, old_bindings[[name]], envir = server_env)
      } else if (exists(name, envir = server_env, inherits = FALSE)) {
        rm(list = name, envir = server_env)
      }
    }
  }, add = TRUE)

  assign(
    "renderMetabDuplicateFilterPlot",
    function(filterPlot, ...) {
      captured$filter_plot_call <- list(
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
        list(
          detect_duplicates = FALSE,
          resolve_duplicates = FALSE,
          revert_duplicates = FALSE
        ),
        output,
        list(ns = function(value) value)
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
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    observe = function(expr, ...) invisible(NULL),
    req = function(...) list(...)[[1]],
    showNotification = function(ui, ...) invisible(NULL),
    removeNotification = function(id, ...) invisible(NULL),
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) eval(substitute(expr), parent.frame()),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    .package = "logger"
  )

  workflow_data <- list(
    state_manager = list(
      getState = function() NULL,
      saveState = function(...) invisible(NULL),
      getHistory = function() "baseline",
      revertToState = function(state_name) invisible(state_name)
    ),
    config_list = list(config = "value")
  )

  mod_metab_qc_duplicates_server(
    id = "duplicates",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "duplicates")
  expect_type(captured$filter_plot_call$filter_plot, "closure")
  expect_identical(captured$output$filter_plot, "filter-plot-render-token")
})
