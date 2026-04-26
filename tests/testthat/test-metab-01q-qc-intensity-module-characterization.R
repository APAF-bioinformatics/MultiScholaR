# fidelity-coverage-compare: shared
library(testthat)

source("helpers-scoped-mocked-bindings.R")

register_binding_teardown(
  asNamespace("shiny"),
  c(
    "moduleServer",
    "reactiveVal",
    "observeEvent",
    "renderText",
    "renderUI",
    "renderPlot",
    "req",
    "showNotification",
    "removeNotification"
  ),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("logger"),
  c("log_info", "log_error", "log_warn"),
  .local_envir = environment()
)
register_binding_teardown(
  asNamespace("MultiScholaR"),
  c("metaboliteIntensityFiltering", "updateMetaboliteFiltering"),
  .local_envir = environment()
)

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
    file.path(repo_root, "R", "mod_metab_qc_intensity_ui_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_intensity.R")
  ),
  symbols = c(
    "buildMetabIntensityAssayTabsUi",
    "renderMetabIntensityFilterPlot",
    "buildMetabIntensityFilterStats",
    "updateMetabIntensityFilterQcPlot",
    "saveMetabIntensityFilterState",
    "buildMetabIntensityFilterSummary",
    "reportMetabIntensityFilterSuccess",
    "reportMetabIntensityFilterRevertSuccess",
    "mod_metab_qc_intensity_server"
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

.metab_intensity_test_env <- environment()

makeMetabIntensityAssayData <- function(metabolite_data, metabolite_id_column = "metabolite_id") {
  sample_columns <- unique(unlist(lapply(
    metabolite_data,
    function(assay) setdiff(names(assay), metabolite_id_column)
  )))
  sample_columns <- sort(sample_columns)

  args <- list(
    Class = "MetaboliteAssayData",
    metabolite_data = metabolite_data,
    metabolite_id_column = metabolite_id_column
  )

  slots <- methods::slotNames("MetaboliteAssayData")
  if ("annotation_id_column" %in% slots) {
    args$annotation_id_column <- "annotation_id"
  }
  if ("database_identifier_type" %in% slots) {
    args$database_identifier_type <- "test_id"
  }
  if ("internal_standard_regex" %in% slots) {
    args$internal_standard_regex <- NA_character_
  }
  if ("design_matrix" %in% slots) {
    args$design_matrix <- data.frame(
      Sample_ID = sample_columns,
      group = "A",
      stringsAsFactors = FALSE
    )
  }
  if ("sample_id" %in% slots) {
    args$sample_id <- "Sample_ID"
  }
  if ("group_id" %in% slots) {
    args$group_id <- "group"
  }
  if ("technical_replicate_id" %in% slots) {
    args$technical_replicate_id <- NA_character_
  }
  if ("args" %in% slots) {
    args$args <- list()
  }

  do.call(methods::new, args)
}

metabIntensityHasExtractedSeams <- function() {
  all(vapply(
    c(
      "buildMetabIntensityAssayTabsUi",
      "renderMetabIntensityFilterPlot",
      "buildMetabIntensityFilterStats",
      "updateMetabIntensityFilterQcPlot",
      "saveMetabIntensityFilterState",
      "buildMetabIntensityFilterSummary",
      "reportMetabIntensityFilterSuccess",
      "reportMetabIntensityFilterRevertSuccess"
    ),
    exists,
    logical(1),
    envir = .metab_intensity_test_env,
    inherits = FALSE
  ))
}

test_that("metabolomics intensity module preserves apply filter public behavior", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  package_ns <- asNamespace("MultiScholaR")
  server_fn <- if (exists(
    "mod_metab_qc_intensity_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    get("mod_metab_qc_intensity_server", envir = package_ns, inherits = FALSE)
  } else {
    mod_metab_qc_intensity_server
  }
  server_env <- environment(server_fn)

  current_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m2"),
        Sample1 = c(10, 20),
        Sample2 = c(5, 25),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = "m3",
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        Sample2 = 25,
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = "m3",
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )

  scoped_mocked_bindings(
    metaboliteIntensityFiltering = function(
      theObject,
      metabolites_intensity_cutoff_percentile,
      metabolites_proportion_of_samples_below_cutoff
    ) {
      captured$filter_call <- list(
        the_object = theObject,
        intensity_cutoff_percentile = metabolites_intensity_cutoff_percentile,
        proportion_below_cutoff = metabolites_proportion_of_samples_below_cutoff
      )
      filtered_s4
    },
    updateMetaboliteFiltering = function(theObject, step_name, omics_type, return_grid, overwrite) {
      captured$qc_update <- list(
        the_object = theObject,
        step_name = step_name,
        omics_type = omics_type,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot-token"
    },
    .env = server_env
  )

  state_manager <- list(
    getState = function() current_s4,
    saveState = function(...) {
      captured$save_state_call <- list(...)
      invisible(NULL)
    },
    getHistory = function() "baseline",
    revertToState = function(state_name) invisible(state_name)
  )

  scoped_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      module(
        list(
          apply_filter = TRUE,
          revert_filter = FALSE,
          intensity_cutoff_percentile = 10,
          proportion_below_cutoff = 0.5
        ),
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    req = function(...) list(...)[[1]],
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(message) {
      captured$log_message <- message
      invisible(NULL)
    },
    log_error = function(message) {
      captured$error_log_message <- message
      invisible(NULL)
    },
    log_warn = function(message) {
      captured$warning_log_message <- message
      invisible(NULL)
    },
    .env = asNamespace("logger")
  )

  server_fn(
    id = "intensity",
    workflow_data = list(
      state_manager = state_manager,
      config_list = list(config = "value")
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$filter_call$the_object, current_s4)
  expect_identical(captured$filter_call$intensity_cutoff_percentile, 10)
  expect_identical(captured$filter_call$proportion_below_cutoff, 0.5)
  expect_identical(captured$qc_update$the_object, filtered_s4)
  expect_identical(captured$qc_update$step_name, "2_Intensity_Filtered")
  expect_identical(captured$qc_update$omics_type, "metabolomics")
  expect_true(isTRUE(captured$qc_update$return_grid))
  expect_true(isTRUE(captured$qc_update$overwrite))
  expect_identical(captured$save_state_call$state_name, "metab_intensity_filtered")
  expect_identical(captured$save_state_call$s4_data_object, filtered_s4)
  expect_identical(captured$save_state_call$config_object, list(config = "value"))
  expect_match(captured$output$filter_results, "Metabolite Intensity Filter Applied Successfully", fixed = TRUE)
  expect_match(captured$output$filter_results, "Plasma: 2 -> 1", fixed = TRUE)
  expect_false(is.null(captured$output$assay_results_tabs))
  expect_identical(captured$removed_notification, "metab_intensity_filter_working")
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
})

test_that("metabolomics intensity module preserves revert public behavior", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  package_ns <- asNamespace("MultiScholaR")
  server_fn <- if (exists(
    "mod_metab_qc_intensity_server",
    envir = package_ns,
    inherits = FALSE
  )) {
    get("mod_metab_qc_intensity_server", envir = package_ns, inherits = FALSE)
  } else {
    mod_metab_qc_intensity_server
  }

  scoped_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      module(
        list(
          apply_filter = FALSE,
          revert_filter = TRUE,
          intensity_cutoff_percentile = 10,
          proportion_below_cutoff = 0.5
        ),
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
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    req = function(...) list(...)[[1]],
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(message) {
      captured$log_message <- message
      invisible(NULL)
    },
    log_error = function(message) {
      captured$error_log_message <- message
      invisible(NULL)
    },
    .env = asNamespace("logger")
  )

  state_manager <- list(
    getState = function() NULL,
    saveState = function(...) invisible(NULL),
    getHistory = function() c("baseline_metab_state", "filtered_metab_state"),
    revertToState = function(state_name) {
      captured$reverted_state_name <- state_name
      invisible(state_name)
    }
  )

  server_fn(
    id = "intensity",
    workflow_data = list(
      state_manager = state_manager,
      config_list = list(config = "value")
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$reverted_state_name, "baseline_metab_state")
  expect_match(
    captured$output$filter_results,
    "Reverted to previous state: baseline_metab_state",
    fixed = TRUE
  )
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
})

test_that("metabolomics intensity assay tabs seam returns null for empty stats and builds namespaced assay tabs", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  expect_null(
    buildMetabIntensityAssayTabsUi(
      stats = NULL,
      ns = function(value) value
    )
  )

  stats <- data.frame(
    Assay = c("Plasma", "Serum"),
    Original = c(10, 6),
    Filtered = c(8, 5),
    Removed = c(2, 1),
    Percent_Retained = c(80, 83.3),
    stringsAsFactors = FALSE
  )

  tabs_ui <- buildMetabIntensityAssayTabsUi(
    stats = stats,
    ns = function(value) paste0("ns-", value),
    tabPanelFn = function(title, ...) {
      list(type = "tab", title = title, children = list(...))
    },
    brFn = function() {
      list(type = "br")
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
    },
    tabsetPanelFn = function(id, ...) {
      list(type = "tabset", id = id, tabs = list(...))
    }
  )

  expect_identical(tabs_ui$type, "tabset")
  expect_identical(tabs_ui$id, "ns-assay_stats_tabs")
  expect_identical(
    vapply(tabs_ui$tabs, `[[`, character(1), "title"),
    c("Plasma", "Serum")
  )
  expect_identical(
    tabs_ui$tabs[[1L]]$children[[2L]]$body[[1L]]$rows[[1L]]$cells[[1L]]$children[[1L]]$text,
    "Original metabolites:"
  )
  expect_identical(
    tabs_ui$tabs[[1L]]$children[[2L]]$body[[1L]]$rows[[1L]]$cells[[2L]]$children[[1L]],
    10
  )
  expect_identical(
    tabs_ui$tabs[[1L]]$children[[2L]]$body[[1L]]$rows[[4L]]$cells[[2L]]$children[[1L]],
    "80%"
  )
})

test_that("metabolomics intensity module server delegates the assay tabs render path through the seam", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_intensity_server)
  binding_name <- "buildMetabIntensityAssayTabsUi"
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

  stats <- data.frame(
    Assay = "Plasma",
    Original = 10,
    Filtered = 8,
    Removed = 2,
    Percent_Retained = 80,
    stringsAsFactors = FALSE
  )
  reactive_values <- list(NULL, stats)
  reactive_index <- 0L

  assign(
    binding_name,
    function(stats, ns, ...) {
      captured$stats <- stats
      captured$tabset_id <- ns("assay_stats_tabs")
      "assay-tabs-token"
    },
    envir = server_env
  )

  scoped_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      stored <- reactive_values[[reactive_index]]

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
      invisible(NULL)
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      substitute(expr)
    },
    req = function(...) {
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .env = asNamespace("logger")
  )

  mod_metab_qc_intensity_server(
    id = "intensity",
    workflow_data = list(state_manager = NULL),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "intensity")
  expect_identical(captured$stats, stats)
  expect_identical(captured$tabset_id, "intensity-assay_stats_tabs")
  expect_identical(output$assay_results_tabs, "assay-tabs-token")
})

test_that("metabolomics intensity filter-plot seam dispatches grid and ggplot payloads", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  captured$req_values <- list()
  grob_object <- structure(list(label = "grid"), class = c("gtable", "grob"))
  ggplot_object <- structure(list(label = "plot"), class = "ggplot")

  grob_visible <- withVisible(
    renderMetabIntensityFilterPlot(
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
    renderMetabIntensityFilterPlot(
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

test_that("metabolomics intensity stats seam counts per-assay metabolites and retention", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  current_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m1", "m2"),
        Sample1 = c(10, 10, 20),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = c("s1", "s2"),
        Sample1 = c(7, 9),
        stringsAsFactors = FALSE
      )
    )
  )
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = "s2",
        Sample1 = 9,
        stringsAsFactors = FALSE
      )
    )
  )

  stats <- buildMetabIntensityFilterStats(
    currentS4 = current_s4,
    filteredS4 = filtered_s4
  )

  expect_equal(
    stats,
    data.frame(
      Assay = c("Plasma", "Serum"),
      Original = c(2, 2),
      Filtered = c(1, 1),
      Removed = c(1, 1),
      Percent_Retained = c(50, 50),
      stringsAsFactors = FALSE
    )
  )
})

test_that("metabolomics intensity QC plot seam refreshes the plot and stores it", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  captured$plot_updates <- list()
  captured$warnings <- character()
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        stringsAsFactors = FALSE
      )
    )
  )

  visible <- withVisible(
    updateMetabIntensityFilterQcPlot(
      filteredS4 = filtered_s4,
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
        "plot-token"
      },
      logWarnFn = function(message) {
        captured$warnings <- c(captured$warnings, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$qc_update$the_object, filtered_s4)
  expect_identical(captured$qc_update$step_name, "2_Intensity_Filtered")
  expect_identical(captured$qc_update$omics_type, "metabolomics")
  expect_true(isTRUE(captured$qc_update$return_grid))
  expect_true(isTRUE(captured$qc_update$overwrite))
  expect_identical(captured$plot_updates, list("plot-token"))
  expect_identical(captured$warnings, character())
  expect_identical(visible$value, "plot-token")
})

test_that("metabolomics intensity QC plot seam falls back to a null plot on refresh errors", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  captured$plot_updates <- list()
  captured$warnings <- character()
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        stringsAsFactors = FALSE
      )
    )
  )

  visible <- withVisible(
    updateMetabIntensityFilterQcPlot(
      filteredS4 = filtered_s4,
      omicType = "metabolomics",
      setFilterPlotFn = function(value) {
        captured$plot_updates <- c(captured$plot_updates, list(value))
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

  expect_false(visible$visible)
  expect_identical(
    captured$warnings,
    "Could not generate QC plot: plot refresh failed"
  )
  expect_length(captured$plot_updates, 1L)
  expect_null(captured$plot_updates[[1L]])
  expect_null(visible$value)
})

test_that("metabolomics intensity state-save seam persists filtered state and returns the saved state name", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        stringsAsFactors = FALSE
      )
    )
  )
  state_manager <- list(
    saveState = function(...) {
      captured$save_state_call <- list(...)
      invisible(NULL)
    }
  )
  config_object <- list(config = "value")

  visible <- withVisible(
    saveMetabIntensityFilterState(
      stateManager = state_manager,
      filteredS4 = filtered_s4,
      configObject = config_object,
      intensityCutoffPercentile = 10,
      proportionBelowCutoff = 0.5
    )
  )

  expect_false(visible$visible)
  expect_identical(visible$value, "metab_intensity_filtered")
  expect_identical(captured$save_state_call$state_name, "metab_intensity_filtered")
  expect_identical(captured$save_state_call$s4_data_object, filtered_s4)
  expect_identical(captured$save_state_call$config_object, config_object)
  expect_identical(
    captured$save_state_call$description,
    "Applied metabolite intensity filter (percentile: 10%, proportion: 0.50)"
  )
})

test_that("metabolomics intensity summary seam assembles totals and result text", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  stats <- data.frame(
    Assay = c("Plasma", "Serum"),
    Original = c(10, 6),
    Filtered = c(8, 5),
    Removed = c(2, 1),
    Percent_Retained = c(80, 83.3),
    stringsAsFactors = FALSE
  )

  summary <- buildMetabIntensityFilterSummary(
    stats = stats,
    intensityCutoffPercentile = 10,
    proportionBelowCutoff = 0.5,
    stateName = "custom_metab_state"
  )

  expect_equal(summary$totalOriginal, 16)
  expect_equal(summary$totalFiltered, 13)
  expect_equal(summary$totalRemoved, 3)
  expect_match(
    summary$resultText,
    "Metabolite Intensity Filter Applied Successfully",
    fixed = TRUE
  )
  expect_match(summary$resultText, "Intensity cutoff percentile: 10%", fixed = TRUE)
  expect_match(summary$resultText, "Max proportion below cutoff: 0.50", fixed = TRUE)
  expect_match(
    summary$resultText,
    "  Plasma: 10 -> 8 \\(removed 2, 80\\.0% retained\\)"
  )
  expect_match(
    summary$resultText,
    "  Serum: 6 -> 5 \\(removed 1, 83\\.3% retained\\)"
  )
  expect_match(
    summary$resultText,
    "Total: 16 -> 13 metabolites \\(removed 3\\)"
  )
  expect_match(summary$resultText, "State saved as: 'custom_metab_state'", fixed = TRUE)
})

test_that("metabolomics intensity success-reporting seam renders results and reports completion", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  stats <- data.frame(
    Assay = c("Plasma", "Serum"),
    Original = c(10, 6),
    Filtered = c(8, 5),
    Removed = c(2, 1),
    Percent_Retained = c(80, 83.3),
    stringsAsFactors = FALSE
  )

  visible <- withVisible(
    reportMetabIntensityFilterSuccess(
      stats = stats,
      intensityCutoffPercentile = 10,
      proportionBelowCutoff = 0.5,
      savedStateName = "saved-metab-state",
      output = output,
      buildSummaryFn = function(stats, intensityCutoffPercentile, proportionBelowCutoff, stateName) {
        captured$summary_call <- list(
          stats = stats,
          intensity_cutoff_percentile = intensityCutoffPercentile,
          proportion_below_cutoff = proportionBelowCutoff,
          state_name = stateName
        )
        list(
          totalOriginal = sum(stats$Original),
          totalFiltered = sum(stats$Filtered),
          totalRemoved = sum(stats$Removed),
          resultText = "intensity-summary-token"
        )
      },
      renderTextFn = function(text) {
        captured$render_text_value <- text
        paste0("rendered::", text)
      },
      logInfoFn = function(message) {
        captured$log_message <- message
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        captured$removed_notification <- id
        invisible(NULL)
      },
      showNotificationFn = function(message, ...) {
        captured$notification <- c(list(message = message), list(...))
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$summary_call$stats, stats)
  expect_identical(captured$summary_call$intensity_cutoff_percentile, 10)
  expect_identical(captured$summary_call$proportion_below_cutoff, 0.5)
  expect_identical(captured$summary_call$state_name, "saved-metab-state")
  expect_identical(captured$render_text_value, "intensity-summary-token")
  expect_identical(output$filter_results, "rendered::intensity-summary-token")
  expect_identical(
    captured$log_message,
    "Metabolite intensity filter applied: removed 3 metabolites"
  )
  expect_identical(captured$removed_notification, "metab_intensity_filter_working")
  expect_identical(
    captured$notification,
    list(
      message = "Intensity filter applied: 13 metabolites retained",
      type = "message"
    )
  )
  expect_equal(visible$value$totalOriginal, 16)
  expect_equal(visible$value$totalFiltered, 13)
  expect_equal(visible$value$totalRemoved, 3)
})

test_that("metabolomics intensity revert-reporting seam resets outputs and reports completion", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    reportMetabIntensityFilterRevertSuccess(
      prevStateName = "baseline_metab_state",
      output = output,
      setFilterStatsFn = function(value) {
        captured$stats_reset_value <- value
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$plot_reset_value <- value
        invisible(NULL)
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
  expect_identical(
    visible$value,
    "Reverted to previous state: baseline_metab_state"
  )
  expect_identical(
    captured$render_text_value,
    "Reverted to previous state: baseline_metab_state"
  )
  expect_identical(
    output$filter_results,
    "rendered::Reverted to previous state: baseline_metab_state"
  )
  expect_null(captured$stats_reset_value)
  expect_null(captured$plot_reset_value)
  expect_identical(
    captured$log_message,
    "Reverted metabolite intensity filter to baseline_metab_state"
  )
  expect_identical(
    captured$notification,
    list(
      message = "Reverted successfully",
      type = "message"
    )
  )
})

test_that("metabolomics intensity module server delegates the filter-plot render path through the seam", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_intensity_server)
  binding_name <- "renderMetabIntensityFilterPlot"
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
      captured$filter_plot_call <- list(
        filter_plot = filterPlot
      )
      "filter-plot-render-token"
    },
    envir = server_env
  )

  scoped_mocked_bindings(
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
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    renderUI = function(expr) substitute(expr),
    renderText = function(expr) substitute(expr),
    renderPlot = function(expr) eval(substitute(expr), parent.frame()),
    req = function(...) list(...)[[1]],
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .env = asNamespace("logger")
  )

  mod_metab_qc_intensity_server(
    id = "intensity",
    workflow_data = list(
      state_manager = list(
        getState = function() NULL,
        saveState = function(...) invisible(NULL),
        getHistory = function() "baseline",
        revertToState = function(state_name) invisible(state_name)
      ),
      config_list = list(config = "value")
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "intensity")
  expect_type(captured$filter_plot_call$filter_plot, "closure")
  expect_identical(captured$output$filter_plot, "filter-plot-render-token")
})

test_that("metabolomics intensity module apply observer delegates QC plot refresh, state save, and success reporting through seams", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$reactive_values <- list()
  server_env <- environment(mod_metab_qc_intensity_server)
  binding_names <- c(
    "buildMetabIntensityFilterStats",
    "updateMetabIntensityFilterQcPlot",
    "saveMetabIntensityFilterState",
    "reportMetabIntensityFilterSuccess",
    "metaboliteIntensityFiltering"
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

  current_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = c("m1", "m2"),
        Sample1 = c(10, 20),
        Sample2 = c(5, 25),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = c("m3"),
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )
  filtered_s4 <- makeMetabIntensityAssayData(
    metabolite_data = list(
      Plasma = data.frame(
        metabolite_id = "m2",
        Sample1 = 20,
        Sample2 = 25,
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        metabolite_id = "m3",
        Sample1 = 7,
        Sample2 = 9,
        stringsAsFactors = FALSE
      )
    )
  )
  expected_stats <- data.frame(
    Assay = c("Plasma", "Serum"),
    Original = c(2, 1),
    Filtered = c(1, 1),
    Removed = c(1, 0),
    Percent_Retained = c(50, 100),
    stringsAsFactors = FALSE
  )

  assign(
    "metaboliteIntensityFiltering",
    function(
      theObject,
      metabolites_intensity_cutoff_percentile,
      metabolites_proportion_of_samples_below_cutoff
    ) {
      captured$filter_call <- list(
        the_object = theObject,
        intensity_cutoff_percentile = metabolites_intensity_cutoff_percentile,
        proportion_below_cutoff = metabolites_proportion_of_samples_below_cutoff
      )
      filtered_s4
    },
    envir = server_env
  )
  assign(
    "updateMetabIntensityFilterQcPlot",
    function(filteredS4, omicType, setFilterPlotFn, ...) {
      captured$qc_plot_helper_call <- list(
        filtered_s4 = filteredS4,
        omic_type = omicType,
        set_filter_plot_fn = setFilterPlotFn
      )
      setFilterPlotFn("plot-token")
      invisible("plot-token")
    },
    envir = server_env
  )
  assign(
    "saveMetabIntensityFilterState",
    function(
      stateManager,
      filteredS4,
      configObject,
      intensityCutoffPercentile,
      proportionBelowCutoff,
      ...
    ) {
      captured$save_state_helper_call <- list(
        state_manager = stateManager,
        filtered_s4 = filteredS4,
        config_object = configObject,
        intensity_cutoff_percentile = intensityCutoffPercentile,
        proportion_below_cutoff = proportionBelowCutoff
      )
      invisible("saved-metab-state")
    },
    envir = server_env
  )
  assign(
    "reportMetabIntensityFilterSuccess",
    function(
      stats,
      intensityCutoffPercentile,
      proportionBelowCutoff,
      savedStateName,
      output,
      ...
    ) {
      captured$success_report_call <- list(
        stats = stats,
        intensity_cutoff_percentile = intensityCutoffPercentile,
        proportion_below_cutoff = proportionBelowCutoff,
        state_name = savedStateName,
        output = output
      )
      output$filter_results <- "success-report-token"
      invisible("success-report-token")
    },
    envir = server_env
  )
  assign(
    "buildMetabIntensityFilterStats",
    function(currentS4, filteredS4) {
      captured$stats_call <- list(
        current_s4 = currentS4,
        filtered_s4 = filteredS4
      )
      expected_stats
    },
    envir = server_env
  )

  scoped_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          apply_filter = TRUE,
          revert_filter = FALSE,
          intensity_cutoff_percentile = 10,
          proportion_below_cutoff = 0.5
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
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) substitute(expr),
    req = function(...) list(...)[[1]],
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(message) invisible(NULL),
    log_error = function(message) {
      captured$error_log_message <- message
      invisible(NULL)
    },
    log_warn = function(message) {
      captured$warn_log_message <- message
      invisible(NULL)
    },
    .env = asNamespace("logger")
  )

  workflow_data <- list(
    state_manager = list(
      getState = function() current_s4,
      saveState = function(...) {
        captured$save_state_call <- list(...)
        invisible(NULL)
      },
      getHistory = function() "baseline",
      revertToState = function(state_name) invisible(state_name)
    ),
    config_list = list(config = "value")
  )

  mod_metab_qc_intensity_server(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "intensity")
  expect_identical(captured$filter_call$the_object, current_s4)
  expect_identical(captured$filter_call$intensity_cutoff_percentile, 10)
  expect_identical(captured$filter_call$proportion_below_cutoff, 0.5)
  expect_identical(captured$stats_call$current_s4, current_s4)
  expect_identical(captured$stats_call$filtered_s4, filtered_s4)
  expect_identical(captured$success_report_call$stats, expected_stats)
  expect_identical(captured$success_report_call$intensity_cutoff_percentile, 10)
  expect_identical(captured$success_report_call$proportion_below_cutoff, 0.5)
  expect_identical(captured$success_report_call$state_name, "saved-metab-state")
  expect_identical(captured$success_report_call$output, captured$output)
  expect_identical(captured$qc_plot_helper_call$filtered_s4, filtered_s4)
  expect_identical(captured$qc_plot_helper_call$omic_type, "metabolomics")
  expect_type(captured$qc_plot_helper_call$set_filter_plot_fn, "closure")
  expect_identical(captured$save_state_helper_call$state_manager, workflow_data$state_manager)
  expect_identical(captured$save_state_helper_call$filtered_s4, filtered_s4)
  expect_identical(captured$save_state_helper_call$config_object, workflow_data$config_list)
  expect_identical(captured$save_state_helper_call$intensity_cutoff_percentile, 10)
  expect_identical(captured$save_state_helper_call$proportion_below_cutoff, 0.5)
  expect_identical(length(captured$reactive_values), 2L)
  expect_identical(captured$reactive_values[[1L]]$set_count, 1L)
  expect_identical(captured$reactive_values[[1L]]$last_value, "plot-token")
  expect_identical(captured$reactive_values[[2L]]$set_count, 1L)
  expect_identical(captured$reactive_values[[2L]]$last_value, expected_stats)
  expect_identical(captured$output$filter_results, "success-report-token")
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
  expect_false(exists("warn_log_message", envir = captured, inherits = FALSE))
  expect_identical(
    captured$notifications[[1L]],
    list(
      message = "Applying metabolite intensity filter...",
      id = "metab_intensity_filter_working",
      duration = NULL
    )
  )
  expect_identical(length(captured$notifications), 1L)
  expect_false(exists("removed_notification", envir = captured, inherits = FALSE))
})

test_that("metabolomics intensity module revert observer delegates reset and reporting through the seam", {
  skip_if_not(
    metabIntensityHasExtractedSeams(),
    "Extracted metabolomics intensity seams are target-only."
  )

  captured <- new.env(parent = emptyenv())
  captured$reactive_values <- list()
  server_env <- environment(mod_metab_qc_intensity_server)
  binding_name <- "reportMetabIntensityFilterRevertSuccess"
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
    function(prevStateName, output, setFilterStatsFn, setFilterPlotFn, ...) {
      captured$revert_report_call <- list(
        prev_state_name = prevStateName,
        output = output,
        set_filter_stats_fn = setFilterStatsFn,
        set_filter_plot_fn = setFilterPlotFn
      )
      output$filter_results <- "revert-report-token"
      setFilterStatsFn(NULL)
      setFilterPlotFn(NULL)
      invisible("revert-report-token")
    },
    envir = server_env
  )

  scoped_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          apply_filter = FALSE,
          revert_filter = TRUE,
          intensity_cutoff_percentile = 10,
          proportion_below_cutoff = 0.5
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
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) substitute(expr),
    req = function(...) list(...)[[1]],
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )
  scoped_mocked_bindings(
    log_info = function(message) {
      captured$log_message <- message
      invisible(NULL)
    },
    log_error = function(message) {
      captured$error_log_message <- message
      invisible(NULL)
    },
    .env = asNamespace("logger")
  )

  workflow_data <- list(
    state_manager = list(
      getState = function() NULL,
      saveState = function(...) invisible(NULL),
      getHistory = function() c("baseline_metab_state", "filtered_metab_state"),
      revertToState = function(state_name) {
        captured$reverted_state_name <- state_name
        invisible(state_name)
      }
    ),
    config_list = list(config = "value")
  )

  mod_metab_qc_intensity_server(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "intensity")
  expect_identical(captured$reverted_state_name, "baseline_metab_state")
  expect_identical(
    captured$revert_report_call$prev_state_name,
    "baseline_metab_state"
  )
  expect_identical(captured$revert_report_call$output, captured$output)
  expect_type(captured$revert_report_call$set_filter_stats_fn, "closure")
  expect_type(captured$revert_report_call$set_filter_plot_fn, "closure")
  expect_identical(length(captured$reactive_values), 2L)
  expect_identical(captured$reactive_values[[1L]]$set_count, 1L)
  expect_null(captured$reactive_values[[1L]]$last_value)
  expect_identical(captured$reactive_values[[2L]]$set_count, 1L)
  expect_null(captured$reactive_values[[2L]]$last_value)
  expect_identical(captured$output$filter_results, "revert-report-token")
  expect_false(exists("log_message", envir = captured, inherits = FALSE))
  expect_false(exists("error_log_message", envir = captured, inherits = FALSE))
  expect_false(exists("notifications", envir = captured, inherits = FALSE))
  expect_false(exists("removed_notification", envir = captured, inherits = FALSE))
})
