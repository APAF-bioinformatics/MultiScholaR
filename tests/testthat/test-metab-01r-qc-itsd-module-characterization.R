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
    file.path(repo_root, "R", "mod_metab_qc_itsd.R"),
    file.path(repo_root, "R", "mod_metab_qc_itsd_ui.R"),
    file.path(repo_root, "R", "mod_metab_qc_itsd_server.R")
  ),
  symbols = c(
    "analyzeMetabQcItsdData",
    "buildMetabQcItsdSummaryUi",
    "buildMetabQcItsdVizTabsUi",
    "buildMetabQcItsdCvPlot",
    "buildMetabQcItsdIntensityPlot",
    "runMetabQcItsdServerBody",
    "mod_metab_qc_itsd_ui",
    "mod_metab_qc_itsd_server"
  ),
  env = environment()
)

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      internal_standard_regex = "character",
      metabolite_id_column = "character",
      annotation_id_column = "character",
      sample_id = "character",
      metabolite_data = "list"
    ),
    prototype = list(
      internal_standard_regex = NA_character_,
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id",
      sample_id = "sample_id",
      metabolite_data = list()
    )
  )
}

test_that("metabolomics ITSD analysis seam preserves pattern fallback, column search order, and long-data assembly", {
  capture <- new.env(parent = emptyenv())
  capture$metric_calls <- list()
  capture$logs <- character()

  current_s4 <- methods::new(
    "MetaboliteAssayData",
    internal_standard_regex = "^IS_",
    metabolite_id_column = "feature_id",
    annotation_id_column = "annotation_id",
    sample_id = "sample_id",
    design_matrix = data.frame(
      sample_id = c("Sample1", "Sample2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    metabolite_data = list(
      Plasma = data.frame(
        feature_id = c("IS_A", "M1"),
        annotation_id = c("ann_a", "ann_b"),
        Sample1 = c(10, 100),
        Sample2 = c(12, 120),
        stringsAsFactors = FALSE
      ),
      Serum = data.frame(
        feature_id = c("x1", "x2"),
        annotation_id = c("ann_c", "ann_d"),
        Name = c("IS_B", "M2"),
        Sample1 = c(30, 300),
        Sample2 = c(36, 360),
        stringsAsFactors = FALSE
      )
    )
  )

  visible <- withVisible(
    analyzeMetabQcItsdData(
      currentS4 = current_s4,
      inputPattern = "",
      getInternalStandardMetricsFn = function(assay_data, is_pattern, metabolite_id_col, sample_id_col) {
        capture$metric_calls[[length(capture$metric_calls) + 1L]] <- list(
          metabolite_id_col = metabolite_id_col,
          sample_id_col = sample_id_col,
          assay_cols = colnames(assay_data),
          pattern = is_pattern
        )

        matched <- grepl("^IS_", assay_data[[metabolite_id_col]])
        if (!any(matched)) {
          return(data.frame(
            is_id = character(),
            mean_intensity = numeric(),
            cv = numeric()
          ))
        }

        data.frame(
          is_id = assay_data[[metabolite_id_col]][matched],
          mean_intensity = rep(100, sum(matched)),
          cv = rep(if (identical(metabolite_id_col, "Name")) 22 else 12, sum(matched))
        )
      },
      getMetaboliteQuantDataFn = function(assay_data) {
        list(sample_names = c("Sample1", "Sample2"))
      },
      pivotLongerFn = function(data, cols, names_to, values_to) {
        id_col <- setdiff(names(data), cols)
        out <- data.frame(
          value_id = rep(data[[id_col]], each = length(cols)),
          sample_name = rep(cols, times = nrow(data)),
          intensity_value = as.vector(t(as.matrix(data[cols]))),
          stringsAsFactors = FALSE
        )
        names(out) <- c(id_col, names_to, values_to)
        out
      },
      allOfFn = function(values) values,
      logInfoFn = function(message) {
        capture$logs <- c(capture$logs, message)
        invisible(NULL)
      }
    )
  )

  expect_true(visible$visible)
  expect_identical(visible$value$pattern, "^IS_")
  expect_identical(visible$value$nIsTotal, 2L)
  expect_identical(
    vapply(capture$metric_calls, `[[`, character(1), "metabolite_id_col"),
    c("feature_id", "feature_id", "annotation_id", "Name")
  )
  expect_true(all(vapply(capture$metric_calls, `[[`, character(1), "sample_id_col") == "sample_id"))
  expect_identical(
    capture$logs,
    c(
      "Analyzing internal standards with pattern: ^IS_",
      "Found IS in column: feature_id",
      "Found IS in column: Name"
    )
  )
  expect_identical(unique(visible$value$metrics$assay), c("Plasma", "Serum"))
  expect_identical(
    unique(visible$value$longData$IS_ID),
    c("IS_A", "IS_B")
  )
  expect_match(visible$value$resultText, "Pattern used: \\^IS_")
  expect_match(visible$value$resultText, "Plasma: 1 internal standards", fixed = TRUE)
  expect_match(visible$value$resultText, "Serum: 1 internal standards", fixed = TRUE)
})

test_that("metabolomics ITSD summary UI seam assembles placeholder and assay status rows", {
  placeholder <- buildMetabQcItsdSummaryUi(
    metrics = NULL,
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
  expect_match(placeholder$children[[2]], "Click 'Analyze'", fixed = TRUE)

  summary_ui <- buildMetabQcItsdSummaryUi(
    metrics = data.frame(
      assay = c("Plasma", "Serum", "Urine"),
      cv = c(10, 22, 35),
      stringsAsFactors = FALSE
    ),
    iconFn = function(name, style = NULL) {
      list(type = "icon", name = name, style = style)
    },
    listTagFn = function(items, style = NULL) {
      list(type = "ul", items = items, style = style)
    },
    itemTagFn = function(...) {
      list(type = "li", children = list(...))
    },
    horizontalRuleFn = function() {
      list(type = "hr")
    },
    smallTagFn = function(..., style = NULL) {
      list(type = "small", children = list(...), style = style)
    },
    tagListFn = function(...) {
      list(type = "tagList", children = list(...))
    }
  )

  expect_identical(summary_ui$children[[1]]$type, "ul")
  expect_identical(summary_ui$children[[1]]$style, "list-style: none; padding-left: 0;")
  expect_length(summary_ui$children[[1]]$items, 3L)
  expect_identical(summary_ui$children[[1]]$items[[1]]$children[[1]]$name, "check-circle")
  expect_identical(summary_ui$children[[1]]$items[[2]]$children[[1]]$name, "exclamation-circle")
  expect_identical(summary_ui$children[[1]]$items[[3]]$children[[1]]$name, "times-circle")
  expect_identical(summary_ui$children[[1]]$items[[1]]$children[[1]]$style, "color: green")
  expect_identical(summary_ui$children[[1]]$items[[2]]$children[[1]]$style, "color: orange")
  expect_identical(summary_ui$children[[1]]$items[[3]]$children[[1]]$style, "color: red")
  expect_identical(summary_ui$children[[1]]$items[[1]]$children[[2]], " Plasma: 1 IS (median CV: 10.0%)")
  expect_identical(summary_ui$children[[1]]$items[[2]]$children[[2]], " Serum: 1 IS (median CV: 22.0%)")
  expect_identical(summary_ui$children[[1]]$items[[3]]$children[[2]], " Urine: 1 IS (median CV: 35.0%)")
  expect_identical(summary_ui$children[[2]]$type, "hr")
  expect_identical(summary_ui$children[[3]]$type, "small")
  expect_identical(summary_ui$children[[3]]$children[[1]]$name, "info-circle")
  expect_match(summary_ui$children[[3]]$children[[2]], "CV < 15%: Good", fixed = TRUE)
})

test_that("metabolomics ITSD visualization tabs seam returns null without metrics and builds namespaced plot tabs", {
  expect_null(
    buildMetabQcItsdVizTabsUi(
      metrics = NULL,
      ns = function(value) value
    )
  )

  tabs_ui <- buildMetabQcItsdVizTabsUi(
    metrics = data.frame(
      assay = "Plasma",
      cv = 12,
      stringsAsFactors = FALSE
    ),
    ns = function(value) paste0("ns-", value),
    tabPanelFn = function(title, ...) {
      list(type = "tab", title = title, children = list(...))
    },
    brFn = function() {
      list(type = "br")
    },
    jquiResizableFn = function(content) {
      list(type = "resizable", child = content)
    },
    plotOutputFn = function(outputId, height = NULL) {
      list(type = "plot", outputId = outputId, height = height)
    },
    tabsetPanelFn = function(id, ...) {
      list(type = "tabset", id = id, tabs = list(...))
    }
  )

  expect_identical(tabs_ui$type, "tabset")
  expect_identical(tabs_ui$id, "ns-is_viz_tabset")
  expect_identical(
    vapply(tabs_ui$tabs, `[[`, character(1), "title"),
    c("CV Distribution", "Intensity Trends")
  )
  expect_identical(tabs_ui$tabs[[1L]]$children[[2L]]$type, "resizable")
  expect_identical(tabs_ui$tabs[[1L]]$children[[2L]]$child$outputId, "ns-cv_plot")
  expect_identical(tabs_ui$tabs[[1L]]$children[[2L]]$child$height, "500px")
  expect_identical(tabs_ui$tabs[[2L]]$children[[2L]]$child$outputId, "ns-intensity_plot")
  expect_identical(tabs_ui$tabs[[2L]]$children[[2L]]$child$height, "500px")
})

test_that("metabolomics ITSD CV plot seam preserves reorder, status buckets, and lollipop plot metadata", {
  plot <- buildMetabQcItsdCvPlot(
    metrics = data.frame(
      is_id = c("IS_high", "IS_low", "IS_mid"),
      cv = c(35, 10, 22),
      assay = c("Plasma", "Plasma", "Serum"),
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(
    levels(plot$data$is_id),
    c("IS_low", "IS_mid", "IS_high")
  )
  expect_identical(
    as.character(plot$data$cv_status),
    c("Review (>30%)", "Good (<15%)", "Acceptable (15-30%)")
  )
  expect_identical(
    levels(plot$data$cv_status),
    c("Good (<15%)", "Acceptable (15-30%)", "Review (>30%)")
  )
  expect_identical(plot$labels$title, "Internal Standard CV")
  expect_identical(plot$labels$subtitle, "Dashed lines: 15% (green) and 30% (orange) thresholds")
  expect_identical(plot$labels$y, "CV (%)")
  expect_identical(plot$scales$get_scales("colour")$name, "CV Status")
  expect_length(plot$layers, 4L)
  expect_true(inherits(plot$facet, "FacetWrap"))
  expect_true(inherits(plot$coordinates, "CoordFlip"))
})

test_that("metabolomics ITSD intensity plot seam preserves labels, layers, and theme metadata", {
  plot <- buildMetabQcItsdIntensityPlot(
    longData = data.frame(
      Sample = c("S1", "S2", "S1", "S2"),
      Intensity = c(100, 120, 80, 85),
      IS_ID = c("IS_A", "IS_A", "IS_B", "IS_B"),
      assay = c("Plasma", "Plasma", "Serum", "Serum"),
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$title, "Internal Standard Intensity Across Samples")
  expect_identical(plot$labels$x, "Sample")
  expect_identical(plot$labels$y, "Intensity")
  expect_identical(plot$labels$colour, "Internal Standard")
  expect_length(plot$layers, 2L)
  expect_true(inherits(plot$facet, "FacetWrap"))
  expect_identical(plot$theme$legend.position, "bottom")
  expect_identical(plot$theme$axis.text.x$angle, 90)
  expect_identical(plot$theme$axis.text.x$hjust, 1)
  expect_identical(plot$theme$axis.text.x$vjust, 0.5)
  expect_identical(plot$theme$axis.text.x$size, 6)
})

test_that("metabolomics ITSD server body seam preserves helper handoffs and reactive wiring", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$reactive_updates <- list(metrics = list(), long_data = list())
  output <- new.env(parent = emptyenv())
  reactive_names <- c("metrics", "long_data")
  reactive_initial_values <- list(NULL, NULL)
  reactive_index <- 0L

  visible <- withVisible(
    runMetabQcItsdServerBody(
      input = list(analyze_is = TRUE, is_pattern = "manual-regex"),
      output = output,
      session = list(ns = function(value) paste("itsd", value, sep = "-")),
      workflowData = list(
        state_manager = list(
          getState = function() "current-s4"
        )
      ),
      omicType = "metabolomics",
      experimentLabel = "Experiment A",
      reactiveValFn = function(value = NULL) {
        reactive_index <<- reactive_index + 1L
        name <- reactive_names[[reactive_index]]
        stored <- reactive_initial_values[[reactive_index]]

        function(new_value) {
          if (missing(new_value)) {
            stored
          } else {
            stored <<- new_value
            updates <- captured$reactive_updates[[name]]
            updates[[length(updates) + 1L]] <- list(value = new_value)
            captured$reactive_updates[[name]] <- updates
            invisible(NULL)
          }
        }
      },
      observeEventFn = function(eventExpr, handlerExpr, ...) {
        if (isTRUE(eval(substitute(eventExpr), parent.frame()))) {
          eval(substitute(handlerExpr), parent.frame())
        }
        invisible(NULL)
      },
      reqFn = function(...) {
        args <- list(...)
        if (length(args) == 0) {
          return(invisible(NULL))
        }
        args[[1L]]
      },
      showNotificationFn = function(ui, id = NULL, duration = NULL, type = NULL, ...) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          ui = ui,
          id = id,
          duration = duration,
          type = type
        )
        invisible(NULL)
      },
      removeNotificationFn = function(id) {
        captured$removed_notification <- id
        invisible(NULL)
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      renderUiFn = function(expr) {
        list(expr = substitute(expr), env = parent.frame())
      },
      renderPlotFn = function(expr) {
        list(expr = substitute(expr), env = parent.frame())
      },
      analyzeMetabQcItsdDataFn = function(currentS4, inputPattern, ...) {
        captured$analysis_call <- list(
          current_s4 = currentS4,
          input_pattern = inputPattern
        )
        list(
          metrics = "metrics-token",
          longData = "long-data-token",
          resultText = "summary-text",
          nIsTotal = 2L
        )
      },
      buildMetabQcItsdSummaryUiFn = function(metrics, ...) {
        captured$summary_call <- list(metrics = metrics)
        "summary-ui-token"
      },
      buildMetabQcItsdVizTabsUiFn = function(metrics, ns, ...) {
        captured$viz_call <- list(
          metrics = metrics,
          tabset_id = ns("is_viz_tabset")
        )
        "viz-ui-token"
      },
      buildMetabQcItsdCvPlotFn = function(metrics, ...) {
        captured$cv_plot_call <- list(metrics = metrics)
        "cv-plot-token"
      },
      buildMetabQcItsdIntensityPlotFn = function(longData, ...) {
        captured$intensity_plot_call <- list(long_data = longData)
        "intensity-plot-token"
      },
      logInfoFn = function(message) {
        captured$logs <- c(captured$logs, message)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$errors <- c(captured$errors, message)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_null(visible$value)
  expect_identical(
    captured$analysis_call,
    list(
      current_s4 = "current-s4",
      input_pattern = "manual-regex"
    )
  )
  expect_identical(captured$reactive_updates$metrics[[1L]]$value, "metrics-token")
  expect_identical(captured$reactive_updates$long_data[[1L]]$value, "long-data-token")
  expect_identical(output$is_results, "summary-text")
  expect_identical(captured$removed_notification, "is_analysis_working")
  expect_identical(captured$logs, "IS analysis complete: 2 standards found")
  expect_length(captured$notifications, 2L)
  expect_identical(captured$notifications[[1L]]$ui, "Analyzing internal standards...")
  expect_identical(captured$notifications[[2L]]$ui, "Found 2 internal standards")

  expect_identical(eval(output$is_summary$expr, output$is_summary$env), "summary-ui-token")
  expect_identical(captured$summary_call$metrics, "metrics-token")

  expect_identical(eval(output$is_viz_tabs$expr, output$is_viz_tabs$env), "viz-ui-token")
  expect_identical(captured$viz_call$metrics, "metrics-token")
  expect_identical(captured$viz_call$tabset_id, "itsd-is_viz_tabset")

  expect_identical(eval(output$cv_plot$expr, output$cv_plot$env), "cv-plot-token")
  expect_identical(captured$cv_plot_call$metrics, "metrics-token")

  expect_identical(eval(output$intensity_plot$expr, output$intensity_plot$env), "intensity-plot-token")
  expect_identical(captured$intensity_plot_call$long_data, "long-data-token")
})

test_that("metabolomics ITSD module server delegates the analyze observer through the seam and clears stale long data", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$reactive_updates <- list(metrics = list(), long_data = list())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_name <- "analyzeMetabQcItsdData"
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
    function(currentS4, inputPattern, ...) {
      captured$helper_call <- list(
        current_s4 = currentS4,
        input_pattern = inputPattern
      )
      list(
        metrics = "metrics-token",
        longData = NULL,
        resultText = "summary-text",
        nIsTotal = 2L
      )
    },
    envir = server_env
  )

  reactive_names <- c("metrics", "long_data")
  reactive_initial_values <- list(NULL, "stale-long-data")
  reactive_index <- 0L

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(analyze_is = TRUE, is_pattern = "manual-regex"),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      name <- reactive_names[[reactive_index]]
      stored <- reactive_initial_values[[reactive_index]]

      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          updates <- captured$reactive_updates[[name]]
          updates[[length(updates) + 1L]] <- list(value = new_value)
          captured$reactive_updates[[name]] <- updates
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
    req = function(...) {
      args <- list(...)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotification = function(ui, id = NULL, duration = NULL, type = NULL, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        ui = ui,
        id = id,
        duration = duration,
        type = type
      )
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured$removed_notification <- id
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = list(
      state_manager = list(
        getState = function() "current-s4"
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_identical(
    captured$helper_call,
    list(
      current_s4 = "current-s4",
      input_pattern = "manual-regex"
    )
  )
  expect_identical(captured$reactive_updates$metrics[[1L]]$value, "metrics-token")
  expect_null(captured$reactive_updates$long_data[[1L]]$value)
  expect_identical(output$is_results, "summary-text")
  expect_identical(captured$removed_notification, "is_analysis_working")
  expect_length(captured$notifications, 2L)
  expect_identical(captured$notifications[[1L]]$ui, "Analyzing internal standards...")
  expect_identical(captured$notifications[[1L]]$id, "is_analysis_working")
  expect_identical(captured$notifications[[2L]]$ui, "Found 2 internal standards")
  expect_identical(captured$notifications[[2L]]$type, "message")
})

test_that("metabolomics ITSD module server delegates through the server body seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_name <- "runMetabQcItsdServerBody"
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
        list(analyze_is = FALSE, is_pattern = "manual-regex"),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
    },
    .package = "shiny"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = expected_workflow_data,
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_identical(captured$server_body_call$input$is_pattern, "manual-regex")
  expect_identical(captured$server_body_call$output, output)
  expect_identical(captured$server_body_call$session$ns("cv_plot"), "itsd-cv_plot")
  expect_identical(captured$server_body_call$workflow_data, expected_workflow_data)
  expect_identical(captured$server_body_call$omic_type, "metabolomics")
  expect_identical(captured$server_body_call$experiment_label, "Experiment A")
})

test_that("metabolomics ITSD module server delegates the summary render through the seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_names <- c("analyzeMetabQcItsdData", "buildMetabQcItsdSummaryUi")
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

  expected_metrics <- data.frame(
    assay = "Plasma",
    cv = 12,
    is_id = "IS_A",
    stringsAsFactors = FALSE
  )

  assign(
    "analyzeMetabQcItsdData",
    function(...) {
      stop("analyzeMetabQcItsdData should not be called directly")
    },
    envir = server_env
  )
  assign(
    "buildMetabQcItsdSummaryUi",
    function(metrics, ...) {
      captured$summary_call <- list(metrics = metrics)
      "summary-ui-token"
    },
    envir = server_env
  )

  reactive_initial_values <- list(expected_metrics, NULL)
  reactive_index <- 0L

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(analyze_is = FALSE, is_pattern = ""),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      stored <- reactive_initial_values[[reactive_index]]

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
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) {
      list(expr = substitute(expr), env = parent.frame())
    },
    renderPlot = function(expr) substitute(expr),
    req = function(...) {
      args <- list(...)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = list(
      state_manager = list(
        getState = function() "current-s4"
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_match(
    paste(deparse(output$is_summary$expr), collapse = " "),
    "buildMetabQcItsdSummaryUi",
    fixed = TRUE
  )
  expect_identical(eval(output$is_summary$expr, output$is_summary$env), "summary-ui-token")
  expect_identical(captured$summary_call$metrics, expected_metrics)
})

test_that("metabolomics ITSD module server delegates the visualization tabs render through the seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_name <- "buildMetabQcItsdVizTabsUi"
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

  expected_metrics <- data.frame(
    assay = "Plasma",
    cv = 12,
    is_id = "IS_A",
    stringsAsFactors = FALSE
  )

  assign(
    binding_name,
    function(metrics, ns, ...) {
      captured$viz_call <- list(
        metrics = metrics,
        tabset_id = ns("is_viz_tabset")
      )
      "viz-tabs-token"
    },
    envir = server_env
  )

  reactive_initial_values <- list(expected_metrics, NULL)
  reactive_index <- 0L

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(analyze_is = FALSE, is_pattern = ""),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      stored <- reactive_initial_values[[reactive_index]]

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
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) {
      list(expr = substitute(expr), env = parent.frame())
    },
    renderPlot = function(expr) substitute(expr),
    req = function(...) {
      args <- list(...)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = list(
      state_manager = list(
        getState = function() "current-s4"
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_match(
    paste(deparse(output$is_viz_tabs$expr), collapse = " "),
    "buildMetabQcItsdVizTabsUi",
    fixed = TRUE
  )
  expect_identical(eval(output$is_viz_tabs$expr, output$is_viz_tabs$env), "viz-tabs-token")
  expect_identical(captured$viz_call$metrics, expected_metrics)
  expect_identical(captured$viz_call$tabset_id, "itsd-is_viz_tabset")
})

test_that("metabolomics ITSD module server delegates the CV plot render through the seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_name <- "buildMetabQcItsdCvPlot"
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

  expected_metrics <- data.frame(
    assay = "Plasma",
    cv = 12,
    is_id = "IS_A",
    stringsAsFactors = FALSE
  )

  assign(
    binding_name,
    function(metrics, ...) {
      captured$cv_plot_call <- list(metrics = metrics)
      "cv-plot-token"
    },
    envir = server_env
  )

  reactive_initial_values <- list(expected_metrics, NULL)
  reactive_index <- 0L

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(analyze_is = FALSE, is_pattern = ""),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      stored <- reactive_initial_values[[reactive_index]]

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
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) {
      list(expr = substitute(expr), env = parent.frame())
    },
    req = function(...) {
      args <- list(...)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = list(
      state_manager = list(
        getState = function() "current-s4"
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_match(
    paste(deparse(output$cv_plot$expr), collapse = " "),
    "buildMetabQcItsdCvPlot",
    fixed = TRUE
  )
  expect_identical(eval(output$cv_plot$expr, output$cv_plot$env), "cv-plot-token")
  expect_identical(captured$cv_plot_call$metrics, expected_metrics)
})

test_that("metabolomics ITSD module server delegates the intensity plot render through the seam", {
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_qc_itsd_server)
  binding_name <- "buildMetabQcItsdIntensityPlot"
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

  expected_long_data <- data.frame(
    Sample = c("S1", "S2"),
    Intensity = c(100, 120),
    IS_ID = c("IS_A", "IS_A"),
    assay = c("Plasma", "Plasma"),
    stringsAsFactors = FALSE
  )

  assign(
    binding_name,
    function(longData, ...) {
      captured$intensity_plot_call <- list(longData = longData)
      "intensity-plot-token"
    },
    envir = server_env
  )

  reactive_initial_values <- list(NULL, expected_long_data)
  reactive_index <- 0L

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      module(
        list(analyze_is = FALSE, is_pattern = ""),
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      reactive_index <<- reactive_index + 1L
      stored <- reactive_initial_values[[reactive_index]]

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
    renderText = function(expr) substitute(expr),
    renderUI = function(expr) substitute(expr),
    renderPlot = function(expr) {
      list(expr = substitute(expr), env = parent.frame())
    },
    req = function(...) {
      args <- list(...)
      if (length(args) == 0) {
        return(invisible(NULL))
      }
      args[[1L]]
    },
    showNotification = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_qc_itsd_server(
    id = "itsd",
    workflow_data = list(
      state_manager = list(
        getState = function() "current-s4"
      )
    ),
    omic_type = "metabolomics",
    experiment_label = "Experiment A"
  )

  expect_identical(captured$module_id, "itsd")
  expect_match(
    paste(deparse(output$intensity_plot$expr), collapse = " "),
    "buildMetabQcItsdIntensityPlot",
    fixed = TRUE
  )
  expect_identical(eval(output$intensity_plot$expr, output$intensity_plot$env), "intensity-plot-token")
  expect_identical(captured$intensity_plot_call$longData, expected_long_data)
})
