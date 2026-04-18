library(testthat)
suppressPackageStartupMessages(library(shiny))

source(test_path("..", "..", "R", "mod_lipid_qc_itsd.R"), local = environment())

renderUiToHtml <- function(expr) {
    value <- eval.parent(substitute(expr))
    htmltools::renderTags(value)$html
}

renderUiValue <- function(expr) {
    eval.parent(substitute(expr))
}

renderPlotValue <- function(expr) {
    eval.parent(substitute(expr))
}

test_that("registerLipidItsdSummaryOutput keeps the empty-state render shell stable", {
    output_env <- new.env(parent = emptyenv())
    metrics_calls <- 0L

    registerLipidItsdSummaryOutput(
        output = output_env,
        isMetricsVal = function() {
            metrics_calls <<- metrics_calls + 1L
            NULL
        },
        renderUiFn = renderUiToHtml
    )

    expect_identical(metrics_calls, 1L)
    expect_match(output_env$is_summary, "Click 'Analyze' to detect internal standards.", fixed = TRUE)
    expect_match(output_env$is_summary, "color: #666;", fixed = TRUE)
})

test_that("registerLipidItsdSummaryOutput keeps assay summary rendering stable", {
    output_env <- new.env(parent = emptyenv())
    metrics <- data.frame(
        assay = c("AssayA", "AssayA", "AssayB", "AssayB", "AssayC"),
        is_id = c("IS1", "IS2", "IS3", "IS4", "IS5"),
        cv = c(10, 12, 18, 22, 32),
        stringsAsFactors = FALSE
    )

    registerLipidItsdSummaryOutput(
        output = output_env,
        isMetricsVal = function() metrics,
        renderUiFn = renderUiToHtml
    )

    expect_match(output_env$is_summary, "AssayA: 2 IS (median CV: 11.0%)", fixed = TRUE)
    expect_match(output_env$is_summary, "AssayB: 2 IS (median CV: 20.0%)", fixed = TRUE)
    expect_match(output_env$is_summary, "AssayC: 1 IS (median CV: 32.0%)", fixed = TRUE)
    expect_match(output_env$is_summary, "color: green", fixed = TRUE)
    expect_match(output_env$is_summary, "color: orange", fixed = TRUE)
    expect_match(output_env$is_summary, "color: red", fixed = TRUE)
    expect_match(output_env$is_summary, "Acceptable", fixed = TRUE)
    expect_match(output_env$is_summary, "Review", fixed = TRUE)
})

test_that("registerLipidItsdVisualizationTabsOutput keeps the empty-state render shell stable", {
    output_env <- new.env(parent = emptyenv())
    metrics_calls <- 0L

    registerLipidItsdVisualizationTabsOutput(
        output = output_env,
        ns = function(id) paste("itsd", id, sep = "-"),
        isMetricsVal = function() {
            metrics_calls <<- metrics_calls + 1L
            NULL
        },
        renderUiFn = renderUiValue
    )

    expect_identical(metrics_calls, 1L)
    expect_null(output_env$is_viz_tabs)
})

test_that("registerLipidItsdVisualizationTabsOutput keeps visualization tabs stable", {
    output_env <- new.env(parent = emptyenv())
    metrics <- data.frame(
        assay = c("AssayA", "AssayB"),
        is_id = c("IS1", "IS2"),
        cv = c(10, 20),
        stringsAsFactors = FALSE
    )

    registerLipidItsdVisualizationTabsOutput(
        output = output_env,
        ns = function(id) paste("itsd", id, sep = "-"),
        isMetricsVal = function() metrics,
        renderUiFn = renderUiToHtml
    )

    expect_match(output_env$is_viz_tabs, "CV Distribution", fixed = TRUE)
    expect_match(output_env$is_viz_tabs, "Intensity Trends", fixed = TRUE)
    expect_match(output_env$is_viz_tabs, "itsd-is_viz_tabset", fixed = TRUE)
    expect_match(output_env$is_viz_tabs, "itsd-cv_plot", fixed = TRUE)
    expect_match(output_env$is_viz_tabs, "itsd-intensity_plot", fixed = TRUE)
})

test_that("registerLipidItsdCvPlotOutput keeps the CV plot shell stable", {
    output_env <- new.env(parent = emptyenv())
    metrics <- data.frame(
        assay = c("AssayA", "AssayA", "AssayB"),
        is_id = c("IS-high", "IS-low", "IS-mid"),
        cv = c(35, 10, 22),
        stringsAsFactors = FALSE
    )

    registerLipidItsdCvPlotOutput(
        output = output_env,
        isMetricsVal = function() metrics,
        renderPlotFn = renderPlotValue
    )

    expect_s3_class(output_env$cv_plot, "ggplot")
    expect_identical(
        levels(output_env$cv_plot$data$is_id),
        c("IS-low", "IS-mid", "IS-high")
    )
    expect_identical(
        levels(output_env$cv_plot$data$cv_status),
        c("Good (<15%)", "Acceptable (15-30%)", "Review (>30%)")
    )
    expect_identical(
        as.character(output_env$cv_plot$data$cv_status),
        c("Review (>30%)", "Good (<15%)", "Acceptable (15-30%)")
    )
    expect_identical(output_env$cv_plot$labels$title, "Internal Standard CV")
    expect_identical(output_env$cv_plot$labels$subtitle, "Dashed lines: 15% (green) and 30% (orange) thresholds")
    expect_identical(output_env$cv_plot$labels$y, "CV (%)")
    expect_length(output_env$cv_plot$layers, 4L)
    expect_s3_class(output_env$cv_plot$facet, "FacetWrap")
    expect_s3_class(output_env$cv_plot$coordinates, "CoordFlip")
})

loadLipidQcItsdModuleHarness <- function() {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc_itsd.R"),
        warn = FALSE
    )
    source_lines <- sub(
        "shiny::moduleServer",
        "moduleServer",
        source_lines,
        fixed = TRUE
    )

    module_env <- new.env(parent = globalenv())
    module_env$moduleServer <- function(id, module, session = shiny::getDefaultReactiveDomain()) {
        assign("capturedModule", module, envir = module_env)
        invisible(NULL)
    }

    eval(parse(text = source_lines), envir = module_env)

    module_env$itsdSummaryOutputCalls <- list()
    module_env$itsdVisualizationTabsOutputCalls <- list()
    module_env$itsdCvPlotOutputCalls <- list()
    module_env$registerLipidItsdSummaryOutput <- function(output, isMetricsVal, renderUiFn = shiny::renderUI) {
        module_env$itsdSummaryOutputCalls <- c(
            module_env$itsdSummaryOutputCalls,
            list(list(
                output = output,
                isMetricsVal = isMetricsVal,
                renderUiFn = renderUiFn
            ))
        )
        output
    }
    module_env$registerLipidItsdVisualizationTabsOutput <- function(output, ns, isMetricsVal, renderUiFn = shiny::renderUI) {
        module_env$itsdVisualizationTabsOutputCalls <- c(
            module_env$itsdVisualizationTabsOutputCalls,
            list(list(
                output = output,
                ns = ns,
                isMetricsVal = isMetricsVal,
                renderUiFn = renderUiFn
            ))
        )
        output
    }
    module_env$registerLipidItsdCvPlotOutput <- function(output, isMetricsVal, renderPlotFn = shiny::renderPlot) {
        module_env$itsdCvPlotOutputCalls <- c(
            module_env$itsdCvPlotOutputCalls,
            list(list(
                output = output,
                isMetricsVal = isMetricsVal,
                renderPlotFn = renderPlotFn
            ))
        )
        output
    }

    workflow_data <- list(
        state_manager = list(
            getState = function(...) NULL
        )
    )

    module_env$mod_lipid_qc_itsd_server(
        "itsd",
        workflow_data = workflow_data,
        omic_type = "lipidomics",
        experiment_label = "Lipidomics"
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "itsd", session = session)
    session$flushReact()

    list(moduleEnv = module_env, session = session)
}

test_that("mod_lipid_qc_itsd_server delegates summary output registration through the seam", {
    harness <- loadLipidQcItsdModuleHarness()

    expect_length(harness$moduleEnv$itsdSummaryOutputCalls, 1L)
    expect_true(is.function(harness$moduleEnv$itsdSummaryOutputCalls[[1]]$isMetricsVal))
    expect_true(is.function(harness$moduleEnv$itsdSummaryOutputCalls[[1]]$renderUiFn))
})

test_that("mod_lipid_qc_itsd_server delegates visualization-tab registration through the seam", {
    harness <- loadLipidQcItsdModuleHarness()

    expect_length(harness$moduleEnv$itsdVisualizationTabsOutputCalls, 1L)
    expect_true(is.function(harness$moduleEnv$itsdVisualizationTabsOutputCalls[[1]]$ns))
    expect_true(is.function(harness$moduleEnv$itsdVisualizationTabsOutputCalls[[1]]$isMetricsVal))
    expect_true(is.function(harness$moduleEnv$itsdVisualizationTabsOutputCalls[[1]]$renderUiFn))
})

test_that("mod_lipid_qc_itsd_server delegates CV plot registration through the seam", {
    harness <- loadLipidQcItsdModuleHarness()

    expect_length(harness$moduleEnv$itsdCvPlotOutputCalls, 1L)
    expect_true(is.function(harness$moduleEnv$itsdCvPlotOutputCalls[[1]]$isMetricsVal))
    expect_true(is.function(harness$moduleEnv$itsdCvPlotOutputCalls[[1]]$renderPlotFn))
})
