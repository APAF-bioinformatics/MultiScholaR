library(testthat)
suppressPackageStartupMessages(library(shiny))

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingLipidNormTargetFiles <- function() {
  required_paths <- c(
    "R/mod_lipid_norm_support_helpers.R",
    "R/mod_lipid_norm_observer_helpers.R",
    "R/mod_lipid_norm_workflow_helpers.R",
    "R/mod_lipid_norm_runtime_helpers.R",
    "R/mod_lipid_norm_server_helpers.R",
    "R/mod_lipid_norm_ui_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only lipid norm helper file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingLipidNormTargetFiles()

source(test_path("..", "..", "R", "mod_lipid_norm_support_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_observer_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_workflow_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_runtime_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_ui_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_norm.R"), local = environment())

test_that("buildLipidNormOptionsControlPanel keeps the left control-panel contract stable", {
    rendered <- htmltools::renderTags(
        buildLipidNormOptionsControlPanel(
            ns = function(id) paste0("norm-", id)
        )
    )$html

    expect_match(rendered, "Normalization Options", fixed = TRUE)
    expect_match(rendered, "Plot Aesthetics", fixed = TRUE)
    expect_match(rendered, "RUV-III Batch Correction", fixed = TRUE)
    expect_match(rendered, "Run Normalization Pipeline", fixed = TRUE)
    expect_match(rendered, "norm-color_variable", fixed = TRUE)
    expect_match(rendered, "norm-itsd_aggregation", fixed = TRUE)
    expect_match(rendered, "norm-k_penalty_weight", fixed = TRUE)
    expect_match(rendered, "norm-export_session", fixed = TRUE)
})

test_that("buildLipidNormQcTabsetPanel keeps the right QC-tabset contract stable", {
    rendered <- htmltools::renderTags(
        buildLipidNormQcTabsetPanel(
            ns = function(id) paste0("norm-", id)
        )
    )$html

    expect_match(rendered, "ITSD Selection", fixed = TRUE)
    expect_match(rendered, "RUV QC", fixed = TRUE)
    expect_match(rendered, "Correlation Filtering", fixed = TRUE)
    expect_match(rendered, "Normalization Pipeline Log", fixed = TRUE)
    expect_match(rendered, "norm-norm_qc_tabs", fixed = TRUE)
    expect_match(rendered, "norm-itsd_selection_ui", fixed = TRUE)
    expect_match(rendered, "norm-correlation_filter_summary", fixed = TRUE)
    expect_match(rendered, "norm-final_qc_plot", fixed = TRUE)
})

test_that("mod_lipid_norm_ui delegates the left control panel through the helper seam", {
    helper_env <- environment(mod_lipid_norm_ui)
    original_helper <- get("buildLipidNormOptionsControlPanel", envir = helper_env)
    on.exit(
        assign("buildLipidNormOptionsControlPanel", original_helper, envir = helper_env),
        add = TRUE
    )

    assign(
        "buildLipidNormOptionsControlPanel",
        function(ns) shiny::tags$div("stub-control-panel", id = ns("stub_panel")),
        envir = helper_env
    )

    rendered <- htmltools::renderTags(
        mod_lipid_norm_ui("norm")
    )$html

    expect_match(rendered, "stub-control-panel", fixed = TRUE)
    expect_match(rendered, "norm-stub_panel", fixed = TRUE)
    expect_match(rendered, "ITSD Selection", fixed = TRUE)
})

test_that("mod_lipid_norm_ui delegates the right QC tabset through the helper seam", {
    helper_env <- environment(mod_lipid_norm_ui)
    original_helper <- get("buildLipidNormQcTabsetPanel", envir = helper_env)
    on.exit(
        assign("buildLipidNormQcTabsetPanel", original_helper, envir = helper_env),
        add = TRUE
    )

    assign(
        "buildLipidNormQcTabsetPanel",
        function(ns) shiny::tags$div("stub-qc-tabset", id = ns("stub_qc_tabset")),
        envir = helper_env
    )

    rendered <- htmltools::renderTags(
        mod_lipid_norm_ui("norm")
    )$html

    expect_match(rendered, "stub-qc-tabset", fixed = TRUE)
    expect_match(rendered, "norm-stub_qc_tabset", fixed = TRUE)
    expect_match(rendered, "Normalization Options", fixed = TRUE)
})

loadLipidNormModuleHarness <- function(
    triggerExportSession = FALSE,
    triggerApplyCorrelationFilter = FALSE,
    triggerSkipCorrelationFilter = FALSE,
    triggerResetNormalization = FALSE,
    triggerRunNormalization = FALSE,
    triggerSelectedTab = FALSE,
    selectedTabValue = NULL,
    stubModuleServerPublicWrapper = FALSE,
    stubModuleServerEntryShell = FALSE,
    stubModuleServerShell = FALSE,
    stubStartupRuntime = FALSE,
    stubServerRuntime = FALSE,
    stubStartupObserverRuntime = FALSE,
    stubPrimaryStartupOutputs = FALSE,
    stubItsdSelectionRuntime = FALSE,
    stubPostNormalizationOutputs = FALSE
) {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_norm.R"),
        warn = FALSE
    )
    source_lines <- sub(
        "shiny::moduleServer",
        "moduleServer",
        source_lines,
        fixed = TRUE
    )
    source_lines <- gsub(
        "shiny::observeEvent",
        "observeEvent",
        source_lines,
        fixed = TRUE
    )

    module_env <- new.env(parent = globalenv())
    module_env$moduleServer <- function(id, module, session = shiny::getDefaultReactiveDomain()) {
        assign("capturedModule", module, envir = module_env)
        invisible(NULL)
    }
    module_env$capturedObserverEvents <- character()
    module_env$observeEvent <- function(eventExpr, handlerExpr, ...) {
        event_label <- paste(deparse(substitute(eventExpr)), collapse = "")
        module_env$capturedObserverEvents <- c(
            module_env$capturedObserverEvents,
            event_label
        )

        if (isTRUE(triggerExportSession) && identical(event_label, "input$export_session")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        if (isTRUE(triggerApplyCorrelationFilter) && identical(event_label, "input$apply_correlation_filter")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        if (isTRUE(triggerSkipCorrelationFilter) && identical(event_label, "input$skip_correlation_filter")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        if (isTRUE(triggerResetNormalization) && identical(event_label, "input$reset_normalization")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        if (isTRUE(triggerRunNormalization) && identical(event_label, "input$run_normalization")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        if (isTRUE(triggerSelectedTab) && identical(event_label, "selected_tab()")) {
            shiny::isolate(eval(substitute(handlerExpr), envir = parent.frame()))
        }

        invisible(NULL)
    }

    helper_files <- c(
        "mod_lipid_norm_support_helpers.R",
        "mod_lipid_norm_observer_helpers.R",
        "mod_lipid_norm_workflow_helpers.R",
        "mod_lipid_norm_runtime_helpers.R",
        "mod_lipid_norm_server_helpers.R",
        "mod_lipid_norm_ui_helpers.R",
        "mod_lipid_norm_ui.R",
        "mod_lipid_norm_server.R"
    )
    for (helper_file in helper_files) {
        helper_lines <- readLines(
            test_path("..", "..", "R", helper_file),
            warn = FALSE
        )
        helper_lines <- sub(
            "shiny::moduleServer",
            "moduleServer",
            helper_lines,
            fixed = TRUE
        )
        helper_lines <- gsub(
            "shiny::observeEvent",
            "observeEvent",
            helper_lines,
            fixed = TRUE
        )
        eval(parse(text = helper_lines), envir = module_env)
    }

    eval(parse(text = source_lines), envir = module_env)

    module_env$staticQcImageHelperCalls <- list()
    module_env$assayLabelHelperCalls <- list()
    module_env$qcImageRendererBuilderCalls <- list()
    module_env$assayLabelBuilderCalls <- list()
    module_env$normLogHelperCalls <- list()
    module_env$normLogBuilderCalls <- list()
    module_env$itsdSelectionHelperCalls <- list()
    module_env$ruvQcHelperCalls <- list()
    module_env$correlationFilterSummaryHelperCalls <- list()
    module_env$finalQcPlotHelperCalls <- list()
    module_env$designDrivenChoiceObserverHelperCalls <- list()
    module_env$assayNameInitializationHelperCalls <- list()
    module_env$ruvCancorOutputsHelperCalls <- list()
    module_env$itsdTableHelperCalls <- list()
    module_env$itsdSelectionTrackingHelperCalls <- list()
    module_env$runNormalizationObserverHelperCalls <- list()
    module_env$resetNormalizationObserverHelperCalls <- list()
    module_env$applyCorrelationFilterObserverHelperCalls <- list()
    module_env$skipCorrelationFilterObserverHelperCalls <- list()
    module_env$exportSessionObserverHelperCalls <- list()
    module_env$selectedTabPreNormalizationObserverHelperCalls <- list()
    module_env$runNormalizationHelperCalls <- list()
    module_env$preNormalizationQcHelperCalls <- list()
    module_env$selectedTabAutoTriggerHelperCalls <- list()
    module_env$exportSessionHelperCalls <- list()
    module_env$applyCorrelationFilterHelperCalls <- list()
    module_env$skipCorrelationFilterHelperCalls <- list()
    module_env$resetNormalizationHelperCalls <- list()
    module_env$reactiveStateHelperCalls <- list()
    module_env$startupRuntimeHelperCalls <- list()
    module_env$serverRuntimeHelperCalls <- list()
    module_env$moduleServerPublicWrapperHelperCalls <- list()
    module_env$moduleServerEntryShellHelperCalls <- list()
    module_env$moduleServerShellHelperCalls <- list()
    module_env$startupObserverRuntimeHelperCalls <- list()
    module_env$primaryStartupOutputsHelperCalls <- list()
    module_env$itsdSelectionRuntimeHelperCalls <- list()
    module_env$postNormalizationOutputsHelperCalls <- list()
    module_env$addLogBuilderCalls <- list()
    module_env$getPlotAestheticsBuilderCalls <- list()
    module_env$compositeFromFilesBuilderCalls <- list()
    module_env$itsdSelectionBuilderCalls <- list()
    module_env$ruvQcUiBuilderCalls <- list()
    module_env$correlationFilterSummaryBuilderCalls <- list()
    module_env$finalQcPlotBuilderCalls <- list()
    module_env$reactiveStateStub <- new.env(parent = emptyenv())
    module_env$addLogStub <- function(message) {
        paste("add log stub", message)
    }
    module_env$getPlotAestheticsStub <- function() {
        list(
            color_var = "builder-color",
            shape_var = "builder-shape"
        )
    }
    module_env$generateCompositeFromFilesStub <- function(plot_files, ncol = 3, row_labels = NULL, column_labels = NULL) {
        list(
            plot_files = plot_files,
            ncol = ncol,
            row_labels = row_labels,
            column_labels = column_labels
        )
    }
    module_env$assayLabelRendererStub <- function(assaySlot) paste("assay label renderer stub", assaySlot)
    module_env$normLogRendererStub <- function() "norm log renderer stub"
    module_env$itsdSelectionRendererStub <- function() "itsd selection renderer stub"
    module_env$ruvQcUiRendererStub <- function() "ruv qc ui renderer stub"
    module_env$correlationFilterSummaryRendererStub <- function() "correlation summary renderer stub"
    module_env$finalQcPlotRendererStub <- function() "final qc plot renderer stub"
    module_env$qcImageRendererStub <- function(assaySlot, plotType, stagePrefix) {
        sprintf("qc image renderer stub %s|%s|%s", assaySlot, plotType, stagePrefix)
    }
    module_env$registerLipidNormStaticQcImageOutputs <- function(output, renderQcImageForAssay) {
        module_env$staticQcImageHelperCalls <- c(
            module_env$staticQcImageHelperCalls,
            list(list(
                output = output,
                renderQcImageForAssay = renderQcImageForAssay
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormAssayLabelOutputs <- function(output, renderAssayLabel) {
        module_env$assayLabelHelperCalls <- c(
            module_env$assayLabelHelperCalls,
            list(list(
                output = output,
                renderAssayLabel = renderAssayLabel
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormLogOutput <- function(output, renderNormLog) {
        module_env$normLogHelperCalls <- c(
            module_env$normLogHelperCalls,
            list(list(
                output = output,
                renderNormLog = renderNormLog
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormItsdSelectionOutput <- function(output, renderItsdSelectionUi) {
        module_env$itsdSelectionHelperCalls <- c(
            module_env$itsdSelectionHelperCalls,
            list(list(
                output = output,
                renderItsdSelectionUi = renderItsdSelectionUi
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormRuvQcOutput <- function(output, renderRuvQcUi) {
        module_env$ruvQcHelperCalls <- c(
            module_env$ruvQcHelperCalls,
            list(list(
                output = output,
                renderRuvQcUi = renderRuvQcUi
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormCorrelationFilterSummaryOutput <- function(output, renderCorrelationFilterSummary) {
        module_env$correlationFilterSummaryHelperCalls <- c(
            module_env$correlationFilterSummaryHelperCalls,
            list(list(
                output = output,
                renderCorrelationFilterSummary = renderCorrelationFilterSummary
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormFinalQcPlotOutput <- function(output, renderFinalQcPlot) {
        module_env$finalQcPlotHelperCalls <- c(
            module_env$finalQcPlotHelperCalls,
            list(list(
                output = output,
                renderFinalQcPlot = renderFinalQcPlot
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormDesignDrivenChoiceObserver <- function(session, workflowData, ...) {
        module_env$designDrivenChoiceObserverHelperCalls <- c(
            module_env$designDrivenChoiceObserverHelperCalls,
            list(list(
                session = session,
                workflowData = workflowData
            ))
        )
        invisible(session)
    }
    module_env$registerLipidNormAssayNameInitializationObserver <- function(workflowData, normData, ...) {
        module_env$assayNameInitializationHelperCalls <- c(
            module_env$assayNameInitializationHelperCalls,
            list(list(
                workflowData = workflowData,
                normData = normData
            ))
        )
        invisible(normData)
    }
    module_env$registerLipidNormRuvCancorOutputs <- function(output, normData) {
        module_env$ruvCancorOutputsHelperCalls <- c(
            module_env$ruvCancorOutputsHelperCalls,
            list(list(
                output = output,
                normData = normData
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormItsdTableOutputs <- function(output, workflowData, normData) {
        module_env$itsdTableHelperCalls <- c(
            module_env$itsdTableHelperCalls,
            list(list(
                output = output,
                workflowData = workflowData,
                normData = normData
            ))
        )
        invisible(output)
    }
    module_env$registerLipidNormItsdSelectionTracking <- function(input, normData) {
        module_env$itsdSelectionTrackingHelperCalls <- c(
            module_env$itsdSelectionTrackingHelperCalls,
            list(list(
                input = input,
                normData = normData
            ))
        )
        invisible(normData)
    }
    module_env$registerLipidNormRunNormalizationObserver <- function(...) {
        args <- list(...)
        module_env$runNormalizationObserverHelperCalls <- c(
            module_env$runNormalizationObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerRunNormalization)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "input$run_normalization"
            )
            do.call(module_env$handleLipidNormRunNormalization, args)
        }

        invisible(args$input)
    }
    module_env$registerLipidNormResetNormalizationObserver <- function(...) {
        args <- list(...)
        module_env$resetNormalizationObserverHelperCalls <- c(
            module_env$resetNormalizationObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerResetNormalization)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "input$reset_normalization"
            )
            do.call(module_env$handleLipidNormResetNormalization, args[names(args) != "input"])
        }

        invisible(args$input)
    }
    module_env$registerLipidNormApplyCorrelationFilterObserver <- function(...) {
        args <- list(...)
        module_env$applyCorrelationFilterObserverHelperCalls <- c(
            module_env$applyCorrelationFilterObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerApplyCorrelationFilter)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "input$apply_correlation_filter"
            )
            do.call(module_env$handleLipidNormApplyCorrelationFilter, args)
        }

        invisible(args$input)
    }
    module_env$registerLipidNormSkipCorrelationFilterObserver <- function(...) {
        args <- list(...)
        module_env$skipCorrelationFilterObserverHelperCalls <- c(
            module_env$skipCorrelationFilterObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerSkipCorrelationFilter)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "input$skip_correlation_filter"
            )
            do.call(module_env$handleLipidNormSkipCorrelationFilter, args[names(args) != "input"])
        }

        invisible(args$input)
    }
    module_env$registerLipidNormExportSessionObserver <- function(...) {
        args <- list(...)
        module_env$exportSessionObserverHelperCalls <- c(
            module_env$exportSessionObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerExportSession)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "input$export_session"
            )
            do.call(module_env$handleLipidNormExportSession, args)
        }

        invisible(args$input)
    }
    module_env$registerLipidNormSelectedTabPreNormalizationObserver <- function(...) {
        args <- list(...)
        module_env$selectedTabPreNormalizationObserverHelperCalls <- c(
            module_env$selectedTabPreNormalizationObserverHelperCalls,
            list(args)
        )

        if (isTRUE(triggerSelectedTab)) {
            module_env$capturedObserverEvents <- c(
                module_env$capturedObserverEvents,
                "selected_tab()"
            )
            do.call(module_env$handleLipidNormSelectedTabPreNormalizationTrigger, c(
                list(selectedTabValue = args$selectedTab()),
                args[names(args) != "selectedTab"]
            ))
        }

        invisible(args$selectedTab)
    }
    module_env$handleLipidNormRunNormalization <- function(...) {
        module_env$runNormalizationHelperCalls <- c(
            module_env$runNormalizationHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormPreNormalizationQc <- function(...) {
        module_env$preNormalizationQcHelperCalls <- c(
            module_env$preNormalizationQcHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormSelectedTabPreNormalizationTrigger <- function(...) {
        module_env$selectedTabAutoTriggerHelperCalls <- c(
            module_env$selectedTabAutoTriggerHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormExportSession <- function(...) {
        module_env$exportSessionHelperCalls <- c(
            module_env$exportSessionHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormApplyCorrelationFilter <- function(...) {
        module_env$applyCorrelationFilterHelperCalls <- c(
            module_env$applyCorrelationFilterHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormSkipCorrelationFilter <- function(...) {
        module_env$skipCorrelationFilterHelperCalls <- c(
            module_env$skipCorrelationFilterHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$handleLipidNormResetNormalization <- function(...) {
        module_env$resetNormalizationHelperCalls <- c(
            module_env$resetNormalizationHelperCalls,
            list(list(...))
        )
        invisible(TRUE)
    }
    module_env$createLipidNormReactiveState <- function(...) {
        module_env$reactiveStateHelperCalls <- c(
            module_env$reactiveStateHelperCalls,
            list(list(...))
        )
        module_env$reactiveStateStub
    }
    module_env$createLipidNormStartupRuntimeStub <- function(...) {
        args <- list(...)
        module_env$startupRuntimeHelperCalls <- c(
            module_env$startupRuntimeHelperCalls,
            list(args)
        )
        list(
            addLog = module_env$addLogStub,
            getPlotAesthetics = module_env$getPlotAestheticsStub,
            generateCompositeFromFiles = module_env$generateCompositeFromFilesStub,
            renderQcImageForAssay = module_env$qcImageRendererStub,
            renderNormLog = module_env$normLogRendererStub,
            renderItsdSelectionUi = module_env$itsdSelectionRendererStub,
            renderRuvQcUi = module_env$ruvQcUiRendererStub,
            renderAssayLabel = module_env$assayLabelRendererStub,
            renderCorrelationFilterSummary = module_env$correlationFilterSummaryRendererStub,
            renderFinalQcPlot = module_env$finalQcPlotRendererStub
        )
    }
    module_env$registerLipidNormStartupObserverRuntimeStub <- function(...) {
        args <- list(...)
        module_env$startupObserverRuntimeHelperCalls <- c(
            module_env$startupObserverRuntimeHelperCalls,
            list(args)
        )
        invisible(args$session)
    }
    module_env$registerLipidNormServerRuntimeStub <- function(...) {
        args <- list(...)
        module_env$serverRuntimeHelperCalls <- c(
            module_env$serverRuntimeHelperCalls,
            list(args)
        )
        invisible(args$input)
    }
    module_env$runLipidNormModuleServerPublicWrapperStub <- function(...) {
        args <- list(...)
        module_env$moduleServerPublicWrapperHelperCalls <- c(
            module_env$moduleServerPublicWrapperHelperCalls,
            list(args)
        )
        invisible(NULL)
    }
    module_env$runLipidNormModuleServerEntryShellStub <- function(...) {
        args <- list(...)
        module_env$moduleServerEntryShellHelperCalls <- c(
            module_env$moduleServerEntryShellHelperCalls,
            list(args)
        )
        invisible(NULL)
    }
    module_env$runLipidNormModuleServerShellStub <- function(...) {
        args <- list(...)
        module_env$moduleServerShellHelperCalls <- c(
            module_env$moduleServerShellHelperCalls,
            list(args)
        )
        invisible(module_env$reactiveStateStub)
    }
    module_env$registerLipidNormPrimaryStartupOutputsStub <- function(...) {
        args <- list(...)
        module_env$primaryStartupOutputsHelperCalls <- c(
            module_env$primaryStartupOutputsHelperCalls,
            list(args)
        )
        invisible(args$output)
    }
    module_env$registerLipidNormItsdSelectionRuntimeStub <- function(...) {
        args <- list(...)
        module_env$itsdSelectionRuntimeHelperCalls <- c(
            module_env$itsdSelectionRuntimeHelperCalls,
            list(args)
        )
        invisible(args$normData)
    }
    module_env$registerLipidNormPostNormalizationOutputsStub <- function(...) {
        args <- list(...)
        module_env$postNormalizationOutputsHelperCalls <- c(
            module_env$postNormalizationOutputsHelperCalls,
            list(args)
        )
        invisible(args$output)
    }
    module_env$buildLipidNormAssayLabelRenderer <- function(normData, ...) {
        module_env$assayLabelBuilderCalls <- c(
            module_env$assayLabelBuilderCalls,
            list(list(normData = normData))
        )
        function(assaySlot) module_env$assayLabelRendererStub(assaySlot)
    }
    module_env$buildLipidNormQcImageRenderer <- function(normData, experimentPaths, ...) {
        module_env$qcImageRendererBuilderCalls <- c(
            module_env$qcImageRendererBuilderCalls,
            list(list(
                normData = normData,
                experimentPaths = experimentPaths
            ))
        )
        module_env$qcImageRendererStub
    }
    module_env$buildLipidNormLogRenderer <- function(normData, ...) {
        module_env$normLogBuilderCalls <- c(
            module_env$normLogBuilderCalls,
            list(list(normData = normData))
        )
        module_env$normLogRendererStub
    }
    module_env$buildLipidNormAddLog <- function(normData, ...) {
        module_env$addLogBuilderCalls <- c(
            module_env$addLogBuilderCalls,
            list(list(normData = normData))
        )
        module_env$addLogStub
    }
    module_env$buildLipidNormPlotAestheticsGetter <- function(input, ...) {
        module_env$getPlotAestheticsBuilderCalls <- c(
            module_env$getPlotAestheticsBuilderCalls,
            list(list(input = input))
        )
        module_env$getPlotAestheticsStub
    }
    module_env$buildLipidNormCompositeFromFilesGenerator <- function(...) {
        module_env$compositeFromFilesBuilderCalls <- c(
            module_env$compositeFromFilesBuilderCalls,
            list(list(...))
        )
        module_env$generateCompositeFromFilesStub
    }
    module_env$buildLipidNormItsdSelectionUiRenderer <- function(normData, ns, ...) {
        module_env$itsdSelectionBuilderCalls <- c(
            module_env$itsdSelectionBuilderCalls,
            list(list(
                normData = normData,
                ns = ns
            ))
        )
        module_env$itsdSelectionRendererStub
    }
    module_env$buildLipidNormRuvQcUiRenderer <- function(normData, ns, ...) {
        module_env$ruvQcUiBuilderCalls <- c(
            module_env$ruvQcUiBuilderCalls,
            list(list(
                normData = normData,
                ns = ns
            ))
        )
        module_env$ruvQcUiRendererStub
    }
    module_env$buildLipidNormCorrelationFilterSummaryRenderer <- function(normData, ...) {
        module_env$correlationFilterSummaryBuilderCalls <- c(
            module_env$correlationFilterSummaryBuilderCalls,
            list(list(normData = normData))
        )
        module_env$correlationFilterSummaryRendererStub
    }
    module_env$buildLipidNormFinalQcPlotRenderer <- function(normData, getPlotAestheticsFn, ...) {
        module_env$finalQcPlotBuilderCalls <- c(
            module_env$finalQcPlotBuilderCalls,
            list(list(
                normData = normData,
                getPlotAestheticsFn = getPlotAestheticsFn
            ))
        )
        module_env$finalQcPlotRendererStub
    }

    selected_tab <- NULL
    mock_selected_tab_state <- NULL
    if (!is.null(selectedTabValue)) {
        selected_tab <- function() selectedTabValue
        mock_selected_tab_state <- structure(list(), class = "LipidomicsAssayData")
    }

    workflow_data <- list(
        state_manager = list(getState = function() mock_selected_tab_state),
        design_matrix = NULL,
        qc_params = NULL
    )
    experiment_paths <- list(lipid_qc_dir = tempdir())

    if (isTRUE(stubModuleServerEntryShell)) {
        module_env$runLipidNormModuleServerEntryShell <- module_env$runLipidNormModuleServerEntryShellStub
    }
    if (isTRUE(stubModuleServerPublicWrapper)) {
        module_env$runLipidNormModuleServerPublicWrapper <- module_env$runLipidNormModuleServerPublicWrapperStub
    }
    if (isTRUE(stubModuleServerShell)) {
        module_env$runLipidNormModuleServerShell <- module_env$runLipidNormModuleServerShellStub
    }
    if (isTRUE(stubStartupRuntime)) {
        module_env$createLipidNormStartupRuntime <- module_env$createLipidNormStartupRuntimeStub
    }
    if (isTRUE(stubServerRuntime)) {
        module_env$registerLipidNormServerRuntime <- module_env$registerLipidNormServerRuntimeStub
    }
    if (isTRUE(stubStartupObserverRuntime)) {
        module_env$registerLipidNormStartupObserverRuntime <- module_env$registerLipidNormStartupObserverRuntimeStub
    }
    if (isTRUE(stubPrimaryStartupOutputs)) {
        module_env$registerLipidNormPrimaryStartupOutputs <- module_env$registerLipidNormPrimaryStartupOutputsStub
    }
    if (isTRUE(stubItsdSelectionRuntime)) {
        module_env$registerLipidNormItsdSelectionRuntime <- module_env$registerLipidNormItsdSelectionRuntimeStub
    }
    if (isTRUE(stubPostNormalizationOutputs)) {
        module_env$registerLipidNormPostNormalizationOutputs <- module_env$registerLipidNormPostNormalizationOutputsStub
    }

    module_env$mod_lipid_norm_server(
        "norm",
        workflow_data = workflow_data,
        experiment_paths = experiment_paths,
        omic_type = "lipidomics",
        experiment_label = "Lipidomics",
        selected_tab = selected_tab
    )

    if (exists("capturedModule", envir = module_env, inherits = FALSE)) {
        session <- shiny:::MockShinySession$new()
        shiny::callModule(module_env$capturedModule, "norm", session = session)
        session$flushReact()
    }

    module_env
}

test_that("registerLipidNormStaticQcImageOutputs keeps the 24 static QC bindings stable", {
    output <- new.env(parent = emptyenv())
    calls <- list()

    render_stub <- function(assay_slot, plot_type, stage_prefix) {
        calls[[length(calls) + 1]] <<- list(
            assay_slot = assay_slot,
            plot_type = plot_type,
            stage_prefix = stage_prefix
        )
        sprintf("%s|%s|%s", assay_slot, plot_type, stage_prefix)
    }

    result <- registerLipidNormStaticQcImageOutputs(output, render_stub)

    expected_ids <- c(
        "pca_post_filter_assay1",
        "pca_post_norm_assay1",
        "pca_ruv_corrected_assay1",
        "pca_post_filter_assay2",
        "pca_post_norm_assay2",
        "pca_ruv_corrected_assay2",
        "density_post_filter_assay1",
        "density_post_norm_assay1",
        "density_ruv_corrected_assay1",
        "density_post_filter_assay2",
        "density_post_norm_assay2",
        "density_ruv_corrected_assay2",
        "rle_post_filter_assay1",
        "rle_post_norm_assay1",
        "rle_ruv_corrected_assay1",
        "rle_post_filter_assay2",
        "rle_post_norm_assay2",
        "rle_ruv_corrected_assay2",
        "correlation_post_filter_assay1",
        "correlation_post_norm_assay1",
        "correlation_ruv_corrected_assay1",
        "correlation_post_filter_assay2",
        "correlation_post_norm_assay2",
        "correlation_ruv_corrected_assay2"
    )
    expected_calls <- list(
        list(assay_slot = 1, plot_type = "pca", stage_prefix = "pre_norm"),
        list(assay_slot = 1, plot_type = "pca", stage_prefix = "post_norm"),
        list(assay_slot = 1, plot_type = "pca", stage_prefix = "ruv_corrected"),
        list(assay_slot = 2, plot_type = "pca", stage_prefix = "pre_norm"),
        list(assay_slot = 2, plot_type = "pca", stage_prefix = "post_norm"),
        list(assay_slot = 2, plot_type = "pca", stage_prefix = "ruv_corrected"),
        list(assay_slot = 1, plot_type = "density", stage_prefix = "pre_norm"),
        list(assay_slot = 1, plot_type = "density", stage_prefix = "post_norm"),
        list(assay_slot = 1, plot_type = "density", stage_prefix = "ruv_corrected"),
        list(assay_slot = 2, plot_type = "density", stage_prefix = "pre_norm"),
        list(assay_slot = 2, plot_type = "density", stage_prefix = "post_norm"),
        list(assay_slot = 2, plot_type = "density", stage_prefix = "ruv_corrected"),
        list(assay_slot = 1, plot_type = "rle", stage_prefix = "pre_norm"),
        list(assay_slot = 1, plot_type = "rle", stage_prefix = "post_norm"),
        list(assay_slot = 1, plot_type = "rle", stage_prefix = "ruv_corrected"),
        list(assay_slot = 2, plot_type = "rle", stage_prefix = "pre_norm"),
        list(assay_slot = 2, plot_type = "rle", stage_prefix = "post_norm"),
        list(assay_slot = 2, plot_type = "rle", stage_prefix = "ruv_corrected"),
        list(assay_slot = 1, plot_type = "correlation", stage_prefix = "pre_norm"),
        list(assay_slot = 1, plot_type = "correlation", stage_prefix = "post_norm"),
        list(assay_slot = 1, plot_type = "correlation", stage_prefix = "ruv_corrected"),
        list(assay_slot = 2, plot_type = "correlation", stage_prefix = "pre_norm"),
        list(assay_slot = 2, plot_type = "correlation", stage_prefix = "post_norm"),
        list(assay_slot = 2, plot_type = "correlation", stage_prefix = "ruv_corrected")
    )

    expect_identical(result, output)
    expect_setequal(ls(output), expected_ids)
    expect_identical(output$pca_post_filter_assay1, "1|pca|pre_norm")
    expect_identical(output$correlation_ruv_corrected_assay2, "2|correlation|ruv_corrected")
    expect_identical(calls, expected_calls)
})

test_that("registerLipidNormAssayLabelOutputs keeps the 8 assay label bindings stable", {
    output <- new.env(parent = emptyenv())
    calls <- list()

    render_stub <- function(assay_slot) {
        calls[[length(calls) + 1]] <<- assay_slot
        paste("Assay slot", assay_slot)
    }

    result <- registerLipidNormAssayLabelOutputs(output, render_stub)

    expected_ids <- c(
        "assay1_label_pca",
        "assay2_label_pca",
        "assay1_label_density",
        "assay2_label_density",
        "assay1_label_rle",
        "assay2_label_rle",
        "assay1_label_correlation",
        "assay2_label_correlation"
    )

    expect_identical(result, output)
    expect_setequal(ls(output), expected_ids)
    expect_identical(output$assay1_label_pca, "Assay slot 1")
    expect_identical(output$assay2_label_correlation, "Assay slot 2")
    expect_identical(unlist(calls), c(1, 2, 1, 2, 1, 2, 1, 2))
})

test_that("buildLipidNormAssayLabelRenderer keeps the assay-label render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- c("Positive Mode", "Negative-Mode")

    render_assay_label <- buildLipidNormAssayLabelRenderer(
        normData = norm_data,
        renderTextFn = function(expr) {
            force(expr)
        }
    )

    expect_true(is.function(render_assay_label))
    expect_identical(render_assay_label(1), "Assay: Positive Mode")
    expect_identical(render_assay_label(2), "Assay: Negative-Mode")
    expect_identical(render_assay_label(3), "Assay 3: (detecting...)")

    norm_data$assay_names <- NULL
    expect_identical(render_assay_label(1), "Assay 1: (detecting...)")
})

test_that("buildLipidNormQcImageRenderer keeps the QC image render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$plot_refresh_trigger <- 1
    norm_data$assay_names <- c("Positive Mode", "Negative-Mode")
    resolved_paths <- character()

    render_qc_image_for_assay <- buildLipidNormQcImageRenderer(
        normData = norm_data,
        experimentPaths = list(lipid_qc_dir = "/tmp/lipid-qc"),
        renderImageFn = function(expr, deleteFile) {
            list(
                image = force(expr),
                deleteFile = deleteFile
            )
        },
        filePathFn = function(...) {
            path <- file.path(...)
            resolved_paths <<- c(resolved_paths, path)
            path
        },
        fileExistsFn = function(path) {
            identical(path, "/tmp/lipid-qc/positive_mode_post_norm_pca.png")
        }
    )

    expect_true(is.function(render_qc_image_for_assay))
    expect_identical(
        render_qc_image_for_assay(1, "pca", "post_norm"),
        list(
            image = list(
                src = "/tmp/lipid-qc/positive_mode_post_norm_pca.png",
                contentType = "image/png",
                width = "100%",
                height = "auto",
                alt = "pca - Positive Mode"
            ),
            deleteFile = FALSE
        )
    )
    expect_identical(
        render_qc_image_for_assay(2, "density", "pre_norm"),
        list(
            image = list(
                src = "",
                alt = "Plot not generated yet: negative_mode_pre_norm_density.png"
            ),
            deleteFile = FALSE
        )
    )

    norm_data$assay_names <- NULL
    expect_identical(
        render_qc_image_for_assay(1, "rle", "post_norm"),
        list(
            image = list(src = "", alt = "Assay not detected yet"),
            deleteFile = FALSE
        )
    )

    expect_identical(
        resolved_paths,
        c(
            "/tmp/lipid-qc/positive_mode_post_norm_pca.png",
            "/tmp/lipid-qc/negative_mode_pre_norm_density.png"
        )
    )
})

test_that("registerLipidNormLogOutput keeps the normalization log binding stable", {
    output <- new.env(parent = emptyenv())
    calls <- 0

    render_stub <- function() {
        calls <<- calls + 1
        "norm log renderer"
    }

    result <- registerLipidNormLogOutput(output, render_stub)

    expect_identical(result, output)
    expect_identical(ls(output), "norm_log")
    expect_identical(output$norm_log, "norm log renderer")
    expect_identical(calls, 1)
})

test_that("buildLipidNormLogRenderer keeps the normalization-log render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$normalization_log <- character()

    render_norm_log <- buildLipidNormLogRenderer(
        normData = norm_data,
        renderTextFn = function(expr) {
            force(expr)
        }
    )

    expect_true(is.function(render_norm_log))
    expect_identical(
        render_norm_log(),
        "Normalization log will appear here as you apply steps..."
    )

    norm_data$normalization_log <- c("Step 1 complete", "Step 2 complete")
    expect_identical(
        render_norm_log(),
        "Step 1 complete\nStep 2 complete"
    )
})

test_that("buildLipidNormAddLog keeps the normalization-log append contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$normalization_log <- "Existing entry"

    add_log <- buildLipidNormAddLog(
        normData = norm_data,
        timeFn = function() as.POSIXct("2026-04-15 12:34:56", tz = "UTC"),
        formatTimeFn = function(x, format) base::format(x, format = format, tz = "UTC")
    )

    expect_true(is.function(add_log))
    add_log("Started normalization")
    add_log("Completed RUV")

    expect_identical(
        norm_data$normalization_log,
        c(
            "Existing entry",
            "[12:34:56] Started normalization",
            "[12:34:56] Completed RUV"
        )
    )
})

test_that("buildLipidNormPlotAestheticsGetter keeps the plot-aesthetics fallback contract stable", {
    input <- new.env(parent = emptyenv())
    input$color_variable <- NULL
    input$shape_variable <- ""

    get_plot_aesthetics <- buildLipidNormPlotAestheticsGetter(
        input = input
    )

    expect_true(is.function(get_plot_aesthetics))
    expect_identical(
        get_plot_aesthetics(),
        list(
            color_var = "group",
            shape_var = "group"
        )
    )

    input$color_variable <- "condition"
    input$shape_variable <- "batch"
    expect_identical(
        get_plot_aesthetics(),
        list(
            color_var = "condition",
            shape_var = "batch"
        )
    )
})

test_that("buildLipidNormCompositeFromFilesGenerator keeps the package-gate contract stable", {
    warnings <- character()
    info_messages <- character()

    generate_composite_from_files <- buildLipidNormCompositeFromFilesGenerator(
        requireNamespaceFn = function(package, quietly = TRUE) {
            expect_identical(package, "patchwork")
            FALSE
        },
        warningFn = function(message) {
            warnings <<- c(warnings, message)
        },
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        }
    )

    expect_true(is.function(generate_composite_from_files))
    expect_null(generate_composite_from_files(character()))
    expect_identical(
        info_messages,
        "[generateCompositeFromFiles] Generating composite from 0 files..."
    )
    expect_identical(
        warnings,
        "patchwork package required for composite generation"
    )
})

test_that("registerLipidNormItsdSelectionOutput keeps the ITSD UI binding stable", {
    output <- new.env(parent = emptyenv())
    calls <- 0

    render_stub <- function() {
        calls <<- calls + 1
        "itsd selection renderer"
    }

    result <- registerLipidNormItsdSelectionOutput(output, render_stub)

    expect_identical(result, output)
    expect_identical(ls(output), "itsd_selection_ui")
    expect_identical(output$itsd_selection_ui, "itsd selection renderer")
    expect_identical(calls, 1)
})

test_that("registerLipidNormRuvQcOutput keeps the RUV QC UI binding stable", {
    output <- new.env(parent = emptyenv())
    calls <- 0

    render_stub <- function() {
        calls <<- calls + 1
        "ruv qc renderer"
    }

    result <- registerLipidNormRuvQcOutput(output, render_stub)

    expect_identical(result, output)
    expect_identical(ls(output), "ruv_qc_ui")
    expect_identical(output$ruv_qc_ui, "ruv qc renderer")
    expect_identical(calls, 1)
})

test_that("registerLipidNormDesignDrivenChoiceObserver keeps the design-driven choice update contract stable", {
    update_calls <- list()
    session <- list(id = "mock-session")
    workflow_data <- list(
        design_matrix = data.frame(
            sample_id = c("s1", "s2"),
            group = c("A", "B"),
            factor1 = c("x", "y"),
            batch = c("b1", "b2"),
            stringsAsFactors = FALSE
        )
    )

    result <- registerLipidNormDesignDrivenChoiceObserver(
        session = session,
        workflowData = workflow_data,
        observeFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        updateSelectInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
            invisible(NULL)
        }
    )

    expect_identical(result, session)
    expect_identical(update_calls, list(
        list(
            session = session,
            inputId = "color_variable",
            choices = c("sample_id", "group", "factor1", "batch"),
            selected = "group"
        ),
        list(
            session = session,
            inputId = "shape_variable",
            choices = c("sample_id", "group", "factor1", "batch"),
            selected = "group"
        ),
        list(
            session = session,
            inputId = "ruv_grouping_variable",
            choices = c("group", "factor1", "batch"),
            selected = "group"
        )
    ))
})

test_that("registerLipidNormAssayNameInitializationObserver keeps the assay-name initialization contract stable", {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                lipid_id_column = "character",
                annotation_id_column = "character"
            )
        )
    }

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            `Positive Mode` = data.frame(feature = 1),
            `Negative Mode` = data.frame(feature = 2)
        )
    )
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- NULL
    norm_data$itsd_selections <- list(`Positive Mode` = 3L)
    info_messages <- character()
    warn_messages <- character()

    result <- registerLipidNormAssayNameInitializationObserver(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        normData = norm_data,
        observeFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(value) {
            expect_false(is.null(value))
            invisible(value)
        },
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarnFn = function(message) {
            warn_messages <<- c(warn_messages, message)
        }
    )

    expect_identical(result, norm_data)
    expect_identical(norm_data$assay_names, c("Positive Mode", "Negative Mode"))
    expect_identical(norm_data$itsd_selections$`Positive Mode`, 3L)
    expect_false("Negative Mode" %in% names(norm_data$itsd_selections))
    expect_identical(info_messages, "Detected assays: Positive Mode, Negative Mode")
    expect_identical(warn_messages, character())
})

test_that("registerLipidNormRunNormalizationObserver keeps the observer shell handoff stable", {
    input <- list(run_normalization = 7L)
    observer_event <- NULL
    handler_args <- NULL

    result <- registerLipidNormRunNormalizationObserver(
        input = input,
        workflowData = list(state_manager = "state-manager"),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        omicType = "lipidomics",
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        getPlotAestheticsFn = function() list(color_var = "group"),
        generateCompositeFromFilesFn = function(...) list(...),
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        handleRunNormalizationFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 7L)
    expect_setequal(
        names(handler_args),
        c(
            "input",
            "workflowData",
            "experimentPaths",
            "omicType",
            "normData",
            "addLog",
            "getPlotAestheticsFn",
            "generateCompositeFromFilesFn"
        )
    )
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(handler_args$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(handler_args$omicType, "lipidomics")
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
    expect_true(is.function(handler_args$getPlotAestheticsFn))
    expect_true(is.function(handler_args$generateCompositeFromFilesFn))
})

test_that("registerLipidNormResetNormalizationObserver keeps the observer shell handoff stable", {
    input <- list(reset_normalization = 9L)
    observer_event <- NULL
    handler_args <- NULL

    result <- registerLipidNormResetNormalizationObserver(
        input = input,
        workflowData = list(state_manager = "state-manager"),
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        handleResetNormalizationFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 9L)
    expect_setequal(
        names(handler_args),
        c("workflowData", "normData", "addLog")
    )
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
})

test_that("registerLipidNormApplyCorrelationFilterObserver keeps the observer shell handoff stable", {
    input <- list(apply_correlation_filter = 11L)
    observer_event <- NULL
    handler_args <- NULL

    result <- registerLipidNormApplyCorrelationFilterObserver(
        input = input,
        workflowData = list(state_manager = "state-manager"),
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        handleApplyCorrelationFilterFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 11L)
    expect_setequal(
        names(handler_args),
        c("input", "workflowData", "normData", "addLog")
    )
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
})

test_that("registerLipidNormSkipCorrelationFilterObserver keeps the observer shell handoff stable", {
    input <- list(skip_correlation_filter = 13L)
    observer_event <- NULL
    handler_args <- NULL

    result <- registerLipidNormSkipCorrelationFilterObserver(
        input = input,
        workflowData = list(state_manager = "state-manager"),
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        handleSkipCorrelationFilterFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 13L)
    expect_setequal(
        names(handler_args),
        c("workflowData", "normData", "addLog")
    )
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
})

test_that("registerLipidNormExportSessionObserver keeps the observer shell handoff stable", {
    input <- list(export_session = 15L)
    observer_event <- NULL
    handler_args <- NULL

    result <- registerLipidNormExportSessionObserver(
        input = input,
        workflowData = list(state_manager = "state-manager"),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        experimentLabel = "Lipidomics",
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        },
        handleExportSessionFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, input)
    expect_identical(observer_event, 15L)
    expect_setequal(
        names(handler_args),
        c("input", "workflowData", "experimentPaths", "experimentLabel", "normData", "addLog")
    )
    expect_identical(handler_args$input, input)
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(handler_args$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(handler_args$experimentLabel, "Lipidomics")
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
})

test_that("registerLipidNormSelectedTabPreNormalizationObserver keeps the observer shell handoff stable", {
    selected_tab <- function() "norm"
    observer_event <- NULL
    observer_ignore_init <- NULL
    handler_args <- NULL

    result <- registerLipidNormSelectedTabPreNormalizationObserver(
        selectedTab = selected_tab,
        workflowData = list(state_manager = "state-manager"),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        normData = new.env(parent = emptyenv()),
        addLog = function(message) message,
        getPlotAestheticsFn = function() list(color_var = "group"),
        observeEventFn = function(eventExpr, handlerExpr, ignoreInit) {
            observer_event <<- eventExpr
            observer_ignore_init <<- ignoreInit
            force(handlerExpr)
            "registered"
        },
        handleSelectedTabPreNormalizationTriggerFn = function(...) {
            handler_args <<- list(...)
            invisible(TRUE)
        }
    )

    expect_identical(result, selected_tab)
    expect_identical(observer_event, "norm")
    expect_false(observer_ignore_init)
    expect_setequal(
        names(handler_args),
        c("selectedTabValue", "workflowData", "experimentPaths", "normData", "addLog", "getPlotAestheticsFn")
    )
    expect_identical(handler_args$selectedTabValue, "norm")
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(handler_args$experimentPaths$lipid_qc_dir, tempdir())
    expect_true(is.environment(handler_args$normData))
    expect_true(is.function(handler_args$addLog))
    expect_true(is.function(handler_args$getPlotAestheticsFn))
})

test_that("buildLipidNormItsdSelectionUiRenderer keeps the ITSD UI render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- c("Positive Mode", "Negative-Mode")
    req_calls <- list()

    render_itsd_selection_ui <- buildLipidNormItsdSelectionUiRenderer(
        normData = norm_data,
        ns = function(id) paste0("norm-", id),
        renderUiFn = function(expr) {
            force(expr)
        },
        reqFn = function(value) {
            req_calls[[length(req_calls) + 1]] <<- value
            invisible(value)
        },
        mapFn = function(.x, .f) {
            lapply(.x, .f)
        },
        wellPanelFn = function(...) {
            list(type = "wellPanel", children = list(...))
        },
        headerFn = function(label) {
            list(type = "h5", label = label)
        },
        dataTableOutputFn = function(id) {
            list(type = "dataTableOutput", id = id)
        },
        brFn = function() {
            list(type = "br")
        },
        tagListFn = function(x) {
            list(type = "tagList", items = x)
        }
    )

    expect_true(is.function(render_itsd_selection_ui))
    expect_identical(
        render_itsd_selection_ui(),
        list(
            type = "tagList",
            items = list(
                list(
                    type = "wellPanel",
                    children = list(
                        list(type = "h5", label = "Assay: Positive Mode"),
                        list(type = "dataTableOutput", id = "norm-itsd_table_positive_mode"),
                        list(type = "br")
                    )
                ),
                list(
                    type = "wellPanel",
                    children = list(
                        list(type = "h5", label = "Assay: Negative-Mode"),
                        list(type = "dataTableOutput", id = "norm-itsd_table_negative_mode"),
                        list(type = "br")
                    )
                )
            )
        )
    )
    expect_identical(req_calls, list(c("Positive Mode", "Negative-Mode")))
})

test_that("buildLipidNormRuvQcUiRenderer keeps the RUV QC UI render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- c("Positive Mode", "Negative-Mode")
    req_calls <- list()

    render_ruv_qc_ui <- buildLipidNormRuvQcUiRenderer(
        normData = norm_data,
        ns = function(id) paste0("norm-", id),
        renderUiFn = function(expr) {
            force(expr)
        },
        reqFn = function(value) {
            req_calls[[length(req_calls) + 1]] <<- value
            invisible(value)
        },
        mapFn = function(.x, .f) {
            lapply(.x, .f)
        },
        tagListFn = function(...) {
            items <- list(...)
            if (length(items) == 1 && is.list(items[[1]]) && is.null(names(items))) {
                items <- items[[1]]
            }
            list(type = "tagList", items = items)
        },
        fluidRowFn = function(...) {
            list(type = "fluidRow", children = list(...))
        },
        columnFn = function(width, ...) {
            list(type = "column", width = width, children = list(...))
        },
        headerFn = function(label, style = NULL) {
            list(type = "h5", label = label, style = style)
        },
        subHeaderFn = function(label) {
            list(type = "h6", label = label)
        },
        plotOutputFn = function(id, height = NULL) {
            list(type = "plotOutput", id = id, height = height)
        },
        resizableFn = function(child) {
            list(type = "resizable", child = child)
        },
        wellPanelFn = function(...) {
            list(type = "wellPanel", children = list(...))
        },
        verbatimTextOutputFn = function(id) {
            list(type = "verbatimTextOutput", id = id)
        },
        brFn = function() {
            list(type = "br")
        },
        dataTableOutputFn = function(id) {
            list(type = "dataTableOutput", id = id)
        },
        hrFn = function() {
            list(type = "hr")
        }
    )

    expect_true(is.function(render_ruv_qc_ui))
    expect_identical(
        render_ruv_qc_ui(),
        list(
            type = "tagList",
            items = list(
                list(
                    type = "tagList",
                    items = list(
                        list(
                            type = "fluidRow",
                            children = list(
                                list(
                                    type = "column",
                                    width = 12,
                                    children = list(
                                        list(
                                            type = "h5",
                                            label = "Assay: Positive Mode",
                                            style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
                                        )
                                    )
                                )
                            )
                        ),
                        list(
                            type = "fluidRow",
                            children = list(
                                list(
                                    type = "column",
                                    width = 8,
                                    children = list(
                                        list(
                                            type = "resizable",
                                            child = list(
                                                type = "plotOutput",
                                                id = "norm-cancor_plot_positive_mode",
                                                height = "400px"
                                            )
                                        )
                                    )
                                ),
                                list(
                                    type = "column",
                                    width = 4,
                                    children = list(
                                        list(
                                            type = "wellPanel",
                                            children = list(
                                                list(type = "h6", label = "Optimization Summary"),
                                                list(type = "verbatimTextOutput", id = "norm-ruv_summary_positive_mode"),
                                                list(type = "br"),
                                                list(type = "h6", label = "Results Table"),
                                                list(type = "dataTableOutput", id = "norm-ruv_table_positive_mode")
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        list(type = "hr")
                    )
                ),
                list(
                    type = "tagList",
                    items = list(
                        list(
                            type = "fluidRow",
                            children = list(
                                list(
                                    type = "column",
                                    width = 12,
                                    children = list(
                                        list(
                                            type = "h5",
                                            label = "Assay: Negative-Mode",
                                            style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
                                        )
                                    )
                                )
                            )
                        ),
                        list(
                            type = "fluidRow",
                            children = list(
                                list(
                                    type = "column",
                                    width = 8,
                                    children = list(
                                        list(
                                            type = "resizable",
                                            child = list(
                                                type = "plotOutput",
                                                id = "norm-cancor_plot_negative_mode",
                                                height = "400px"
                                            )
                                        )
                                    )
                                ),
                                list(
                                    type = "column",
                                    width = 4,
                                    children = list(
                                        list(
                                            type = "wellPanel",
                                            children = list(
                                                list(type = "h6", label = "Optimization Summary"),
                                                list(type = "verbatimTextOutput", id = "norm-ruv_summary_negative_mode"),
                                                list(type = "br"),
                                                list(type = "h6", label = "Results Table"),
                                                list(type = "dataTableOutput", id = "norm-ruv_table_negative_mode")
                                            )
                                        )
                                    )
                                )
                            )
                        ),
                        list(type = "hr")
                    )
                )
            )
        )
    )
    expect_identical(req_calls, list(c("Positive Mode", "Negative-Mode")))
})

test_that("registerLipidNormCorrelationFilterSummaryOutput keeps the correlation summary binding stable", {
    output <- new.env(parent = emptyenv())
    calls <- 0

    render_stub <- function() {
        calls <<- calls + 1
        "correlation summary renderer"
    }

    result <- registerLipidNormCorrelationFilterSummaryOutput(output, render_stub)

    expect_identical(result, output)
    expect_identical(ls(output), "correlation_filter_summary")
    expect_identical(output$correlation_filter_summary, "correlation summary renderer")
    expect_identical(calls, 1)
})

test_that("buildLipidNormCorrelationFilterSummaryRenderer keeps the summary contract stable", {
    if (!methods::isClass("MockLipidNormSummaryObject")) {
        methods::setClass(
            "MockLipidNormSummaryObject",
            slots = c(design_matrix = "data.frame")
        )
    }

    norm_data <- new.env(parent = emptyenv())
    norm_data$correlation_filtering_complete <- TRUE
    norm_data$correlation_results <- list(
        `Positive Mode` = data.frame(pearson_correlation = c(0.8, 0.6))
    )
    norm_data$ruv_corrected_obj <- methods::new(
        "MockLipidNormSummaryObject",
        design_matrix = data.frame(sample = c("s1", "s2", "s3", "s4"))
    )
    norm_data$post_norm_obj <- NULL
    norm_data$correlation_filtered_obj <- methods::new(
        "MockLipidNormSummaryObject",
        design_matrix = data.frame(sample = c("s1", "s2", "s3"))
    )

    render_summary <- buildLipidNormCorrelationFilterSummaryRenderer(
        normData = norm_data,
        renderTextFn = function(expr) {
            force(expr)
        }
    )

    expect_true(is.function(render_summary))
    expect_identical(
        render_summary(),
        paste0(
            "=== Correlation Filtering Summary ===\n",
            "\n[Positive Mode]\n  Sample pairs: 2\n  Correlation: mean=0.700, min=0.600, max=0.800",
            "\n\n[Sample Filtering]\n  Original: 4 samples\n  After filtering: 3 samples\n  Removed: 1 samples"
        )
    )
})

test_that("buildLipidNormFinalQcPlotRenderer keeps the final QC render contract stable", {
    norm_data <- new.env(parent = emptyenv())
    norm_data$correlation_filtering_complete <- TRUE
    norm_data$ruv_complete <- FALSE
    norm_data$post_norm_obj <- structure(list(id = "post-norm"), class = "mockS4")
    norm_data$ruv_corrected_obj <- structure(list(id = "ruv"), class = "mockS4")
    norm_data$correlation_filtered_obj <- structure(list(id = "filtered"), class = "mockS4")

    plot_pca_call <- NULL
    empty_plot_calls <- list()
    aesthetics_calls <- 0

    render_final_qc <- buildLipidNormFinalQcPlotRenderer(
        normData = norm_data,
        getPlotAestheticsFn = function() {
            aesthetics_calls <<- aesthetics_calls + 1
            list(color_var = "batch", shape_var = "group")
        },
        renderPlotFn = function(expr) {
            force(expr)
        },
        reqFn = function(value) {
            expect_true(isTRUE(value))
            invisible(value)
        },
        plotPcaFn = function(theObject, grouping_variable, shape_variable, title) {
            plot_pca_call <<- list(
                theObject = theObject,
                grouping_variable = grouping_variable,
                shape_variable = shape_variable,
                title = title
            )
            list("plot-a", "plot-b")
        },
        wrapPlotsFn = function(plots, ncol) {
            list(plots = plots, ncol = ncol)
        },
        emptyPlotFn = function(label, size) {
            empty_plot_calls[[length(empty_plot_calls) + 1]] <<- list(
                label = label,
                size = size
            )
            list(label = label, size = size)
        }
    )

    expect_true(is.function(render_final_qc))
    expect_identical(
        render_final_qc(),
        list(plots = list("plot-a", "plot-b"), ncol = 1)
    )
    expect_identical(plot_pca_call$theObject, norm_data$correlation_filtered_obj)
    expect_identical(plot_pca_call$grouping_variable, "batch")
    expect_identical(plot_pca_call$shape_variable, "group")
    expect_identical(plot_pca_call$title, "Final QC - PCA")
    expect_identical(aesthetics_calls, 1)
    expect_identical(empty_plot_calls, list())
})

test_that("registerLipidNormFinalQcPlotOutput keeps the final QC plot binding stable", {
    output <- new.env(parent = emptyenv())
    calls <- 0

    render_stub <- function() {
        calls <<- calls + 1
        "final qc plot renderer"
    }

    result <- registerLipidNormFinalQcPlotOutput(output, render_stub)

    expect_identical(result, output)
    expect_identical(ls(output), "final_qc_plot")
    expect_identical(output$final_qc_plot, "final qc plot renderer")
    expect_identical(calls, 1)
})

test_that("registerLipidNormRuvCancorOutputs keeps the per-assay RUV cancor bindings stable", {
    output <- new.env(parent = emptyenv())
    datatable_calls <- list()
    success_table <- data.frame(
        k = c(1, 2)
        , separation_score = c(0.12, 0.34)
    )
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- c("Positive Mode", "Negative Mode")
    norm_data$ruv_optimization_results <- list(
        `Positive Mode` = list(
            success = TRUE
            , cancor_plot = "positive-cancor-plot"
            , best_k = 3
            , best_percentage = 12.5
            , separation_score = 0.87654
            , control_genes_index = c(TRUE, FALSE, TRUE)
            , optimization_results = success_table
        )
        , `Negative Mode` = list(
            success = FALSE
            , error = "optimizer failed"
        )
    )

    result <- registerLipidNormRuvCancorOutputs(
        output = output
        , normData = norm_data
        , observeFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        }
        , reqFn = function(value) {
            expect_true(!is.null(value))
            invisible(value)
        }
        , renderPlotFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        }
        , renderTextFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        }
        , renderDataTableFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        }
        , datatableFn = function(data, options, rownames) {
            datatable_calls[[length(datatable_calls) + 1]] <<- list(
                data = data
                , options = options
                , rownames = rownames
            )
            list(
                data = data
                , options = options
                , rownames = rownames
            )
        }
    )

    expected_ids <- c(
        "cancor_plot_positive_mode",
        "ruv_summary_positive_mode",
        "ruv_table_positive_mode",
        "cancor_plot_negative_mode",
        "ruv_summary_negative_mode",
        "ruv_table_negative_mode"
    )

    expect_identical(result, output)
    expect_setequal(ls(output), expected_ids)
    expect_identical(output$cancor_plot_positive_mode, "positive-cancor-plot")
    expect_s3_class(output$cancor_plot_negative_mode, "ggplot")
    expect_identical(
        output$ruv_summary_positive_mode,
        paste(
            "Best k: 3",
            "Best %: 12.5",
            "Separation: 0.8765",
            "Controls: 2",
            sep = "\n"
        )
    )
    expect_identical(output$ruv_summary_negative_mode, "Failed: optimizer failed")
    expect_identical(output$ruv_table_positive_mode, list(
        data = success_table
        , options = list(pageLength = 5, dom = "t")
        , rownames = FALSE
    ))
    expect_null(output$ruv_table_negative_mode)
    expect_identical(datatable_calls, list(list(
        data = success_table
        , options = list(pageLength = 5, dom = "t")
        , rownames = FALSE
    )))
})

test_that("registerLipidNormItsdTableOutputs keeps the per-assay ITSD table bindings stable", {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                lipid_id_column = "character",
                annotation_id_column = "character"
            )
        )
    }

    output <- new.env(parent = emptyenv())
    build_calls <- list()
    datatable_calls <- list()
    style_calls <- list()
    round_calls <- list()
    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            `Positive Mode` = data.frame(raw = 1:2),
            `Negative Mode` = data.frame(raw = 3:4)
        ),
        lipid_id_column = "lipid_id",
        annotation_id_column = "annotation"
    )
    selection_table <- data.frame(
        lipid_id = c("L1", "L2"),
        annotation = c("A", "B"),
        mean_intensity = c(12.345, 67.89),
        cv_percent = c(4.567, 8.91),
        is_candidate = c(TRUE, FALSE)
    )

    result <- registerLipidNormItsdTableOutputs(
        output = output,
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        normData = list(
            assay_names = c("Positive Mode", "Negative Mode")
        ),
        observeFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(value) {
            expect_true(!is.null(value))
            invisible(value)
        },
        renderDataTableFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
        },
        datatableFn = function(data, selection, filter, options, rownames) {
            datatable_calls[[length(datatable_calls) + 1]] <<- list(
                data = data,
                selection = selection,
                filter = filter,
                options = options,
                rownames = rownames
            )
            list(
                data = data,
                selection = selection,
                filter = filter,
                options = options,
                rownames = rownames
            )
        },
        formatStyleFn = function(widget, column, backgroundColor) {
            style_calls[[length(style_calls) + 1]] <<- list(
                column = column,
                backgroundColor = backgroundColor
            )
            widget$style <- list(column = column, backgroundColor = backgroundColor)
            widget
        },
        styleEqualFn = function(levels, values) {
            list(levels = levels, values = values)
        },
        formatRoundFn = function(widget, columns, digits) {
            round_calls[[length(round_calls) + 1]] <<- list(
                columns = columns,
                digits = digits
            )
            widget$round <- list(columns = columns, digits = digits)
            widget
        },
        buildLipidItsdSelectionTableFn = function(assay_data, lipid_id_col, annotation_cols) {
            build_calls[[length(build_calls) + 1]] <<- list(
                assay_data = assay_data,
                lipid_id_col = lipid_id_col,
                annotation_cols = annotation_cols
            )
            selection_table
        }
    )

    expected_ids <- c("itsd_table_positive_mode", "itsd_table_negative_mode")
    expected_style <- list(levels = TRUE, values = "#d4edda")

    expect_identical(result, output)
    expect_setequal(ls(output), expected_ids)
    expect_identical(output$itsd_table_positive_mode$data, selection_table)
    expect_identical(
        output$itsd_table_positive_mode$selection,
        list(mode = "multiple", selected = 1L)
    )
    expect_identical(output$itsd_table_positive_mode$filter, "top")
    expect_identical(
        output$itsd_table_positive_mode$options,
        list(pageLength = 10, scrollX = TRUE, order = list(list(4, "desc"), list(3, "asc")))
    )
    expect_false(output$itsd_table_positive_mode$rownames)
    expect_identical(output$itsd_table_positive_mode$style, list(
        column = "is_candidate",
        backgroundColor = expected_style
    ))
    expect_identical(output$itsd_table_positive_mode$round, list(
        columns = c("mean_intensity", "cv_percent"),
        digits = 2
    ))
    expect_identical(output$itsd_table_negative_mode$data, selection_table)
    expect_length(build_calls, 2)
    expect_identical(build_calls[[1]]$assay_data, current_s4@lipid_data[["Positive Mode"]])
    expect_identical(build_calls[[2]]$assay_data, current_s4@lipid_data[["Negative Mode"]])
    expect_identical(build_calls[[1]]$lipid_id_col, "lipid_id")
    expect_identical(build_calls[[1]]$annotation_cols, "annotation")
    expect_identical(style_calls, rep(list(list(
        column = "is_candidate",
        backgroundColor = expected_style
    )), 2))
    expect_identical(round_calls, rep(list(list(
        columns = c("mean_intensity", "cv_percent"),
        digits = 2
    )), 2))
})

test_that("registerLipidNormItsdSelectionTracking keeps the per-assay ITSD selection tracking stable", {
    observed_events <- list()
    log_messages <- character()
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- c("Positive Mode", "Negative Mode")
    norm_data$itsd_selections <- list()
    input <- list(
        itsd_table_positive_mode_rows_selected = c(1L, 3L),
        itsd_table_negative_mode_rows_selected = integer()
    )

    result <- registerLipidNormItsdSelectionTracking(
        input = input,
        normData = norm_data,
        observeFn = function(expr) {
            eval(substitute(expr), envir = parent.frame())
            invisible(NULL)
        },
        reqFn = function(value) {
            expect_identical(value, norm_data$assay_names)
            invisible(value)
        },
        walkFn = function(.x, .f) purrr::walk(.x, .f),
        observeEventFn = function(eventExpr, handlerExpr, ignoreNULL) {
            observed_events[[length(observed_events) + 1]] <<- list(
                value = eventExpr,
                ignoreNULL = ignoreNULL
            )
            eval(substitute(handlerExpr), envir = parent.frame())
            invisible(NULL)
        },
        logInfoFn = function(message) {
            log_messages <<- c(log_messages, message)
        }
    )

    expect_identical(result, norm_data)
    expect_identical(norm_data$itsd_selections[["Positive Mode"]], c(1L, 3L))
    expect_identical(norm_data$itsd_selections[["Negative Mode"]], integer())
    expect_identical(observed_events, list(
        list(value = c(1L, 3L), ignoreNULL = FALSE),
        list(value = integer(), ignoreNULL = FALSE)
    ))
    expect_identical(log_messages, c(
        "ITSD selection updated for Positive Mode : 2 features selected",
        "ITSD selection updated for Negative Mode : 0 features selected"
    ))
})

test_that("handleLipidNormRunNormalization keeps the skip-path contract stable", {
    progress_calls <- list()
    qc_calls <- list()
    saved_states <- list()
    composite_call <- NULL
    saved_plot <- NULL
    notifications <- list()
    log_messages <- character()
    log_transform_call <- NULL
    current_s4 <- structure(list(id = "current"), class = "mockS4")
    log2_s4 <- structure(list(id = "log2"), class = "mockS4")
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- list(
        getState = function() current_s4,
        saveState = function(state_name, s4_data_object, config_object, description) {
            saved_states[[length(saved_states) + 1]] <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
        }
    )
    workflow_data$config_list <- list(mode = "test")
    norm_data <- new.env(parent = emptyenv())
    norm_data$itsd_selections <- list()
    norm_data$assay_names <- "Positive Mode"
    norm_data$plot_refresh_trigger <- 0
    norm_data$normalization_complete <- FALSE
    norm_data$ruv_complete <- FALSE
    norm_data$ruv_optimization_results <- list()
    norm_data$ruv_corrected_obj <- NULL
    norm_data$post_norm_obj <- NULL
    norm_data$post_log2_obj <- NULL
    norm_data$post_filter_obj <- NULL
    experiment_paths <- list(lipid_qc_dir = tempdir())

    result <- handleLipidNormRunNormalization(
        input = list(
            apply_itsd = FALSE,
            itsd_aggregation = "mean",
            log_offset = 1,
            norm_method = "none",
            ruv_mode = "skip",
            auto_percentage_min = 5,
            auto_percentage_max = 25,
            ruv_grouping_variable = "group",
            separation_metric = "metric",
            k_penalty_weight = 2,
            adaptive_k_penalty = TRUE,
            ruv_k = 3,
            ruv_percentage = 10
        ),
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        omicType = "lipidomics",
        normData = norm_data,
        addLog = function(message) {
            log_messages <<- c(log_messages, message)
        },
        getPlotAestheticsFn = function() {
            list(color_var = "condition", shape_var = "batch")
        },
        reqFn = function(value) {
            expect_false(is.null(value))
            invisible(value)
        },
        withProgressFn = function(message, value, expr) {
            progress_calls[[length(progress_calls) + 1]] <<- list(
                message = message,
                value = value
            )
            force(expr)
        },
        incProgressFn = function(value, detail) {
            progress_calls[[length(progress_calls) + 1]] <<- list(
                value = value,
                detail = detail
            )
        },
        generateLipidQcPlotsFn = function(theObject, experiment_paths, stage, grouping_variable, shape_variable) {
            qc_calls[[length(qc_calls) + 1]] <<- list(
                theObject = theObject,
                experiment_paths = experiment_paths,
                stage = stage,
                grouping_variable = grouping_variable,
                shape_variable = shape_variable
            )
            invisible(NULL)
        },
        normaliseUntransformedDataFn = function(...) {
            stop("unexpected ITSD normalization")
        },
        logTransformAssaysFn = function(theObject, offset) {
            log_transform_call <<- list(theObject = theObject, offset = offset)
            log2_s4
        },
        normaliseBetweenSamplesFn = function(...) {
            stop("unexpected between-sample normalization")
        },
        runLipidPerAssayRuvOptimizationFn = function(...) {
            stop("unexpected RUV optimization")
        },
        extractLipidBestKPerAssayFn = function(...) {
            stop("unexpected best-k extraction")
        },
        extractLipidCtrlPerAssayFn = function(...) {
            stop("unexpected control extraction")
        },
        ruvIII_C_VaryingFn = function(...) {
            stop("unexpected RUV correction")
        },
        generateCompositeFromFilesFn = function(plot_files, ncol, row_labels, column_labels) {
            composite_call <<- list(
                plot_files = plot_files,
                ncol = ncol,
                row_labels = row_labels,
                column_labels = column_labels
            )
            list(plot = "composite-plot", width = 8, height = 10)
        },
        savePlotFn = function(plot, dir, filename, width, height, dpi, limitsize) {
            saved_plot <<- list(
                plot = plot,
                dir = dir,
                filename = filename,
                width = width,
                height = height,
                dpi = dpi,
                limitsize = limitsize
            )
            invisible(NULL)
        },
        dirExistsFn = function(path) {
            expect_identical(path, tempdir())
            TRUE
        },
        showNotificationFn = function(...) {
            notifications[[length(notifications) + 1]] <<- list(...)
        },
        logWarnFn = function(message) {
            stop(sprintf("unexpected logWarn call: %s", message))
        },
        logErrorFn = function(message) {
            stop(sprintf("unexpected logError call: %s", message))
        }
    )

    expect_true(result)
    expect_identical(norm_data$post_filter_obj, current_s4)
    expect_identical(norm_data$post_log2_obj, log2_s4)
    expect_identical(norm_data$post_norm_obj, log2_s4)
    expect_identical(norm_data$ruv_corrected_obj, log2_s4)
    expect_true(norm_data$normalization_complete)
    expect_true(norm_data$ruv_complete)
    expect_identical(norm_data$plot_refresh_trigger, 1)
    expect_identical(log_transform_call, list(theObject = current_s4, offset = 1))
    expect_length(qc_calls, 2)
    expect_identical(qc_calls[[1]]$stage, "post_filter")
    expect_identical(qc_calls[[2]]$stage, "post_norm")
    expect_identical(qc_calls[[1]]$grouping_variable, "condition")
    expect_identical(qc_calls[[1]]$shape_variable, "batch")
    expect_identical(saved_states, list(
        list(
            state_name = "lipid_log2",
            s4_data_object = log2_s4,
            config_object = workflow_data$config_list,
            description = "Log2 transformation (offset: 1 )"
        ),
        list(
            state_name = "lipid_normalized",
            s4_data_object = log2_s4,
            config_object = workflow_data$config_list,
            description = "Between-sample normalization (method: none )"
        )
    ))
    expect_identical(progress_calls[[1]], list(
        message = "Running normalization pipeline...",
        value = 0
    ))
    expect_identical(
        vapply(progress_calls[-1], `[[`, character(1), "detail"),
        c(
            "Capturing pre-normalization state...",
            "Applying ITSD normalization...",
            "Applying log2 transformation...",
            "Applying between-sample normalization...",
            "Running RUV-III batch correction..."
        )
    )
    expect_identical(composite_call$ncol, 2)
    expect_identical(composite_call$column_labels, c("Pre-Normalisation", "Post-Normalisation"))
    expect_length(composite_call$plot_files, 8)
    expect_identical(names(composite_call$row_labels), c(
        "positive_mode_pca",
        "positive_mode_density",
        "positive_mode_rle",
        "positive_mode_correlation"
    ))
    expect_identical(saved_plot, list(
        plot = "composite-plot",
        dir = tempdir(),
        filename = "lipidomics_composite_QC_figure",
        width = 8,
        height = 10,
        dpi = 150,
        limitsize = FALSE
    ))
    expect_identical(log_messages, c(
        "Starting normalization pipeline...",
        "Post-filtering state captured",
        "Pre-normalization QC plots generated",
        "ITSD normalization skipped",
        "Applying log2 transformation (offset: 1 )",
        "Log2 transformation complete",
        "Between-sample normalization complete",
        "Post-normalization QC plots generated",
        "RUV-III skipped",
        "Generating composite QC figure...",
        sprintf("Composite QC figure saved to: %s", file.path(tempdir(), "composite_QC_figure")),
        "Normalization pipeline complete!"
    ))
    expect_identical(notifications, list(list(
        "Normalization pipeline complete!",
        type = "message"
    )))
})

test_that("handleLipidNormPreNormalizationQc keeps the pre-QC generation contract stable", {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(lipid_data = "list")
        )
    }

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            `Positive Mode` = data.frame(feature = 1),
            `Negative Mode` = data.frame(feature = 2)
        )
    )
    norm_data <- new.env(parent = emptyenv())
    norm_data$assay_names <- NULL
    norm_data$plot_refresh_trigger <- 0
    norm_data$pre_norm_qc_generated <- FALSE
    qc_plot_call <- NULL
    log_messages <- character()
    info_messages <- character()
    warn_messages <- character()
    error_messages <- character()

    result <- handleLipidNormPreNormalizationQc(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        normData = norm_data,
        addLog = function(message) {
            log_messages <<- c(log_messages, message)
        },
        getPlotAestheticsFn = function() {
            list(color_var = "group", shape_var = "batch")
        },
        reqFn = function(value) {
            expect_false(is.null(value))
            invisible(value)
        },
        generateLipidQcPlotsFn = function(theObject, experiment_paths, stage, grouping_variable, shape_variable) {
            qc_plot_call <<- list(
                theObject = theObject,
                experiment_paths = experiment_paths,
                stage = stage,
                grouping_variable = grouping_variable,
                shape_variable = shape_variable
            )
            invisible(TRUE)
        },
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarnFn = function(message) {
            warn_messages <<- c(warn_messages, message)
        },
        logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        }
    )

    expect_true(result)
    expect_identical(norm_data$assay_names, c("Positive Mode", "Negative Mode"))
    expect_identical(norm_data$plot_refresh_trigger, 1)
    expect_true(norm_data$pre_norm_qc_generated)
    expect_identical(qc_plot_call$theObject, current_s4)
    expect_identical(qc_plot_call$experiment_paths, list(lipid_qc_dir = tempdir()))
    expect_identical(qc_plot_call$stage, "post_filter")
    expect_identical(qc_plot_call$grouping_variable, "group")
    expect_identical(qc_plot_call$shape_variable, "batch")
    expect_identical(log_messages, character())
    expect_true(any(grepl("GENERATING PRE-NORMALIZATION QC PLOTS", info_messages, fixed = TRUE)))
    expect_true(any(grepl("Set assay names:", info_messages, fixed = TRUE)))
    expect_true(any(grepl("generated successfully", info_messages, fixed = TRUE)))
    expect_identical(warn_messages, character())
    expect_identical(error_messages, character())
})

test_that("handleLipidNormSelectedTabPreNormalizationTrigger keeps the auto-trigger contract stable", {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(lipid_data = "list")
        )
    }

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(`Positive Mode` = data.frame(feature = 1))
    )
    norm_data <- new.env(parent = emptyenv())
    norm_data$pre_norm_qc_generated <- FALSE
    helper_call <- NULL
    progress_calls <- list()
    info_messages <- character()

    result <- handleLipidNormSelectedTabPreNormalizationTrigger(
        selectedTabValue = "norm",
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        normData = norm_data,
        addLog = function(message) {
            invisible(message)
        },
        getPlotAestheticsFn = function() {
            list(color_var = "group", shape_var = "batch")
        },
        reqFn = function(value) {
            expect_false(is.null(value))
            invisible(value)
        },
        withProgressFn = function(message, value, expr) {
            progress_calls[[length(progress_calls) + 1]] <<- list(
                message = message,
                value = value
            )
            force(expr)
        },
        handlePreNormalizationQcFn = function(...) {
            helper_call <<- list(...)
            invisible(FALSE)
        },
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        }
    )

    expect_true(result)
    expect_true(norm_data$pre_norm_qc_generated)
    expect_identical(progress_calls, list(list(
        message = "Generating Pre-Normalization QC...",
        value = 0.5
    )))
    expect_identical(info_messages, c(
        "Normalization tab selected - checking if pre-QC needed",
        "Auto-triggering pre-normalization QC plots"
    ))
    expect_setequal(
        names(helper_call),
        c("workflowData", "experimentPaths", "normData", "addLog", "getPlotAestheticsFn")
    )
    expect_identical(helper_call$workflowData$state_manager$getState(), current_s4)
    expect_identical(helper_call$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(helper_call$normData, norm_data)
    expect_true(is.function(helper_call$addLog))
    expect_true(is.function(helper_call$getPlotAestheticsFn))
})

test_that("handleLipidNormExportSession keeps the incomplete-normalization warning contract stable", {
    notification <- NULL
    logged_message <- NULL
    state_manager <- list(marker = "state-manager")

    result <- handleLipidNormExportSession(
        input = list(),
        workflowData = list(state_manager = state_manager),
        experimentPaths = list(),
        experimentLabel = "Lipidomics",
        normData = list(normalization_complete = FALSE),
        addLog = function(message) stop(sprintf("unexpected addLog call: %s", message)),
        reqFn = function(value) {
            expect_identical(value, state_manager)
            invisible(value)
        },
        showNotificationFn = function(message, type, duration) {
            notification <<- list(
                message = message,
                type = type,
                duration = duration
            )
        },
        logInfoFn = function(message) {
            logged_message <<- message
        }
    )

    expect_false(result)
    expect_identical(logged_message, "=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
    expect_identical(notification, list(
        message = "Please complete normalization before exporting session data.",
        type = "warning",
        duration = 5
    ))
})

test_that("handleLipidNormSkipCorrelationFilter keeps the completion contract stable", {
    saved_state <- NULL
    notification <- NULL
    log_message <- NULL
    post_norm_obj <- structure(list(id = "post-norm"), class = "mockS4")
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- list(
        saveState = function(state_name, s4_data_object, config_object, description) {
            saved_state <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
        }
    )
    workflow_data$config_list <- list(mode = "test")
    workflow_data$tab_status <- list(quality_control = "pending", normalization = "pending")

    result <- handleLipidNormSkipCorrelationFilter(
        workflowData = workflow_data,
        normData = list(
            ruv_complete = FALSE,
            normalization_complete = TRUE,
            ruv_corrected_obj = NULL,
            post_norm_obj = post_norm_obj
        ),
        addLog = function(message) {
            log_message <<- message
        },
        reqFn = function(value) {
            expect_true(!is.null(value))
            invisible(value)
        },
        showNotificationFn = function(message, type) {
            notification <<- list(message = message, type = type)
        }
    )

    expect_true(result)
    expect_identical(saved_state, list(
        state_name = "lipid_norm_complete",
        s4_data_object = post_norm_obj,
        config_object = workflow_data$config_list,
        description = "Normalization complete (correlation filtering skipped)"
    ))
    expect_identical(workflow_data$tab_status$quality_control, "complete")
    expect_identical(workflow_data$tab_status$normalization, "complete")
    expect_identical(log_message, "Correlation filtering skipped - ready for DE analysis")
    expect_identical(notification, list(
        message = "Normalization complete! Proceeding to DE analysis.",
        type = "message"
    ))
})

test_that("handleLipidNormApplyCorrelationFilter keeps the completion contract stable", {
    saved_state <- NULL
    correlation_call <- NULL
    filter_call <- NULL
    notifications <- list()
    removed_notification <- NULL
    log_messages <- character()
    info_message <- NULL
    current_s4 <- structure(list(id = "ruv"), class = "mockS4")
    filtered_s4 <- structure(list(id = "filtered"), class = "mockS4")
    corr_results <- list(
        pos = data.frame(pearson_correlation = c(0.91, 0.84))
    )
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- list(
        saveState = function(state_name, s4_data_object, config_object, description) {
            saved_state <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
        }
    )
    workflow_data$config_list <- list(mode = "test")
    workflow_data$tab_status <- list(quality_control = "pending", normalization = "pending")
    norm_data <- new.env(parent = emptyenv())
    norm_data$ruv_complete <- TRUE
    norm_data$normalization_complete <- TRUE
    norm_data$ruv_corrected_obj <- current_s4
    norm_data$post_norm_obj <- NULL
    norm_data$correlation_results <- list()
    norm_data$correlation_filtered_obj <- NULL
    norm_data$correlation_filtering_complete <- FALSE

    result <- handleLipidNormApplyCorrelationFilter(
        input = list(
            min_pearson_correlation_threshold = 0.65,
            ruv_grouping_variable = "group"
        ),
        workflowData = workflow_data,
        normData = norm_data,
        addLog = function(message) {
            log_messages <<- c(log_messages, message)
        },
        reqFn = function(value) {
            expect_true(!is.null(value))
            invisible(value)
        },
        showNotificationFn = function(...) {
            notifications[[length(notifications) + 1]] <<- list(...)
        },
        removeNotificationFn = function(id) {
            removed_notification <<- id
        },
        pearsonCorForSamplePairsFn = function(theObject, correlation_group) {
            correlation_call <<- list(
                theObject = theObject,
                correlation_group = correlation_group
            )
            corr_results
        },
        filterSamplesByLipidCorrelationThresholdFn = function(theObject, pearson_correlation_per_pair, min_pearson_correlation_threshold) {
            filter_call <<- list(
                theObject = theObject,
                pearson_correlation_per_pair = pearson_correlation_per_pair,
                min_pearson_correlation_threshold = min_pearson_correlation_threshold
            )
            filtered_s4
        },
        logInfoFn = function(message) {
            info_message <<- message
        },
        logErrorFn = function(message) {
            stop(sprintf("unexpected logError call: %s", message))
        }
    )

    expect_true(result)
    expect_identical(correlation_call, list(
        theObject = current_s4,
        correlation_group = "group"
    ))
    expect_identical(filter_call, list(
        theObject = current_s4,
        pearson_correlation_per_pair = corr_results,
        min_pearson_correlation_threshold = 0.65
    ))
    expect_identical(norm_data$correlation_results, corr_results)
    expect_identical(norm_data$correlation_filtered_obj, filtered_s4)
    expect_true(norm_data$correlation_filtering_complete)
    expect_identical(saved_state, list(
        state_name = "lipid_correlation_filtered",
        s4_data_object = filtered_s4,
        config_object = workflow_data$config_list,
        description = "Correlation filtering (threshold: 0.65 )"
    ))
    expect_identical(workflow_data$tab_status$quality_control, "complete")
    expect_identical(workflow_data$tab_status$normalization, "complete")
    expect_identical(log_messages, c(
        "Applying correlation filter (threshold: 0.65 )",
        "Correlation filtering complete"
    ))
    expect_identical(removed_notification, "corr_working")
    expect_identical(info_message, "Calculating Pearson correlations per sample pair...")
    expect_length(notifications, 2)
    expect_identical(notifications[[1]], list(
        "Applying correlation filter...",
        id = "corr_working",
        duration = NULL
    ))
    expect_identical(notifications[[2]], list(
        "Correlation filtering complete! Ready for DE analysis.",
        type = "message"
    ))
})

test_that("handleLipidNormResetNormalization keeps the reset contract stable", {
    saved_state <- NULL
    notification <- NULL
    log_message <- NULL
    post_filter_obj <- structure(list(id = "post-filter"), class = "mockS4")
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$state_manager <- list(
        saveState = function(state_name, s4_data_object, config_object, description) {
            saved_state <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
        }
    )
    workflow_data$config_list <- list(mode = "test")
    norm_data <- new.env(parent = emptyenv())
    norm_data$post_filter_obj <- post_filter_obj
    norm_data$normalization_complete <- TRUE
    norm_data$ruv_complete <- TRUE
    norm_data$correlation_filtering_complete <- TRUE
    norm_data$post_norm_obj <- structure(list(id = "post-norm"), class = "mockS4")
    norm_data$ruv_corrected_obj <- structure(list(id = "ruv"), class = "mockS4")
    norm_data$correlation_filtered_obj <- structure(list(id = "corr"), class = "mockS4")
    norm_data$ruv_optimization_results <- list(pos = list(success = TRUE))

    result <- handleLipidNormResetNormalization(
        workflowData = workflow_data,
        normData = norm_data,
        addLog = function(message) {
            log_message <<- message
        },
        reqFn = function(value) {
            expect_true(!is.null(value))
            invisible(value)
        },
        showNotificationFn = function(message, type) {
            notification <<- list(message = message, type = type)
        }
    )

    expect_true(result)
    expect_identical(saved_state, list(
        state_name = "lipid_reset",
        s4_data_object = post_filter_obj,
        config_object = workflow_data$config_list,
        description = "Reset to pre-normalization state"
    ))
    expect_false(norm_data$normalization_complete)
    expect_false(norm_data$ruv_complete)
    expect_false(norm_data$correlation_filtering_complete)
    expect_null(norm_data$post_norm_obj)
    expect_null(norm_data$ruv_corrected_obj)
    expect_null(norm_data$correlation_filtered_obj)
    expect_identical(norm_data$ruv_optimization_results, list())
    expect_identical(log_message, "Reset to pre-normalization state")
    expect_identical(notification, list(
        message = "Reset to pre-normalization state",
        type = "message"
    ))
})

test_that("mod_lipid_norm_server routes startup output registration through the top-level seams", {
    module_env <- loadLipidNormModuleHarness()

    expect_length(module_env$reactiveStateHelperCalls, 1)
    expect_identical(module_env$reactiveStateHelperCalls[[1]], list())
    expect_length(module_env$addLogBuilderCalls, 1)
    expect_setequal(
        names(module_env$addLogBuilderCalls[[1]]),
        "normData"
    )
    expect_identical(module_env$addLogBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_length(module_env$getPlotAestheticsBuilderCalls, 1)
    expect_setequal(
        names(module_env$getPlotAestheticsBuilderCalls[[1]]),
        "input"
    )
    expect_false(is.null(module_env$getPlotAestheticsBuilderCalls[[1]]$input))
    expect_length(module_env$compositeFromFilesBuilderCalls, 1)
    expect_identical(module_env$compositeFromFilesBuilderCalls[[1]], list())
    expect_length(module_env$staticQcImageHelperCalls, 1)
    expect_length(module_env$qcImageRendererBuilderCalls, 1)
    expect_setequal(
        names(module_env$qcImageRendererBuilderCalls[[1]]),
        c("normData", "experimentPaths")
    )
    expect_identical(module_env$qcImageRendererBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_identical(module_env$qcImageRendererBuilderCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(
        module_env$staticQcImageHelperCalls[[1]]$renderQcImageForAssay(2, "rle", "post_norm"),
        "qc image renderer stub 2|rle|post_norm"
    )
    expect_length(module_env$assayLabelBuilderCalls, 1)
    expect_setequal(
        names(module_env$assayLabelBuilderCalls[[1]]),
        "normData"
    )
    expect_identical(module_env$assayLabelBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_length(module_env$assayLabelHelperCalls, 1)
    expect_identical(
        module_env$assayLabelHelperCalls[[1]]$renderAssayLabel(2),
        "assay label renderer stub 2"
    )
    expect_length(module_env$normLogBuilderCalls, 1)
    expect_setequal(
        names(module_env$normLogBuilderCalls[[1]]),
        "normData"
    )
    expect_identical(module_env$normLogBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_length(module_env$normLogHelperCalls, 1)
    expect_identical(
        module_env$normLogHelperCalls[[1]]$renderNormLog,
        module_env$normLogRendererStub
    )
    expect_length(module_env$itsdSelectionBuilderCalls, 1)
    expect_setequal(
        names(module_env$itsdSelectionBuilderCalls[[1]]),
        c("normData", "ns")
    )
    expect_identical(module_env$itsdSelectionBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_true(is.function(module_env$itsdSelectionBuilderCalls[[1]]$ns))
    expect_length(module_env$itsdSelectionHelperCalls, 1)
    expect_identical(
        module_env$itsdSelectionHelperCalls[[1]]$renderItsdSelectionUi,
        module_env$itsdSelectionRendererStub
    )
    expect_length(module_env$ruvQcUiBuilderCalls, 1)
    expect_setequal(
        names(module_env$ruvQcUiBuilderCalls[[1]]),
        c("normData", "ns")
    )
    expect_identical(module_env$ruvQcUiBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_true(is.function(module_env$ruvQcUiBuilderCalls[[1]]$ns))
    expect_length(module_env$ruvQcHelperCalls, 1)
    expect_identical(
        module_env$ruvQcHelperCalls[[1]]$renderRuvQcUi,
        module_env$ruvQcUiRendererStub
    )
    expect_length(module_env$correlationFilterSummaryBuilderCalls, 1)
    expect_setequal(
        names(module_env$correlationFilterSummaryBuilderCalls[[1]]),
        "normData"
    )
    expect_identical(module_env$correlationFilterSummaryBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_length(module_env$correlationFilterSummaryHelperCalls, 1)
    expect_identical(
        module_env$correlationFilterSummaryHelperCalls[[1]]$renderCorrelationFilterSummary,
        module_env$correlationFilterSummaryRendererStub
    )
    expect_length(module_env$finalQcPlotBuilderCalls, 1)
    expect_setequal(
        names(module_env$finalQcPlotBuilderCalls[[1]]),
        c("normData", "getPlotAestheticsFn")
    )
    expect_identical(module_env$finalQcPlotBuilderCalls[[1]]$normData, module_env$reactiveStateStub)
    expect_identical(
        module_env$finalQcPlotBuilderCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_length(module_env$finalQcPlotHelperCalls, 1)
    expect_identical(
        module_env$finalQcPlotHelperCalls[[1]]$renderFinalQcPlot,
        module_env$finalQcPlotRendererStub
    )
    expect_length(module_env$designDrivenChoiceObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$designDrivenChoiceObserverHelperCalls[[1]]),
        c("session", "workflowData")
    )
    expect_false(is.null(module_env$designDrivenChoiceObserverHelperCalls[[1]]$session))
    expect_identical(module_env$designDrivenChoiceObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
})

test_that("createLipidNormReactiveState keeps the normalization state defaults stable", {
    args <- NULL
    reactive_state <- createLipidNormReactiveState(
        reactiveValuesFn = function(...) {
            args <<- list(...)
            args
        }
    )

    expect_identical(reactive_state, args)
    expect_named(
        reactive_state,
        c(
            "assay_names",
            "itsd_selections",
            "pre_norm_qc_generated",
            "normalization_complete",
            "ruv_complete",
            "correlation_filtering_complete",
            "qc_plot_paths",
            "post_filter_obj",
            "post_itsd_obj",
            "post_log2_obj",
            "post_norm_obj",
            "ruv_corrected_obj",
            "correlation_filtered_obj",
            "ruv_optimization_results",
            "correlation_results",
            "plot_refresh_trigger",
            "normalization_log"
        ),
        ignore.order = FALSE
    )
    expect_null(reactive_state$assay_names)
    expect_identical(reactive_state$itsd_selections, list())
    expect_false(reactive_state$pre_norm_qc_generated)
    expect_false(reactive_state$normalization_complete)
    expect_false(reactive_state$ruv_complete)
    expect_false(reactive_state$correlation_filtering_complete)
    expect_identical(reactive_state$qc_plot_paths, list())
    expect_null(reactive_state$post_filter_obj)
    expect_null(reactive_state$post_itsd_obj)
    expect_null(reactive_state$post_log2_obj)
    expect_null(reactive_state$post_norm_obj)
    expect_null(reactive_state$ruv_corrected_obj)
    expect_null(reactive_state$correlation_filtered_obj)
    expect_identical(reactive_state$ruv_optimization_results, list())
    expect_identical(reactive_state$correlation_results, list())
    expect_identical(reactive_state$plot_refresh_trigger, 0)
    expect_identical(reactive_state$normalization_log, character(0))
})

test_that("createLipidNormStartupRuntime keeps the startup builder bundle stable", {
    add_log_stub <- function(message) paste("log", message)
    get_plot_aesthetics_stub <- function() list(color_var = "condition")
    generate_composite_stub <- function(...) list(...)
    qc_image_renderer_stub <- function(...) "qc renderer"
    norm_log_renderer_stub <- function() "norm log"
    itsd_selection_renderer_stub <- function() "itsd selection"
    ruv_qc_ui_renderer_stub <- function() "ruv qc"
    assay_label_renderer_stub <- function(...) "assay label"
    correlation_summary_renderer_stub <- function() "correlation summary"
    final_qc_plot_renderer_stub <- function() "final qc"

    add_log_args <- NULL
    get_plot_aesthetics_args <- NULL
    composite_generator_args <- NULL
    qc_image_renderer_args <- NULL
    norm_log_renderer_args <- NULL
    itsd_selection_renderer_args <- NULL
    ruv_qc_ui_renderer_args <- NULL
    assay_label_renderer_args <- NULL
    correlation_summary_renderer_args <- NULL
    final_qc_plot_renderer_args <- NULL

    runtime <- createLipidNormStartupRuntime(
        input = "input-stub",
        experimentPaths = list(lipid_qc_dir = "qc-dir"),
        normData = "norm-state",
        ns = function(id) paste("ns", id),
        buildAddLogFn = function(normData, ...) {
            add_log_args <<- list(normData = normData)
            add_log_stub
        },
        buildPlotAestheticsGetterFn = function(input, ...) {
            get_plot_aesthetics_args <<- list(input = input)
            get_plot_aesthetics_stub
        },
        buildCompositeFromFilesGeneratorFn = function(...) {
            composite_generator_args <<- list(...)
            generate_composite_stub
        },
        buildQcImageRendererFn = function(normData, experimentPaths, ...) {
            qc_image_renderer_args <<- list(
                normData = normData,
                experimentPaths = experimentPaths
            )
            qc_image_renderer_stub
        },
        buildLogRendererFn = function(normData, ...) {
            norm_log_renderer_args <<- list(normData = normData)
            norm_log_renderer_stub
        },
        buildItsdSelectionUiRendererFn = function(normData, ns, ...) {
            itsd_selection_renderer_args <<- list(normData = normData, ns = ns)
            itsd_selection_renderer_stub
        },
        buildRuvQcUiRendererFn = function(normData, ns, ...) {
            ruv_qc_ui_renderer_args <<- list(normData = normData, ns = ns)
            ruv_qc_ui_renderer_stub
        },
        buildAssayLabelRendererFn = function(normData, ...) {
            assay_label_renderer_args <<- list(normData = normData)
            assay_label_renderer_stub
        },
        buildCorrelationFilterSummaryRendererFn = function(normData, ...) {
            correlation_summary_renderer_args <<- list(normData = normData)
            correlation_summary_renderer_stub
        },
        buildFinalQcPlotRendererFn = function(normData, getPlotAestheticsFn, ...) {
            final_qc_plot_renderer_args <<- list(
                normData = normData,
                getPlotAestheticsFn = getPlotAestheticsFn
            )
            final_qc_plot_renderer_stub
        }
    )

    expect_named(
        runtime,
        c(
            "addLog",
            "getPlotAesthetics",
            "generateCompositeFromFiles",
            "renderQcImageForAssay",
            "renderNormLog",
            "renderItsdSelectionUi",
            "renderRuvQcUi",
            "renderAssayLabel",
            "renderCorrelationFilterSummary",
            "renderFinalQcPlot"
        ),
        ignore.order = FALSE
    )
    expect_identical(runtime$addLog, add_log_stub)
    expect_identical(runtime$getPlotAesthetics, get_plot_aesthetics_stub)
    expect_identical(runtime$generateCompositeFromFiles, generate_composite_stub)
    expect_identical(runtime$renderQcImageForAssay, qc_image_renderer_stub)
    expect_identical(runtime$renderNormLog, norm_log_renderer_stub)
    expect_identical(runtime$renderItsdSelectionUi, itsd_selection_renderer_stub)
    expect_identical(runtime$renderRuvQcUi, ruv_qc_ui_renderer_stub)
    expect_identical(runtime$renderAssayLabel, assay_label_renderer_stub)
    expect_identical(runtime$renderCorrelationFilterSummary, correlation_summary_renderer_stub)
    expect_identical(runtime$renderFinalQcPlot, final_qc_plot_renderer_stub)
    expect_identical(add_log_args$normData, "norm-state")
    expect_identical(get_plot_aesthetics_args$input, "input-stub")
    expect_identical(composite_generator_args, list())
    expect_identical(qc_image_renderer_args$normData, "norm-state")
    expect_identical(qc_image_renderer_args$experimentPaths$lipid_qc_dir, "qc-dir")
    expect_identical(norm_log_renderer_args$normData, "norm-state")
    expect_identical(itsd_selection_renderer_args$normData, "norm-state")
    expect_identical(itsd_selection_renderer_args$ns("probe"), "ns probe")
    expect_identical(ruv_qc_ui_renderer_args$normData, "norm-state")
    expect_identical(ruv_qc_ui_renderer_args$ns("probe"), "ns probe")
    expect_identical(assay_label_renderer_args$normData, "norm-state")
    expect_identical(correlation_summary_renderer_args$normData, "norm-state")
    expect_identical(final_qc_plot_renderer_args$normData, "norm-state")
    expect_identical(
        final_qc_plot_renderer_args$getPlotAestheticsFn,
        get_plot_aesthetics_stub
    )
})

test_that("registerLipidNormPrimaryStartupOutputs keeps the startup output handoff stable", {
    output <- new.env(parent = emptyenv())
    startup_runtime <- list(
        renderNormLog = function() "norm log",
        renderItsdSelectionUi = function() "itsd selection",
        renderRuvQcUi = function() "ruv qc",
        renderQcImageForAssay = function(...) "qc image",
        renderAssayLabel = function(...) "assay label"
    )
    log_calls <- list()
    itsd_calls <- list()
    ruv_calls <- list()
    static_qc_calls <- list()
    assay_label_calls <- list()

    result <- registerLipidNormPrimaryStartupOutputs(
        output = output,
        startupRuntime = startup_runtime,
        registerLogOutputFn = function(output, renderNormLog) {
            log_calls[[length(log_calls) + 1]] <<- list(
                output = output,
                renderNormLog = renderNormLog
            )
            invisible(output)
        },
        registerItsdSelectionOutputFn = function(output, renderItsdSelectionUi) {
            itsd_calls[[length(itsd_calls) + 1]] <<- list(
                output = output,
                renderItsdSelectionUi = renderItsdSelectionUi
            )
            invisible(output)
        },
        registerRuvQcOutputFn = function(output, renderRuvQcUi) {
            ruv_calls[[length(ruv_calls) + 1]] <<- list(
                output = output,
                renderRuvQcUi = renderRuvQcUi
            )
            invisible(output)
        },
        registerStaticQcImageOutputsFn = function(output, renderQcImageForAssay) {
            static_qc_calls[[length(static_qc_calls) + 1]] <<- list(
                output = output,
                renderQcImageForAssay = renderQcImageForAssay
            )
            invisible(output)
        },
        registerAssayLabelOutputsFn = function(output, renderAssayLabel) {
            assay_label_calls[[length(assay_label_calls) + 1]] <<- list(
                output = output,
                renderAssayLabel = renderAssayLabel
            )
            invisible(output)
        }
    )

    expect_identical(result, output)
    expect_identical(log_calls[[1]]$output, output)
    expect_identical(log_calls[[1]]$renderNormLog, startup_runtime$renderNormLog)
    expect_identical(itsd_calls[[1]]$output, output)
    expect_identical(itsd_calls[[1]]$renderItsdSelectionUi, startup_runtime$renderItsdSelectionUi)
    expect_identical(ruv_calls[[1]]$output, output)
    expect_identical(ruv_calls[[1]]$renderRuvQcUi, startup_runtime$renderRuvQcUi)
    expect_identical(static_qc_calls[[1]]$output, output)
    expect_identical(static_qc_calls[[1]]$renderQcImageForAssay, startup_runtime$renderQcImageForAssay)
    expect_identical(assay_label_calls[[1]]$output, output)
    expect_identical(assay_label_calls[[1]]$renderAssayLabel, startup_runtime$renderAssayLabel)
})

test_that("registerLipidNormItsdSelectionRuntime keeps the ITSD registration handoff stable", {
    input <- new.env(parent = emptyenv())
    output <- new.env(parent = emptyenv())
    workflow_data <- list(state_manager = "state-manager")
    norm_data <- new.env(parent = emptyenv())
    itsd_table_calls <- list()
    itsd_tracking_calls <- list()

    result <- registerLipidNormItsdSelectionRuntime(
        input = input,
        output = output,
        workflowData = workflow_data,
        normData = norm_data,
        registerItsdTableOutputsFn = function(output, workflowData, normData) {
            itsd_table_calls[[length(itsd_table_calls) + 1]] <<- list(
                output = output,
                workflowData = workflowData,
                normData = normData
            )
            invisible(output)
        },
        registerItsdSelectionTrackingFn = function(input, normData) {
            itsd_tracking_calls[[length(itsd_tracking_calls) + 1]] <<- list(
                input = input,
                normData = normData
            )
            invisible(normData)
        }
    )

    expect_identical(result, norm_data)
    expect_identical(itsd_table_calls[[1]]$output, output)
    expect_identical(itsd_table_calls[[1]]$workflowData, workflow_data)
    expect_identical(itsd_table_calls[[1]]$normData, norm_data)
    expect_identical(itsd_tracking_calls[[1]]$input, input)
    expect_identical(itsd_tracking_calls[[1]]$normData, norm_data)
})

test_that("registerLipidNormPostNormalizationOutputs keeps the post-normalization output handoff stable", {
    output <- new.env(parent = emptyenv())
    norm_data <- new.env(parent = emptyenv())
    startup_runtime <- list(
        renderCorrelationFilterSummary = function() "correlation summary",
        renderFinalQcPlot = function() "final qc plot"
    )
    ruv_calls <- list()
    correlation_summary_calls <- list()
    final_qc_calls <- list()

    result <- registerLipidNormPostNormalizationOutputs(
        output = output,
        normData = norm_data,
        startupRuntime = startup_runtime,
        registerRuvCancorOutputsFn = function(output, normData) {
            ruv_calls[[length(ruv_calls) + 1]] <<- list(
                output = output,
                normData = normData
            )
            invisible(output)
        },
        registerCorrelationFilterSummaryOutputFn = function(output, renderCorrelationFilterSummary) {
            correlation_summary_calls[[length(correlation_summary_calls) + 1]] <<- list(
                output = output,
                renderCorrelationFilterSummary = renderCorrelationFilterSummary
            )
            invisible(output)
        },
        registerFinalQcPlotOutputFn = function(output, renderFinalQcPlot) {
            final_qc_calls[[length(final_qc_calls) + 1]] <<- list(
                output = output,
                renderFinalQcPlot = renderFinalQcPlot
            )
            invisible(output)
        }
    )

    expect_identical(result, output)
    expect_identical(ruv_calls[[1]]$output, output)
    expect_identical(ruv_calls[[1]]$normData, norm_data)
    expect_identical(correlation_summary_calls[[1]]$output, output)
    expect_identical(
        correlation_summary_calls[[1]]$renderCorrelationFilterSummary,
        startup_runtime$renderCorrelationFilterSummary
    )
    expect_identical(final_qc_calls[[1]]$output, output)
    expect_identical(
        final_qc_calls[[1]]$renderFinalQcPlot,
        startup_runtime$renderFinalQcPlot
    )
})

test_that("registerLipidNormStartupObserverRuntime keeps the startup observer handoff stable", {
    session <- structure(list(id = "norm-session"), class = "mock-session")
    workflow_data <- list(state_manager = list(getState = function() NULL))
    experiment_paths <- list(lipid_qc_dir = tempdir())
    norm_data <- new.env(parent = emptyenv())
    add_log <- function(message) paste("log", message)
    get_plot_aesthetics <- function() list(color_var = "group")
    selected_tab <- function() "norm"
    assay_init_calls <- list()
    selected_tab_calls <- list()
    design_choice_calls <- list()

    result <- registerLipidNormStartupObserverRuntime(
        session = session,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        normData = norm_data,
        addLog = add_log,
        getPlotAestheticsFn = get_plot_aesthetics,
        selectedTab = selected_tab,
        registerAssayNameInitializationObserverFn = function(workflowData, normData) {
            assay_init_calls[[length(assay_init_calls) + 1]] <<- list(
                workflowData = workflowData,
                normData = normData
            )
            invisible(normData)
        },
        registerSelectedTabPreNormalizationObserverFn = function(selectedTab, workflowData, experimentPaths, normData, addLog, getPlotAestheticsFn) {
            selected_tab_calls[[length(selected_tab_calls) + 1]] <<- list(
                selectedTab = selectedTab,
                workflowData = workflowData,
                experimentPaths = experimentPaths,
                normData = normData,
                addLog = addLog,
                getPlotAestheticsFn = getPlotAestheticsFn
            )
            invisible(selectedTab)
        },
        registerDesignDrivenChoiceObserverFn = function(session, workflowData) {
            design_choice_calls[[length(design_choice_calls) + 1]] <<- list(
                session = session,
                workflowData = workflowData
            )
            invisible(session)
        }
    )

    expect_identical(result, session)
    expect_identical(assay_init_calls[[1]]$workflowData, workflow_data)
    expect_identical(assay_init_calls[[1]]$normData, norm_data)
    expect_identical(selected_tab_calls[[1]]$selectedTab, selected_tab)
    expect_identical(selected_tab_calls[[1]]$workflowData, workflow_data)
    expect_identical(selected_tab_calls[[1]]$experimentPaths, experiment_paths)
    expect_identical(selected_tab_calls[[1]]$normData, norm_data)
    expect_identical(selected_tab_calls[[1]]$addLog, add_log)
    expect_identical(selected_tab_calls[[1]]$getPlotAestheticsFn, get_plot_aesthetics)
    expect_identical(design_choice_calls[[1]]$session, session)
    expect_identical(design_choice_calls[[1]]$workflowData, workflow_data)
})

test_that("registerLipidNormServerRuntime keeps the server-runtime handoff stable", {
    input <- new.env(parent = emptyenv())
    output <- new.env(parent = emptyenv())
    session <- structure(list(id = "norm-session"), class = "mock-session")
    workflow_data <- list(state_manager = list(getState = function() NULL))
    experiment_paths <- list(lipid_qc_dir = tempdir())
    norm_data <- new.env(parent = emptyenv())
    startup_runtime <- list(
        addLog = function(message) paste("log", message),
        getPlotAesthetics = function() list(color_var = "group", shape_var = "batch"),
        generateCompositeFromFiles = function() "composite",
        renderNormLog = function() "norm log renderer",
        renderItsdSelectionUi = function() "itsd selection renderer",
        renderRuvQcUi = function() "ruv qc renderer",
        renderQcImageForAssay = function(...) "qc image renderer",
        renderAssayLabel = function(...) "assay label renderer",
        renderCorrelationFilterSummary = function() "correlation summary renderer",
        renderFinalQcPlot = function() "final qc plot renderer"
    )
    add_log <- function(message) paste("add log", message)
    get_plot_aesthetics <- function() list(color_var = "group", shape_var = "batch")
    generate_composite <- function() "composite"
    selected_tab <- function() "norm"
    startup_observer_calls <- list()
    primary_output_calls <- list()
    itsd_runtime_calls <- list()
    run_observer_calls <- list()
    reset_observer_calls <- list()
    post_output_calls <- list()
    apply_observer_calls <- list()
    skip_observer_calls <- list()
    export_observer_calls <- list()

    result <- registerLipidNormServerRuntime(
        input = input,
        output = output,
        session = session,
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        omicType = "lipidomics",
        experimentLabel = "Lipidomics",
        normData = norm_data,
        startupRuntime = startup_runtime,
        addLog = add_log,
        getPlotAestheticsFn = get_plot_aesthetics,
        generateCompositeFromFilesFn = generate_composite,
        selectedTab = selected_tab,
        registerStartupObserverRuntimeFn = function(...) {
            startup_observer_calls[[length(startup_observer_calls) + 1]] <<- list(...)
            invisible(session)
        },
        registerPrimaryStartupOutputsFn = function(...) {
            primary_output_calls[[length(primary_output_calls) + 1]] <<- list(...)
            invisible(output)
        },
        registerItsdSelectionRuntimeFn = function(...) {
            itsd_runtime_calls[[length(itsd_runtime_calls) + 1]] <<- list(...)
            invisible(norm_data)
        },
        registerRunNormalizationObserverFn = function(...) {
            run_observer_calls[[length(run_observer_calls) + 1]] <<- list(...)
            invisible(input)
        },
        registerResetNormalizationObserverFn = function(...) {
            reset_observer_calls[[length(reset_observer_calls) + 1]] <<- list(...)
            invisible(input)
        },
        registerPostNormalizationOutputsFn = function(...) {
            post_output_calls[[length(post_output_calls) + 1]] <<- list(...)
            invisible(output)
        },
        registerApplyCorrelationFilterObserverFn = function(...) {
            apply_observer_calls[[length(apply_observer_calls) + 1]] <<- list(...)
            invisible(input)
        },
        registerSkipCorrelationFilterObserverFn = function(...) {
            skip_observer_calls[[length(skip_observer_calls) + 1]] <<- list(...)
            invisible(input)
        },
        registerExportSessionObserverFn = function(...) {
            export_observer_calls[[length(export_observer_calls) + 1]] <<- list(...)
            invisible(input)
        }
    )

    expect_identical(result, input)
    expect_identical(startup_observer_calls[[1]]$session, session)
    expect_identical(startup_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(startup_observer_calls[[1]]$experimentPaths, experiment_paths)
    expect_identical(startup_observer_calls[[1]]$normData, norm_data)
    expect_identical(startup_observer_calls[[1]]$addLog, add_log)
    expect_identical(startup_observer_calls[[1]]$getPlotAestheticsFn, get_plot_aesthetics)
    expect_identical(startup_observer_calls[[1]]$selectedTab, selected_tab)
    expect_identical(primary_output_calls[[1]]$output, output)
    expect_identical(primary_output_calls[[1]]$startupRuntime, startup_runtime)
    expect_identical(itsd_runtime_calls[[1]]$input, input)
    expect_identical(itsd_runtime_calls[[1]]$output, output)
    expect_identical(itsd_runtime_calls[[1]]$workflowData, workflow_data)
    expect_identical(itsd_runtime_calls[[1]]$normData, norm_data)
    expect_identical(run_observer_calls[[1]]$input, input)
    expect_identical(run_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(run_observer_calls[[1]]$experimentPaths, experiment_paths)
    expect_identical(run_observer_calls[[1]]$omicType, "lipidomics")
    expect_identical(run_observer_calls[[1]]$normData, norm_data)
    expect_identical(run_observer_calls[[1]]$addLog, add_log)
    expect_identical(run_observer_calls[[1]]$getPlotAestheticsFn, get_plot_aesthetics)
    expect_identical(run_observer_calls[[1]]$generateCompositeFromFilesFn, generate_composite)
    expect_identical(reset_observer_calls[[1]]$input, input)
    expect_identical(reset_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(reset_observer_calls[[1]]$normData, norm_data)
    expect_identical(reset_observer_calls[[1]]$addLog, add_log)
    expect_identical(post_output_calls[[1]]$output, output)
    expect_identical(post_output_calls[[1]]$normData, norm_data)
    expect_identical(post_output_calls[[1]]$startupRuntime, startup_runtime)
    expect_identical(apply_observer_calls[[1]]$input, input)
    expect_identical(apply_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(apply_observer_calls[[1]]$normData, norm_data)
    expect_identical(apply_observer_calls[[1]]$addLog, add_log)
    expect_identical(skip_observer_calls[[1]]$input, input)
    expect_identical(skip_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(skip_observer_calls[[1]]$normData, norm_data)
    expect_identical(skip_observer_calls[[1]]$addLog, add_log)
    expect_identical(export_observer_calls[[1]]$input, input)
    expect_identical(export_observer_calls[[1]]$workflowData, workflow_data)
    expect_identical(export_observer_calls[[1]]$experimentPaths, experiment_paths)
    expect_identical(export_observer_calls[[1]]$experimentLabel, "Lipidomics")
    expect_identical(export_observer_calls[[1]]$normData, norm_data)
    expect_identical(export_observer_calls[[1]]$addLog, add_log)
})

test_that("runLipidNormModuleServerShell keeps the wrapper shell handoff stable", {
    input <- new.env(parent = emptyenv())
    output <- new.env(parent = emptyenv())
    session <- list(ns = function(id) paste("ns", id))
    reactive_state <- new.env(parent = emptyenv())
    workflow_data <- list(state_manager = list(getState = function() NULL))
    experiment_paths <- list(lipid_qc_dir = tempdir())
    startup_runtime <- list(
        addLog = function(message) paste("log", message),
        getPlotAesthetics = function() list(color_var = "group", shape_var = "batch"),
        generateCompositeFromFiles = function() "composite"
    )
    log_messages <- character()
    startup_runtime_args <- NULL
    server_runtime_args <- NULL

    result <- runLipidNormModuleServerShell(
        input = input,
        output = output,
        session = session,
        id = "norm",
        workflowData = workflow_data,
        experimentPaths = experiment_paths,
        omicType = "lipidomics",
        experimentLabel = "Lipidomics",
        selectedTab = function() "norm",
        logInfoFn = function(message) {
            log_messages <<- c(log_messages, message)
            invisible(NULL)
        },
        createReactiveStateFn = function() reactive_state,
        createStartupRuntimeFn = function(...) {
            startup_runtime_args <<- list(...)
            startup_runtime
        },
        registerServerRuntimeFn = function(...) {
            server_runtime_args <<- list(...)
            invisible(input)
        }
    )

    expect_identical(result, reactive_state)
    expect_identical(
        log_messages,
        c(
            "=== METABOLOMICS NORMALIZATION MODULE STARTED ===",
            "Module ID: norm"
        )
    )
    expect_setequal(
        names(startup_runtime_args),
        c("input", "experimentPaths", "normData", "ns")
    )
    expect_identical(startup_runtime_args$input, input)
    expect_identical(startup_runtime_args$experimentPaths, experiment_paths)
    expect_identical(startup_runtime_args$normData, reactive_state)
    expect_identical(startup_runtime_args$ns("probe"), "ns probe")
    expect_setequal(
        names(server_runtime_args),
        c(
            "input", "output", "session", "workflowData", "experimentPaths",
            "omicType", "experimentLabel", "normData", "startupRuntime",
            "addLog", "getPlotAestheticsFn", "generateCompositeFromFilesFn",
            "selectedTab"
        )
    )
    expect_identical(server_runtime_args$input, input)
    expect_identical(server_runtime_args$output, output)
    expect_identical(server_runtime_args$session, session)
    expect_identical(server_runtime_args$workflowData, workflow_data)
    expect_identical(server_runtime_args$experimentPaths, experiment_paths)
    expect_identical(server_runtime_args$omicType, "lipidomics")
    expect_identical(server_runtime_args$experimentLabel, "Lipidomics")
    expect_identical(server_runtime_args$normData, reactive_state)
    expect_identical(server_runtime_args$startupRuntime, startup_runtime)
    expect_identical(server_runtime_args$addLog, startup_runtime$addLog)
    expect_identical(
        server_runtime_args$getPlotAestheticsFn,
        startup_runtime$getPlotAesthetics
    )
    expect_identical(
        server_runtime_args$generateCompositeFromFilesFn,
        startup_runtime$generateCompositeFromFiles
    )
    expect_identical(server_runtime_args$selectedTab(), "norm")
})

test_that("runLipidNormModuleServerEntryShell keeps the public entry-shell handoff stable", {
    module_server_id <- NULL
    captured_module <- NULL
    module_server_shell_args <- NULL

    result <- runLipidNormModuleServerEntryShell(
        id = "norm",
        workflowData = list(state_manager = list(getState = function() NULL)),
        experimentPaths = list(lipid_qc_dir = tempdir()),
        omicType = "lipidomics",
        experimentLabel = "Lipidomics",
        selectedTab = function() "norm",
        moduleServerFn = function(id, module) {
            module_server_id <<- id
            captured_module <<- module
            invisible("module-server-stub")
        },
        runModuleServerShellFn = function(...) {
            module_server_shell_args <<- list(...)
            invisible("module-shell-stub")
        }
    )

    expect_identical(result, "module-server-stub")
    expect_identical(module_server_id, "norm")
    expect_true(is.function(captured_module))

    captured_module(
        input = structure(list(run_normalization = 0), class = "reactivevalues"),
        output = new.env(parent = emptyenv()),
        session = list(ns = identity)
    )

    expect_setequal(
        names(module_server_shell_args),
        c(
            "input", "output", "session", "id", "workflowData",
            "experimentPaths", "omicType", "experimentLabel", "selectedTab"
        )
    )
    expect_identical(module_server_shell_args$id, "norm")
    expect_identical(module_server_shell_args$omicType, "lipidomics")
    expect_identical(module_server_shell_args$experimentLabel, "Lipidomics")
    expect_identical(module_server_shell_args$selectedTab(), "norm")
})

test_that("runLipidNormModuleServerPublicWrapper keeps the breadcrumb public-wrapper handoff stable", {
    module_server_entry_shell_args <- NULL

    result <- runLipidNormModuleServerPublicWrapper(
        id = "norm",
        workflow_data = list(state_manager = list(getState = function() NULL)),
        experiment_paths = list(lipid_qc_dir = tempdir()),
        omic_type = "lipidomics",
        experiment_label = "Lipidomics",
        selected_tab = function() "norm",
        runModuleServerEntryShellFn = function(...) {
            module_server_entry_shell_args <<- list(...)
            invisible("entry-shell-stub")
        }
    )

    expect_identical(result, "entry-shell-stub")
    expect_setequal(
        names(module_server_entry_shell_args),
        c(
            "id", "workflowData", "experimentPaths",
            "omicType", "experimentLabel", "selectedTab"
        )
    )
    expect_identical(module_server_entry_shell_args$id, "norm")
    expect_identical(module_server_entry_shell_args$omicType, "lipidomics")
    expect_identical(module_server_entry_shell_args$experimentLabel, "Lipidomics")
    expect_identical(module_server_entry_shell_args$selectedTab(), "norm")
})

test_that("mod_lipid_norm_server routes the breadcrumb public wrapper through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubModuleServerPublicWrapper = TRUE,
        selectedTabValue = "norm"
    )

    expect_length(module_env$moduleServerPublicWrapperHelperCalls, 1)
    expect_setequal(
        names(module_env$moduleServerPublicWrapperHelperCalls[[1]]),
        c(
            "id", "workflow_data", "experiment_paths",
            "omic_type", "experiment_label", "selected_tab"
        )
    )
    expect_identical(module_env$moduleServerPublicWrapperHelperCalls[[1]]$id, "norm")
    expect_true(inherits(
        module_env$moduleServerPublicWrapperHelperCalls[[1]]$workflow_data$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$moduleServerPublicWrapperHelperCalls[[1]]$experiment_paths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(module_env$moduleServerPublicWrapperHelperCalls[[1]]$omic_type, "lipidomics")
    expect_identical(module_env$moduleServerPublicWrapperHelperCalls[[1]]$experiment_label, "Lipidomics")
    expect_identical(module_env$moduleServerPublicWrapperHelperCalls[[1]]$selected_tab(), "norm")
    expect_length(module_env$moduleServerEntryShellHelperCalls, 0)
    expect_length(module_env$moduleServerShellHelperCalls, 0)
    expect_length(module_env$reactiveStateHelperCalls, 0)
    expect_length(module_env$startupRuntimeHelperCalls, 0)
    expect_length(module_env$serverRuntimeHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the public module wrapper through the entry-shell seam", {
    module_env <- loadLipidNormModuleHarness(
        stubModuleServerEntryShell = TRUE,
        selectedTabValue = "norm"
    )

    expect_length(module_env$moduleServerEntryShellHelperCalls, 1)
    expect_setequal(
        names(module_env$moduleServerEntryShellHelperCalls[[1]]),
        c(
            "id", "workflowData", "experimentPaths",
            "omicType", "experimentLabel", "selectedTab"
        )
    )
    expect_identical(module_env$moduleServerEntryShellHelperCalls[[1]]$id, "norm")
    expect_true(inherits(
        module_env$moduleServerEntryShellHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$moduleServerEntryShellHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(module_env$moduleServerEntryShellHelperCalls[[1]]$omicType, "lipidomics")
    expect_identical(module_env$moduleServerEntryShellHelperCalls[[1]]$experimentLabel, "Lipidomics")
    expect_identical(module_env$moduleServerEntryShellHelperCalls[[1]]$selectedTab(), "norm")
    expect_length(module_env$moduleServerShellHelperCalls, 0)
    expect_length(module_env$reactiveStateHelperCalls, 0)
    expect_length(module_env$startupRuntimeHelperCalls, 0)
    expect_length(module_env$serverRuntimeHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the module wrapper through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubModuleServerShell = TRUE,
        selectedTabValue = "norm"
    )

    expect_length(module_env$moduleServerShellHelperCalls, 1)
    expect_setequal(
        names(module_env$moduleServerShellHelperCalls[[1]]),
        c(
            "input", "output", "session", "id", "workflowData",
            "experimentPaths", "omicType", "experimentLabel", "selectedTab"
        )
    )
    expect_false(is.null(module_env$moduleServerShellHelperCalls[[1]]$input))
    expect_false(is.null(module_env$moduleServerShellHelperCalls[[1]]$output))
    expect_false(is.null(module_env$moduleServerShellHelperCalls[[1]]$session))
    expect_identical(module_env$moduleServerShellHelperCalls[[1]]$id, "norm")
    expect_true(inherits(
        module_env$moduleServerShellHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$moduleServerShellHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(module_env$moduleServerShellHelperCalls[[1]]$omicType, "lipidomics")
    expect_identical(module_env$moduleServerShellHelperCalls[[1]]$experimentLabel, "Lipidomics")
    expect_identical(module_env$moduleServerShellHelperCalls[[1]]$selectedTab(), "norm")
    expect_length(module_env$reactiveStateHelperCalls, 0)
    expect_length(module_env$startupRuntimeHelperCalls, 0)
    expect_length(module_env$serverRuntimeHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the startup runtime through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(stubStartupRuntime = TRUE)

    expect_length(module_env$startupRuntimeHelperCalls, 1)
    expect_setequal(
        names(module_env$startupRuntimeHelperCalls[[1]]),
        c("input", "experimentPaths", "normData", "ns")
    )
    expect_false(is.null(module_env$startupRuntimeHelperCalls[[1]]$input))
    expect_identical(
        module_env$startupRuntimeHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(
        module_env$startupRuntimeHelperCalls[[1]]$normData,
        module_env$reactiveStateStub
    )
    expect_identical(
        module_env$startupRuntimeHelperCalls[[1]]$ns("probe"),
        "norm-probe"
    )
    expect_length(module_env$normLogHelperCalls, 1)
    expect_identical(
        module_env$normLogHelperCalls[[1]]$renderNormLog,
        module_env$normLogRendererStub
    )
    expect_length(module_env$staticQcImageHelperCalls, 1)
    expect_identical(
        module_env$staticQcImageHelperCalls[[1]]$renderQcImageForAssay,
        module_env$qcImageRendererStub
    )
    expect_length(module_env$finalQcPlotHelperCalls, 1)
    expect_identical(
        module_env$finalQcPlotHelperCalls[[1]]$renderFinalQcPlot,
        module_env$finalQcPlotRendererStub
    )
    expect_length(module_env$runNormalizationObserverHelperCalls, 1)
    expect_identical(
        module_env$runNormalizationObserverHelperCalls[[1]]$addLog,
        module_env$addLogStub
    )
    expect_identical(
        module_env$runNormalizationObserverHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_identical(
        module_env$runNormalizationObserverHelperCalls[[1]]$generateCompositeFromFilesFn,
        module_env$generateCompositeFromFilesStub
    )
})

test_that("mod_lipid_norm_server routes the server runtime through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubStartupRuntime = TRUE,
        stubServerRuntime = TRUE,
        selectedTabValue = "norm"
    )

    expect_length(module_env$serverRuntimeHelperCalls, 1)
    expect_setequal(
        names(module_env$serverRuntimeHelperCalls[[1]]),
        c(
            "input", "output", "session", "workflowData", "experimentPaths",
            "omicType", "experimentLabel", "normData", "startupRuntime",
            "addLog", "getPlotAestheticsFn", "generateCompositeFromFilesFn",
            "selectedTab"
        )
    )
    expect_false(is.null(module_env$serverRuntimeHelperCalls[[1]]$input))
    expect_false(is.null(module_env$serverRuntimeHelperCalls[[1]]$output))
    expect_false(is.null(module_env$serverRuntimeHelperCalls[[1]]$session))
    expect_true(inherits(
        module_env$serverRuntimeHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(module_env$serverRuntimeHelperCalls[[1]]$omicType, "lipidomics")
    expect_identical(module_env$serverRuntimeHelperCalls[[1]]$experimentLabel, "Lipidomics")
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$normData,
        module_env$reactiveStateStub
    )
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$startupRuntime$renderNormLog,
        module_env$normLogRendererStub
    )
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$addLog,
        module_env$addLogStub
    )
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$generateCompositeFromFilesFn,
        module_env$generateCompositeFromFilesStub
    )
    expect_identical(
        module_env$serverRuntimeHelperCalls[[1]]$selectedTab(),
        "norm"
    )
    expect_length(module_env$startupObserverRuntimeHelperCalls, 0)
    expect_length(module_env$primaryStartupOutputsHelperCalls, 0)
    expect_length(module_env$itsdSelectionRuntimeHelperCalls, 0)
    expect_length(module_env$runNormalizationObserverHelperCalls, 0)
    expect_length(module_env$resetNormalizationObserverHelperCalls, 0)
    expect_length(module_env$postNormalizationOutputsHelperCalls, 0)
    expect_length(module_env$applyCorrelationFilterObserverHelperCalls, 0)
    expect_length(module_env$skipCorrelationFilterObserverHelperCalls, 0)
    expect_length(module_env$exportSessionObserverHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the primary startup outputs through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubStartupRuntime = TRUE,
        stubPrimaryStartupOutputs = TRUE
    )

    expect_length(module_env$primaryStartupOutputsHelperCalls, 1)
    expect_setequal(
        names(module_env$primaryStartupOutputsHelperCalls[[1]]),
        c("output", "startupRuntime")
    )
    expect_false(is.null(module_env$primaryStartupOutputsHelperCalls[[1]]$output))
    expect_identical(
        module_env$primaryStartupOutputsHelperCalls[[1]]$startupRuntime$renderNormLog,
        module_env$normLogRendererStub
    )
    expect_identical(
        module_env$primaryStartupOutputsHelperCalls[[1]]$startupRuntime$renderItsdSelectionUi,
        module_env$itsdSelectionRendererStub
    )
    expect_identical(
        module_env$primaryStartupOutputsHelperCalls[[1]]$startupRuntime$renderRuvQcUi,
        module_env$ruvQcUiRendererStub
    )
    expect_identical(
        module_env$primaryStartupOutputsHelperCalls[[1]]$startupRuntime$renderQcImageForAssay,
        module_env$qcImageRendererStub
    )
    expect_identical(
        module_env$primaryStartupOutputsHelperCalls[[1]]$startupRuntime$renderAssayLabel,
        module_env$assayLabelRendererStub
    )
})

test_that("mod_lipid_norm_server routes the ITSD registration runtime through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(stubItsdSelectionRuntime = TRUE)

    expect_length(module_env$itsdSelectionRuntimeHelperCalls, 1)
    expect_setequal(
        names(module_env$itsdSelectionRuntimeHelperCalls[[1]]),
        c("input", "output", "workflowData", "normData")
    )
    expect_false(is.null(module_env$itsdSelectionRuntimeHelperCalls[[1]]$input))
    expect_false(is.null(module_env$itsdSelectionRuntimeHelperCalls[[1]]$output))
    expect_identical(
        module_env$itsdSelectionRuntimeHelperCalls[[1]]$workflowData$state_manager$getState(),
        NULL
    )
    expect_false(is.null(module_env$itsdSelectionRuntimeHelperCalls[[1]]$normData))
    expect_length(module_env$itsdTableHelperCalls, 0)
    expect_length(module_env$itsdSelectionTrackingHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the post-normalization outputs through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubStartupRuntime = TRUE,
        stubPostNormalizationOutputs = TRUE
    )

    expect_length(module_env$postNormalizationOutputsHelperCalls, 1)
    expect_setequal(
        names(module_env$postNormalizationOutputsHelperCalls[[1]]),
        c("output", "normData", "startupRuntime")
    )
    expect_false(is.null(module_env$postNormalizationOutputsHelperCalls[[1]]$output))
    expect_identical(
        module_env$postNormalizationOutputsHelperCalls[[1]]$normData,
        module_env$reactiveStateStub
    )
    expect_identical(
        module_env$postNormalizationOutputsHelperCalls[[1]]$startupRuntime$renderCorrelationFilterSummary,
        module_env$correlationFilterSummaryRendererStub
    )
    expect_identical(
        module_env$postNormalizationOutputsHelperCalls[[1]]$startupRuntime$renderFinalQcPlot,
        module_env$finalQcPlotRendererStub
    )
    expect_length(module_env$ruvCancorOutputsHelperCalls, 0)
    expect_length(module_env$correlationFilterSummaryHelperCalls, 0)
    expect_length(module_env$finalQcPlotHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes the startup observer wiring through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        stubStartupRuntime = TRUE,
        stubStartupObserverRuntime = TRUE,
        selectedTabValue = "norm"
    )

    expect_length(module_env$startupObserverRuntimeHelperCalls, 1)
    expect_setequal(
        names(module_env$startupObserverRuntimeHelperCalls[[1]]),
        c("session", "workflowData", "experimentPaths", "normData", "addLog", "getPlotAestheticsFn", "selectedTab")
    )
    expect_false(is.null(module_env$startupObserverRuntimeHelperCalls[[1]]$session))
    expect_true(inherits(
        module_env$startupObserverRuntimeHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$startupObserverRuntimeHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(
        module_env$startupObserverRuntimeHelperCalls[[1]]$normData,
        module_env$reactiveStateStub
    )
    expect_identical(
        module_env$startupObserverRuntimeHelperCalls[[1]]$addLog,
        module_env$addLogStub
    )
    expect_identical(
        module_env$startupObserverRuntimeHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_identical(
        module_env$startupObserverRuntimeHelperCalls[[1]]$selectedTab(),
        "norm"
    )
    expect_length(module_env$assayNameInitializationHelperCalls, 0)
    expect_length(module_env$selectedTabPreNormalizationObserverHelperCalls, 0)
    expect_length(module_env$designDrivenChoiceObserverHelperCalls, 0)
})

test_that("mod_lipid_norm_server routes assay-name initialization through the top-level seam", {
    module_env <- loadLipidNormModuleHarness()

    expect_length(module_env$assayNameInitializationHelperCalls, 1)
    expect_setequal(
        names(module_env$assayNameInitializationHelperCalls[[1]]),
        c("workflowData", "normData")
    )
    expect_identical(module_env$assayNameInitializationHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_false(is.null(module_env$assayNameInitializationHelperCalls[[1]]$normData))
})

test_that("mod_lipid_norm_server routes the per-assay RUV cancor outputs through the top-level seam", {
    module_env <- loadLipidNormModuleHarness()

    expect_length(module_env$ruvCancorOutputsHelperCalls, 1)
    expect_setequal(
        names(module_env$ruvCancorOutputsHelperCalls[[1]]),
        c("output", "normData")
    )
    expect_false(is.null(module_env$ruvCancorOutputsHelperCalls[[1]]$normData))
})

test_that("mod_lipid_norm_server routes the per-assay ITSD tables through the top-level seam", {
    module_env <- loadLipidNormModuleHarness()

    expect_length(module_env$itsdTableHelperCalls, 1)
    expect_setequal(
        names(module_env$itsdTableHelperCalls[[1]]),
        c("output", "workflowData", "normData")
    )
    expect_identical(module_env$itsdTableHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_false(is.null(module_env$itsdTableHelperCalls[[1]]$normData))
})

test_that("mod_lipid_norm_server routes the ITSD-selection tracking through the top-level seam", {
    module_env <- loadLipidNormModuleHarness()

    expect_length(module_env$itsdSelectionTrackingHelperCalls, 1)
    expect_setequal(
        names(module_env$itsdSelectionTrackingHelperCalls[[1]]),
        c("input", "normData")
    )
    expect_false(is.null(module_env$itsdSelectionTrackingHelperCalls[[1]]$input))
    expect_false(is.null(module_env$itsdSelectionTrackingHelperCalls[[1]]$normData))
})

test_that("mod_lipid_norm_server routes the run-normalization observer through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(triggerRunNormalization = TRUE)

    expect_true("input$run_normalization" %in% module_env$capturedObserverEvents)
    expect_length(module_env$runNormalizationObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$runNormalizationObserverHelperCalls[[1]]),
        c("input", "workflowData", "experimentPaths", "omicType", "normData", "addLog", "getPlotAestheticsFn", "generateCompositeFromFilesFn")
    )
    expect_identical(module_env$runNormalizationObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_identical(module_env$runNormalizationObserverHelperCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(module_env$runNormalizationObserverHelperCalls[[1]]$omicType, "lipidomics")
    expect_identical(module_env$runNormalizationObserverHelperCalls[[1]]$addLog, module_env$addLogStub)
    expect_identical(
        module_env$runNormalizationObserverHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_identical(
        module_env$runNormalizationObserverHelperCalls[[1]]$generateCompositeFromFilesFn,
        module_env$generateCompositeFromFilesStub
    )
    expect_length(module_env$runNormalizationHelperCalls, 1)
    expect_setequal(
        names(module_env$runNormalizationHelperCalls[[1]]),
        c("input", "workflowData", "experimentPaths", "omicType", "normData", "addLog", "getPlotAestheticsFn", "generateCompositeFromFilesFn")
    )
    expect_identical(module_env$runNormalizationHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_identical(module_env$runNormalizationHelperCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(module_env$runNormalizationHelperCalls[[1]]$omicType, "lipidomics")
    expect_identical(module_env$runNormalizationHelperCalls[[1]]$addLog, module_env$addLogStub)
    expect_identical(
        module_env$runNormalizationHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_identical(
        module_env$runNormalizationHelperCalls[[1]]$generateCompositeFromFilesFn,
        module_env$generateCompositeFromFilesStub
    )
})

test_that("mod_lipid_norm_server routes the selected-tab pre-normalization trigger through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(
        triggerSelectedTab = TRUE,
        selectedTabValue = "norm"
    )

    expect_true("selected_tab()" %in% module_env$capturedObserverEvents)
    expect_length(module_env$selectedTabPreNormalizationObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]),
        c("selectedTab", "workflowData", "experimentPaths", "normData", "addLog", "getPlotAestheticsFn")
    )
    expect_identical(
        module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]$selectedTab(),
        "norm"
    )
    expect_true(inherits(
        module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(
        module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]$experimentPaths$lipid_qc_dir,
        tempdir()
    )
    expect_identical(
        module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]$addLog,
        module_env$addLogStub
    )
    expect_identical(
        module_env$selectedTabPreNormalizationObserverHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
    expect_length(module_env$selectedTabAutoTriggerHelperCalls, 1)
    expect_setequal(
        names(module_env$selectedTabAutoTriggerHelperCalls[[1]]),
        c("selectedTabValue", "workflowData", "experimentPaths", "normData", "addLog", "getPlotAestheticsFn")
    )
    expect_identical(module_env$selectedTabAutoTriggerHelperCalls[[1]]$selectedTabValue, "norm")
    expect_true(inherits(
        module_env$selectedTabAutoTriggerHelperCalls[[1]]$workflowData$state_manager$getState(),
        "LipidomicsAssayData"
    ))
    expect_identical(module_env$selectedTabAutoTriggerHelperCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(module_env$selectedTabAutoTriggerHelperCalls[[1]]$addLog, module_env$addLogStub)
    expect_identical(
        module_env$selectedTabAutoTriggerHelperCalls[[1]]$getPlotAestheticsFn,
        module_env$getPlotAestheticsStub
    )
})

test_that("mod_lipid_norm_server routes the export-session observer through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(triggerExportSession = TRUE)

    expect_true("input$export_session" %in% module_env$capturedObserverEvents)
    expect_length(module_env$exportSessionObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$exportSessionObserverHelperCalls[[1]]),
        c("input", "workflowData", "experimentPaths", "experimentLabel", "normData", "addLog")
    )
    expect_false(is.null(module_env$exportSessionObserverHelperCalls[[1]]$input))
    expect_identical(module_env$exportSessionObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_identical(module_env$exportSessionObserverHelperCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(module_env$exportSessionObserverHelperCalls[[1]]$experimentLabel, "Lipidomics")
    expect_true(is.function(module_env$exportSessionObserverHelperCalls[[1]]$addLog))
    expect_length(module_env$exportSessionHelperCalls, 1)
    expect_identical(
        module_env$exportSessionHelperCalls[[1]]$input,
        module_env$exportSessionObserverHelperCalls[[1]]$input
    )
    expect_identical(
        module_env$exportSessionHelperCalls[[1]]$workflowData,
        module_env$exportSessionObserverHelperCalls[[1]]$workflowData
    )
    expect_identical(
        module_env$exportSessionHelperCalls[[1]]$experimentPaths,
        module_env$exportSessionObserverHelperCalls[[1]]$experimentPaths
    )
    expect_identical(module_env$exportSessionHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_identical(module_env$exportSessionHelperCalls[[1]]$experimentPaths$lipid_qc_dir, tempdir())
    expect_identical(module_env$exportSessionHelperCalls[[1]]$experimentLabel, "Lipidomics")
    expect_true(is.function(module_env$exportSessionHelperCalls[[1]]$addLog))
})

test_that("mod_lipid_norm_server routes the correlation-filter observer through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(triggerApplyCorrelationFilter = TRUE)

    expect_true("input$apply_correlation_filter" %in% module_env$capturedObserverEvents)
    expect_length(module_env$applyCorrelationFilterObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$applyCorrelationFilterObserverHelperCalls[[1]]),
        c("input", "workflowData", "normData", "addLog")
    )
    expect_false(is.null(module_env$applyCorrelationFilterObserverHelperCalls[[1]]$input))
    expect_identical(module_env$applyCorrelationFilterObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_true(is.function(module_env$applyCorrelationFilterObserverHelperCalls[[1]]$addLog))
    expect_length(module_env$applyCorrelationFilterHelperCalls, 1)
    expect_identical(module_env$applyCorrelationFilterHelperCalls[[1]]$input, module_env$applyCorrelationFilterObserverHelperCalls[[1]]$input)
    expect_identical(module_env$applyCorrelationFilterHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_true(is.function(module_env$applyCorrelationFilterHelperCalls[[1]]$addLog))
})

test_that("mod_lipid_norm_server routes the skip-correlation observer through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(triggerSkipCorrelationFilter = TRUE)

    expect_true("input$skip_correlation_filter" %in% module_env$capturedObserverEvents)
    expect_length(module_env$skipCorrelationFilterObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$skipCorrelationFilterObserverHelperCalls[[1]]),
        c("input", "workflowData", "normData", "addLog")
    )
    expect_false(is.null(module_env$skipCorrelationFilterObserverHelperCalls[[1]]$input))
    expect_identical(module_env$skipCorrelationFilterObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_true(is.function(module_env$skipCorrelationFilterObserverHelperCalls[[1]]$addLog))
    expect_length(module_env$skipCorrelationFilterHelperCalls, 1)
    expect_identical(
        module_env$skipCorrelationFilterHelperCalls[[1]]$workflowData,
        module_env$skipCorrelationFilterObserverHelperCalls[[1]]$workflowData
    )
    expect_identical(
        module_env$skipCorrelationFilterHelperCalls[[1]]$normData,
        module_env$skipCorrelationFilterObserverHelperCalls[[1]]$normData
    )
    expect_identical(module_env$skipCorrelationFilterHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_true(is.function(module_env$skipCorrelationFilterHelperCalls[[1]]$addLog))
})

test_that("mod_lipid_norm_server routes the reset observer through the top-level seam", {
    module_env <- loadLipidNormModuleHarness(triggerResetNormalization = TRUE)

    expect_true("input$reset_normalization" %in% module_env$capturedObserverEvents)
    expect_length(module_env$resetNormalizationObserverHelperCalls, 1)
    expect_setequal(
        names(module_env$resetNormalizationObserverHelperCalls[[1]]),
        c("input", "workflowData", "normData", "addLog")
    )
    expect_identical(module_env$resetNormalizationObserverHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_false(is.null(module_env$resetNormalizationObserverHelperCalls[[1]]$input))
    expect_true(is.function(module_env$resetNormalizationObserverHelperCalls[[1]]$addLog))
    expect_length(module_env$resetNormalizationHelperCalls, 1)
    expect_identical(module_env$resetNormalizationHelperCalls[[1]]$workflowData$state_manager$getState(), NULL)
    expect_true(is.function(module_env$resetNormalizationHelperCalls[[1]]$addLog))
})
