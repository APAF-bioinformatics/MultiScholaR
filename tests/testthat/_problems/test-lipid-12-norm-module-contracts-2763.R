# Extracted from test-lipid-12-norm-module-contracts.R:2763

# prequel ----------------------------------------------------------------------
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
makeValidLipidNormContractS4 <- function(
    lipid_data,
    lipidIdColumn = "lipid_id",
    annotationIdColumn = "annotation"
) {
    lipidData <- lapply(lipid_data, function(assayData) {
        assayData <- as.data.frame(assayData)
        if (!lipidIdColumn %in% names(assayData)) {
            assayData[[lipidIdColumn]] <- paste0("L", seq_len(nrow(assayData)))
        }
        if (!"Sample_1" %in% names(assayData)) {
            assayData[["Sample_1"]] <- seq_len(nrow(assayData))
        }
        assayData
    })

    methods::new(
        "LipidomicsAssayData",
        lipid_data = lipidData,
        lipid_id_column = lipidIdColumn,
        annotation_id_column = annotationIdColumn,
        database_identifier_type = "Unknown",
        internal_standard_regex = NA_character_,
        design_matrix = data.frame(
            Sample_ID = "Sample_1",
            group = "WT",
            batch = "B1",
            stringsAsFactors = FALSE
        ),
        sample_id = "Sample_ID",
        group_id = "group",
        technical_replicate_id = NA_character_,
        args = list()
    )
}
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

# test -------------------------------------------------------------------------
saved_state <- NULL
notification <- NULL
log_message <- NULL
post_filter_obj <- structure(list(id = "post-filter"), class = "mockS4")
workflow_data <- new.env(parent = emptyenv())
workflow_data$config_list <- list(mode = "test")
omicsTestOwner(post_filter_obj, "lipid_normalized", workflow = workflow_data, omic_type = "lipidomics")
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
expect_identical(workflow_data$state_manager$current_state, "lipid_reset")
expect_identical(workflow_data$state_manager$getState("lipid_reset"), post_filter_obj)
expect_false(norm_data$normalization_complete)
expect_false(norm_data$ruv_complete)
expect_false(norm_data$correlation_filtering_complete)
expect_null(norm_data$post_norm_obj)
expect_null(norm_data$ruv_corrected_obj)
expect_null(norm_data$correlation_filtered_obj)
expect_identical(norm_data$ruv_optimization_results, list())
expect_identical(log_message, "Reset to pre-normalization state")
