library(testthat)
suppressPackageStartupMessages(library(shiny))

source(test_path("..", "..", "R", "mod_lipid_qc_intensity_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_qc_intensity.R"), local = environment())

test_that("registerLipidIntensityAssayResultsOutput keeps the empty-state render shell stable", {
    output_env <- new.env(parent = emptyenv())
    stats_calls <- 0L

    registerLipidIntensityAssayResultsOutput(
        output = output_env,
        filterStatsVal = function() {
            stats_calls <<- stats_calls + 1L
            NULL
        },
        ns = function(id) paste0("intensity-", id),
        renderUiFn = function(expr) force(expr)
    )

    expect_identical(stats_calls, 1L)
    expect_false(exists("assay_results_tabs", envir = output_env, inherits = FALSE))
})

test_that("registerLipidIntensityAssayResultsOutput keeps per-assay summary tabs stable", {
    output_env <- new.env(parent = emptyenv())
    stats_df <- data.frame(
        Assay = c("AssayA", "Assay 1/A"),
        Original = c(10, 8),
        Filtered = c(8, 5),
        Removed = c(2, 3),
        Percent_Retained = c(80, 62.5),
        stringsAsFactors = FALSE
    )

    registerLipidIntensityAssayResultsOutput(
        output = output_env,
        filterStatsVal = function() stats_df,
        ns = function(id) paste0("intensity-", id),
        renderUiFn = function(expr) htmltools::renderTags(expr)$html
    )

    expect_match(output_env$assay_results_tabs, "intensity-assay_stats_tabs", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "AssayA", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "Assay 1/A", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "Original lipids:", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "After filtering:", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "Percent retained:", fixed = TRUE)
    expect_match(output_env$assay_results_tabs, "62.5%", fixed = TRUE)
})

test_that("registerLipidIntensityFilterPlotOutput keeps the empty-state render shell stable", {
    output_env <- new.env(parent = emptyenv())
    requested_plot <- "unset"
    grid_draw_calls <- 0L
    print_calls <- 0L

    registerLipidIntensityFilterPlotOutput(
        output = output_env,
        filterPlotVal = function() NULL,
        renderPlotFn = function(expr) {
            try(force(expr), silent = TRUE)
            "rendered-plot"
        },
        reqFn = function(value) {
            requested_plot <<- value
            stop("plot required")
        },
        gridDrawFn = function(plot_obj) {
            grid_draw_calls <<- grid_draw_calls + 1L
            invisible(plot_obj)
        },
        printFn = function(plot_obj) {
            print_calls <<- print_calls + 1L
            invisible(plot_obj)
        }
    )

    expect_identical(output_env$filter_plot, "rendered-plot")
    expect_null(requested_plot)
    expect_identical(grid_draw_calls, 0L)
    expect_identical(print_calls, 0L)
})

test_that("registerLipidIntensityFilterPlotOutput keeps plot dispatch stable", {
    grob_output_env <- new.env(parent = emptyenv())
    grob_plot <- structure(list(id = "grid"), class = "grob")
    grob_draws <- list()
    grob_print_calls <- 0L

    registerLipidIntensityFilterPlotOutput(
        output = grob_output_env,
        filterPlotVal = function() grob_plot,
        renderPlotFn = function(expr) force(expr),
        reqFn = identity,
        gridDrawFn = function(plot_obj) {
            grob_draws <<- c(grob_draws, list(plot_obj))
            invisible(plot_obj)
        },
        printFn = function(plot_obj) {
            grob_print_calls <<- grob_print_calls + 1L
            invisible(plot_obj)
        }
    )

    ggplot_output_env <- new.env(parent = emptyenv())
    ggplot_plot <- structure(list(id = "gg"), class = "ggplot")
    ggplot_draw_calls <- 0L
    ggplot_prints <- list()

    registerLipidIntensityFilterPlotOutput(
        output = ggplot_output_env,
        filterPlotVal = function() ggplot_plot,
        renderPlotFn = function(expr) force(expr),
        reqFn = identity,
        gridDrawFn = function(plot_obj) {
            ggplot_draw_calls <<- ggplot_draw_calls + 1L
            invisible(plot_obj)
        },
        printFn = function(plot_obj) {
            ggplot_prints <<- c(ggplot_prints, list(plot_obj))
            invisible(plot_obj)
        }
    )

    expect_length(grob_draws, 1L)
    expect_identical(grob_draws[[1]], grob_plot)
    expect_identical(grob_print_calls, 0L)
    expect_identical(ggplot_draw_calls, 0L)
    expect_length(ggplot_prints, 1L)
    expect_identical(ggplot_prints[[1]], ggplot_plot)
})

ensureLipidIntensityWorkflowTestClass <- function() {
    if (!methods::isClass("LipidomicsAssayData")) {
        methods::setClass(
            "LipidomicsAssayData",
            slots = c(
                lipid_data = "list",
                lipid_id_column = "character"
            )
        )
    }
}

test_that("registerLipidIntensityApplyFilterObserver keeps the apply-filter workflow shell stable", {
    ensureLipidIntensityWorkflowTestClass()

    input <- list(
        apply_filter = 7L,
        intensity_cutoff_percentile = 10L,
        proportion_below_cutoff = 0.5
    )
    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = c("L1", "L1", "L2"),
                intensity = c(12, 9, 7),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                feature = c("F1", "F2"),
                intensity = c(4, 2),
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id"
    )
    filtered_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            AssayA = data.frame(
                lipid_id = "L1",
                intensity = 12,
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(
                feature = "F1",
                intensity = 4,
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id"
    )
    filter_stats_value <- NULL
    filter_plot_value <- "stale-plot"
    info_messages <- character()
    warn_messages <- character()
    notifications <- list()
    removed_notifications <- character()
    save_state_args <- NULL
    filtering_args <- NULL
    qc_plot_args <- NULL
    output_env <- new.env(parent = emptyenv())
    qc_plot <- structure(list(id = "qc-grid"), class = "grob")

    result <- registerLipidIntensityApplyFilterObserver(
        input = input,
        workflowData = list(
            state_manager = list(
                getState = function() current_s4,
                saveState = function(state_name, s4_data_object, config_object, description) {
                    save_state_args <<- list(
                        state_name = state_name,
                        s4_data_object = s4_data_object,
                        config_object = config_object,
                        description = description
                    )
                    invisible(NULL)
                }
            ),
            config_list = list(cache = "cfg")
        ),
        output = output_env,
        filterStatsVal = function(value) {
            filter_stats_value <<- value
        },
        filterPlotVal = function(value) {
            filter_plot_value <<- value
        },
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            expect_identical(event, 7L)
            force(handler)
            invisible(NULL)
        },
        reqFn = function(value) {
            if (is.null(value)) {
                stop("required value missing")
            }
            invisible(value)
        },
        showNotificationFn = function(message, type = NULL, id = NULL, duration = NULL) {
            notifications <<- c(notifications, list(list(
                message = message,
                type = type,
                id = id,
                duration = duration
            )))
        },
        removeNotificationFn = function(id) {
            removed_notifications <<- c(removed_notifications, id)
        },
        renderTextFn = function(expr) force(expr),
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarnFn = function(message) {
            warn_messages <<- c(warn_messages, message)
        },
        logErrorFn = function(...) {
            stop("logErrorFn should not be called on success")
        },
        lipidIntensityFilteringFn = function(theObject, lipids_intensity_cutoff_percentile, lipids_proportion_of_samples_below_cutoff) {
            filtering_args <<- list(
                theObject = theObject,
                percentile = lipids_intensity_cutoff_percentile,
                proportion = lipids_proportion_of_samples_below_cutoff
            )
            filtered_s4
        },
        updateLipidFilteringFn = function(theObject, step_name, omics_type, return_grid, overwrite) {
            qc_plot_args <<- list(
                theObject = theObject,
                step_name = step_name,
                omics_type = omics_type,
                return_grid = return_grid,
                overwrite = overwrite
            )
            qc_plot
        }
    )

    expect_identical(result, input)
    expect_identical(filtering_args$theObject, current_s4)
    expect_identical(filtering_args$percentile, 10L)
    expect_identical(filtering_args$proportion, 0.5)
    expect_identical(qc_plot_args$theObject, filtered_s4)
    expect_identical(qc_plot_args$step_name, "2_Intensity_Filtered")
    expect_identical(qc_plot_args$omics_type, "lipidomics")
    expect_identical(qc_plot_args$return_grid, TRUE)
    expect_identical(qc_plot_args$overwrite, TRUE)
    expect_equal(
        filter_stats_value,
        data.frame(
            Assay = c("AssayA", "AssayB"),
            Original = c(2, 2),
            Filtered = c(1, 1),
            Removed = c(1, 1),
            Percent_Retained = c(50, 50),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(filter_plot_value, qc_plot)
    expect_identical(save_state_args$state_name, "lipid_intensity_filtered")
    expect_identical(save_state_args$s4_data_object, filtered_s4)
    expect_identical(save_state_args$config_object, list(cache = "cfg"))
    expect_identical(
        save_state_args$description,
        "Applied lipid intensity filter (percentile: 10%, proportion: 0.50)"
    )
    expect_match(output_env$filter_results, "Lipid Intensity Filter Applied Successfully", fixed = TRUE)
    expect_match(output_env$filter_results, "AssayA: 2 -> 1", fixed = TRUE)
    expect_match(output_env$filter_results, "AssayB: 2 -> 1", fixed = TRUE)
    expect_match(output_env$filter_results, "Total: 4 -> 2 lipids (removed 2)", fixed = TRUE)
    expect_match(output_env$filter_results, "State saved as: 'lipid_intensity_filtered'", fixed = TRUE)
    expect_true(any(grepl(
        "Applying lipid intensity filter: percentile = 10 , proportion = 0.5",
        info_messages,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "Lipid intensity filter applied: removed 2 lipids",
        info_messages,
        fixed = TRUE
    )))
    expect_identical(warn_messages, character())
    expect_identical(removed_notifications, "lipid_intensity_filter_working")
    expect_identical(notifications[[1]], list(
        message = "Applying lipid intensity filter...",
        type = NULL,
        id = "lipid_intensity_filter_working",
        duration = NULL
    ))
    expect_identical(notifications[[2]], list(
        message = "Intensity filter applied: 2 lipids retained",
        type = "message",
        id = NULL,
        duration = NULL
    ))
})

test_that("registerLipidIntensityApplyFilterObserver keeps the apply-filter failure shell stable", {
    input <- list(
        apply_filter = 8L,
        intensity_cutoff_percentile = 10L,
        proportion_below_cutoff = 0.5
    )
    filter_stats_called <- FALSE
    filter_plot_called <- FALSE
    error_messages <- character()
    notifications <- list()
    removed_notifications <- character()

    result <- registerLipidIntensityApplyFilterObserver(
        input = input,
        workflowData = list(
            state_manager = list(
                getState = function() structure(list(marker = TRUE), class = "not_lipid_s4"),
                saveState = function(...) {
                    stop("saveState should not be called on failure")
                }
            ),
            config_list = list()
        ),
        output = new.env(parent = emptyenv()),
        filterStatsVal = function(...) {
            filter_stats_called <<- TRUE
        },
        filterPlotVal = function(...) {
            filter_plot_called <<- TRUE
        },
        omicType = "lipidomics",
        observeEventFn = function(event, handler) {
            expect_identical(event, 8L)
            force(handler)
            invisible(NULL)
        },
        reqFn = function(value) {
            if (is.null(value)) {
                stop("required value missing")
            }
            invisible(value)
        },
        showNotificationFn = function(message, type = NULL, id = NULL, duration = NULL) {
            notifications <<- c(notifications, list(list(
                message = message,
                type = type,
                id = id,
                duration = duration
            )))
        },
        removeNotificationFn = function(id) {
            removed_notifications <<- c(removed_notifications, id)
        },
        renderTextFn = function(...) {
            stop("renderTextFn should not be called on failure")
        },
        logInfoFn = function(...) {
            stop("logInfoFn should not be called on failure")
        },
        logWarnFn = function(...) {
            stop("logWarnFn should not be called on failure")
        },
        logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        },
        lipidIntensityFilteringFn = function(...) {
            stop("lipidIntensityFilteringFn should not be called on failure")
        },
        updateLipidFilteringFn = function(...) {
            stop("updateLipidFilteringFn should not be called on failure")
        }
    )

    expect_identical(result, input)
    expect_false(filter_stats_called)
    expect_false(filter_plot_called)
    expect_identical(
        error_messages,
        "Error applying lipid intensity filter: Current state is not a LipidomicsAssayData object"
    )
    expect_identical(removed_notifications, "lipid_intensity_filter_working")
    expect_identical(notifications[[1]], list(
        message = "Applying lipid intensity filter...",
        type = NULL,
        id = "lipid_intensity_filter_working",
        duration = NULL
    ))
    expect_identical(notifications[[2]], list(
        message = "Error applying lipid intensity filter: Current state is not a LipidomicsAssayData object",
        type = "error",
        id = NULL,
        duration = 15
    ))
})

test_that("registerLipidIntensityRevertObserver keeps the revert observer shell handoff stable", {
    input <- list(revert_filter = 9L)
    reverted_state <- NULL
    filter_stats_value <- "stale-stats"
    filter_plot_value <- "stale-plot"
    info_messages <- character()
    notifications <- list()
    output_env <- new.env(parent = emptyenv())

    result <- registerLipidIntensityRevertObserver(
        input = input,
        workflowData = list(
            state_manager = list(
                getHistory = function() c("initial", "lipid_intensity_filtered", "lipid_norm"),
                revertToState = function(state_name) {
                    reverted_state <<- state_name
                    invisible(state_name)
                }
            )
        ),
        output = output_env,
        filterStatsVal = function(value) {
            filter_stats_value <<- value
        },
        filterPlotVal = function(value) {
            filter_plot_value <<- value
        },
        observeEventFn = function(event, handler) {
            expect_identical(event, 9L)
            force(handler)
            invisible(NULL)
        },
        renderTextFn = function(expr) force(expr),
        logInfoFn = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logErrorFn = function(...) {
            stop("logErrorFn should not be called on success")
        },
        showNotificationFn = function(message, type = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type)))
        }
    )

    expect_identical(result, input)
    expect_identical(reverted_state, "lipid_intensity_filtered")
    expect_identical(output_env$filter_results, "Reverted to previous state: lipid_intensity_filtered")
    expect_null(filter_stats_value)
    expect_null(filter_plot_value)
    expect_identical(info_messages, "Reverted lipid intensity filter to lipid_intensity_filtered")
    expect_identical(notifications, list(list(
        message = "Reverted successfully",
        type = "message"
    )))
})

test_that("registerLipidIntensityRevertObserver keeps the revert failure notification stable", {
    input <- list(revert_filter = 10L)
    notifications <- list()
    error_messages <- character()

    result <- registerLipidIntensityRevertObserver(
        input = input,
        workflowData = list(
            state_manager = list(
                getHistory = function() "initial",
                revertToState = function(...) {
                    stop("revertToState should not be called without prior history")
                }
            )
        ),
        output = new.env(parent = emptyenv()),
        filterStatsVal = function(...) {
            stop("filterStatsVal should not be called on failure")
        },
        filterPlotVal = function(...) {
            stop("filterPlotVal should not be called on failure")
        },
        observeEventFn = function(event, handler) {
            expect_identical(event, 10L)
            force(handler)
            invisible(NULL)
        },
        renderTextFn = function(...) {
            stop("renderTextFn should not be called on failure")
        },
        logInfoFn = function(...) {
            stop("logInfoFn should not be called on failure")
        },
        logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        },
        showNotificationFn = function(message, type = NULL) {
            notifications <<- c(notifications, list(list(message = message, type = type)))
        }
    )

    expect_identical(result, input)
    expect_identical(error_messages, "Error reverting: No previous state to revert to.")
    expect_identical(notifications, list(list(
        message = "Error reverting: No previous state to revert to.",
        type = "error"
    )))
})

loadLipidQcIntensityModuleHarness <- function() {
    source_lines <- readLines(
        test_path("..", "..", "R", "mod_lipid_qc_intensity.R"),
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

    module_env$assayResultsOutputCalls <- list()
    module_env$registerLipidIntensityAssayResultsOutput <- function(output, filterStatsVal, ns, renderUiFn = shiny::renderUI) {
        module_env$assayResultsOutputCalls <- c(
            module_env$assayResultsOutputCalls,
            list(list(
                output = output,
                filterStatsVal = filterStatsVal,
                assayStatsTabId = ns("assay_stats_tabs"),
                renderUiFn = renderUiFn
            ))
        )
        output
    }
    module_env$filterPlotOutputCalls <- list()
    module_env$registerLipidIntensityFilterPlotOutput <- function(output, filterPlotVal, renderPlotFn = shiny::renderPlot, reqFn = shiny::req, gridDrawFn = grid::grid.draw, printFn = print) {
        module_env$filterPlotOutputCalls <- c(
            module_env$filterPlotOutputCalls,
            list(list(
                output = output,
                filterPlotVal = filterPlotVal,
                renderPlotFn = renderPlotFn,
                reqFn = reqFn,
                gridDrawFn = gridDrawFn,
                printFn = printFn
            ))
        )
        output
    }
    module_env$revertObserverCalls <- list()
    module_env$registerLipidIntensityRevertObserver <- function(input, workflowData, output, filterStatsVal, filterPlotVal, observeEventFn = shiny::observeEvent, renderTextFn = shiny::renderText, logInfoFn = logger::log_info, logErrorFn = logger::log_error, showNotificationFn = shiny::showNotification) {
        module_env$revertObserverCalls <- c(
            module_env$revertObserverCalls,
            list(list(
                input = input,
                workflowData = workflowData,
                output = output,
                filterStatsVal = filterStatsVal,
                filterPlotVal = filterPlotVal,
                observeEventFn = observeEventFn,
                renderTextFn = renderTextFn,
                logInfoFn = logInfoFn,
                logErrorFn = logErrorFn,
                showNotificationFn = showNotificationFn
            ))
        )
        invisible(input)
    }
    module_env$applyFilterObserverCalls <- list()
    module_env$registerLipidIntensityApplyFilterObserver <- function(input, workflowData, output, filterStatsVal, filterPlotVal, omicType, observeEventFn = shiny::observeEvent, reqFn = shiny::req, showNotificationFn = shiny::showNotification, removeNotificationFn = shiny::removeNotification, renderTextFn = shiny::renderText, logInfoFn = logger::log_info, logWarnFn = logger::log_warn, logErrorFn = logger::log_error, lipidIntensityFilteringFn = NULL, updateLipidFilteringFn = NULL) {
        module_env$applyFilterObserverCalls <- c(
            module_env$applyFilterObserverCalls,
            list(list(
                input = input,
                workflowData = workflowData,
                output = output,
                filterStatsVal = filterStatsVal,
                filterPlotVal = filterPlotVal,
                omicType = omicType,
                observeEventFn = observeEventFn,
                reqFn = reqFn,
                showNotificationFn = showNotificationFn,
                removeNotificationFn = removeNotificationFn,
                renderTextFn = renderTextFn,
                logInfoFn = logInfoFn,
                logWarnFn = logWarnFn,
                logErrorFn = logErrorFn,
                lipidIntensityFilteringFn = lipidIntensityFilteringFn,
                updateLipidFilteringFn = updateLipidFilteringFn
            ))
        )
        invisible(input)
    }

    workflow_data <- list(
        state_manager = list(
            getState = function(...) NULL,
            getHistory = function(...) character(),
            revertToState = function(...) invisible(NULL),
            saveState = function(...) invisible(NULL)
        ),
        config_list = list()
    )

    module_env$mod_lipid_qc_intensity_server(
        "intensity",
        workflow_data = workflow_data,
        omic_type = "lipidomics",
        experiment_label = "Lipidomics"
    )

    session <- shiny:::MockShinySession$new()
    shiny::callModule(module_env$capturedModule, "intensity", session = session)
    session$flushReact()

    list(moduleEnv = module_env, session = session)
}

test_that("mod_lipid_qc_intensity_server delegates apply-filter observer registration through the seam", {
    harness <- loadLipidQcIntensityModuleHarness()

    expect_length(harness$moduleEnv$applyFilterObserverCalls, 1L)
    expect_s3_class(harness$moduleEnv$applyFilterObserverCalls[[1]]$input, "reactivevalues")
    expect_true(is.function(harness$moduleEnv$applyFilterObserverCalls[[1]]$filterStatsVal))
    expect_true(is.function(harness$moduleEnv$applyFilterObserverCalls[[1]]$filterPlotVal))
    expect_identical(harness$moduleEnv$applyFilterObserverCalls[[1]]$omicType, "lipidomics")
    expect_true(is.function(harness$moduleEnv$applyFilterObserverCalls[[1]]$observeEventFn))
    expect_identical(harness$moduleEnv$applyFilterObserverCalls[[1]]$reqFn, shiny::req)
    expect_true(is.function(harness$moduleEnv$applyFilterObserverCalls[[1]]$renderTextFn))
})

test_that("mod_lipid_qc_intensity_server delegates assay-results output registration through the seam", {
    harness <- loadLipidQcIntensityModuleHarness()

    expect_length(harness$moduleEnv$assayResultsOutputCalls, 1L)
    expect_true(is.function(harness$moduleEnv$assayResultsOutputCalls[[1]]$filterStatsVal))
    expect_identical(
        harness$moduleEnv$assayResultsOutputCalls[[1]]$assayStatsTabId,
        "intensity-assay_stats_tabs"
    )
    expect_true(is.function(harness$moduleEnv$assayResultsOutputCalls[[1]]$renderUiFn))
})

test_that("mod_lipid_qc_intensity_server delegates filter-plot output registration through the seam", {
    harness <- loadLipidQcIntensityModuleHarness()

    expect_length(harness$moduleEnv$filterPlotOutputCalls, 1L)
    expect_true(is.function(harness$moduleEnv$filterPlotOutputCalls[[1]]$filterPlotVal))
    expect_true(is.function(harness$moduleEnv$filterPlotOutputCalls[[1]]$renderPlotFn))
    expect_identical(harness$moduleEnv$filterPlotOutputCalls[[1]]$reqFn, shiny::req)
})

test_that("mod_lipid_qc_intensity_server delegates revert observer registration through the seam", {
    harness <- loadLipidQcIntensityModuleHarness()

    expect_length(harness$moduleEnv$revertObserverCalls, 1L)
    expect_s3_class(harness$moduleEnv$revertObserverCalls[[1]]$input, "reactivevalues")
    expect_true(is.function(harness$moduleEnv$revertObserverCalls[[1]]$filterStatsVal))
    expect_true(is.function(harness$moduleEnv$revertObserverCalls[[1]]$filterPlotVal))
    expect_true(is.function(harness$moduleEnv$revertObserverCalls[[1]]$observeEventFn))
    expect_true(is.function(harness$moduleEnv$revertObserverCalls[[1]]$renderTextFn))
})
