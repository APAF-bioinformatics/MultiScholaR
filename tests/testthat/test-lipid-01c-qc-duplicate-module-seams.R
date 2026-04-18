library(testthat)

source(test_path("..", "..", "R", "mod_lipid_qc_duplicates_server_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_qc_duplicates.R"), local = environment())

test_that("initializeLipidDuplicateServerState keeps reactive-value defaults stable", {
    reactive_val_defaults <- logical()
    reactive_val_calls <- 0L

    state <- initializeLipidDuplicateServerState(
        reactiveValFn = function(value) {
            reactive_val_calls <<- reactive_val_calls + 1L
            reactive_val_defaults <<- c(reactive_val_defaults, is.null(value))

            structure(
                list(index = reactive_val_calls, value = value),
                class = "mock_lipid_duplicate_reactive"
            )
        }
    )

    expect_identical(names(state), c("duplicateInfo", "resolutionStats", "filterPlot"))
    expect_identical(reactive_val_calls, 3L)
    expect_identical(reactive_val_defaults, c(TRUE, TRUE, TRUE))
    expect_s3_class(state$duplicateInfo, "mock_lipid_duplicate_reactive")
    expect_s3_class(state$resolutionStats, "mock_lipid_duplicate_reactive")
    expect_s3_class(state$filterPlot, "mock_lipid_duplicate_reactive")
    expect_identical(state$duplicateInfo$index, 1L)
    expect_identical(state$resolutionStats$index, 2L)
    expect_identical(state$filterPlot$index, 3L)
})

test_that("buildLipidDuplicateSummaryUi keeps the pre-detection prompt stable", {
    rendered <- htmltools::renderTags(
        buildLipidDuplicateSummaryUi(NULL)
    )$html

    expect_match(rendered, "Detect Duplicates", fixed = TRUE)
    expect_match(rendered, "circle-info", fixed = TRUE)
})

test_that("buildLipidDuplicateSummaryUi keeps per-assay counts and status icons stable", {
    rendered <- htmltools::renderTags(
        buildLipidDuplicateSummaryUi(list(
            AssayA = NULL,
            "Assay 1/A" = data.frame(
                lipid_id = c("L1", "L2"),
                count = c(2L, 3L),
                stringsAsFactors = FALSE
            )
        ))
    )$html

    expect_match(rendered, "AssayA: 0 duplicate IDs", fixed = TRUE)
    expect_match(rendered, "Assay 1/A: 2 duplicate IDs", fixed = TRUE)
    expect_match(rendered, "circle-check", fixed = TRUE)
    expect_match(rendered, "triangle-exclamation", fixed = TRUE)
})

test_that("registerLipidDuplicateSummaryOutput keeps the summary render registration shell stable", {
    output_env <- new.env(parent = emptyenv())
    duplicate_info_calls <- 0L
    build_arg <- NULL

    result <- registerLipidDuplicateSummaryOutput(
        output = output_env,
        duplicateInfoVal = function() {
            duplicate_info_calls <<- duplicate_info_calls + 1L

            list(
                AssayA = NULL,
                `Assay 1` = data.frame(
                    lipid_id = c("L1", "L2"),
                    count = c(2L, 3L),
                    stringsAsFactors = FALSE
                )
            )
        },
        renderUiFn = function(expr) {
            htmltools::renderTags(expr)$html
        },
        buildSummaryUiFn = function(dupList) {
            build_arg <<- dupList
            shiny::tags$span("summary-shell-ok")
        }
    )

    expect_identical(result, output_env)
    expect_identical(duplicate_info_calls, 1L)
    expect_identical(names(build_arg), c("AssayA", "Assay 1"))
    expect_equal(nrow(build_arg[["Assay 1"]]), 2)
    expect_match(as.character(output_env$duplicate_summary), "summary-shell-ok", fixed = TRUE)
})

test_that("registerLipidDuplicateTablesOutput keeps the duplicate-table render registration shell stable", {
    output_env <- new.env(parent = emptyenv())
    duplicate_info_calls <- 0L
    build_args <- NULL

    result <- registerLipidDuplicateTablesOutput(
        output = output_env,
        duplicateInfoVal = function() {
            duplicate_info_calls <<- duplicate_info_calls + 1L

            list(
                AssayA = NULL,
                `Assay 1` = data.frame(
                    lipid_id = c("L1", "L2"),
                    count = c(2L, 3L),
                    stringsAsFactors = FALSE
                )
            )
        },
        ns = function(id) paste0("dup-", id),
        renderUiFn = function(expr) {
            htmltools::renderTags(expr)$html
        },
        buildTablesUiFn = function(dupList, ns) {
            build_args <<- list(
                dupList = dupList,
                ns_value = ns("tables-shell")
            )

            shiny::tags$span("tables-shell-ok")
        }
    )

    expect_identical(result, output_env)
    expect_identical(duplicate_info_calls, 1L)
    expect_identical(names(build_args$dupList), c("AssayA", "Assay 1"))
    expect_equal(nrow(build_args$dupList[["Assay 1"]]), 2)
    expect_identical(build_args$ns_value, "dup-tables-shell")
    expect_match(as.character(output_env$duplicate_tables), "tables-shell-ok", fixed = TRUE)
})

test_that("registerLipidDuplicateFilterPlotOutput keeps the QC-progress grob render shell stable", {
    output_env <- new.env(parent = emptyenv())
    req_values <- list()
    drawn_plot <- NULL
    plot_object <- structure(list(token = "plot-grob"), class = "grob")

    result <- registerLipidDuplicateFilterPlotOutput(
        output = output_env,
        filterPlotVal = function() plot_object,
        renderPlotFn = function(expr) {
            force(expr)
            "plot-shell-ok"
        },
        reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            invisible(value)
        },
        gridDrawFn = function(value) {
            drawn_plot <<- value
            invisible(NULL)
        },
        printFn = function(...) {
            stop("printFn should not be called for grob plots")
        }
    )

    expect_identical(result, output_env)
    expect_identical(req_values, list(plot_object))
    expect_identical(drawn_plot, plot_object)
    expect_identical(output_env$filter_plot, "plot-shell-ok")
})

test_that("registerLipidDuplicateFilterPlotOutput keeps the QC-progress ggplot render shell stable", {
    output_env <- new.env(parent = emptyenv())
    req_values <- list()
    printed_plot <- NULL
    plot_object <- structure(list(token = "plot-ggplot"), class = "ggplot")

    result <- registerLipidDuplicateFilterPlotOutput(
        output = output_env,
        filterPlotVal = function() plot_object,
        renderPlotFn = function(expr) {
            force(expr)
            "plot-shell-ok"
        },
        reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            invisible(value)
        },
        gridDrawFn = function(...) {
            stop("gridDrawFn should not be called for ggplot plots")
        },
        printFn = function(value) {
            printed_plot <<- value
            invisible(NULL)
        }
    )

    expect_identical(result, output_env)
    expect_identical(req_values, list(plot_object))
    expect_identical(printed_plot, plot_object)
    expect_identical(output_env$filter_plot, "plot-shell-ok")
})

test_that("buildLipidDuplicateTablesUi keeps the pre-detection empty state stable", {
    expect_null(
        buildLipidDuplicateTablesUi(
            dupList = NULL,
            ns = function(id) paste0("dup-", id)
        )
    )
})

test_that("buildLipidDuplicateTablesUi keeps the no-duplicates panel stable", {
    rendered <- htmltools::renderTags(
        buildLipidDuplicateTablesUi(
            dupList = list(
                AssayA = NULL,
                AssayB = data.frame(lipid_id = character(), count = integer())
            ),
            ns = function(id) paste0("dup-", id)
        )
    )$html

    expect_match(rendered, "No duplicates found in any assay!", fixed = TRUE)
    expect_match(rendered, "All lipid IDs are unique. No resolution needed.", fixed = TRUE)
})

test_that("buildLipidDuplicateTablesUi keeps duplicate tab ids and assay labels stable", {
    rendered <- htmltools::renderTags(
        buildLipidDuplicateTablesUi(
            dupList = list(
                "Assay 1/A" = data.frame(
                    lipid_id = c("L1", "L2"),
                    count = c(2L, 3L),
                    stringsAsFactors = FALSE
                ),
                AssayB = NULL
            ),
            ns = function(id) paste0("dup-", id)
        )
    )$html

    expect_match(rendered, "Assay 1/A", fixed = TRUE)
    expect_match(rendered, "dup-dup_tables_tabs", fixed = TRUE)
    expect_match(rendered, "dup-dup_table_Assay_1_A", fixed = TRUE)
    expect_false(grepl("AssayB", rendered, fixed = TRUE))
})

test_that("registerLipidDuplicateTableOutputs keeps sanitized ids and duplicate-only registration stable", {
    output_env <- new.env(parent = emptyenv())

    registered_ids <- registerLipidDuplicateTableOutputs(
        output = output_env,
        dupList = list(
            AssayA = NULL,
            "Assay 1/A" = data.frame(
                lipid_id = c("L1", "L2"),
                count = c(2L, 3L),
                stringsAsFactors = FALSE
            ),
            AssayB = data.frame(lipid_id = character(), count = integer())
        ),
        renderTable = function(dupDf) {
            list(
                row_count = nrow(dupDf),
                lipid_ids = dupDf$lipid_id
            )
        }
    )

    expect_equal(registered_ids, "dup_table_Assay_1_A")
    expect_equal(ls(output_env), "dup_table_Assay_1_A")
    expect_equal(output_env[["dup_table_Assay_1_A"]]$row_count, 2)
    expect_equal(output_env[["dup_table_Assay_1_A"]]$lipid_ids, c("L1", "L2"))
})

test_that("registerLipidDuplicateTableObserver keeps the table-registration observer shell stable", {
    output_env <- new.env(parent = emptyenv())
    register_args <- NULL
    req_values <- list()

    result <- registerLipidDuplicateTableObserver(
        output = output_env,
        duplicateInfoVal = function() {
            list(
                AssayA = NULL,
                `Assay 1` = data.frame(
                    lipid_id = c("L1", "L2"),
                    count = c(2L, 3L),
                    stringsAsFactors = FALSE
                )
            )
        },
        observeFn = function(handler) {
            force(handler)
            invisible(NULL)
        },
        reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            invisible(value)
        },
        registerTableOutputsFn = function(output, dupList) {
            register_args <<- list(
                output = output,
                dupList = dupList
            )

            "dup_table_Assay_1"
        }
    )

    expect_identical(result, output_env)
    expect_length(req_values, 1)
    expect_identical(req_values[[1]], register_args$dupList)
    expect_identical(register_args$output, output_env)
    expect_identical(names(register_args$dupList), c("AssayA", "Assay 1"))
    expect_equal(nrow(register_args$dupList[["Assay 1"]]), 2)
})

ensureLipidDuplicateWorkflowTestClass <- function() {
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

test_that("handleLipidDuplicateDetection keeps the detect-summary workflow stable", {
    ensureLipidDuplicateWorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        lipid_id_column = "lipid_id"
    )
    info_messages <- character()

    result <- handleLipidDuplicateDetection(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        reqFn = function(...) invisible(NULL),
        findDuplicatesFn = function(theObject) {
            expect_identical(theObject, current_s4)
            list(
                AssayA = NULL,
                `Assay 1` = data.frame(
                    lipid_id = c("L1", "L2"),
                    count = c(2L, 3L),
                    stringsAsFactors = FALSE
                )
            )
        },
        logInfoFn = function(msg) {
            info_messages <<- c(info_messages, msg)
        }
    )

    expect_equal(names(result$duplicatesList), c("AssayA", "Assay 1"))
    expect_equal(result$totalDuplicates, 2)
    expect_equal(result$notificationType, "warning")
    expect_equal(result$notificationMessage, "Detection complete: 2 duplicate IDs found")
    expect_true(any(grepl("Detected 2 duplicate feature IDs across assays", info_messages, fixed = TRUE)))
})

test_that("handleLipidDuplicateDetection keeps the zero-duplicate notification stable", {
    ensureLipidDuplicateWorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(),
        lipid_id_column = "lipid_id"
    )

    result <- handleLipidDuplicateDetection(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4
            )
        ),
        reqFn = function(...) invisible(NULL),
        findDuplicatesFn = function(...) {
            list(
                AssayA = NULL,
                AssayB = data.frame(lipid_id = character(), count = integer())
            )
        },
        logInfoFn = function(...) invisible(NULL)
    )

    expect_equal(result$totalDuplicates, 0)
    expect_equal(result$notificationType, "message")
    expect_equal(result$notificationMessage, "Detection complete: 0 duplicate IDs found")
})

test_that("applyLipidDuplicateDetectionResult keeps the post-detection apply contract stable", {
    duplicate_info_value <- NULL
    shown_notifications <- list()

    detection_result <- list(
        duplicatesList = list(
            AssayA = NULL,
            `Assay 1` = data.frame(
                lipid_id = c("L1", "L2"),
                count = c(2L, 3L),
                stringsAsFactors = FALSE
            )
        )
        , totalDuplicates = 2
        , notificationMessage = "Detection complete: 2 duplicate IDs found"
        , notificationType = "warning"
    )

    result <- applyLipidDuplicateDetectionResult(
        detectionResult = detection_result
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , showNotificationFn = function(message, type) {
            shown_notifications <<- c(shown_notifications, list(list(
                message = message,
                type = type
            )))
        }
    )

    expect_identical(duplicate_info_value, detection_result$duplicatesList)
    expect_identical(result$totalDuplicates, 2)
    expect_identical(result$notificationMessage, "Detection complete: 2 duplicate IDs found")
    expect_identical(result$notificationType, "warning")
    expect_identical(shown_notifications, list(list(
        message = "Detection complete: 2 duplicate IDs found",
        type = "warning"
    )))
})

test_that("registerLipidDuplicateDetectObserver keeps the detect observer shell handoff stable", {
    input <- list(detect_duplicates = 7L)
    req_values <- list()
    handler_args <- NULL
    apply_args <- NULL
    duplicate_info_value <- NULL
    detection_result <- list(
        duplicatesList = list(
            AssayA = NULL,
            `Assay 1` = data.frame(
                lipid_id = c("L1", "L2"),
                count = c(2L, 3L),
                stringsAsFactors = FALSE
            )
        )
        , totalDuplicates = 2
        , notificationMessage = "Detection complete: 2 duplicate IDs found"
        , notificationType = "warning"
    )

    result <- registerLipidDuplicateDetectObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , observeEventFn = function(event, handler) {
            expect_identical(event, 7L)
            force(handler)
            invisible(NULL)
        }
        , reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            invisible(value)
        }
        , handleDetectionFn = function(workflowData) {
            handler_args <<- list(workflowData = workflowData)
            detection_result
        }
        , applyDetectionResultFn = function(detectionResult, duplicateInfoVal, showNotificationFn) {
            apply_args <<- list(
                detectionResult = detectionResult,
                showNotificationFn = showNotificationFn
            )

            duplicateInfoVal(detectionResult$duplicatesList)
        }
        , logErrorFn = function(...) {
            stop("logErrorFn should not be called")
        }
        , showNotificationFn = function(...) {
            stop("showNotificationFn should be delegated to applyDetectionResultFn")
        }
    )

    expect_identical(result, input)
    expect_identical(req_values, list("state-manager"))
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(apply_args$detectionResult, detection_result)
    expect_true(is.function(apply_args$showNotificationFn))
    expect_identical(duplicate_info_value, detection_result$duplicatesList)
})

test_that("registerLipidDuplicateDetectObserver keeps the detect failure notification stable", {
    input <- list(detect_duplicates = 3L)
    notification_calls <- list()
    error_messages <- character()

    result <- registerLipidDuplicateDetectObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , duplicateInfoVal = function(...) {
            stop("duplicateInfoVal should not be called on failure path")
        }
        , observeEventFn = function(event, handler) {
            expect_identical(event, 3L)
            force(handler)
            invisible(NULL)
        }
        , reqFn = function(value) invisible(value)
        , handleDetectionFn = function(...) {
            stop("detect exploded")
        }
        , applyDetectionResultFn = function(...) {
            stop("applyDetectionResultFn should not be called")
        }
        , logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        }
        , showNotificationFn = function(message, type = NULL) {
            notification_calls <<- c(notification_calls, list(list(
                message = message,
                type = type
            )))
        }
    )

    expect_identical(result, input)
    expect_identical(notification_calls, list(list(
        message = "Error detecting duplicates: detect exploded",
        type = "error"
    )))
    expect_identical(error_messages, "Error detecting duplicates: detect exploded")
})

test_that("handleLipidDuplicateResolution keeps the resolve-save-summary workflow stable", {
    ensureLipidDuplicateWorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            `Assay 1` = data.frame(
                lipid_id = c("L1", "L1", "L2"),
                sample_a = c(10, 4, 7),
                sample_b = c(9, 3, 6),
                label = c("keep", "drop", "single"),
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id"
    )
    saved_state <- NULL
    filtering_call <- NULL
    warn_messages <- character()

    result <- handleLipidDuplicateResolution(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4,
                saveState = function(state_name, s4_data_object, config_object, description) {
                    saved_state <<- list(
                        state_name = state_name,
                        s4_data_object = s4_data_object,
                        config_object = config_object,
                        description = description
                    )
                }
            ),
            config_list = list(config_marker = TRUE)
        ),
        omicType = "lipidomics",
        reqFn = function(...) invisible(NULL),
        resolveDuplicatesFn = function(assay_tibble, id_col, sample_cols) {
            expect_equal(id_col, "lipid_id")
            expect_equal(sample_cols, c("sample_a", "sample_b"))
            assay_tibble[c(1, 3), , drop = FALSE]
        },
        updateFilteringFn = function(theObject, step_name, omics_type, return_grid, overwrite) {
            filtering_call <<- list(
                rows = nrow(theObject@lipid_data[[1]]),
                step_name = step_name,
                omics_type = omics_type,
                return_grid = return_grid,
                overwrite = overwrite
            )
            "plot-token"
        },
        logWarnFn = function(msg) {
            warn_messages <<- c(warn_messages, msg)
        }
    )

    expect_equal(result$totalRemoved, 1)
    expect_equal(result$statsList[["Assay 1"]]$original, 3)
    expect_equal(result$statsList[["Assay 1"]]$resolved, 2)
    expect_equal(result$statsList[["Assay 1"]]$removed, 1)
    expect_equal(result$qcPlot, "plot-token")
    expect_match(result$resultText, "Assay 1: 3 -> 2 rows \\(removed 1 duplicates\\)")
    expect_match(result$resultText, "Total duplicate rows removed: 1", fixed = TRUE)
    expect_equal(saved_state$state_name, "lipid_duplicates_resolved")
    expect_equal(saved_state$config_object, list(config_marker = TRUE))
    expect_equal(saved_state$description, "Resolved duplicate lipid features by keeping highest intensity")
    expect_equal(nrow(saved_state$s4_data_object@lipid_data[[1]]), 2)
    expect_equal(filtering_call$rows, 2)
    expect_equal(filtering_call$step_name, "3_Duplicates_Resolved")
    expect_equal(filtering_call$omics_type, "lipidomics")
    expect_true(filtering_call$return_grid)
    expect_true(filtering_call$overwrite)
    expect_length(warn_messages, 0)
})

test_that("handleLipidDuplicateResolution keeps the no-numeric fallback and QC-plot warning stable", {
    ensureLipidDuplicateWorkflowTestClass()

    current_s4 <- methods::new(
        "LipidomicsAssayData",
        lipid_data = list(
            data.frame(
                lipid_id = c("L1", "L1"),
                label = c("x", "y"),
                stringsAsFactors = FALSE
            )
        ),
        lipid_id_column = "lipid_id"
    )
    saved_state <- NULL
    warn_messages <- character()

    result <- handleLipidDuplicateResolution(
        workflowData = list(
            state_manager = list(
                getState = function() current_s4,
                saveState = function(state_name, s4_data_object, config_object, description) {
                    saved_state <<- list(
                        state_name = state_name,
                        s4_data_object = s4_data_object,
                        config_object = config_object,
                        description = description
                    )
                }
            ),
            config_list = list()
        ),
        omicType = "lipidomics",
        reqFn = function(...) invisible(NULL),
        resolveDuplicatesFn = function(...) {
            stop("resolveDuplicatesFn should not be called")
        },
        updateFilteringFn = function(...) {
            stop("plot exploded")
        },
        logWarnFn = function(msg) {
            warn_messages <<- c(warn_messages, msg)
        }
    )

    expect_equal(names(result$statsList), "Assay_1")
    expect_equal(result$statsList[["Assay_1"]]$original, 2)
    expect_equal(result$statsList[["Assay_1"]]$resolved, 2)
    expect_equal(result$statsList[["Assay_1"]]$removed, 0)
    expect_null(result$qcPlot)
    expect_equal(result$totalRemoved, 0)
    expect_match(result$resultText, "Total duplicate rows removed: 0", fixed = TRUE)
    expect_equal(nrow(saved_state$s4_data_object@lipid_data[[1]]), 2)
    expect_true(any(grepl("No numeric columns found in assay: Assay_1", warn_messages, fixed = TRUE)))
    expect_true(any(grepl("Could not generate QC plot: plot exploded", warn_messages, fixed = TRUE)))
})

test_that("applyLipidDuplicateResolutionResult keeps the post-resolution reset and notification contract stable", {
    resolution_stats_value <- "stale-stats"
    duplicate_info_value <- list(stale = TRUE)
    filter_plot_value <- "stale-plot"
    rendered_text <- NULL
    info_messages <- character()
    removed_notification_ids <- character()
    shown_notifications <- list()
    output_env <- new.env(parent = emptyenv())

    resolution_result <- list(
        statsList = list(AssayA = list(original = 3, resolved = 2, removed = 1))
        , qcPlot = "plot-token"
        , totalRemoved = 1
        , resultText = "Duplicate Resolution Complete"
    )

    result <- applyLipidDuplicateResolutionResult(
        resolutionResult = resolution_result
        , resolutionStatsVal = function(value) {
            resolution_stats_value <<- value
        }
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , filterPlotVal = function(value) {
            filter_plot_value <<- value
        }
        , output = output_env
        , renderTextFn = function(text) {
            rendered_text <<- text
            list(rendered = text)
        }
        , logInfoFn = function(msg) {
            info_messages <<- c(info_messages, msg)
        }
        , removeNotificationFn = function(id) {
            removed_notification_ids <<- c(removed_notification_ids, id)
        }
        , showNotificationFn = function(message, type) {
            shown_notifications <<- c(shown_notifications, list(list(
                message = message,
                type = type
            )))
        }
    )

    expect_identical(resolution_stats_value, resolution_result$statsList)
    expect_null(duplicate_info_value)
    expect_identical(filter_plot_value, "plot-token")
    expect_identical(rendered_text, "Duplicate Resolution Complete")
    expect_identical(output_env$resolution_results, list(rendered = "Duplicate Resolution Complete"))
    expect_identical(result$successMessage, "Duplicates resolved: 1 rows removed")
    expect_identical(result$notificationType, "message")
    expect_identical(result$workingNotificationId, "lipid_dup_resolve_working")
    expect_identical(removed_notification_ids, "lipid_dup_resolve_working")
    expect_identical(shown_notifications, list(list(
        message = "Duplicates resolved: 1 rows removed",
        type = "message"
    )))
    expect_true(any(grepl("Resolved duplicates: removed 1 rows", info_messages, fixed = TRUE)))
})

test_that("registerLipidDuplicateResolveObserver keeps the resolve observer shell handoff stable", {
    input <- list(resolve_duplicates = 9L)
    req_values <- list()
    notification_calls <- list()
    handler_args <- NULL
    apply_args <- NULL
    output_env <- new.env(parent = emptyenv())
    resolution_stats_value <- NULL
    duplicate_info_value <- list(stale = TRUE)
    filter_plot_value <- NULL
    resolution_result <- list(
        statsList = list(AssayA = list(original = 3, resolved = 2, removed = 1))
        , qcPlot = "plot-token"
        , totalRemoved = 1
        , resultText = "Duplicate Resolution Complete"
    )

    result <- registerLipidDuplicateResolveObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , omicType = "lipidomics"
        , resolutionStatsVal = function(value) {
            resolution_stats_value <<- value
        }
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , filterPlotVal = function(value) {
            filter_plot_value <<- value
        }
        , output = output_env
        , observeEventFn = function(event, handler) {
            expect_identical(event, 9L)
            force(handler)
            invisible(NULL)
        }
        , reqFn = function(value) {
            req_values <<- c(req_values, list(value))
            invisible(value)
        }
        , showNotificationFn = function(message, type = NULL, id = NULL, duration = NULL) {
            notification_calls <<- c(notification_calls, list(list(
                message = message,
                type = type,
                id = id,
                duration = duration
            )))
        }
        , handleResolutionFn = function(workflowData, omicType) {
            handler_args <<- list(
                workflowData = workflowData,
                omicType = omicType
            )
            resolution_result
        }
        , applyResolutionResultFn = function(
            resolutionResult,
            resolutionStatsVal,
            duplicateInfoVal,
            filterPlotVal,
            output,
            removeNotificationFn,
            showNotificationFn,
            workingNotificationId
        ) {
            apply_args <<- list(
                resolutionResult = resolutionResult,
                output = output,
                workingNotificationId = workingNotificationId,
                removeNotificationFn = removeNotificationFn,
                showNotificationFn = showNotificationFn
            )

            resolutionStatsVal(resolutionResult$statsList)
            duplicateInfoVal(NULL)
            filterPlotVal(resolutionResult$qcPlot)
            output$resolution_results <- resolutionResult$resultText
        }
        , logErrorFn = function(...) {
            stop("logErrorFn should not be called")
        }
        , removeNotificationFn = function(...) {
            stop("removeNotificationFn should not be called on success path")
        }
    )

    expect_identical(result, input)
    expect_identical(req_values, list("state-manager"))
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(handler_args$omicType, "lipidomics")
    expect_identical(notification_calls, list(list(
        message = "Resolving duplicate features...",
        type = NULL,
        id = "lipid_dup_resolve_working",
        duration = NULL
    )))
    expect_identical(apply_args$resolutionResult, resolution_result)
    expect_identical(apply_args$output, output_env)
    expect_identical(apply_args$workingNotificationId, "lipid_dup_resolve_working")
    expect_true(is.function(apply_args$removeNotificationFn))
    expect_true(is.function(apply_args$showNotificationFn))
    expect_identical(resolution_stats_value, resolution_result$statsList)
    expect_null(duplicate_info_value)
    expect_identical(filter_plot_value, "plot-token")
    expect_identical(output_env$resolution_results, "Duplicate Resolution Complete")
})

test_that("registerLipidDuplicateResolveObserver keeps failure cleanup stable", {
    input <- list(resolve_duplicates = 4L)
    notification_calls <- list()
    error_messages <- character()
    removed_notification_ids <- character()

    result <- registerLipidDuplicateResolveObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , omicType = "lipidomics"
        , resolutionStatsVal = function(...) {
            stop("resolutionStatsVal should not be called on failure path")
        }
        , duplicateInfoVal = function(...) {
            stop("duplicateInfoVal should not be called on failure path")
        }
        , filterPlotVal = function(...) {
            stop("filterPlotVal should not be called on failure path")
        }
        , output = new.env(parent = emptyenv())
        , observeEventFn = function(event, handler) {
            expect_identical(event, 4L)
            force(handler)
            invisible(NULL)
        }
        , reqFn = function(value) invisible(value)
        , showNotificationFn = function(message, type = NULL, id = NULL, duration = NULL) {
            notification_calls <<- c(notification_calls, list(list(
                message = message,
                type = type,
                id = id,
                duration = duration
            )))
        }
        , handleResolutionFn = function(...) {
            stop("resolution exploded")
        }
        , applyResolutionResultFn = function(...) {
            stop("applyResolutionResultFn should not be called")
        }
        , logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        }
        , removeNotificationFn = function(id) {
            removed_notification_ids <<- c(removed_notification_ids, id)
        }
    )

    expect_identical(result, input)
    expect_identical(notification_calls, list(
        list(
            message = "Resolving duplicate features...",
            type = NULL,
            id = "lipid_dup_resolve_working",
            duration = NULL
        ),
        list(
            message = "Error resolving duplicates: resolution exploded",
            type = "error",
            id = NULL,
            duration = 15
        )
    ))
    expect_identical(error_messages, "Error resolving duplicates: resolution exploded")
    expect_identical(removed_notification_ids, "lipid_dup_resolve_working")
})

test_that("applyLipidDuplicateRevertResult keeps the post-revert reset and notification contract stable", {
    resolution_stats_value <- "stale-stats"
    duplicate_info_value <- list(stale = TRUE)
    filter_plot_value <- "stale-plot"
    rendered_text <- NULL
    shown_notifications <- list()
    output_env <- new.env(parent = emptyenv())

    revert_result <- list(
        resultText = "Reverted to previous state: lipid_duplicates_resolved"
        , notificationMessage = "Reverted successfully"
        , notificationType = "message"
    )

    result <- applyLipidDuplicateRevertResult(
        revertResult = revert_result
        , resolutionStatsVal = function(value) {
            resolution_stats_value <<- value
        }
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , filterPlotVal = function(value) {
            filter_plot_value <<- value
        }
        , output = output_env
        , renderTextFn = function(text) {
            rendered_text <<- text
            list(rendered = text)
        }
        , showNotificationFn = function(message, type) {
            shown_notifications <<- c(shown_notifications, list(list(
                message = message,
                type = type
            )))
        }
    )

    expect_null(resolution_stats_value)
    expect_null(duplicate_info_value)
    expect_null(filter_plot_value)
    expect_identical(rendered_text, revert_result$resultText)
    expect_identical(output_env$resolution_results, list(rendered = revert_result$resultText))
    expect_identical(result$notificationMessage, "Reverted successfully")
    expect_identical(result$notificationType, "message")
    expect_identical(shown_notifications, list(list(
        message = "Reverted successfully",
        type = "message"
    )))
})

test_that("registerLipidDuplicateRevertObserver keeps the revert observer shell handoff stable", {
    input <- list(revert_duplicates = 5L)
    handler_args <- NULL
    apply_args <- NULL
    output_env <- new.env(parent = emptyenv())
    resolution_stats_value <- "stale-stats"
    duplicate_info_value <- list(stale = TRUE)
    filter_plot_value <- "stale-plot"
    revert_result <- list(
        resultText = "Reverted to previous state: lipid_duplicates_resolved"
        , notificationMessage = "Reverted successfully"
        , notificationType = "message"
    )

    result <- registerLipidDuplicateRevertObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , resolutionStatsVal = function(value) {
            resolution_stats_value <<- value
        }
        , duplicateInfoVal = function(value) {
            duplicate_info_value <<- value
        }
        , filterPlotVal = function(value) {
            filter_plot_value <<- value
        }
        , output = output_env
        , observeEventFn = function(event, handler) {
            expect_identical(event, 5L)
            force(handler)
            invisible(NULL)
        }
        , handleRevertFn = function(workflowData) {
            handler_args <<- list(workflowData = workflowData)
            revert_result
        }
        , applyRevertResultFn = function(
            revertResult,
            resolutionStatsVal,
            duplicateInfoVal,
            filterPlotVal,
            output,
            showNotificationFn
        ) {
            apply_args <<- list(
                revertResult = revertResult,
                output = output,
                showNotificationFn = showNotificationFn
            )

            resolutionStatsVal(NULL)
            duplicateInfoVal(NULL)
            filterPlotVal(NULL)
            output$resolution_results <- revertResult$resultText
        }
        , logErrorFn = function(...) {
            stop("logErrorFn should not be called")
        }
        , showNotificationFn = function(...) {
            stop("showNotificationFn should be delegated to applyRevertResultFn")
        }
    )

    expect_identical(result, input)
    expect_identical(handler_args$workflowData$state_manager, "state-manager")
    expect_identical(apply_args$revertResult, revert_result)
    expect_identical(apply_args$output, output_env)
    expect_true(is.function(apply_args$showNotificationFn))
    expect_null(resolution_stats_value)
    expect_null(duplicate_info_value)
    expect_null(filter_plot_value)
    expect_identical(output_env$resolution_results, revert_result$resultText)
})

test_that("registerLipidDuplicateRevertObserver keeps the revert failure notification stable", {
    input <- list(revert_duplicates = 6L)
    notification_calls <- list()
    error_messages <- character()

    result <- registerLipidDuplicateRevertObserver(
        input = input
        , workflowData = list(state_manager = "state-manager")
        , resolutionStatsVal = function(...) {
            stop("resolutionStatsVal should not be called on failure path")
        }
        , duplicateInfoVal = function(...) {
            stop("duplicateInfoVal should not be called on failure path")
        }
        , filterPlotVal = function(...) {
            stop("filterPlotVal should not be called on failure path")
        }
        , output = new.env(parent = emptyenv())
        , observeEventFn = function(event, handler) {
            expect_identical(event, 6L)
            force(handler)
            invisible(NULL)
        }
        , handleRevertFn = function(...) {
            stop("revert exploded")
        }
        , applyRevertResultFn = function(...) {
            stop("applyRevertResultFn should not be called")
        }
        , logErrorFn = function(message) {
            error_messages <<- c(error_messages, message)
        }
        , showNotificationFn = function(message, type = NULL) {
            notification_calls <<- c(notification_calls, list(list(
                message = message,
                type = type
            )))
        }
    )

    expect_identical(result, input)
    expect_identical(notification_calls, list(list(
        message = "Error reverting: revert exploded",
        type = "error"
    )))
    expect_identical(error_messages, "Error reverting: revert exploded")
})

test_that("handleLipidDuplicateRevert keeps the revert workflow contract stable", {
    reverted_state <- NULL
    info_messages <- character()

    result <- handleLipidDuplicateRevert(
        workflowData = list(
            state_manager = list(
                getHistory = function() c("initial", "lipid_duplicates_resolved", "lipid_filtered"),
                revertToState = function(state_name) {
                    reverted_state <<- state_name
                    list(reverted = state_name)
                }
            )
        ),
        reqFn = function(...) invisible(NULL),
        logInfoFn = function(msg) {
            info_messages <<- c(info_messages, msg)
        }
    )

    expect_equal(reverted_state, "lipid_duplicates_resolved")
    expect_equal(result$previousStateName, "lipid_duplicates_resolved")
    expect_equal(result$resultText, "Reverted to previous state: lipid_duplicates_resolved")
    expect_equal(result$notificationMessage, "Reverted successfully")
    expect_equal(result$notificationType, "message")
    expect_true(any(grepl(
        "Reverted duplicate resolution to lipid_duplicates_resolved",
        info_messages,
        fixed = TRUE
    )))
})

test_that("handleLipidDuplicateRevert keeps the no-history failure stable", {
    expect_error(
        handleLipidDuplicateRevert(
            workflowData = list(
                state_manager = list(
                    getHistory = function() "initial",
                    revertToState = function(...) {
                        stop("revertToState should not be called")
                    }
                )
            ),
            reqFn = function(...) invisible(NULL),
            logInfoFn = function(...) invisible(NULL)
        ),
        "No previous state to revert to\\."
    )
})

test_that("registerLipidDuplicateServerBindings keeps the live server registration shell stable", {
    input <- list(
        detect_duplicates = 1L,
        resolve_duplicates = 2L,
        revert_duplicates = 3L
    )
    output_env <- new.env(parent = emptyenv())
    workflow_data <- list(state_manager = "state-manager")
    duplicate_info_val <- function(...) NULL
    resolution_stats_val <- function(...) NULL
    filter_plot_val <- function(...) NULL
    call_log <- list()

    result <- registerLipidDuplicateServerBindings(
        input = input
        , output = output_env
        , workflowData = workflow_data
        , omicType = "lipidomics"
        , ns = function(id) paste0("dup-", id)
        , duplicateInfoVal = duplicate_info_val
        , resolutionStatsVal = resolution_stats_val
        , filterPlotVal = filter_plot_val
        , registerDetectObserverFn = function(input, workflowData, duplicateInfoVal) {
            call_log <<- c(call_log, list(list(
                step = "detect",
                input = input,
                workflowData = workflowData,
                duplicateInfoVal = duplicateInfoVal
            )))
        }
        , registerSummaryOutputFn = function(output, duplicateInfoVal) {
            call_log <<- c(call_log, list(list(
                step = "summary",
                output = output,
                duplicateInfoVal = duplicateInfoVal
            )))
        }
        , registerTablesOutputFn = function(output, duplicateInfoVal, ns) {
            call_log <<- c(call_log, list(list(
                step = "tables",
                output = output,
                duplicateInfoVal = duplicateInfoVal,
                ns_value = ns("dup_tables_tabs")
            )))
        }
        , registerTableObserverFn = function(output, duplicateInfoVal) {
            call_log <<- c(call_log, list(list(
                step = "table_observer",
                output = output,
                duplicateInfoVal = duplicateInfoVal
            )))
        }
        , registerResolveObserverFn = function(
            input,
            workflowData,
            omicType,
            resolutionStatsVal,
            duplicateInfoVal,
            filterPlotVal,
            output
        ) {
            call_log <<- c(call_log, list(list(
                step = "resolve",
                input = input,
                workflowData = workflowData,
                omicType = omicType,
                resolutionStatsVal = resolutionStatsVal,
                duplicateInfoVal = duplicateInfoVal,
                filterPlotVal = filterPlotVal,
                output = output
            )))
        }
        , registerRevertObserverFn = function(
            input,
            workflowData,
            resolutionStatsVal,
            duplicateInfoVal,
            filterPlotVal,
            output
        ) {
            call_log <<- c(call_log, list(list(
                step = "revert",
                input = input,
                workflowData = workflowData,
                resolutionStatsVal = resolutionStatsVal,
                duplicateInfoVal = duplicateInfoVal,
                filterPlotVal = filterPlotVal,
                output = output
            )))
        }
        , registerFilterPlotOutputFn = function(output, filterPlotVal) {
            call_log <<- c(call_log, list(list(
                step = "filter_plot",
                output = output,
                filterPlotVal = filterPlotVal
            )))
        }
    )

    expect_identical(result, output_env)
    expect_identical(vapply(call_log, `[[`, character(1), "step"), c(
        "detect",
        "summary",
        "tables",
        "table_observer",
        "resolve",
        "revert",
        "filter_plot"
    ))
    expect_identical(call_log[[1]]$input, input)
    expect_identical(call_log[[1]]$workflowData, workflow_data)
    expect_identical(call_log[[1]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[2]]$output, output_env)
    expect_identical(call_log[[2]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[3]]$output, output_env)
    expect_identical(call_log[[3]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[3]]$ns_value, "dup-dup_tables_tabs")
    expect_identical(call_log[[4]]$output, output_env)
    expect_identical(call_log[[4]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[5]]$input, input)
    expect_identical(call_log[[5]]$workflowData, workflow_data)
    expect_identical(call_log[[5]]$omicType, "lipidomics")
    expect_identical(call_log[[5]]$resolutionStatsVal, resolution_stats_val)
    expect_identical(call_log[[5]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[5]]$filterPlotVal, filter_plot_val)
    expect_identical(call_log[[5]]$output, output_env)
    expect_identical(call_log[[6]]$input, input)
    expect_identical(call_log[[6]]$workflowData, workflow_data)
    expect_identical(call_log[[6]]$resolutionStatsVal, resolution_stats_val)
    expect_identical(call_log[[6]]$duplicateInfoVal, duplicate_info_val)
    expect_identical(call_log[[6]]$filterPlotVal, filter_plot_val)
    expect_identical(call_log[[6]]$output, output_env)
    expect_identical(call_log[[7]]$output, output_env)
    expect_identical(call_log[[7]]$filterPlotVal, filter_plot_val)
})

test_that("runLipidDuplicateModuleServerShell keeps the moduleServer body handoff stable", {
    input <- list(
        detect_duplicates = 1L,
        resolve_duplicates = 2L,
        revert_duplicates = 3L
    )
    output_env <- new.env(parent = emptyenv())
    session <- list(ns = function(id) paste0("dup-", id))
    workflow_data <- list(state_manager = "state-manager")
    duplicate_state <- list(
        duplicateInfo = function(...) NULL,
        resolutionStats = function(...) NULL,
        filterPlot = function(...) NULL
    )
    initialize_state_calls <- 0L
    register_args <- NULL

    result <- runLipidDuplicateModuleServerShell(
        input = input
        , output = output_env
        , session = session
        , workflowData = workflow_data
        , omicType = "lipidomics"
        , initializeServerStateFn = function() {
            initialize_state_calls <<- initialize_state_calls + 1L
            duplicate_state
        }
        , registerServerBindingsFn = function(...) {
            register_args <<- list(...)
            invisible(output_env)
        }
    )

    expect_identical(result, duplicate_state)
    expect_identical(initialize_state_calls, 1L)
    expect_setequal(
        names(register_args),
        c(
            "input", "output", "workflowData", "omicType",
            "ns", "duplicateInfoVal", "resolutionStatsVal", "filterPlotVal"
        )
    )
    expect_identical(register_args$input, input)
    expect_identical(register_args$output, output_env)
    expect_identical(register_args$workflowData, workflow_data)
    expect_identical(register_args$omicType, "lipidomics")
    expect_identical(register_args$ns("dup_tables_tabs"), "dup-dup_tables_tabs")
    expect_identical(register_args$duplicateInfoVal, duplicate_state$duplicateInfo)
    expect_identical(register_args$resolutionStatsVal, duplicate_state$resolutionStats)
    expect_identical(register_args$filterPlotVal, duplicate_state$filterPlot)
})
