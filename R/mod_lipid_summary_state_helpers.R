# ============================================================================
# mod_lipid_summary_state_helpers.R
# ============================================================================
# Purpose: Session state and parameter helpers for the lipid summary module
# ============================================================================

buildLipidSummaryTemplateStatus <- function(projectDirs, omicType, fileExistsFn = file.exists) {
    if (!omicType %in% names(projectDirs)) {
        return("[WARNING] Project directories not available")
    }

    template_dir <- file.path(
        projectDirs[[omicType]]$base_dir,
        "scripts",
        omicType
    )
    lipid_template <- file.path(template_dir, "lipidomics_report.rmd")

    if (fileExistsFn(lipid_template)) {
        "Template: Lipidomics Report [OK]"
    } else {
        "[WARNING] Report template will be downloaded when generating report"
    }
}

registerLipidSummaryTemplateStatusOutput <- function(
    output,
    projectDirs,
    omicType,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    buildStatusFn = buildLipidSummaryTemplateStatus
) {
    output$template_status <- renderTextFn({
        reqFn(projectDirs)
        buildStatusFn(projectDirs, omicType)
    })

    output
}

buildLipidSummarySessionState <- function(
    input,
    values,
    projectDirs,
    omicType,
    timeFn = Sys.time
) {
    list(
        experiment_label = input$experiment_label,
        description = input$description,
        timestamp = timeFn(),
        omic_type = omicType,
        workflow_args_saved = values$workflow_args_saved,
        files_copied = values$files_copied,
        report_generated = values$report_generated,
        report_path = values$report_path,
        project_dirs = projectDirs
    )
}

lipidSummaryFirstNonNull <- function(...) {
    values <- list(...)
    for (value in values) {
        if (!is.null(value)) {
            return(value)
        }
    }
    NULL
}

getLipidSummaryWorkflowValue <- function(workflowData, name) {
    if (is.null(workflowData)) {
        return(NULL)
    }

    tryCatch(workflowData[[name]], error = function(e) NULL)
}

getLipidSummaryStateObject <- function(workflowData) {
    stateManager <- getLipidSummaryWorkflowValue(workflowData, "state_manager")
    if (is.null(stateManager)) {
        return(NULL)
    }

    dataStates <- c(
        "lipid_correlation_filtered",
        "lipid_norm_complete",
        "lipid_ruv_corrected",
        "lipid_normalized",
        "loaded_for_de",
        "lipid_qc_complete"
    )

    history <- tryCatch(stateManager$getHistory(), error = function(e) character())
    matchingStates <- dataStates[dataStates %in% history]
    stateName <- if (length(matchingStates) > 0L) matchingStates[[1]] else NULL
    if (is.null(stateName) || is.na(stateName)) {
        stateName <- tryCatch(stateManager$current_state, error = function(e) NULL)
    }
    if (is.null(stateName) || is.na(stateName)) {
        return(NULL)
    }

    tryCatch(stateManager$getState(stateName), error = function(e) NULL)
}

getLipidSummaryObjectArgs <- function(object) {
    if (is.null(object) || !isS4(object) || !"args" %in% methods::slotNames(object)) {
        return(NULL)
    }

    tryCatch(object@args, error = function(e) NULL)
}

buildLipidSummaryParameterPayload <- function(workflowData, finalS4Object = NULL) {
    payload <- list(
        config_list = getLipidSummaryWorkflowValue(workflowData, "config_list"),
        contrasts_tbl = getLipidSummaryWorkflowValue(workflowData, "contrasts_tbl"),
        da_ui_params = getLipidSummaryWorkflowValue(workflowData, "da_ui_params"),
        normalization_ui_params = getLipidSummaryWorkflowValue(workflowData, "normalization_ui_params"),
        itsd_ui_params = getLipidSummaryWorkflowValue(workflowData, "itsd_ui_params"),
        ruv_optimization_result = getLipidSummaryWorkflowValue(workflowData, "ruv_optimization_result"),
        s4_args = getLipidSummaryObjectArgs(finalS4Object)
    )

    payload[!vapply(payload, is.null, logical(1))]
}

prepareLipidSummarySessionState <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    timestamp,
    workflowData = NULL
) {
    sessionState <- list(
        experiment_label = inputValues$experiment_label,
        description = inputValues$description,
        timestamp = timestamp,
        omic_type = omicType,
        workflow_args_saved = values$workflow_args_saved,
        files_copied = values$files_copied,
        report_generated = values$report_generated,
        report_path = values$report_path,
        project_dirs = projectDirs
    )

    if (is.null(workflowData)) {
        return(sessionState)
    }

    finalS4Object <- getLipidSummaryStateObject(workflowData)
    parameterPayload <- buildLipidSummaryParameterPayload(workflowData, finalS4Object)
    globalParameters <- lipidSummaryFirstNonNull(
        parameterPayload$config_list$globalParameters,
        parameterPayload$s4_args$globalParameters,
        list()
    )

    sessionState$workflow_type <- lipidSummaryFirstNonNull(
        globalParameters$workflow_type,
        "lipidomics"
    )
    sessionState$report_template <- lipidSummaryFirstNonNull(
        globalParameters$report_template,
        "lipidomics_report.rmd"
    )
    sessionState$parameter_payload <- parameterPayload

    if (!is.null(finalS4Object) && isS4(finalS4Object)) {
        slotNames <- methods::slotNames(finalS4Object)
        if ("lipid_data" %in% slotNames) {
            assayData <- tryCatch(finalS4Object@lipid_data, error = function(e) NULL)
            if (is.list(assayData)) {
                sessionState$assay_names <- names(assayData)
                sessionState$feature_counts <- vapply(assayData, nrow, integer(1))
            }
        }
        if (all(c("design_matrix", "sample_id") %in% slotNames)) {
            designMatrix <- tryCatch(finalS4Object@design_matrix, error = function(e) NULL)
            sampleIdCol <- tryCatch(finalS4Object@sample_id, error = function(e) NULL)
            if (is.data.frame(designMatrix) && length(sampleIdCol) == 1L && sampleIdCol %in% names(designMatrix)) {
                sessionState$sample_count <- length(unique(designMatrix[[sampleIdCol]]))
            }
        }
    }

    contrastsTbl <- getLipidSummaryWorkflowValue(workflowData, "contrasts_tbl")
    if (is.data.frame(contrastsTbl)) {
        sessionState$contrast_count <- nrow(contrastsTbl)
    }

    sessionState
}

classifyLipidSummaryCopyFailures <- function(
    copyFailures,
    requiredDisplayNames = c("Contrasts Table", "Design Matrix", "Study Parameters")
) {
    if (is.null(copyFailures) || length(copyFailures) == 0L || isTRUE(copyFailures)) {
        return(list(required = list(), optional = list()))
    }

    isRequired <- vapply(copyFailures, function(failure) {
        !is.null(failure$display_name) && failure$display_name %in% requiredDisplayNames
    }, logical(1))

    list(
        required = copyFailures[isRequired],
        optional = copyFailures[!isRequired]
    )
}

handleLipidSummaryExportSessionState <- function(
    input,
    values,
    projectDirs,
    omicType,
    workflowData = NULL,
    saveRDSFn = saveRDS,
    showNotificationFn = shiny::showNotification,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    dateFn = Sys.Date,
    buildSessionStateFn = buildLipidSummarySessionState
) {
    tryCatch({
        sessionExportPath <- file.path(
            projectDirs[[omicType]]$source_dir,
            paste0("session_state_", dateFn(), ".RDS")
        )

        sessionState <- buildSessionStateFn(
            input = input,
            values = values,
            projectDirs = projectDirs,
            omicType = omicType
        )
        if (!is.null(workflowData) && identical(buildSessionStateFn, buildLipidSummarySessionState)) {
            sessionState <- prepareLipidSummarySessionState(
                inputValues = input,
                projectDirs = projectDirs,
                omicType = omicType,
                values = values,
                timestamp = Sys.time(),
                workflowData = workflowData
            )
        }

        saveRDSFn(sessionState, sessionExportPath)

        showNotificationFn(
            paste("Session state exported to:", sessionExportPath),
            type = "message"
        )
        logInfoFn(paste("Session state exported to:", sessionExportPath))

        sessionExportPath
    }, error = function(e) {
        showNotificationFn(
            paste("Export failed:", e$message),
            type = "error",
            duration = 10
        )
        logErrorFn(skipFormatterFn(paste("Failed to export session state:", e$message)))
        NULL
    })
}

registerLipidSummaryExportSessionObserver <- function(
    input,
    values,
    projectDirs,
    omicType,
    workflowData = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handleExportFn = handleLipidSummaryExportSessionState
) {
    observeEventFn(input$export_session_state, {
        reqFn(input$experiment_label)
        handleExportFn(
            input = input,
            values = values,
            projectDirs = projectDirs,
            omicType = omicType,
            workflowData = workflowData
        )
    })

    invisible(input)
}

collectLipidSummaryWorkflowArgsContext <- function(
    workflowData,
    detectFn = purrr::detect,
    assignFn = assign,
    globalEnv = .GlobalEnv,
    catFn = cat
) {
    finalS4Object <- NULL

    if (!is.null(workflowData$state_manager)) {
        dataStates <- c(
            "lipid_correlation_filtered",
            "lipid_norm_complete",
            "lipid_ruv_corrected",
            "lipid_normalized",
            "loaded_for_de",
            "lipid_qc_complete"
        )

        availableStates <- workflowData$state_manager$getHistory()
        dataStateUsed <- detectFn(dataStates, ~ .x %in% availableStates)

        if (!is.null(dataStateUsed)) {
            finalS4Object <- workflowData$state_manager$getState(dataStateUsed)
        }

        if (!is.null(finalS4Object)) {
            catFn(sprintf("SESSION SUMMARY: Retrieved DATA S4 object from state '%s'\n", dataStateUsed))
            catFn(sprintf("SESSION SUMMARY: S4 object class: %s\n", class(finalS4Object)))

            argsAvailable <- tryCatch({
                !is.null(finalS4Object@args)
            }, error = function(e) {
                FALSE
            })

            if (argsAvailable) {
                catFn(sprintf("SESSION SUMMARY: S4 @args contains %d function groups\n", length(finalS4Object@args)))
            } else {
                catFn("SESSION SUMMARY: S4 @args is NULL or slot doesn't exist\n")
            }
        } else {
            catFn("SESSION SUMMARY: No data S4 object found in any valid state\n")
        }
    } else {
        catFn("SESSION SUMMARY: No state manager available\n")
    }

    contrastsTbl <- NULL
    if (!is.null(workflowData) && !is.null(workflowData$contrasts_tbl)) {
        contrastsTbl <- workflowData$contrasts_tbl
        catFn("SESSION SUMMARY: Using contrasts_tbl from workflow_data\n")
    }

    if (!is.null(workflowData) && !is.null(workflowData$config_list)) {
        assignFn("config_list", workflowData$config_list, envir = globalEnv)
        catFn("SESSION SUMMARY: Config list available with", length(workflowData$config_list), "items\n")
    }

    list(
        finalS4Object = finalS4Object,
        contrastsTbl = contrastsTbl
    )
}
