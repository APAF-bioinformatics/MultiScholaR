# ============================================================================
# mod_lipid_summary_state_helpers.R
# ============================================================================
# Purpose: Session state and parameter helpers for the lipid summary module
# ============================================================================

buildLipidSummaryTemplateStatus <- function(projectDirs, omicType, fileExistsFn = file.exists, profile = NULL) {
    if (!omicType %in% names(projectDirs)) {
        return("[WARNING] Project directories not available")
    }

    template_dir <- file.path(
        projectDirs[[omicType]]$base_dir,
        "scripts",
        omicType
    )
    lipid_template <- file.path(template_dir, profile$report_template %||% "lipidomics_report.rmd")

    if (fileExistsFn(lipid_template)) {
        sprintf("Template: %s Report [OK]", profile$label %||% "Lipidomics")
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
    buildSessionStateFn = buildLipidSummarySessionState,
    profile = NULL
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
                workflowData = workflowData,
                profile = profile
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

        # KIT-071: the final state is the head's (the owner's current state when it is a data state).
        dataStateUsed <- workflowSummaryFinalStateName(workflowData$state_manager, dataStates, detectFn = detectFn)

        if (!is.null(dataStateUsed)) {
            # KIT-082: on a managed session the worker writes the integration object, and the parameters here read the id column.
            finalS4Object <- if (!is.null(workflowOmicsManagedActions("lipidomics"))) readLipidWorkflowIds(workflowData$state_manager, dataStateUsed)
                else workflowData$state_manager$getState(dataStateUsed)
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
        # KIT-030: the configuration is read from workflow_data (createWorkflowArgsFromConfig takes it); no global copy.
        catFn("SESSION SUMMARY: Config list available with", length(workflowData$config_list), "items\n")
    }

    list(
        finalS4Object = finalS4Object,
        contrastsTbl = contrastsTbl,
        dataStateUsed = if (!is.null(workflowData$state_manager)) dataStateUsed
    )
}
