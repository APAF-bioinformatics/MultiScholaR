initializeLipidSummarySessionBootstrap <- function(
    session,
    experimentLabel = NULL,
    updateTextInputFn = shiny::updateTextInput,
    reactiveValuesFn = shiny::reactiveValues
) {
    if (!is.null(experimentLabel)) {
        updateTextInputFn(session, "experiment_label", value = experimentLabel)
    }

    reactiveValuesFn(
        workflow_args_saved = FALSE,
        files_copied = FALSE,
        report_generated = FALSE,
        report_path = NULL
    )
}

registerLipidSummaryInitialOutputs <- function(
    output,
    renderTextFn = shiny::renderText,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions
) {
    output$session_summary <- renderTextFn({
        "Ready to save workflow parameters and generate report"
    })

    output$report_ready <- reactiveFn({ FALSE })
    outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)

    output
}

runLipidSummarySaveWorkflowArgsObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    workflowData,
    values,
    output,
    reqFn = shiny::req,
    handleSaveFn = handleLipidSummarySaveWorkflowArgs,
    ...
) {
    reqFn(inputValues$experiment_label)
    result <- handleSaveFn(
        input = inputValues,
        values = values,
        output = output,
        projectDirs = projectDirs,
        omicType = omicType,
        workflowData = workflowData,
        ...
    )

    if (is.list(result)) {
        c(list(status = "success"), result)
    } else {
        list(status = "fallback", studyParamsFile = result)
    }
}

runLipidSummaryCopyToPublicationObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    workflowData,
    values,
    output,
    reqFn = shiny::req,
    handleCopyFn = handleLipidSummaryCopyToPublication,
    ...
) {
    reqFn(inputValues$experiment_label)
    result <- handleCopyFn(
        input = inputValues,
        values = values,
        output = output,
        projectDirs = projectDirs,
        omicType = omicType,
        workflowData = workflowData,
        ...
    )

    if (is.null(result)) {
        errorMessage <- tryCatch(as.character(output$copy_status), error = function(e) "Copy failed")
        return(list(status = "error", errorMessage = errorMessage))
    }

    c(list(status = "success"), result)
}

runLipidSummaryGenerateReportObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    output,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    incProgressFn = shiny::incProgress,
    showNotificationFn = shiny::showNotification,
    dirCreateFn = dir.create,
    fileExistsFn = file.exists,
    systemFileFn = system.file,
    fileCopyFn = file.copy,
    downloadFileFn = utils::download.file,
    renderReportFn = NULL,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions,
    downloadHandlerFn = shiny::downloadHandler,
    renderTextFn = shiny::renderText,
    sysTimeFn = Sys.time,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    catFn = cat,
    printFn = print
) {
    reqFn(inputValues$experiment_label)
    reqFn(values$files_copied)

    if (!omicType %in% names(projectDirs) || is.null(projectDirs[[omicType]]$base_dir)) {
        showNotificationFn(
            "Error: Project directories not properly initialized",
            type = "error",
            duration = 10
        )
        return(list(status = "invalid_project_dirs", errorMessage = "Project directories not properly initialized"))
    }

    withProgressFn(message = "Generating report...", {
        templateFilename <- "lipidomics_report.rmd"
        reportTemplatePath <- file.path(
            projectDirs[[omicType]]$base_dir,
            "scripts",
            omicType,
            templateFilename
        )
        templateSource <- NULL

        if (!fileExistsFn(reportTemplatePath)) {
            incProgressFn(0.2, detail = paste("Checking for", templateFilename, "template..."))
            dirCreateFn(dirname(reportTemplatePath), recursive = TRUE, showWarnings = FALSE)
            templateSource <- tryCatch({
                pkgFile <- systemFileFn("reports", "lipidomics", templateFilename, package = "MultiScholaR")
                if (fileExistsFn(pkgFile) && pkgFile != "") {
                    fileCopyFn(pkgFile, reportTemplatePath)
                    showNotificationFn(paste(templateFilename, "template copied from package"), type = "message")
                    "package"
                } else {
                    templateUrl <- paste0(
                        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/lipidomics/",
                        templateFilename
                    )
                    downloadFileFn(templateUrl, destfile = reportTemplatePath, quiet = TRUE)
                    showNotificationFn(paste(templateFilename, "template downloaded"), type = "message")
                    "github"
                }
            }, error = function(e) {
                errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
                showNotificationFn(paste("Template retrieval failed:", errorMessage), type = "error", duration = 10)
                logErrorFn(skipFormatterFn(paste("Failed to retrieve template:", errorMessage)))
                NULL
            })
        } else {
            templateSource <- "existing"
        }

        incProgressFn(0.5, detail = "Rendering report...")
        tryCatch({
            if (is.null(renderReportFn)) {
                if (!exists("RenderReport", inherits = TRUE)) {
                    showNotificationFn(
                        "Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.",
                        type = "error",
                        duration = 15
                    )
                    return(list(
                        status = "missing_render_report",
                        templateFilename = templateFilename,
                        reportTemplatePath = reportTemplatePath,
                        templateSource = templateSource
                    ))
                }
                renderReportFn <- get("RenderReport", inherits = TRUE)
            }

            renderedPath <- renderReportFn(
                omic_type = omicType,
                experiment_label = inputValues$experiment_label,
                rmd_filename = templateFilename
            )

            if (!is.null(renderedPath) && fileExistsFn(renderedPath)) {
                values$report_generated <- TRUE
                values$report_path <- renderedPath
                output$report_ready <- reactiveFn({ TRUE })
                outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)
                output$download_report <- downloadHandlerFn(
                    filename = function() basename(renderedPath),
                    content = function(file) fileCopyFn(renderedPath, file)
                )
                showNotificationFn("Report generated successfully!", type = "message")
                output$session_summary <- renderTextFn({
                    paste(
                        "Workflow args created for:", inputValues$experiment_label,
                        "\nDescription:", inputValues$description,
                        "\nTimestamp:", sysTimeFn(),
                        "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
                        "\nReport location:", renderedPath
                    )
                })
                list(
                    status = "success",
                    templateFilename = templateFilename,
                    reportTemplatePath = reportTemplatePath,
                    templateSource = templateSource,
                    renderedPath = renderedPath
                )
            } else {
                showNotificationFn("Report generation failed - no output file created", type = "error", duration = 10)
                list(
                    status = "missing_output",
                    templateFilename = templateFilename,
                    reportTemplatePath = reportTemplatePath,
                    templateSource = templateSource,
                    renderedPath = renderedPath
                )
            }
        }, error = function(e) {
            errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
            catFn("DEBUG66: REPORT GENERATION ERROR\n")
            catFn(sprintf("Error class: %s\n", class(e)[1]))
            catFn(sprintf("Error message: %s\n", errorMessage))
            catFn("Full error object:\n")
            printFn(e)
            showNotificationFn(paste("Report generation failed:", errorMessage), type = "error", duration = 15)
            list(
                status = "error",
                templateFilename = templateFilename,
                reportTemplatePath = reportTemplatePath,
                templateSource = templateSource,
                errorMessage = errorMessage
            )
        })
    })
}

runLipidSummaryExportSessionObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    workflowData = NULL,
    reqFn = shiny::req,
    saveRdsFn = saveRDS,
    sysDateFn = Sys.Date,
    sysTimeFn = Sys.time,
    showNotificationFn = shiny::showNotification,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter
) {
    reqFn(inputValues$experiment_label)

    tryCatch({
        sessionExportPath <- file.path(
            projectDirs[[omicType]]$source_dir,
            paste0("session_state_", sysDateFn(), ".RDS")
        )
        sessionState <- prepareLipidSummarySessionState(
            inputValues = inputValues,
            projectDirs = projectDirs,
            omicType = omicType,
            values = values,
            timestamp = sysTimeFn(),
            workflowData = workflowData
        )
        saveRdsFn(sessionState, sessionExportPath)
        showNotificationFn(paste("Session state exported to:", sessionExportPath), type = "message")
        logInfoFn(paste("Session state exported to:", sessionExportPath))

        list(
            status = "success",
            sessionExportPath = sessionExportPath,
            sessionState = sessionState
        )
    }, error = function(e) {
        errorMessage <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
        showNotificationFn(paste("Export failed:", errorMessage), type = "error", duration = 10)
        logErrorFn(skipFormatterFn(paste("Failed to export session state:", errorMessage)))
        list(status = "error", errorMessage = errorMessage)
    })
}

