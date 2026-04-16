# ============================================================================
# mod_lipid_summary_server_helpers.R
# ============================================================================
# Purpose: Extracted helper surface for the lipid summary module server wrapper
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

handleLipidSummaryExportSessionState <- function(
    input,
    values,
    projectDirs,
    omicType,
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
            omicType = omicType
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

handleLipidSummarySaveWorkflowArgs <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    workflowData = NULL,
    collectContextFn = collectLipidSummaryWorkflowArgsContext,
    createWorkflowArgsFn = createWorkflowArgsFromConfig,
    saveRDSFn = saveRDS,
    showNotificationFn = shiny::showNotification,
    renderTextFn = shiny::renderText,
    fileExistsFn = file.exists,
    dirExistsFn = dir.exists,
    dirCreateFn = dir.create,
    writeLinesFn = writeLines,
    filePathFn = file.path,
    timeFn = Sys.time,
    catFn = cat
) {
    catFn("SESSION SUMMARY: Starting workflow args save process\n")

    tryCatch({
        context <- collectContextFn(workflowData = workflowData, catFn = catFn)

        catFn("SESSION SUMMARY: Creating study_parameters.txt file with S4 parameters\n")
        studyParamsFile <- createWorkflowArgsFn(
            workflow_name = input$experiment_label,
            description = input$description,
            source_dir_path = projectDirs[[omicType]]$source_dir,
            final_s4_object = context$finalS4Object,
            contrasts_tbl = context$contrastsTbl,
            workflow_data = workflowData
        )

        catFn("SESSION SUMMARY: Successfully created study_parameters.txt at:", studyParamsFile, "\n")
        catFn("SESSION SUMMARY: Saving Integration S4 Object...\n")

        integrationDir <- projectDirs[[omicType]]$integration_dir
        if (is.null(integrationDir)) {
            integrationDir <- filePathFn(projectDirs[[omicType]]$base_dir, "integration")
        }

        if (!dirExistsFn(integrationDir)) {
            dirCreateFn(integrationDir, recursive = TRUE, showWarnings = FALSE)
        }

        s4Filename <- sprintf("%s_%s_final_s4.RDS", omicType, input$experiment_label)
        s4Filepath <- filePathFn(integrationDir, s4Filename)

        saveRDSFn(context$finalS4Object, s4Filepath)
        catFn(sprintf("SESSION SUMMARY: Saved Integration S4 object to: %s\n", s4Filepath))
        showNotificationFn("Saved Integration S4 Object", type = "message")

        values$workflow_args_saved <- TRUE
        showNotificationFn("Study parameters saved successfully", type = "message")

        output$session_summary <- renderTextFn({
            paste(
                "Study parameters created for:", input$experiment_label,
                "\nDescription:", input$description,
                "\nTimestamp:", timeFn(),
                "\nFile:", studyParamsFile,
                "\nSource: Final S4 object @args + config_list",
                "\nIntegration Object:", s4Filename,
                "\nStatus: Parameters saved [OK]"
            )
        })

        list(
            studyParamsFile = studyParamsFile,
            s4Filepath = s4Filepath,
            s4Filename = s4Filename
        )
    }, error = function(e) {
        catFn("SESSION SUMMARY ERROR:", e$message, "\n")

        basicParamsFile <- filePathFn(projectDirs[[omicType]]$source_dir, "study_parameters.txt")
        if (!fileExistsFn(basicParamsFile)) {
            writeLinesFn(c(
                "Study Parameters",
                "================",
                "",
                paste("Workflow Name:", input$experiment_label),
                paste("Description:", input$description),
                paste("Timestamp:", timeFn()),
                paste("Error:", e$message)
            ), basicParamsFile)
            catFn("SESSION SUMMARY: Created basic fallback file at:", basicParamsFile, "\n")
        }

        values$workflow_args_saved <- TRUE
        showNotificationFn("Study parameters saved with warnings", type = "warning")

        basicParamsFile
    })
}

registerLipidSummarySaveWorkflowArgsObserver <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    workflowData = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handleSaveFn = handleLipidSummarySaveWorkflowArgs
) {
    observeEventFn(input$save_workflow_args, {
        reqFn(input$experiment_label)
        handleSaveFn(
            input = input,
            values = values,
            output = output,
            projectDirs = projectDirs,
            omicType = omicType,
            workflowData = workflowData
        )
    })

    invisible(input)
}

handleLipidSummaryCopyToPublication <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    workflowData = NULL,
    withProgressFn = shiny::withProgress,
    filePathFn = file.path,
    fileExistsFn = file.exists,
    writeLinesFn = writeLines,
    timeFn = Sys.time,
    readTsvFn = readr::read_tsv,
    existsFn = exists,
    assignFn = assign,
    globalEnv = .GlobalEnv,
    copyResultsSummaryFn = copyToResultsSummary,
    showNotificationFn = shiny::showNotification,
    renderTextFn = shiny::renderText,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    catFn = cat,
    tracebackFn = traceback
) {
    if (!values$workflow_args_saved) {
        basicParamsFile <- filePathFn(projectDirs[[omicType]]$source_dir, "study_parameters.txt")
        if (!fileExistsFn(basicParamsFile)) {
            catFn("SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")
            basicContent <- paste(
                "Study Parameters",
                "================",
                "",
                paste("Workflow Name:", input$experiment_label),
                paste("Description:", input$description),
                paste("Timestamp:", timeFn()),
                paste("Note: Some parameters could not be saved due to serialization issues"),
                sep = "\n"
            )
            tryCatch({
                writeLinesFn(basicContent, basicParamsFile)
                values$workflow_args_saved <- TRUE
            }, error = function(e) {
                logErrorFn(skipFormatterFn(paste("Failed to create basic study_parameters.txt:", e$message)))
            })
        }
    }

    catFn("SESSION SUMMARY: Copy to Publication button clicked\n")
    catFn("SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
    catFn("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")

    withProgressFn(message = "Copying files to publication directory...", {
        tryCatch({
            contrastsTbl <- NULL
            designMatrix <- NULL

            if (!is.null(workflowData)) {
                if (!is.null(workflowData$contrasts_tbl)) {
                    contrastsTbl <- workflowData$contrasts_tbl
                    catFn("SESSION SUMMARY: Got contrasts_tbl from workflow_data\n")
                }
                if (!is.null(workflowData$design_matrix)) {
                    designMatrix <- workflowData$design_matrix
                    catFn("SESSION SUMMARY: Got design_matrix from workflow_data\n")
                }
            }

            if (is.null(designMatrix)) {
                designMatrixFile <- filePathFn(projectDirs[[omicType]]$source_dir, "design_matrix.tab")
                if (fileExistsFn(designMatrixFile)) {
                    designMatrix <- readTsvFn(designMatrixFile, show_col_types = FALSE)
                    catFn("SESSION SUMMARY: Loaded design_matrix from file\n")
                }
            }

            if (is.null(contrastsTbl)) {
                contrastsFile <- filePathFn(projectDirs[[omicType]]$source_dir, "contrasts_tbl.tab")
                if (fileExistsFn(contrastsFile)) {
                    contrastsTbl <- readTsvFn(contrastsFile, show_col_types = FALSE)
                    catFn("SESSION SUMMARY: Loaded contrasts_tbl from file\n")
                }
            }

            catFn("SESSION SUMMARY: About to call copyToResultsSummary\n")
            catFn("SESSION SUMMARY: omic_type =", omicType, "\n")
            catFn("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
            catFn("SESSION SUMMARY: contrasts_tbl is", ifelse(is.null(contrastsTbl), "NULL", "available"), "\n")
            catFn("SESSION SUMMARY: design_matrix is", ifelse(is.null(designMatrix), "NULL", "available"), "\n")

            if (!existsFn("project_dirs", envir = globalEnv)) {
                catFn("SESSION SUMMARY: Setting project_dirs in global environment\n")
                assignFn("project_dirs", projectDirs, envir = globalEnv)
            }

            copyResultsSummaryFn(
                omic_type = omicType,
                experiment_label = input$experiment_label,
                contrasts_tbl = contrastsTbl,
                design_matrix = designMatrix,
                force = TRUE
            )

            values$files_copied <- TRUE

            output$copy_status <- renderTextFn("Files copied to publication directory successfully [OK]")
            showNotificationFn("Publication files copied", type = "message")

            output$session_summary <- renderTextFn({
                paste(
                    "Workflow args created for:", input$experiment_label,
                    "\nDescription:", input$description,
                    "\nTimestamp:", timeFn(),
                    "\nStatus: Arguments saved [OK], Files copied [OK]"
                )
            })

            list(
                contrastsTbl = contrastsTbl,
                designMatrix = designMatrix
            )
        }, error = function(e) {
            output$copy_status <- renderTextFn(paste("Error:", e$message))
            showNotificationFn(paste("Copy error:", e$message), type = "error", duration = 10)
            logErrorFn(skipFormatterFn(paste("Failed to copy files:", e$message)))
            catFn("SESSION SUMMARY ERROR:", e$message, "\n")
            catFn("SESSION SUMMARY Traceback:\n")
            tracebackFn()
            NULL
        })
    })
}

registerLipidSummaryCopyToPublicationObserver <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    workflowData = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handleCopyFn = handleLipidSummaryCopyToPublication
) {
    observeEventFn(input$copy_to_publication, {
        reqFn(input$experiment_label)
        handleCopyFn(
            input = input,
            values = values,
            output = output,
            projectDirs = projectDirs,
            omicType = omicType,
            workflowData = workflowData
        )
    }, ignoreInit = TRUE)

    invisible(input)
}

handleLipidSummaryGenerateReport <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    templateFilename = "lipidomics_report.rmd",
    withProgressFn = shiny::withProgress,
    showNotificationFn = shiny::showNotification,
    incProgressFn = shiny::incProgress,
    filePathFn = file.path,
    fileExistsFn = file.exists,
    dirCreateFn = dir.create,
    dirnameFn = dirname,
    systemFileFn = system.file,
    fileCopyFn = file.copy,
    downloadFileFn = download.file,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    renderReportExistsFn = function() exists("RenderReport"),
    renderReportFn = RenderReport,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions,
    downloadHandlerFn = shiny::downloadHandler,
    renderTextFn = shiny::renderText,
    timeFn = Sys.time,
    basenameFn = basename,
    catFn = cat,
    printFn = print
) {
    if (!omicType %in% names(projectDirs) || is.null(projectDirs[[omicType]]$base_dir)) {
        showNotificationFn(
            "Error: Project directories not properly initialized",
            type = "error",
            duration = 10
        )
        return(invisible(NULL))
    }

    withProgressFn(message = "Generating report...", {
        catFn(sprintf("REPORT: Selected template: %s\n", templateFilename))

        reportTemplatePath <- filePathFn(
            projectDirs[[omicType]]$base_dir,
            "scripts",
            omicType,
            templateFilename
        )

        if (!fileExistsFn(reportTemplatePath)) {
            incProgressFn(0.2, detail = paste("Checking for", templateFilename, "template..."))
            dirCreateFn(dirnameFn(reportTemplatePath), recursive = TRUE, showWarnings = FALSE)

            templateReady <- tryCatch({
                pkgFile <- systemFileFn("reports", "lipidomics", templateFilename, package = "MultiScholaR")

                if (fileExistsFn(pkgFile) && pkgFile != "") {
                    catFn(sprintf("REPORT: Found template in package at: %s\n", pkgFile))
                    fileCopyFn(pkgFile, reportTemplatePath)
                    logInfoFn("Template copied from package to: {reportTemplatePath}")
                    showNotificationFn(paste(templateFilename, "template copied from package"), type = "message")
                } else {
                    catFn("REPORT: Template not found in package, downloading from GitHub...\n")

                    templateUrl <- paste0(
                        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/lipidomics/",
                        templateFilename
                    )

                    catFn(sprintf("REPORT: Downloading template from: %s\n", templateUrl))
                    downloadFileFn(templateUrl, destfile = reportTemplatePath, quiet = TRUE)
                    logInfoFn("Template downloaded to: {reportTemplatePath}")
                    showNotificationFn(paste(templateFilename, "template downloaded"), type = "message")
                }

                TRUE
            }, error = function(e) {
                showNotificationFn(
                    paste("Template retrieval failed:", e$message),
                    type = "error",
                    duration = 10
                )
                logErrorFn(skipFormatterFn(paste("Failed to retrieve template:", e$message)))
                FALSE
            })

            if (!isTRUE(templateReady)) {
                return(invisible(NULL))
            }
        } else {
            catFn(sprintf("REPORT: Template already exists at: %s\n", reportTemplatePath))
        }

        incProgressFn(0.5, detail = "Rendering report...")

        tryCatch({
            if (!renderReportExistsFn()) {
                showNotificationFn(
                    "Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.",
                    type = "error",
                    duration = 15
                )
                return(invisible(NULL))
            }

            logInfoFn("Calling RenderReport with omic_type: {omicType}, experiment_label: {input$experiment_label}")

            renderedPath <- renderReportFn(
                omic_type = omicType,
                experiment_label = input$experiment_label,
                rmd_filename = templateFilename
            )

            logInfoFn("RenderReport returned path: {renderedPath}")

            if (!is.null(renderedPath) && fileExistsFn(renderedPath)) {
                values$report_generated <- TRUE
                values$report_path <- renderedPath

                output$report_ready <- reactiveFn({ TRUE })
                outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)

                output$download_report <- downloadHandlerFn(
                    filename = function() basenameFn(renderedPath),
                    content = function(file) fileCopyFn(renderedPath, file)
                )

                showNotificationFn("Report generated successfully!", type = "message")

                output$session_summary <- renderTextFn({
                    paste(
                        "Workflow args created for:", input$experiment_label,
                        "\nDescription:", input$description,
                        "\nTimestamp:", timeFn(),
                        "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
                        "\nReport location:", renderedPath
                    )
                })
            } else {
                showNotificationFn(
                    "Report generation failed - no output file created",
                    type = "error",
                    duration = 10
                )
            }
        }, error = function(e) {
            errorMsg <- paste("Report generation failed:", e$message)

            catFn("DEBUG66: REPORT GENERATION ERROR\n")
            catFn(sprintf("Error class: %s\n", class(e)[1]))
            catFn(sprintf("Error message: %s\n", e$message))
            catFn("Full error object:\n")
            printFn(e)

            logErrorFn("Failed to generate report: {e$message}")
            logErrorFn("Error class: {class(e)[1]}")

            showNotificationFn(errorMsg, type = "error", duration = 15)
            showNotificationFn(
                "Debug info: Check R console for detailed error trace",
                type = "warning",
                duration = 10
            )
        })
    })

    invisible(NULL)
}

registerLipidSummaryGenerateReportObserver <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handleGenerateReportFn = handleLipidSummaryGenerateReport
) {
    observeEventFn(input$generate_report, {
        reqFn(input$experiment_label)
        reqFn(values$files_copied)
        handleGenerateReportFn(
            input = input,
            values = values,
            output = output,
            projectDirs = projectDirs,
            omicType = omicType
        )
    })

    invisible(input)
}

handleLipidSummaryPushToGithub <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    withProgressFn = shiny::withProgress,
    optionsFn = options,
    pushProjectFn = pushProjectToGithubFromDirs,
    showNotificationFn = shiny::showNotification,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    renderTextFn = shiny::renderText,
    timeFn = Sys.time
) {
    withProgressFn(message = "Pushing to GitHub...", {
        tryCatch({
            optionsFn(
                github_org = input$github_org,
                github_user_email = input$github_email,
                github_user_name = input$github_username
            )

            pushProjectFn(
                project_dirs = projectDirs,
                omic_type = omicType,
                experiment_label = input$experiment_label,
                project_id = input$project_id
            )

            showNotificationFn("Successfully pushed to GitHub", type = "message")

            output$session_summary <- renderTextFn({
                paste(
                    "Workflow args created for:", input$experiment_label,
                    "\nDescription:", input$description,
                    "\nTimestamp:", timeFn(),
                    "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]"
                )
            })

            invisible(TRUE)
        }, error = function(e) {
            showNotificationFn(
                paste("GitHub push failed:", e$message),
                type = "error",
                duration = 10
            )
            logErrorFn(skipFormatterFn(paste("Failed to push to GitHub:", e$message)))
            invisible(NULL)
        })
    })
}

registerLipidSummaryPushToGithubObserver <- function(
    input,
    values,
    output,
    projectDirs,
    omicType,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    handlePushFn = handleLipidSummaryPushToGithub
) {
    observeEventFn(input$push_to_github, {
        reqFn(
            input$enable_github,
            input$github_org,
            input$github_email,
            input$github_username,
            input$project_id
        )
        reqFn(values$report_generated)

        handlePushFn(
            input = input,
            values = values,
            output = output,
            projectDirs = projectDirs,
            omicType = omicType
        )
    })

    invisible(input)
}

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
