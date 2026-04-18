runMetabSummaryExportSessionObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
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
        session_export_path <- file.path(
            projectDirs[[omicType]]$source_dir,
            paste0("session_state_", sysDateFn(), ".RDS")
        )

        session_state <- list(
            experiment_label = inputValues$experiment_label,
            description = inputValues$description,
            timestamp = sysTimeFn(),
            omic_type = omicType,
            workflow_args_saved = values$workflow_args_saved,
            files_copied = values$files_copied,
            report_generated = values$report_generated,
            report_path = values$report_path,
            project_dirs = projectDirs
        )

        saveRdsFn(session_state, session_export_path)

        showNotificationFn(
            paste("Session state exported to:", session_export_path),
            type = "message"
        )
        logInfoFn(paste("Session state exported to:", session_export_path))

        invisible(list(
            status = "success",
            sessionExportPath = session_export_path,
            sessionState = session_state
        ))
    }, error = function(e) {
        error_message <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

        showNotificationFn(
            paste("Export failed:", error_message),
            type = "error",
            duration = 10
        )
        logErrorFn(
            skipFormatterFn(
                paste("Failed to export session state:", error_message)
            )
        )

        invisible(list(
            status = "error",
            errorMessage = error_message
        ))
    })
}

runMetabSummarySaveWorkflowArgsObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    workflowData,
    values,
    output,
    reqFn = shiny::req,
    detectFn = purrr::detect,
    assignFn = assign,
    createWorkflowArgsFromConfigFn = createWorkflowArgsFromConfig,
    dirExistsFn = dir.exists,
    dirCreateFn = dir.create,
    saveRdsFn = saveRDS,
    fileExistsFn = file.exists,
    writeLinesFn = writeLines,
    renderTextFn = shiny::renderText,
    showNotificationFn = shiny::showNotification,
    sysTimeFn = Sys.time,
    catFn = cat,
    globalEnv = .GlobalEnv
) {
    input <- inputValues
    project_dirs <- projectDirs
    omic_type <- omicType
    workflow_data <- workflowData

    reqFn(input$experiment_label)

    catFn("SESSION SUMMARY: Starting workflow args save process\n")

    tryCatch({
        final_s4_object <- NULL
        data_state_used <- NULL

        if (!is.null(workflow_data$state_manager)) {
            data_states <- c(
                "metab_correlation_filtered", "metab_norm_complete",
                "metab_ruv_corrected", "metab_normalized",
                "loaded_for_de", "metab_qc_complete"
            )

            available_states <- workflow_data$state_manager$getHistory()
            data_state_used <- detectFn(data_states, ~ .x %in% available_states)

            if (!is.null(data_state_used)) {
                final_s4_object <- workflow_data$state_manager$getState(data_state_used)
            }

            if (!is.null(final_s4_object)) {
                catFn(sprintf(
                    "SESSION SUMMARY: Retrieved DATA S4 object from state '%s'\n",
                    data_state_used
                ))
                catFn(sprintf(
                    "SESSION SUMMARY: S4 object class: %s\n",
                    class(final_s4_object)
                ))

                args_available <- tryCatch({
                    !is.null(final_s4_object@args)
                }, error = function(e) {
                    FALSE
                })

                if (args_available) {
                    catFn(sprintf(
                        "SESSION SUMMARY: S4 @args contains %d function groups\n",
                        length(final_s4_object@args)
                    ))
                } else {
                    catFn("SESSION SUMMARY: S4 @args is NULL or slot doesn't exist\n")
                }
            } else {
                catFn("SESSION SUMMARY: No data S4 object found in any valid state\n")
            }
        } else {
            catFn("SESSION SUMMARY: No state manager available\n")
        }

        contrasts_tbl <- NULL
        if (!is.null(workflow_data) && !is.null(workflow_data$contrasts_tbl)) {
            contrasts_tbl <- workflow_data$contrasts_tbl
            catFn("SESSION SUMMARY: Using contrasts_tbl from workflow_data\n")
        }

        config_list_assigned <- FALSE
        if (!is.null(workflow_data) && !is.null(workflow_data$config_list)) {
            assignFn("config_list", workflow_data$config_list, envir = globalEnv)
            config_list_assigned <- TRUE
            catFn(
                "SESSION SUMMARY: Config list available with",
                length(workflow_data$config_list),
                "items\n"
            )
        }

        catFn("SESSION SUMMARY: Creating study_parameters.txt file with S4 parameters\n")
        study_params_file <- createWorkflowArgsFromConfigFn(
            workflow_name = input$experiment_label,
            description = input$description,
            source_dir_path = project_dirs[[omic_type]]$source_dir,
            final_s4_object = final_s4_object,
            contrasts_tbl = contrasts_tbl,
            workflow_data = workflow_data
        )

        catFn(
            "SESSION SUMMARY: Successfully created study_parameters.txt at:",
            study_params_file,
            "\n"
        )

        catFn("SESSION SUMMARY: Saving Integration S4 Object...\n")
        integration_dir <- project_dirs[[omic_type]]$integration_dir
        if (is.null(integration_dir)) {
            integration_dir <- file.path(
                project_dirs[[omic_type]]$base_dir,
                "integration"
            )
        }

        if (!dirExistsFn(integration_dir)) {
            dirCreateFn(integration_dir, recursive = TRUE, showWarnings = FALSE)
        }

        s4_filename <- sprintf(
            "%s_%s_final_s4.RDS",
            omic_type,
            input$experiment_label
        )
        s4_filepath <- file.path(integration_dir, s4_filename)

        saveRdsFn(final_s4_object, s4_filepath)
        catFn(sprintf(
            "SESSION SUMMARY: Saved Integration S4 object to: %s\n",
            s4_filepath
        ))
        showNotificationFn("Saved Integration S4 Object", type = "message")

        values$workflow_args_saved <- TRUE
        showNotificationFn("Study parameters saved successfully", type = "message")

        output$session_summary <- renderTextFn({
            paste(
                "Study parameters created for:", input$experiment_label,
                "\nDescription:", input$description,
                "\nTimestamp:", sysTimeFn(),
                "\nFile:", study_params_file,
                "\nSource: Final S4 object @args + config_list",
                "\nIntegration Object:", s4_filename,
                "\nStatus: Parameters saved [OK]"
            )
        })

        invisible(list(
            status = "success",
            studyParamsFile = study_params_file,
            s4Filepath = s4_filepath,
            dataStateUsed = data_state_used,
            configListAssigned = config_list_assigned
        ))
    }, error = function(e) {
        error_message <- if (inherits(e, "condition")) {
            conditionMessage(e)
        } else {
            as.character(e)
        }

        catFn("SESSION SUMMARY ERROR:", error_message, "\n")

        basic_params_file <- file.path(
            project_dirs[[omic_type]]$source_dir,
            "study_parameters.txt"
        )
        fallback_created <- FALSE

        if (!fileExistsFn(basic_params_file)) {
            writeLinesFn(c(
                "Study Parameters",
                "================",
                "",
                paste("Workflow Name:", input$experiment_label),
                paste("Description:", input$description),
                paste("Timestamp:", sysTimeFn()),
                paste("Error:", error_message)
            ), basic_params_file)
            catFn(
                "SESSION SUMMARY: Created basic fallback file at:",
                basic_params_file,
                "\n"
            )
            fallback_created <- TRUE
        }

        values$workflow_args_saved <- TRUE
        showNotificationFn(
            "Study parameters saved with warnings",
            type = "warning"
        )

        invisible(list(
            status = "warning",
            basicParamsFile = basic_params_file,
            fallbackCreated = fallback_created,
            errorMessage = error_message
        ))
    })
}

runMetabSummaryCopyToPublicationObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    workflowData,
    values,
    output,
    reqFn = shiny::req,
    fileExistsFn = file.exists,
    sysTimeFn = Sys.time,
    writeLinesFn = writeLines,
    withProgressFn = shiny::withProgress,
    readTsvFn = readr::read_tsv,
    existsFn = exists,
    assignFn = assign,
    copyToResultsSummaryFn = copyToResultsSummary,
    renderTextFn = shiny::renderText,
    showNotificationFn = shiny::showNotification,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter,
    catFn = cat,
    tracebackFn = traceback,
    globalEnv = .GlobalEnv
) {
    input <- inputValues
    project_dirs <- projectDirs
    omic_type <- omicType
    workflow_data <- workflowData

    reqFn(input$experiment_label)

    basic_params_file <- file.path(
        project_dirs[[omic_type]]$source_dir,
        "study_parameters.txt"
    )
    fallback_created <- FALSE

    if (!values$workflow_args_saved) {
        if (!fileExistsFn(basic_params_file)) {
            catFn("SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")
            basic_content <- paste(
                "Study Parameters",
                "================",
                "",
                paste("Workflow Name:", input$experiment_label),
                paste("Description:", input$description),
                paste("Timestamp:", sysTimeFn()),
                paste("Note: Some parameters could not be saved due to serialization issues"),
                sep = "\n"
            )
            tryCatch({
                writeLinesFn(basic_content, basic_params_file)
                values$workflow_args_saved <- TRUE
                fallback_created <- TRUE
            }, error = function(e) {
                logErrorFn(
                    skipFormatterFn(
                        paste(
                            "Failed to create basic study_parameters.txt:",
                            e$message
                        )
                    )
                )
            })
        }
    }

    catFn("SESSION SUMMARY: Copy to Publication button clicked\n")
    catFn("SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
    catFn("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")

    invisible(withProgressFn(message = "Copying files to publication directory...", {
        tryCatch({
            contrasts_tbl <- NULL
            design_matrix <- NULL

            if (!is.null(workflow_data)) {
                if (!is.null(workflow_data$contrasts_tbl)) {
                    contrasts_tbl <- workflow_data$contrasts_tbl
                    catFn("SESSION SUMMARY: Got contrasts_tbl from workflow_data\n")
                }
                if (!is.null(workflow_data$design_matrix)) {
                    design_matrix <- workflow_data$design_matrix
                    catFn("SESSION SUMMARY: Got design_matrix from workflow_data\n")
                }
            }

            if (is.null(design_matrix)) {
                design_matrix_file <- file.path(
                    project_dirs[[omic_type]]$source_dir,
                    "design_matrix.tab"
                )
                if (fileExistsFn(design_matrix_file)) {
                    design_matrix <- readTsvFn(
                        design_matrix_file,
                        show_col_types = FALSE
                    )
                    catFn("SESSION SUMMARY: Loaded design_matrix from file\n")
                }
            }

            if (is.null(contrasts_tbl)) {
                contrasts_file <- file.path(
                    project_dirs[[omic_type]]$source_dir,
                    "contrasts_tbl.tab"
                )
                if (fileExistsFn(contrasts_file)) {
                    contrasts_tbl <- readTsvFn(
                        contrasts_file,
                        show_col_types = FALSE
                    )
                    catFn("SESSION SUMMARY: Loaded contrasts_tbl from file\n")
                }
            }

            catFn("SESSION SUMMARY: About to call copyToResultsSummary\n")
            catFn("SESSION SUMMARY: omic_type =", omic_type, "\n")
            catFn("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
            catFn(
                "SESSION SUMMARY: contrasts_tbl is",
                ifelse(is.null(contrasts_tbl), "NULL", "available"),
                "\n"
            )
            catFn(
                "SESSION SUMMARY: design_matrix is",
                ifelse(is.null(design_matrix), "NULL", "available"),
                "\n"
            )

            project_dirs_assigned <- FALSE
            if (!existsFn("project_dirs", envir = globalEnv)) {
                catFn("SESSION SUMMARY: Setting project_dirs in global environment\n")
                assignFn("project_dirs", project_dirs, envir = globalEnv)
                project_dirs_assigned <- TRUE
            }

            copyToResultsSummaryFn(
                omic_type = omic_type,
                experiment_label = input$experiment_label,
                contrasts_tbl = contrasts_tbl,
                design_matrix = design_matrix,
                force = TRUE
            )

            values$files_copied <- TRUE

            output$copy_status <- renderTextFn(
                "Files copied to publication directory successfully [OK]"
            )
            showNotificationFn("Publication files copied", type = "message")

            output$session_summary <- renderTextFn({
                paste(
                    "Workflow args created for:", input$experiment_label,
                    "\nDescription:", input$description,
                    "\nTimestamp:", sysTimeFn(),
                    "\nStatus: Arguments saved [OK], Files copied [OK]"
                )
            })

            list(
                status = "success",
                basicParamsFile = basic_params_file,
                fallbackCreated = fallback_created,
                projectDirsAssigned = project_dirs_assigned,
                contrastsTbl = contrasts_tbl,
                designMatrix = design_matrix
            )
        }, error = function(e) {
            error_message <- if (inherits(e, "condition")) {
                conditionMessage(e)
            } else {
                as.character(e)
            }

            output$copy_status <- renderTextFn(paste("Error:", error_message))
            showNotificationFn(
                paste("Copy error:", error_message),
                type = "error",
                duration = 10
            )
            logErrorFn(
                skipFormatterFn(
                    paste("Failed to copy files:", error_message)
                )
            )
            catFn("SESSION SUMMARY ERROR:", error_message, "\n")
            catFn("SESSION SUMMARY Traceback:\n")
            tracebackFn()

            list(
                status = "error",
                basicParamsFile = basic_params_file,
                fallbackCreated = fallback_created,
                errorMessage = error_message
            )
        })
    }))
}

runMetabSummaryGenerateReportObserverShell <- function(
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
    existsFn = exists,
    getFn = get,
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
    input <- inputValues
    project_dirs <- projectDirs
    omic_type <- omicType

    reqFn(input$experiment_label)
    reqFn(values$files_copied)

    if (!omic_type %in% names(project_dirs) || is.null(project_dirs[[omic_type]]$base_dir)) {
        showNotificationFn(
            "Error: Project directories not properly initialized",
            type = "error",
            duration = 10
        )
        return(invisible(list(
            status = "invalid_project_dirs",
            errorMessage = "Project directories not properly initialized"
        )))
    }

    invisible(withProgressFn(message = "Generating report...", {
        template_filename <- "metabolomics_report.rmd"

        catFn(sprintf("REPORT: Selected template: %s\n", template_filename))

        report_template_path <- file.path(
            project_dirs[[omic_type]]$base_dir,
            "scripts",
            omic_type,
            template_filename
        )
        templateSource <- NULL

        if (!fileExistsFn(report_template_path)) {
            incProgressFn(
                0.2,
                detail = paste("Checking for", template_filename, "template...")
            )

            dirCreateFn(
                dirname(report_template_path),
                recursive = TRUE,
                showWarnings = FALSE
            )

            templateSource <- tryCatch({
                pkg_file <- systemFileFn(
                    "reports",
                    "metabolomics",
                    template_filename,
                    package = "MultiScholaR"
                )

                if (fileExistsFn(pkg_file) && pkg_file != "") {
                    catFn(sprintf("REPORT: Found template in package at: %s\n", pkg_file))
                    fileCopyFn(pkg_file, report_template_path)
                    logInfoFn("Template copied from package to: {report_template_path}")
                    showNotificationFn(
                        paste(template_filename, "template copied from package"),
                        type = "message"
                    )
                    "package"
                } else {
                    catFn("REPORT: Template not found in package, downloading from GitHub...\n")

                    template_url <- paste0(
                        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/metabolomics/",
                        template_filename
                    )

                    catFn(sprintf("REPORT: Downloading template from: %s\n", template_url))
                    downloadFileFn(
                        template_url,
                        destfile = report_template_path,
                        quiet = TRUE
                    )
                    logInfoFn("Template downloaded to: {report_template_path}")
                    showNotificationFn(
                        paste(template_filename, "template downloaded"),
                        type = "message"
                    )
                    "github"
                }
            }, error = function(e) {
                error_message <- if (inherits(e, "condition")) {
                    conditionMessage(e)
                } else {
                    as.character(e)
                }

                showNotificationFn(
                    paste("Template retrieval failed:", error_message),
                    type = "error",
                    duration = 10
                )
                logErrorFn(
                    skipFormatterFn(
                        paste("Failed to retrieve template:", error_message)
                    )
                )

                NULL
            })
        } else {
            catFn(sprintf("REPORT: Template already exists at: %s\n", report_template_path))
            templateSource <- "existing"
        }

        incProgressFn(0.5, detail = "Rendering report...")

        tryCatch({
            renderReportImpl <- renderReportFn

            if (is.null(renderReportImpl)) {
                if (!existsFn("RenderReport", inherits = TRUE)) {
                    showNotificationFn(
                        "Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.",
                        type = "error",
                        duration = 15
                    )
                    return(invisible(list(
                        status = "missing_render_report",
                        templateFilename = template_filename,
                        reportTemplatePath = report_template_path,
                        templateSource = templateSource
                    )))
                }

                renderReportImpl <- getFn("RenderReport", inherits = TRUE)
            }

            logInfoFn(
                "Calling RenderReport with omic_type: {omic_type}, experiment_label: {input$experiment_label}"
            )

            rendered_path <- renderReportImpl(
                omic_type = omic_type,
                experiment_label = input$experiment_label,
                rmd_filename = template_filename
            )

            logInfoFn("RenderReport returned path: {rendered_path}")

            if (!is.null(rendered_path) && fileExistsFn(rendered_path)) {
                values$report_generated <- TRUE
                values$report_path <- rendered_path

                output$report_ready <- reactiveFn({ TRUE })
                outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)

                output$download_report <- downloadHandlerFn(
                    filename = function() basename(rendered_path),
                    content = function(file) fileCopyFn(rendered_path, file)
                )

                showNotificationFn("Report generated successfully!", type = "message")

                output$session_summary <- renderTextFn({
                    paste("Workflow args created for:", input$experiment_label,
                          "\nDescription:", input$description,
                          "\nTimestamp:", sysTimeFn(),
                          "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
                          "\nReport location:", rendered_path)
                })

                list(
                    status = "success",
                    templateFilename = template_filename,
                    reportTemplatePath = report_template_path,
                    templateSource = templateSource,
                    renderedPath = rendered_path
                )
            } else {
                showNotificationFn(
                    "Report generation failed - no output file created",
                    type = "error",
                    duration = 10
                )

                list(
                    status = "missing_output",
                    templateFilename = template_filename,
                    reportTemplatePath = report_template_path,
                    templateSource = templateSource,
                    renderedPath = rendered_path
                )
            }
        }, error = function(e) {
            error_message <- if (inherits(e, "condition")) {
                conditionMessage(e)
            } else {
                as.character(e)
            }

            error_msg <- paste("Report generation failed:", error_message)

            catFn("DEBUG66: REPORT GENERATION ERROR\n")
            catFn(sprintf("Error class: %s\n", class(e)[1]))
            catFn(sprintf("Error message: %s\n", error_message))
            catFn("Full error object:\n")
            printFn(e)

            logErrorFn("Failed to generate report: {e$message}")
            logErrorFn("Error class: {class(e)[1]}")

            showNotificationFn(error_msg, type = "error", duration = 15)

            showNotificationFn(
                "Debug info: Check R console for detailed error trace",
                type = "warning",
                duration = 10
            )

            list(
                status = "error",
                templateFilename = template_filename,
                reportTemplatePath = report_template_path,
                templateSource = templateSource,
                errorMessage = error_message
            )
        })
    }))
}

runMetabSummaryGithubPushObserverShell <- function(
    inputValues,
    projectDirs,
    omicType,
    values,
    output,
    reqFn = shiny::req,
    withProgressFn = shiny::withProgress,
    optionsFn = options,
    pushProjectToGithubFromDirsFn = pushProjectToGithubFromDirs,
    showNotificationFn = shiny::showNotification,
    renderTextFn = shiny::renderText,
    sysTimeFn = Sys.time,
    logErrorFn = logger::log_error,
    skipFormatterFn = logger::skip_formatter
) {
    reqFn(
        inputValues$enable_github,
        inputValues$github_org,
        inputValues$github_email,
        inputValues$github_username,
        inputValues$project_id
    )
    reqFn(values$report_generated)

    invisible(withProgressFn(message = "Pushing to GitHub...", {
        tryCatch({
            optionsFn(
                github_org = inputValues$github_org,
                github_user_email = inputValues$github_email,
                github_user_name = inputValues$github_username
            )

            pushProjectToGithubFromDirsFn(
                project_dirs = projectDirs,
                omic_type = omicType,
                experiment_label = inputValues$experiment_label,
                project_id = inputValues$project_id
            )

            showNotificationFn("Successfully pushed to GitHub", type = "message")

            output$session_summary <- renderTextFn({
                paste("Workflow args created for:", inputValues$experiment_label,
                      "\nDescription:", inputValues$description,
                      "\nTimestamp:", sysTimeFn(),
                      "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]")
            })

            list(status = "success")
        }, error = function(e) {
            error_message <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

            showNotificationFn(
                paste("GitHub push failed:", error_message),
                type = "error",
                duration = 10
            )
            logErrorFn(
                skipFormatterFn(
                    paste("Failed to push to GitHub:", error_message)
                )
            )

            list(
                status = "error",
                errorMessage = error_message
            )
        })
    }))
}

