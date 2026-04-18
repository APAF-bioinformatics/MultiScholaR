setupMetabSummaryTemplateStatusOutput <- function(output, projectDirs, omicType,
                                                  renderText = shiny::renderText,
                                                  req = shiny::req) {
    output$template_status <- renderText({
        req(projectDirs)

        if (!omicType %in% names(projectDirs)) {
            return("[WARNING] Project directories not available")
        }

        template_dir <- file.path(
            projectDirs[[omicType]]$base_dir,
            "scripts",
            omicType
        )

        metab_template <- file.path(template_dir, "metabolomics_report.rmd")

        if (file.exists(metab_template)) {
            "Template: Metabolomics Report [OK]"
        } else {
            "[WARNING] Report template will be downloaded when generating report"
        }
    })

    invisible(NULL)
}

setupMetabSummaryBootstrapOutputs <- function(output,
                                              renderText = shiny::renderText,
                                              reactive = shiny::reactive,
                                              outputOptions = shiny::outputOptions) {
    output$session_summary <- renderText({
        "Ready to save workflow parameters and generate report"
    })

    output$report_ready <- reactive({ FALSE })
    outputOptions(output, "report_ready", suspendWhenHidden = FALSE)

    invisible(NULL)
}

setupMetabSummaryServerBootstrapState <- function(
    session,
    experimentLabel,
    reactiveValuesFn = shiny::reactiveValues,
    updateTextInputFn = shiny::updateTextInput
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

registerMetabSummaryServerObservers <- function(
    input,
    output,
    projectDirs,
    omicType,
    workflowData,
    values,
    observeEventFn = shiny::observeEvent,
    runSaveWorkflowArgsObserverShellFn = runMetabSummarySaveWorkflowArgsObserverShell,
    runCopyToPublicationObserverShellFn = runMetabSummaryCopyToPublicationObserverShell,
    runGenerateReportObserverShellFn = runMetabSummaryGenerateReportObserverShell,
    runGithubPushObserverShellFn = runMetabSummaryGithubPushObserverShell,
    runExportSessionObserverShellFn = runMetabSummaryExportSessionObserverShell
) {
    registrations <- list(
        saveWorkflowArgsObserver = observeEventFn(input$save_workflow_args, {
            runSaveWorkflowArgsObserverShellFn(
                inputValues = list(
                    experiment_label = input$experiment_label,
                    description = input$description
                ),
                projectDirs = projectDirs,
                omicType = omicType,
                workflowData = workflowData,
                values = values,
                output = output
            )
        }),
        copyToPublicationObserver = observeEventFn(input$copy_to_publication, {
            runCopyToPublicationObserverShellFn(
                inputValues = list(
                    experiment_label = input$experiment_label,
                    description = input$description
                ),
                projectDirs = projectDirs,
                omicType = omicType,
                workflowData = workflowData,
                values = values,
                output = output
            )
        }, ignoreInit = TRUE),
        generateReportObserver = observeEventFn(input$generate_report, {
            runGenerateReportObserverShellFn(
                inputValues = list(
                    experiment_label = input$experiment_label,
                    description = input$description
                ),
                projectDirs = projectDirs,
                omicType = omicType,
                values = values,
                output = output
            )
        }),
        githubPushObserver = observeEventFn(input$push_to_github, {
            runGithubPushObserverShellFn(
                inputValues = list(
                    enable_github = input$enable_github,
                    github_org = input$github_org,
                    github_email = input$github_email,
                    github_username = input$github_username,
                    project_id = input$project_id,
                    experiment_label = input$experiment_label,
                    description = input$description
                ),
                projectDirs = projectDirs,
                omicType = omicType,
                values = values,
                output = output
            )
        }),
        exportSessionObserver = observeEventFn(input$export_session_state, {
            runExportSessionObserverShellFn(
                inputValues = list(
                    experiment_label = input$experiment_label,
                    description = input$description
                ),
                projectDirs = projectDirs,
                omicType = omicType,
                values = values
            )
        })
    )

    invisible(registrations)
}

