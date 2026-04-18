initializeProtSummaryDefaultOutputs <- function(output,
                                                initialSessionSummary = "Ready to save workflow parameters and generate report",
                                                renderTextFn = shiny::renderText,
                                                reactiveFn = shiny::reactive,
                                                outputOptionsFn = shiny::outputOptions) {
  output$session_summary <- renderTextFn({
    initialSessionSummary
  })

  output$report_ready <- reactiveFn({ FALSE })
  outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)

  invisible(TRUE)
}

observeProtSummaryWorkflowArgsSave <- function(input,
                                               output,
                                               values,
                                               projectDirs,
                                               workflowData = NULL,
                                               omicType = "proteomics",
                                               completeWorkflowArgsSaveFn = completeProtSummaryWorkflowArgsSave,
                                               resolveFinalS4StateFn = resolveProtSummaryFinalS4State,
                                               assignFn = assign,
                                               createWorkflowArgsFn = createWorkflowArgsFromConfig,
                                               saveRDSFn = saveRDS,
                                               renderTextFn = shiny::renderText,
                                               showNotificationFn = shiny::showNotification,
                                               dirExistsFn = dir.exists,
                                               dirCreateFn = dir.create,
                                               fileExistsFn = file.exists,
                                               writeLinesFn = writeLines,
                                               timestampFn = Sys.time,
                                               catFn = cat) {
  shiny::observeEvent(input$save_workflow_args, {
    shiny::req(input$experiment_label)

    catFn("SESSION SUMMARY: Starting workflow args save process\n")

    completeWorkflowArgsSaveFn(
      output = output,
      values = values,
      projectDirs = projectDirs,
      workflowData = workflowData,
      omicType = omicType,
      experimentLabel = input$experiment_label,
      description = input$description,
      resolveFinalS4StateFn = resolveFinalS4StateFn,
      assignFn = assignFn,
      createWorkflowArgsFn = createWorkflowArgsFn,
      saveRDSFn = saveRDSFn,
      renderTextFn = renderTextFn,
      showNotificationFn = showNotificationFn,
      dirExistsFn = dirExistsFn,
      dirCreateFn = dirCreateFn,
      fileExistsFn = fileExistsFn,
      writeLinesFn = writeLinesFn,
      timestampFn = timestampFn,
      catFn = catFn
    )
  })
}

observeProtSummaryPublicationCopy <- function(input,
                                              output,
                                              values,
                                              projectDirs,
                                              workflowData = NULL,
                                              omicType = "proteomics",
                                              fallbackBootstrapFn = bootstrapProtSummaryCopyFallbackStudyParams,
                                              prepareCopyInputsFn = prepareProtSummaryCopyInputs,
                                              runPublicationCopyFn = runProtSummaryPublicationCopy,
                                              handleCopyErrorFn = handleProtSummaryPublicationCopyError,
                                              withProgressFn = shiny::withProgress,
                                              catFn = cat) {
  shiny::observeEvent(input$copy_to_publication, {
    shiny::req(input$experiment_label)

    if (!values$workflow_args_saved) {
      fallbackBootstrapFn(
        values = values,
        projectDirs = projectDirs,
        omicType = omicType,
        experimentLabel = input$experiment_label,
        description = input$description
      )
    }

    catFn("   SESSION SUMMARY: Copy to Publication button clicked\n")
    catFn("   SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
    catFn("   SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")

    withProgressFn(message = "Copying files to publication directory...", {
      tryCatch({
        copyInputs <- prepareCopyInputsFn(
          workflowData = workflowData,
          projectDirs = projectDirs,
          omicType = omicType
        )
        contrastsTbl <- copyInputs$contrastsTbl
        designMatrix <- copyInputs$designMatrix

        catFn("   SESSION SUMMARY: About to call copyToResultsSummary\n")
        catFn("   SESSION SUMMARY: omic_type =", omicType, "\n")
        catFn("   SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
        catFn(
          "   SESSION SUMMARY: contrasts_tbl is",
          ifelse(is.null(contrastsTbl), "NULL", "available"),
          "\n"
        )
        catFn(
          "   SESSION SUMMARY: design_matrix is",
          ifelse(is.null(designMatrix), "NULL", "available"),
          "\n"
        )

        runPublicationCopyFn(
          output = output,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = input$experiment_label,
          description = input$description,
          contrastsTbl = contrastsTbl,
          designMatrix = designMatrix
        )
      }, error = function(e) {
        handleCopyErrorFn(
          output = output,
          error = e
        )
      })
    })
  }, ignoreInit = TRUE)
}

observeProtSummarySessionStateExport <- function(input,
                                                 values,
                                                 projectDirs,
                                                 omicType = "proteomics",
                                                 completeSessionStateExportFn = completeProtSummarySessionStateExport) {
  shiny::observeEvent(input$export_session_state, {
    shiny::req(input$experiment_label)

    completeSessionStateExportFn(
      projectDirs = projectDirs,
      omicType = omicType,
      experimentLabel = input$experiment_label,
      description = input$description,
      workflowArgsSaved = values$workflow_args_saved,
      filesCopied = values$files_copied,
      reportGenerated = values$report_generated,
      reportPath = values$report_path
    )
  })
}

observeProtSummaryReportGeneration <- function(input,
                                               output,
                                               values,
                                               projectDirs,
                                               workflowData = NULL,
                                               omicType = "proteomics",
                                               validateProjectDirsFn = validateProtSummaryProjectDirs,
                                               runReportProgressFn = runProtSummaryReportProgress,
                                               withProgressFn = shiny::withProgress) {
  shiny::observeEvent(input$generate_report, {
    shiny::req(input$experiment_label)
    shiny::req(values$files_copied)

    if (!validateProjectDirsFn(
      projectDirs = projectDirs,
      omicType = omicType
    )) {
      return()
    }

    withProgressFn(message = "Generating report...", {
      runReportProgressFn(
        output = output,
        values = values,
        workflowData = workflowData,
        projectDirs = projectDirs,
        omicType = omicType,
        experimentLabel = input$experiment_label,
        description = input$description
      )
    })
  })
}

observeProtSummaryGithubPush <- function(input,
                                         output,
                                         values,
                                         projectDirs,
                                         omicType = "proteomics",
                                         completeGithubPushFn = completeProtSummaryGithubPush,
                                         withProgressFn = shiny::withProgress) {
  shiny::observeEvent(input$push_to_github, {
    shiny::req(
      input$enable_github,
      input$github_org,
      input$github_email,
      input$github_username,
      input$project_id
    )
    shiny::req(values$report_generated)

    withProgressFn(message = "Pushing to GitHub...", {
      completeGithubPushFn(
        output = output,
        projectDirs = projectDirs,
        omicType = omicType,
        experimentLabel = input$experiment_label,
        description = input$description,
        githubOrg = input$github_org,
        githubEmail = input$github_email,
        githubUsername = input$github_username,
        projectId = input$project_id
      )
    })
  })
}

registerProtSummaryTemplateStatusOutput <- function(output,
                                                    projectDirs,
                                                    omicType = "proteomics",
                                                    renderTextFn = shiny::renderText,
                                                    reqFn = shiny::req,
                                                    buildTemplateStatusFn = buildProtSummaryTemplateStatus) {
  output$template_status <- renderTextFn({
    reqFn(projectDirs)

    buildTemplateStatusFn(
      projectDirs = projectDirs,
      omicType = omicType
    )
  })

  invisible(TRUE)
}

