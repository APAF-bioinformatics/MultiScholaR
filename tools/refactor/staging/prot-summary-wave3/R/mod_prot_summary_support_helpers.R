completeProtSummaryWorkflowArgsSave <- function(output,
                                                values,
                                                projectDirs,
                                                workflowData,
                                                omicType = "proteomics",
                                                experimentLabel,
                                                description,
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
  tryCatch({
    finalS4State <- resolveFinalS4StateFn(workflowData)
    finalS4Object <- finalS4State$finalS4Object

    contrastsTbl <- NULL
    if (!is.null(workflowData) && !is.null(workflowData$contrasts_tbl)) {
      contrastsTbl <- workflowData$contrasts_tbl
      catFn("SESSION SUMMARY: Using contrasts_tbl from workflow_data\n")
    }

    if (!is.null(workflowData) && !is.null(workflowData$config_list)) {
      assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
      catFn(
        "SESSION SUMMARY: Config list available with",
        length(workflowData$config_list),
        "items\n"
      )
    }

    catFn(
      "SESSION SUMMARY: Creating study_parameters.txt file with S4 parameters and RUV results\n"
    )
    studyParamsFile <- createWorkflowArgsFn(
      workflow_name = experimentLabel,
      description = description,
      source_dir_path = projectDirs[[omicType]]$source_dir,
      final_s4_object = finalS4Object,
      contrasts_tbl = contrastsTbl,
      workflow_data = workflowData
    )

    catFn(
      "SESSION SUMMARY: Successfully created study_parameters.txt at:",
      studyParamsFile,
      "\n"
    )

    catFn("SESSION SUMMARY: Saving Integration S4 Object...\n")
    integrationDir <- projectDirs[[omicType]]$integration_dir
    if (is.null(integrationDir)) {
      integrationDir <- file.path(projectDirs[[omicType]]$base_dir, "integration")
    }

    if (!dirExistsFn(integrationDir)) {
      dirCreateFn(integrationDir, recursive = TRUE, showWarnings = FALSE)
    }

    s4Filename <- sprintf("%s_%s_final_s4.RDS", omicType, experimentLabel)
    s4Filepath <- file.path(integrationDir, s4Filename)

    saveRDSFn(finalS4Object, s4Filepath)
    catFn(
      sprintf("SESSION SUMMARY: Saved Integration S4 object to: %s\n", s4Filepath)
    )
    showNotificationFn("Saved Integration S4 Object", type = "message")

    values$workflow_args_saved <- TRUE
    showNotificationFn("Study parameters saved successfully", type = "message")

    output$session_summary <- renderTextFn({
      paste(
        "Study parameters created for:",
        experimentLabel,
        "\nDescription:",
        description,
        "\nTimestamp:",
        timestampFn(),
        "\nFile:",
        studyParamsFile,
        "\nSource: Final S4 object @args + config_list",
        "\nIntegration Object:",
        s4Filename,
        "\nStatus: Parameters saved [OK]"
      )
    })

    TRUE
  }, error = function(e) {
    catFn("SESSION SUMMARY ERROR:", e$message, "\n")

    basicParamsFile <- file.path(projectDirs[[omicType]]$source_dir, "study_parameters.txt")
    if (!fileExistsFn(basicParamsFile)) {
      writeLinesFn(c(
        "Study Parameters",
        "================",
        "",
        paste("Workflow Name:", experimentLabel),
        paste("Description:", description),
        paste("Timestamp:", timestampFn()),
        paste("Error:", e$message)
      ), basicParamsFile)
      catFn("SESSION SUMMARY: Created basic fallback file at:", basicParamsFile, "\n")
    }

    values$workflow_args_saved <- TRUE
    showNotificationFn("Study parameters saved with warnings", type = "warning")
    FALSE
  })
}

prepareProtSummarySessionStateExport <- function(projectDirs,
                                                 omicType = "proteomics",
                                                 experimentLabel,
                                                 description,
                                                 workflowArgsSaved,
                                                 filesCopied,
                                                 reportGenerated,
                                                 reportPath,
                                                 exportDate = Sys.Date(),
                                                 timestamp = Sys.time()) {
  sessionExportPath <- file.path(
    projectDirs[[omicType]]$source_dir,
    paste0("session_state_", exportDate, ".RDS")
  )

  sessionState <- list(
    experiment_label = experimentLabel,
    description = description,
    timestamp = timestamp,
    omic_type = omicType,
    workflow_args_saved = workflowArgsSaved,
    files_copied = filesCopied,
    report_generated = reportGenerated,
    report_path = reportPath,
    project_dirs = projectDirs
  )

  list(sessionExportPath = sessionExportPath, sessionState = sessionState)
}

completeProtSummarySessionStateExport <- function(projectDirs,
                                                  omicType = "proteomics",
                                                  experimentLabel,
                                                  description,
                                                  workflowArgsSaved,
                                                  filesCopied,
                                                  reportGenerated,
                                                  reportPath,
                                                  prepareExportFn = prepareProtSummarySessionStateExport,
                                                  saveRDSFn = saveRDS,
                                                  showNotificationFn = shiny::showNotification,
                                                  logInfoFn = function(message) logger::log_info(message),
                                                  logErrorFn = function(message) {
                                                    logger::log_error(
                                                      logger::skip_formatter(message)
                                                    )
                                                  }) {
  tryCatch({
    sessionStateExport <- prepareExportFn(
      projectDirs = projectDirs,
      omicType = omicType,
      experimentLabel = experimentLabel,
      description = description,
      workflowArgsSaved = workflowArgsSaved,
      filesCopied = filesCopied,
      reportGenerated = reportGenerated,
      reportPath = reportPath
    )

    saveRDSFn(
      sessionStateExport$sessionState,
      sessionStateExport$sessionExportPath
    )

    exportMessage <- paste(
      "Session state exported to:",
      sessionStateExport$sessionExportPath
    )
    showNotificationFn(exportMessage, type = "message")
    logInfoFn(exportMessage)

    TRUE
  }, error = function(e) {
    showNotificationFn(
      paste("Export failed:", e$message),
      type = "error",
      duration = 10
    )
    logErrorFn(paste("Failed to export session state:", e$message))
    FALSE
  })
}

bootstrapProtSummaryCopyFallbackStudyParams <- function(values,
                                                        projectDirs,
                                                        omicType = "proteomics",
                                                        experimentLabel,
                                                        description,
                                                        fileExistsFn = file.exists,
                                                        writeLinesFn = writeLines,
                                                        timestampFn = Sys.time,
                                                        logErrorFn = function(message) {
                                                          logger::log_error(
                                                            logger::skip_formatter(message)
                                                          )
                                                        },
                                                        catFn = cat) {
  if (isTRUE(values$workflow_args_saved)) {
    return(FALSE)
  }

  basicParamsFile <- file.path(
    projectDirs[[omicType]]$source_dir,
    "study_parameters.txt"
  )

  if (fileExistsFn(basicParamsFile)) {
    return(FALSE)
  }

  catFn("   SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")

  basicContent <- paste(
    "Study Parameters",
    "================",
    "",
    paste("Workflow Name:", experimentLabel),
    paste("Description:", description),
    paste("Timestamp:", timestampFn()),
    paste("Note: Some parameters could not be saved due to serialization issues"),
    sep = "\n"
  )

  tryCatch({
    writeLinesFn(basicContent, basicParamsFile)
    values$workflow_args_saved <- TRUE
    TRUE
  }, error = function(e) {
    logErrorFn(paste("Failed to create basic study_parameters.txt:", e$message))
    FALSE
  })
}

prepareProtSummaryCopyInputs <- function(workflowData,
                                         projectDirs,
                                         omicType = "proteomics",
                                         readTsvFn = readr::read_tsv,
                                         catFn = cat) {
  contrastsTbl <- NULL
  designMatrix <- NULL

  if (!is.null(workflowData)) {
    if (!is.null(workflowData$contrasts_tbl)) {
      contrastsTbl <- workflowData$contrasts_tbl
      catFn("   SESSION SUMMARY Step: Got contrasts_tbl from workflow_data\n")
    }
    if (!is.null(workflowData$design_matrix)) {
      designMatrix <- workflowData$design_matrix
      catFn("   SESSION SUMMARY Step: Got design_matrix from workflow_data\n")
    }
  }

  if (is.null(designMatrix)) {
    designMatrixFile <- file.path(projectDirs[[omicType]]$source_dir, "design_matrix.tab")
    if (file.exists(designMatrixFile)) {
      designMatrix <- readTsvFn(designMatrixFile, show_col_types = FALSE)
      catFn("   SESSION SUMMARY Step: Loaded design_matrix from file\n")
    }
  }

  if (is.null(contrastsTbl)) {
    contrastsFile <- file.path(projectDirs[[omicType]]$source_dir, "contrasts_tbl.tab")
    if (file.exists(contrastsFile)) {
      contrastsTbl <- readTsvFn(contrastsFile, show_col_types = FALSE)
      catFn("   SESSION SUMMARY Step: Loaded contrasts_tbl from file\n")
    }
  }

  list(contrastsTbl = contrastsTbl, designMatrix = designMatrix)
}

runProtSummaryPublicationCopy <- function(output,
                                          values,
                                          projectDirs,
                                          omicType = "proteomics",
                                          experimentLabel,
                                          description,
                                          contrastsTbl = NULL,
                                          designMatrix = NULL,
                                          existsFn = exists,
                                          assignFn = assign,
                                          copyFn = copyToResultsSummary,
                                          renderTextFn = shiny::renderText,
                                          showNotificationFn = shiny::showNotification,
                                          timestampFn = Sys.time,
                                          catFn = cat) {
  if (!existsFn("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
    catFn("   SESSION SUMMARY: Setting project_dirs in global environment\n")
    assignFn("project_dirs", projectDirs, envir = .GlobalEnv)
  }

  catFn(
    "   SESSION SUMMARY: project_dirs keys:",
    paste(names(projectDirs), collapse = ", "),
    "\n"
  )
  catFn(
    "   SESSION SUMMARY: Using omic_type =",
    omicType,
    "experiment_label =",
    experimentLabel,
    "\n"
  )

  copyArgs <- list(
    omic_type = omicType,
    experiment_label = experimentLabel,
    contrasts_tbl = contrastsTbl,
    design_matrix = designMatrix,
    force = TRUE
  )
  do.call(copyFn, copyArgs)

  values$files_copied <- TRUE

  output$copy_status <- renderTextFn(
    "Files copied to publication directory successfully [OK]"
  )
  showNotificationFn("Publication files copied", type = "message")

  output$session_summary <- renderTextFn({
    paste(
      "Workflow args created for:",
      experimentLabel,
      "\nDescription:",
      description,
      "\nTimestamp:",
      timestampFn(),
      "\nStatus: Arguments saved [OK], Files copied [OK]"
    )
  })

  copyArgs
}

handleProtSummaryPublicationCopyError <- function(output,
                                                  error,
                                                  renderTextFn = shiny::renderText,
                                                  showNotificationFn = shiny::showNotification,
                                                  logErrorFn = function(message) {
                                                    logger::log_error(
                                                      logger::skip_formatter(message)
                                                    )
                                                  },
                                                  catFn = cat,
                                                  tracebackFn = traceback) {
  output$copy_status <- renderTextFn(paste("Error:", error$message))
  showNotificationFn(
    paste("Copy error:", error$message),
    type = "error",
    duration = 10
  )
  logErrorFn(paste("Failed to copy files:", error$message))
  catFn("   SESSION SUMMARY ERROR:", error$message, "\n")
  catFn("   SESSION SUMMARY Traceback:\n")
  tracebackFn()

  FALSE
}

runProtSummaryGithubPush <- function(projectDirs,
                                     omicType = "proteomics",
                                     experimentLabel,
                                     githubOrg,
                                     githubEmail,
                                     githubUsername,
                                     projectId,
                                     optionsFn = options,
                                     pushFn = pushProjectToGithubFromDirs) {
  githubOptions <- list(
    github_org = githubOrg,
    github_user_email = githubEmail,
    github_user_name = githubUsername
  )
  do.call(optionsFn, githubOptions)

  pushArgs <- list(
    project_dirs = projectDirs,
    omic_type = omicType,
    experiment_label = experimentLabel,
    project_id = projectId
  )
  do.call(pushFn, pushArgs)

  list(githubOptions = githubOptions, pushArgs = pushArgs)
}

completeProtSummaryGithubPush <- function(output,
                                          projectDirs,
                                          omicType = "proteomics",
                                          experimentLabel,
                                          description,
                                          githubOrg,
                                          githubEmail,
                                          githubUsername,
                                          projectId,
                                          pushGithubFn = runProtSummaryGithubPush,
                                          renderTextFn = shiny::renderText,
                                          showNotificationFn = shiny::showNotification,
                                          timestampFn = Sys.time,
                                          logErrorFn = function(message) {
                                            logger::log_error(
                                              logger::skip_formatter(message)
                                            )
                                          }) {
  tryCatch({
    pushGithubFn(
      projectDirs = projectDirs,
      omicType = omicType,
      experimentLabel = experimentLabel,
      githubOrg = githubOrg,
      githubEmail = githubEmail,
      githubUsername = githubUsername,
      projectId = projectId
    )

    showNotificationFn("Successfully pushed to GitHub", type = "message")

    output$session_summary <- renderTextFn({
      paste("Workflow args created for:", experimentLabel,
            "\nDescription:", description,
            "\nTimestamp:", timestampFn(),
            "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]")
    })

    TRUE
  }, error = function(e) {
    showNotificationFn(
      paste("GitHub push failed:", e$message),
      type = "error",
      duration = 10
    )
    logErrorFn(paste("Failed to push to GitHub:", e$message))
    FALSE
  })
}

