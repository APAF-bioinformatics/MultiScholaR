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

summariseProtSummaryWorkflowMetadata <- function(workflowData,
                                                resolveReportTemplateFn = resolveProtSummaryReportTemplate) {
  globalParameters <- tryCatch(
    workflowData$config_list$globalParameters,
    error = function(e) NULL
  )
  workflowType <- globalParameters$workflow_type
  reportTemplate <- globalParameters$report_template

  if (!is.null(workflowData)) {
    templateInfo <- tryCatch(
      resolveReportTemplateFn(
        workflowData = workflowData,
        catFn = function(...) invisible(NULL)
      ),
      error = function(e) NULL
    )
    if (!is.null(templateInfo)) {
      workflowType <- templateInfo$workflowTypeDetected
      reportTemplate <- templateInfo$templateFilename
    }
  }

  list(
    workflow_type = workflowType,
    report_template = reportTemplate,
    parameters = list(
      globalParameters = globalParameters,
      da_ui_params = tryCatch(workflowData$da_ui_params, error = function(e) NULL),
      enrichment_ui_params = tryCatch(workflowData$enrichment_ui_params, error = function(e) NULL),
      ruv_optimization_result = tryCatch(workflowData$ruv_optimization_result, error = function(e) NULL)
    )
  )
}

prepareProtSummarySessionStateExport <- function(projectDirs,
                                                 omicType = "proteomics",
                                                 experimentLabel,
                                                 description,
                                                 workflowArgsSaved,
                                                 filesCopied,
                                                 reportGenerated,
                                                 reportPath,
                                                 workflowData = NULL,
                                                 exportDate = Sys.Date(),
                                                 timestamp = Sys.time(),
                                                 summariseWorkflowMetadataFn = summariseProtSummaryWorkflowMetadata) {
  sessionExportPath <- file.path(
    projectDirs[[omicType]]$source_dir,
    paste0("session_state_", exportDate, ".RDS")
  )
  workflowMetadata <- summariseWorkflowMetadataFn(workflowData)

  sessionState <- list(
    experiment_label = experimentLabel,
    description = description,
    timestamp = timestamp,
    omic_type = omicType,
    workflow_args_saved = workflowArgsSaved,
    files_copied = filesCopied,
    report_generated = reportGenerated,
    report_path = reportPath,
    workflow_type = workflowMetadata$workflow_type,
    report_template = workflowMetadata$report_template,
    parameter_payload = workflowMetadata$parameters,
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
                                                  workflowData = NULL,
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
    exportArgs <- list(
      projectDirs = projectDirs,
      omicType = omicType,
      experimentLabel = experimentLabel,
      description = description,
      workflowArgsSaved = workflowArgsSaved,
      filesCopied = filesCopied,
      reportGenerated = reportGenerated,
      reportPath = reportPath
    )
    if ("workflowData" %in% names(formals(prepareExportFn))) {
      exportArgs$workflowData <- workflowData
    }
    sessionStateExport <- do.call(prepareExportFn, exportArgs)

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

