resolveProtSummaryFinalS4State <- function(workflowData,
                                           logPrefix = "SESSION SUMMARY",
                                           catFn = cat) {
  finalS4Object <- NULL
  dataStateUsed <- NULL
  stateManager <- workflowData$state_manager

  if (is.null(stateManager)) {
    catFn(sprintf("%s: No state manager available\n", logPrefix))
    return(list(finalS4Object = NULL, dataStateUsed = NULL))
  }

  dataStates <- c(
    "correlation_filtered",
    "ruv_corrected",
    "protein_replicate_filtered",
    "imputed"
  )
  availableStates <- stateManager$getHistory()
  dataStateUsed <- purrr::detect(dataStates, ~ .x %in% availableStates)

  if (!is.null(dataStateUsed)) {
    finalS4Object <- stateManager$getState(dataStateUsed)
  }

  if (is.null(finalS4Object)) {
    catFn(sprintf("%s: No data S4 object found in any valid state\n", logPrefix))
    return(list(finalS4Object = NULL, dataStateUsed = dataStateUsed))
  }

  catFn(sprintf(
    "%s: Retrieved DATA S4 object from state '%s'\n",
    logPrefix,
    dataStateUsed
  ))
  catFn(sprintf("%s: S4 object class: %s\n", logPrefix, class(finalS4Object)))

  argsAvailable <- tryCatch({
    !is.null(finalS4Object@args)
  }, error = function(e) {
    FALSE
  })

  if (argsAvailable) {
    catFn(sprintf(
      "%s: S4 @args contains %d function groups\n",
      logPrefix,
      length(finalS4Object@args)
    ))
  } else {
    catFn(sprintf("%s: S4 @args is NULL or slot doesn't exist\n", logPrefix))
  }

  list(finalS4Object = finalS4Object, dataStateUsed = dataStateUsed)
}

buildProtSummaryTemplateStatus <- function(projectDirs,
                                           omicType = "proteomics") {
  if (!omicType %in% names(projectDirs)) {
    return("[WARNING] Project directories not available")
  }

  templateDir <- file.path(
    projectDirs[[omicType]]$base_dir,
    "scripts",
    omicType
  )

  diannTemplate <- file.path(templateDir, "DIANN_report.rmd")
  tmtTemplate <- file.path(templateDir, "TMT_report.rmd")

  templatesStatus <- character()
  if (file.exists(diannTemplate)) {
    templatesStatus <- c(templatesStatus, "DIA-NN [OK]")
  }
  if (file.exists(tmtTemplate)) {
    templatesStatus <- c(templatesStatus, "TMT [OK]")
  }

  if (length(templatesStatus) > 0) {
    paste("Templates:", paste(templatesStatus, collapse = ", "))
  } else {
    "[WARNING] Report templates will be downloaded when generating report"
  }
}

resolveProtSummaryReportTemplate <- function(workflowData,
                                             logPrefix = "REPORT",
                                             catFn = cat) {
  workflowTypeDetected <- NULL
  dataStateUsed <- NULL

  if (!is.null(workflowData) &&
      !is.null(workflowData$config_list) &&
      !is.null(workflowData$config_list$globalParameters$workflow_type)) {
    workflowTypeDetected <- workflowData$config_list$globalParameters$workflow_type
    catFn(sprintf(
      "   %s: Detected workflow_type from workflow_data: %s\n",
      logPrefix,
      workflowTypeDetected
    ))
  }

  if ((is.null(workflowTypeDetected) || !nzchar(workflowTypeDetected)) &&
      !is.null(workflowData) &&
      !is.null(workflowData$state_manager)) {
    finalS4State <- resolveProtSummaryFinalS4State(
      workflowData,
      logPrefix = logPrefix,
      catFn = catFn
    )
    currentS4 <- finalS4State$finalS4Object
    dataStateUsed <- finalS4State$dataStateUsed
    workflowTypeFromS4 <- tryCatch(
      currentS4@args$globalParameters$workflow_type,
      error = function(e) NULL
    )

    if (!is.null(workflowTypeFromS4) && nzchar(workflowTypeFromS4)) {
      workflowTypeDetected <- workflowTypeFromS4
      catFn(sprintf(
        "   %s: Detected workflow_type from S4 object (state: %s): %s\n",
        logPrefix,
        dataStateUsed,
        workflowTypeDetected
      ))
    }
  }

  if ((is.null(workflowTypeDetected) || !nzchar(workflowTypeDetected)) &&
      exists("config_list", envir = .GlobalEnv)) {
    configList <- get("config_list", envir = .GlobalEnv)
    workflowTypeFromGlobal <- configList$globalParameters$workflow_type

    if (!is.null(workflowTypeFromGlobal) && nzchar(workflowTypeFromGlobal)) {
      workflowTypeDetected <- workflowTypeFromGlobal
      catFn(sprintf(
        "   %s: Detected workflow_type from global config_list: %s\n",
        logPrefix,
        workflowTypeDetected
      ))
    }
  }

  if (is.null(workflowTypeDetected) || !nzchar(workflowTypeDetected)) {
    workflowTypeDetected <- "DIA"
    catFn(sprintf(
      "   %s: Using default workflow_type: DIA (no workflow type found)\n",
      logPrefix
    ))
  }

  templateFilename <- if (tolower(workflowTypeDetected) %in% c("tmt", "tmt_pd")) {
    "TMT_report.rmd"
  } else if (tolower(workflowTypeDetected) == "lfq") {
    "LFQ_report.rmd"
  } else {
    "DIANN_report.rmd"
  }

  catFn(sprintf(
    "   %s: Selected template: %s for workflow_type: %s\n",
    logPrefix,
    templateFilename,
    workflowTypeDetected
  ))

  list(
    workflowTypeDetected = workflowTypeDetected,
    templateFilename = templateFilename,
    dataStateUsed = dataStateUsed
  )
}

ensureProtSummaryReportTemplate <- function(projectDirs,
                                            templateFilename,
                                            omicType = "proteomics",
                                            packageName = "MultiScholaR",
                                            packageReportSubdir = "proteomics",
                                            systemFileFn = system.file,
                                            fileExistsFn = file.exists,
                                            dirCreateFn = dir.create,
                                            fileCopyFn = file.copy,
                                            downloadFileFn = download.file,
                                            showNotificationFn = shiny::showNotification,
                                            logInfoFn = function(message) logger::log_info(message),
                                            catFn = cat) {
  reportTemplatePath <- file.path(
    projectDirs[[omicType]]$base_dir,
    "scripts",
    omicType,
    templateFilename
  )

  if (fileExistsFn(reportTemplatePath)) {
    catFn(sprintf("   REPORT: Template already exists at: %s\n", reportTemplatePath))
    return(list(
      reportTemplatePath = reportTemplatePath,
      templateSource = "existing",
      templateUrl = NULL
    ))
  }

  dirCreateFn(dirname(reportTemplatePath), recursive = TRUE, showWarnings = FALSE)

  pkgFile <- systemFileFn(
    "reports",
    packageReportSubdir,
    templateFilename,
    package = packageName
  )

  if (nzchar(pkgFile) && fileExistsFn(pkgFile)) {
    catFn(sprintf("   REPORT: Found template in package at: %s\n", pkgFile))
    fileCopyFn(pkgFile, reportTemplatePath)
    logInfoFn(sprintf("Template copied from package to: %s", reportTemplatePath))
    showNotificationFn(
      paste(templateFilename, "template copied from package"),
      type = "message"
    )

    return(list(
      reportTemplatePath = reportTemplatePath,
      templateSource = "package",
      templateUrl = NULL
    ))
  }

  templateUrl <- paste0(
    "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/",
    packageReportSubdir,
    "/",
    templateFilename
  )

  catFn("   REPORT: Template not found in package, downloading from GitHub...\n")
  catFn(sprintf("   REPORT: Downloading template from: %s\n", templateUrl))
  downloadFileFn(templateUrl, destfile = reportTemplatePath, quiet = TRUE)
  logInfoFn(sprintf("Template downloaded to: %s", reportTemplatePath))
  showNotificationFn(
    paste(templateFilename, "template downloaded"),
    type = "message"
  )

  list(
    reportTemplatePath = reportTemplatePath,
    templateSource = "github",
    templateUrl = templateUrl
  )
}

retrieveProtSummaryReportTemplateAsset <- function(projectDirs,
                                                   templateFilename,
                                                   omicType = "proteomics",
                                                   packageName = "MultiScholaR",
                                                   packageReportSubdir = "proteomics",
                                                   ensureTemplateFn = ensureProtSummaryReportTemplate,
                                                   showNotificationFn = shiny::showNotification,
                                                   logErrorFn = function(message) {
                                                     logger::log_error(logger::skip_formatter(message))
                                                   }) {
  tryCatch({
    ensureTemplateFn(
      projectDirs = projectDirs,
      templateFilename = templateFilename,
      omicType = omicType,
      packageName = packageName,
      packageReportSubdir = packageReportSubdir
    )
  }, error = function(e) {
    showNotificationFn(
      paste("Template retrieval failed:", e$message),
      type = "error",
      duration = 10
    )
    logErrorFn(paste("Failed to retrieve template:", e$message))
    NULL
  })
}

validateProtSummaryProjectDirs <- function(projectDirs,
                                           omicType = "proteomics",
                                           showNotificationFn = shiny::showNotification) {
  if (!omicType %in% names(projectDirs) || is.null(projectDirs[[omicType]]$base_dir)) {
    showNotificationFn(
      "Error: Project directories not properly initialized",
      type = "error",
      duration = 10
    )
    return(FALSE)
  }

  TRUE
}

activateProtSummaryRenderedReport <- function(output,
                                              values,
                                              renderedPath,
                                              experimentLabel,
                                              description,
                                              fileExistsFn = file.exists,
                                              reactiveFn = shiny::reactive,
                                              outputOptionsFn = shiny::outputOptions,
                                              downloadHandlerFn = shiny::downloadHandler,
                                              basenameFn = basename,
                                              fileCopyFn = file.copy,
                                              renderTextFn = shiny::renderText,
                                              showNotificationFn = shiny::showNotification,
                                              timestampFn = Sys.time) {
  if (is.null(renderedPath) || !fileExistsFn(renderedPath)) {
    return(FALSE)
  }

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
    paste("Workflow args created for:", experimentLabel,
          "\nDescription:", description,
          "\nTimestamp:", timestampFn(),
          "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
          "\nReport location:", renderedPath)
  })

  TRUE
}

runProtSummaryReportGeneration <- function(output,
                                           values,
                                           omicType = "proteomics",
                                           experimentLabel,
                                           description,
                                           templateFilename,
                                           renderReportAvailableFn = function() {
                                             exists("RenderReport", mode = "function", inherits = TRUE)
                                           },
                                           renderReportFn = NULL,
                                           activateReportFn = activateProtSummaryRenderedReport,
                                           showNotificationFn = shiny::showNotification,
                                           logInfoFn = function(message) logger::log_info(message),
                                           logErrorFn = function(message) logger::log_error(message),
                                           catFn = cat,
                                           printFn = print,
                                           tracebackFn = traceback) {
  if (!renderReportAvailableFn()) {
    showNotificationFn(
      "Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.",
      type = "error",
      duration = 15
    )
    return(FALSE)
  }

  if (is.null(renderReportFn)) {
    renderReportFn <- get("RenderReport", mode = "function", inherits = TRUE)
  }

  tryCatch({
    logInfoFn(sprintf(
      "Calling RenderReport with omic_type: %s, experiment_label: %s",
      omicType,
      experimentLabel
    ))

    renderedPath <- renderReportFn(
      omic_type = omicType,
      experiment_label = experimentLabel,
      rmd_filename = templateFilename
    )

    logInfoFn(sprintf("RenderReport returned path: %s", renderedPath))

    reportActivated <- activateReportFn(
      output = output,
      values = values,
      renderedPath = renderedPath,
      experimentLabel = experimentLabel,
      description = description
    )

    if (!isTRUE(reportActivated)) {
      showNotificationFn(
        "Report generation failed - no output file created",
        type = "error",
        duration = 10
      )
    }

    isTRUE(reportActivated)
  }, error = function(e) {
    errorMsg <- paste("Report generation failed:", e$message)

    catFn("   DABUG66: REPORT GENERATION ERROR\n")
    catFn(sprintf("      Error class: %s\n", class(e)[1]))
    catFn(sprintf("      Error message: %s\n", e$message))
    catFn("      Full error object:\n")
    printFn(e)

    logErrorFn(sprintf("Failed to generate report: %s", e$message))
    logErrorFn(sprintf("Error class: %s", class(e)[1]))

    showNotificationFn(errorMsg, type = "error", duration = 15)
    showNotificationFn(
      "Debug info: Check R console for detailed error trace",
      type = "warning",
      duration = 10
    )

    tracebackFn()
    FALSE
  })
}

runProtSummaryReportProgress <- function(output,
                                         values,
                                         workflowData,
                                         projectDirs,
                                         omicType = "proteomics",
                                         experimentLabel,
                                         description,
                                         resolveReportTemplateFn = resolveProtSummaryReportTemplate,
                                         retrieveTemplateAssetFn = retrieveProtSummaryReportTemplateAsset,
                                         runReportGenerationFn = runProtSummaryReportGeneration,
                                         incProgressFn = shiny::incProgress) {
  incProgressFn(0.1, detail = "Detecting workflow type...")

  reportTemplate <- resolveReportTemplateFn(
    workflowData,
    logPrefix = "REPORT"
  )
  templateFilename <- reportTemplate$templateFilename

  incProgressFn(0.2, detail = paste("Checking for", templateFilename, "template..."))

  templateAsset <- retrieveTemplateAssetFn(
    projectDirs = projectDirs,
    templateFilename = templateFilename,
    omicType = omicType
  )

  if (is.null(templateAsset)) {
    return(FALSE)
  }

  incProgressFn(0.5, detail = "Rendering report...")

  runReportGenerationFn(
    output = output,
    values = values,
    omicType = omicType,
    experimentLabel = experimentLabel,
    description = description,
    templateFilename = templateFilename
  )
}

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

