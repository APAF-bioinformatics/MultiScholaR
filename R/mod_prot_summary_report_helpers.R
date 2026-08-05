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

  templateMap <- c(
    "DIA-NN" = "DIANN_report.rmd",
    "DIA-NN limpa" = "DIANN_limpa_report.rmd",
    "TMT" = "TMT_report.rmd",
    "LFQ" = "LFQ_report.rmd"
  )
  templatesStatus <- vapply(names(templateMap), function(label) {
    templatePath <- file.path(templateDir, templateMap[[label]])
    if (file.exists(templatePath)) {
      paste(label, "[OK]")
    } else {
      paste(label, "[missing]")
    }
  }, character(1))

  if (any(grepl("\\[OK\\]", templatesStatus))) {
    paste("Templates:", paste(templatesStatus, collapse = ", "))
  } else {
    "[WARNING] Report templates will be downloaded when generating report"
  }
}

resolveProtSummaryExpectedTemplate <- function(workflowType,
                                               limpaRequested = FALSE) {
  if (is.null(workflowType) || !nzchar(as.character(workflowType))) {
    workflowType <- "DIA"
  }
  workflowType <- tolower(as.character(workflowType))

  if (isTRUE(limpaRequested)) {
    return("DIANN_limpa_report.rmd")
  }

  if (workflowType %in% c("tmt", "tmt_pd")) {
    return("TMT_report.rmd")
  }

  if (workflowType == "lfq") {
    return("LFQ_report.rmd")
  }

  "DIANN_report.rmd"
}

resolveProtSummaryTemplateSelection <- function(workflowType,
                                                requestedTemplate = NULL,
                                                limpaRequested = FALSE) {
  expectedTemplate <- resolveProtSummaryExpectedTemplate(
    workflowType = workflowType,
    limpaRequested = limpaRequested
  )
  requestedTemplate <- if (!is.null(requestedTemplate) && nzchar(requestedTemplate)) {
    basename(requestedTemplate)
  } else {
    NULL
  }

  if (is.null(workflowType) || !nzchar(as.character(workflowType))) {
    workflowType <- "DIA"
  }
  compatibleTemplates <- switch(tolower(as.character(workflowType)),
    tmt = "TMT_report.rmd",
    tmt_pd = "TMT_report.rmd",
    lfq = "LFQ_report.rmd",
    dia = if (isTRUE(limpaRequested)) "DIANN_limpa_report.rmd" else "DIANN_report.rmd",
    expectedTemplate
  )

  if (!is.null(requestedTemplate) && requestedTemplate %in% compatibleTemplates) {
    return(list(
      templateFilename = requestedTemplate,
      expectedTemplate = expectedTemplate,
      requestedTemplate = requestedTemplate,
      requestHonoured = TRUE,
      staleTemplateIgnored = FALSE
    ))
  }

  list(
    templateFilename = expectedTemplate,
    expectedTemplate = expectedTemplate,
    requestedTemplate = requestedTemplate,
    requestHonoured = FALSE,
    staleTemplateIgnored = !is.null(requestedTemplate)
  )
}

resolveProtSummaryReportTemplate <- function(workflowData,
                                             logPrefix = "REPORT",
                                             catFn = cat) {
  workflowTypeDetected <- NULL
  dataStateUsed <- NULL
  reportTemplateRequested <- NULL
  limpaRequested <- FALSE

  if (!is.null(workflowData) &&
      !is.null(workflowData$config_list) &&
      !is.null(workflowData$config_list$globalParameters)) {
    workflowTypeDetected <- workflowData$config_list$globalParameters$workflow_type
    reportTemplateRequested <- workflowData$config_list$globalParameters$report_template
    limpaRequested <- isTRUE(workflowData$config_list$globalParameters$use_limpa)
    if (!is.null(reportTemplateRequested) && nzchar(reportTemplateRequested)) {
      reportTemplateRequested <- basename(reportTemplateRequested)
      catFn(sprintf(
        "   %s: Detected report_template from workflow_data: %s\n",
        logPrefix,
        reportTemplateRequested
      ))
    }
  }

  if (!is.null(workflowTypeDetected) && nzchar(workflowTypeDetected)) {
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
    reportTemplateFromS4 <- tryCatch(
      currentS4@args$globalParameters$report_template,
      error = function(e) NULL
    )
    useLimpaFromS4 <- tryCatch(
      isTRUE(currentS4@args$globalParameters$use_limpa) ||
        !is.null(currentS4@args$limpa_dpc_quant_results) ||
        !is.null(currentS4@args$proteinMissingValueImputationLimpa),
      error = function(e) FALSE
    )

    if ((is.null(reportTemplateRequested) || !nzchar(reportTemplateRequested)) &&
        !is.null(reportTemplateFromS4) &&
        nzchar(reportTemplateFromS4)) {
      reportTemplateRequested <- basename(reportTemplateFromS4)
      catFn(sprintf(
        "   %s: Detected report_template from S4 object (state: %s): %s\n",
        logPrefix,
        dataStateUsed,
        reportTemplateRequested
      ))
    }
    limpaRequested <- isTRUE(limpaRequested) ||
      isTRUE(useLimpaFromS4)

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
    reportTemplateFromGlobal <- configList$globalParameters$report_template
    useLimpaFromGlobal <- isTRUE(configList$globalParameters$use_limpa)

    if ((is.null(reportTemplateRequested) || !nzchar(reportTemplateRequested)) &&
        !is.null(reportTemplateFromGlobal) &&
        nzchar(reportTemplateFromGlobal)) {
      reportTemplateRequested <- basename(reportTemplateFromGlobal)
      catFn(sprintf(
        "   %s: Detected report_template from global config_list: %s\n",
        logPrefix,
        reportTemplateRequested
      ))
    }
    limpaRequested <- isTRUE(limpaRequested) ||
      isTRUE(useLimpaFromGlobal)

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

  templateSelection <- resolveProtSummaryTemplateSelection(
    workflowType = workflowTypeDetected,
    requestedTemplate = reportTemplateRequested,
    limpaRequested = limpaRequested
  )
  templateFilename <- templateSelection$templateFilename

  if (isTRUE(templateSelection$staleTemplateIgnored)) {
    catFn(sprintf(
      "   %s: Ignoring stale report_template '%s' for workflow_type: %s; using %s\n",
      logPrefix,
      templateSelection$requestedTemplate,
      workflowTypeDetected,
      templateFilename
    ))
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
    dataStateUsed = dataStateUsed,
    requestedTemplate = templateSelection$requestedTemplate,
    staleTemplateIgnored = templateSelection$staleTemplateIgnored
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

