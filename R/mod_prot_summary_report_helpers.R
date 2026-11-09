
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
                                         incProgressFn = shiny::incProgress,
                                         exportRecordFn = exportWorkflowAnalysisRecordBundle,
                                         showNotificationFn = shiny::showNotification) {
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

  generated <- runReportGenerationFn(
    output = output,
    values = values,
    omicType = omicType,
    experimentLabel = experimentLabel,
    description = description,
    templateFilename = templateFilename
  )

  # KIT-029: the analysis record bundle (JSON, Word, Excel) beside the report, for an on-disk owner; a failure
  # is shown and leaves the report as generated.
  if (isTRUE(generated) && !is.null(workflowData) && inherits(workflowData$state_manager, "DIAStateAdapter")) {
    incProgressFn(0.9, detail = "Exporting the analysis record...")
    recordDir <- projectDirs[[omicType]]$study_report_dir %||% projectDirs[[omicType]]$results_summary_dir
    tryCatch(exportRecordFn(workflowData$state_manager, recordDir), error = function(e) showNotificationFn(
      paste("The analysis record bundle was not exported:", conditionMessage(e)), type = "warning", duration = NULL))
  }

  generated
}

