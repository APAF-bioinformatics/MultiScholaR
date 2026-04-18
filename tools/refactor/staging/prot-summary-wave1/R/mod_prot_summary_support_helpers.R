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

