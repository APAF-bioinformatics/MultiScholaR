# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

skipIfMissingMultiScholaRBindings <- function(...) {
  missing <- setdiff(c(...), ls(envir = asNamespace("MultiScholaR")))
  if (length(missing) > 0) {
    testthat::skip(sprintf("Target-only extracted helper(s) not present: %s", paste(missing, collapse = ", ")))
  }
}

skipIfMissingMultiScholaRBindings(
  "resolveProtSummaryFinalS4State",
  "buildProtSummaryTemplateStatus",
  "observeProtSummaryReportGeneration"
)

makeFunctionWithOverrides <- function(fun, replacements) {
  funOverride <- fun
  environment(funOverride) <- list2env(replacements, parent = environment(fun))
  funOverride
}

withCleanGlobalObjects <- function(objectNames, code) {
  hadExisting <- vapply(
    objectNames,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
  oldValues <- lapply(seq_along(objectNames), function(i) {
    if (hadExisting[[i]]) {
      get(objectNames[[i]], envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
  })
  names(oldValues) <- objectNames

  on.exit({
    for (name in rev(objectNames)) {
      if (hadExisting[[name]]) {
        assign(name, oldValues[[name]], envir = .GlobalEnv)
      } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = name, envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  force(code)
}

if (!methods::isClass("mockProtSummaryArgsCarrier")) {
  methods::setClass("mockProtSummaryArgsCarrier", slots = c(args = "list"))
}

makeProtSummaryStateManager <- function(history, objectsByState) {
  requested <- character()

  list(
    manager = list(
      getHistory = function() history,
      getState = function(stateName) {
        requested <<- c(requested, stateName)
        objectsByState[[stateName]]
      }
    ),
    requestedStates = function() requested
  )
}

captureProtSummaryMessages <- function() {
  messages <- character()

  list(
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    },
    getMessages = function() messages
  )
}

test_that("resolveProtSummaryFinalS4State prefers the highest-priority saved data state", {
  stateObject <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(), normalization = list())
  )
  stateFixture <- makeProtSummaryStateManager(
    history = c("imputed", "correlation_filtered"),
    objectsByState = list(
      correlation_filtered = stateObject,
      imputed = methods::new("mockProtSummaryArgsCarrier", args = list(legacy = list()))
    )
  )
  messages <- captureProtSummaryMessages()

  resolved <- resolveProtSummaryFinalS4State(
    workflowData = list(state_manager = stateFixture$manager),
    catFn = messages$catFn
  )

  expect_identical(resolved$finalS4Object, stateObject)
  expect_identical(resolved$dataStateUsed, "correlation_filtered")
  expect_identical(stateFixture$requestedStates(), "correlation_filtered")
  expect_true(any(grepl(
    "Retrieved DATA S4 object from state 'correlation_filtered'",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "S4 @args contains 2 function groups",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("resolveProtSummaryFinalS4State reports when no valid state manager data is available", {
  messages <- captureProtSummaryMessages()

  missingManager <- resolveProtSummaryFinalS4State(
    workflowData = list(state_manager = NULL),
    catFn = messages$catFn
  )
  missingState <- resolveProtSummaryFinalS4State(
    workflowData = list(
      state_manager = list(
        getHistory = function() c("analysis_complete"),
        getState = function(stateName) stop(sprintf("unexpected state request: %s", stateName))
      )
    ),
    catFn = messages$catFn
  )

  expect_null(missingManager$finalS4Object)
  expect_null(missingManager$dataStateUsed)
  expect_null(missingState$finalS4Object)
  expect_null(missingState$dataStateUsed)
  expect_true(any(grepl(
    "No state manager available",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "No data S4 object found in any valid state",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("buildProtSummaryTemplateStatus reports available report templates", {
  fixtureDir <- tempfile("prot-summary-templates-")
  templateDir <- file.path(fixtureDir, "scripts", "proteomics")
  dir.create(templateDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines("---", file.path(templateDir, "DIANN_report.rmd"))
  writeLines("---", file.path(templateDir, "TMT_report.rmd"))

  status <- buildProtSummaryTemplateStatus(
    projectDirs = list(proteomics = list(base_dir = fixtureDir)),
    omicType = "proteomics"
  )

  expect_identical(status, "Templates: DIA-NN [OK], TMT [OK]")
})

test_that("buildProtSummaryTemplateStatus warns when project directories are unavailable", {
  expect_identical(
    buildProtSummaryTemplateStatus(projectDirs = list(), omicType = "proteomics"),
    "[WARNING] Project directories not available"
  )
})

test_that("resolveProtSummaryReportTemplate prefers workflow_data over fallback sources", {
  messages <- captureProtSummaryMessages()
  workflowData <- list(
    config_list = list(globalParameters = list(workflow_type = "TMT_PD")),
    state_manager = list(
      getHistory = function() stop("state manager should not be consulted"),
      getState = function(stateName) stop(sprintf("unexpected state request: %s", stateName))
    )
  )

  withCleanGlobalObjects("config_list", {
    assign(
      "config_list",
      list(globalParameters = list(workflow_type = "LFQ")),
      envir = .GlobalEnv
    )

    resolved <- resolveProtSummaryReportTemplate(
      workflowData = workflowData,
      catFn = messages$catFn
    )

    expect_identical(resolved$workflowTypeDetected, "TMT_PD")
    expect_identical(resolved$templateFilename, "TMT_report.rmd")
    expect_null(resolved$dataStateUsed)
  })

  expect_true(any(grepl(
    "Detected workflow_type from workflow_data: TMT_PD",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Selected template: TMT_report.rmd for workflow_type: TMT_PD",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("resolveProtSummaryReportTemplate falls back to final S4 workflow type mapping", {
  stateObject <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(workflow_type = "LFQ"))
  )
  stateFixture <- makeProtSummaryStateManager(
    history = c("ruv_corrected"),
    objectsByState = list(ruv_corrected = stateObject)
  )
  messages <- captureProtSummaryMessages()

  resolved <- resolveProtSummaryReportTemplate(
    workflowData = list(
      config_list = NULL,
      state_manager = stateFixture$manager
    ),
    catFn = messages$catFn
  )

  expect_identical(resolved$workflowTypeDetected, "LFQ")
  expect_identical(resolved$templateFilename, "LFQ_report.rmd")
  expect_identical(resolved$dataStateUsed, "ruv_corrected")
  expect_identical(stateFixture$requestedStates(), "ruv_corrected")
  expect_true(any(grepl(
    "Detected workflow_type from S4 object (state: ruv_corrected): LFQ",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Selected template: LFQ_report.rmd for workflow_type: LFQ",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("ensureProtSummaryReportTemplate copies the packaged template into the project", {
  fixtureDir <- tempfile("prot-summary-template-helper-")
  projectBaseDir <- file.path(fixtureDir, "project")
  packageDir <- file.path(fixtureDir, "package")
  templateFilename <- "TMT_report.rmd"
  packageTemplate <- file.path(packageDir, templateFilename)
  dir.create(projectBaseDir, recursive = TRUE)
  dir.create(packageDir, recursive = TRUE)
  writeLines("---", packageTemplate)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  notifications <- list()
  logs <- character()
  messages <- captureProtSummaryMessages()

  templateInfo <- ensureProtSummaryReportTemplate(
    projectDirs = list(proteomics = list(base_dir = projectBaseDir)),
    templateFilename = templateFilename,
    systemFileFn = function(...) packageTemplate,
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      logs <<- c(logs, message)
      invisible(NULL)
    },
    catFn = messages$catFn
  )

  expect_identical(templateInfo$templateSource, "package")
  expect_identical(
    templateInfo$reportTemplatePath,
    file.path(projectBaseDir, "scripts", "proteomics", templateFilename)
  )
  expect_true(file.exists(templateInfo$reportTemplatePath))
  expect_identical(readLines(templateInfo$reportTemplatePath), "---")
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "TMT_report.rmd template copied from package")
  expect_identical(notifications[[1]]$type, "message")
  expect_true(any(grepl("Template copied from package to:", logs, fixed = TRUE)))
  expect_true(any(grepl(
    "Found template in package",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("retrieveProtSummaryReportTemplateAsset delegates template lookup through the seam", {
  projectDirs <- list(proteomics = list(base_dir = tempfile("prot-summary-template-project-")))
  recorded <- new.env(parent = emptyenv())
  templateAsset <- list(
    reportTemplatePath = "/tmp/prot-summary-template.rmd",
    templateSource = "existing",
    templateUrl = NULL
  )

  retrieved <- retrieveProtSummaryReportTemplateAsset(
    projectDirs = projectDirs,
    templateFilename = "TMT_report.rmd",
    ensureTemplateFn = function(projectDirs,
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
      recorded$templateArgs <- list(
        projectDirs = projectDirs,
        templateFilename = templateFilename,
        omicType = omicType,
        packageName = packageName,
        packageReportSubdir = packageReportSubdir
      )
      templateAsset
    },
    showNotificationFn = function(...) stop("helper should not notify on success"),
    logErrorFn = function(...) stop("helper should not log errors on success")
  )

  expect_identical(retrieved, templateAsset)
  expect_identical(recorded$templateArgs$projectDirs, projectDirs)
  expect_identical(recorded$templateArgs$templateFilename, "TMT_report.rmd")
  expect_identical(recorded$templateArgs$omicType, "proteomics")
  expect_identical(recorded$templateArgs$packageName, "MultiScholaR")
  expect_identical(recorded$templateArgs$packageReportSubdir, "proteomics")
})

test_that("retrieveProtSummaryReportTemplateAsset preserves the template-retrieval error shell", {
  notifications <- list()
  logErrors <- character()

  retrieved <- retrieveProtSummaryReportTemplateAsset(
    projectDirs = list(proteomics = list(base_dir = tempfile("prot-summary-template-project-"))),
    templateFilename = "TMT_report.rmd",
    ensureTemplateFn = function(...) stop("download boom"),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    }
  )

  expect_null(retrieved)
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "Template retrieval failed: download boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
  expect_identical(logErrors, "Failed to retrieve template: download boom")
})

test_that("validateProtSummaryProjectDirs accepts initialized project directories", {
  notifications <- list()

  validated <- validateProtSummaryProjectDirs(
    projectDirs = list(proteomics = list(base_dir = tempfile("prot-summary-base-"))),
    omicType = "proteomics",
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    }
  )

  expect_true(validated)
  expect_identical(notifications, list())
})

test_that("initializeProtSummaryDefaultOutputs wires the default summary outputs", {
  output <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())

  initialized <- initializeProtSummaryDefaultOutputs(
    output = output,
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    reactiveFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    outputOptionsFn = function(output, name, suspendWhenHidden = TRUE) {
      recorded$outputOptions <- list(name = name, suspendWhenHidden = suspendWhenHidden)
      invisible(NULL)
    }
  )

  expect_true(initialized)
  expect_identical(
    output$session_summary,
    "Ready to save workflow parameters and generate report"
  )
  expect_false(output$report_ready)
  expect_identical(recorded$outputOptions$name, "report_ready")
  expect_false(recorded$outputOptions$suspendWhenHidden)
})

test_that("validateProtSummaryProjectDirs preserves the project-dir validation error shell", {
  notifications <- list()

  validated <- validateProtSummaryProjectDirs(
    projectDirs = list(proteomics = list(base_dir = NULL)),
    omicType = "proteomics",
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    }
  )

  expect_false(validated)
  expect_identical(length(notifications), 1L)
  expect_identical(
    notifications[[1]]$message,
    "Error: Project directories not properly initialized"
  )
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
})

test_that("activateProtSummaryRenderedReport wires report-ready outputs from a rendered file", {
  fixtureDir <- tempfile("prot-summary-render-helper-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  renderedPath <- file.path(fixtureDir, "summary-report.html")
  writeLines("<html><body>report</body></html>", renderedPath)

  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$report_generated <- FALSE
  values$report_path <- NULL
  recorded <- new.env(parent = emptyenv())
  notifications <- list()
  fixedTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")

  activated <- activateProtSummaryRenderedReport(
    output = output,
    values = values,
    renderedPath = renderedPath,
    experimentLabel = "summary-demo",
    description = "demo run",
    reactiveFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    outputOptionsFn = function(output, name, suspendWhenHidden = TRUE) {
      recorded$outputOptions <- list(name = name, suspendWhenHidden = suspendWhenHidden)
      invisible(NULL)
    },
    downloadHandlerFn = function(filename, content) {
      list(filename = filename, content = content)
    },
    fileCopyFn = function(from, to) {
      recorded$fileCopy <- list(from = from, to = to)
      TRUE
    },
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    timestampFn = function() fixedTimestamp
  )

  expect_true(activated)
  expect_true(values$report_generated)
  expect_identical(values$report_path, renderedPath)
  expect_true(output$report_ready)
  expect_identical(recorded$outputOptions$name, "report_ready")
  expect_false(recorded$outputOptions$suspendWhenHidden)
  expect_identical(output$download_report$filename(), "summary-report.html")

  downloadCopyTarget <- tempfile("prot-summary-download-copy-")
  output$download_report$content(downloadCopyTarget)
  expect_identical(recorded$fileCopy$from, renderedPath)
  expect_identical(recorded$fileCopy$to, downloadCopyTarget)

  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "Report generated successfully!")
  expect_identical(notifications[[1]]$type, "message")
  expect_identical(
    output$session_summary,
    paste("Workflow args created for:", "summary-demo",
          "\nDescription:", "demo run",
          "\nTimestamp:", fixedTimestamp,
          "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
          "\nReport location:", renderedPath)
  )
})

test_that("activateProtSummaryRenderedReport returns FALSE when the rendered file is missing", {
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$report_generated <- FALSE
  values$report_path <- NULL

  activated <- activateProtSummaryRenderedReport(
    output = output,
    values = values,
    renderedPath = tempfile("prot-summary-missing-render-"),
    showNotificationFn = function(...) stop("notification should not be called"),
    outputOptionsFn = function(...) stop("outputOptions should not be called"),
    downloadHandlerFn = function(...) stop("downloadHandler should not be called"),
    renderTextFn = function(...) stop("renderText should not be called")
  )

  expect_false(activated)
  expect_false(values$report_generated)
  expect_null(values$report_path)
  expect_false(exists("report_ready", envir = output, inherits = FALSE))
  expect_false(exists("download_report", envir = output, inherits = FALSE))
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
})

test_that("runProtSummaryReportGeneration delegates report rendering and activation", {
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())

  rendered <- runProtSummaryReportGeneration(
    output = output,
    values = values,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    templateFilename = "TMT_report.rmd",
    renderReportAvailableFn = function() TRUE,
    renderReportFn = function(...) {
      recorded$renderArgs <- list(...)
      "/tmp/summary-report.html"
    },
    activateReportFn = function(output,
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
      recorded$activateArgs <- list(
        output = output,
        values = values,
        renderedPath = renderedPath,
        experimentLabel = experimentLabel,
        description = description
      )
      TRUE
    },
    showNotificationFn = function(...) stop("helper should not notify on success"),
    logInfoFn = function(message) {
      recorded$logInfo <- c(recorded$logInfo, message)
      invisible(NULL)
    },
    logErrorFn = function(...) stop("helper should not log errors on success"),
    catFn = function(...) stop("helper should not emit debug output on success"),
    printFn = function(...) stop("helper should not print on success"),
    tracebackFn = function(...) stop("helper should not traceback on success")
  )

  expect_true(rendered)
  expect_identical(
    recorded$renderArgs,
    list(
      omic_type = "proteomics",
      experiment_label = "summary-demo",
      rmd_filename = "TMT_report.rmd"
    )
  )
  expect_identical(recorded$activateArgs$output, output)
  expect_identical(recorded$activateArgs$values, values)
  expect_identical(recorded$activateArgs$renderedPath, "/tmp/summary-report.html")
  expect_identical(recorded$activateArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$activateArgs$description, "demo run")
  expect_identical(
    recorded$logInfo,
    c(
      "Calling RenderReport with omic_type: proteomics, experiment_label: summary-demo",
      "RenderReport returned path: /tmp/summary-report.html"
    )
  )
})

test_that("runProtSummaryReportGeneration preserves the report-generation error shell", {
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  notifications <- list()
  logErrors <- character()
  messages <- character()
  printed <- list()
  tracebackCalled <- FALSE

  rendered <- runProtSummaryReportGeneration(
    output = output,
    values = values,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    templateFilename = "TMT_report.rmd",
    renderReportAvailableFn = function() TRUE,
    renderReportFn = function(...) stop("render boom"),
    activateReportFn = function(...) stop("activation should not be called"),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    },
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    },
    printFn = function(object) {
      printed <<- c(printed, list(conditionMessage(object)))
      invisible(NULL)
    },
    tracebackFn = function() {
      tracebackCalled <<- TRUE
      invisible(NULL)
    }
  )

  expect_false(rendered)
  expect_true(any(grepl("REPORT GENERATION ERROR", messages, fixed = TRUE)))
  expect_true(any(grepl("Error class: simpleError", messages, fixed = TRUE)))
  expect_true(any(grepl("Error message: render boom", messages, fixed = TRUE)))
  expect_identical(printed, list("render boom"))
  expect_identical(
    logErrors,
    c(
      "Failed to generate report: render boom",
      "Error class: simpleError"
    )
  )
  expect_true(tracebackCalled)
  expect_identical(length(notifications), 2L)
  expect_identical(notifications[[1]]$message, "Report generation failed: render boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 15)
  expect_identical(
    notifications[[2]]$message,
    "Debug info: Check R console for detailed error trace"
  )
  expect_identical(notifications[[2]]$type, "warning")
  expect_identical(notifications[[2]]$duration, 10)
})

test_that("runProtSummaryReportProgress delegates report progress orchestration through seams", {
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  progress <- list()

  progressed <- runProtSummaryReportProgress(
    output = output,
    values = values,
    workflowData = list(config_list = list(globalParameters = list(workflow_type = "TMT"))),
    projectDirs = list(proteomics = list(base_dir = tempfile("prot-summary-progress-"))),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    resolveReportTemplateFn = function(workflowData,
                                       logPrefix = "REPORT",
                                       catFn = cat) {
      recorded$resolveArgs <- list(
        workflowData = workflowData,
        logPrefix = logPrefix
      )
      list(
        workflowTypeDetected = "TMT",
        templateFilename = "TMT_report.rmd",
        dataStateUsed = "ruv_corrected"
      )
    },
    retrieveTemplateAssetFn = function(projectDirs,
                                       templateFilename,
                                       omicType = "proteomics",
                                       packageName = "MultiScholaR",
                                       packageReportSubdir = "proteomics",
                                       ensureTemplateFn = ensureProtSummaryReportTemplate,
                                       showNotificationFn = shiny::showNotification,
                                       logErrorFn = function(message) {
                                         logger::log_error(logger::skip_formatter(message))
                                       }) {
      recorded$retrieveArgs <- list(
        projectDirs = projectDirs,
        templateFilename = templateFilename,
        omicType = omicType
      )
      list(
        reportTemplatePath = tempfile("prot-summary-template-"),
        templateSource = "existing",
        templateUrl = NULL
      )
    },
    runReportGenerationFn = function(output,
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
      recorded$runArgs <- list(
        output = output,
        values = values,
        omicType = omicType,
        experimentLabel = experimentLabel,
        description = description,
        templateFilename = templateFilename
      )
      TRUE
    },
    incProgressFn = function(amount, detail = NULL, message = NULL, value = NULL) {
      progress <<- c(
        progress,
        list(list(amount = amount, detail = detail, message = message, value = value))
      )
      invisible(NULL)
    }
  )

  expect_true(progressed)
  expect_identical(recorded$resolveArgs$workflowData$config_list$globalParameters$workflow_type, "TMT")
  expect_identical(recorded$resolveArgs$logPrefix, "REPORT")
  expect_identical(recorded$retrieveArgs$templateFilename, "TMT_report.rmd")
  expect_identical(recorded$retrieveArgs$omicType, "proteomics")
  expect_identical(recorded$runArgs$output, output)
  expect_identical(recorded$runArgs$values, values)
  expect_identical(recorded$runArgs$omicType, "proteomics")
  expect_identical(recorded$runArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$runArgs$description, "demo run")
  expect_identical(recorded$runArgs$templateFilename, "TMT_report.rmd")
  expect_identical(
    progress,
    list(
      list(amount = 0.1, detail = "Detecting workflow type...", message = NULL, value = NULL),
      list(amount = 0.2, detail = "Checking for TMT_report.rmd template...", message = NULL, value = NULL),
      list(amount = 0.5, detail = "Rendering report...", message = NULL, value = NULL)
    )
  )
})

test_that("runProtSummaryReportProgress stops before rendering when template retrieval fails", {
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  progress <- list()
  renderCalled <- FALSE

  progressed <- runProtSummaryReportProgress(
    output = output,
    values = values,
    workflowData = list(config_list = list(globalParameters = list(workflow_type = "LFQ"))),
    projectDirs = list(proteomics = list(base_dir = tempfile("prot-summary-progress-fail-"))),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    resolveReportTemplateFn = function(...) {
      list(
        workflowTypeDetected = "LFQ",
        templateFilename = "LFQ_report.rmd",
        dataStateUsed = NULL
      )
    },
    retrieveTemplateAssetFn = function(...) NULL,
    runReportGenerationFn = function(...) {
      renderCalled <<- TRUE
      TRUE
    },
    incProgressFn = function(amount, detail = NULL, message = NULL, value = NULL) {
      progress <<- c(
        progress,
        list(list(amount = amount, detail = detail, message = message, value = value))
      )
      invisible(NULL)
    }
  )

  expect_false(progressed)
  expect_false(renderCalled)
  expect_identical(
    progress,
    list(
      list(amount = 0.1, detail = "Detecting workflow type...", message = NULL, value = NULL),
      list(amount = 0.2, detail = "Checking for LFQ_report.rmd template...", message = NULL, value = NULL)
    )
  )
})

test_that("completeProtSummaryWorkflowArgsSave preserves the save success shell", {
  fixtureDir <- tempfile("prot-summary-save-helper-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  recorded <- new.env(parent = emptyenv())
  recorded$notifications <- list()
  fixedTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
  finalS4Object <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )
  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )

  withCleanGlobalObjects("config_list", {
    completed <- completeProtSummaryWorkflowArgsSave(
      output = output,
      values = values,
      projectDirs = projectDirs,
      workflowData = workflowData,
      omicType = "proteomics",
      experimentLabel = "summary-demo",
      description = "demo run",
      resolveFinalS4StateFn = function(workflowData,
                                       logPrefix = "SESSION SUMMARY",
                                       catFn = cat) {
        recorded$resolveArgs <- list(workflowData = workflowData, logPrefix = logPrefix)
        list(finalS4Object = finalS4Object, dataStateUsed = "correlation_filtered")
      },
      assignFn = function(x, value, envir = .GlobalEnv) {
        recorded$assigned <- list(name = x, value = value, envir = envir)
        base::assign(x, value, envir = envir)
        invisible(value)
      },
      createWorkflowArgsFn = function(...) {
        recorded$workflowArgs <- list(...)
        file.path(sourceDir, "study_parameters.txt")
      },
      saveRDSFn = function(object, file) {
        recorded$saveRDS <- list(object = object, file = file)
        invisible(NULL)
      },
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, type = "default", duration = NULL) {
        recorded$notifications <- c(
          recorded$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      dirExistsFn = function(path) {
        recorded$dirExists <- c(recorded$dirExists, path)
        FALSE
      },
      dirCreateFn = function(path, recursive = FALSE, showWarnings = TRUE) {
        recorded$dirCreate <- list(
          path = path,
          recursive = recursive,
          showWarnings = showWarnings
        )
        invisible(TRUE)
      },
      timestampFn = function() fixedTimestamp,
      catFn = function(...) invisible(NULL)
    )

    expect_true(completed)
    expect_true(exists("config_list", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(get("config_list", envir = .GlobalEnv), workflowData$config_list)
  })

  expect_true(values$workflow_args_saved)
  expect_identical(recorded$resolveArgs$workflowData, workflowData)
  expect_identical(recorded$resolveArgs$logPrefix, "SESSION SUMMARY")
  expect_identical(recorded$assigned$name, "config_list")
  expect_identical(recorded$assigned$value, workflowData$config_list)
  expect_identical(recorded$workflowArgs$workflow_name, "summary-demo")
  expect_identical(recorded$workflowArgs$description, "demo run")
  expect_identical(recorded$workflowArgs$source_dir_path, sourceDir)
  expect_identical(recorded$workflowArgs$final_s4_object, finalS4Object)
  expect_identical(recorded$workflowArgs$contrasts_tbl, workflowData$contrasts_tbl)
  expect_identical(recorded$workflowArgs$workflow_data, workflowData)
  expect_identical(recorded$dirExists, integrationDir)
  expect_identical(recorded$dirCreate$path, integrationDir)
  expect_true(recorded$dirCreate$recursive)
  expect_false(recorded$dirCreate$showWarnings)
  expect_identical(recorded$saveRDS$object, finalS4Object)
  expect_identical(
    recorded$saveRDS$file,
    file.path(integrationDir, "proteomics_summary-demo_final_s4.RDS")
  )
  expect_identical(length(recorded$notifications), 2L)
  expect_identical(recorded$notifications[[1]]$message, "Saved Integration S4 Object")
  expect_identical(recorded$notifications[[1]]$type, "message")
  expect_identical(recorded$notifications[[2]]$message, "Study parameters saved successfully")
  expect_identical(recorded$notifications[[2]]$type, "message")
  expect_match(output$session_summary, "Study parameters created for: summary-demo", fixed = TRUE)
  expect_match(output$session_summary, "Description: demo run", fixed = TRUE)
  expect_match(output$session_summary, paste("Timestamp:", fixedTimestamp), fixed = TRUE)
  expect_match(
    output$session_summary,
    file.path(sourceDir, "study_parameters.txt"),
    fixed = TRUE
  )
  expect_match(
    output$session_summary,
    "Integration Object: proteomics_summary-demo_final_s4.RDS",
    fixed = TRUE
  )
  expect_match(
    output$session_summary,
    "Status: Parameters saved [OK]",
    fixed = TRUE
  )
})

test_that("completeProtSummaryWorkflowArgsSave preserves the save warning shell", {
  fixtureDir <- tempfile("prot-summary-save-helper-error-")
  sourceDir <- file.path(fixtureDir, "source")
  dir.create(sourceDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  notifications <- list()
  messages <- character()
  recorded <- new.env(parent = emptyenv())
  fixedTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
  finalS4Object <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )

  completed <- completeProtSummaryWorkflowArgsSave(
    output = output,
    values = values,
    projectDirs = projectDirs,
    workflowData = list(
      state_manager = list(dummy = TRUE),
      contrasts_tbl = NULL,
      config_list = NULL
    ),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    resolveFinalS4StateFn = function(...) {
      list(finalS4Object = finalS4Object, dataStateUsed = "correlation_filtered")
    },
    createWorkflowArgsFn = function(...) stop("save boom"),
    saveRDSFn = function(...) stop("saveRDSFn should not be called on failure"),
    renderTextFn = function(...) stop("renderTextFn should not be called on failure"),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    fileExistsFn = function(path) {
      recorded$fileExists <- c(recorded$fileExists, path)
      FALSE
    },
    writeLinesFn = function(text, con) {
      recorded$fallback <- list(text = text, con = con)
      invisible(NULL)
    },
    timestampFn = function() fixedTimestamp,
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    }
  )

  expect_false(completed)
  expect_true(values$workflow_args_saved)
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
  expect_identical(
    recorded$fileExists,
    file.path(sourceDir, "study_parameters.txt")
  )
  expect_identical(recorded$fallback$con, file.path(sourceDir, "study_parameters.txt"))
  expect_identical(
    recorded$fallback$text,
    c(
      "Study Parameters",
      "================",
      "",
      "Workflow Name: summary-demo",
      "Description: demo run",
      paste("Timestamp:", fixedTimestamp),
      "Error: save boom"
    )
  )
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "Study parameters saved with warnings")
  expect_identical(notifications[[1]]$type, "warning")
  expect_true(any(grepl("SESSION SUMMARY ERROR:save boom", messages, fixed = TRUE)))
  expect_true(any(grepl("Created basic fallback file", messages, fixed = TRUE)))
})

test_that("prepareProtSummarySessionStateExport builds export path and payload", {
  exportTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
  projectDirs <- list(
    proteomics = list(
      source_dir = tempfile("prot-summary-export-source-"),
      integration_dir = tempfile("prot-summary-export-integration-"),
      base_dir = tempfile("prot-summary-export-base-")
    )
  )

  exportPayload <- prepareProtSummarySessionStateExport(
    projectDirs = projectDirs,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = FALSE,
    reportPath = NULL,
    exportDate = as.Date("2026-04-16"),
    timestamp = exportTimestamp
  )

  expect_identical(
    exportPayload$sessionExportPath,
    file.path(projectDirs$proteomics$source_dir, "session_state_2026-04-16.RDS")
  )
  expect_identical(exportPayload$sessionState$experiment_label, "summary-demo")
  expect_identical(exportPayload$sessionState$description, "demo run")
  expect_identical(exportPayload$sessionState$timestamp, exportTimestamp)
  expect_true(exportPayload$sessionState$workflow_args_saved)
  expect_true(exportPayload$sessionState$files_copied)
  expect_false(exportPayload$sessionState$report_generated)
  expect_null(exportPayload$sessionState$report_path)
  expect_identical(exportPayload$sessionState$project_dirs, projectDirs)
})

test_that("completeProtSummarySessionStateExport preserves the export success shell", {
  projectDirs <- list(
    proteomics = list(
      source_dir = tempfile("prot-summary-export-complete-source-"),
      integration_dir = tempfile("prot-summary-export-complete-integration-"),
      base_dir = tempfile("prot-summary-export-complete-base-")
    )
  )
  recorded <- new.env(parent = emptyenv())
  exportPayload <- list(
    sessionExportPath = file.path(projectDirs$proteomics$source_dir, "session_state_2026-04-16.RDS"),
    sessionState = list(exported = TRUE)
  )

  completed <- completeProtSummarySessionStateExport(
    projectDirs = projectDirs,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = TRUE,
    reportPath = "/tmp/summary-report.html",
    prepareExportFn = function(projectDirs,
                               omicType = "proteomics",
                               experimentLabel,
                               description,
                               workflowArgsSaved,
                               filesCopied,
                               reportGenerated,
                               reportPath,
                               exportDate = Sys.Date(),
                               timestamp = Sys.time()) {
      recorded$prepareArgs <- list(
        projectDirs = projectDirs,
        omicType = omicType,
        experimentLabel = experimentLabel,
        description = description,
        workflowArgsSaved = workflowArgsSaved,
        filesCopied = filesCopied,
        reportGenerated = reportGenerated,
        reportPath = reportPath
      )
      exportPayload
    },
    saveRDSFn = function(object, file) {
      recorded$saveRDS <- list(object = object, file = file)
      invisible(NULL)
    },
    showNotificationFn = function(message, type = "default", duration = NULL) {
      recorded$notifications <- c(
        recorded$notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      recorded$logInfo <- c(recorded$logInfo, message)
      invisible(NULL)
    },
    logErrorFn = function(...) stop("logErrorFn should not be called on success")
  )

  expect_true(completed)
  expect_identical(recorded$prepareArgs$projectDirs, projectDirs)
  expect_identical(recorded$prepareArgs$omicType, "proteomics")
  expect_identical(recorded$prepareArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$prepareArgs$description, "demo run")
  expect_true(recorded$prepareArgs$workflowArgsSaved)
  expect_true(recorded$prepareArgs$filesCopied)
  expect_true(recorded$prepareArgs$reportGenerated)
  expect_identical(recorded$prepareArgs$reportPath, "/tmp/summary-report.html")
  expect_identical(recorded$saveRDS$object, exportPayload$sessionState)
  expect_identical(recorded$saveRDS$file, exportPayload$sessionExportPath)
  expect_identical(length(recorded$notifications), 1L)
  expect_identical(
    recorded$notifications[[1]]$message,
    paste("Session state exported to:", exportPayload$sessionExportPath)
  )
  expect_identical(recorded$notifications[[1]]$type, "message")
  expect_null(recorded$notifications[[1]]$duration)
  expect_identical(
    recorded$logInfo,
    paste("Session state exported to:", exportPayload$sessionExportPath)
  )
})

test_that("completeProtSummarySessionStateExport preserves the export error shell", {
  notifications <- list()
  logErrors <- character()

  completed <- completeProtSummarySessionStateExport(
    projectDirs = list(proteomics = list(source_dir = tempfile("prot-summary-export-error-"))),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = TRUE,
    reportPath = "/tmp/summary-report.html",
    prepareExportFn = function(...) {
      list(
        sessionExportPath = tempfile("prot-summary-export-error-target-"),
        sessionState = list(exported = TRUE)
      )
    },
    saveRDSFn = function(...) stop("save boom"),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logInfoFn = function(...) stop("logInfoFn should not be called on error"),
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    }
  )

  expect_false(completed)
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "Export failed: save boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
  expect_identical(logErrors, "Failed to export session state: save boom")
})

test_that("bootstrapProtSummaryCopyFallbackStudyParams creates a basic study-parameters fallback file", {
  fixtureDir <- tempfile("prot-summary-copy-bootstrap-")
  sourceDir <- file.path(fixtureDir, "source")
  dir.create(sourceDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  messages <- captureProtSummaryMessages()
  fixedTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")

  created <- bootstrapProtSummaryCopyFallbackStudyParams(
    values = values,
    projectDirs = list(proteomics = list(source_dir = sourceDir)),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    timestampFn = function() fixedTimestamp,
    logErrorFn = function(...) stop("helper should not log errors on success"),
    catFn = messages$catFn
  )

  expect_true(created)
  expect_true(values$workflow_args_saved)
  expect_identical(
    readLines(file.path(sourceDir, "study_parameters.txt")),
    c(
      "Study Parameters",
      "================",
      "",
      "Workflow Name: summary-demo",
      "Description: demo run",
      paste("Timestamp:", fixedTimestamp),
      "Note: Some parameters could not be saved due to serialization issues"
    )
  )
  expect_true(any(grepl(
    "Creating basic study_parameters.txt as fallback",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("bootstrapProtSummaryCopyFallbackStudyParams preserves the fallback-write error shell", {
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  logErrors <- character()

  created <- bootstrapProtSummaryCopyFallbackStudyParams(
    values = values,
    projectDirs = list(proteomics = list(source_dir = tempfile("prot-summary-copy-bootstrap-fail-"))),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    writeLinesFn = function(...) stop("disk full"),
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    },
    catFn = function(...) invisible(NULL)
  )

  expect_false(created)
  expect_false(values$workflow_args_saved)
  expect_identical(logErrors, "Failed to create basic study_parameters.txt: disk full")
})

test_that("prepareProtSummaryCopyInputs prefers workflow_data tables over file fallbacks", {
  messages <- captureProtSummaryMessages()
  contrastsTbl <- data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
  designMatrix <- data.frame(sample = "S1", stringsAsFactors = FALSE)

  copyInputs <- prepareProtSummaryCopyInputs(
    workflowData = list(
      contrasts_tbl = contrastsTbl,
      design_matrix = designMatrix
    ),
    projectDirs = list(proteomics = list(source_dir = tempfile("prot-summary-copy-source-"))),
    readTsvFn = function(...) stop("file fallback should not be used"),
    catFn = messages$catFn
  )

  expect_identical(copyInputs$contrastsTbl, contrastsTbl)
  expect_identical(copyInputs$designMatrix, designMatrix)
  expect_true(any(grepl(
    "Got contrasts_tbl from workflow_data",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Got design_matrix from workflow_data",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("prepareProtSummaryCopyInputs falls back to source files when workflow_data tables are missing", {
  fixtureDir <- tempfile("prot-summary-copy-fallback-")
  sourceDir <- file.path(fixtureDir, "source")
  dir.create(sourceDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  designMatrixFile <- file.path(sourceDir, "design_matrix.tab")
  contrastsFile <- file.path(sourceDir, "contrasts_tbl.tab")
  file.create(designMatrixFile)
  file.create(contrastsFile)

  messages <- captureProtSummaryMessages()
  requestedPaths <- character()

  copyInputs <- prepareProtSummaryCopyInputs(
    workflowData = list(),
    projectDirs = list(proteomics = list(source_dir = sourceDir)),
    readTsvFn = function(path, show_col_types = FALSE) {
      requestedPaths <<- c(requestedPaths, basename(path))

      if (identical(basename(path), "design_matrix.tab")) {
        return(data.frame(sample = "S1", stringsAsFactors = FALSE))
      }

      data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
    },
    catFn = messages$catFn
  )

  expect_identical(requestedPaths, c("design_matrix.tab", "contrasts_tbl.tab"))
  expect_identical(copyInputs$designMatrix$sample, "S1")
  expect_identical(copyInputs$contrastsTbl$contrast, "A_vs_B")
  expect_true(any(grepl(
    "Loaded design_matrix from file",
    messages$getMessages(),
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "Loaded contrasts_tbl from file",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("runProtSummaryPublicationCopy bootstraps copy execution and applies success updates", {
  fixtureDir <- tempfile("prot-summary-copy-helper-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$files_copied <- FALSE
  recorded <- new.env(parent = emptyenv())
  fixedTimestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
  delegatedContrasts <- data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
  delegatedDesign <- data.frame(sample = "S1", stringsAsFactors = FALSE)

  withCleanGlobalObjects("project_dirs", {
    if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "project_dirs", envir = .GlobalEnv)
    }

    recorded$copyArgs <- runProtSummaryPublicationCopy(
      output = output,
      values = values,
      projectDirs = projectDirs,
      omicType = "proteomics",
      experimentLabel = "summary-demo",
      description = "demo run",
      contrastsTbl = delegatedContrasts,
      designMatrix = delegatedDesign,
      assignFn = function(x, value, envir = .GlobalEnv) {
        recorded$assigned <- list(name = x, value = value, envir = envir)
        base::assign(x, value, envir = envir)
        invisible(value)
      },
      copyFn = function(...) {
        recorded$delegated <- list(...)
        invisible(NULL)
      },
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, type = "default", duration = NULL) {
        recorded$notifications <- c(
          recorded$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      timestampFn = function() fixedTimestamp,
      catFn = function(...) invisible(NULL)
    )

    expect_true(exists("project_dirs", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(get("project_dirs", envir = .GlobalEnv), projectDirs)
  })

  expect_identical(recorded$assigned$name, "project_dirs")
  expect_identical(recorded$assigned$value, projectDirs)
  expect_identical(
    recorded$copyArgs,
    list(
      omic_type = "proteomics",
      experiment_label = "summary-demo",
      contrasts_tbl = delegatedContrasts,
      design_matrix = delegatedDesign,
      force = TRUE
    )
  )
  expect_identical(recorded$delegated, recorded$copyArgs)
  expect_true(values$files_copied)
  expect_identical(
    output$copy_status,
    "Files copied to publication directory successfully [OK]"
  )
  expect_identical(length(recorded$notifications), 1L)
  expect_identical(recorded$notifications[[1]]$message, "Publication files copied")
  expect_identical(recorded$notifications[[1]]$type, "message")
  expect_match(output$session_summary, "Workflow args created for: summary-demo", fixed = TRUE)
  expect_match(output$session_summary, "Description: demo run", fixed = TRUE)
  expect_match(output$session_summary, paste("Timestamp:", fixedTimestamp), fixed = TRUE)
  expect_match(
    output$session_summary,
    "Status: Arguments saved [OK], Files copied [OK]",
    fixed = TRUE
  )
})

test_that("observeProtSummaryPublicationCopy wires the copy observer through the existing seams", {
  fixtureDir <- tempfile("prot-summary-copy-observer-helper-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(state_manager = list(dummy = TRUE))
  values <- shiny::reactiveValues(
    workflow_args_saved = FALSE,
    files_copied = FALSE,
    report_generated = FALSE,
    report_path = NULL
  )
  recorded <- new.env(parent = emptyenv())
  delegatedContrasts <- data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
  delegatedDesign <- data.frame(sample = "S1", stringsAsFactors = FALSE)

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummaryPublicationCopy,
    list(
      bootstrapProtSummaryCopyFallbackStudyParams = function(values,
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
        recorded$events <- c(recorded$events, "bootstrap")
        recorded$bootstrapArgs <- list(
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description
        )
        values$workflow_args_saved <- TRUE
        TRUE
      },
      prepareProtSummaryCopyInputs = function(workflowData,
                                              projectDirs,
                                              omicType = "proteomics",
                                              readTsvFn = readr::read_tsv,
                                              catFn = cat) {
        recorded$events <- c(recorded$events, "prepare")
        recorded$prepareArgs <- list(
          workflowData = workflowData,
          projectDirs = projectDirs,
          omicType = omicType
        )
        list(contrastsTbl = delegatedContrasts, designMatrix = delegatedDesign)
      },
      runProtSummaryPublicationCopy = function(output,
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
        recorded$events <- c(recorded$events, "copy")
        recorded$copyArgs <- list(
          output = output,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          contrastsTbl = contrastsTbl,
          designMatrix = designMatrix
        )
        invisible(NULL)
      },
      handleProtSummaryPublicationCopyError = function(...) {
        stop("copy observer helper should not enter the error seam on success")
      }
    )
  )

  testServer(
    function(input, output, session) {
      helperUnderTest(
        input = input,
        output = output,
        values = values,
        projectDirs = projectDirs,
        workflowData = workflowData
      )
    },
    {
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(copy_to_publication = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$events, c("bootstrap", "prepare", "copy"))
  expect_identical(recorded$bootstrapArgs$values, values)
  expect_identical(recorded$bootstrapArgs$projectDirs, projectDirs)
  expect_identical(recorded$bootstrapArgs$omicType, "proteomics")
  expect_identical(recorded$bootstrapArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$bootstrapArgs$description, "demo run")
  expect_identical(recorded$prepareArgs$workflowData, workflowData)
  expect_identical(recorded$prepareArgs$projectDirs, projectDirs)
  expect_identical(recorded$prepareArgs$omicType, "proteomics")
  expect_identical(recorded$copyArgs$output, recorded$expectedOutput)
  expect_identical(recorded$copyArgs$values, values)
  expect_identical(recorded$copyArgs$projectDirs, projectDirs)
  expect_identical(recorded$copyArgs$omicType, "proteomics")
  expect_identical(recorded$copyArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$copyArgs$description, "demo run")
  expect_identical(recorded$copyArgs$contrastsTbl, delegatedContrasts)
  expect_identical(recorded$copyArgs$designMatrix, delegatedDesign)
})

test_that("handleProtSummaryPublicationCopyError preserves the copy error shell", {
  output <- new.env(parent = emptyenv())
  notifications <- list()
  logErrors <- character()
  messages <- character()
  tracebackCalled <- FALSE

  handled <- handleProtSummaryPublicationCopyError(
    output = output,
    error = simpleError("copy boom"),
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    },
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    },
    tracebackFn = function() {
      tracebackCalled <<- TRUE
      invisible(NULL)
    }
  )

  expect_false(handled)
  expect_identical(output$copy_status, "Error: copy boom")
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "Copy error: copy boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
  expect_identical(logErrors, "Failed to copy files: copy boom")
  expect_true(any(grepl("SESSION SUMMARY ERROR:", messages, fixed = TRUE)))
  expect_true(any(grepl("copy boom", messages, fixed = TRUE)))
  expect_true(any(grepl("SESSION SUMMARY Traceback:", messages, fixed = TRUE)))
  expect_true(tracebackCalled)
})

test_that("runProtSummaryGithubPush configures GitHub options and delegates the push call", {
  fixtureDir <- tempfile("prot-summary-github-helper-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  githubPush <- runProtSummaryGithubPush(
    projectDirs = projectDirs,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    githubOrg = "APAF-bioinformatics",
    githubEmail = "summary@example.org",
    githubUsername = "summary-user",
    projectId = "MSR-42",
    optionsFn = function(...) {
      recorded$options <- list(...)
      invisible(NULL)
    },
    pushFn = function(...) {
      recorded$pushArgs <- list(...)
      invisible(NULL)
    }
  )

  expect_identical(
    githubPush$githubOptions,
    list(
      github_org = "APAF-bioinformatics",
      github_user_email = "summary@example.org",
      github_user_name = "summary-user"
    )
  )
  expect_identical(recorded$options, githubPush$githubOptions)
  expect_identical(
    githubPush$pushArgs,
    list(
      project_dirs = projectDirs,
      omic_type = "proteomics",
      experiment_label = "summary-demo",
      project_id = "MSR-42"
    )
  )
  expect_identical(recorded$pushArgs, githubPush$pushArgs)
})

test_that("completeProtSummaryGithubPush preserves the GitHub success shell", {
  fixtureDir <- tempfile("prot-summary-github-complete-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  output <- new.env(parent = emptyenv())
  recorded <- new.env(parent = emptyenv())
  fixedTimestamp <- as.POSIXct("2026-04-16 14:23:45", tz = "UTC")

  completed <- completeProtSummaryGithubPush(
    output = output,
    projectDirs = projectDirs,
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    githubOrg = "APAF-bioinformatics",
    githubEmail = "summary@example.org",
    githubUsername = "summary-user",
    projectId = "MSR-42",
    pushGithubFn = function(projectDirs,
                            omicType = "proteomics",
                            experimentLabel,
                            githubOrg,
                            githubEmail,
                            githubUsername,
                            projectId) {
      recorded$pushArgs <- list(
        projectDirs = projectDirs,
        omicType = omicType,
        experimentLabel = experimentLabel,
        githubOrg = githubOrg,
        githubEmail = githubEmail,
        githubUsername = githubUsername,
        projectId = projectId
      )
      invisible(list(status = "ok"))
    },
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      recorded$notifications <- c(
        recorded$notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    timestampFn = function() fixedTimestamp,
    logErrorFn = function(...) stop("logErrorFn should not be called on success")
  )

  expect_true(completed)
  expect_identical(recorded$pushArgs$projectDirs, projectDirs)
  expect_identical(recorded$pushArgs$omicType, "proteomics")
  expect_identical(recorded$pushArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$pushArgs$githubOrg, "APAF-bioinformatics")
  expect_identical(recorded$pushArgs$githubEmail, "summary@example.org")
  expect_identical(recorded$pushArgs$githubUsername, "summary-user")
  expect_identical(recorded$pushArgs$projectId, "MSR-42")
  expect_identical(length(recorded$notifications), 1L)
  expect_identical(recorded$notifications[[1]]$message, "Successfully pushed to GitHub")
  expect_identical(recorded$notifications[[1]]$type, "message")
  expect_match(output$session_summary, "Workflow args created for: summary-demo", fixed = TRUE)
  expect_match(output$session_summary, "Description: demo run", fixed = TRUE)
  expect_match(output$session_summary, paste("Timestamp:", fixedTimestamp), fixed = TRUE)
  expect_match(
    output$session_summary,
    "Status: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]",
    fixed = TRUE
  )
})

test_that("observeProtSummaryGithubPush wires the GitHub observer through the completion seam", {
  fixtureDir <- tempfile("prot-summary-github-observer-helper-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummaryGithubPush,
    list(
      completeProtSummaryGithubPush = function(output,
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
        recorded$completionArgs <- list(
          output = output,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          githubOrg = githubOrg,
          githubEmail = githubEmail,
          githubUsername = githubUsername,
          projectId = projectId
        )
        invisible(NULL)
      }
    )
  )

  serverUnderTest <- function(input, output, session) {
    values <- shiny::reactiveValues(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
    recorded$values <- values

    helperUnderTest(
      input = input,
      output = output,
      values = values,
      projectDirs = projectDirs
    )
  }

  testServer(
    serverUnderTest,
    {
      recorded$values$report_generated <- TRUE
      session$setInputs(
        experiment_label = "summary-demo",
        description = "demo run",
        enable_github = TRUE,
        github_org = "APAF-bioinformatics",
        github_email = "summary@example.org",
        github_username = "summary-user",
        project_id = "MSR-42"
      )
      session$flushReact()

      recorded$expectedOutput <- output
      session$setInputs(push_to_github = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$completionArgs$output, recorded$expectedOutput)
  expect_identical(recorded$completionArgs$projectDirs, projectDirs)
  expect_identical(recorded$completionArgs$omicType, "proteomics")
  expect_identical(recorded$completionArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$completionArgs$description, "demo run")
  expect_identical(recorded$completionArgs$githubOrg, "APAF-bioinformatics")
  expect_identical(recorded$completionArgs$githubEmail, "summary@example.org")
  expect_identical(recorded$completionArgs$githubUsername, "summary-user")
  expect_identical(recorded$completionArgs$projectId, "MSR-42")
})

test_that("observeProtSummarySessionStateExport wires the export observer through the completion seam", {
  fixtureDir <- tempfile("prot-summary-export-observer-helper-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummarySessionStateExport,
    list(
      completeProtSummarySessionStateExport = function(projectDirs,
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
        recorded$completionArgs <- list(
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          workflowArgsSaved = workflowArgsSaved,
          filesCopied = filesCopied,
          reportGenerated = reportGenerated,
          reportPath = reportPath
        )
        invisible(NULL)
      }
    )
  )

  serverUnderTest <- function(input, output, session) {
    values <- shiny::reactiveValues(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
    recorded$values <- values

    helperUnderTest(
      input = input,
      values = values,
      projectDirs = projectDirs
    )
  }

  testServer(
    serverUnderTest,
    {
      recorded$values$workflow_args_saved <- TRUE
      recorded$values$files_copied <- TRUE
      recorded$values$report_generated <- TRUE
      recorded$values$report_path <- file.path(fixtureDir, "summary-report.html")
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(export_session_state = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$completionArgs$projectDirs, projectDirs)
  expect_identical(recorded$completionArgs$omicType, "proteomics")
  expect_identical(recorded$completionArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$completionArgs$description, "demo run")
  expect_true(recorded$completionArgs$workflowArgsSaved)
  expect_true(recorded$completionArgs$filesCopied)
  expect_true(recorded$completionArgs$reportGenerated)
  expect_identical(
    recorded$completionArgs$reportPath,
    file.path(fixtureDir, "summary-report.html")
  )
})

test_that("completeProtSummaryGithubPush preserves the GitHub error shell", {
  output <- new.env(parent = emptyenv())
  notifications <- list()
  logErrors <- character()

  completed <- completeProtSummaryGithubPush(
    output = output,
    projectDirs = list(proteomics = list(base_dir = tempfile("prot-summary-github-error-"))),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    githubOrg = "APAF-bioinformatics",
    githubEmail = "summary@example.org",
    githubUsername = "summary-user",
    projectId = "MSR-42",
    pushGithubFn = function(...) stop("push boom"),
    renderTextFn = function(...) stop("renderTextFn should not be called on error"),
    showNotificationFn = function(message, type = "default", duration = NULL) {
      notifications <<- c(
        notifications,
        list(list(message = message, type = type, duration = duration))
      )
      invisible(NULL)
    },
    logErrorFn = function(message) {
      logErrors <<- c(logErrors, message)
      invisible(NULL)
    }
  )

  expect_false(completed)
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
  expect_identical(length(notifications), 1L)
  expect_identical(notifications[[1]]$message, "GitHub push failed: push boom")
  expect_identical(notifications[[1]]$type, "error")
  expect_identical(notifications[[1]]$duration, 10)
  expect_identical(logErrors, "Failed to push to GitHub: push boom")
})

test_that("observeProtSummaryWorkflowArgsSave wires the save observer through the completion seam", {
  fixtureDir <- tempfile("prot-summary-save-observer-helper-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )
  values <- shiny::reactiveValues(
    workflow_args_saved = FALSE,
    files_copied = FALSE,
    report_generated = FALSE,
    report_path = NULL
  )
  recorded <- new.env(parent = emptyenv())
  messages <- captureProtSummaryMessages()

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummaryWorkflowArgsSave,
    list(
      completeProtSummaryWorkflowArgsSave = function(output,
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
        recorded$completionArgs <- list(
          output = output,
          values = values,
          projectDirs = projectDirs,
          workflowData = workflowData,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
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
        invisible(TRUE)
      }
    )
  )

  serverUnderTest <- function(input, output, session) {
    recorded$values <- values

    helperUnderTest(
      input = input,
      output = output,
      values = values,
      projectDirs = projectDirs,
      workflowData = workflowData,
      catFn = messages$catFn
    )
  }

  testServer(
    serverUnderTest,
    {
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$setInputs(save_workflow_args = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$completionArgs$output, recorded$expectedOutput)
  expect_identical(recorded$completionArgs$values, recorded$expectedValues)
  expect_identical(recorded$completionArgs$projectDirs, projectDirs)
  expect_identical(recorded$completionArgs$workflowData, workflowData)
  expect_identical(recorded$completionArgs$omicType, "proteomics")
  expect_identical(recorded$completionArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$completionArgs$description, "demo run")
  expect_true(is.function(recorded$completionArgs$resolveFinalS4StateFn))
  expect_true(is.function(recorded$completionArgs$assignFn))
  expect_true(is.function(recorded$completionArgs$createWorkflowArgsFn))
  expect_true(is.function(recorded$completionArgs$saveRDSFn))
  expect_true(is.function(recorded$completionArgs$renderTextFn))
  expect_true(is.function(recorded$completionArgs$showNotificationFn))
  expect_true(is.function(recorded$completionArgs$dirExistsFn))
  expect_true(is.function(recorded$completionArgs$dirCreateFn))
  expect_true(is.function(recorded$completionArgs$fileExistsFn))
  expect_true(is.function(recorded$completionArgs$writeLinesFn))
  expect_true(is.function(recorded$completionArgs$timestampFn))
  expect_identical(recorded$completionArgs$catFn, messages$catFn)
  expect_true(any(grepl(
    "SESSION SUMMARY: Starting workflow args save process",
    messages$getMessages(),
    fixed = TRUE
  )))
})

test_that("mod_prot_summary_server save observer delegates final S4 lookup through the seam", {
  fixtureDir <- tempfile("prot-summary-save-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  finalS4Object <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )
  recorded <- new.env(parent = emptyenv())

  withCleanGlobalObjects("config_list", {
    serverUnderTest <- makeFunctionWithOverrides(
      mod_prot_summary_server,
      list(
        resolveProtSummaryFinalS4State = function(workflowData,
                                                  logPrefix = "SESSION SUMMARY",
                                                  catFn = cat) {
          recorded$workflowData <- workflowData
          recorded$logPrefix <- logPrefix
          list(finalS4Object = finalS4Object, dataStateUsed = "correlation_filtered")
        },
        createWorkflowArgsFromConfig = function(...) {
          recorded$workflowArgs <- list(...)
          file.path(sourceDir, "study_parameters.txt")
        },
        saveRDS = function(object, file) {
          recorded$saveRDS <- list(object = object, file = file)
          invisible(NULL)
        }
      )
    )

    testServer(
      serverUnderTest,
      args = list(
        project_dirs = projectDirs,
        workflow_data = workflowData
      ),
      {
        session$setInputs(experiment_label = "summary-demo", description = "demo run")
        session$flushReact()

        session$setInputs(save_workflow_args = 1)
        session$flushReact()
      }
    )
  })

  expect_identical(recorded$workflowData, workflowData)
  expect_identical(recorded$logPrefix, "SESSION SUMMARY")
  expect_identical(recorded$workflowArgs$workflow_name, "summary-demo")
  expect_identical(recorded$workflowArgs$description, "demo run")
  expect_identical(recorded$workflowArgs$source_dir_path, sourceDir)
  expect_identical(recorded$workflowArgs$final_s4_object, finalS4Object)
  expect_identical(recorded$saveRDS$object, finalS4Object)
  expect_identical(
    recorded$saveRDS$file,
    file.path(integrationDir, "proteomics_summary-demo_final_s4.RDS")
  )
})

test_that("mod_prot_summary_server save observer delegates completion through the seam", {
  fixtureDir <- tempfile("prot-summary-save-delegate-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      resolveProtSummaryFinalS4State = function(...) {
        stop("resolveProtSummaryFinalS4State should be delegated through the helper")
      },
      createWorkflowArgsFromConfig = function(...) {
        stop("createWorkflowArgsFromConfig should be delegated through the helper")
      },
      saveRDS = function(...) {
        stop("saveRDS should be delegated through the helper")
      },
      completeProtSummaryWorkflowArgsSave = function(output,
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
        recorded$args <- list(
          output = output,
          values = values,
          projectDirs = projectDirs,
          workflowData = workflowData,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          resolveFinalS4StateFn = resolveFinalS4StateFn,
          assignFn = assignFn,
          createWorkflowArgsFn = createWorkflowArgsFn,
          saveRDSFn = saveRDSFn
        )
        values$workflow_args_saved <- TRUE
        TRUE
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(save_workflow_args = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$workflowData, workflowData)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(recorded$args$experimentLabel, "summary-demo")
  expect_identical(recorded$args$description, "demo run")
  expect_true(is.function(recorded$args$resolveFinalS4StateFn))
  expect_true(is.function(recorded$args$assignFn))
  expect_true(is.function(recorded$args$createWorkflowArgsFn))
  expect_true(is.function(recorded$args$saveRDSFn))
})

test_that("mod_prot_summary_server delegates save observer registration through the seam", {
  fixtureDir <- tempfile("prot-summary-save-registration-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
    config_list = list(globalParameters = list(workflow_type = "DIA"))
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      observeProtSummaryWorkflowArgsSave = function(input,
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
        recorded$args <- list(
          input = input,
          output = output,
          values = values,
          projectDirs = projectDirs,
          workflowData = workflowData,
          omicType = omicType,
          completeWorkflowArgsSaveFn = completeWorkflowArgsSaveFn,
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
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedInput <- input
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$flushReact()
    }
  )

  expect_identical(recorded$args$input, recorded$expectedInput)
  expect_identical(recorded$args$output, recorded$expectedOutput)
  expect_identical(recorded$args$values, recorded$expectedValues)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$workflowData, workflowData)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_true(is.function(recorded$args$completeWorkflowArgsSaveFn))
  expect_true(is.function(recorded$args$resolveFinalS4StateFn))
  expect_true(is.function(recorded$args$assignFn))
  expect_true(is.function(recorded$args$createWorkflowArgsFn))
  expect_true(is.function(recorded$args$saveRDSFn))
  expect_true(is.function(recorded$args$renderTextFn))
  expect_true(is.function(recorded$args$showNotificationFn))
  expect_true(is.function(recorded$args$dirExistsFn))
  expect_true(is.function(recorded$args$dirCreateFn))
  expect_true(is.function(recorded$args$fileExistsFn))
  expect_true(is.function(recorded$args$writeLinesFn))
  expect_true(is.function(recorded$args$timestampFn))
  expect_true(is.function(recorded$args$catFn))
})

test_that("registerProtSummaryTemplateStatusOutput renders template status through the seam", {
  fixtureDir <- tempfile("prot-summary-status-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  helperUnderTest <- makeFunctionWithOverrides(
    registerProtSummaryTemplateStatusOutput,
    list(
      buildProtSummaryTemplateStatus = function(projectDirs,
                                                omicType = "proteomics") {
        recorded$projectDirs <- projectDirs
        recorded$omicType <- omicType
        "Templates: delegated status"
      }
    )
  )

  testServer(
    function(input, output, session) {
      helperUnderTest(
        output = output,
        projectDirs = projectDirs
      )
    },
    {
      session$flushReact()

      expect_identical(output$template_status, "Templates: delegated status")
    }
  )

  expect_identical(recorded$projectDirs, projectDirs)
  expect_identical(recorded$omicType, "proteomics")
})

test_that("mod_prot_summary_server delegates template status registration through the seam", {
  fixtureDir <- tempfile("prot-summary-status-registration-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      registerProtSummaryTemplateStatusOutput = function(output,
                                                         projectDirs,
                                                         omicType = "proteomics",
                                                         renderTextFn = shiny::renderText,
                                                         reqFn = shiny::req,
                                                         buildTemplateStatusFn = buildProtSummaryTemplateStatus) {
        recorded$args <- list(
          output = output,
          projectDirs = projectDirs,
          omicType = omicType,
          renderTextFn = renderTextFn,
          reqFn = reqFn,
          buildTemplateStatusFn = buildTemplateStatusFn
        )
        output$template_status <- renderTextFn({
          "Templates: delegated registration"
        })
        invisible(TRUE)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(project_dirs = projectDirs),
    {
      recorded$expectedOutput <- output
      session$flushReact()

      expect_identical(output$template_status, "Templates: delegated registration")
    }
  )

  expect_identical(recorded$args$output, recorded$expectedOutput)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(recorded$args$renderTextFn, shiny::renderText)
  expect_identical(recorded$args$reqFn, shiny::req)
  expect_identical(recorded$args$buildTemplateStatusFn, buildProtSummaryTemplateStatus)
})

test_that("observeProtSummaryReportGeneration wires the report observer through the existing seams", {
  fixtureDir <- tempfile("prot-summary-report-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    config_list = NULL
  )
  values <- shiny::reactiveValues(
    workflow_args_saved = FALSE,
    files_copied = TRUE,
    report_generated = FALSE,
    report_path = NULL
  )
  recorded <- new.env(parent = emptyenv())
  recorded$expectedValues <- values

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummaryReportGeneration,
    list(
      validateProtSummaryProjectDirs = function(projectDirs,
                                                omicType = "proteomics",
                                                showNotificationFn = shiny::showNotification) {
        recorded$validationArgs <- list(
          projectDirs = projectDirs,
          omicType = omicType
        )
        TRUE
      },
      runProtSummaryReportProgress = function(output,
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
        recorded$progressArgs <- list(
          output = output,
          values = values,
          workflowData = workflowData,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description
        )
        TRUE
      }
    )
  )

  testServer(
    function(input, output, session) {
      helperUnderTest(
        input = input,
        output = output,
        values = values,
        projectDirs = projectDirs,
        workflowData = workflowData
      )
    },
    {
      recorded$expectedOutput <- output
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(generate_report = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$validationArgs$projectDirs, projectDirs)
  expect_identical(recorded$validationArgs$omicType, "proteomics")
  expect_identical(recorded$progressArgs$output, recorded$expectedOutput)
  expect_identical(recorded$progressArgs$values, recorded$expectedValues)
  expect_identical(recorded$progressArgs$workflowData, workflowData)
  expect_identical(recorded$progressArgs$projectDirs, projectDirs)
  expect_identical(recorded$progressArgs$omicType, "proteomics")
  expect_identical(recorded$progressArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$progressArgs$description, "demo run")
})

test_that("observeProtSummaryReportGeneration stops when project-dir validation fails", {
  fixtureDir <- tempfile("prot-summary-report-invalid-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  values <- shiny::reactiveValues(
    workflow_args_saved = FALSE,
    files_copied = TRUE,
    report_generated = FALSE,
    report_path = NULL
  )
  recorded <- new.env(parent = emptyenv())

  helperUnderTest <- makeFunctionWithOverrides(
    observeProtSummaryReportGeneration,
    list(
      validateProtSummaryProjectDirs = function(projectDirs,
                                                omicType = "proteomics",
                                                showNotificationFn = shiny::showNotification) {
        recorded$validationArgs <- list(
          projectDirs = projectDirs,
          omicType = omicType
        )
        FALSE
      },
      runProtSummaryReportProgress = function(...) {
        recorded$progressed <- TRUE
        stop("report progress helper should not run")
      }
    )
  )

  testServer(
    function(input, output, session) {
      helperUnderTest(
        input = input,
        output = output,
        values = values,
        projectDirs = projectDirs
      )
    },
    {
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(generate_report = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$validationArgs$projectDirs, projectDirs)
  expect_identical(recorded$validationArgs$omicType, "proteomics")
  expect_false(isTRUE(recorded$progressed))
})

test_that("mod_prot_summary_server delegates report observer registration through the seam", {
  fixtureDir <- tempfile("prot-summary-report-registration-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  workflowData <- list(
    state_manager = list(dummy = TRUE),
    config_list = NULL
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      observeProtSummaryReportGeneration = function(input,
                                                    output,
                                                    values,
                                                    projectDirs,
                                                    workflowData = NULL,
                                                    omicType = "proteomics",
                                                    validateProjectDirsFn = validateProtSummaryProjectDirs,
                                                    runReportProgressFn = runProtSummaryReportProgress,
                                                    withProgressFn = shiny::withProgress) {
        recorded$args <- list(
          input = input,
          output = output,
          values = values,
          projectDirs = projectDirs,
          workflowData = workflowData,
          omicType = omicType,
          validateProjectDirsFn = validateProjectDirsFn,
          runReportProgressFn = runReportProgressFn,
          withProgressFn = withProgressFn
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedInput <- input
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$flushReact()
    }
  )

  expect_identical(recorded$args$input, recorded$expectedInput)
  expect_identical(recorded$args$output, recorded$expectedOutput)
  expect_identical(recorded$args$values, recorded$expectedValues)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$workflowData, workflowData)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(recorded$args$validateProjectDirsFn, validateProtSummaryProjectDirs)
  expect_identical(recorded$args$runReportProgressFn, runProtSummaryReportProgress)
  expect_identical(recorded$args$withProgressFn, shiny::withProgress)
})

test_that("mod_prot_summary_server copy observer delegates fallback bootstrap through the seam", {
  fixtureDir <- tempfile("prot-summary-copy-bootstrap-observer-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(state_manager = list(dummy = TRUE))
  recorded <- new.env(parent = emptyenv())
  delegatedContrasts <- data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
  delegatedDesign <- data.frame(sample = "S1", stringsAsFactors = FALSE)

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      bootstrapProtSummaryCopyFallbackStudyParams = function(values,
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
        recorded$bootstrapArgs <- list(
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description
        )
        recorded$events <- c(recorded$events, "bootstrap")
        values$workflow_args_saved <- TRUE
        TRUE
      },
      prepareProtSummaryCopyInputs = function(workflowData,
                                              projectDirs,
                                              omicType = "proteomics",
                                              readTsvFn = readr::read_tsv,
                                              catFn = cat) {
        list(contrastsTbl = delegatedContrasts, designMatrix = delegatedDesign)
      },
      runProtSummaryPublicationCopy = function(output,
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
        recorded$events <- c(recorded$events, "copy")
        recorded$copyArgs <- list(
          output = output,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          contrastsTbl = contrastsTbl,
          designMatrix = designMatrix
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(copy_to_publication = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$bootstrapArgs$values, recorded$expectedValues)
  expect_identical(recorded$bootstrapArgs$projectDirs, projectDirs)
  expect_identical(recorded$bootstrapArgs$omicType, "proteomics")
  expect_identical(recorded$bootstrapArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$bootstrapArgs$description, "demo run")
  expect_identical(recorded$events, c("bootstrap", "copy"))
  expect_identical(recorded$copyArgs$output, recorded$expectedOutput)
  expect_identical(recorded$copyArgs$values, recorded$expectedValues)
  expect_identical(recorded$copyArgs$projectDirs, projectDirs)
  expect_identical(recorded$copyArgs$omicType, "proteomics")
  expect_identical(recorded$copyArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$copyArgs$description, "demo run")
  expect_identical(recorded$copyArgs$contrastsTbl, delegatedContrasts)
  expect_identical(recorded$copyArgs$designMatrix, delegatedDesign)
})

test_that("mod_prot_summary_server copy observer delegates table resolution into the execution seam", {
  fixtureDir <- tempfile("prot-summary-copy-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(state_manager = list(dummy = TRUE))
  recorded <- new.env(parent = emptyenv())
  delegatedContrasts <- data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
  delegatedDesign <- data.frame(sample = "S1", stringsAsFactors = FALSE)

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      prepareProtSummaryCopyInputs = function(workflowData,
                                              projectDirs,
                                              omicType = "proteomics",
                                              readTsvFn = readr::read_tsv,
                                              catFn = cat) {
        recorded$workflowData <- workflowData
        recorded$projectDirs <- projectDirs
        recorded$omicType <- omicType
        list(contrastsTbl = delegatedContrasts, designMatrix = delegatedDesign)
      },
      runProtSummaryPublicationCopy = function(output,
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
        recorded$copyArgs <- list(
          output = output,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          contrastsTbl = contrastsTbl,
          designMatrix = designMatrix
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      values$workflow_args_saved <- TRUE
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(copy_to_publication = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$workflowData, workflowData)
  expect_identical(recorded$projectDirs, projectDirs)
  expect_identical(recorded$omicType, "proteomics")
  expect_identical(recorded$copyArgs$output, recorded$expectedOutput)
  expect_identical(recorded$copyArgs$values, recorded$expectedValues)
  expect_identical(recorded$copyArgs$projectDirs, projectDirs)
  expect_identical(recorded$copyArgs$omicType, "proteomics")
  expect_identical(recorded$copyArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$copyArgs$description, "demo run")
  expect_identical(recorded$copyArgs$contrastsTbl, delegatedContrasts)
  expect_identical(recorded$copyArgs$designMatrix, delegatedDesign)
})

test_that("mod_prot_summary_server copy observer delegates error handling through the seam", {
  fixtureDir <- tempfile("prot-summary-copy-error-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  workflowData <- list(state_manager = list(dummy = TRUE))
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      prepareProtSummaryCopyInputs = function(workflowData,
                                              projectDirs,
                                              omicType = "proteomics",
                                              readTsvFn = readr::read_tsv,
                                              catFn = cat) {
        list(
          contrastsTbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
          designMatrix = data.frame(sample = "S1", stringsAsFactors = FALSE)
        )
      },
      runProtSummaryPublicationCopy = function(output,
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
        stop("copy boom")
      },
      handleProtSummaryPublicationCopyError = function(output,
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
        recorded$errorArgs <- list(
          output = output,
          errorMessage = conditionMessage(error)
        )
        invisible(FALSE)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedOutput <- output
      values$workflow_args_saved <- TRUE
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(copy_to_publication = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$errorArgs$output, recorded$expectedOutput)
  expect_identical(recorded$errorArgs$errorMessage, "copy boom")
})

test_that("mod_prot_summary_server export observer delegates export completion through the seam", {
  fixtureDir <- tempfile("prot-summary-export-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      completeProtSummarySessionStateExport = function(projectDirs,
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
        recorded$completionArgs <- list(
          projectDirs = projectDirs,
          omicType = omicType,
          experimentLabel = experimentLabel,
          description = description,
          workflowArgsSaved = workflowArgsSaved,
          filesCopied = filesCopied,
          reportGenerated = reportGenerated,
          reportPath = reportPath
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(project_dirs = projectDirs),
    {
      values$workflow_args_saved <- TRUE
      values$files_copied <- TRUE
      values$report_generated <- TRUE
      values$report_path <- file.path(fixtureDir, "summary-report.html")
      session$setInputs(experiment_label = "summary-demo", description = "demo run")
      session$flushReact()

      session$setInputs(export_session_state = 1)
      session$flushReact()
    }
  )

  expect_identical(recorded$completionArgs$projectDirs, projectDirs)
  expect_identical(recorded$completionArgs$omicType, "proteomics")
  expect_identical(recorded$completionArgs$experimentLabel, "summary-demo")
  expect_identical(recorded$completionArgs$description, "demo run")
  expect_true(recorded$completionArgs$workflowArgsSaved)
  expect_true(recorded$completionArgs$filesCopied)
  expect_true(recorded$completionArgs$reportGenerated)
  expect_identical(
    recorded$completionArgs$reportPath,
    file.path(fixtureDir, "summary-report.html")
  )
})

test_that("mod_prot_summary_server delegates GitHub observer registration through the seam", {
  fixtureDir <- tempfile("prot-summary-github-")
  sourceDir <- file.path(fixtureDir, "source")
  integrationDir <- file.path(fixtureDir, "integration")
  dir.create(sourceDir, recursive = TRUE)
  dir.create(integrationDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = sourceDir,
      integration_dir = integrationDir,
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      observeProtSummaryGithubPush = function(input,
                                              output,
                                              values,
                                              projectDirs,
                                              omicType = "proteomics",
                                              completeGithubPushFn = completeProtSummaryGithubPush,
                                              withProgressFn = shiny::withProgress) {
        recorded$args <- list(
          input = input,
          output = output,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          completeGithubPushFn = completeGithubPushFn,
          withProgressFn = withProgressFn
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(project_dirs = projectDirs),
    {
      recorded$expectedInput <- input
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$flushReact()
    }
  )

  expect_identical(recorded$args$input, recorded$expectedInput)
  expect_identical(recorded$args$output, recorded$expectedOutput)
  expect_identical(recorded$args$values, recorded$expectedValues)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(recorded$args$completeGithubPushFn, completeProtSummaryGithubPush)
  expect_identical(recorded$args$withProgressFn, shiny::withProgress)
})

test_that("mod_prot_summary_server delegates copy observer registration through the seam", {
  fixtureDir <- tempfile("prot-summary-copy-registration-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  workflowData <- list(state_manager = list(dummy = TRUE))
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      observeProtSummaryPublicationCopy = function(input,
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
        recorded$args <- list(
          input = input,
          output = output,
          values = values,
          projectDirs = projectDirs,
          workflowData = workflowData,
          omicType = omicType,
          fallbackBootstrapFn = fallbackBootstrapFn,
          prepareCopyInputsFn = prepareCopyInputsFn,
          runPublicationCopyFn = runPublicationCopyFn,
          handleCopyErrorFn = handleCopyErrorFn,
          withProgressFn = withProgressFn,
          catFn = catFn
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(
      project_dirs = projectDirs,
      workflow_data = workflowData
    ),
    {
      recorded$expectedInput <- input
      recorded$expectedOutput <- output
      recorded$expectedValues <- values
      session$flushReact()
    }
  )

  expect_identical(recorded$args$input, recorded$expectedInput)
  expect_identical(recorded$args$output, recorded$expectedOutput)
  expect_identical(recorded$args$values, recorded$expectedValues)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$workflowData, workflowData)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(recorded$args$fallbackBootstrapFn, bootstrapProtSummaryCopyFallbackStudyParams)
  expect_identical(recorded$args$prepareCopyInputsFn, prepareProtSummaryCopyInputs)
  expect_identical(recorded$args$runPublicationCopyFn, runProtSummaryPublicationCopy)
  expect_identical(recorded$args$handleCopyErrorFn, handleProtSummaryPublicationCopyError)
  expect_identical(recorded$args$withProgressFn, shiny::withProgress)
  expect_identical(recorded$args$catFn, base::cat)
})

test_that("mod_prot_summary_server initializes default outputs through the seam", {
  fixtureDir <- tempfile("prot-summary-init-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      initializeProtSummaryDefaultOutputs = function(output,
                                                     initialSessionSummary = "Ready to save workflow parameters and generate report",
                                                     renderTextFn = shiny::renderText,
                                                     reactiveFn = shiny::reactive,
                                                     outputOptionsFn = shiny::outputOptions) {
        recorded$output <- output
        recorded$initialSessionSummary <- initialSessionSummary
        output$session_summary <- renderTextFn({
          "Delegated default summary"
        })
        output$report_ready <- reactiveFn({ FALSE })
        outputOptionsFn(output, "report_ready", suspendWhenHidden = FALSE)
        invisible(TRUE)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(project_dirs = projectDirs),
    {
      session$flushReact()
    }
  )

  expect_true(!is.null(recorded$output))
  expect_identical(
    recorded$initialSessionSummary,
    "Ready to save workflow parameters and generate report"
  )
})

test_that("mod_prot_summary_server delegates export observer registration through the seam", {
  fixtureDir <- tempfile("prot-summary-export-registration-")
  dir.create(fixtureDir, recursive = TRUE)
  on.exit(unlink(fixtureDir, recursive = TRUE, force = TRUE), add = TRUE)

  projectDirs <- list(
    proteomics = list(
      source_dir = file.path(fixtureDir, "source"),
      integration_dir = file.path(fixtureDir, "integration"),
      base_dir = fixtureDir
    )
  )
  recorded <- new.env(parent = emptyenv())

  serverUnderTest <- makeFunctionWithOverrides(
    mod_prot_summary_server,
    list(
      observeProtSummarySessionStateExport = function(input,
                                                      values,
                                                      projectDirs,
                                                      omicType = "proteomics",
                                                      completeSessionStateExportFn = completeProtSummarySessionStateExport) {
        recorded$args <- list(
          input = input,
          values = values,
          projectDirs = projectDirs,
          omicType = omicType,
          completeSessionStateExportFn = completeSessionStateExportFn
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    serverUnderTest,
    args = list(project_dirs = projectDirs),
    {
      recorded$expectedInput <- input
      recorded$expectedValues <- values
      session$flushReact()
    }
  )

  expect_identical(recorded$args$input, recorded$expectedInput)
  expect_identical(recorded$args$values, recorded$expectedValues)
  expect_identical(recorded$args$projectDirs, projectDirs)
  expect_identical(recorded$args$omicType, "proteomics")
  expect_identical(
    recorded$args$completeSessionStateExportFn,
    completeProtSummarySessionStateExport
  )
})
