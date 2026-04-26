# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

multiScholaRNamespace <- function() {
  asNamespace("MultiScholaR")
}

hasMultiScholaRBinding <- function(name) {
  exists(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

getMultiScholaRBinding <- function(name) {
  get(name, envir = multiScholaRNamespace(), inherits = FALSE)
}

assignSelectedBindings <- function(symbols, env) {
  for (symbol in symbols) {
    if (hasMultiScholaRBinding(symbol)) {
      assign(symbol, getMultiScholaRBinding(symbol), envir = env)
    } else {
      assign(
        symbol,
        eval(bquote(function(...) {
          skip(.(paste("requires extracted helper binding:", symbol)))
        })),
        envir = env
      )
    }
  }
}

assignSelectedBindings(
  symbols = c(
    "buildLipidSummaryTemplateStatus",
    "registerLipidSummaryTemplateStatusOutput",
    "buildLipidSummarySessionState",
    "handleLipidSummaryExportSessionState",
    "registerLipidSummaryExportSessionObserver",
    "collectLipidSummaryWorkflowArgsContext",
    "handleLipidSummarySaveWorkflowArgs",
    "registerLipidSummarySaveWorkflowArgsObserver",
    "handleLipidSummaryCopyToPublication",
    "registerLipidSummaryCopyToPublicationObserver",
    "handleLipidSummaryGenerateReport",
    "registerLipidSummaryGenerateReportObserver",
    "handleLipidSummaryPushToGithub",
    "registerLipidSummaryPushToGithubObserver",
    "initializeLipidSummarySessionBootstrap",
    "registerLipidSummaryInitialOutputs",
    "mod_lipid_summary_server"
  ),
  env = environment()
)

if (!methods::isClass("LipidSummaryCharacterizationS4")) {
  methods::setClass("LipidSummaryCharacterizationS4", slots = c(args = "list"))
}

makeLipidSummaryPublicHarness <- function() {
  base_dir <- tempfile("lipid-summary-public-")
  source_dir <- file.path(base_dir, "source")
  integration_dir <- file.path(base_dir, "integration")
  dir.create(file.path(base_dir, "scripts", "lipidomics"), recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)
  dir.create(integration_dir, recursive = TRUE)

  list(
    project_dirs = list(
      lipidomics = list(
        base_dir = base_dir,
        source_dir = source_dir,
        integration_dir = integration_dir
      )
    ),
    workflow_data = list(
      design_matrix = data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE),
      contrasts_tbl = data.frame(contrast = "A-B", stringsAsFactors = FALSE)
    )
  )
}

test_that("mod_lipid_summary_server preserves public wrapper initialization behavior", {
  harness <- makeLipidSummaryPublicHarness()

  testServer(
    mod_lipid_summary_server,
    args = list(
      project_dirs = harness$project_dirs,
      omic_type = "lipidomics",
      experiment_label = "Lipid Session",
      workflow_data = harness$workflow_data
    ),
    {
      session$flushReact()

      expect_true(exists("output"))
      expect_true(exists("session"))
      expect_true(exists("input"))
      expect_s3_class(session, "session_proxy")
      expect_match(
        output$session_summary,
        "Ready to save workflow parameters and generate report",
        fixed = TRUE
      )
      expect_match(
        output$template_status,
        "downloaded when generating report",
        fixed = TRUE
      )
    }
  )
})

test_that("lipid summary template and output helpers preserve current text variants", {
  base_dir <- tempfile("lipid-summary-template-")
  template_dir <- file.path(base_dir, "scripts", "lipidomics")
  dir.create(template_dir, recursive = TRUE)

  project_dirs <- list(lipidomics = list(base_dir = base_dir))

  expect_identical(
    buildLipidSummaryTemplateStatus(
      projectDirs = list(other = list(base_dir = base_dir)),
      omicType = "lipidomics"
    ),
    "[WARNING] Project directories not available"
  )
  expect_identical(
    buildLipidSummaryTemplateStatus(
      projectDirs = project_dirs,
      omicType = "lipidomics"
    ),
    "[WARNING] Report template will be downloaded when generating report"
  )

  writeLines("report", file.path(template_dir, "lipidomics_report.rmd"))
  expect_identical(
    buildLipidSummaryTemplateStatus(
      projectDirs = project_dirs,
      omicType = "lipidomics"
    ),
    "Template: Lipidomics Report [OK]"
  )

  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())

  result <- registerLipidSummaryTemplateStatusOutput(
    output = output,
    projectDirs = project_dirs,
    omicType = "lipidomics",
    renderTextFn = function(expr) {
      captured$rendered <- force(expr)
      structure(list(kind = "render-text"), class = "render-text")
    },
    reqFn = function(value) {
      captured$req <- value
      invisible(value)
    },
    buildStatusFn = function(projectDirsArg, omicTypeArg) {
      captured$project_dirs <- projectDirsArg
      captured$omic_type <- omicTypeArg
      "template-status"
    }
  )

  expect_identical(result, output)
  expect_identical(captured$req, project_dirs)
  expect_identical(captured$project_dirs, project_dirs)
  expect_identical(captured$omic_type, "lipidomics")
  expect_identical(captured$rendered, "template-status")
  expect_s3_class(output$template_status, "render-text")
})

test_that("lipid summary bootstrap helpers preserve output and state initialization", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$options <- list()

  result <- registerLipidSummaryInitialOutputs(
    output = output,
    renderTextFn = function(expr) force(expr),
    reactiveFn = function(expr) force(expr),
    outputOptionsFn = function(outputArg, name, suspendWhenHidden) {
      captured$options[[length(captured$options) + 1L]] <- list(
        output = outputArg,
        name = name,
        suspendWhenHidden = suspendWhenHidden
      )
      invisible(NULL)
    }
  )

  expect_identical(result, output)
  expect_identical(output$session_summary, "Ready to save workflow parameters and generate report")
  expect_false(output$report_ready)
  expect_identical(captured$options[[1L]]$output, output)
  expect_identical(captured$options[[1L]]$name, "report_ready")
  expect_false(captured$options[[1L]]$suspendWhenHidden)

  captured$update <- NULL
  captured$reactive_values <- NULL
  reactive_values <- list(kind = "reactive-values")

  bootstrap <- initializeLipidSummarySessionBootstrap(
    session = "session-token",
    experimentLabel = "Lipid Session",
    updateTextInputFn = function(session, inputId, value) {
      captured$update <- list(session = session, inputId = inputId, value = value)
      invisible(NULL)
    },
    reactiveValuesFn = function(...) {
      captured$reactive_values <- list(...)
      reactive_values
    }
  )

  expect_identical(bootstrap, reactive_values)
  expect_identical(
    captured$update,
    list(session = "session-token", inputId = "experiment_label", value = "Lipid Session")
  )
  expect_identical(
    captured$reactive_values,
    list(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
  )

  captured$update <- NULL
  initializeLipidSummarySessionBootstrap(
    session = "session-token",
    experimentLabel = NULL,
    updateTextInputFn = function(...) {
      captured$update <- list(...)
      invisible(NULL)
    },
    reactiveValuesFn = function(...) list(...)
  )
  expect_null(captured$update)
})

test_that("lipid summary workflow context preserves state, contrasts, and config selection", {
  final_s4 <- methods::new(
    "LipidSummaryCharacterizationS4",
    args = list(normalization = list(method = "median"))
  )
  assigned_env <- new.env(parent = emptyenv())
  messages <- character()
  workflow_data <- list(
    state_manager = list(
      getHistory = function() c("lipid_qc_complete", "lipid_normalized"),
      getState = function(state_name) {
        expect_identical(state_name, "lipid_normalized")
        final_s4
      }
    ),
    contrasts_tbl = data.frame(contrast = "A-B", stringsAsFactors = FALSE),
    config_list = list(alpha = 0.05, method = "median")
  )

  context <- collectLipidSummaryWorkflowArgsContext(
    workflowData = workflow_data,
    assignFn = function(x, value, envir) {
      assign(x, value, envir = envir)
      invisible(NULL)
    },
    globalEnv = assigned_env,
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    }
  )

  expect_identical(context$finalS4Object, final_s4)
  expect_identical(context$contrastsTbl, workflow_data$contrasts_tbl)
  expect_identical(get("config_list", envir = assigned_env), workflow_data$config_list)
  expect_true(any(grepl("Retrieved DATA S4 object from state 'lipid_normalized'", messages, fixed = TRUE)))
  expect_true(any(grepl("Using contrasts_tbl from workflow_data", messages, fixed = TRUE)))

  messages <- character()
  empty_context <- collectLipidSummaryWorkflowArgsContext(
    workflowData = list(state_manager = NULL),
    catFn = function(...) {
      messages <<- c(messages, paste0(..., collapse = ""))
      invisible(NULL)
    }
  )
  expect_null(empty_context$finalS4Object)
  expect_null(empty_context$contrastsTbl)
  expect_true(any(grepl("No state manager available", messages, fixed = TRUE)))
})

test_that("lipid summary save-workflow helper preserves success and fallback side effects", {
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  output <- new.env(parent = emptyenv())
  notifications <- list()
  captured <- new.env(parent = emptyenv())
  captured$dirs <- character()

  result <- handleLipidSummarySaveWorkflowArgs(
    input = list(experiment_label = "lipid-study", description = "workflow summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(source_dir = "/tmp/lipid-source", base_dir = "/tmp/lipid-base")),
    omicType = "lipidomics",
    workflowData = list(marker = "workflow-data"),
    collectContextFn = function(workflowData, catFn) {
      expect_identical(workflowData$marker, "workflow-data")
      list(finalS4Object = list(kind = "final-s4"), contrastsTbl = list(kind = "contrasts"))
    },
    createWorkflowArgsFn = function(...) {
      captured$create_args <- list(...)
      "/tmp/lipid-source/study_parameters.txt"
    },
    saveRDSFn = function(object, file) {
      captured$save <- list(object = object, file = file)
      invisible(file)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      notifications[[length(notifications) + 1L]] <<- list(message = message, type = type, duration = duration)
      invisible(NULL)
    },
    renderTextFn = function(expr) force(expr),
    dirExistsFn = function(path) FALSE,
    dirCreateFn = function(path, recursive, showWarnings) {
      captured$dirs <- c(captured$dirs, path)
      invisible(TRUE)
    },
    timeFn = function() as.POSIXct("2026-04-16 13:00:00", tz = "UTC"),
    catFn = function(...) invisible(NULL)
  )

  expect_true(values$workflow_args_saved)
  expect_identical(result$studyParamsFile, "/tmp/lipid-source/study_parameters.txt")
  expect_identical(result$s4Filename, "lipidomics_lipid-study_final_s4.RDS")
  expect_identical(result$s4Filepath, "/tmp/lipid-base/integration/lipidomics_lipid-study_final_s4.RDS")
  expect_identical(captured$dirs, "/tmp/lipid-base/integration")
  expect_identical(captured$create_args$workflow_name, "lipid-study")
  expect_identical(captured$create_args$source_dir_path, "/tmp/lipid-source")
  expect_identical(captured$save$object, list(kind = "final-s4"))
  expect_match(output$session_summary, "Parameters saved [OK]", fixed = TRUE)
  expect_identical(
    notifications,
    list(
      list(message = "Saved Integration S4 Object", type = "message", duration = NULL),
      list(message = "Study parameters saved successfully", type = "message", duration = NULL)
    )
  )

  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  output <- new.env(parent = emptyenv())
  fallback <- new.env(parent = emptyenv())
  fallback$notifications <- list()

  fallback_result <- handleLipidSummarySaveWorkflowArgs(
    input = list(experiment_label = "lipid-study", description = "fallback summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(source_dir = "/tmp/lipid-source", base_dir = "/tmp/lipid-base")),
    omicType = "lipidomics",
    collectContextFn = function(...) stop("context unavailable"),
    fileExistsFn = function(path) FALSE,
    writeLinesFn = function(text, con) {
      fallback$written <- list(text = text, con = con)
      invisible(NULL)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      fallback$notifications[[length(fallback$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    timeFn = function() as.POSIXct("2026-04-16 13:30:00", tz = "UTC"),
    catFn = function(...) invisible(NULL)
  )

  expect_true(values$workflow_args_saved)
  expect_identical(fallback_result, "/tmp/lipid-source/study_parameters.txt")
  expect_identical(fallback$written$con, fallback_result)
  expect_true(any(grepl("Error: context unavailable", fallback$written$text, fixed = TRUE)))
  expect_identical(
    fallback$notifications,
    list(list(message = "Study parameters saved with warnings", type = "warning", duration = NULL))
  )
})

test_that("lipid summary copy-publication helper preserves file fallback, success, and error behavior", {
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  values$files_copied <- FALSE
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$read_paths <- character()
  captured$notifications <- list()
  assigned_env <- new.env(parent = emptyenv())

  source_dir <- file.path(tempdir(), "lipid-summary-copy-source")
  design_path <- file.path(source_dir, "design_matrix.tab")
  contrasts_path <- file.path(source_dir, "contrasts_tbl.tab")
  basic_path <- file.path(source_dir, "study_parameters.txt")

  result <- handleLipidSummaryCopyToPublication(
    input = list(experiment_label = "lipid-study", description = "copy summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(source_dir = source_dir)),
    omicType = "lipidomics",
    workflowData = list(),
    withProgressFn = function(message, expr) {
      captured$progress <- message
      force(expr)
    },
    fileExistsFn = function(path) path %in% c(design_path, contrasts_path),
    writeLinesFn = function(text, con) {
      captured$written <- list(text = text, con = con)
      invisible(NULL)
    },
    timeFn = function() as.POSIXct("2026-04-16 16:00:00", tz = "UTC"),
    readTsvFn = function(path, show_col_types = FALSE) {
      captured$read_paths <- c(captured$read_paths, path)
      if (identical(path, design_path)) {
        return(data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE))
      }
      if (identical(path, contrasts_path)) {
        return(data.frame(contrast = "A-B", stringsAsFactors = FALSE))
      }
      stop("unexpected path")
    },
    existsFn = function(name, envir) FALSE,
    assignFn = function(x, value, envir) {
      assign(x, value, envir = envir)
      invisible(NULL)
    },
    globalEnv = assigned_env,
    copyResultsSummaryFn = function(...) {
      captured$copy_args <- list(...)
      invisible(TRUE)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    renderTextFn = function(expr) force(expr),
    logErrorFn = function(message) stop(message),
    skipFormatterFn = identity,
    catFn = function(...) invisible(NULL),
    tracebackFn = function() invisible(NULL)
  )

  expect_identical(captured$progress, "Copying files to publication directory...")
  expect_true(values$workflow_args_saved)
  expect_true(values$files_copied)
  expect_identical(captured$written$con, basic_path)
  expect_identical(captured$read_paths, c(design_path, contrasts_path))
  expect_identical(get("project_dirs", envir = assigned_env), list(lipidomics = list(source_dir = source_dir)))
  expect_identical(captured$copy_args$omic_type, "lipidomics")
  expect_identical(captured$copy_args$experiment_label, "lipid-study")
  expect_true(captured$copy_args$force)
  expect_identical(result$designMatrix, data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE))
  expect_identical(result$contrastsTbl, data.frame(contrast = "A-B", stringsAsFactors = FALSE))
  expect_identical(output$copy_status, "Files copied to publication directory successfully [OK]")
  expect_match(output$session_summary, "Files copied [OK]", fixed = TRUE)
  expect_identical(
    captured$notifications,
    list(list(message = "Publication files copied", type = "message", duration = NULL))
  )

  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- TRUE
  values$files_copied <- FALSE
  output <- new.env(parent = emptyenv())
  error_capture <- new.env(parent = emptyenv())
  error_capture$notifications <- list()

  error_result <- handleLipidSummaryCopyToPublication(
    input = list(experiment_label = "lipid-study", description = "copy summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(source_dir = source_dir)),
    omicType = "lipidomics",
    workflowData = list(design_matrix = data.frame(sample = "S1"), contrasts_tbl = data.frame(contrast = "A-B")),
    withProgressFn = function(message, expr) force(expr),
    copyResultsSummaryFn = function(...) stop("copy failed"),
    showNotificationFn = function(message, type, duration = NULL) {
      error_capture$notifications[[length(error_capture$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    renderTextFn = function(expr) force(expr),
    logErrorFn = function(message) {
      error_capture$error <- message
      invisible(NULL)
    },
    skipFormatterFn = function(message) paste0("skip:", message),
    catFn = function(...) invisible(NULL),
    tracebackFn = function() {
      error_capture$traceback <- TRUE
      invisible(NULL)
    }
  )

  expect_null(error_result)
  expect_identical(output$copy_status, "Error: copy failed")
  expect_identical(
    error_capture$notifications,
    list(list(message = "Copy error: copy failed", type = "error", duration = 10))
  )
  expect_identical(error_capture$error, "skip:Failed to copy files: copy failed")
  expect_true(error_capture$traceback)
})

test_that("lipid summary report helper preserves template, render, and guard behavior", {
  values <- new.env(parent = emptyenv())
  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$progress <- list()
  captured$options <- list()
  captured$info <- character()

  base_dir <- tempfile("lipid-summary-report-")
  template_path <- file.path(base_dir, "scripts", "lipidomics", "lipidomics_report.rmd")
  dir.create(dirname(template_path), recursive = TRUE)
  writeLines("report", template_path)
  rendered_path <- file.path(base_dir, "reports", "lipid-report.html")
  dir.create(dirname(rendered_path), recursive = TRUE)
  writeLines("<html></html>", rendered_path)

  handleLipidSummaryGenerateReport(
    input = list(experiment_label = "lipid-study", description = "report summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(base_dir = base_dir)),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) {
      captured$progress_message <- message
      force(expr)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    incProgressFn = function(value, detail = NULL) {
      captured$progress[[length(captured$progress) + 1L]] <- list(value = value, detail = detail)
      invisible(NULL)
    },
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
      invisible(NULL)
    },
    logErrorFn = function(message) stop(message),
    renderReportExistsFn = function() TRUE,
    renderReportFn = function(...) {
      captured$render_args <- list(...)
      rendered_path
    },
    reactiveFn = function(expr) force(expr),
    outputOptionsFn = function(outputArg, name, suspendWhenHidden) {
      captured$options[[length(captured$options) + 1L]] <- list(
        output = outputArg,
        name = name,
        suspendWhenHidden = suspendWhenHidden
      )
      invisible(NULL)
    },
    downloadHandlerFn = function(filename, content) {
      captured$download <- list(filename = filename, content = content)
      "download-handler"
    },
    renderTextFn = function(expr) force(expr),
    timeFn = function() as.POSIXct("2026-04-16 14:30:00", tz = "UTC"),
    catFn = function(...) invisible(NULL),
    printFn = function(...) invisible(NULL)
  )

  expect_true(values$report_generated)
  expect_identical(values$report_path, rendered_path)
  expect_identical(captured$progress_message, "Generating report...")
  expect_identical(captured$render_args$rmd_filename, "lipidomics_report.rmd")
  expect_identical(output$report_ready, TRUE)
  expect_identical(output$download_report, "download-handler")
  expect_identical(captured$download$filename(), "lipid-report.html")
  expect_match(output$session_summary, "Report generated [OK]", fixed = TRUE)
  expect_identical(captured$options[[1L]]$name, "report_ready")
  expect_false(captured$options[[1L]]$suspendWhenHidden)
  expect_identical(
    captured$notifications,
    list(list(message = "Report generated successfully!", type = "message", duration = NULL))
  )

  guard_notifications <- list()
  guard_result <- withVisible(handleLipidSummaryGenerateReport(
    input = list(experiment_label = "lipid-study"),
    values = values,
    output = output,
    projectDirs = list(),
    omicType = "lipidomics",
    showNotificationFn = function(message, type, duration = NULL) {
      guard_notifications[[length(guard_notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    }
  ))
  expect_false(guard_result$visible)
  expect_null(guard_result$value)
  expect_identical(
    guard_notifications,
    list(list(message = "Error: Project directories not properly initialized", type = "error", duration = 10))
  )
})

test_that("lipid summary report helper preserves template retrieval and render failure branches", {
  values <- new.env(parent = emptyenv())
  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()

  base_dir <- tempfile("lipid-summary-report-template-")
  package_template <- tempfile("lipid-package-template-", fileext = ".rmd")
  writeLines("report", package_template)
  rendered_path <- file.path(base_dir, "reports", "missing-output.html")

  handleLipidSummaryGenerateReport(
    input = list(experiment_label = "lipid-study", description = "report summary"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(base_dir = base_dir)),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) force(expr),
    showNotificationFn = function(message, type, duration = NULL) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    incProgressFn = function(...) invisible(NULL),
    fileExistsFn = function(path) identical(path, package_template),
    dirCreateFn = function(path, recursive = TRUE, showWarnings = FALSE) {
      captured$dir <- list(path = path, recursive = recursive, showWarnings = showWarnings)
      invisible(TRUE)
    },
    systemFileFn = function(...) package_template,
    fileCopyFn = function(from, to, ...) {
      captured$copy <- list(from = from, to = to)
      TRUE
    },
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
      invisible(NULL)
    },
    logErrorFn = function(message) {
      captured$error <- c(captured$error, message)
      invisible(NULL)
    },
    renderReportExistsFn = function() TRUE,
    renderReportFn = function(...) rendered_path,
    catFn = function(...) invisible(NULL),
    printFn = function(...) invisible(NULL)
  )

  expect_identical(captured$copy$from, package_template)
  expect_match(captured$copy$to, "lipidomics_report.rmd", fixed = TRUE)
  expect_false(values$report_generated)
  expect_identical(
    captured$notifications,
    list(
      list(message = "lipidomics_report.rmd template copied from package", type = "message", duration = NULL),
      list(message = "Report generation failed - no output file created", type = "error", duration = 10)
    )
  )

  missing_render_notifications <- list()
  handleLipidSummaryGenerateReport(
    input = list(experiment_label = "lipid-study"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(base_dir = base_dir)),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) force(expr),
    showNotificationFn = function(message, type, duration = NULL) {
      missing_render_notifications[[length(missing_render_notifications) + 1L]] <<- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    incProgressFn = function(...) invisible(NULL),
    fileExistsFn = function(path) TRUE,
    renderReportExistsFn = function() FALSE,
    catFn = function(...) invisible(NULL),
    printFn = function(...) invisible(NULL)
  )
  expect_identical(
    missing_render_notifications,
    list(list(
      message = "Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.",
      type = "error",
      duration = 15
    ))
  )

  render_error <- new.env(parent = emptyenv())
  render_error$notifications <- list()
  handleLipidSummaryGenerateReport(
    input = list(experiment_label = "lipid-study"),
    values = values,
    output = output,
    projectDirs = list(lipidomics = list(base_dir = base_dir)),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) force(expr),
    showNotificationFn = function(message, type, duration = NULL) {
      render_error$notifications[[length(render_error$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    incProgressFn = function(...) invisible(NULL),
    fileExistsFn = function(path) TRUE,
    renderReportExistsFn = function() TRUE,
    renderReportFn = function(...) stop("render exploded"),
    logInfoFn = function(...) invisible(NULL),
    logErrorFn = function(message) {
      render_error$log <- c(render_error$log, message)
      invisible(NULL)
    },
    catFn = function(...) invisible(NULL),
    printFn = function(...) invisible(NULL)
  )

  expect_identical(render_error$notifications[[1L]]$message, "Report generation failed: render exploded")
  expect_identical(render_error$notifications[[1L]]$type, "error")
  expect_identical(render_error$notifications[[2L]]$type, "warning")
  expect_true(any(grepl("Failed to generate report", render_error$log, fixed = TRUE)))
})

test_that("lipid summary GitHub helper preserves success and error side effects", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()

  result <- handleLipidSummaryPushToGithub(
    input = list(
      enable_github = TRUE,
      github_org = "apaf-bioinformatics",
      github_email = "team@example.org",
      github_username = "lipid-user",
      project_id = "proj-001",
      experiment_label = "lipid-study",
      description = "github push"
    ),
    values = list(report_generated = TRUE),
    output = output,
    projectDirs = list(lipidomics = list(base_dir = tempdir())),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) {
      captured$progress <- message
      force(expr)
    },
    optionsFn = function(...) {
      captured$options <- list(...)
      invisible(NULL)
    },
    pushProjectFn = function(...) {
      captured$push <- list(...)
      invisible(TRUE)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    logErrorFn = function(message) stop(message),
    skipFormatterFn = identity,
    renderTextFn = function(expr) force(expr),
    timeFn = function() as.POSIXct("2026-04-16 15:00:00", tz = "UTC")
  )

  expect_true(result)
  expect_identical(captured$progress, "Pushing to GitHub...")
  expect_identical(captured$options$github_org, "apaf-bioinformatics")
  expect_identical(captured$push$project_id, "proj-001")
  expect_match(output$session_summary, "GitHub pushed [OK]", fixed = TRUE)
  expect_identical(
    captured$notifications,
    list(list(message = "Successfully pushed to GitHub", type = "message", duration = NULL))
  )

  error_capture <- new.env(parent = emptyenv())
  error_capture$notifications <- list()
  error_result <- handleLipidSummaryPushToGithub(
    input = list(
      github_org = "apaf-bioinformatics",
      github_email = "team@example.org",
      github_username = "lipid-user",
      project_id = "proj-001",
      experiment_label = "lipid-study",
      description = "github push"
    ),
    values = list(report_generated = TRUE),
    output = output,
    projectDirs = list(lipidomics = list(base_dir = tempdir())),
    omicType = "lipidomics",
    withProgressFn = function(message, expr) force(expr),
    pushProjectFn = function(...) stop("push failed"),
    showNotificationFn = function(message, type, duration = NULL) {
      error_capture$notifications[[length(error_capture$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    logErrorFn = function(message) {
      error_capture$error <- message
      invisible(NULL)
    },
    skipFormatterFn = function(message) paste0("skip:", message)
  )

  expect_null(error_result)
  expect_identical(
    error_capture$notifications,
    list(list(message = "GitHub push failed: push failed", type = "error", duration = 10))
  )
  expect_identical(error_capture$error, "skip:Failed to push to GitHub: push failed")
})

test_that("lipid summary export-session helpers preserve payload, success, and error behavior", {
  fixed_time <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")
  input <- list(experiment_label = "lipid-study", description = "session export")
  values <- list(
    workflow_args_saved = TRUE,
    files_copied = FALSE,
    report_generated = TRUE,
    report_path = "/tmp/report.html"
  )
  project_dirs <- list(lipidomics = list(source_dir = "/tmp/lipid-source"))

  state <- buildLipidSummarySessionState(
    input = input,
    values = values,
    projectDirs = project_dirs,
    omicType = "lipidomics",
    timeFn = function() fixed_time
  )
  expect_identical(
    state,
    list(
      experiment_label = "lipid-study",
      description = "session export",
      timestamp = fixed_time,
      omic_type = "lipidomics",
      workflow_args_saved = TRUE,
      files_copied = FALSE,
      report_generated = TRUE,
      report_path = "/tmp/report.html",
      project_dirs = project_dirs
    )
  )

  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  result <- handleLipidSummaryExportSessionState(
    input = input,
    values = values,
    projectDirs = project_dirs,
    omicType = "lipidomics",
    saveRDSFn = function(object, file) {
      captured$save <- list(object = object, file = file)
      invisible(file)
    },
    showNotificationFn = function(message, type, duration = NULL) {
      captured$notifications[[length(captured$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    logInfoFn = function(message) {
      captured$info <- message
      invisible(NULL)
    },
    logErrorFn = function(message) stop(message),
    dateFn = function() as.Date("2026-04-16"),
    buildSessionStateFn = function(...) list(kind = "session-state")
  )

  expect_identical(result, "/tmp/lipid-source/session_state_2026-04-16.RDS")
  expect_identical(captured$save$object, list(kind = "session-state"))
  expect_identical(captured$save$file, result)
  expect_identical(captured$info, paste("Session state exported to:", result))
  expect_identical(
    captured$notifications,
    list(list(message = paste("Session state exported to:", result), type = "message", duration = NULL))
  )

  error_capture <- new.env(parent = emptyenv())
  error_capture$notifications <- list()
  error_result <- handleLipidSummaryExportSessionState(
    input = input,
    values = values,
    projectDirs = project_dirs,
    omicType = "lipidomics",
    saveRDSFn = function(...) stop("disk full"),
    showNotificationFn = function(message, type, duration = NULL) {
      error_capture$notifications[[length(error_capture$notifications) + 1L]] <- list(
        message = message,
        type = type,
        duration = duration
      )
      invisible(NULL)
    },
    logInfoFn = function(message) stop(message),
    logErrorFn = function(message) {
      error_capture$error <- message
      invisible(NULL)
    },
    skipFormatterFn = function(message) paste0("skip:", message),
    dateFn = function() as.Date("2026-04-16"),
    buildSessionStateFn = function(...) list(kind = "session-state")
  )

  expect_null(error_result)
  expect_identical(
    error_capture$notifications,
    list(list(message = "Export failed: disk full", type = "error", duration = 10))
  )
  expect_identical(error_capture$error, "skip:Failed to export session state: disk full")
})

test_that("lipid summary observer registration helpers preserve event handoffs", {
  input <- list(
    save_workflow_args = 11L,
    copy_to_publication = 12L,
    generate_report = 13L,
    push_to_github = 14L,
    export_session_state = 15L,
    experiment_label = "lipid-study",
    enable_github = TRUE,
    github_org = "apaf-bioinformatics",
    github_email = "team@example.org",
    github_username = "lipid-user",
    project_id = "proj-001"
  )
  values <- list(
    workflow_args_saved = TRUE,
    files_copied = TRUE,
    report_generated = TRUE
  )
  output <- new.env(parent = emptyenv())
  project_dirs <- list(lipidomics = list(source_dir = tempdir(), base_dir = tempdir()))
  workflow_data <- list(marker = "workflow-data")
  events <- list()
  req_values <- list()
  handlers <- list()

  observeEventFn <- function(event, handler, ignoreInit = FALSE) {
    events[[length(events) + 1L]] <<- list(event = event, ignoreInit = ignoreInit)
    force(handler)
    invisible("registered")
  }
  reqFn <- function(...) {
    args <- list(...)
    req_values[[length(req_values) + 1L]] <<- args
    if (length(args) == 1L) {
      args[[1L]]
    } else {
      args
    }
  }

  expect_identical(
    registerLipidSummarySaveWorkflowArgsObserver(
      input = input,
      values = values,
      output = output,
      projectDirs = project_dirs,
      omicType = "lipidomics",
      workflowData = workflow_data,
      observeEventFn = observeEventFn,
      reqFn = reqFn,
      handleSaveFn = function(...) {
        handlers$save <<- list(...)
        invisible(TRUE)
      }
    ),
    input
  )
  expect_identical(
    registerLipidSummaryCopyToPublicationObserver(
      input = input,
      values = values,
      output = output,
      projectDirs = project_dirs,
      omicType = "lipidomics",
      workflowData = workflow_data,
      observeEventFn = observeEventFn,
      reqFn = reqFn,
      handleCopyFn = function(...) {
        handlers$copy <<- list(...)
        invisible(TRUE)
      }
    ),
    input
  )
  expect_identical(
    registerLipidSummaryGenerateReportObserver(
      input = input,
      values = values,
      output = output,
      projectDirs = project_dirs,
      omicType = "lipidomics",
      observeEventFn = observeEventFn,
      reqFn = reqFn,
      handleGenerateReportFn = function(...) {
        handlers$generate <<- list(...)
        invisible(TRUE)
      }
    ),
    input
  )
  expect_identical(
    registerLipidSummaryPushToGithubObserver(
      input = input,
      values = values,
      output = output,
      projectDirs = project_dirs,
      omicType = "lipidomics",
      observeEventFn = observeEventFn,
      reqFn = reqFn,
      handlePushFn = function(...) {
        handlers$push <<- list(...)
        invisible(TRUE)
      }
    ),
    input
  )
  expect_identical(
    registerLipidSummaryExportSessionObserver(
      input = input,
      values = values,
      projectDirs = project_dirs,
      omicType = "lipidomics",
      observeEventFn = observeEventFn,
      reqFn = reqFn,
      handleExportFn = function(...) {
        handlers$export <<- list(...)
        invisible(TRUE)
      }
    ),
    input
  )

  expect_identical(vapply(events, `[[`, numeric(1), "event"), c(11, 12, 13, 14, 15))
  expect_identical(vapply(events, `[[`, logical(1), "ignoreInit"), c(FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(req_values[[1L]], list("lipid-study"))
  expect_identical(req_values[[2L]], list("lipid-study"))
  expect_identical(req_values[[3L]], list("lipid-study"))
  expect_identical(req_values[[4L]], list(TRUE))
  expect_identical(req_values[[5L]], list(TRUE, "apaf-bioinformatics", "team@example.org", "lipid-user", "proj-001"))
  expect_identical(req_values[[6L]], list(TRUE))
  expect_identical(req_values[[7L]], list("lipid-study"))
  expect_identical(handlers$save$workflowData, workflow_data)
  expect_identical(handlers$copy$workflowData, workflow_data)
  expect_identical(handlers$generate$omicType, "lipidomics")
  expect_identical(handlers$push$projectDirs, project_dirs)
  expect_identical(handlers$export$values, values)
})
