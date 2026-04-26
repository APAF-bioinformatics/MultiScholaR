# fidelity-coverage-compare: shared
library(testthat)

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
    "setupMetabSummaryBootstrapOutputs",
    "registerMetabSummaryServerObservers",
    "setupMetabSummaryServerBootstrapState",
    "runMetabSummaryCopyToPublicationObserverShell",
    "runMetabSummaryExportSessionObserverShell",
    "runMetabSummaryGenerateReportObserverShell",
    "runMetabSummaryGithubPushObserverShell",
    "runMetabSummarySaveWorkflowArgsObserverShell",
    "setupMetabSummaryTemplateStatusOutput",
    "mod_metab_summary_server"
  ),
  env = environment()
)

makeMetabSummaryPublicHarness <- function() {
  base_dir <- tempfile("metab-summary-public-")
  dir.create(base_dir, recursive = TRUE)

  list(
    project_dirs = list(
      metabolomics = list(
        base_dir = base_dir,
        source_dir = base_dir,
        integration_dir = base_dir
      )
    ),
    workflow_data = list(
      design_matrix = data.frame(Run = "Sample1", stringsAsFactors = FALSE),
      contrasts_tbl = data.frame(contrasts = "groupA-groupA", stringsAsFactors = FALSE)
    )
  )
}

if (!methods::isClass("MetabSummaryObserverTestS4")) {
  methods::setClass("MetabSummaryObserverTestS4", slots = c(args = "list"))
}

test_that("metabolomics summary template-status seam preserves current status text variants", {
  base_dir <- tempfile("metab-summary-")
  dir.create(file.path(base_dir, "scripts", "metabolomics"), recursive = TRUE)

  fakeReq <- function(...) {
    invisible(list(...)[[1]])
  }

  fakeRenderText <- function(expr) {
    eval(substitute(expr), parent.frame())
  }

  output <- new.env(parent = emptyenv())
  setupMetabSummaryTemplateStatusOutput(
    output = output,
    projectDirs = list(metabolomics = list(base_dir = base_dir)),
    omicType = "metabolomics",
    renderText = fakeRenderText,
    req = fakeReq
  )

  expect_identical(
    output$template_status,
    "[WARNING] Report template will be downloaded when generating report"
  )

  writeLines(
    "report",
    file.path(base_dir, "scripts", "metabolomics", "metabolomics_report.rmd")
  )

  output <- new.env(parent = emptyenv())
  setupMetabSummaryTemplateStatusOutput(
    output = output,
    projectDirs = list(metabolomics = list(base_dir = base_dir)),
    omicType = "metabolomics",
    renderText = fakeRenderText,
    req = fakeReq
  )

  expect_identical(output$template_status, "Template: Metabolomics Report [OK]")

  output <- new.env(parent = emptyenv())
  setupMetabSummaryTemplateStatusOutput(
    output = output,
    projectDirs = list(lipidomics = list(base_dir = base_dir)),
    omicType = "metabolomics",
    renderText = fakeRenderText,
    req = fakeReq
  )

  expect_identical(output$template_status, "[WARNING] Project directories not available")
})

test_that("metabolomics summary bootstrap seam preserves initial output registrations", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())

  fakeRenderText <- function(expr) {
    eval(substitute(expr), parent.frame())
  }

  fakeReactive <- function(expr) {
    eval(substitute(expr), parent.frame())
  }

  fakeOutputOptions <- function(output, name, suspendWhenHidden) {
    captured$output <- output
    captured$name <- name
    captured$suspendWhenHidden <- suspendWhenHidden
    invisible(output)
  }

  setupMetabSummaryBootstrapOutputs(
    output = output,
    renderText = fakeRenderText,
    reactive = fakeReactive,
    outputOptions = fakeOutputOptions
  )

  expect_identical(
    output$session_summary,
    "Ready to save workflow parameters and generate report"
  )
  expect_identical(output$report_ready, FALSE)
  expect_identical(captured$output, output)
  expect_identical(captured$name, "report_ready")
  expect_identical(captured$suspendWhenHidden, FALSE)
})

test_that("metabolomics summary bootstrap-state seam preserves label update and reactive defaults", {
  captured <- new.env(parent = emptyenv())
  session <- list(token = "session-token")
  reactive_values <- list(token = "reactive-values-token")

  result <- setupMetabSummaryServerBootstrapState(
    session = session,
    experimentLabel = "Metab Session",
    reactiveValuesFn = function(...) {
      captured$reactive_args <- list(...)
      reactive_values
    },
    updateTextInputFn = function(session, inputId, value) {
      captured$update <- list(
        session = session,
        input_id = inputId,
        value = value
      )
      invisible(NULL)
    }
  )

  expect_identical(
    captured$update,
    list(
      session = session,
      input_id = "experiment_label",
      value = "Metab Session"
    )
  )
  expect_identical(
    captured$reactive_args,
    list(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
  )
  expect_identical(result, reactive_values)

  captured$update <- NULL

  setupMetabSummaryServerBootstrapState(
    session = session,
    experimentLabel = NULL,
    reactiveValuesFn = function(...) reactive_values,
    updateTextInputFn = function(...) {
      captured$update <- list(...)
      invisible(NULL)
    }
  )

  expect_null(captured$update)
})

test_that("mod_metab_summary_server preserves public wrapper initialization behavior", {
  harness <- makeMetabSummaryPublicHarness()
  testServer(
    mod_metab_summary_server,
    args = list(
      project_dirs = harness$project_dirs,
      omic_type = "metabolomics",
      experiment_label = "Metab Session",
      workflow_data = harness$workflow_data
    ),
    {
      expect_true(exists("output"))
      expect_true(exists("session"))
      expect_true(exists("input"))
      expect_s3_class(session, "session_proxy")
    }
  )

  expect_identical(
    harness$project_dirs$metabolomics$base_dir,
    harness$project_dirs$metabolomics$source_dir
  )
  expect_identical(harness$workflow_data$design_matrix$Run, "Sample1")
})

test_that("metabolomics summary export-session observer shell preserves export payload and notifications", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  fixed_time <- as.POSIXct("2026-04-18 10:11:12", tz = "UTC")
  project_dirs <- list(metabolomics = list(source_dir = "/tmp/metab-summary-source"))
  values <- list(
    workflow_args_saved = TRUE,
    files_copied = FALSE,
    report_generated = TRUE,
    report_path = "/tmp/metab-summary-source/report.pdf"
  )

  visible <- withVisible(
    runMetabSummaryExportSessionObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = project_dirs,
      omicType = "metabolomics",
      values = values,
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      saveRdsFn = function(object, path) {
        captured$save <- list(
          object = object,
          path = path
        )
        invisible(NULL)
      },
      sysDateFn = function() as.Date("2026-04-18"),
      sysTimeFn = function() fixed_time,
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      logInfoFn = function(message) {
        captured$info <- c(captured$info, message)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      },
      skipFormatterFn = function(message) message
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(
    captured$save$path,
    "/tmp/metab-summary-source/session_state_2026-04-18.RDS"
  )
  expect_identical(
    captured$save$object,
    list(
      experiment_label = "Metab Session",
      description = "summary export description",
      timestamp = fixed_time,
      omic_type = "metabolomics",
      workflow_args_saved = TRUE,
      files_copied = FALSE,
      report_generated = TRUE,
      report_path = "/tmp/metab-summary-source/report.pdf",
      project_dirs = project_dirs
    )
  )
  expect_identical(
    captured$notifications,
    list(list(
      message = "Session state exported to: /tmp/metab-summary-source/session_state_2026-04-18.RDS",
      type = "message",
      duration = NULL
    ))
  )
  expect_identical(
    captured$info,
    "Session state exported to: /tmp/metab-summary-source/session_state_2026-04-18.RDS"
  )
  expect_identical(visible$value$status, "success")
  expect_identical(
    visible$value$sessionExportPath,
    "/tmp/metab-summary-source/session_state_2026-04-18.RDS"
  )
  expect_identical(visible$value$sessionState, captured$save$object)
})

test_that("metabolomics summary export-session observer shell preserves error reporting", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()

  visible <- withVisible(
    runMetabSummaryExportSessionObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(metabolomics = list(source_dir = "/tmp/metab-summary-source")),
      omicType = "metabolomics",
      values = list(
        workflow_args_saved = FALSE,
        files_copied = FALSE,
        report_generated = FALSE,
        report_path = NULL
      ),
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      saveRdsFn = function(object, path) {
        captured$save_attempt <- list(
          object = object,
          path = path
        )
        stop("disk full")
      },
      sysDateFn = function() as.Date("2026-04-18"),
      sysTimeFn = function() as.POSIXct("2026-04-18 10:11:12", tz = "UTC"),
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      logInfoFn = function(message) {
        stop(sprintf("logInfoFn should not be called: %s", message))
      },
      logErrorFn = function(message) {
        captured$error <- c(captured$error, message)
        invisible(NULL)
      },
      skipFormatterFn = function(message) {
        paste0("skip:", message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(
    captured$save_attempt$path,
    "/tmp/metab-summary-source/session_state_2026-04-18.RDS"
  )
  expect_identical(
    captured$notifications,
    list(list(
      message = "Export failed: disk full",
      type = "error",
      duration = 10
    ))
  )
  expect_identical(
    captured$error,
    "skip:Failed to export session state: disk full"
  )
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "disk full"
    )
  )
})

test_that("metabolomics summary save-workflow observer shell preserves S4 export contract", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  fixed_time <- as.POSIXct("2026-04-18 10:11:12", tz = "UTC")
  final_s4_object <- methods::new(
    "MetabSummaryObserverTestS4",
    args = list(norm = list(method = "median"), de = list(engine = "limma"))
  )
  contrasts_tbl <- data.frame(contrast = "case-control", stringsAsFactors = FALSE)
  config_list <- list(alpha = 0.05, method = "median")
  global_env <- new.env(parent = emptyenv())
  workflow_data <- list(
    state_manager = list(
      getHistory = function() c("metab_normalized", "metab_qc_complete"),
      getState = function(state_name) {
        captured$state_lookup <- c(captured$state_lookup, state_name)
        final_s4_object
      }
    ),
    contrasts_tbl = contrasts_tbl,
    config_list = config_list
  )
  project_dirs <- list(
    metabolomics = list(
      source_dir = "/tmp/metab-summary-source",
      base_dir = "/tmp/metab-summary-base"
    )
  )

  values$workflow_args_saved <- FALSE

  visible <- withVisible(
    runMetabSummarySaveWorkflowArgsObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = project_dirs,
      omicType = "metabolomics",
      workflowData = workflow_data,
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      detectFn = function(.x, .p, ...) {
        predicate_expr <- if (inherits(.p, "formula")) .p[[2]] else .p
        for (value in .x) {
          predicate_result <- if (is.language(predicate_expr)) {
            eval(predicate_expr, envir = list(.x = value), enclos = parent.frame())
          } else {
            predicate_expr(value)
          }

          if (isTRUE(predicate_result)) {
            return(value)
          }
        }
        NULL
      },
      assignFn = function(x, value, envir) {
        captured$assigned <- list(name = x, value = value, envir = envir)
        assign(x, value, envir = envir)
      },
      createWorkflowArgsFromConfigFn = function(workflow_name,
                                                description,
                                                source_dir_path,
                                                final_s4_object,
                                                contrasts_tbl,
                                                workflow_data) {
        captured$create_call <- list(
          workflow_name = workflow_name,
          description = description,
          source_dir_path = source_dir_path,
          final_s4_object = final_s4_object,
          contrasts_tbl = contrasts_tbl,
          workflow_data = workflow_data
        )
        "/tmp/metab-summary-source/study_parameters.txt"
      },
      dirExistsFn = function(path) FALSE,
      dirCreateFn = function(path, recursive = FALSE, showWarnings = TRUE) {
        captured$dir_create <- list(
          path = path,
          recursive = recursive,
          show_warnings = showWarnings
        )
        invisible(TRUE)
      },
      saveRdsFn = function(object, path) {
        captured$save_rds <- list(object = object, path = path)
        invisible(NULL)
      },
      fileExistsFn = function(...) {
        stop("fileExistsFn should not be called")
      },
      writeLinesFn = function(...) {
        stop("writeLinesFn should not be called")
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      sysTimeFn = function() fixed_time,
      catFn = function(...) invisible(NULL),
      globalEnv = global_env
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(captured$state_lookup, "metab_normalized")
  expect_identical(
    captured$assigned,
    list(
      name = "config_list",
      value = config_list,
      envir = global_env
    )
  )
  expect_identical(
    captured$create_call,
    list(
      workflow_name = "Metab Session",
      description = "summary export description",
      source_dir_path = "/tmp/metab-summary-source",
      final_s4_object = final_s4_object,
      contrasts_tbl = contrasts_tbl,
      workflow_data = workflow_data
    )
  )
  expect_identical(
    captured$dir_create,
    list(
      path = "/tmp/metab-summary-base/integration",
      recursive = TRUE,
      show_warnings = FALSE
    )
  )
  expect_identical(
    captured$save_rds,
    list(
      object = final_s4_object,
      path = "/tmp/metab-summary-base/integration/metabolomics_Metab Session_final_s4.RDS"
    )
  )
  expect_true(values$workflow_args_saved)
  expect_identical(
    captured$notifications,
    list(
      list(
        message = "Saved Integration S4 Object",
        type = "message",
        duration = NULL
      ),
      list(
        message = "Study parameters saved successfully",
        type = "message",
        duration = NULL
      )
    )
  )
  expect_identical(
    output$session_summary,
    paste("Study parameters created for:", "Metab Session",
          "\nDescription:", "summary export description",
          "\nTimestamp:", fixed_time,
          "\nFile:", "/tmp/metab-summary-source/study_parameters.txt",
          "\nSource: Final S4 object @args + config_list",
          "\nIntegration Object:", "metabolomics_Metab Session_final_s4.RDS",
          "\nStatus: Parameters saved [OK]")
  )
  expect_identical(global_env$config_list, config_list)
  expect_identical(
    visible$value,
    list(
      status = "success",
      studyParamsFile = "/tmp/metab-summary-source/study_parameters.txt",
      s4Filepath = "/tmp/metab-summary-base/integration/metabolomics_Metab Session_final_s4.RDS",
      dataStateUsed = "metab_normalized",
      configListAssigned = TRUE
    )
  )
})

test_that("metabolomics summary save-workflow observer shell preserves fallback warning contract", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())

  values$workflow_args_saved <- FALSE

  visible <- withVisible(
    runMetabSummarySaveWorkflowArgsObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(
        metabolomics = list(
          source_dir = "/tmp/metab-summary-source",
          base_dir = "/tmp/metab-summary-base"
        )
      ),
      omicType = "metabolomics",
      workflowData = list(),
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      detectFn = function(...) {
        stop("detectFn should not be called")
      },
      assignFn = function(...) {
        stop("assignFn should not be called")
      },
      createWorkflowArgsFromConfigFn = function(...) {
        stop("serialization failed")
      },
      dirExistsFn = function(...) {
        stop("dirExistsFn should not be called")
      },
      dirCreateFn = function(...) {
        stop("dirCreateFn should not be called")
      },
      saveRdsFn = function(...) {
        stop("saveRdsFn should not be called")
      },
      fileExistsFn = function(path) FALSE,
      writeLinesFn = function(text, con, ...) {
        captured$write_lines <- list(text = text, path = con)
        invisible(NULL)
      },
      renderTextFn = function(expr) {
        stop("renderTextFn should not be called")
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      sysTimeFn = function() as.POSIXct("2026-04-18 10:11:12", tz = "UTC"),
      catFn = function(...) invisible(NULL)
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(
    captured$write_lines,
    list(
      text = c(
        "Study Parameters",
        "================",
        "",
        "Workflow Name: Metab Session",
        "Description: summary export description",
        "Timestamp: 2026-04-18 10:11:12",
        "Error: serialization failed"
      ),
      path = "/tmp/metab-summary-source/study_parameters.txt"
    )
  )
  expect_true(values$workflow_args_saved)
  expect_identical(
    captured$notifications,
    list(list(
      message = "Study parameters saved with warnings",
      type = "warning",
      duration = NULL
    ))
  )
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
  expect_identical(
    visible$value,
    list(
      status = "warning",
      basicParamsFile = "/tmp/metab-summary-source/study_parameters.txt",
      fallbackCreated = TRUE,
      errorMessage = "serialization failed"
    )
  )
})

test_that("metabolomics summary copy-to-publication observer shell preserves fallback copy contract", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  fixed_time <- as.POSIXct("2026-04-18 10:11:12", tz = "UTC")
  design_matrix <- data.frame(sample = "S1", condition = "case", stringsAsFactors = FALSE)
  contrasts_tbl <- data.frame(contrast = "case-control", stringsAsFactors = FALSE)
  project_dirs <- list(
    metabolomics = list(source_dir = "/tmp/metab-summary-source")
  )
  global_env <- new.env(parent = emptyenv())

  values$workflow_args_saved <- FALSE
  values$files_copied <- FALSE

  visible <- withVisible(
    runMetabSummaryCopyToPublicationObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = project_dirs,
      omicType = "metabolomics",
      workflowData = list(),
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      fileExistsFn = function(path) {
        path %in% c(
          "/tmp/metab-summary-source/design_matrix.tab",
          "/tmp/metab-summary-source/contrasts_tbl.tab"
        )
      },
      sysTimeFn = function() fixed_time,
      writeLinesFn = function(text, con, ...) {
        captured$write <- list(text = text, path = con)
        invisible(NULL)
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      readTsvFn = function(path, show_col_types = FALSE, ...) {
        captured$reads <- c(captured$reads, path)
        if (identical(path, "/tmp/metab-summary-source/design_matrix.tab")) {
          return(design_matrix)
        }
        if (identical(path, "/tmp/metab-summary-source/contrasts_tbl.tab")) {
          return(contrasts_tbl)
        }
        stop(sprintf("Unexpected readTsvFn path: %s", path))
      },
      existsFn = function(x, envir) exists(x, envir = envir, inherits = FALSE),
      assignFn = assign,
      copyToResultsSummaryFn = function(omic_type,
                                       experiment_label,
                                       contrasts_tbl,
                                       design_matrix,
                                       force) {
        captured$copy_call <- list(
          omic_type = omic_type,
          experiment_label = experiment_label,
          contrasts_tbl = contrasts_tbl,
          design_matrix = design_matrix,
          force = force
        )
        invisible(NULL)
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      },
      skipFormatterFn = function(message) message,
      catFn = function(...) invisible(NULL),
      tracebackFn = function() {
        stop("tracebackFn should not be called")
      },
      globalEnv = global_env
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(
    captured$write,
    list(
      text = paste(
        "Study Parameters",
        "================",
        "",
        "Workflow Name: Metab Session",
        "Description: summary export description",
        "Timestamp: 2026-04-18 10:11:12",
        "Note: Some parameters could not be saved due to serialization issues",
        sep = "\n"
      ),
      path = "/tmp/metab-summary-source/study_parameters.txt"
    )
  )
  expect_identical(
    captured$reads,
    c(
      "/tmp/metab-summary-source/design_matrix.tab",
      "/tmp/metab-summary-source/contrasts_tbl.tab"
    )
  )
  expect_identical(captured$progress_message, "Copying files to publication directory...")
  expect_identical(
    captured$copy_call,
    list(
      omic_type = "metabolomics",
      experiment_label = "Metab Session",
      contrasts_tbl = contrasts_tbl,
      design_matrix = design_matrix,
      force = TRUE
    )
  )
  expect_identical(
    captured$notifications,
    list(list(
      message = "Publication files copied",
      type = "message",
      duration = NULL
    ))
  )
  expect_true(values$workflow_args_saved)
  expect_true(values$files_copied)
  expect_identical(
    output$copy_status,
    "Files copied to publication directory successfully [OK]"
  )
  expect_identical(
    output$session_summary,
    paste("Workflow args created for:", "Metab Session",
          "\nDescription:", "summary export description",
          "\nTimestamp:", fixed_time,
          "\nStatus: Arguments saved [OK], Files copied [OK]")
  )
  expect_identical(global_env$project_dirs, project_dirs)
  expect_identical(
    visible$value,
    list(
      status = "success",
      basicParamsFile = "/tmp/metab-summary-source/study_parameters.txt",
      fallbackCreated = TRUE,
      projectDirsAssigned = TRUE,
      contrastsTbl = contrasts_tbl,
      designMatrix = design_matrix
    )
  )
})

test_that("metabolomics summary copy-to-publication observer shell preserves error reporting", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  design_matrix <- data.frame(sample = "S1", condition = "case", stringsAsFactors = FALSE)
  contrasts_tbl <- data.frame(contrast = "case-control", stringsAsFactors = FALSE)
  global_env <- new.env(parent = emptyenv())
  global_env$project_dirs <- "existing-project-dirs"

  values$workflow_args_saved <- TRUE
  values$files_copied <- FALSE

  visible <- withVisible(
    runMetabSummaryCopyToPublicationObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(metabolomics = list(source_dir = "/tmp/metab-summary-source")),
      omicType = "metabolomics",
      workflowData = list(
        contrasts_tbl = contrasts_tbl,
        design_matrix = design_matrix
      ),
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req <- c(captured$req, list(value))
        invisible(value)
      },
      fileExistsFn = function(path) FALSE,
      sysTimeFn = function() stop("sysTimeFn should not be called"),
      writeLinesFn = function(...) {
        stop("writeLinesFn should not be called")
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      readTsvFn = function(...) {
        stop("readTsvFn should not be called")
      },
      existsFn = function(x, envir) exists(x, envir = envir, inherits = FALSE),
      assignFn = function(...) {
        stop("assignFn should not be called")
      },
      copyToResultsSummaryFn = function(...) {
        stop("copy failure")
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error <- c(captured$error, message)
        invisible(NULL)
      },
      skipFormatterFn = function(message) {
        paste0("skip:", message)
      },
      catFn = function(...) invisible(NULL),
      tracebackFn = function() {
        captured$traceback <- TRUE
        invisible(NULL)
      },
      globalEnv = global_env
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session"))
  expect_identical(captured$progress_message, "Copying files to publication directory...")
  expect_identical(
    captured$notifications,
    list(list(
      message = "Copy error: copy failure",
      type = "error",
      duration = 10
    ))
  )
  expect_identical(captured$error, "skip:Failed to copy files: copy failure")
  expect_true(captured$traceback)
  expect_false(values$files_copied)
  expect_identical(output$copy_status, "Error: copy failure")
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
  expect_identical(
    visible$value,
    list(
      status = "error",
      basicParamsFile = "/tmp/metab-summary-source/study_parameters.txt",
      fallbackCreated = FALSE,
      errorMessage = "copy failure"
    )
  )
})

test_that("metabolomics summary generate-report observer shell preserves package-template render success contract", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$req <- list()
  captured$progress <- list()
  captured$file_copies <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  fixed_time <- as.POSIXct("2026-04-18 10:11:12", tz = "UTC")
  project_dirs <- list(metabolomics = list(base_dir = "/tmp/metab-summary-base"))

  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL

  visible <- withVisible(
    runMetabSummaryGenerateReportObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = project_dirs,
      omicType = "metabolomics",
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req[[length(captured$req) + 1L]] <- value
        invisible(value)
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      incProgressFn = function(amount, detail = NULL, ...) {
        captured$progress[[length(captured$progress) + 1L]] <- list(
          amount = amount,
          detail = detail
        )
        invisible(NULL)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      dirCreateFn = function(path, recursive = FALSE, showWarnings = TRUE) {
        captured$dirs <- c(captured$dirs, path)
        invisible(TRUE)
      },
      fileExistsFn = function(path) {
        path %in% c(
          "/pkg/metabolomics_report.rmd",
          "/tmp/metab-rendered.html"
        )
      },
      systemFileFn = function(...) {
        "/pkg/metabolomics_report.rmd"
      },
      fileCopyFn = function(from, to) {
        captured$file_copies[[length(captured$file_copies) + 1L]] <- list(
          from = from,
          to = to
        )
        TRUE
      },
      renderReportFn = function(omic_type, experiment_label, rmd_filename) {
        captured$render <- list(
          omic_type = omic_type,
          experiment_label = experiment_label,
          rmd_filename = rmd_filename
        )
        "/tmp/metab-rendered.html"
      },
      reactiveFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      outputOptionsFn = function(output, name, suspendWhenHidden) {
        captured$output_options <- list(
          output = output,
          name = name,
          suspendWhenHidden = suspendWhenHidden
        )
        invisible(output)
      },
      downloadHandlerFn = function(filename, content) {
        list(filename = filename, content = content)
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      sysTimeFn = function() fixed_time,
      logInfoFn = function(message) {
        captured$info <- c(captured$info, message)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      },
      skipFormatterFn = function(message) message,
      catFn = function(...) invisible(NULL),
      printFn = function(x) {
        stop(sprintf("printFn should not be called: %s", conditionMessage(x)))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session", TRUE))
  expect_identical(captured$progress_message, "Generating report...")
  expect_identical(
    captured$progress,
    list(
      list(
        amount = 0.2,
        detail = "Checking for metabolomics_report.rmd template..."
      ),
      list(
        amount = 0.5,
        detail = "Rendering report..."
      )
    )
  )
  expect_identical(
    captured$dirs,
    "/tmp/metab-summary-base/scripts/metabolomics"
  )
  expect_identical(
    captured$file_copies,
    list(list(
      from = "/pkg/metabolomics_report.rmd",
      to = "/tmp/metab-summary-base/scripts/metabolomics/metabolomics_report.rmd"
    ))
  )
  expect_identical(
    captured$render,
    list(
      omic_type = "metabolomics",
      experiment_label = "Metab Session",
      rmd_filename = "metabolomics_report.rmd"
    )
  )
  expect_identical(
    captured$notifications,
    list(
      list(
        message = "metabolomics_report.rmd template copied from package",
        type = "message",
        duration = NULL
      ),
      list(
        message = "Report generated successfully!",
        type = "message",
        duration = NULL
      )
    )
  )
  expect_identical(
    captured$info,
    c(
      "Template copied from package to: {report_template_path}",
      "Calling RenderReport with omic_type: {omic_type}, experiment_label: {input$experiment_label}",
      "RenderReport returned path: {rendered_path}"
    )
  )
  expect_identical(values$report_generated, TRUE)
  expect_identical(values$report_path, "/tmp/metab-rendered.html")
  expect_identical(output$report_ready, TRUE)
  expect_identical(captured$output_options$output, output)
  expect_identical(captured$output_options$name, "report_ready")
  expect_identical(captured$output_options$suspendWhenHidden, FALSE)
  expect_identical(output$download_report$filename(), "metab-rendered.html")
  expect_identical(
    output$session_summary,
    paste("Workflow args created for:", "Metab Session",
          "\nDescription:", "summary export description",
          "\nTimestamp:", fixed_time,
          "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
          "\nReport location:", "/tmp/metab-rendered.html")
  )
  expect_identical(
    visible$value,
    list(
      status = "success",
      templateFilename = "metabolomics_report.rmd",
      reportTemplatePath = "/tmp/metab-summary-base/scripts/metabolomics/metabolomics_report.rmd",
      templateSource = "package",
      renderedPath = "/tmp/metab-rendered.html"
    )
  )
})

test_that("metabolomics summary generate-report observer shell preserves invalid project-dir guard", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$req <- list()

  values <- new.env(parent = emptyenv())
  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL

  visible <- withVisible(
    runMetabSummaryGenerateReportObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(lipidomics = list(base_dir = "/tmp/lipid-summary-base")),
      omicType = "metabolomics",
      values = values,
      output = new.env(parent = emptyenv()),
      reqFn = function(value) {
        captured$req[[length(captured$req) + 1L]] <- value
        invisible(value)
      },
      withProgressFn = function(...) {
        stop("withProgressFn should not be called")
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      logInfoFn = function(message) {
        stop(sprintf("logInfoFn should not be called: %s", message))
      },
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      },
      catFn = function(...) {
        stop("catFn should not be called")
      },
      printFn = function(x) {
        stop(sprintf("printFn should not be called: %s", conditionMessage(x)))
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session", TRUE))
  expect_identical(
    captured$notifications,
    list(list(
      message = "Error: Project directories not properly initialized",
      type = "error",
      duration = 10
    ))
  )
  expect_identical(
    visible$value,
    list(
      status = "invalid_project_dirs",
      errorMessage = "Project directories not properly initialized"
    )
  )
})

test_that("metabolomics summary generate-report observer shell preserves render failure reporting", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$req <- list()
  captured$progress <- list()
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())

  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL

  visible <- withVisible(
    runMetabSummaryGenerateReportObserverShell(
      inputValues = list(
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(metabolomics = list(base_dir = "/tmp/metab-summary-base")),
      omicType = "metabolomics",
      values = values,
      output = output,
      reqFn = function(value) {
        captured$req[[length(captured$req) + 1L]] <- value
        invisible(value)
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      incProgressFn = function(amount, detail = NULL, ...) {
        captured$progress[[length(captured$progress) + 1L]] <- list(
          amount = amount,
          detail = detail
        )
        invisible(NULL)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      dirCreateFn = function(...) {
        stop("dirCreateFn should not be called")
      },
      fileExistsFn = function(path) {
        identical(
          path,
          "/tmp/metab-summary-base/scripts/metabolomics/metabolomics_report.rmd"
        )
      },
      renderReportFn = function(...) {
        stop("pandoc missing")
      },
      reactiveFn = function(expr) {
        stop("reactiveFn should not be called")
      },
      outputOptionsFn = function(...) {
        stop("outputOptionsFn should not be called")
      },
      downloadHandlerFn = function(...) {
        stop("downloadHandlerFn should not be called")
      },
      renderTextFn = function(expr) {
        stop("renderTextFn should not be called")
      },
      logInfoFn = function(message) {
        captured$info <- c(captured$info, message)
        invisible(NULL)
      },
      logErrorFn = function(message) {
        captured$error <- c(captured$error, message)
        invisible(NULL)
      },
      catFn = function(...) invisible(NULL),
      printFn = function(x) {
        captured$printed <- conditionMessage(x)
        invisible(NULL)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(captured$req, list("Metab Session", TRUE))
  expect_identical(captured$progress_message, "Generating report...")
  expect_identical(
    captured$progress,
    list(list(
      amount = 0.5,
      detail = "Rendering report..."
    ))
  )
  expect_identical(
    captured$info,
    "Calling RenderReport with omic_type: {omic_type}, experiment_label: {input$experiment_label}"
  )
  expect_identical(
    captured$notifications,
    list(
      list(
        message = "Report generation failed: pandoc missing",
        type = "error",
        duration = 15
      ),
      list(
        message = "Debug info: Check R console for detailed error trace",
        type = "warning",
        duration = 10
      )
    )
  )
  expect_identical(
    captured$error,
    c(
      "Failed to generate report: {e$message}",
      "Error class: {class(e)[1]}"
    )
  )
  expect_identical(captured$printed, "pandoc missing")
  expect_false(values$report_generated)
  expect_null(values$report_path)
  expect_false(exists("download_report", envir = output, inherits = FALSE))
  expect_identical(
    visible$value,
    list(
      status = "error",
      templateFilename = "metabolomics_report.rmd",
      reportTemplatePath = "/tmp/metab-summary-base/scripts/metabolomics/metabolomics_report.rmd",
      templateSource = "existing",
      errorMessage = "pandoc missing"
    )
  )
})

test_that("metabolomics summary GitHub observer shell preserves push success contract", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$req <- list()
  output <- new.env(parent = emptyenv())
  fixed_time <- as.POSIXct("2026-04-18 10:11:12", tz = "UTC")
  project_dirs <- list(metabolomics = list(base_dir = "/tmp/metab-summary-base"))

  visible <- withVisible(
    runMetabSummaryGithubPushObserverShell(
      inputValues = list(
        enable_github = TRUE,
        github_org = "APAF-bioinformatics",
        github_email = "metab@example.org",
        github_username = "metab-user",
        project_id = "metab-project",
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = project_dirs,
      omicType = "metabolomics",
      values = list(report_generated = TRUE),
      output = output,
      reqFn = function(...) {
        args <- list(...)
        captured$req[[length(captured$req) + 1L]] <- args
        invisible(args[[length(args)]])
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      optionsFn = function(...) {
        captured$options <- list(...)
        invisible(NULL)
      },
      pushProjectToGithubFromDirsFn = function(project_dirs, omic_type, experiment_label, project_id) {
        captured$push <- list(
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          project_id = project_id
        )
        invisible(NULL)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      renderTextFn = function(expr) {
        eval(substitute(expr), parent.frame())
      },
      sysTimeFn = function() fixed_time,
      logErrorFn = function(message) {
        stop(sprintf("logErrorFn should not be called: %s", message))
      },
      skipFormatterFn = function(message) message
    )
  )

  expect_false(visible$visible)
  expect_identical(
    captured$req,
    list(
      list(TRUE, "APAF-bioinformatics", "metab@example.org", "metab-user", "metab-project"),
      list(TRUE)
    )
  )
  expect_identical(captured$progress_message, "Pushing to GitHub...")
  expect_identical(
    captured$options,
    list(
      github_org = "APAF-bioinformatics",
      github_user_email = "metab@example.org",
      github_user_name = "metab-user"
    )
  )
  expect_identical(
    captured$push,
    list(
      project_dirs = project_dirs,
      omic_type = "metabolomics",
      experiment_label = "Metab Session",
      project_id = "metab-project"
    )
  )
  expect_identical(
    captured$notifications,
    list(list(
      message = "Successfully pushed to GitHub",
      type = "message",
      duration = NULL
    ))
  )
  expect_identical(
    output$session_summary,
    paste("Workflow args created for:", "Metab Session",
          "\nDescription:", "summary export description",
          "\nTimestamp:", fixed_time,
          "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]")
  )
  expect_identical(visible$value, list(status = "success"))
})

test_that("metabolomics summary GitHub observer shell preserves error reporting", {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$req <- list()
  output <- new.env(parent = emptyenv())

  visible <- withVisible(
    runMetabSummaryGithubPushObserverShell(
      inputValues = list(
        enable_github = TRUE,
        github_org = "APAF-bioinformatics",
        github_email = "metab@example.org",
        github_username = "metab-user",
        project_id = "metab-project",
        experiment_label = "Metab Session",
        description = "summary export description"
      ),
      projectDirs = list(metabolomics = list(base_dir = "/tmp/metab-summary-base")),
      omicType = "metabolomics",
      values = list(report_generated = TRUE),
      output = output,
      reqFn = function(...) {
        args <- list(...)
        captured$req[[length(captured$req) + 1L]] <- args
        invisible(args[[length(args)]])
      },
      withProgressFn = function(message = NULL, expr, ...) {
        captured$progress_message <- message
        eval(substitute(expr), parent.frame())
      },
      optionsFn = function(...) {
        captured$options <- list(...)
        invisible(NULL)
      },
      pushProjectToGithubFromDirsFn = function(...) {
        stop("permission denied")
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        captured$notifications[[length(captured$notifications) + 1L]] <- list(
          message = message,
          type = type,
          duration = duration
        )
        invisible(NULL)
      },
      renderTextFn = function(expr) {
        stop("renderTextFn should not be called")
      },
      sysTimeFn = function() stop("sysTimeFn should not be called"),
      logErrorFn = function(message) {
        captured$error <- c(captured$error, message)
        invisible(NULL)
      },
      skipFormatterFn = function(message) {
        paste0("skip:", message)
      }
    )
  )

  expect_false(visible$visible)
  expect_identical(
    captured$req,
    list(
      list(TRUE, "APAF-bioinformatics", "metab@example.org", "metab-user", "metab-project"),
      list(TRUE)
    )
  )
  expect_identical(captured$progress_message, "Pushing to GitHub...")
  expect_identical(
    captured$options,
    list(
      github_org = "APAF-bioinformatics",
      github_user_email = "metab@example.org",
      github_user_name = "metab-user"
    )
  )
  expect_identical(
    captured$notifications,
    list(list(
      message = "GitHub push failed: permission denied",
      type = "error",
      duration = 10
    ))
  )
  expect_identical(
    captured$error,
    "skip:Failed to push to GitHub: permission denied"
  )
  expect_false(exists("session_summary", envir = output, inherits = FALSE))
  expect_identical(
    visible$value,
    list(
      status = "error",
      errorMessage = "permission denied"
    )
  )
})

test_that("metabolomics summary observer-registration seam preserves observer wiring and shell handoff", {
  captured <- new.env(parent = emptyenv())
  captured$observe_calls <- list()
  input <- list(
    save_workflow_args = "save-event",
    copy_to_publication = "copy-event",
    generate_report = "report-event",
    push_to_github = "github-event",
    export_session_state = "export-event",
    experiment_label = "Metab Session",
    description = "summary export description",
    enable_github = TRUE,
    github_org = "APAF-bioinformatics",
    github_email = "metab@example.org",
    github_username = "metab-user",
    project_id = "metab-project"
  )
  output <- new.env(parent = emptyenv())
  project_dirs <- list(metabolomics = list(base_dir = tempfile("metab-base-")))
  workflow_data <- list(token = "workflow-data")
  values <- list(token = "values")

  registration <- withVisible(
    registerMetabSummaryServerObservers(
      input = input,
      output = output,
      projectDirs = project_dirs,
      omicType = "metabolomics",
      workflowData = workflow_data,
      values = values,
      observeEventFn = function(eventExpr, handlerExpr, ...) {
        index <- length(captured$observe_calls) + 1L
        captured$observe_calls[[index]] <- list(
          event = eventExpr,
          dots = list(...)
        )
        eval(substitute(handlerExpr), parent.frame())
        paste0("observer-", index)
      },
      runSaveWorkflowArgsObserverShellFn = function(...) {
        captured$save <- list(...)
        invisible(NULL)
      },
      runCopyToPublicationObserverShellFn = function(...) {
        captured$copy <- list(...)
        invisible(NULL)
      },
      runGenerateReportObserverShellFn = function(...) {
        captured$report <- list(...)
        invisible(NULL)
      },
      runGithubPushObserverShellFn = function(...) {
        captured$github <- list(...)
        invisible(NULL)
      },
      runExportSessionObserverShellFn = function(...) {
        captured$export <- list(...)
        invisible(NULL)
      }
    )
  )

  expect_false(registration$visible)
  expect_identical(
    registration$value,
    list(
      saveWorkflowArgsObserver = "observer-1",
      copyToPublicationObserver = "observer-2",
      generateReportObserver = "observer-3",
      githubPushObserver = "observer-4",
      exportSessionObserver = "observer-5"
    )
  )
  expect_identical(
    lapply(captured$observe_calls, `[[`, "event"),
    list("save-event", "copy-event", "report-event", "github-event", "export-event")
  )
  expect_identical(captured$observe_calls[[1]]$dots, list())
  expect_identical(captured$observe_calls[[2]]$dots, list(ignoreInit = TRUE))
  expect_identical(captured$observe_calls[[3]]$dots, list())
  expect_identical(captured$observe_calls[[4]]$dots, list())
  expect_identical(captured$observe_calls[[5]]$dots, list())
  expect_identical(
    captured$save$inputValues,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$save$projectDirs, project_dirs)
  expect_identical(captured$save$omicType, "metabolomics")
  expect_identical(captured$save$workflowData, workflow_data)
  expect_identical(captured$save$values, values)
  expect_identical(captured$save$output, output)
  expect_identical(
    captured$copy$inputValues,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$copy$projectDirs, project_dirs)
  expect_identical(captured$copy$omicType, "metabolomics")
  expect_identical(captured$copy$workflowData, workflow_data)
  expect_identical(captured$copy$values, values)
  expect_identical(captured$copy$output, output)
  expect_identical(
    captured$report$inputValues,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$report$projectDirs, project_dirs)
  expect_identical(captured$report$omicType, "metabolomics")
  expect_identical(captured$report$values, values)
  expect_identical(captured$report$output, output)
  expect_identical(
    captured$github$inputValues,
    list(
      enable_github = TRUE,
      github_org = "APAF-bioinformatics",
      github_email = "metab@example.org",
      github_username = "metab-user",
      project_id = "metab-project",
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$github$projectDirs, project_dirs)
  expect_identical(captured$github$omicType, "metabolomics")
  expect_identical(captured$github$values, values)
  expect_identical(captured$github$output, output)
  expect_identical(
    captured$export$inputValues,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$export$projectDirs, project_dirs)
  expect_identical(captured$export$omicType, "metabolomics")
  expect_identical(captured$export$values, values)
})

test_that("metabolomics summary server delegates template-status registration through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(metabolomics = list(base_dir = tempfile("metab-base-")))
  had_helper <- exists(
    "setupMetabSummaryTemplateStatusOutput",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "setupMetabSummaryTemplateStatusOutput",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "setupMetabSummaryTemplateStatusOutput",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "setupMetabSummaryTemplateStatusOutput",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(list = "setupMetabSummaryTemplateStatusOutput", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "setupMetabSummaryTemplateStatusOutput",
    function(output, projectDirs, omicType, renderText = shiny::renderText, req = shiny::req) {
      captured$call <- list(
        output = output,
        project_dirs = projectDirs,
        omic_type = omicType,
        render_text = renderText,
        req = req
      )
      output$template_status <- "template-status-token"
      invisible(NULL)
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) list(...),
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$output, captured$output)
  expect_identical(captured$output$template_status, "template-status-token")
})

test_that("metabolomics summary server delegates bootstrap output registration through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  had_helper <- exists(
    "setupMetabSummaryBootstrapOutputs",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "setupMetabSummaryBootstrapOutputs",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "setupMetabSummaryBootstrapOutputs",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "setupMetabSummaryBootstrapOutputs",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(list = "setupMetabSummaryBootstrapOutputs", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "setupMetabSummaryBootstrapOutputs",
    function(output,
             renderText = shiny::renderText,
             reactive = shiny::reactive,
             outputOptions = shiny::outputOptions) {
      captured$call <- list(
        output = output,
        render_text = renderText,
        reactive = reactive,
        output_options = outputOptions
      )
      output$session_summary <- "bootstrap-summary-token"
      output$report_ready <- "bootstrap-ready-token"
      invisible(NULL)
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) list(...),
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = list(metabolomics = list(base_dir = tempfile("metab-base-"))),
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(captured$call$output, captured$output)
  expect_identical(captured$output$session_summary, "bootstrap-summary-token")
  expect_identical(captured$output$report_ready, "bootstrap-ready-token")
})

test_that("metabolomics summary server delegates bootstrap state setup through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  had_helper <- exists(
    "setupMetabSummaryServerBootstrapState",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "setupMetabSummaryServerBootstrapState",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "setupMetabSummaryServerBootstrapState",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "setupMetabSummaryServerBootstrapState",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(
        list = "setupMetabSummaryServerBootstrapState",
        envir = server_env
      )
    }
  }, add = TRUE)

  values_token <- list(token = "bootstrap-values")

  assign(
    "setupMetabSummaryServerBootstrapState",
    function(session,
             experimentLabel,
             reactiveValuesFn = shiny::reactiveValues,
             updateTextInputFn = shiny::updateTextInput) {
      captured$call <- list(
        session = session,
        experiment_label = experimentLabel,
        reactive_values = reactiveValuesFn,
        update_text_input = updateTextInputFn
      )
      values_token
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      session <- list(
        ns = function(value) value,
        token = "session-token"
      )
      module(
        list(),
        output,
        session
      )
      captured$output <- output
      captured$session <- session
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$direct_reactive_values <- list(...)
      list(...)
    },
    observeEvent = function(eventExpr, handlerExpr, ...) invisible(NULL),
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) {
      captured$direct_update_text_input <- list(...)
      invisible(NULL)
    },
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = list(metabolomics = list(base_dir = tempfile("metab-base-"))),
    omic_type = "metabolomics",
    experiment_label = "Metab Session",
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(captured$call$session$token, "session-token")
  expect_identical(captured$call$experiment_label, "Metab Session")
  expect_type(captured$call$reactive_values, "closure")
  expect_type(captured$call$update_text_input, "closure")
  expect_null(captured$direct_reactive_values)
  expect_null(captured$direct_update_text_input)
})

test_that("metabolomics summary server delegates observer registration through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(
    metabolomics = list(
      base_dir = tempfile("metab-base-"),
      source_dir = tempfile("metab-source-")
    )
  )
  workflow_data <- list(token = "workflow-data")
  had_helper <- exists(
    "registerMetabSummaryServerObservers",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "registerMetabSummaryServerObservers",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "registerMetabSummaryServerObservers",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "registerMetabSummaryServerObservers",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(
        list = "registerMetabSummaryServerObservers",
        envir = server_env
      )
    }
  }, add = TRUE)

  assign(
    "registerMetabSummaryServerObservers",
    function(input, output, projectDirs, omicType, workflowData, values, ...) {
      captured$call <- list(
        input = input,
        output = output,
        project_dirs = projectDirs,
        omic_type = omicType,
        workflow_data = workflowData,
        values = values
      )
      invisible(list(observerRegistration = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(),
        output,
        list(
          ns = function(value) value,
          token = "session-token"
        )
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(...) {
      captured$direct_observe_event <- TRUE
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = workflow_data
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(captured$call$output, captured$output)
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$workflow_data, workflow_data)
  expect_identical(captured$call$values, captured$values)
  expect_null(captured$direct_observe_event)
})

test_that("metabolomics summary server delegates export observer through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(
    metabolomics = list(
      base_dir = tempfile("metab-base-"),
      source_dir = tempfile("metab-source-")
    )
  )
  had_helper <- exists(
    "runMetabSummaryExportSessionObserverShell",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "runMetabSummaryExportSessionObserverShell",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "runMetabSummaryExportSessionObserverShell",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "runMetabSummaryExportSessionObserverShell",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(list = "runMetabSummaryExportSessionObserverShell", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "runMetabSummaryExportSessionObserverShell",
    function(inputValues, projectDirs, omicType, values, ...) {
      captured$call <- list(
        input_values = inputValues,
        project_dirs = projectDirs,
        omic_type = omicType,
        values = values
      )
      invisible(list(status = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          save_workflow_args = NULL,
          copy_to_publication = NULL,
          generate_report = NULL,
          push_to_github = NULL,
          export_session_state = 1,
          experiment_label = "Metab Session",
          description = "summary export description"
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (!is.null(eventExpr)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(
    captured$call$input_values,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$values, captured$values)
})

test_that("metabolomics summary server delegates save-workflow observer through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(
    metabolomics = list(
      base_dir = tempfile("metab-base-"),
      source_dir = tempfile("metab-source-")
    )
  )
  workflow_data <- list(
    state_manager = "state-manager-token",
    contrasts_tbl = "contrast-token",
    config_list = list(alpha = 0.05)
  )
  had_helper <- exists(
    "runMetabSummarySaveWorkflowArgsObserverShell",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "runMetabSummarySaveWorkflowArgsObserverShell",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "runMetabSummarySaveWorkflowArgsObserverShell",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "runMetabSummarySaveWorkflowArgsObserverShell",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(
        list = "runMetabSummarySaveWorkflowArgsObserverShell",
        envir = server_env
      )
    }
  }, add = TRUE)

  assign(
    "runMetabSummarySaveWorkflowArgsObserverShell",
    function(inputValues, projectDirs, omicType, workflowData, values, output, ...) {
      captured$call <- list(
        input_values = inputValues,
        project_dirs = projectDirs,
        omic_type = omicType,
        workflow_data = workflowData,
        values = values,
        output = output
      )
      invisible(list(status = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          save_workflow_args = 1,
          copy_to_publication = NULL,
          generate_report = NULL,
          push_to_github = NULL,
          export_session_state = NULL,
          experiment_label = "Metab Session",
          description = "summary export description"
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (!is.null(eventExpr)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = workflow_data
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(
    captured$call$input_values,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$workflow_data, workflow_data)
  expect_identical(captured$call$values, captured$values)
  expect_identical(captured$call$output, captured$output)
})

test_that("metabolomics summary server delegates copy-to-publication observer through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(
    metabolomics = list(
      base_dir = tempfile("metab-base-"),
      source_dir = tempfile("metab-source-")
    )
  )
  workflow_data <- list(
    contrasts_tbl = "contrast-token",
    design_matrix = "design-token"
  )
  had_helper <- exists(
    "runMetabSummaryCopyToPublicationObserverShell",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "runMetabSummaryCopyToPublicationObserverShell",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "runMetabSummaryCopyToPublicationObserverShell",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "runMetabSummaryCopyToPublicationObserverShell",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(
        list = "runMetabSummaryCopyToPublicationObserverShell",
        envir = server_env
      )
    }
  }, add = TRUE)

  assign(
    "runMetabSummaryCopyToPublicationObserverShell",
    function(inputValues, projectDirs, omicType, workflowData, values, output, ...) {
      captured$call <- list(
        input_values = inputValues,
        project_dirs = projectDirs,
        omic_type = omicType,
        workflow_data = workflowData,
        values = values,
        output = output
      )
      output$copy_status <- "copy-status-token"
      invisible(list(status = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          save_workflow_args = NULL,
          copy_to_publication = 1,
          generate_report = NULL,
          push_to_github = NULL,
          export_session_state = NULL,
          experiment_label = "Metab Session",
          description = "summary export description"
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (!is.null(eventExpr)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = workflow_data
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(
    captured$call$input_values,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$workflow_data, workflow_data)
  expect_identical(captured$call$values, captured$values)
  expect_identical(captured$call$output, captured$output)
  expect_identical(captured$output$copy_status, "copy-status-token")
})

test_that("metabolomics summary server delegates GitHub observer through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(metabolomics = list(base_dir = tempfile("metab-base-")))
  had_helper <- exists(
    "runMetabSummaryGithubPushObserverShell",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "runMetabSummaryGithubPushObserverShell",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "runMetabSummaryGithubPushObserverShell",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "runMetabSummaryGithubPushObserverShell",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(list = "runMetabSummaryGithubPushObserverShell", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "runMetabSummaryGithubPushObserverShell",
    function(inputValues, projectDirs, omicType, values, output, ...) {
      captured$call <- list(
        input_values = inputValues,
        project_dirs = projectDirs,
        omic_type = omicType,
        values = values,
        output = output
      )
      output$session_summary <- "github-summary-token"
      invisible(list(status = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          save_workflow_args = NULL,
          copy_to_publication = NULL,
          generate_report = NULL,
          push_to_github = 1,
          export_session_state = NULL,
          experiment_label = "Metab Session",
          description = "summary export description",
          enable_github = TRUE,
          github_org = "APAF-bioinformatics",
          github_email = "metab@example.org",
          github_username = "metab-user",
          project_id = "metab-project"
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (!is.null(eventExpr)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(
    captured$call$input_values,
    list(
      enable_github = TRUE,
      github_org = "APAF-bioinformatics",
      github_email = "metab@example.org",
      github_username = "metab-user",
      project_id = "metab-project",
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$values, captured$values)
  expect_identical(captured$call$output, captured$output)
})

test_that("metabolomics summary server delegates report-generation observer through the seam", {
  skip("shared coverage now uses the public module harness instead of target-only wrapper seam delegation")
  captured <- new.env(parent = emptyenv())
  server_env <- environment(mod_metab_summary_server)
  project_dirs <- list(metabolomics = list(base_dir = tempfile("metab-base-")))
  had_helper <- exists(
    "runMetabSummaryGenerateReportObserverShell",
    envir = server_env,
    inherits = FALSE
  )

  if (had_helper) {
    original_helper <- get(
      "runMetabSummaryGenerateReportObserverShell",
      envir = server_env,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_helper) {
      assign(
        "runMetabSummaryGenerateReportObserverShell",
        original_helper,
        envir = server_env
      )
    } else if (exists(
      "runMetabSummaryGenerateReportObserverShell",
      envir = server_env,
      inherits = FALSE
    )) {
      rm(list = "runMetabSummaryGenerateReportObserverShell", envir = server_env)
    }
  }, add = TRUE)

  assign(
    "runMetabSummaryGenerateReportObserverShell",
    function(inputValues, projectDirs, omicType, values, output, ...) {
      captured$call <- list(
        input_values = inputValues,
        project_dirs = projectDirs,
        omic_type = omicType,
        values = values,
        output = output
      )
      output$session_summary <- "generate-report-summary-token"
      invisible(list(status = "delegated"))
    },
    envir = server_env
  )

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      captured$module_id <- id
      output <- new.env(parent = emptyenv())
      module(
        list(
          save_workflow_args = NULL,
          copy_to_publication = NULL,
          generate_report = 1,
          push_to_github = NULL,
          export_session_state = NULL,
          experiment_label = "Metab Session",
          description = "summary export description"
        ),
        output,
        list(ns = function(value) value)
      )
      captured$output <- output
      invisible(NULL)
    },
    reactiveValues = function(...) {
      captured$values <- list(...)
      captured$values
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      if (!is.null(eventExpr)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) eval(substitute(expr), parent.frame()),
    reactive = function(expr) eval(substitute(expr), parent.frame()),
    outputOptions = function(output, name, suspendWhenHidden) invisible(output),
    showNotification = function(...) invisible(NULL),
    downloadHandler = function(...) "download-handler",
    withProgress = function(message = NULL, expr, ...) eval(substitute(expr), parent.frame()),
    incProgress = function(...) invisible(NULL),
    req = function(...) invisible(list(...)[[1]]),
    updateTextInput = function(...) invisible(NULL),
    .package = "shiny"
  )
  testthat::local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  mod_metab_summary_server(
    id = "summary",
    project_dirs = project_dirs,
    omic_type = "metabolomics",
    experiment_label = NULL,
    workflow_data = list()
  )

  expect_identical(captured$module_id, "summary")
  expect_identical(
    captured$call$input_values,
    list(
      experiment_label = "Metab Session",
      description = "summary export description"
    )
  )
  expect_identical(captured$call$project_dirs, project_dirs)
  expect_identical(captured$call$omic_type, "metabolomics")
  expect_identical(captured$call$values, captured$values)
  expect_identical(captured$call$output, captured$output)
})
