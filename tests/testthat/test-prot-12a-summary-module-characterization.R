# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

if (!methods::isClass("mockProtSummaryArgsCarrier")) {
  methods::setClass("mockProtSummaryArgsCarrier", slots = c(args = "list"))
}

test_that("proteomics summary module preserves startup outputs and template warning", {
  fixture_dir <- tempfile("prot-summary-module-")
  source_dir <- file.path(fixture_dir, "source")
  integration_dir <- file.path(fixture_dir, "integration")
  dir.create(source_dir, recursive = TRUE)
  dir.create(integration_dir, recursive = TRUE)
  dir.create(file.path(fixture_dir, "scripts", "proteomics"), recursive = TRUE)

  testServer(
    mod_prot_summary_server,
    args = list(
      project_dirs = list(
        proteomics = list(
          source_dir = source_dir,
          integration_dir = integration_dir,
          base_dir = fixture_dir
        )
      ),
      experiment_label = "summary-demo",
      workflow_data = list()
    ),
    {
      session$flushReact()

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

test_that("proteomics summary state and template helpers preserve selection and directory behavior", {
  skip_if_not(
    exists("resolveProtSummaryFinalS4State", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("buildProtSummaryTemplateStatus", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("validateProtSummaryProjectDirs", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  resolve_final_state <- get("resolveProtSummaryFinalS4State", envir = asNamespace("MultiScholaR"))
  build_template_status <- get("buildProtSummaryTemplateStatus", envir = asNamespace("MultiScholaR"))
  validate_project_dirs <- get("validateProtSummaryProjectDirs", envir = asNamespace("MultiScholaR"))

  recorded_messages <- character()
  final_s4 <- methods::new(
    "mockProtSummaryArgsCarrier",
    args = list(globalParameters = list(workflow_type = "DIA"))
  )
  workflow_data <- list(
    state_manager = list(
      getHistory = function() c("imputed", "protein_replicate_filtered"),
      getState = function(state_name) {
        if (identical(state_name, "protein_replicate_filtered")) {
          return(final_s4)
        }
        NULL
      }
    )
  )

  result <- resolve_final_state(
    workflowData = workflow_data,
    catFn = function(...) {
      recorded_messages <<- c(recorded_messages, paste0(..., collapse = ""))
      invisible(NULL)
    }
  )

  expect_identical(result$dataStateUsed, "protein_replicate_filtered")
  expect_identical(result$finalS4Object, final_s4)
  expect_true(any(grepl("Retrieved DATA S4 object", recorded_messages, fixed = TRUE)))

  missing_result <- resolve_final_state(
    workflowData = list(state_manager = NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_null(missing_result$finalS4Object)
  expect_null(missing_result$dataStateUsed)

  fixture_dir <- tempfile("prot-summary-template-")
  template_dir <- file.path(fixture_dir, "scripts", "proteomics")
  dir.create(template_dir, recursive = TRUE)

  expect_identical(
    build_template_status(
      projectDirs = list(proteomics = list(base_dir = fixture_dir)),
      omicType = "proteomics"
    ),
    "[WARNING] Report templates will be downloaded when generating report"
  )

  writeLines("report", file.path(template_dir, "DIANN_report.rmd"))
  writeLines("report", file.path(template_dir, "TMT_report.rmd"))
  status_text <- build_template_status(
    projectDirs = list(proteomics = list(base_dir = fixture_dir)),
    omicType = "proteomics"
  )
  expect_match(status_text, "DIA-NN \\[OK\\]")
  expect_match(status_text, "TMT \\[OK\\]")

  notifications <- list()
  expect_true(
    validate_project_dirs(
      projectDirs = list(proteomics = list(base_dir = fixture_dir)),
      omicType = "proteomics",
      showNotificationFn = function(...) {
        notifications <<- append(notifications, list(list(...)))
        invisible(NULL)
      }
    )
  )
  expect_false(
    validate_project_dirs(
      projectDirs = list(),
      omicType = "proteomics",
      showNotificationFn = function(...) {
        notifications <<- append(notifications, list(list(...)))
        invisible(NULL)
      }
    )
  )
  expect_identical(notifications[[1]]$type, "error")
})

test_that("proteomics summary report helpers preserve template retrieval and activation behavior", {
  skip_if_not(
    exists("resolveProtSummaryReportTemplate", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("ensureProtSummaryReportTemplate", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("retrieveProtSummaryReportTemplateAsset", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("activateProtSummaryRenderedReport", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  resolve_report_template <- get("resolveProtSummaryReportTemplate", envir = asNamespace("MultiScholaR"))
  ensure_report_template <- get("ensureProtSummaryReportTemplate", envir = asNamespace("MultiScholaR"))
  retrieve_report_template <- get("retrieveProtSummaryReportTemplateAsset", envir = asNamespace("MultiScholaR"))
  activate_rendered_report <- get("activateProtSummaryRenderedReport", envir = asNamespace("MultiScholaR"))

  workflow_data <- list(config_list = list(globalParameters = list(workflow_type = "tmt")))
  template_choice <- resolve_report_template(
    workflowData = workflow_data,
    catFn = function(...) invisible(NULL)
  )
  expect_identical(template_choice$templateFilename, "TMT_report.rmd")

  default_choice <- resolve_report_template(
    workflowData = list(state_manager = NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(default_choice$templateFilename, "DIANN_report.rmd")

  s4_choice <- resolve_report_template(
    workflowData = list(
      state_manager = list(
        getHistory = function() "imputed",
        getState = function(state_name) {
          methods::new(
            "mockProtSummaryArgsCarrier",
            args = list(globalParameters = list(workflow_type = "lfq"))
          )
        }
      )
    ),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(s4_choice$templateFilename, "LFQ_report.rmd")

  config_list <- list(globalParameters = list(workflow_type = "tmt_pd"))
  assign("config_list", config_list, envir = .GlobalEnv)
  on.exit(rm("config_list", envir = .GlobalEnv), add = TRUE)
  global_choice <- resolve_report_template(
    workflowData = list(state_manager = NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(global_choice$templateFilename, "TMT_report.rmd")

  fixture_dir <- tempfile("prot-summary-report-template-")
  template_dir <- file.path(fixture_dir, "scripts", "proteomics")
  dir.create(template_dir, recursive = TRUE)
  report_template_path <- file.path(template_dir, "DIANN_report.rmd")
  package_template <- tempfile("pkg-template-", fileext = ".rmd")
  writeLines("report", package_template)
  notifications <- list()
  copied_paths <- list()

  ensured <- ensure_report_template(
    projectDirs = list(proteomics = list(base_dir = fixture_dir)),
    templateFilename = "DIANN_report.rmd",
    systemFileFn = function(...) package_template,
    fileExistsFn = function(path) file.exists(path),
    dirCreateFn = function(path, recursive = TRUE, showWarnings = FALSE) {
      dir.create(path, recursive = recursive, showWarnings = showWarnings)
    },
    fileCopyFn = function(from, to, ...) {
      copied_paths <<- append(copied_paths, list(list(from = from, to = to)))
      file.copy(from, to, overwrite = TRUE)
    },
    showNotificationFn = function(...) {
      notifications <<- append(notifications, list(list(...)))
      invisible(NULL)
    },
    logInfoFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(ensured$templateSource, "package")
  expect_true(file.exists(report_template_path))
  expect_identical(copied_paths[[1]]$from, package_template)

  existing <- ensure_report_template(
    projectDirs = list(proteomics = list(base_dir = fixture_dir)),
    templateFilename = "DIANN_report.rmd",
    systemFileFn = function(...) "",
    fileExistsFn = function(path) file.exists(path),
    dirCreateFn = function(path, recursive = TRUE, showWarnings = FALSE) {
      dir.create(path, recursive = recursive, showWarnings = showWarnings)
    },
    fileCopyFn = function(...) stop("should not copy when template already exists"),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(existing$templateSource, "existing")

  github_fixture_dir <- tempfile("prot-summary-report-github-")
  github_template_dir <- file.path(github_fixture_dir, "scripts", "proteomics")
  dir.create(github_template_dir, recursive = TRUE)
  github_template_path <- file.path(github_template_dir, "LFQ_report.rmd")
  github_downloads <- character()
  github_ensured <- ensure_report_template(
    projectDirs = list(proteomics = list(base_dir = github_fixture_dir)),
    templateFilename = "LFQ_report.rmd",
    systemFileFn = function(...) "",
    fileExistsFn = function(path) file.exists(path),
    dirCreateFn = function(path, recursive = TRUE, showWarnings = FALSE) {
      dir.create(path, recursive = recursive, showWarnings = showWarnings)
    },
    fileCopyFn = function(...) stop("should not copy package template"),
    downloadFileFn = function(url, destfile, quiet = TRUE) {
      github_downloads <<- c(github_downloads, url)
      writeLines("report", destfile)
      invisible(NULL)
    },
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_identical(github_ensured$templateSource, "github")
  expect_true(file.exists(github_template_path))
  expect_length(github_downloads, 1L)

  retrieved <- retrieve_report_template(
    projectDirs = list(proteomics = list(base_dir = fixture_dir)),
    templateFilename = "DIANN_report.rmd",
    ensureTemplateFn = function(...) ensured,
    showNotificationFn = function(...) invisible(NULL),
    logErrorFn = function(...) stop("unexpected retrieval error")
  )
  expect_identical(retrieved$templateSource, "package")

  error_notifications <- list()
  expect_null(
    retrieve_report_template(
      projectDirs = list(proteomics = list(base_dir = fixture_dir)),
      templateFilename = "broken.rmd",
      ensureTemplateFn = function(...) stop("template boom"),
      showNotificationFn = function(...) {
        error_notifications <<- append(error_notifications, list(list(...)))
        invisible(NULL)
      },
      logErrorFn = function(...) invisible(NULL)
    )
  )
  expect_identical(error_notifications[[1]]$type, "error")

  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$report_generated <- FALSE
  values$report_path <- NULL
  rendered_path <- tempfile("prot-summary-rendered-", fileext = ".html")
  writeLines("<html>report</html>", rendered_path)
  activation_notifications <- list()

  expect_true(
    activate_rendered_report(
      output = output,
      values = values,
      renderedPath = rendered_path,
      experimentLabel = "summary-demo",
      description = "demo run",
      reactiveFn = function(expr) eval(substitute(expr), parent.frame()),
      outputOptionsFn = function(...) invisible(NULL),
      downloadHandlerFn = function(filename, content) list(filename = filename, content = content),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) {
        activation_notifications <<- append(activation_notifications, list(list(...)))
        invisible(NULL)
      },
      timestampFn = function() as.POSIXct("2026-04-20 12:00:00", tz = "UTC")
    )
  )
  expect_true(values$report_generated)
  expect_identical(values$report_path, rendered_path)
  expect_identical(output$report_ready, TRUE)
  expect_match(output$session_summary, "Report generated \\[OK\\]")

  expect_false(
    activate_rendered_report(
      output = new.env(parent = emptyenv()),
      values = shiny::reactiveValues(report_generated = FALSE, report_path = NULL),
      renderedPath = tempfile("prot-summary-missing-", fileext = ".html"),
      experimentLabel = "summary-demo",
      description = "demo run"
    )
  )
})

test_that("proteomics summary report execution helpers preserve success and failure orchestration", {
  skip_if_not(
    exists("runProtSummaryReportGeneration", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("runProtSummaryReportProgress", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  run_report_generation <- get("runProtSummaryReportGeneration", envir = asNamespace("MultiScholaR"))
  run_report_progress <- get("runProtSummaryReportProgress", envir = asNamespace("MultiScholaR"))
  output <- new.env(parent = emptyenv())
  values <- shiny::reactiveValues(report_generated = FALSE, report_path = NULL)
  generated_path <- tempfile("prot-summary-generated-", fileext = ".html")
  writeLines("<html>report</html>", generated_path)
  notifications <- list()
  progress_details <- character()

  expect_true(
    run_report_generation(
      output = output,
      values = values,
      omicType = "proteomics",
      experimentLabel = "summary-demo",
      description = "demo run",
      templateFilename = "DIANN_report.rmd",
      renderReportAvailableFn = function() TRUE,
      renderReportFn = function(...) generated_path,
      activateReportFn = function(...) TRUE,
      showNotificationFn = function(...) {
        notifications <<- append(notifications, list(list(...)))
        invisible(NULL)
      },
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) stop("unexpected report error"),
      catFn = function(...) invisible(NULL),
      printFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL)
    )
  )

  expect_false(
    run_report_generation(
      output = output,
      values = values,
      omicType = "proteomics",
      experimentLabel = "summary-demo",
      description = "demo run",
      templateFilename = "DIANN_report.rmd",
      renderReportAvailableFn = function() TRUE,
      renderReportFn = function(...) stop("render boom"),
      activateReportFn = function(...) TRUE,
      showNotificationFn = function(...) {
        notifications <<- append(notifications, list(list(...)))
        invisible(NULL)
      },
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      printFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL)
    )
  )

  progress_result <- run_report_progress(
    output = output,
    values = values,
    workflowData = list(config_list = list(globalParameters = list(workflow_type = "lfq"))),
    projectDirs = list(proteomics = list(base_dir = tempdir())),
    omicType = "proteomics",
    experimentLabel = "summary-demo",
    description = "demo run",
    resolveReportTemplateFn = function(...) list(templateFilename = "LFQ_report.rmd"),
    retrieveTemplateAssetFn = function(...) list(reportTemplatePath = "/tmp/LFQ_report.rmd"),
    runReportGenerationFn = function(...) "generated",
    incProgressFn = function(value, detail = NULL) {
      progress_details <<- c(progress_details, paste(value, detail))
      invisible(NULL)
    }
  )
  expect_identical(progress_result, "generated")
  expect_length(progress_details, 3L)

  expect_false(
    run_report_progress(
      output = output,
      values = values,
      workflowData = list(),
      projectDirs = list(proteomics = list(base_dir = tempdir())),
      omicType = "proteomics",
      experimentLabel = "summary-demo",
      description = "demo run",
      resolveReportTemplateFn = function(...) list(templateFilename = "DIANN_report.rmd"),
      retrieveTemplateAssetFn = function(...) NULL,
      runReportGenerationFn = function(...) stop("should not run"),
      incProgressFn = function(...) invisible(NULL)
    )
  )
})

test_that("proteomics summary save and export helpers preserve success and fallback behavior", {
  skip_if_not(
    exists("completeProtSummaryWorkflowArgsSave", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("prepareProtSummarySessionStateExport", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("completeProtSummarySessionStateExport", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  complete_workflow_save <- get("completeProtSummaryWorkflowArgsSave", envir = asNamespace("MultiScholaR"))
  prepare_session_export <- get("prepareProtSummarySessionStateExport", envir = asNamespace("MultiScholaR"))
  complete_session_export <- get("completeProtSummarySessionStateExport", envir = asNamespace("MultiScholaR"))

  fixture_dir <- tempfile("prot-summary-save-")
  source_dir <- file.path(fixture_dir, "source")
  integration_dir <- file.path(fixture_dir, "integration")
  dir.create(source_dir, recursive = TRUE)
  dir.create(integration_dir, recursive = TRUE)

  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  final_s4 <- methods::new("mockProtSummaryArgsCarrier", args = list(globalParameters = list(workflow_type = "DIA")))
  notifications <- list()

  expect_true(
    complete_workflow_save(
      output = output,
      values = values,
      projectDirs = list(proteomics = list(source_dir = source_dir, integration_dir = integration_dir, base_dir = fixture_dir)),
      workflowData = list(
        contrasts_tbl = data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE),
        config_list = list(globalParameters = list(workflow_type = "DIA"))
      ),
      experimentLabel = "summary-demo",
      description = "demo run",
      resolveFinalS4StateFn = function(...) list(finalS4Object = final_s4),
      assignFn = function(...) invisible(NULL),
      createWorkflowArgsFn = function(...) file.path(source_dir, "study_parameters.txt"),
      saveRDSFn = function(object, file) invisible(NULL),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) {
        notifications <<- append(notifications, list(list(...)))
        invisible(NULL)
      },
      dirExistsFn = dir.exists,
      dirCreateFn = dir.create,
      fileExistsFn = file.exists,
      writeLinesFn = writeLines,
      timestampFn = function() as.POSIXct("2026-04-20 12:30:00", tz = "UTC"),
      catFn = function(...) invisible(NULL)
    )
  )
  expect_true(values$workflow_args_saved)
  expect_match(output$session_summary, "Parameters saved \\[OK\\]")

  fallback_values <- new.env(parent = emptyenv())
  fallback_values$workflow_args_saved <- FALSE
  expect_false(
    complete_workflow_save(
      output = new.env(parent = emptyenv()),
      values = fallback_values,
      projectDirs = list(proteomics = list(source_dir = source_dir, integration_dir = integration_dir, base_dir = fixture_dir)),
      workflowData = list(),
      experimentLabel = "summary-demo",
      description = "demo run",
      resolveFinalS4StateFn = function(...) stop("serialize boom"),
      assignFn = function(...) invisible(NULL),
      createWorkflowArgsFn = function(...) stop("should not build"),
      saveRDSFn = function(...) stop("should not save"),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      dirExistsFn = dir.exists,
      dirCreateFn = dir.create,
      fileExistsFn = file.exists,
      writeLinesFn = writeLines,
      timestampFn = function() "2026-04-20 12:31:00 UTC",
      catFn = function(...) invisible(NULL)
    )
  )
  expect_true(file.exists(file.path(source_dir, "study_parameters.txt")))

  prepared_export <- prepare_session_export(
    projectDirs = list(proteomics = list(source_dir = source_dir)),
    experimentLabel = "summary-demo",
    description = "demo run",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = TRUE,
    reportPath = "/tmp/report.html",
    exportDate = as.Date("2026-04-20"),
    timestamp = as.POSIXct("2026-04-20 12:32:00", tz = "UTC")
  )
  expect_match(prepared_export$sessionExportPath, "session_state_2026-04-20.RDS$", perl = TRUE)
  expect_identical(prepared_export$sessionState$experiment_label, "summary-demo")

  expect_true(
    complete_session_export(
      projectDirs = list(proteomics = list(source_dir = source_dir)),
      experimentLabel = "summary-demo",
      description = "demo run",
      workflowArgsSaved = TRUE,
      filesCopied = TRUE,
      reportGenerated = TRUE,
      reportPath = "/tmp/report.html",
      prepareExportFn = function(...) prepared_export,
      saveRDSFn = function(object, path) invisible(NULL),
      showNotificationFn = function(...) invisible(NULL),
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) stop("unexpected export error")
    )
  )
  expect_false(
    complete_session_export(
      projectDirs = list(proteomics = list(source_dir = source_dir)),
      experimentLabel = "summary-demo",
      description = "demo run",
      workflowArgsSaved = TRUE,
      filesCopied = TRUE,
      reportGenerated = TRUE,
      reportPath = "/tmp/report.html",
      prepareExportFn = function(...) prepared_export,
      saveRDSFn = function(...) stop("export boom"),
      showNotificationFn = function(...) invisible(NULL),
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL)
    )
  )
})

test_that("proteomics summary copy and github helpers preserve copy, error, and push behavior", {
  skip_if_not(
    exists("bootstrapProtSummaryCopyFallbackStudyParams", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("prepareProtSummaryCopyInputs", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("runProtSummaryPublicationCopy", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("handleProtSummaryPublicationCopyError", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("runProtSummaryGithubPush", envir = asNamespace("MultiScholaR"), inherits = FALSE) &&
      exists("completeProtSummaryGithubPush", envir = asNamespace("MultiScholaR"), inherits = FALSE)
  )

  bootstrap_copy_fallback <- get("bootstrapProtSummaryCopyFallbackStudyParams", envir = asNamespace("MultiScholaR"))
  prepare_copy_inputs <- get("prepareProtSummaryCopyInputs", envir = asNamespace("MultiScholaR"))
  run_publication_copy <- get("runProtSummaryPublicationCopy", envir = asNamespace("MultiScholaR"))
  handle_copy_error <- get("handleProtSummaryPublicationCopyError", envir = asNamespace("MultiScholaR"))
  run_github_push <- get("runProtSummaryGithubPush", envir = asNamespace("MultiScholaR"))
  complete_github_push <- get("completeProtSummaryGithubPush", envir = asNamespace("MultiScholaR"))

  fixture_dir <- tempfile("prot-summary-copy-")
  source_dir <- file.path(fixture_dir, "source")
  dir.create(source_dir, recursive = TRUE)

  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- FALSE
  values$files_copied <- FALSE
  expect_true(
    bootstrap_copy_fallback(
      values = values,
      projectDirs = list(proteomics = list(source_dir = source_dir)),
      experimentLabel = "summary-demo",
      description = "demo run",
      fileExistsFn = file.exists,
      writeLinesFn = writeLines,
      timestampFn = function() "2026-04-20 12:40:00 UTC",
      logErrorFn = function(...) stop("unexpected fallback error"),
      catFn = function(...) invisible(NULL)
    )
  )
  expect_true(values$workflow_args_saved)
  expect_true(file.exists(file.path(source_dir, "study_parameters.txt")))

  writeLines("sample\tgroup\nS1\tA\n", file.path(source_dir, "design_matrix.tab"))
  writeLines("contrast\nA_vs_B\n", file.path(source_dir, "contrasts_tbl.tab"))

  inputs <- prepare_copy_inputs(
    workflowData = list(),
    projectDirs = list(proteomics = list(source_dir = source_dir)),
    readTsvFn = function(path, show_col_types = FALSE) {
      if (grepl("design_matrix", path, fixed = TRUE)) {
        return(data.frame(sample = "S1", group = "A", stringsAsFactors = FALSE))
      }
      data.frame(contrast = "A_vs_B", stringsAsFactors = FALSE)
    },
    catFn = function(...) invisible(NULL)
  )
  expect_identical(inputs$designMatrix$sample, "S1")
  expect_identical(inputs$contrastsTbl$contrast, "A_vs_B")

  output <- new.env(parent = emptyenv())
  copy_args <- run_publication_copy(
    output = output,
    values = values,
    projectDirs = list(proteomics = list(source_dir = source_dir)),
    experimentLabel = "summary-demo",
    description = "demo run",
    contrastsTbl = inputs$contrastsTbl,
    designMatrix = inputs$designMatrix,
    existsFn = function(...) FALSE,
    assignFn = function(...) invisible(NULL),
    copyFn = function(...) invisible(NULL),
    renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
    showNotificationFn = function(...) invisible(NULL),
    timestampFn = function() as.POSIXct("2026-04-20 12:41:00", tz = "UTC"),
    catFn = function(...) invisible(NULL)
  )
  expect_true(values$files_copied)
  expect_identical(copy_args$omic_type, "proteomics")
  expect_match(output$session_summary, "Files copied \\[OK\\]")

  output <- new.env(parent = emptyenv())
  expect_false(
    handle_copy_error(
      output = output,
      error = simpleError("copy boom"),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL)
    )
  )
  expect_match(output$copy_status, "copy boom", fixed = TRUE)

  github_result <- run_github_push(
    projectDirs = list(proteomics = list(source_dir = source_dir)),
    experimentLabel = "summary-demo",
    githubOrg = "apaf",
    githubEmail = "summary@example.org",
    githubUsername = "summary-user",
    projectId = "MSR-001",
    optionsFn = function(...) invisible(NULL),
    pushFn = function(...) invisible(NULL)
  )
  expect_identical(github_result$pushArgs$project_id, "MSR-001")

  output <- new.env(parent = emptyenv())
  expect_true(
    complete_github_push(
      output = output,
      projectDirs = list(proteomics = list(source_dir = source_dir)),
      experimentLabel = "summary-demo",
      description = "demo run",
      githubOrg = "apaf",
      githubEmail = "summary@example.org",
      githubUsername = "summary-user",
      projectId = "MSR-001",
      pushGithubFn = function(...) invisible(NULL),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      timestampFn = function() as.POSIXct("2026-04-20 12:42:00", tz = "UTC"),
      logErrorFn = function(...) stop("unexpected github error")
    )
  )
  expect_match(output$session_summary, "GitHub pushed \\[OK\\]")

  expect_false(
    complete_github_push(
      output = new.env(parent = emptyenv()),
      projectDirs = list(proteomics = list(source_dir = source_dir)),
      experimentLabel = "summary-demo",
      description = "demo run",
      githubOrg = "apaf",
      githubEmail = "summary@example.org",
      githubUsername = "summary-user",
      projectId = "MSR-001",
      pushGithubFn = function(...) stop("github boom"),
      renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      timestampFn = function() as.POSIXct("2026-04-20 12:43:00", tz = "UTC"),
      logErrorFn = function(...) invisible(NULL)
    )
  )
})
