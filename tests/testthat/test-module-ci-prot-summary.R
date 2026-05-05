test_that("MCI-011.1 workflow parameters roundtrip for every proteomics summary branch", {
  branches <- list(
    DIA = module_ci_prot_summary_branch_config("DIA", enrichment = TRUE, ruv = "run", da_q_cutoff = 0.05),
    DIA_limpa = module_ci_prot_summary_branch_config("DIA_limpa", enrichment = TRUE, ruv = "run", da_q_cutoff = 0.01),
    TMT = module_ci_prot_summary_branch_config("TMT", enrichment = FALSE, ruv = "skip", da_q_cutoff = 0.1),
    MaxQuant_LFQ = module_ci_prot_summary_branch_config("MaxQuant_LFQ", enrichment = TRUE, ruv = "skip", da_q_cutoff = 0.2),
    FragPipe_LFQ = module_ci_prot_summary_branch_config("FragPipe_LFQ", enrichment = FALSE, ruv = "run", da_q_cutoff = 0.05)
  )

  for (branch_name in names(branches)) {
    config <- branches[[branch_name]]
    project_dirs <- module_ci_prot_summary_paths()
    workflow <- module_ci_prot_summary_workflow(config)
    output <- new.env(parent = emptyenv())
    values <- new.env(parent = emptyenv())
    notifications <- character()

    ok <- suppressMessages(suppressWarnings(completeProtSummaryWorkflowArgsSave(
      output = output,
      values = values,
      projectDirs = project_dirs,
      workflowData = workflow,
      experimentLabel = paste0("MCI011-", branch_name),
      description = paste("module-ci", branch_name),
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, ...) {
        notifications <<- c(notifications, message)
      },
      timestampFn = function() as.POSIXct("2026-05-05 10:00:00", tz = "UTC"),
      catFn = function(...) invisible(NULL)
    )))

    expect_true(isTRUE(ok), info = branch_name)
    expect_true(isTRUE(values$workflow_args_saved), info = branch_name)
    study_file <- file.path(project_dirs$proteomics$source_dir, "study_parameters.txt")
    integration_file <- file.path(
      project_dirs$proteomics$integration_dir,
      paste0("proteomics_MCI011-", branch_name, "_final_s4.RDS")
    )
    expect_true(file.exists(study_file), info = branch_name)
    expect_true(file.exists(integration_file), info = branch_name)
    saved_object <- readRDS(integration_file)
    expect_equal(saved_object@args$globalParameters$workflow_type, config$workflow_type, info = branch_name)
    expect_equal(saved_object@args$globalParameters$report_template, config$report_template, info = branch_name)
    expect_equal(saved_object@args$deAnalysisParameters$q_cutoff, config$da_q_cutoff, info = branch_name)
    expect_true("Saved Integration S4 Object" %in% notifications, info = branch_name)
  }
})

test_that("MCI-011.2 publication copy preserves required outputs and tolerates missing optional artifacts", {
  skip_if_not_installed("openxlsx")
  complete_dirs <- module_ci_prot_summary_paths()
  module_ci_prot_summary_write_publication_inputs(complete_dirs, include_optional = TRUE)
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  workflow <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
  notifications <- character()

  copy_args <- suppressMessages(suppressWarnings(runProtSummaryPublicationCopy(
    output = output,
    values = values,
    projectDirs = complete_dirs,
    experimentLabel = "MCI011-copy",
    description = "complete artifacts",
    contrastsTbl = workflow$contrasts_tbl,
    designMatrix = workflow$design_matrix,
    existsFn = function(...) FALSE,
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(message, ...) {
      notifications <<- c(notifications, message)
    },
    catFn = function(...) invisible(NULL)
  )))

  expect_true(isTRUE(values$files_copied))
  expect_equal(copy_args$omic_type, "proteomics")
  expect_true("Publication files copied" %in% notifications)
  expect_true(file.exists(file.path(complete_dirs$proteomics$results_summary_dir, "Publication_tables", "DA_results_proteomics.xlsx")))
  expect_true(file.exists(file.path(complete_dirs$proteomics$results_summary_dir, "Publication_tables", "Pathway_enrichment_results_proteomics.xlsx")))
  expect_true(file.exists(file.path(complete_dirs$proteomics$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv")))
  expect_true(file.exists(file.path(complete_dirs$proteomics$results_summary_dir, "Publication_tables", "ruv_normalised_results.RDS")))
  copied_ruv <- readr::read_tsv(
    file.path(complete_dirs$proteomics$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv"),
    show_col_types = FALSE
  )
  expect_setequal(names(copied_ruv), c("Protein.Group", "Run", "normalised_intensity"))

  missing_optional_dirs <- module_ci_prot_summary_paths()
  module_ci_prot_summary_write_publication_inputs(missing_optional_dirs, include_optional = FALSE, ruv_applied = FALSE)
  missing_values <- new.env(parent = emptyenv())
  missing_output <- new.env(parent = emptyenv())
  expect_no_error(suppressMessages(suppressWarnings(runProtSummaryPublicationCopy(
    output = missing_output,
    values = missing_values,
    projectDirs = missing_optional_dirs,
    experimentLabel = "MCI011-missing-optional",
    description = "missing optional artifacts",
    contrastsTbl = workflow$contrasts_tbl,
    designMatrix = workflow$design_matrix,
    existsFn = function(...) FALSE,
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  ))))
  expect_true(isTRUE(missing_values$files_copied))
  expect_true(file.exists(file.path(missing_optional_dirs$proteomics$results_summary_dir, "Publication_tables", "Pathway_enrichment_results_proteomics.xlsx")))
  expect_true(file.exists(file.path(missing_optional_dirs$proteomics$results_summary_dir, "Publication_tables", "normalised_results.tsv")))
  expect_true(file.exists(file.path(missing_optional_dirs$proteomics$results_summary_dir, "Publication_tables", "normalised_results.RDS")))
})

test_that("MCI-011.3 report-template selection is deterministic and branch-correct", {
  project_dirs <- module_ci_prot_summary_paths()
  module_ci_prot_summary_write_template_set(project_dirs)
  status <- buildProtSummaryTemplateStatus(project_dirs, "proteomics")
  expect_match(status, "DIA-NN \\[OK\\]")
  expect_match(status, "DIA-NN limpa \\[OK\\]")
  expect_match(status, "TMT \\[OK\\]")
  expect_match(status, "LFQ \\[OK\\]")

  expected <- c(
    DIA = "DIANN_report.rmd",
    DIA_limpa = "DIANN_limpa_report.rmd",
    TMT = "TMT_report.rmd",
    MaxQuant_LFQ = "LFQ_report.rmd",
    FragPipe_LFQ = "LFQ_report.rmd"
  )

  for (branch_name in names(expected)) {
    workflow <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config(branch_name))
    resolved <- resolveProtSummaryReportTemplate(workflow, catFn = function(...) invisible(NULL))
    expect_equal(resolved$templateFilename, expected[[branch_name]], info = branch_name)
    expect_false(isTRUE(resolved$staleTemplateIgnored), info = branch_name)
  }

  stale_tmt <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("TMT"))
  stale_tmt$config_list$globalParameters$report_template <- "DIANN_report.rmd"
  resolved_stale <- resolveProtSummaryReportTemplate(stale_tmt, catFn = function(...) invisible(NULL))
  expect_equal(resolved_stale$templateFilename, "TMT_report.rmd")
  expect_true(isTRUE(resolved_stale$staleTemplateIgnored))

  stale_dia <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
  stale_dia$config_list$globalParameters$report_template <- "DIANN_limpa_report.rmd"
  stale_dia$config_list$globalParameters$use_limpa <- FALSE
  resolved_stale_dia <- resolveProtSummaryReportTemplate(stale_dia, catFn = function(...) invisible(NULL))
  expect_equal(resolved_stale_dia$templateFilename, "DIANN_report.rmd")
  expect_true(isTRUE(resolved_stale_dia$staleTemplateIgnored))

  stale_limpa <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA_limpa"))
  stale_limpa$config_list$globalParameters$report_template <- "DIANN_report.rmd"
  stale_limpa$config_list$globalParameters$use_limpa <- TRUE
  resolved_stale_limpa <- resolveProtSummaryReportTemplate(stale_limpa, catFn = function(...) invisible(NULL))
  expect_equal(resolved_stale_limpa$templateFilename, "DIANN_limpa_report.rmd")
  expect_true(isTRUE(resolved_stale_limpa$staleTemplateIgnored))
})

test_that("MCI-011.4 rendered-report stubs activate every proteomics template", {
  templates <- c("DIANN_report.rmd", "DIANN_limpa_report.rmd", "TMT_report.rmd", "LFQ_report.rmd")

  for (template in templates) {
    output <- new.env(parent = emptyenv())
    values <- new.env(parent = emptyenv())
    render_dir <- tempfile("module-ci-prot-summary-render-")

    rendered <- runProtSummaryReportGeneration(
      output = output,
      values = values,
      experimentLabel = paste0("MCI011-", template),
      description = "deterministic render stub",
      templateFilename = template,
      renderReportAvailableFn = function() TRUE,
      renderReportFn = function(omic_type, experiment_label, rmd_filename) {
        expect_equal(omic_type, "proteomics")
        expect_equal(rmd_filename, template)
        module_ci_prot_summary_render_stub(render_dir, rmd_filename)
      },
      activateReportFn = function(output, values, renderedPath, experimentLabel, description) {
        activateProtSummaryRenderedReport(
          output = output,
          values = values,
          renderedPath = renderedPath,
          experimentLabel = experimentLabel,
          description = description,
          reactiveFn = function(expr) eval(substitute(expr), envir = parent.frame()),
          outputOptionsFn = function(...) invisible(NULL),
          downloadHandlerFn = function(filename, content) list(filename = filename, content = content),
          renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
          showNotificationFn = function(...) invisible(NULL)
        )
      },
      showNotificationFn = function(...) invisible(NULL),
      logInfoFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      printFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL)
    )

    expect_true(isTRUE(rendered), info = template)
    expect_true(isTRUE(values$report_generated), info = template)
    expect_true(file.exists(values$report_path), info = template)
    expect_equal(output$download_report$filename(), sub("\\.rmd$", ".html", template, ignore.case = TRUE), info = template)
  }
})

test_that("MCI-011.5 session-state export preserves workflow metadata and parameter payload", {
  config <- module_ci_prot_summary_branch_config("DIA_limpa", enrichment = TRUE, ruv = "run", da_q_cutoff = 0.01)
  workflow <- module_ci_prot_summary_workflow(config)
  project_dirs <- module_ci_prot_summary_paths()
  notifications <- character()
  logs <- character()

  ok <- completeProtSummarySessionStateExport(
    projectDirs = project_dirs,
    experimentLabel = "MCI011-session",
    description = "session export",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = TRUE,
    reportPath = "/tmp/MCI011-session.html",
    workflowData = workflow,
    showNotificationFn = function(message, ...) {
      notifications <<- c(notifications, message)
    },
    logInfoFn = function(message) {
      logs <<- c(logs, message)
    },
    logErrorFn = function(...) stop("export should not log errors")
  )

  expect_true(isTRUE(ok))
  export_path <- file.path(project_dirs$proteomics$source_dir, paste0("session_state_", Sys.Date(), ".RDS"))
  expect_true(file.exists(export_path))
  loaded <- readRDS(export_path)
  expect_equal(loaded$omic_type, "proteomics")
  expect_equal(loaded$workflow_type, "DIA")
  expect_equal(loaded$report_template, "DIANN_limpa_report.rmd")
  expect_equal(loaded$parameter_payload$globalParameters$workflow_type, "DIA")
  expect_equal(loaded$parameter_payload$da_ui_params$q_cutoff, 0.01)
  expect_equal(loaded$parameter_payload$enrichment_ui_params$database_source, "gprofiler2")
  expect_equal(loaded$parameter_payload$ruv_optimization_result$best_k, 2L)
  expect_true(any(grepl("Session state exported to:", notifications, fixed = TRUE)))
  expect_true(any(grepl("Session state exported to:", logs, fixed = TRUE)))
})

test_that("MCI-011.6 artifact scorecard entries cover report-facing outputs", {
  skip_if_not_installed("openxlsx")
  project_dirs <- module_ci_prot_summary_paths()
  module_ci_prot_summary_write_publication_inputs(project_dirs, include_optional = TRUE)
  workflow <- module_ci_prot_summary_workflow(module_ci_prot_summary_branch_config("DIA"))
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())

  suppressMessages(suppressWarnings(runProtSummaryPublicationCopy(
    output = output,
    values = values,
    projectDirs = project_dirs,
    experimentLabel = "MCI011-scorecard",
    description = "scorecard",
    contrastsTbl = workflow$contrasts_tbl,
    designMatrix = workflow$design_matrix,
    existsFn = function(...) FALSE,
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    showNotificationFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )))

  session_export <- prepareProtSummarySessionStateExport(
    projectDirs = project_dirs,
    experimentLabel = "MCI011-scorecard",
    description = "scorecard",
    workflowArgsSaved = TRUE,
    filesCopied = TRUE,
    reportGenerated = TRUE,
    reportPath = file.path(project_dirs$proteomics$results_summary_dir, "Study_report", "report.html"),
    workflowData = workflow
  )
  saveRDS(session_export$sessionState, session_export$sessionExportPath)
  report_path <- module_ci_prot_summary_render_stub(
    file.path(project_dirs$proteomics$results_summary_dir, "Study_report"),
    "DIANN_report.rmd"
  )

  scorecard <- module_ci_prot_summary_scorecard(
    project_dirs = project_dirs,
    session_export_path = session_export$sessionExportPath,
    report_path = report_path
  )

  expect_setequal(
    scorecard$artifact,
    c("study_parameters", "da_workbook", "enrichment_workbook", "session_state", "report_file")
  )
  expect_true(all(scorecard$exists))
})
