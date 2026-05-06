test_that("MCI-023.1 summary parameters roundtrip for lipid assay branches", {
  configs <- list(
    single_pos = module_ci_lipid_summary_config("single_pos", normalization_method = "none", ruv_mode = "skip", itsd_method = "skip", da_q_cutoff = 0.1),
    pos_gcms = module_ci_lipid_summary_config("pos_gcms", normalization_method = "scale", ruv_mode = "manual", itsd_method = "manual_assay_specific", da_q_cutoff = 0.01),
    pos_neg = module_ci_lipid_summary_config("pos_neg", normalization_method = "quantile", ruv_mode = "automatic", itsd_method = "regex_global", da_q_cutoff = 0.05)
  )

  for (branch_name in names(configs)) {
    config <- configs[[branch_name]]
    project_dirs <- module_ci_lipid_summary_paths()
    workflow <- module_ci_lipid_summary_workflow(config)
    output <- new.env(parent = emptyenv())
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- FALSE
    notifications <- character()

    result <- suppressMessages(suppressWarnings(runLipidSummarySaveWorkflowArgsObserverShell(
      inputValues = list(
        experiment_label = paste0("MCI023-", branch_name),
        description = paste("module-ci", branch_name)
      ),
      projectDirs = project_dirs,
      omicType = "lipidomics",
      workflowData = workflow,
      values = values,
      output = output,
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, ...) {
        notifications <<- c(notifications, message)
      },
      timeFn = function() as.POSIXct("2026-05-06 10:00:00", tz = "UTC"),
      catFn = function(...) invisible(NULL)
    )))

    expect_equal(result$status, "success", info = branch_name)
    expect_true(isTRUE(values$workflow_args_saved), info = branch_name)
    expect_true(file.exists(result$studyParamsFile), info = branch_name)
    expect_true(file.exists(result$s4Filepath), info = branch_name)

    saved_object <- readRDS(result$s4Filepath)
    expect_equal(saved_object@args$globalParameters$assay_layout, config$layout, info = branch_name)
    expect_equal(saved_object@args$globalParameters$workflow_type, "lipidomics", info = branch_name)
    expect_equal(saved_object@args$globalParameters$report_template, "lipidomics_report.rmd", info = branch_name)
    expect_equal(saved_object@args$normalisation_method, config$normalization_method, info = branch_name)
    expect_equal(saved_object@args$daAnalysisParameters$q_cutoff, config$da_q_cutoff, info = branch_name)
    expect_equal(saved_object@args$ITSDNormalization$method_type, config$itsd_method, info = branch_name)
    expect_true("Saved Integration S4 Object" %in% notifications, info = branch_name)

    study_text <- paste(readLines(result$studyParamsFile, warn = FALSE), collapse = "\n")
    expect_match(study_text, "Assay Information", info = branch_name)
    expect_match(study_text, paste(names(saved_object@lipid_data), collapse = "|"), info = branch_name)
    expect_match(study_text, "ITSD Normalization", info = branch_name)
    expect_match(study_text, "Between-Sample Normalization", info = branch_name)
    expect_match(study_text, "RUV-III Batch Correction", info = branch_name)
    expect_match(study_text, "KO_vs_WT", info = branch_name)
    expect_match(study_text, "groupKO-groupWT", info = branch_name)
  }
})

test_that("MCI-023.2 publication copy preserves required outputs and reports optional gaps", {
  skip_if_not_installed("openxlsx")

  complete_dirs <- module_ci_lipid_summary_paths()
  complete_config <- module_ci_lipid_summary_config("pos_neg", ruv_mode = "automatic")
  module_ci_lipid_summary_write_publication_inputs(complete_dirs, complete_config, include_optional = TRUE)
  complete_workflow <- module_ci_lipid_summary_workflow(complete_config)
  complete_values <- new.env(parent = emptyenv())
  complete_values$workflow_args_saved <- TRUE
  complete_values$files_copied <- FALSE
  complete_output <- new.env(parent = emptyenv())
  complete_notifications <- character()

  complete_result <- module_ci_lipid_summary_with_project_dirs(complete_dirs, {
    suppressWarnings(suppressMessages(runLipidSummaryCopyToPublicationObserverShell(
      inputValues = list(experiment_label = "MCI023-copy-complete", description = "complete artifacts"),
      projectDirs = complete_dirs,
      omicType = "lipidomics",
      workflowData = complete_workflow,
      values = complete_values,
      output = complete_output,
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, ...) {
        complete_notifications <<- c(complete_notifications, message)
      },
      catFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL),
      globalEnv = .GlobalEnv
    )))
  })

  expect_equal(complete_result$status, "success")
  expect_false("optionalMissing" %in% names(complete_result))
  expect_true(isTRUE(complete_values$files_copied))
  expect_true("Publication files copied" %in% complete_notifications)
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_tables", "DA_results_lipidomics.xlsx")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_tables", "ruv_normalised_results.RDS")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_tables", "da_lipidomics_num_sig_de_molecules.tab")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "QC_figures", "lcms_pos_pre_norm_pca.png")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_figures", "Heatmaps", "artifact.txt")))
  expect_true(file.exists(file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_figures", "NumSigDeMolecules", "artifact.txt")))

  copied_norm <- readr::read_tsv(
    file.path(complete_dirs$lipidomics$results_summary_dir, "Publication_tables", "RUV_normalised_results.tsv"),
    show_col_types = FALSE
  )
  expect_setequal(names(copied_norm), c("lipid_id", "assay", "Run", "normalised_intensity"))

  optional_dirs <- module_ci_lipid_summary_paths()
  module_ci_lipid_summary_write_publication_inputs(optional_dirs, complete_config, include_optional = FALSE)
  optional_values <- new.env(parent = emptyenv())
  optional_values$workflow_args_saved <- TRUE
  optional_values$files_copied <- FALSE
  optional_output <- new.env(parent = emptyenv())
  optional_notifications <- character()
  optional_result <- module_ci_lipid_summary_with_project_dirs(optional_dirs, {
    suppressWarnings(suppressMessages(runLipidSummaryCopyToPublicationObserverShell(
      inputValues = list(experiment_label = "MCI023-copy-optional", description = "missing optional artifacts"),
      projectDirs = optional_dirs,
      omicType = "lipidomics",
      workflowData = complete_workflow,
      values = optional_values,
      output = optional_output,
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(message, ...) {
        optional_notifications <<- c(optional_notifications, message)
      },
      catFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL),
      globalEnv = .GlobalEnv
    )))
  })

  expect_equal(optional_result$status, "success")
  expect_true(length(optional_result$optionalMissing) > 0L)
  expect_true("Publication files copied with optional missing artifacts" %in% optional_notifications)
  expect_true(file.exists(file.path(optional_dirs$lipidomics$results_summary_dir, "Publication_tables", "DA_results_lipidomics.xlsx")))

  required_dirs <- module_ci_lipid_summary_paths()
  required_values <- new.env(parent = emptyenv())
  required_values$workflow_args_saved <- TRUE
  required_values$files_copied <- FALSE
  required_output <- new.env(parent = emptyenv())
  required_result <- module_ci_lipid_summary_with_project_dirs(required_dirs, {
    suppressWarnings(suppressMessages(runLipidSummaryCopyToPublicationObserverShell(
      inputValues = list(experiment_label = "MCI023-copy-required", description = "missing required artifacts"),
      projectDirs = required_dirs,
      omicType = "lipidomics",
      workflowData = list(),
      values = required_values,
      output = required_output,
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      logErrorFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL),
      globalEnv = .GlobalEnv
    )))
  })
  expect_equal(required_result$status, "error")
  expect_match(required_result$errorMessage, "Required publication artifacts missing")
  expect_false(isTRUE(required_values$files_copied))
})

test_that("MCI-023.3 report rendering produces deterministic push-safe outputs", {
  project_dirs <- module_ci_lipid_summary_paths()
  module_ci_lipid_summary_write_template(project_dirs)
  output <- new.env(parent = emptyenv())
  values <- new.env(parent = emptyenv())
  values$files_copied <- TRUE
  values$report_generated <- FALSE
  values$report_path <- NULL

  rendered <- runLipidSummaryGenerateReportObserverShell(
    inputValues = list(experiment_label = "MCI023-render", description = "render stub"),
    projectDirs = project_dirs,
    omicType = "lipidomics",
    values = values,
    output = output,
    renderReportFn = function(omic_type, experiment_label, rmd_filename) {
      expect_equal(omic_type, "lipidomics")
      expect_equal(experiment_label, "MCI023-render")
      expect_equal(rmd_filename, "lipidomics_report.rmd")
      module_ci_lipid_summary_render_stub(project_dirs$lipidomics$results_summary_dir, rmd_filename)
    },
    reactiveFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    outputOptionsFn = function(...) invisible(NULL),
    downloadHandlerFn = function(filename, content) list(filename = filename, content = content),
    renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
    withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
    incProgressFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )

  expect_equal(rendered$status, "success")
  expect_true(isTRUE(values$report_generated))
  expect_true(file.exists(values$report_path))
  expect_equal(output$download_report$filename(), "lipidomics_report.html")
  expect_equal(rendered$templateFilename, "lipidomics_report.rmd")
  expect_match(rendered$reportTemplatePath, "lipidomics_report[.]rmd$")

  missing_output <- runLipidSummaryGenerateReportObserverShell(
    inputValues = list(experiment_label = "MCI023-render-missing", description = "missing render"),
    projectDirs = project_dirs,
    omicType = "lipidomics",
    values = list(files_copied = TRUE),
    output = new.env(parent = emptyenv()),
    renderReportFn = function(...) file.path(project_dirs$lipidomics$results_summary_dir, "missing.html"),
    withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
    incProgressFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL),
    catFn = function(...) invisible(NULL)
  )
  expect_equal(missing_output$status, "missing_output")
})

test_that("MCI-023.4 study/session parameter payload preserves lipid metadata", {
  config <- module_ci_lipid_summary_config("pos_neg", normalization_method = "cyclicloess", ruv_mode = "manual", da_q_cutoff = 0.01)
  workflow <- module_ci_lipid_summary_workflow(config)
  session_state <- prepareLipidSummarySessionState(
    inputValues = list(experiment_label = "MCI023-params", description = "payload"),
    projectDirs = module_ci_lipid_summary_paths(),
    omicType = "lipidomics",
    values = list(
      workflow_args_saved = TRUE,
      files_copied = TRUE,
      report_generated = FALSE,
      report_path = NULL
    ),
    timestamp = as.POSIXct("2026-05-06 10:00:00", tz = "UTC"),
    workflowData = workflow
  )

  expect_equal(session_state$workflow_type, "lipidomics")
  expect_equal(session_state$report_template, "lipidomics_report.rmd")
  expect_equal(session_state$assay_names, c("LCMS_Pos", "LCMS_Neg"))
  expect_equal(session_state$sample_count, 6L)
  expect_equal(session_state$feature_counts, c(LCMS_Pos = 6L, LCMS_Neg = 6L))
  expect_equal(session_state$contrast_count, 2L)
  expect_equal(session_state$parameter_payload$contrasts_tbl$friendly_names, c("KO_vs_WT", "RES_vs_WT"))
  expect_equal(session_state$parameter_payload$contrasts_tbl$contrasts, c("groupKO-groupWT", "groupRES-groupWT"))
  expect_equal(session_state$parameter_payload$normalization_ui_params$normalisation_method, "cyclicloess")
  expect_equal(session_state$parameter_payload$da_ui_params$q_cutoff, 0.01)
  expect_equal(session_state$parameter_payload$itsd_ui_params$method_type, "regex_global")
  expect_equal(session_state$parameter_payload$ruv_optimization_result$best_k, 2L)
  expect_equal(session_state$parameter_payload$s4_args$globalParameters$assay_layout, "pos_neg")
})

test_that("MCI-023.5 summary session-state export writes a reloadable schema-valid RDS", {
  config <- module_ci_lipid_summary_config("pos_gcms", ruv_mode = "automatic")
  workflow <- module_ci_lipid_summary_workflow(config)
  project_dirs <- module_ci_lipid_summary_paths()
  values <- list(
    workflow_args_saved = TRUE,
    files_copied = TRUE,
    report_generated = TRUE,
    report_path = file.path(project_dirs$lipidomics$results_summary_dir, "lipidomics_report.html")
  )

  export <- runLipidSummaryExportSessionObserverShell(
    inputValues = list(experiment_label = "MCI023-session", description = "session export"),
    projectDirs = project_dirs,
    omicType = "lipidomics",
    values = values,
    workflowData = workflow,
    sysDateFn = function() as.Date("2026-05-06"),
    sysTimeFn = function() as.POSIXct("2026-05-06 10:00:00", tz = "UTC"),
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL)
  )

  expect_equal(export$status, "success")
  expect_true(file.exists(export$sessionExportPath))
  loaded <- readRDS(export$sessionExportPath)
  expect_equal(loaded$omic_type, "lipidomics")
  expect_equal(loaded$workflow_type, "lipidomics")
  expect_equal(loaded$report_template, "lipidomics_report.rmd")
  expect_equal(loaded$assay_names, c("LCMS_Pos", "GCMS"))
  expect_equal(loaded$parameter_payload$da_ui_params$formula_string, "~ 0 + group + batch")
  expect_equal(loaded$parameter_payload$s4_args$globalParameters$report_template, "lipidomics_report.rmd")
})

test_that("MCI-023.6 artifact scorecard covers report-facing outputs", {
  skip_if_not_installed("openxlsx")
  project_dirs <- module_ci_lipid_summary_paths()
  config <- module_ci_lipid_summary_config("pos_neg", ruv_mode = "automatic")
  workflow <- module_ci_lipid_summary_workflow(config)
  module_ci_lipid_summary_write_publication_inputs(project_dirs, config, include_optional = TRUE)
  values <- new.env(parent = emptyenv())
  values$workflow_args_saved <- TRUE
  values$files_copied <- FALSE
  output <- new.env(parent = emptyenv())

  module_ci_lipid_summary_with_project_dirs(project_dirs, {
    suppressWarnings(suppressMessages(runLipidSummaryCopyToPublicationObserverShell(
      inputValues = list(experiment_label = "MCI023-scorecard", description = "scorecard"),
      projectDirs = project_dirs,
      omicType = "lipidomics",
      workflowData = workflow,
      values = values,
      output = output,
      renderTextFn = function(expr) eval(substitute(expr), envir = parent.frame()),
      withProgressFn = function(message = NULL, expr, ...) eval(substitute(expr), envir = parent.frame()),
      showNotificationFn = function(...) invisible(NULL),
      catFn = function(...) invisible(NULL),
      tracebackFn = function() invisible(NULL),
      globalEnv = .GlobalEnv
    )))
  })

  report_path <- module_ci_lipid_summary_render_stub(
    file.path(project_dirs$lipidomics$results_summary_dir, "Study_report")
  )
  session_export <- runLipidSummaryExportSessionObserverShell(
    inputValues = list(experiment_label = "MCI023-scorecard", description = "scorecard"),
    projectDirs = project_dirs,
    omicType = "lipidomics",
    values = list(
      workflow_args_saved = TRUE,
      files_copied = TRUE,
      report_generated = TRUE,
      report_path = report_path
    ),
    workflowData = workflow,
    showNotificationFn = function(...) invisible(NULL),
    logInfoFn = function(...) invisible(NULL)
  )

  scorecard <- module_ci_lipid_summary_scorecard(
    project_dirs = project_dirs,
    session_export_path = session_export$sessionExportPath,
    report_path = report_path
  )

  expect_setequal(
    scorecard$artifact,
    c(
      "study_parameters", "da_workbook", "normalized_tsv",
      "normalized_rds", "qc_figure", "num_sig",
      "session_state", "report_file"
    )
  )
  expect_true(all(scorecard$exists))
})
