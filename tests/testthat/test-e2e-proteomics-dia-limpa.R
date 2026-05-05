library(testthat)

test_that("E2E-006 DIA limpa workflow preserves imputation metadata through report generation", {
  skip_if_no_e2e_browser()

  case_id <- "E2E-006-proteomics-dia-limpa"
  lane <- read_e2e_manifest()[["prot_dia_limpa"]]
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  project_base <- e2e_project_base_dir(project_root, case_id)
  paths <- e2e_proteomics_project_paths(project_root, case_id)

  driver <- e2e_new_app_driver(
    case_id = case_id,
    launch_surface = "run_app",
    timeout = 60000L,
    load_timeout = 180000L
  )
  on.exit(e2e_stop_driver(driver), add = TRUE)

  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-proteomics", timeout = 60000L)
    e2e_select_omics(driver, "proteomics", assert_digest = TRUE, timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    e2e_wait_for_selector(driver, "prot-tab-setup", timeout = 60000L)
    e2e_assert_project_dirs_exist(paths)

    e2e_upload_lane_inputs(driver, lane)
    e2e_click_process_import(driver, "proteomics")
    e2e_wait_for_step_status(driver, "proteomics", "setup_import", "complete", timeout = 120000L)
    e2e_assert_state_digest(
      driver,
      selected_omics = "proteomics",
      initialized_omics = "proteomics",
      required_project_dir_keys = "proteomics",
      workflow_types = list(proteomics = "DIA"),
      timeout = 60000L
    )

    e2e_complete_prot_design_from_manifest(driver, lane, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "design_matrix", "complete", timeout = 120000L)

    e2e_run_prot_dia_limpa_qc(driver, timeout = 180000L)
    e2e_wait_for_step_status(driver, "proteomics", "quality_control", "complete", timeout = 180000L)

    e2e_run_prot_dia_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "proteomics", "normalization", "complete", timeout = 240000L)
    session_artifacts <- e2e_assert_prot_filtered_session_artifacts(project_base)
    session_data <- session_artifacts$session_data
    expect_identical(session_data$workflow_type, "DIA")
    expect_identical(session_data$report_template, "DIANN_limpa_report.rmd")
    expect_true(isTRUE(session_data$limpa_applied))
    expect_true(isTRUE(session_data$use_limpa))
    expect_equal(session_data$final_sample_count, lane$sample_count)
    expect_gt(session_data$final_protein_count, 0L)

    limpa_args <- session_data$current_s4_object@args
    expect_true(isTRUE(limpa_args$globalParameters$use_limpa))
    expect_identical(limpa_args$globalParameters$report_template, "DIANN_limpa_report.rmd")
    expect_false(is.null(limpa_args$proteinMissingValueImputationLimpa))
    expect_false(is.null(limpa_args$limpa_dpc_quant_results))
    expect_identical(limpa_args$limpa_dpc_quant_results$dpc_method, "limpa_dpc_quant")
    expect_match(
      limpa_args$limpa_dpc_quant_results$quantification_method,
      "^limpa_dpc_quant",
      perl = TRUE
    )

    e2e_assert_file_nonempty(file.path(paths$source_dir, "limpa_parameters.RDS"))
    e2e_assert_file_nonempty(file.path(paths$source_dir, "limpa_dpc_quant_results.RDS"))
    filtered_summary <- paste(readLines(session_artifacts$summary, warn = FALSE), collapse = "\n")
    expect_match(filtered_summary, "Report Template: DIANN_limpa_report.rmd", fixed = TRUE)
    expect_match(filtered_summary, "limpa DPC-Quant: yes", fixed = TRUE)

    report_download <- e2e_run_prot_da_summary_report(
      driver,
      lane = lane,
      project_base_dir = project_base,
      experiment_label = case_id,
      case_id = paste0(case_id, "-report"),
      description = "E2E-006 DIA limpa GUI workflow with DPC-Quant imputation branch",
      timeout = 240000L
    )
    e2e_wait_for_step_status(driver, "proteomics", "differential_expression", "complete", timeout = 240000L)
    e2e_assert_file_nonempty(report_download)
    e2e_assert_file_nonempty(file.path(paths$source_dir, "study_parameters.txt"))
    e2e_assert_file_nonempty(file.path(paths$source_dir, lane$report_template))
    expect_false(
      file.exists(file.path(paths$source_dir, "DIANN_report.rmd")),
      info = "DIA limpa workflow must not route through the canonical DIA report template"
    )

    study_parameters <- paste(readLines(file.path(paths$source_dir, "study_parameters.txt"), warn = FALSE), collapse = "\n")
    expect_match(study_parameters, "DIANN_limpa_report.rmd", fixed = TRUE)
    expect_match(study_parameters, "proteinMissingValueImputationLimpa", fixed = TRUE)
    expect_match(study_parameters, "limpa_dpc_quant", fixed = TRUE)

    report_outputs <- list.files(
      paths$results_summary_dir,
      pattern = "^DIANN_limpa_report_proteomics_.*\\.(docx|html|pdf)$",
      full.names = TRUE
    )
    expect_gt(length(report_outputs), 0L)
    invisible(lapply(report_outputs, e2e_assert_file_nonempty))

    session_state_exports <- list.files(
      paths$source_dir,
      pattern = "^session_state_.*\\.RDS$",
      full.names = TRUE
    )
    expect_gt(length(session_state_exports), 0L)
    invisible(lapply(session_state_exports, e2e_assert_file_nonempty))

    e2e_assert_no_console_errors(driver)
    e2e_assert_no_error_notifications(driver)
  })
})
