library(testthat)

.e2e_run_protein_level_proteomics_lane <- function(lane_id) {
  skip_if_no_e2e_browser()

  manifest <- read_e2e_manifest()
  lane <- manifest[[lane_id]]
  case_id <- paste("E2E-005", lane_id, lane$import_tool, sep = "-")
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
    e2e_assert_step_statuses(
      driver,
      "proteomics",
      list(setup_import = "complete", design_matrix = "pending"),
      timeout = 60000L
    )
    e2e_assert_state_digest(
      driver,
      selected_omics = "proteomics",
      initialized_omics = "proteomics",
      required_project_dir_keys = "proteomics",
      workflow_types = list(proteomics = lane$workflow_type),
      timeout = 60000L
    )

    e2e_complete_prot_design_from_manifest(driver, lane, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "design_matrix", "complete", timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "quality_control", "complete", timeout = 120000L)
    e2e_assert_step_statuses(
      driver,
      "proteomics",
      list(
        design_matrix = "complete",
        quality_control = "complete",
        normalization = "pending"
      ),
      timeout = 60000L
    )

    e2e_run_prot_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "proteomics", "normalization", "complete", timeout = 240000L)
    session_artifacts <- e2e_assert_prot_filtered_session_artifacts(project_base)
    expect_identical(session_artifacts$session_data$workflow_type, lane$workflow_type)
    expect_equal(session_artifacts$session_data$final_sample_count, lane$sample_count)
    expect_gt(session_artifacts$session_data$final_protein_count, 0L)

    report_download <- e2e_run_prot_da_summary_report(
      driver,
      lane = lane,
      project_base_dir = project_base,
      experiment_label = case_id,
      case_id = paste0(case_id, "-report"),
      description = sprintf(
        "E2E-005 %s proteomics workflow via %s import",
        lane$workflow_type,
        lane$import_tool
      ),
      timeout = 240000L
    )
    e2e_wait_for_step_status(driver, "proteomics", "differential_expression", "complete", timeout = 240000L)
    e2e_assert_file_nonempty(report_download)
    e2e_assert_file_nonempty(file.path(paths$source_dir, "study_parameters.txt"))
    e2e_assert_file_nonempty(file.path(paths$source_dir, lane$report_template))

    if (!identical(lane$report_template, "DIANN_report.rmd")) {
      expect_false(
        file.exists(file.path(paths$source_dir, "DIANN_report.rmd")),
        info = sprintf("%s must not route through the DIA report template", lane_id)
      )
    }

    report_stem <- tools::file_path_sans_ext(lane$report_template)
    report_outputs <- list.files(
      paths$results_summary_dir,
      pattern = sprintf("^%s_proteomics_.*\\.(docx|html|pdf)$", report_stem),
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
}

test_that("E2E-005 Proteome Discoverer TMT reaches branch-specific report/export", {
  .e2e_run_protein_level_proteomics_lane("prot_tmt")
})

test_that("E2E-005 MaxQuant LFQ reaches branch-specific report/export", {
  .e2e_run_protein_level_proteomics_lane("prot_lfq")
})

test_that("E2E-005 FragPipe LFQ reaches branch-specific report/export", {
  .e2e_run_protein_level_proteomics_lane("prot_lfq_fragpipe")
})
