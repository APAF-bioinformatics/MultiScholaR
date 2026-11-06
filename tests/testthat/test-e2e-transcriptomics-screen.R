library(testthat)

# KIT-101 W5: the transcriptomics e2e lane `trx_screen` (tests/testdata/e2e/manifest.json), in the form of the other
# omics' lanes: the committed fixture (the KIT-040 trx_screen recipe, 12 samples, 300 genes, sha256 in the manifest),
# the browser from the tile through import, design, QC, normalisation, DA and the summary's RUVIII export, save, copy,
# report (the lane's stub template, seeded as the other lanes seed theirs) and session export. The KIT-097 journey
# (test-e2e-trx-journey.R) stays the trx_edge_12 proof.

.e2e_trx_screen <- function(step, input) e2e_input_id("transcriptomics", step, input)

test_that("E2E-KIT101 transcriptomics trx_screen lane reaches the summary, the RUVIII export and the report through the browser", {
  skip_if_no_e2e_browser()
  started <- Sys.time()
  lane <- read_e2e_manifest()[["trx_screen"]]
  design <- e2e_read_lane_design(lane)
  groups <- unlist(lane$groups)
  contrast <- strsplit(lane$expected_contrasts[[1L]], "_vs_", fixed = TRUE)[[1L]]
  case_id <- "E2E-KIT101-trx_screen"
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  base <- e2e_project_base_dir(project_root, case_id)
  results_summary <- file.path(base, "results_summary", "transcriptomics")
  source_dir <- file.path(base, "scripts", "transcriptomics")
  timeout <- 300000L
  driver <- e2e_new_app_driver(case_id = case_id, launch_surface = "run_app", timeout = 60000L, load_timeout = 180000L)
  on.exit(e2e_stop_driver(driver), add = TRUE)
  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-transcriptomics", timeout = 60000L)
    e2e_select_omics(driver, "transcriptomics", timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    # Import: the lane's count table through the browser file input.
    e2e_wait_for_selector(driver, "trx-tab-import", timeout = timeout)
    e2e_upload_file(driver, .e2e_trx_screen("import", "assay1_file_std"), e2e_lane_seed_path(lane))
    e2e_wait_for_idle(driver, timeout = timeout)
    e2e_click_testid(driver, "trx-import-process", timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "setup_import", "complete", timeout = timeout)
    # Design: the lane's groups, every sample assigned, the manifest's contrast.
    builder <- function(input) .e2e_trx_screen("design", paste0("builder-", input))
    e2e_switch_workflow_tab(driver, "transcriptomics", "design", timeout = timeout)
    e2e_wait_for_selector(driver, "trx-tab-design", timeout = timeout)
    e2e_set_input_and_idle(driver, builder("main_tabs"), "Factors", timeout = timeout)
    for (group in groups) {
      e2e_set_input_and_idle(driver, builder("new_factor"), group, timeout = timeout)
      e2e_click_input_id(driver, builder("add_factor"), timeout = timeout)
    }
    e2e_set_input_and_idle(driver, builder("main_tabs"), "Assign Metadata", timeout = timeout)
    available <- e2e_selectize_option_values(driver, builder("selected_runs"), timeout = timeout)
    for (group in groups) {
      e2e_set_selectize_values(driver, builder("selected_runs"), e2e_match_available_values(design$sample[design$group == group], available),
        timeout = timeout)
      e2e_wait_for_input_id(driver, builder("replicate_start"), timeout = timeout)
      e2e_set_selectize_values(driver, builder("factor1_select"), group, timeout = timeout, multiple = FALSE)
      e2e_click_input_id(driver, builder("assign_metadata"), timeout = timeout)
    }
    e2e_set_input_and_idle(driver, builder("main_tabs"), "Contrasts", timeout = timeout)
    e2e_set_selectize_values(driver, builder("contrast_group1"), contrast[[1L]], timeout = timeout, multiple = FALSE)
    e2e_set_selectize_values(driver, builder("contrast_group2"), contrast[[2L]], timeout = timeout, multiple = FALSE)
    e2e_click_input_id(driver, builder("add_contrast"), timeout = timeout)
    e2e_click_input_id(driver, builder("save_results"), timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "design_matrix", "complete", timeout = timeout)
    # QC: the expression filter and the finalise.
    e2e_switch_workflow_tab(driver, "transcriptomics", "qc", timeout = timeout)
    e2e_wait_for_selector(driver, "trx-tab-qc", timeout = timeout)
    e2e_click_testid(driver, "trx-qc-apply-filter", timeout = timeout)
    e2e_wait_for_r6_state(driver, "transcriptomics", "trx_expression_filtered", timeout = timeout)
    e2e_set_input_and_idle(driver, .e2e_trx_screen("qc", "trx_qc_tabs"), "Finalize QC", timeout = timeout)
    e2e_click_input_id(driver, .e2e_trx_screen("qc", "s4_finalize-finalize_qc"), timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "quality_control", "complete", timeout = timeout)
    # Normalisation and DA.
    e2e_switch_workflow_tab(driver, "transcriptomics", "norm", timeout = timeout)
    e2e_click_testid(driver, "trx-norm-run", timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "normalization", "complete", timeout = timeout)
    e2e_switch_workflow_tab(driver, "transcriptomics", "de", timeout = timeout)
    e2e_click_testid(driver, "trx-da-run", timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "differential_analysis", "complete", timeout = timeout)
    # Summary: the RUVIII export, save, copy, the report with the lane's template, the session export.
    e2e_switch_workflow_tab(driver, "transcriptomics", "summary", timeout = timeout)
    e2e_wait_for_selector(driver, "trx-tab-summary", timeout = timeout)
    e2e_click_testid(driver, "trx-summary-ruviii-export", timeout = timeout)
    e2e_wait_for_output_text(driver, .e2e_trx_screen("summary", "ruviii_status"), "M = group", timeout = timeout)
    e2e_set_input_and_idle(driver, .e2e_trx_screen("summary", "experiment_label"), case_id, timeout = timeout)
    e2e_set_input_and_idle(driver, .e2e_trx_screen("summary", "description"), "KIT-101 transcriptomics trx_screen lane", timeout = timeout)
    e2e_seed_report_template(lane, base)
    e2e_click_testid(driver, "trx-summary-save-args", timeout = timeout)
    e2e_click_testid(driver, "trx-summary-copy-publication", timeout = timeout)
    e2e_click_testid(driver, "trx-summary-generate-report", timeout = timeout)
    e2e_wait_for_selector(driver, "trx-summary-download-report", timeout = timeout)
    report <- e2e_get_download(driver, .e2e_trx_screen("summary", "download_report"), artifact_dir = e2e_case_artifact_dir(case_id),
      filename = paste0("downloaded-", lane$report_template))
    e2e_assert_file_nonempty(report)
    e2e_click_testid(driver, "trx-summary-export-state", timeout = timeout)
    e2e_assert_no_error_notifications(driver)
  })
  # Every state committed, the DA on the RUVs head with the edgeR QL model, the export's M, the lane's expected exports.
  record <- jsonlite::read_json(file.path(source_dir, "analysis_record.json"), simplifyVector = FALSE)
  operations <- vapply(record$steps, function(step) step$operation, "")
  expect_true(all(c("trx_expression_filter", "trx_qc_finalize", "trx_ruvs") %in% operations), info = paste(operations, collapse = " "))
  expect_identical(record$analyses$differential_abundance$source_state, "trx_ruvs")
  expect_identical(record$analyses$differential_abundance$model$model, "edgeR_QL")
  expect_identical(record$analyses$trx_ruviii_export$outcomes$ruviii_m_source, "group")
  exports <- list.files(base, recursive = TRUE)
  for (name in unlist(lane$expected_exports)) expect_true(any(basename(exports) == name), info = name)
  expect_gt(length(list.files(source_dir, pattern = "^session_state_.*\\.RDS$")), 0L)
  expect_lt(as.numeric(difftime(Sys.time(), started, units = "secs")), 1800)
})
