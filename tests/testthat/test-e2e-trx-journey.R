library(testthat)

# KIT-097: the one transcriptomics browser journey (F13: one attempt, 1,800 s, the kit library; any skip is a failure).
# The 12-sample corpus is the KIT-040 generator's trx_edge_12, regenerated and hash-checked (trxDeclReference). The app
# is driven from the tile through import, design, QC, normalisation, DA and the summary's RUVIII export, save, copy and
# report; the report is the package's transcriptomics template (no stub is seeded).

.e2e_trx <- function(step, input) e2e_input_id("transcriptomics", step, input)

.e2e_trx_project_paths <- function(project_dir, experiment_label) {
  base_dir <- e2e_project_base_dir(project_dir, experiment_label)
  results <- file.path(base_dir, "results", "transcriptomics")
  list(base_dir = base_dir, source_dir = file.path(base_dir, "scripts", "transcriptomics"),
    results_summary_dir = file.path(base_dir, "results_summary", "transcriptomics"), da_output_dir = file.path(results, "da_genes"),
    publication_graphs_dir = file.path(results, "publication_graphs"), gene_qc_dir = file.path(results, "gene_qc"),
    integration_dir = file.path(base_dir, "integration"))
}

.e2e_trx_import <- function(driver, counts_file, timeout) {
  e2e_wait_for_selector(driver, "trx-tab-import", timeout = timeout)
  e2e_upload_file(driver, .e2e_trx("import", "assay1_file_std"), counts_file)
  e2e_wait_for_idle(driver, timeout = timeout)
  e2e_click_testid(driver, "trx-import-process", timeout = timeout)
}

.e2e_trx_design <- function(driver, design, timeout) {
  builder <- function(input) .e2e_trx("design", paste0("builder-", input))
  groups <- c("RPMI", "Sera")
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
    e2e_set_selectize_values(driver, builder("selected_runs"), e2e_match_available_values(design$Run[design$group == group], available), timeout = timeout)
    e2e_wait_for_input_id(driver, builder("replicate_start"), timeout = timeout)
    e2e_set_selectize_values(driver, builder("factor1_select"), group, timeout = timeout, multiple = FALSE)
    e2e_click_input_id(driver, builder("assign_metadata"), timeout = timeout)
  }
  e2e_set_input_and_idle(driver, builder("main_tabs"), "Contrasts", timeout = timeout)
  e2e_set_selectize_values(driver, builder("contrast_group1"), "Sera", timeout = timeout, multiple = FALSE)
  e2e_set_selectize_values(driver, builder("contrast_group2"), "RPMI", timeout = timeout, multiple = FALSE)
  e2e_click_input_id(driver, builder("add_contrast"), timeout = timeout)
  e2e_click_input_id(driver, builder("save_results"), timeout = timeout)
}

.e2e_trx_qc <- function(driver, timeout) {
  e2e_switch_workflow_tab(driver, "transcriptomics", "qc", timeout = timeout)
  e2e_wait_for_selector(driver, "trx-tab-qc", timeout = timeout)
  e2e_click_testid(driver, "trx-qc-apply-filter", timeout = timeout)
  e2e_wait_for_r6_state(driver, "transcriptomics", "trx_expression_filtered", timeout = timeout)
  e2e_set_input_and_idle(driver, .e2e_trx("qc", "trx_qc_tabs"), "Finalize QC", timeout = timeout)
  e2e_click_input_id(driver, .e2e_trx("qc", "s4_finalize-finalize_qc"), timeout = timeout)
}

.e2e_trx_summary <- function(driver, experiment_label, case_id, timeout) {
  e2e_switch_workflow_tab(driver, "transcriptomics", "summary", timeout = timeout)
  e2e_wait_for_selector(driver, "trx-tab-summary", timeout = timeout)
  e2e_click_testid(driver, "trx-summary-ruviii-export", timeout = timeout)
  e2e_wait_for_output_text(driver, .e2e_trx("summary", "ruviii_status"), "M = group", timeout = timeout)
  e2e_set_input_and_idle(driver, .e2e_trx("summary", "experiment_label"), experiment_label, timeout = timeout)
  e2e_set_input_and_idle(driver, .e2e_trx("summary", "description"), "KIT-097 transcriptomics browser journey", timeout = timeout)
  e2e_click_testid(driver, "trx-summary-save-args", timeout = timeout)
  e2e_click_testid(driver, "trx-summary-copy-publication", timeout = timeout)
  e2e_click_testid(driver, "trx-summary-generate-report", timeout = timeout)
  e2e_wait_for_selector(driver, "trx-summary-download-report", timeout = timeout)
  report <- e2e_get_download(driver, .e2e_trx("summary", "download_report"), artifact_dir = e2e_case_artifact_dir(case_id),
    filename = "transcriptomics_report.docx")
  e2e_click_testid(driver, "trx-summary-export-state", timeout = timeout)
  report
}

test_that("E2E-KIT097 transcriptomics: the GUI journey commits every state, runs the RUVIII export and renders the report", {
  skip_if_no_e2e_browser()
  skipWithoutTrxDeclReference()
  started <- Sys.time()
  reference <- trxDeclReference("trx_edge_12")
  design <- utils::read.delim(file.path(reference$dir, "design_matrix.tab"), stringsAsFactors = FALSE)
  case_id <- "E2E-KIT097-trx_edge_12"
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  paths <- .e2e_trx_project_paths(project_root, case_id)
  timeout <- 300000L
  driver <- e2e_new_app_driver(case_id = case_id, launch_surface = "run_app", timeout = 60000L, load_timeout = 180000L)
  on.exit(e2e_stop_driver(driver), add = TRUE)
  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-transcriptomics", timeout = 60000L)
    e2e_select_omics(driver, "transcriptomics", timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    .e2e_trx_import(driver, file.path(reference$dir, "counts.tsv"), timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "setup_import", "complete", timeout = timeout)
    .e2e_trx_design(driver, design, timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "design_matrix", "complete", timeout = timeout)
    .e2e_trx_qc(driver, timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "quality_control", "complete", timeout = timeout)
    e2e_switch_workflow_tab(driver, "transcriptomics", "norm", timeout = timeout)
    e2e_click_testid(driver, "trx-norm-run", timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "normalization", "complete", timeout = timeout)
    e2e_switch_workflow_tab(driver, "transcriptomics", "de", timeout = timeout)
    e2e_click_testid(driver, "trx-da-run", timeout = timeout)
    e2e_wait_for_step_status(driver, "transcriptomics", "differential_analysis", "complete", timeout = timeout)
    report <- .e2e_trx_summary(driver, case_id, case_id, timeout)
    e2e_assert_file_nonempty(report)
    e2e_assert_no_error_notifications(driver)
  })
  # Every state committed, the DA on the RUVs head, the export's M, the model as it ran.
  record <- jsonlite::read_json(file.path(paths$source_dir, "analysis_record.json"), simplifyVector = FALSE)
  operations <- vapply(record$steps, function(step) step$operation, "")
  expect_true(all(c("trx_expression_filter", "trx_qc_finalize", "trx_ruvs") %in% operations), info = paste(operations, collapse = " "))
  expect_identical(record$analyses$differential_abundance$source_state, "trx_ruvs")
  expect_identical(record$analyses$differential_abundance$model$model, "edgeR_QL")
  expect_false(isTRUE(record$analyses$differential_abundance$model$treat_applied))
  expect_identical(record$analyses$trx_ruviii_export$outcomes$ruviii_m_source, "group")
  expect_true(file.exists(file.path(paths$integration_dir, "transcriptomics_normalized_data.tsv")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "Publication_tables", "DA_results_transcriptomics.xlsx")))
  expect_true(file.exists(file.path(paths$results_summary_dir, "QC_figures", "ruviii_before_and_after.png")))
  # The rendered report is kept beside the journey's artifacts for the wording grep (F05).
  expect_lt(as.numeric(difftime(Sys.time(), started, units = "secs")), 1800)
})
