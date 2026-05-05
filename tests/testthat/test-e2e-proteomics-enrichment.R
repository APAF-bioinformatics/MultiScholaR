library(testthat)

.e2e_assert_enrichment_archive <- function(archive_path, backend) {
  e2e_assert_file_nonempty(archive_path)
  listing <- utils::unzip(archive_path, list = TRUE)
  expected_result <- if (identical(backend, "gprofiler2")) {
    "gprofiler2_results.tsv"
  } else {
    "clusterProfileR_results.tsv"
  }
  unexpected_result <- if (identical(backend, "gprofiler2")) {
    "clusterProfileR_results.tsv"
  } else {
    "gprofiler2_results.tsv"
  }

  expect_true(expected_result %in% listing$Name)
  expect_false(unexpected_result %in% listing$Name)
  expect_true("enrichment_analysis_summary.txt" %in% listing$Name)

  extract_dir <- tempfile("e2e-enrichment-archive-")
  dir.create(extract_dir)
  utils::unzip(archive_path, exdir = extract_dir)
  summary_text <- paste(readLines(file.path(extract_dir, "enrichment_analysis_summary.txt"), warn = FALSE), collapse = "\n")
  expect_match(summary_text, sprintf("Analysis Method: %s", backend), fixed = TRUE)

  result_table <- read.delim(file.path(extract_dir, expected_result), stringsAsFactors = FALSE)
  expect_true("analysis_method" %in% names(result_table))
  expect_true(all(result_table$analysis_method == backend))

  invisible(result_table)
}

.e2e_assert_enrichment_publication_workbook <- function(paths, backend) {
  workbook <- file.path(
    paths$results_summary_dir,
    "Publication_tables",
    "Pathway_enrichment_results_proteomics.xlsx"
  )
  e2e_assert_file_nonempty(workbook)

  if (requireNamespace("openxlsx", quietly = TRUE)) {
    index <- openxlsx::read.xlsx(workbook, sheet = "Enrichment_Index")
    expect_gt(nrow(index), 0L)
    first_sheet <- index$Sheet[[1]]
    data <- openxlsx::read.xlsx(workbook, sheet = first_sheet)
    expect_true("analysis_method" %in% names(data))
    expect_true(all(data$analysis_method == backend))
  }

  invisible(workbook)
}

.e2e_run_prot_enrichment_backend_lane <- function(backend) {
  skip_if_no_e2e_browser()

  lane <- read_e2e_manifest()[["prot_dia"]]
  case_id <- paste("E2E-007-proteomics-enrichment", backend, sep = "-")
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  project_base <- e2e_project_base_dir(project_root, case_id)
  paths <- e2e_proteomics_project_paths(project_root, case_id)
  organism_taxid <- if (identical(backend, "gprofiler2")) "9606" else "999999"

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
    e2e_complete_prot_design_from_manifest(driver, lane, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "design_matrix", "complete", timeout = 120000L)

    e2e_run_prot_dia_qc(driver, timeout = 180000L)
    e2e_wait_for_step_status(driver, "proteomics", "quality_control", "complete", timeout = 180000L)

    e2e_run_prot_dia_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "proteomics", "normalization", "complete", timeout = 240000L)
    session_artifacts <- e2e_assert_prot_filtered_session_artifacts(project_base)
    expect_identical(session_artifacts$session_data$workflow_type, "DIA")

    e2e_run_prot_da(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "proteomics", "differential_expression", "complete", timeout = 240000L)

    enrich_result <- e2e_run_prot_enrichment_backend(
      driver,
      backend = backend,
      organism_taxid = organism_taxid,
      expected_contrast = lane$expected_contrasts[[1]],
      timeout = 240000L
    )
    e2e_wait_for_step_status(driver, "proteomics", "enrichment_analysis", "complete", timeout = 240000L)
    expect_match(enrich_result$status_text, sprintf("Method: %s", backend), fixed = TRUE)

    pathway_files <- list.files(
      paths$pathway_dir,
      pattern = "_enrichment_results\\.tsv$",
      full.names = TRUE
    )
    expect_gt(length(pathway_files), 0L)
    pathway_table <- read.delim(pathway_files[[1]], stringsAsFactors = FALSE)
    expect_true("analysis_method" %in% names(pathway_table))
    expect_true(all(pathway_table$analysis_method == backend))

    archive <- e2e_get_download(
      driver,
      e2e_input_id("proteomics", "enrich", "download_enrichment_results"),
      artifact_dir = e2e_case_artifact_dir(case_id),
      filename = paste0("downloaded-enrichment-", backend, ".zip")
    )
    .e2e_assert_enrichment_archive(archive, backend)

    report_download <- e2e_run_prot_summary_report(
      driver,
      lane = lane,
      project_base_dir = project_base,
      experiment_label = case_id,
      case_id = paste0(case_id, "-report"),
      description = sprintf("E2E-007 proteomics enrichment branch via %s", backend),
      timeout = 240000L
    )
    e2e_assert_file_nonempty(report_download)

    study_parameters <- paste(readLines(file.path(paths$source_dir, "study_parameters.txt"), warn = FALSE), collapse = "\n")
    expect_match(study_parameters, "Enrichment Analysis UI Parameters", fixed = TRUE)
    expect_match(study_parameters, sprintf("database_source = %s", backend), fixed = TRUE)
    expect_match(study_parameters, sprintf("organism_selected = %s", organism_taxid), fixed = TRUE)

    .e2e_assert_enrichment_publication_workbook(paths, backend)

    e2e_assert_no_console_errors(driver)
    e2e_assert_no_error_notifications(driver)
  })
}

test_that("E2E-007 gprofiler2 enrichment branch reaches export/report artifacts", {
  .e2e_run_prot_enrichment_backend_lane("gprofiler2")
})

test_that("E2E-007 clusterProfiler enrichment branch reaches export/report artifacts", {
  .e2e_run_prot_enrichment_backend_lane("clusterprofiler")
})
