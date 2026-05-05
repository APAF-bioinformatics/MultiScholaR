library(testthat)

.e2e_assert_metab_import_state <- function(driver, lane, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  digest <- e2e_state_digest(driver, timeout = timeout)
  setup_log <- digest$export_paths$metabolomics$setup_import
  expected_assays <- unlist(lane$assays)

  expect_false(is.null(setup_log))
  expect_identical(setup_log$n_assays, length(expected_assays))
  expect_setequal(setup_log$assay_names, expected_assays)
  expect_identical(setup_log$n_samples, lane$sample_count)
  expect_identical(setup_log$detected_format, "custom")

  invisible(setup_log)
}

.e2e_assert_metab_design_artifacts <- function(paths, expected_assays) {
  e2e_assert_file_nonempty(file.path(paths$source_dir, "design_matrix.tab"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "contrast_strings.tab"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "assay_manifest.txt"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "column_mapping.json"))

  assay_manifest <- readLines(file.path(paths$source_dir, "assay_manifest.txt"), warn = FALSE)
  expect_setequal(assay_manifest, expected_assays)

  for (assay in expected_assays) {
    e2e_assert_file_nonempty(file.path(paths$source_dir, paste0("data_cln_", assay, ".tab")))
  }

  invisible(paths)
}

.e2e_assert_metab_da_artifacts <- function(paths, lane, expected_assays) {
  da_files <- list.files(
    paths$da_output_dir,
    pattern = "^de_.+_metabolites_.+_long_annot\\.(tsv|xlsx)$",
    full.names = TRUE
  )
  expect_gt(length(da_files), 0L)

  tsv_files <- da_files[grepl("\\.tsv$", da_files)]
  expect_gt(length(tsv_files), 0L)
  da_table <- do.call(
    rbind,
    lapply(tsv_files, read.delim, stringsAsFactors = FALSE, check.names = FALSE)
  )
  expected_friendly <- unlist(lane$expected_contrasts)
  expected_raw <- paste0(
    "group",
    sub("_vs_.*$", "", expected_friendly),
    "-group",
    sub("^.*_vs_", "", expected_friendly)
  )

  expect_true("comparison" %in% names(da_table))
  expect_true("friendly_name" %in% names(da_table))
  expect_true("fdr_qvalue" %in% names(da_table))
  expect_true(any(da_table$comparison %in% expected_raw))
  expect_true(any(da_table$friendly_name %in% expected_friendly))

  expected_prefixes <- vapply(expected_assays, function(assay) {
    if (grepl("pos", assay, ignore.case = TRUE)) {
      "posmode"
    } else if (grepl("neg", assay, ignore.case = TRUE)) {
      "negmode"
    } else {
      gsub("[^A-Za-z0-9]", "", tolower(assay))
    }
  }, character(1))
  present_prefixes <- unique(sub("^de_([^_]+)_metabolites_.*$", "\\1", basename(da_files)))
  expect_true(all(expected_prefixes %in% present_prefixes))

  invisible(da_files)
}

.e2e_assert_metab_summary_artifacts <- function(paths, lane, expected_assays) {
  e2e_assert_file_nonempty(file.path(paths$source_dir, "study_parameters.txt"))
  study_parameters <- paste(readLines(file.path(paths$source_dir, "study_parameters.txt"), warn = FALSE), collapse = "\n")
  expected_friendly <- unlist(lane$expected_contrasts)
  expected_raw <- paste0(
    "group",
    sub("_vs_.*$", "", expected_friendly),
    "-group",
    sub("^.*_vs_", "", expected_friendly)
  )

  expect_match(study_parameters, "metabolomics", ignore.case = TRUE)
  expect_match(study_parameters, expected_friendly[[1]], fixed = TRUE)
  expect_match(study_parameters, expected_raw[[1]], fixed = TRUE)

  publication_table <- file.path(
    paths$results_summary_dir,
    "Publication_tables",
    "DA_results_metabolomics.xlsx"
  )
  e2e_assert_file_nonempty(publication_table)

  session_exports <- list.files(paths$source_dir, pattern = "^session_state_.*\\.RDS$", full.names = TRUE)
  expect_gt(length(session_exports), 0L)

  invisible(list(
    study_parameters = study_parameters,
    publication_table = publication_table,
    session_exports = session_exports,
    expected_assays = expected_assays,
    report_template = lane$report_template
  ))
}

.e2e_run_metab_lane <- function(lane_id) {
  skip_if_no_e2e_browser()

  lane <- read_e2e_manifest()[[lane_id]]
  expected_assays <- unlist(lane$assays)
  case_id <- paste("E2E-008", lane_id, sep = "-")
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  project_base <- e2e_project_base_dir(project_root, case_id)
  paths <- e2e_metabolomics_project_paths(project_root, case_id)

  driver <- e2e_new_app_driver(
    case_id = case_id,
    launch_surface = "run_app",
    timeout = 60000L,
    load_timeout = 180000L
  )
  on.exit(e2e_stop_driver(driver), add = TRUE)

  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-metabolomics", timeout = 60000L)
    e2e_select_omics(driver, "metabolomics", assert_digest = TRUE, timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    e2e_wait_for_selector(driver, "metab-tab-import", timeout = 60000L)
    e2e_assert_project_dirs_exist(paths)

    e2e_upload_lane_inputs(driver, lane)
    e2e_click_process_import(driver, "metabolomics")
    e2e_wait_for_step_status(driver, "metabolomics", "setup_import", "complete", timeout = 120000L)
    .e2e_assert_metab_import_state(driver, lane, timeout = 120000L)

    e2e_complete_metab_design_from_manifest(driver, lane, timeout = 180000L)
    e2e_wait_for_step_status(driver, "metabolomics", "design_matrix", "complete", timeout = 180000L)
    .e2e_assert_metab_design_artifacts(paths, expected_assays)

    e2e_run_metab_qc_finalize(driver, timeout = 120000L)
    e2e_wait_for_step_status(driver, "metabolomics", "quality_control", "complete", timeout = 120000L)

    e2e_run_metab_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "metabolomics", "normalization", "complete", timeout = 240000L)
    session_artifacts <- e2e_assert_metab_filtered_session_artifacts(project_base, expected_assays = expected_assays)
    expect_setequal(session_artifacts$session_data$feature_counts |> names(), expected_assays)

    e2e_run_metab_da(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "metabolomics", "differential_analysis", "complete", timeout = 240000L)
    .e2e_assert_metab_da_artifacts(paths, lane, expected_assays)

    report_download <- e2e_run_metab_summary_report(
      driver,
      lane = lane,
      project_base_dir = project_base,
      experiment_label = case_id,
      case_id = paste0(case_id, "-report"),
      description = sprintf("E2E-008 metabolomics browser workflow for %s", lane_id),
      timeout = 240000L
    )
    e2e_assert_file_nonempty(report_download)
    .e2e_assert_metab_summary_artifacts(paths, lane, expected_assays)

    e2e_assert_no_console_errors(driver)
    e2e_assert_no_error_notifications(driver)
  })
}

test_that("E2E-008 metabolomics LC lane reaches DA and report artifacts", {
  .e2e_run_metab_lane("metab_lc")
})

test_that("E2E-008 metabolomics GC lane reaches DA and report artifacts", {
  .e2e_run_metab_lane("metab_gc")
})

test_that("E2E-008 metabolomics combined LC plus GC lane preserves assay partitioning", {
  .e2e_run_metab_lane("metab_combined")
})
