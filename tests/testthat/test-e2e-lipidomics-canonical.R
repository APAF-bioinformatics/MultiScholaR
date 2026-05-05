library(testthat)

.e2e_assert_lipid_import_state <- function(driver, lane, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  digest <- e2e_state_digest(driver, timeout = timeout)
  setup_log <- digest$export_paths$lipidomics$setup_import
  expected_assays <- unlist(lane$assays)

  expect_false(is.null(setup_log))
  expect_identical(setup_log$n_assays, length(expected_assays))
  expect_setequal(setup_log$assay_names, expected_assays)
  expect_identical(setup_log$n_samples, lane$sample_count)
  expect_identical(setup_log$detected_format, "lipidsearch")

  invisible(setup_log)
}

.e2e_assert_lipid_design_artifacts <- function(paths, expected_assays) {
  e2e_assert_file_nonempty(file.path(paths$source_dir, "design_matrix.tab"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "contrast_strings.tab"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "assay_manifest.txt"))
  e2e_assert_file_nonempty(file.path(paths$source_dir, "column_mapping.json"))

  assay_manifest <- readLines(file.path(paths$source_dir, "assay_manifest.txt"), warn = FALSE)
  expect_setequal(assay_manifest, expected_assays)

  column_mapping <- jsonlite::read_json(file.path(paths$source_dir, "column_mapping.json"), simplifyVector = TRUE)
  expect_identical(column_mapping$lipid_id_col, "LipidName")
  expect_identical(column_mapping$annotation_col, "LipidClass")
  expect_setequal(column_mapping$sample_columns, c("WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"))

  for (assay in expected_assays) {
    e2e_assert_file_nonempty(file.path(paths$source_dir, paste0("data_cln_", assay, ".tab")))
  }

  invisible(paths)
}

.e2e_lipid_expected_prefixes <- function(expected_assays) {
  vapply(expected_assays, function(assay) {
    if (grepl("pos", assay, ignore.case = TRUE)) {
      "posmode"
    } else if (grepl("neg", assay, ignore.case = TRUE)) {
      "negmode"
    } else {
      gsub("[^A-Za-z0-9]", "", tolower(assay))
    }
  }, character(1))
}

.e2e_assert_lipid_da_artifacts <- function(paths, lane, expected_assays) {
  da_files <- list.files(
    paths$da_output_dir,
    pattern = "^de_.+_lipids_.+_long_annot\\.(tsv|xlsx)$",
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
  expect_true("assay" %in% names(da_table))
  expect_true(any(da_table$comparison %in% expected_raw))
  expect_true(any(da_table$friendly_name %in% expected_friendly))
  expect_setequal(unique(da_table$assay), expected_assays)

  present_prefixes <- unique(sub("^de_([^_]+)_lipids_.*$", "\\1", basename(da_files)))
  expect_true(all(.e2e_lipid_expected_prefixes(expected_assays) %in% present_prefixes))

  invisible(da_files)
}

.e2e_assert_lipid_summary_artifacts <- function(paths, lane, expected_assays) {
  e2e_assert_file_nonempty(file.path(paths$source_dir, "study_parameters.txt"))
  study_parameters <- paste(readLines(file.path(paths$source_dir, "study_parameters.txt"), warn = FALSE), collapse = "\n")
  expected_friendly <- unlist(lane$expected_contrasts)
  expected_raw <- paste0(
    "group",
    sub("_vs_.*$", "", expected_friendly),
    "-group",
    sub("^.*_vs_", "", expected_friendly)
  )

  expect_match(study_parameters, "lipidomics", ignore.case = TRUE)
  expect_match(study_parameters, expected_friendly[[1]], fixed = TRUE)
  expect_match(study_parameters, expected_raw[[1]], fixed = TRUE)
  for (assay in expected_assays) {
    expect_match(study_parameters, assay, fixed = TRUE)
  }

  publication_table <- file.path(
    paths$results_summary_dir,
    "Publication_tables",
    "DA_results_lipidomics.xlsx"
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

.e2e_run_lipid_import_only <- function(lane, case_id) {
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  paths <- e2e_lipidomics_project_paths(project_root, case_id)

  driver <- e2e_new_app_driver(
    case_id = case_id,
    launch_surface = "run_app",
    timeout = 60000L,
    load_timeout = 180000L
  )
  on.exit(e2e_stop_driver(driver), add = TRUE)

  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-lipidomics", timeout = 60000L)
    e2e_select_omics(driver, "lipidomics", assert_digest = TRUE, timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    e2e_wait_for_selector(driver, "lipid-tab-import", timeout = 60000L)
    e2e_assert_project_dirs_exist(paths)

    e2e_upload_lane_inputs(driver, lane)
    e2e_click_process_import(driver, "lipidomics")
    e2e_wait_for_step_status(driver, "lipidomics", "setup_import", "complete", timeout = 120000L)
    .e2e_assert_lipid_import_state(driver, lane, timeout = 120000L)

    e2e_assert_no_console_errors(driver)
    e2e_assert_no_error_notifications(driver)
  })
}

test_that("E2E-009 lipidomics LCMS positive plus negative reaches DA and report artifacts", {
  skip_if_no_e2e_browser()

  lane <- read_e2e_manifest()[["lipid_canonical"]]
  expected_assays <- unlist(lane$assays)
  case_id <- "E2E-009-lipid_canonical"
  project_root <- file.path(tempdir(), paste0(case_id, "-", Sys.getpid()))
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  project_base <- e2e_project_base_dir(project_root, case_id)
  paths <- e2e_lipidomics_project_paths(project_root, case_id)

  driver <- e2e_new_app_driver(
    case_id = case_id,
    launch_surface = "run_app",
    timeout = 60000L,
    load_timeout = 180000L
  )
  on.exit(e2e_stop_driver(driver), add = TRUE)

  e2e_with_failure_artifacts(driver, case_id, {
    e2e_wait_for_selector(driver, "tile-lipidomics", timeout = 60000L)
    e2e_select_omics(driver, "lipidomics", assert_digest = TRUE, timeout = 60000L)
    e2e_complete_project_setup(driver, case_id, project_root, timeout = 60000L)
    e2e_wait_for_selector(driver, "lipid-tab-import", timeout = 60000L)
    e2e_assert_project_dirs_exist(paths)

    e2e_upload_lane_inputs(driver, lane)
    e2e_click_process_import(driver, "lipidomics")
    e2e_wait_for_step_status(driver, "lipidomics", "setup_import", "complete", timeout = 120000L)
    .e2e_assert_lipid_import_state(driver, lane, timeout = 120000L)

    e2e_complete_lipid_design_from_manifest(driver, lane, timeout = 180000L)
    e2e_wait_for_step_status(driver, "lipidomics", "design_matrix", "complete", timeout = 180000L)
    .e2e_assert_lipid_design_artifacts(paths, expected_assays)

    e2e_run_lipid_qc_finalize(driver, timeout = 120000L)
    e2e_wait_for_step_status(driver, "lipidomics", "quality_control", "complete", timeout = 120000L)

    e2e_run_lipid_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "lipidomics", "normalization", "complete", timeout = 240000L)
    session_artifacts <- e2e_assert_lipid_filtered_session_artifacts(project_base, expected_assays = expected_assays)
    expect_setequal(session_artifacts$session_data$feature_counts |> names(), expected_assays)

    e2e_run_lipid_da(driver, expected_assays = expected_assays, timeout = 240000L)
    e2e_wait_for_step_status(driver, "lipidomics", "differential_analysis", "complete", timeout = 240000L)
    .e2e_assert_lipid_da_artifacts(paths, lane, expected_assays)

    report_download <- e2e_run_lipid_summary_report(
      driver,
      lane = lane,
      project_base_dir = project_base,
      experiment_label = case_id,
      case_id = paste0(case_id, "-report"),
      description = "E2E-009 lipidomics canonical multi-assay browser workflow",
      timeout = 240000L
    )
    e2e_assert_file_nonempty(report_download)
    .e2e_assert_lipid_summary_artifacts(paths, lane, expected_assays)

    e2e_assert_no_console_errors(driver)
    e2e_assert_no_error_notifications(driver)
  })
})

test_that("E2E-009 lipidomics GCMS-named assay imports without LC-only assumptions", {
  skip_if_no_e2e_browser()

  lane <- read_e2e_manifest()[["lipid_canonical"]]
  lane$assay2_file <- "lipidsearch_gcms.txt"
  lane$assays <- list("LCMS_Pos", "GCMS")

  .e2e_run_lipid_import_only(lane, "E2E-009-lipid_gcms_import")
})
