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
    guiEntryInstallUniprotCache(paths)

    e2e_upload_lane_inputs(driver, lane)
    e2e_set_input_and_idle(driver, e2e_input_id("proteomics", "import", "dia_route"), "DIA_LIMPA", timeout = 60000L)
    e2e_click_process_import(driver, "proteomics")
    e2e_wait_for_step_status(driver, "proteomics", "setup_import", "complete", timeout = 120000L)
    e2e_assert_state_digest(
      driver,
      selected_omics = "proteomics",
      initialized_omics = "proteomics",
      required_project_dir_keys = "proteomics",
      workflow_types = list(proteomics = "DIA_LIMPA"),
      timeout = 60000L
    )

    e2e_complete_prot_design_from_manifest(driver, lane, timeout = 120000L)
    e2e_wait_for_step_status(driver, "proteomics", "design_matrix", "complete", timeout = 120000L)
    e2e_assert_no_uniprot_network(driver)

    # KIT-104 F-25 (D12): the route's QC ends at peptide_count_filtered; its QC tab status is not waited on (P-16).
    e2e_run_prot_dia_limpa_qc(driver, timeout = 180000L)
    e2e_wait_for_r6_state(driver, "proteomics", "peptide_count_filtered", timeout = 60000L)

    e2e_run_prot_dia_limpa_normalization_export(driver, timeout = 240000L)
    e2e_wait_for_step_status(driver, "proteomics", "normalization", "complete", timeout = 240000L)

    # KIT-104 F-25: the native export's checks, expected values from the prot_dia_limpa entry oracle's session_export
    # step (the head and history it exports, and the files the export adds to the source directory: the archive, its
    # pointer, completion record and summary, and the metadata sidecars; the legacy limpa sidecars are not among them).
    oracle <- jsonlite::read_json(guiEntryOraclePath("prot_dia_limpa"), simplifyVector = FALSE)$steps
    exported <- oracle$session_export
    digest <- e2e_state_digest(driver, timeout = 60000L)
    expect_identical(e2e_digest_r6_current_state(digest, "proteomics"), exported$head)
    expect_identical(e2e_digest_r6_state_history(digest, "proteomics"), as.character(unlist(exported$history)))

    session_artifacts <- e2e_assert_prot_filtered_session_artifacts(project_base)
    session_data <- session_artifacts$session_data
    # The export summary's Report Template is not checked: the native export reads it from globalParameters, which no
    # route sets since KIT-046 (the report step chooses the template by workflow type), so it reads NA on every route
    # (F104-01, queued to KIT-108). The template is checked at the report below.
    expect_identical(session_data$workflow_type, "DIA_LIMPA")
    expect_true(isTRUE(session_data$limpa_applied))
    expect_equal(session_data$final_sample_count, lane$sample_count)
    expect_gt(session_data$final_protein_count, 0L)

    added <- setdiff(names(exported$artifacts), names(oracle$ruv_covariates$artifacts))
    added <- sub("^scripts/proteomics/", "", grep("^scripts/proteomics/[^/]+$", added, value = TRUE))
    expect_length(added, 8L)
    source_files <- list.files(paths$source_dir)
    for (name in added) {
      pattern <- paste0("^", gsub("<TIME>", "[0-9]{8}_[0-9]{6}", gsub("<HEX>", "[A-Za-z0-9_-]+",
        gsub(".", "[.]", name, fixed = TRUE), fixed = TRUE), fixed = TRUE), "$")
      matched <- grep(pattern, source_files, value = TRUE)
      expect_length(matched, 1L)
      invisible(lapply(file.path(paths$source_dir, matched), e2e_assert_file_nonempty))
    }
    expect_false(any(c("limpa_parameters.RDS", "limpa_dpc_quant_results.RDS") %in% source_files))
    filtered_summary <- paste(readLines(session_artifacts$summary, warn = FALSE), collapse = "\n")
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

    # KIT-104 W4 (F-25): the native path's study parameters: the DPC-Quant args groups (limpaDpcQuant,
    # limpa_dpc_quant_results), not the legacy proteinMissingValueImputationLimpa group or a template in
    # globalParameters (the template is checked above by the copied report file; F104-01).
    study_parameters <- paste(readLines(file.path(paths$source_dir, "study_parameters.txt"), warn = FALSE), collapse = "\n")
    expect_match(study_parameters, "[limpaDpcQuant]", fixed = TRUE)
    expect_match(study_parameters, "[limpa_dpc_quant_results]", fixed = TRUE)
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
