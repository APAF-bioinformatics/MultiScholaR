library(testthat)

new_fake_e2e_driver <- function(digest = NULL, logs = list()) {
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  if (is.null(digest)) {
    digest <- jsonlite::toJSON(
      list(
        selected_omics = c("proteomics"),
        initialized_omics = c("proteomics"),
        project_dir_keys = c("proteomics", "scripts"),
        workflow_type_per_omic = list(proteomics = "DIA"),
        step_status_per_omic = list(proteomics = list(setup_import = "complete")),
        active_tab_per_omic = list(proteomics = "setup"),
        export_paths = list(proteomics = list()),
        report_fingerprints = list()
      ),
      auto_unbox = TRUE
    )
  }
  state$digest <- digest
  state$logs <- logs
  state$stopped <- FALSE

  record <- function(method, args) {
    state$calls <- c(state$calls, list(list(method = method, args = args)))
  }

  list(
    .state = state,
    click = function(...) {
      record("click", list(...))
      invisible(NULL)
    },
    wait_for_idle = function(...) {
      record("wait_for_idle", list(...))
      invisible(NULL)
    },
    wait_for_js = function(...) {
      record("wait_for_js", list(...))
      TRUE
    },
    set_inputs = function(...) {
      record("set_inputs", list(...))
      invisible(NULL)
    },
    upload_file = function(...) {
      record("upload_file", list(...))
      invisible(NULL)
    },
    get_value = function(output) {
      record("get_value", list(output = output))
      if (identical(output, "test_state_digest")) state$digest else NULL
    },
    get_logs = function() {
      record("get_logs", list())
      state$logs
    },
    get_js = function(script) {
      record("get_js", list(script = script))
      list()
    },
    get_html = function(...) {
      record("get_html", list(...))
      "<html><body>debug</body></html>"
    },
    get_screenshot = function(path, ...) {
      record("get_screenshot", c(list(path = path), list(...)))
      writeBin(as.raw(c(1, 2, 3, 4)), path)
      path
    },
    get_download = function(output_id, ...) {
      record("get_download", c(list(output_id = output_id), list(...)))
      path <- tempfile(fileext = ".html")
      writeLines("<html>report</html>", path)
      path
    },
    stop = function() {
      state$stopped <- TRUE
      invisible(NULL)
    }
  )
}

test_that("e2e_selector builds stable data-testid selectors", {
  expect_identical(e2e_selector("tile-proteomics"), "[data-testid=\"tile-proteomics\"]")
  expect_error(e2e_selector(""), "non-empty")
  expect_identical(
    e2e_input_selector("proteomics_workflow-quality_control-protein_qc-replicate_filter-apply_protein_replicate_filter"),
    "[id=\"proteomics_workflow-quality_control-protein_qc-replicate_filter-apply_protein_replicate_filter\"]"
  )
  expect_error(e2e_input_selector(""), "non-empty")
})

test_that("e2e_input_id maps all current omic module namespaces", {
  expect_identical(
    e2e_input_id("proteomics", "import", "search_results_standard"),
    "proteomics_workflow-setup_import-search_results_standard"
  )
  expect_identical(
    e2e_input_id("metabolomics", "import", "assay1_file_std"),
    "metabolomics_workflow-import-assay1_file_std"
  )
  expect_identical(
    e2e_input_id("lipidomics", "summary", "generate_report"),
    "lipidomics_workflow-summary-generate_report"
  )
  expect_identical(
    e2e_input_id("proteomics", "qc", "peptide_qc-qvalue_filter-apply_qvalue_filter"),
    "proteomics_workflow-quality_control-peptide_qc-qvalue_filter-apply_qvalue_filter"
  )
  expect_error(e2e_input_id("transcriptomics", "import", "file"), "Unsupported omic")
})

test_that("e2e_launch_spec supports package and app-dir launch surfaces", {
  expect_true(is.function(e2e_launch_spec("run_app")))
  expect_true(dir.exists(e2e_launch_spec("app_dir")))
})

test_that("click, input, upload, and setup helpers call the driver through stable contracts", {
  driver <- new_fake_e2e_driver()
  upload <- tempfile(fileext = ".tsv")
  writeLines("a\tb", upload)
  project_dir <- tempdir()

  e2e_click_testid(driver, "btn-start-analysis")
  e2e_set_input(driver, "experiment_label", "e2e_case")
  e2e_upload_file(driver, "proteomics_workflow-setup_import-search_results_standard", upload)
  e2e_click_input_id(driver, "proteomics_workflow-design_matrix-builder-save_results")
  e2e_trigger_action_input_id(driver, "proteomics_workflow-quality_control-protein_qc-replicate_filter-apply_protein_replicate_filter")
  e2e_complete_project_setup(driver, "e2e_case", project_dir)
  e2e_click_summary_action(driver, "proteomics", "copy_publication")

  methods <- vapply(driver$.state$calls, `[[`, character(1), "method")
  expect_true("click" %in% methods)
  expect_true("set_inputs" %in% methods)
  expect_true("upload_file" %in% methods)
  expect_true("wait_for_idle" %in% methods)
  expect_true(any(vapply(
    driver$.state$calls,
    function(call) {
      identical(call$method, "get_js") &&
        grepl(
          "proteomics_workflow-quality_control-protein_qc-replicate_filter-apply_protein_replicate_filter",
          call$args$script,
          fixed = TRUE
        ) &&
        grepl("Shiny.setInputValue", call$args$script, fixed = TRUE) &&
        grepl("priority: 'event'", call$args$script, fixed = TRUE)
    },
    logical(1)
  )))
})

test_that("selectize helpers resolve imported run case without inventing fixture-only IDs", {
  expect_identical(
    e2e_match_available_values(
      c("WT_1", "KO_1"),
      c("wt_1", "wt_2", "ko_1", "ko_2")
    ),
    c("wt_1", "ko_1")
  )

  expect_identical(
    e2e_match_available_values("WT_1", character()),
    "WT_1"
  )
})

test_that("workflow tab helper maps proteomics step keys to UI tab values", {
  driver <- new_fake_e2e_driver()

  e2e_switch_workflow_tab(driver, "proteomics", "design_matrix")
  e2e_switch_workflow_tab(driver, "proteomics", "quality_control")
  e2e_switch_workflow_tab(driver, "proteomics", "differential_expression")

  js_calls <- Filter(function(call) identical(call$method, "get_js"), driver$.state$calls)
  scripts <- paste(vapply(js_calls, function(call) call$args$script, character(1L)), collapse = "\n")
  expect_match(scripts, '"design"', fixed = TRUE)
  expect_match(scripts, '"qc"', fixed = TRUE)
  expect_match(scripts, '"da"', fixed = TRUE)
})

test_that("state digest assertions parse test_state_digest and validate invariants", {
  driver <- new_fake_e2e_driver()

  digest <- e2e_assert_state_digest(
    driver,
    selected_omics = "proteomics",
    initialized_omics = "proteomics",
    required_project_dir_keys = "proteomics",
    workflow_types = list(proteomics = "DIA")
  )

  expect_identical(unlist(digest$selected_omics), "proteomics")
})

test_that("step-status assertions and waits inspect per-omic digest state", {
  digest <- jsonlite::toJSON(
    list(
      selected_omics = c("proteomics"),
      initialized_omics = c("proteomics"),
      project_dir_keys = c("proteomics"),
      workflow_type_per_omic = list(proteomics = "DIA"),
      step_status_per_omic = list(
        proteomics = list(
          setup_import = "complete",
          design_matrix = "complete",
          quality_control = "complete",
          normalization = "complete",
          differential_expression = "complete"
        )
      ),
      active_tab_per_omic = list(proteomics = "session_summary"),
      export_paths = list(proteomics = list()),
      report_fingerprints = list()
    ),
    auto_unbox = TRUE
  )
  driver <- new_fake_e2e_driver(digest = digest)

  expect_no_error(e2e_assert_step_statuses(
    driver,
    "proteomics",
    list(setup_import = "complete", differential_expression = "complete")
  ))
  expect_no_error(e2e_wait_for_step_status(
    driver,
    "proteomics",
    "normalization",
    "complete",
    timeout = 1000
  ))
})

test_that("project path and filtered-session artifact assertions validate DIA handoff files", {
  project_root <- tempfile("e2e-project-")
  experiment_label <- "E2E-004"
  dir.create(file.path(project_root, experiment_label, "scripts", "proteomics"), recursive = TRUE)
  dir.create(file.path(project_root, experiment_label, "data", "proteomics"), recursive = TRUE)
  dir.create(file.path(project_root, experiment_label, "results", "proteomics", "da_proteins"), recursive = TRUE)
  dir.create(file.path(project_root, experiment_label, "results_summary", "proteomics"), recursive = TRUE)

  paths <- e2e_proteomics_project_paths(project_root, experiment_label)
  expect_no_error(e2e_assert_project_dirs_exist(paths))

  lane <- read_e2e_manifest()[["prot_dia"]]
  seeded <- e2e_seed_report_template(lane, paths$base_dir)
  expect_true(file.exists(seeded))

  saveRDS(
    list(
      r6_current_state_name = "correlation_filtered",
      r6_complete_states = list(correlation_filtered = list(marker = TRUE)),
      r6_state_history = c("protein_replicate_filtered", "correlation_filtered")
    ),
    file.path(paths$source_dir, "filtered_session_data_latest.rds")
  )
  writeLines("Filtered Session Data Export Summary", file.path(paths$source_dir, "filtered_session_summary.txt"))
  expect_no_error(e2e_assert_prot_filtered_session_artifacts(paths$base_dir))
})

test_that("proteomics lane upload helper supplies search results and required FASTA", {
  driver <- new_fake_e2e_driver()
  lane <- read_e2e_manifest()[["prot_dia"]]

  e2e_upload_lane_inputs(driver, lane)

  upload_calls <- Filter(
    function(call) identical(call$method, "upload_file"),
    driver$.state$calls
  )
  expect_equal(length(upload_calls), 2L)
  uploaded_ids <- vapply(upload_calls, function(call) names(call$args)[[1L]], character(1L))
  expect_setequal(
    uploaded_ids,
    c(
      e2e_input_id("proteomics", "import", "search_results_standard"),
      e2e_input_id("proteomics", "import", "fasta_file_standard")
    )
  )
})

test_that("proteomics DIA scenario helpers dispatch design, QC, normalization, and summary actions", {
  driver <- new_fake_e2e_driver()
  lane <- read_e2e_manifest()[["prot_dia"]]
  project_root <- tempfile("e2e-dia-project-")
  experiment_label <- "E2E-004"
  dir.create(project_root, recursive = TRUE)
  paths <- e2e_proteomics_project_paths(project_root, experiment_label)
  dir.create(paths$source_dir, recursive = TRUE)

  e2e_complete_prot_design_from_manifest(driver, lane)
  e2e_run_prot_dia_qc(driver)
  e2e_run_prot_dia_normalization_export(driver)
  report_download <- e2e_run_prot_da_summary_report(
    driver,
    lane = lane,
    project_base_dir = paths$base_dir,
    experiment_label = experiment_label
  )

  methods <- vapply(driver$.state$calls, `[[`, character(1), "method")
  expect_true("set_inputs" %in% methods)
  expect_true("get_js" %in% methods)
  expect_true("get_download" %in% methods)
  expect_true(file.exists(report_download))
})

test_that("console and notification assertions flag driver-visible failures", {
  clean_driver <- new_fake_e2e_driver(logs = list(list(level = "INFO", text = "ok")))
  expect_no_error(e2e_assert_no_console_errors(clean_driver))

  datatable_lifecycle_driver <- new_fake_e2e_driver(logs = list(list(
    location = "chromote",
    level = "log",
    message = paste(
      "Couldn't find table with id proteomics_workflow-design_matrix-builder-data_table",
      "datatables-binding-0.34.0/datatables.js:1526:14"
    )
  )))
  expect_no_error(e2e_assert_no_console_errors(datatable_lifecycle_driver))

  error_driver <- new_fake_e2e_driver(logs = list(list(level = "SEVERE", text = "ReferenceError: x")))
  expect_failure(e2e_assert_no_console_errors(error_driver))
})

test_that("failure artifact capture writes actionable debug files", {
  withr::with_envvar(
    c(MULTISCHOLAR_E2E_ARTIFACT_DIR = tempfile("e2e-artifacts-"))
    , {
      driver <- new_fake_e2e_driver()
      artifact_dir <- e2e_capture_failure_artifacts(
        driver,
        case_id = "E2E-003 fake failure",
        reason = "synthetic failure"
      )

      expect_true(file.exists(file.path(artifact_dir, "failure-reason.txt")))
      expect_true(file.exists(file.path(artifact_dir, "dom.html")))
      expect_true(file.exists(file.path(artifact_dir, "browser-logs.json")))
      expect_true(file.exists(file.path(artifact_dir, "screenshot.png")))
    }
  )
})

test_that("run_app shinytest2 smoke reaches the home selector when browser dependencies exist", {
  skip_if_no_e2e_browser()

  driver <- e2e_new_app_driver(
    case_id = "E2E-003-run-app-smoke",
    launch_surface = "run_app"
  )
  on.exit(e2e_stop_driver(driver), add = TRUE)

  e2e_wait_for_selector(driver, "tile-proteomics")
  e2e_assert_no_console_errors(driver)
})
