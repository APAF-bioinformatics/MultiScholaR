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
    stop = function() {
      state$stopped <- TRUE
      invisible(NULL)
    }
  )
}

test_that("e2e_selector builds stable data-testid selectors", {
  expect_identical(e2e_selector("tile-proteomics"), "[data-testid=\"tile-proteomics\"]")
  expect_error(e2e_selector(""), "non-empty")
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
  e2e_complete_project_setup(driver, "e2e_case", project_dir)
  e2e_click_summary_action(driver, "proteomics", "copy_publication")

  methods <- vapply(driver$.state$calls, `[[`, character(1), "method")
  expect_true("click" %in% methods)
  expect_true("set_inputs" %in% methods)
  expect_true("upload_file" %in% methods)
  expect_true("wait_for_idle" %in% methods)
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

test_that("console and notification assertions flag driver-visible failures", {
  clean_driver <- new_fake_e2e_driver(logs = list(list(level = "INFO", text = "ok")))
  expect_no_error(e2e_assert_no_console_errors(clean_driver))

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
