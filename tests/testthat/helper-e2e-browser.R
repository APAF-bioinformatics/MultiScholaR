# Shared shinytest2 browser harness for GUI E2E tests.
# The helpers keep workflow tests readable and centralize failure artifacts.

.E2E_BROWSER_DEFAULT_TIMEOUT <- 30000L
.E2E_BROWSER_DEFAULT_LOAD_TIMEOUT <- 120000L
.E2E_VALID_OMICS <- c("proteomics", "metabolomics", "lipidomics")

e2e_browser_dependencies_available <- function() {
  if (!requireNamespace("shinytest2", quietly = TRUE)) {
    return(FALSE)
  }
  if (!requireNamespace("chromote", quietly = TRUE)) {
    return(FALSE)
  }

  chrome_path <- tryCatch(
    chromote::find_chrome(),
    error = function(e) ""
  )
  any(nzchar(chrome_path))
}

skip_if_no_e2e_browser <- function() {
  testthat::skip_if_not(
    e2e_browser_dependencies_available(),
    "shinytest2/chromote browser dependencies are not available"
  )
}

e2e_sanitize_id <- function(value) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(value))
}

e2e_artifact_root <- function() {
  root <- Sys.getenv("MULTISCHOLAR_E2E_ARTIFACT_DIR", unset = "")
  if (!nzchar(root)) {
    root <- testthat::test_path("_e2e_artifacts")
  }
  root
}

e2e_case_artifact_dir <- function(case_id, create = TRUE) {
  path <- file.path(e2e_artifact_root(), e2e_sanitize_id(case_id))
  if (isTRUE(create) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

e2e_selector <- function(test_id) {
  if (length(test_id) != 1L || !nzchar(test_id)) {
    stop("test_id must be a non-empty character scalar", call. = FALSE)
  }
  escaped <- gsub("\"", "\\\\\"", test_id, fixed = TRUE)
  sprintf("[data-testid=\"%s\"]", escaped)
}

e2e_launch_spec <- function(launch_surface = c("run_app", "app_dir")) {
  launch_surface <- match.arg(launch_surface)

  switch(
    launch_surface,
    run_app = function() {
      options(multischolar.test_mode = TRUE)
      MultiScholaR::run_app(test_mode = TRUE)
    },
    app_dir = normalizePath(
      testthat::test_path("..", "..", "inst", "app"),
      mustWork = TRUE
    )
  )
}

e2e_platform_variant <- function() {
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  paste(
    e2e_sanitize_id(Sys.info()[["sysname"]]),
    paste0("r", e2e_sanitize_id(r_version)),
    sep = "-"
  )
}

e2e_new_app_driver <- function(
    case_id,
    launch_surface = c("run_app", "app_dir"),
    seed = 20260428L,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    load_timeout = .E2E_BROWSER_DEFAULT_LOAD_TIMEOUT,
    width = 1440L,
    height = 1100L,
    artifact_dir = e2e_case_artifact_dir(case_id),
    ...
) {
  skip_if_no_e2e_browser()
  launch_surface <- match.arg(launch_surface)
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }

  AppDriver <- shinytest2::AppDriver
  do.call(
    AppDriver$new,
    c(
      list(
        app_dir = e2e_launch_spec(launch_surface),
        name = e2e_sanitize_id(case_id),
        variant = e2e_platform_variant(),
        seed = seed,
        timeout = timeout,
        load_timeout = load_timeout,
        width = width,
        height = height,
        screenshot_args = list(selector = "viewport"),
        expect_values_screenshot_args = FALSE,
        options = list(
          multischolar.test_mode = TRUE,
          shiny.testmode = TRUE
        )
      ),
      list(...)
    )
  )
}

e2e_driver_has_method <- function(driver, method) {
  is.function(driver[[method]])
}

e2e_call_driver_method <- function(driver, method, ..., required = TRUE) {
  if (!e2e_driver_has_method(driver, method)) {
    if (isTRUE(required)) {
      stop(
        sprintf("E2E driver does not expose method '%s'", method),
        call. = FALSE
      )
    }
    return(NULL)
  }

  driver[[method]](...)
}

e2e_wait_for_idle <- function(
    driver,
    duration = 200L,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_call_driver_method(
    driver,
    "wait_for_idle",
    duration = duration,
    timeout = timeout,
    required = FALSE
  )
  invisible(driver)
}

e2e_wait_for_js <- function(
    driver,
    script,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_call_driver_method(
    driver,
    "wait_for_js",
    script = script,
    timeout = timeout
  )
}

e2e_wait_for_selector <- function(
    driver,
    test_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  selector <- e2e_selector(test_id)
  script <- sprintf(
    "document.querySelector(%s) !== null",
    jsonlite::toJSON(selector, auto_unbox = TRUE)
  )
  e2e_wait_for_js(driver, script = script, timeout = timeout)
}

e2e_click_testid <- function(
    driver,
    test_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_call_driver_method(driver, "click", selector = e2e_selector(test_id))
  e2e_wait_for_idle(driver, timeout = timeout)
  invisible(driver)
}

e2e_set_input <- function(
    driver,
    input_id,
    value,
    wait = TRUE,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    allow_no_input_binding = FALSE,
    priority = c("input", "event")
) {
  priority <- match.arg(priority)
  args <- list(value)
  names(args) <- input_id
  args$wait_ <- wait
  args$timeout_ <- timeout
  args$allow_no_input_binding_ <- allow_no_input_binding
  args$priority_ <- priority
  do.call(driver$set_inputs, args)
  invisible(driver)
}

e2e_upload_file <- function(
    driver,
    input_id,
    path,
    wait = TRUE,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  if (!file.exists(path)) {
    stop(sprintf("Upload file does not exist: %s", path), call. = FALSE)
  }
  args <- list(path)
  names(args) <- input_id
  args$wait_ <- wait
  args$timeout_ <- timeout
  do.call(driver$upload_file, args)
  invisible(driver)
}

e2e_input_id <- function(omic, step, input) {
  if (!omic %in% .E2E_VALID_OMICS) {
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  }

  step_id <- switch(
    step,
    import = switch(
      omic,
      proteomics = "setup_import",
      metabolomics = "import",
      lipidomics = "import"
    ),
    design = switch(
      omic,
      proteomics = "design_matrix",
      metabolomics = "design",
      lipidomics = "design"
    ),
    norm = switch(
      omic,
      proteomics = "normalization",
      metabolomics = "norm",
      lipidomics = "norm"
    ),
    da = switch(
      omic,
      proteomics = "differential_expression",
      metabolomics = "de",
      lipidomics = "de"
    ),
    summary = switch(
      omic,
      proteomics = "session_summary",
      metabolomics = "summary",
      lipidomics = "summary"
    ),
    stop(sprintf("Unsupported workflow step: %s", step), call. = FALSE)
  )

  paste0(omic, "_workflow-", step_id, "-", input)
}

e2e_lane_seed_path <- function(lane, fixture_root = .e2e_fixture_root()) {
  normalizePath(
    file.path(fixture_root, lane$fixture_dir, lane$seed_file),
    mustWork = TRUE
  )
}

e2e_lane_design_path <- function(lane, fixture_root = .e2e_fixture_root()) {
  normalizePath(
    file.path(fixture_root, lane$fixture_dir, "design.tsv"),
    mustWork = TRUE
  )
}

e2e_upload_lane_inputs <- function(driver, lane) {
  seed_path <- e2e_lane_seed_path(lane)
  omic <- lane$omic_type

  if (identical(omic, "proteomics")) {
    e2e_upload_file(
      driver,
      e2e_input_id("proteomics", "import", "search_results_standard"),
      seed_path
    )
    if (!is.null(lane$fasta_file)) {
      fasta_path <- normalizePath(
        file.path(.e2e_fixture_root(), lane$fixture_dir, lane$fasta_file),
        mustWork = TRUE
      )
      e2e_upload_file(
        driver,
        e2e_input_id("proteomics", "import", "fasta_file_standard"),
        fasta_path
      )
    }
    e2e_set_input(
      driver,
      e2e_input_id("proteomics", "import", "format_override"),
      lane$import_tool
    )
    return(invisible(driver))
  }

  if (identical(omic, "metabolomics")) {
    e2e_upload_file(
      driver,
      e2e_input_id("metabolomics", "import", "assay1_file_std"),
      seed_path
    )
    e2e_set_input(
      driver,
      e2e_input_id("metabolomics", "import", "vendor_format"),
      lane$import_tool
    )
    return(invisible(driver))
  }

  if (identical(omic, "lipidomics")) {
    e2e_upload_file(
      driver,
      e2e_input_id("lipidomics", "import", "assay1_file_std"),
      seed_path
    )
    e2e_set_input(
      driver,
      e2e_input_id("lipidomics", "import", "vendor_format"),
      lane$import_tool
    )
    return(invisible(driver))
  }

  stop(sprintf("Unsupported lane omic_type: %s", omic), call. = FALSE)
}

e2e_select_omics <- function(
    driver,
    omics,
    assert_digest = FALSE,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  bad_omics <- setdiff(omics, .E2E_VALID_OMICS)
  if (length(bad_omics) > 0L) {
    stop(
      sprintf("Unsupported omics requested: %s", paste(bad_omics, collapse = ", ")),
      call. = FALSE
    )
  }

  for (omic in omics) {
    e2e_click_testid(driver, paste0("tile-", omic), timeout = timeout)
  }

  if (isTRUE(assert_digest)) {
    e2e_assert_state_digest(
      driver,
      selected_omics = omics,
      timeout = timeout
    )
  }

  invisible(driver)
}

e2e_complete_project_setup <- function(
    driver,
    experiment_label,
    project_dir,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  if (!dir.exists(project_dir)) {
    stop(sprintf("project_dir does not exist: %s", project_dir), call. = FALSE)
  }

  e2e_click_testid(driver, "btn-start-analysis", timeout = timeout)
  e2e_set_input(driver, "experiment_label", experiment_label, timeout = timeout)
  e2e_set_input(driver, "project_dir", normalizePath(project_dir), timeout = timeout)
  e2e_click_testid(driver, "btn-confirm-setup", timeout = timeout)
  invisible(driver)
}

e2e_switch_workflow_tab <- function(driver, omic, tab_value) {
  input <- switch(
    omic,
    proteomics = "proteomics_workflow-workflow_tabs",
    metabolomics = "metabolomics_workflow-metabolomics_tabs",
    lipidomics = "lipidomics_workflow-lipidomics_tabs",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_set_input(driver, input, tab_value)
}

e2e_click_process_import <- function(driver, omic) {
  test_id <- switch(
    omic,
    proteomics = "prot-import-process",
    metabolomics = "metab-import-process",
    lipidomics = "lipid-import-process",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id)
}

e2e_click_norm_export <- function(driver, omic) {
  test_id <- switch(
    omic,
    proteomics = "prot-norm-export-session",
    metabolomics = "metab-norm-export-session",
    lipidomics = "lipid-norm-export-session",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id)
}

e2e_click_da_load_session <- function(driver, omic) {
  test_id <- switch(
    omic,
    proteomics = "prot-da-load-session",
    metabolomics = "metab-da-load-session",
    lipidomics = "lipid-da-load-session",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id)
}

e2e_click_da_run <- function(driver, omic) {
  test_id <- switch(
    omic,
    proteomics = "prot-da-run",
    metabolomics = "metab-da-run",
    lipidomics = "lipid-da-run",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id)
}

e2e_click_summary_action <- function(driver, omic, action) {
  prefix <- switch(
    omic,
    proteomics = "prot",
    metabolomics = "metab",
    lipidomics = "lipid",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  suffix <- switch(
    action,
    save_workflow_args = "summary-save-args",
    copy_publication = "summary-copy-publication",
    generate_report = "summary-generate-report",
    download_report = "summary-download-report",
    export_session = "summary-export-state",
    stop(sprintf("Unsupported summary action: %s", action), call. = FALSE)
  )
  e2e_click_testid(driver, paste(prefix, suffix, sep = "-"))
}

e2e_state_digest <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  e2e_wait_for_idle(driver, timeout = timeout)
  raw <- e2e_call_driver_method(
    driver,
    "get_value",
    output = "test_state_digest"
  )

  if (is.list(raw)) {
    return(raw)
  }
  if (is.null(raw) || !nzchar(raw)) {
    stop("test_state_digest output was empty", call. = FALSE)
  }

  jsonlite::fromJSON(raw, simplifyVector = FALSE)
}

e2e_assert_state_digest <- function(
    driver,
    selected_omics = NULL,
    initialized_omics = NULL,
    required_project_dir_keys = NULL,
    workflow_types = NULL,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  digest <- e2e_state_digest(driver, timeout = timeout)

  if (!is.null(selected_omics)) {
    testthat::expect_setequal(unlist(digest$selected_omics), selected_omics)
  }
  if (!is.null(initialized_omics)) {
    testthat::expect_setequal(unlist(digest$initialized_omics), initialized_omics)
  }
  if (!is.null(required_project_dir_keys)) {
    testthat::expect_true(all(required_project_dir_keys %in% unlist(digest$project_dir_keys)))
  }
  if (!is.null(workflow_types)) {
    for (omic in names(workflow_types)) {
      testthat::expect_identical(
        digest$workflow_type_per_omic[[omic]],
        workflow_types[[omic]]
      )
    }
  }

  invisible(digest)
}

e2e_collect_console_errors <- function(driver, allow = character()) {
  logs <- e2e_call_driver_method(driver, "get_logs", required = FALSE)
  if (is.null(logs) || length(logs) == 0L) {
    return(character())
  }

  log_text <- unlist(lapply(logs, function(item) {
    paste(unlist(item), collapse = " ")
  }), use.names = FALSE)

  error_pattern <- "error|severe|uncaught|typeerror|referenceerror|syntaxerror"
  errors <- log_text[grepl(error_pattern, log_text, ignore.case = TRUE)]
  if (length(allow) > 0L) {
    for (pattern in allow) {
      errors <- errors[!grepl(pattern, errors, ignore.case = TRUE)]
    }
  }

  errors
}

e2e_assert_no_console_errors <- function(driver, allow = character()) {
  errors <- e2e_collect_console_errors(driver, allow = allow)
  testthat::expect_equal(errors, character())
  invisible(driver)
}

e2e_collect_notifications <- function(driver) {
  script <- paste(
    "Array.from(document.querySelectorAll('.shiny-notification')).map(function(node) {",
    "return {",
    "text: (node.textContent || '').trim(),",
    "className: node.className || ''",
    "};",
    "}).filter(function(item) { return item.text.length > 0; })"
  )

  result <- e2e_call_driver_method(
    driver,
    "get_js",
    script = script,
    required = FALSE
  )
  if (is.null(result)) list() else result
}

e2e_assert_no_error_notifications <- function(driver) {
  notifications <- e2e_collect_notifications(driver)
  if (length(notifications) == 0L) {
    return(invisible(driver))
  }

  text <- unlist(lapply(notifications, function(item) {
    paste(unlist(item), collapse = " ")
  }), use.names = FALSE)
  errors <- text[grepl("error|failed|exception", text, ignore.case = TRUE)]
  testthat::expect_equal(errors, character())
  invisible(driver)
}

e2e_get_download <- function(
    driver,
    output_id,
    artifact_dir = e2e_case_artifact_dir(output_id),
    filename = NULL,
    ...
) {
  path <- e2e_call_driver_method(driver, "get_download", output_id, ...)
  if (is.null(filename)) {
    return(path)
  }

  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dest <- file.path(artifact_dir, filename)
  file.copy(path, dest, overwrite = TRUE)
  dest
}

e2e_assert_file_nonempty <- function(path) {
  testthat::expect_true(file.exists(path))
  testthat::expect_gt(file.info(path)$size, 0)
  invisible(path)
}

e2e_capture_failure_artifacts <- function(
    driver,
    case_id,
    reason = "failure",
    artifact_dir = e2e_case_artifact_dir(case_id)
) {
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }

  reason_path <- file.path(artifact_dir, "failure-reason.txt")
  writeLines(as.character(reason), reason_path)

  tryCatch(
    {
      html <- e2e_call_driver_method(driver, "get_html", "html", required = FALSE)
      if (!is.null(html)) {
        writeLines(as.character(html), file.path(artifact_dir, "dom.html"))
      }
    },
    error = function(e) NULL
  )

  tryCatch(
    {
      logs <- e2e_call_driver_method(driver, "get_logs", required = FALSE)
      if (!is.null(logs)) {
        jsonlite::write_json(
          logs,
          file.path(artifact_dir, "browser-logs.json"),
          auto_unbox = TRUE,
          pretty = TRUE
        )
      }
    },
    error = function(e) NULL
  )

  tryCatch(
    {
      e2e_call_driver_method(
        driver,
        "get_screenshot",
        file.path(artifact_dir, "screenshot.png"),
        selector = "viewport",
        required = FALSE
      )
    },
    error = function(e) NULL
  )

  invisible(artifact_dir)
}

e2e_with_failure_artifacts <- function(
    driver,
    case_id,
    expr,
    artifact_dir = e2e_case_artifact_dir(case_id)
) {
  tryCatch(
    force(expr),
    error = function(e) {
      e2e_capture_failure_artifacts(
        driver = driver,
        case_id = case_id,
        reason = conditionMessage(e),
        artifact_dir = artifact_dir
      )
      stop(e)
    }
  )
}

e2e_stop_driver <- function(driver) {
  e2e_call_driver_method(driver, "stop", required = FALSE)
  invisible(NULL)
}
