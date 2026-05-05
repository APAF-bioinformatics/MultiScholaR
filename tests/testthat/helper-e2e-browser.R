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

e2e_input_selector <- function(input_id) {
  if (length(input_id) != 1L || !nzchar(input_id)) {
    stop("input_id must be a non-empty character scalar", call. = FALSE)
  }
  escaped <- gsub("\"", "\\\\\"", input_id, fixed = TRUE)
  sprintf("[id=\"%s\"]", escaped)
}

e2e_launch_spec <- function(launch_surface = c("run_app", "app_dir")) {
  launch_surface <- match.arg(launch_surface)
  project_dir <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)

  switch(
    launch_surface,
    run_app = local({
      app_project_dir <- project_dir
      function() {
        setwd(app_project_dir)
        options(multischolar.test_mode = TRUE)
        pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
        run_app(test_mode = TRUE)
      }
    }),
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

e2e_wait_for_input_id <- function(
    driver,
    input_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  script <- sprintf(
    "document.getElementById(%s) !== null",
    jsonlite::toJSON(input_id, auto_unbox = TRUE)
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

e2e_click_input_id <- function(
    driver,
    input_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_wait_for_input_id(driver, input_id, timeout = timeout)
  script <- sprintf(
    "document.getElementById(%s).click(); true",
    jsonlite::toJSON(input_id, auto_unbox = TRUE)
  )
  e2e_call_driver_method(driver, "get_js", script = script)
  e2e_wait_for_idle(driver, timeout = timeout)
  invisible(driver)
}

e2e_trigger_action_input_id <- function(
    driver,
    input_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_wait_for_input_id(driver, input_id, timeout = timeout)
  script <- sprintf(
    paste(
      "(function(){",
      "var inputId = %s;",
      "var node = document.getElementById(inputId);",
      "if (!node) { return false; }",
      "node.click();",
      "if (window.Shiny && Shiny.setInputValue) {",
      "Shiny.setInputValue(inputId, Date.now() + Math.random(), { priority: 'event' });",
      "}",
      "return true;",
      "})()"
    ),
    jsonlite::toJSON(input_id, auto_unbox = TRUE)
  )
  e2e_call_driver_method(driver, "get_js", script = script)
  e2e_wait_for_idle(driver, timeout = timeout)
  invisible(driver)
}

e2e_click_optional_input_id <- function(
    driver,
    input_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  script <- sprintf(
    paste(
      "var node = document.getElementById(%s);",
      "if (node) { node.click(); true; } else { false; }"
    ),
    jsonlite::toJSON(input_id, auto_unbox = TRUE)
  )
  clicked <- e2e_call_driver_method(driver, "get_js", script = script, required = FALSE)
  if (isTRUE(clicked)) {
    e2e_wait_for_idle(driver, timeout = timeout)
  }
  invisible(clicked)
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

e2e_set_selectize_values <- function(
    driver,
    input_id,
    values,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    multiple = TRUE
) {
  input_js <- jsonlite::toJSON(input_id, auto_unbox = TRUE)
  values_js <- jsonlite::toJSON(as.character(values), auto_unbox = FALSE)
  multiple_js <- if (isTRUE(multiple)) "true" else "false"
  set_script <- sprintf(
    paste(
      "(function(){",
      "var inputId = %s;",
      "var values = %s.map(String);",
      "var multiple = %s;",
      "var payload = multiple ? values : (values.length > 0 ? values[0] : '');",
      "var el = document.getElementById(inputId);",
      "if (!el) { return false; }",
      "if (el.selectize) {",
      "values.forEach(function(value) {",
      "if (!Object.prototype.hasOwnProperty.call(el.selectize.options, value)) {",
      "el.selectize.addOption({ value: value, label: value, text: value });",
      "}",
      "});",
      "el.selectize.setValue(payload, false);",
      "el.selectize.refreshOptions(false);",
      "} else {",
      "Array.from(el.options).forEach(function(option) {",
      "option.selected = values.indexOf(String(option.value)) !== -1;",
      "});",
      "values.forEach(function(value) {",
      "if (!Array.from(el.options).some(function(option) { return String(option.value) === value; })) {",
      "var option = new Option(value, value, true, true);",
      "el.add(option);",
      "}",
      "});",
      "el.dispatchEvent(new Event('change', { bubbles: true }));",
      "}",
      "if (window.Shiny && Shiny.setInputValue) {",
      "Shiny.setInputValue(inputId, payload, { priority: 'event' });",
      "}",
      "return true;",
      "})()"
    ),
    input_js,
    values_js,
    multiple_js
  )
  value_script <- sprintf(
    paste(
      "(function(){",
      "var inputId = %s;",
      "var expected = %s.map(String);",
      "var el = document.getElementById(inputId);",
      "if (!el) { return false; }",
      "var actual = [];",
      "if (el.selectize) {",
      "actual = el.selectize.getValue();",
      "actual = Array.isArray(actual) ? actual.map(String) : [String(actual)];",
      "} else {",
      "actual = Array.from(el.selectedOptions).map(function(option) { return String(option.value); });",
      "}",
      "return actual.length === expected.length && expected.every(function(value) {",
      "return actual.indexOf(value) !== -1;",
      "});",
      "})()"
    ),
    input_js,
    values_js
  )

  e2e_wait_for_input_id(driver, input_id, timeout = timeout)
  e2e_call_driver_method(driver, "get_js", script = set_script)
  e2e_wait_for_js(driver, script = value_script, timeout = timeout)
  e2e_wait_for_idle(driver, timeout = timeout)
  invisible(driver)
}

e2e_selectize_option_values <- function(
    driver,
    input_id,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  input_js <- jsonlite::toJSON(input_id, auto_unbox = TRUE)
  script <- sprintf(
    paste(
      "(function(){",
      "var el = document.getElementById(%s);",
      "if (!el) { return []; }",
      "if (el.selectize) { return Object.keys(el.selectize.options || {}); }",
      "return Array.from(el.options).map(function(option) { return String(option.value); });",
      "})()"
    ),
    input_js
  )

  e2e_wait_for_input_id(driver, input_id, timeout = timeout)
  result <- e2e_call_driver_method(driver, "get_js", script = script, required = FALSE)
  as.character(unlist(result, use.names = FALSE))
}

e2e_match_available_values <- function(values, available) {
  values <- as.character(values)
  available <- as.character(available)
  if (length(available) == 0L) {
    return(values)
  }

  vapply(values, function(value) {
    exact <- available[available == value]
    if (length(exact) > 0L) {
      return(exact[[1L]])
    }

    case_insensitive <- available[tolower(available) == tolower(value)]
    if (length(case_insensitive) > 0L) {
      return(case_insensitive[[1L]])
    }

    suffix_match <- available[endsWith(tolower(available), paste0("_", tolower(value)))]
    if (length(suffix_match) > 0L) {
      return(suffix_match[[1L]])
    }

    value
  }, character(1), USE.NAMES = FALSE)
}

e2e_set_input_and_idle <- function(
    driver,
    input_id,
    value,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    priority = c("input", "event")
) {
  priority <- match.arg(priority)
  e2e_set_input(
    driver = driver,
    input_id = input_id,
    value = value,
    wait = FALSE,
    timeout = timeout,
    priority = priority
  )
  e2e_wait_for_idle(driver, timeout = timeout)
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
    qc = switch(
      omic,
      proteomics = "quality_control",
      metabolomics = "qc",
      lipidomics = "qc"
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
    enrich = switch(
      omic,
      proteomics = "enrichment_analysis",
      metabolomics = "enrichment",
      lipidomics = "enrichment"
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

e2e_project_base_dir <- function(project_dir, experiment_label) {
  file.path(normalizePath(project_dir, mustWork = TRUE), experiment_label)
}

e2e_proteomics_project_paths <- function(project_dir, experiment_label) {
  base_dir <- e2e_project_base_dir(project_dir, experiment_label)
  list(
    base_dir = base_dir,
    source_dir = file.path(base_dir, "scripts", "proteomics"),
    data_dir = file.path(base_dir, "data", "proteomics"),
    results_dir = file.path(base_dir, "results", "proteomics"),
    results_summary_dir = file.path(base_dir, "results_summary", "proteomics"),
    da_output_dir = file.path(base_dir, "results", "proteomics", "da_proteins"),
    pathway_dir = file.path(base_dir, "results", "proteomics", "pathway_enrichment")
  )
}

e2e_metabolomics_project_paths <- function(project_dir, experiment_label) {
  base_dir <- e2e_project_base_dir(project_dir, experiment_label)
  list(
    base_dir = base_dir,
    source_dir = file.path(base_dir, "scripts", "metabolomics"),
    data_dir = file.path(base_dir, "data", "metabolomics"),
    results_dir = file.path(base_dir, "results", "metabolomics"),
    results_summary_dir = file.path(base_dir, "results_summary", "metabolomics"),
    da_output_dir = file.path(base_dir, "results", "metabolomics", "da_metabolites"),
    publication_graphs_dir = file.path(base_dir, "results", "metabolomics", "publication_graphs"),
    metabolite_qc_dir = file.path(base_dir, "results", "metabolomics", "metabolite_qc"),
    pathway_dir = file.path(base_dir, "results", "metabolomics", "pathway_enrichment")
  )
}

e2e_assert_project_dirs_exist <- function(
    paths,
    required = c(
      "base_dir",
      "source_dir",
      "data_dir",
      "results_dir",
      "results_summary_dir",
      "da_output_dir"
    )
) {
  for (name in required) {
    testthat::expect_true(
      dir.exists(paths[[name]]),
      info = sprintf("Expected project directory '%s' to exist at %s", name, paths[[name]])
    )
  }
  invisible(paths)
}

e2e_seed_report_template <- function(lane, project_base_dir) {
  source <- normalizePath(
    file.path(.e2e_fixture_root(), "report_templates", lane$report_template),
    mustWork = TRUE
  )
  dest_dir <- file.path(project_base_dir, "scripts", lane$omic_type)
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  }
  dest <- file.path(dest_dir, lane$report_template)
  file.copy(source, dest, overwrite = TRUE)
  e2e_assert_file_nonempty(dest)
}

e2e_assert_prot_filtered_session_artifacts <- function(project_base_dir) {
  source_dir <- file.path(project_base_dir, "scripts", "proteomics")
  latest <- file.path(source_dir, "filtered_session_data_latest.rds")
  summary <- file.path(source_dir, "filtered_session_summary.txt")

  e2e_assert_file_nonempty(latest)
  e2e_assert_file_nonempty(summary)

  session_data <- readRDS(latest)
  testthat::expect_true("r6_current_state_name" %in% names(session_data))
  testthat::expect_true("r6_complete_states" %in% names(session_data))
  testthat::expect_true("r6_state_history" %in% names(session_data))
  testthat::expect_gt(length(session_data$r6_complete_states), 0)
  testthat::expect_true(session_data$r6_current_state_name %in% names(session_data$r6_complete_states))

  invisible(list(latest = latest, summary = summary, session_data = session_data))
}

e2e_assert_metab_filtered_session_artifacts <- function(project_base_dir, expected_assays = NULL) {
  source_dir <- file.path(project_base_dir, "scripts", "metabolomics")
  latest <- file.path(source_dir, "metab_filtered_session_data_latest.rds")
  summary <- file.path(source_dir, "metab_filtered_session_summary.txt")

  e2e_assert_file_nonempty(latest)
  e2e_assert_file_nonempty(summary)

  session_data <- readRDS(latest)
  testthat::expect_identical(session_data$omic_type, "metabolomics")
  testthat::expect_true("r6_current_state_name" %in% names(session_data))
  testthat::expect_true("current_s4_object" %in% names(session_data))
  testthat::expect_true("assay_names" %in% names(session_data))
  testthat::expect_true(isTRUE(session_data$normalization_complete))
  testthat::expect_true(isTRUE(session_data$correlation_filtering_complete))
  testthat::expect_identical(session_data$normalization_method, "none")
  testthat::expect_identical(session_data$ruv_mode, "skip")
  testthat::expect_false(isTRUE(session_data$itsd_applied))
  if (!is.null(expected_assays)) {
    testthat::expect_setequal(session_data$assay_names, unlist(expected_assays))
  }

  invisible(list(latest = latest, summary = summary, session_data = session_data))
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
    metab_import <- function(input) e2e_input_id("metabolomics", "import", input)
    assay_names <- unlist(lane$assays)
    assay_names <- assay_names[!is.na(assay_names) & nzchar(assay_names)]

    if (length(assay_names) > 0L) {
      e2e_set_input_and_idle(driver, metab_import("assay1_name"), assay_names[[1]])
    }
    e2e_set_input_and_idle(
      driver,
      metab_import("vendor_format"),
      lane$import_tool
    )
    e2e_upload_file(
      driver,
      metab_import("assay1_file_std"),
      seed_path
    )

    if (!is.null(lane$assay2_file) && length(assay_names) >= 2L) {
      assay2_path <- normalizePath(
        file.path(.e2e_fixture_root(), lane$fixture_dir, lane$assay2_file),
        mustWork = TRUE
      )
      e2e_set_input_and_idle(driver, metab_import("assay2_name"), assay_names[[2]])
      e2e_upload_file(
        driver,
        metab_import("assay2_file_std"),
        assay2_path
      )
    }

    if (identical(lane$import_tool, "custom")) {
      e2e_set_input_and_idle(driver, metab_import("metabolite_id_col_custom"), "Feature.Name")
      e2e_set_input_and_idle(driver, metab_import("annotation_col_custom"), "")
      e2e_set_input_and_idle(driver, metab_import("sample_cols_pattern"), "^(WT|KO)_")
      e2e_set_input_and_idle(driver, metab_import("is_pattern"), "^ITSD|_d[0-9]+$")
      e2e_set_input_and_idle(driver, metab_import("sanitize_names"), FALSE)
    }

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

e2e_switch_workflow_tab <- function(driver, omic, tab_value, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  input <- switch(
    omic,
    proteomics = "proteomics_workflow-workflow_tabs",
    metabolomics = "metabolomics_workflow-metabolomics_tabs",
    lipidomics = "lipidomics_workflow-lipidomics_tabs",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  ui_value <- switch(
    omic,
    proteomics = c(
      setup_import = "setup",
      design_matrix = "design",
      quality_control = "qc",
      normalization = "normalization",
      differential_expression = "da",
      enrichment_analysis = "enrichment",
      session_summary = "session_summary"
    )[[tab_value]] %||% tab_value,
    tab_value
  )
  input_js <- jsonlite::toJSON(input, auto_unbox = TRUE)
  value_js <- jsonlite::toJSON(ui_value, auto_unbox = TRUE)
  exists_script <- sprintf(
    "(function(){var tabs=document.getElementById(%s);return !!(tabs && tabs.querySelector('a[data-value=\"' + %s + '\"]'));})()",
    input_js,
    value_js
  )
  click_script <- sprintf(
    "(function(){var tabs=document.getElementById(%s);var tab=tabs && tabs.querySelector('a[data-value=\"' + %s + '\"]');if(!tab){return false;}tab.click();return true;})()",
    input_js,
    value_js
  )
  e2e_wait_for_js(driver, script = exists_script, timeout = timeout)
  e2e_call_driver_method(driver, "get_js", script = click_script, required = FALSE)
  e2e_wait_for_idle(driver, timeout = timeout)
  invisible(driver)
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

e2e_click_da_load_session <- function(driver, omic, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  test_id <- switch(
    omic,
    proteomics = "prot-da-load-session",
    metabolomics = "metab-da-load-session",
    lipidomics = "lipid-da-load-session",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id, timeout = timeout)
}

e2e_click_da_run <- function(driver, omic, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  if (identical(omic, "metabolomics")) {
    return(e2e_trigger_action_input_id(
      driver,
      e2e_input_id("metabolomics", "da", "run_da_analysis"),
      timeout = timeout
    ))
  }

  test_id <- switch(
    omic,
    proteomics = "prot-da-run",
    metabolomics = "metab-da-run",
    lipidomics = "lipid-da-run",
    stop(sprintf("Unsupported omic: %s", omic), call. = FALSE)
  )
  e2e_click_testid(driver, test_id, timeout = timeout)
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

e2e_read_lane_design <- function(lane) {
  utils::read.table(
    e2e_lane_design_path(lane),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
}

e2e_complete_prot_design_from_manifest <- function(
    driver,
    lane,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  design <- e2e_read_lane_design(lane)
  groups <- unique(design$group)
  builder <- function(input) e2e_input_id("proteomics", "design", paste0("builder-", input))

  e2e_switch_workflow_tab(driver, "proteomics", "design_matrix")
  e2e_wait_for_selector(driver, "prot-tab-design", timeout = timeout)

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Factors", timeout = timeout)
  for (group in groups) {
    e2e_set_input_and_idle(driver, builder("new_factor"), group, timeout = timeout)
    e2e_click_input_id(driver, builder("add_factor"), timeout = timeout)
  }

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Assign Metadata", timeout = timeout)
  available_runs <- e2e_selectize_option_values(driver, builder("selected_runs"), timeout = timeout)
  for (group in groups) {
    group_design <- design[design$group == group, , drop = FALSE]
    selected_runs <- e2e_match_available_values(group_design$sample, available_runs)
    e2e_set_selectize_values(driver, builder("selected_runs"), selected_runs, timeout = timeout)
    e2e_wait_for_input_id(driver, builder("replicate_start"), timeout = timeout)
    e2e_set_selectize_values(
      driver,
      builder("factor1_select"),
      group,
      timeout = timeout,
      multiple = FALSE
    )
    e2e_click_input_id(driver, builder("assign_metadata"), timeout = timeout)
  }

  contrast <- lane$expected_contrasts[[1]]
  contrast_parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  if (length(contrast_parts) != 2L) {
    stop(sprintf("Unsupported expected_contrasts format: %s", contrast), call. = FALSE)
  }

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Contrasts", timeout = timeout)
  contrast_choices <- e2e_selectize_option_values(driver, builder("contrast_group1"), timeout = timeout)
  contrast_group1 <- e2e_match_available_values(contrast_parts[[1]], contrast_choices)
  contrast_group2 <- e2e_match_available_values(contrast_parts[[2]], contrast_choices)
  e2e_set_selectize_values(
    driver,
    builder("contrast_group1"),
    contrast_group1,
    timeout = timeout,
    multiple = FALSE
  )
  e2e_set_selectize_values(
    driver,
    builder("contrast_group2"),
    contrast_group2,
    timeout = timeout,
    multiple = FALSE
  )
  e2e_click_input_id(driver, builder("add_contrast"), timeout = timeout)
  e2e_click_input_id(driver, builder("save_results"), timeout = timeout)
  invisible(driver)
}

e2e_complete_metab_design_from_manifest <- function(
    driver,
    lane,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  design <- e2e_read_lane_design(lane)
  groups <- unique(design$group)
  builder <- function(input) e2e_input_id("metabolomics", "design", paste0("builder-", input))

  e2e_switch_workflow_tab(driver, "metabolomics", "design", timeout = timeout)
  e2e_wait_for_selector(driver, "metab-tab-design", timeout = timeout)

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Factors", timeout = timeout)
  for (group in groups) {
    e2e_set_input_and_idle(driver, builder("new_factor"), group, timeout = timeout)
    e2e_click_input_id(driver, builder("add_factor"), timeout = timeout)
  }

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Assign Metadata", timeout = timeout)
  available_runs <- e2e_selectize_option_values(driver, builder("selected_runs"), timeout = timeout)
  for (group in groups) {
    group_design <- design[design$group == group, , drop = FALSE]
    selected_runs <- e2e_match_available_values(group_design$sample, available_runs)
    e2e_set_selectize_values(driver, builder("selected_runs"), selected_runs, timeout = timeout)
    e2e_wait_for_input_id(driver, builder("replicate_start"), timeout = timeout)
    e2e_set_selectize_values(
      driver,
      builder("factor1_select"),
      group,
      timeout = timeout,
      multiple = FALSE
    )
    e2e_click_input_id(driver, builder("assign_metadata"), timeout = timeout)
  }

  contrast <- lane$expected_contrasts[[1]]
  contrast_parts <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  if (length(contrast_parts) != 2L) {
    stop(sprintf("Unsupported expected_contrasts format: %s", contrast), call. = FALSE)
  }

  e2e_set_input_and_idle(driver, builder("main_tabs"), "Contrasts", timeout = timeout)
  contrast_choices <- e2e_selectize_option_values(driver, builder("contrast_group1"), timeout = timeout)
  contrast_group1 <- e2e_match_available_values(contrast_parts[[1]], contrast_choices)
  contrast_group2 <- e2e_match_available_values(contrast_parts[[2]], contrast_choices)
  e2e_set_selectize_values(
    driver,
    builder("contrast_group1"),
    contrast_group1,
    timeout = timeout,
    multiple = FALSE
  )
  e2e_set_selectize_values(
    driver,
    builder("contrast_group2"),
    contrast_group2,
    timeout = timeout,
    multiple = FALSE
  )
  e2e_click_input_id(driver, builder("add_contrast"), timeout = timeout)
  e2e_click_input_id(driver, builder("save_results"), timeout = timeout)
  invisible(driver)
}

e2e_run_prot_dia_qc <- function(
    driver,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    rollup_method = "iq"
) {
  qc <- function(input) e2e_input_id("proteomics", "qc", input)

  peptide_steps <- list(
    list(tab = "Q-Value Filter", action = "peptide_qc-qvalue_filter-apply_qvalue_filter"),
    list(tab = "Precursor Rollup", action = "peptide_qc-rollup-apply_rollup"),
    list(tab = "Intensity Filter", action = "peptide_qc-intensity_filter-apply_intensity_filter"),
    list(tab = "Protein Peptides", action = "peptide_qc-protein_peptide_filter-apply_protein_peptide_filter"),
    list(tab = "Sample Quality", action = "peptide_qc-sample_filter-apply_sample_filter"),
    list(tab = "Replicate Filter", action = "peptide_qc-replicate_filter-apply_replicate_filter"),
    list(tab = "Imputation", action = "peptide_qc-imputation-apply_imputation")
  )

  protein_steps <- list(
    list(tab = "IQ Protein Rollup", action = "protein_qc-rollup-apply_iq_rollup"),
    list(tab = "Accession Cleanup", action = "protein_qc-cleanup-apply_accession_cleanup"),
    list(tab = "Protein Intensity Filter", action = "protein_qc-intensity_filter-apply_protein_intensity_filter"),
    list(tab = "Duplicate Removal", action = "protein_qc-duplicate_removal-apply_duplicate_removal"),
    list(tab = "Protein Replicate Filter", action = "protein_qc-replicate_filter-apply_protein_replicate_filter")
  )

  e2e_switch_workflow_tab(driver, "proteomics", "quality_control")
  e2e_wait_for_selector(driver, "prot-tab-qc", timeout = timeout)
  e2e_set_input_and_idle(driver, qc("qc_tabs_lfq"), "Peptide QC", timeout = timeout)
  for (step in peptide_steps) {
    e2e_set_input_and_idle(driver, qc("peptide_qc-peptide_filter_tabs"), step$tab, timeout = timeout)
    e2e_click_input_id(driver, qc(step$action), timeout = timeout)
  }

  e2e_set_input_and_idle(driver, qc("qc_tabs_lfq"), "Protein QC", timeout = timeout)
  for (step in protein_steps) {
    e2e_set_input_and_idle(driver, qc("protein_qc-protein_filter_tabs"), step$tab, timeout = timeout)
    if (identical(step$action, "protein_qc-rollup-apply_iq_rollup")) {
      e2e_set_input_and_idle(
        driver,
        qc("protein_qc-rollup-rollup_method"),
        rollup_method,
        timeout = timeout
      )
    }
    click_fn <- if (identical(step$action, "protein_qc-replicate_filter-apply_protein_replicate_filter")) {
      e2e_trigger_action_input_id
    } else {
      e2e_click_input_id
    }
    click_fn(driver, qc(step$action), timeout = timeout)
  }

  invisible(driver)
}

e2e_run_prot_dia_limpa_qc <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  e2e_run_prot_dia_qc(driver, timeout = timeout, rollup_method = "limpa")
}

e2e_run_prot_normalization_export <- function(
    driver,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_switch_workflow_tab(driver, "proteomics", "normalization")
  e2e_wait_for_selector(driver, "prot-tab-normalization", timeout = timeout)
  e2e_set_input_and_idle(driver, e2e_input_id("proteomics", "norm", "ruv_mode"), "skip", timeout = timeout)
  e2e_click_testid(driver, "prot-norm-run", timeout = timeout)
  e2e_set_input_and_idle(
    driver,
    e2e_input_id("proteomics", "norm", "norm_qc_tabs"),
    "Correlation Filtering",
    timeout = timeout
  )
  e2e_click_input_id(driver, e2e_input_id("proteomics", "norm", "skip_correlation_filter"), timeout = timeout)
  e2e_click_norm_export(driver, "proteomics")
  invisible(driver)
}

e2e_run_prot_dia_normalization_export <- e2e_run_prot_normalization_export

e2e_run_metab_qc_finalize <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  qc <- function(input) e2e_input_id("metabolomics", "qc", input)

  e2e_switch_workflow_tab(driver, "metabolomics", "qc", timeout = timeout)
  e2e_wait_for_selector(driver, "metab-tab-qc", timeout = timeout)
  e2e_set_input_and_idle(driver, qc("metab_qc_tabs"), "Finalize QC", timeout = timeout)
  e2e_click_input_id(driver, qc("s4_finalize-finalize_qc"), timeout = timeout)
  invisible(driver)
}

e2e_run_metab_normalization_export <- function(
    driver,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  norm <- function(input) e2e_input_id("metabolomics", "norm", input)

  e2e_switch_workflow_tab(driver, "metabolomics", "norm", timeout = timeout)
  e2e_wait_for_selector(driver, "metab-tab-norm", timeout = timeout)
  e2e_set_input_and_idle(driver, norm("apply_itsd"), FALSE, timeout = timeout)
  e2e_set_input_and_idle(driver, norm("norm_method"), "none", timeout = timeout)
  e2e_set_input_and_idle(driver, norm("ruv_mode"), "skip", timeout = timeout)
  e2e_click_testid(driver, "metab-norm-run", timeout = timeout)
  e2e_set_input_and_idle(driver, norm("norm_qc_tabs"), "Correlation Filtering", timeout = timeout)
  e2e_click_input_id(driver, norm("skip_correlation_filter"), timeout = timeout)
  e2e_click_norm_export(driver, "metabolomics")
  invisible(driver)
}

e2e_run_metab_da <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  e2e_switch_workflow_tab(driver, "metabolomics", "de", timeout = timeout)
  e2e_wait_for_selector(driver, "metab-tab-de", timeout = timeout)
  e2e_click_da_load_session(driver, "metabolomics", timeout = timeout)
  e2e_click_da_run(driver, "metabolomics", timeout = timeout)
  invisible(driver)
}

e2e_run_metab_summary_report <- function(
    driver,
    lane,
    project_base_dir,
    experiment_label,
    case_id,
    description,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_switch_workflow_tab(driver, "metabolomics", "summary", timeout = timeout)
  e2e_wait_for_selector(driver, "metab-tab-summary", timeout = timeout)
  e2e_set_input_and_idle(
    driver,
    e2e_input_id("metabolomics", "summary", "experiment_label"),
    experiment_label,
    timeout = timeout
  )
  e2e_set_input_and_idle(
    driver,
    e2e_input_id("metabolomics", "summary", "description"),
    description,
    timeout = timeout
  )
  e2e_seed_report_template(lane, project_base_dir)
  e2e_click_summary_action(driver, "metabolomics", "save_workflow_args")
  e2e_click_summary_action(driver, "metabolomics", "copy_publication")
  e2e_click_summary_action(driver, "metabolomics", "generate_report")
  e2e_wait_for_selector(driver, "metab-summary-download-report", timeout = timeout)
  report_download <- e2e_get_download(
    driver,
    e2e_input_id("metabolomics", "summary", "download_report"),
    artifact_dir = e2e_case_artifact_dir(case_id),
    filename = paste0("downloaded-", lane$report_template)
  )
  e2e_assert_file_nonempty(report_download)
  e2e_click_summary_action(driver, "metabolomics", "export_session")

  invisible(report_download)
}

e2e_run_prot_da <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  e2e_switch_workflow_tab(driver, "proteomics", "differential_expression")
  e2e_wait_for_selector(driver, "prot-tab-da", timeout = timeout)
  e2e_click_da_load_session(driver, "proteomics")
  e2e_click_da_run(driver, "proteomics")
  e2e_click_optional_input_id(
    driver,
    e2e_input_id("proteomics", "da", "acknowledge_qvalue_warning"),
    timeout = timeout
  )
  invisible(driver)
}

e2e_wait_for_output_text <- function(
    driver,
    output_id,
    pattern,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    fixed = TRUE
) {
  deadline <- Sys.time() + (timeout / 1000)
  last_value <- NULL

  repeat {
    e2e_wait_for_idle(driver, timeout = timeout)
    last_value <- e2e_call_driver_method(
      driver,
      "get_value",
      output = output_id,
      required = FALSE
    )
    if (!is.null(last_value) && !is.list(last_value)) {
      text <- paste(as.character(last_value), collapse = "\n")
      if (grepl(pattern, text, fixed = fixed)) {
        return(text)
      }
    }
    if (Sys.time() >= deadline) {
      break
    }
    Sys.sleep(0.25)
  }

  stop(
    sprintf("Output %s did not match %s. Last value: %s", output_id, pattern, paste(last_value, collapse = "\n")),
    call. = FALSE
  )
}

e2e_run_prot_enrichment_backend <- function(
    driver,
    backend = c("gprofiler2", "clusterprofiler"),
    organism_taxid = NULL,
    expected_contrast = "KO_vs_WT",
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  backend <- match.arg(backend)
  enrich <- function(input) e2e_input_id("proteomics", "enrich", input)
  organism_taxid <- organism_taxid %||% if (identical(backend, "gprofiler2")) "9606" else "999999"

  e2e_switch_workflow_tab(driver, "proteomics", "enrichment_analysis", timeout = timeout)
  e2e_wait_for_selector(driver, "prot-tab-enrichment", timeout = timeout)
  e2e_set_input_and_idle(driver, enrich("organism_taxid"), organism_taxid, timeout = timeout)
  e2e_set_input_and_idle(driver, enrich("enrichment_method_tabs"), backend, timeout = timeout)
  contrast_input_js <- jsonlite::toJSON(enrich("selected_contrast"), auto_unbox = TRUE)
  contrast_value_js <- jsonlite::toJSON(expected_contrast, auto_unbox = TRUE)
  e2e_wait_for_js(
    driver,
    script = sprintf(
      "(function(){var el=document.getElementById(%s);return !!(el && Array.from(el.options || []).some(function(o){return o.value === %s;}));})()",
      contrast_input_js,
      contrast_value_js
    ),
    timeout = timeout
  )
  e2e_set_input_and_idle(driver, enrich("selected_contrast"), expected_contrast, timeout = timeout)
  e2e_set_input_and_idle(driver, enrich("up_cutoff"), 0, timeout = timeout)
  e2e_set_input_and_idle(driver, enrich("down_cutoff"), 0, timeout = timeout)
  e2e_set_input_and_idle(driver, enrich("q_cutoff"), 0.05, timeout = timeout)
  e2e_click_input_id(driver, enrich("run_enrichment_analysis"), timeout = timeout)
  status_text <- e2e_wait_for_output_text(
    driver,
    enrich("enrichment_status"),
    sprintf("Method: %s", backend),
    timeout = timeout
  )

  invisible(list(status_text = status_text, backend = backend, organism_taxid = organism_taxid))
}

e2e_run_prot_summary_report <- function(
    driver,
    lane,
    project_base_dir,
    experiment_label,
    case_id = "E2E-004-proteomics-dia-report",
    description = "E2E-004 canonical DIA GUI workflow",
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_switch_workflow_tab(driver, "proteomics", "session_summary")
  e2e_wait_for_selector(driver, "prot-tab-session-summary", timeout = timeout)
  e2e_set_input_and_idle(
    driver,
    e2e_input_id("proteomics", "summary", "experiment_label"),
    experiment_label,
    timeout = timeout
  )
  e2e_set_input_and_idle(
    driver,
    e2e_input_id("proteomics", "summary", "description"),
    description,
    timeout = timeout
  )
  e2e_seed_report_template(lane, project_base_dir)
  e2e_click_summary_action(driver, "proteomics", "save_workflow_args")
  e2e_click_summary_action(driver, "proteomics", "copy_publication")
  e2e_click_summary_action(driver, "proteomics", "generate_report")
  e2e_wait_for_selector(driver, "prot-summary-download-report", timeout = timeout)
  report_download <- e2e_get_download(
    driver,
    e2e_input_id("proteomics", "summary", "download_report"),
    artifact_dir = e2e_case_artifact_dir(case_id),
    filename = paste0("downloaded-", lane$report_template)
  )
  e2e_assert_file_nonempty(report_download)
  e2e_click_summary_action(driver, "proteomics", "export_session")

  invisible(report_download)
}

e2e_run_prot_da_summary_report <- function(
    driver,
    lane,
    project_base_dir,
    experiment_label,
    case_id = "E2E-004-proteomics-dia-report",
    description = "E2E-004 canonical DIA GUI workflow",
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT
) {
  e2e_run_prot_da(driver, timeout = timeout)
  e2e_run_prot_summary_report(
    driver = driver,
    lane = lane,
    project_base_dir = project_base_dir,
    experiment_label = experiment_label,
    case_id = case_id,
    description = description,
    timeout = timeout
  )
}

e2e_state_digest <- function(driver, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  deadline <- Sys.time() + (timeout / 1000)
  raw <- NULL

  repeat {
    e2e_wait_for_idle(driver, timeout = timeout)
    raw <- e2e_call_driver_method(
      driver,
      "get_value",
      output = "test_state_digest"
    )

    if (is.list(raw)) {
      return(raw)
    }
    if (!is.null(raw) && nzchar(raw)) {
      return(jsonlite::fromJSON(raw, simplifyVector = FALSE))
    }
    if (Sys.time() >= deadline) {
      break
    }
    Sys.sleep(0.25)
  }

  stop("test_state_digest output was empty", call. = FALSE)
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

e2e_digest_step_status <- function(digest, omic, step) {
  statuses <- digest$step_status_per_omic[[omic]]
  if (is.null(statuses)) {
    return(NULL)
  }
  statuses[[step]]
}

e2e_assert_step_statuses <- function(driver, omic, expected, timeout = .E2E_BROWSER_DEFAULT_TIMEOUT) {
  digest <- e2e_state_digest(driver, timeout = timeout)
  for (step in names(expected)) {
    testthat::expect_identical(
      e2e_digest_step_status(digest, omic, step),
      expected[[step]],
      info = sprintf("Expected %s/%s to be %s", omic, step, expected[[step]])
    )
  }
  invisible(digest)
}

e2e_wait_for_step_status <- function(
    driver,
    omic,
    step,
    status,
    timeout = .E2E_BROWSER_DEFAULT_TIMEOUT,
    interval = 0.25
) {
  deadline <- Sys.time() + (timeout / 1000)
  last_status <- NULL

  repeat {
    digest <- tryCatch(e2e_state_digest(driver, timeout = timeout), error = function(e) NULL)
    if (!is.null(digest)) {
      last_status <- e2e_digest_step_status(digest, omic, step)
      if (identical(last_status, status)) {
        return(invisible(digest))
      }
    }
    if (Sys.time() >= deadline) {
      break
    }
    Sys.sleep(interval)
  }

  testthat::fail(sprintf(
    "Timed out waiting for %s/%s to become %s; last status was %s",
    omic,
    step,
    status,
    last_status %||% "<NULL>"
  ))
}

e2e_collect_console_errors <- function(driver, allow = character()) {
  logs <- e2e_call_driver_method(driver, "get_logs", required = FALSE)
  if (is.null(logs) || length(logs) == 0L) {
    return(character())
  }

  log_entries <- lapply(logs, function(item) {
    if (is.list(item)) {
      level <- unlist(item$level %||% "", use.names = FALSE)
      level <- if (length(level) > 0L) tolower(as.character(level[[1L]])) else ""
      list(
        level = level,
        text = paste(unlist(item), collapse = " ")
      )
    } else {
      list(
        level = "",
        text = as.character(item)
      )
    }
  })

  javascript_error_pattern <- "\\b(uncaught|typeerror|referenceerror|syntaxerror|rangeerror|securityerror)\\b"
  errors <- unlist(lapply(log_entries, function(item) {
    high_severity <- item$level %in% c("error", "severe")
    javascript_exception <- any(grepl(javascript_error_pattern, item$text, ignore.case = TRUE))
    if (high_severity || javascript_exception) item$text else character()
  }), use.names = FALSE)

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
      if (inherits(e, "expectation")) {
        stop(conditionMessage(e), call. = FALSE)
      }
      stop(e)
    }
  )
}

e2e_stop_driver <- function(driver) {
  e2e_call_driver_method(driver, "stop", required = FALSE)
  invisible(NULL)
}
