library(testthat)

# --- is_test_mode ---

test_that("is_test_mode() returns FALSE by default", {
  withr::with_options(
    list(multischolar.test_mode = NULL)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , expect_false(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns TRUE when option is TRUE", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , expect_true(is_test_mode())
  )
})

test_that("is_test_mode() returns FALSE when option is FALSE", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , expect_false(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns TRUE when env var is 'true'", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "true")
      , expect_true(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns TRUE when env var is 'TRUE' (case insensitive)", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "TRUE")
      , expect_true(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns TRUE when env var is '1'", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "1")
      , expect_true(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns FALSE when env var is 'false'", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "false")
      , expect_false(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns TRUE when env var is 'True' (mixed case)", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "True")
      , expect_true(is_test_mode())
    )
  )
})

test_that("is_test_mode() returns FALSE when env var is '0'", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "0")
      , expect_false(is_test_mode())
    )
  )
})

test_that("is_test_mode() option TRUE takes precedence over env var 'false'", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "false")
      , expect_true(is_test_mode())
    )
  )
})

test_that("is_test_mode() is not cached — withr::with_options toggle works", {
  # Verify reading fresh each call so toggling in the same R session works
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , expect_false(is_test_mode())
  )
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , expect_true(is_test_mode())
  )
})

# --- testid ---

test_that("testid() appends data-testid attribute to a shiny tag", {
  tag <- shiny::div("hello")
  result <- testid(tag, "my-test-id")
  html <- as.character(result)
  expect_match(html, 'data-testid="my-test-id"')
})

test_that("testid() preserves the original tag type and children", {
  tag <- shiny::span(shiny::strong("text"))
  result <- testid(tag, "span-id")
  html <- as.character(result)
  expect_match(html, "<span")
  expect_match(html, "data-testid=\"span-id\"")
})

test_that("testid() works on actionButton (btn-start-analysis pattern)", {
  btn <- shiny::actionButton("start_analysis", "Start Analysis", class = "btn-primary")
  result <- testid(btn, "btn-start-analysis")
  html <- as.character(result)
  expect_match(html, 'data-testid="btn-start-analysis"')
  expect_match(html, 'id="start_analysis"')
})

test_that("testid() works on textInput (input-experiment-label pattern)", {
  inp <- shiny::textInput("experiment_label", "Experiment Label:", value = "test")
  result <- testid(inp, "input-experiment-label")
  html <- as.character(result)
  expect_match(html, 'data-testid="input-experiment-label"')
  expect_match(html, 'id="experiment_label"')
})

test_that("testid() works on a nested div (omic tile pattern)", {
  tile <- shiny::div(
    id = "proteomics_box"
    , class = "omic-selection-box"
    , shiny::h3("Proteomics")
  )
  result <- testid(tile, "tile-proteomics")
  html <- as.character(result)
  expect_match(html, 'data-testid="tile-proteomics"')
  expect_match(html, 'id="proteomics_box"')
})

# --- test_mode_dir_input ---

test_that("test_mode_dir_input() returns a textInput with data-testid", {
  ns <- shiny::NS("mymodule")
  result <- test_mode_dir_input(ns, "project_dir", "Project Directory", "/tmp")
  html <- as.character(result)
  expect_match(html, 'data-testid="project_dir"')
  expect_match(html, 'id="mymodule-project_dir"')
  expect_match(html, 'value="/tmp"')
})

# --- test_mode_digest_ui ---

test_that("test_mode_digest_ui() returns NULL when test mode is off", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , expect_null(test_mode_digest_ui())
    )
  )
})

test_that("test_mode_digest_ui() returns a hidden div when test mode is on", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      result <- test_mode_digest_ui(ns = shiny::NS("mod"))
      expect_false(is.null(result))
    }
  )
})

test_that("test_mode_digest_ui() uses output id 'test_state_digest' (no double underscores)", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      ns <- shiny::NS("mod")
      result <- test_mode_digest_ui(ns = ns)
      html <- as.character(result)
      expect_match(html, "test_state_digest")
      expect_false(grepl("__test_state_digest__", html))
    }
  )
})

test_that("test_mode_digest_ui(NULL) works via identity NS default", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      result <- test_mode_digest_ui(ns = NULL)
      expect_false(is.null(result))
      html <- as.character(result)
      expect_match(html, "test_state_digest")
    }
  )
})

# --- collect_state_digest ---

test_that("collect_state_digest() returns complete list with empty inputs", {
  result <- collect_state_digest()
  expected_keys <- c(
    "selected_omics", "initialized_omics", "project_dir_keys"
    , "experiment_label", "workflow_type_per_omic", "step_status_per_omic"
    , "active_tab_per_omic", "export_paths", "report_fingerprints"
  )
  expect_true(all(expected_keys %in% names(result)))
})

test_that("collect_state_digest() returns correct selected_omics and initialized_omics", {
  vals <- list(
    selected_omics = c("proteomics", "metabolomics")
    , initialized_omics = c("proteomics")
  )
  result <- collect_state_digest(values = vals)
  expect_identical(result$selected_omics, c("proteomics", "metabolomics"))
  expect_identical(result$initialized_omics, c("proteomics"))
})

test_that("collect_state_digest() extracts non-NULL project_dir_keys", {
  vals <- list(
    project_dirs = list(proteomics = "/tmp/prot", metabolomics = NULL, lipidomics = "/tmp/lipid")
  )
  result <- collect_state_digest(values = vals)
  expect_setequal(result$project_dir_keys, c("proteomics", "lipidomics"))
})

test_that("collect_state_digest() handles partial workflow_states", {
  ws <- list(
    proteomics = list(
      workflow_type = "LFQ"
      , current_state = "step2"
      , states = list(step1 = list(), step2 = list())
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_identical(result$workflow_type_per_omic[["proteomics"]], "LFQ")
  expect_identical(result$active_tab_per_omic[["proteomics"]], "step2")
  expect_setequal(result$step_status_per_omic[["proteomics"]], c("step1", "step2"))
})

test_that("collect_state_digest() returns NULL (not NA) for missing report files", {
  tmp <- tempfile(fileext = ".html")
  # tmp does not exist
  vals <- list(report_files = c(missing = tmp))
  result <- collect_state_digest(values = vals)
  expect_null(result$report_fingerprints[[tmp]])
})

test_that("collect_state_digest() returns md5 hash for existing report files", {
  tmp <- tempfile(fileext = ".html")
  writeLines("<html><body>report</body></html>", tmp)
  on.exit(unlink(tmp), add = TRUE)

  vals <- list(report_files = c(report = tmp))
  result <- collect_state_digest(values = vals)
  # Should be a non-null, non-NA character hash
  fp <- result$report_fingerprints[[tmp]]
  expect_false(is.null(fp))
  expect_false(is.na(fp))
  expect_match(fp, "^[a-f0-9]{32}$")
})

test_that("collect_state_digest() calls active_tab function (reactive path)", {
  # Simulate a module server that exposes active_tab as a function (reactive)
  active_tab_fn <- function() "normalization"
  ws <- list(
    proteomics = list(
      active_tab = active_tab_fn
      , tab_status = list(setup_import = "complete", normalization = "pending")
      , processing_log = list()
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_identical(result$active_tab_per_omic[["proteomics"]], "normalization")
})

test_that("collect_state_digest() reads workflow_type from state_manager field", {
  sm <- WorkflowState$new()
  sm$setWorkflowType("TMT")
  ws <- list(
    proteomics = list(
      state_manager = sm
      , tab_status = list(setup_import = "complete")
      , processing_log = list()
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_identical(result$workflow_type_per_omic[["proteomics"]], "TMT")
})

test_that("collect_state_digest() reads tab_status as step_status_per_omic", {
  expected_status <- list(
    setup_import = "complete"
    , design_matrix = "pending"
    , quality_control = "disabled"
  )
  ws <- list(
    metabolomics = list(
      tab_status = expected_status
      , processing_log = list()
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_identical(result$step_status_per_omic[["metabolomics"]], expected_status)
})

test_that("collect_state_digest() returns export_paths from processing_log", {
  log_entry <- list(step = "normalization", path = "/tmp/out.rds")
  ws <- list(
    lipidomics = list(
      processing_log = list(log_entry)
      , tab_status = list()
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_identical(result$export_paths[["lipidomics"]], list(log_entry))
})

test_that("collect_state_digest() returns NULL for active_tab when no active_tab key exists", {
  ws <- list(
    proteomics = list(
      tab_status = list(setup_import = "complete")
      , processing_log = list()
    )
  )
  result <- collect_state_digest(workflow_states = ws)
  expect_null(result$active_tab_per_omic[["proteomics"]])
})

test_that("collect_state_digest() handles multiple omics with partial initialization", {
  # Only proteomics initialized; metabolomics not in workflow_states
  vals <- list(
    selected_omics = c("proteomics", "metabolomics")
    , initialized_omics = "proteomics"
  )
  ws <- list(
    proteomics = list(
      active_tab = function() "design"
      , tab_status = list(setup_import = "complete", design_matrix = "pending")
      , processing_log = list()
    )
  )
  result <- collect_state_digest(values = vals, workflow_states = ws)
  expect_identical(result$selected_omics, c("proteomics", "metabolomics"))
  expect_identical(result$initialized_omics, "proteomics")
  expect_identical(result$active_tab_per_omic[["proteomics"]], "design")
  expect_null(result$active_tab_per_omic[["metabolomics"]])
})

test_that("collect_state_digest() returns empty lists when workflow_states is empty", {
  result <- collect_state_digest(workflow_states = list())
  expect_identical(result$workflow_type_per_omic, list())
  expect_identical(result$step_status_per_omic, list())
  expect_identical(result$active_tab_per_omic, list())
  expect_identical(result$export_paths, list())
})
