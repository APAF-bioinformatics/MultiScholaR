library(testthat)

# Helper: render a UI function to HTML string
ui_html <- function(fn, ...) {
  as.character(fn(...))
}

# ── app_ui ─────────────────────────────────────────────────────────────────────

test_that("app_ui in test mode contains tile-proteomics data-testid", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(app_ui, NULL)
      expect_match(html, 'data-testid="tile-proteomics"', fixed = TRUE)
    }
  )
})

test_that("app_ui in test mode exposes test_state_digest output", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(app_ui, NULL)
      expect_true(grepl("test_state_digest", html, fixed = TRUE))
    }
  )
})

test_that("app_ui in production mode does not expose test_state_digest", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        html <- ui_html(app_ui, NULL)
        expect_false(grepl("test_state_digest", html, fixed = TRUE))
      }
    )
  )
})

# ── mod_prot_import_ui: test mode ──────────────────────────────────────────────

test_that("mod_prot_import_ui in test mode: fileInput elements present", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(mod_prot_import_ui, "test")
      # fileInput renders <input type="file"> elements
      expect_true(grepl('type="file"', html, fixed = TRUE))
      expect_match(html, 'data-testid="prot-import-search-results-file"', fixed = TRUE)
      expect_match(html, 'data-testid="prot-import-fasta-file"', fixed = TRUE)
    }
  )
})

test_that("mod_prot_import_ui in test mode: shinyFilesButton absent", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(mod_prot_import_ui, "test")
      # When test mode is active, shinyFiles is disabled; standard IDs appear
      expect_true(grepl("search_results_standard", html, fixed = TRUE))
      # shinyFilesButton-only ID should not appear (it's replaced by fileInput)
      expect_false(grepl('"test-search_results"', html, fixed = TRUE))
    }
  )
})

test_that("mod_prot_import_ui in test mode: prot-import-process data-testid present", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(mod_prot_import_ui, "test")
      expect_match(html, 'data-testid="prot-import-process"', fixed = TRUE)
    }
  )
})

# ── metabolomics/lipidomics import UI: test mode upload contract ──────────────

test_that("mod_metab_import_ui in test mode exposes browser-upload controls", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(mod_metab_import_ui, "test")
      expect_match(html, "assay1_file_std", fixed = TRUE)
      expect_match(html, "assay2_file_std", fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-assay1-file"', fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-assay2-file"', fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-process"', fixed = TRUE)
    }
  )
})

test_that("mod_lipid_import_ui in test mode exposes browser-upload controls", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      html <- ui_html(mod_lipid_import_ui, "test")
      expect_match(html, "assay1_file_std", fixed = TRUE)
      expect_match(html, "assay2_file_std", fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-assay1-file"', fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-assay2-file"', fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-process"', fixed = TRUE)
    }
  )
})

test_that("standard metabolomics file input assigns upload datapath and triggers preview", {
  local_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  callback_path <- NULL

  result <- handleMetabImportStandardAssayFileSelection(
    fileInput = data.frame(datapath = "/tmp/metab.tsv"),
    localData = local_data,
    localField = "assay1_file",
    output = output,
    outputId = "assay1_path",
    onSelected = function(path) callback_path <<- path,
    renderTextFn = function(value) value
  )

  expect_identical(result, "/tmp/metab.tsv")
  expect_identical(local_data$assay1_file, "/tmp/metab.tsv")
  expect_identical(output$assay1_path, "/tmp/metab.tsv")
  expect_identical(callback_path, "/tmp/metab.tsv")
})

test_that("standard lipidomics file input assigns upload datapath and triggers preview", {
  selected_path <- NULL
  rendered_path <- NULL
  callback_called <- FALSE

  result <- handleLipidImportStandardFileSelection(
    fileInput = data.frame(datapath = "/tmp/lipid.tsv"),
    assignSelectedPath = function(path) selected_path <<- path,
    setRenderedPath = function(value) rendered_path <<- value,
    onSelected = function() callback_called <<- TRUE,
    assignAssayPath = function(selectedPath, assignSelectedPath, setRenderedPath, onSelected = NULL) {
      assignSelectedPath(selectedPath)
      setRenderedPath(selectedPath)
      onSelected()
      selectedPath
    }
  )

  expect_identical(result, "/tmp/lipid.tsv")
  expect_identical(selected_path, "/tmp/lipid.tsv")
  expect_identical(rendered_path, "/tmp/lipid.tsv")
  expect_true(callback_called)
})

test_that("lipidomics shinyFiles probe is disabled in test mode", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      result <- probeLipidImportShinyFilesAvailability(
        requireNamespaceFn = function(...) TRUE,
        emitMessage = function(...) NULL
      )
      expect_false(result)
    }
  )
})

# ── mod_prot_import_ui: production mode ───────────────────────────────────────

test_that("mod_prot_import_ui in production mode: shinyFilesButton present when shinyFiles available", {
  skip_if_not(
    requireNamespace("shinyFiles", quietly = TRUE)
    , "shinyFiles not installed"
  )
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        html <- ui_html(mod_prot_import_ui, "test")
        # shinyFilesButton uses the un-suffixed ID; fileInput uses _standard suffix
        expect_false(grepl("search_results_standard", html, fixed = TRUE))
      }
    )
  )
})

test_that("mod_prot_import_ui production mode still renders prot-import-process testid", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        html <- ui_html(mod_prot_import_ui, "test")
        expect_match(html, 'data-testid="prot-import-process"', fixed = TRUE)
      }
    )
  )
})

# ── test_state_digest ID check ─────────────────────────────────────────────────

test_that("test_state_digest present in test mode, absent in production (round-trip check)", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , expect_true(grepl("test_state_digest", ui_html(app_ui, NULL), fixed = TRUE))
  )
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , expect_false(grepl("test_state_digest", ui_html(app_ui, NULL), fixed = TRUE))
    )
  )
})

# ── mod_prot_norm_ui ───────────────────────────────────────────────────────────

test_that("mod_prot_norm_ui contains prot-norm-run data-testid", {
  html <- ui_html(mod_prot_norm_ui, "test")
  expect_match(html, 'data-testid="prot-norm-run"', fixed = TRUE)
})

test_that("mod_prot_norm_ui contains prot-norm-export-session data-testid", {
  html <- ui_html(mod_prot_norm_ui, "test")
  expect_match(html, 'data-testid="prot-norm-export-session"', fixed = TRUE)
})

# ── mod_prot_da_ui ────────────────────────────────────────────────────────────

test_that("mod_prot_da_ui contains prot-da-run data-testid", {
  html <- ui_html(mod_prot_da_ui, "test")
  expect_match(html, 'data-testid="prot-da-run"', fixed = TRUE)
})

test_that("mod_prot_da_ui contains prot-da-load-session data-testid", {
  html <- ui_html(mod_prot_da_ui, "test")
  expect_match(html, 'data-testid="prot-da-load-session"', fixed = TRUE)
})

# ── mod_prot_summary_ui ───────────────────────────────────────────────────────

test_that("mod_prot_summary_ui contains prot-summary-generate-report data-testid", {
  html <- ui_html(mod_prot_summary_ui, "test")
  expect_match(html, 'data-testid="prot-summary-generate-report"', fixed = TRUE)
})

test_that("all summary UIs expose publication-copy and save-args test IDs", {
  expect_match(ui_html(mod_prot_summary_ui, "test"), 'data-testid="prot-summary-save-args"', fixed = TRUE)
  expect_match(ui_html(mod_prot_summary_ui, "test"), 'data-testid="prot-summary-copy-publication"', fixed = TRUE)

  expect_match(ui_html(mod_metab_summary_ui, "test"), 'data-testid="metab-summary-save-args"', fixed = TRUE)
  expect_match(ui_html(mod_metab_summary_ui, "test"), 'data-testid="metab-summary-copy-publication"', fixed = TRUE)

  expect_match(ui_html(mod_lipid_summary_ui, "test"), 'data-testid="lipid-summary-save-args"', fixed = TRUE)
  expect_match(ui_html(mod_lipid_summary_ui, "test"), 'data-testid="lipid-summary-copy-publication"', fixed = TRUE)
})

test_that("mod_prot_summary_ui contains prot-summary-export-state data-testid", {
  html <- ui_html(mod_prot_summary_ui, "test")
  expect_match(html, 'data-testid="prot-summary-export-state"', fixed = TRUE)
})

# ── No regression: all import UIs render without error ────────────────────────

test_that("all import UIs render without error in production mode", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        expect_no_error(mod_prot_import_ui("test"))
        expect_no_error(mod_metab_import_ui("test"))
        expect_no_error(mod_lipid_import_ui("test"))
      }
    )
  )
})

test_that("all import UIs render without error in test mode", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      expect_no_error(mod_prot_import_ui("test"))
      expect_no_error(mod_metab_import_ui("test"))
      expect_no_error(mod_lipid_import_ui("test"))
    }
  )
})

# ── No regression: all orchestrator UIs render without error ──────────────────

test_that("orchestrator UIs render without error in production mode", {
  withr::with_options(
    list(multischolar.test_mode = FALSE)
    , withr::with_envvar(
      c(MULTISCHOLAR_TEST_MODE = "")
      , {
        expect_no_error(mod_prot_norm_ui("test"))
        expect_no_error(mod_prot_da_ui("test"))
        expect_no_error(mod_prot_summary_ui("test"))
      }
    )
  )
})

test_that("orchestrator UIs render without error in test mode", {
  withr::with_options(
    list(multischolar.test_mode = TRUE)
    , {
      expect_no_error(mod_prot_norm_ui("test"))
      expect_no_error(mod_prot_da_ui("test"))
      expect_no_error(mod_prot_summary_ui("test"))
    }
  )
})
