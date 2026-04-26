# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

makeMetabImportFixtureData <- function(label = "import") {
  data.frame(
    "Peak ID" = c(paste0(label, "_m1"), paste0(label, "_m2")),
    "Sample A" = c(10, 11),
    "Sample B" = c(20, 21),
    Annotation = c("A", "B"),
    check.names = FALSE
  )
}

makeMetabImportHarness <- function() {
  capture <- new.env(parent = emptyenv())
  capture$notifications <- list()
  capture$removed_notifications <- character()
  capture$update_select_calls <- list()
  capture$update_text_calls <- list()
  capture$workflow_types <- character()

  root_dir <- tempfile("metab-import-characterization-")
  dir.create(root_dir, recursive = TRUE)
  writeLines(
    c(
      "Peak ID,Sample A,Sample B,Annotation",
      "m1,10,20,A",
      "m2,11,21,B"
    ),
    file.path(root_dir, "assay1.csv")
  )
  writeLines(
    c(
      "Peak ID,Sample A,Sample B,Annotation",
      "n1,5,7,A",
      "n2,6,8,B"
    ),
    file.path(root_dir, "assay2.csv")
  )

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_tbl <- NULL
  workflow_data$data_format <- NULL
  workflow_data$data_type <- NULL
  workflow_data$column_mapping <- NULL
  workflow_data$processing_log <- list(setup_import = NULL)
  workflow_data$tab_status <- list(setup_import = "incomplete", downstream = "pending")
  workflow_data$state_manager <- list(
    setWorkflowType = function(type) {
      capture$workflow_types <<- c(capture$workflow_types, type)
      invisible(type)
    }
  )

  list(
    capture = capture,
    workflow_data = workflow_data,
    experiment_paths = list(source_dir = tempdir()),
    volumes = c(home = root_dir),
    session = shiny::MockShinySession$new()
  )
}

launchMetabImportModule <- function(harness) {
  shiny::withReactiveDomain(harness$session, {
    mod_metab_import_server(
      id = "import",
      workflow_data = harness$workflow_data,
      experiment_paths = harness$experiment_paths,
      volumes = harness$volumes
    )
  })

  harness$session$getReturned()
}

makeShinyFilesSelection <- function(filename, root = "home") {
  list(root = root, files = list(c(filename)))
}

setMetabImportInputs <- function(session, ...) {
  values <- list(...)
  names(values) <- paste0("import-", names(values))
  do.call(session$setInputs, values)
  shiny:::flushReact()
}

installMetabImportCommonMocks <- function(harness, importer_fn) {
  local_mocked_bindings(
    detectMetabolomicsFormat = function(...) {
      list(format = "msdial", confidence = 0.95)
    },
    importMetabMSDIALData = importer_fn,
    validateMetabColumnMapping = function(...) {
      list(
        valid = TRUE,
        summary = list(n_metabolites = 2L, n_samples = 2L, pct_missing = 0),
        warnings = character(),
        errors = character()
      )
    },
    .env = asNamespace("MultiScholaR")
  )

  local_mocked_bindings(
    updateSelectInput = function(session, inputId, choices = NULL, selected = NULL) {
      harness$capture$update_select_calls[[length(harness$capture$update_select_calls) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        choices = choices,
        selected = selected
      )
      invisible(NULL)
    },
    updateTextInput = function(session, inputId, value = NULL, label = NULL, placeholder = NULL) {
      harness$capture$update_text_calls[[length(harness$capture$update_text_calls) + 1L]] <<- list(
        session = session,
        inputId = inputId,
        value = value,
        label = label,
        placeholder = placeholder
      )
      invisible(NULL)
    },
    showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL, type = "default", session = getDefaultReactiveDomain()) {
      harness$capture$notifications[[length(harness$capture$notifications) + 1L]] <<- list(
        ui = ui,
        id = id,
        type = type,
        duration = duration
      )
      if (is.null(id)) {
        invisible(paste0("notification-", length(harness$capture$notifications)))
      } else {
        invisible(id)
      }
    },
    removeNotification = function(id, session = getDefaultReactiveDomain()) {
      harness$capture$removed_notifications <<- c(harness$capture$removed_notifications, id)
      invisible(NULL)
    },
    .env = asNamespace("shiny")
  )

  local_mocked_bindings(
    make_clean_names = function(x, ...) c("sample_a", "sample_b"),
    .env = asNamespace("janitor")
  )
}

test_that("mod_metab_import_server preserves assay-selection hydration behavior", {
  harness <- makeMetabImportHarness()
  assay_data <- makeMetabImportFixtureData()
  import_calls <- character()

  installMetabImportCommonMocks(
    harness,
    importer_fn = function(path) {
      import_calls <<- c(import_calls, path)
      expect_identical(path, file.path(unname(harness$volumes["home"]), "assay1.csv"))
      list(
        data = assay_data,
        detected_columns = list(
          metabolite_id = "Peak ID",
          annotation = "Annotation"
        ),
        sample_columns = c("Sample A", "Sample B"),
        is_pattern = "^IS_"
      )
    }
  )

  launchMetabImportModule(harness)

  setMetabImportInputs(
    harness$session,
    vendor_format = "custom",
    assay1_name = "LCMS_Pos",
    assay1_file = makeShinyFilesSelection("assay1.csv")
  )

  expect_identical(import_calls, file.path(unname(harness$volumes["home"]), "assay1.csv"))
})

test_that("mod_metab_import_server preserves process-import success behavior", {
  harness <- makeMetabImportHarness()
  assay_data <- makeMetabImportFixtureData()

  installMetabImportCommonMocks(
    harness,
    importer_fn = function(path) {
      expect_identical(path, file.path(unname(harness$volumes["home"]), "assay1.csv"))
      list(
        data = assay_data,
        detected_columns = list(
          metabolite_id = "Peak ID",
          annotation = "Annotation"
        ),
        sample_columns = c("Sample A", "Sample B"),
        is_pattern = "^IS_"
      )
    }
  )

  launchMetabImportModule(harness)

  setMetabImportInputs(
    harness$session,
    vendor_format = "custom",
    assay1_name = "LCMS_Pos",
    metabolite_id_col_custom = "Peak ID",
    annotation_col_custom = "Annotation",
    sample_cols_pattern = "^Sample",
    sanitize_names = TRUE,
    is_pattern = "^IS_",
    assay1_file = makeShinyFilesSelection("assay1.csv")
  )

  setMetabImportInputs(harness$session, process_import = 1)

  expect_identical(names(harness$workflow_data$data_tbl), "LCMS_Pos")
  expect_identical(
    names(harness$workflow_data$data_tbl$LCMS_Pos),
    c("Peak ID", "sample_a", "sample_b", "Annotation")
  )
  expect_identical(harness$workflow_data$data_format, "custom")
  expect_identical(harness$workflow_data$data_type, "metabolite")
  expect_identical(harness$workflow_data$column_mapping$metabolite_id_col, "Peak ID")
  expect_identical(harness$workflow_data$column_mapping$annotation_col, "Annotation")
  expect_identical(harness$workflow_data$column_mapping$sample_columns, c("sample_a", "sample_b"))
  expect_identical(harness$workflow_data$tab_status$setup_import, "complete")
  expect_identical(harness$capture$workflow_types, "metabolomics_standard")
  expect_identical(harness$workflow_data$processing_log$setup_import$assay_names, "LCMS_Pos")
  expect_identical(harness$workflow_data$processing_log$setup_import$n_assays, 1L)
  expect_identical(harness$workflow_data$processing_log$setup_import$n_samples, 2L)
})

test_that("mod_metab_import_server preserves assay-selection error behavior", {
  harness <- makeMetabImportHarness()
  import_calls <- character()

  installMetabImportCommonMocks(
    harness,
    importer_fn = function(path) {
      import_calls <<- c(import_calls, path)
      stop("assay1 import boom")
    }
  )

  launchMetabImportModule(harness)

  expect_no_error(
    setMetabImportInputs(
      harness$session,
      vendor_format = "custom",
      assay1_name = "LCMS_Pos",
      assay1_file = makeShinyFilesSelection("assay1.csv")
    )
  )

  expect_identical(import_calls, file.path(unname(harness$volumes["home"]), "assay1.csv"))
  expect_null(harness$workflow_data$data_tbl)
  expect_null(harness$workflow_data$processing_log$setup_import)
})

test_that("mod_metab_import_server preserves process-import error behavior", {
  harness <- makeMetabImportHarness()
  assay_data <- makeMetabImportFixtureData()

  installMetabImportCommonMocks(
    harness,
    importer_fn = function(path) {
      if (identical(path, file.path(unname(harness$volumes["home"]), "assay1.csv"))) {
        return(list(
          data = assay_data,
          detected_columns = list(
            metabolite_id = "Peak ID",
            annotation = "Annotation"
          ),
          sample_columns = c("Sample A", "Sample B"),
          is_pattern = "^IS_"
        ))
      }

      stop("assay2 boom")
    }
  )

  launchMetabImportModule(harness)

  setMetabImportInputs(
    harness$session,
    vendor_format = "custom",
    assay1_name = "LCMS_Pos",
    assay2_name = "LCMS_Neg",
    metabolite_id_col_custom = "Peak ID",
    annotation_col_custom = "Annotation",
    sample_cols_pattern = "^Sample",
    sanitize_names = FALSE,
    assay1_file = makeShinyFilesSelection("assay1.csv")
  )
  setMetabImportInputs(harness$session, assay2_file = makeShinyFilesSelection("assay2.csv"))
  setMetabImportInputs(harness$session, process_import = 1)

  expect_null(harness$workflow_data$data_tbl)
  expect_null(harness$workflow_data$column_mapping)
  expect_identical(harness$workflow_data$tab_status$setup_import, "incomplete")
  expect_identical(harness$capture$workflow_types, character())
  expect_null(harness$workflow_data$processing_log$setup_import)
})
