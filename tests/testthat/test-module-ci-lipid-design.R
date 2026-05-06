library(testthat)

test_that("MCI-019.1 builder design matrix covers accepted and invalid lipidomics design shapes", {
  accepted <- c("two_group", "three_group", "unbalanced", "batch_aware", "extra_metadata")
  for (kind in accepted) {
    payload <- module_ci_lipid_design_payload(
      kind = kind,
      layout = if (identical(kind, "three_group")) "combined" else "lc_pair",
      contrasts = if (identical(kind, "three_group")) {
        data.frame(contrasts = c("groupKO-groupWT", "groupRES-groupWT"), stringsAsFactors = FALSE)
      } else if (identical(kind, "batch_aware")) {
        data.frame(contrasts = "groupKO_B1-groupWT_B1", stringsAsFactors = FALSE)
      } else {
        module_ci_lipid_design_contrasts("friendly")
      }
    )

    preflight <- module_ci_lipid_design_assert_preflight_passes(payload)
    module_ci_lipid_design_assert_samples_once(payload$design_matrix, payload$column_mapping$sample_columns)
    expect_false(any(is.na(payload$design_matrix$group) | payload$design_matrix$group == ""), info = kind)
    expect_false(any(is.na(payload$design_matrix$replicates)), info = kind)
    expect_true(length(preflight$model_terms) >= 2L, info = kind)

    if (identical(kind, "extra_metadata")) {
      expect_true(all(c("site", "instrument_batch") %in% names(payload$design_matrix)))
    }
  }

  invalid_payload <- module_ci_lipid_design_payload(kind = "missing_group")
  invalid <- validateLipidDesignDaPreflight(
    designMatrix = invalid_payload$design_matrix,
    assayList = invalid_payload$data_cln,
    contrastsTbl = invalid_payload$contrasts_tbl,
    formulaString = "~ 0 + group",
    columnMapping = invalid_payload$column_mapping
  )
  expect_false(invalid$valid)
  expect_true(any(grepl("blank assignments", invalid$errors, fixed = TRUE)))
})

test_that("MCI-019.2 multi-assay alignment accepts valid lipid assays and rejects corrupt sample sets", {
  design <- module_ci_lipid_design_matrix("two_group")
  mapping <- module_ci_lipid_design_column_mapping(design$Run)

  valid_cases <- list(
    lc_pair = module_ci_lipid_design_assays("lc_pair", samples = design$Run),
    gc_only = module_ci_lipid_design_assays("gc", samples = design$Run),
    combined_reordered = module_ci_lipid_design_assays("combined", samples = design$Run, reorder_secondary = TRUE)
  )
  for (case_name in names(valid_cases)) {
    alignment <- validateLipidDesignAlignment(
      designMatrix = design,
      assayList = valid_cases[[case_name]],
      columnMapping = mapping
    )
    expect_true(alignment$valid, info = paste(case_name, paste(alignment$errors, collapse = "; ")))
    expect_setequal(
      unlist(lapply(alignment$assay_summaries, `[[`, "sample_columns"), use.names = FALSE),
      design$Run
    )
  }

  invalid_cases <- list(
    missing_sample = module_ci_lipid_design_assays("lc_pair", samples = design$Run, missing_sample = "WT_2"),
    extra_sample = module_ci_lipid_design_assays("lc_pair", samples = design$Run, extra_sample = "EXTRA_1"),
    case_varied = module_ci_lipid_design_assays("lc_pair", samples = design$Run, case_varied = TRUE)
  )
  for (case_name in names(invalid_cases)) {
    alignment <- validateLipidDesignAlignment(
      designMatrix = design,
      assayList = invalid_cases[[case_name]],
      columnMapping = mapping
    )
    expect_false(alignment$valid, info = case_name)
  }
  expect_true(any(grepl("missing mapped sample", validateLipidDesignAlignment(design, invalid_cases$missing_sample, mapping)$errors)))
  expect_true(any(grepl("sample-like column", validateLipidDesignAlignment(design, invalid_cases$extra_sample, mapping)$errors)))
  expect_true(any(grepl("case-varied sample names", validateLipidDesignAlignment(design, invalid_cases$case_varied, mapping)$errors)))
})

test_that("MCI-019.3 contrast matrix preserves lipid friendly/raw serialization and rejects unsafe contrasts", {
  design <- module_ci_lipid_design_matrix("two_group")

  raw <- validateLipidDesignContrasts(design, module_ci_lipid_design_contrasts("raw_terms"), "~ 0 + group")
  expect_true(raw$valid)
  expect_identical(raw$contrasts_tbl$contrasts, "groupKO-groupWT")
  expect_identical(raw$contrasts_tbl$friendly_names, "KO_vs_WT")
  expect_identical(raw$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")

  friendly <- validateLipidDesignContrasts(design, module_ci_lipid_design_contrasts("friendly"), "~ 0 + group")
  expect_true(friendly$valid)
  expect_identical(friendly$contrasts_tbl$friendly_names, "KO_vs_WT")
  expect_identical(friendly$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")

  reversed <- validateLipidDesignContrasts(design, module_ci_lipid_design_contrasts("reversed"), "~ 0 + group")
  expect_true(reversed$valid)
  expect_identical(reversed$contrasts_tbl$contrasts, "groupWT-groupKO")

  duplicate <- validateLipidDesignContrasts(design, module_ci_lipid_design_contrasts("duplicate"), "~ 0 + group")
  expect_false(duplicate$valid)
  expect_true(any(grepl("duplicate contrast definitions", duplicate$errors, fixed = TRUE)))

  invalid <- validateLipidDesignContrasts(design, module_ci_lipid_design_contrasts("invalid_term"), "~ 0 + group")
  expect_false(invalid$valid)
  expect_true(any(grepl("terms absent from the model matrix", invalid$errors, fixed = TRUE)))

  no_contrasts <- validateLipidDesignContrasts(design, NULL, "~ 0 + group")
  expect_false(no_contrasts$valid)
  expect_true(any(grepl("non-empty data frame", no_contrasts$errors, fixed = TRUE)))
})

test_that("MCI-019.4 existing-design import and browser controls preserve lipid builder equivalence", {
  payload <- module_ci_lipid_design_payload(layout = "combined")
  import_root <- module_ci_lipid_design_write_import_pack(tempfile("module-ci-lipid-design-import-"), payload)
  workflow <- module_ci_lipid_design_workflow(payload)
  workflow$state_manager <- list(saveState = function(...) invisible(NULL))
  qc_value <- FALSE
  assigned_names <- character()
  s4_token <- structure(list(kind = "lipid-s4"), class = "mock_lipid_s4")

  runLipidDesignImportConfirmationShell(
    workflowData = workflow,
    experimentPaths = list(source_dir = import_root),
    importPath = import_root,
    qcTrigger = function(value) qc_value <<- value,
    removeModalFn = function() invisible(NULL),
    removeNotificationFn = function(...) invisible(NULL),
    showNotificationFn = function(...) invisible(NULL),
    readConfigFileFn = function(file) module_ci_lipid_design_config(),
    vroomFn = function(path, show_col_types = FALSE) {
      utils::read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    },
    readJsonFn = function(path, simplifyVector = TRUE) {
      jsonlite::read_json(path, simplifyVector = simplifyVector)
    },
    assignFn = function(name, value, envir) assigned_names <<- c(assigned_names, name),
    createLipidomicsAssayDataFn = function(...) s4_token,
    workflowStateClass = list(new = function(type) stop("existing state manager should be reused")),
    updateLipidFilteringFn = function(...) invisible(NULL)
  )

  expect_identical(workflow$design_matrix$Run, payload$design_matrix$Run)
  expect_identical(workflow$design_matrix$group, payload$design_matrix$group)
  expect_identical(names(workflow$data_cln), names(payload$data_cln))
  expect_identical(workflow$column_mapping$sample_columns, payload$column_mapping$sample_columns)
  expect_identical(workflow$contrasts_tbl$contrasts, payload$contrasts_tbl$contrasts)
  expect_identical(workflow$contrasts_tbl$friendly_names, payload$contrasts_tbl$friendly_names)
  expect_true(qc_value)
  expect_true(all(c("config_list", "contrasts_tbl") %in% assigned_names))

  captured <- new.env(parent = emptyenv())
  modal <- buildLipidDesignImportModal(
    ns = function(id) paste0("lipid-design-", id),
    modalDialogFn = function(..., title = NULL, footer = NULL, easyClose = TRUE) list(title = title, body = list(...), footer = footer, easyClose = easyClose),
    paragraphFn = function(...) paste(...),
    helpTextFn = function(...) paste(...),
    dirButtonFn = function(id, label, title) {
      captured$dir_button <- list(id = id, label = label, title = title)
      list(kind = "dirButton", id = id, label = label, title = title)
    },
    verbatimTextOutputFn = function(id, placeholder = FALSE) list(kind = "verbatim", id = id, placeholder = placeholder),
    tagListFn = function(...) list(...),
    modalButtonFn = function(label) list(kind = "modalButton", label = label),
    actionButtonFn = function(inputId, label, class = NULL) {
      captured$action_button <- list(inputId = inputId, label = label, class = class)
      list(kind = "actionButton", inputId = inputId, label = label, class = class)
    }
  )
  expect_identical(modal$title, "Import Existing Design Matrix")
  expect_identical(captured$dir_button$id, "lipid-design-import_dir")
  expect_identical(captured$action_button$inputId, "lipid-design-confirm_import")
})

test_that("MCI-019.5 design artifacts preserve expected lipid filenames and serialized state", {
  payload <- module_ci_lipid_design_payload(layout = "lc_pair")
  paths <- module_ci_lipid_design_paths()
  workflow <- module_ci_lipid_design_workflow(payload)
  saved_states <- list()
  workflow$state_manager <- list(saveState = function(...) saved_states[[length(saved_states) + 1L]] <<- list(...))
  s4_token <- structure(list(kind = "lipid-s4"), class = "mock_lipid_s4")
  qc_value <- FALSE

  local_mocked_bindings(
    createLipidomicsAssayData = function(...) s4_token,
    .env = asNamespace("MultiScholaR")
  )

  shiny::withReactiveDomain(shiny::MockShinySession$new(), {
    runLipidDesignBuilderObserverShell(
      results = payload,
      workflowData = workflow,
      experimentPaths = paths,
      qcTrigger = function(value) qc_value <<- value
    )
  })

  expected_files <- c(
    "design_matrix.tab",
    "contrast_strings.tab",
    "assay_manifest.txt",
    "column_mapping.json",
    "manifest.json",
    "config.ini",
    paste0("data_cln_", names(payload$data_cln), ".tab")
  )
  for (file_name in expected_files) {
    module_ci_lipid_design_file_nonempty(file.path(paths$source_dir, file_name))
  }

  design <- utils::read.delim(file.path(paths$source_dir, "design_matrix.tab"), sep = "\t", check.names = FALSE)
  contrasts <- readLines(file.path(paths$source_dir, "contrast_strings.tab"), warn = FALSE)
  manifest <- readLines(file.path(paths$source_dir, "assay_manifest.txt"), warn = FALSE)
  column_mapping <- jsonlite::read_json(file.path(paths$source_dir, "column_mapping.json"), simplifyVector = TRUE)

  expect_identical(design$Run, payload$design_matrix$Run)
  expect_identical(contrasts, payload$contrasts_tbl$contrasts)
  expect_identical(manifest, names(payload$data_cln))
  expect_identical(column_mapping$sample_columns, payload$column_mapping$sample_columns)
  expect_identical(saved_states[[1L]]$state_name, "lipid_raw_data_s4")
  expect_true(qc_value)
  expect_identical(workflow$tab_status$design_matrix, "complete")
})

test_that("MCI-019.6 DA preflight catches lipid formula, contrast, and assay/design drift before state advances", {
  payload <- module_ci_lipid_design_payload(layout = "combined")
  valid <- module_ci_lipid_design_assert_preflight_passes(payload)
  expect_true(all(c("groupWT", "groupKO") %in% valid$model_terms))

  invalid_formula <- validateLipidDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = payload$contrasts_tbl,
    formulaString = "~ 0 + missing_factor",
    columnMapping = payload$column_mapping
  )
  expect_false(invalid_formula$valid)
  expect_true(any(grepl("cannot produce a model matrix", invalid_formula$errors, fixed = TRUE)))

  invalid_contrast <- validateLipidDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = module_ci_lipid_design_contrasts("invalid_term"),
    formulaString = "~ 0 + group",
    columnMapping = payload$column_mapping
  )
  expect_false(invalid_contrast$valid)
  expect_true(any(grepl("terms absent from the model matrix", invalid_contrast$errors, fixed = TRUE)))

  drifted_payload <- payload
  drifted_payload$data_cln <- module_ci_lipid_design_assays("combined", samples = payload$design_matrix$Run, missing_sample = "KO_2")
  drifted <- validateLipidDesignDaPreflight(
    designMatrix = drifted_payload$design_matrix,
    assayList = drifted_payload$data_cln,
    contrastsTbl = drifted_payload$contrasts_tbl,
    formulaString = "~ 0 + group",
    columnMapping = drifted_payload$column_mapping
  )
  expect_false(drifted$valid)
  expect_true(any(grepl("missing mapped sample", drifted$errors)))

  no_contrast_design_save <- validateLipidDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = NULL,
    formulaString = "~ 0 + group",
    columnMapping = payload$column_mapping,
    requireContrasts = FALSE
  )
  expect_true(no_contrast_design_save$valid, info = paste(no_contrast_design_save$errors, collapse = "; "))

  invalid_workflow <- module_ci_lipid_design_workflow(payload)
  invalid_workflow$tab_status <- list(design_matrix = "pending")
  invalid_workflow$state_manager <- list(saveState = function(...) stop("invalid design should not be saved"))
  paths <- module_ci_lipid_design_paths()
  invalid_results <- payload
  invalid_results$data_cln <- drifted_payload$data_cln
  local_mocked_bindings(
    createLipidomicsAssayData = function(...) stop("invalid design should not create S4 state"),
    .env = asNamespace("MultiScholaR")
  )

  shiny::withReactiveDomain(shiny::MockShinySession$new(), {
    expect_null(runLipidDesignBuilderObserverShell(
      results = invalid_results,
      workflowData = invalid_workflow,
      experimentPaths = paths
    ))
  })
  expect_null(invalid_workflow$design_matrix)
  expect_identical(invalid_workflow$tab_status$design_matrix, "pending")
  expect_false(file.exists(file.path(paths$source_dir, "design_matrix.tab")))
})
