test_that("MCI-013.1 builder design matrix covers accepted and invalid design shapes", {
  accepted <- c("two_group", "three_group", "unbalanced", "batch_aware", "extra_metadata")
  for (kind in accepted) {
    payload <- module_ci_metab_design_payload(
      kind = kind,
      layout = if (identical(kind, "three_group")) "combined" else "lc_pair",
      contrasts = if (identical(kind, "three_group")) {
        data.frame(contrasts = c("groupKO-groupWT", "groupRES-groupWT"), stringsAsFactors = FALSE)
      } else if (identical(kind, "batch_aware")) {
        data.frame(contrasts = "groupKO_B1-groupWT_B1", stringsAsFactors = FALSE)
      } else {
        module_ci_metab_design_contrasts("friendly")
      }
    )

    preflight <- module_ci_metab_design_assert_preflight_passes(payload)
    module_ci_metab_design_assert_samples_once(payload$design_matrix, payload$column_mapping$sample_columns)
    expect_false(any(is.na(payload$design_matrix$group) | payload$design_matrix$group == ""), info = kind)
    expect_false(any(is.na(payload$design_matrix$replicates)), info = kind)
    expect_true(length(preflight$model_terms) >= 2L, info = kind)

    if (identical(kind, "extra_metadata")) {
      expect_true(all(c("site", "age_days", "instrument_batch") %in% names(payload$design_matrix)))
    }
  }

  invalid_payload <- module_ci_metab_design_payload(kind = "missing_group")
  invalid <- validateMetabDesignDaPreflight(
    designMatrix = invalid_payload$design_matrix,
    assayList = invalid_payload$data_cln,
    contrastsTbl = invalid_payload$contrasts_tbl,
    formulaString = "~ 0 + group",
    columnMapping = invalid_payload$column_mapping
  )
  expect_false(invalid$valid)
  expect_true(any(grepl("blank assignments", invalid$errors, fixed = TRUE)))
})

test_that("MCI-013.2 multi-assay alignment accepts identical/reordered assays and rejects corrupt sample sets", {
  design <- module_ci_metab_design_matrix("two_group")
  mapping <- module_ci_metab_design_column_mapping(design$Run)

  valid_cases <- list(
    lc_pair = module_ci_metab_design_assays("lc_pair", samples = design$Run),
    gc_only = module_ci_metab_design_assays("gc", samples = design$Run),
    combined_reordered = module_ci_metab_design_assays("combined", samples = design$Run, reorder_secondary = TRUE)
  )
  for (case_name in names(valid_cases)) {
    alignment <- validateMetabDesignAlignment(
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
    missing_sample = module_ci_metab_design_assays("lc_pair", samples = design$Run, missing_sample = "WT_2"),
    extra_sample = module_ci_metab_design_assays("lc_pair", samples = design$Run, extra_sample = "EXTRA_1"),
    case_varied = module_ci_metab_design_assays("lc_pair", samples = design$Run, case_varied = TRUE)
  )
  for (case_name in names(invalid_cases)) {
    alignment <- validateMetabDesignAlignment(
      designMatrix = design,
      assayList = invalid_cases[[case_name]],
      columnMapping = mapping
    )
    expect_false(alignment$valid, info = case_name)
  }
  expect_true(any(grepl("missing mapped sample", validateMetabDesignAlignment(design, invalid_cases$missing_sample, mapping)$errors)))
  expect_true(any(grepl("sample-like column", validateMetabDesignAlignment(design, invalid_cases$extra_sample, mapping)$errors)))
  expect_true(any(grepl("case-varied sample names", validateMetabDesignAlignment(design, invalid_cases$case_varied, mapping)$errors)))
})

test_that("MCI-013.3 contrast matrix preserves friendly/raw serialization and rejects unsafe contrasts", {
  design <- module_ci_metab_design_matrix("two_group")

  raw <- validateMetabDesignContrasts(design, module_ci_metab_design_contrasts("raw_terms"), "~ 0 + group")
  expect_true(raw$valid)
  expect_identical(raw$contrasts_tbl$contrasts, "groupKO-groupWT")
  expect_identical(raw$contrasts_tbl$friendly_names, "KO_vs_WT")
  expect_identical(raw$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")

  friendly <- validateMetabDesignContrasts(design, module_ci_metab_design_contrasts("friendly"), "~ 0 + group")
  expect_true(friendly$valid)
  expect_identical(friendly$contrasts_tbl$friendly_names, "KO_vs_WT")
  expect_identical(friendly$contrasts_tbl$full_format, "KO_vs_WT=groupKO-groupWT")

  reversed <- validateMetabDesignContrasts(design, module_ci_metab_design_contrasts("reversed"), "~ 0 + group")
  expect_true(reversed$valid)
  expect_identical(reversed$contrasts_tbl$contrasts, "groupWT-groupKO")

  duplicate <- validateMetabDesignContrasts(design, module_ci_metab_design_contrasts("duplicate"), "~ 0 + group")
  expect_false(duplicate$valid)
  expect_true(any(grepl("duplicate contrast definitions", duplicate$errors, fixed = TRUE)))

  invalid <- validateMetabDesignContrasts(design, module_ci_metab_design_contrasts("invalid_term"), "~ 0 + group")
  expect_false(invalid$valid)
  expect_true(any(grepl("terms absent from the model matrix", invalid$errors, fixed = TRUE)))

  no_contrasts <- validateMetabDesignContrasts(design, NULL, "~ 0 + group")
  expect_false(no_contrasts$valid)
  expect_true(any(grepl("non-empty data frame", no_contrasts$errors, fixed = TRUE)))
})

test_that("MCI-013.4 existing-design import and browser controls preserve builder equivalence", {
  payload <- module_ci_metab_design_payload(layout = "combined")
  import_root <- module_ci_metab_design_write_import_pack(tempfile("module-ci-metab-design-import-"), payload)
  workflow <- module_ci_metab_design_workflow(payload)

  # KIT-103 W5: the shared design tab's import (runOmicsDesignImport over the core's file check, reader and commit,
  # C10) under the metab configuration; the QC baseline is quiet and the notifications silent.
  config <- omicGuiConfig("metabolomics")
  config$design$design_import$qc_baseline <- function(...) invisible(NULL)
  testthat::local_mocked_bindings(removeModal = function(...) invisible(NULL),
    removeNotification = function(...) invisible(NULL), showNotification = function(...) invisible(NULL),
    .package = "shiny")
  expect_true(checkOmicsDesignImportFiles(import_root)$ok)
  # KIT-075: the commit binds the session's owner before holding the tables, which attaches to the reactive domain.
  shiny::withReactiveDomain(omicsOwnerSession(), runOmicsDesignImport(
    workflowData = workflow,
    experimentPaths = list(source_dir = import_root),
    importPath = import_root,
    qcTrigger = NULL,
    config = config
  ))

  expect_identical(workflow$design_matrix$Run, payload$design_matrix$Run)
  expect_identical(workflow$design_matrix$group, payload$design_matrix$group)
  expect_identical(names(getWorkflowOmicsBusTables(workflow, "data_cln")), names(payload$data_cln))
  expect_identical(workflow$column_mapping$sample_columns, payload$column_mapping$sample_columns)
  expect_identical(workflow$contrasts_tbl$contrasts, payload$contrasts_tbl$contrasts)

  # KIT-103 W5: the shared frame's import modal in test mode; its directory input and confirm button carry the
  # metab test ids.
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(testid = function(tag, id) {
    html <- as.character(tag)
    if (grepl("metab-design-import_dir", html, fixed = TRUE)) captured$dir_button <- list(testid = id)
    if (grepl("metab-design-confirm_import", html, fixed = TRUE)) captured$action_button <- list(testid = id)
    shiny::tagAppendAttributes(tag, "data-testid" = id)
  })
  buildOmicsDesignImportModal(ns = function(id) paste0("metab-design-", id), config = config, testMode = TRUE)
  expect_identical(captured$dir_button$testid, "metab-design-import-dir")
  expect_identical(captured$action_button$testid, "metab-design-import-confirm")
})

test_that("MCI-013.5 design artifacts preserve expected filenames and serialized state", {
  payload <- module_ci_metab_design_payload(layout = "lc_pair")
  paths <- module_ci_metab_design_paths()
  workflow <- module_ci_metab_design_workflow(payload)
  # KIT-103 W5: the artifacts are the shared design tab's commit's (runOmicsDesignCommit over commitOmicsDesign,
  # whose files are writeOmicsDesignArtifacts') under the metab configuration, as KIT-103 W4 did for MCI-019.5.
  shiny::withReactiveDomain(shiny::MockShinySession$new(), {
    runOmicsDesignCommit(
      results = payload,
      workflowData = workflow,
      experimentPaths = paths,
      qcTrigger = NULL,
      config = omicGuiConfig("metabolomics")
    )
  })
  source_dir <- paths$source_dir

  expected_files <- c(
    "design_matrix.tab",
    "contrast_strings.tab",
    "assay_manifest.txt",
    "column_mapping.json",
    paste0("data_cln_", names(payload$data_cln), ".tab")
  )
  for (file_name in expected_files) {
    module_ci_metab_design_file_nonempty(file.path(source_dir, file_name))
  }

  design <- utils::read.delim(file.path(source_dir, "design_matrix.tab"), sep = "\t", check.names = FALSE)
  contrasts <- readLines(file.path(source_dir, "contrast_strings.tab"), warn = FALSE)
  manifest <- readLines(file.path(source_dir, "assay_manifest.txt"), warn = FALSE)
  column_mapping <- jsonlite::read_json(file.path(source_dir, "column_mapping.json"), simplifyVector = TRUE)

  expect_identical(design$Run, payload$design_matrix$Run)
  expect_identical(contrasts, payload$contrasts_tbl$contrasts)
  expect_identical(manifest, names(payload$data_cln))
  expect_identical(column_mapping$sample_columns, payload$column_mapping$sample_columns)

  digest <- list(
    design_matrix = design[c("Run", "group", "replicates")],
    contrasts = contrasts,
    assay_manifest = manifest,
    column_mapping = column_mapping
  )
  expect_identical(digest$design_matrix$Run, payload$design_matrix$Run)
  expect_identical(digest$assay_manifest, c("LCMS_Pos", "LCMS_Neg"))
})

test_that("MCI-013.6 DA preflight catches invalid formulas, contrasts, and assay/design drift", {
  payload <- module_ci_metab_design_payload(layout = "combined")
  valid <- module_ci_metab_design_assert_preflight_passes(payload)
  expect_true(all(c("groupWT", "groupKO") %in% valid$model_terms))

  invalid_formula <- validateMetabDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = payload$contrasts_tbl,
    formulaString = "~ 0 + missing_factor",
    columnMapping = payload$column_mapping
  )
  expect_false(invalid_formula$valid)
  expect_true(any(grepl("cannot produce a model matrix", invalid_formula$errors, fixed = TRUE)))

  invalid_contrast <- validateMetabDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = module_ci_metab_design_contrasts("invalid_term"),
    formulaString = "~ 0 + group",
    columnMapping = payload$column_mapping
  )
  expect_false(invalid_contrast$valid)
  expect_true(any(grepl("terms absent from the model matrix", invalid_contrast$errors, fixed = TRUE)))

  drifted_payload <- payload
  drifted_payload$data_cln <- module_ci_metab_design_assays("combined", samples = payload$design_matrix$Run, missing_sample = "KO_2")
  drifted <- validateMetabDesignDaPreflight(
    designMatrix = drifted_payload$design_matrix,
    assayList = drifted_payload$data_cln,
    contrastsTbl = drifted_payload$contrasts_tbl,
    formulaString = "~ 0 + group",
    columnMapping = drifted_payload$column_mapping
  )
  expect_false(drifted$valid)
  expect_true(any(grepl("missing mapped sample", drifted$errors)))

  no_contrast_strict <- validateMetabDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = NULL,
    formulaString = "~ 0 + group",
    columnMapping = payload$column_mapping,
    requireContrasts = TRUE
  )
  expect_false(no_contrast_strict$valid)

  no_contrast_design_save <- validateMetabDesignDaPreflight(
    designMatrix = payload$design_matrix,
    assayList = payload$data_cln,
    contrastsTbl = NULL,
    formulaString = "~ 0 + group",
    columnMapping = payload$column_mapping,
    requireContrasts = FALSE
  )
  expect_true(no_contrast_design_save$valid, info = paste(no_contrast_design_save$errors, collapse = "; "))
})
