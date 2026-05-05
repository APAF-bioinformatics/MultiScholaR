library(testthat)

test_that("MCI-005.1 builder matrix preserves every imported sample exactly once", {
  for (kind in c("two_group", "three_group", "unbalanced", "multi_factor", "interaction_ready")) {
    design <- module_ci_prot_design_matrix(kind)
    groups <- unique(design$group)
    contrast <- module_ci_prot_contrast_data(c(groups[[min(2L, length(groups))]], groups[[1L]]))
    payload <- buildProtDesignSaveResultsPayload(
      designMatrix = design,
      currentRemovedSamples = character(),
      dataCln = module_ci_prot_design_data(design$Run),
      contrastData = contrast,
      configList = list(deAnalysisParameters = list(formula_string = "~ old")),
      formulaString = "~ 0 + group"
    )

    expect_false(is.null(payload), info = kind)
    module_ci_prot_assert_design_samples_once(payload$design_matrix, design$Run)
    expect_true("tech_rep_group" %in% names(payload$design_matrix))
    expect_false(any(is.na(payload$design_matrix$group) | payload$design_matrix$group == ""))
    expect_setequal(unique(payload$data_cln$Run), payload$design_matrix$Run)
    expect_identical(payload$config_list$deAnalysisParameters$formula_string, "~ 0 + group")
  }

  removed_payload <- module_ci_prot_build_design_payload(
    kind = "two_group",
    removed = "WT_2"
  )
  module_ci_prot_assert_design_samples_once(removed_payload$design_matrix, c("WT_1", "KO_1", "KO_2"))
  expect_false("WT_2" %in% removed_payload$data_cln$Run)
})

test_that("MCI-005.2 contrast matrix covers valid, reversed, invalid, duplicate, and empty contrast states", {
  contrasts <- NULL
  contrasts <- applyProtDesignContrastAppend(contrasts, "KO", "WT")
  contrasts <- applyProtDesignContrastAppend(contrasts, "WT", "KO")
  contrasts_after_duplicate <- applyProtDesignContrastAppend(contrasts, "KO", "WT")
  contrasts_after_invalid <- applyProtDesignContrastAppend(contrasts_after_duplicate, "KO", "KO")
  contrasts_after_empty <- applyProtDesignContrastAppend(contrasts_after_invalid, "", "WT")

  expect_equal(nrow(contrasts), 2L)
  expect_identical(contrasts_after_duplicate, contrasts)
  expect_identical(contrasts_after_invalid, contrasts)
  expect_identical(contrasts_after_empty, contrasts)

  contrast_table <- buildProtDesignSaveResultsContrastsTable(contrasts, "~ 0 + group")
  expect_identical(contrast_table$contrasts, c("groupKO-groupWT", "groupWT-groupKO"))
  expect_identical(contrast_table$friendly_names, c("KO_vs_WT", "WT_vs_KO"))
  expect_false(anyDuplicated(contrast_table$friendly_names) > 0L)

  design <- module_ci_prot_design_matrix("two_group")
  module_ci_prot_assert_contrast_terms_in_model(design, contrast_table, formula = "~ 0 + group")
  expect_null(buildProtDesignSaveResultsContrastsTable(NULL, "~ 0 + group"))

  bad_contrast <- data.frame(contrasts = "groupMISSING-groupWT", stringsAsFactors = FALSE)
  expect_failure(module_ci_prot_assert_contrast_terms_in_model(design, bad_contrast, formula = "~ 0 + group"))
})

test_that("MCI-005.3 sample assignment edge cases accept exact/sanitized/reordered samples and reject unsafe rows", {
  sanitized <- module_ci_prot_design_matrix("two_group", samples = module_ci_prot_design_samples("sanitized"))
  module_ci_prot_assert_design_samples_once(sanitized, module_ci_prot_design_samples("sanitized"))

  reordered <- sanitized[rev(seq_len(nrow(sanitized))), , drop = FALSE]
  module_ci_prot_assert_design_samples_once(reordered, module_ci_prot_design_samples("sanitized"))

  renamed <- applyProtDesignSingleRenameUpdate(
    designMatrix = sanitized,
    dataCln = module_ci_prot_design_data(sanitized$Run),
    originalName = "x123_sample_a",
    newName = "WT_1"
  )
  expect_true("WT_1" %in% renamed$designMatrix$Run)
  expect_true("WT_1" %in% renamed$dataCln$Run)

  suffix_design <- module_ci_prot_design_matrix("two_group", samples = module_ci_prot_design_samples("suffix"))
  suffix_update <- applyProtDesignBulkRenameUpdates(
    designMatrix = suffix_design,
    dataCln = module_ci_prot_design_data(suffix_design$Run),
    selectedSamples = suffix_design$Run,
    newNames = sub("\\.raw$", "", suffix_design$Run)
  )
  module_ci_prot_assert_design_samples_once(suffix_update$designMatrix, sub("\\.raw$", "", suffix_design$Run))

  expected <- module_ci_prot_design_samples("balanced")
  module_ci_prot_assert_design_rejected(
    data.frame(Run = c("WT_1", "wt_2", "KO_1", "KO_2"), stringsAsFactors = FALSE),
    expected
  )
  module_ci_prot_assert_design_rejected(
    data.frame(Run = c("WT_1", "WT_1", "KO_1", "KO_2"), stringsAsFactors = FALSE),
    expected
  )
  module_ci_prot_assert_design_rejected(
    data.frame(Run = c(expected, "EXTRA_1"), stringsAsFactors = FALSE),
    expected
  )
})

test_that("MCI-005.4 existing-design import and builder-created designs produce equivalent artifacts", {
  payload <- module_ci_prot_build_design_payload()
  import_root <- module_ci_prot_write_design_import_pack(tempfile("module-ci-prot-design-import-"), payload)
  workflow_data <- module_ci_prot_design_workflow_data("DIA")
  experiment_paths <- list(source_dir = import_root, results_dir = tempfile("module-ci-prot-design-results-"))
  dir.create(experiment_paths$results_dir)

  artifacts <- resolveProtDesignImportArtifacts(import_root)
  expect_true(artifacts$ok)
  expect_identical(basename(artifacts$designFile), "design_matrix.tab")
  expect_identical(basename(artifacts$dataClnFile), "data_cln.tab")
  expect_identical(basename(artifacts$contrastFile), "contrast_strings.tab")

  imported <- loadProtDesignImportedConfigAndTables(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    designFile = artifacts$designFile,
    dataClnFile = artifacts$dataClnFile,
    contrastFile = artifacts$contrastFile,
    readConfig = function(file) list(
      globalParameters = list(workflow_type = "DIA"),
      deAnalysisParameters = list(formula_string = "~ 0 + group")
    ),
    readTabular = module_ci_prot_read_tsv,
    systemFileFn = function(...) "",
    fileExists = file.exists,
    fileCopy = file.copy,
    downloadFile = function(...) stop("network not allowed in module-CI design import"),
    showNotification = function(...) NULL,
    assignFn = function(...) NULL
  )

  expect_identical(
    imported$importedDesign[c("Run", "factor1", "factor2", "group", "batch", "tech_rep_group")],
    payload$design_matrix[c("Run", "factor1", "factor2", "group", "batch", "tech_rep_group")]
  )
  expect_identical(as.integer(imported$importedDesign$replicates), as.integer(payload$design_matrix$replicates))
  expect_true(all(is.na(imported$importedDesign$factor3)))
  expect_true(all(is.na(payload$design_matrix$factor3)))
  expect_true(all(is.na(imported$importedDesign$tech_reps)))
  expect_true(all(is.na(payload$design_matrix$tech_reps)))
  expect_equal(imported$importedDataCln, payload$data_cln, ignore_attr = TRUE)
  expect_identical(imported$importedContrasts$contrasts, payload$contrasts_tbl$contrasts)

  html <- htmltools::renderTags(mod_prot_design_ui("prot-design"))$html
  expect_match(html, "prot-design-import-existing", fixed = TRUE)
  expect_match(html, "show_import_modal", fixed = TRUE)
})

test_that("MCI-005.5 design artifact writes preserve required filenames, config, manifest, and state digest", {
  payload <- module_ci_prot_build_design_payload()
  source_dir <- tempfile("module-ci-prot-design-source-")
  dir.create(source_dir)
  workflow_data <- module_ci_prot_design_workflow_data("DIA")
  workflow_data$design_matrix <- payload$design_matrix
  workflow_data$data_cln <- payload$data_cln
  workflow_data$contrasts_tbl <- payload$contrasts_tbl

  persistProtDesignBuilderArtifacts(payload, workflow_data, source_dir)

  design_path <- file.path(source_dir, "design_matrix.tab")
  contrast_path <- file.path(source_dir, "contrast_strings.tab")
  manifest_path <- file.path(source_dir, "manifest.json")
  config_path <- file.path(source_dir, "config.ini")

  module_ci_assert_file_nonempty(design_path)
  module_ci_assert_file_nonempty(file.path(source_dir, "data_cln.tab"))
  module_ci_assert_file_nonempty(contrast_path)
  module_ci_assert_file_nonempty(manifest_path)
  module_ci_assert_file_nonempty(config_path)

  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  expect_identical(manifest$design_matrix_path, "design_matrix.tab")
  expect_identical(manifest$contrast_strings_path, "contrast_strings.tab")
  expect_identical(readLines(contrast_path), payload$contrasts_tbl$contrasts)

  digest <- collect_state_digest(
    values = list(selected_omics = "proteomics"),
    workflow_states = list(proteomics = workflow_data)
  )
  expect_identical(digest$workflow_type_per_omic$proteomics, "DIA")
  expect_identical(digest$step_status_per_omic$proteomics$design_matrix, "pending")
})

test_that("MCI-005.6 design artifacts pass lightweight downstream DA setup validation before workflow advancement", {
  payload <- module_ci_prot_build_design_payload()
  module_ci_prot_validate_da_setup(
    data = payload$data_cln,
    design = payload$design_matrix,
    contrasts = payload$contrasts_tbl,
    formula = "~ 0 + group"
  )

  invalid_design <- payload$design_matrix[-1, , drop = FALSE]
  module_ci_prot_assert_design_rejected(invalid_design, unique(as.character(payload$data_cln$Run)))

  invalid_contrast <- data.frame(contrasts = "groupKO-groupMISSING", stringsAsFactors = FALSE)
  model_cols <- colnames(stats::model.matrix(~ 0 + group, data = payload$design_matrix))
  expect_false(all(module_ci_prot_parse_contrast_terms(invalid_contrast$contrasts) %in% model_cols))
})
