test_that("MCI-012.1 custom import matrix preserves mappings and numeric coercion", {
  cases <- list(
    annotation_present = module_ci_metab_custom_assay(annotation = TRUE, numeric_as_character = FALSE, sample_prefix = "WT"),
    annotation_absent = module_ci_metab_custom_assay(annotation = FALSE, numeric_as_character = FALSE, sample_prefix = "WT"),
    regex_variation = module_ci_metab_custom_assay(annotation = TRUE, numeric_as_character = TRUE, sample_prefix = "KO")
  )

  for (case_name in names(cases)) {
    assay <- cases[[case_name]]
    annotation_col <- if ("Annotation" %in% names(assay)) "Annotation" else ""
    sample_pattern <- if (identical(case_name, "regex_variation")) "^(Sample|QC)-" else "^(WT|KO)_"
    sample_cols <- grep(sample_pattern, names(assay), value = TRUE)

    run <- module_ci_metab_import_run(
      assay = assay,
      assay_name = paste0("LCMS_Pos_", case_name),
      metabolite_col = "Feature.Name",
      annotation_col = annotation_col,
      sample_cols = sample_cols,
      sanitize_names = FALSE
    )

    expect_identical(run$result$status, "success", info = case_name)
    expect_identical(run$workflow$data_format, "custom", info = case_name)
    expect_identical(run$workflow$column_mapping$metabolite_id_col, "Feature.Name", info = case_name)
    expect_identical(run$workflow$column_mapping$sample_columns, sample_cols, info = case_name)
    expect_true(all(vapply(run$workflow$data_tbl[[1]][sample_cols], is.numeric, logical(1))), info = case_name)
    expect_true(run$result$artifactResult$written, info = case_name)
  }
})

test_that("MCI-012.2 MS-DIAL-style import matrix preserves metadata and edge intensities", {
  assay <- module_ci_metab_msdial_assay(include_name = TRUE, duplicate_feature = TRUE, zero_missing = TRUE)
  paths <- module_ci_metab_import_paths()
  assay_file <- module_ci_metab_write_table(assay, file.path(paths$raw_dir, "msdial_positive.tsv"))

  selected <- prepareMetabImportAssaySelectionState(assay_file)
  expect_identical(selected$formatInfo$format, "msdial")
  expect_identical(selected$selectedMetaboliteId, "Peak ID")
  expect_identical(selected$selectedAnnotation, "Name")
  expect_setequal(selected$importResult$sample_columns, c("WT_1", "WT_2", "KO_1", "KO_2"))
  expect_false("Total Score" %in% selected$importResult$sample_columns)

  validation <- validateMetabColumnMapping(
    data = selected$importResult$data,
    metabolite_id_column = "Peak ID",
    sample_columns = selected$importResult$sample_columns
  )
  expect_true(validation$valid)
  expect_true(any(grepl("duplicate metabolite IDs", validation$warnings, fixed = TRUE)))

  run <- module_ci_metab_import_run(
    assay = selected$importResult$data,
    paths = paths,
    assay_name = "LCMS_Pos",
    assay_file = assay_file,
    vendor_format = "msdial",
    detected_format = selected$formatInfo$format,
    metabolite_col = "Peak ID",
    annotation_col = "Name",
    sample_cols = selected$importResult$sample_columns
  )

  expect_identical(run$result$status, "success")
  expect_identical(run$workflow$data_format, "msdial")
  expect_true(any(is.na(run$workflow$data_tbl$LCMS_Pos$WT_2)))
  expect_true(any(run$workflow$data_tbl$LCMS_Pos$WT_1 == 0))
  expect_identical(run$workflow$processing_log$setup_import$n_assays, 1L)
})

test_that("MCI-012.3 multi-assay imports preserve assay names and reject ambiguous names", {
  primary <- module_ci_metab_custom_assay()
  secondary <- module_ci_metab_custom_assay(annotation = TRUE)
  paths <- module_ci_metab_import_paths()
  secondary_file <- module_ci_metab_write_table(secondary, file.path(paths$raw_dir, "lcms_neg.tsv"))

  lcms_pair <- module_ci_metab_import_run(
    assay = primary,
    paths = paths,
    assay_name = "LCMS_Pos",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Neg"
  )
  expect_identical(lcms_pair$result$status, "success")
  expect_identical(names(lcms_pair$workflow$data_tbl), c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(unname(lcms_pair$workflow$processing_log$setup_import$n_metabolites), c(3L, 3L))

  gc_only <- module_ci_metab_import_run(
    assay = primary,
    assay_name = "GCMS",
    assay2_file = secondary_file,
    assay2_name = ""
  )
  expect_identical(gc_only$result$status, "success")
  expect_identical(names(gc_only$workflow$data_tbl), "GCMS")

  unusual <- module_ci_metab_import_run(
    assay = primary,
    assay_name = "GCMS-2026_A",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Pos_2026"
  )
  expect_identical(unusual$result$status, "success")
  expect_identical(names(unusual$workflow$data_tbl), c("GCMS-2026_A", "LCMS_Pos_2026"))

  duplicate <- module_ci_metab_import_run(
    assay = primary,
    assay_name = "LCMS_Pos",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Pos"
  )
  expect_identical(duplicate$result$status, "error")
  expect_match(duplicate$result$error$message, "Duplicate metabolomics assay names", fixed = TRUE)
  expect_null(duplicate$workflow$data_tbl)
  expect_identical(duplicate$workflow$tab_status$setup_import, "incomplete")

  unsafe <- module_ci_metab_import_run(
    assay = primary,
    assay_name = "LCMS/Pos"
  )
  expect_identical(unsafe$result$status, "error")
  expect_match(unsafe$result$error$message, "Assay names must not contain path separators", fixed = TRUE)
})

test_that("MCI-012.4 invalid import matrix fails before workflow state advances", {
  valid_assay <- module_ci_metab_custom_assay()
  invalid_cases <- list(
    missing_id_column = list(assay = valid_assay, metabolite_col = "Missing.ID", sample_cols = c("WT_1", "WT_2")),
    no_sample_columns = list(assay = valid_assay, metabolite_col = "Feature.Name", sample_cols = character()),
    non_numeric_samples = list(
      assay = within(valid_assay, WT_1 <- c("ok", "bad", "values")),
      metabolite_col = "Feature.Name",
      sample_cols = c("WT_1", "WT_2")
    ),
    mismatched_pattern = list(assay = valid_assay, metabolite_col = "Feature.Name", sample_cols = c("Missing_1")),
    malformed_payload = list(assay = data.frame(Bad = c("{", "}")), metabolite_col = "Feature.Name", sample_cols = "WT_1")
  )

  for (case_name in names(invalid_cases)) {
    case <- invalid_cases[[case_name]]
    workflow <- module_ci_metab_import_workflow()
    run <- module_ci_metab_import_run(
      assay = case$assay,
      workflow = workflow,
      metabolite_col = case$metabolite_col,
      sample_cols = case$sample_cols
    )

    expect_identical(run$result$status, "error", info = case_name)
    expect_null(workflow$data_tbl, info = case_name)
    expect_null(workflow$column_mapping, info = case_name)
    expect_identical(workflow$tab_status$setup_import, "incomplete", info = case_name)
    expect_null(workflow$processing_log$setup_import, info = case_name)
  }
})

test_that("MCI-012.5 browser import UI exposes standard upload controls in test mode", {
  withr::with_options(
    list(multischolar.test_mode = TRUE),
    {
      html <- htmltools::renderTags(mod_metab_import_ui("module-ci"))$html
      expect_match(html, "assay1_file_std", fixed = TRUE)
      expect_match(html, "assay2_file_std", fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-assay1-file"', fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-assay2-file"', fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-process"', fixed = TRUE)
      expect_match(html, 'data-testid="metab-import-vendor-format"', fixed = TRUE)
    }
  )
})

test_that("MCI-012.6 import artifacts and state digest payload are complete", {
  assay <- module_ci_metab_custom_assay(numeric_as_character = TRUE)
  paths <- module_ci_metab_import_paths()
  run <- module_ci_metab_import_run(
    assay = assay,
    paths = paths,
    assay_name = "LCMS_Pos",
    sanitize_names = FALSE
  )

  expect_identical(run$result$status, "success")
  artifacts <- module_ci_metab_import_artifacts(paths, "LCMS_Pos")
  expect_true(all(file.exists(artifacts)))

  cleaned <- read.delim(artifacts[["data_cln_LCMS_Pos"]], check.names = FALSE)
  expect_setequal(names(cleaned), names(run$workflow$data_tbl$LCMS_Pos))
  expect_true(is.numeric(cleaned$WT_1))

  manifest <- readLines(artifacts[["assay_manifest"]], warn = FALSE)
  expect_identical(manifest, "LCMS_Pos")

  column_mapping <- jsonlite::read_json(artifacts[["column_mapping"]], simplifyVector = TRUE)
  expect_identical(column_mapping$metabolite_id_col, "Feature.Name")
  expect_identical(column_mapping$sample_columns, c("WT_1", "WT_2", "KO_1", "KO_2"))

  summary <- read.delim(artifacts[["import_summary"]], check.names = FALSE)
  expect_identical(summary$assay_name, "LCMS_Pos")
  expect_identical(summary$feature_count, 3L)
  expect_identical(summary$sample_count, 4L)

  digest_payload <- list(
    setup_import = run$workflow$tab_status$setup_import,
    column_mapping = run$workflow$column_mapping,
    processing_log = run$workflow$processing_log$setup_import
  )
  expect_identical(digest_payload$setup_import, "complete")
  expect_true(digest_payload$processing_log$artifacts$written)
  expect_true(length(digest_payload$processing_log$artifacts$source_copies) >= 1L)
  expect_true(file.exists(digest_payload$processing_log$artifacts$source_copies[[1]]))
})
