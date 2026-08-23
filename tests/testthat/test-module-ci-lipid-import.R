test_that("MCI-018.1 LipidSearch matrix preserves routing, mappings, and boundaries", {
  lipidsearch_cases <- list(
    LCMS_Pos = list(assay = module_ci_lipid_lipidsearch_assay(), assay_name = "LCMS_Pos", annotation_col = "LipidClass"),
    LCMS_Neg = list(assay = module_ci_lipid_lipidsearch_assay(zero_values = TRUE), assay_name = "LCMS_Neg", annotation_col = "LipidClass"),
    GCMS_named = list(assay = module_ci_lipid_lipidsearch_assay(extra_metadata = TRUE), assay_name = "GCMS", annotation_col = "LipidClass"),
    missing_lipid_class = list(assay = module_ci_lipid_lipidsearch_assay(include_lipid_class = FALSE), assay_name = "LCMS_Pos_NoClass", annotation_col = ""),
    duplicate_lipid_ids = list(assay = module_ci_lipid_lipidsearch_assay(duplicate_lipid = TRUE), assay_name = "LCMS_Pos_Duplicate", annotation_col = "LipidClass")
  )

  for (case_name in names(lipidsearch_cases)) {
    case <- lipidsearch_cases[[case_name]]
    paths <- module_ci_lipid_import_paths()
    assay_file <- module_ci_lipid_write_table(case$assay, file.path(paths$raw_dir, paste0(case_name, ".txt")))
    preview <- module_ci_lipid_preview(assay_file)

    expect_identical(preview$detectedFormat, "lipidsearch", info = case_name)
    expect_setequal(preview$importResult$sample_columns, module_ci_lipid_sample_cols())
    expect_false("BatchFlag" %in% preview$importResult$sample_columns, info = case_name)

    run <- module_ci_lipid_import_run(
      assay = preview$importResult$data,
      paths = paths,
      assay_name = case$assay_name,
      assay_file = assay_file,
      vendor_format = "lipidsearch",
      detected_format = preview$detectedFormat,
      lipid_id_col = "LipidName",
      annotation_col = case$annotation_col,
      sample_cols = preview$importResult$sample_columns
    )

    expect_identical(run$result$status, "success", info = case_name)
    expect_identical(run$workflow$data_format, "lipidsearch", info = case_name)
    expect_identical(names(run$workflow$data_tbl), case$assay_name, info = case_name)
    expect_true(run$result$artifactResult$written, info = case_name)
    module_ci_assert_numeric_finite(run$workflow$data_tbl[[case$assay_name]], module_ci_lipid_sample_cols(), allow_na = TRUE)
  }

  missing_id <- module_ci_lipid_import_run(
    assay = module_ci_lipid_lipidsearch_assay(include_lipid_name = FALSE),
    lipid_id_col = "LipidName"
  )
  expect_null(missing_id$result)
  expect_match(missing_id$error_messages[[1]], "Lipid ID column not found: LipidName", fixed = TRUE)
  expect_null(missing_id$workflow$data_tbl)
  expect_identical(missing_id$workflow$tab_status$setup_import, "incomplete")
})

test_that("MCI-018.2 MS-DIAL matrix preserves annotation, duplicate, and intensity semantics", {
  msdial <- module_ci_lipid_msdial_assay(include_annotation = TRUE, duplicate_feature = TRUE, zero_missing = TRUE)
  paths <- module_ci_lipid_import_paths()
  assay_file <- module_ci_lipid_write_table(msdial, file.path(paths$raw_dir, "msdial_positive.tsv"))

  preview <- module_ci_lipid_preview(assay_file)
  expect_identical(preview$detectedFormat, "msdial")
  expect_identical(preview$updates$lipidId$selected, "Peak ID")
  expect_identical(preview$updates$annotation$selected, "Name")
  expect_setequal(preview$importResult$sample_columns, module_ci_lipid_sample_cols())
  expect_false("Total Score" %in% preview$importResult$sample_columns)

  run <- module_ci_lipid_import_run(
    assay = preview$importResult$data,
    paths = paths,
    assay_file = assay_file,
    vendor_format = "msdial",
    detected_format = preview$detectedFormat,
    lipid_id_col = "Peak ID",
    annotation_col = "Name",
    sample_cols = preview$importResult$sample_columns
  )

  expect_identical(run$result$status, "success")
  expect_identical(run$workflow$data_format, "msdial")
  expect_true(any(run$workflow$data_tbl$LCMS_Pos$WT_1 == 0))
  expect_true(any(is.na(run$workflow$data_tbl$LCMS_Pos$WT_2)))
  expect_true(any(grepl("duplicate lipid IDs", run$result$validationResult$warnings, fixed = TRUE)))

  no_annotation <- module_ci_lipid_msdial_assay(include_annotation = FALSE, duplicate_feature = FALSE)
  no_annotation_file <- module_ci_lipid_write_table(no_annotation, file.path(paths$raw_dir, "msdial_no_annotation.tsv"))
  no_annotation_preview <- module_ci_lipid_preview(no_annotation_file)
  no_annotation_run <- module_ci_lipid_import_run(
    assay = no_annotation_preview$importResult$data,
    paths = module_ci_lipid_import_paths(),
    assay_file = no_annotation_file,
    vendor_format = "msdial",
    detected_format = no_annotation_preview$detectedFormat,
    lipid_id_col = "Peak ID",
    annotation_col = "",
    sample_cols = no_annotation_preview$importResult$sample_columns
  )

  expect_identical(no_annotation_run$result$status, "success")
  expect_null(no_annotation_run$workflow$column_mapping$annotation_col)
})

test_that("MCI-018.3 custom and vendor overrides fail closed without fallback", {
  custom <- module_ci_lipid_custom_assay(sample_prefix = "sample")
  expect_identical(
    resolveLipidImportEffectiveColumn(custom, "custom", "ignored", "LIPID_ID"),
    "lipid_id"
  )
  expect_identical(
    resolveLipidImportSampleColumns(custom, vendorFormat = "custom", sampleColsPattern = "^sample", excludeNormalized = TRUE),
    c("Sample A", "Sample B")
  )

  fallback_called <- character()
  expect_error(
    loadLipidImportAssayPreview(
      assay1File = "unknown_vendor.tsv",
      readHeaders = function(path) c("lipid_id", "annotation", "WT_1"),
      detectFormat = function(headers, filename) {
        list(format = "unsupported_vendor", confidence = 0.11)
      },
      importMsdial = function(path) {
        fallback_called <<- c(fallback_called, path)
      },
      importLipidSearch = function(...) stop("LipidSearch fallback should not be used")
    ),
    class = "multischolar_format_unknown"
  )
  expect_identical(fallback_called, character())

  expect_error(
    module_ci_lipid_import_run(
      assay = module_ci_lipid_msdial_assay(),
      vendor_format = "lipidsearch",
      detected_format = "msdial",
      lipid_id_col = "Peak ID",
      annotation_col = "Name",
      sample_cols = module_ci_lipid_sample_cols()
    ),
    class = "multischolar_format_mismatch"
  )

  custom_override <- module_ci_lipid_import_run(
    assay = module_ci_lipid_lipidsearch_assay(),
    vendor_format = "custom",
    detected_format = "lipidsearch",
    lipid_id_col = "LipidName",
    annotation_col = "LipidClass",
    sample_cols = module_ci_lipid_sample_cols()
  )
  expect_identical(custom_override$result$status, "success")
  expect_identical(custom_override$workflow$data_format, "custom")
})

test_that("MCI-018.4 multi-assay routing matrix preserves assay2 parity and rejects ambiguous routes", {
  primary <- module_ci_lipid_lipidsearch_assay()
  secondary <- module_ci_lipid_lipidsearch_assay(zero_values = TRUE)
  paths <- module_ci_lipid_import_paths()
  secondary_file <- module_ci_lipid_write_table(secondary, file.path(paths$raw_dir, "lcms_neg.txt"))

  pair <- module_ci_lipid_import_run(
    assay = primary,
    paths = paths,
    assay_name = "LCMS_Pos",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Neg"
  )
  expect_identical(pair$result$status, "success")
  expect_identical(names(pair$workflow$data_tbl), c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(unname(pair$workflow$processing_log$setup_import$n_lipids), c(3L, 3L))

  absent <- module_ci_lipid_import_run(
    assay = primary,
    assay_name = "GCMS",
    assay2_file = NULL,
    assay2_name = ""
  )
  expect_identical(absent$result$status, "success")
  expect_identical(names(absent$workflow$data_tbl), "GCMS")

  unusual <- module_ci_lipid_import_run(
    assay = primary,
    assay_name = "GCMS-2026_A",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Pos_2026"
  )
  expect_identical(unusual$result$status, "success")
  expect_identical(names(unusual$workflow$data_tbl), c("GCMS-2026_A", "LCMS_Pos_2026"))

  duplicate <- module_ci_lipid_import_run(
    assay = primary,
    assay_name = "LCMS_Pos",
    assay2_file = secondary_file,
    assay2_name = "LCMS_Pos"
  )
  expect_null(duplicate$result)
  expect_match(duplicate$error_messages[[1]], "Duplicate lipidomics assay names", fixed = TRUE)
  expect_null(duplicate$workflow$data_tbl)

  empty_assay2 <- file.path(paths$raw_dir, "empty_assay2.txt")
  writeLines(character(), empty_assay2)
  empty <- module_ci_lipid_import_run(
    assay = primary,
    assay2_file = empty_assay2,
    assay2_name = "LCMS_Neg"
  )
  expect_null(empty$result)
  expect_null(empty$workflow$data_tbl)

  reader_calls <- character()
  reader <- resolveLipidImportSecondAssayReader(
    "lipidsearch",
    importMsdial = function(...) {
      reader_calls <<- c(reader_calls, "msdial")
      list(data = data.frame())
    },
    importLipidSearch = function(...) {
      reader_calls <<- c(reader_calls, "lipidsearch")
      list(data = data.frame())
    }
  )
  callLipidImportSecondAssayReader(reader, "assay2.txt", lipidIdCol = "LipidName", annotationCol = "LipidClass")
  expect_identical(reader_calls, "lipidsearch")
})

test_that("MCI-018.5 invalid input matrix fails before downstream state advances", {
  valid <- module_ci_lipid_lipidsearch_assay()
  invalid_cases <- list(
    missing_id_column = list(assay = valid, lipid_id_col = "MissingLipid", sample_cols = module_ci_lipid_sample_cols()),
    no_sample_columns = list(assay = valid, lipid_id_col = "LipidName", sample_cols = character()),
    non_numeric_samples = list(assay = module_ci_lipid_lipidsearch_assay(non_numeric_samples = TRUE), lipid_id_col = "LipidName", sample_cols = module_ci_lipid_sample_cols()),
    mismatched_sample_names = list(assay = valid, lipid_id_col = "LipidName", sample_cols = c("Missing_1", "Missing_2")),
    no_numeric_data = list(assay = data.frame(LipidName = c("L1", "L2"), WT_1 = c("low", "high")), lipid_id_col = "LipidName", sample_cols = "WT_1")
  )

  for (case_name in names(invalid_cases)) {
    case <- invalid_cases[[case_name]]
    workflow <- module_ci_lipid_import_workflow()
    run <- module_ci_lipid_import_run(
      assay = case$assay,
      workflow = workflow,
      lipid_id_col = case$lipid_id_col,
      sample_cols = case$sample_cols
    )

    expect_null(run$result, info = case_name)
    expect_null(workflow$data_tbl, info = case_name)
    expect_null(workflow$column_mapping, info = case_name)
    expect_identical(workflow$tab_status$setup_import, "incomplete", info = case_name)
    expect_identical(workflow$tab_status$design_matrix, "disabled", info = case_name)
    expect_null(workflow$processing_log$setup_import, info = case_name)
  }

  malformed_paths <- module_ci_lipid_import_paths()
  malformed_path <- module_ci_lipid_write_malformed(file.path(malformed_paths$raw_dir, "malformed.tsv"))
  malformed_preview <- loadLipidImportAssayPreview(
    malformed_path,
    vendorFormat = "custom"
  )
  malformed_run <- module_ci_lipid_import_run(
    assay = malformed_preview$importResult$data,
    paths = malformed_paths,
    assay_file = malformed_path,
    vendor_format = "custom",
    detected_format = malformed_preview$detectedFormat,
    lipid_id_col = "LipidName",
    sample_cols = malformed_preview$importResult$sample_columns
  )
  expect_null(malformed_run$result)
  expect_null(malformed_run$workflow$data_tbl)
})

test_that("MCI-018.6 artifact and state digest matrix preserves downstream design handoff", {
  paths <- module_ci_lipid_import_paths()
  primary <- module_ci_lipid_lipidsearch_assay()
  secondary <- module_ci_lipid_lipidsearch_assay(zero_values = TRUE)
  primary_file <- module_ci_lipid_write_table(primary, file.path(paths$raw_dir, "lcms_pos.txt"))
  secondary_file <- module_ci_lipid_write_table(secondary, file.path(paths$raw_dir, "lcms_neg.txt"))

  run <- module_ci_lipid_import_run(
    assay = primary,
    paths = paths,
    assay_name = "LCMS_Pos",
    assay_file = primary_file,
    assay2_file = secondary_file,
    assay2_name = "LCMS_Neg"
  )
  expect_identical(run$result$status, "success")

  artifacts <- module_ci_lipid_import_artifacts(paths, c("LCMS_Pos", "LCMS_Neg"))
  expect_true(all(file.exists(artifacts)))
  cleaned <- module_ci_assert_table_schema(
    artifacts[["data_cln_LCMS_Pos"]],
    required_columns = c("LipidName", module_ci_lipid_sample_cols())
  )
  module_ci_assert_checksum_unchanged(primary[c("LipidName", module_ci_lipid_sample_cols())], cleaned[c("LipidName", module_ci_lipid_sample_cols())])

  manifest <- readLines(artifacts[["assay_manifest"]], warn = FALSE)
  expect_identical(manifest, c("LCMS_Pos", "LCMS_Neg"))

  column_mapping <- jsonlite::read_json(artifacts[["column_mapping"]], simplifyVector = TRUE)
  expect_identical(column_mapping$lipid_id_col, "LipidName")
  expect_identical(column_mapping$annotation_col, "LipidClass")
  expect_identical(column_mapping$sample_columns, module_ci_lipid_sample_cols())

  summary <- module_ci_assert_table_schema(
    artifacts[["import_summary"]],
    required_columns = c("assay_name", "feature_count", "sample_count", "sample_columns"),
    min_rows = 2L
  )
  expect_identical(summary$assay_name, c("LCMS_Pos", "LCMS_Neg"))
  expect_identical(summary$sample_count, c(4L, 4L))

  digest <- module_ci_lipid_import_state_digest(run$workflow)
  expect_identical(digest$setup_import, "complete")
  expect_identical(digest$design_matrix, "pending")
  expect_identical(digest$workflow_types, "lipidomics_standard")
  expect_true(digest$processing_log$artifacts$written)
  expect_true(length(digest$processing_log$artifacts$source_copies) >= 2L)
  expect_true(all(file.exists(digest$processing_log$artifacts$source_copies)))

  withr::with_options(
    list(multischolar.test_mode = TRUE),
    {
      html <- htmltools::renderTags(mod_lipid_import_ui("module-ci"))$html
      expect_match(html, "assay1_file_std", fixed = TRUE)
      expect_match(html, "assay2_file_std", fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-assay1-file"', fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-assay2-file"', fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-process"', fixed = TRUE)
      expect_match(html, 'data-testid="lipid-import-vendor-format"', fixed = TRUE)
    }
  )
})
