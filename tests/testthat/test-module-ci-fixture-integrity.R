library(testthat)

test_that("module CI fixture pack catalog validates required omics, formats, and classes", {
  catalog <- read_module_ci_fixture_packs()
  expected_formats <- module_ci_expected_import_formats()
  expected_pack_count <- sum(lengths(expected_formats)) * length(.MODULE_CI_FIXTURE_CLASSES)

  expect_identical(catalog$schema_version, "1.0.0")
  expect_identical(length(catalog$packs), expected_pack_count)
  expect_setequal(
    vapply(catalog$fixture_classes, `[[`, character(1), "fixture_class"),
    .MODULE_CI_FIXTURE_CLASSES
  )

  for (omic in names(expected_formats)) {
    omic_packs <- module_ci_fixture_packs(catalog, omic = omic)
    expect_setequal(
      unique(vapply(omic_packs, `[[`, character(1), "fixture_class")),
      .MODULE_CI_FIXTURE_CLASSES
    )
    expect_true(any(vapply(omic_packs, function(pack) {
      pack$case_category %in% c("multi_assay", "branch_specific")
    }, logical(1))))

    for (import_format in expected_formats[[omic]]) {
      format_packs <- module_ci_fixture_packs(catalog, omic = omic, import_format = import_format)
      expect_setequal(
        vapply(format_packs, `[[`, character(1), "fixture_class"),
        .MODULE_CI_FIXTURE_CLASSES
      )
    }
  }
})

test_that("module CI fixture files match oracle sidecar contracts", {
  catalog <- read_module_ci_fixture_packs()

  for (pack in catalog$packs) {
    data <- module_ci_read_pack_table(pack, "data_path")
    design <- module_ci_read_pack_table(pack, "design_path")
    oracle <- module_ci_read_pack_oracle(pack)
    feature_col <- oracle$schema$feature_id_col
    sample_columns <- module_ci_as_character(oracle$schema$sample_columns)
    data_checksum_columns <- module_ci_as_character(oracle$input_checksums$data_columns)
    design_checksum_columns <- module_ci_as_character(oracle$input_checksums$design_columns)

    expect_identical(oracle$fixture_class, pack$fixture_class)
    expect_identical(oracle$omic, pack$omic)
    expect_identical(oracle$import_format, pack$import_format)
    expect_identical(oracle$ci_lane, pack$ci_lane)
    expect_identical(oracle$push_safe, pack$push_safe)

    module_ci_assert_table_schema(
      module_ci_pack_path(pack, "data_path"),
      required_columns = module_ci_as_character(oracle$schema$required_columns),
      min_rows = 1L
    )
    module_ci_assert_table_schema(
      module_ci_pack_path(pack, "design_path"),
      required_columns = c("sample", "group", "batch"),
      min_rows = 1L
    )
    module_ci_assert_checksum_equal(
      data,
      oracle$input_checksums$data_selected_columns_md5,
      data_checksum_columns
    )
    module_ci_assert_checksum_equal(
      design,
      oracle$input_checksums$design_selected_columns_md5,
      design_checksum_columns
    )
    module_ci_assert_sample_identity(sample_columns, sample_columns)
    module_ci_assert_feature_identity(
      module_ci_as_character(oracle$feature_ids),
      as.character(data[[feature_col]])
    )
    module_ci_assert_assay_provenance(data, module_ci_as_character(oracle$assay_names))
    module_ci_assert_expected_missingness(data, oracle$expected_missingness_by_column)
    module_ci_assert_nonfinite_columns(
      data,
      oracle$expected_nonfinite_columns,
      columns = sample_columns
    )

    if (identical(pack$fixture_class, "invalid_design")) {
      expect_true(anyDuplicated(as.character(design$sample)) > 0L)
      expect_false(setequal(as.character(design$sample), sample_columns))
      expect_true(all(c("design_sample_mismatch", "duplicate_design_samples") %in%
        module_ci_as_character(oracle$expected_validation_errors)))
    } else {
      module_ci_assert_design_alignment(design, sample_columns)
    }

    if (identical(pack$fixture_class, "duplicates")) {
      expect_failure(module_ci_assert_no_duplicate_keys(data, feature_col))
      expect_true(length(module_ci_as_character(oracle$duplicate_feature_ids)) > 0L)
    } else {
      module_ci_assert_no_duplicate_keys(data, feature_col)
    }

    if (!identical(pack$fixture_class, "nonfinite")) {
      module_ci_assert_numeric_finite(data, columns = sample_columns, allow_na = TRUE)
    }
  }
})

test_that("module CI corruption sentinels fail when fixtures drift", {
  catalog <- read_module_ci_fixture_packs()
  pack <- module_ci_fixture_packs(
    catalog,
    omic = "proteomics",
    import_format = "dia",
    fixture_class = "happy_path"
  )[[1]]
  oracle <- module_ci_read_pack_oracle(pack)
  data <- module_ci_read_pack_table(pack, "data_path")
  data_checksum_columns <- module_ci_as_character(oracle$input_checksums$data_columns)

  expect_no_error(module_ci_assert_checksum_equal(
    data,
    oracle$input_checksums$data_selected_columns_md5,
    data_checksum_columns
  ))

  drifted <- data
  drifted[[data_checksum_columns[[2]]]][[1]] <- drifted[[data_checksum_columns[[2]]]][[1]] + 1
  expect_failure(module_ci_assert_checksum_equal(
    drifted,
    oracle$input_checksums$data_selected_columns_md5,
    data_checksum_columns
  ))

  expect_failure(module_ci_assert_sample_identity(
    before = module_ci_as_character(oracle$sample_order),
    after = module_ci_as_character(oracle$sample_order)[-1]
  ))
  expect_failure(module_ci_assert_expected_missingness(
    data,
    stats::setNames(as.list(rep(99L, length(data_checksum_columns) - 1L)), data_checksum_columns[-1])
  ))
})

test_that("module CI fixture generator reproduces catalog and sidecar checksums", {
  skip_if(Sys.which("Rscript") == "", "Rscript not available")
  root <- tempfile("module-ci-regenerated-")
  dir.create(root)

  generator <- file.path(module_ci_fixture_root(), "generate-fixtures.R")
  output <- system2("Rscript", c(generator, "--root", root), stdout = TRUE, stderr = TRUE)
  exit_status <- attr(output, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  expect_identical(as.integer(exit_status), 0L, info = paste(output, collapse = "\n"))
  expect_true(file.exists(file.path(root, "fixture_packs.json")))

  regenerated <- read_module_ci_fixture_packs(file.path(root, "fixture_packs.json"), validate = TRUE)
  current <- read_module_ci_fixture_packs()
  expect_identical(module_ci_pack_ids(regenerated$packs), module_ci_pack_ids(current$packs))

  pack <- regenerated$packs[["MCI-002-proteomics-dia-happy-path"]]
  oracle <- module_ci_read_pack_oracle(pack, fixture_root = root)
  data <- module_ci_read_pack_table(pack, "data_path", fixture_root = root)
  module_ci_assert_checksum_equal(
    data,
    oracle$input_checksums$data_selected_columns_md5,
    module_ci_as_character(oracle$input_checksums$data_columns)
  )
})

test_that("module CI push and nightly fixture labels are explicit and bounded", {
  catalog <- read_module_ci_fixture_packs()
  push_packs <- module_ci_fixture_packs(catalog, ci_lane = "push", push_safe = TRUE)
  nightly_packs <- module_ci_fixture_packs(catalog, ci_lane = "nightly", push_safe = FALSE)

  expect_true(length(push_packs) > 0L)
  expect_true(length(nightly_packs) > 0L)
  expect_false(any(vapply(push_packs, function(pack) {
    identical(pack$fixture_class, "large_enough_for_plots")
  }, logical(1))))
  expect_true(all(vapply(push_packs, function(pack) {
    pack$estimated_rows <= 8L && pack$estimated_columns <= 12L
  }, logical(1))))
  expect_true(all(vapply(nightly_packs, function(pack) {
    identical(pack$fixture_class, "large_enough_for_plots") &&
      identical(pack$ci_lane, "nightly") &&
      identical(pack$push_safe, FALSE)
  }, logical(1))))
})
