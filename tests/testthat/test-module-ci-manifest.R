library(testthat)

test_that("module CI manifest validates and returns named scenarios", {
  manifest <- read_module_ci_manifest()

  expect_identical(manifest$schema_version, "1.0.0")
  expect_true(length(manifest$scenarios) >= 4L)
  expect_setequal(
    names(manifest$scenarios),
    module_ci_scenario_ids(manifest$scenarios)
  )
  expect_true(all(vapply(manifest$scenarios, function(scenario) {
    file.exists(file.path(module_ci_fixture_root(), scenario$fixture$path))
  }, logical(1))))
})

test_that("module CI manifest filters scenarios by omic, module, runtime, and ticket", {
  manifest <- read_module_ci_manifest()

  prot <- module_ci_scenarios(manifest, omic = "proteomics")
  expect_setequal(
    module_ci_scenario_ids(prot),
    c(
      "MCI-004.1-proteomics-import-schema-smoke",
      "MCI-005.1-proteomics-design-schema-smoke",
      "MCI-006.1-proteomics-peptide-qc-schema-smoke",
      "MCI-007.1-proteomics-protein-qc-schema-smoke",
      "MCI-008.1-proteomics-normalization-schema-smoke",
      "MCI-009.1-proteomics-da-schema-smoke",
      "MCI-010.1-proteomics-enrichment-schema-smoke",
      "MCI-011.1-proteomics-summary-report-schema-smoke"
    )
  )

  import_unit <- module_ci_scenarios(
    manifest,
    module_family = "import",
    runtime = "unit-contract"
  )
  expect_setequal(
    module_ci_scenario_ids(import_unit),
    c(
      "MCI-004.1-proteomics-import-schema-smoke",
      "MCI-012.1-metabolomics-import-schema-smoke",
      "MCI-018.1-lipidomics-import-schema-smoke"
    )
  )

  ticket <- module_ci_scenarios(manifest, ticket_id = "MCI-001")
  expect_identical(module_ci_scenario_ids(ticket), "MCI-001.1-foundation-manifest-smoke")

  metab_qc <- module_ci_scenarios(
    manifest,
    omic = "metabolomics",
    module_family = "qc",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(metab_qc),
    "MCI-014.1-metabolomics-qc-schema-smoke"
  )

  metab_norm <- module_ci_scenarios(
    manifest,
    omic = "metabolomics",
    module_family = "normalization",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(metab_norm),
    "MCI-015.1-metabolomics-normalization-schema-smoke"
  )

  metab_da <- module_ci_scenarios(
    manifest,
    omic = "metabolomics",
    module_family = "differential_abundance",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(metab_da),
    "MCI-016.1-metabolomics-da-schema-smoke"
  )

  metab_summary <- module_ci_scenarios(
    manifest,
    omic = "metabolomics",
    module_family = "summary_report",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(metab_summary),
    "MCI-017.1-metabolomics-summary-report-schema-smoke"
  )

  lipid_import <- module_ci_scenarios(
    manifest,
    omic = "lipidomics",
    module_family = "import",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(lipid_import),
    "MCI-018.1-lipidomics-import-schema-smoke"
  )
  expect_true("module-ci-lipid-import" %in% lipid_import[[1]]$required_tests)

  lipid_design <- module_ci_scenarios(
    manifest,
    omic = "lipidomics",
    module_family = "design",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(lipid_design),
    "MCI-019.1-lipidomics-design-schema-smoke"
  )
  expect_true("module-ci-lipid-design" %in% lipid_design[[1]]$required_tests)

  lipid_qc <- module_ci_scenarios(
    manifest,
    omic = "lipidomics",
    module_family = "qc",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(lipid_qc),
    "MCI-020.1-lipidomics-qc-schema-smoke"
  )
  expect_true("module-ci-lipid-qc" %in% lipid_qc[[1]]$required_tests)

  lipid_norm <- module_ci_scenarios(
    manifest,
    omic = "lipidomics",
    module_family = "normalization",
    runtime = "unit-contract"
  )
  expect_identical(
    module_ci_scenario_ids(lipid_norm),
    "MCI-021.1-lipidomics-normalization-schema-smoke"
  )
  expect_true("module-ci-lipid-norm" %in% lipid_norm[[1]]$required_tests)
})

test_that("module CI manifest rejects duplicate scenario IDs and missing fixtures", {
  manifest <- read_module_ci_manifest(validate = FALSE)
  manifest$scenarios[[2]]$scenario_id <- manifest$scenarios[[1]]$scenario_id
  expect_error(validate_module_ci_manifest(manifest), "duplicate scenario_id")

  manifest <- read_module_ci_manifest(validate = FALSE)
  manifest$scenarios[[1]]$fixture$path <- "fixtures/missing.tsv"
  expect_error(validate_module_ci_manifest(manifest), "fixture does not exist")
})
