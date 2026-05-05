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
      "MCI-009.1-proteomics-da-schema-smoke"
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
})

test_that("module CI manifest rejects duplicate scenario IDs and missing fixtures", {
  manifest <- read_module_ci_manifest(validate = FALSE)
  manifest$scenarios[[2]]$scenario_id <- manifest$scenarios[[1]]$scenario_id
  expect_error(validate_module_ci_manifest(manifest), "duplicate scenario_id")

  manifest <- read_module_ci_manifest(validate = FALSE)
  manifest$scenarios[[1]]$fixture$path <- "fixtures/missing.tsv"
  expect_error(validate_module_ci_manifest(manifest), "fixture does not exist")
})
