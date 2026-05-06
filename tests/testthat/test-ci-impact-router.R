library(testthat)

test_that("ci impact router selects targeted lipid normalization lanes", {
  payload <- ci_impact_run("R/mod_lipid_norm_server.R")

  expect_identical(payload$impact_level, "targeted")
  expect_true("lipidomics/normalization" %in% ci_impact_modules(payload))
  expect_true("^e2e-lipidomics-canonical$" %in% ci_impact_e2e_filters(payload))
  expect_true(isTRUE(payload$run_cross_omic))
  expect_true("lipidomics-normalization" %in% ci_impact_rule_ids(payload))
})

test_that("ci impact router allows docs-only no-op", {
  payload <- ci_impact_run("docs/qa/module-ci.md")

  expect_identical(payload$impact_level, "none")
  expect_identical(payload$module_count, 0L)
  expect_identical(payload$e2e_count, 0L)
  expect_true("docs/qa/module-ci.md" %in% unlist(payload$no_op_files, use.names = FALSE))
  expect_true("docs-noop" %in% ci_impact_rule_ids(payload))
})

test_that("ci impact router fails closed for unknown R source", {
  payload <- ci_impact_run("R/unmapped_new_file.R")

  expect_identical(payload$impact_level, "broad")
  expect_true("proteomics/all" %in% ci_impact_modules(payload))
  expect_true("metabolomics/all" %in% ci_impact_modules(payload))
  expect_true("lipidomics/all" %in% ci_impact_modules(payload))
  expect_true("^e2e-proteomics-dia$" %in% ci_impact_e2e_filters(payload))
})

test_that("ci impact router does not let docs suppress code impact", {
  payload <- ci_impact_run(c("docs/qa/module-ci.md", "R/mod_lipid_da.R"))

  expect_true(payload$impact_level %in% c("targeted", "omic", "broad"))
  expect_true("lipidomics/differential_abundance" %in% ci_impact_modules(payload))
  expect_true("docs/qa/module-ci.md" %in% unlist(payload$no_op_files, use.names = FALSE))
})

test_that("ci impact router output is deterministic and deduplicated", {
  payload_a <- ci_impact_run(c("R/mod_lipid_da.R", "R/func_lipid_da.R"))
  payload_b <- ci_impact_run(c("R/func_lipid_da.R", "R/mod_lipid_da.R"))

  expect_identical(ci_impact_modules(payload_a), ci_impact_modules(payload_b))
  expect_identical(ci_impact_e2e_filters(payload_a), ci_impact_e2e_filters(payload_b))
  expect_identical(length(ci_impact_e2e_filters(payload_a)), length(unique(ci_impact_e2e_filters(payload_a))))
})
