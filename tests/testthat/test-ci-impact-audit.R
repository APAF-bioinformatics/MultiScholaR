library(testthat)

test_that("ci impact audit covers representative routing bundles", {
  cases <- list(
    lipid_norm = list(path = "R/mod_lipid_norm_server.R", impact = "targeted", module = "lipidomics/normalization", e2e = "^e2e-lipidomics-canonical$"),
    metab_da = list(path = "R/mod_metab_da_server.R", impact = "targeted", module = "metabolomics/differential_abundance", e2e = "^e2e-metabolomics-lc-gc$"),
    prot_limpa = list(path = "R/func_prot_limpa_qc_helpers.R", impact = "targeted", module = "proteomics/normalization", e2e = "^e2e-proteomics-dia-limpa$"),
    state = list(path = "R/utils_workflow_state.R", impact = "broad", module = "proteomics/all", e2e = "^e2e-proteomics-dia$"),
    ci = list(path = ".github/workflows/module-ci.yml", impact = "broad", module = "all/ci_release_gate", e2e = "^e2e-lipidomics-canonical$")
  )

  for (case in cases) {
    payload <- ci_impact_run(case$path)
    expect_identical(payload$impact_level, case$impact, info = case$path)
    expect_true(case$module %in% ci_impact_modules(payload), info = case$path)
    expect_true(case$e2e %in% ci_impact_e2e_filters(payload), info = case$path)
  }
})

test_that("ci impact audit treats docs-only and unknown source as opposite boundaries", {
  docs <- ci_impact_run("docs/qa/module-ci.md")
  expect_identical(docs$impact_level, "none")
  expect_identical(docs$module_count, 0L)

  unknown <- ci_impact_run("R/new_unmapped_surface.R")
  expect_identical(unknown$impact_level, "broad")
  expect_true("proteomics/all" %in% ci_impact_modules(unknown))
})
