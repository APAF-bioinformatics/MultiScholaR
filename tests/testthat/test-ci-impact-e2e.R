library(testthat)

test_that("ci impact E2E metadata selects current omic workflow filters", {
  payload <- ci_impact_run("R/mod_lipid_norm_server.R")
  expect_true("^e2e-lipidomics-canonical$" %in% ci_impact_e2e_filters(payload))

  payload <- ci_impact_run("R/mod_metab_import_server.R")
  expect_true("^e2e-metabolomics-lc-gc$" %in% ci_impact_e2e_filters(payload))

  payload <- ci_impact_run("R/func_prot_limpa_qc_helpers.R")
  expect_true("^e2e-proteomics-dia-limpa$" %in% ci_impact_e2e_filters(payload))
})

test_that("run-e2e-ci accepts lane metadata and writes a run manifest", {
  artifact_dir <- tempfile("e2e-impact-run-")
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      testthat::test_path("..", "..", "tools", "ci", "run-e2e-ci.R"),
      "--lane", "lipid_canonical",
      "--reporter", "summary",
      "--artifact-dir", artifact_dir,
      "--dry-run", "true"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))

  manifest <- jsonlite::read_json(file.path(artifact_dir, "e2e-run-manifest.json"), simplifyVector = FALSE)
  expect_identical(manifest$result, "passed")
  expect_identical(manifest$lanes$lipid_canonical$test_filter, "^e2e-lipidomics-canonical$")
  expect_true("^e2e-lipidomics-canonical$" %in% unlist(manifest$filters, use.names = FALSE))
})
