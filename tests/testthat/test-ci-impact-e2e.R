library(testthat)

test_that("ci impact E2E metadata selects current omic workflow filters", {
  payload <- ci_impact_run("R/mod_lipid_norm_server.R")
  expect_true("^e2e-lipidomics-canonical$" %in% ci_impact_e2e_filters(payload))

  payload <- ci_impact_run("R/mod_metab_import_server.R")
  expect_true("^e2e-metabolomics-lc-gc$" %in% ci_impact_e2e_filters(payload))

  payload <- ci_impact_run("R/func_prot_limpa_qc_helpers.R")
  expect_true("^e2e-proteomics-dia-limpa$" %in% ci_impact_e2e_filters(payload))
})

test_that("ci impact marks only canonical DIA as browser-required", {
  proteomics <- ci_impact_e2e_rows(ci_impact_run("R/func_prot_qc_exclusion_helpers.R"))
  required <- Filter(function(row) isTRUE(row$browser_required), proteomics)
  expect_length(required, 1L)
  expect_identical(required[[1L]]$lane, "prot_dia")
  expect_identical(required[[1L]]$filter, "^e2e-proteomics-dia$")

  lipidomics <- ci_impact_e2e_rows(ci_impact_run("R/mod_lipid_norm_server.R"))
  expect_false(any(vapply(lipidomics, function(row) isTRUE(row$browser_required), logical(1L))))
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
  expect_false(manifest$browser$required)
  expect_true(is.logical(manifest$browser$preflight$available))
})

test_that("run-e2e-ci fails before tests when required Chromium is unavailable", {
  artifact_dir <- tempfile("e2e-required-browser-")
  missing_chrome <- file.path(artifact_dir, "missing-chromium")
  output <- suppressWarnings(
    system2(
      file.path(R.home("bin"), "Rscript"),
      c(
        testthat::test_path("..", "..", "tools", "ci", "run-e2e-ci.R"),
        "--lane", "prot_dia",
        "--reporter", "summary",
        "--artifact-dir", artifact_dir,
        "--browser-required", "true",
        "--dry-run", "true"
      ),
      stdout = TRUE,
      stderr = TRUE,
      env = c(sprintf("CHROMOTE_CHROME=%s", missing_chrome))
    )
  )
  expect_identical(attr(output, "status"), 1L, info = paste(output, collapse = "\n"))

  manifest <- jsonlite::read_json(
    file.path(artifact_dir, "e2e-run-manifest.json"),
    simplifyVector = FALSE
  )
  expect_identical(manifest$result, "failed")
  expect_true(manifest$browser$required)
  expect_false(manifest$browser$preflight$available)
  expect_match(manifest$failure_reason, "Required E2E browser preflight failed")
})
