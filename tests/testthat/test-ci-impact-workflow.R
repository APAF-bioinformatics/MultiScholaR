library(testthat)

test_that("ci impact workflow exposes detect-impact and dynamic matrix jobs", {
  workflow <- ci_impact_workflow_text()

  expect_match(workflow, "detect-impact:", fixed = TRUE)
  expect_match(workflow, "Rscript tools/ci/detect-impact.R", fixed = TRUE)
  expect_match(workflow, "module-ci-impacted:", fixed = TRUE)
  expect_match(workflow, "matrix: ${{ fromJSON(needs.detect-impact.outputs.module_matrix) }}", fixed = TRUE)
  expect_match(workflow, "module-ci-browser-impacted:", fixed = TRUE)
  expect_match(workflow, "e2e-impacted:", fixed = TRUE)
  expect_match(workflow, "Rscript tools/ci/run-e2e-ci.R", fixed = TRUE)
  expect_match(workflow, "MULTISCHOLAR_E2E_BROWSER_REQUIRED", fixed = TRUE)
  expect_match(workflow, "matrix.browser_required", fixed = TRUE)
  expect_match(workflow, "test-output.log", fixed = TRUE)
  expect_match(workflow, "impact-routing-summary:", fixed = TRUE)
})

test_that("ci impact workflow keeps full nightly and release gates conservative", {
  workflow <- ci_impact_workflow_text()

  expect_match(workflow, "full-gates:", fixed = TRUE)
  expect_match(workflow, "filter = \"^e2e-\"", fixed = TRUE)
  expect_match(workflow, "peptide-qc-browser-required", fixed = TRUE)
  expect_match(workflow, "Run mandatory DIA peptide-QC browser evidence", fixed = TRUE)
  expect_match(workflow, "release-candidate-gate:", fixed = TRUE)
  expect_match(workflow, "test \"${{ needs.full-gates.result }}\" = \"success\"", fixed = TRUE)
})

test_that("ci impact workflow uploads only allowlisted artifact roots", {
  workflow <- ci_impact_workflow_text()

  expected <- c(
    "tests/testthat/_module_ci_artifacts/impact/",
    "tests/testthat/_module_ci_artifacts/full/",
    "tests/testthat/_module_ci_artifacts/release-gate/",
    "tests/testthat/_e2e_artifacts/impact/",
    "tests/testthat/_e2e_artifacts/full/"
  )
  invisible(lapply(expected, function(path) expect_match(workflow, path, fixed = TRUE)))

  forbidden <- c("dev/", "renv/", "renv.lock", ".Rprofile", ".ticket-config.json", "Workbooks/")
  expect_false(any(vapply(forbidden, grepl, logical(1), x = workflow, fixed = TRUE)))
})
