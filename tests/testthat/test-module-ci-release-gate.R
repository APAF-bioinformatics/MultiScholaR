library(testthat)

module_ci_workflow_path <- function() {
  testthat::test_path("..", "..", ".github", "workflows", "module-ci.yml")
}

module_ci_workflow_text <- function() {
  paste(readLines(module_ci_workflow_path(), warn = FALSE), collapse = "\n")
}

expect_workflow_contains <- function(text, needle) {
  expect_true(grepl(needle, text, fixed = TRUE), info = sprintf("Missing workflow text: %s", needle))
}

test_that("MCI-025 push and pull-request jobs enforce all unit-contract module matrices", {
  workflow <- module_ci_workflow_text()

  expect_workflow_contains(workflow, "push:")
  expect_workflow_contains(workflow, "pull_request:")
  expect_workflow_contains(workflow, "workflow_dispatch:")
  expect_workflow_contains(workflow, "schedule:")
  expect_workflow_contains(workflow, "module-ci-unit:")
  expect_workflow_contains(workflow, "github.event_name == 'push'")
  expect_workflow_contains(workflow, "github.event_name == 'pull_request'")

  expected_unit_jobs <- list(
    foundation = c("name: foundation", "omic: all", "module: foundation"),
    cross_module = c("name: cross-module-integrity", "omic: all", "module: cross_module_integrity"),
    proteomics = c("name: proteomics", "omic: proteomics", "module: all"),
    metabolomics = c("name: metabolomics", "omic: metabolomics", "module: all"),
    lipidomics = c("name: lipidomics", "omic: lipidomics", "module: all")
  )

  for (job in expected_unit_jobs) {
    invisible(lapply(job, function(needle) expect_workflow_contains(workflow, needle)))
  }

  expect_workflow_contains(workflow, "Rscript tools/ci/run-module-ci.R")
  expect_workflow_contains(workflow, "--runtime ${{ matrix.runtime }}")
  expect_workflow_contains(workflow, "Rscript tools/ci/module-ci-scorecard.R")
  expect_workflow_contains(workflow, "$GITHUB_STEP_SUMMARY")
})

test_that("MCI-025 browser and representative E2E smoke cover every current omic", {
  workflow <- module_ci_workflow_text()

  expected_browser_jobs <- list(
    harness = c("name: harness", "filter: \"^(module-ci-harness|e2e-browser-harness)$\""),
    proteomics = c("name: proteomics-import-module", "filter: \"^module-ci-prot-import$\""),
    metabolomics = c("name: metabolomics-import-module", "filter: \"^module-ci-metab-import$\""),
    lipidomics = c("name: lipidomics-import-module", "filter: \"^module-ci-lipid-import$\"")
  )

  for (job in expected_browser_jobs) {
    invisible(lapply(job, function(needle) expect_workflow_contains(workflow, needle)))
  }

  expected_e2e_jobs <- list(
    proteomics = c("name: proteomics-dia", "filter: \"^e2e-proteomics-dia$\""),
    metabolomics = c("name: metabolomics-lc-gc", "filter: \"^e2e-metabolomics-lc-gc$\""),
    lipidomics = c("name: lipidomics-canonical", "filter: \"^e2e-lipidomics-canonical$\"")
  )

  for (job in expected_e2e_jobs) {
    invisible(lapply(job, function(needle) expect_workflow_contains(workflow, needle)))
  }

  expect_workflow_contains(workflow, "any::shinytest2")
  expect_workflow_contains(workflow, "any::chromote")
  expect_workflow_contains(workflow, "MULTISCHOLAR_E2E_ARTIFACT_DIR")
})

test_that("MCI-025 nightly and release gates run the full corpus before promotion", {
  workflow <- module_ci_workflow_text()

  expect_workflow_contains(workflow, "full-gates:")
  expect_workflow_contains(workflow, "github.event_name == 'schedule'")
  expect_workflow_contains(workflow, "startsWith(github.ref, 'refs/tags/v')")
  expect_workflow_contains(workflow, "startsWith(github.ref, 'refs/heads/release/')")
  expect_workflow_contains(workflow, "github.event.inputs.gate == 'nightly'")
  expect_workflow_contains(workflow, "github.event.inputs.gate == 'release'")

  expected_full_gates <- c(
    "name: all-module-unit-contracts",
    "name: all-current-e2e-workflows",
    "name: alternate-launch-browser-smoke",
    "name: report-export-fidelity",
    "name: cross-omic-coexistence",
    "filter = \"^e2e-\"",
    "filter = \"^(launch-test-mode|e2e-test-mode-integration|e2e-browser-harness)$\"",
    "filter = \"^(e2e-cross-omic|workflow-server-shared|module-ci-sentinels)$\""
  )
  invisible(lapply(expected_full_gates, function(needle) expect_workflow_contains(workflow, needle)))

  expect_workflow_contains(workflow, "release-candidate-gate:")
  expect_workflow_contains(workflow, "needs.full-gates.result")
  expect_workflow_contains(workflow, "final release tagging is blocked unless full gates pass")
  expect_workflow_contains(workflow, "test \"${{ needs.full-gates.result }}\" = \"success\"")
})

test_that("MCI-025 artifact upload paths exclude ignored local files and expose debug outputs", {
  workflow <- module_ci_workflow_text()

  forbidden_paths <- c("renv/", "renv.lock", ".Rprofile", ".ticket-config.json", "dev/", "Workbooks/")
  expect_false(any(vapply(forbidden_paths, grepl, logical(1), x = workflow, fixed = TRUE)))

  expected_artifact_paths <- c(
    "tests/testthat/_module_ci_artifacts/push/",
    "tests/testthat/_module_ci_artifacts/browser-smoke/",
    "tests/testthat/_module_ci_artifacts/full/",
    "tests/testthat/_module_ci_artifacts/release-gate/",
    "tests/testthat/_e2e_artifacts/browser-smoke/",
    "tests/testthat/_e2e_artifacts/representative/",
    "tests/testthat/_e2e_artifacts/full/"
  )
  invisible(lapply(expected_artifact_paths, function(needle) expect_workflow_contains(workflow, needle)))

  expect_workflow_contains(workflow, "actions/upload-artifact@v4")
  expect_workflow_contains(workflow, "retention-days: 14")
  expect_workflow_contains(workflow, "retention-days: 30")
})

test_that("MCI-025 scorecards retain scenario status and artifact paths", {
  artifact_dir <- tempfile("module-ci-scorecard-")
  dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(artifact_dir, "module-ci-run-manifest.json")

  jsonlite::write_json(
    list(
      schema_version = "1.0.0",
      generated_at = "2026-05-06T00:00:00Z",
      started_at = "2026-05-06T00:00:00Z",
      finished_at = "2026-05-06T00:00:01Z",
      invocation = list(omic = "all", module = "ci_release_gate", runtime = "unit-contract"),
      result = "passed",
      test_filter = "^(module-ci-release-gate)$",
      test_names = list("module-ci-release-gate"),
      artifact_dir = artifact_dir,
      scenarios = list(list(
        scenario_id = "MCI-025.1-ci-release-gate-enforcement",
        ticket_id = "MCI-025",
        item_id = "MCI-025.1",
        omic = "all",
        module = "ci_release_gate",
        runtime = "unit-contract",
        ci_lane = "release",
        result = "passed",
        artifacts = list("release-gate-summary.md"),
        run_artifacts = list(
          file.path(artifact_dir, "module-ci-run-manifest.json"),
          file.path(artifact_dir, "module-ci-scorecard.json"),
          file.path(artifact_dir, "module-ci-scorecard.md")
        )
      )),
      artifact_retention = list(
        scorecards = list("module-ci-scorecard.json", "module-ci-scorecard.md"),
        module_ci_artifacts = "tests/testthat/_module_ci_artifacts/**",
        e2e_artifacts = "tests/testthat/_e2e_artifacts/**"
      )
    ),
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(testthat::test_path("..", "..", "tools", "ci", "module-ci-scorecard.R"), "--artifact-dir", artifact_dir),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))

  scorecard <- jsonlite::read_json(file.path(artifact_dir, "module-ci-scorecard.json"), simplifyVector = FALSE)
  row <- scorecard$scenarios[[1]]
  expect_identical(row$scenario_id, "MCI-025.1-ci-release-gate-enforcement")
  expect_identical(row$omic, "all")
  expect_identical(row$module, "ci_release_gate")
  expect_identical(row$runtime, "unit-contract")
  expect_identical(row$status, "passed")
  expect_identical(row$result, "passed")
  expect_identical(row$ci_lane, "release")
  expect_true("release-gate-summary.md" %in% unlist(row$artifact_paths, use.names = FALSE))

  scorecard_md <- paste(readLines(file.path(artifact_dir, "module-ci-scorecard.md"), warn = FALSE), collapse = "\n")
  expect_match(scorecard_md, "Artifact paths", fixed = TRUE)
  expect_match(scorecard_md, "release-gate-summary.md", fixed = TRUE)
  expect_match(scorecard_md, "ci_release_gate", fixed = TRUE)
})
