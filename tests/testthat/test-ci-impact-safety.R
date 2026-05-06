library(testthat)

test_that("ci impact safety escalates shared state and CI self changes", {
  state_payload <- ci_impact_run("R/utils_workflow_state.R")
  expect_identical(state_payload$impact_level, "broad")
  expect_true(isTRUE(state_payload$run_cross_omic))
  expect_true("shared-state-and-workflow" %in% ci_impact_rule_ids(state_payload))

  ci_payload <- ci_impact_run(".github/workflows/module-ci.yml")
  expect_identical(ci_payload$impact_level, "broad")
  expect_true("all/ci_release_gate" %in% ci_impact_modules(ci_payload))
  expect_true("ci-self-protection" %in% ci_impact_rule_ids(ci_payload))
})

test_that("ci impact safety reports broad fallback for git diff errors", {
  payload <- ci_impact_run(
    character(),
    extra_args = c("--base-ref", "0000000000000000000000000000000000000000", "--head-ref", "HEAD")
  )

  expect_identical(payload$impact_level, "broad")
  expect_gt(length(payload$errors), 0L)
  expect_true("proteomics/all" %in% ci_impact_modules(payload))
})

test_that("ci impact safety output excludes ignored local path exposure", {
  payload <- ci_impact_run(".github/workflows/module-ci.yml")
  encoded <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null")

  forbidden <- c("dev/", "renv/", "renv.lock", ".Rprofile", ".ticket-config.json", "Workbooks/")
  expect_false(any(vapply(forbidden, grepl, logical(1), x = encoded, fixed = TRUE)))
})
