library(testthat)

test_that("ci impact scorecard summary contains routing explanation", {
  artifact_dir <- tempfile("impact-scorecard-")
  dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(artifact_dir, "impact.json")
  summary <- file.path(artifact_dir, "impact.md")
  payload <- ci_impact_run("R/mod_lipid_da.R", output = output, summary = summary)

  expect_true(file.exists(output))
  expect_true(file.exists(summary))
  expect_identical(payload$impact_level, "targeted")
  expect_true("lipidomics/differential_abundance" %in% ci_impact_modules(payload))

  summary_text <- paste(readLines(summary, warn = FALSE), collapse = "\n")
  expect_match(summary_text, "CI Impact Routing Summary", fixed = TRUE)
  expect_match(summary_text, "Selected Module CI", fixed = TRUE)
  expect_match(summary_text, "Selected E2E", fixed = TRUE)
  expect_match(summary_text, "lipidomics-da", fixed = TRUE)
})

test_that("ci impact scorecard includes stable machine fields", {
  payload <- ci_impact_run("R/mod_lipid_da.R")
  expected_fields <- c(
    "schema_version", "resolver_version", "generated_at", "map_schema_version",
    "impact_level", "changed_files", "matched_rules", "module_matrix",
    "browser_matrix", "e2e_matrix", "run_cross_omic", "run_report_export",
    "impact_summary", "errors"
  )

  expect_true(all(expected_fields %in% names(payload)))
})
