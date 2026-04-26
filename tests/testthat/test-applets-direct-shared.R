# fidelity-coverage-compare: shared
library(testthat)

renderAppletUi <- function(ui) {
  htmltools::renderTags(ui)$html
}

test_that("workflow stepper preserves current status progression and markup", {
  steps <- list(
    list(name = "Import", key = "import", icon = "file-import"),
    list(name = "Design", key = "design", icon = "table"),
    list(name = "Summary", key = "summary", icon = "clipboard")
  )

  mixed_html <- renderAppletUi(
    render_workflow_stepper(
      steps = steps,
      tab_status = list(import = "complete", design = "pending", summary = "disabled")
    )
  )

  expect_match(mixed_html, "workflow-stepper", fixed = TRUE)
  expect_match(mixed_html, "step-completed", fixed = TRUE)
  expect_match(mixed_html, "step-current", fixed = TRUE)
  expect_match(mixed_html, "step-pending", fixed = TRUE)
  expect_match(mixed_html, "Import", fixed = TRUE)
  expect_match(mixed_html, "Design", fixed = TRUE)

  all_complete_html <- renderAppletUi(
    render_workflow_stepper(
      steps = steps,
      tab_status = list(import = "complete", design = "complete", summary = "complete")
    )
  )

  expect_false(grepl("step-current", all_complete_html, fixed = TRUE))
})

test_that("RunApplet preserves current validation and empty-source errors", {
  expect_error(
    RunApplet("bad-applet", "proteomics", "exp"),
    "Invalid applet_type specified",
    fixed = TRUE
  )

  expect_error(
    RunApplet("designMatrix", "bad-omic", "exp"),
    "Invalid omic_type specified",
    fixed = TRUE
  )

  expect_error(
    RunApplet("designMatrix", "proteomics", ""),
    "'experiment_label' must be a single non-empty character string",
    fixed = TRUE
  )

  local({
    expect_error(
      RunApplet(
        "designMatrix",
        "proteomics",
        "exp",
        project_dirs_object_name = "missing_project_dirs"
      ),
      "not found in the calling environment",
      fixed = TRUE
    )
  })

  local({
    broken_project_dirs <- "not-a-list"
    expect_error(
      RunApplet(
        "designMatrix",
        "proteomics",
        "exp",
        project_dirs_object_name = "broken_project_dirs"
      ),
      "is not a list as expected",
      fixed = TRUE
    )
  })

  local({
    missing_key_project_dirs <- list(
      proteomics_other = list(source_dir = tempdir(), output_dir = tempdir())
    )
    expect_error(
      RunApplet(
        "designMatrix",
        "proteomics",
        "exp",
        project_dirs_object_name = "missing_key_project_dirs"
      ),
      "No directory information found",
      fixed = TRUE
    )
  })

  local({
    missing_source_project_dirs <- list(proteomics_exp = list(output_dir = tempdir()))
    expect_error(
      RunApplet(
        "designMatrix",
        "proteomics",
        "exp",
        project_dirs_object_name = "missing_source_project_dirs"
      ),
      "does not contain 'source_dir'",
      fixed = TRUE
    )
  })

  source_dir <- tempfile("design-applet-source-")
  output_dir <- tempfile("design-applet-output-")
  dir.create(source_dir)
  dir.create(output_dir)

  expect_error(
    .runDesignMatrixApplet(
      source_dir = source_dir,
      output_dir = output_dir,
      omic_type = "proteomics",
      force = TRUE
    ),
    "No sample files",
    fixed = TRUE
  )
})
