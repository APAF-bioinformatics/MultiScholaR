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

# KIT-103 W6: the design applet's builder (back-references #4 and #5) as .runDesignMatrixApplet
# calls it (func_general_applets.R:260-277): the proteomics builder UI, and its server with plain
# reactiveVals and no DIA owner. The fresh design is mixedsort; saving needs an owner since KIT-048
# (F103-03, parked P-15), so the save is the refusal only. Recorded before W6 (de2d891c).
test_that("the design applet's builder renders and starts a mixedsort fresh design without an owner", {
  html <- htmltools::renderTags(mod_prot_design_builder_ui("designBuilder"))$html
  html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", html)
  html <- gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", html)
  expect_identical(digest::digest(as.character(html), algo = "sha256", serialize = FALSE),
    "7961d7eb971fa7b26cb1c9a14ac6dfd5b3ec283b0bd6694594ddf89835790bbe")

  applet_seen <- new.env(parent = emptyenv())
  applet_seen$choices <- list()
  applet_seen$notes <- character()
  local_mocked_bindings(
    updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...) {
      applet_seen$choices[[length(applet_seen$choices) + 1L]] <- list(id = inputId, choices = choices)
      invisible(NULL)
    },
    showNotification = function(ui, ..., type = "default") {
      applet_seen$notes <- c(applet_seen$notes, paste0(type, ": ", ui))
      invisible(NULL)
    },
    .package = "shiny"
  )
  applet_table <- data.frame(Run = c("S10", "S2", "S1", "S2"), Intensity = c("1", "2", "3", "4"))
  shiny::testServer(mod_prot_design_builder_server, args = list(
    data_tbl = shiny::reactiveVal(applet_table),
    config_list = shiny::reactiveVal(list(deAnalysisParameters = list(formula_string = "~ 0 + group"))),
    column_mapping = shiny::reactiveVal(list(quantity_col = "Intensity"))
  ), {
    session$flushReact()
    session$setInputs(save_results = 1L)
    expect_null(session$returned())
  })
  # The initial state fills the sample selections in the fresh design's (mixedsort) order.
  expect_identical(applet_seen$choices[[1L]], list(id = "sample_to_rename", choices = c("S1", "S2", "S10")))
  expect_identical(applet_seen$notes,
    "warning: No samples have been assigned to groups. Please assign metadata before saving.")
})
