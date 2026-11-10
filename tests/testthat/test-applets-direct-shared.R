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

test_that("RunApplet preserves current validation errors", {
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

})

# GUI audit R9 (P-15): what the applet cannot serve is refused by class before anything is
# launched. The applet is proteomics' (as before KIT-048: no other omic ever had one), and it
# works on the caller's `data_tbl` and `config_list`, as the workbooks define them.
test_that("the design applet refuses an omic it cannot serve and missing or run-less inputs", {
  source_dir <- withr::local_tempdir("design-applet-source-")
  output_dir <- withr::local_tempdir("design-applet-output-")
  runApplet <- function(omic_type = "proteomics") {
    .runDesignMatrixApplet(source_dir = source_dir, output_dir = output_dir, omic_type = omic_type, force = TRUE)
  }

  for (omic_type in c("metabolomics", "lipidomics", "transcriptomics", "integration")) {
    expect_error(runApplet(omic_type), class = "design_applet_omic_unsupported")
  }
  local({
    project_dirs <- list(metabolomics_exp = list(source_dir = source_dir, output_dir = output_dir))
    expect_error(RunApplet("designMatrix", "metabolomics", "exp", force = TRUE),
      class = "design_applet_omic_unsupported")
  })

  expect_error(runApplet(), class = "design_applet_input_missing")
  # The table the applet made from the file names before this repair has no runs.
  data_tbl <- data.frame(File = "S1.tsv", FullPath = file.path(source_dir, "S1.tsv"), SampleName = "S1")
  expect_error(runApplet(), class = "design_applet_input_invalid")
  data_tbl <- data.frame(Run = "S1")
  expect_error(runApplet(), "config_list", class = "design_applet_input_missing")
  expect_identical(list.files(source_dir), character())
})

# GUI audit R9 (P-15), replacing the KIT-103 W6 block that fed the builder a made-up table and
# pinned the owner-less save's refusal (F103-03). The applet as a workbook runs it: `RunApplet()`
# from an environment holding the report read by vroom and the configuration. Only the browser is
# replaced: `runApp` drives the applet's own app object (its UI and server) through testServer.
test_that("the design applet builds, renames and saves a design on its own temporary owner", {
  skip_if_not(workflowOwnerCapable())
  html <- htmltools::renderTags(mod_prot_design_builder_ui("designBuilder"))$html
  html <- gsub("data-tabsetid=\"[0-9]+\"", "data-tabsetid=\"<ID>\"", html)
  html <- gsub("tab-[0-9]+-([0-9]+)", "tab-<ID>-\\1", html)
  expect_identical(digest::digest(as.character(html), algo = "sha256", serialize = FALSE),
    "7961d7eb971fa7b26cb1c9a14ac6dfd5b3ec283b0bd6694594ddf89835790bbe")

  source_dir <- withr::local_tempdir("design-applet-source-")
  output_dir <- withr::local_tempdir("design-applet-output-")
  report_path <- file.path(withr::local_tempdir("design-applet-data-"), "report.tsv")
  utils::write.table(data.frame(
    Run = rep(c("S10", "S2", "S1"), each = 2L), Protein.Group = rep(c("P1", "P2"), 3L),
    Stripped.Sequence = rep(c("AAK", "CCK"), 3L), Precursor.Quantity = c(10, 20, 30, 40, 50, 60),
    Q.Value = 0.001), report_path, sep = "\t", row.names = FALSE, quote = FALSE)

  applet_seen <- new.env(parent = emptyenv())
  applet_seen$choices <- list()
  setBuilderInputs <- function(session, ...) {
    values <- list(...)
    names(values) <- paste0("designBuilder-", names(values))
    do.call(session$setInputs, values)
  }
  local_mocked_bindings(
    runApp = function(appDir, ...) {
      shiny::testServer(appDir, {
        session$flushReact()
        applet_seen$owner <- session$userData$dia_capacity_manager
        applet_seen$store <- applet_seen$owner$.context$store$directory
        applet_seen$store_existed <- dir.exists(applet_seen$store)
        setBuilderInputs(session, sample_to_rename = "S10", new_sample_name = "T10")
        setBuilderInputs(session, rename_sample = 1L)
        for (factor_name in c("Ctl", "Trt")) {
          setBuilderInputs(session, new_factor = factor_name)
          setBuilderInputs(session, add_factor = which(c("Ctl", "Trt") == factor_name))
        }
        # S2 is left unassigned: the saved design and its table drop it.
        assignments <- c(S1 = "Ctl", T10 = "Trt")
        for (index in seq_along(assignments)) {
          setBuilderInputs(session, selected_runs = names(assignments)[[index]],
            factor1_select = assignments[[index]], factor2_select = "", factor3_select = "", replicate_start = 1)
          setBuilderInputs(session, assign_metadata = index)
        }
        setBuilderInputs(session, contrast_group1 = "Trt", contrast_group2 = "Ctl")
        setBuilderInputs(session, add_contrast = 1L)
        setBuilderInputs(session, formula_string = "~ 0 + group")
        applet_seen$before_save <- applet_seen$stopped
        setBuilderInputs(session, save_results = 1L)
      })
      applet_seen$stopped
    },
    stopApp = function(returnValue = invisible()) {
      applet_seen$stopped <- returnValue
      invisible(NULL)
    },
    updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL, ...) {
      applet_seen$choices[[length(applet_seen$choices) + 1L]] <- list(id = inputId, choices = choices)
      invisible(NULL)
    },
    showNotification = function(ui, ...) invisible(NULL),
    .package = "shiny"
  )

  workbook <- local({
    project_dirs <- list(proteomics_exp = list(source_dir = source_dir, output_dir = output_dir))
    data_tbl <- vroom::vroom(report_path, show_col_types = FALSE)
    config_list <- list(deAnalysisParameters = list(formula_string = "~ 0 + 1"), kept = list(value = "1"))
    returned <- RunApplet("designMatrix", "proteomics", "exp", force = TRUE)
    environment()
  })
  returned <- workbook$returned

  # The applet ran on an owner of its own, and the fresh design is in mixedsort order.
  expect_s3_class(applet_seen$owner, "DIAStateAdapter")
  expect_true(applet_seen$store_existed)
  expect_identical(applet_seen$choices[[1L]], list(id = "sample_to_rename", choices = c("S1", "S2", "S10")))
  expect_null(applet_seen$before_save)

  # The save result: the renamed, assigned design, its table, the contrast and the configuration.
  expect_identical(names(returned), c("design_matrix", "data_cln", "contrasts_tbl", "config_list"))
  expect_identical(returned$design_matrix$Run, c("S1", "T10"))
  expect_identical(returned$design_matrix$group, c("Ctl", "Trt"))
  expect_identical(returned$design_matrix$replicates, c(1L, 1L))
  expect_s3_class(returned$data_cln, "data.frame")
  expect_identical(returned$data_cln$Run, c("T10", "T10", "S1", "S1"))
  expect_identical(returned$data_cln$Precursor.Quantity, c(10, 20, 50, 60))
  expect_identical(returned$contrasts_tbl$contrasts, "groupTrt-groupCtl")
  expect_identical(returned$config_list,
    list(deAnalysisParameters = list(formula_string = "~ 0 + group"), kept = list(value = "1")))

  # What the workbooks' later chunks read: the variables where the applet was run, and the files.
  for (name in names(returned)) expect_identical(workbook[[name]], returned[[name]])
  expect_setequal(list.files(source_dir), c("design_matrix.tab", "data_cln.tab", "contrast_strings.tab"))
  design_file <- utils::read.table(file.path(source_dir, "design_matrix.tab"), sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)
  expect_identical(design_file[c("Run", "group", "replicates")],
    data.frame(Run = c("S1", "T10"), group = c("Ctl", "Trt"), replicates = c(1L, 1L)))
  data_file <- utils::read.delim(file.path(source_dir, "data_cln.tab"), stringsAsFactors = FALSE)
  expect_equal(data_file, as.data.frame(returned$data_cln))
  expect_identical(readLines(file.path(source_dir, "contrast_strings.tab")), "groupTrt-groupCtl")
  expect_identical(readRDS(file.path(output_dir, "design_matrix_results.rds")), returned)

  # Nothing of the applet's store is left behind.
  expect_true(applet_seen$owner$.context$store$closed)
  expect_false(dir.exists(applet_seen$store))
  expect_identical(list.files(tempdir(), pattern = "^design-applet-store-"), character())
})
