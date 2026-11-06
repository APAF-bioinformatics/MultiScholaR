# KIT-096 MCI-029: the transcriptomics module in the app, and the stable ids KIT-098's driver uses.

test_that("MCI-029.1 the app's dispatch finds the transcriptomics module, not the 'Module not yet implemented' fallback", {
  # app_server looks the module up by name in its own environment (mod_<omic>_ui / mod_<omic>_server).
  env <- environment(app_server)
  expect_true(exists("mod_transcriptomics_ui", envir = env))
  expect_true(exists("mod_transcriptomics_server", envir = env))
  ui <- get("mod_transcriptomics_ui", envir = env)
  html <- htmltools::renderTags(ui(paste0("transcriptomics", "_workflow")))$html
  expect_false(grepl("Module not yet implemented", html, fixed = TRUE))
  expect_match(html, "Transcriptomics Workflow", fixed = TRUE)
})

test_that("MCI-029.1 the stable ids: module, tabs, sub-modules and test ids", {
  html <- htmltools::renderTags(mod_transcriptomics_ui("transcriptomics_workflow"))$html
  expect_match(html, "id=\"transcriptomics_workflow-transcriptomics_tabs\"", fixed = TRUE)
  for (value in c("import", "design", "qc")) expect_match(html, sprintf("data-value=\"%s\"", value), fixed = TRUE)
  for (id in c("trx-workflow-tabs", "trx-tab-import", "trx-tab-design", "trx-tab-qc")) {
    expect_match(html, sprintf("data-testid=\"%s\"", id), fixed = TRUE)
  }
  expect_match(html, "transcriptomics_workflow-import-process_import", fixed = TRUE)
  expect_match(html, "transcriptomics_workflow-design-builder-save_results", fixed = TRUE)
  expect_match(html, "transcriptomics_workflow-qc-dynamic_qc_tabs", fixed = TRUE)
  expect_false(grepl("data-testid=\"lipid-", html, fixed = TRUE))
})

test_that("MCI-029.1 the module server holds the bus, the session's owner and the stepper's statuses", {
  skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir()
  project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project"), omic_types = "transcriptomics", label = NULL,
    force = TRUE, reuse_existing = TRUE))
  suppressMessages(shiny::testServer(mod_transcriptomics_server, args = list(id = "transcriptomics_workflow", project_dirs = project_dirs,
    omic_type = "transcriptomics", experiment_label = "kit096"), {
    owner <- shiny::isolate(workflow_data$state_manager)
    expect_true(inherits(owner, "DIAStateAdapter"))
    expect_identical(owner$.context$omic_type, "transcriptomics")
    expect_identical(session$userData$omics_workflow_data$transcriptomics, workflow_data)
    expect_identical(shiny::isolate(workflow_data$tab_status),
      list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled"))
    expect_match(paste(as.character(output$workflow_progress$html), collapse = ""), "Design")
    # The design's completion opens the QC step.
    status <- shiny::isolate(workflow_data$tab_status)
    status$design_matrix <- "complete"
    workflow_data$tab_status <- status
    session$flushReact()
    expect_identical(shiny::isolate(workflow_data$tab_status$quality_control), "pending")
    shiny::isolate(closeDIAStateAdapter(owner))
  }))
})
