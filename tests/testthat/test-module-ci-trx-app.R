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
  for (value in c("import", "design", "qc", "norm", "de", "summary")) expect_match(html, sprintf("data-value=\"%s\"", value), fixed = TRUE)
  for (id in c("trx-workflow-tabs", "trx-tab-import", "trx-tab-design", "trx-tab-qc", "trx-tab-norm", "trx-tab-de", "trx-tab-summary",
    "trx-norm-run", "trx-da-run", "trx-da-download", "trx-summary-ruviii-export", "trx-summary-generate-report")) {
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
      list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled", normalization = "disabled",
        differential_analysis = "disabled", session_summary = "disabled"))
    expect_match(paste(as.character(output$workflow_progress$html), collapse = ""), "Design")
    # The design's completion opens the QC step.
    status <- shiny::isolate(workflow_data$tab_status)
    status$design_matrix <- "complete"
    workflow_data$tab_status <- status
    session$flushReact()
    expect_identical(shiny::isolate(workflow_data$tab_status$quality_control), "pending")
    # KIT-097: each later tab opens when the one before completes.
    for (link in list(c("quality_control", "normalization"), c("normalization", "differential_analysis"), c("differential_analysis", "session_summary"))) {
      status <- shiny::isolate(workflow_data$tab_status)
      status[[link[[1L]]]] <- "complete"
      workflow_data$tab_status <- status
      session$flushReact()
      expect_identical(shiny::isolate(workflow_data$tab_status[[link[[2L]]]]), "pending", info = link[[2L]])
    }
    shiny::isolate(closeDIAStateAdapter(owner))
  }))
})

test_that("MCI-029.2 the module in process: import, design, QC, normalisation and DA each complete; the DA's completion keeps every status, the digest reads it, and the interactive volcano is built only on demand", {
  skip_if_not(workflowOwnerCapable())
  skip_if_not_installed("edgeR")
  root <- withr::local_tempdir()
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true")
  invisible(utils::capture.output(project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project"),
    omic_types = "transcriptomics", label = NULL, force = TRUE, reuse_existing = TRUE))))
  dir <- module_ci_trx_corpus()
  design <- module_ci_trx_corpus_design(dir)
  seen <- new.env(parent = emptyenv())
  seen$widgets <- 0L
  quiet <- function(...) invisible(NULL)
  # KIT-097 e2e defect: the browser hung after DA while the Glimma volcano of every gene was built and sent on completion.
  # The DA tab's interactive volcano is counted here; the browser sends the checkbox's declared default. KIT-106 W3 (D3):
  # the volcano is the binding's `trx_volcano_page`, built in process by its page action (the binding holds the widget
  # by value, so the action is counted).
  testthat::local_mocked_bindings(computeWorkflowOmicsVolcanoPageAction = function(manager, request, variant) {
    seen$widgets <- seen$widgets + 1L
    list(display = list(page = NULL))
  }, .package = "MultiScholaR")
  html <- htmltools::renderTags(mod_lipid_da_ui("de", profile = trxOmicProfile()))$html
  interactive_default <- grepl("<input id=\"de-volcano_interactive\"[^>]*checked", html)
  testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    if (identical(type, "error")) seen$errors <- c(seen$errors, paste(as.character(ui), collapse = " "))
    invisible(NULL)
  }, removeNotification = quiet, .package = "shiny")
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(mod_transcriptomics_server, args = list(id = "transcriptomics_workflow",
    project_dirs = project_dirs, omic_type = "transcriptomics", experiment_label = "kit097"), {
    set <- function(values) do.call(session$setInputs, values)
    # KIT-102 W5: the shared feature-table panel's preview and processing, then the frame's completion.
    trx_config <- omicGuiConfig("transcriptomics")
    preview <- omicImportFeatureTablePreview(file.path(dir, "counts.tsv"), trx_config)
    processOmicImportFeatureTable(workflow_data, list(assay1_name = "counts", assay1_data = preview$assay_data,
      assay1_file = file.path(dir, "counts.tsv"), vendor_format = "trx_counts", detected_format = preview$detected_format,
      id_col = preview$updates$id$selected, annotation_col = preview$updates$annotation$selected,
      experiment_paths = project_dirs$transcriptomics), trx_config,
      complete = function(...) completeOmicImport(workflow_data, trx_config, notify = FALSE),
      showNotification = quiet, removeNotification = quiet)
    session$flushReact()
    builder <- function(name) paste0("design-builder-", name)
    for (i in 1:2) set(stats::setNames(list(c("RPMI", "Sera")[[i]], i), builder(c("new_factor", "add_factor"))))
    for (i in 1:2) {
      group <- c("RPMI", "Sera")[[i]]
      set(stats::setNames(list(design$Run[design$group == group], group, "", "", 1L),
        builder(c("selected_runs", "factor1_select", "factor2_select", "factor3_select", "replicate_start"))))
      set(stats::setNames(list(i), builder("assign_metadata")))
    }
    set(stats::setNames(list("Sera", "RPMI", "~ 0 + group"), builder(c("contrast_group1", "contrast_group2", "formula_string"))))
    set(stats::setNames(list(1L, 1L), builder(c("add_contrast", "save_results"))))
    set(list(`qc-filter-apply_filter` = 1L))
    set(list(`qc-s4_finalize-finalize_qc` = 1L))
    session$flushReact()
    seen$after_qc <- shiny::isolate(workflow_data$tab_status)
    set(list(`norm-run_normalization` = 1L))
    session$flushReact()
    seen$after_norm <- shiny::isolate(workflow_data$tab_status)
    set(list(`de-volcano_interactive` = interactive_default))
    set(list(`de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 1, `de-run_da_analysis` = 1L))
    session$flushReact()
    seen$after_da <- shiny::isolate(workflow_data$tab_status)
    # The selectors as the browser holds them after the finaliser's updates.
    set(list(`de-volcano_contrast` = "Sera_vs_RPMI", `de-volcano_assay` = "counts"))
    try(output[["de-volcano_glimma"]], silent = TRUE)
    seen$widgets_on_completion <- seen$widgets
    seen$static <- tryCatch({ output[["de-volcano_static"]]; "drawn" }, error = function(e) conditionMessage(e))
    set(list(`de-volcano_interactive` = TRUE))
    try(output[["de-volcano_glimma"]], silent = TRUE)
    seen$widgets_on_demand <- seen$widgets
    seen$head <- shiny::isolate(workflow_data$state_manager$current_state)
    seen$digest <- tryCatch(jsonlite::fromJSON(jsonlite::toJSON(collect_state_digest(workflow_states = list(transcriptomics = workflow_data)),
      auto_unbox = TRUE, null = "null"), simplifyVector = FALSE), error = function(e) conditionMessage(e))
    shiny::isolate(closeDIAStateAdapter(workflow_data$state_manager))
  }))))
  expect_null(seen$errors)
  expect_identical(seen$after_qc$quality_control, "complete")
  expect_identical(seen$after_norm$normalization, "complete")
  expect_identical(seen$after_norm$differential_analysis, "pending")
  expect_identical(seen$head, "trx_ruvs")
  expect_identical(names(seen$after_da), c("setup_import", "design_matrix", "quality_control", "normalization", "differential_analysis", "session_summary"))
  expect_identical(seen$after_da$differential_analysis, "complete")
  expect_identical(seen$after_da$session_summary, "pending")
  expect_identical(seen$digest$step_status_per_omic$transcriptomics$differential_analysis, "complete")
  # The interactive volcano is not built when the DA completes; the static one is, and the interactive one on demand.
  expect_identical(seen$widgets_on_completion, 0L)
  expect_identical(seen$static, "drawn")
  expect_identical(seen$widgets_on_demand, 1L)
})
