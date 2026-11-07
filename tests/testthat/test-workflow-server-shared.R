# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeWorkflowStateStub <- function(default_type = "DIA") {
  list(
    new = function(audit_enabled = TRUE) {
      list(workflow_type = default_type, audit_enabled = audit_enabled)
    }
  )
}

makeRecordingServerStub <- function(captured, key) {
  force(captured)
  force(key)

  function(id, ...) {
    captured$calls[[key]] <- list(id = id, args = list(...))
    invisible(NULL)
  }
}

makeWorkflowStepperStub <- function(steps, tab_status) {
  shiny::tags$div(
    class = "workflow-stepper-stub",
    paste(vapply(steps, `[[`, character(1), "name"), collapse = "|"),
    paste(names(tab_status), collapse = "|")
  )
}

renderSharedOutput <- function(value) {
  paste(capture.output(print(value)), collapse = "\n")
}

expectWorkflowGatesAreMonotonic <- function(workflow_data, session, differential_key) {
  reset_status <- workflow_data$tab_status
  reset_status$design_matrix <- "disabled"
  reset_status$quality_control <- "complete"
  reset_status$normalization <- "pending"
  reset_status[[differential_key]] <- "disabled"
  reset_status$session_summary <- "disabled"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$design_matrix <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$quality_control, "complete")

  reset_status <- workflow_data$tab_status
  reset_status$quality_control <- "pending"
  reset_status$normalization <- "complete"
  reset_status[[differential_key]] <- "pending"
  reset_status$session_summary <- "disabled"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$quality_control <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$normalization, "complete")

  reset_status <- workflow_data$tab_status
  reset_status$quality_control <- "complete"
  reset_status$normalization <- "pending"
  reset_status[[differential_key]] <- "complete"
  reset_status$session_summary <- "pending"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs$normalization <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status[[differential_key]], "complete")

  reset_status <- workflow_data$tab_status
  reset_status$normalization <- "complete"
  reset_status[[differential_key]] <- "pending"
  reset_status$session_summary <- "complete"
  workflow_data$tab_status <- reset_status
  session$flushReact()

  updated_tabs <- workflow_data$tab_status
  updated_tabs[[differential_key]] <- "complete"
  workflow_data$tab_status <- updated_tabs
  session$flushReact()
  expect_equal(workflow_data$tab_status$session_summary, "complete")
}

test_that("proteomics workflow server wires modules and advances workflow state", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  fixture_dir <- tempfile("prot-workflow-server-")
  source_dir <- file.path(fixture_dir, "source")
  results_dir <- file.path(fixture_dir, "results")
  dir.create(source_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE)
  saveRDS(
    data.frame(sequence = c("PEPTIDE_1", "PEPTIDE_2"), stringsAsFactors = FALSE),
    file.path(source_dir, "aa_seq_tbl_final.RDS")
  )

  had_global_aa <- exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE)
  old_global_aa <- if (had_global_aa) get("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE) else NULL
  withr::defer({
    if (had_global_aa) {
      assign("aa_seq_tbl_final", old_global_aa, envir = .GlobalEnv)
    } else if (exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE)) {
      rm(list = "aa_seq_tbl_final", envir = .GlobalEnv)
    }
    unlink(fixture_dir, recursive = TRUE, force = TRUE)
  })

  localNamespaceBindings(
    package_ns,
    list(
      createDIASessionOwner = function(workflow_data, experiment_paths, session, ...) makeWorkflowStateStub("DIA")$new(FALSE), # KIT-048: the owner is created at session start
      mod_prot_import_server = makeRecordingServerStub(captured, "import"),
      mod_prot_design_server = makeRecordingServerStub(captured, "design"),
      mod_prot_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_prot_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_prot_da_server = makeRecordingServerStub(captured, "da"),
      mod_prot_enrich_server = makeRecordingServerStub(captured, "enrich"),
      mod_prot_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_proteomics_server,
    args = list(
      project_dirs = list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = source_dir,
          results_dir = results_dir
        )
      ),
      omic_type = "proteomics",
      experiment_label = "workflow-demo"
    ),
    {
      session$setInputs(workflow_tabs = "normalization")
      session$flushReact()

      expect_equal(captured$calls$import$id, "setup_import")
      expect_equal(captured$calls$design$id, "design_matrix")
      expect_equal(captured$calls$qc$id, "quality_control")
      expect_false(workflow_data$state_manager$audit_enabled)
      expect_equal(captured$calls$norm$id, "normalization")
      expect_equal(captured$calls$da$id, "differential_expression")
      expect_equal(captured$calls$enrich$id, "enrichment_analysis")
      expect_equal(captured$calls$summary$id, "session_summary")
      expect_equal(captured$calls$summary$args[[1]], list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = source_dir,
          results_dir = results_dir
        )
      ))
      expect_true(is.function(captured$calls$norm$args[[5]]))
      expect_equal(captured$calls$norm$args[[5]](), "normalization")
      # KIT-013 moved the proteomics globals onto the workflow bus (KIT-034, quarantine Q-15).
      expect_false(exists("aa_seq_tbl_final", envir = .GlobalEnv, inherits = FALSE))
      expect_false(is.null(shiny::isolate(workflow_data$aa_seq_tbl_final)))

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Import", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_expression, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_expression <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$enrichment_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$enrichment_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")
    }
  )
})

test_that("proteomics workflow server bypasses QC for TMT workflows", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  fixture_dir <- tempfile("prot-workflow-tmt-")
  dir.create(file.path(fixture_dir, "source"), recursive = TRUE)
  withr::defer(unlink(fixture_dir, recursive = TRUE, force = TRUE))

  localNamespaceBindings(
    package_ns,
    list(
      createDIASessionOwner = function(workflow_data, experiment_paths, session, ...) makeWorkflowStateStub("TMT")$new(FALSE), # KIT-048: the owner is created at session start
      mod_prot_import_server = makeRecordingServerStub(captured, "import"),
      mod_prot_design_server = makeRecordingServerStub(captured, "design"),
      mod_prot_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_prot_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_prot_da_server = makeRecordingServerStub(captured, "da"),
      mod_prot_enrich_server = makeRecordingServerStub(captured, "enrich"),
      mod_prot_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_proteomics_server,
    args = list(
      project_dirs = list(
        proteomics = list(
          base_dir = fixture_dir,
          source_dir = file.path(fixture_dir, "source"),
          results_dir = file.path(fixture_dir, "results")
        )
      ),
      omic_type = "proteomics",
      experiment_label = "workflow-demo"
    ),
    {
      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()

      expect_equal(workflow_data$tab_status$quality_control, "complete")
      expect_equal(workflow_data$tab_status$normalization, "pending")
    }
  )
})

test_that("lipidomics workflow server wires modules, status progression, and launcher shell", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      createWorkflowSessionOwner = function(workflow_data, experiment_paths, session, ...) makeWorkflowStateStub("DIA")$new(FALSE), # KIT-048: the owner is created at session start
      mod_lipid_import_server = makeRecordingServerStub(captured, "import"),
      mod_lipid_design_server = makeRecordingServerStub(captured, "design"),
      mod_lipid_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_lipid_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_lipid_da_server = makeRecordingServerStub(captured, "da"),
      mod_lipid_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_lipidomics_server,
    args = list(
      project_dirs = list(
        lipidomics = list(
          base_dir = tempfile("lipid-workflow-base-"),
          source_dir = tempfile("lipid-workflow-source-"),
          results_dir = tempfile("lipid-workflow-results-")
        )
      ),
      omic_type = "lipidomics",
      experiment_label = "lipid-demo"
    ),
    {
      session$setInputs(lipidomics_tabs = "norm")
      session$flushReact()

      expect_equal(captured$calls$import$id, "import")
      expect_equal(captured$calls$design$id, "design")
      expect_equal(captured$calls$qc$id, "qc")
      expect_equal(captured$calls$norm$id, "norm")
      expect_equal(captured$calls$da$id, "de")
      expect_equal(captured$calls$summary$id, "summary")
      expect_equal(captured$calls$norm$args$selected_tab(), "norm")

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Normalize", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")

      expectWorkflowGatesAreMonotonic(
        workflow_data = workflow_data,
        session = session,
        differential_key = "differential_analysis"
      )
    }
  )

  launcher_capture <- new.env(parent = emptyenv())
  launcher_capture$server_calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      mod_lipidomics_ui = function(id) {
        launcher_capture$ui_id <- id
        shiny::div("lipidomics launcher ui")
      },
      mod_lipidomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        launcher_capture$server_calls[[length(launcher_capture$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        invisible(NULL)
      }
    )
  )

  local_mocked_bindings(
    shinyApp = function(ui, server, ...) {
      list(ui = ui, server = server, dots = list(...))
    },
    .package = "shiny"
  )

  launcher_dir <- tempfile("lipidomics-launcher-")
  withr::defer(unlink(launcher_dir, recursive = TRUE, force = TRUE))
  app <- run_lipidomics_app(base_dir = launcher_dir)

  expect_equal(launcher_capture$ui_id, "lipidomics_app")
  expect_true(dir.exists(file.path(launcher_dir, "data")))
  expect_true(dir.exists(file.path(launcher_dir, "results")))
  expect_true(dir.exists(file.path(launcher_dir, "scripts")))

  app$server(list(), new.env(parent = emptyenv()), list())

  expect_equal(launcher_capture$server_calls[[1]]$id, "lipidomics_app")
  expect_equal(launcher_capture$server_calls[[1]]$omic_type, "lipidomics")
  expect_equal(launcher_capture$server_calls[[1]]$experiment_label, "Lipidomics Test")
})

test_that("metabolomics workflow server wires modules, status progression, and launcher shell", {
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      createWorkflowSessionOwner = function(workflow_data, experiment_paths, session, ...) makeWorkflowStateStub("DIA")$new(FALSE), # KIT-048: the owner is created at session start
      mod_metab_import_server = makeRecordingServerStub(captured, "import"),
      mod_metab_design_server = makeRecordingServerStub(captured, "design"),
      mod_metab_qc_server = makeRecordingServerStub(captured, "qc"),
      mod_metab_norm_server = makeRecordingServerStub(captured, "norm"),
      mod_metab_da_server = makeRecordingServerStub(captured, "da"),
      mod_metab_summary_server = makeRecordingServerStub(captured, "summary"),
      render_workflow_stepper = makeWorkflowStepperStub
    )
  )

  local_mocked_bindings(
    log_info = function(...) invisible(NULL),
    log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL),
    .package = "logger"
  )

  testServer(
    mod_metabolomics_server,
    args = list(
      project_dirs = list(
        metabolomics = list(
          base_dir = tempfile("metab-workflow-base-"),
          source_dir = tempfile("metab-workflow-source-"),
          results_dir = tempfile("metab-workflow-results-")
        )
      ),
      omic_type = "metabolomics",
      experiment_label = "metab-demo"
    ),
    {
      session$setInputs(metabolomics_tabs = "norm")
      session$flushReact()

      expect_equal(captured$calls$import$id, "import")
      expect_equal(captured$calls$design$id, "design")
      expect_equal(captured$calls$qc$id, "qc")
      expect_equal(captured$calls$norm$id, "norm")
      expect_equal(captured$calls$da$id, "de")
      expect_equal(captured$calls$summary$id, "summary")
      expect_equal(captured$calls$norm$args$selected_tab(), "norm")

      progress_markup <- renderSharedOutput(output$workflow_progress)
      expect_match(progress_markup, "Normalize", fixed = TRUE)
      expect_match(progress_markup, "Summary", fixed = TRUE)

      updated_tabs <- workflow_data$tab_status
      updated_tabs$design_matrix <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$quality_control, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$quality_control <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$normalization, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$normalization <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$differential_analysis, "pending")

      updated_tabs <- workflow_data$tab_status
      updated_tabs$differential_analysis <- "complete"
      workflow_data$tab_status <- updated_tabs
      session$flushReact()
      expect_equal(workflow_data$tab_status$session_summary, "pending")

      expectWorkflowGatesAreMonotonic(
        workflow_data = workflow_data,
        session = session,
        differential_key = "differential_analysis"
      )
    }
  )

  launcher_capture <- new.env(parent = emptyenv())
  launcher_capture$server_calls <- list()

  localNamespaceBindings(
    package_ns,
    list(
      mod_metabolomics_ui = function(id) {
        launcher_capture$ui_id <- id
        shiny::div("metabolomics launcher ui")
      },
      mod_metabolomics_server = function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
        launcher_capture$server_calls[[length(launcher_capture$server_calls) + 1L]] <- list(
          id = id,
          project_dirs = project_dirs,
          omic_type = omic_type,
          experiment_label = experiment_label,
          volumes = volumes
        )
        invisible(NULL)
      }
    )
  )

  local_mocked_bindings(
    shinyApp = function(ui, server, ...) {
      list(ui = ui, server = server, dots = list(...))
    },
    .package = "shiny"
  )

  launcher_dir <- tempfile("metabolomics-launcher-")
  withr::defer(unlink(launcher_dir, recursive = TRUE, force = TRUE))
  app <- run_metabolomics_app(base_dir = launcher_dir)

  expect_equal(launcher_capture$ui_id, "metabolomics_app")
  expect_true(dir.exists(file.path(launcher_dir, "data")))
  expect_true(dir.exists(file.path(launcher_dir, "results")))
  expect_true(dir.exists(file.path(launcher_dir, "scripts")))

  app$server(list(), new.env(parent = emptyenv()), list())

  expect_equal(launcher_capture$server_calls[[1]]$id, "metabolomics_app")
  expect_equal(launcher_capture$server_calls[[1]]$omic_type, "metabolomics")
  expect_equal(launcher_capture$server_calls[[1]]$experiment_label, "Metabolomics Test")
})

# KIT-102 W4: the shared shell (`mod_omic_workflow_server`) and the import frame (`mod_omic_import_server`).

# Each omic's tab servers, in today's call order: the stub name's tab and the module id.
workflowTabServers <- function(omic) {
  switch(omic,
    proteomics = c(import = "setup_import", design = "design_matrix", qc = "quality_control", norm = "normalization",
      da = "differential_expression", enrich = "enrichment_analysis", summary = "session_summary"),
    c(import = "import", design = "design", qc = "qc", norm = "norm", da = "de", summary = "summary"))
}

# Stub every tab server of an omic (the recording stub, with the call order), the owner and the stepper.
localWorkflowServerStubs <- function(omic, captured, .local_envir = parent.frame()) {
  module <- c(proteomics = "prot", lipidomics = "lipid", metabolomics = "metab", transcriptomics = "lipid")[[omic]]
  captured$calls <- list()
  captured$order <- character()
  stubs <- lapply(names(workflowTabServers(omic)), function(tab) {
    force(tab)
    function(id, ...) {
      captured$calls[[tab]] <- list(id = id, args = list(...))
      captured$order <- c(captured$order, tab)
      invisible(NULL)
    }
  })
  names(stubs) <- paste0("mod_", module, "_", names(workflowTabServers(omic)), "_server")
  owner <- function(workflow_data, experiment_paths, session, ...) {
    captured$owner_args <- list(...)
    makeWorkflowStateStub("DIA")$new(FALSE)
  }
  localNamespaceBindings(asNamespace("MultiScholaR"), c(stubs, list(createWorkflowSessionOwner = owner,
    createDIASessionOwner = owner, render_workflow_stepper = makeWorkflowStepperStub)), .local_envir = .local_envir)
  testthat::local_mocked_bindings(log_info = function(...) invisible(NULL), log_warn = function(...) invisible(NULL),
    log_error = function(...) invisible(NULL), .package = "logger", .env = .local_envir)
}

workflowProjectDirs <- function(omic) {
  stats::setNames(list(list(base_dir = tempfile("w4-base-"), source_dir = tempfile("w4-source-"),
    results_dir = tempfile("w4-results-"))), omic)
}

argNames <- function(args) if (is.null(names(args))) rep("", length(args)) else names(args)

test_that("KIT-102 W4: each shell calls its tab servers in today's order with today's exact arguments", {
  # The argument names by position, as each pre-KIT-102 shell passed them (mod_lipidomics.R:207-263,
  # mod_proteomics.R:328-374, mod_transcriptomics.R:258-264); "" is a positional argument.
  feature <- list(import = c("", "", "volumes"), design = c("", "", "volumes", "qc_trigger"), qc = c("", "", "", "", ""),
    norm = c("", "", "", "", "selected_tab"), da = c("", "", "", ""), summary = c("", "", "", ""))
  shapes <- list(lipidomics = feature, metabolomics = feature,
    transcriptomics = c(lapply(feature[c("import", "design", "qc")], c, "profile"), list(norm = c("", "", "", "", "profile")),
      lapply(feature[c("da", "summary")], c, "profile")),
    proteomics = list(import = c("", "", "volumes"), design = c("", "", "volumes", "qc_trigger"),
      qc = c("", "", "", "", "qc_trigger"), norm = c("", "", "", "", ""), da = c("", "", "", "", ""),
      enrich = c("", "", "", "", ""), summary = c("", "", "", "")))
  for (omic in names(shapes)) {
    captured <- new.env(parent = emptyenv())
    localWorkflowServerStubs(omic, captured)
    project_dirs <- workflowProjectDirs(omic)
    server <- get(paste0("mod_", omic, "_server"), envir = asNamespace("MultiScholaR"))
    testServer(server, args = list(project_dirs = project_dirs, omic_type = omic, experiment_label = "w4-demo"), {
      tabs <- workflowTabServers(omic)
      # testServer masks a caller's expression with the module frame and its parent: the shell holds no name the
      # pre-KIT-102 shells did not (a new local shadowed the entry-point suite's `steps`).
      before <- c("input", "output", "session", "ns", "workflow_data", "paths_key", "experiment_paths", "qc_trigger", "selected_tab", "update_tab_status",
        "checkpoint_status", "checkpoint_volumes", "aa_seq_file_path", "aa_seq_tbl_final", "profile", "chain", "link")
      expect_identical(setdiff(ls(session$env, all.names = TRUE), before), character(), info = omic)
      expect_identical(setdiff(ls(parent.env(session$env), all.names = TRUE),
        c("id", "project_dirs", "omic_type", "experiment_label", "volumes", "config")), character(), info = omic)
      expect_identical(captured$order, names(tabs), info = omic)
      expect_identical(vapply(captured$calls[names(tabs)], `[[`, character(1), "id"), tabs, info = omic)
      expect_identical(lapply(captured$calls[names(tabs)], function(call) argNames(call$args)), shapes[[omic]], info = omic)
      for (tab in setdiff(names(tabs), "summary")) {
        expect_identical(captured$calls[[tab]]$args[[1L]], workflow_data, info = paste(omic, tab))
        expect_identical(captured$calls[[tab]]$args[[2L]], project_dirs[[omic]], info = paste(omic, tab))
      }
      summary <- captured$calls$summary$args
      expect_identical(unname(summary[1:4]), list(project_dirs, omic, "w4-demo", workflow_data), info = omic)
      expect_null(captured$calls$import$args$volumes)
      expect_identical(captured$calls$qc$args[[3L]], omic)
      expect_identical(captured$calls$qc$args[[4L]], "w4-demo")
      qc_trigger <- captured$calls$qc$args$qc_trigger %||% captured$calls$qc$args[[5L]]
      expect_identical(captured$calls$design$args$qc_trigger, qc_trigger, info = omic)
      expect_s3_class(qc_trigger, "reactiveVal")
      if (identical(omic, "transcriptomics")) {
        for (tab in names(tabs)) expect_identical(captured$calls[[tab]]$args$profile, trxOmicProfile(), info = tab)
      }
      selected <- switch(omic, proteomics = captured$calls$norm$args[[5L]], transcriptomics = NULL,
        captured$calls$norm$args$selected_tab)
      do.call(session$setInputs, stats::setNames(list("qc"), omicGuiConfig(omic)$shell$tabset_id))
      if (!is.null(selected)) expect_identical(selected(), "qc", info = omic)
      expect_identical(workflow_data$active_tab(), "qc", info = omic)
      if (identical(omic, "proteomics")) {
        expect_identical(captured$calls$da$args[[5L]], selected)
        expect_identical(captured$calls$enrich$args[[5L]], selected)
        # The owner is the DIA alias with its default omic; the native-restore QC trigger is kept.
        expect_identical(captured$owner_args, list())
        expect_null(qc_trigger())
        workflow_data$native_session_restored <- Sys.time()
        session$flushReact()
        expect_true(qc_trigger())
      } else {
        expect_identical(captured$owner_args, list(omic_type = omic), info = omic)
      }
    })
  }
})

test_that("KIT-102 W4: the shared shell builds the bus, the statuses, the registration and the chain from the configuration", {
  captured <- new.env(parent = emptyenv())
  localWorkflowServerStubs("lipidomics", captured)
  config <- omicGuiConfig("lipidomics")
  # Three tabs only: the chain is derived from the configured tabs.
  config$shell$tabs <- config$shell$tabs[1:3]
  testServer(mod_omic_workflow_server, args = list(project_dirs = workflowProjectDirs("lipidomics"),
    omic_type = "lipidomics", experiment_label = "w4-demo", config = config), {
    expect_identical(captured$order, c("import", "design", "qc"))
    expect_identical(shiny::isolate(workflow_data$tab_status),
      list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled"))
    expect_identical(shiny::isolate(workflow_data$data_type), "lipid")
    expect_true("lipid_assay_data" %in% names(shiny::isolate(shiny::reactiveValuesToList(workflow_data))))
    expect_identical(session$userData$omics_workflow_data$lipidomics, workflow_data)
    expect_null(session$userData$workflow_data)
    markup <- renderSharedOutput(output$workflow_progress)
    expect_match(markup, "Import|Design|QC", fixed = TRUE)
    status <- workflow_data$tab_status
    status$design_matrix <- "complete"
    workflow_data$tab_status <- status
    session$flushReact()
    expect_identical(shiny::isolate(workflow_data$tab_status$quality_control), "pending")
    # A later tab is not reopened.
    status <- workflow_data$tab_status
    status$quality_control <- "complete"
    workflow_data$tab_status <- status
    session$flushReact()
    expect_identical(names(shiny::isolate(workflow_data$tab_status)), c("setup_import", "design_matrix", "quality_control"))
  })

  captured <- new.env(parent = emptyenv())
  localWorkflowServerStubs("proteomics", captured)
  testServer(mod_omic_workflow_server, args = list(project_dirs = workflowProjectDirs("proteomics"),
    omic_type = "proteomics", experiment_label = "w4-demo", config = omicGuiConfig("proteomics")), {
    expect_identical(session$userData$workflow_data, workflow_data)
    expect_null(session$userData$omics_workflow_data)
    expect_identical(shiny::isolate(workflow_data$tab_status), list(setup_import = "pending", design_matrix = "disabled",
      quality_control = "disabled", normalization = "disabled", differential_expression = "disabled",
      enrichment_analysis = "disabled", session_summary = "disabled"))
  })

  # Missing project paths: the bus is registered, nothing else starts (today's early return).
  captured <- new.env(parent = emptyenv())
  localWorkflowServerStubs("metabolomics", captured)
  testServer(mod_omic_workflow_server, args = list(project_dirs = list(), omic_type = "metabolomics",
    experiment_label = "w4-demo", config = omicGuiConfig("metabolomics")), {
    expect_identical(session$userData$omics_workflow_data$metabolomics, workflow_data)
    expect_identical(captured$order, character())
    expect_null(shiny::isolate(workflow_data$state_manager))
  })
})

test_that("KIT-102 W4: the import frame runs the panel's server in its own session and offers the completion helper", {
  seen <- new.env(parent = emptyenv())
  notes <- new.env(parent = emptyenv())
  notes$shown <- list()
  panel_server <- function(input, output, session, workflow_data, experiment_paths, volumes, config, complete) {
    seen$ns <- session$ns("vendor_format")
    seen$args <- list(workflow_data = workflow_data, experiment_paths = experiment_paths, volumes = volumes,
      omic = config$shell$omic_type)
    seen$complete <- complete
    "panel-result"
  }
  # The panels are phase C's (KIT-102 W5, W6): their servers are injected through the panel registry.
  panels <- omicImportSourcePanel
  local_mocked_bindings(omicImportSourcePanel = function(source_panel) {
    utils::modifyList(panels(source_panel), list(server = panel_server))
  }, .package = "MultiScholaR")
  local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
    type = "default", session = NULL) {
    notes$shown[[length(notes$shown) + 1L]] <- list(ui = ui, type = type)
    invisible(id)
  }, .package = "shiny")
  newBus <- function(design) shiny::reactiveValues(tab_status = list(setup_import = "pending", design_matrix = design,
    quality_control = "disabled"))

  # Feature table (metab): import completes and Design opens in the same write (the KIT-102 unification).
  bus <- newBus("disabled")
  seen$writes <- 0L
  testServer(mod_omic_import_server, args = list(id = "import", workflow_data = bus, experiment_paths = list(base_dir = "b"),
    volumes = NULL, config = omicGuiConfig("metabolomics")), {
    expect_identical(seen$ns, "import-vendor_format")
    expect_identical(seen$args, list(workflow_data = bus, experiment_paths = list(base_dir = "b"), volumes = NULL,
      omic = "metabolomics"))
    expect_identical(session$getReturned(), "panel-result")
    shiny::observe({
      bus$tab_status
      seen$writes <- seen$writes + 1L
    })
    session$flushReact()
    before <- seen$writes
    seen$complete()
    session$flushReact()
    expect_identical(shiny::isolate(bus$tab_status),
      list(setup_import = "complete", design_matrix = "pending", quality_control = "disabled"))
    expect_identical(seen$writes - before, 1L)
    expect_identical(notes$shown, list(list(ui = "Data imported successfully!", type = "message")))
  })

  # Feature table: a Design step already past `disabled` is left as it is (lipid's rule today).
  notes$shown <- list()
  bus <- newBus("complete")
  testServer(mod_omic_import_server, args = list(id = "import", workflow_data = bus, experiment_paths = list(),
    volumes = NULL, config = omicGuiConfig("lipidomics")), {
    seen$complete(notify = FALSE)
    expect_identical(shiny::isolate(bus$tab_status$design_matrix), "complete")
    expect_identical(shiny::isolate(bus$tab_status$setup_import), "complete")
    expect_identical(notes$shown, list())
  })

  # Proteomics: Design is set pending unconditionally (completeProtImportSuccessState today), its own success text.
  bus <- newBus("complete")
  testServer(mod_omic_import_server, args = list(id = "setup_import", workflow_data = bus, experiment_paths = list(),
    volumes = NULL, config = omicGuiConfig("proteomics")), {
    expect_identical(seen$ns, "setup_import-vendor_format")
    seen$complete()
    expect_identical(shiny::isolate(bus$tab_status$design_matrix), "pending")
    expect_identical(notes$shown, list(list(ui = "Data import successful!", type = "message")))
  })
})
