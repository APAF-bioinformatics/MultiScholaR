# KIT-096: the transcriptomics tabs' module CI (MCI-026 import, MCI-027 design, MCI-028 QC, MCI-029 the module and its
# stable ids). The fixtures are the KIT-040 generator's (tests/testdata/mci/fx/transcriptomics, written by
# generate-fixtures.R): a trx_screen corpus and its one-value fixtures. Each test has its own session owner; the tabs'
# own processing, builder and commit functions run under the transcriptomics profile.

module_ci_trx_dir <- function(...) file.path(module_ci_fixture_root(), "fx", "transcriptomics", ...)

module_ci_trx_corpus <- function() module_ci_trx_dir("trx_screen")

module_ci_trx_fixture <- function(name) module_ci_trx_dir("fixtures", name)

module_ci_trx_truth <- function(dir = module_ci_trx_corpus()) jsonlite::read_json(file.path(dir, "truth.json"), simplifyVector = TRUE)

#' A pending transcriptomics owner on a fresh project, with the tab bus the module server keeps.
module_ci_trx_fixture_owner <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  paths <- list(base_dir = root, source_dir = file.path(root, "source"), gene_qc_dir = file.path(root, "gene_qc"),
    da_output_dir = file.path(root, "da_genes"), publication_graphs_dir = file.path(root, "publication_graphs"),
    integration_dir = file.path(root, "integration"))
  for (p in paths[-1]) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  workflow <- new.env(parent = emptyenv())
  workflow$config_list <- list()
  workflow$processing_log <- list()
  workflow$tab_status <- list(setup_import = "pending", design_matrix = "disabled", quality_control = "disabled")
  omicsTestOwner(workflow = workflow, experiment_paths = paths, omic_type = "transcriptomics", bind = FALSE, env = env)
  list(workflow = workflow, paths = paths, session = omicsOwnerSession())
}

#' The import tab's processing on a count file, its notifications and errors captured.
module_ci_trx_import <- function(fixture, counts_file, id_col = NULL, name_col = NULL) {
  # KIT-102 W5: the shared feature-table panel's preview and processing under the transcriptomics configuration. The
  # result is the pre-KIT-102 processing's: NULL on a refusal, else its success list.
  config <- omicGuiConfig("transcriptomics")
  preview <- omicImportFeatureTablePreview(counts_file, config)
  seen <- new.env(parent = emptyenv())
  seen$notes <- list()
  seen$errors <- character()
  notify <- function(message, type = NULL, ...) seen$notes <- c(seen$notes, list(list(message = message, type = type)))
  run <- shiny::withReactiveDomain(fixture$session, processOmicImportFeatureTable(fixture$workflow,
    request = list(assay1_name = "counts", assay1_data = preview$assay_data, assay1_file = counts_file,
      id_col = id_col %||% preview$updates$id$selected, annotation_col = name_col %||% preview$updates$annotation$selected,
      sample_columns = preview$import_result$sample_columns, vendor_format = "trx_counts",
      detected_format = preview$detected_format, is_pattern = NULL, sanitize_names = NULL, assay2_file = NULL,
      assay2_name = NULL, experiment_paths = fixture$paths),
    config = config, complete = function(...) completeOmicImport(fixture$workflow, config, showNotification = notify),
    showNotification = notify, removeNotification = function(...) invisible(NULL),
    logError = function(message) seen$errors <- c(seen$errors, message), logInfo = function(...) invisible(NULL)))
  result <- if (identical(run$status, "success")) {
    list(status = "success", assayList = run$payload$assays, sampleColumns = run$payload$sample_columns,
      validationResult = run$payload$validation, artifactResult = run$payload$artifacts)
  }
  list(result = result, preview = preview, notes = seen$notes, errors = seen$errors)
}

#' The corpus's design (design_matrix.tab) with its technical-replicate sets from design.tsv.
module_ci_trx_corpus_design <- function(dir = module_ci_trx_corpus()) {
  design <- utils::read.delim(file.path(dir, "design_matrix.tab"), stringsAsFactors = FALSE)
  technical <- utils::read.delim(file.path(dir, "design.tsv"), stringsAsFactors = FALSE)
  design$tech_rep_group <- technical$tech_rep_group[match(design$Run, technical$sample)]
  design
}

#' Drive the builder (the design tab's shared input ids) on the imported tables, then save. Groups are added as factors in
#' `levels` order and assigned in that order, with each run's biological replicate; technical replicates are assigned
#' when `technical` is TRUE.
#' @return The builder's result (what the design tab commits).
module_ci_trx_build_design <- function(fixture, design, levels, technical = TRUE, remove = character(), contrast = TRUE) {
  workflow <- fixture$workflow
  out <- new.env(parent = emptyenv())
  shiny::testServer(mod_lipid_design_builder_server, args = list(
    data_tbl = function() getWorkflowOmicsBusTables(workflow, "data_tbl"), config_list = shiny::reactive(workflow$config_list),
    column_mapping = shiny::reactive(workflow$column_mapping), state_owner = function() workflow$state_manager,
    profile = trxOmicProfile()), {
    for (i in seq_along(levels)) {
      session$setInputs(new_factor = levels[[i]])
      session$setInputs(add_factor = i)
    }
    out$initial_runs <- shiny::isolate(design_matrix()$Run)
    clicks <- 0L
    for (group in levels) {
      runs <- design$Run[design$group == group]
      if (!length(runs)) next
      clicks <- clicks + 1L
      session$setInputs(selected_runs = runs, factor1_select = group, factor2_select = "", factor3_select = "", replicate_start = 1L)
      session$setInputs(assign_metadata = clicks)
    }
    # Each corpus biological sample's runs as technical replicates, numbered by the corpus's replicate (the driver's way).
    if (isTRUE(technical)) {
      tech <- 0L
      for (set in unique(design$tech_rep_group)) {
        runs <- design$Run[design$tech_rep_group == set]
        if (length(runs) < 2L) next
        tech <- tech + 1L
        session$setInputs(tech_rep_samples = runs, tech_rep_assignment_mode = "manual",
          manual_replicate_number = design$replicates[match(runs[[1L]], design$Run)])
        session$setInputs(assign_tech_reps = tech)
      }
    }
    if (length(remove)) {
      session$setInputs(samples_to_remove = remove)
      session$setInputs(remove_samples = 1L)
    }
    if (isTRUE(contrast) && length(levels) >= 2L) {
      session$setInputs(contrast_group1 = levels[[2L]], contrast_group2 = levels[[1L]], formula_string = "~ 0 + group")
      session$setInputs(add_contrast = 1L)
    } else {
      session$setInputs(formula_string = "~ 0 + group")
    }
    session$setInputs(save_results = 1L)
    out$result <- shiny::isolate(session$returned())
  })
  result <- out$result
  result$data_cln <- resolveWorkflowOwnedValue(result$data_cln, workflow$state_manager)
  attr(result, "initial_runs") <- out$initial_runs
  result
}

#' The design tab's commit of a builder result under the transcriptomics profile, its notifications captured.
module_ci_trx_commit <- function(fixture, result) {
  # A builder run passed as the argument runs here, outside the commit's reactive domain (testServer refuses one).
  force(result)
  seen <- new.env(parent = emptyenv())
  seen$notes <- list()
  testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    seen$notes <- c(seen$notes, list(list(message = paste(as.character(ui), collapse = " "), type = type)))
    invisible(NULL)
  }, showModal = function(...) invisible(NULL), removeModal = function(...) invisible(NULL), .package = "shiny")
  qc <- shiny::reactiveVal(NULL)
  shiny::withReactiveDomain(fixture$session, runLipidDesignBuilderObserverShell(results = result,
    workflowData = fixture$workflow, experimentPaths = fixture$paths, qcTrigger = function(value) seen$qc <- value,
    profile = trxOmicProfile()))
  list(notes = seen$notes, qc = seen$qc)
}

#' Import the corpus and commit its design through the builder: the design tab's state, ready for QC.
module_ci_trx_designed <- function(fixture, dir = module_ci_trx_corpus(), technical = TRUE) {
  module_ci_trx_import(fixture, file.path(dir, "counts.tsv"))
  design <- module_ci_trx_corpus_design(dir)
  result <- module_ci_trx_build_design(fixture, design, c("RPMI", "Sera"), technical = technical)
  module_ci_trx_commit(fixture, result)
  invisible(result)
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-097: the normalisation, DA and summary tabs (MCI-030 norm, MCI-031 DA, MCI-032 summary).
# ---------------------------------------------------------------------------------------------------------------

#' The builder's contrasts table (the lipid shape KIT-096's design tab saves) from contrast entries.
module_ci_trx_builder_contrasts <- function(entries) {
  parsed <- lapply(entries, trxParseContrastEntry)
  data.frame(contrasts = vapply(parsed, `[[`, "", "expr"), friendly_names = vapply(parsed, `[[`, "", "name"),
    full_format = entries, stringsAsFactors = FALSE)
}

#' A transcriptomics owner at `trx_qc_complete` on a setupDirectories project (the paths the application passes): the
#' corpus's counts and design (with its technical replicates when `technical`), the design state with its provenance,
#' the filter and the QC finalise through the wrappers.
#' @param contrasts Contrast entries; the corpus's contrast file when NULL.
module_ci_trx_project <- function(dir, technical = FALSE, contrasts = NULL, env = parent.frame()) {
  testthat::skip_if_not(workflowOwnerCapable())
  root <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", .local_envir = env)
  invisible(utils::capture.output(project_dirs <- suppressMessages(setupDirectories(base_dir = file.path(root, "project", "kit097"),
    omic_types = "transcriptomics", label = NULL, force = TRUE, reuse_existing = TRUE))))
  paths <- project_dirs$transcriptomics
  dir.create(file.path(paths$source_dir, "import_sources"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(dir, "counts.tsv"), file.path(paths$source_dir, "import_sources", "counts.tsv"))
  counts <- vroom::vroom(file.path(dir, "counts.tsv"), show_col_types = FALSE)
  design <- if (technical) trxDeclTechnicalDesign(dir) else vroom::vroom(file.path(dir, "design_matrix.tab"), show_col_types = FALSE)
  utils::write.table(design, file.path(paths$source_dir, "design_matrix.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
  entries <- contrasts %||% trxParseContrastFile(file.path(dir, "contrast_strings.tab"))
  session <- omicsOwnerSession()
  workflow <- list2env(list(state_update_trigger = NULL, config_list = list(), contrasts_tbl = module_ci_trx_builder_contrasts(entries),
    design_matrix = as.data.frame(design), processing_log = list(),
    tab_status = list(quality_control = "complete", normalization = "pending", differential_analysis = "disabled", session_summary = "disabled")))
  workflow$state_manager <- createDIASessionOwner(workflow, paths, session, omic_type = "transcriptomics")
  shiny::withReactiveDomain(session, bindWorkflowOwnerTypeAtDesign(workflow, "TRANSCRIPTOMICS", "trx_counts"))
  object <- createTranscriptomicsAssayData(list(counts = counts), design, c("RPMI", "Sera"),
    roles = list(technical_replicate_id = if (technical) "tech_rep_group" else NA_character_))
  owner <- workflow$state_manager
  recordTrxDesignProvenance(owner, object, source_dir = paths$source_dir)
  owner$saveState("trx_raw_data_s4", object, list(), "Initial TranscriptomicsAssayData S4 object created after design matrix")
  runTrxReferenceExpressionFilter(owner, list())
  commitWorkflowOmicsStoreStep(owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
  withr::defer(closeDIAStateAdapter(owner), envir = env)
  list(owner = owner, workflow = workflow, session = session, paths = paths, project_dirs = project_dirs, source_dir = paths$source_dir,
    object = object, dir = dir)
}

#' A quiet progress for a managed runner under test.
module_ci_trx_quiet_progress <- function(message) {
  list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL), close = function(...) invisible(NULL))
}

#' The normalisation, DA and summary tabs under the transcriptomics profile, as the module starts them, with the runner when
#' managed. The DA tab's values are kept on the session for the test.
module_ci_trx_tabs_server <- function(fixture, managed = FALSE) {
  force(fixture)
  force(managed)
  function(input, output, session) {
    if (managed) {
      runner <- installWorkflowOmicsActions(session, fixture$workflow, "transcriptomics", "kit097", fixture$paths)
      runner$progress <- module_ci_trx_quiet_progress
    }
    profile <- trxOmicProfile()
    # A tab that cannot start is recorded, so that the journey's assertions report what is missing.
    start <- function(name, expr) tryCatch(expr, error = function(error) session$userData$start_errors[[name]] <- conditionMessage(error))
    start("norm", mod_lipid_norm_server("norm", fixture$workflow, fixture$paths, "transcriptomics", "kit097", profile = profile))
    start("de", session$userData$da <- mod_lipid_da_server("de", fixture$workflow, fixture$paths, "transcriptomics", "kit097", profile = profile))
    start("summary", mod_lipid_summary_server("summary", fixture$project_dirs, "transcriptomics", "kit097", fixture$workflow, profile = profile))
  }
}

#' Drive the tabs: the normalisation with `norm` inputs, then (unless `da = FALSE`) the DA, then each RUVIII export in
#' `ruviii` (the M sources, in order). A failure to drive a step is kept (`out$errors`), so that the test's assertions on
#' the record and the tab report what is missing.
#' @return A list: `norm_log`, `ruv_summary`, `da` (the DA tab's results list), `summary_stats`, `ruviii_status`, `errors`,
#'   `notes` (the notifications shown).
module_ci_trx_tabs_journey <- function(fixture, managed = FALSE, norm = list(ruv_mode = "automatic"), da = TRUE, ruviii = character(),
                                       summary = NULL) {
  out <- new.env(parent = emptyenv())
  out$errors <- character()
  out$notes <- list()
  step <- function(name, expr) tryCatch(expr, error = function(error) out$errors[[name]] <- conditionMessage(error))
  testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
    out$notes <- c(out$notes, list(list(message = paste(as.character(ui), collapse = " "), type = type)))
    invisible(NULL)
  }, .package = "shiny")
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(module_ci_trx_tabs_server(fixture, managed), {
    step("norm", {
      do.call(session$setInputs, stats::setNames(norm, paste0("norm-", names(norm))))
      session$setInputs(`norm-run_normalization` = 1L)
      session$flushReact()
      out$norm_log <- output$`norm-norm_log`
      out$ruv_summary <- output$`norm-ruv_summary_counts`
    })
    if (isTRUE(da)) step("da", {
      session$setInputs(`de-da_q_val_thresh` = 0.05, `de-treat_lfc_cutoff` = 1, `de-run_da_analysis` = 1L)
      session$flushReact()
      out$da <- shiny::isolate(session$userData$da$da_results_list)
      out$summary_stats <- output$`de-da_summary_stats`
      out$da_reference <- fixture$owner$.context$snapshot$da_results$entity_id
    })
    for (i in seq_along(ruviii)) step(paste0("ruviii_", i), {
      session$setInputs(`summary-ruviii_m_source` = ruviii[[i]])
      session$setInputs(`summary-run_ruviii_export` = i)
      session$flushReact()
      out$ruviii_status <- output$`summary-ruviii_status`
    })
    if (is.function(summary)) step("summary", summary(session, output, out))
    out$errors <- c(unlist(session$userData$start_errors), out$errors)
  }))))
  as.list(out)
}

#' The DA tab's long table as the topTags table it was annotated from (KIT-093's columns, rows named by gene id).
module_ci_trx_top_table <- function(long) {
  table <- as.data.frame(long[, c("logFC", "logCPM", "F", "PValue", "FDR")])
  rownames(table) <- long$GeneID
  table
}

#' The first contrast's rows of the DA tab's long table.
module_ci_trx_first_contrast <- function(results) {
  long <- results$da_genes_long
  long[long$comparison == unique(long$comparison)[[1L]], , drop = FALSE]
}
