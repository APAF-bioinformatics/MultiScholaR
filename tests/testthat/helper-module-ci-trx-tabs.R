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
  preview <- trxImportPreview(counts_file)
  seen <- new.env(parent = emptyenv())
  seen$notes <- list()
  seen$errors <- character()
  result <- shiny::withReactiveDomain(fixture$session, runTrxImportProcessing(
    workflowData = fixture$workflow, assay1Name = "counts", assay1Data = preview$assayData, assay1File = counts_file,
    lipidIdCol = id_col %||% preview$updates$lipidId$selected, annotationCol = name_col %||% preview$updates$annotation$selected,
    sampleColumns = preview$importResult$sample_columns, vendorFormat = "trx_counts", detectedFormat = preview$detectedFormat,
    isPattern = NULL, sanitizeNames = NULL, assay2File = NULL, assay2Name = NULL, experimentPaths = fixture$paths,
    notify = function(message, type = NULL, ...) seen$notes <- c(seen$notes, list(list(message = message, type = type))),
    removeNotify = function(...) invisible(NULL), logError = function(message) seen$errors <- c(seen$errors, message)))
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
