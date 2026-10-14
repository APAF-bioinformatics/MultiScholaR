# KIT-019: the DA publication outside the application process, recorded on the KIT-008 screen member's project
# (8 samples, the on-disk default). The application-side process opens a reflinked copy as the DA tab's Open does,
# runs the DA kernel in process with the publication deferred (`publish = FALSE`), then a separate R process runs
# the registered `da_publication` worker handler on the same store (the request directory as the worker gets it);
# the application-side process completes and publishes the worker's files into the project (the controller's two
# phases) and records whether Glimma's namespace was loaded in it, its anonymous memory before and after the DA
# step, the published files and the saved page assets. Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-019/da-publication-probe.R
suppressMessages(devtools::load_all(quiet = TRUE))
main <- function() {
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
member <- file.path(b$root, "kit-gui-008", "member")
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit019-da-publication-probe")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project", "dia_baseline_smoke"), work)) == 0L)
project <- file.path(work, "dia_baseline_smoke")
results <- file.path(project, "results", "proteomics")
paths <- list(base_dir = project, results_dir = results, source_dir = file.path(project, "scripts", "proteomics"),
    protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
    da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
    pathway_dir = file.path(results, "pathway_enrichment"),
    time_dir = file.path(results, "publication_graphs"), integration_dir = file.path(project, "integration"))
for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
source_dir <- normalizePath(paths$source_dir, mustWork = TRUE)
archive <- list.files(source_dir, "[.]msdia[.]tar$", full.names = TRUE)
stopifnot(length(archive) == 1L)
roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
stopifnot(length(roots) == 1L)
anon <- function() {
    line <- grep("^RssAnon:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
    as.numeric(strsplit(trimws(sub("^[^:]+:", "", line)), " +")[[1L]][[1L]]) * 1024
}
timed <- function(expr) { started <- proc.time()[["elapsed"]]; value <- force(expr); list(value = value, seconds = round(proc.time()[["elapsed"]] - started, 1)) }
identity <- readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
record <- list(schema = "kit.da_publication_probe.v1", ticket = "KIT-019", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), project_label = "kit-gui-008 member (screen, 8 samples)")
workflow <- shiny::reactiveValues(state_manager = manager, data_type = "peptide", data_format = "diann",
    tab_status = list(), contrasts_tbl = NULL, uniprot_dat_cln = NULL, design_matrix = NULL, config_list = NULL,
    da_analysis_results_list = NULL, da_ui_params = NULL, qc_params = list())
shiny::testServer(function(input, output, session) {}, {
    workflow$state_manager <- attachDIAStateSession(manager, session, workflow)
    imported <- openDIANativeSession(archive, source_dir)
    view <- adoptDIANativeSession(imported, workflow, session, source_dir, paths)
    workflow$contrasts_tbl <- view$contrasts_tbl
    workflow$design_matrix <- view$design_matrix
    workflow$config_list <- view$config_list
    uniprot_path <- file.path(source_dir, "uniprot_dat_cln.RDS")
    workflow$uniprot_dat_cln <- if (is.data.frame(view$uniprot_dat_cln)) view$uniprot_dat_cln else
        if (file.exists(uniprot_path)) readRDS(uniprot_path) else NULL
    da_data <- shiny::reactiveValues(current_s4_object = view$current_s4_object,
        contrasts_available = view$contrasts_tbl$contrasts, analysis_complete = FALSE, da_results_list = NULL)
    manager <- workflow$state_manager
    namespaces_before <- loadedNamespaces()
    record$application$glimma_loaded_before_da <<- "Glimma" %in% namespaces_before
    record$application$anon_before_da_bytes <<- anon()
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    # The DA kernel in process, the publication deferred (the module's managed path).
    kernel <- timed(runDIAReferenceDaAnalysis(settings, da_data, workflow, paths, publish = FALSE))
    manifest <- kernel$value
    stopifnot(isDIAResultManifest(manifest))
    record$application$kernel_seconds <<- kernel$seconds
    record$application$anon_after_kernel_bytes <<- anon()
    head <- getDIAStateReference(manager)
    annotations <- shiny::isolate(workflow$uniprot_dat_cln)
    annotations_reference <- if (is.data.frame(annotations) && nrow(annotations)) saveDIADisplayTable(manager, annotations) else NULL
    # The worker: a separate R process running the registered handler on the same store from a request directory.
    request_dir <- file.path(work, "request")
    dir.create(request_dir)
    saveRDS(list(experiment_label = "dia_baseline_smoke", experiment_paths = paths, omic_type = "proteomics"),
        file.path(request_dir, "publication-context.rds"), version = 3)
    request <- list(operation = "da_publication", input_references = list(state = head),
        parameters = list(annotations = annotations_reference), directory = request_dir, worker_cap = 12 * 1024^3,
        project_root = dirname(manager$.context$store$directory), owner = manager$.context$store$identity)
    saveRDS(request, file.path(request_dir, "probe-request.rds"), version = 3)
    worker_script <- file.path(work, "worker.R")
    writeLines(c('suppressMessages(devtools::load_all(quiet = TRUE))',
        'request <- readRDS(commandArgs(TRUE)[[1L]])',
        'store <- newWorkflowEntityStore(request$project_root, request$owner$project_id, request$owner$workflow_id)',
        'manager <- newDIAStateAdapter(store, reopen = TRUE)',
        'anon <- function() { line <- grep("^RssAnon:", readLines("/proc/self/status", warn = FALSE), value = TRUE); as.numeric(strsplit(trimws(sub("^[^:]+:", "", line)), " +")[[1L]][[1L]]) * 1024 }',
        'started <- proc.time()[["elapsed"]]',
        'value <- computeWorkflowDaPublicationAction(manager, request, "whole")',
        'saveRDS(list(value = value, seconds = round(proc.time()[["elapsed"]] - started, 1), glimma_loaded = "Glimma" %in% loadedNamespaces(), anon_bytes = anon()), file.path(request$directory, "probe-worker-result.rds"), version = 3)',
        'closeDIAStateAdapter(manager)'), worker_script)
    store_root <- dirname(manager$.context$store$directory)
    store_identity <- manager$.context$store$identity
    closeDIAStateAdapter(manager)
    status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(worker_script), shQuote(file.path(request_dir, "probe-request.rds"))),
        stdout = file.path(work, "worker.log"), stderr = file.path(work, "worker.log"))
    if (!identical(status, 0L)) {
        cat(tail(readLines(file.path(work, "worker.log"), warn = FALSE), 40), sep = "\n")
        stop("the worker process failed (status ", status, ")")
    }
    worker <- readRDS(file.path(request_dir, "probe-worker-result.rds"))
    record$worker <<- list(seconds = worker$seconds, glimma_loaded = worker$glimma_loaded, anon_bytes = worker$anon_bytes,
        files = length(worker$value$files), bytes = worker$value$summary$bytes)
    # The application side again: completion and publication (the controller's phases), then the display.
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, store_identity$project_id, store_identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    manager$.context$session_format <- "native_v1"
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    workflow$state_manager <- attachDIAStateSession(manager, session, workflow)
    session$userData$dia_capacity_manager <- manager
    completed <- withDIAStateScope(manager, function() completeDIADaPublication(list(state_manager = manager), worker$value))
    stopifnot(identical(completed$stateReference$entity_id, head$entity_id))
    owner <- list(manager = manager, control = list(generation = 1L))
    published <- timed(withDIAStateScope(manager, function() publishWorkflowDaResult(owner,
        list(request_dir = request_dir, root = work), list(operation = "da_publication", state_reference = head))))
    read <- readManagedDaResult(manager, list(state_reference = head), workflow)
    completeDIADaDisplay(da_data, workflow, read$manifest)
    contrast <- read$manifest$contrasts$friendly_names[[1L]]
    shown <- withDIAStateScope(manager, function() timed(generateProtDAVolcanoPlotGlimma(read$manifest, contrast,
        da_q_val_thresh = 0.05, args_row_id = read$manifest$protein_id_column)))
    record$application$publication_seconds <<- published$seconds
    record$application$anon_after_publication_bytes <<- anon()
    record$application$glimma_loaded_after_da <<- "Glimma" %in% loadedNamespaces()
    record$application$namespaces_added <<- sort(setdiff(loadedNamespaces(), namespaces_before))
    record$application$status_text <<- read$resultText
    record$application$interactive_shown <<- inherits(shown$value, "shiny.tag") && identical(shown$value$name, "iframe")
    record$application$interactive_seconds <<- shown$seconds
    record$project$da_exports <<- list.files(paths$da_output_dir)
    record$project$volcano_files <<- list.files(file.path(paths$publication_graphs_dir, "Volcano_Plots"))
    record$project$interactive_files <<- list.files(file.path(paths$publication_graphs_dir, "Interactive_Volcano_Plots"))
    record$project$numsig_files <<- list.files(file.path(paths$publication_graphs_dir, "NumSigDaMolecules"))
    record$project$completed_outputs <<- manager$.context$snapshot$completed_outputs[[head$entity_id]]
    closeDIAStateAdapter(manager)
})
record$result <- list(glimma_kept_out_of_application = !record$application$glimma_loaded_after_da && record$worker$glimma_loaded,
    anon_delta_da_step_bytes = record$application$anon_after_publication_bytes - record$application$anon_before_da_bytes)
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-019/da-publication-probe.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT glimma in application:", record$application$glimma_loaded_after_da, " in worker:", record$worker$glimma_loaded,
    " anon delta GiB:", round(record$result$anon_delta_da_step_bytes / 2^30, 3), " shown:", record$application$interactive_shown, "\n")
}
main()
