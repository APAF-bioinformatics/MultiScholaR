# KIT-016: the enrichment declaration and the DA display products exercised in process on the KIT-008 screen
# member's project (8 samples, the on-disk default), recorded. A reflinked copy of the project is opened as the
# DA tab's "Load filtered session" does (the KIT-012 probe's path), the DA step is rerun through the KIT-004
# declaration (its publication now also saves each contrast's interactive volcano page as a rendered asset),
# then: the static volcano display product and the interactive asset are resolved for the first contrast; the
# enrichment module body runs on the manifest in test mode (the deterministic double, never the live service),
# through the enrichment declaration; the declaration's results are compared with the pre-ticket kernel path
# (the double called directly on the module's process arguments) by hash; the completed state is the KIT-003
# result class on the adapter. Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-016/in-process-enrichment-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
# The body runs inside a function so `on.exit` removes the working copy (audit PH-04).
main <- function() {
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
member <- file.path(b$root, "kit-gui-008", "member")
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit016-enrichment-check")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project", "dia_baseline_smoke"), work)) == 0L)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
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
identity <- readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
timed <- function(expr) { started <- proc.time()[["elapsed"]]; value <- force(expr); list(value = value, seconds = round(proc.time()[["elapsed"]] - started, 1)) }
record <- list(schema = "kit.in_process_enrichment_check.v1", ticket = "KIT-016", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), project = "kit-gui-008 member (screen, 8 samples)")
options(multischolar.test_mode = TRUE)
workflow <- shiny::reactiveValues(state_manager = manager, data_type = "peptide", data_format = "diann",
    tab_status = list(), contrasts_tbl = NULL, uniprot_dat_cln = NULL, design_matrix = NULL, config_list = NULL,
    da_analysis_results_list = NULL, da_ui_params = NULL, qc_params = list())
shiny::testServer(function(input, output, session) {}, {
    workflow$state_manager <- attachDIAStateSession(manager, session, workflow)
    imported <- openDIANativeSession(archive, source_dir)
    view <- adoptDIANativeSession(imported, workflow, session, source_dir, paths)
    stopifnot(isSavedWorkflowReference(view$current_s4_object))
    workflow$contrasts_tbl <- view$contrasts_tbl
    workflow$design_matrix <- view$design_matrix
    workflow$config_list <- view$config_list
    uniprot_path <- file.path(source_dir, "uniprot_dat_cln.RDS")
    workflow$uniprot_dat_cln <- if (is.data.frame(view$uniprot_dat_cln)) view$uniprot_dat_cln else
        if (file.exists(uniprot_path)) readRDS(uniprot_path) else NULL
    da_data <- shiny::reactiveValues(current_s4_object = view$current_s4_object,
        contrasts_available = view$contrasts_tbl$contrasts, analysis_complete = FALSE, da_results_list = NULL)
    record$glimma_loaded_before_da <<- "Glimma" %in% loadedNamespaces()
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    analysed <- timed(runDIAReferenceDaAnalysis(settings, da_data, workflow, paths))
    manifest <- analysed$value
    stopifnot(isDIAResultManifest(manifest))
    record$glimma_loaded_after_da_publication <<- "Glimma" %in% loadedNamespaces()
    contrast <- manifest$contrasts$friendly_names[[1L]]
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        products_before <- length(manager$.context$snapshot$products)
        product <- timed(generateDIAReferenceVolcanoProduct(manifest, contrast))
        again <- timed(generateDIAReferenceVolcanoProduct(manifest, contrast))
        plot <- timed(generateProtDAVolcanoStatic(manifest, contrast, da_q_val_thresh = 0.05, lfc_threshold = 0))
        page <- timed(resolveDIAInteractiveVolcanoAsset(manifest, contrast, 0.05))
        shown <- timed(generateProtDAVolcanoPlotGlimma(manifest, contrast, da_q_val_thresh = 0.05, args_row_id = manifest$protein_id_column))
        record$volcano <<- list(contrast = contrast, product_rows = product$value$rows, product_columns = product$value$columns,
            product_seconds = product$seconds, hit_seconds = again$seconds, same_product = identical(digest::digest(product$value$data), digest::digest(again$value$data)),
            products_added = length(manager$.context$snapshot$products) - products_before,
            static_plot = inherits(plot$value, "ggplot"), static_seconds = plot$seconds,
            interactive_asset_bytes = if (is.null(page$value)) NA_integer_ else nchar(page$value, type = "bytes"),
            interactive_shown = inherits(shown$value, "shiny.tag") && identical(shown$value$name, "iframe"),
            interactive_seconds = shown$seconds)
    })
    # Enrichment on the manifest through the declaration (test mode: the deterministic double).
    enrichmentData <- createProtEnrichReactiveValues()
    enrichmentData$current_s4_object <- manifest$source
    enrichmentData$da_results_data <- manifest
    enrichmentData$contrasts_available <- manifest$contrasts$friendly_names
    workflow$da_analysis_results_list <- manifest
    input <- list(selected_contrast = contrast, organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05,
        enable_organism_filter = FALSE, exclude_iea = FALSE, correction_method = "gSCS")
    supported <- buildProtEnrichSupportedOrganisms()
    method_fn <- function() resolveProtEnrichAnalysisMethod(input$organism_taxid, supported)
    captured <- new.env()
    prepare <- function(...) { config <- prepareProtEnrichProcessExecution(...); captured$args <- config$enrichmentArgs; config }
    completed <- timed(runProtEnrichAnalysisBody(input = input, enrichmentData = enrichmentData, workflowData = workflow,
        experimentPaths = paths, currentAnalysisMethodFn = method_fn, prepareProcessExecutionFn = prepare,
        showNotificationFn = function(...) NULL, incProgressFn = function(...) NULL, catFn = function(...) NULL))
    declared <- shiny::isolate(enrichmentData$enrichment_results_full)
    stopifnot(methods::is(declared, "EnrichmentResults"))
    direct_args <- captured$args$processArgs
    direct_args$pathway_dir <- file.path(work, "pathway-direct"); dir.create(direct_args$pathway_dir)
    direct <- do.call(runProtEnrichDeterministicProcessEnrichments, direct_args)
    product <- hydrateSavedWorkflowEntity(manager$.context$store, shiny::isolate(enrichmentData$enrichment_product))$value
    state <- manager$current_state
    restored <- manager$getState("enrichment_completed")
    record$enrichment <<- list(contrast = contrast, seconds = completed$seconds, analysis_method = completed$value$analysisMethod,
        service = shiny::isolate(enrichmentData$enrichment_service), product_service = product$service,
        product_results = unname(product$results),
        exact_against_direct_double = identical(digest::digest(declared@enrichment_data), digest::digest(direct@enrichment_data)) &&
            identical(digest::digest(declared@enrichment_summaries), digest::digest(direct@enrichment_summaries)),
        holder_is_reference = isSavedWorkflowReference(shiny::isolate(enrichmentData$current_s4_object)),
        current_state = state, completed_state_class = class(restored)[[1L]],
        restored_exact = identical(digest::digest(restored@enrichment_data), digest::digest(declared@enrichment_data)),
        contrasts = names(declared@enrichment_data))
})
closeDIAStateAdapter(shiny::isolate(workflow$state_manager))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-016/in-process-enrichment-check.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT enrichment exact:", record$enrichment$exact_against_direct_double, " volcano shown:", record$volcano$interactive_shown,
    " state:", record$enrichment$current_state, "\n")
}
main()
