# KIT-022 item 1: the DA contrast identities before and after the argument commit moves ahead of the contrast
# loop. On a reflinked copy of the KIT-008 screen member's project the DA kernel runs (publication deferred) and
# the manifest's source, the head, each contrast's result-table and product entity ids and the state's products
# are recorded. Run once from a worktree at the pre-ticket commit (`before`) and once from the ticket (`after`);
# identical table and product ids mean the contrast identity is unchanged and the products would be cache hits.
# KIT-035 form of the KIT-022 identity check: the same in-process DA run on the KIT-008 screen project, plus the limma
# warnings the run emits and the base fit's recorded factors (block reason, duplicate correlation), recorded before and
# after the technical-replicate declaration change so the contrast tables are compared by content.
# KIT-025 form: the KIT-035 in-process DA run on the KIT-008 screen project, with every state fingerprinted before
# the store is reopened on the provenance-recording build and again after the DA step, and the head's step
# provenance read for every name in state_history. Record invariance: every state fingerprint identical; the head
# sha differs (additive fields) and is recorded, not judged.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-025/provenance-check.R <label>
suppressMessages(devtools::load_all(quiet = TRUE))
main <- function(label) {
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
member <- file.path(b$root, "kit-gui-008", "member")
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit022-da-identity-", label))
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project", "dia_baseline_smoke"), work)) == 0L)
project <- file.path(work, "dia_baseline_smoke")
results <- file.path(project, "results", "proteomics")
paths <- list(base_dir = project, results_dir = results, source_dir = file.path(project, "scripts", "proteomics"),
    protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
    da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
    pathway_dir = file.path(results, "pathway_enrichment"), time_dir = file.path(results, "publication_graphs"),
    integration_dir = file.path(project, "integration"))
for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
source_dir <- normalizePath(paths$source_dir, mustWork = TRUE)
archive <- list.files(source_dir, "[.]msdia[.]tar$", full.names = TRUE)
roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
stopifnot(length(archive) == 1L, length(roots) == 1L)
api <- function(name) get(name, asNamespace("MultiScholaR"))
store_dir <- file.path(roots, ".workflow-entities")
fingerprints_before <- h$kitFingerprintStore(store_dir, api)
identity <- readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
record <- list(schema = "kit.provenance_check.v1", ticket = "KIT-025", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
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
    da_data <- shiny::reactiveValues(current_s4_object = view$current_s4_object,
        contrasts_available = view$contrasts_tbl$contrasts, analysis_complete = FALSE, da_results_list = NULL)
    manager <- workflow$state_manager
    source_before <- getDIAStateReference(manager)
    settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
    limma_warnings <- character()
    manifest <- withCallingHandlers(runDIAReferenceDaAnalysis(settings, da_data, workflow, paths, publish = FALSE),
        warning = function(w) { limma_warnings <<- c(limma_warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    record$limma_warnings <<- limma_warnings
    snapshot <- manager$.context$snapshot
    history <- unlist(snapshot$state_history)
    provenance <- snapshot$step_provenance
    record$state_history <<- history
    record$provenance_complete <<- all(history %in% names(provenance))
    record$step_provenance <<- lapply(provenance[history], function(p) p[c("operation", "variant", "parameter_names", "build", "source_commit", "defaults_sha256", "seed", "committed_at")])
    registered <- setdiff(unique(vapply(provenance[history], `[[`, "", "operation")), c("legacy", "design_builder", "import"))
    record$registered_operations_resolve <<- all(vapply(registered, function(op) !is.null(tryCatch(workflowActionSpec(op), error = function(e) NULL)), logical(1)))
    record$parameter_names_match_registry <<- all(vapply(registered, function(op) {
        entries <- Filter(function(p) identical(p$operation, op), provenance[history])
        all(vapply(entries, function(p) identical(p$parameter_names, workflowActionSpec(op)$parameters), logical(1)))
    }, logical(1)))
    build <- workflowSoftwareIdentity()$build
    record$software_environment <<- snapshot$software_environments[[build]][c("version", "source_commit", "r_version")]
    record$software_packages <<- snapshot$software_environments[[build]]$packages[c("MultiScholaR", "limma", "iq", "arrow", "duckdb")]
    fit_reference <- manager$.context$da_base_fit_reference
    record$base_fit_factors <<- fit_reference$metadata$factors
    record$manifest_model <<- manifest$model
    head <- getDIAStateReference(manager)
    state <- readDIAStateMetadata(manager, head)$data
    record$state_before_da <<- source_before$entity_id
    record$manifest_source <<- manifest$source$entity_id
    record$head_after_da <<- head$entity_id
    record$source_is_head <<- identical(manifest$source$entity_id, head$entity_id)
    # Entity ids are per save (a fresh copy of the store assigns new ones), so the comparison is by content: the
    # result table's rows, each product's value and the contrast operation's identity (the stripped arguments).
    record$contrasts <<- withDIAStateScope(manager, function() stats::setNames(lapply(names(manifest$results), function(name) {
        result <- manifest$results[[name]]
        table <- materializeDIATable(manager$.context, result$tables$da_proteins_long)
        products <- lapply(result$products, function(p) {
            value <- hydrateSavedWorkflowEntity(manager$.context$store, p)$value
            digest::digest(value$table %||% value, algo = "sha256")
        })
        i <- match(name, manifest$contrasts$contrasts)
        identity <- diaComponentScientificIdentity(manager, manifest$source, "DA-contrast-output", c(manifest$parameters,
            list(contrast = manifest$contrasts$full_format[[i]], annotations = state$components$protein_id_table$entity_id,
                name = name, comparison = manifest$contrasts$friendly_names[[i]])))
        list(rows = nrow(table), table_sha256 = digest::digest(table, algo = "sha256"), products = products,
            identity_sha256 = digest::digest(identity, algo = "sha256"))
    }), names(manifest$results)))
    record$state_products <<- sort(unname(vapply(state$products, function(p) p$entity_id, character(1))))
    # Diagnostic: what the commit changed on the state's slots (the identity reads them with the DA args stripped).
    before_slots <- readDIAStateMetadata(manager, source_before)$data$slots
    after_slots <- state$slots
    record$slot_names_changed <<- sort(unique(c(setdiff(names(after_slots), names(before_slots)),
        names(before_slots)[!vapply(names(before_slots), function(n) identical(before_slots[[n]], after_slots[[n]]), logical(1))])))
    record$arg_names_changed <<- sort(unique(c(setdiff(names(after_slots$args), names(before_slots$args)),
        names(before_slots$args)[!vapply(names(before_slots$args), function(n) identical(before_slots$args[[n]], after_slots$args[[n]]), logical(1))])))
    # A second run with the same settings on the same store: identical product entity ids mean cache hits.
    first_ids <- lapply(manifest$results, function(r) c(r$tables$da_proteins_long$entity_id, unname(vapply(r$products, function(p) p$entity_id, character(1)))))
    again <- runDIAReferenceDaAnalysis(settings, da_data, workflow, paths, publish = FALSE)
    second_ids <- lapply(again$results, function(r) c(r$tables$da_proteins_long$entity_id, unname(vapply(r$products, function(p) p$entity_id, character(1)))))
    head_after_second <- getDIAStateReference(manager)$entity_id
    second <- list(hits = identical(first_ids, second_ids), source_is_head = identical(again$source$entity_id, head_after_second),
        holder_is_head = identical(shiny::isolate(da_data$current_s4_object)$entity_id, head_after_second))
    record$second_run <<- second
    # The identity of the contrast operation on the pre-commit source and on the head, in one process (the
    # software identity differs between source trees under load_all, so cross-run identity hashes are not comparable).
    identity_of <- function(reference) withDIAStateScope(manager, function() {
        i <- 1L
        digest::digest(diaComponentScientificIdentity(manager, reference, "DA-contrast-output", c(manifest$parameters,
            list(contrast = manifest$contrasts$full_format[[i]], annotations = state$components$protein_id_table$entity_id,
                name = manifest$contrasts$contrasts[[i]], comparison = manifest$contrasts$friendly_names[[i]]))), algo = "sha256")
    })
    record$identity_equal_in_process <<- identical(identity_of(source_before), identity_of(head))
    record$da_args <<- state$slots$args$daAnalysisUI[c("q_value_threshold", "log_fold_change_cutoff", "formula_string")]
    closeDIAStateAdapter(manager)
})
# The adoption replaces the store folder (KIT-013), so the after-fingerprint reads the live owner's directory.
live_dir <- shiny::isolate(workflow$state_manager)$.context$store$directory
fingerprints_after <- h$kitFingerprintStore(live_dir, api)
same <- vapply(names(fingerprints_before$states), function(name) {
    a <- fingerprints_before$states[[name]]; b <- fingerprints_after$states[[name]]
    !is.null(b) && identical(a$components, b$components) && identical(a$support, b$support) && identical(a$summary_sha256, b$summary_sha256)
}, logical(1))
record$state_fingerprints_identical <- all(same)
record$states_compared <- names(same)
# The DA step commits its arguments onto the head state (KIT-022), so that one state's summary moves by the run itself.
record$states_differing <- names(same)[!same]
record$head_sha256 <- list(before = fingerprints_before$head_sha256, after = fingerprints_after$head_sha256)
cat("PROVENANCE complete:", record$provenance_complete, "| registered resolve:", record$registered_operations_resolve,
    "| parameter names match:", record$parameter_names_match_registry, "| state fingerprints identical:", record$state_fingerprints_identical,
    "(", length(same), "states; differing:", paste(record$states_differing, collapse = ","), ") | head differs:", !identical(record$head_sha256$before, record$head_sha256$after), "\n")
print(t(sapply(record$step_provenance, function(p) c(operation = p$operation, variant = as.character(p$variant), params = paste(p$parameter_names, collapse = ",")))))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-025", paste0("provenance-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", label, "source is head:", record$source_is_head, " table:", substr(record$contrasts[[1]]$table_sha256, 1, 12), " identity:", substr(record$contrasts[[1]]$identity_sha256, 1, 12), " second run hits:", record$second_run$hits, " identity equal in process:", record$identity_equal_in_process, "\n")
}
main(commandArgs(TRUE)[[1L]])
