# KIT-010: the protein replicate filter on the row partitioner checked against the KIT-005 oracle on the real
# 300-sample protein population before the counted members. Input: the KIT-005 control journey's state after
# duplicate removal (`journey_duplicates`, 2,500 proteins x 297 samples), materialised as a snapshot store per
# variant by the harness so the evidence store is never written. Each variant (`whole`: the GUI's in-process
# path; `rows`: the partitioner) runs in process from the working tree with the journey's recorded parameters
# (grouping `group`, 4 cores recorded), the state is committed through the protein completion, and its
# fingerprint under the current canonical form is compared with the `journey_protein_replicate` oracle as the
# member verifier compares it (components, support entries, args summary; products reported). Not a counted
# run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-010/in-process-replicate-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
input <- b$populations$journey_duplicates
oracle <- b$oracles$journey_protein_replicate
stopifnot(identical(input$state, "duplicates_removed"), is.null(h$kitPopulationInvalid(input)), !is.null(oracle$class))
parameters <- list(grouping_variable = "group", parallel_cores = 4L)
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit010-replicate-check")
unlink(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
run <- function(variant) {
    snapshot <- h$kitPopulationStore(input, file.path(work, variant, "population-store"), api)
    identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    reference <- getDIAStateReference(manager)
    stopifnot(identical(manager$current_state, "duplicates_removed"), identical(reference$entity_id, input$entity_id))
    started <- proc.time()[["elapsed"]]
    value <- computeDIAReferenceProteinReplicateFilter(manager, reference, parameters, variant = variant)
    compute_seconds <- round(proc.time()[["elapsed"]] - started, 1)
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    completeDIAProteinAction(workflow, value)
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    stopifnot(identical(manager$current_state, "protein_replicate_filtered"))
    actual <- h$kitFingerprintState(manager$.context$store, getDIAStateReference(manager), api)
    mismatches <- h$kitOracleStateMismatches(actual, oracle)
    list(variant = variant, compute_seconds = compute_seconds, seconds = seconds, rows = actual$rows,
        proteins = value$protein_count, execution = value$execution, mismatches = mismatches,
        products_reported = names(actual$products), exact = !length(mismatches))
}
whole <- run("whole")
cat(sprintf("whole: %.1f s, rows=%d exact=%s %s\n", whole$seconds, whole$rows, whole$exact, paste(whole$mismatches, collapse = ",")))
rows <- run("rows")
cat(sprintf("rows: %.1f s, rows=%d pages=%s exact=%s %s\n", rows$seconds, rows$rows, rows$execution$batches, rows$exact, paste(rows$mismatches, collapse = ",")))
record <- list(schema = "kit.in_process_replicate_check.v1", ticket = "KIT-010", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    population = list(name = "journey_duplicates", state = input$state, rows = input$rows, entity_id = input$entity_id),
    oracle = list(name = "journey_protein_replicate", rows = oracle$rows), parameters = parameters,
    whole = whole, rows = rows, both_exact = isTRUE(whole$exact) && isTRUE(rows$exact))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-010/in-process-replicate-check.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", if (record$both_exact) "both variants reproduce the journey_protein_replicate oracle" else "MISMATCH", "\n")
