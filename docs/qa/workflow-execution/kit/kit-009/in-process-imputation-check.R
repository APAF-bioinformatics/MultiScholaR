# KIT-009: the imputation step on the feature-group partitioner checked against the KIT-005 oracle on the
# real 300-sample population before the counted members. Input: the KIT-005 control journey's state after
# the replicate filter (`journey_replicate`, 5,515,154 rows), materialised as a snapshot store per variant
# by the harness so the evidence store is never written. Each variant (`grouped`: the per-group loop as it
# stood before KIT-009; `partitioned`: the partitioner) runs in process from the working tree with the
# journey's recorded parameters (proportion 0.5, no exclusion column), the state is committed through the
# shared finishing step, and its fingerprint under the current canonical form is compared with the
# `journey_imputation` oracle as the member verifier compares it (components, support entries, args
# summary; products reported). Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-009/in-process-imputation-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
input <- b$populations$journey_replicate
oracle <- b$oracles$journey_imputation
stopifnot(identical(input$state, "replicate_filtered"), is.null(h$kitPopulationInvalid(input)), !is.null(oracle$class))
parameters <- list(proportion_missing_values = 0.5, exclusion_column = NULL)
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit009-imputation-check")
unlink(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
anon <- function() snapshotWorkflowMemory()$process$RssAnon
run <- function(variant) {
    snapshot <- h$kitPopulationStore(input, file.path(work, variant, "population-store"), api)
    identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    reference <- getDIAStateReference(manager)
    stopifnot(identical(manager$current_state, "replicate_filtered"), identical(reference$entity_id, input$entity_id))
    gc(full = TRUE)
    before <- anon()
    started <- proc.time()[["elapsed"]]
    value <- computeDIAReferencePeptideGroups(manager, reference, "peptideMissingValueImputation", parameters,
        variant = variant)
    compute_seconds <- round(proc.time()[["elapsed"]] - started, 1)
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
        value$section, value$counts)
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    after <- anon()
    stopifnot(identical(manager$current_state, "imputed"))
    actual <- h$kitFingerprintState(manager$.context$store, getDIAStateReference(manager), api)
    mismatches <- h$kitOracleStateMismatches(actual, oracle)
    list(variant = variant, compute_seconds = compute_seconds, seconds = seconds,
        anon_before = before, anon_after = after,
        rows = actual$rows, batches = value$execution$batches, groups = value$execution$groups,
        engine_limit_bytes = value$execution$engine_limit_bytes, admission = value$execution$admission,
        summary = value$summary, mismatches = mismatches,
        products_reported = names(actual$products), exact = !length(mismatches))
}
grouped <- run("grouped")
cat(sprintf("grouped: %.1f s, rows=%d batches=%d groups=%d exact=%s %s\n", grouped$seconds, grouped$rows, grouped$batches,
    grouped$groups, grouped$exact, paste(grouped$mismatches, collapse = ",")))
partitioned <- run("partitioned")
cat(sprintf("partitioned: %.1f s, rows=%d batches=%d groups=%d exact=%s %s\n", partitioned$seconds, partitioned$rows,
    partitioned$batches, partitioned$groups, partitioned$exact, paste(partitioned$mismatches, collapse = ",")))
record <- list(schema = "kit.in_process_imputation_check.v1", ticket = "KIT-009", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    population = list(name = "journey_replicate", state = input$state, rows = input$rows, entity_id = input$entity_id),
    oracle = list(name = "journey_imputation", rows = oracle$rows), parameters = parameters,
    grouped = grouped, partitioned = partitioned,
    same_batches = identical(grouped$batches, partitioned$batches),
    both_exact = isTRUE(grouped$exact) && isTRUE(partitioned$exact))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-009/in-process-imputation-check.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", if (record$both_exact) "both variants reproduce the journey_imputation oracle" else "MISMATCH", "\n")
