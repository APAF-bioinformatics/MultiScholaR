# KIT-014: the four peptide QC declarations checked against the KIT-005 oracles on the real 300-sample populations
# before the counted single members. Each step runs in process from the working tree on the KIT-005 control
# journey's state that fed it, materialised as a snapshot store by the harness so the evidence store is never
# written: the Q-value filter on the raw design state (`journey_design`), the protein identification evidence
# filter on `journey_qvalue`, the intensity filter on `journey_rollup` and the replicate filter on `journey_sample`,
# with the journey's recorded parameters and the registry's production variant; the state is committed through
# the step's finishing step and its fingerprint under the current canonical form is compared with the journey's
# next state as the member verifier compares it (components, support entries, args summary; products reported).
# Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-014/in-process-peptide-qc-check.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
steps <- list(
    list(name = "qvalue", operation = "peptide_qvalue_filter", population = "journey_design", oracle = "journey_qvalue",
        parameters = list(qvalue_threshold = 0.01, global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
            proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report")),
    list(name = "protein_peptide", operation = "peptide_protein_evidence_filter", population = "journey_qvalue",
        oracle = "journey_protein_peptide", parameters = list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 1)),
    list(name = "intensity", operation = "peptide_intensity_filter", population = "journey_rollup", oracle = "journey_intensity",
        parameters = list(strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1, peptides_intensity_cutoff_percentile = 0.1,
            grouping_variable = "group")),
    list(name = "replicate", operation = "peptide_replicate_filter", population = "journey_sample", oracle = "journey_replicate",
        parameters = list(replicate_group_column = "replicates")))
only <- commandArgs(TRUE)
if (length(only)) steps <- Filter(function(step) step$name %in% only, steps)
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit014-peptide-qc-check")
unlink(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
anon <- function() snapshotWorkflowMemory()$process$RssAnon
run <- function(step) {
    input <- b$populations[[step$population]]
    oracle <- b$oracles[[step$oracle]]
    stopifnot(is.null(h$kitPopulationInvalid(input)), !is.null(oracle$class))
    entry <- workflowActionSpec(step$operation)
    validateWorkflowActionParameters(entry$name, step$parameters, workflowExecutionPolicy(FALSE))
    snapshot <- h$kitPopulationStore(input, file.path(work, step$name, "population-store"), api)
    identity <- readSavedStoreIdentity(file.path(snapshot, "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(snapshot), identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    reference <- getDIAStateReference(manager)
    stopifnot(identical(manager$current_state, input$state), identical(reference$entity_id, input$entity_id))
    gc(full = TRUE)
    before <- anon()
    started <- proc.time()[["elapsed"]]
    value <- runProjected(manager, entry$declaration, reference, c(step$parameters, list(variant = entry$default_variant)))
    compute_seconds <- round(proc.time()[["elapsed"]] - started, 1)
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    if (identical(value$section, "srlQvalueProteotypicPeptideClean")) finishDIAReferenceQvalueStep(workflow, value) else
        finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
            value$section, value$counts)
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    after <- anon()
    stopifnot(identical(manager$current_state, value$state_name))
    actual <- h$kitFingerprintState(manager$.context$store, getDIAStateReference(manager), api)
    mismatches <- h$kitOracleStateMismatches(actual, oracle)
    list(step = step$name, operation = step$operation, population = step$population, oracle = step$oracle,
        variant = entry$default_variant, parameters = step$parameters, input_rows = input$rows, oracle_rows = oracle$rows,
        compute_seconds = compute_seconds, seconds = seconds, anon_before = before, anon_after = after,
        rows = actual$rows, batches = value$execution$batches, groups = value$execution$groups,
        engine_limit_bytes = value$execution$engine_limit_bytes, admission = value$execution$admission,
        mismatches = mismatches, products_reported = names(actual$products), exact = !length(mismatches))
}
results <- list()
for (step in steps) {
    result <- run(step)
    results[[step$name]] <- result
    cat(sprintf("%s: %.1f s (compute %.1f s), rows=%d exact=%s %s\n", step$name, result$seconds, result$compute_seconds,
        result$rows, result$exact, paste(result$mismatches, collapse = ",")))
    gc(full = TRUE)
}
record <- list(schema = "kit.in_process_peptide_qc_check.v1", ticket = "KIT-014", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), steps = unname(results),
    all_exact = all(vapply(results, `[[`, logical(1), "exact")))
out <- file.path(repo, "docs/qa/workflow-execution/kit/kit-014",
    if (length(only)) paste0("in-process-peptide-qc-check-", paste(only, collapse = "-"), ".json") else "in-process-peptide-qc-check.json")
jsonlite::write_json(record, out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", if (record$all_exact) "every step reproduces its journey oracle" else "MISMATCH", "\n")
