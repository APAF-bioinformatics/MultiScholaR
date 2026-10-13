# KIT-015: the four protein QC declarations checked against the KIT-005 oracles on the real 300-sample populations
# before the counted single members. Each step runs in process from the working tree on the KIT-005 control
# journey's state that fed it, materialised as a snapshot store by the harness so the evidence store is never
# written: the IQ rollup on `journey_imputation`, the accession cleanup on `journey_iq` (with the journey's FASTA
# table, read from its project), the protein intensity filter on `journey_accession` and the duplicate removal on
# `journey_protein_intensity`, with the parameters the oracle states record and the registry's production
# variant; the state is committed through the protein finishing step and its fingerprint under the current
# canonical form is compared with the journey's next state as the member verifier compares it (components,
# support entries, args summary; products reported). Not a counted run; the record is written beside this script.
# Usage: R_LIBS=<proof>/core/browser-smoke.X9aj9P/library Rscript docs/qa/workflow-execution/kit/kit-015/in-process-protein-qc-check.R [step ...]
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
api <- function(name) get(name, asNamespace("MultiScholaR"))
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
journey <- dirname(dirname(dirname(dirname(dirname(dirname(b$populations$journey_iq$store))))))
fasta <- file.path(journey, "project/dia_baseline_smoke/scripts/proteomics/aa_seq_tbl_final.RDS")
stopifnot(file.exists(fasta))
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit015-protein-qc-check")
unlink(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
steps <- list(
    list(name = "iq", operation = "protein_iq_rollup", population = "journey_imputation", oracle = "journey_iq",
        parameters = list(peptide_qc_dir = file.path(work, "iq", "peptide"), protein_qc_dir = file.path(work, "iq", "protein"))),
    list(name = "accession", operation = "protein_accession_cleanup", population = "journey_iq", oracle = "journey_accession",
        parameters = list(delimiter = ";", aggregation_method = "mean", fasta_table = fasta)),
    list(name = "protein_intensity", operation = "protein_intensity_filter", population = "journey_accession",
        oracle = "journey_protein_intensity", parameters = list(strict_mode = FALSE, min_reps_per_group = 2L, min_groups = 2L,
            proteins_intensity_cutoff_percentile = 1L)),
    list(name = "duplicates", operation = "protein_duplicate_removal", population = "journey_protein_intensity",
        oracle = "journey_duplicates", parameters = list(aggregation_method = "mean")))
only <- commandArgs(TRUE)
if (length(only)) steps <- Filter(function(step) step$name %in% only, steps)
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
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
    seconds <- round(proc.time()[["elapsed"]] - started, 1)
    after <- anon()
    stopifnot(identical(manager$current_state, value$state_name))
    actual <- h$kitFingerprintState(manager$.context$store, getDIAStateReference(manager), api)
    mismatches <- h$kitOracleStateMismatches(actual, oracle)
    list(step = step$name, operation = step$operation, population = step$population, oracle = step$oracle,
        variant = entry$default_variant, parameters = step$parameters[setdiff(names(step$parameters), c("peptide_qc_dir", "protein_qc_dir"))],
        input_rows = input$rows, oracle_rows = oracle$rows, compute_seconds = compute_seconds, seconds = seconds,
        anon_before = before, anon_after = after, rows = actual$rows, protein_count = value$protein_count,
        execution = value$execution, mismatches = mismatches, products_reported = names(actual$products), exact = !length(mismatches))
}
results <- list()
for (step in steps) {
    result <- run(step)
    results[[step$name]] <- result
    cat(sprintf("%s: %.1f s (compute %.1f s), rows=%d exact=%s %s\n", step$name, result$seconds, result$compute_seconds,
        result$rows, result$exact, paste(result$mismatches, collapse = ",")))
    gc(full = TRUE)
}
record <- list(schema = "kit.in_process_protein_qc_check.v1", ticket = "KIT-015", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), steps = unname(results),
    all_exact = all(vapply(results, `[[`, logical(1), "exact")))
out <- file.path(repo, "docs/qa/workflow-execution/kit/kit-015",
    if (length(only)) paste0("in-process-protein-qc-check-", paste(only, collapse = "-"), ".json") else "in-process-protein-qc-check.json")
jsonlite::write_json(record, out, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", if (record$all_exact) "every step reproduces its journey oracle" else "MISMATCH", "\n")
