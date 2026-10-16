# KIT-027 in-process check on the screen fixture (the KIT-026 check extended) (the KIT-008 screen member's project, reflinked): every state
# fingerprinted before the store is reopened on the inputs-recording build; the screen corpus report and FASTA
# recorded by role (the function the import's completion calls); the IQ rollup rerun in process from `imputed`
# (its commit writes the `iq_output` and `peptide_input` roles) and the accession cleanup rerun with the bound
# state's parameters (its config carries `fasta_sha256` from the head). Record invariance: every state fingerprint
# identical before and after, the two recomputed states included (their tables compared, `fasta_sha256` checked
# against the head and the corpus hash).
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-026/inputs-check.R <label>
suppressMessages(devtools::load_all(quiet = TRUE))
main <- function(label) {
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
corpus <- b$corpora$screen$files
member <- file.path(b$root, "kit-gui-008", "member")
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit027-record-", label))
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
roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
stopifnot(length(roots) == 1L)
api <- function(name) get(name, asNamespace("MultiScholaR"))
store_dir <- file.path(roots, ".workflow-entities")
before <- h$kitFingerprintStore(store_dir, api)
identity <- readSavedStoreIdentity(file.path(store_dir, "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
on.exit(closeDIAStateAdapter(manager), add = TRUE, after = FALSE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
record <- list(schema = "kit.record_check.v1", ticket = "KIT-027", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
# Inputs recorded before roles existed open with role NA.
record$reopened_roles <- vapply(manager$.context$snapshot$file_dependencies, function(d) as.character(d$role), character(1))
accession_before <- readDIAStateMetadata(manager, getDIAStateReference(manager, "protein_accession_cleaned"))
withDIAStateScope(manager, function()
    recordDIAWorkflowInputs(manager, c(diann_report = corpus[["report.tsv"]]$path, fasta = corpus[["proteins.fasta"]]$path)))
# Every state after the inputs are recorded (the rerun below truncates the history after `imputed`).
closeDIAStateAdapter(manager)
recorded <- h$kitFingerprintStore(store_dir, api)
manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"
manager$.context$reference_states <- TRUE
manager$.context$experiment_paths <- paths
withDIAStateScope(manager, function() {
    selectDIAStateReference(manager, "imputed")
    workflow <- list2env(list(state_manager = manager), parent = emptyenv())
    runProteinIqRollupApplyStep(workflow, paths, captureCheckpointFn = function(...) NULL, showNotificationFn = function(...) NULL)
    value <- computeDIAReferenceAccessionCleanup(manager, getDIAStateReference(manager),
        accession_before$config[c("delimiter", "aggregation_method")])
    owner <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    completeDIAProteinAction(owner, value)
})
# KIT-027: the analysis record built twice from the store, validated, written, and every atomic argument of every
# state restored identical() from the JSON; the full-object read counter unchanged by the build.
reads_before <- manager$.context$full_s4_reads
hydrates_before <- manager$.context$hydrates
json_a <- file.path(work, "analysis_record-a.json"); json_b <- file.path(work, "analysis_record-b.json")
fixed_time <- as.POSIXct("2026-10-17 10:00:00", tz = "UTC")
writeWorkflowAnalysisRecord(manager, json_a, now = fixed_time)
writeWorkflowAnalysisRecord(manager, json_b, now = fixed_time)
record$record_full_s4_reads_unchanged <- identical(manager$.context$full_s4_reads, reads_before)
record$record_hydrates_delta <- manager$.context$hydrates - hydrates_before
record$record_byte_identical <- identical(readLines(json_a), readLines(json_b))
parsed <- jsonlite::read_json(json_a, simplifyVector = FALSE)
validation <- validateWorkflowAnalysisRecord(parsed)
record$record_valid <- validation$valid
record$record_errors <- validation$errors
record$record_incomplete <- validation$incomplete
history_now <- unlist(manager$.context$snapshot$state_history)
record$record_steps_equal_history <- identical(vapply(parsed$steps, function(x) x$state, ""), history_now)
record$record_operations <- vapply(parsed$steps, function(x) x$operation, "")
names(record$record_operations) <- history_now
mismatched <- character(); compared <- 0L
compare <- function(original, encoded, path) {
    if (is.data.frame(original) || is.null(original) || is.function(original) || is.environment(original)) return(invisible())
    if (is.list(original)) { for (name in names(original)) compare(original[[name]], encoded$value[[match(name, unlist(encoded$names))]], c(path, name)); return(invisible()) }
    compared <<- compared + 1L
    if (!identical(decodeWorkflowRecordValue(encoded), original)) mismatched <<- c(mismatched, paste(path, collapse = "/"))
}
for (step in parsed$steps) {
    meta <- tryCatch(readDIAStateMetadata(manager, getDIAStateReference(manager, step$state)), error = function(e) NULL)
    args <- meta$data$slots$args
    if (is.null(args)) next
    by_name <- stats::setNames(step$parameters, vapply(step$parameters, function(p) p$name, ""))
    if (identical(step$operation, "legacy")) {
        for (name in names(by_name)) compare(args[[name]], by_name[[name]]$value, c(step$state, name))
    } else for (name in names(by_name)) if (identical(by_name[[name]]$source, "args"))
        compare(args[[names(workflowSectionOperations)[match(step$operation, workflowSectionOperations)]]][[name]], by_name[[name]]$value, c(step$state, name))
}
record$round_trip_values_compared <- compared
record$round_trip_mismatches <- mismatched
accession_step <- Filter(function(x) identical(x$state, "protein_accession_cleaned"), parsed$steps)[[1]]
fasta_parameter <- Filter(function(p) identical(p$name, "fasta_sha256"), accession_step$parameters)
fasta_input <- Filter(function(x) identical(x$role, "fasta"), parsed$inputs)
record$record_fasta <- list(source = if (length(fasta_parameter)) fasta_parameter[[1]]$source else NA,
    value_equals_input = length(fasta_parameter) && length(fasta_input) &&
        identical(decodeWorkflowRecordValue(fasta_parameter[[1]]$value), fasta_input[[1]]$sha256))
record$record_inputs <- vapply(parsed$inputs, function(x) as.character(x$role %||% NA), "")
record$record_external_resources <- length(parsed$external_resources)
record$record_limits <- parsed$limits
file.copy(json_a, file.path(repo, "docs/qa/workflow-execution/kit/kit-027", paste0("analysis_record-", label, ".json")), overwrite = TRUE)
snapshot <- manager$.context$snapshot
inputs <- snapshot$file_dependencies
record$roles <- vapply(inputs, function(d) as.character(d$role), character(1))
record$corpus_inputs <- lapply(c(diann_report = "report.tsv", fasta = "proteins.fasta"), function(file) {
    role <- if (file == "report.tsv") "diann_report" else "fasta"
    list(recorded = inputs[[role]]$sha256, corpus_binding = corpus[[file]]$sha256,
        sha256sum = substr(system2("sha256sum", corpus[[file]]$path, stdout = TRUE), 1L, 64L),
        bytes_equal = isTRUE(inputs[[role]]$bytes == corpus[[file]]$bytes))
})
record$execution_inputs <- lapply(c("iq_output", "peptide_input"), function(role) {
    entry <- inputs[[paste0(role, "_file")]]
    list(role = entry$role, recorded = entry$sha256, sha256sum = substr(system2("sha256sum", entry$source, stdout = TRUE), 1L, 64L))
})
names(record$execution_inputs) <- c("iq_output", "peptide_input")
accession <- readDIAStateMetadata(manager, getDIAStateReference(manager, "protein_accession_cleaned"))
record$accession <- list(fasta_sha256 = accession$config$fasta_sha256, head_fasta = inputs$fasta$sha256,
    cleanup_applied = accession$config$cleanup_applied, slots_args_has_fasta = "fasta_sha256" %in%
        names(unlist(accession$data$slots$args, recursive = FALSE)))
entry <- snapshot$step_provenance$protein_accession_cleaned
record$accession_provenance <- list(operation = entry$operation, parameter_names = entry$parameter_names,
    names_in_config = entry$parameter_names %in% names(accession$config))
closeDIAStateAdapter(manager)
after <- h$kitFingerprintStore(store_dir, api)
compare <- function(x, y) vapply(names(y$states), function(name) {
    a <- x$states[[name]]; b <- y$states[[name]]
    !is.null(a) && identical(a$components, b$components) && identical(a$support, b$support) && identical(a$summary_sha256, b$summary_sha256)
}, logical(1))
same <- c(compare(before, recorded), compare(before, after)[c("protein_s4_created", "protein_accession_cleaned")])
record$states_compared <- names(same)
record$state_fingerprints_identical <- all(same)
record$states_differing <- names(same)[!same]
record$recomputed_states_identical <- same[c("protein_s4_created", "protein_accession_cleaned")]
record$fasta_matches_input <- after$states$protein_accession_cleaned$fasta_matches_input
oracle_like <- lapply(before$states, function(x) list(fingerprint = x))
after_like <- lapply(after$states, function(x) list(fingerprint = x))
record$verifier_mismatches <- as.vector(h$kitReproduceMismatches(oracle_like[c("protein_s4_created", "protein_accession_cleaned")],
    after_like[c("protein_s4_created", "protein_accession_cleaned")]))
record$head_sha256 <- list(before = before$head_sha256, inputs_recorded = recorded$head_sha256, after = after$head_sha256)
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-027", paste0("record-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("ROLES:", paste(names(record$roles), record$roles, sep = "=", collapse = " "), "\n")
cat("CORPUS:", paste(names(record$corpus_inputs), vapply(record$corpus_inputs, function(x) identical(x$recorded, x$sha256sum) &&
    identical(x$recorded, x$corpus_binding) && x$bytes_equal, logical(1)), sep = "=", collapse = " "), "\n")
cat("EXECUTION:", paste(names(record$execution_inputs), vapply(record$execution_inputs, function(x) identical(x$recorded, x$sha256sum), logical(1)), sep = "=", collapse = " "), "\n")
cat("ACCESSION fasta_sha256 equals head:", identical(record$accession$fasta_sha256, record$accession$head_fasta),
    "| applied:", record$accession$cleanup_applied, "| in slots args:", record$accession$slots_args_has_fasta,
    "| provenance names in config:", all(record$accession_provenance$names_in_config), "\n")
cat("RECORD valid:", record$record_valid, "| steps = history:", record$record_steps_equal_history, "| byte identical:", record$record_byte_identical,
    "| full_s4_reads unchanged:", record$record_full_s4_reads_unchanged, "| round trip:", record$round_trip_values_compared, "values,",
    length(record$round_trip_mismatches), "mismatches | fasta:", record$record_fasta$source, record$record_fasta$value_equals_input,
    "| inputs:", paste(record$record_inputs, collapse = ","), "| incomplete:", length(record$record_incomplete), "\n")
cat("VERIFIER fasta_matches_input:", record$fasta_matches_input, "| mismatches:", paste(record$verifier_mismatches, collapse = ","), "\n")
cat("FINGERPRINTS identical:", record$state_fingerprints_identical, "(", length(same), "comparisons, the two recomputed states last; differing:",
    paste(record$states_differing, collapse = ","), ")\n")
}
main(commandArgs(TRUE)[[1L]])
