# KIT-028 replay check (in process, ceiling 15 min). The origin is the kit-journey-028 member's project (kit_027,
# screen corpus, import through DA). Its record is written and read back as JSON; each replay runs on its own
# reflinked copy of the project whose store is cut to the import and design states (kitPrepareReplayProject).
#   exact:  the record replayed; every state's fingerprint equal to the origin's (kitReproduceMismatches).
#   edited: the protein evidence filter's minimum peptides per protein raised 1 -> 16 (ticket rev 1.3; the fixture has at least 14 peptides per protein, so 16 removes the proteins with 14 or 15: the DIA-NN
#           q-value thresholds are fixed at 0.01 by the science contract, so the Q-value edit is refused); expected
#           rows from a direct kernel call on the Q-value state first; equal before the filter, different at it
#           with the expected rows, different after it.
#   inputs: a FASTA that differs from the record's sha256 is refused before any step (history unchanged).
#   legacy: a record whose Q-value step is `legacy` stops with workflow_replay_unsupported naming it.
# Usage: R_LIBS=<kit_027 library> Rscript docs/qa/workflow-execution/kit/kit-028/replay-check.R <label>
suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
label <- commandArgs(TRUE)[[1L]]
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
api <- function(name) get(name, asNamespace("MultiScholaR"))
member <- file.path(b$root, "kit-journey-028", "candidate-1")
journey <- readRDS(file.path(member, "journey-result.rds"))
origin_project <- journey$project
store_relative <- substring(dirname(journey$store), nchar(origin_project) + 2L)
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit028-replay-", label))
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
copy_project <- function(name) {
    target <- file.path(work, name)
    dir.create(target)
    stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(origin_project, "."), target)) == 0L)
    target
}
project_paths <- function(project) {
    results <- file.path(project, "dia_baseline_smoke", "results", "proteomics")
    list(base_dir = file.path(project, "dia_baseline_smoke"), results_dir = results,
        source_dir = file.path(project, "dia_baseline_smoke", "scripts", "proteomics"),
        protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
        da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
        pathway_dir = file.path(results, "pathway_enrichment"), time_dir = file.path(results, "publication_graphs"),
        integration_dir = file.path(project, "dia_baseline_smoke", "integration"))
}
record <- list(schema = "kit.replay_check.v1", ticket = "KIT-028", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), origin = member)
# The record, from a read-only copy of the origin, written and read back.
origin_copy <- copy_project("origin")
origin_root <- file.path(origin_copy, store_relative)
identity <- readSavedStoreIdentity(file.path(origin_root, ".workflow-entities", "store.json"))
manager <- newDIAStateAdapter(newWorkflowEntityStore(origin_root, identity$project_id, identity$workflow_id), reopen = TRUE)
json <- file.path(work, "analysis_record.json")
writeWorkflowAnalysisRecord(manager, json)
edit_input <- manager$getState("qvalue_filtered")
closeDIAStateAdapter(manager)
file.copy(json, file.path(repo, "docs/qa/workflow-execution/kit/kit-028", paste0("analysis_record-origin-", label, ".json")), overwrite = TRUE)
analysis <- jsonlite::read_json(json, simplifyVector = FALSE)
states <- vapply(analysis$steps, function(s) s$state, "")
operations <- vapply(analysis$steps, function(s) s$operation, "")
preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
record$states <- data.frame(state = states, operation = operations,
    role = ifelse(states %in% preconditions, "precondition", "replayed"), stringsAsFactors = FALSE)
corpus <- b$corpora$screen$files
inputs <- c(diann_report = corpus[["report.tsv"]]$path, fasta = corpus[["proteins.fasta"]]$path)
run_replay <- function(name, edited_record) {
    project <- copy_project(name)
    root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
    scratch <- file.path(work, paste0(name, "-scratch")); dir.create(scratch)
    t0 <- Sys.time()
    result <- h$kitReplayJourney(edited_record, root, origin_root, inputs, project_paths(project), scratch, api)
    result$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    result
}
# 1. Exact replay.
exact <- run_replay("exact", analysis)
record$exact <- list(seconds = exact$seconds, states = exact$states, mismatches = exact$mismatches,
    replayed_steps = vapply(exact$per_step, function(x) x$state, ""))
cat("EXACT:", length(exact$states), "states,", length(exact$mismatches), "mismatches",
    if (length(exact$mismatches)) paste0("(", paste(exact$mismatches, collapse = ", "), ")") else "", "|", round(exact$seconds), "s\n")
# 2. Edited record: expected rows from a direct kernel call before replay.
qstep <- which(states == "protein_peptide_filtered")
recorded <- stats::setNames(lapply(analysis$steps[[qstep]]$parameters, function(p) decodeWorkflowRecordValue(p$value)),
    vapply(analysis$steps[[qstep]]$parameters, function(p) p$name, ""))
kernel <- tryCatch(filterMinNumPeptidesPerProtein(edit_input, num_peptides_per_protein_thresh = 16,
    num_peptidoforms_per_protein_thresh = recorded$num_peptidoforms_per_protein_thresh), error = function(e) e)
expected_rows <- if (inherits(kernel, "error")) NA_integer_ else nrow(kernel@peptide_data)
record$edited_parameter <- list(state = "protein_peptide_filtered", name = "num_peptides_per_protein_thresh", from = recorded$num_peptides_per_protein_thresh, to = 16)
record$edited_expected_rows <- expected_rows
record$edited_kernel_error <- if (inherits(kernel, "error")) conditionMessage(kernel) else NULL
edit_input <- kernel <- NULL
edited <- analysis
for (i in seq_along(edited$steps[[qstep]]$parameters)) if (identical(edited$steps[[qstep]]$parameters[[i]]$name, "num_peptides_per_protein_thresh"))
    edited$steps[[qstep]]$parameters[[i]]$value$value <- 16
changed <- run_replay("edited", edited)
mismatched_states <- unique(sub("/.*$", "", changed$mismatches))
before <- states[seq_len(qstep - 1L)]; after <- states[seq(qstep + 1L, length(states))]
record$edited <- list(seconds = changed$seconds, mismatched_states = mismatched_states,
    before_identical = !any(before %in% mismatched_states), at_differs = "protein_peptide_filtered" %in% mismatched_states,
    after_all_differ = all(after %in% mismatched_states),
    replay_rows = changed$replay$states$protein_peptide_filtered$rows, origin_rows = changed$origin$states$protein_peptide_filtered$rows,
    rows_equal_expected = identical(as.integer(changed$replay$states$protein_peptide_filtered$rows), as.integer(expected_rows)))
cat("EDITED: before identical", record$edited$before_identical, "| at differs", record$edited$at_differs,
    "| after all differ", record$edited$after_all_differ, "| rows", record$edited$replay_rows, "expected", expected_rows,
    "origin", record$edited$origin_rows, "|", round(changed$seconds), "s\n")
# 3. Input refusal before any step.
project <- copy_project("inputs")
root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
wrong <- file.path(work, "other.fasta"); writeLines(c(">sp|X|OTHER", "MKV"), wrong)
refusal <- tryCatch(replayWorkflowAnalysisRecord(analysis, root, c(inputs[["diann_report"]], fasta = wrong) |> stats::setNames(c("diann_report", "fasta")),
    project_paths(project), file.path(work, "inputs-scratch")), workflow_replay_inputs = function(e) e)
history_after <- unlist(h$kitReadHead(file.path(root, ".workflow-entities"), api)$history %||% NA)
record$inputs_refused <- inherits(refusal, "workflow_replay_inputs")
record$inputs_message <- if (inherits(refusal, "condition")) conditionMessage(refusal) else NA
cat("INPUTS refused:", record$inputs_refused, "|", record$inputs_message, "\n")
# 4. Legacy step refusal (the Q-value step marked legacy).
legacy <- analysis
legacy$steps[[which(states == "qvalue_filtered")]]$operation <- "legacy"
project <- copy_project("legacy")
root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
stopped <- tryCatch(replayWorkflowAnalysisRecord(legacy, root, inputs, project_paths(project), {d <- file.path(work, "legacy-scratch"); dir.create(d); d}),
    workflow_replay_unsupported = function(e) e)
record$legacy_refused <- inherits(stopped, "workflow_replay_unsupported")
record$legacy_message <- if (inherits(stopped, "condition")) conditionMessage(stopped) else NA
cat("LEGACY refused:", record$legacy_refused, "|", record$legacy_message, "\n")
record$total_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
record$within_ceiling <- record$total_seconds <= 900
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-028", paste0("replay-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("TOTAL", round(record$total_seconds), "s; within 15 min:", record$within_ceiling, "\n")
