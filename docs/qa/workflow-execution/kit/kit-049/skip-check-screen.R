# KIT-049 skip check (in process). Origin: the kit-journey-028 member's project (screen corpus, import through DA,
# recorded provenance). Its record is written and read back; each replay runs on a reflinked copy whose store is cut
# to the import and design states (kitPrepareReplayProject, as KIT-028's check does).
#   default: the record replayed as recorded on this build; every state's fingerprint equal to the origin's, so the
#            default path (skip absent) is byte-identical to the build the origin ran on.
#   skip:    the imputation step's parameters gain `skip = TRUE` (user_set); the states before `imputed` equal the
#            origin's, `imputed` differs, and on the replay store `Peptide.Imputed` equals the raw quantity column on
#            every row, `is_imputed` is FALSE everywhere, the state's config says skipped, the record written from the
#            replay store carries `skip = TRUE` with `outcomes$skipped`, and the replay ran every later step through DA.
# Usage: KIT049_INSTALLED=1 R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-049/skip-check-screen.R <label>
if (identical(Sys.getenv("KIT049_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
started <- Sys.time()
label <- commandArgs(TRUE)[[1L]]
repo <- normalizePath(".")
`%||%` <- function(a, b) if (is.null(a)) b else a
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
member <- file.path(b$root, "kit-journey-028", "candidate-1")
journey <- readRDS(file.path(member, "journey-result.rds"))
origin_project <- journey$project
store_relative <- substring(dirname(journey$store), nchar(origin_project) + 2L)
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit049-", label))
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
copy_project <- function(name) {
    target <- file.path(work, name); dir.create(target)
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
out <- list(schema = "kit.skip_check.v1", ticket = "KIT-049", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), origin = member,
    library = list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_))
origin_copy <- copy_project("origin")
origin_root <- file.path(origin_copy, store_relative)
identity <- api("readSavedStoreIdentity")(file.path(origin_root, ".workflow-entities", "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(origin_root, identity$project_id, identity$workflow_id), reopen = TRUE)
json <- file.path(work, "analysis_record.json")
api("writeWorkflowAnalysisRecord")(manager, json)
api("closeDIAStateAdapter")(manager)
analysis <- jsonlite::read_json(json, simplifyVector = FALSE)
states <- vapply(analysis$steps, function(s) s$state, "")
operations <- vapply(analysis$steps, function(s) s$operation, "")
preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
corpus <- b$corpora$screen$files
inputs <- c(diann_report = corpus[["report.tsv"]]$path, fasta = corpus[["proteins.fasta"]]$path)
run_replay <- function(name, edited_record) {
    project <- copy_project(name)
    root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
    scratch <- file.path(work, paste0(name, "-scratch")); dir.create(scratch)
    t0 <- Sys.time()
    result <- h$kitReplayJourney(edited_record, root, origin_root, inputs, project_paths(project), scratch, api)
    result$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    result$root <- root
    result
}
# 1. The default path: the record as recorded, replayed on this build.
exact <- run_replay("default", analysis)
out$default <- list(seconds = exact$seconds, states = exact$states, mismatches = exact$mismatches,
    replayed = vapply(exact$per_step, function(x) paste(x$operation, x$state), ""))
cat("DEFAULT:", length(exact$states), "states,", length(exact$mismatches), "mismatches",
    if (length(exact$mismatches)) paste0("(", paste(exact$mismatches, collapse = ", "), ")") else "", "|", round(exact$seconds), "s\n")
# 2. The skip path: the imputation step's recorded parameters gain skip = TRUE.
istep <- which(states == "imputed")
stopifnot(length(istep) == 1L, identical(operations[[istep]], "peptide_imputation"))
edited <- analysis
edited$steps[[istep]]$parameters[[length(edited$steps[[istep]]$parameters) + 1L]] <- list(name = "skip", source = "provenance", provenance = "user_set",
    value = list(type = "logical", value = TRUE, length = 1L))
skip <- run_replay("skip", edited)
mismatched_states <- unique(sub("/.*$", "", skip$mismatches))
before <- states[seq_len(istep - 1L)]; after <- states[seq(istep + 1L, length(states))]
out$skip <- list(seconds = skip$seconds, mismatched_states = mismatched_states,
    before_identical = !any(before %in% mismatched_states), at_differs = "imputed" %in% mismatched_states,
    after_all_differ = all(after %in% mismatched_states),
    replayed = vapply(skip$per_step, function(x) paste(x$operation, x$state), ""),
    da_replayed = any(vapply(skip$per_step, function(x) identical(x$operation, "DA-contrast-output"), logical(1))),
    replay_rows = skip$replay$states$imputed$rows, origin_rows = skip$origin$states$imputed$rows)
# The skipped state on the replay store: the column identity, the flags, the config, the record.
identity <- api("readSavedStoreIdentity")(file.path(skip$root, ".workflow-entities", "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(skip$root, identity$project_id, identity$workflow_id), reopen = TRUE)
on.exit(api("closeDIAStateAdapter")(manager), add = TRUE)
api("withDIAStateScope")(manager, function() {
    head <- h$kitReadHead(file.path(skip$root, ".workflow-entities"), api)
    record <- api("readDIAStateMetadata")(manager, head$references$imputed)
    spec <- record$data
    raw <- spec$slots$raw_quantity_column
    imputed_column <- spec$slots$args$peptideMissingValueImputation$imputed_value_column %||% "Peptide.Imputed"
    table <- api("materializeDIATable")(manager$.context, spec$components$peptide_data, c(raw, imputed_column, "is_imputed"))
    out$skip_state <<- list(rows = nrow(table), raw_column = raw, imputed_column = imputed_column,
        imputed_equals_raw = isTRUE(all.equal(table[[imputed_column]], table[[raw]])) && identical(is.na(table[[imputed_column]]), is.na(table[[raw]])),
        raw_missing = sum(is.na(table[[raw]])), imputed_missing = sum(is.na(table[[imputed_column]])),
        is_imputed_any = any(table$is_imputed %in% TRUE),
        args_skip = isTRUE(spec$slots$args$peptideMissingValueImputation$skip),
        summary_status = spec$slots$args$peptideMissingValueImputation$imputation_summary$status,
        summary_reason = spec$slots$args$peptideMissingValueImputation$imputation_summary$skip_reason,
        config_skip = isTRUE(record$config$skip), config_skipped = isTRUE(record$config$skipped),
        description = record$description, head_state = head$current_state, history = head$history)
    replay_json <- file.path(work, "analysis_record-skip.json")
    api("writeWorkflowAnalysisRecord")(manager, replay_json)
    replayed_record <- jsonlite::read_json(replay_json, simplifyVector = FALSE)
    step <- Filter(function(s) identical(s$state, "imputed"), replayed_record$steps)[[1L]]
    skip_parameter <- Filter(function(p) identical(p$name, "skip"), step$parameters)
    out$skip_record <<- list(parameters = vapply(step$parameters, function(p) paste(p$name, p$provenance), ""),
        skip_present = length(skip_parameter) == 1L,
        skip_value = if (length(skip_parameter)) api("decodeWorkflowRecordValue")(skip_parameter[[1L]]$value) else NA,
        skip_provenance = if (length(skip_parameter)) skip_parameter[[1L]]$provenance else NA_character_,
        outcomes = step$outcomes, da_present = !is.null(replayed_record$analyses$differential_abundance),
        valid = api("validateWorkflowAnalysisRecord")(replayed_record)$valid)
    file.copy(replay_json, file.path(repo, "docs/qa/workflow-execution/kit/kit-049", paste0("analysis_record-skip-", label, ".json")), overwrite = TRUE)
})
cat("SKIP: before identical", out$skip$before_identical, "| at differs", out$skip$at_differs, "| after all differ", out$skip$after_all_differ,
    "| DA replayed", out$skip$da_replayed, "| rows", out$skip$replay_rows, "origin", out$skip$origin_rows, "|", round(skip$seconds), "s\n")
cat("SKIP STATE: imputed==raw", out$skip_state$imputed_equals_raw, "| missing raw/imputed", out$skip_state$raw_missing, "/", out$skip_state$imputed_missing,
    "| any is_imputed", out$skip_state$is_imputed_any, "| args skip", out$skip_state$args_skip, "| summary", out$skip_state$summary_status, out$skip_state$summary_reason,
    "| config skip/skipped", out$skip_state$config_skip, out$skip_state$config_skipped, "\n")
cat("SKIP RECORD: skip present", out$skip_record$skip_present, "value", out$skip_record$skip_value, "provenance", out$skip_record$skip_provenance,
    "| outcomes", paste(names(out$skip_record$outcomes), unlist(out$skip_record$outcomes), collapse = " "), "| DA", out$skip_record$da_present, "| valid", out$skip_record$valid, "\n")
out$total_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-049", paste0("skip-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
cat("TOTAL", round(out$total_seconds), "s\n")
