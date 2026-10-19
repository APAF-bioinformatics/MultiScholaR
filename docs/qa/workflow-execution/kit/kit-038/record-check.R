# KIT-038 record check (in process). The origin is the kit-journey-028 member's project (kit_027, screen corpus,
# import through DA), opened read-only from a reflinked copy:
#   rebuilt: the record built by this tree, diffed against kit-028/analysis_record-origin-screen.json (fixed rows
#            appear, no @args dump on the IQ/P14/design steps, no `no_default`, outcomes and build present, the
#            environment fields "not recorded" because that store's environment entry predates KIT-038);
#   valid:   the origin record (16 `no_default`) and the rebuilt record both validate; a bogus string does not;
#   replay:  the rebuilt record replays exactly (kitPrepareReplayProject + kitReplayJourney, as KIT-028);
#   env:     the record of the replayed store (states written by this build) carries the environment fields;
#   study:   study_parameters.txt written from the rebuilt record: no truncation marks; diff against the
#            origin project's file (written by the KIT-033 build) confined to the Analysis Record section and the
#            Timestamp line.
# Usage: R_LIBS=<proof library> Rscript docs/qa/workflow-execution/kit/kit-038/record-check.R <label>
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
work <- file.path(Sys.getenv("HOME"), ".cache", paste0("multischolar-kit038-", label))
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
out <- list(schema = "kit.record_check.v1", ticket = "KIT-038", label = label, date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), origin = member)
open_store <- function(root) {
    identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
    newDIAStateAdapter(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
}
step_summary <- function(record) lapply(record$steps, function(step) {
    provenance <- vapply(step$parameters, function(p) p$provenance, "")
    list(state = step$state, operation = step$operation, parameters = length(step$parameters),
        fixed = sum(provenance == "fixed"), unknown = sum(provenance == "unknown"),
        no_default = sum(provenance == "no_default"), not_set = sum(provenance == "not_set"),
        no_shipped_default = sum(provenance == "no_shipped_default"),
        names = vapply(step$parameters, function(p) p$name, ""), outcomes = names(step$outcomes),
        build = step$software$build %||% NA_character_)
})
# 1. The rebuilt record against the origin record.
origin_copy <- copy_project("origin")
origin_root <- file.path(origin_copy, store_relative)
manager <- open_store(origin_root)
json <- file.path(work, "analysis_record.json")
writeWorkflowAnalysisRecord(manager, json)
head_build <- manager$.context$snapshot$step_provenance[[manager$current_state]]$build
paths <- project_paths(origin_copy)
baseline_study <- file.path(paths$source_dir, "study_parameters.txt")
baseline_lines <- if (file.exists(baseline_study)) readLines(baseline_study) else NULL
rebuilt_for_study <- jsonlite::read_json(json, simplifyVector = FALSE)
contrasts <- do.call(rbind, lapply(rebuilt_for_study$analyses$differential_abundance$contrasts, function(row) as.data.frame(lapply(row, unlist))))
invisible(capture.output(createWorkflowArgsFromConfig("dia_baseline_smoke", source_dir_path = paths$source_dir,
    final_s4_object = getDIAStateReference(manager), contrasts_tbl = contrasts,
    workflow_data = list(state_manager = manager, config_list = list()))))
study <- readLines(baseline_study)
closeDIAStateAdapter(manager)
rebuilt <- jsonlite::read_json(json, simplifyVector = FALSE)
origin <- jsonlite::read_json(file.path(repo, "docs/qa/workflow-execution/kit/kit-028/analysis_record-origin-screen.json"), simplifyVector = FALSE)
file.copy(json, file.path(repo, "docs/qa/workflow-execution/kit/kit-038", paste0("analysis_record-rebuilt-", label, ".json")), overwrite = TRUE)
out$origin_steps <- step_summary(origin)
out$rebuilt_steps <- step_summary(rebuilt)
by_state <- function(record) stats::setNames(record$steps, vapply(record$steps, function(s) s$state, ""))
rs <- by_state(rebuilt); os <- by_state(origin)
named <- function(step) vapply(step$parameters, function(p) p$name, "")
out$iq <- list(origin_parameters = named(os$protein_s4_created), rebuilt_parameters = named(rs$protein_s4_created),
    rebuilt_fixed = vapply(Filter(function(p) identical(p$provenance, "fixed"), rs$protein_s4_created$parameters), function(p) p$name, ""))
p14_state <- names(rs)[vapply(rs, function(s) identical(s$operation, "precursor_rollup"), logical(1))]
out$p14 <- list(state = p14_state, parameters = lapply(rs[[p14_state]]$parameters, function(p) list(name = p$name, provenance = p$provenance, value = p$value$value)))
out$design <- list(origin_parameters = length(os$raw_data_s4$parameters), rebuilt_parameters = named(rs$raw_data_s4))
all_prov <- function(record) unlist(lapply(record$steps, function(s) vapply(s$parameters, function(p) p$provenance, "")))
out$provenance_counts <- list(origin = as.list(table(all_prov(origin))), rebuilt = as.list(table(all_prov(rebuilt))))
out$fixed_rows_rebuilt <- sum(all_prov(rebuilt) == "fixed")
out$no_default_rebuilt <- sum(all_prov(rebuilt) == "no_default")
ruv_state <- names(rs)[vapply(rs, function(s) identical(s$operation, "normalization-ruv"), logical(1))]
out$ruv <- list(state = ruv_state, outcomes = rs[[ruv_state]]$outcomes,
    requested_k = Filter(function(p) identical(p$name, "ruv_k"), rs[[ruv_state]]$parameters)[[1]][c("source", "provenance")],
    ruv_mode_provenance = Filter(function(p) identical(p$name, "ruv_mode"), rs[[ruv_state]]$parameters)[[1]]$provenance)
out$software <- rebuilt$software[c("version", "sha", "build", "blas", "lapack", "os_release", "locale", "duckdb_threads")]
out$build_matches_head <- identical(rebuilt$software$build, head_build)
out$step_builds_present <- all(vapply(rebuilt$steps, function(s) is.character(s$software$build) && nzchar(s$software$build), logical(1)))
out$outcomes_present <- all(vapply(rebuilt$steps, function(s) is.list(s$outcomes) && !is.null(s$outcomes$rows), logical(1)))
# 2. Validation: origin (legacy spelling) and rebuilt valid; a bogus string invalid.
out$validation <- list(origin = validateWorkflowAnalysisRecord(origin)[c("valid", "errors")],
    rebuilt = validateWorkflowAnalysisRecord(rebuilt)[c("valid", "errors")])
bogus <- rebuilt; bogus$steps[[3]]$parameters[[1]]$provenance <- "bogus"
out$validation$bogus <- validateWorkflowAnalysisRecord(bogus)[c("valid", "errors")]
# 3. study_parameters.txt from the rebuilt record.
out$study <- list(lines = length(study), truncation_marks = sum(grepl("omitted|\\[list object\\]|\\.\\.\\.", study)),
    record_section = any(grepl("^Analysis Record", study)), baseline_present = !is.null(baseline_lines))
if (!is.null(baseline_lines)) {
    strip <- function(x) { x <- x[!grepl("^Timestamp:", x)]; i <- grep("^Analysis Record", x); if (length(i)) x[seq_len(i[[1]] - 1L)] else x }
    out$study$outside_record_section_identical <- identical(strip(baseline_lines), strip(study))
    out$study$baseline_lines <- length(baseline_lines)
}
# 4. Exact replay of the rebuilt record, then the replayed store's own record (written by this build).
states <- vapply(rebuilt$steps, function(s) s$state, "")
operations <- vapply(rebuilt$steps, function(s) s$operation, "")
preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
corpus <- b$corpora$screen$files
inputs <- c(diann_report = corpus[["report.tsv"]]$path, fasta = corpus[["proteins.fasta"]]$path)
project <- copy_project("exact")
root <- h$kitPrepareReplayProject(file.path(project, store_relative), preconditions, api)
scratch <- file.path(work, "exact-scratch"); dir.create(scratch)
t0 <- Sys.time()
exact <- h$kitReplayJourney(rebuilt, root, origin_root, inputs, project_paths(project), scratch, api)
out$replay <- list(seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")), states = length(exact$states), mismatches = exact$mismatches)
cat("REPLAY:", length(exact$states), "states,", length(exact$mismatches), "mismatches |", round(out$replay$seconds), "s\n")
replayed <- open_store(root)
env_record <- buildWorkflowAnalysisRecord(replayed)
closeDIAStateAdapter(replayed)
out$replayed_store_software <- env_record$software[c("build", "blas", "lapack", "omp_threads", "blas_threads", "os_release", "locale")]
out$replayed_store_env_recorded <- !any(vapply(env_record$software[c("blas", "lapack", "os_release", "locale")],
    function(v) identical(v, "not recorded"), logical(1)))
out$seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
jsonlite::write_json(out, file.path(repo, "docs/qa/workflow-execution/kit/kit-038", paste0("record-check-", label, ".json")),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("IQ rebuilt parameters:", paste(out$iq$rebuilt_parameters, collapse = ","), "| P14:", paste(vapply(out$p14$parameters, function(p) paste0(p$name, "=", p$value, "[", p$provenance, "]"), ""), collapse = ","),
    "| design:", paste(out$design$rebuilt_parameters, collapse = ","), "\n")
cat("fixed rows", out$fixed_rows_rebuilt, "| no_default", out$no_default_rebuilt, "| ruv outcomes", paste(names(out$ruv$outcomes), unlist(lapply(out$ruv$outcomes, function(v) if (is.null(v)) "NULL" else as.character(v))), collapse = " "), "| requested ruv_k", out$ruv$requested_k$provenance, "| ruv_mode", out$ruv$ruv_mode_provenance, "\n")
cat("valid origin", out$validation$origin$valid, "rebuilt", out$validation$rebuilt$valid, "bogus", out$validation$bogus$valid, "| build matches head", out$build_matches_head,
    "| blas", out$software$blas, "| replayed-store env recorded", out$replayed_store_env_recorded, "\n")
cat("STUDY lines", out$study$lines, "truncation", out$study$truncation_marks, "record section", out$study$record_section, "baseline", out$study$baseline_present,
    "outside-section identical", out$study$outside_record_section_identical, "|", round(out$seconds), "s total\n")
