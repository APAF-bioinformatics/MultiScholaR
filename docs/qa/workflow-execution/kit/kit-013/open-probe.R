# KIT-013: the CORE-014 read-only store probe before and after two Opens on the screen fixture.
# The session is the KIT-011 screen lifecycle member's project (kit_011 build, on-disk default; 18 states,
# session archive exported by the journey) copied with reflinks into a scratch directory, so the evidence
# project is never written. Each Open goes through the DA tab's native-session path (openDIANativeSession
# + adoptDIANativeSession); the store probe classifies the live store's entities against the closure Save
# archives, as core/core-014/store-closure-probe.R does, and the folder listing shows what an Open leaves
# beside the project (hand-off item 5). Not a counted run; the record is written beside this script.
# Usage: Rscript docs/qa/workflow-execution/kit/kit-013/open-probe.R
suppressMessages(devtools::load_all(quiet = TRUE))
repo <- normalizePath(".")
h <- new.env(parent = globalenv())
sys.source(file.path(repo, "tools/profiling/workflow_kit.R"), h)
b <- h$kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
member <- file.path(b$root, "kit-lifecycle-complete", "member")
stopifnot(identical(readRDS(file.path(member, "journey-result.rds"))$status, "complete"))
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit013-open-probe")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)
status <- system2("cp", c("-a", "--reflink=auto", file.path(member, "project", "dia_baseline_smoke"), work))
stopifnot(status == 0L)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
project <- file.path(work, "dia_baseline_smoke")
# The member's experiment paths (setupDirectories layout; setupDirectories itself prompts on an existing project).
results <- file.path(project, "results", "proteomics")
paths <- list(base_dir = project, results_dir = results, source_dir = file.path(project, "scripts", "proteomics"),
    protein_qc_dir = file.path(results, "protein_qc"), peptide_qc_dir = file.path(results, "peptide_qc"),
    da_output_dir = file.path(results, "da_proteins"), publication_graphs_dir = file.path(results, "publication_graphs"),
    time_dir = file.path(results, "publication_graphs"), integration_dir = file.path(project, "integration"))
stopifnot(all(dir.exists(unlist(paths))))
source_dir <- normalizePath(paths$source_dir, mustWork = TRUE)
archive <- list.files(source_dir, "[.]msdia[.]tar$", full.names = TRUE)
stopifnot(length(archive) == 1L)

private_folders <- function() {
    folders <- list.files(source_dir, "^[.]dia-", all.files = TRUE, include.dirs = TRUE)
    bytes <- vapply(folders, function(folder) sum(file.size(list.files(file.path(source_dir, folder),
        recursive = TRUE, full.names = TRUE, all.files = TRUE))), numeric(1))
    list(names = folders, bytes = unname(bytes))
}
store_probe <- function(manager) {
    # CORE-014 classification on the live owner's store: entities referenced from the head/snapshot
    # closure (what Save archives) against the rest.
    ctx <- manager$.context
    store <- ctx$store
    refs <- expandWorkflowTableReferences(store, collectDIANativeReferences(list(ctx$snapshot, ctx$head)))
    live <- unique(vapply(refs, function(r) basename(dirname(r$payload$relative_path)), character(1)))
    ents <- list.dirs(file.path(store$directory, "entities"), recursive = FALSE)
    size <- vapply(ents, function(d) sum(file.size(list.files(d, recursive = TRUE, full.names = TRUE))), numeric(1))
    dead <- !basename(ents) %in% live
    list(store = basename(dirname(store$directory)), entities = length(ents), closure = sum(!dead),
        dead = sum(dead), bytes = sum(size), dead_bytes = sum(size[dead]),
        current_state = ctx$snapshot$current_state, states = length(ctx$snapshot$states))
}
open_store <- function() {
    roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
    stopifnot(length(roots) == 1L)
    identity <- readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id,
        max_cache_bytes = 0), reopen = TRUE)
    manager$.context$session_format <- "native_v1"
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    manager
}

record <- list(schema = "kit.open_probe.v1", ticket = "KIT-013", date = format(Sys.Date()),
    head_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    session = list(source = "kit-lifecycle-complete member project (kit_011 screen journey)",
        build = b$lifecycle$build, archive = basename(archive)),
    steps = list())
step <- function(label, manager, seconds = NA_real_) {
    folders <- private_folders()
    entry <- list(label = label, seconds = seconds, private_folders = folders$names,
        private_folder_bytes = folders$bytes, live_store = store_probe(manager))
    record$steps[[length(record$steps) + 1L]] <<- entry
    cat(sprintf("%s: folders=%s live=%s entities=%d closure=%d dead=%d dead_bytes=%.1fMB state=%s\n", label,
        paste(folders$names, collapse = ","), entry$live_store$store, entry$live_store$entities,
        entry$live_store$closure, entry$live_store$dead, entry$live_store$dead_bytes / 1e6, entry$live_store$current_state))
    invisible(entry)
}

manager <- open_store()
workflow <- shiny::reactiveValues(state_manager = manager, data_type = "peptide", data_format = "diann")
shiny::testServer(function(input, output, session) {}, {
    workflow$state_manager <- attachDIAStateSession(manager, session, workflow)
    step("before", workflow$state_manager)
    for (i in 1:2) {
        started <- proc.time()[["elapsed"]]
        imported <- openDIANativeSession(archive, source_dir)
        view <- adoptDIANativeSession(imported, workflow, session, source_dir, paths)
        stopifnot(!is.null(view$current_s4_object) || !is.null(view$final_state_name) || is.list(view))
        step(paste0("after_open_", i), workflow$state_manager, round(proc.time()[["elapsed"]] - started, 1))
    }
    # Audit PX-03: a Save from the opened session after the replaced stores are gone, and a reopen of it.
    manager <- workflow$state_manager
    ctx <- manager$.context
    metadata <- resolveDIAStateValue(ctx, ctx$native_metadata_reference)
    # The GUI rebuilds metadata from the bus (collectDIANativeSessionMetadata); the stored view still
    # carries the legacy history aliases the exporter refuses, so they are dropped here.
    metadata[intersect(names(metadata), c("r6_complete_states", "current_s4_object",
        "correlation_filtered_s4", "state_manager"))] <- NULL
    export_dir <- file.path(work, "export-after-open")
    dir.create(export_dir)
    started <- proc.time()[["elapsed"]]
    exported <- runDIANativeSessionExport(manager, metadata, export_dir,
        extra_references = ctx$native_extra_references, dependencies = resolveDIANativeDependencies(manager))
    reopen_dir <- file.path(work, "reopen-after-export")
    dir.create(reopen_dir)
    reopened <- openDIANativeSession(exported$publication$archive_path, reopen_dir)
    record$save_after_open <<- list(seconds = round(proc.time()[["elapsed"]] - started, 1),
        archive = basename(exported$publication$archive_path),
        archive_bytes = file.size(exported$publication$archive_path),
        reopened_state = reopened$manager$current_state,
        reopened_states = length(reopened$manager$.context$snapshot$states),
        states_identical = identical(names(reopened$manager$.context$snapshot$states), names(ctx$snapshot$states)))
    closeDIAStateAdapter(reopened$manager)
    cat(sprintf("save_after_open: %s (%.1f MB) reopened at %s with %d states\n", record$save_after_open$archive,
        record$save_after_open$archive_bytes / 1e6, record$save_after_open$reopened_state, record$save_after_open$reopened_states))
    closeDIAStateAdapter(workflow$state_manager)
})
folders <- vapply(record$steps, function(s) paste(s$private_folders, collapse = ","), character(1))
record$result <- list(
    before_folders = record$steps[[1L]]$private_folders,
    after_first_open = record$steps[[2L]]$private_folders,
    after_second_open = record$steps[[3L]]$private_folders,
    replaced_folder_removed_each_time = length(record$steps[[2L]]$private_folders) == 1L &&
        length(record$steps[[3L]]$private_folders) == 1L &&
        !identical(record$steps[[2L]]$private_folders, record$steps[[3L]]$private_folders) &&
        !record$steps[[1L]]$private_folders %in% record$steps[[2L]]$private_folders,
    removed_bytes = c(record$steps[[1L]]$private_folder_bytes, record$steps[[2L]]$private_folder_bytes))
jsonlite::write_json(record, file.path(repo, "docs/qa/workflow-execution/kit/kit-013/open-probe.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
cat("RESULT", if (isTRUE(record$result$replaced_folder_removed_each_time)) "replaced folder removed after each Open" else
    "orphaned folder remains", "\n")
