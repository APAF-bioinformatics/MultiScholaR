# KIT-031 acceptance: the analysis record of the GUI member's store (kit-journey-031l, kit_031) validates with no legacy
# step and no missing parameter, and replays exact: a copy of the member's project, its store cut to the import and
# design states (kitPrepareReplayProject), the record replayed in process (replayWorkflowAnalysisRecord through
# kitReplayJourney) and every state compared with the member's (kitReproduceMismatches), the DA re-run from the record.
# Inputs are hash-checked against the member's import copies and design file. Read-only on the member.
# Usage: KIT031_INSTALLED=1 R_LIBS=<kit_031 library> Rscript docs/qa/workflow-execution/kit/kit-031/member-record-replay.R <member dir> <out json>
if (identical(Sys.getenv("KIT031_INSTALLED"), "1")) suppressMessages(library(MultiScholaR)) else suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(TRUE); member <- args[[1L]]; out <- args[[2L]]
started <- Sys.time()
`%||%` <- function(a, b) if (is.null(a)) b else a
api <- function(name) get(name, asNamespace("MultiScholaR"))
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
project <- file.path(member, "project")
store_dir <- h$kitFindStore(project)
stopifnot(!is.null(store_dir))
store_root <- dirname(store_dir)
source_dir <- dirname(store_root)
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(store_root, identity$project_id, identity$workflow_id), reopen = TRUE)
record <- api("buildWorkflowAnalysisRecord")(manager)
api("closeDIAStateAdapter")(manager)
work <- file.path(Sys.getenv("HOME"), ".cache", "multischolar-kit031-member-replay")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
json <- file.path(work, "analysis_record.json")
jsonlite::write_json(record, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
read_back <- jsonlite::read_json(json, simplifyVector = FALSE)
check <- api("validateWorkflowAnalysisRecord")(read_back)
operations <- vapply(read_back$steps, function(s) s$operation, "")
cat(sprintf("record: valid %s; steps %d; legacy %d; incomplete: %s\n", check$valid, length(operations), sum(operations == "legacy"),
    paste(check$incomplete, collapse = " | ")))
# The replay project: a reflinked copy of the member's source directory, its store cut to the preconditions.
copy <- file.path(work, "copy"); dir.create(copy)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(source_dir, "."), copy)) == 0L)
copy_root <- file.path(copy, basename(store_root))
preconditions <- vapply(read_back$steps[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))], function(s) s$state, "")
prepared <- h$kitPrepareReplayProject(copy_root, preconditions, api)
roles <- vapply(read_back$inputs, function(i) i$role, "")
basenames <- vapply(read_back$inputs, function(i) i$basename, "")
inputs <- stats::setNames(ifelse(roles == "design_file", file.path(source_dir, basenames), file.path(source_dir, "import_sources", basenames)), roles)
paths <- list(source_dir = file.path(work, "paths", "source"), lipid_qc_dir = file.path(work, "paths", "qc"),
    da_output_dir = file.path(work, "paths", "da"), publication_graphs_dir = file.path(work, "paths", "pub"))
for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
scratch <- file.path(work, "scratch"); dir.create(scratch)
replay <- suppressWarnings(suppressMessages(h$kitReplayJourney(read_back, prepared, store_root, inputs, paths, scratch, api)))
mismatches <- as.character(replay$mismatches)
additional <- attr(replay$mismatches, "additional_products") %||% character()
# A tampered input is refused before any step.
tampered <- inputs; tampered[["design_file"]] <- inputs[[which(roles != "design_file")[[1L]]]]
refused <- tryCatch({ api("replayWorkflowAnalysisRecord")(read_back, prepared, tampered, paths, scratch); FALSE },
    workflow_replay_inputs = function(e) TRUE)
result <- list(schema = "kit.member_record_replay.v1", ticket = "KIT-031", member = member, store = store_dir,
    record = list(valid = check$valid, errors = check$errors, incomplete = check$incomplete, operations = operations,
        inputs = lapply(read_back$inputs, function(i) i[c("role", "basename", "sha256")]),
        da = read_back$analyses$differential_abundance[c("operation", "source_state", "model")]),
    replay = list(preconditions = preconditions, states = replay$states, mismatches = mismatches, additional_products = additional,
        per_step = replay$per_step, exact = !length(mismatches)),
    inputs_refused = refused, seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
jsonlite::write_json(result, out, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
cat(sprintf("replay: %s; states %d; mismatches %d; additional products %d; tampered input refused %s; %.0f s\n",
    if (!length(mismatches)) "EXACT" else "DIFFERS", length(replay$states), length(mismatches), length(additional), refused, result$seconds))
