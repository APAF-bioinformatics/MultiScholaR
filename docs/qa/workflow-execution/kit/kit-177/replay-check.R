# KIT-177 gate F3 (ticket section 8): does a label-free and a TMT record replay on this commit, state by state?
#
# Usage, from the package root on a run-tests runner (the commit under test installed, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-177/replay-check.R [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-177/replay
#
# For the e2e fixtures prot_lfq_maxquant (MaxQuant, LFQ) and prot_tmt_pd (Proteome Discoverer, TMT): the package's own
# import and design step in process, as test-workflow-native-session.R's KIT-020 test drives them (an on-disk owner, the
# inputs registered by role, the reader, `applyProtImportResultToWorkflow`, `buildProtDesignStateCheckpoint`), then the
# route's declared operations (route.R: Accession Cleanup, Protein Intensity Filter, Duplicate Removal, Protein Replicate
# Filter, normalisation with RUV skipped, correlation filter 0.5, DA). The record is written, the store copied and cut
# to the record's import and design states (tools/profiling/workflow_kit.R `kitPrepareReplayProject`) and the record
# replayed into it (`replayWorkflowAnalysisRecord`, inputs checked by sha256). Compared state by state: each replayed
# state with the first run's by identical(), and where not identical the slots that differ (for `args`, the entries that
# differ, flagged when either value names a path); the DA table alike; and the stores by the kit fingerprint
# (`kitFingerprintStore`, `kitReproduceMismatches`), the project's established replay criterion.
# Written: <out>/replay.tsv (lane, commit, measure, value) and <out>/replay-slots.tsv (lane, state, slot, entry, path_bearing).
# The outcome per lane (ticket F3): `identical` (outcome 1), `not_supported` (replay stops: outcome 2), or
# `in_order_not_identical` (the replayed history equals the record's, a state or the DA table not identical: outcome 3).

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) args[[1L]] else "docs/qa/workflow-execution/kit/kit-177/replay"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-177 F3 at", commit, "\n")

ns <- asNamespace("MultiScholaR")
env <- new.env(parent = ns)
sys.source("docs/qa/workflow-execution/kit/kit-177/route.R", env)
kit <- new.env(parent = globalenv())
sys.source(file.path("tools", "profiling", "workflow_kit.R"), kit)
api <- function(name) get(name, envir = ns)

rows <- list()
slot_rows <- list()
put <- function(lane, measure, value) {
    value <- if (is.null(value) || !length(value)) NA_character_ else paste(as.character(value), collapse = ";")
    rows[[length(rows) + 1L]] <<- data.frame(lane = lane, commit = commit, measure = measure,
        value = gsub("[\t\r\n]+", " ", value), stringsAsFactors = FALSE)
    utils::write.table(do.call(rbind, rows), file.path(out_dir, "replay.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
isPath <- function(x) is.character(x) && any(grepl("/", x, fixed = TRUE))
slotDifferences <- function(lane, state, a, b) {
    if (!isS4(a) || !isS4(b) || !identical(class(a), class(b))) {
        slot_rows[[length(slot_rows) + 1L]] <<- data.frame(lane = lane, state = state, slot = "<class>", entry = "",
            path_bearing = FALSE, stringsAsFactors = FALSE)
        return("<class>")
    }
    differing <- character()
    for (slot in methods::slotNames(a)) {
        x <- methods::slot(a, slot)
        y <- methods::slot(b, slot)
        if (identical(x, y)) next
        differing <- c(differing, slot)
        entries <- if (is.list(x) && is.list(y) && !is.data.frame(x)) {
            names <- union(names(x), names(y))
            names[!vapply(names, function(n) identical(x[[n]], y[[n]]), logical(1))]
        } else ""
        for (entry in entries) {
            flat <- c(unlist(if (nzchar(entry)) x[[entry]] else NULL), unlist(if (nzchar(entry)) y[[entry]] else NULL))
            slot_rows[[length(slot_rows) + 1L]] <<- data.frame(lane = lane, state = state, slot = slot, entry = entry,
                path_bearing = isPath(flat), stringsAsFactors = FALSE)
        }
    }
    differing
}
readDa <- function(manager) tryCatch(withDIAStateScope(manager, function() {
    reference <- manager$.context$snapshot$da_results
    if (is.null(reference)) return(NULL)
    # As the record reads it (`buildWorkflowAnalysisRecord`): the head's DA manifest, then the contrast's long table.
    manifest <- resolveWorkflowStateValue(manager$.context, reference)
    long <- readDIAContrastOutput(manifest$results[[1L]], "da_proteins_long")$da_proteins_long
    if (is.null(long)) NULL else as.data.frame(long)
}), error = function(e) NULL)
readDa <- (function(f) { environment(f) <- env; f })(readDa)

lanes <- list(
    prot_lfq_maxquant = list(format = "maxquant", type = "LFQ", file = "proteinGroups.txt", reader = "importMaxQuantData"),
    prot_tmt_pd = list(format = "pd_tmt", type = "TMT", file = "pd_tmt_peptides.tsv", reader = "importProteomeDiscovererTMTData"))
scratch_root <- file.path(tempdir(), "kit177-replay")
dir.create(scratch_root, showWarnings = FALSE)

for (lane in names(lanes)) {
    spec <- lanes[[lane]]
    started <- Sys.time()
    where <- character()
    outcome <- tryCatch(withCallingHandlers(local({
        fixture <- normalizePath(file.path("tests", "testdata", "e2e", lane))
        root <- file.path(scratch_root, lane, "origin")
        paths <- list(source_dir = file.path(root, "scripts", "proteomics"), protein_qc_dir = file.path(root, "results", "pqc"),
            peptide_qc_dir = file.path(root, "results", "peptide_qc"), da_output_dir = file.path(root, "results", "da"),
            publication_graphs_dir = file.path(root, "results", "pub"), time_dir = file.path(root, "results", "time"),
            results_dir = file.path(root, "results"))
        for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
        session <- new.env(parent = emptyenv())
        session$userData <- new.env(parent = emptyenv())
        session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- api("WorkflowState")$new(audit_enabled = FALSE)
        workflow$config_list <- env$k177Config()
        workflow$config_list$globalParameters$workflow_type <- spec$type
        api("initializeDIANativeImport")(workflow, paths, session, spec$format, "iq")
        manager <- workflow$state_manager
        put(lane, "owner_class", class(manager)[[1L]])
        manager$.context$experiment_paths <- paths
        report <- file.path(fixture, spec$file)
        design_file <- utils::read.delim(file.path(fixture, "design.tsv"), stringsAsFactors = FALSE)
        api("withDIAStateScope")(manager, function() {
            api("registerDIAInputFiles")(manager, stats::setNames(c(report, file.path(fixture, "proteins.fasta")),
                c(paste0(spec$format, "_report"), "fasta")))
            imported <- api(spec$reader)(report)
            api("applyProtImportResultToWorkflow")(workflow, imported, spec$format, file.path(fixture, "proteins.fasta"),
                logInfo = function(...) NULL, showNotification = function(...) NULL)
            api("completeDIAInputFiles")(manager)
            runs <- unique(as.character(imported$data[[imported$column_mapping$run_col]]))
            sample <- sub("^.*((WT|KO)_[0-9]+)$", "\\1", runs)
            workflow$design_matrix <- data.frame(Run = runs, group = design_file$group[match(sample, design_file$sample)],
                replicates = as.integer(sub("^.*_([0-9]+)$", "\\1", sample)), stringsAsFactors = FALSE)
            api("buildProtDesignStateCheckpoint")(workflow, spec$type, "Save")
        })
        put(lane, "design_transform", manager$.context$snapshot$step_provenance$protein_s4_initial$fixed$transform)
        first <- env$k177RunPath(manager, env$k177NormParameters(ruv_mode = "skip"),
            list(name = "KO.vs.WT", full = "KO.vs.WT=groupKO-groupWT"))
        put(lane, "first_run_states", paste(first$states, collapse = ","))
        json <- file.path(scratch_root, lane, "analysis_record.json")
        api("writeWorkflowAnalysisRecord")(manager, json)
        states_first <- lapply(stats::setNames(first$states, first$states), function(name)
            api("withDIAStateScope")(manager, function() manager$getState(name)))
        # The first run's DA table read from its store the same way as the replay's (the in-memory read kept beside it).
        da_memory <- first$da
        da_first <- readDa(manager)
        put(lane, "da_first_store_read_rows", if (is.null(da_first)) "not read" else nrow(da_first))
        put(lane, "da_first_store_equals_memory_by_name", !is.null(da_first) && setequal(names(da_first), names(da_memory)) &&
            nrow(da_first) == nrow(da_memory))
        store_root <- dirname(manager$.context$store$directory)
        put(lane, "store_root_is_root", identical(normalizePath(store_root), normalizePath(root)))
        api("closeDIAStateAdapter")(manager)
        record <- jsonlite::read_json(json, simplifyVector = FALSE)
        states <- vapply(record$steps, function(step) step$state, "")
        operations <- vapply(record$steps, function(step) step$operation, "")
        put(lane, "record_workflow", record$workflow)
        put(lane, "record_steps", paste(operations, collapse = ","))
        put(lane, "record_states", paste(states, collapse = ","))
        put(lane, "record_valid", api("validateWorkflowAnalysisRecord")(record)$valid)
        preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
        put(lane, "preconditions", paste(preconditions, collapse = ","))
        copy <- file.path(scratch_root, lane, "replay")
        dir.create(copy, showWarnings = FALSE)
        stopifnot(system2("cp", c("-a", shQuote(file.path(store_root, ".")), shQuote(copy))) == 0L)
        prepared <- kit$kitPrepareReplayProject(copy, preconditions, api)
        inputs <- character()
        for (input in record$inputs) {
            if (is.null(input$role) || is.null(input$basename) || is.na(input$basename)) next
            path <- file.path(fixture, input$basename)
            if (file.exists(path)) inputs[[input$role]] <- path
        }
        put(lane, "inputs_found", paste(names(inputs), collapse = ","))
        replay_paths <- lapply(paths, function(p) sub(root, copy, p, fixed = TRUE))
        for (p in replay_paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
        replayed <- list()
        operations_seen <- character()
        da_replayed <- NULL
        utils::capture.output(replay <- suppressMessages(suppressWarnings(api("replayWorkflowAnalysisRecord")(record, prepared,
            inputs, replay_paths, file.path(scratch_root, lane, "replay-scratch"), on_step = function(step, replay_manager) {
                operations_seen <<- c(operations_seen, step$operation %||% "")
                if (!is.null(step$state) && step$state %in% states) replayed[[step$state]] <<-
                    api("withDIAStateScope")(replay_manager, function() replay_manager$getState(step$state))
                da <- readDa(replay_manager)
                if (!is.null(da)) da_replayed <<- da
            }))))
        put(lane, "replayed_operations", paste(operations_seen, collapse = ","))
        put(lane, "replayed_states", paste(replay$states, collapse = ","))
        history_equal <- identical(as.character(replay$states), states)
        put(lane, "history_equal", history_equal)
        later <- setdiff(states, preconditions)
        # VAL-79 diagnostic: every replayed state read after the replayed DA, from the replay store reopened at its end
        # (the first run's states were read after its DA).
        identity <- api("readSavedStoreIdentity")(file.path(prepared, ".workflow-entities", "store.json"))
        reopened <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(prepared, identity$project_id, identity$workflow_id),
            reopen = TRUE)
        replayed_end <- lapply(stats::setNames(later, later), function(name)
            api("withDIAStateScope")(reopened, function() reopened$getState(name)))
        da_end <- readDa(reopened)
        api("closeDIAStateAdapter")(reopened)
        same <- vapply(later, function(state) identical(replayed_end[[state]], states_first[[state]]), logical(1))
        same_on_step <- vapply(later, function(state) identical(replayed[[state]], states_first[[state]]), logical(1))
        put(lane, "states_compared", length(later))
        put(lane, "states_identical_after_da", sum(same))
        put(lane, "states_identical_on_step", sum(same_on_step))
        put(lane, "states_not_identical", paste(later[!same], collapse = ","))
        allowed <- TRUE
        for (state in later[!same]) {
            slots <- slotDifferences(lane, state, states_first[[state]], replayed_end[[state]])
            put(lane, paste0("slots_differing:", state), paste(slots, collapse = ","))
            a <- states_first[[state]]@args
            b <- replayed_end[[state]]@args
            entries <- union(names(a), names(b))
            entries <- entries[!vapply(entries, function(n) identical(a[[n]], b[[n]]), logical(1))]
            put(lane, paste0("args_entries_differing:", state), paste(entries, collapse = ","))
            inner <- if (identical(entries, "daAnalysisUI")) {
                x <- a$daAnalysisUI
                y <- b$daAnalysisUI
                fields <- union(names(x), names(y))
                fields[!vapply(fields, function(n) identical(x[[n]], y[[n]]), logical(1))]
            } else NA_character_
            put(lane, paste0("daAnalysisUI_fields_differing:", state), paste(inner, collapse = ","))
            if (!identical(slots, "args") || !identical(entries, "daAnalysisUI") || !identical(inner, "timestamp")) allowed <- FALSE
        }
        da_replayed <- da_end %||% da_replayed
        put(lane, "da_replayed_read", !is.null(da_replayed))
        put(lane, "da_identical", !is.null(da_replayed) && identical(da_replayed, da_first))
        da_by_name <- FALSE
        if (!is.null(da_replayed)) {
            put(lane, "da_columns_first_order", paste(names(da_first), collapse = ","))
            put(lane, "da_columns_replay_order", paste(names(da_replayed), collapse = ","))
            put(lane, "da_column_sets_equal", setequal(names(da_first), names(da_replayed)))
            put(lane, "da_rows", paste(nrow(da_first), nrow(da_replayed)))
            id <- if ("Protein.Ids" %in% names(da_first)) "Protein.Ids" else names(da_first)[[1L]]
            key <- function(table) paste(table[[id]], table[["comparison"]] %||% "")
            first <- da_first[order(key(da_first)), , drop = FALSE]
            again <- da_replayed[order(key(da_replayed)), , drop = FALSE]
            rows_same <- nrow(first) == nrow(again) && identical(key(first), key(again))
            put(lane, "da_row_keys_equal", rows_same)
            put(lane, "da_row_order_equal", identical(key(da_first), key(da_replayed)))
            columns <- intersect(names(da_first), names(da_replayed))
            differing <- if (rows_same) columns[!vapply(columns, function(n) identical(unname(first[[n]]), unname(again[[n]])),
                logical(1))] else columns
            put(lane, "da_columns_differing_by_name", paste(differing, collapse = ","))
            da_by_name <- rows_same && setequal(names(da_first), names(da_replayed)) && !length(differing)
        }
        put(lane, "da_equal_by_name", da_by_name)
        put(lane, "only_settings_timestamp_and_column_order", history_equal && allowed && da_by_name)
        origin <- kit$kitFingerprintStore(file.path(store_root, ".workflow-entities"), api)
        again <- kit$kitFingerprintStore(file.path(prepared, ".workflow-entities"), api)
        as_populations <- function(fingerprints) lapply(fingerprints$states, function(x) list(fingerprint = x))
        mismatches <- as.vector(kit$kitReproduceMismatches(as_populations(origin), as_populations(again)))
        put(lane, "kit_fingerprint_mismatches", if (length(mismatches)) paste(mismatches, collapse = ",") else "none")
        if (history_equal && all(same) && identical(da_replayed, da_first)) "identical" else
            if (history_equal && allowed && da_by_name) "in_order_settings_timestamp_and_column_order_only" else
            if (history_equal) "in_order_not_identical" else "history_differs"
    }), error = function(e) {
        calls <- vapply(utils::tail(sys.calls(), 8L), function(call) substr(paste(deparse(call, nlines = 1L), collapse = " "),
            1L, 120L), "")
        where <<- calls
    }), error = function(e) {
        put(lane, "stopped", paste(class(e)[[1L]], conditionMessage(e)))
        put(lane, "stopped_in", paste(where, collapse = " <- "))
        "not_supported_or_stopped"
    })
    put(lane, "outcome", outcome)
    put(lane, "seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
if (length(slot_rows)) utils::write.table(do.call(rbind, slot_rows), file.path(out_dir, "replay-slots.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
cat("KIT-177 F3 done\n")
