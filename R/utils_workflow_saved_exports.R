#' Export one DIA generation, retaining only descriptors for retry and owned reuse.
#' @param manager Active DIA adapter.
#' @param workflowData,normData,input Existing module inputs.
#' @param sourceDir Existing session destination.
#' @param collectFn,saveFn Existing collection and artifact producers.
#' @param serializeFn Writer for the complete legacy RDS, before owner release.
#' @return Legacy-shaped output references; no resident session history.
#' @noRd
runDIAOwnedSessionExport <- function(manager, workflowData, normData, input, sourceDir,
    collectFn, saveFn, serializeFn = saveRDS) {
    context <- manager$.context
    root <- normalizePath(sourceDir, mustWork = TRUE)
    request <- context$export_pending
    if (is.null(request) || !identical(request$generation, context$head$entity_id) ||
        !identical(request$root, root)) {
        stamp <- Sys.time()
        # Keep the existing timestamp filename convention without overwriting a
        # completed same-second export. This affects filenames, not scientific time.
        while (file.exists(file.path(root, paste0("filtered_session_data_",
            format(stamp, "%Y%m%d_%H%M%S"), ".rds")))) stamp <- stamp + 1
        request <- list(request_id = basename(tempfile("dia-export-")),
            timestamp = format(stamp, "%Y%m%d_%H%M%S"), root = root,
            generation = context$head$entity_id, identity = context$store$identity,
            state_reference = context$snapshot$states[[context$snapshot$current_state]])
        context$export_pending <- request
    }
    paths <- workflowExportPaths(root, request$request_id)
    if (!is.null(request$metadata_reference) && file.exists(file.path(paths$request_dir, "journal.json"))) {
        session_data <- resolveDIAStateValue(context, request$metadata_reference)
    } else {
        session_data <- serializeDIAOwnedSessionExport(manager, request,
            workflowData, normData, input, collectFn, serializeFn)
        request <- context$export_pending
        # The child frame and working-map roots are gone; the full session is
        # durable on disk. Reclaim that unreachable history before file copying.
        gc(full = TRUE)
        workflowExportBoundary("after_history_release")
    }
    artifacts <- saveFn(sessionData = session_data, sourceDir = root, publicationRequest = request)
    if (!identical(artifacts$publication$status, "complete")) {
        cli::cli_abort("session publication did not complete", class = "saved_entity_outputs_pending")
    }
    request$sha256 <- artifacts$publication$session_sha256
    request$bytes <- artifacts$publication$session_bytes
    request$path <- artifacts$latestFilepath
    request$completion <- file.path(root, paste0("export_complete_", request$request_id, ".json"))
    context$owned_export <- request
    context$export_pending <- NULL
    list(sessionData = NULL, exportArtifacts = artifacts)
}

#' Serialize a complete DIA history without returning it to publication callers.
#' @param manager Active DIA adapter.
#' @param request Owned generation and destination metadata.
#' @param workflowData,normData,input Existing module inputs.
#' @param collectFn,serializeFn Existing session collector and legacy RDS writer.
#' @return Metadata-only legacy view; all full-history call-frame roots expire here.
#' @noRd
serializeDIAOwnedSessionExport <- function(manager, request, workflowData, normData,
    input, collectFn, serializeFn) {
    context <- manager$.context
    on.exit({ context$working <- new.env(parent = emptyenv()) }, add = TRUE)
    session_data <- collectFn(workflowData = workflowData, normData = normData, input = input)
    if (!is.list(session_data) ||
        !identical(names(session_data$r6_complete_states), names(context$snapshot$states))) {
        cli::cli_abort("DIA export history is incomplete", class = "saved_entity_corrupt")
    }
    request$correlation_reference <- saveDIAOwnedValue(session_data$correlation_filtered_s4,
        "export-correlation", manager)
    metadata <- session_data
    metadata[c("r6_complete_states", "current_s4_object", "correlation_filtered_s4")] <- NULL
    metadata$config_list <- prepareDIAStateRecord(list(data = NULL, config = metadata$config_list))$config
    request$metadata_reference <- saveDIAStateValue(context, metadata, "export-metadata")
    context$export_pending <- request
    paths <- workflowExportPaths(request$root, request$request_id)
    stageWorkflowLegacyExport(session_data, paths, request$request_id, request$timestamp,
        serializeFn = serializeFn, serialize_only = TRUE)
    metadata
}

#' Refuse dependent work until the current export has finished publication.
#' @param manager Existing state owner, possibly NULL or native.
#' @return NULL invisibly, or a pending-output condition.
#' @noRd
assertDIAExportReady <- function(manager, allow_action = NULL) {
    if (usesDIAReferenceStates(manager)) {
        reference <- getDIAStateReference(manager)
        request <- manager$.context$snapshot$action_outputs[[reference$entity_id]]
        completed <- !is.null(request) && isTRUE(manager$.context$snapshot$completed_outputs[[reference$entity_id]][[request$label]])
        if (!is.null(request) && !completed && !identical(request$label, allow_action)) {
            paths <- workflowExportPaths(request$root, request$request_id)
            journal <- if (file.exists(file.path(paths$request_dir, "journal.json")))
                jsonlite::read_json(file.path(paths$request_dir, "journal.json")) else NULL
            if (!identical(journal$status, "complete")) cli::cli_abort(
                paste0(request$label, " outputs are pending; retry that action"),
                class = "saved_entity_outputs_pending")
        }
    }
    if (usesDIAReferenceStates(manager) && manager$current_state %in% c("qvalue_filtered",
        "protein_peptide_filtered", "precursor_rollup", "intensity_filtered", "sample_filtered",
        "replicate_filtered", "imputed", "protein_s4_created", "protein_accession_cleaned",
        "protein_intensity_filtered", "duplicates_removed", "protein_replicate_filtered") &&
        !is.null(resumeDIAQvalueOutput(manager, readDIAStateMetadata(manager), hydrate = FALSE))) {
        cli::cli_abort("QC outputs are pending; retry Apply to publish them",
            class = "saved_entity_outputs_pending")
    }
    if (inherits(manager, "DIAStateAdapter") && !is.null(manager$.context$export_pending) &&
        !identical(allow_action, "session-export")) {
        cli::cli_abort("session outputs are pending; retry export before dependent actions",
            class = "saved_entity_outputs_pending")
    }
    invisible(NULL)
}

#' Publish included DIA action files separately from their accepted science.
#' @param manager,reference Owner and accepted scientific checkpoint.
#' @param label Existing action name used by its retry boundary.
#' @param destinations Named existing output directories.
#' @param produce Function writing only to corresponding staged directories.
#' @return Bounded display handoff after all files have been acknowledged.
#' @noRd
publishDIAActionOutputs <- function(manager, reference, label, destinations, produce, identity = NULL) {
    context <- manager$.context
    destinations <- destinations[!vapply(destinations, is.null, logical(1))]
    root <- context$experiment_paths$source_dir %||% destinations[[1L]]
    request_key <- digest::digest(list(state = reference$entity_id, action = label, identity = identity), algo = "sha256")
    request_id <- paste0("dia-", label, "-", substr(request_key, 1L, 32L))
    paths <- workflowExportPaths(root, request_id)
    journal_path <- file.path(paths$request_dir, "journal.json")
    owner <- function() {
        checkWorkflowEntityStore(context$store)
        if (!identical(getDIAStateReference(manager), reference)) {
            cli::cli_abort("action publication was superseded", class = "saved_entity_scope")
        }
    }
    owner()
    snapshot <- context$snapshot
    snapshot$action_outputs[[reference$entity_id]] <- list(root = root, request_id = request_id, label = label)
    snapshot$completed_outputs[[reference$entity_id]][[label]] <- FALSE
    publishDIAStateHead(context, snapshot)
    journal <- if (file.exists(journal_path)) jsonlite::read_json(journal_path) else list(
        schema = "multischolar.dia_qc_publication.v1", status = "outputs_pending",
        request_id = request_id, state_id = reference$entity_id,
        incarnation = context$store$identity$incarnation, files = NULL)
    if (!identical(journal$state_id, reference$entity_id) ||
        !identical(journal$incarnation, context$store$identity$incarnation)) {
        cli::cli_abort("action publication identity differs", class = "saved_entity_scope")
    }
    was_complete <- identical(journal$status, "complete")
    writeWorkflowExportJournal(paths, journal)
    staged_root <- file.path(paths$request_dir, "files")
    # Until an inventory is journaled, staged files belong to no accepted attempt.
    if (is.null(journal$files)) unlink(staged_root, recursive = TRUE)
    staged <- stats::setNames(lapply(names(destinations), function(role) {
        checkSavedEntityString(role, "output role")
        path <- file.path(staged_root, role)
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        path
    }), names(destinations))
    display_path <- file.path(paths$request_dir, "display.rds")
    if (is.null(journal$files)) {
        workflowExportBoundary("before_action_outputs")
        display <- produce(staged)
        if (as.numeric(object.size(display)) > 67108864) cli::cli_abort(
            "action display handoff exceeds limits", class = "saved_entity_limit")
        saveRDS(display, display_path, version = 3)
        members <- list.files(staged_root, recursive = TRUE, full.names = FALSE)
        if (!length(members) || length(members) > 10000L) cli::cli_abort(
            "action output inventory is empty or oversized", class = "saved_entity_outputs_pending")
        journal$files <- stats::setNames(lapply(members, function(member) {
            checkSavedEntityPath(context$store, member)
            path <- file.path(staged_root, member)
            if (file.size(path) <= 0) cli::cli_abort("action output is empty", class = "saved_entity_outputs_pending")
            list(sha256 = digest::digest(file = path, algo = "sha256"), bytes = file.size(path))
        }), members)
        journal$display_sha256 <- digest::digest(file = display_path, algo = "sha256")
        syncWorkflowExportPaths(c(file.path(staged_root, members), display_path, staged_root), paths$root)
        journal$status <- "staged"
        writeWorkflowExportJournal(paths, journal)
    }
    previous <- file.path(paths$request_dir, "previous")
    for (member in names(journal$files)) {
        owner()
        fields <- strsplit(member, "/", fixed = TRUE)[[1L]]
        if (length(fields) < 2L || !fields[[1L]] %in% names(destinations)) {
            cli::cli_abort("action output path differs", class = "saved_entity_corrupt")
        }
        source <- file.path(staged_root, member)
        expected <- journal$files[[member]]
        if (!file.exists(source) || !isTRUE(file.size(source) == expected$bytes) ||
            !identical(digest::digest(file = source, algo = "sha256"), expected$sha256)) {
            cli::cli_abort("staged action output changed", class = "saved_entity_corrupt")
        }
        destination <- file.path(destinations[[fields[[1L]]]], paste(fields[-1L], collapse = "/"))
        if (!dir.exists(dirname(destination))) dir.create(dirname(destination), recursive = TRUE)
        backup <- file.path(previous, member)
        if (file.exists(destination) && !file.exists(backup)) {
            dir.create(dirname(backup), recursive = TRUE, showWarnings = FALSE)
            if (!file.copy(destination, backup)) cli::cli_abort("cannot retain previous action output")
        }
        temporary <- tempfile(".publish-", tmpdir = dirname(destination))
        if (!file.copy(source, temporary)) cli::cli_abort("cannot stage action output publication")
        syncWorkflowExportPaths(temporary, dirname(destination))
        if (!file.rename(temporary, destination)) cli::cli_abort("cannot publish action output")
        syncWorkflowExportPaths(c(destination, dirname(destination)), dirname(destination))
        workflowExportBoundary("between_action_publications")
    }
    if (!identical(digest::digest(file = display_path, algo = "sha256"), journal$display_sha256)) {
        cli::cli_abort("action display changed", class = "saved_entity_corrupt")
    }
    workflowExportBoundary("before_action_acknowledgement")
    journal$status <- "complete"
    writeWorkflowExportJournal(paths, journal)
    acknowledgeDIAOutputs(manager, reference, label)
    if (!was_complete) requestDIAReclamation(manager, reference, request_id, label)
    readRDS(display_path)
}

#' Record portable completion independently of local publication-journal paths.
#' @param manager,reference Accepted state owner and descriptor.
#' @param label Completed output group.
#' @return NULL invisibly; this publishes metadata, never a scientific checkpoint.
#' @noRd
acknowledgeDIAOutputs <- function(manager, reference, label) {
    snapshot <- manager$.context$snapshot
    snapshot$completed_outputs[[reference$entity_id]][[label]] <- TRUE
    publishDIAStateHead(manager$.context, snapshot)
    invisible(NULL)
}

#' Reuse only the unchanged completed export belonging to the current DIA head.
#' @param manager Current state owner, possibly the ordinary native manager.
#' @param path Selected session file; external or changed files return NULL.
#' @return Minimal legacy session view, or NULL for the existing full restore path.
#' @noRd
resolveDIAOwnedExport <- function(manager, path) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    context <- manager$.context
    checkWorkflowEntityStore(context$store)
    assertDIAExportReady(manager)
    request <- context$owned_export
    if (is.null(request) || !file.exists(path) ||
        !identical(request$identity, context$store$identity) ||
        !identical(request$generation, context$head$entity_id) ||
        !identical(normalizePath(path), normalizePath(request$path)) ||
        !identical(as.numeric(file.size(path)), as.numeric(request$bytes)) ||
        !file.exists(request$completion)) return(NULL)
    completion <- jsonlite::read_json(request$completion)
    if (!identical(completion$status, "complete") ||
        !identical(completion$request_id, request$request_id) ||
        !identical(completion$session_sha256, request$sha256) ||
        !identical(digest::digest(file = path, algo = "sha256"), request$sha256)) return(NULL)
    session <- resolveDIAStateValue(context, request$metadata_reference)
    session$current_s4_object <- resolveDIAStateValue(context, request$state_reference)$data
    session$correlation_filtered_s4 <- resolveDIAOwnedValue(request$correlation_reference, manager)
    session
}

#' Identify one immutable descriptor accepted by the shared saved-entity store.
#' @param value Candidate value.
#' @return TRUE for descriptor-shaped references.
#' @noRd
isSavedWorkflowReference <- function(value) {
    is.list(value) && !is.data.frame(value) && identical(value[["schema", exact = TRUE]], 1L) &&
        is.character(value$entity_id) && length(value$entity_id) == 1L &&
        grepl("^[a-f0-9]{64}$", value$entity_id) && is.list(value$payload)
}

#' Create one export-scoped decoder owner with identity de-duplication.
#' @param store Open shared saved-entity store.
#' @return Mutable scope; call releaseWorkflowExportHydrationScope when finished.
#' @noRd
newWorkflowExportHydrationScope <- function(store) {
    checkWorkflowEntityStore(store)
    scope <- new.env(parent = emptyenv())
    scope$store <- store
    scope$values <- new.env(hash = TRUE, parent = emptyenv())
    scope$unique_decodes <- 0L
    scope$reference_uses <- 0L
    scope$mutable_clones <- 0L
    scope$released <- FALSE
    scope
}

#' Decode one reference once within an active export scope.
#' @param scope Export hydration scope.
#' @param reference Accepted saved descriptor.
#' @return Hydrated value; mutable-reference payloads receive an isolated clone.
#' @noRd
hydrateWorkflowExportReference <- function(scope, reference) {
    if (!is.environment(scope) || isTRUE(scope$released) ||
        !isSavedWorkflowReference(reference)) {
        cli::cli_abort("invalid export hydration request", class = "saved_entity_invalid")
    }
    scope$reference_uses <- scope$reference_uses + 1L
    id <- reference$entity_id
    if (!exists(id, envir = scope$values, inherits = FALSE)) {
        value <- hydrateSavedWorkflowEntity(scope$store, reference)
        assign(id, value, envir = scope$values)
        scope$unique_decodes <- scope$unique_decodes + 1L
    }
    value <- get(id, envir = scope$values, inherits = FALSE)
    if (hasMutableWorkflowEntityOwner(value)) {
        scope$mutable_clones <- scope$mutable_clones + 1L
        value <- unserialize(serialize(value, NULL, version = 3))
    }
    value
}

#' Hydrate the existing complete legacy session fields within one scope.
#' @param scope Active export hydration scope.
#' @param session_template Existing self-contained session shape with optional refs.
#' @return Complete resident legacy session; no reference-only substitution.
#' @noRd
hydrateWorkflowExportSession <- function(scope, session_template) {
    if (!is.list(session_template) || !is.list(session_template$r6_complete_states)) {
        cli::cli_abort("invalid legacy session template", class = "saved_entity_invalid")
    }
    session <- session_template
    states <- session$r6_complete_states
    for (name in names(states)) {
        state <- states[[name]]
        if (!is.list(state) || is.null(names(state))) {
            cli::cli_abort("invalid legacy history record", class = "saved_entity_invalid")
        }
        if (isSavedWorkflowReference(state$data)) {
            state$data <- hydrateWorkflowExportReference(scope, state$data)
        }
        states[[name]] <- state
    }
    session$r6_complete_states <- states
    for (field in c("current_s4_object", "correlation_filtered_s4")) {
        if (isSavedWorkflowReference(session[[field]])) {
            session[[field]] <- hydrateWorkflowExportReference(scope, session[[field]])
        }
    }
    session
}

#' Report export-scope counters without hydrating another value.
#' @param scope Active export hydration scope.
#' @return Decode/use/clone and resident-entry counters.
#' @noRd
workflowExportHydrationInfo <- function(scope) {
    if (!is.environment(scope)) {
        cli::cli_abort("invalid export hydration scope", class = "saved_entity_invalid")
    }
    list(unique_decodes = scope$unique_decodes,
        reference_uses = scope$reference_uses,
        mutable_clones = scope$mutable_clones,
        resident_entries = if (isTRUE(scope$released)) 0L else
            length(ls(scope$values, all.names = TRUE)))
}

#' Release every export-owned hydrated value.
#' @param scope Export hydration scope.
#' @return TRUE invisibly; repeated release is harmless.
#' @noRd
releaseWorkflowExportHydrationScope <- function(scope) {
    if (!is.environment(scope)) {
        cli::cli_abort("invalid export hydration scope", class = "saved_entity_invalid")
    }
    if (isTRUE(scope$released)) return(invisible(TRUE))
    scope$values <- new.env(hash = TRUE, parent = emptyenv())
    scope$store <- NULL
    scope$released <- TRUE
    invisible(TRUE)
}

#' Named publication boundaries for deterministic failure injection.
#' @param phase Publication phase.
#' @return NULL invisibly.
#' @noRd
workflowExportBoundary <- function(phase) invisible(NULL)

#' Validate a request-owned export directory and token.
#' @param output_dir Existing destination directory.
#' @param request_id Idempotent request token.
#' @return Canonical root and request staging paths.
#' @noRd
workflowExportPaths <- function(output_dir, request_id) {
    if (!is.character(output_dir) || length(output_dir) != 1L ||
        !dir.exists(output_dir) || !is.character(request_id) ||
        length(request_id) != 1L || is.na(request_id) ||
        !grepl("^[A-Za-z0-9_-]{1,80}$", request_id)) {
        cli::cli_abort("invalid export publication request", class = "saved_entity_invalid")
    }
    root <- normalizePath(output_dir, mustWork = TRUE)
    staging_root <- file.path(root, ".workflow-export")
    if (!dir.exists(staging_root) && !dir.create(staging_root)) {
        cli::cli_abort("cannot create export staging root", class = "saved_entity_write")
    }
    request_dir <- file.path(staging_root, request_id)
    if (!dir.exists(request_dir) && !dir.create(request_dir)) {
        cli::cli_abort("cannot create export request directory", class = "saved_entity_write")
    }
    list(root = root, staging_root = normalizePath(staging_root),
        request_dir = normalizePath(request_dir))
}

#' Flush explicit export-owned files and directories.
#' @param paths Existing paths.
#' @param root Export destination root.
#' @return TRUE invisibly.
#' @noRd
syncWorkflowExportPaths <- function(paths, root) {
    command <- unname(Sys.which("sync"))
    if (!identical(Sys.info()[["sysname"]], "Linux") || !nzchar(command) ||
        !length(paths) || any(!file.exists(paths) & !dir.exists(paths))) {
        cli::cli_abort("durable export sync is unavailable", class = "saved_entity_unsupported")
    }
    root <- normalizePath(root, mustWork = TRUE)
    paths <- unique(normalizePath(paths, mustWork = TRUE))
    if (!all(paths == root | startsWith(paths, paste0(root, "/")))) {
        cli::cli_abort("export sync path escapes destination", class = "saved_entity_path")
    }
    if (!identical(system2(command, c("--", shQuote(paths))), 0L)) {
        cli::cli_abort("durable export sync failed", class = "saved_entity_write")
    }
    invisible(TRUE)
}

#' Atomically persist metadata for an already durable owned export stage.
#' @param paths Existing request paths.
#' @param journal Existing v1 publication metadata; never scientific objects.
#' @return Journal invisibly, after its file and directory are durable.
#' @noRd
writeWorkflowExportJournal <- function(paths, journal) {
    temporary <- tempfile("journal-", tmpdir = paths$request_dir)
    tryCatch(jsonlite::write_json(journal, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null"),
        error = function(error) {
            unlink(temporary)
            stop(error)
        })
    syncWorkflowExportPaths(temporary, paths$root)
    journal_path <- file.path(paths$request_dir, "journal.json")
    if (!file.rename(temporary, journal_path)) {
        cli::cli_abort("export journal publication failed", class = "saved_entity_write")
    }
    syncWorkflowExportPaths(c(journal_path, paths$request_dir, paths$staging_root, paths$root), paths$root)
    invisible(journal)
}

#' Finish a DIA latest-file copy using only a durable serialized-session journal.
#' @param paths Existing request paths.
#' @param copyFn,boundaryFn Existing copy operation and failure-injection boundary.
#' @return Outputs-pending journal with an independently verified latest checksum.
#' @noRd
completeWorkflowLegacyCopy <- function(paths, copyFn = file.copy,
    boundaryFn = workflowExportBoundary) {
    journal <- jsonlite::read_json(file.path(paths$request_dir, "journal.json"))
    if (!identical(journal$schema, "multischolar.workflow_export.v1") ||
        !identical(journal$status, "outputs_pending") ||
        !identical(journal$request_id, basename(paths$request_dir)) ||
        !is.null(journal$latest_sha256) || length(journal$additional) ||
        !identical(journal$staged_session, file.path(paths$request_dir, "session.rds")) ||
        !identical(journal$staged_latest, file.path(paths$request_dir, "latest.rds")) ||
        !isTRUE(as.numeric(file.size(journal$staged_session)) == journal$session_bytes)) {
        cli::cli_abort("serialized export journal is invalid", class = "saved_entity_corrupt")
    }
    boundaryFn("before_latest_copy")
    if (!isTRUE(copyFn(journal$staged_session, journal$staged_latest, overwrite = TRUE))) {
        cli::cli_abort("latest legacy copy failed", class = "saved_entity_write")
    }
    latest_hash <- digest::digest(file = journal$staged_latest, algo = "sha256")
    if (!identical(latest_hash, journal$session_sha256)) {
        cli::cli_abort("latest legacy copy differs", class = "saved_entity_corrupt")
    }
    syncWorkflowExportPaths(journal$staged_latest, paths$root)
    journal$latest_sha256 <- latest_hash
    writeWorkflowExportJournal(paths, journal)
    boundaryFn("after_latest_copy")
    journal
}

#' Stage one complete legacy serialization and independently verified copies.
#' @param session_data Complete self-contained legacy session.
#' @param paths Request paths from workflowExportPaths.
#' @param request_id Request token.
#' @param timestamp Fixed YYYYmmdd_HHMMSS string.
#' @param additional_files Named list of source and relative destination entries.
#' @param serializeFn,copyFn Injectable file operations for failure tests.
#' @param boundaryFn Failure-injection boundary.
#' @param serialize_only Internal DIA phase: journal serialization before any copies.
#' @return Outputs-pending journal.
#' @noRd
stageWorkflowLegacyExport <- function(session_data, paths, request_id, timestamp,
    additional_files = list(), serializeFn = saveRDS, copyFn = file.copy,
    boundaryFn = workflowExportBoundary, serialize_only = FALSE) {
    if (!is.list(session_data) || !grepl("^[0-9]{8}_[0-9]{6}$", timestamp) ||
        !is.logical(serialize_only) || length(serialize_only) != 1L || is.na(serialize_only) ||
        isTRUE(serialize_only) && length(additional_files) ||
        !is.list(additional_files) ||
        length(additional_files) && (is.null(names(additional_files)) ||
            anyDuplicated(names(additional_files)))) {
        cli::cli_abort("invalid staged legacy export", class = "saved_entity_invalid")
    }
    session_name <- paste0("filtered_session_data_", timestamp, ".rds")
    staged_session <- file.path(paths$request_dir, "session.rds")
    staged_latest <- file.path(paths$request_dir, "latest.rds")
    boundaryFn("before_serialization")
    serializeFn(session_data, staged_session)
    boundaryFn("after_serialization")
    if (!file.exists(staged_session) || file.size(staged_session) <= 0) {
        cli::cli_abort("legacy session serialization produced no file",
            class = "saved_entity_write")
    }
    session_hash <- digest::digest(file = staged_session, algo = "sha256")
    if (serialize_only) {
        journal <- list(schema = "multischolar.workflow_export.v1",
            request_id = request_id, status = "outputs_pending",
            session_name = session_name, session_sha256 = session_hash,
            session_bytes = file.size(staged_session), staged_session = staged_session,
            staged_latest = staged_latest, latest_sha256 = NULL,
            additional = list(), serialization_count = 1L,
            scientific_checkpoint_count = 0L, promotion_authorized = FALSE)
        syncWorkflowExportPaths(c(staged_session, paths$request_dir), paths$root)
        writeWorkflowExportJournal(paths, journal)
        boundaryFn("after_serialized_checkpoint")
        return(journal)
    }
    if (!isTRUE(copyFn(staged_session, staged_latest, overwrite = TRUE))) {
        cli::cli_abort("latest legacy copy failed", class = "saved_entity_write")
    }
    latest_hash <- digest::digest(file = staged_latest, algo = "sha256")
    if (!identical(latest_hash, session_hash)) {
        cli::cli_abort("latest legacy copy differs", class = "saved_entity_corrupt")
    }
    additions <- list()
    if (length(additional_files)) {
        additions_dir <- file.path(paths$request_dir, "additional")
        if (!dir.exists(additions_dir) && !dir.create(additions_dir)) {
            cli::cli_abort("cannot create additional-output staging",
                class = "saved_entity_write")
        }
        for (index in seq_along(additional_files)) {
            spec <- additional_files[[index]]
            destination <- spec$destination
            if (!is.list(spec) || !is.character(spec$source) ||
                length(spec$source) != 1L || !file.exists(spec$source) ||
                !is.character(destination) || length(destination) != 1L ||
                grepl("(^/|\\\\|(^|/)\\.\\.(/|$))", destination)) {
                cli::cli_abort("invalid additional export file", class = "saved_entity_path")
            }
            staged <- file.path(additions_dir, sprintf("%04d", index))
            if (!isTRUE(copyFn(spec$source, staged, overwrite = TRUE))) {
                cli::cli_abort("additional export copy failed", class = "saved_entity_write")
            }
            source_hash <- digest::digest(file = spec$source, algo = "sha256")
            copy_hash <- digest::digest(file = staged, algo = "sha256")
            if (!identical(source_hash, copy_hash)) {
                cli::cli_abort("additional export copy differs", class = "saved_entity_corrupt")
            }
            additions[[names(additional_files)[[index]]]] <- list(staged = staged,
                destination = destination, source_sha256 = source_hash,
                copy_sha256 = copy_hash, bytes = file.size(staged))
        }
    }
    journal <- list(schema = "multischolar.workflow_export.v1",
        request_id = request_id, status = "outputs_pending",
        session_name = session_name, session_sha256 = session_hash,
        session_bytes = file.size(staged_session), staged_session = staged_session,
        staged_latest = staged_latest, latest_sha256 = latest_hash,
        additional = additions, serialization_count = 1L,
        scientific_checkpoint_count = 0L, promotion_authorized = FALSE)
    # Staged bytes are durable before the journal that names them is published.
    syncWorkflowExportPaths(c(staged_session, staged_latest,
        paths$request_dir, paths$staging_root, paths$root,
        unlist(lapply(additions, `[[`, "staged"), use.names = FALSE)), paths$root)
    writeWorkflowExportJournal(paths, journal)
}

#' Publish or recover one staged export, with latest and completion last.
#' @param paths Request paths.
#' @param boundaryFn Failure-injection boundary.
#' @param renameFn Injectable owned-file rename.
#' @return Complete publication receipt.
#' @noRd
publishWorkflowLegacyExport <- function(paths, boundaryFn = workflowExportBoundary,
    renameFn = file.rename) {
    journal_path <- file.path(paths$request_dir, "journal.json")
    if (!file.exists(journal_path)) {
        cli::cli_abort("export journal is missing", class = "saved_entity_missing")
    }
    journal <- jsonlite::read_json(journal_path, simplifyVector = FALSE)
    if (!identical(journal$schema, "multischolar.workflow_export.v1") ||
        !journal$status %in% c("outputs_pending", "complete")) {
        cli::cli_abort("export journal is invalid", class = "saved_entity_corrupt")
    }
    if (is.null(journal$latest_sha256)) {
        cli::cli_abort("latest export copy is pending", class = "saved_entity_outputs_pending")
    }
    session_destination <- file.path(paths$root, journal$session_name)
    if (!file.exists(session_destination)) {
        if (!file.exists(journal$staged_session) ||
            !isTRUE(renameFn(journal$staged_session, session_destination))) {
            cli::cli_abort("timestamped session publication failed",
                class = "saved_entity_write")
        }
    }
    if (!identical(digest::digest(file = session_destination, algo = "sha256"),
        journal$session_sha256)) {
        cli::cli_abort("timestamped session publication differs",
            class = "saved_entity_corrupt")
    }
    boundaryFn("after_timestamped_publication")
    for (name in names(journal$additional)) {
        spec <- journal$additional[[name]]
        destination <- file.path(paths$root, spec$destination)
        parent <- dirname(destination)
        if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
            cli::cli_abort("cannot create additional-output destination",
                class = "saved_entity_write")
        }
        if (!file.exists(destination)) {
            if (!file.exists(spec$staged) || !isTRUE(renameFn(spec$staged, destination))) {
                cli::cli_abort("additional-output publication failed",
                    class = "saved_entity_write")
            }
        }
        if (!identical(digest::digest(file = destination, algo = "sha256"),
            spec$copy_sha256)) {
            cli::cli_abort("published additional output differs",
                class = "saved_entity_corrupt")
        }
    }
    boundaryFn("before_latest_publication")
    latest_destination <- file.path(paths$root, "filtered_session_data_latest.rds")
    replacement <- file.path(paths$request_dir, "latest.replacement")
    if (file.exists(journal$staged_latest) &&
        !isTRUE(renameFn(journal$staged_latest, replacement))) {
        cli::cli_abort("latest session publication failed", class = "saved_entity_write")
    }
    # A retry after loss between the two renames finds only the replacement.
    if (file.exists(replacement) && !isTRUE(renameFn(replacement, latest_destination))) {
        cli::cli_abort("latest session publication failed", class = "saved_entity_write")
    }
    if (!file.exists(latest_destination) ||
        !identical(digest::digest(file = latest_destination, algo = "sha256"),
            journal$latest_sha256)) {
        cli::cli_abort("published latest session differs", class = "saved_entity_corrupt")
    }
    boundaryFn("before_completion_acknowledgement")
    journal$status <- "complete"
    journal$timestamped_path <- session_destination
    journal$latest_path <- latest_destination
    journal$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    syncWorkflowExportPaths(c(session_destination, latest_destination), paths$root)
    writeWorkflowExportJournal(paths, journal)
    completion <- file.path(paths$root, paste0("export_complete_",
        journal$request_id, ".json"))
    temporary <- tempfile(".export-complete-", tmpdir = paths$root)
    tryCatch(jsonlite::write_json(journal, temporary, auto_unbox = TRUE, pretty = TRUE, null = "null"),
        error = function(error) {
            unlink(temporary)
            stop(error)
        })
    syncWorkflowExportPaths(temporary, paths$root)
    if (!file.rename(temporary, completion)) {
        cli::cli_abort("export completion publication failed", class = "saved_entity_write")
    }
    syncWorkflowExportPaths(c(completion, paths$root), paths$root)
    journal
}

#' Serialize once or resume output-only publication for one request.
#' @param session_data Complete legacy session, required only before journal creation.
#' @param output_dir Existing output directory.
#' @param request_id Idempotent request token.
#' @param timestamp Fixed export timestamp.
#' @param additional_files Required owned copies.
#' @param serializeFn,copyFn,renameFn,boundaryFn Injectable test operations.
#' @return Complete publication receipt.
#' @noRd
runWorkflowLegacyExport <- function(session_data = NULL, output_dir, request_id,
    timestamp, additional_files = list(), serializeFn = saveRDS,
    copyFn = file.copy, renameFn = file.rename,
    boundaryFn = workflowExportBoundary) {
    paths <- workflowExportPaths(output_dir, request_id)
    journal_path <- file.path(paths$request_dir, "journal.json")
    if (!file.exists(journal_path)) {
        if (is.null(session_data)) {
            cli::cli_abort("session data is required before staging",
                class = "saved_entity_missing")
        }
        stageWorkflowLegacyExport(session_data, paths, request_id, timestamp,
            additional_files, serializeFn, copyFn, boundaryFn)
    }
    journal <- jsonlite::read_json(journal_path)
    if (is.null(journal$latest_sha256)) completeWorkflowLegacyCopy(paths, copyFn, boundaryFn)
    publishWorkflowLegacyExport(paths, boundaryFn, renameFn)
}

#' Fully validate an externally restored legacy session file.
#' @param path Existing session RDS.
#' @param expected_state_names Ordered expected states.
#' @return Decoded externally validated session.
#' @noRd
validateExternalWorkflowLegacySession <- function(path, expected_state_names) {
    if (!file.exists(path) || is.na(file.size(path)) || file.size(path) <= 0) {
        cli::cli_abort("external legacy session is missing", class = "saved_entity_missing")
    }
    session <- tryCatch(readRDS(path), error = function(error) {
        cli::cli_abort("external legacy session cannot be decoded",
            class = "saved_entity_corrupt", parent = error)
    })
    history <- if (is.list(session))
        unlist(session$r6_state_history, use.names = FALSE) else NULL
    if (!is.list(session) || !is.list(session$r6_complete_states) ||
        !identical(names(session$r6_complete_states), expected_state_names) ||
        !is.character(history) || !length(history) || anyNA(history) ||
        !all(history %in% expected_state_names)) {
        cli::cli_abort("external legacy session structure is invalid",
            class = "saved_entity_corrupt")
    }
    for (state in session$r6_complete_states) {
        if (!is.list(state) || isSavedWorkflowReference(state$data)) {
            cli::cli_abort("external legacy session is not self-contained",
                class = "saved_entity_corrupt")
        }
        if (isS4(state$data)) methods::validObject(state$data)
    }
    for (field in c("current_s4_object", "correlation_filtered_s4")) {
        if (isSavedWorkflowReference(session[[field]])) {
            cli::cli_abort("external legacy session contains a live reference",
                class = "saved_entity_corrupt")
        }
        if (isS4(session[[field]])) methods::validObject(session[[field]])
    }
    session
}

#' Resolve a committed Q-value result whose required figures are still pending.
#' @param manager Current component-state owner.
#' @param record Its already-read bounded checkpoint metadata.
#' @param hydrate FALSE checks pending status without decoding the display product.
#' @return Saved output inputs or NULL; never evaluates the scientific producer.
#' @noRd
resumeDIAQvalueOutput <- function(manager, record, hydrate = TRUE) {
    if (!manager$current_state %in% c("qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
        "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed", "protein_s4_created",
        "protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed",
        "protein_replicate_filtered")) return(NULL)
    selected_reference <- getDIAStateReference(manager)
    if (isTRUE(manager$.context$snapshot$completed_outputs[[selected_reference$entity_id]]$qc)) return(NULL)
    product <- if (manager$current_state == "qvalue_filtered") record$data$products$qvalue_progress else
        record$data$products$qc_progress
    reference <- getDIAStateReference(manager)
    product <- manager$.context$snapshot$output_products[[reference$entity_id]] %||% product
    root <- manager$.context$experiment_paths$time_dir
    if (is.null(product) || is.null(root)) return(NULL)
    path <- file.path(root, ".workflow-export", paste0("qvalue-", substr(reference$entity_id, 1L, 40L)),
        "journal.json")
    if (file.exists(path)) {
        if (file.size(path) > 65536) cli::cli_abort("oversized Q-value journal", class = "saved_entity_corrupt")
        journal <- jsonlite::read_json(path)
        if (!identical(journal$schema, "multischolar.dia_qc_publication.v1") ||
            !identical(journal$state_id, reference$entity_id) ||
            !identical(journal$incarnation, manager$.context$store$identity$incarnation)) {
            cli::cli_abort("Q-value publication identity differs", class = "saved_entity_scope")
        }
        if (identical(journal$status, "complete")) return(NULL)
    }
    if (!hydrate) return(TRUE)
    value <- hydrateSavedWorkflowEntity(manager$.context$store, product)
    checkSavedWorkflowConsumerPayload("scientific_product", value)
    c(list(stateReference = reference), value$value[setdiff(names(value$value), "stateReference")])
}

#' Publish required Q-value figures through the existing durable journal paths.
#' @param manager Current state owner.
#' @param result Accepted references and bounded display inputs.
#' @param arguments Existing filtering-plot arguments with concrete project paths.
#' @param render_fn Existing plot producer; numerical QC is already committed.
#' @return Bounded plot grid after publication acknowledgement.
#' @noRd
publishDIAQvalueOutput <- function(manager, result, arguments, render_fn) {
    reference <- result$stateReference
    context <- manager$.context
    check_owner <- function() {
        checkWorkflowEntityStore(context$store)
        if (!identical(getDIAStateReference(manager), reference)) {
            cli::cli_abort("Q-value publication was superseded", class = "saved_entity_scope")
        }
    }
    check_owner()
    key <- paste(arguments$omic_type, arguments$experiment_label, sep = "_")
    directories <- arguments$project_dirs[[key]]
    root <- directories$time_dir
    paths <- workflowExportPaths(root, paste0("qvalue-", substr(reference$entity_id, 1L, 40L)))
    journal_path <- file.path(paths$request_dir, "journal.json")
    journal <- if (file.exists(journal_path)) {
        if (file.size(journal_path) > 65536) cli::cli_abort("oversized Q-value journal", class = "saved_entity_corrupt")
        jsonlite::read_json(journal_path)
    } else {
        journal <- list(schema = "multischolar.dia_qc_publication.v1", status = "outputs_pending",
            request_id = basename(paths$request_dir), state_id = reference$entity_id,
            incarnation = context$store$identity$incarnation, files = NULL,
            scientific_checkpoint_count = 0L)
        writeWorkflowExportJournal(paths, journal)
        journal
    }
    if (!identical(journal$schema, "multischolar.dia_qc_publication.v1") ||
        !identical(journal$state_id, reference$entity_id) ||
        !identical(journal$incarnation, context$store$identity$incarnation)) {
        cli::cli_abort("Q-value publication identity differs", class = "saved_entity_scope")
    }
    was_complete <- identical(journal$status, "complete")
    staged <- file.path(paths$request_dir, "figures")
    grid_path <- file.path(paths$request_dir, "display.rds")
    expected <- as.vector(outer(paste0(arguments$step_name, "_", c("proteins_total",
        "proteins_per_run", "peptides_total", "peptides_per_protein", "peptides_per_run",
        "combined_plots")), c("png", "pdf"), paste, sep = "."))
    if (is.null(journal$files)) {
        if (!dir.exists(staged)) dir.create(staged)
        arguments$project_dirs[[key]]$time_dir <- staged
        previous <- if (exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE))
            get("filtering_progress", envir = .GlobalEnv) else NULL
        succeeded <- FALSE
        on.exit(if (!succeeded) {
            if (is.null(previous)) {
                if (exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE))
                    rm("filtering_progress", envir = .GlobalEnv)
            } else assign("filtering_progress", previous, envir = .GlobalEnv)
        }, add = TRUE)
        prior <- if (is.null(result$priorProgress)) methods::new("FilteringProgress") else
            do.call(methods::new, c(list(Class = "FilteringProgress"), result$priorProgress))
        assign("filtering_progress", prior, envir = .GlobalEnv)
        workflowExportBoundary("before_qvalue_figures")
        grid <- do.call(render_fn, arguments)
        files <- file.path(staged, expected)
        if (is.null(grid) || any(!file.exists(files)) || any(file.size(files) <= 0)) {
            cli::cli_abort("required Q-value figures are incomplete", class = "saved_entity_outputs_pending")
        }
        saveRDS(list(grid = grid, progress = get("filtering_progress", envir = .GlobalEnv)), grid_path)
        journal$files <- stats::setNames(lapply(files, \(path) list(
            sha256 = digest::digest(file = path, algo = "sha256"), bytes = file.size(path))), expected)
        journal$display_sha256 <- digest::digest(file = grid_path, algo = "sha256")
        syncWorkflowExportPaths(c(files, grid_path, staged), paths$root)
        writeWorkflowExportJournal(paths, journal)
        succeeded <- TRUE
    }
    if (!identical(sort(names(journal$files)), sort(expected)) ||
        any(grepl("[/\\\\]", names(journal$files)))) {
        cli::cli_abort("Q-value publication inventory differs", class = "saved_entity_corrupt")
    }
    previous_dir <- file.path(paths$request_dir, "previous")
    if (!dir.exists(previous_dir)) dir.create(previous_dir)
    for (name in names(journal$files)) {
        check_owner()
        source <- file.path(staged, name)
        destination <- file.path(paths$root, name)
        spec <- journal$files[[name]]
        valid <- file.exists(destination) && identical(digest::digest(file = destination, algo = "sha256"), spec$sha256)
        if (!valid) {
            if (!file.exists(source) || !identical(digest::digest(file = source, algo = "sha256"), spec$sha256)) {
                cli::cli_abort("staged Q-value figure differs", class = "saved_entity_corrupt")
            }
            backup <- file.path(previous_dir, name)
            if (file.exists(destination) && !file.exists(backup)) {
                if (!file.copy(destination, backup)) cli::cli_abort("cannot retain prior Q-value figure", class = "saved_entity_write")
                syncWorkflowExportPaths(c(backup, previous_dir), paths$root)
            }
            if (!file.rename(source, destination)) cli::cli_abort("cannot publish Q-value figure", class = "saved_entity_write")
            syncWorkflowExportPaths(c(destination, paths$root), paths$root)
        }
        workflowExportBoundary("after_qvalue_figure_publication")
    }
    check_owner()
    workflowExportBoundary("before_qvalue_acknowledgement")
    if (!file.exists(grid_path) || !identical(digest::digest(file = grid_path, algo = "sha256"), journal$display_sha256)) {
        cli::cli_abort("saved Q-value display differs", class = "saved_entity_corrupt")
    }
    display <- readRDS(grid_path)
    assign("filtering_progress", display$progress, envir = .GlobalEnv)
    journal$status <- "complete"
    writeWorkflowExportJournal(paths, journal)
    acknowledgeDIAOutputs(manager, reference, "qc")
    if (!was_complete) requestDIAReclamation(manager, reference, journal$request_id, "qc")
    display$grid
}
