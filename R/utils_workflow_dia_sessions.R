# Native DIA sessions reuse committed checkpoint files; they never serialize history.

#' Run a bounded native-export I/O pipeline, preserving every process failure.
#' @param commands List of argument vectors; paths are quoted, never shell code.
#' @param log_path Request-owned diagnostic file.
#' @return Small command output (a SHA-256 line or no output).
#' @noRd
runDIANativeIO <- function(commands, log_path) {
    shell <- unname(Sys.which("bash"))
    if (!identical(Sys.info()[["sysname"]], "Linux") || !nzchar(shell)) {
        cli::cli_abort("native export I/O requires bash on Linux", class = "saved_entity_unsupported")
    }
    pipeline <- vapply(commands, \(arguments) {
        command <- unname(Sys.which(arguments[[1L]]))
        if (!nzchar(command)) {
            cli::cli_abort("native export I/O command is unavailable", class = "saved_entity_unsupported")
        }
        paste(shQuote(c(command, arguments[-1L])), collapse = " ")
    }, character(1))
    # pipefail is essential: a successful hasher/writer cannot hide a failed reader.
    output <- suppressWarnings(system2(shell, c("--noprofile", "--norc", "-o", "pipefail",
        "-c", shQuote(paste(pipeline, collapse = " | "))), stdout = TRUE, stderr = log_path,
        env = c("LC_ALL=C", "TAR_OPTIONS=")))
    if (!file.exists(log_path) || file.size(log_path) > 67108864) {
        cli::cli_abort("native I/O log exceeds metadata limits", class = "saved_entity_limit")
    }
    if (!is.null(attr(output, "status"))) {
        cli::cli_abort(paste("native export I/O failed:",
            paste(utils::head(readLines(log_path, warn = FALSE), 10L), collapse = "\n")),
            class = "saved_entity_write")
    }
    output
}

#' Snapshot an owned file using copy-on-write when the filesystem supports it.
#' @param source Validated regular input file.
#' @param destination Fresh path in the caller's private staging directory.
#' @param log_path Request-owned bounded command log.
#' @return NULL invisibly; callers independently verify content before publication.
#' @noRd
snapshotWorkflowFile <- function(source, destination, log_path) {
    before <- workflowTableFileIdentity(source)
    destination_link <- Sys.readlink(destination)
    if (file.exists(destination) || !is.na(destination_link) && nzchar(destination_link)) {
        cli::cli_abort("snapshot destination exists", class = "saved_entity_conflict")
    }
    # GNU cp falls back to a regular copy; a reflink is independent of subsequent
    # source writes, unlike a hard link. Existing publication hash checks remain.
    runDIANativeIO(list(c("cp", "--reflink=auto", "--no-clobber", "--",
        source, destination)), log_path)
    if (!identical(before, workflowTableFileIdentity(source))) {
        cli::cli_abort("snapshot source changed", class = "saved_entity_corrupt")
    }
    syncWorkflowExportPaths(destination, dirname(destination))
    invisible(NULL)
}

#' Release clean cached pages of explicit validated export inputs, not their files.
#' @param paths Regular files already validated by the export owner.
#' @param log_path Request-owned diagnostic file.
#' @return NULL invisibly; no global cache operation or file deletion.
#' @noRd
releaseDIANativeFileCache <- function(paths, log_path) {
    releaseWorkflowFileCache(paths, log_path)
}

#' Release clean pages of explicit owned files across workflow adapters.
#' @param paths Closed regular files; never a global cache target.
#' @param log_path Existing owner-controlled bounded I/O log.
#' @return NULL invisibly.
#' @noRd
releaseWorkflowFileCache <- function(paths, log_path) {
    for (path in unique(paths)) {
        info <- file.info(path)
        if (is.na(info$size) || isTRUE(info$isdir) || nzchar(Sys.readlink(path))) {
            cli::cli_abort("invalid native cache-release file", class = "saved_entity_path")
        }
        if (info$size < 4194304) next
        runDIANativeIO(list(c("dd", paste0("if=", path), "iflag=nocache", "count=0",
            "status=none")), log_path)
    }
    invisible(NULL)
}

#' Run the supported streaming archive tool with bounded diagnostic output.
#' @param arguments Fixed tar arguments, each quoted separately.
#' @param log_path Request-owned command log.
#' @param output_path Optional owned archive written through a bounded cache-dropping stream.
#' @return Command output, bounded by the existing metadata ceiling.
#' @noRd
runDIANativeTar <- function(arguments, log_path, output_path = NULL) {
    command <- unname(Sys.which("tar"))
    if (!identical(Sys.info()[["sysname"]], "Linux") || !nzchar(command)) {
        cli::cli_abort("native DIA archives require GNU tar on Linux", class = "saved_entity_unsupported")
    }
    version <- system2(command, "--version", stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(version, "status")) || !any(grepl("GNU tar", version, fixed = TRUE))) {
        cli::cli_abort("native DIA archives require GNU tar", class = "saved_entity_unsupported")
    }
    if (!is.null(output_path)) {
        # Synchronous 4 MiB writes make pages clean before dd advises eviction.
        # This is an uncompressed byte stream: the native PAX format is unchanged.
        return(runDIANativeIO(list(c(command, arguments), c("dd", paste0("of=", output_path),
            "bs=4M", "iflag=fullblock", "oflag=dsync,nocache", "status=none")), log_path))
    }
    status <- system2(command, shQuote(arguments), stdout = log_path, stderr = log_path,
        env = c("LC_ALL=C", "TAR_OPTIONS="))
    if (!file.exists(log_path) || file.size(log_path) > 67108864) {
        cli::cli_abort("native archive log exceeds metadata limits", class = "saved_entity_limit")
    }
    output <- readLines(log_path, warn = FALSE)
    if (!identical(status, 0L)) {
        cli::cli_abort(paste("native archive command failed:",
            paste(utils::head(output, 10L), collapse = "\n")), class = "saved_entity_write")
    }
    output
}

#' Accept only the native archive's fixed regular-file layout.
#' @param names Archive-relative member paths, never shell expressions.
#' @return TRUE invisibly; aliases/traversal/duplicates are rejected.
#' @noRd
checkDIANativeMemberNames <- function(names) {
    pattern <- paste0("^(manifest\\.json|session-references\\.rds|",
        "\\.workflow-entities/(store\\.json|dia-state-head\\.rds|",
        "entities/[A-Za-z0-9_-]{1,80}/(commit\\.json|descriptor\\.rds|",
        "payload\\.(rds|qs2|parquet)))|files/[a-f0-9]{64})$")
    if (!is.character(names) || !length(names) || length(names) > 40004L ||
        anyNA(names) || anyDuplicated(names) || any(!grepl(pattern, names))) {
        cli::cli_abort("invalid native archive member paths", class = "saved_entity_path")
    }
    invisible(TRUE)
}

#' Inspect every archive member before extraction into a private directory.
#' @param path Native uncompressed archive.
#' @param directory Existing request-owned inspection directory.
#' @return Member names and declared sizes; no payload is decoded.
#' @noRd
inspectDIANativeArchive <- function(path, directory) {
    info <- file.info(path)
    if (!file.exists(path) || isTRUE(info$isdir) || is.na(info$size) || info$size < 1024) {
        cli::cli_abort("native archive is missing or incomplete", class = "saved_entity_missing")
    }
    lines <- runDIANativeTar(c("--list", "--verbose", "--full-time", "--numeric-owner",
        "--quoting-style=escape", "--force-local", "--file", path),
        file.path(directory, "archive-list.log"))
    # Only ordinary files are supported: links, devices and directories are not
    # emitted by our writer and cannot be smuggled into a project on import.
    pattern <- paste0("^-[rwxStTs-]{9}[[:space:]]+[0-9]+/[0-9]+[[:space:]]+",
        "([0-9]+)[[:space:]]+[0-9]{4}-[0-9]{2}-[0-9]{2}[[:space:]]+",
        "[0-9:.]+[[:space:]]+(.+)$")
    matches <- regmatches(lines, regexec(pattern, lines))
    if (!length(matches) || any(lengths(matches) != 3L)) {
        cli::cli_abort("unsupported native archive entry type", class = "saved_entity_path")
    }
    names <- vapply(matches, `[[`, character(1), 3L)
    sizes <- as.numeric(vapply(matches, `[[`, character(1), 2L))
    checkDIANativeMemberNames(names)
    if (any(!is.finite(sizes) | sizes < 0) || sum(sizes) > info$size ||
        !all(c("manifest.json", "session-references.rds", ".workflow-entities/store.json",
            ".workflow-entities/dia-state-head.rds") %in% names) ||
        any(sizes[names %in% c("manifest.json", "session-references.rds")] > 67108864)) {
        cli::cli_abort("native archive sizes or required entries differ", class = "saved_entity_corrupt")
    }
    stats::setNames(sizes, names)
}

#' Collect explicit descriptor roots without visiting any scientific payload.
#' @param value Snapshot/reference metadata only.
#' @param references Existing identity-keyed descriptor list.
#' @return Distinct descriptors, preserving aliases by immutable identity.
#' @noRd
collectDIANativeReferences <- function(value, references = list()) {
    found <- new.env(hash = TRUE, parent = emptyenv())
    count <- 0L
    visit <- function(node, depth = 0L) {
        if (depth > 16L) cli::cli_abort("native reference nesting exceeds limits", class = "saved_entity_limit")
        if (isSavedWorkflowReference(node)) {
            previous <- get0(node$entity_id, found, inherits = FALSE)
            if (!is.null(previous) && !identical(previous, node)) {
                cli::cli_abort("conflicting native session references", class = "saved_entity_conflict")
            }
            if (is.null(previous)) count <<- count + 1L
            if (count > 10000L) cli::cli_abort("native references exceed registry limits", class = "saved_entity_limit")
            assign(node$entity_id, node, found)
        } else if (is.list(node) && !is.object(node)) {
            for (item in node) visit(item, depth + 1L)
        } else if (is.object(node) || !is.atomic(node) && !is.null(node)) {
            cli::cli_abort("native reference roots contain a live payload", class = "saved_entity_invalid")
        }
    }
    visit(references)
    visit(value)
    references <- as.list(found, all.names = TRUE)
    if (as.numeric(object.size(references)) > 67108864) {
        cli::cli_abort("native references exceed registry limits", class = "saved_entity_limit")
    }
    references
}

#' Bind one regular owned file's bytes and fingerprint to an archive path.
#' @param path Existing file.
#' @param log_path Optional export-owned log; enables drop-behind streaming hashing.
#' @return Bounded inventory entry.
#' @noRd
describeDIANativeFile <- function(path, log_path = NULL) {
    identity <- workflowTableFileIdentity(path)
    link <- Sys.readlink(path)
    if (isTRUE(file.info(path)$isdir) || !is.na(link) && nzchar(link)) {
        cli::cli_abort("native archive sources must be regular owned files", class = "saved_entity_path")
    }
    hash <- hashWorkflowFile(path, log_path)
    if (!identical(identity, workflowTableFileIdentity(path))) {
        cli::cli_abort("native archive source changed during validation", class = "saved_entity_corrupt")
    }
    list(bytes = as.numeric(identity$size), sha256 = hash)
}

#' Build one native archive from a frozen head and committed immutable files.
#' @param manager Current DIA adapter.
#' @param metadata_reference Previously persisted metadata-only session value.
#' @param extra_references Required GUI-owned references outside the history head.
#' @param dependencies Named source/relative-destination file specifications.
#' @param paths Existing export request paths.
#' @param request Fixed request identity, head and timestamp.
#' @return Durable outputs-pending journal; no scientific object is hydrated.
#' @noRd
stageDIANativeSession <- function(manager, metadata_reference, extra_references,
    dependencies, paths, request) {
    context <- manager$.context
    store <- context$store
    checkWorkflowEntityStore(store)
    if (store$open_resources != 0L || !identical(context$head$entity_id, request$generation)) {
        cli::cli_abort("native export head changed or resources remain open", class = "saved_entity_scope")
    }
    references <- expandWorkflowTableReferences(store, collectDIANativeReferences(list(context$snapshot, context$head,
        metadata_reference, extra_references)))
    inventory <- list()
    store_members <- character()
    identities <- list()
    io_log <- file.path(paths$request_dir, "native-io.log")
    for (reference in references) {
        validateSavedWorkflowReference(store, reference)
        releaseDIANativeFileCache(file.path(store$directory, reference$payload$relative_path), io_log)
        relative <- c(file.path(dirname(reference$payload$relative_path), c("descriptor.rds", "commit.json")),
            reference$payload$relative_path)
        for (member in relative) {
            path <- file.path(store$directory, member)
            archive_name <- paste0(".workflow-entities/", member)
            # Payload validation already computed this fingerprint. Descriptor
            # and commit files are small, independent archived files.
            inventory[[archive_name]] <- if (identical(member, reference$payload$relative_path))
                list(bytes = reference$payload$bytes, sha256 = reference$payload$sha256) else
                    describeDIANativeFile(path)
            identities[[path]] <- workflowTableFileIdentity(path)
            store_members <- c(store_members, archive_name)
        }
    }
    staged_store <- file.path(paths$request_dir, ".workflow-entities")
    if (!dir.exists(staged_store) && !dir.create(staged_store)) {
        cli::cli_abort("cannot stage native session metadata", class = "saved_entity_write")
    }
    jsonlite::write_json(store$identity, file.path(staged_store, "store.json"), auto_unbox = TRUE)
    saveRDS(context$head, file.path(staged_store, "dia-state-head.rds"), version = 3)
    roots <- list(metadata = metadata_reference, extra = extra_references)
    saveRDS(roots, file.path(paths$request_dir, "session-references.rds"), version = 3)
    stage_members <- c(".workflow-entities/store.json", ".workflow-entities/dia-state-head.rds",
        "session-references.rds")
    if (!is.list(dependencies) || length(dependencies) > store$max_entries ||
        length(dependencies) && (is.null(names(dependencies)) || anyDuplicated(names(dependencies)))) {
        cli::cli_abort("invalid native dependency map", class = "saved_entity_invalid")
    }
    dependency_map <- list()
    for (name in names(dependencies)) {
        spec <- dependencies[[name]]
        checkSavedEntityString(name, "dependency name")
        checkSavedEntityPath(store, spec$destination)
        description <- describeDIANativeFile(spec$source, io_log)
        if (!is.null(spec$sha256) && (!identical(description$sha256, spec$sha256) ||
            !isTRUE(description$bytes == spec$bytes))) {
            cli::cli_abort("scientific execution dependency changed", class = "saved_entity_corrupt")
        }
        member <- paste0("files/", description$sha256)
        destination <- file.path(paths$request_dir, member)
        if (!dir.exists(dirname(destination))) dir.create(dirname(destination))
        if (!file.exists(destination)) {
            snapshotWorkflowFile(spec$source, destination, io_log)
        }
        if (!identical(describeDIANativeFile(destination, io_log), description)) {
            cli::cli_abort("native dependency copy differs", class = "saved_entity_corrupt")
        }
        inventory[[member]] <- description
        identities[[destination]] <- workflowTableFileIdentity(destination)
        dependency_map[[name]] <- list(member = member, destination = spec$destination,
            source_provenance = spec$source_provenance %||% spec$source)
        stage_members <- union(stage_members, member)
    }
    for (member in setdiff(stage_members, names(inventory))) inventory[[member]] <-
        describeDIANativeFile(file.path(paths$request_dir, member), io_log)
    checkDIANativeMemberNames(c("manifest.json", names(inventory)))
    # KIT-020 audit F-07: the store's workflow identity is the one source of the manifest's type (the manager's
    # own field is the fresh default until an import sets it, and is not consulted).
    manifest <- list(schema = if (usesDIAReferenceStates(manager)) "multischolar.dia_native_session.v3" else
        "multischolar.dia_native_session.v1", request_id = request$request_id,
        generation = request$generation, workflow = store$identity$workflow_id,
        route = workflowRouteForType(store$identity$workflow_id), identity = store$identity,
        software = as.character(utils::packageVersion("MultiScholaR")),
        dependencies = dependency_map, files = inventory)
    jsonlite::write_json(manifest, file.path(paths$request_dir, "manifest.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    writeBin(store_members, file.path(paths$request_dir, "store-members.nul"))
    writeBin(c("manifest.json", stage_members), file.path(paths$request_dir, "stage-members.nul"))
    archive <- file.path(paths$request_dir, "session.msdia.tar")
    workflowExportBoundary("before_native_archive")
    runDIANativeTar(c("--create", "--format=pax", "--no-recursion", "--hard-dereference",
        "--force-local", "--file", "-", "--null", "--directory", dirname(store$directory),
        "--files-from", file.path(paths$request_dir, "store-members.nul"),
        "--directory", paths$request_dir,
        "--files-from", file.path(paths$request_dir, "stage-members.nul")),
        file.path(paths$request_dir, "archive-create.log"), output_path = archive)
    for (path in names(identities)) if (!identical(identities[[path]], workflowTableFileIdentity(path))) {
        cli::cli_abort("native source changed during packaging", class = "saved_entity_corrupt")
    }
    releaseDIANativeFileCache(names(identities), io_log)
    checkWorkflowEntityStore(store)
    if (!identical(context$head$entity_id, request$generation)) {
        cli::cli_abort("native export generation changed", class = "saved_entity_scope")
    }
    listed <- inspectDIANativeArchive(archive, paths$request_dir)
    if (!setequal(names(listed), c("manifest.json", names(inventory)))) {
        cli::cli_abort("native archive inventory differs", class = "saved_entity_corrupt")
    }
    description <- describeDIANativeFile(archive, io_log)
    journal <- request
    journal$archive_bytes <- description$bytes
    journal$archive_sha256 <- description$sha256
    journal$staged_archive <- archive
    syncWorkflowExportPaths(c(archive, paths$request_dir), paths$root)
    writeWorkflowExportJournal(paths, journal)
    workflowExportBoundary("after_native_archive")
    journal
}

#' Read a bounded journal for native-only publication/recovery.
#' @param paths Request paths produced by workflowExportPaths.
#' @return Validated outputs-pending or complete native journal.
#' @noRd
readDIANativeJournal <- function(paths) {
    path <- file.path(paths$request_dir, "journal.json")
    if (!file.exists(path) || file.size(path) > 67108864) {
        cli::cli_abort("native export journal is missing or oversized", class = "saved_entity_corrupt")
    }
    journal <- jsonlite::read_json(path)
    if (!identical(journal$schema, "multischolar.dia_native_export.v1") ||
        !identical(journal$request_id, basename(paths$request_dir)) ||
        !journal$status %in% c("outputs_pending", "complete") ||
        !is.character(journal$generation) || !grepl("^[a-f0-9]{64}$", journal$generation) ||
        !identical(journal$archive_name, paste0("filtered_session_", journal$timestamp, "_",
            journal$request_id, ".msdia.tar"))) {
        cli::cli_abort("native export journal differs from request", class = "saved_entity_corrupt")
    }
    journal
}

#' Publish one immutable native archive and its small latest pointer last.
#' @param paths Existing request paths.
#' @param owned_staged TRUE only immediately after staging in the same call.
#' @return Completed native publication receipt.
#' @noRd
publishDIANativeSession <- function(paths, owned_staged = FALSE) {
    journal <- readDIANativeJournal(paths)
    staged <- file.path(paths$request_dir, "session.msdia.tar")
    destination <- file.path(paths$root, journal$archive_name)
    io_log <- file.path(paths$request_dir, "native-io.log")
    if (is.null(journal$archive_sha256) || is.null(journal$archive_bytes)) {
        cli::cli_abort("native archive packaging is pending", class = "saved_entity_outputs_pending")
    }
    expected <- list(bytes = as.numeric(journal$archive_bytes), sha256 = journal$archive_sha256)
    if (!file.exists(destination)) {
        if (!file.exists(staged) || !isTRUE(as.numeric(file.size(staged)) == expected$bytes) ||
            !owned_staged && !identical(describeDIANativeFile(staged, io_log), expected)) {
            cli::cli_abort("staged native archive differs", class = "saved_entity_corrupt")
        }
        if (!file.rename(staged, destination)) {
            cli::cli_abort("native archive publication failed", class = "saved_entity_write")
        }
    } else if (!identical(describeDIANativeFile(destination, io_log), expected)) {
        cli::cli_abort("published native archive differs", class = "saved_entity_corrupt")
    }
    syncWorkflowExportPaths(c(destination, paths$root), paths$root)
    workflowExportBoundary("after_native_archive_publication")
    pointer <- list(schema = "multischolar.dia_native_latest.v1", status = "complete",
        request_id = journal$request_id, generation = journal$generation,
        archive = journal$archive_name, bytes = expected$bytes, sha256 = expected$sha256)
    completion <- file.path(paths$root, paste0("native_export_complete_", journal$request_id, ".json"))
    temporary <- file.path(paths$request_dir, "completion.json")
    jsonlite::write_json(pointer, temporary, auto_unbox = TRUE, pretty = TRUE)
    syncWorkflowExportPaths(temporary, paths$root)
    if (!file.rename(temporary, completion)) {
        cli::cli_abort("native completion publication failed", class = "saved_entity_write")
    }
    workflowExportBoundary("before_native_latest_publication")
    latest <- file.path(paths$root, "filtered_session_latest.json")
    temporary <- file.path(paths$request_dir, "latest.json")
    jsonlite::write_json(pointer, temporary, auto_unbox = TRUE, pretty = TRUE)
    syncWorkflowExportPaths(c(temporary, completion, paths$root), paths$root)
    if (!file.rename(temporary, latest)) {
        cli::cli_abort("native latest publication failed", class = "saved_entity_write")
    }
    syncWorkflowExportPaths(c(latest, paths$root), paths$root)
    workflowExportBoundary("before_native_completion_acknowledgement")
    journal$status <- "complete"
    journal$archive_path <- destination
    journal$latest_path <- latest
    journal$completion <- completion
    writeWorkflowExportJournal(paths, journal)
    journal
}

#' Find one interrupted publication for the unchanged live generation.
#' @param manager Current adapter.
#' @param output_dir Existing source/export directory.
#' @return Matching journal, or NULL; ambiguous recovery is refused.
#' @noRd
findDIANativePendingExport <- function(manager, output_dir) {
    root <- file.path(output_dir, ".workflow-export")
    if (!dir.exists(root)) return(NULL)
    entries <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    if (length(entries) > 10000L) cli::cli_abort("export recovery index exceeds limits", class = "saved_entity_limit")
    matches <- list()
    for (directory in entries) {
        path <- file.path(directory, "journal.json")
        if (!file.exists(path)) next
        if (file.size(path) > 67108864) cli::cli_abort("export journal exceeds limits", class = "saved_entity_limit")
        journal <- jsonlite::read_json(path)
        if (identical(journal$schema, "multischolar.dia_native_export.v1") &&
            identical(journal$status, "outputs_pending") &&
            identical(journal$generation, manager$.context$head$entity_id) &&
            identical(journal$store_incarnation, manager$.context$store$identity$incarnation)) {
            matches[[length(matches) + 1L]] <- journal
        }
    }
    if (length(matches) > 1L) cli::cli_abort("ambiguous pending native exports", class = "saved_entity_conflict")
    if (length(matches)) matches[[1L]] else NULL
}

#' Save or resume a native DIA session without assembling its history.
#' @param manager Live DIA adapter.
#' @param metadata Explicit non-history session metadata; NULL on recovery.
#' @param output_dir Existing export destination.
#' @param extra_references Named additional GUI-owned descriptors.
#' @param dependencies Named source/destination file dependencies.
#' @param prepare_outputs Optional existing metadata-output callback.
#' @return Native publication receipt and metadata reference only.
#' @noRd
runDIANativeSessionExport <- function(manager, metadata, output_dir,
    extra_references = list(), dependencies = list(), prepare_outputs = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort("native export requires a DIA state owner", class = "saved_entity_scope")
    }
    assertDIAExportReady(manager, allow_action = "session-export")
    context <- manager$.context
    checkWorkflowEntityStore(context$store)
    root <- normalizePath(output_dir, mustWork = TRUE)
    request <- context$export_pending
    if (is.null(request)) request <- findDIANativePendingExport(manager, root)
    if (!is.null(request) && (!identical(request$schema, "multischolar.dia_native_export.v1") ||
        !identical(request$generation, context$head$entity_id) || !identical(request$root, root))) {
        cli::cli_abort("another export generation is pending", class = "saved_entity_outputs_pending")
    }
    if (is.null(request)) {
        if (!is.list(metadata) || any(c("r6_complete_states", "current_s4_object",
            "correlation_filtered_s4", "state_manager") %in% names(metadata))) {
            cli::cli_abort("native session metadata contains legacy history fields", class = "saved_entity_invalid")
        }
        if (!is.null(metadata$config_list)) metadata$config_list <-
            prepareDIAStateRecord(list(data = NULL, config = metadata$config_list))$config
        reference <- saveDIAStateValue(context, metadata, "native-session-metadata", codec = "rds")
        request <- list(schema = "multischolar.dia_native_export.v1", status = "outputs_pending",
            request_id = basename(tempfile("native-export-")), generation = context$head$entity_id,
            root = root, store_incarnation = context$store$identity$incarnation,
            timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"), scientific_checkpoint_count = 0L)
        request$archive_name <- paste0("filtered_session_", request$timestamp, "_", request$request_id, ".msdia.tar")
        paths <- workflowExportPaths(root, request$request_id)
        # Keep descriptors in RDS so factors/integers retain their exact types;
        # the recovery journal contains only scalar control metadata.
        saveRDS(list(metadata = reference, extra = extra_references, dependencies = dependencies),
            file.path(paths$request_dir, "request-references.rds"), version = 3)
        request$references_sha256 <- digest::digest(file = file.path(paths$request_dir,
            "request-references.rds"), algo = "sha256")
        syncWorkflowExportPaths(file.path(paths$request_dir, "request-references.rds"), root)
        writeWorkflowExportJournal(paths, request)
    } else paths <- workflowExportPaths(root, request$request_id)
    context$export_pending <- request
    request <- readDIANativeJournal(paths)
    reference_path <- file.path(paths$request_dir, "request-references.rds")
    if (!file.exists(reference_path) || file.size(reference_path) > 67108864 ||
        !identical(digest::digest(file = reference_path, algo = "sha256"), request$references_sha256)) {
        cli::cli_abort("native recovery references differ", class = "saved_entity_corrupt")
    }
    roots <- readRDS(reference_path)
    if (is.function(prepare_outputs)) {
        prepare_outputs(resolveDIAStateValue(context, roots$metadata), request)
    }
    freshly_staged <- is.null(request$archive_sha256)
    if (freshly_staged) stageDIANativeSession(manager, roots$metadata, roots$extra,
        roots$dependencies, paths, request)
    publication <- publishDIANativeSession(paths, owned_staged = freshly_staged)
    context$owned_export <- c(publication, list(format = "native_v1", identity = context$store$identity,
        metadata_reference = roots$metadata, extra_references = roots$extra,
        archive_identity = workflowTableFileIdentity(publication$archive_path)))
    context$export_pending <- NULL
    list(publication = publication, metadata_reference = roots$metadata)
}

#' Convert old whole-state records at the explicit session-import boundary.
#' @param manager Imported native owner; the previous live project is still intact.
#' @return TRUE after one-record-at-a-time conversion and atomic head publication.
#' @noRd
convertDIANativeLegacyCheckpoints <- function(manager) {
    if (usesDIAReferenceStates(manager)) return(invisible(TRUE))
    context <- manager$.context
    snapshot <- context$snapshot
    for (name in names(snapshot$states)) {
        record <- materializeDIAStateRecord(context, snapshot$states[[name]])
        if (!is.null(record$data)) {
            if (!isS4(record$data)) cli::cli_abort("legacy checkpoint is not a scientific state",
                class = "saved_entity_corrupt")
            checkpoint <- checkpointDIAExecutionResult(manager, record$data)
            record$data <- checkpoint
            snapshot$states[[name]] <- saveDIAStateValue(context, record, name, codec = "rds")
            snapshot$component_dependencies[[name]] <- collectDIANativeReferences(list(
                checkpoint$components, checkpoint$support, checkpoint$products, checkpoint$dependencies))
            snapshot$sizes[[name]] <- list(rows = checkpoint$dimensions$rows,
                columns = length(checkpoint$dimensions$columns), note = "converted component dimensions")
            if (!is.null(checkpoint$metadata)) snapshot$state_metadata[[name]] <- checkpoint$metadata
        }
        record <- NULL
        context$working <- new.env(parent = emptyenv())
        gc(full = TRUE)
    }
    snapshot$reference_states <- TRUE
    publishDIAStateHead(context, snapshot)
    context$reference_states <- TRUE
    invisible(TRUE)
}

#' Open a validated native archive in a new private project directory.
#' @param path User-selected native archive; never an owned-generation shortcut.
#' @param project_root Existing destination project directory.
#' @return New manager, metadata and dependency map; previous session is untouched.
#' @noRd
openDIANativeSession <- function(path, project_root) {
    path <- normalizePath(path, mustWork = TRUE)
    project_root <- normalizePath(project_root, mustWork = TRUE)
    directory <- tempfile(".dia-native-import-", tmpdir = project_root)
    if (!dir.create(directory, mode = "0700")) {
        cli::cli_abort("cannot stage native import", class = "saved_entity_write")
    }
    listed <- inspectDIANativeArchive(path, directory)
    # Inspection forbids links/traversal before any extraction. Extraction is
    # confined to a newly created private directory, never an existing project.
    runDIANativeTar(c("--extract", "--force-local", "--no-same-owner", "--no-same-permissions",
        "--keep-old-files", "--file", path, "--directory", directory),
        file.path(directory, "archive-extract.log"))
    manifest <- jsonlite::read_json(file.path(directory, "manifest.json"))
    if (!manifest$schema %in% c("multischolar.dia_native_session.v1", "multischolar.dia_native_session.v2", "multischolar.dia_native_session.v3") ||
        !is.character(manifest$workflow) || !manifest$workflow %in% names(workflowRoutesByType) ||
        !identical(manifest$route, workflowRoutesByType[[manifest$workflow]]) ||
        !is.list(manifest$files) || !setequal(names(listed), c("manifest.json", names(manifest$files)))) {
        cli::cli_abort("native session manifest is invalid", class = "saved_entity_corrupt")
    }
    for (member in names(manifest$files)) {
        actual <- describeDIANativeFile(file.path(directory, member))
        expected <- manifest$files[[member]]
        if (!identical(actual$sha256, expected$sha256) ||
            !isTRUE(actual$bytes == expected$bytes) || !isTRUE(actual$bytes == listed[[member]])) {
            cli::cli_abort("native session file integrity differs", class = "saved_entity_corrupt")
        }
    }
    identity <- readSavedStoreIdentity(file.path(directory, ".workflow-entities", "store.json"))
    if (!identical(identity, manifest$identity) || !identical(manifest$workflow, identity$workflow_id)) {
        cli::cli_abort("native session store identity differs", class = "saved_entity_scope")
    }
    store <- newWorkflowEntityStore(directory, identity$project_id, identity$workflow_id, max_cache_bytes = 0)
    accepted <- FALSE
    on.exit(if (!accepted) closeWorkflowEntityStore(store), add = TRUE)
    descriptor_paths <- grep("/descriptor[.]rds$", names(manifest$files), value = TRUE)
    for (member in descriptor_paths) {
        descriptor_path <- file.path(directory, member)
        commit_path <- file.path(dirname(descriptor_path), "commit.json")
        if (file.size(descriptor_path) > 65536 || !file.exists(commit_path) || file.size(commit_path) > 8192) {
            cli::cli_abort("native descriptor bounds differ", class = "saved_entity_corrupt")
        }
        descriptor <- readRDS(descriptor_path)
        commit <- jsonlite::read_json(commit_path)
        expected <- newSavedWorkflowDescriptor(store, descriptor$kind, descriptor$identity,
            descriptor$payload, descriptor$metadata, descriptor$parent_id, descriptor$revision)
        payload_member <- paste0(".workflow-entities/", descriptor$payload$relative_path)
        inventory <- manifest$files[[payload_member]]
        if (!identical(descriptor, expected) || !identical(commit$schema, 1L) ||
            !identical(commit$entity_id, descriptor$entity_id) ||
            !identical(commit$descriptor_sha256, manifest$files[[member]]$sha256) ||
            is.null(inventory) || !identical(inventory$sha256, descriptor$payload$sha256) ||
            !isTRUE(inventory$bytes == descriptor$payload$bytes)) {
            cli::cli_abort("native committed descriptor differs", class = "saved_entity_corrupt")
        }
        registerSavedWorkflowMetadata(store, descriptor)
    }
    roots <- readRDS(file.path(directory, "session-references.rds"))
    manager <- newDIAStateAdapter(store, reopen = TRUE)
    if (manifest$schema %in% c("multischolar.dia_native_session.v2", "multischolar.dia_native_session.v3") &&
        !usesDIAReferenceStates(manager)) {
        cli::cli_abort("component archive has no component state head", class = "saved_entity_corrupt")
    }
    if (!identical(manager$.context$head$entity_id, manifest$generation)) {
        cli::cli_abort("native session generation differs", class = "saved_entity_scope")
    }
    references <- expandWorkflowTableReferences(store, collectDIANativeReferences(list(manager$.context$snapshot, manager$.context$head, roots)))
    for (reference in references) if (!identical(reference,
        get0(reference$entity_id, envir = store$records, inherits = FALSE))) {
        cli::cli_abort("native session reference graph is incomplete", class = "saved_entity_corrupt")
    }
    if (manifest$schema %in% c("multischolar.dia_native_session.v2", "multischolar.dia_native_session.v3")) {
        for (state in setdiff(names(manager$.context$snapshot$states), "initial")) {
            record <- readDIAStateMetadata(manager, getDIAStateReference(manager, state))
            children <- expandWorkflowTableReferences(store, collectDIANativeReferences(list(record$data$components, record$data$support, record$data$products, record$data$dependencies)))
            for (child in children) if (!identical(child,
                get0(child$entity_id, envir = store$records, inherits = FALSE))) {
                cli::cli_abort("component checkpoint dependency is absent", class = "saved_entity_corrupt")
            }
        }
    }
    metadata <- resolveDIAStateValue(manager$.context, roots$metadata)
    if (!is.list(metadata) || any(c("r6_complete_states", "current_s4_object",
        "correlation_filtered_s4", "state_manager") %in% names(metadata))) {
        cli::cli_abort("native metadata contains legacy history", class = "saved_entity_corrupt")
    }
    dependencies <- lapply(manifest$dependencies, function(dependency) {
        checkDIANativeMemberNames(dependency$member)
        checkSavedEntityPath(store, dependency$destination)
        if (!dependency$member %in% names(manifest$files)) {
            cli::cli_abort("native dependency is absent", class = "saved_entity_corrupt")
        }
        path <- file.path(directory, dependency$member)
        c(dependency, list(path = path, file_identity = workflowTableFileIdentity(path)))
    })
    manager$.context$sql_rollup <- TRUE
    manager$.context$session_format <- "native_v1"
    manager$.context$native_dependencies <- dependencies
    manager$.context$native_metadata_reference <- roots$metadata
    manager$.context$native_extra_references <- roots$extra
    manager$.context$owned_export <- manager$.context$export_pending <- NULL
    # All extracted files and directory entries become durable before a GUI can
    # replace its previous manager with this one.
    durable <- c(file.path(directory, names(manifest$files)), file.path(directory, "manifest.json"),
        list.dirs(directory, recursive = TRUE, full.names = TRUE))
    for (chunk in split(durable, ceiling(seq_along(durable) / 128L))) {
        syncWorkflowExportPaths(chunk, project_root)
    }
    syncWorkflowExportPaths(c(directory, project_root), project_root)
    workflowExportBoundary("before_native_import_adoption")
    # KIT-020 audit F-07: the store's workflow identity is authoritative for the reopened owner's type (the
    # saved field may hold the fresh default of a manager whose import never set it); set after the generation
    # check, which compares the archive's head with the store's.
    if (!identical(manager$workflow_type, identity$workflow_id)) manager$setWorkflowType(identity$workflow_id)
    accepted <- TRUE
    list(manager = manager, metadata = metadata, extra_references = roots$extra,
        dependencies = dependencies, directory = directory, manifest = manifest)
}

#' The on-disk native DIA path is the production default (KIT-011, user decision 2026-10-11): no
#' option gate remains. The platform limits stay in `initializeDIANativeImport` (Linux; the DIA-NN IQ route and
#' the protein-level formats, KIT-020).
#' @return TRUE.
#' @noRd
isDIANativeDefaultEnabled <- function() TRUE

#' The proteomics workflow type an import format produces (KIT-020): DIA-NN precursors are `DIA`, Proteome
#' Discoverer TMT abundances are `TMT`, the other protein-level readers are `LFQ`; an unknown format gives NA
#' and never attaches the on-disk owner. The owner serves peptide-level DIA-NN tables and protein-level tables
#' of the other readers only (audit F-11): a Spectronaut peptide export (`quantity_type = "pep"`) gives NA.
#' @param format Detected import format.
#' @param data_type The reader's `data_type` (`protein` or `peptide`) when known; NULL before the file is read.
#' @return One of `DIA`, `TMT`, `LFQ` or NA.
#' @noRd
workflowTypeForImportFormat <- function(format, data_type = NULL) {
    if (!is.character(format) || length(format) != 1L) return(NA_character_)
    type <- switch(format, diann = "DIA", pd_tmt = "TMT", maxquant = , fragpipe = , spectronaut = "LFQ", NA_character_)
    if (is.na(type) || is.null(data_type)) return(type)
    expected <- if (identical(format, "diann")) "peptide" else "protein"
    if (identical(data_type, expected)) type else NA_character_
}

#' The native session manifest's `route` for each workflow type (KIT-020): the DIA route is `iq`; the
#' protein-level types have no quantification route and are `protein`.
#' @noRd
workflowRoutesByType <- c(DIA = "iq", LFQ = "protein", TMT = "protein")

#' The manifest route of a workflow type.
#' @param type One of `names(workflowRoutesByType)`.
#' @return The route string; aborts with `saved_entity_scope` for any other type.
#' @noRd
workflowRouteForType <- function(type) {
    if (!is.character(type) || length(type) != 1L || !type %in% names(workflowRoutesByType))
        cli::cli_abort("workflow type {.val {type}} has no native session route", class = "saved_entity_scope")
    workflowRoutesByType[[type]]
}

#' The store label of the owner an import format attaches (KIT-020): `DIA-IQ-native` for DIA-NN, otherwise
#' `<FORMAT>-native`.
#' @param format Detected import format.
#' @return The store's project id.
#' @noRd
workflowNativeStoreLabel <- function(format) {
    if (identical(format, "diann")) "DIA-IQ-native" else paste0(toupper(format), "-native")
}

#' Select the on-disk owner before import or configuration creates heavy values.
#' @param workflow_data Existing GUI bus.
#' @param paths Existing experiment paths.
#' @param session Existing Shiny session.
#' @param format Detected import format.
#' @param method Setup IQ/Limpa choice.
#' @param data_type The reader's `data_type` when known before the read (the Spectronaut quantity choice).
#' @return Current or newly attached state manager, invisibly.
#' @noRd
initializeDIANativeImport <- function(workflow_data, paths, session, format, method = "iq", data_type = NULL) {
    manager <- workflow_data$state_manager
    workflow_type <- workflowTypeForImportFormat(format, data_type)
    if (inherits(manager, "DIAStateAdapter")) {
        # KIT-020: an owner serves one route, fixed by its store's identity (the manager's own type is the fresh
        # default until the import sets it): the type must match, and an owner the import attached (a `-native`
        # label) accepts only its own format (audit F-09); for DIA-NN only the IQ method.
        identity <- manager$.context$store$identity
        if (is.na(workflow_type) || !identical(identity$workflow_id, workflow_type) ||
            (endsWith(identity$project_id, "-native") && !identical(identity$project_id, workflowNativeStoreLabel(format))) ||
            (identical(format, "diann") && !identical(method, "iq"))) {
            cli::cli_abort("reset the project before changing the route of an on-disk session", class = "saved_entity_scope")
        }
        if (identical(format, "diann")) workflow_data$dia_quantification_method <- method
        return(invisible(manager))
    }
    if (!isDIANativeDefaultEnabled() || !identical(Sys.info()[["sysname"]], "Linux") || is.na(workflow_type))
        return(invisible(manager))
    if (identical(format, "diann")) {
        workflow_data$dia_quantification_method <- method
        if (identical(method, "limpa")) return(invisible(manager))
        if (!identical(method, "iq")) cli::cli_abort("reset the project before activating native DIA-IQ", class = "saved_entity_scope")
    }
    if (!identical(manager$current_state, "initial") || length(manager$getHistory()) != 1L) {
        cli::cli_abort("reset the project before activating an on-disk session", class = "saved_entity_scope")
    }
    if (!requireNamespace("qs2", quietly = TRUE) || !requireNamespace("duckdb", quietly = TRUE) ||
        !requireNamespace("arrow", quietly = TRUE)) {
        cli::cli_abort("native DIA requires qs2, DuckDB and Arrow before import", class = "saved_entity_unsupported")
    }
    root <- normalizePath(paths$source_dir, mustWork = TRUE)
    runDIANativeTar("--version", tempfile("native-preflight-", tmpdir = root))
    for (command in c("dd", "sha256sum")) {
        version <- runDIANativeIO(list(c(command, "--version")),
            tempfile("native-io-preflight-", tmpdir = root))
        if (!length(version) || !grepl(paste0("^", command,
            " \\((GNU )?coreutils\\) "), version[[1L]])) {
            cli::cli_abort("native export requires GNU coreutils", class = "saved_entity_unsupported")
        }
    }
    state_root <- tempfile(".dia-native-", tmpdir = root)
    if (!dir.create(state_root, mode = "0700")) {
        cli::cli_abort("cannot create native DIA storage", class = "saved_entity_write")
    }
    store <- newWorkflowEntityStore(state_root, workflowNativeStoreLabel(format), workflow_type, max_cache_bytes = 0)
    manager <- newDIAStateAdapter(store, manager)
    manager$.context$sql_rollup <- TRUE
    manager$.context$session_format <- "native_v1"
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    workflow_data$state_manager <- attachDIAStateSession(manager, session, workflow_data)
    invisible(manager)
}

#' Stage an ordinary DIA history into a durable owner without mutating its source.
#' @param native Existing ordinary WorkflowState.
#' @param experiment_paths Existing project paths containing source_dir.
#' @return Reopened reference-backed owner with the same accepted topology.
#' @noRd
stageDIAGUIReferenceOwner <- function(native, experiment_paths) {
    if (!inherits(native, "WorkflowState") || !is.list(experiment_paths) ||
        !is.character(experiment_paths$source_dir) ||
        length(experiment_paths$source_dir) != 1L ||
        !dir.exists(experiment_paths$source_dir)) {
        cli::cli_abort("ordinary DIA owner cannot be staged", class = "saved_entity_scope")
    }
    before <- list(states = names(native$states), history = native$state_history,
        current = native$current_state, audits = names(native$audit_records))
    shadow <- native$clone(deep = FALSE)
    root <- tempfile(".dia-gui-adoption-", tmpdir = experiment_paths$source_dir)
    if (!dir.create(root, mode = "0700")) {
        cli::cli_abort("cannot create GUI adoption storage", class = "saved_entity_write")
    }
    store <- newWorkflowEntityStore(root, "DIA-GUI-pilot", "DIA", max_cache_bytes = 0)
    manager <- newDIAStateAdapter(store, shadow)
    accepted <- FALSE
    on.exit(if (!accepted) closeDIAStateAdapter(manager), add = TRUE)
    convertDIANativeLegacyCheckpoints(manager)
    manager$.context$session_format <- "native_v1"
    manager$.context$experiment_paths <- experiment_paths
    identity <- store$identity
    closeDIAStateAdapter(manager)
    reopened_store <- newWorkflowEntityStore(root, identity$project_id,
        identity$workflow_id, max_cache_bytes = 0)
    manager <- newDIAStateAdapter(reopened_store, reopen = TRUE)
    manager$.context$reference_states <- TRUE
    manager$.context$session_format <- "native_v1"
    manager$.context$experiment_paths <- experiment_paths
    after <- list(states = names(native$states), history = native$state_history,
        current = native$current_state, audits = names(native$audit_records))
    if (!identical(before, after) ||
        !identical(as.character(unlist(manager$getHistory())),
            as.character(unlist(native$getHistory()))) ||
        !identical(manager$current_state, native$current_state)) {
        closeDIAStateAdapter(manager)
        cli::cli_abort("GUI adoption changed accepted history", class = "saved_entity_corrupt")
    }
    accepted <- TRUE
    manager
}

#' Prepare immutable GUI roots before adoption or a resident checkpoint export.
#' @param manager,workflow_data Staged owner and current GUI values.
#' @return Named saved GUI roots.
#' @noRd
prepareDIAGUIReferences <- function(manager, workflow_data) {
    # Stage the known import/design assay holders before changing the live bus.
    # Their existing consumers already resolve saved values through the DIA codec.
    holders <- c("data_tbl", "data_cln", "peptide_data")
    references <- stats::setNames(lapply(holders, function(field) {
        value <- workflow_data[[field]]
        if (is.null(value)) return(NULL)
        if (!is.data.frame(value) && !isS4(value) && !isSavedWorkflowReference(value)) {
            cli::cli_abort("unsupported resident GUI assay holder", class = "saved_entity_unsupported")
        }
        saveDIAOwnedValue(value, field, manager)
    }), holders)
    if (!is.null(references$data_tbl)) {
        mapping <- workflow_data$column_mapping
        if (!is.list(mapping) || is.null(mapping$run_col)) {
            cli::cli_abort("resident import column mapping is required for adoption",
                class = "saved_entity_unsupported")
        }
        recordDIAImportMetadata(manager, workflow_data$data_tbl,
            list(data = workflow_data$data_tbl, data_type = workflow_data$data_type %||% "peptide",
                column_mapping = mapping), references$data_tbl)
        # Settle the Design observer's durable roots before dispatch freezes the
        # head. Its later reactive refresh must not invalidate the running job.
        buildDIAReferenceDesignState(manager, references$data_tbl,
            workflow_data$config_list, mapping)
        saveDIAOwnedValue(references$data_cln, "design-cleaned-data", manager)
        snapshot <- manager$.context$snapshot
        snapshot$design_current <- references$data_cln
        publishDIAStateHead(manager$.context, snapshot)
        # Notify the existing design observer that the source representation
        # changed without changing the user's design or removing any samples.
        manager$.context$design_adoption_pending <- references$data_tbl$entity_id
    }
    references
}

#' Swap a fully reopened pilot owner into one existing Shiny session.
#' @param manager Staged reference-backed owner.
#' @param workflow_data Existing reactive workflow bus.
#' @param session Current Shiny session.
#' @return Attached manager.
#' @noRd
adoptDIAGUIReferenceOwner <- function(manager, workflow_data, session) {
    if (!inherits(manager, "DIAStateAdapter") ||
        !isTRUE(manager$.context$reference_states) ||
        isTRUE(manager$.context$store$closed)) {
        cli::cli_abort("staged GUI owner is unavailable", class = "saved_entity_scope")
    }
    previous <- workflow_data$state_manager
    adopted <- FALSE
    previous_session_owner <- session$userData$dia_capacity_manager
    on.exit(if (!adopted && !identical(previous, manager)) {
        closeDIAStateAdapter(manager)
        session$userData$dia_capacity_manager <- previous_session_owner
    }, add = TRUE)
    if (!identical(as.character(unlist(previous$getHistory())),
        as.character(unlist(manager$getHistory()))) ||
        !identical(previous$current_state, manager$current_state)) {
        cli::cli_abort("staged GUI history was superseded", class = "workflow_execution_stale")
    }
    references <- prepareDIAGUIReferences(manager, workflow_data)
    holders <- names(references)
    workflowExecutionBoundary("before_gui_owner_swap")
    if (!identical(workflow_data$state_manager, previous)) {
        cli::cli_abort("GUI owner was replaced during adoption", class = "workflow_execution_stale")
    }
    manager$.context$gui_pilot_adopted <- TRUE
    workflow_data$state_manager <- attachDIAStateSession(manager, session, workflow_data,
        scope = "perf-107-sample-filter")
    for (field in holders) workflow_data[[field]] <- references[[field]]
    adopted <- TRUE
    if (inherits(previous, "DIAStateAdapter") && !identical(previous, manager)) {
        closeDIAStateAdapter(previous)
    }
    manager
}

#' Attach an externally validated native session only after its DA view is valid.
#' @param imported Result from openDIANativeSession.
#' @param workflow_data Existing GUI bus, left untouched on validation failure.
#' @param session Existing Shiny session.
#' @param source_dir Existing metadata output directory.
#' @return Validated current-state DA view.
#' @noRd
adoptDIANativeSession <- function(imported, workflow_data, session, source_dir, experiment_paths = NULL) {
    manager <- imported$manager
    adopted <- FALSE
    on.exit(if (!adopted) {
        extracted <- dirname(manager$.context$store$directory)
        closeDIAStateAdapter(manager)
        removeDIAReplacedSessionFolder(extracted, NULL, list(source_dir = source_dir))
    }, add = TRUE)
    convertDIANativeLegacyCheckpoints(manager)
    for (field in names(imported$extra_references)) {
        reference <- imported$extra_references[[field]]
        if (!identical(reference$payload$codec, "qs2")) next
        value <- resolveDIAStateValue(manager$.context, reference)$data
        if (is.data.frame(value) || isS4(value)) {
            imported$extra_references[[field]] <- saveDIAOwnedValue(value, field, manager)
        }
        value <- NULL
        manager$.context$working <- new.env(parent = emptyenv())
        gc(full = TRUE)
    }
    manager$.context$native_extra_references <- imported$extra_references
    view <- resolveDIANativeSessionView(manager, manager$.context$native_metadata_reference,
        imported$extra_references)
    validateProtDaFilteredSession(view, source_dir, manager)
    if (!is.null(experiment_paths)) restoreDIANativeGuiFiles(imported, experiment_paths)
    manager$.context$experiment_paths <- experiment_paths
    previous <- workflow_data$state_manager
    replaced_root <- if (inherits(previous, "DIAStateAdapter")) dirname(previous$.context$store$directory) else NULL
    if (inherits(previous, "DIAStateAdapter")) closeDIAStateAdapter(previous)
    previous <- NULL
    workflow_data$state_manager <- attachDIAStateSession(manager, session, workflow_data)
    for (field in c("data_tbl", "data_cln")) workflow_data[[field]] <- imported$extra_references[[field]]
    for (field in c("aa_seq_tbl_final", "uniprot_dat_cln", "data_format", "data_type", "column_mapping")) {
        workflow_data[[field]] <- imported$metadata[[field]]
    }
    workflow_data$dia_quantification_method <- "iq"
    workflow_data$native_session_restored <- Sys.time()
    adopted <- TRUE
    # The replaced private folder goes once the swap is complete (KIT-013 item 5).
    if (!is.null(replaced_root) && !is.null(experiment_paths)) {
        removeDIAReplacedSessionFolder(replaced_root, dirname(manager$.context$store$directory), experiment_paths)
    }
    view
}

#' Ingest a supported old proteomics envelope (`filtered_session_data`) into a fresh native owner after validation.
#' KIT-042: every workflow type the owner serves (`DIA`, `LFQ`, `TMT`) is adopted; a session from the retired
#' in-memory limpa route is refused by name.
#' @param value Decoded legacy session; its one-off full read remains explicit.
#' @param workflow_data,session Existing GUI owners.
#' @param source_dir Existing project storage directory.
#' @param experiment_paths Optional destination paths for current table/metadata files.
#' @return TRUE when adopted, FALSE for an unchanged non-native workflow (another omics type, or not Linux).
#' @noRd
ingestDIANativeLegacySession <- function(value, workflow_data, session, source_dir,
    experiment_paths = NULL) {
    enabled <- isDIANativeDefaultEnabled() || isDIANativeSession(workflow_data$state_manager)
    workflow_type <- value$workflow_type
    if (!enabled || !identical(Sys.info()[["sysname"]], "Linux") || !is.character(workflow_type) ||
        length(workflow_type) != 1L || !workflow_type %in% names(workflowRoutesByType)) return(FALSE)
    if (isTRUE(value$use_limpa) || isTRUE(value$limpa_applied) || !is.null(value$limpa_dpc_quant_results)) {
        cli::cli_abort("sessions from the retired in-memory limpa route cannot be adopted; rerun on the DIA with limpa workflow",
            class = "workflow_session_unsupported")
    }
    records <- value$r6_complete_states
    history <- value$r6_state_history
    current <- value$r6_current_state_name
    if (!is.list(records) || is.null(names(records)) || anyDuplicated(names(records)) ||
        !is.character(current) || length(current) != 1L || !current %in% names(records) ||
        !all(unlist(history) %in% names(records))) {
        cli::cli_abort("legacy DIA session history is incomplete", class = "saved_entity_corrupt")
    }
    for (record in records) {
        if (!is.list(record) || isSavedWorkflowReference(record$data) ||
            !is.null(record$data) && !methods::is(record$data, "PeptideQuantitativeData") &&
                !methods::is(record$data, "ProteinQuantitativeData")) {
            cli::cli_abort("legacy DIA history is not self-contained", class = "saved_entity_corrupt")
        }
        if (isS4(record$data)) methods::validObject(record$data)
    }
    checkSavedWorkflowPayload(value)
    validateProtDaFilteredSession(value, source_dir)
    root <- tempfile(".dia-native-legacy-", tmpdir = source_dir)
    if (!dir.create(root, mode = "0700")) cli::cli_abort("cannot stage legacy import", class = "saved_entity_write")
    # The store label follows the import's (KIT-020): the DIA-NN IQ owner, or the session's recorded format;
    # a session without a recorded format is labelled by its type.
    label <- if (identical(workflow_type, "DIA")) "DIA-IQ-native" else
        workflowNativeStoreLabel(if (is.character(value$data_format) && length(value$data_format) == 1L &&
            nzchar(value$data_format)) value$data_format else tolower(workflow_type))
    store <- newWorkflowEntityStore(root, label, workflow_type, max_cache_bytes = 0)
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    adopted <- FALSE
    on.exit(if (!adopted) closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$sql_rollup <- TRUE
    manager$.context$session_format <- "native_v1"
    manager$.context$experiment_paths <- experiment_paths
    restoreDIAStateSession(manager, records, history, current)
    convertDIANativeLegacyCheckpoints(manager)
    snapshot <- manager$.context$snapshot
    snapshot$workflow_type <- workflow_type
    snapshot$state_metadata[[current]] <- describeDIAProteinState(value$current_s4_object)
    # KIT-042: every adopted state's provenance is `legacy`, written at ingest by the adopting build (the reopen
    # backfill of KIT-025 then finds nothing to add).
    for (name in setdiff(unlist(history), names(snapshot$step_provenance))) {
        legacy <- diaStepProvenanceEntry(snapshot, list(operation = "legacy", variant = NA_character_,
            parameter_names = character()))
        snapshot <- legacy$snapshot
        snapshot$step_provenance[[name]] <- legacy$entry
    }
    publishDIAStateHead(manager$.context, snapshot)
    if (!is.null(experiment_paths)) {
        restoreDIANativeLegacyArtifacts(value, manager, experiment_paths, root)
    }
    metadata <- value
    metadata[c("r6_complete_states", "current_s4_object", "correlation_filtered_s4", "state_manager")] <- NULL
    if (!is.null(metadata$ruv_optimization_result$best_cancor_plot) &&
        !inherits(metadata$ruv_optimization_result$best_cancor_plot, "grob")) {
        metadata$ruv_optimization_result$best_cancor_plot <- ggplot2::ggplotGrob(
            stripProtNormQcPlotOwners(metadata$ruv_optimization_result$best_cancor_plot))
    }
    manager$.context$native_metadata_reference <- saveDIAStateValue(manager$.context,
        metadata, "legacy-import-metadata", codec = "rds")
    manager$.context$native_extra_references <- list()
    previous <- workflow_data$state_manager
    if (inherits(previous, "DIAStateAdapter")) {
        # KIT-042: the bus payloads the replaced owner held (the import's table references) would dangle once
        # that store closes; the legacy envelope carries no import table, so they are released.
        for (field in c("data_tbl", "data_cln")) {
            if (isSavedWorkflowReference(workflow_data[[field]])) workflow_data[[field]] <- NULL
        }
        closeDIAStateAdapter(previous)
    }
    previous <- NULL
    workflow_data$state_manager <- attachDIAStateSession(manager, session, workflow_data)
    if (identical(workflow_type, "DIA")) workflow_data$dia_quantification_method <- "iq"
    workflow_data$native_session_restored <- Sys.time()
    adopted <- TRUE
    TRUE
}

#' Recreate current legacy table/metadata artifacts without rerunning science.
#' @param value Validated old envelope.
#' @param manager Fresh native destination owner.
#' @param paths Existing destination experiment paths.
#' @param root Owned import staging directory.
#' @return TRUE invisibly; old conflicting project files remain untouched.
#' @noRd
restoreDIANativeLegacyArtifacts <- function(value, manager, paths, root) {
    source <- file.path(root, "legacy-source")
    protein <- file.path(root, "legacy-protein-qc")
    dir.create(source)
    dir.create(protein)
    ruv <- value$ruv_optimization_result %||% list(ruv_skipped = !isTRUE(value$ruv_applied))
    saveProtNormCorrelationResults(finalS4ForDe = value$current_s4_object,
        experimentPaths = list(protein_qc_dir = protein), normRuvOptimizationResult = ruv,
        compatibilityExport = FALSE)
    saveProtNormExportMetadataFiles(value, source, strict = TRUE)
    dependencies <- list()
    for (area in c("source", "protein_qc")) {
        directory <- if (area == "source") source else protein
        for (file in list.files(directory, full.names = TRUE)) {
            dependencies[[length(dependencies) + 1L]] <- list(path = file,
                destination = paste(area, basename(file), sep = "/"))
        }
    }
    restoreDIANativeGuiFiles(list(manager = manager, dependencies = dependencies), paths)
}

#' Identify the explicitly selected native-session route without hydrating state.
#' @param manager Existing state owner.
#' @return TRUE only for the native DIA route.
#' @noRd
isDIANativeSession <- function(manager) {
    inherits(manager, "DIAStateAdapter") && identical(manager$.context$session_format, "native_v1")
}

#' Resolve declared dependency copies after a project has moved.
#' @param manager Live native manager.
#' @return Export specifications preserving original source provenance.
#' @noRd
resolveDIANativeDependencies <- function(manager) {
    dependencies <- manager$.context$snapshot$file_dependencies %||% list()
    # KIT-026: recorded inputs (the report, the FASTA, the design) carry no source; they are not archive members.
    dependencies <- Filter(function(dependency) !is.null(dependency$source), dependencies)
    for (name in names(dependencies)) {
        restored <- manager$.context$native_dependencies[[name]]
        if (!is.null(restored) && identical(dependencies[[name]]$source, restored$source_provenance)) {
            if (!identical(workflowTableFileIdentity(restored$path), restored$file_identity)) {
                cli::cli_abort("restored dependency changed", class = "saved_entity_corrupt")
            }
            dependencies[[name]]$source_provenance <- dependencies[[name]]$source
            dependencies[[name]]$source <- restored$path
        }
    }
    dependencies
}

#' Map only existing project output areas used by the DIA GUI consumers.
#' @param paths Existing experiment path list.
#' @return Named native destination roots; no directory is created.
#' @noRd
diaNativeGuiRoots <- function(paths) {
    roots <- list(source = paths$source_dir, protein_qc = paths$protein_qc_dir,
        peptide_qc = paths$peptide_qc_dir, da = paths$da_output_dir)
    if (!is.null(paths$protein_qc_dir)) roots$publication_graphs <- paths$publication_graphs_dir %||%
        file.path(dirname(paths$protein_qc_dir), "publication_graphs")
    roots[!vapply(roots, is.null, logical(1))]
}

#' Include declared scientific/report artifacts, not a blanket project backup.
#' @param manager Native owner with experiment paths.
#' @param metadata Frozen metadata for the export.
#' @return File specifications including execution dependencies and existing QC assets.
#' @noRd
collectDIANativeGuiFiles <- function(manager, metadata) {
    paths <- manager$.context$experiment_paths
    if (is.null(paths$source_dir) || is.null(paths$protein_qc_dir)) {
        cli::cli_abort("native export lacks its project output paths", class = "saved_entity_scope")
    }
    roots <- diaNativeGuiRoots(paths)
    specs <- resolveDIANativeDependencies(manager)
    add <- function(area, relative) {
        source <- file.path(roots[[area]], relative)
        destination <- paste(area, relative, sep = "/")
        key <- paste0("gui_", sprintf("%05d", length(specs) + 1L))
        specs[[key]] <<- list(source = source, destination = destination)
    }
    final_stub <- paste0(if (!identical(metadata$ruv_mode, "skip")) "ruv_" else "",
        "normalised_results_cln_with_replicates")
    add("protein_qc", paste0(final_stub, ".tsv"))
    protein_names <- c("normalized_protein_matrix_pre_ruv.tsv", "remove_proteins_with_only_one_rep.tsv",
        paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")),
        paste0(rep(c("pre_norm", "post_norm", "ruv_corrected"), each = 4L), "_",
            rep(c("pca", "density", "rle", "correlation"), 3L), ".png"), "ruv_corrected_cancor.png")
    for (name in protein_names[file.exists(file.path(roots$protein_qc, protein_names))]) add("protein_qc", name)
    source_files <- c(accession_cleanup_results = "accession_cleanup_results.RDS",
        qc_params = "qc_params.RDS", fasta_metadata = "fasta_metadata.RDS",
        ruv_optimization_result = "ruv_optimization_results.RDS",
        mixed_species_analysis = "mixed_species_analysis.RDS",
        aa_seq_tbl_final = "aa_seq_tbl_final.RDS", uniprot_dat_cln = "uniprot_dat_cln.RDS")
    for (field in names(source_files)) {
        value <- metadata[[field]]
        if (!is.null(value) && (!is.data.frame(value) || nrow(value) > 0L)) add("source", source_files[[field]])
    }
    add("source", "filtered_session_summary.txt")
    if (file.exists(file.path(roots$source, "config.ini"))) add("source", "config.ini")
    if (file.exists(file.path(roots$source, "study_parameters.txt"))) add("source", "study_parameters.txt")
    if (file.exists(file.path(roots$source, "analysis_record.json"))) add("source", "analysis_record.json")
    if (!is.null(roots$peptide_qc) && file.exists(file.path(roots$peptide_qc, "peptide_values_imputed.tsv"))) {
        add("peptide_qc", "peptide_values_imputed.tsv")
    }
    graph_root <- file.path(roots$publication_graphs, "filtering_qc")
    if (dir.exists(graph_root)) {
        files <- list.files(graph_root, pattern = "[.](png|pdf|rds)$", recursive = TRUE, full.names = FALSE)
        for (name in files) add("publication_graphs", paste0("filtering_qc/", name))
    }
    if (!is.null(manager$.context$snapshot$da_results)) {
        if (!is.null(roots$da) && dir.exists(roots$da)) {
            files <- list.files(roots$da, recursive = TRUE, full.names = FALSE)
            if (length(files) > 10000L) cli::cli_abort("DA output inventory exceeds limits")
            for (name in files) add("da", name)
        }
        for (area in c("Volcano_Plots", "NumSigDaMolecules", "Heatmap")) {
            directory <- file.path(roots$publication_graphs, area)
            if (dir.exists(directory)) for (name in list.files(directory, recursive = TRUE))
                add("publication_graphs", paste(area, name, sep = "/"))
        }
    }
    specs
}

#' Refresh the automatic native archive for the completed Summary generation.
#' @param workflow Current workflow with accepted state and result references.
#' @param source_dir Existing native publication directory.
#' @return Completed native publication; no assay or S4 is reconstructed.
#' @noRd
refreshDIANativeSummaryArchive <- function(workflow, source_dir) {
    manager <- workflow$state_manager
    context <- manager$.context
    assertDIAExportReady(manager)
    owned <- context$owned_export
    if (is.null(owned$metadata_reference)) cli::cli_abort("export the native session before Summary")
    if (identical(owned$generation, context$head$entity_id)) return(invisible(owned))
    metadata <- resolveDIAStateValue(context, owned$metadata_reference)
    metadata$r6_current_state_name <- manager$current_state
    metadata$r6_state_history <- manager$state_history
    metadata$da_ui_params <- workflow$da_ui_params
    metadata$export_timestamp <- Sys.time()
    runDIANativeSessionExport(manager, metadata, source_dir,
        extra_references = owned$extra_references,
        dependencies = collectDIANativeGuiFiles(manager, metadata),
        prepare_outputs = function(metadata, request) {
            saveProtNormExportMetadataFiles(metadata, source_dir, strict = TRUE)
            summary <- buildProtNormExportSummaryContent(metadata, request$archive_name,
                timeFn = function() metadata$export_timestamp)
            writeLines(summary, file.path(source_dir, "filtered_session_summary.txt"))
        })
}

#' Require the existing proof service guard for an explicit S4 compatibility export.
#' @return TRUE invisibly or an explicit unsupported-export condition.
#' @noRd
assertDIACompatibilityGuard <- function() {
    if (!identical(Sys.info()[["sysname"]], "Linux") || !file.exists("/proc/self/cgroup")) {
        cli::cli_abort("S4 compatibility export requires the existing guarded runner", class = "saved_entity_scope")
    }
    cgroup <- readLines("/proc/self/cgroup", warn = FALSE)
    cgroup <- sub("^0::", "", cgroup[startsWith(cgroup, "0::")])
    seconds <- suppressWarnings(as.numeric(Sys.getenv("MULTISCHOLAR_DIA_DIAGNOSTIC_SECONDS")))
    if (length(cgroup) != 1L || !grepl("multischolar-dia[0-9]+-", basename(cgroup)) ||
        length(seconds) != 1L || !is.finite(seconds) || seconds < 1 || seconds > 3600) {
        cli::cli_abort("S4 compatibility export requires the existing guarded runner", class = "saved_entity_scope")
    }
    root <- file.path("/sys/fs/cgroup", sub("^/", "", cgroup))
    limit <- suppressWarnings(as.numeric(readLines(file.path(root, "memory.max"), warn = FALSE)))
    swap <- readLines(file.path(root, "memory.swap.max"), warn = FALSE)
    if (length(limit) != 1L || !is.finite(limit) || limit > 16 * 1024^3 || !identical(swap, "0")) {
        cli::cli_abort("S4 compatibility export requires the unchanged memory and swap guards",
            class = "saved_entity_scope")
    }
    invisible(TRUE)
}

#' Export a real S4 RDS only in response to an explicit compatibility request.
#' @param manager Current owner; completed native outputs are never modified here.
#' @param path Requested RDS destination.
#' @param state_name Optional selected state.
#' @return Destination path after atomic publication; failures leave the native session valid.
#' @noRd
exportDIAS4Compatibility <- function(manager, path, state_name = NULL) {
    assertDIACompatibilityGuard()
    temporary <- tempfile(".s4-compatibility-", tmpdir = dirname(path))
    on.exit(unlink(temporary), add = TRUE)
    withDIAStateScope(manager, function() {
        object <- manager$getState(state_name)
        if (!isS4(object)) cli::cli_abort("compatibility export requires a scientific state")
        saveRDS(object, temporary, version = 3)
        invisible(NULL)
    })
    syncWorkflowExportPaths(temporary, dirname(path))
    if (!file.rename(temporary, path)) cli::cli_abort("cannot publish S4 compatibility export")
    syncWorkflowExportPaths(c(path, dirname(path)), dirname(path))
    invisible(path)
}

#' Restore archived GUI files only into vacant or byte-identical destinations.
#' @param imported Fully validated native session.
#' @param paths Existing destination experiment paths.
#' @return The files this call created (character, possibly empty), invisibly, so a failed or
#'   stale Open can remove them again (KIT-013); conflicting projects are not overwritten.
#' @noRd
restoreDIANativeGuiFiles <- function(imported, paths) {
    roots <- diaNativeGuiRoots(paths)
    copies <- list()
    for (dependency in imported$dependencies) {
        parts <- strsplit(dependency$destination, "/", fixed = TRUE)[[1L]]
        if (!parts[[1L]] %in% names(roots)) next
        root <- normalizePath(roots[[parts[[1L]]]], mustWork = TRUE)
        relative <- paste(parts[-1L], collapse = "/")
        checkSavedEntityPath(imported$manager$.context$store, relative)
        target <- file.path(root, relative)
        parent <- dirname(target)
        while (!dir.exists(parent)) parent <- dirname(parent)
        canonical <- normalizePath(parent, mustWork = TRUE)
        if (!identical(canonical, root) && !startsWith(canonical, paste0(root, "/"))) {
            cli::cli_abort("native artifact destination escapes project", class = "saved_entity_path")
        }
        source <- describeDIANativeFile(dependency$path)
        if (file.exists(target)) {
            if (!identical(describeDIANativeFile(target), source)) {
                cli::cli_abort("native restore conflicts with project files; choose a fresh project",
                    class = "saved_entity_conflict")
            }
        } else copies[[target]] <- list(source = dependency$path, description = source, root = root)
    }
    created <- character()
    completed <- FALSE
    on.exit(if (!completed) discardDIANativeGuiFiles(created), add = TRUE)
    for (target in names(copies)) {
        spec <- copies[[target]]
        if (!dir.exists(dirname(target))) dir.create(dirname(target), recursive = TRUE)
        created <- c(created, target)
        if (!file.copy(spec$source, target, overwrite = FALSE) ||
            !identical(describeDIANativeFile(target), spec$description)) {
            cli::cli_abort("native artifact restore failed", class = "saved_entity_write")
        }
        syncWorkflowExportPaths(c(target, dirname(target), spec$root), spec$root)
    }
    completed <- TRUE
    invisible(created)
}

#' Remove the GUI files a restore created, after the Open that restored them failed.
#' @param created Paths returned by `restoreDIANativeGuiFiles`.
#' @return NULL invisibly.
#' @noRd
discardDIANativeGuiFiles <- function(created) {
    for (target in created) if (file.exists(target)) unlink(target)
    invisible(NULL)
}

#' Remove the private session folder an Open just replaced (hand-off item 5, CORE-014).
#' @param root The replaced manager's project root (its store directory's parent).
#' @param keep The adopted manager's project root; never removed.
#' @param paths Experiment paths whose `source_dir` must contain the folder.
#' @return TRUE when the folder was removed, FALSE when the rule did not apply, invisibly. Only the
#'   application's own dot-prefixed folders directly under the project's source directory qualify
#'   (`.dia-native-*`, `.dia-native-import-*`, `.dia-gui-adoption-*`); a user's folder is never
#'   touched.
#' @noRd
removeDIAReplacedSessionFolder <- function(root, keep, paths) {
    if (!is.character(root) || length(root) != 1L || is.na(root) || !dir.exists(root) ||
        !is.list(paths) || !is.character(paths$source_dir) || !dir.exists(paths$source_dir)) {
        return(invisible(FALSE))
    }
    root <- normalizePath(root, mustWork = TRUE)
    source <- normalizePath(paths$source_dir, mustWork = TRUE)
    kept <- if (is.character(keep) && length(keep) == 1L && dir.exists(keep)) normalizePath(keep) else ""
    if (identical(root, kept) || !identical(dirname(root), source) ||
        !grepl("^[.]dia-(native|native-import|gui-adoption)-", basename(root))) {
        return(invisible(FALSE))
    }
    unlink(root, recursive = TRUE)
    invisible(!dir.exists(root))
}

#' Gather the existing export metadata without materializing historical states.
#' @param workflow_data,norm_data,input Existing GUI data owners and settings.
#' @return Metadata value and explicitly named extra descriptor roots.
#' @noRd
collectDIANativeSessionMetadata <- function(workflow_data, norm_data, input) {
    manager <- workflow_data$state_manager
    if (!isDIANativeSession(manager)) cli::cli_abort("native DIA session is not active", class = "saved_entity_scope")
    context <- manager$.context
    state <- context$snapshot$current_state
    summary <- context$snapshot$state_metadata[[state]]
    if (is.null(summary)) {
        cli::cli_abort("native export needs component metadata; explicitly convert the legacy session first",
            class = "saved_entity_unsupported")
    }
    config <- workflow_data$config_list
    global <- config$globalParameters
    # KIT-042: no refusal by `use_limpa` here; the limpa route is a workflow type of its own (KIT-046), and the
    # in-memory limpa route never attaches an owner.
    extra <- list()
    for (field in c("data_tbl", "data_cln")) {
        value <- workflow_data[[field]]
        if (!is.null(value)) extra[[field]] <- saveDIAOwnedValue(value, field, manager)
    }
    correlation <- norm_data$correlation_filtered_obj
    if (!is.null(correlation)) extra$correlation_filtered <-
        saveDIAOwnedValue(correlation, "correlation_filtered_obj", manager)
    for (field in c("normalization_input_reference", "normalized_protein_obj", "ruv_normalized_obj")) {
        if (isSavedWorkflowReference(norm_data[[field]])) extra[[field]] <- norm_data[[field]]
    }
    metadata <- list(r6_current_state_name = state, r6_state_history = context$snapshot$state_history,
        contrasts_tbl = workflow_data$contrasts_tbl, design_matrix = workflow_data$design_matrix,
        config_list = config, export_timestamp = Sys.time(), normalization_method = input$norm_method,
        ruv_mode = input$ruv_mode, ruv_applied = input$ruv_mode != "skip", ruv_k = norm_data$best_k,
        correlation_threshold = norm_data$correlation_threshold,
        output_policy = "native_first_v1",
        # KIT-042: the owner's type (the store identity, KIT-020) is authoritative; the config and the state's
        # own record are fallbacks for stores written before the type was set at import.
        workflow_type = context$store$identity$workflow_id %||% global$workflow_type %||% summary$workflow_type %||% "DIA",
        report_template = global$report_template %||% summary$report_template,
        use_limpa = FALSE, limpa_applied = FALSE, limpa_parameters = NULL, limpa_dpc_quant_results = NULL,
        final_protein_count = summary$proteins, final_sample_count = summary$samples)
    for (field in c("fasta_metadata", "accession_cleanup_results", "ruv_optimization_result",
        "qc_params", "protein_counts", "mixed_species_analysis", "aa_seq_tbl_final", "uniprot_dat_cln",
        "data_format", "data_type", "column_mapping")) metadata[field] <- list(workflow_data[[field]])
    list(metadata = metadata, references = extra)
}

#' Build the existing DA view from one native head, never its full history.
#' @param manager Live native adapter.
#' @param metadata_reference Saved session metadata.
#' @param extra_references Named GUI-owned references.
#' @return Legacy-shaped current-state view for the existing DA validator.
#' @noRd
resolveDIANativeSessionView <- function(manager, metadata_reference, extra_references) {
    context <- manager$.context
    metadata <- resolveDIAStateValue(context, metadata_reference)
    state <- context$snapshot$current_state
    if (usesDIAReferenceStates(manager)) {
        metadata$current_s4_object <- getDIAStateReference(manager, state)
        metadata$correlation_filtered_s4 <- if (hasDIAState(manager, "correlation_filtered"))
            getDIAStateReference(manager, "correlation_filtered") else metadata$current_s4_object
        metadata$state_metadata <- readDIAStateMetadata(manager, metadata$current_s4_object)$data$slots
    } else {
        metadata$current_s4_object <- materializeDIAStateRecord(context, context$snapshot$states[[state]])$data
        metadata$correlation_filtered_s4 <- if (is.null(extra_references$correlation_filtered))
            metadata$current_s4_object else resolveDIAOwnedValue(extra_references$correlation_filtered, manager)
    }
    metadata$r6_current_state_name <- state
    metadata$r6_state_history <- context$snapshot$state_history
    metadata
}

#' Resolve only an unchanged completed native publication owned by this session.
#' @param manager Current native adapter.
#' @param pointer_path Expected latest pointer in the active project.
#' @return Minimal DA view, or NULL for a genuine external restore.
#' @noRd
resolveDIANativeOwnedExport <- function(manager, pointer_path) {
    if (!isDIANativeSession(manager)) return(NULL)
    context <- manager$.context
    checkWorkflowEntityStore(context$store)
    assertDIAExportReady(manager)
    request <- context$owned_export
    if (!identical(request$format, "native_v1") ||
        !identical(request$identity, context$store$identity) ||
        !identical(request$generation, context$head$entity_id) ||
        !file.exists(pointer_path) || file.size(pointer_path) > 8192 ||
        !identical(normalizePath(pointer_path), request$latest_path)) return(NULL)
    pointer <- jsonlite::read_json(pointer_path)
    if (!identical(pointer$status, "complete") || !identical(pointer$request_id, request$request_id) ||
        !identical(pointer$generation, request$generation) ||
        !identical(pointer$archive, basename(request$archive_path)) ||
        !identical(pointer$sha256, request$archive_sha256) || !file.exists(request$completion) ||
        !identical(workflowTableFileIdentity(request$archive_path), request$archive_identity)) {
        cli::cli_abort("owned native publication changed", class = "saved_entity_corrupt")
    }
    resolveDIANativeSessionView(manager, request$metadata_reference, request$extra_references)
}

#' Resolve the single archive named by a bounded project latest pointer.
#' @param pointer_path Native latest pointer.
#' @return Verified archive path; external loading still validates its contents.
#' @noRd
readDIANativeLatestPath <- function(pointer_path) {
    readDIANativeLatestArtifact(pointer_path)$path
}

#' Resolve a native publication and retain its validated digest for this copy operation.
#' @param pointer_path Existing bounded latest pointer.
#' @return Path and already-verified byte description.
#' @noRd
readDIANativeLatestArtifact <- function(pointer_path) {
    if (!file.exists(pointer_path) || file.size(pointer_path) > 8192) {
        cli::cli_abort("native latest pointer is missing or oversized", class = "saved_entity_corrupt")
    }
    pointer <- jsonlite::read_json(pointer_path)
    if (!identical(pointer$schema, "multischolar.dia_native_latest.v1") ||
        !identical(pointer$status, "complete") || !is.character(pointer$archive) ||
        length(pointer$archive) != 1L ||
        !grepl("^filtered_session_[0-9]{8}_[0-9]{6}_native-export-[A-Za-z0-9_-]+[.]msdia[.]tar$", pointer$archive)) {
        cli::cli_abort("native latest pointer is invalid", class = "saved_entity_corrupt")
    }
    path <- file.path(dirname(pointer_path), pointer$archive)
    actual <- describeDIANativeFile(path)
    if (!isTRUE(actual$bytes == pointer$bytes) || !identical(actual$sha256, pointer$sha256)) {
        cli::cli_abort("latest native archive differs", class = "saved_entity_corrupt")
    }
    list(path = path, description = actual)
}

#' Copy a completed portable native session for an included file consumer.
#' @param source_dir Directory containing the completed native latest pointer.
#' @param destination Explicit archive filename, never an S4 RDS filename.
#' @return Published archive path after independent byte verification.
#' @noRd
copyDIANativeSessionArtifact <- function(source_dir, destination) {
    if (!endsWith(destination, ".msdia.tar")) {
        cli::cli_abort("native artifact requires an archive filename", class = "saved_entity_path")
    }
    artifact <- readDIANativeLatestArtifact(file.path(source_dir, "filtered_session_latest.json"))
    source <- artifact$path
    root <- normalizePath(dirname(destination), mustWork = TRUE)
    temporary <- tempfile(".native-copy-", tmpdir = root)
    log_path <- tempfile(".native-copy-log-", tmpdir = root)
    on.exit(unlink(c(temporary, log_path)), add = TRUE)
    expected <- artifact$description
    if (file.exists(destination) && identical(describeDIANativeFile(destination, log_path), expected)) {
        return(destination)
    }
    snapshotWorkflowFile(source, temporary, log_path)
    if (!identical(describeDIANativeFile(temporary, log_path), expected)) {
        cli::cli_abort("native artifact copy differs", class = "saved_entity_corrupt")
    }
    syncWorkflowExportPaths(temporary, root)
    if (!file.rename(temporary, destination)) {
        cli::cli_abort("native artifact publication failed", class = "saved_entity_write")
    }
    syncWorkflowExportPaths(c(destination, root), root)
    destination
}

#' Recognize the supported pre-protein DIA checkpoint boundary.
#' @param manager Existing workflow owner.
#' @return Whether its selected state is a supported peptide checkpoint.
#' @noRd
isDIAPeptideCheckpoint <- function(manager) {
    !is.null(manager) && identical(manager$workflow_type, "DIA") &&
        manager$current_state %in% c("raw_data_s4", "qvalue_filtered",
            "protein_peptide_filtered", "precursor_rollup", "intensity_filtered",
            "sample_filtered", "replicate_filtered", "imputed")
}

#' Validate the supported peptide checkpoint route before staging or replacement.
#' @param metadata Workflow/checkpoint metadata; no assay hydration.
#' @return NULL invisibly.
#' @noRd
validateDIAPeptideCheckpointRoute <- function(metadata) {
    method <- metadata$dia_quantification_method
    global <- metadata$config_list$globalParameters
    if ((!is.null(method) && !identical(method, "iq")) || isTRUE(global$use_limpa) ||
        isTRUE(metadata$use_limpa) || isTRUE(metadata$limpa_applied) ||
        !is.null(metadata$limpa_parameters) || !is.null(metadata$limpa_dpc_quant_results)) {
        cli::cli_abort("peptide checkpoint archives support the DIA IQ route only",
            class = "saved_entity_unsupported")
    }
    invisible(NULL)
}

#' Validate every metadata field consumed during peptide checkpoint restoration.
#' @param metadata Native archive metadata.
#' @return NULL invisibly; invalid metadata cannot replace a live owner.
#' @noRd
validateDIAPeptideCheckpointMetadata <- function(metadata) {
    checkSavedWorkflowPayload(metadata)
    validateDIAPeptideCheckpointRoute(metadata)
    design <- metadata$peptide_design
    valid_design <- function(value) is.data.frame(value) && nrow(value) > 0L &&
        all(c("Run", "group") %in% names(value)) && !anyNA(value$Run) &&
        !anyDuplicated(value$Run) && all(nzchar(as.character(value$Run)))
    if (!identical(metadata$checkpoint_kind, "peptide-qc") ||
        !identical(metadata$workflow_type, "DIA") || !is.list(metadata$config_list) ||
        !valid_design(metadata$design_matrix) || !is.list(design) ||
        !valid_design(design$design_matrix) ||
        !(is.character(design$groups) || is.factor(design$groups)) ||
        !(is.character(design$factors) || is.factor(design$factors)) ||
        !is.character(design$removed_samples) || !is.data.frame(design$contrasts) ||
        !is.character(design$formula) || length(design$formula) != 1L ||
        is.na(design$formula) || !nzchar(design$formula) ||
        !is.list(metadata$column_mapping) ||
        !is.character(metadata$column_mapping$run_col) ||
        length(metadata$column_mapping$run_col) != 1L ||
        is.na(metadata$column_mapping$run_col) || !nzchar(metadata$column_mapping$run_col)) {
        cli::cli_abort("peptide checkpoint restore metadata is incomplete or invalid",
            class = "saved_entity_invalid")
    }
    invisible(NULL)
}

#' Require a quiescent workflow before checkpoint save or replacement.
#' @param workflow_data,session Existing GUI owners.
#' @return NULL invisibly.
#' @noRd
checkDIAPeptideCheckpointIdle <- function(workflow_data, session, saving = FALSE) {
    app <- session$userData$workflow_execution
    if (inherits(app, "WorkflowExecutionApplication") && !is.null(app$active)) {
        cli::cli_abort("wait for the running action before saving or opening a checkpoint",
            class = "workflow_execution_busy")
    }
    manager <- workflow_data$state_manager
    if (inherits(manager, "DIAStateAdapter")) assertDIAExportReady(manager,
        allow_action = if (saving) "session-export" else NULL)
    invisible(NULL)
}

#' Save peptide history with the existing portable native archive and publisher.
#' @param workflow_data,session,paths Existing GUI state, session and project paths.
#' @return Existing publication receipt. No assay is returned.
#' @noRd
saveDIAPeptideCheckpoint <- function(workflow_data, session, paths) {
    validateDIAPeptideCheckpointRoute(workflow_data)
    checkDIAPeptideCheckpointIdle(workflow_data, session, saving = TRUE)
    manager <- workflow_data$state_manager
    if (!isDIAPeptideCheckpoint(manager)) {
        cli::cli_abort("select a DIA peptide state before saving a checkpoint",
            class = "saved_entity_unsupported")
    }
    design_reader <- session$userData$peptide_checkpoint_design
    if (!is.function(design_reader)) cli::cli_abort("Design checkpoint metadata is unavailable")
    if (!isDIANativeSession(manager)) {
        # Saving a resident workflow exports a staged shadow. It never replaces
        # the live owner, including when export fails.
        manager <- stageDIAGUIReferenceOwner(manager, paths)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        prepared <- prepareDIAGUIReferences(manager, workflow_data)
        snapshot <- new.env(parent = emptyenv())
        for (field in c("config_list", "contrasts_tbl", "design_matrix", "fasta_metadata",
            "accession_cleanup_results", "ruv_optimization_result", "qc_params", "protein_counts",
            "mixed_species_analysis", "aa_seq_tbl_final", "uniprot_dat_cln", "data_format",
            "data_type", "column_mapping")) snapshot[[field]] <- workflow_data[[field]]
        snapshot$state_manager <- manager
        snapshot$data_tbl <- prepared$data_tbl
        snapshot$data_cln <- prepared$data_cln
        workflow_data <- snapshot
    }
    # Early checkpoints have no normalized-protein summary. The existing archive
    # permits metadata records; keep the later-stage collector/DA contract intact.
    metadata <- list(workflow_type = "DIA", checkpoint_kind = "peptide-qc",
        r6_current_state_name = manager$current_state, r6_state_history = manager$getHistory(),
        export_timestamp = Sys.time(), peptide_design = design_reader(), ruv_mode = "skip",
        dia_quantification_method = "iq")
    for (field in c("config_list", "design_matrix", "contrasts_tbl", "qc_params",
        "aa_seq_tbl_final", "uniprot_dat_cln", "data_format", "data_type", "column_mapping",
        "fasta_metadata", "accession_cleanup_results")) metadata[field] <- list(workflow_data[[field]])
    references <- lapply(c("data_tbl", "data_cln"), function(field) {
        saveDIAOwnedValue(workflow_data[[field]], field, manager)
    })
    names(references) <- c("data_tbl", "data_cln")
    collected <- list(metadata = metadata, references = references)
    collected$metadata$config_list <- prepareDIAStateRecord(
        list(data = NULL, config = collected$metadata$config_list))$config
    validateDIAPeptideCheckpointMetadata(collected$metadata)
    output <- file.path(paths$source_dir, "peptide-checkpoints")
    if (!dir.exists(output) && !dir.create(output)) cli::cli_abort("cannot create checkpoint directory")
    dependencies <- resolveDIANativeDependencies(manager)
    roots <- diaNativeGuiRoots(paths)
    # Only existing peptide outputs belong to an early checkpoint. The later
    # normalization exporter still requires its own protein outputs.
    for (area in intersect(c("peptide_qc", "publication_graphs"), names(roots))) {
        source <- if (area == "publication_graphs") file.path(roots[[area]], "filtering_qc") else roots[[area]]
        if (!dir.exists(source)) next
        files <- list.files(source, "[.](tsv|csv|png|pdf|rds)$", recursive = TRUE, full.names = TRUE)
        if (length(files) > 10000L) cli::cli_abort("checkpoint output inventory exceeds limits")
        for (file in files) {
            relative <- substring(file, nchar(roots[[area]]) + 2L)
            dependencies[[paste0("peptide_gui_", length(dependencies) + 1L)]] <-
                list(source = file, destination = paste(area, relative, sep = "/"))
        }
    }
    runDIANativeSessionExport(manager, collected$metadata, output,
        extra_references = collected$references, dependencies = dependencies)$publication
}

#' Open a validated peptide archive without going through the protein DA loader.
#' @param path Explicit local native archive.
#' @param workflow_data,session,paths Current GUI ownership and project paths.
#' @return Restored state name after a successful owner replacement.
#' @noRd
openDIAPeptideCheckpoint <- function(path, workflow_data, session, paths) {
    checkDIAPeptideCheckpointIdle(workflow_data, session)
    previous <- workflow_data$state_manager
    imported <- openDIANativeSession(path, paths$source_dir)
    manager <- imported$manager
    adopted <- FALSE
    # A failed Open closes the staged owner and removes the private folder it extracted into.
    on.exit(if (!adopted) {
        extracted <- dirname(manager$.context$store$directory)
        closeDIAStateAdapter(manager)
        removeDIAReplacedSessionFolder(extracted, NULL, paths)
    }, add = TRUE)
    metadata <- imported$metadata
    validateDIAPeptideCheckpointMetadata(metadata)
    if (!identical(metadata$checkpoint_kind, "peptide-qc") || !isDIAPeptideCheckpoint(manager) ||
        !is.data.frame(metadata$peptide_design$design_matrix) ||
        !is.list(metadata$column_mapping) || is.null(metadata$column_mapping$run_col) ||
        !isSavedWorkflowReference(imported$extra_references$data_tbl) ||
        !isSavedWorkflowReference(imported$extra_references$data_cln)) {
        cli::cli_abort("archive is not a supported peptide checkpoint", class = "saved_entity_invalid")
    }
    convertDIANativeLegacyCheckpoints(manager)
    for (reference in imported$extra_references) validateSavedWorkflowReference(manager$.context$store, reference)
    if (!identical(workflow_data$state_manager, previous)) {
        cli::cli_abort("project changed during checkpoint open", class = "workflow_execution_stale")
    }
    manager$.context$experiment_paths <- paths
    manager$.context$gui_pilot_adopted <- inherits(previous, "DIAStateAdapter") &&
        isTRUE(previous$.context$gui_pilot_adopted)
    manager$.context$native_extra_references <- imported$extra_references
    manager$.context$peptide_design_restore <- c(metadata$peptide_design,
        list(source_reference = imported$extra_references$data_tbl$entity_id))
    fields <- c("state_manager", "data_tbl", "data_cln", "config_list", "design_matrix",
        "contrasts_tbl", "qc_params", "aa_seq_tbl_final", "uniprot_dat_cln", "data_format",
        "data_type", "column_mapping", "fasta_metadata", "accession_cleanup_results",
        "dia_quantification_method", "native_session_restored", "state_update_trigger")
    previous_values <- stats::setNames(lapply(fields, function(field) workflow_data[[field]]), fields)
    previous_session_owner <- session$userData$dia_capacity_manager
    restored_files <- character()
    # Register rollback before the first live assignment and before any project file is
    # restored (KIT-013 item 2): a stale or failed Open puts the bus back and removes the
    # files this Open created. Assays remain shared; no deep clone is introduced here.
    on.exit(if (!adopted) {
        for (field in fields) workflow_data[[field]] <- previous_values[[field]]
        session$userData$dia_capacity_manager <- previous_session_owner
        discardDIANativeGuiFiles(restored_files)
    }, add = TRUE, after = FALSE)
    restored_files <- restoreDIANativeGuiFiles(imported, paths)
    if (!identical(workflow_data$state_manager, previous)) {
        cli::cli_abort("project changed during checkpoint open", class = "workflow_execution_stale")
    }
    workflow_data$state_manager <- attachDIAStateSession(manager, session, workflow_data)
    workflowExecutionBoundary("peptide_checkpoint_after_owner_swap")
    for (field in c("data_tbl", "data_cln")) workflow_data[[field]] <- imported$extra_references[[field]]
    # The design and annotation payloads are owned by the bus only (KIT-013 item 1).
    for (field in c("config_list", "design_matrix", "contrasts_tbl", "qc_params",
        "aa_seq_tbl_final", "uniprot_dat_cln", "data_format", "data_type", "column_mapping",
        "fasta_metadata", "accession_cleanup_results")) workflow_data[[field]] <- metadata[[field]]
    workflowExecutionBoundary("peptide_checkpoint_after_metadata")
    workflow_data$dia_quantification_method <- "iq"
    workflow_data$native_session_restored <- Sys.time()
    workflow_data$state_update_trigger <- Sys.time()
    adopted <- TRUE
    replaced_root <- if (inherits(previous, "DIAStateAdapter")) dirname(previous$.context$store$directory) else NULL
    if (inherits(previous, "DIAStateAdapter")) closeDIAStateAdapter(previous)
    # The replaced private folder is removed only after the previous owner is released and the
    # rollback handler above has exited (KIT-013 item 5); it runs after that handler.
    if (!is.null(replaced_root)) on.exit(if (adopted) removeDIAReplacedSessionFolder(replaced_root,
        dirname(manager$.context$store$directory), paths), add = TRUE, after = TRUE)
    manager$current_state
}
