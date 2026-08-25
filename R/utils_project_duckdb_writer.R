# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

PROJECT_REGISTRY_OWNER_SCHEMA <- "multischolar.project_registry_owner"
PROJECT_REGISTRY_OWNER_VERSION <- 1L
.PROJECT_REGISTRY_ACTIVE_WRITERS <- new.env(parent = emptyenv())

#' Build the process-local writer registry key
#'
#' @param registry A validated project registry.
#'
#' @return A normalized writer-lock path.
#' @noRd
projectRegistryWriterKey <- function(registry) {
    normalizePath(
        projectRegistryPath(registry, "lock"),
        winslash = "/",
        mustWork = FALSE
    )
}

projectRegistryEnsureDirectory <- function(registry, path_name) {
    path <- projectRegistryPath(registry, path_name)
    if (file.exists(path) && !dir.exists(path)) {
        projectRegistryAbort(
            sprintf("project registry directory '%s' is occupied", path_name),
            "multischolar_registry_path_conflict"
        )
    }
    if (!dir.exists(path)) {
        created <- dir.create(path, recursive = TRUE, mode = "0700", showWarnings = FALSE)
        if (!isTRUE(created) && !dir.exists(path)) {
            projectRegistryAbort(
                sprintf("project registry directory '%s' could not be created", path_name),
                "multischolar_registry_write_failed"
            )
        }
    }
    projectRegistryPath(registry, path_name, must_exist = TRUE)
}

projectRegistryEnsureParent <- function(registry, path_name) {
    registry <- validateProjectRegistry(registry)
    relative_parent <- artifactNormalizeRelativePath(dirname(
        registry$relative_paths[[path_name]]
    ))
    path <- artifactResolveContainedPath(registry$project_root, relative_parent)
    if (file.exists(path) && !dir.exists(path)) {
        projectRegistryAbort(
            sprintf("project registry parent for '%s' is occupied", path_name),
            "multischolar_registry_path_conflict"
        )
    }
    if (!dir.exists(path)) {
        created <- dir.create(path, recursive = TRUE, mode = "0700", showWarnings = FALSE)
        if (!isTRUE(created) && !dir.exists(path)) {
            projectRegistryAbort(
                sprintf("project registry parent for '%s' could not be created", path_name),
                "multischolar_registry_write_failed"
            )
        }
    }
    artifactResolveContainedPath(registry$project_root, relative_parent, must_exist = TRUE)
}

projectRegistryHost <- function() {
    host <- unname(Sys.info()[["nodename"]])
    if (!workflowCapabilityScalarString(host)) host <- "unknown-host"
    host
}

projectRegistryPackageVersion <- function() {
    tryCatch(
        as.character(utils::packageVersion("MultiScholaR")),
        error = \(error) "development"
    )
}

projectRegistryProcessStartedAt <- function() {
    if (!requireNamespace("ps", quietly = TRUE)) return(NA_character_)
    started <- tryCatch(
        ps::ps_create_time(ps::ps_handle()),
        error = \(error) NA_real_
    )
    if (is.na(started)) return(NA_character_)
    format(as.POSIXct(started, origin = "1970-01-01", tz = "UTC"),
        "%Y-%m-%dT%H:%M:%OS6Z",
        tz = "UTC"
    )
}

newProjectRegistryOwner <- function(registry) {
    registry <- validateProjectRegistry(registry)
    list(
        schema = PROJECT_REGISTRY_OWNER_SCHEMA,
        schema_version = PROJECT_REGISTRY_OWNER_VERSION,
        project_id = registry$project_id,
        host = projectRegistryHost(),
        pid = as.integer(Sys.getpid()),
        process_started_at = projectRegistryProcessStartedAt(),
        session_id = artifactOpaqueId("session"),
        owner_token = artifactOpaqueId("owner"),
        acquired_at = artifactRefUtcNow(),
        package_version = projectRegistryPackageVersion()
    )
}

validateProjectRegistryOwner <- function(owner) {
    required <- c(
        "schema", "schema_version", "project_id", "host", "pid",
        "process_started_at", "session_id", "owner_token", "acquired_at",
        "package_version"
    )
    if (is.list(owner) && "process_started_at" %in% names(owner) &&
        is.null(owner$process_started_at)) {
        owner$process_started_at <- NA_character_
    }
    valid <- is.list(owner) && identical(names(owner), required) &&
        identical(owner$schema, PROJECT_REGISTRY_OWNER_SCHEMA) &&
        identical(as.integer(owner$schema_version), PROJECT_REGISTRY_OWNER_VERSION) &&
        workflowCapabilityScalarString(owner$project_id) &&
        workflowCapabilityScalarString(owner$host) &&
        length(owner$pid) == 1L && is.numeric(owner$pid) && !is.na(owner$pid) &&
        owner$pid >= 1L && artifactRefValidUtc(owner$acquired_at) &&
        workflowCapabilityScalarString(owner$package_version)
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry writer owner metadata is malformed",
            "multischolar_malformed_registry_owner"
        )
    }
    valid_process_start <- length(owner$process_started_at) == 1L &&
        (is.na(owner$process_started_at) ||
            artifactRefValidUtc(owner$process_started_at))
    if (!isTRUE(valid_process_start)) {
        projectRegistryAbort(
            "project registry process start timestamp is malformed",
            "multischolar_malformed_registry_owner"
        )
    }
    if (!is.na(owner$process_started_at) &&
        !artifactRefValidUtc(owner$process_started_at)) {
        projectRegistryAbort(
            "project registry process start timestamp is malformed",
            "multischolar_malformed_registry_owner"
        )
    }
    artifactRefValidateId(owner$session_id, "session_id", "session")
    artifactRefValidateId(owner$owner_token, "owner_token", "owner")
    owner$schema_version <- as.integer(owner$schema_version)
    owner$pid <- as.integer(owner$pid)
    owner
}

projectRegistryReadOwner <- function(registry, missing_ok = FALSE) {
    path <- projectRegistryPath(registry, "owner")
    if (!file.exists(path)) {
        if (isTRUE(missing_ok)) return(NULL)
        projectRegistryAbort(
            "project registry writer owner metadata is missing",
            "multischolar_missing_registry_owner"
        )
    }
    owner <- tryCatch(
        jsonlite::read_json(path, simplifyVector = TRUE),
        error = \(error) projectRegistryAbort(
            "project registry writer owner metadata cannot be read",
            "multischolar_malformed_registry_owner",
            parent = error
        )
    )
    validateProjectRegistryOwner(owner)
}

projectRegistryWriteOwner <- function(registry, owner) {
    owner <- validateProjectRegistryOwner(owner)
    owner_path <- projectRegistryPath(registry, "owner")
    owner_directory <- dirname(owner_path)
    temporary <- tempfile("project-registry-owner-", tmpdir = owner_directory)
    on.exit(unlink(temporary, force = FALSE), add = TRUE)
    encoded <- jsonlite::toJSON(
        owner,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA,
        pretty = TRUE
    )
    connection <- file(temporary, open = "wb")
    tryCatch(
        writeBin(charToRaw(as.character(encoded)), connection),
        error = \(error) projectRegistryAbort(
            "project registry writer owner metadata could not be written",
            "multischolar_registry_write_failed",
            parent = error
        ),
        finally = close(connection)
    )
    if (file.exists(owner_path) && unlink(owner_path, force = FALSE) != 0L) {
        projectRegistryAbort(
            "stale project registry owner metadata could not be replaced",
            "multischolar_registry_write_failed"
        )
    }
    if (!file.rename(temporary, owner_path)) {
        projectRegistryAbort(
            "project registry writer owner metadata could not be published",
            "multischolar_registry_write_failed"
        )
    }
    invisible(owner_path)
}

projectRegistryAcquireWriter <- function(registry, timeout_ms = NULL) {
    registry <- validateProjectRegistry(registry)
    projectRegistryRequireDependencies()
    projectRegistryAssertLocalFilesystem(registry$filesystem_type)
    projectRegistryEnsureParent(registry, "lock")
    lock_path <- projectRegistryPath(registry, "lock")
    if (is.null(timeout_ms)) timeout_ms <- registry$resource_policy$writer_timeout_ms
    timeout_ms <- as.integer(projectRegistryNonnegativeInteger(
        timeout_ms,
        "writer_timeout_ms",
        maximum = 3600000L
    ))
    writer_key <- projectRegistryWriterKey(registry)
    if (exists(
        writer_key,
        envir = .PROJECT_REGISTRY_ACTIVE_WRITERS,
        inherits = FALSE
    )) {
        projectRegistryAbort(
            "project registry already has an active in-process writer",
            "multischolar_registry_writer_busy"
        )
    }
    lock_handle <- tryCatch(
        filelock::lock(lock_path, exclusive = TRUE, timeout = timeout_ms),
        error = \(error) projectRegistryAbort(
            "project registry writer guard could not be acquired",
            "multischolar_registry_writer_busy",
            parent = error
        )
    )
    if (is.null(lock_handle)) {
        projectRegistryAbort(
            "project registry already has an active writer",
            "multischolar_registry_writer_busy"
        )
    }
    release_on_error <- TRUE
    on.exit({
        if (isTRUE(release_on_error)) filelock::unlock(lock_handle)
    }, add = TRUE)
    previous <- projectRegistryReadOwner(registry, missing_ok = TRUE)
    if (!is.null(previous) && !identical(previous$host, projectRegistryHost())) {
        projectRegistryAbort(
            "cross-host project registry writer recovery is unsupported",
            "multischolar_cross_host_registry_owner",
            owner_host = previous$host,
            current_host = projectRegistryHost()
        )
    }
    if (!is.null(previous) && !identical(previous$project_id, registry$project_id)) {
        projectRegistryAbort(
            "project registry owner belongs to another project",
            "multischolar_cross_project_registry_owner"
        )
    }
    owner <- newProjectRegistryOwner(registry)
    projectRegistryWriteOwner(registry, owner)
    assign(
        writer_key,
        owner$owner_token,
        envir = .PROJECT_REGISTRY_ACTIVE_WRITERS
    )
    release_on_error <- FALSE
    structure(
        list(
            registry = registry,
            lock_handle = lock_handle,
            owner = owner
        ),
        class = c("MultiScholaRProjectRegistryWriter", "list")
    )
}

validateProjectRegistryWriter <- function(writer) {
    valid <- inherits(writer, "MultiScholaRProjectRegistryWriter") &&
        identical(names(writer), c("registry", "lock_handle", "owner"))
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry writer guard is invalid",
            "multischolar_invalid_registry_writer"
        )
    }
    validateProjectRegistry(writer$registry)
    validateProjectRegistryOwner(writer$owner)
    writer
}

projectRegistryReleaseWriter <- function(writer, owner_token) {
    writer <- validateProjectRegistryWriter(writer)
    artifactRefValidateId(owner_token, "owner_token", "owner")
    current <- projectRegistryReadOwner(writer$registry)
    expected <- writer$owner$owner_token
    if (!identical(owner_token, expected) ||
        !identical(current$owner_token, expected)) {
        projectRegistryAbort(
            "project registry writer release token does not match its owner",
            "multischolar_registry_owner_token_mismatch"
        )
    }
    owner_path <- projectRegistryPath(writer$registry, "owner", must_exist = TRUE)
    if (unlink(owner_path, force = FALSE) != 0L) {
        projectRegistryAbort(
            "project registry writer owner metadata could not be removed",
            "multischolar_registry_write_failed"
        )
    }
    unlock_error <- NULL
    unlocked <- tryCatch(
        filelock::unlock(writer$lock_handle),
        error = \(error) {
            unlock_error <<- error
            FALSE
        }
    )
    if (!isTRUE(unlocked)) {
        try(projectRegistryWriteOwner(writer$registry, writer$owner), silent = TRUE)
        projectRegistryAbort(
            "project registry writer guard could not be released",
            "multischolar_registry_writer_release_failed",
            parent = unlock_error
        )
    }
    writer_key <- projectRegistryWriterKey(writer$registry)
    active_token <- get0(
        writer_key,
        envir = .PROJECT_REGISTRY_ACTIVE_WRITERS,
        inherits = FALSE
    )
    if (identical(active_token, expected)) {
        rm(list = writer_key, envir = .PROJECT_REGISTRY_ACTIVE_WRITERS)
    }
    invisible(TRUE)
}
