# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

PROJECT_REGISTRY_SCHEMA_ID <- "multischolar.project_registry"
PROJECT_REGISTRY_SCHEMA_VERSION <- 1L
PROJECT_REGISTRY_REQUIRED_PACKAGES <- c("DBI", "duckdb", "filelock")
PROJECT_REGISTRY_UNSUPPORTED_FILESYSTEMS <- c(
    "9p", "afs", "cifs", "fuse.sshfs", "nfs", "nfs4", "smbfs", "sshfs"
)

projectRegistryAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_project_registry_error"),
        ...
    )
}

projectRegistryAssertDependencies <- function(availability) {
    valid <- is.logical(availability) && !is.null(names(availability)) &&
        identical(names(availability), PROJECT_REGISTRY_REQUIRED_PACKAGES) &&
        !anyNA(availability)
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry dependency availability is malformed",
            "multischolar_invalid_registry_dependency_check"
        )
    }
    missing <- names(availability)[!availability]
    if (length(missing) > 0L) {
        projectRegistryAbort(
            sprintf(
                paste0(
                    "artifact mode requires optional packages: %s; ",
                    "install them before initializing artifact storage"
                ),
                paste(missing, collapse = ", ")
            ),
            "multischolar_missing_artifact_dependencies",
            missing_packages = missing
        )
    }
    invisible(TRUE)
}

projectRegistryRequireDependencies <- function() {
    availability <- vapply(
        PROJECT_REGISTRY_REQUIRED_PACKAGES,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
    projectRegistryAssertDependencies(availability)
}

projectRegistryPositiveInteger <- function(value, field, maximum = .Machine$integer.max) {
    valid <- length(value) == 1L && is.numeric(value) && !is.na(value) &&
        is.finite(value) && value >= 1 && value == floor(value) && value <= maximum
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            sprintf("project registry policy '%s' must be a bounded positive integer", field),
            "multischolar_invalid_registry_resource_policy"
        )
    }
    as.numeric(value)
}

projectRegistryNonnegativeInteger <- function(value, field, maximum) {
    valid <- length(value) == 1L && is.numeric(value) && !is.na(value) &&
        is.finite(value) && value >= 0 && value == floor(value) && value <= maximum
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            sprintf(
                "project registry policy '%s' must be a bounded non-negative integer",
                field
            ),
            "multischolar_invalid_registry_resource_policy"
        )
    }
    as.numeric(value)
}

projectRegistryDefaultThreads <- function() {
    detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (length(detected) != 1L || is.na(detected) || detected < 1L) detected <- 1L
    as.integer(min(4L, detected))
}

normalizeProjectRegistryPolicy <- function(resource_policy = NULL) {
    if (is.null(resource_policy)) resource_policy <- list()
    if (!is.list(resource_policy) ||
        (length(resource_policy) > 0L && is.null(names(resource_policy)))) {
        projectRegistryAbort(
            "project registry resource policy must be a named list",
            "multischolar_invalid_registry_resource_policy"
        )
    }
    allowed <- c(
        "threads", "duckdb_memory_limit_bytes", "process_rss_limit_bytes",
        "temp_size_limit_bytes", "max_result_rows", "max_result_bytes",
        "preserve_insertion_order", "writer_timeout_ms"
    )
    unknown <- setdiff(names(resource_policy), allowed)
    if (length(unknown) > 0L || any(!nzchar(names(resource_policy)))) {
        projectRegistryAbort(
            "project registry resource policy contains unsupported settings",
            "multischolar_invalid_registry_resource_policy",
            unsupported_settings = unknown
        )
    }
    defaults <- list(
        threads = projectRegistryDefaultThreads(),
        duckdb_memory_limit_bytes = 512 * 1024^2,
        process_rss_limit_bytes = 2 * 1024^3,
        temp_size_limit_bytes = 2 * 1024^3,
        max_result_rows = 10000L,
        max_result_bytes = 64 * 1024^2,
        preserve_insertion_order = TRUE,
        writer_timeout_ms = 0L
    )
    policy <- utils::modifyList(defaults, resource_policy, keep.null = FALSE)
    policy$threads <- as.integer(projectRegistryPositiveInteger(
        policy$threads,
        "threads",
        maximum = 64L
    ))
    byte_fields <- c(
        "duckdb_memory_limit_bytes", "process_rss_limit_bytes",
        "temp_size_limit_bytes", "max_result_bytes"
    )
    for (field in byte_fields) {
        policy[[field]] <- projectRegistryPositiveInteger(
            policy[[field]],
            field,
            maximum = 2^53
        )
    }
    policy$max_result_rows <- as.integer(projectRegistryPositiveInteger(
        policy$max_result_rows,
        "max_result_rows",
        maximum = 1000000L
    ))
    policy$writer_timeout_ms <- projectRegistryNonnegativeInteger(
        policy$writer_timeout_ms,
        "writer_timeout_ms",
        maximum = 3600000L
    )
    if (!is.logical(policy$preserve_insertion_order) ||
        length(policy$preserve_insertion_order) != 1L ||
        is.na(policy$preserve_insertion_order)) {
        projectRegistryAbort(
            "project registry insertion-order policy must be true or false",
            "multischolar_invalid_registry_resource_policy"
        )
    }
    if (policy$process_rss_limit_bytes <= policy$duckdb_memory_limit_bytes ||
        policy$max_result_bytes >= policy$process_rss_limit_bytes) {
        projectRegistryAbort(
            "process RSS must exceed DuckDB memory and bounded result limits",
            "multischolar_invalid_registry_resource_policy"
        )
    }
    policy
}

projectRegistryFilesystemType <- function(project_root) {
    root <- artifactNormalizeProjectRoot(project_root)
    mount_info <- "/proc/self/mountinfo"
    if (!identical(Sys.info()[["sysname"]], "Linux") || !file.exists(mount_info)) {
        return("unknown_local")
    }
    lines <- readLines(mount_info, warn = FALSE)
    candidates <- lapply(lines, function(line) {
        halves <- strsplit(line, " - ", fixed = TRUE)[[1L]]
        if (length(halves) != 2L) return(NULL)
        left <- strsplit(halves[[1L]], " ", fixed = TRUE)[[1L]]
        right <- strsplit(halves[[2L]], " ", fixed = TRUE)[[1L]]
        if (length(left) < 5L || length(right) < 1L) return(NULL)
        mount_point <- gsub("\\\\040", " ", left[[5L]], fixed = TRUE)
        mount_point <- normalizePath(mount_point, winslash = "/", mustWork = FALSE)
        contained <- identical(mount_point, "/") ||
            artifactPathIsContained(root, mount_point)
        if (!isTRUE(contained)) return(NULL)
        list(path = mount_point, type = tolower(right[[1L]]))
    })
    candidates <- Filter(Negate(is.null), candidates)
    if (length(candidates) == 0L) return("unknown_local")
    lengths <- vapply(candidates, \(candidate) nchar(candidate$path), integer(1))
    candidates[[which.max(lengths)]]$type
}

projectRegistryAssertLocalFilesystem <- function(filesystem_type) {
    if (!workflowCapabilityScalarString(filesystem_type)) {
        projectRegistryAbort(
            "project registry filesystem type could not be determined",
            "multischolar_unsupported_registry_filesystem"
        )
    }
    normalized <- tolower(filesystem_type)
    unsupported <- normalized %in% PROJECT_REGISTRY_UNSUPPORTED_FILESYSTEMS ||
        grepl("^(fuse[.])?(nfs|smb|sshfs)", normalized)
    if (isTRUE(unsupported)) {
        projectRegistryAbort(
            sprintf("project registry writer locking is unsupported on '%s'", normalized),
            "multischolar_unsupported_registry_filesystem",
            filesystem_type = normalized
        )
    }
    invisible(normalized)
}

newProjectRegistry <- function(paths, project_id, resource_policy = NULL) {
    if (!inherits(paths, "MultiScholaRArtifactPaths") ||
        !workflowCapabilityScalarString(project_id)) {
        projectRegistryAbort(
            "project registry requires validated artifact paths and one project ID",
            "multischolar_invalid_project_registry"
        )
    }
    filesystem_type <- projectRegistryFilesystemType(paths$project_root)
    projectRegistryAssertLocalFilesystem(filesystem_type)
    registry_relative <- paths$relative_paths$registry
    state_relative <- paths$relative_paths$state_root
    structure(
        list(
            schema = PROJECT_REGISTRY_SCHEMA_ID,
            schema_version = PROJECT_REGISTRY_SCHEMA_VERSION,
            project_id = project_id,
            project_root = paths$project_root,
            filesystem_type = filesystem_type,
            relative_paths = list(
                database = registry_relative,
                lock = artifactNormalizeRelativePath(
                    file.path(paths$relative_paths$locks, "project-registry.lock")
                ),
                owner = artifactNormalizeRelativePath(
                    file.path(paths$relative_paths$locks, "project-registry.owner.json")
                ),
                temporary = paths$relative_paths$duckdb_tmp,
                backups = artifactNormalizeRelativePath(file.path(state_relative, "backups"))
            ),
            resource_policy = normalizeProjectRegistryPolicy(resource_policy)
        ),
        class = c("MultiScholaRProjectRegistry", "list")
    )
}

validateProjectRegistry <- function(registry) {
    required <- c(
        "schema", "schema_version", "project_id", "project_root",
        "filesystem_type", "relative_paths", "resource_policy"
    )
    valid <- inherits(registry, "MultiScholaRProjectRegistry") &&
        identical(names(registry), required) &&
        identical(registry$schema, PROJECT_REGISTRY_SCHEMA_ID) &&
        identical(registry$schema_version, PROJECT_REGISTRY_SCHEMA_VERSION) &&
        workflowCapabilityScalarString(registry$project_id) &&
        identical(
            artifactNormalizeProjectRoot(registry$project_root),
            registry$project_root
        )
    if (!isTRUE(valid)) {
        projectRegistryAbort(
            "project registry structure or version is invalid",
            "multischolar_invalid_project_registry"
        )
    }
    projectRegistryAssertLocalFilesystem(registry$filesystem_type)
    lapply(registry$relative_paths, artifactNormalizeRelativePath)
    normalizeProjectRegistryPolicy(registry$resource_policy)
    registry
}

projectRegistryPath <- function(registry, path_name, must_exist = FALSE) {
    registry <- validateProjectRegistry(registry)
    if (!workflowCapabilityScalarString(path_name) ||
        !path_name %in% names(registry$relative_paths)) {
        projectRegistryAbort(
            "requested project registry path is not package-owned",
            "multischolar_unknown_registry_path"
        )
    }
    artifactResolveContainedPath(
        registry$project_root,
        registry$relative_paths[[path_name]],
        must_exist = must_exist
    )
}

projectRegistryForContext <- function(context, resource_policy = NULL) {
    if (!inherits(context, "WorkflowContext") || !context$isBound()) {
        projectRegistryAbort(
            "project registry requires a bound workflow context",
            "multischolar_invalid_registry_context"
        )
    }
    decision <- context$getStorageDecision()
    if (identical(decision$effective_backend, "memory")) return(NULL)
    identity <- context$getIdentity()
    newProjectRegistry(
        context$getPaths(),
        identity$project_id,
        resource_policy = resource_policy
    )
}

projectRegistryProcessRss <- function() {
    status_path <- "/proc/self/status"
    if (file.exists(status_path)) {
        lines <- readLines(status_path, warn = FALSE)
        rss_line <- grep("^VmRSS:[[:space:]]+[0-9]+[[:space:]]+kB$", lines, value = TRUE)
        if (length(rss_line) == 1L) {
            kibibytes <- as.numeric(gsub("[^0-9]", "", rss_line))
            return(kibibytes * 1024)
        }
    }
    if (requireNamespace("ps", quietly = TRUE)) {
        return(as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]))
    }
    NA_real_
}

projectRegistryAssertRss <- function(registry, stage) {
    registry <- validateProjectRegistry(registry)
    rss <- projectRegistryProcessRss()
    if (!is.na(rss) && rss > registry$resource_policy$process_rss_limit_bytes) {
        projectRegistryAbort(
            sprintf("process RSS exceeds the project registry limit at '%s'", stage),
            "multischolar_registry_rss_limit_exceeded",
            rss_bytes = rss,
            limit_bytes = registry$resource_policy$process_rss_limit_bytes
        )
    }
    invisible(rss)
}
