# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_PATH_SCHEMA_VERSION <- 1L
.ARTIFACT_MANIFEST_SCHEMA <- "multischolar.artifact_manifest"
.ARTIFACT_MANIFEST_SCHEMA_VERSION <- 1L
.ARTIFACT_WINDOWS_RESERVED <- c(
    "CON", "PRN", "AUX", "NUL",
    paste0("COM", 1:9),
    paste0("LPT", 1:9)
)

artifactPathAbort <- function(message, class) {
    rlang::abort(message, class = c(class, "multischolar_artifact_path_error"))
}

artifactNormalizeProjectRoot <- function(project_root) {
    if (!workflowCapabilityScalarString(project_root) ||
        !dir.exists(project_root)) {
        artifactPathAbort(
            "artifact project root must be an existing directory",
            "multischolar_invalid_project_root"
        )
    }
    normalizePath(project_root, winslash = "/", mustWork = TRUE)
}

artifactValidatePathComponent <- function(value, field) {
    invalid <- !workflowCapabilityScalarString(value) ||
        !identical(value, trimws(value)) ||
        value %in% c(".", "..") ||
        grepl("[/\\\\<>:\"|?*]", value) ||
        grepl("[[:cntrl:]]", value) ||
        toupper(sub("[.].*$", "", value)) %in% .ARTIFACT_WINDOWS_RESERVED
    if (isTRUE(invalid)) {
        artifactPathAbort(
            sprintf("artifact path component '%s' is invalid", field),
            "multischolar_invalid_artifact_path_component"
        )
    }
    value
}

artifactNormalizeRelativePath <- function(relative_path) {
    if (!workflowCapabilityScalarString(relative_path) ||
        grepl("^[A-Za-z]:", relative_path) ||
        startsWith(relative_path, "/") ||
        grepl("[[:cntrl:]]", relative_path)) {
        artifactPathAbort(
            "artifact metadata path must be project-relative",
            "multischolar_invalid_relative_artifact_path"
        )
    }
    normalized <- gsub("\\\\", "/", relative_path)
    components <- strsplit(normalized, "/", fixed = TRUE)[[1L]]
    if (length(components) == 0L || any(!nzchar(components)) ||
        any(components %in% c(".", ".."))) {
        artifactPathAbort(
            "artifact metadata path contains traversal or empty components",
            "multischolar_invalid_relative_artifact_path"
        )
    }
    for (index in seq_along(components)) {
        artifactValidatePathComponent(components[[index]], paste0("path[", index, "]"))
    }
    paste(components, collapse = "/")
}

artifactPathIsContained <- function(path, project_root) {
    comparable_path <- path
    comparable_root <- project_root
    if (.Platform$OS.type == "windows") {
        comparable_path <- tolower(comparable_path)
        comparable_root <- tolower(comparable_root)
    }
    identical(comparable_path, comparable_root) ||
        startsWith(comparable_path, paste0(comparable_root, "/"))
}

artifactResolveContainedPath <- function(
    project_root,
    relative_path,
    must_exist = FALSE
) {
    root <- artifactNormalizeProjectRoot(project_root)
    relative <- artifactNormalizeRelativePath(relative_path)
    components <- strsplit(relative, "/", fixed = TRUE)[[1L]]
    candidate <- do.call(file.path, as.list(c(root, components)))
    candidate <- gsub("\\\\", "/", candidate)
    if (isTRUE(must_exist) && !file.exists(candidate) && !dir.exists(candidate)) {
        artifactPathAbort(
            sprintf("artifact path does not exist: %s", relative),
            "multischolar_missing_artifact_path"
        )
    }
    ancestor <- candidate
    suffix <- character()
    while (!file.exists(ancestor) && !dir.exists(ancestor)) {
        parent <- dirname(ancestor)
        if (identical(parent, ancestor)) break
        suffix <- c(basename(ancestor), suffix)
        ancestor <- parent
    }
    resolved_ancestor <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
    if (!artifactPathIsContained(resolved_ancestor, root)) {
        artifactPathAbort(
            "artifact path escapes the project root through a symlink",
            "multischolar_artifact_path_escape"
        )
    }
    resolved <- if (length(suffix) == 0L) {
        resolved_ancestor
    } else {
        do.call(file.path, as.list(c(resolved_ancestor, suffix)))
    }
    resolved <- gsub("\\\\", "/", resolved)
    if (!artifactPathIsContained(resolved, root)) {
        artifactPathAbort(
            "artifact path escapes the project root",
            "multischolar_artifact_path_escape"
        )
    }
    resolved
}

artifactWorkflowRelativePaths <- function(identity) {
    omic_label <- artifactValidatePathComponent(identity$omic_label, "omic_label")
    workflow_slug <- artifactValidatePathComponent(
        identity$workflow_slug,
        "workflow_slug"
    )
    data_root <- file.path("data", omic_label, workflow_slug)
    workflow_state_root <- file.path("state", omic_label, workflow_slug)
    paths <- list(
        data_root = data_root,
        tables = file.path(data_root, "tables"),
        cache = file.path(data_root, "cache"),
        staging = file.path(data_root, "staging"),
        trash = file.path(data_root, "trash"),
        state_root = "state",
        registry = file.path("state", "multischolar.duckdb"),
        locks = file.path("state", "locks"),
        duckdb_tmp = file.path("state", "tmp", "duckdb"),
        workflow_state_root = workflow_state_root,
        generations = file.path(workflow_state_root, "generations"),
        artifact_manifest = file.path(workflow_state_root, "artifact_manifest.json"),
        workflow_state = file.path(workflow_state_root, "workflow_state.json")
    )
    lapply(paths, artifactNormalizeRelativePath)
}

buildArtifactPaths <- function(project_root, identity) {
    required <- c("omic_type", "omic_label", "workflow_slug")
    if (!all(required %in% names(identity)) ||
        !all(vapply(identity[required], workflowCapabilityScalarString, logical(1)))) {
        artifactPathAbort(
            "artifact path identity is incomplete",
            "multischolar_invalid_workflow_identity"
        )
    }
    root <- artifactNormalizeProjectRoot(project_root)
    relative_paths <- artifactWorkflowRelativePaths(identity)
    legacy_relative <- artifactNormalizeRelativePath(
        file.path("data", artifactValidatePathComponent(identity$omic_type, "omic_type"))
    )
    structure(
        list(
            schema_version = .ARTIFACT_PATH_SCHEMA_VERSION,
            project_root = root,
            labels = list(
                omic_type = identity$omic_type,
                omic_label = identity$omic_label,
                workflow_slug = identity$workflow_slug
            ),
            relative_paths = relative_paths,
            legacy_read_fallbacks = list(
                data_root = list(
                    relative_path = legacy_relative,
                    access = "read_only",
                    migration = "explicit_required"
                )
            )
        ),
        class = c("MultiScholaRArtifactPaths", "list")
    )
}

artifactPath <- function(paths, path_name, must_exist = FALSE) {
    if (!inherits(paths, "MultiScholaRArtifactPaths") ||
        !workflowCapabilityScalarString(path_name) ||
        !path_name %in% names(paths$relative_paths)) {
        artifactPathAbort(
            "requested artifact path label is invalid",
            "multischolar_unknown_artifact_path"
        )
    }
    artifactResolveContainedPath(
        paths$project_root,
        paths$relative_paths[[path_name]],
        must_exist = must_exist
    )
}

artifactPathMetadata <- function(paths) {
    if (!inherits(paths, "MultiScholaRArtifactPaths")) {
        artifactPathAbort(
            "artifact paths object is invalid",
            "multischolar_invalid_artifact_paths"
        )
    }
    list(
        schema_version = paths$schema_version,
        labels = paths$labels,
        relative_paths = paths$relative_paths,
        legacy_read_fallbacks = paths$legacy_read_fallbacks
    )
}

artifactDirectoryHasContent <- function(path) {
    dir.exists(path) && length(list.files(
        path,
        all.files = TRUE,
        no.. = TRUE,
        recursive = FALSE
    )) > 0L
}

detectWorkflowProjectState <- function(project_root, identity) {
    root <- artifactNormalizeProjectRoot(project_root)
    relative <- artifactWorkflowRelativePaths(identity)
    manifest_path <- artifactResolveContainedPath(root, relative$artifact_manifest)
    workflow_root <- artifactResolveContainedPath(root, relative$workflow_state_root)
    data_root <- artifactResolveContainedPath(root, relative$data_root)
    if (file.exists(manifest_path)) {
        parse_error <- NULL
        manifest <- tryCatch(
            jsonlite::read_json(manifest_path, simplifyVector = TRUE),
            error = \(error) {
                parse_error <<- conditionMessage(error)
                NULL
            }
        )
        if (is.null(manifest) || !is.list(manifest)) {
            return(list(status = "artifact_corrupt", reason = parse_error))
        }
        version <- workflowStateVersionValue(manifest$schema_version)
        if (is.na(version) || !identical(manifest$schema, .ARTIFACT_MANIFEST_SCHEMA)) {
            return(list(status = "artifact_corrupt", reason = "invalid_manifest_schema"))
        }
        if (version > .ARTIFACT_MANIFEST_SCHEMA_VERSION) {
            return(list(
                status = "artifact_future_schema",
                schema_version = version,
                reason = "future_manifest_schema"
            ))
        }
        if (!identical(version, .ARTIFACT_MANIFEST_SCHEMA_VERSION)) {
            return(list(status = "artifact_corrupt", reason = "unsupported_manifest_schema"))
        }
        return(list(
            status = "artifact_valid",
            schema_version = version,
            reason = "valid_manifest"
        ))
    }
    # The registry is project-wide; only workflow-scoped evidence can make this
    # exact context corrupt when its manifest is absent.
    if (artifactDirectoryHasContent(workflow_root) ||
        artifactDirectoryHasContent(data_root)) {
        return(list(status = "artifact_corrupt", reason = "artifact_evidence_without_manifest"))
    }
    legacy_relative <- artifactNormalizeRelativePath(file.path("data", identity$omic_type))
    legacy_path <- artifactResolveContainedPath(root, legacy_relative)
    if (artifactDirectoryHasContent(legacy_path)) {
        return(list(status = "legacy_memory", reason = "legacy_data_root_has_content"))
    }
    list(status = "new", reason = "no_project_state")
}
