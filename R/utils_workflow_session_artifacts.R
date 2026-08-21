# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.WORKFLOW_SESSION_SCHEMA <- "multischolar.workflow_session"
.WORKFLOW_SESSION_VERSION <- 1L
.WORKFLOW_SESSION_METADATA_LIMIT <- 1024L * 1024L

workflowSessionAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_workflow_session_error"),
        ...
    )
}

workflowSessionProjectRelativePath <- function(path, project_root) {
    root <- artifactNormalizeProjectRoot(project_root)
    if (!workflowCapabilityScalarString(path)) {
        workflowSessionAbort(
            "workflow session path must be one string",
            "multischolar_invalid_workflow_session_path"
        )
    }
    expanded <- path.expand(path)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
        expanded <- file.path(root, expanded)
    }
    absolute <- normalizePath(expanded, winslash = "/", mustWork = FALSE)
    comparable_root <- root
    comparable_path <- absolute
    if (.Platform$OS.type == "windows") {
        comparable_root <- tolower(comparable_root)
        comparable_path <- tolower(comparable_path)
    }
    prefix <- paste0(comparable_root, "/")
    if (!startsWith(comparable_path, prefix)) {
        workflowSessionAbort(
            "workflow session path is outside the project root",
            "multischolar_workflow_session_path_escape",
            path = path
        )
    }
    relative <- substring(absolute, nchar(root) + 2L)
    artifactNormalizeRelativePath(relative)
}

workflowSessionEncodeMetadata <- function(value, owner = "session metadata") {
    artifactWorkflowStateAssertSafeMetadata(value, owner)
    encoded <- as.character(jsonlite::serializeJSON(value, digits = NA))
    if (nchar(encoded, type = "bytes") > .WORKFLOW_SESSION_METADATA_LIMIT) {
        workflowSessionAbort(
            sprintf("%s exceeds the portable session metadata limit", owner),
            "multischolar_workflow_session_metadata_too_large",
            owner = owner
        )
    }
    encoded
}

workflowSessionDecodeMetadata <- function(value, owner = "session metadata") {
    if (!workflowCapabilityScalarString(value) || !jsonlite::validate(value)) {
        workflowSessionAbort(
            sprintf("%s is not valid serialized metadata", owner),
            "multischolar_malformed_workflow_session",
            owner = owner
        )
    }
    decoded <- tryCatch(
        jsonlite::unserializeJSON(value),
        error = function(error) workflowSessionAbort(
            sprintf("%s could not be decoded", owner),
            "multischolar_malformed_workflow_session",
            owner = owner,
            parent = error
        )
    )
    artifactWorkflowStateAssertSafeMetadata(decoded, owner)
    decoded
}

workflowSessionContentDigest <- function(manifest) {
    candidate <- manifest
    candidate$manifest_digest <- NULL
    canonical <- jsonlite::fromJSON(
        workflowSessionJson(unclass(candidate)),
        simplifyVector = FALSE
    )
    digest::digest(
        workflowSessionJson(canonical),
        algo = .artifactHashAlgorithm,
        serialize = FALSE
    )
}

workflowSessionValidateLineage <- function(lineage, current_generation_id) {
    if (!is.list(lineage) || length(lineage) == 0L) {
        workflowSessionAbort(
            "workflow session active lineage is empty",
            "multischolar_invalid_workflow_session_lineage"
        )
    }
    required <- c(
        "generation_id", "logical_name", "manifest_relative_path",
        "manifest_digest"
    )
    for (entry in lineage) {
        valid <- is.list(entry) && identical(names(entry), required) &&
            workflowCapabilityScalarString(entry$generation_id) &&
            workflowCapabilityScalarString(entry$logical_name) &&
            workflowCapabilityScalarString(entry$manifest_relative_path)
        if (!isTRUE(valid)) {
            workflowSessionAbort(
                "workflow session active lineage entry is invalid",
                "multischolar_invalid_workflow_session_lineage"
            )
        }
        artifactNormalizeRelativePath(entry$manifest_relative_path)
        artifactRefValidateDigest(entry$manifest_digest, "manifest_digest")
    }
    generation_ids <- vapply(lineage, `[[`, character(1), "generation_id")
    if (anyDuplicated(generation_ids) > 0L ||
        !identical(tail(generation_ids, 1L), current_generation_id)) {
        workflowSessionAbort(
            "workflow session current generation is not the unique lineage tail",
            "multischolar_invalid_workflow_session_lineage"
        )
    }
    invisible(TRUE)
}

workflowSessionValidateFingerprint <- function(fingerprint) {
    required <- c(
        "relative_path", "byte_digest", "generation_id",
        "state_semantic_digest"
    )
    valid <- is.list(fingerprint) && identical(names(fingerprint), required) &&
        workflowCapabilityScalarString(fingerprint$relative_path) &&
        workflowCapabilityScalarString(fingerprint$generation_id)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "workflow session compatibility fingerprint is invalid",
            "multischolar_invalid_workflow_session_fingerprint"
        )
    }
    artifactNormalizeRelativePath(fingerprint$relative_path)
    artifactRefValidateDigest(fingerprint$byte_digest, "byte_digest")
    artifactRefValidateDigest(
        fingerprint$state_semantic_digest,
        "state_semantic_digest"
    )
    invisible(TRUE)
}

workflowSessionValidateState <- function(state) {
    required <- c(
        "schema", "schema_version", "backend", "current_generation_id",
        "current_state", "active_lineage", "workflow_type", "audit_enabled",
        "state_semantic_digest"
    )
    version <- workflowStateVersionValue(state$schema_version)
    valid <- is.list(state) && identical(names(state), required) &&
        identical(state$schema, ARTIFACT_WORKFLOW_STATE_SCHEMA) &&
        identical(version, ARTIFACT_WORKFLOW_STATE_VERSION) &&
        identical(state$backend, "artifact") &&
        workflowCapabilityScalarString(state$current_generation_id) &&
        workflowCapabilityScalarString(state$current_state) &&
        is.logical(state$audit_enabled) && length(state$audit_enabled) == 1L &&
        !is.na(state$audit_enabled)
    if (!isTRUE(valid)) {
        workflowSessionAbort(
            "workflow session state snapshot is invalid",
            "multischolar_invalid_workflow_session_state"
        )
    }
    artifactRefValidateDigest(
        state$state_semantic_digest,
        "state_semantic_digest"
    )
    invisible(TRUE)
}

validateWorkflowSessionManifest <- function(manifest) {
    required <- c(
        "schema", "schema_version", "backend", "identity", "descriptor",
        "workflow_state", "stage_refs", "plot_refs", "metadata_json",
        "compatibility", "exported_at", "manifest_digest"
    )
    version <- workflowStateVersionValue(manifest$schema_version)
    if (!is.list(manifest) || !identical(names(manifest), required) ||
        !identical(manifest$schema, .WORKFLOW_SESSION_SCHEMA) ||
        is.na(version) || version > .WORKFLOW_SESSION_VERSION ||
        !identical(version, .WORKFLOW_SESSION_VERSION) ||
        !identical(manifest$backend, "artifact")) {
        workflowSessionAbort(
            "workflow session schema or version is unsupported",
            "multischolar_unsupported_workflow_session"
        )
    }
    identity_fields <- c(
        "project_id", "workflow_id", "omic_type", "omic_label",
        "workflow_slug"
    )
    descriptor_fields <- c(
        "descriptor_id", "descriptor_version", "descriptor_digest"
    )
    valid_identity <- is.list(manifest$identity) &&
        identical(names(manifest$identity), identity_fields) &&
        all(vapply(
            manifest$identity,
            workflowCapabilityScalarString,
            logical(1)
        ))
    valid_descriptor <- is.list(manifest$descriptor) &&
        identical(names(manifest$descriptor), descriptor_fields) &&
        workflowCapabilityScalarString(manifest$descriptor$descriptor_id) &&
        workflowCapabilityScalarString(
            manifest$descriptor$descriptor_version
        )
    if (!isTRUE(valid_identity) || !isTRUE(valid_descriptor)) {
        workflowSessionAbort(
            "workflow session identity or descriptor is invalid",
            "multischolar_invalid_workflow_session_identity"
        )
    }
    artifactRefValidateDigest(
        manifest$descriptor$descriptor_digest,
        "descriptor_digest"
    )
    workflowSessionValidateState(manifest$workflow_state)
    workflowSessionValidateLineage(
        manifest$workflow_state$active_lineage,
        manifest$workflow_state$current_generation_id
    )
    if (!is.list(manifest$stage_refs) || !is.list(manifest$plot_refs)) {
        workflowSessionAbort(
            "workflow session references are invalid",
            "multischolar_invalid_workflow_session_refs"
        )
    }
    workflowSessionDecodeMetadata(manifest$metadata_json)
    workflowSessionValidateFingerprint(manifest$compatibility)
    if (!artifactRefValidUtc(manifest$exported_at)) {
        workflowSessionAbort(
            "workflow session export timestamp is invalid",
            "multischolar_invalid_workflow_session_timestamp"
        )
    }
    artifactRefValidateDigest(manifest$manifest_digest, "manifest_digest")
    if (!identical(
        manifest$manifest_digest,
        workflowSessionContentDigest(manifest)
    )) {
        workflowSessionAbort(
            "workflow session manifest digest does not match its contents",
            "multischolar_workflow_session_digest_mismatch"
        )
    }
    manifest$schema_version <- version
    manifest$workflow_state$schema_version <- workflowStateVersionValue(
        manifest$workflow_state$schema_version
    )
    manifest
}

workflowSessionJson <- function(value) {
    as.character(jsonlite::toJSON(
        value,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA,
        pretty = TRUE
    ))
}

workflowSessionWriteCandidate <- function(
    manifest,
    path,
    failure_injector = NULL
) {
    manifest <- validateWorkflowSessionManifest(manifest)
    parent <- dirname(path)
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
        workflowSessionAbort(
            "workflow session directory could not be created",
            "multischolar_workflow_session_write_failed"
        )
    }
    temporary <- file.path(
        parent,
        paste0(".", basename(path), ".", artifactOpaqueId("tmp"), ".tmp")
    )
    keep_temporary <- FALSE
    on.exit({
        if (!keep_temporary && file.exists(temporary)) {
            unlink(temporary, force = FALSE)
        }
    }, add = TRUE)
    artifactStoreInvokeFailure(
        failure_injector,
        "before_session_write",
        list(path = path)
    )
    connection <- file(temporary, open = "wb")
    tryCatch(
        writeBin(
            charToRaw(paste0(
                workflowSessionJson(unclass(manifest)),
                "\n"
            )),
            connection
        ),
        finally = close(connection)
    )
    decoded <- jsonlite::read_json(temporary, simplifyVector = FALSE)
    validateWorkflowSessionManifest(decoded)
    artifactStoreInvokeFailure(
        failure_injector,
        "after_session_validation",
        list(path = path)
    )
    keep_temporary <- TRUE
    temporary
}

workflowSessionPublishImmutable <- function(
    manifest,
    path,
    failure_injector = NULL
) {
    if (file.exists(path) || dir.exists(path)) {
        workflowSessionAbort(
            "immutable workflow session already exists",
            "multischolar_workflow_session_already_exists",
            path = path
        )
    }
    temporary <- workflowSessionWriteCandidate(
        manifest,
        path,
        failure_injector
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "immutable workflow session could not be published",
            "multischolar_workflow_session_publish_failed",
            path = path
        )
    }
    invisible(path)
}

workflowSessionReplaceLatest <- function(
    manifest,
    path,
    failure_injector = NULL
) {
    temporary <- workflowSessionWriteCandidate(
        manifest,
        path,
        failure_injector
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = FALSE), add = TRUE)
    artifactStoreInvokeFailure(
        failure_injector,
        "before_latest_publication",
        list(path = path)
    )
    if (!isTRUE(file.rename(temporary, path))) {
        workflowSessionAbort(
            "latest workflow session could not be atomically replaced",
            "multischolar_workflow_session_publish_failed",
            path = path
        )
    }
    invisible(path)
}

readWorkflowSessionManifest <- function(path) {
    if (!workflowCapabilityScalarString(path) || !file.exists(path) ||
        dir.exists(path)) {
        workflowSessionAbort(
            "workflow session manifest does not exist",
            "multischolar_missing_workflow_session",
            path = path
        )
    }
    manifest <- tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(error) workflowSessionAbort(
            "workflow session manifest could not be parsed",
            "multischolar_malformed_workflow_session",
            path = path,
            parent = error
        )
    )
    validateWorkflowSessionManifest(manifest)
}
