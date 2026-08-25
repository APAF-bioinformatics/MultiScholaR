# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

ARTIFACT_REGISTRY_SCOPE_SCHEMA <- "multischolar.artifact_registry_scope"
ARTIFACT_REGISTRY_SCOPE_VERSION <- 1L

#' Build the exact registry identifier for one workflow tuple
#'
#' @param identity A complete bound workflow identity.
#'
#' @return A deterministic registry workflow identifier.
#' @noRd
artifactRegistryScopedWorkflowId <- function(identity) {
    key_digest <- digest::digest(
        workflowCapabilityKey(identity),
        algo = "sha256",
        serialize = FALSE
    )
    paste(
        identity$workflow_id,
        identity$workflow_slug,
        substr(key_digest, 1L, 16L),
        sep = "::"
    )
}

#' Resolve the registry-scope sidecar path
#'
#' @param store A validated artifact store.
#'
#' @return A project-relative sidecar path.
#' @noRd
artifactRegistryScopePath <- function(store) {
    store <- validateArtifactStore(store)
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$workflow_state_root,
        "registry-scope.json"
    ))
}

#' Validate immutable registry-scope metadata
#'
#' @param scope Decoded registry-scope metadata.
#' @param identity The public bound workflow identity.
#'
#' @return The validated scope metadata.
#' @noRd
validateArtifactRegistryScope <- function(scope, identity) {
    required <- c(
        "schema", "schema_version", "project_id", "public_workflow_id",
        "registry_workflow_id", "omic_type", "omic_label", "workflow_slug",
        "capability_key_digest", "created_at"
    )
    expected_digest <- digest::digest(
        workflowCapabilityKey(identity),
        algo = "sha256",
        serialize = FALSE
    )
    valid <- is.list(scope) && identical(names(scope), required) &&
        identical(scope$schema, ARTIFACT_REGISTRY_SCOPE_SCHEMA) &&
        identical(
            as.integer(scope$schema_version),
            ARTIFACT_REGISTRY_SCOPE_VERSION
        ) && identical(scope$project_id, identity$project_id) &&
        identical(scope$public_workflow_id, identity$workflow_id) &&
        identical(
            scope$registry_workflow_id,
            artifactRegistryScopedWorkflowId(identity)
        ) && identical(scope$omic_type, identity$omic_type) &&
        identical(scope$omic_label, identity$omic_label) &&
        identical(scope$workflow_slug, identity$workflow_slug) &&
        identical(scope$capability_key_digest, expected_digest) &&
        artifactRefValidUtc(scope$created_at)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact registry scope is incompatible with its workflow",
            "multischolar_incompatible_artifact_registry_scope"
        )
    }
    scope$schema_version <- as.integer(scope$schema_version)
    scope
}

#' Create immutable registry-scope metadata
#'
#' @param identity The public bound workflow identity.
#'
#' @return A validated registry-scope metadata list.
#' @noRd
newArtifactRegistryScope <- function(identity) {
    validateArtifactRegistryScope(
        list(
            schema = ARTIFACT_REGISTRY_SCOPE_SCHEMA,
            schema_version = ARTIFACT_REGISTRY_SCOPE_VERSION,
            project_id = identity$project_id,
            public_workflow_id = identity$workflow_id,
            registry_workflow_id = artifactRegistryScopedWorkflowId(identity),
            omic_type = identity$omic_type,
            omic_label = identity$omic_label,
            workflow_slug = identity$workflow_slug,
            capability_key_digest = digest::digest(
                workflowCapabilityKey(identity),
                algo = "sha256",
                serialize = FALSE
            ),
            created_at = artifactRefUtcNow()
        ),
        identity
    )
}

#' Resolve the registry identity for an artifact workflow
#'
#' Existing roots without a scope sidecar retain their legacy workflow ID.
#'
#' @param store A validated artifact store.
#' @param identity The public bound workflow identity.
#' @param create_scope Whether a new artifact root must publish scoped metadata.
#'
#' @return A workflow identity whose workflow ID is registry-specific.
#' @noRd
artifactRegistryIdentity <- function(store, identity, create_scope = FALSE) {
    store <- validateArtifactStore(store)
    if (!is.logical(create_scope) || length(create_scope) != 1L ||
        is.na(create_scope)) {
        artifactWorkflowStateAbort(
            "artifact registry scope creation flag must be true or false",
            "multischolar_invalid_artifact_registry_scope"
        )
    }
    relative_path <- artifactRegistryScopePath(store)
    path <- artifactStoreResolveFile(store, relative_path)
    scope <- NULL
    if (file.exists(path)) {
        scope <- validateArtifactRegistryScope(
            artifactStoreReadJson(store, relative_path),
            identity
        )
    } else if (isTRUE(create_scope)) {
        scope <- newArtifactRegistryScope(identity)
        temporary_path <- artifactNormalizeRelativePath(paste0(
            relative_path,
            ".",
            artifactOpaqueId("tmp"),
            ".tmp"
        ))
        artifactStoreWriteJson(store, scope, temporary_path)
        on.exit({
            temporary <- artifactStoreResolveFile(store, temporary_path)
            if (file.exists(temporary)) unlink(temporary, force = FALSE)
        }, add = TRUE)
        artifactStorePublishFile(store, temporary_path, relative_path)
    }
    registry_identity <- identity
    if (!is.null(scope)) {
        registry_identity$workflow_id <- scope$registry_workflow_id
    }
    registry_identity
}
