# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_STATE_REFERENCE_HINT_SCHEMA <-
    "multischolar.artifact_state_reference_hint"
.ARTIFACT_STATE_REFERENCE_HINT_VERSION <- 1L
.ARTIFACT_STATE_REFERENCE_SCHEMA <- "multischolar.artifact_state_reference"
.ARTIFACT_STATE_REFERENCE_VERSION <- 1L

newArtifactStateReferenceHint <- function(
    method,
    normalized_parameters,
    software,
    lineage,
    compaction = list(
        enabled = FALSE,
        reason = "metadata_only_generation"
    )
) {
    hint <- list(
        schema = .ARTIFACT_STATE_REFERENCE_HINT_SCHEMA,
        schema_version = .ARTIFACT_STATE_REFERENCE_HINT_VERSION,
        representation = "state_reference",
        method = method,
        normalized_parameters = normalized_parameters,
        software = software,
        lineage = lineage,
        compaction = compaction
    )
    artifactWorkflowStateValidateStateReferenceHint(hint)
}

artifactWorkflowStateValidateStateReferenceHint <- function(hint) {
    required <- c(
        "schema", "schema_version", "representation", "method",
        "normalized_parameters", "software", "lineage", "compaction"
    )
    valid <- is.list(hint) && identical(names(hint), required) &&
        identical(hint$schema, .ARTIFACT_STATE_REFERENCE_HINT_SCHEMA) &&
        identical(
            as.integer(hint$schema_version),
            .ARTIFACT_STATE_REFERENCE_HINT_VERSION
        ) && identical(hint$representation, "state_reference") &&
        workflowStateScalarString(hint$method) &&
        is.list(hint$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(hint$software) &&
        artifactWorkflowStateValidSelectionLineage(hint$lineage) &&
        is.list(hint$compaction) && identical(hint$compaction$enabled, FALSE) &&
        workflowStateScalarString(hint$compaction$reason)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact state-reference hint is malformed",
            "multischolar_invalid_artifact_state_reference"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(hint, "state-reference hint")
    hint$schema_version <- as.integer(hint$schema_version)
    hint
}

artifactWorkflowStateWriteStateReferenceData <- function(
    previous_manifest,
    parent_object,
    state_object,
    hint
) {
    hint <- artifactWorkflowStateValidateStateReferenceHint(hint)
    if (!identical(parent_object, state_object)) {
        artifactWorkflowStateAbort(
            "artifact state reference requires an unchanged parent object",
            "multischolar_artifact_state_reference_requires_materialization"
        )
    }
    lineage <- hint$lineage
    lineage$parent_generation_id <- previous_manifest$generation_id
    lineage$parent_content_identity <- previous_manifest$data$semantic_digest
    metadata <- list(
        schema = .ARTIFACT_STATE_REFERENCE_SCHEMA,
        schema_version = .ARTIFACT_STATE_REFERENCE_VERSION,
        representation = "state_reference",
        parent_generation_id = previous_manifest$generation_id,
        parent_semantic_digest = previous_manifest$data$semantic_digest,
        parent_manifest_digest = previous_manifest$manifest_digest,
        method = hint$method,
        normalized_parameters = hint$normalized_parameters,
        software = hint$software,
        lineage = lineage,
        compaction = hint$compaction
    )
    list(
        s4_class = previous_manifest$s4_class,
        semantic_digest = previous_manifest$data$semantic_digest,
        codec = previous_manifest$data$codec,
        metadata_json = artifactWorkflowStateSerializeMetadata(
            metadata,
            "state-reference metadata"
        ),
        artifact_refs = list(),
        reused = TRUE,
        new_reference_names = character(),
        dependencies = list()
    )
}

artifactWorkflowStateValidateStateReference <- function(metadata, manifest) {
    required <- c(
        "schema", "schema_version", "representation",
        "parent_generation_id", "parent_semantic_digest",
        "parent_manifest_digest", "method", "normalized_parameters",
        "software", "lineage", "compaction"
    )
    valid <- is.list(metadata) && identical(names(metadata), required) &&
        identical(metadata$schema, .ARTIFACT_STATE_REFERENCE_SCHEMA) &&
        identical(
            as.integer(metadata$schema_version),
            .ARTIFACT_STATE_REFERENCE_VERSION
        ) && identical(metadata$representation, "state_reference") &&
        identical(metadata$parent_generation_id, manifest$parent_generation_id) &&
        identical(metadata$parent_semantic_digest, manifest$data$semantic_digest) &&
        workflowStateScalarString(metadata$method) &&
        is.list(metadata$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(metadata$software) &&
        artifactWorkflowStateValidSelectionLineage(metadata$lineage) &&
        is.list(metadata$compaction) &&
        identical(metadata$compaction$enabled, FALSE) &&
        workflowStateScalarString(metadata$compaction$reason) &&
        length(manifest$data$artifact_refs) == 0L
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact state-reference metadata is malformed",
            "multischolar_invalid_artifact_state_reference"
        )
    }
    artifactRefValidateDigest(
        metadata$parent_semantic_digest,
        "parent_semantic_digest"
    )
    artifactRefValidateDigest(
        metadata$parent_manifest_digest,
        "parent_manifest_digest"
    )
    artifactWorkflowStateAssertSafeMetadata(metadata, "state-reference metadata")
    metadata
}

artifactWorkflowStateHydrateStateReference <- function(
    store,
    manifest,
    metadata,
    hydrate_fn,
    visited_generation_ids
) {
    metadata <- artifactWorkflowStateValidateStateReference(metadata, manifest)
    parent_manifest <- artifactWorkflowStateReadManifest(
        store,
        artifactWorkflowStateManifestRelativePath(
            store,
            metadata$parent_generation_id
        )
    )
    valid_parent <- identical(
        parent_manifest$data$semantic_digest,
        metadata$parent_semantic_digest
    ) && identical(
        parent_manifest$manifest_digest,
        metadata$parent_manifest_digest
    ) && identical(parent_manifest$s4_class, manifest$s4_class)
    if (!isTRUE(valid_parent)) {
        artifactWorkflowStateAbort(
            "artifact state reference differs from its immutable parent",
            "multischolar_artifact_state_reference_parent_mismatch"
        )
    }
    artifactWorkflowStateHydrateData(
        store,
        parent_manifest,
        hydrate_fn,
        visited_generation_ids
    )
}

artifactWorkflowStateIsStateReference <- function(metadata) {
    is.list(metadata) &&
        identical(metadata$schema, .ARTIFACT_STATE_REFERENCE_SCHEMA)
}
