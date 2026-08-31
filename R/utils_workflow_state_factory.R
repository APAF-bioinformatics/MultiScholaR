# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

newWorkflowState <- function(
    audit_enabled = TRUE,
    workflow_context = NULL,
    resource_policy = NULL,
    dehydrate_fn = dehydrateDiaS4Artifact,
    validate_bundle_fn = validateDiaS4Bundle,
    hydrate_fn = hydrateDiaS4Artifact,
    verify_hydration_fn = artifactWorkflowStateVerifyHydrationInline,
    workflow_descriptor = NULL,
    descriptor_catalogue = NULL,
    codec_catalogue = NULL,
    settled_bootstrap = NULL
) {
    memory <- is.null(workflow_context) ||
        !inherits(workflow_context, "WorkflowContext") ||
        !workflow_context$isBound() ||
        identical(
            workflow_context$getStorageDecision()$effective_backend,
            "memory"
        )
    if (memory) return(WorkflowState$new(audit_enabled = audit_enabled))
    descriptor_contract <- NULL
    if (!is.null(descriptor_catalogue)) {
        candidate <- findArtifactWorkflowDescriptor(
            workflow_context$getIdentity(),
            descriptor_catalogue
        )
        if (is.null(candidate)) {
            artifactWorkflowStateAbort(
                "artifact context has no exact descriptor in the injected catalogue",
                "multischolar_artifact_state_descriptor_mismatch"
            )
        }
        if (!is.null(workflow_descriptor) && !identical(candidate, workflow_descriptor)) {
            artifactWorkflowStateAbort(
                "injected workflow descriptor disagrees with its immutable catalogue",
                "multischolar_artifact_state_descriptor_mismatch"
            )
        }
        workflow_descriptor <- candidate
    }
    if (!is.null(workflow_descriptor)) {
        workflow_descriptor <- validateArtifactWorkflowDescriptor(workflow_descriptor)
        identity <- workflow_context$getIdentity()
        decision <- workflow_context$getStorageDecision()
        matching_identity <- identical(
            workflowCapabilityKey(identity),
            workflowCapabilityKey(workflow_descriptor$identity)
        )
        matching_capability <- identical(
            decision$capability_id,
            workflow_descriptor$descriptor_id
        ) && identical(
            decision$capability_version,
            workflow_descriptor$descriptor_version
        )
        if (!isTRUE(matching_identity) || !isTRUE(matching_capability)) {
            artifactWorkflowStateAbort(
                "artifact WorkflowState descriptor does not match its context decision",
                "multischolar_artifact_state_descriptor_mismatch"
            )
        }
        if (is.null(codec_catalogue)) codec_catalogue <- artifactS4CodecCatalogue()
        adapter <- artifactCodecAdapter(workflow_descriptor, codec_catalogue)
        dehydrate_fn <- adapter$dehydrate
        validate_bundle_fn <- adapter$validate
        hydrate_fn <- adapter$hydrate
        descriptor_contract <- list(
            descriptor_id = workflow_descriptor$descriptor_id,
            descriptor_version = workflow_descriptor$descriptor_version,
            descriptor_digest = workflow_descriptor$descriptor_digest
        )
    } else if (!is.null(codec_catalogue)) {
        artifactWorkflowStateAbort(
            "artifact codec catalogue requires an exact workflow descriptor",
            "multischolar_artifact_state_descriptor_mismatch"
        )
    }
    ArtifactWorkflowState$new(
        workflow_context = workflow_context,
        audit_enabled = audit_enabled,
        resource_policy = resource_policy,
        dehydrate_fn = dehydrate_fn,
        validate_bundle_fn = validate_bundle_fn,
        hydrate_fn = hydrate_fn,
        verify_hydration_fn = verify_hydration_fn,
        descriptor_contract = descriptor_contract,
        settled_bootstrap = settled_bootstrap
    )
}
