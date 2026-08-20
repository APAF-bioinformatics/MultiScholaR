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
    hydrate_fn = hydrateDiaS4Artifact
) {
    memory <- is.null(workflow_context) ||
        !inherits(workflow_context, "WorkflowContext") ||
        !workflow_context$isBound() ||
        identical(
            workflow_context$getStorageDecision()$effective_backend,
            "memory"
        )
    if (memory) return(WorkflowState$new(audit_enabled = audit_enabled))
    ArtifactWorkflowState$new(
        workflow_context = workflow_context,
        audit_enabled = audit_enabled,
        resource_policy = resource_policy,
        dehydrate_fn = dehydrate_fn,
        validate_bundle_fn = validate_bundle_fn,
        hydrate_fn = hydrate_fn
    )
}
