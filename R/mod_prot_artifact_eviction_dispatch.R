# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Resolve the exact descriptor bound to a proteomics workflow
#' @param workflow_data Mutable proteomics workflow state.
#' @return The exact bound descriptor, or `NULL` for memory/unbound workflows.
#' @noRd
protArtifactBoundDescriptor <- function(workflow_data) {
    context <- tryCatch(workflow_data$workflow_context, error = \(...) NULL)
    if (!inherits(context, "WorkflowContext") || !context$isBound()) return(NULL)
    findArtifactWorkflowDescriptor(
        context$getIdentity(),
        artifactWorkflowDescriptorCatalogue()
    )
}

#' Safely dispatch current-session proteomics artifact settlement
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @param rollout_fn Effective rollout provider.
#' @param log_warn Warning logger used for additive settlement failures.
#' @return The recorded exact-tuple settlement result, invisibly.
#' @noRd
settleProtArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn
) {
    descriptor <- protArtifactBoundDescriptor(workflow_data)
    if (is.null(descriptor)) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            evicted = FALSE,
            reason = "artifact_descriptor_unavailable"
        ))
    }
    if (identical(
        descriptor$descriptor_id,
        artifactDiaWorkflowDescriptor()$descriptor_id
    )) {
        return(settleProtDiaArtifactWorkflowSafely(
            workflow_data,
            experiment_paths,
            experiment_label,
            storage_policy,
            resource_policy,
            rollout_fn,
            log_warn
        ))
    }
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    settleProtNonDiaArtifactWorkflowSafely(
        workflow_data,
        experiment_paths,
        experiment_label,
        descriptor,
        storage_policy,
        resource_policy,
        rollout_fn,
        log_warn
    )
}

#' Safely dispatch resumed proteomics payload eviction
#' @param workflow_data Mutable proteomics workflow state.
#' @param rollout_fn Effective rollout provider.
#' @param log_warn Warning logger used for additive eviction failures.
#' @param ... Additional exact-tuple eviction controls used by tests.
#' @return The recorded exact-tuple eviction result, invisibly.
#' @noRd
evictProtArtifactWorkflowPayloadsSafely <- function(
    workflow_data,
    rollout_fn = \(context) context$getStorageDecision()$effective_rollout,
    log_warn = logger::log_warn,
    ...
) {
    descriptor <- protArtifactBoundDescriptor(workflow_data)
    if (is.null(descriptor)) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            evicted = FALSE,
            reason = "artifact_descriptor_unavailable"
        ))
    }
    if (identical(
        descriptor$descriptor_id,
        artifactDiaWorkflowDescriptor()$descriptor_id
    )) {
        return(evictProtDiaWorkflowPayloadsSafely(
            workflow_data,
            rollout_fn,
            log_warn,
            ...
        ))
    }
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    evictProtNonDiaWorkflowPayloadsSafely(
        workflow_data,
        descriptor,
        rollout_fn,
        log_warn,
        ...
    )
}
