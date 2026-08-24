# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Probe exact proteomics artifact resume candidates
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param storage_policy Optional workflow storage policy.
#' @return A list containing zero or one exact resume candidate.
#' @noRd
probeProtArtifactResumeCandidates <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    if (!artifactStageReadthroughEnabled(storage_policy)) {
        return(list(reason = "readthrough_disabled", candidates = list()))
    }
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) {
        return(list(reason = "explicit_memory_backend", candidates = list()))
    }
    dia <- protDiaArtifactProbeManifest(experiment_paths, experiment_label)
    nondia <- probeProtNonDiaArtifactManifests(
        experiment_paths,
        experiment_label
    )$found
    candidates <- list()
    if (isTRUE(dia$found)) {
        candidates[[artifactDiaWorkflowDescriptor()$descriptor_id]] <- list(
            kind = "dia",
            probe = dia
        )
    }
    if (length(nondia) == 1L) {
        capability_id <- names(nondia)[[1L]]
        candidates[[capability_id]] <- list(
            kind = "nondia",
            probe = nondia[[1L]]
        )
    }
    if (length(candidates) > 1L) {
        rlang::abort(
            "multiple proteomics artifact manifests make resume ambiguous",
            class = c(
                "multischolar_ambiguous_prot_artifact_project",
                "multischolar_prot_artifact_resume_error"
            ),
            capability_ids = names(candidates)
        )
    }
    list(
        reason = if (length(candidates) == 0L) {
            "artifact_manifest_absent"
        } else {
            "artifact_manifest_found"
        },
        candidates = candidates
    )
}

#' Safely dispatch exact proteomics artifact resume
#' @param workflow_data Mutable proteomics workflow state.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param storage_policy Optional workflow storage policy.
#' @param resource_policy Optional project registry resource policy.
#' @param log_warn Warning logger used for additive resume failures.
#' @return The recorded exact-tuple resume result, invisibly.
#' @noRd
resumeProtArtifactWorkflowSafely <- function(
    workflow_data,
    experiment_paths,
    experiment_label,
    storage_policy = NULL,
    resource_policy = NULL,
    log_warn = logger::log_warn
) {
    candidates <- tryCatch(
        probeProtArtifactResumeCandidates(
            experiment_paths,
            experiment_label,
            storage_policy
        ),
        error = \(error) error
    )
    if (inherits(candidates, "error")) {
        log_warn(paste(
            "proteomics artifact resume probe failed without changing workflow state:",
            conditionMessage(candidates)
        ))
        result <- list(
            enabled = TRUE,
            ok = FALSE,
            resumed = FALSE,
            reason = "artifact_resume_probe_failed",
            artifact_project = TRUE,
            error_class = class(candidates),
            error_message = conditionMessage(candidates)
        )
        return(recordArtifactStageResult(workflow_data, "resume", result))
    }
    if (length(candidates$candidates) == 0L) {
        result <- list(
            enabled = FALSE,
            ok = TRUE,
            resumed = FALSE,
            reason = candidates$reason,
            artifact_project = FALSE
        )
        return(recordArtifactStageResult(workflow_data, "resume", result))
    }
    capability_id <- names(candidates$candidates)[[1L]]
    candidate <- candidates$candidates[[1L]]
    if (identical(candidate$kind, "dia")) {
        result <- resumeProtDiaArtifactWorkflowSafely(
            workflow_data,
            experiment_paths,
            experiment_label,
            storage_policy,
            resource_policy,
            log_warn
        )
    } else {
        result <- resumeProtNonDiaArtifactWorkflowSafely(
            workflow_data,
            experiment_paths,
            experiment_label,
            capability_id,
            storage_policy,
            resource_policy,
            log_warn = log_warn
        )
    }
    result$capability_id <- capability_id
    recordArtifactStageResult(workflow_data, "resume", result)
}
