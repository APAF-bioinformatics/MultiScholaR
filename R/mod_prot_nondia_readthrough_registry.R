# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_NONDIA_READTHROUGH_VERSION <- 1L
.PROT_NONDIA_IMPORT_ROLES <- "canonical_data"
.PROT_NONDIA_DESIGN_ROLES <- c(
    "cleaned_data", "design_matrix", "contrasts", "args", "annotations",
    "sequences"
)

#' Define supported non-DIA read-through descriptors
#' @return A named list containing the three exact certified descriptors.
#' @noRd
protNonDiaReadthroughDescriptors <- function() {
    ids <- names(artifactProteomicsNonDiaDescriptorSpecs())
    catalogue <- artifactWorkflowDescriptorCatalogue()
    descriptors <- artifactDescriptorCatalogueValues(catalogue)
    descriptors <- Filter(
        \(descriptor) descriptor$descriptor_id %in% ids,
        descriptors
    )
    stats::setNames(
        descriptors,
        vapply(descriptors, `[[`, character(1), "descriptor_id")
    )[ids]
}

#' Resolve one supported non-DIA read-through descriptor
#' @param capability_id Exact workflow capability identifier.
#' @return A validated workflow descriptor.
#' @noRd
protNonDiaReadthroughDescriptor <- function(capability_id) {
    if (!workflowCapabilityScalarString(capability_id)) {
        protNonDiaArtifactAbort(
            "non-DIA read-through requires one exact capability identifier",
            "multischolar_invalid_prot_nondia_resume_descriptor"
        )
    }
    descriptor <- protNonDiaReadthroughDescriptors()[[capability_id]]
    if (is.null(descriptor)) {
        protNonDiaArtifactAbort(
            "non-DIA read-through capability is not certified",
            "multischolar_invalid_prot_nondia_resume_descriptor",
            capability_id = capability_id
        )
    }
    descriptor
}

#' Define an exact non-DIA proteomics read-through adapter
#' @param descriptor Exact supported workflow descriptor.
#' @return A descriptor-bound artifact read-through adapter.
#' @noRd
protNonDiaArtifactReadthroughAdapter <- function(descriptor) {
    descriptor <- protNonDiaReadthroughDescriptor(descriptor$descriptor_id)
    newArtifactStageReadthroughAdapter(
        descriptor = descriptor,
        owner_label = "non-DIA proteomics",
        default_omic_label = "proteomics",
        abort_fn = protNonDiaArtifactAbort,
        conditions = c(
            corrupt_manifest =
                "multischolar_corrupt_prot_nondia_artifact_manifest",
            invalid_context =
                "multischolar_invalid_prot_nondia_artifact_context",
            descriptor_pin_missing =
                "multischolar_prot_nondia_descriptor_pin_missing",
            descriptor_pin_mismatch =
                "multischolar_prot_nondia_descriptor_pin_mismatch",
            registry_ref_mismatch =
                "multischolar_prot_nondia_registry_ref_mismatch",
            missing_committed_stage =
                "multischolar_missing_prot_nondia_committed_stage",
            invalid_run_parameters =
                "multischolar_invalid_prot_nondia_run_parameters",
            missing_current_state =
                "multischolar_missing_prot_nondia_current_state",
            incomplete_contract =
                "multischolar_incomplete_prot_nondia_readthrough_contract",
            mixed_generation =
                "multischolar_mixed_prot_nondia_artifact_generation",
            invalid_optional_payload =
                "multischolar_invalid_prot_nondia_run_parameters",
            invalid_contrasts =
                "multischolar_inexact_prot_nondia_artifact_hydration",
            invalid_column_mapping =
                "multischolar_invalid_prot_nondia_run_parameters"
        )
    )
}

#' Probe all certified non-DIA proteomics manifest paths
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @return A list of manifest probes, including at most one found manifest.
#' @noRd
probeProtNonDiaArtifactManifests <- function(
    experiment_paths,
    experiment_label
) {
    descriptors <- protNonDiaReadthroughDescriptors()
    probes <- lapply(descriptors, \(descriptor) {
        artifactStageProbeManifest(
            protNonDiaArtifactReadthroughAdapter(descriptor),
            experiment_paths,
            experiment_label
        )
    })
    names(probes) <- names(descriptors)
    found <- probes[vapply(probes, \(probe) isTRUE(probe$found), logical(1))]
    if (length(found) > 1L) {
        protNonDiaArtifactAbort(
            "multiple non-DIA proteomics artifact manifests are ambiguous",
            "multischolar_ambiguous_prot_nondia_artifact_project",
            capability_ids = names(found)
        )
    }
    list(probes = probes, found = found)
}

#' Probe one exact non-DIA proteomics artifact manifest
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param capability_id Exact workflow capability identifier.
#' @return An exact manifest probe result.
#' @noRd
protNonDiaArtifactProbeManifest <- function(
    experiment_paths,
    experiment_label,
    capability_id
) {
    descriptor <- protNonDiaReadthroughDescriptor(capability_id)
    artifactStageProbeManifest(
        protNonDiaArtifactReadthroughAdapter(descriptor),
        experiment_paths,
        experiment_label
    )
}

#' Create an exact non-DIA proteomics resume context
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param capability_id Optional exact workflow capability identifier.
#' @param storage_policy Optional workflow storage policy.
#' @return A prepared resume context with an explicit enabled flag and reason.
#' @noRd
createProtNonDiaResumeContext <- function(
    experiment_paths,
    experiment_label,
    capability_id = NULL,
    storage_policy = NULL
) {
    if (!artifactStageReadthroughEnabled(storage_policy)) {
        return(list(enabled = FALSE, reason = "readthrough_disabled"))
    }
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) {
        return(list(enabled = FALSE, reason = "explicit_memory_backend"))
    }
    if (is.null(capability_id)) {
        found <- probeProtNonDiaArtifactManifests(
            experiment_paths,
            experiment_label
        )$found
        if (length(found) == 0L) {
            return(list(enabled = FALSE, reason = "artifact_manifest_absent"))
        }
        capability_id <- names(found)[[1L]]
    }
    descriptor <- protNonDiaReadthroughDescriptor(capability_id)
    createArtifactStageResumeContext(
        protNonDiaArtifactReadthroughAdapter(descriptor),
        experiment_paths,
        experiment_label,
        storage_policy
    )
}

#' Test whether one certified non-DIA artifact project is present
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param storage_policy Optional workflow storage policy.
#' @return A scalar logical.
#' @noRd
protNonDiaArtifactProjectPresent <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    if (!artifactStageReadthroughEnabled(storage_policy)) return(FALSE)
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) return(FALSE)
    tryCatch(
        length(probeProtNonDiaArtifactManifests(
            experiment_paths,
            experiment_label
        )$found) == 1L,
        error = \(...) TRUE
    )
}

#' Collect exact non-DIA proteomics resume evidence
#' @param context Exact bound workflow context.
#' @param descriptor Exact supported workflow descriptor.
#' @param resource_policy Optional project registry resource policy.
#' @param payload_validation Payload validation mode.
#' @return Immutable parent-linked import, design, and state evidence.
#' @noRd
collectProtNonDiaResumeEvidence <- function(
    context,
    descriptor,
    resource_policy = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    collectArtifactStageResumeEvidence(
        adapter = protNonDiaArtifactReadthroughAdapter(descriptor),
        context = context,
        import_roles = .PROT_NONDIA_IMPORT_ROLES,
        design_roles = .PROT_NONDIA_DESIGN_ROLES,
        import_parameter_names = c(
            "capability_id", "column_mapping", "column_mapping_serialized",
            "readthrough_contract_version", "workflow_type", "input_format",
            "data_level", "parser_parameters", "sanitize_names"
        ),
        design_parameter_names = c(
            "capability_id", "state_name", "workflow_type",
            "readthrough_contract_version", "parent_import_run_id",
            "parent_import_generation_id", "parent_import_artifact_id",
            "parent_import_semantic_digest", "contrasts_kind",
            "annotation_available", "sequence_data_available"
        ),
        resource_policy = resource_policy,
        payload_validation = payload_validation
    )
}
