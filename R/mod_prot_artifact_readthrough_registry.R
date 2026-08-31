# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_READTHROUGH_VERSION <- 1L
.PROT_DIA_IMPORT_ROLES <- "canonical_data"
.PROT_DIA_DESIGN_ROLES <- c(
    "cleaned_data", "design_matrix", "contrasts", "args", "annotations",
    "sequences"
)

#' Define the DIA-NN read-through compatibility adapter
#' @return A descriptor-bound artifact read-through adapter.
#' @noRd
protDiaArtifactReadthroughAdapter <- function() {
    newArtifactStageReadthroughAdapter(
        descriptor = artifactDiaWorkflowDescriptor(),
        owner_label = "DIA-NN",
        default_omic_label = "proteomics",
        abort_fn = protDiaArtifactAbort,
        compatible_descriptor_contracts = list(
            artifactStageDescriptorContract(artifactDiaWorkflowDescriptorV1())
        ),
        conditions = c(
            corrupt_manifest = "multischolar_corrupt_prot_dia_artifact_manifest",
            invalid_context = "multischolar_invalid_prot_dia_artifact_context",
            descriptor_pin_missing =
                "multischolar_prot_dia_descriptor_pin_missing",
            descriptor_pin_mismatch =
                "multischolar_prot_dia_descriptor_pin_mismatch",
            registry_ref_mismatch =
                "multischolar_prot_dia_registry_ref_mismatch",
            missing_committed_stage =
                "multischolar_missing_prot_dia_committed_stage",
            invalid_run_parameters =
                "multischolar_invalid_prot_dia_run_parameters",
            missing_current_state =
                "multischolar_missing_prot_dia_current_state",
            incomplete_contract =
                "multischolar_incomplete_prot_dia_readthrough_contract",
            mixed_generation =
                "multischolar_mixed_prot_dia_artifact_generation",
            invalid_optional_payload =
                "multischolar_invalid_prot_dia_run_parameters",
            invalid_contrasts =
                "multischolar_inexact_prot_dia_artifact_hydration",
            invalid_column_mapping =
                "multischolar_invalid_prot_dia_run_parameters"
        )
    )
}

protDiaArtifactProbeIdentity <- function(omic_label) {
    artifactStageProbeIdentity(protDiaArtifactReadthroughAdapter(), omic_label)
}

protDiaArtifactProbeManifest <- function(experiment_paths, experiment_label) {
    artifactStageProbeManifest(
        protDiaArtifactReadthroughAdapter(),
        experiment_paths,
        experiment_label
    )
}

createProtDiaResumeContext <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    createArtifactStageResumeContext(
        protDiaArtifactReadthroughAdapter(),
        experiment_paths,
        experiment_label,
        storage_policy
    )
}

protDiaArtifactValidateDescriptorPin <- function(store) {
    artifactStageValidateDescriptorPin(protDiaArtifactReadthroughAdapter(), store)
}

protDiaArtifactReadResources <- function(context, resource_policy = NULL) {
    artifactStageReadResources(
        protDiaArtifactReadthroughAdapter(),
        context,
        resource_policy
    )
}

protDiaArtifactReadRegistryRef <- function(
    resources,
    row,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    artifactStageReadRegistryRef(
        protDiaArtifactReadthroughAdapter(),
        resources,
        row,
        payload_validation
    )
}

protDiaArtifactDataFrameRow <- function(value, index) {
    artifactStageDataFrameRow(value, index)
}

protDiaArtifactResolveRun <- function(
    resources,
    run,
    stage_id,
    required_roles,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    artifactStageResolveRun(
        protDiaArtifactReadthroughAdapter(),
        resources,
        run,
        stage_id,
        required_roles,
        payload_validation
    )
}

resolveProtDiaCommittedStage <- function(
    resources,
    stage_id,
    required_roles,
    run_id = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    resolveArtifactStageCommittedStage(
        protDiaArtifactReadthroughAdapter(),
        resources,
        stage_id,
        required_roles,
        run_id,
        payload_validation
    )
}

protDiaArtifactDecodeParameters <- function(resources, run_id, parameter_names) {
    artifactStageDecodeParameters(
        protDiaArtifactReadthroughAdapter(),
        resources,
        run_id,
        parameter_names
    )
}

protDiaArtifactCurrentStateRow <- function(resources) {
    artifactStageCurrentStateRow(protDiaArtifactReadthroughAdapter(), resources)
}

collectProtDiaResumeEvidence <- function(
    context,
    resource_policy = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    collectArtifactStageResumeEvidence(
        adapter = protDiaArtifactReadthroughAdapter(),
        context = context,
        import_roles = .PROT_DIA_IMPORT_ROLES,
        design_roles = .PROT_DIA_DESIGN_ROLES,
        import_parameter_names = c(
            "readthrough_contract_version", "column_mapping_serialized",
            "use_precursor_norm", "sanitize_names"
        ),
        design_parameter_names = c(
            "state_name", "workflow_type", "readthrough_contract_version",
            "parent_import_run_id", "parent_import_generation_id",
            "parent_import_artifact_id", "parent_import_semantic_digest",
            "contrasts_kind", "annotation_available", "sequence_data_available"
        ),
        resource_policy = resource_policy,
        payload_validation = payload_validation
    )
}
