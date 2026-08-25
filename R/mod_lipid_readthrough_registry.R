# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.LIPID_READTHROUGH_VERSION <- 1L
.LIPID_IMPORT_FIXED_ROLES <- c(
    "assay_manifest", "column_mapping", "source_manifest"
)
.LIPID_DESIGN_FIXED_ROLES <- c(
    "design_matrix", "contrasts", "args", "column_mapping",
    "assay_alignment", "raw_s4_dependencies"
)
.LIPID_IMPORT_PARAMETER_NAMES <- c(
    "capability_id", "assay_order", "assay_roles", "assay_identities",
    "column_mapping", "column_mapping_serialized",
    "readthrough_contract_version", "input_format", "data_level",
    "sample_names_sanitized", "evidence_boundary"
)
.LIPID_DESIGN_PARAMETER_NAMES <- c(
    "capability_id", "state_name", "workflow_type", "formula_string",
    "contrasts_kind", "assay_order", "parent_import_run_id",
    "parent_import_generation_id", "readthrough_contract_version"
)

lipidArtifactReadthroughAdapter <- function() {
    newArtifactStageReadthroughAdapter(
        descriptor = artifactLipidomicsWorkflowDescriptor(),
        owner_label = "lipidomics",
        default_omic_label = "lipidomics",
        abort_fn = lipidArtifactAbort,
        conditions = c(
            corrupt_manifest = "multischolar_corrupt_lipidomics_manifest",
            invalid_context = "multischolar_invalid_lipidomics_resume_context",
            descriptor_pin_missing = "multischolar_lipidomics_descriptor_pin_missing",
            descriptor_pin_mismatch = "multischolar_lipidomics_descriptor_pin_mismatch",
            registry_ref_mismatch = "multischolar_lipidomics_registry_ref_mismatch",
            missing_committed_stage = "multischolar_missing_lipidomics_stage",
            invalid_run_parameters = "multischolar_invalid_lipidomics_run_parameters",
            missing_current_state = "multischolar_missing_lipidomics_state",
            incomplete_contract = "multischolar_incomplete_lipidomics_readthrough",
            mixed_generation = "multischolar_mixed_lipidomics_generation",
            invalid_optional_payload = "multischolar_invalid_lipidomics_payload",
            invalid_contrasts = "multischolar_invalid_lipidomics_contrasts",
            invalid_column_mapping = "multischolar_invalid_lipidomics_mapping"
        )
    )
}

lipidArtifactProbeManifest <- function(experiment_paths, experiment_label) {
    artifactStageProbeManifest(
        lipidArtifactReadthroughAdapter(),
        experiment_paths,
        experiment_label
    )
}

createLipidResumeContext <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    createArtifactStageResumeContext(
        lipidArtifactReadthroughAdapter(),
        experiment_paths,
        experiment_label,
        storage_policy
    )
}

lipidArtifactProjectPresent <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    if (!artifactStageReadthroughEnabled(storage_policy)) return(FALSE)
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) return(FALSE)
    tryCatch(
        isTRUE(lipidArtifactProbeManifest(
            experiment_paths,
            experiment_label
        )$found),
        error = \(...) TRUE
    )
}

lipidArtifactRunHasParameters <- function(
    resources,
    run_id,
    parameter_names
) {
    rows <- projectRegistryQuery(
        resources$session,
        "parameters",
        filters = list(
            workflow_id = resources$registry_identity$workflow_id,
            run_id = run_id
        )
    )
    rows <- rows[rows$parameter_name %in% parameter_names, , drop = FALSE]
    nrow(rows) == length(parameter_names) &&
        !anyDuplicated(rows$parameter_name) &&
        setequal(rows$parameter_name, parameter_names)
}

lipidArtifactStageRoles <- function(stage_id, parameters) {
    assay_order <- unlist(parameters$assay_order, use.names = FALSE)
    if (stage_id == "import") {
        assay_roles <- unlist(parameters$assay_roles, use.names = FALSE)
        valid <- length(assay_order) > 0L &&
            length(assay_roles) == length(assay_order) &&
            !anyDuplicated(assay_order) && !anyDuplicated(assay_roles)
        if (!isTRUE(valid)) return(NULL)
        return(c(assay_roles, .LIPID_IMPORT_FIXED_ROLES))
    }
    if (length(assay_order) == 0L || anyDuplicated(assay_order)) return(NULL)
    c(
        sprintf("cleaned_assay_%04d", seq_along(assay_order)),
        .LIPID_DESIGN_FIXED_ROLES
    )
}

resolveLipidCommittedStage <- function(
    resources,
    stage_id,
    parameter_names,
    run_id = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    adapter <- lipidArtifactReadthroughAdapter()
    payload_validation <- match.arg(payload_validation)
    filters <- list(
        workflow_id = resources$registry_identity$workflow_id,
        status = "completed"
    )
    if (!is.null(run_id)) filters$run_id <- run_id
    runs <- projectRegistryQuery(resources$session, "runs", filters = filters)
    if (nrow(runs) > 0L) {
        for (index in rev(seq_len(nrow(runs)))) {
            run <- artifactStageDataFrameRow(runs, index)
            if (!lipidArtifactRunHasParameters(
                resources,
                run$run_id,
                parameter_names
            )) next
            parameters <- artifactStageDecodeParameters(
                adapter,
                resources,
                run$run_id,
                parameter_names
            )
            roles <- lipidArtifactStageRoles(stage_id, parameters)
            if (is.null(roles)) next
            stage <- artifactStageResolveRun(
                adapter,
                resources,
                run,
                stage_id,
                roles,
                payload_validation
            )
            if (!is.null(stage)) {
                stage$parameters <- parameters
                return(stage)
            }
        }
    }
    artifactStageReadthroughAbort(
        adapter,
        "missing_committed_stage",
        sprintf("lipidomics stage '%s' has no complete artifact set", stage_id),
        stage_id = stage_id,
        requested_run_id = run_id
    )
}

collectLipidResumeEvidence <- function(
    context,
    resource_policy = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    adapter <- lipidArtifactReadthroughAdapter()
    payload_validation <- match.arg(payload_validation)
    resources <- artifactStageReadResources(adapter, context, resource_policy)
    on.exit(closeProjectRegistry(resources$session), add = TRUE)
    design <- resolveLipidCommittedStage(
        resources,
        "design",
        .LIPID_DESIGN_PARAMETER_NAMES,
        payload_validation = payload_validation
    )
    parent_run_id <- design$parameters$parent_import_run_id
    if (!workflowCapabilityScalarString(parent_run_id)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "lipidomics design does not declare one parent import run"
        )
    }
    import <- resolveLipidCommittedStage(
        resources,
        "import",
        .LIPID_IMPORT_PARAMETER_NAMES,
        run_id = parent_run_id,
        payload_validation = payload_validation
    )
    list(
        identity = resources$identity,
        descriptor = adapter$descriptor,
        store = resources$store,
        import = import,
        design = design,
        current_state = artifactStageCurrentStateRow(adapter, resources)
    )
}
