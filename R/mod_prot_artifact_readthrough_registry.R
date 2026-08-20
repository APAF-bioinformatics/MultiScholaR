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

protDiaArtifactProbeIdentity <- function(omic_label) {
    identity <- artifactDiaWorkflowDescriptor()$identity
    identity$omic_label <- artifactValidatePathComponent(omic_label, "omic_label")
    identity
}

protDiaArtifactProbeManifest <- function(experiment_paths, experiment_label) {
    root <- workflowNormalizeProjectRoot(
        workflowProjectRootFromPaths(experiment_paths)
    )
    omic_label <- experiment_paths$omic_label
    if (is.null(omic_label)) omic_label <- experiment_label
    if (is.null(omic_label)) omic_label <- "proteomics"
    identity <- protDiaArtifactProbeIdentity(omic_label)
    relative <- artifactWorkflowRelativePaths(identity)$artifact_manifest
    path <- artifactResolveContainedPath(root, relative)
    if (!file.exists(path)) {
        return(list(found = FALSE, root = root, identity = identity))
    }
    manifest <- tryCatch(
        jsonlite::read_json(path, simplifyVector = TRUE),
        error = \(error) protDiaArtifactAbort(
            "DIA-NN artifact root manifest could not be parsed",
            "multischolar_corrupt_prot_dia_artifact_manifest",
            manifest_path = relative,
            parent = error
        )
    )
    required <- c(
        "schema", "schema_version", "project_id", "workflow_id", "omic_type",
        "workflow_slug", "created_at"
    )
    version <- workflowStateVersionValue(manifest$schema_version)
    valid <- is.list(manifest) && identical(names(manifest), required) &&
        identical(manifest$schema, .ARTIFACT_MANIFEST_SCHEMA) &&
        identical(version, .ARTIFACT_MANIFEST_SCHEMA_VERSION) &&
        workflowCapabilityScalarString(manifest$project_id) &&
        identical(manifest$workflow_id, identity$workflow_id) &&
        identical(manifest$omic_type, identity$omic_type) &&
        identical(manifest$workflow_slug, identity$workflow_slug) &&
        artifactRefValidUtc(manifest$created_at)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN artifact root manifest does not match the canary workflow",
            "multischolar_corrupt_prot_dia_artifact_manifest",
            manifest_path = relative
        )
    }
    list(
        found = TRUE,
        root = root,
        identity = identity,
        manifest = manifest,
        manifest_path = relative
    )
}

createProtDiaResumeContext <- function(
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) {
        return(list(enabled = FALSE, reason = "explicit_memory_backend"))
    }
    probe <- protDiaArtifactProbeManifest(experiment_paths, experiment_label)
    if (!isTRUE(probe$found)) {
        return(list(enabled = FALSE, reason = "artifact_manifest_absent"))
    }
    policy$project_id <- probe$manifest$project_id
    context <- createWorkflowContext(
        experiment_paths,
        omic_type = "proteomics",
        experiment_label = experiment_label,
        storage_policy = policy
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "DIA",
        input_format = "diann",
        data_level = "peptide",
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
    )
    decision <- context$getStorageDecision()
    enabled <- identical(decision$effective_backend, "artifact") &&
        identical(decision$effective_rollout, "dual_write") &&
        identical(context$getIdentity()$project_id, probe$manifest$project_id) &&
        protDiaArtifactIdentityMatches(context$getIdentity())
    if (!isTRUE(enabled)) {
        return(list(enabled = FALSE, reason = decision$reason_code))
    }
    list(
        enabled = TRUE,
        reason = "committed_artifact_project",
        context = context,
        manifest = probe$manifest
    )
}

protDiaArtifactValidateDescriptorPin <- function(store) {
    relative <- artifactWorkflowStateDescriptorPinPath(store)
    path <- artifactStoreResolveFile(store, relative)
    if (!file.exists(path)) {
        protDiaArtifactAbort(
            "DIA-NN read-through requires an immutable descriptor pin",
            "multischolar_prot_dia_descriptor_pin_missing"
        )
    }
    pin <- artifactStoreReadJson(store, relative)
    version <- workflowStateVersionValue(pin$schema_version)
    valid <- identical(pin$schema, ARTIFACT_DESCRIPTOR_PIN_SCHEMA) &&
        identical(version, ARTIFACT_DESCRIPTOR_PIN_VERSION) &&
        identical(pin$project_id, store$project_id) &&
        identical(pin$workflow_id, artifactDiaWorkflowDescriptor()$identity$workflow_id) &&
        identical(pin$contract, protDiaArtifactDescriptorContract())
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN read-through descriptor differs from the project pin",
            "multischolar_prot_dia_descriptor_pin_mismatch"
        )
    }
    invisible(pin)
}

protDiaArtifactReadResources <- function(context, resource_policy = NULL) {
    if (!inherits(context, "WorkflowContext") || !context$isBound() ||
        !protDiaArtifactCoordinatorOwned(list(workflow_context = context))) {
        protDiaArtifactAbort(
            "DIA-NN read-through requires its exact artifact context",
            "multischolar_invalid_prot_dia_artifact_context"
        )
    }
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    protDiaArtifactValidateDescriptorPin(store)
    registry <- projectRegistryForContext(context, resource_policy)
    session <- openProjectRegistryReadOnly(registry)
    list(identity = identity, store = store, registry = registry, session = session)
}

protDiaArtifactReadRegistryRef <- function(resources, row) {
    key <- list(
        project_id = resources$identity$project_id,
        omic_type = resources$identity$omic_type,
        workflow_slug = resources$identity$workflow_slug,
        stage_id = row$stage_id,
        state_role = row$state_role,
        generation_id = row$generation_id
    )
    sidecar_path <- artifactStoreManagedPaths(
        resources$store,
        key,
        row$artifact_id
    )$sidecar
    sidecar <- artifactStoreReadSidecar(
        resources$store,
        sidecar_path,
        validate_payload = TRUE
    )
    ref <- artifactStoreNormalizeRef(sidecar$artifact_ref)
    matches <- identical(ref$artifact_id, row$artifact_id) &&
        identical(ref$logical_key, key) &&
        identical(ref$relative_path, row$relative_path) &&
        identical(ref$codec$id, row$codec_id) &&
        identical(ref$payload_schema$id, row$payload_schema_id) &&
        identical(ref$hash_policy$semantic$digest, row$semantic_digest) &&
        identical(ref$hash_policy$byte$digest, row$byte_digest) &&
        identical(as.numeric(ref$shape$bytes), as.numeric(row$payload_bytes)) &&
        identical(ref$status, row$status) && row$run_id %in% ref$provenance_ids
    if (!isTRUE(matches)) {
        protDiaArtifactAbort(
            "DIA-NN registry artifact differs from its immutable sidecar",
            "multischolar_prot_dia_registry_ref_mismatch",
            artifact_id = row$artifact_id,
            run_id = row$run_id
        )
    }
    ref
}

protDiaArtifactDataFrameRow <- function(value, index) {
    row <- as.list(value[index, , drop = FALSE])
    lapply(row, `[[`, 1L)
}

protDiaArtifactResolveRun <- function(resources, run, stage_id, required_roles) {
    links <- projectRegistryQuery(
        resources$session,
        "run_artifacts",
        filters = list(
            workflow_id = resources$identity$workflow_id,
            run_id = run$run_id,
            direction = "output"
        )
    )
    complete <- nrow(links) == length(required_roles) &&
        setequal(links$artifact_role, required_roles) &&
        identical(sort(as.integer(links$ordinal)), seq_along(required_roles) - 1L)
    if (!isTRUE(complete)) return(NULL)
    refs <- lapply(required_roles, function(role) {
        link <- links[links$artifact_role == role, , drop = FALSE]
        if (nrow(link) != 1L) return(NULL)
        rows <- projectRegistryQuery(
            resources$session,
            "artifacts",
            filters = list(
                workflow_id = resources$identity$workflow_id,
                artifact_id = link$artifact_id[[1L]]
            ),
            limit = 1L
        )
        if (nrow(rows) != 1L) return(NULL)
        row <- protDiaArtifactDataFrameRow(rows, 1L)
        valid <- identical(row$run_id, run$run_id) &&
            identical(row$stage_id, stage_id) &&
            identical(row$state_role, role) && identical(row$status, "committed")
        if (!isTRUE(valid)) return(NULL)
        protDiaArtifactReadRegistryRef(resources, row)
    })
    if (any(vapply(refs, is.null, logical(1)))) return(NULL)
    generations <- unique(vapply(
        refs,
        \(ref) ref$logical_key$generation_id,
        character(1)
    ))
    if (length(generations) != 1L) return(NULL)
    names(refs) <- required_roles
    list(
        stage_id = stage_id,
        run_id = run$run_id,
        generation_id = generations[[1L]],
        refs = refs,
        completed_at = run$completed_at
    )
}

resolveProtDiaCommittedStage <- function(
    resources,
    stage_id,
    required_roles,
    run_id = NULL
) {
    filters <- list(
        workflow_id = resources$identity$workflow_id,
        status = "completed"
    )
    if (!is.null(run_id)) filters$run_id <- run_id
    runs <- projectRegistryQuery(resources$session, "runs", filters = filters)
    if (nrow(runs) > 0L) {
        for (index in rev(seq_len(nrow(runs)))) {
            run <- protDiaArtifactDataFrameRow(runs, index)
            resolved <- protDiaArtifactResolveRun(
                resources,
                run,
                stage_id,
                required_roles
            )
            if (!is.null(resolved)) return(resolved)
        }
    }
    protDiaArtifactAbort(
        sprintf("DIA-NN stage '%s' has no complete committed artifact set", stage_id),
        "multischolar_missing_prot_dia_committed_stage",
        project_id = resources$identity$project_id,
        workflow_id = resources$identity$workflow_id,
        stage_id = stage_id,
        requested_run_id = run_id,
        remediation = paste(
            "Complete a new dual-write import and design generation, or use",
            "one explicit whole-workflow legacy adapter."
        )
    )
}

protDiaArtifactDecodeParameters <- function(resources, run_id, parameter_names) {
    rows <- projectRegistryQuery(
        resources$session,
        "parameters",
        filters = list(
            workflow_id = resources$identity$workflow_id,
            run_id = run_id
        )
    )
    rows <- rows[rows$parameter_name %in% parameter_names, , drop = FALSE]
    if (nrow(rows) != length(parameter_names) ||
        anyDuplicated(rows$parameter_name) > 0L ||
        !setequal(rows$parameter_name, parameter_names)) {
        protDiaArtifactAbort(
            "DIA-NN artifact run parameters are missing or ambiguous",
            "multischolar_invalid_prot_dia_run_parameters",
            run_id = run_id
        )
    }
    values <- lapply(seq_len(nrow(rows)), function(index) {
        value <- jsonlite::fromJSON(
            rows$value_json[[index]],
            simplifyVector = FALSE
        )
        if (!identical(artifactSemanticDigest(value), rows$value_digest[[index]])) {
            protDiaArtifactAbort(
                "DIA-NN artifact run parameter digest is invalid",
                "multischolar_invalid_prot_dia_run_parameters",
                run_id = run_id,
                parameter_name = rows$parameter_name[[index]]
            )
        }
        value
    })
    stats::setNames(values, rows$parameter_name)
}

protDiaArtifactCurrentStateRow <- function(resources) {
    rows <- projectRegistryQuery(
        resources$session,
        "states",
        filters = list(
            workflow_id = resources$identity$workflow_id,
            status = "current"
        )
    )
    if (nrow(rows) != 1L || identical(rows$logical_name[[1L]], "initial")) {
        protDiaArtifactAbort(
            "DIA-NN read-through requires one committed current design state",
            "multischolar_missing_prot_dia_current_state"
        )
    }
    protDiaArtifactDataFrameRow(rows, 1L)
}

collectProtDiaResumeEvidence <- function(context, resource_policy = NULL) {
    resources <- protDiaArtifactReadResources(context, resource_policy)
    on.exit(closeProjectRegistry(resources$session), add = TRUE)
    design <- resolveProtDiaCommittedStage(
        resources,
        "design",
        .PROT_DIA_DESIGN_ROLES
    )
    design_parameters <- c(
        "state_name", "workflow_type", "readthrough_contract_version",
        "parent_import_run_id", "parent_import_generation_id",
        "parent_import_artifact_id", "parent_import_semantic_digest",
        "contrasts_kind", "annotation_available", "sequence_data_available"
    )
    design$parameters <- protDiaArtifactDecodeParameters(
        resources,
        design$run_id,
        design_parameters
    )
    parent_run_id <- design$parameters$parent_import_run_id
    if (!workflowCapabilityScalarString(parent_run_id)) {
        protDiaArtifactAbort(
            "DIA-NN design artifacts do not declare one parent import run",
            "multischolar_incomplete_prot_dia_readthrough_contract",
            run_id = design$run_id
        )
    }
    import <- resolveProtDiaCommittedStage(
        resources,
        "import",
        .PROT_DIA_IMPORT_ROLES,
        run_id = parent_run_id
    )
    import$parameters <- protDiaArtifactDecodeParameters(
        resources,
        import$run_id,
        c(
            "readthrough_contract_version", "column_mapping_serialized",
            "use_precursor_norm", "sanitize_names"
        )
    )
    list(
        identity = resources$identity,
        store = resources$store,
        import = import,
        design = design,
        current_state = protDiaArtifactCurrentStateRow(resources)
    )
}
