# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_STAGE_READTHROUGH_CONDITIONS <- c(
    "corrupt_manifest",
    "invalid_context",
    "descriptor_pin_missing",
    "descriptor_pin_mismatch",
    "registry_ref_mismatch",
    "missing_committed_stage",
    "invalid_run_parameters",
    "missing_current_state",
    "incomplete_contract",
    "mixed_generation",
    "invalid_optional_payload",
    "invalid_contrasts",
    "invalid_column_mapping"
)

#' Define a descriptor-bound artifact read-through adapter
#' @param descriptor Exact immutable workflow descriptor.
#' @param owner_label Human-readable workflow owner.
#' @param default_omic_label Default label used when paths contain no label.
#' @param abort_fn Typed abort function supplied by the workflow adapter.
#' @param conditions Named condition classes for shared failure categories.
#' @return A validated artifact read-through adapter.
#' @noRd
newArtifactStageReadthroughAdapter <- function(
    descriptor,
    owner_label,
    default_omic_label,
    abort_fn,
    conditions
) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    valid <- workflowCapabilityScalarString(owner_label) &&
        workflowCapabilityScalarString(default_omic_label) &&
        is.function(abort_fn) && is.character(conditions) &&
        identical(names(conditions), .ARTIFACT_STAGE_READTHROUGH_CONDITIONS) &&
        all(nzchar(conditions))
    if (!isTRUE(valid)) {
        artifactStagePersistenceAbort(
            "artifact read-through adapter is malformed",
            "multischolar_invalid_artifact_stage_readthrough_adapter"
        )
    }
    structure(
        list(
            descriptor = descriptor,
            owner_label = owner_label,
            default_omic_label = default_omic_label,
            abort_fn = abort_fn,
            conditions = conditions
        ),
        class = "ArtifactStageReadthroughAdapter"
    )
}

#' Validate a descriptor-bound artifact read-through adapter
#' @param adapter Candidate read-through adapter.
#' @return The validated adapter.
#' @noRd
validateArtifactStageReadthroughAdapter <- function(adapter) {
    valid <- inherits(adapter, "ArtifactStageReadthroughAdapter") &&
        is.list(adapter) && identical(
            names(adapter),
            c(
                "descriptor", "owner_label", "default_omic_label",
                "abort_fn", "conditions"
            )
        )
    if (!isTRUE(valid)) {
        artifactStagePersistenceAbort(
            "artifact read-through adapter is invalid",
            "multischolar_invalid_artifact_stage_readthrough_adapter"
        )
    }
    validateArtifactWorkflowDescriptor(adapter$descriptor)
    if (!is.function(adapter$abort_fn) ||
        !identical(
            names(adapter$conditions),
            .ARTIFACT_STAGE_READTHROUGH_CONDITIONS
        )) {
        artifactStagePersistenceAbort(
            "artifact read-through adapter contract changed",
            "multischolar_invalid_artifact_stage_readthrough_adapter"
        )
    }
    adapter
}

#' Abort through a workflow-specific read-through adapter
#' @param adapter Validated read-through adapter.
#' @param condition Shared failure category.
#' @param message Human-readable failure message.
#' @param ... Additional condition fields.
#' @return This function does not return; it signals a typed error.
#' @noRd
artifactStageReadthroughAbort <- function(adapter, condition, message, ...) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    class <- adapter$conditions[[condition]]
    if (!workflowCapabilityScalarString(class)) {
        artifactStagePersistenceAbort(
            "artifact read-through failure category is unknown",
            "multischolar_invalid_artifact_stage_readthrough_adapter",
            condition = condition
        )
    }
    adapter$abort_fn(message, class, ...)
}

#' Resolve the project identity used to probe an artifact manifest
#' @param adapter Validated read-through adapter.
#' @param omic_label Project omic label.
#' @return The descriptor identity with the project-specific omic label.
#' @noRd
artifactStageProbeIdentity <- function(adapter, omic_label) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    identity <- adapter$descriptor$identity
    identity$omic_label <- artifactValidatePathComponent(
        omic_label,
        "omic_label"
    )
    identity
}

#' Probe one exact descriptor manifest without creating resources
#' @param adapter Validated read-through adapter.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @return A manifest probe result with a stable project root and identity.
#' @noRd
artifactStageProbeManifest <- function(
    adapter,
    experiment_paths,
    experiment_label
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    root <- workflowNormalizeProjectRoot(
        workflowProjectRootFromPaths(experiment_paths)
    )
    omic_label <- experiment_paths$omic_label
    if (is.null(omic_label)) omic_label <- experiment_label
    if (is.null(omic_label)) omic_label <- adapter$default_omic_label
    identity <- artifactStageProbeIdentity(adapter, omic_label)
    relative <- artifactWorkflowRelativePaths(identity)$artifact_manifest
    path <- artifactResolveContainedPath(root, relative)
    if (!file.exists(path)) {
        return(list(found = FALSE, root = root, identity = identity))
    }
    manifest <- tryCatch(
        jsonlite::read_json(path, simplifyVector = TRUE),
        error = \(error) artifactStageReadthroughAbort(
            adapter,
            "corrupt_manifest",
            sprintf("%s artifact root manifest could not be parsed", adapter$owner_label),
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
        artifactStageReadthroughAbort(
            adapter,
            "corrupt_manifest",
            sprintf(
                "%s artifact root manifest does not match its workflow",
                adapter$owner_label
            ),
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

#' Test whether artifact read-through is enabled by session policy
#' @param storage_policy Optional workflow storage policy.
#' @return A scalar logical.
#' @noRd
artifactStageReadthroughEnabled <- function(storage_policy = NULL) {
    if (is.null(storage_policy)) return(TRUE)
    if (!is.list(storage_policy)) {
        workflowCapabilityAbort(
            "workflow storage policy must be an R list",
            "multischolar_invalid_storage_policy"
        )
    }
    enabled <- storage_policy$readthrough_enabled
    if (is.null(enabled)) return(TRUE)
    if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
        workflowCapabilityAbort(
            "artifact read-through policy must be true or false",
            "multischolar_invalid_storage_policy"
        )
    }
    enabled
}

#' Create an exact context for artifact read-through
#' @param adapter Validated read-through adapter.
#' @param experiment_paths Workflow project paths.
#' @param experiment_label Project experiment label.
#' @param storage_policy Optional workflow storage policy.
#' @return A prepared context result with an explicit enabled flag and reason.
#' @noRd
createArtifactStageResumeContext <- function(
    adapter,
    experiment_paths,
    experiment_label,
    storage_policy = NULL
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    if (!artifactStageReadthroughEnabled(storage_policy)) {
        return(list(enabled = FALSE, reason = "readthrough_disabled"))
    }
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    if (identical(policy$requested_backend, "memory")) {
        return(list(enabled = FALSE, reason = "explicit_memory_backend"))
    }
    probe <- artifactStageProbeManifest(
        adapter,
        experiment_paths,
        experiment_label
    )
    if (!isTRUE(probe$found)) {
        return(list(enabled = FALSE, reason = "artifact_manifest_absent"))
    }
    policy$project_id <- probe$manifest$project_id
    descriptor <- adapter$descriptor
    identity <- descriptor$identity
    context <- createWorkflowContext(
        experiment_paths,
        omic_type = identity$omic_type,
        experiment_label = experiment_label,
        storage_policy = policy
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = identity$workflow_type,
        input_format = identity$input_format,
        data_level = identity$data_level,
        descriptor_catalogue = artifactWorkflowDescriptorCatalogue()
    )
    decision <- context$getStorageDecision()
    enabled <- identical(decision$effective_backend, "artifact") &&
        decision$effective_rollout %in% .WORKFLOW_ARTIFACT_ROLLOUTS &&
        identical(context$getIdentity()$project_id, probe$manifest$project_id) &&
        artifactStageIdentityMatches(context$getIdentity(), descriptor)
    if (!isTRUE(enabled)) {
        return(list(enabled = FALSE, reason = decision$reason_code))
    }
    list(
        enabled = TRUE,
        reason = "committed_artifact_project",
        context = context,
        descriptor = descriptor,
        manifest = probe$manifest
    )
}

#' Validate the immutable descriptor pin for a read-through project
#' @param adapter Validated read-through adapter.
#' @param store Validated artifact store.
#' @return The validated descriptor pin, invisibly.
#' @noRd
artifactStageValidateDescriptorPin <- function(adapter, store) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    relative <- artifactWorkflowStateDescriptorPinPath(store)
    path <- artifactStoreResolveFile(store, relative)
    if (!file.exists(path)) {
        artifactStageReadthroughAbort(
            adapter,
            "descriptor_pin_missing",
            sprintf(
                "%s read-through requires an immutable descriptor pin",
                adapter$owner_label
            )
        )
    }
    pin <- artifactStoreReadJson(store, relative)
    version <- workflowStateVersionValue(pin$schema_version)
    valid <- identical(pin$schema, ARTIFACT_DESCRIPTOR_PIN_SCHEMA) &&
        identical(version, ARTIFACT_DESCRIPTOR_PIN_VERSION) &&
        identical(pin$project_id, store$project_id) &&
        identical(
            pin$workflow_id,
            adapter$descriptor$identity$workflow_id
        ) && identical(
            pin$contract,
            artifactStageDescriptorContract(adapter$descriptor)
        )
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "descriptor_pin_mismatch",
            sprintf(
                "%s read-through descriptor differs from the project pin",
                adapter$owner_label
            )
        )
    }
    invisible(pin)
}

#' Open read-only resources for one descriptor-bound project
#' @param adapter Validated read-through adapter.
#' @param context Exact bound workflow context.
#' @param resource_policy Optional project registry resource policy.
#' @return Read-only store and registry resources.
#' @noRd
artifactStageReadResources <- function(
    adapter,
    context,
    resource_policy = NULL
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    valid <- inherits(context, "WorkflowContext") && context$isBound() &&
        artifactStageCoordinatorOwned(
            list(workflow_context = context),
            adapter$descriptor
        )
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_context",
            sprintf(
                "%s read-through requires its exact artifact context",
                adapter$owner_label
            )
        )
    }
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    artifactStageValidateDescriptorPin(adapter, store)
    registry_identity <- artifactRegistryIdentity(store, identity)
    registry <- projectRegistryForContext(context, resource_policy)
    session <- openProjectRegistryReadOnly(registry)
    list(
        identity = identity,
        registry_identity = registry_identity,
        store = store,
        registry = registry,
        session = session
    )
}

#' Validate a registry artifact against its immutable sidecar
#' @param adapter Validated read-through adapter.
#' @param resources Read-only artifact resources.
#' @param row One artifact registry row.
#' @param payload_validation Payload validation mode.
#' @return A normalized immutable artifact reference.
#' @noRd
artifactStageReadRegistryRef <- function(
    adapter,
    resources,
    row,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    payload_validation <- match.arg(payload_validation)
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
        validate_payload = identical(payload_validation, "materialized")
    )
    ref <- artifactStoreNormalizeRef(sidecar$artifact_ref)
    if (identical(payload_validation, "digest")) {
        payload_path <- artifactStoreResolveFile(
            resources$store,
            ref$relative_path,
            must_exist = TRUE
        )
        actual_shape <- ref$shape
        actual_shape$bytes <- unname(as.numeric(file.info(payload_path)$size))
        validateArtifactRefPayload(
            ref,
            resources$store$project_root,
            actual_shape
        )
    }
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
        artifactStageReadthroughAbort(
            adapter,
            "registry_ref_mismatch",
            sprintf(
                "%s registry artifact differs from its immutable sidecar",
                adapter$owner_label
            ),
            artifact_id = row$artifact_id,
            run_id = row$run_id
        )
    }
    ref
}

#' Convert one registry data-frame row to scalar fields
#' @param value Registry query result.
#' @param index One-based row index.
#' @return A list containing one scalar value per registry column.
#' @noRd
artifactStageDataFrameRow <- function(value, index) {
    row <- as.list(value[index, , drop = FALSE])
    lapply(row, `[[`, 1L)
}

#' Resolve one complete artifact stage run
#' @param adapter Validated read-through adapter.
#' @param resources Read-only artifact resources.
#' @param run One completed run registry row.
#' @param stage_id Required stage identifier.
#' @param required_roles Exact required state roles.
#' @param payload_validation Payload validation mode.
#' @return A complete stage record, or `NULL` when the run is incomplete.
#' @noRd
artifactStageResolveRun <- function(
    adapter,
    resources,
    run,
    stage_id,
    required_roles,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    payload_validation <- match.arg(payload_validation)
    links <- projectRegistryQuery(
        resources$session,
        "run_artifacts",
        filters = list(
            workflow_id = resources$registry_identity$workflow_id,
            run_id = run$run_id,
            direction = "output"
        )
    )
    complete <- nrow(links) == length(required_roles) &&
        setequal(links$artifact_role, required_roles) &&
        identical(
            sort(as.integer(links$ordinal)),
            seq_along(required_roles) - 1L
        )
    if (!isTRUE(complete)) return(NULL)
    refs <- lapply(required_roles, \(role) {
        link <- links[links$artifact_role == role, , drop = FALSE]
        if (nrow(link) != 1L) return(NULL)
        rows <- projectRegistryQuery(
            resources$session,
            "artifacts",
            filters = list(
                workflow_id = resources$registry_identity$workflow_id,
                artifact_id = link$artifact_id[[1L]]
            ),
            limit = 1L
        )
        if (nrow(rows) != 1L) return(NULL)
        row <- artifactStageDataFrameRow(rows, 1L)
        valid <- identical(row$run_id, run$run_id) &&
            identical(row$stage_id, stage_id) &&
            identical(row$state_role, role) &&
            identical(row$status, "committed")
        if (!isTRUE(valid)) return(NULL)
        artifactStageReadRegistryRef(
            adapter,
            resources,
            row,
            payload_validation = payload_validation
        )
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

#' Resolve the newest complete committed stage generation
#' @param adapter Validated read-through adapter.
#' @param resources Read-only artifact resources.
#' @param stage_id Required stage identifier.
#' @param required_roles Exact required state roles.
#' @param run_id Optional exact parent run identifier.
#' @param payload_validation Payload validation mode.
#' @return A complete committed stage record.
#' @noRd
resolveArtifactStageCommittedStage <- function(
    adapter,
    resources,
    stage_id,
    required_roles,
    run_id = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
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
            resolved <- artifactStageResolveRun(
                adapter,
                resources,
                run,
                stage_id,
                required_roles,
                payload_validation = payload_validation
            )
            if (!is.null(resolved)) return(resolved)
        }
    }
    artifactStageReadthroughAbort(
        adapter,
        "missing_committed_stage",
        sprintf(
            "%s stage '%s' has no complete committed artifact set",
            adapter$owner_label,
            stage_id
        ),
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

#' Decode and validate exact run parameters
#' @param adapter Validated read-through adapter.
#' @param resources Read-only artifact resources.
#' @param run_id Exact run identifier.
#' @param parameter_names Exact required parameter names.
#' @return A named list of decoded parameter values.
#' @noRd
artifactStageDecodeParameters <- function(
    adapter,
    resources,
    run_id,
    parameter_names
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    rows <- projectRegistryQuery(
        resources$session,
        "parameters",
        filters = list(
            workflow_id = resources$registry_identity$workflow_id,
            run_id = run_id
        )
    )
    rows <- rows[rows$parameter_name %in% parameter_names, , drop = FALSE]
    if (nrow(rows) != length(parameter_names) ||
        anyDuplicated(rows$parameter_name) > 0L ||
        !setequal(rows$parameter_name, parameter_names)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_run_parameters",
            sprintf(
                "%s artifact run parameters are missing or ambiguous",
                adapter$owner_label
            ),
            run_id = run_id
        )
    }
    values <- lapply(seq_len(nrow(rows)), \(index) {
        value <- jsonlite::fromJSON(
            rows$value_json[[index]],
            simplifyVector = FALSE
        )
        if (!identical(
            artifactSemanticDigest(value),
            rows$value_digest[[index]]
        )) {
            artifactStageReadthroughAbort(
                adapter,
                "invalid_run_parameters",
                sprintf(
                    "%s artifact run parameter digest is invalid",
                    adapter$owner_label
                ),
                run_id = run_id,
                parameter_name = rows$parameter_name[[index]]
            )
        }
        value
    })
    stats::setNames(values, rows$parameter_name)
}

#' Resolve the one current non-initial workflow state
#' @param adapter Validated read-through adapter.
#' @param resources Read-only artifact resources.
#' @return The current state registry row as scalar fields.
#' @noRd
artifactStageCurrentStateRow <- function(adapter, resources) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    rows <- projectRegistryQuery(
        resources$session,
        "states",
        filters = list(
            workflow_id = resources$registry_identity$workflow_id,
            status = "current"
        )
    )
    if (nrow(rows) != 1L || identical(rows$logical_name[[1L]], "initial")) {
        artifactStageReadthroughAbort(
            adapter,
            "missing_current_state",
            sprintf(
                "%s read-through requires one committed current design state",
                adapter$owner_label
            )
        )
    }
    artifactStageDataFrameRow(rows, 1L)
}

#' Collect one parent-linked import/design read-through evidence bundle
#' @param adapter Validated read-through adapter.
#' @param context Exact bound workflow context.
#' @param import_roles Exact required import roles.
#' @param design_roles Exact required design roles excluding the S4 state role.
#' @param import_parameter_names Exact required import parameters.
#' @param design_parameter_names Exact required design parameters.
#' @param resource_policy Optional project registry resource policy.
#' @param payload_validation Payload validation mode.
#' @return Immutable import, design, current-state, store, and identity evidence.
#' @noRd
collectArtifactStageResumeEvidence <- function(
    adapter,
    context,
    import_roles,
    design_roles,
    import_parameter_names,
    design_parameter_names,
    resource_policy = NULL,
    payload_validation = c("materialized", "digest", "sidecar")
) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    payload_validation <- match.arg(payload_validation)
    resources <- artifactStageReadResources(
        adapter,
        context,
        resource_policy
    )
    on.exit(closeProjectRegistry(resources$session), add = TRUE)
    design <- resolveArtifactStageCommittedStage(
        adapter,
        resources,
        "design",
        design_roles,
        payload_validation = payload_validation
    )
    design$parameters <- artifactStageDecodeParameters(
        adapter,
        resources,
        design$run_id,
        design_parameter_names
    )
    parent_run_id <- design$parameters$parent_import_run_id
    if (!workflowCapabilityScalarString(parent_run_id)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            sprintf(
                "%s design artifacts do not declare one parent import run",
                adapter$owner_label
            ),
            run_id = design$run_id
        )
    }
    import <- resolveArtifactStageCommittedStage(
        adapter,
        resources,
        "import",
        import_roles,
        run_id = parent_run_id,
        payload_validation = payload_validation
    )
    import$parameters <- artifactStageDecodeParameters(
        adapter,
        resources,
        import$run_id,
        import_parameter_names
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

#' Validate a stored artifact reference against its sidecar
#' @param adapter Validated read-through adapter.
#' @param store Validated artifact store.
#' @param ref Artifact reference to validate.
#' @return The normalized reference and immutable sidecar.
#' @noRd
artifactStageValidateStoredRef <- function(adapter, store, ref) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    ref <- artifactStoreNormalizeRef(ref)
    sidecar_path <- artifactStoreManagedPaths(
        store,
        ref$logical_key,
        ref$artifact_id
    )$sidecar
    sidecar <- artifactStoreReadSidecar(
        store,
        sidecar_path,
        validate_payload = FALSE
    )
    if (!identical(artifactStoreNormalizeRef(sidecar$artifact_ref), ref)) {
        artifactStageReadthroughAbort(
            adapter,
            "registry_ref_mismatch",
            sprintf(
                "%s artifact ref differs from its immutable sidecar",
                adapter$owner_label
            ),
            artifact_id = ref$artifact_id
        )
    }
    list(ref = ref, sidecar = sidecar)
}

#' Independently hydrate one rectangular artifact reference
#' @param adapter Validated read-through adapter.
#' @param store Validated artifact store.
#' @param ref Artifact reference to hydrate.
#' @return The decoded rectangular value.
#' @noRd
artifactStageReadTable <- function(adapter, store, ref) {
    validated <- artifactStageValidateStoredRef(adapter, store, ref)
    payload_path <- artifactStoreResolveFile(
        store,
        validated$ref$relative_path,
        must_exist = TRUE
    )
    payload <- arrow::read_parquet(payload_path, as_data_frame = FALSE)
    value <- decodeArtifactRectangular(
        payload,
        validated$sidecar$codec_metadata
    )
    validateArtifactRefPayload(
        validated$ref,
        store$project_root,
        list(
            kind = validated$sidecar$codec_metadata$kind,
            rows = nrow(value),
            columns = ncol(value),
            payloads = 1L,
            bytes = unname(as.numeric(file.info(payload_path)$size))
        )
    )
    value
}

#' Hydrate selected roles from one complete stage
#' @param adapter Validated read-through adapter.
#' @param store Validated artifact store.
#' @param stage Complete stage record.
#' @param roles Exact roles to hydrate.
#' @return A named list of independently hydrated tables.
#' @noRd
readArtifactStageRoles <- function(adapter, store, stage, roles) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    missing <- setdiff(roles, names(stage$refs))
    if (length(missing) > 0L) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            sprintf(
                "%s stage '%s' is missing role(s): %s",
                adapter$owner_label,
                stage$stage_id,
                paste(missing, collapse = ", ")
            ),
            stage_id = stage$stage_id,
            missing_roles = missing
        )
    }
    values <- lapply(
        stage$refs[roles],
        \(ref) artifactStageReadTable(adapter, store, ref)
    )
    names(values) <- roles
    values
}

#' Hydrate every role from one complete generation
#' @param adapter Validated read-through adapter.
#' @param store Validated artifact store.
#' @param stage Complete stage record.
#' @return A named list of independently hydrated tables.
#' @noRd
readArtifactStage <- function(adapter, store, stage) {
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    values <- lapply(
        stage$refs,
        \(ref) artifactStageReadTable(adapter, store, ref)
    )
    names(values) <- names(stage$refs)
    generations <- unique(vapply(
        stage$refs,
        \(ref) ref$logical_key$generation_id,
        character(1)
    ))
    if (length(generations) != 1L ||
        !identical(generations[[1L]], stage$generation_id)) {
        artifactStageReadthroughAbort(
            adapter,
            "mixed_generation",
            sprintf(
                "%s stage hydration combined different generations",
                adapter$owner_label
            ),
            stage_id = stage$stage_id,
            run_id = stage$run_id
        )
    }
    values
}

#' Test a serialized read-through contract version
#' @param value Serialized version value.
#' @param expected Expected integer contract version.
#' @return A scalar logical.
#' @noRd
artifactStageReadthroughContractVersion <- function(value, expected = 1L) {
    identical(workflowStateVersionValue(value), as.integer(expected))
}

#' Normalize a named artifact reference list
#' @param refs Named artifact references.
#' @return Named normalized artifact references.
#' @noRd
artifactStageNormalizedRefs <- function(refs) {
    normalized <- lapply(refs, artifactStoreNormalizeRef)
    names(normalized) <- names(refs)
    normalized
}

#' Restore one optional artifact table
#' @param adapter Validated read-through adapter.
#' @param value Hydrated optional table.
#' @param available Provenance availability flag.
#' @param role Artifact role.
#' @return The hydrated table or `NULL` when absence is proven.
#' @noRd
artifactStageRestoreOptional <- function(adapter, value, available, role) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_optional_payload",
            sprintf(
                "%s '%s' availability flag is invalid",
                adapter$owner_label,
                role
            )
        )
    }
    if (isTRUE(available)) return(value)
    expected <- artifactStageOptionalTable(
        NULL,
        role,
        adapter$abort_fn
    )
    if (!identical(value, expected)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_optional_payload",
            sprintf(
                "%s '%s' absence marker is invalid",
                adapter$owner_label,
                role
            )
        )
    }
    NULL
}

#' Restore the exact contrasts representation from provenance
#' @param adapter Validated read-through adapter.
#' @param value Hydrated contrasts table.
#' @param kind Recorded contrasts representation.
#' @return A data frame, character vector, or `NULL`.
#' @noRd
artifactStageRestoreContrasts <- function(adapter, value, kind) {
    if (identical(kind, "data.frame")) return(value)
    if (identical(kind, "character") && identical(names(value), "contrasts")) {
        return(as.character(value$contrasts))
    }
    expected <- artifactStageOptionalTable(
        NULL,
        "contrasts",
        adapter$abort_fn
    )
    if (identical(kind, "null") && identical(value, expected)) return(NULL)
    artifactStageReadthroughAbort(
        adapter,
        "invalid_contrasts",
        sprintf(
            "%s contrast representation is incompatible with provenance",
            adapter$owner_label
        )
    )
}

#' Restore an exact serialized import column mapping
#' @param adapter Validated read-through adapter.
#' @param parameters Decoded import run parameters.
#' @return The restored column mapping list.
#' @noRd
artifactStageColumnMapping <- function(adapter, parameters) {
    encoded <- parameters$column_mapping_serialized
    mapping <- artifactWorkflowStateUnserializeMetadata(
        encoded,
        sprintf("%s import column mapping", adapter$owner_label)
    )
    if (!is.list(mapping)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_column_mapping",
            sprintf(
                "%s import column mapping is not an R list",
                adapter$owner_label
            )
        )
    }
    mapping
}
