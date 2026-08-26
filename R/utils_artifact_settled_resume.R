# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_SETTLED_RESUME_SCHEMA <- "multischolar.artifact_settled_resume"
.ARTIFACT_SETTLED_RESUME_VERSION <- 1L

artifactSettledResumeRelativePath <- function(store) {
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$workflow_state_root,
        "settled-resume.json"
    ))
}

artifactSettledResumePointerRelativePath <- function(store) {
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$workflow_state_root,
        "settled-resume-latest.json"
    ))
}

artifactSettledResumeCacheRelativePath <- function(store, digest) {
    artifactRefValidateDigest(digest, "settled_resume_cache_digest")
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$workflow_state_root,
        paste0("settled-resume-cache-", digest, ".rds")
    ))
}

artifactSettledResumeRefLocator <- function(store, ref) {
    ref <- artifactStoreNormalizeRef(ref)
    structure(
        list(
            sidecar_relative_path = artifactStoreManagedPaths(
                store,
                ref$logical_key,
                ref$artifact_id
            )$sidecar,
            payload_relative_path = ref$relative_path,
            payload_bytes = as.numeric(ref$shape$bytes),
            proof = artifactPayloadReadthroughRefProof(list(ref))[[1L]]
        ),
        class = c("ArtifactSettledResumeRefLocator", "list")
    )
}

artifactSettledResumeLeanSnapshot <- function(snapshot, store) {
    lean <- snapshot
    lean$import$refs <- lapply(
        snapshot$import$refs,
        \(ref) artifactSettledResumeRefLocator(store, ref)
    )
    lean$design$refs <- lapply(
        snapshot$design$refs,
        \(ref) artifactSettledResumeRefLocator(store, ref)
    )
    metadata <- lean$state_bootstrap$state_metadata
    for (generation_id in names(metadata)) {
        entry <- metadata[[generation_id]]
        entry$artifact_refs <- lapply(
            entry$artifact_refs,
            \(ref) artifactSettledResumeRefLocator(store, ref)
        )
        audit <- entry$audit_metadata
        if (is.list(audit) && is.list(audit$stage_artifact_refs)) {
            audit$stage_artifact_refs <- lapply(
                audit$stage_artifact_refs,
                \(ref) artifactSettledResumeRefLocator(store, ref)
            )
            entry$audit_metadata <- audit
        }
        metadata[[generation_id]] <- entry
    }
    lean$state_bootstrap$state_metadata <- metadata
    lean
}

artifactSettledResumeReadLocator <- function(adapter, store, locator) {
    if (!inherits(locator, "ArtifactSettledResumeRefLocator")) {
        return(artifactStageReadTable(adapter, store, locator))
    }
    sidecar <- artifactStoreReadSidecar(
        store,
        locator$sidecar_relative_path,
        validate_payload = FALSE
    )
    ref <- artifactStoreNormalizeRef(sidecar$artifact_ref)
    if (!identical(
        artifactPayloadReadthroughRefProof(list(ref))[[1L]],
        locator$proof
    )) {
        artifactStageReadthroughAbort(
            adapter,
            "registry_ref_mismatch",
            "settled artifact locator differs from its immutable sidecar"
        )
    }
    artifactStageReadTable(adapter, store, ref)
}

artifactWorkflowStateBootstrapFromExport <- function(store, state_export) {
    state_export <- artifactWorkflowStateExportSnapshot(
        identity = list(
            project_id = state_export$project_id,
            workflow_id = state_export$workflow_id
        ),
        current_generation_id = state_export$current_generation_id,
        current_state = state_export$current_state,
        lineage = state_export$active_lineage,
        workflow_type = state_export$workflow_type,
        audit_enabled = state_export$audit_enabled
    )
    lineage <- state_export$active_lineage
    manifests <- lapply(lineage, \(entry) {
        artifactWorkflowStateReadManifest(store, entry$manifest_relative_path)
    })
    rows <- lapply(seq_along(lineage), \(index) {
        entry <- lineage[[index]]
        manifest <- manifests[[index]]
        list(
            generation_id = manifest$generation_id,
            parent_generation_id = manifest$parent_generation_id,
            logical_name = manifest$logical_name,
            manifest_relative_path = entry$manifest_relative_path,
            status = if (index == length(lineage)) "current" else "historical"
        )
    })
    names_by_generation <- stats::setNames(
        vapply(rows, `[[`, character(1), "logical_name"),
        vapply(rows, `[[`, character(1), "generation_id")
    )
    metadata <- lapply(seq_along(rows), \(index) {
        row <- rows[[index]]
        manifest <- manifests[[index]]
        list(
            timestamp = as.POSIXct(manifest$created_at, tz = "UTC"),
            config = artifactWorkflowStateUnserializeMetadata(
                manifest$config_json,
                "config"
            ),
            description = manifest$description,
            previous_state = if (is.null(manifest$parent_generation_id)) {
                NULL
            } else {
                unname(names_by_generation[[manifest$parent_generation_id]])
            },
            s4_class = manifest$s4_class,
            audit_metadata = artifactWorkflowStateUnserializeMetadata(
                manifest$audit_json,
                "audit metadata"
            ),
            generation_id = manifest$generation_id,
            artifact_refs = manifest$data$artifact_refs,
            status = row$status
        )
    })
    names(metadata) <- vapply(rows, `[[`, character(1), "generation_id")
    list(
        state_rows = rows,
        state_metadata = metadata,
        workflow_type = state_export$workflow_type,
        audit_enabled = state_export$audit_enabled
    )
}

artifactWorkflowStateValidateSettledBootstrap <- function(
    bootstrap,
    store,
    identity
) {
    expected_bootstrap_names <- c(
        "state_rows", "state_metadata", "workflow_type", "audit_enabled"
    )
    valid <- is.list(bootstrap) &&
        identical(names(bootstrap), expected_bootstrap_names) &&
        is.list(bootstrap$state_rows) && length(bootstrap$state_rows) > 0L
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "settled artifact WorkflowState bootstrap is malformed",
            "multischolar_invalid_artifact_state_bootstrap"
        )
    }
    expected_names <- c(
        "generation_id", "parent_generation_id", "logical_name",
        "manifest_relative_path", "status"
    )
    generation_ids <- vapply(
        bootstrap$state_rows,
        `[[`,
        character(1),
        "generation_id"
    )
    valid <- is.list(bootstrap$state_metadata) &&
        identical(names(bootstrap$state_metadata), generation_ids) &&
        workflowStateScalarString(bootstrap$workflow_type) &&
        is.logical(bootstrap$audit_enabled) &&
        length(bootstrap$audit_enabled) == 1L && !is.na(bootstrap$audit_enabled)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "settled artifact WorkflowState metadata is malformed",
            "multischolar_invalid_artifact_state_bootstrap"
        )
    }
    current <- 0L
    deep_validation <- !inherits(
        bootstrap,
        "ValidatedArtifactWorkflowStateBootstrap"
    )
    for (row in bootstrap$state_rows) {
        valid <- is.list(row) && identical(names(row), expected_names) &&
            row$status %in% c("current", "historical")
        if (!isTRUE(valid)) {
            artifactWorkflowStateAbort(
                "settled artifact WorkflowState row is malformed",
                "multischolar_invalid_artifact_state_bootstrap"
            )
        }
        metadata <- bootstrap$state_metadata[[row$generation_id]]
        metadata_valid <- is.list(metadata) &&
            identical(metadata$generation_id, row$generation_id) &&
            identical(metadata$status, row$status)
        if (!isTRUE(metadata_valid)) {
            artifactWorkflowStateAbort(
                "settled artifact WorkflowState metadata differs from its row",
                "multischolar_invalid_artifact_state_bootstrap"
            )
        }
        if (isTRUE(deep_validation)) {
            manifest <- artifactWorkflowStateReadManifest(
                store,
                row$manifest_relative_path
            )
            matches <- identical(manifest$project_id, identity$project_id) &&
                identical(manifest$workflow_id, identity$workflow_id) &&
                identical(manifest$generation_id, row$generation_id) &&
                identical(manifest$parent_generation_id, row$parent_generation_id) &&
                identical(manifest$logical_name, row$logical_name) &&
                identical(manifest$data$artifact_refs, metadata$artifact_refs)
            if (!isTRUE(matches)) {
                artifactWorkflowStateAbort(
                    "settled artifact WorkflowState row differs from its manifest",
                    "multischolar_invalid_artifact_state_bootstrap"
                )
            }
        }
        if (identical(row$status, "current")) current <- current + 1L
    }
    if (current != 1L) {
        artifactWorkflowStateAbort(
            "settled artifact WorkflowState requires one current row",
            "multischolar_invalid_artifact_state_bootstrap"
        )
    }
    structure(
        bootstrap,
        class = c("ValidatedArtifactWorkflowStateBootstrap", "list")
    )
}

artifactSettledResumeSnapshot <- function(
    context,
    descriptor,
    evidence,
    state_export,
    adapter
) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    assay_roles <- unlist(
        evidence$import$parameters$assay_roles,
        use.names = FALSE
    )
    cleaned_roles <- sprintf(
        "cleaned_assay_%04d",
        seq_along(evidence$design$parameters$assay_order)
    )
    settled_tables <- list(
        import = readArtifactStageRoles(
            adapter,
            store,
            evidence$import,
            setdiff(names(evidence$import$refs), assay_roles)
        ),
        design = readArtifactStageRoles(
            adapter,
            store,
            evidence$design,
            setdiff(names(evidence$design$refs), cleaned_roles)
        )
    )
    list(
        schema = .ARTIFACT_SETTLED_RESUME_SCHEMA,
        schema_version = .ARTIFACT_SETTLED_RESUME_VERSION,
        descriptor_contract = artifactStageDescriptorContract(descriptor),
        identity = identity,
        import = evidence$import,
        design = evidence$design,
        current_state = evidence$current_state,
        state_bootstrap = artifactWorkflowStateBootstrapFromExport(
            store,
            state_export
        ),
        settled_tables = settled_tables,
        created_at = artifactRefUtcNow()
    )
}

artifactSettledResumeWrite <- function(context, descriptor, snapshot) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    relative_path <- artifactSettledResumeRelativePath(store)
    serialized <- artifactWorkflowStateSerializeMetadata(
        snapshot,
        "settled artifact resume snapshot"
    )
    cache_temporary <- tempfile("settled-resume-cache-", tmpdir = store$project_root)
    on.exit(if (file.exists(cache_temporary)) unlink(cache_temporary), add = TRUE)
    lean_snapshot <- artifactSettledResumeLeanSnapshot(snapshot, store)
    saveRDS(lean_snapshot, cache_temporary, version = 3L, compress = FALSE)
    cache_digest <- artifactByteDigest(cache_temporary)
    cache_relative <- artifactSettledResumeCacheRelativePath(
        store,
        cache_digest
    )
    cache_target <- artifactStoreResolveFile(store, cache_relative)
    artifactStoreEnsureDirectory(
        store,
        artifactNormalizeRelativePath(dirname(cache_relative))
    )
    if (!file.exists(cache_target) &&
        !isTRUE(file.rename(cache_temporary, cache_target))) {
        artifactStagePersistenceAbort(
            "settled resume cache could not be published",
            "multischolar_settled_resume_publish_failed"
        )
    }
    if (!identical(artifactByteDigest(cache_target), cache_digest)) {
        artifactStagePersistenceAbort(
            "settled resume cache digest is invalid",
            "multischolar_invalid_settled_resume_snapshot"
        )
    }
    wrapper <- list(
        schema = .ARTIFACT_SETTLED_RESUME_SCHEMA,
        schema_version = .ARTIFACT_SETTLED_RESUME_VERSION,
        payload = serialized,
        payload_digest = digest::digest(
            serialized,
            algo = "sha256",
            serialize = FALSE
        ),
        cache_relative_path = cache_relative,
        cache_byte_digest = cache_digest
    )
    temporary <- artifactNormalizeRelativePath(paste0(
        relative_path,
        ".",
        artifactOpaqueId("tmp"),
        ".tmp"
    ))
    artifactStoreWriteJson(store, wrapper, temporary)
    on.exit({
        path <- artifactStoreResolveFile(store, temporary)
        if (file.exists(path)) unlink(path, force = FALSE)
    }, add = TRUE)
    candidate <- artifactStoreReadJson(store, temporary)
    if (!identical(candidate$payload_digest, wrapper$payload_digest) ||
        !identical(candidate$cache_byte_digest, wrapper$cache_byte_digest)) {
        artifactStagePersistenceAbort(
            "settled resume snapshot failed pre-publication validation",
            "multischolar_invalid_settled_resume_snapshot"
        )
    }
    source <- artifactStoreResolveFile(store, temporary, must_exist = TRUE)
    target <- artifactStoreResolveFile(store, relative_path)
    artifactStoreEnsureDirectory(
        store,
        artifactNormalizeRelativePath(dirname(relative_path))
    )
    if (!isTRUE(file.rename(source, target))) {
        artifactStagePersistenceAbort(
            "settled resume snapshot could not be atomically published",
            "multischolar_settled_resume_publish_failed"
        )
    }
    pointer <- list(
        schema = .ARTIFACT_SETTLED_RESUME_SCHEMA,
        schema_version = .ARTIFACT_SETTLED_RESUME_VERSION,
        descriptor_contract = artifactStageDescriptorContract(descriptor),
        identity_digest = artifactSemanticDigest(identity),
        authoritative_relative_path = relative_path,
        authoritative_byte_digest = artifactByteDigest(target),
        cache_relative_path = cache_relative,
        cache_byte_digest = cache_digest
    )
    pointer_relative <- artifactSettledResumePointerRelativePath(store)
    pointer_temporary <- artifactNormalizeRelativePath(paste0(
        pointer_relative,
        ".",
        artifactOpaqueId("tmp"),
        ".tmp"
    ))
    artifactStoreWriteJson(store, pointer, pointer_temporary)
    on.exit({
        path <- artifactStoreResolveFile(store, pointer_temporary)
        if (file.exists(path)) unlink(path, force = FALSE)
    }, add = TRUE)
    pointer_source <- artifactStoreResolveFile(
        store,
        pointer_temporary,
        must_exist = TRUE
    )
    pointer_target <- artifactStoreResolveFile(store, pointer_relative)
    if (!isTRUE(file.rename(pointer_source, pointer_target))) {
        artifactStagePersistenceAbort(
            "settled resume pointer could not be atomically published",
            "multischolar_settled_resume_publish_failed"
        )
    }
    invisible(relative_path)
}

artifactSettledResumeValidateRef <- function(adapter, store, ref) {
    validated <- artifactStageValidateStoredRef(adapter, store, ref)
    path <- artifactStoreResolveFile(
        store,
        validated$ref$relative_path,
        must_exist = TRUE
    )
    shape <- validated$ref$shape
    shape$bytes <- unname(as.numeric(file.info(path)$size))
    validateArtifactRefPayload(validated$ref, store$project_root, shape)
    validated$ref
}

artifactSettledResumeValidateRefs <- function(adapter, store, refs) {
    if (!is.list(refs) || length(refs) == 0L) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume snapshot has no immutable artifact refs"
        )
    }
    lapply(refs, \(ref) artifactSettledResumeValidateRef(adapter, store, ref))
}

artifactSettledResumeReadDirect <- function(context, descriptor, adapter) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    relative_path <- artifactSettledResumeRelativePath(store)
    path <- artifactStoreResolveFile(store, relative_path)
    if (!file.exists(path)) return(NULL)
    artifactStageValidateDescriptorPin(adapter, store)
    wrapper <- artifactStoreReadJson(store, relative_path)
    valid <- identical(wrapper$schema, .ARTIFACT_SETTLED_RESUME_SCHEMA) &&
        identical(
            as.integer(wrapper$schema_version),
            .ARTIFACT_SETTLED_RESUME_VERSION
        ) && workflowCapabilityScalarString(wrapper$payload) &&
        workflowCapabilityScalarString(wrapper$payload_digest) &&
        workflowCapabilityScalarString(wrapper$cache_relative_path) &&
        workflowCapabilityScalarString(wrapper$cache_byte_digest)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume snapshot wrapper is malformed"
        )
    }
    snapshot <- artifactWorkflowStateUnserializeMetadata(
        wrapper$payload,
        "settled artifact resume snapshot"
    )
    expected_names <- c(
        "schema", "schema_version", "descriptor_contract", "identity",
        "import", "design", "current_state", "state_bootstrap",
        "settled_tables", "created_at"
    )
    valid <- is.list(snapshot) && identical(names(snapshot), expected_names) &&
        identical(snapshot$schema, .ARTIFACT_SETTLED_RESUME_SCHEMA) &&
        identical(snapshot$schema_version, .ARTIFACT_SETTLED_RESUME_VERSION) &&
        identical(snapshot$descriptor_contract,
            artifactStageDescriptorContract(descriptor)) &&
        identical(snapshot$identity, identity) &&
        identical(
            digest::digest(
                wrapper$payload,
                algo = "sha256",
                serialize = FALSE
            ),
            wrapper$payload_digest
        )
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume snapshot identity or digest is invalid"
        )
    }
    snapshot$import$refs <- artifactSettledResumeValidateRefs(
        adapter,
        store,
        snapshot$import$refs
    )
    snapshot$design$refs <- artifactSettledResumeValidateRefs(
        adapter,
        store,
        snapshot$design$refs
    )
    assay_roles <- unlist(
        snapshot$import$parameters$assay_roles,
        use.names = FALSE
    )
    cleaned_roles <- sprintf(
        "cleaned_assay_%04d",
        seq_along(snapshot$design$parameters$assay_order)
    )
    actual_tables <- list(
        import = readArtifactStageRoles(
            adapter,
            store,
            snapshot$import,
            setdiff(names(snapshot$import$refs), assay_roles)
        ),
        design = readArtifactStageRoles(
            adapter,
            store,
            snapshot$design,
            setdiff(names(snapshot$design$refs), cleaned_roles)
        )
    )
    if (!identical(actual_tables, snapshot$settled_tables)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume metadata tables differ from immutable artifacts"
        )
    }
    bootstrap <- artifactWorkflowStateValidateSettledBootstrap(
        snapshot$state_bootstrap,
        store,
        identity
    )
    for (row in bootstrap$state_rows) {
        manifest <- artifactWorkflowStateReadManifest(
            store,
            row$manifest_relative_path
        )
        if (length(manifest$data$artifact_refs) > 0L) {
            artifactSettledResumeValidateRefs(
                adapter,
                store,
                manifest$data$artifact_refs
            )
        }
    }
    snapshot$state_bootstrap <- bootstrap
    snapshot$store <- store
    snapshot
}

artifactSettledResumeRead <- function(context, descriptor, adapter) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    adapter <- validateArtifactStageReadthroughAdapter(adapter)
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    relative_path <- artifactSettledResumePointerRelativePath(store)
    path <- artifactStoreResolveFile(store, relative_path)
    if (!file.exists(path)) return(NULL)
    artifactStageValidateDescriptorPin(adapter, store)
    pointer <- artifactStoreReadJson(store, relative_path)
    valid <- identical(pointer$schema, .ARTIFACT_SETTLED_RESUME_SCHEMA) &&
        identical(
            as.integer(pointer$schema_version),
            .ARTIFACT_SETTLED_RESUME_VERSION
        ) && identical(
            pointer$descriptor_contract,
            artifactStageDescriptorContract(descriptor)
        ) && identical(
            pointer$identity_digest,
            artifactSemanticDigest(identity)
        ) && workflowCapabilityScalarString(
            pointer$authoritative_relative_path
        ) && workflowCapabilityScalarString(pointer$authoritative_byte_digest) &&
        workflowCapabilityScalarString(pointer$cache_relative_path) &&
        workflowCapabilityScalarString(pointer$cache_byte_digest)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume snapshot wrapper is malformed"
        )
    }
    cache_path <- artifactStoreResolveFile(
        store,
        pointer$cache_relative_path,
        must_exist = TRUE
    )
    authoritative <- artifactStoreResolveFile(
        store,
        pointer$authoritative_relative_path,
        must_exist = TRUE
    )
    if (!identical(
        artifactByteDigest(authoritative),
        pointer$authoritative_byte_digest
    ) || !identical(
        artifactByteDigest(cache_path),
        pointer$cache_byte_digest
    )) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume acceleration cache digest is invalid"
        )
    }
    snapshot <- tryCatch(
        readRDS(cache_path),
        error = \(error) {
            artifactStageReadthroughAbort(
                adapter,
                "incomplete_contract",
                "settled resume acceleration cache is unreadable",
                parent = error
            )
        }
    )
    expected_names <- c(
        "schema", "schema_version", "descriptor_contract", "identity",
        "import", "design", "current_state", "state_bootstrap",
        "settled_tables", "created_at"
    )
    valid <- is.list(snapshot) && identical(names(snapshot), expected_names) &&
        identical(snapshot$schema, .ARTIFACT_SETTLED_RESUME_SCHEMA) &&
        identical(snapshot$schema_version, .ARTIFACT_SETTLED_RESUME_VERSION) &&
        identical(snapshot$descriptor_contract,
            artifactStageDescriptorContract(descriptor)) &&
        identical(snapshot$identity, identity)
    if (!isTRUE(valid)) {
        artifactStageReadthroughAbort(
            adapter,
            "incomplete_contract",
            "settled resume snapshot identity or digest is invalid"
        )
    }
    snapshot$state_bootstrap <- structure(
        snapshot$state_bootstrap,
        class = c("ValidatedArtifactWorkflowStateBootstrap", "list")
    )
    snapshot$state_bootstrap <- artifactWorkflowStateValidateSettledBootstrap(
        snapshot$state_bootstrap,
        store,
        identity
    )
    for (row in snapshot$state_bootstrap$state_rows) {
        artifactStoreResolveFile(
            store,
            row$manifest_relative_path,
            must_exist = TRUE
        )
    }
    snapshot$store <- store
    snapshot
}
