# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

ARTIFACT_WORKFLOW_STATE_SCHEMA <- "multischolar.artifact_workflow_state"
ARTIFACT_WORKFLOW_STATE_VERSION <- 1L
ARTIFACT_WORKFLOW_STATE_METADATA_LIMIT <- 1024L * 1024L
ARTIFACT_DESCRIPTOR_PIN_SCHEMA <- "multischolar.artifact_descriptor_pin"
ARTIFACT_DESCRIPTOR_PIN_VERSION <- 1L

artifactWorkflowStateAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_workflow_state_error"),
        ...
    )
}

artifactWorkflowStateValidateContext <- function(context) {
    if (!inherits(context, "WorkflowContext") || !context$isBound()) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState requires a bound workflow context",
            "multischolar_invalid_artifact_workflow_context"
        )
    }
    decision <- context$getStorageDecision()
    if (!identical(decision$effective_backend, "artifact")) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState requires an artifact backend decision",
            "multischolar_invalid_artifact_workflow_context"
        )
    }
    if (!inherits(context$getPaths(), "MultiScholaRArtifactPaths")) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState requires validated project paths",
            "multischolar_invalid_artifact_workflow_context"
        )
    }
    context
}

artifactWorkflowStateUnsafeName <- function(value) {
    normalized <- tolower(value)
    normalized %in% c("cache", "caches", "connection", "connections") |
        grepl("_(cache|connection)s?$", normalized)
}

artifactWorkflowStateAbsolutePath <- function(value) {
    is.character(value) && any(grepl("^(/|[A-Za-z]:[/\\\\])", value))
}

artifactWorkflowStateAssertSafeMetadata <- function(value, owner = "metadata") {
    unsafe_type <- is.environment(value) || is.function(value) ||
        typeof(value) %in% c("externalptr", "weakref") ||
        inherits(value, "connection") || isS4(value)
    if (isTRUE(unsafe_type)) {
        artifactWorkflowStateAbort(
            sprintf("artifact state %s contains a runtime-only value", owner),
            "multischolar_unsafe_artifact_state_metadata",
            owner = owner
        )
    }
    if (artifactWorkflowStateAbsolutePath(value)) {
        artifactWorkflowStateAbort(
            sprintf("artifact state %s contains an absolute path", owner),
            "multischolar_absolute_path_in_artifact_state",
            owner = owner
        )
    }
    if (!is.list(value)) return(invisible(TRUE))
    value_names <- names(value)
    if (!is.null(value_names) && any(artifactWorkflowStateUnsafeName(value_names))) {
        artifactWorkflowStateAbort(
            sprintf("artifact state %s contains a cache or connection field", owner),
            "multischolar_unsafe_artifact_state_metadata",
            owner = owner
        )
    }
    for (index in seq_along(value)) {
        label <- if (!is.null(value_names) && nzchar(value_names[[index]])) {
            value_names[[index]]
        } else {
            as.character(index)
        }
        artifactWorkflowStateAssertSafeMetadata(
            value[[index]],
            paste0(owner, "[[", label, "]]")
        )
    }
    invisible(TRUE)
}

artifactWorkflowStateSerializeMetadata <- function(value, owner) {
    artifactWorkflowStateAssertSafeMetadata(value, owner)
    encoded <- as.character(jsonlite::serializeJSON(value, digits = NA))
    if (nchar(encoded, type = "bytes") > ARTIFACT_WORKFLOW_STATE_METADATA_LIMIT) {
        artifactWorkflowStateAbort(
            sprintf("artifact state %s exceeds the inline metadata limit", owner),
            "multischolar_artifact_state_metadata_too_large",
            owner = owner
        )
    }
    encoded
}

artifactWorkflowStateValidateSaveMetadata <- function(config, audit_metadata, description) {
    artifactWorkflowStateSerializeMetadata(config, "config")
    artifactWorkflowStateSerializeMetadata(audit_metadata, "audit metadata")
    if (!is.character(description) || length(description) != 1L || is.na(description)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState description must be one string",
            "multischolar_invalid_artifact_state_metadata"
        )
    }
    invisible(TRUE)
}

artifactWorkflowStateUnserializeMetadata <- function(value, owner) {
    if (!workflowStateScalarString(value) || !jsonlite::validate(value)) {
        artifactWorkflowStateAbort(
            sprintf("artifact state %s is not valid serialized metadata", owner),
            "multischolar_malformed_artifact_state_manifest",
            owner = owner
        )
    }
    decoded <- tryCatch(
        jsonlite::unserializeJSON(value),
        error = \(error) artifactWorkflowStateAbort(
            sprintf("artifact state %s could not be decoded", owner),
            "multischolar_malformed_artifact_state_manifest",
            owner = owner,
            parent = error
        )
    )
    artifactWorkflowStateAssertSafeMetadata(decoded, owner)
    decoded
}

artifactWorkflowStateManifestRelativePath <- function(store, generation_id) {
    artifactRefValidateId(generation_id, "generation_id", "gen")
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$generations,
        generation_id,
        "workflow-state.json"
    ))
}

artifactWorkflowStateEnsureRootManifest <- function(store, identity) {
    relative_path <- store$relative_paths$artifact_manifest
    path <- artifactStoreResolveFile(store, relative_path)
    if (file.exists(path)) {
        manifest <- artifactStoreReadJson(store, relative_path)
        version <- workflowStateVersionValue(manifest$schema_version)
        valid <- identical(manifest$schema, .ARTIFACT_MANIFEST_SCHEMA) &&
            identical(version, .ARTIFACT_MANIFEST_SCHEMA_VERSION) &&
            identical(manifest$project_id, identity$project_id) &&
            identical(manifest$workflow_id, identity$workflow_id)
        if (!isTRUE(valid)) {
            artifactWorkflowStateAbort(
                "artifact workflow root manifest is incompatible",
                "multischolar_incompatible_artifact_workflow_manifest"
            )
        }
        return(invisible(FALSE))
    }
    manifest <- list(
        schema = .ARTIFACT_MANIFEST_SCHEMA,
        schema_version = .ARTIFACT_MANIFEST_SCHEMA_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        omic_type = identity$omic_type,
        workflow_slug = identity$workflow_slug,
        created_at = artifactRefUtcNow()
    )
    temporary_path <- artifactNormalizeRelativePath(paste0(
        relative_path,
        ".",
        artifactOpaqueId("tmp"),
        ".tmp"
    ))
    artifactStoreWriteJson(store, manifest, temporary_path)
    on.exit({
        temporary <- artifactStoreResolveFile(store, temporary_path)
        if (file.exists(temporary)) unlink(temporary, force = FALSE)
    }, add = TRUE)
    artifactStorePublishFile(store, temporary_path, relative_path)
    invisible(TRUE)
}

artifactWorkflowStateDescriptorPinPath <- function(store) {
    artifactNormalizeRelativePath(file.path(
        store$relative_paths$workflow_state_root,
        "artifact-descriptor.json"
    ))
}

artifactWorkflowStateValidateDescriptorContract <- function(contract) {
    required <- c("descriptor_id", "descriptor_version", "descriptor_digest")
    if (!is.list(contract) || !identical(names(contract), required) ||
        !all(vapply(contract, workflowCapabilityScalarString, logical(1)))) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState descriptor contract is malformed",
            "multischolar_invalid_artifact_state_descriptor_contract"
        )
    }
    artifactDescriptorSemver(contract$descriptor_version, "descriptor_version")
    artifactRefValidateDigest(contract$descriptor_digest, "descriptor_digest")
    contract
}

artifactWorkflowStateEnsureDescriptorPin <- function(
    store,
    identity,
    contract = NULL,
    allow_create = FALSE
) {
    relative_path <- artifactWorkflowStateDescriptorPinPath(store)
    path <- artifactStoreResolveFile(store, relative_path)
    if (is.null(contract)) {
        if (file.exists(path)) {
            artifactWorkflowStateAbort(
                "artifact workflow is descriptor-pinned but no descriptor was supplied",
                "multischolar_artifact_state_descriptor_required"
            )
        }
        return(invisible(FALSE))
    }
    contract <- artifactWorkflowStateValidateDescriptorContract(contract)
    if (file.exists(path)) {
        pin <- artifactStoreReadJson(store, relative_path)
        version <- workflowStateVersionValue(pin$schema_version)
        valid <- identical(pin$schema, ARTIFACT_DESCRIPTOR_PIN_SCHEMA) &&
            identical(version, ARTIFACT_DESCRIPTOR_PIN_VERSION) &&
            identical(pin$project_id, identity$project_id) &&
            identical(pin$workflow_id, identity$workflow_id) &&
            identical(pin$contract, contract)
        if (!isTRUE(valid)) {
            artifactWorkflowStateAbort(
                "artifact workflow descriptor pin is incompatible with this session",
                "multischolar_artifact_state_descriptor_pin_mismatch"
            )
        }
        return(invisible(FALSE))
    }
    generations <- artifactStoreResolveFile(
        store, store$relative_paths$generations
    )
    recoverable_root <- !dir.exists(generations) ||
        length(list.files(generations, all.files = FALSE, no.. = TRUE)) == 0L
    if (!isTRUE(allow_create) && !isTRUE(recoverable_root)) {
        artifactWorkflowStateAbort(
            "existing artifact workflow has no immutable descriptor pin",
            "multischolar_unpinned_artifact_state_descriptor"
        )
    }
    pin <- list(
        schema = ARTIFACT_DESCRIPTOR_PIN_SCHEMA,
        schema_version = ARTIFACT_DESCRIPTOR_PIN_VERSION,
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        contract = contract,
        created_at = artifactRefUtcNow()
    )
    temporary_path <- artifactNormalizeRelativePath(paste0(
        relative_path,
        ".",
        artifactOpaqueId("tmp"),
        ".tmp"
    ))
    artifactStoreWriteJson(store, pin, temporary_path)
    on.exit({
        temporary <- artifactStoreResolveFile(store, temporary_path)
        if (file.exists(temporary)) unlink(temporary, force = FALSE)
    }, add = TRUE)
    artifactStorePublishFile(store, temporary_path, relative_path)
    invisible(TRUE)
}

artifactWorkflowStateEnsureMetadata <- function(store, identity, contract = NULL) {
    root_created <- artifactWorkflowStateEnsureRootManifest(store, identity)
    registry_identity <- artifactRegistryIdentity(
        store,
        identity,
        create_scope = root_created
    )
    artifactWorkflowStateEnsureDescriptorPin(
        store,
        identity,
        contract = contract,
        allow_create = root_created
    )
    registry_identity
}

artifactWorkflowStateManifestDigest <- function(manifest) {
    candidate <- manifest
    candidate$manifest_digest <- NULL
    artifactSemanticDigest(candidate)
}

artifactWorkflowStateValidateManifest <- function(manifest) {
    required <- c(
        "schema", "schema_version", "backend", "project_id", "workflow_id",
        "generation_id", "parent_generation_id", "logical_name", "created_at",
        "workflow_type", "audit_enabled", "description", "s4_class", "data",
        "config_json", "audit_json", "manifest_digest"
    )
    if (!is.list(manifest) || !identical(names(manifest), required) ||
        !identical(manifest$schema, ARTIFACT_WORKFLOW_STATE_SCHEMA)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState manifest has unknown or missing fields",
            "multischolar_malformed_artifact_state_manifest"
        )
    }
    version <- workflowStateVersionValue(manifest$schema_version)
    if (is.na(version)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState manifest version is missing",
            "multischolar_malformed_artifact_state_manifest"
        )
    }
    if (version > ARTIFACT_WORKFLOW_STATE_VERSION) {
        artifactWorkflowStateAbort(
            sprintf("unsupported future artifact WorkflowState version: %d", version),
            "multischolar_future_artifact_state_version"
        )
    }
    if (!identical(version, ARTIFACT_WORKFLOW_STATE_VERSION) ||
        !identical(manifest$backend, "artifact")) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState manifest version or backend is unsupported",
            "multischolar_unsupported_artifact_state_version"
        )
    }
    scalar_fields <- c(
        "project_id", "workflow_id", "generation_id", "logical_name",
        "created_at", "workflow_type", "manifest_digest"
    )
    if (!all(vapply(manifest[scalar_fields], workflowStateScalarString, logical(1))) ||
        !artifactRefValidUtc(manifest$created_at) ||
        !manifest$workflow_type %in% .WORKFLOW_STATE_TYPES ||
        !is.character(manifest$description) || length(manifest$description) != 1L ||
        is.na(manifest$description) ||
        !is.logical(manifest$audit_enabled) || length(manifest$audit_enabled) != 1L ||
        is.na(manifest$audit_enabled)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState manifest scalar metadata is invalid",
            "multischolar_malformed_artifact_state_manifest"
        )
    }
    artifactRefValidateId(manifest$generation_id, "generation_id", "gen")
    parent_id <- manifest$parent_generation_id
    if (!is.null(parent_id)) {
        artifactRefValidateId(parent_id, "parent_generation_id", "gen")
    }
    data <- manifest$data
    data_required <- c("semantic_digest", "codec", "metadata_json", "artifact_refs")
    if (!is.list(data) || !identical(names(data), data_required) ||
        !is.list(data$artifact_refs)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState data descriptor is malformed",
            "multischolar_malformed_artifact_state_manifest"
        )
    }
    if (is.null(data$semantic_digest)) {
        if (!is.null(data$codec) || !is.null(data$metadata_json) ||
            length(data$artifact_refs) > 0L || !is.null(manifest$s4_class)) {
            artifactWorkflowStateAbort(
                "empty artifact state contains payload metadata",
                "multischolar_malformed_artifact_state_manifest"
            )
        }
    } else {
        artifactRefValidateDigest(data$semantic_digest, "semantic_digest")
        codec_version <- if (is.list(data$codec)) {
            workflowStateVersionValue(data$codec$version)
        } else {
            NA_integer_
        }
        if (!is.list(data$codec) ||
            !workflowStateScalarString(data$codec$id) ||
            is.na(codec_version) || codec_version < 1L ||
            !workflowStateScalarString(data$metadata_json) ||
            !workflowStateScalarString(manifest$s4_class)) {
            artifactWorkflowStateAbort(
                "artifact state codec metadata is malformed",
                "multischolar_malformed_artifact_state_manifest"
            )
        }
        data$artifact_refs <- lapply(data$artifact_refs, artifactStoreNormalizeRef)
        names(data$artifact_refs) <- names(manifest$data$artifact_refs)
        manifest$data <- data
    }
    artifactWorkflowStateUnserializeMetadata(manifest$config_json, "config")
    artifactWorkflowStateUnserializeMetadata(manifest$audit_json, "audit metadata")
    artifactWorkflowStateAssertSafeMetadata(manifest, "manifest")
    expected_digest <- artifactWorkflowStateManifestDigest(manifest)
    if (!identical(manifest$manifest_digest, expected_digest)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState manifest digest does not match its contents",
            "multischolar_artifact_state_manifest_digest_mismatch"
        )
    }
    manifest
}

artifactWorkflowStateNewManifest <- function(
    identity,
    generation_id,
    parent_generation_id,
    logical_name,
    workflow_type,
    audit_enabled,
    description,
    data,
    config,
    audit_metadata
) {
    manifest <- list(
        schema = ARTIFACT_WORKFLOW_STATE_SCHEMA,
        schema_version = ARTIFACT_WORKFLOW_STATE_VERSION,
        backend = "artifact",
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        generation_id = generation_id,
        parent_generation_id = parent_generation_id,
        logical_name = logical_name,
        created_at = artifactRefUtcNow(),
        workflow_type = workflow_type,
        audit_enabled = isTRUE(audit_enabled),
        description = as.character(description),
        s4_class = data$s4_class,
        data = data[c("semantic_digest", "codec", "metadata_json", "artifact_refs")],
        config_json = artifactWorkflowStateSerializeMetadata(config, "config"),
        audit_json = artifactWorkflowStateSerializeMetadata(
            audit_metadata,
            "audit metadata"
        ),
        manifest_digest = NULL
    )
    manifest$manifest_digest <- artifactWorkflowStateManifestDigest(manifest)
    artifactWorkflowStateValidateManifest(manifest)
}

artifactWorkflowStateWriteManifest <- function(store, manifest) {
    manifest <- artifactWorkflowStateValidateManifest(manifest)
    relative_path <- artifactWorkflowStateManifestRelativePath(
        store,
        manifest$generation_id
    )
    temporary_path <- artifactNormalizeRelativePath(paste0(
        relative_path,
        ".",
        artifactOpaqueId("tmp"),
        ".tmp"
    ))
    artifactStoreWriteJson(store, unclass(manifest), temporary_path)
    on.exit({
        temporary <- artifactStoreResolveFile(store, temporary_path)
        if (file.exists(temporary)) unlink(temporary, force = FALSE)
    }, add = TRUE)
    artifactWorkflowStateReadManifest(store, temporary_path)
    artifactStorePublishFile(store, temporary_path, relative_path)
    relative_path
}

artifactWorkflowStateReadManifest <- function(store, relative_path) {
    manifest <- artifactStoreReadJson(store, relative_path)
    manifest$schema_version <- as.integer(manifest$schema_version)
    manifest$audit_enabled <- as.logical(manifest$audit_enabled)
    if (!is.null(manifest$data$codec)) {
        manifest$data$codec$version <- as.integer(manifest$data$codec$version)
    }
    if (is.list(manifest$data$artifact_refs) &&
        length(manifest$data$artifact_refs) == 0L) {
        manifest$data$artifact_refs <- list()
    }
    artifactWorkflowStateValidateManifest(manifest)
}

artifactWorkflowStateWriteData <- function(
    store,
    identity,
    generation_id,
    state_object,
    previous_manifest = NULL,
    parent_object = NULL,
    persistence_hint = NULL,
    failure_injector = NULL,
    dehydrate_fn = dehydrateDiaS4Artifact,
    validate_bundle_fn = validateDiaS4Bundle
) {
    if (is.null(state_object)) {
        return(list(
            s4_class = NULL,
            semantic_digest = NULL,
            codec = NULL,
            metadata_json = NULL,
            artifact_refs = list(),
            reused = FALSE,
            new_reference_names = character(),
            dependencies = list()
        ))
    }
    if (!is.function(dehydrate_fn) || !is.function(validate_bundle_fn)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState dehydrate adapter is invalid",
            "multischolar_invalid_artifact_state_codec"
        )
    }
    if (!is.null(persistence_hint) &&
        identical(persistence_hint$representation, "state_reference")) {
        if (is.null(previous_manifest) || is.null(parent_object)) {
            artifactWorkflowStateAbort(
                "artifact state reference requires one persisted parent object",
                "multischolar_invalid_artifact_state_reference"
            )
        }
        return(artifactWorkflowStateWriteStateReferenceData(
            previous_manifest,
            parent_object,
            state_object,
            persistence_hint
        ))
    }
    bundle <- artifactWorkflowStateDehydrate(dehydrate_fn, state_object, identity, generation_id)
    bundle <- validate_bundle_fn(bundle)
    semantic_digest <- bundle$metadata$semantic_digest
    if (!is.null(persistence_hint)) {
        if (is.null(previous_manifest) || is.null(parent_object)) {
            artifactWorkflowStateAbort(
                "artifact row selection requires one persisted parent object",
                "multischolar_invalid_artifact_row_selection"
            )
        }
        return(artifactWorkflowStateWritePersistenceData(
            store = store,
            identity = identity,
            generation_id = generation_id,
            state_object = state_object,
            parent_object = parent_object,
            previous_manifest = previous_manifest,
            bundle = bundle,
            hint = persistence_hint,
            failure_injector = failure_injector
        ))
    }
    if (!is.null(previous_manifest) &&
        identical(previous_manifest$data$semantic_digest, semantic_digest)) {
        return(list(
            s4_class = previous_manifest$s4_class,
            semantic_digest = semantic_digest,
            codec = previous_manifest$data$codec,
            metadata_json = previous_manifest$data$metadata_json,
            artifact_refs = previous_manifest$data$artifact_refs,
            reused = TRUE,
            new_reference_names = character(),
            dependencies = list()
        ))
    }
    references <- Map(function(payload, metadata, index) {
        encoded <- structure(
            list(payload = payload, metadata = metadata),
            class = c("MultiScholaRArtifactRectangular", "list")
        )
        logical_key <- list(
            project_id = identity$project_id,
            omic_type = identity$omic_type,
            workflow_slug = identity$workflow_slug,
            stage_id = paste0("state_", substr(generation_id, 5L, 20L)),
            state_role = sprintf("payload_%04d", index),
            generation_id = generation_id
        )
        artifactStoreWriteParquet(
            store,
            encoded,
            logical_key,
            failure_injector = failure_injector
        )
    }, bundle$payloads, bundle$metadata$payloads, seq_along(bundle$payloads))
    names(references) <- names(bundle$payloads)
    descriptor <- bundle$metadata$codec
    list(
        s4_class = bundle$metadata$class_name,
        semantic_digest = semantic_digest,
        codec = descriptor,
        metadata_json = artifactWorkflowStateSerializeMetadata(
            bundle$metadata,
            "S4 codec metadata"
        ),
        artifact_refs = references,
        reused = FALSE,
        new_reference_names = names(references),
        dependencies = list()
    )
}

artifactWorkflowStateHydrateData <- function(
    store,
    manifest,
    hydrate_fn = hydrateDiaS4Artifact,
    visited_generation_ids = character()
) {
    manifest <- artifactWorkflowStateValidateManifest(manifest)
    if (manifest$generation_id %in% visited_generation_ids) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState selection lineage contains a cycle",
            "multischolar_invalid_artifact_state_lineage"
        )
    }
    visited_generation_ids <- c(visited_generation_ids, manifest$generation_id)
    if (is.null(manifest$data$semantic_digest)) return(NULL)
    if (!is.function(hydrate_fn)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState hydrate adapter is invalid",
            "multischolar_invalid_artifact_state_codec"
        )
    }
    metadata <- artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "S4 codec metadata"
    )
    if (artifactWorkflowStateIsStateReference(metadata)) {
        return(artifactWorkflowStateHydrateStateReference(
            store,
            manifest,
            metadata,
            hydrate_fn,
            visited_generation_ids
        ))
    }
    references <- manifest$data$artifact_refs
    payload_entries <- lapply(names(references), function(payload_name) {
        ref <- references[[payload_name]]
        sidecar_path <- artifactStoreManagedPaths(
            store,
            ref$logical_key,
            ref$artifact_id
        )$sidecar
        sidecar <- artifactStoreReadSidecar(store, sidecar_path, validate_payload = TRUE)
        if (!identical(sidecar$artifact_ref, ref)) {
            artifactWorkflowStateAbort(
                "artifact state reference does not match its immutable sidecar",
                "multischolar_artifact_state_ref_mismatch"
            )
        }
        list(
            payload = arrow::read_parquet(
                artifactStoreResolveFile(store, ref$relative_path, must_exist = TRUE),
                as_data_frame = FALSE
            ),
            metadata = sidecar$codec_metadata,
            ref = ref
        )
    })
    names(payload_entries) <- names(references)
    payloads <- if (is.list(metadata$derivation)) {
        artifactWorkflowStateHydrateSelectionPayloads(
            store,
            manifest,
            metadata,
            payload_entries,
            hydrate_parent_fn = function(parent_manifest) {
                artifactWorkflowStateHydrateData(
                    store,
                    parent_manifest,
                    hydrate_fn,
                    visited_generation_ids
                )
            }
        )
    } else {
        expected <- names(metadata$payloads)
        if (!identical(names(payload_entries), expected)) {
            artifactWorkflowStateAbort(
                "artifact state refs do not match exact S4 payload metadata",
                "multischolar_artifact_state_ref_mismatch"
            )
        }
        lapply(payload_entries, `[[`, "payload")
    }
    bundle <- structure(
        list(metadata = metadata, payloads = payloads),
        class = c("MultiScholaRArtifactS4Bundle", "list")
    )
    hydrated <- hydrate_fn(bundle)
    if (!identical(class(hydrated)[[1L]], manifest$s4_class)) {
        artifactWorkflowStateAbort(
            "hydrated artifact state class differs from its manifest",
            "multischolar_artifact_state_class_mismatch"
        )
    }
    hydrated
}

artifactWorkflowStateVerifyHydration <- function(
    store,
    manifest,
    expected_object,
    hydrate_fn = hydrateDiaS4Artifact
) {
    hydrated <- artifactWorkflowStateHydrateData(store, manifest, hydrate_fn)
    if (!identical(hydrated, expected_object)) {
        expected_class <- if (is.null(expected_object)) NULL else class(expected_object)[[1L]]
        hydrated_class <- if (is.null(hydrated)) NULL else class(hydrated)[[1L]]
        artifactWorkflowStateAbort(
            "independent artifact hydration changed the pending scientific state",
            "multischolar_inexact_artifact_state_hydration",
            expected_class = expected_class,
            hydrated_class = hydrated_class
        )
    }
    if (isS4(hydrated) && !identical(methods::validObject(hydrated, test = TRUE), TRUE)) {
        artifactWorkflowStateAbort(
            "independently hydrated artifact state is not a valid S4 object",
            "multischolar_invalid_hydrated_s4_object"
        )
    }
    invisible(hydrated)
}

artifactWorkflowStateObjectDigest <- function(value) {
    digest::digest(
        serialize(value, NULL, version = 3L),
        algo = "sha256",
        serialize = FALSE
    )
}

artifactWorkflowStateVerifyHydrationInline <- function(
    store,
    manifest,
    expected_object,
    hydrate_fn
) {
    hydrated <- artifactWorkflowStateVerifyHydration(
        store,
        manifest,
        expected_object,
        hydrate_fn
    )
    expected_digest <- artifactWorkflowStateObjectDigest(expected_object)
    proof <- list(
        valid = TRUE,
        mode = "inline_exact",
        verifier_pid = as.integer(Sys.getpid()),
        expected_digest = expected_digest,
        hydrated_digest = expected_digest,
        manifest_semantic_digest = manifest$data$semantic_digest,
        generation_id = manifest$generation_id,
        complete_payload_returned = FALSE
    )
    hydrated <- NULL
    artifactResourceDataOnly(proof, "artifact hydration verification proof")
    proof
}

artifactWorkflowStateArtifactRecord <- function(identity, ref, ordinal) {
    list(
        workflow_id = identity$workflow_id,
        artifact_id = ref$artifact_id,
        run_id = NULL,
        generation_id = ref$logical_key$generation_id,
        stage_id = ref$logical_key$stage_id,
        state_role = ref$logical_key$state_role,
        hydration_ordinal = as.integer(ordinal),
        relative_path = ref$relative_path,
        codec_id = ref$codec$id,
        codec_version = ref$codec$version,
        payload_schema_id = ref$payload_schema$id,
        payload_schema_version = ref$payload_schema$version,
        semantic_digest = ref$hash_policy$semantic$digest,
        byte_digest = ref$hash_policy$byte$digest,
        row_count = ref$shape$rows,
        column_count = ref$shape$columns,
        payload_bytes = ref$shape$bytes,
        status = "committed",
        created_at = ref$created_at,
        updated_at = ref$updated_at
    )
}

artifactWorkflowStateJson <- function(value) {
    as.character(jsonlite::toJSON(
        value,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA
    ))
}

artifactWorkflowStateEnsureWorkflow <- function(session, identity) {
    existing <- projectRegistryQuery(
        session,
        "workflows",
        filters = list(workflow_id = identity$workflow_id),
        limit = 1L
    )
    if (nrow(existing) == 1L) {
        same_workflow <- identical(existing$omic_type[[1L]], identity$omic_type) &&
            identical(existing$omic_label[[1L]], identity$omic_label) &&
            identical(existing$workflow_slug[[1L]], identity$workflow_slug)
        if (!isTRUE(same_workflow)) {
            artifactWorkflowStateAbort(
                "artifact registry workflow ID belongs to another workflow",
                "multischolar_artifact_registry_workflow_collision"
            )
        }
        return(invisible(FALSE))
    }
    now <- artifactRefUtcNow()
    projectRegistryWrite(session, "workflow", list(
        workflow_id = identity$workflow_id,
        omic_type = identity$omic_type,
        omic_label = identity$omic_label,
        workflow_slug = identity$workflow_slug,
        status = "active",
        created_at = now,
        updated_at = now
    ))
    invisible(TRUE)
}

artifactWorkflowStateUpdateStatus <- function(
    session,
    identity,
    generation_id,
    from_status,
    to_status,
    timestamp
) {
    connection <- projectRegistrySessionConnection(session, write = TRUE)
    affected <- projectRegistryExecuteBound(
        connection,
        paste0(
            "UPDATE workflow_states SET status = ?, updated_at = ? ",
            "WHERE project_id = ? AND workflow_id = ? AND generation_id = ? ",
            "AND status = ?"
        ),
        list(
            to_status,
            timestamp,
            identity$project_id,
            identity$workflow_id,
            generation_id,
            from_status
        )
    )
    if (!identical(as.integer(affected), 1L)) {
        artifactWorkflowStateAbort(
            "artifact state current-pointer compare-and-set failed",
            "multischolar_artifact_state_compare_and_set_failed"
        )
    }
    invisible(TRUE)
}

artifactWorkflowStateTransaction <- function(session, operation) {
    connection <- projectRegistrySessionConnection(session, write = TRUE)
    DBI::dbBegin(connection)
    committed <- FALSE
    on.exit({
        if (!committed) try(DBI::dbRollback(connection), silent = TRUE)
    }, add = TRUE)
    value <- operation()
    DBI::dbCommit(connection)
    committed <- TRUE
    value
}

artifactWorkflowStateEvents <- function(session, workflow_id) {
    rows <- projectRegistryQuery(
        session,
        "events",
        filters = list(workflow_id = workflow_id)
    )
    if (nrow(rows) == 0L) return(list())
    lapply(seq_len(nrow(rows)), function(index) {
        details <- jsonlite::fromJSON(
            rows$details_json[[index]],
            simplifyVector = FALSE
        )
        list(
            revision = as.integer(index),
            event_type = rows$event_type[[index]],
            state_name = details$state_name,
            previous_state = details$previous_state,
            timestamp = as.POSIXct(rows$recorded_at[[index]], tz = "UTC"),
            metadata = details$metadata
        )
    })
}

#' Release the artifact state hydration cache
#'
#' @param workflow_state A workflow state manager.
#' @return `TRUE` when one cached object was released, otherwise `FALSE`.
#' @noRd
workflowStateReleaseHydrationCache <- function(workflow_state) {
    if (!inherits(workflow_state, "ArtifactWorkflowState") ||
        !is.function(workflow_state$releaseCache)) {
        return(invisible(FALSE))
    }
    workflow_state$releaseCache()
}

#' Release an idle artifact state registry session
#'
#' @param workflow_state A workflow state manager.
#' @return `TRUE` when an open registry session was released.
#' @noRd
workflowStateReleaseRegistrySession <- function(workflow_state) {
    if (!inherits(workflow_state, "ArtifactWorkflowState") ||
        !is.function(workflow_state$releaseRegistrySession)) {
        return(invisible(FALSE))
    }
    workflow_state$releaseRegistrySession()
}

artifactWorkflowStateExportSnapshot <- function(
    identity,
    current_generation_id,
    current_state,
    lineage,
    workflow_type,
    audit_enabled
) {
    list(
        schema = ARTIFACT_WORKFLOW_STATE_SCHEMA,
        schema_version = ARTIFACT_WORKFLOW_STATE_VERSION,
        backend = "artifact",
        project_id = identity$project_id,
        workflow_id = identity$workflow_id,
        current_generation_id = current_generation_id,
        current_state = current_state,
        active_lineage = lapply(lineage, function(row) {
            list(
                generation_id = row$generation_id,
                logical_name = row$logical_name,
                manifest_relative_path = row$manifest_relative_path
            )
        }),
        workflow_type = workflow_type,
        audit_enabled = isTRUE(audit_enabled)
    )
}

artifactWorkflowStateValidateRestoreSnapshot <- function(
    manifest,
    identity,
    current_generation_id,
    schema_version = NULL
) {
    if (!is.null(schema_version) || !is.list(manifest)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState cannot import a memory session manifest",
            "multischolar_unsupported_artifact_state_restore"
        )
    }
    version <- workflowStateVersionValue(manifest$schema_version)
    if (is.na(version)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState restore version is missing",
            "multischolar_unsupported_artifact_state_restore"
        )
    }
    if (version > ARTIFACT_WORKFLOW_STATE_VERSION) {
        artifactWorkflowStateAbort(
            sprintf("unsupported future artifact WorkflowState version: %d", version),
            "multischolar_future_artifact_state_version"
        )
    }
    required <- c(
        "schema", "schema_version", "backend", "project_id", "workflow_id",
        "current_generation_id", "current_state", "active_lineage",
        "workflow_type", "audit_enabled"
    )
    valid <- identical(names(manifest), required) &&
        identical(manifest$schema, ARTIFACT_WORKFLOW_STATE_SCHEMA) &&
        identical(version, ARTIFACT_WORKFLOW_STATE_VERSION) &&
        identical(manifest$backend, "artifact") &&
        identical(manifest$project_id, identity$project_id) &&
        identical(manifest$workflow_id, identity$workflow_id) &&
        identical(manifest$current_generation_id, current_generation_id)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState restore manifest is incompatible",
            "multischolar_unsupported_artifact_state_restore"
        )
    }
    invisible(TRUE)
}
