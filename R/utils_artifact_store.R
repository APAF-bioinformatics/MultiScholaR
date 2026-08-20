# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactStoreSchema <- "multischolar.artifact_store"
.artifactStoreSchemaVersion <- 1L
.artifactStoreProtocolVersion <- 1L
.artifactSidecarSchema <- "multischolar.artifact_sidecar"
.artifactSidecarSchemaVersion <- 1L
.artifactIntentSchema <- "multischolar.artifact_write_intent"
.artifactIntentSchemaVersion <- 1L

artifactStoreAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_store_error"),
        ...
    )
}

newArtifactStore <- function(paths, project_id) {
    if (!inherits(paths, "MultiScholaRArtifactPaths") ||
        !workflowCapabilityScalarString(project_id)) {
        artifactStoreAbort(
            "artifact store requires validated paths and one project ID",
            "multischolar_invalid_artifact_store"
        )
    }
    structure(
        list(
            schema = .artifactStoreSchema,
            schema_version = .artifactStoreSchemaVersion,
            project_id = project_id,
            project_root = paths$project_root,
            labels = paths$labels,
            relative_paths = paths$relative_paths
        ),
        class = c("MultiScholaRArtifactStore", "list")
    )
}

validateArtifactStore <- function(store) {
    required <- c(
        "schema", "schema_version", "project_id", "project_root", "labels",
        "relative_paths"
    )
    valid <- inherits(store, "MultiScholaRArtifactStore") &&
        identical(names(store), required) &&
        identical(store$schema, .artifactStoreSchema) &&
        identical(store$schema_version, .artifactStoreSchemaVersion) &&
        workflowCapabilityScalarString(store$project_id) &&
        identical(artifactNormalizeProjectRoot(store$project_root), store$project_root)
    if (!isTRUE(valid)) {
        artifactStoreAbort(
            "artifact store structure or version is invalid",
            "multischolar_invalid_artifact_store"
        )
    }
    store
}

artifactStoreValidateLogicalKey <- function(store, logical_key) {
    store <- validateArtifactStore(store)
    artifactRefValidateLogicalKey(logical_key)
    expected <- c(
        project_id = store$project_id,
        omic_type = store$labels$omic_type,
        workflow_slug = store$labels$workflow_slug
    )
    actual <- unlist(logical_key[names(expected)], use.names = TRUE)
    if (!identical(actual, expected)) {
        artifactStoreAbort(
            "artifact logical key belongs to another project or workflow",
            "multischolar_cross_project_artifact"
        )
    }
    path_fields <- c("stage_id", "state_role", "generation_id")
    for (field in path_fields) {
        artifactValidatePathComponent(logical_key[[field]], field)
    }
    invisible(logical_key)
}

artifactStoreManagedPaths <- function(store, logical_key, artifact_id) {
    artifactStoreValidateLogicalKey(store, logical_key)
    artifactRefValidateId(artifact_id, "artifact_id", "art")
    payload_directory <- file.path(
        store$relative_paths$tables,
        logical_key$stage_id,
        logical_key$state_role,
        logical_key$generation_id,
        artifact_id
    )
    list(
        payload_directory = artifactNormalizeRelativePath(payload_directory),
        payload = artifactNormalizeRelativePath(
            file.path(payload_directory, "payload.parquet")
        ),
        sidecar = artifactNormalizeRelativePath(file.path(
            store$relative_paths$generations,
            logical_key$generation_id,
            paste0(artifact_id, ".artifact.json")
        )),
        intent = artifactNormalizeRelativePath(file.path(
            store$relative_paths$staging,
            paste0(artifact_id, ".intent.json")
        ))
    )
}

artifactStoreEnsureDirectory <- function(store, relative_path) {
    store <- validateArtifactStore(store)
    relative_path <- artifactNormalizeRelativePath(relative_path)
    path <- artifactResolveContainedPath(store$project_root, relative_path)
    if (file.exists(path) && !dir.exists(path)) {
        artifactStoreAbort(
            sprintf("artifact directory path is occupied: %s", relative_path),
            "multischolar_artifact_path_conflict"
        )
    }
    if (!dir.exists(path)) {
        created <- dir.create(path, recursive = TRUE, mode = "0700", showWarnings = FALSE)
        if (!isTRUE(created) && !dir.exists(path)) {
            artifactStoreAbort(
                sprintf("artifact directory could not be created: %s", relative_path),
                "multischolar_artifact_write_failed"
            )
        }
    }
    artifactResolveContainedPath(store$project_root, relative_path, must_exist = TRUE)
}

artifactStoreResolveFile <- function(store, relative_path, must_exist = FALSE) {
    validateArtifactStore(store)
    artifactResolveContainedPath(
        store$project_root,
        artifactNormalizeRelativePath(relative_path),
        must_exist = must_exist
    )
}

artifactStoreWriteJson <- function(store, value, relative_path) {
    artifactStoreEnsureDirectory(
        store,
        artifactNormalizeRelativePath(dirname(relative_path))
    )
    target <- artifactStoreResolveFile(store, relative_path)
    if (file.exists(target) || dir.exists(target)) {
        artifactStoreAbort(
            sprintf("artifact temporary path already exists: %s", relative_path),
            "multischolar_artifact_path_conflict"
        )
    }
    encoded <- jsonlite::toJSON(
        value,
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = NA,
        pretty = TRUE
    )
    connection <- file(target, open = "wb")
    on.exit(close(connection), add = TRUE)
    tryCatch(
        writeBin(charToRaw(as.character(encoded)), connection),
        error = \(error) artifactStoreAbort(
            sprintf("artifact metadata write failed: %s", conditionMessage(error)),
            "multischolar_artifact_write_failed",
            parent = error
        )
    )
    invisible(target)
}

artifactStorePublishFile <- function(store, temporary_path, final_path) {
    source <- artifactStoreResolveFile(store, temporary_path, must_exist = TRUE)
    target <- artifactStoreResolveFile(store, final_path)
    artifactStoreEnsureDirectory(
        store,
        artifactNormalizeRelativePath(dirname(final_path))
    )
    if (file.exists(target) || dir.exists(target)) {
        artifactStoreAbort(
            sprintf("immutable artifact already exists: %s", final_path),
            "multischolar_artifact_already_exists"
        )
    }
    published <- suppressWarnings(file.link(source, target))
    if (!isTRUE(published)) {
        conflict <- file.exists(target) || dir.exists(target)
        artifactStoreAbort(
            sprintf("artifact finalization failed: %s", final_path),
            if (conflict) {
                "multischolar_artifact_already_exists"
            } else {
                "multischolar_artifact_publish_failed"
            }
        )
    }
    unlink(source, force = FALSE)
    invisible(target)
}

artifactStorePublishDirectory <- function(store, temporary_path, final_path) {
    source <- artifactStoreResolveFile(store, temporary_path, must_exist = TRUE)
    target <- artifactStoreResolveFile(store, final_path)
    artifactStoreEnsureDirectory(
        store,
        artifactNormalizeRelativePath(dirname(final_path))
    )
    if (!dir.exists(source) || file.exists(target) || dir.exists(target)) {
        artifactStoreAbort(
            sprintf("immutable artifact directory cannot be finalized: %s", final_path),
            "multischolar_artifact_already_exists"
        )
    }
    published <- suppressWarnings(file.rename(source, target))
    if (!isTRUE(published)) {
        artifactStoreAbort(
            sprintf("artifact directory finalization failed: %s", final_path),
            "multischolar_artifact_publish_failed"
        )
    }
    invisible(target)
}

artifactStorePayloadShape <- function(path, metadata) {
    payload <- arrow::read_parquet(path, as_data_frame = FALSE)
    decodeArtifactRectangular(payload, metadata)
    list(
        kind = metadata$kind,
        rows = as.integer(metadata$dimensions$rows),
        columns = as.integer(metadata$dimensions$columns),
        payloads = 1L,
        bytes = unname(as.numeric(file.info(path)$size))
    )
}

artifactStoreNewSidecar <- function(ref, codec_metadata, operation_id) {
    ref <- validateArtifactRef(ref)
    list(
        schema = .artifactSidecarSchema,
        schema_version = .artifactSidecarSchemaVersion,
        protocol_version = .artifactStoreProtocolVersion,
        artifact_ref = unclass(ref),
        codec_metadata = codec_metadata,
        registration = list(
            artifact_id = ref$artifact_id,
            logical_key = ref$logical_key,
            provenance_ids = ref$provenance_ids,
            status = ref$status
        ),
        operation_id = operation_id,
        guarantees = list(
            visibility = "same_filesystem_atomic_publication",
            durability = "power_loss_durability_not_claimed_without_fsync",
            overwrite = "forbidden"
        ),
        created_at = artifactRefUtcNow()
    )
}

artifactStoreNewIntent <- function(ref, codec_metadata, managed_paths, temporary_paths,
                                   operation_id) {
    list(
        schema = .artifactIntentSchema,
        schema_version = .artifactIntentSchemaVersion,
        protocol_version = .artifactStoreProtocolVersion,
        operation_id = operation_id,
        artifact_ref = unclass(ref),
        codec_metadata = codec_metadata,
        managed_paths = managed_paths,
        temporary_paths = temporary_paths,
        guarantees = list(
            visibility = "same_filesystem_atomic_publication",
            durability = "power_loss_durability_not_claimed_without_fsync",
            overwrite = "forbidden"
        ),
        created_at = artifactRefUtcNow()
    )
}

artifactStoreInvokeFailure <- function(failure_injector, stage, context) {
    if (is.null(failure_injector)) return(invisible(NULL))
    if (!is.function(failure_injector)) {
        artifactStoreAbort(
            "artifact failure injector must be a function",
            "multischolar_invalid_artifact_failure_injector"
        )
    }
    failure_injector(stage, context)
    invisible(NULL)
}

artifactStoreWriteParquet <- function(store, encoded, logical_key,
                                      provenance_ids = character(),
                                      failure_injector = NULL,
                                      write_parquet_fn = arrow::write_parquet) {
    store <- validateArtifactStore(store)
    if (!inherits(encoded, "MultiScholaRArtifactRectangular") ||
        !is.function(write_parquet_fn)) {
        artifactStoreAbort(
            "artifact store requires one encoded rectangular payload",
            "multischolar_invalid_artifact_store_payload"
        )
    }
    metadata <- artifactValidateRectangularMetadata(encoded$metadata)
    artifactStoreValidateLogicalKey(store, logical_key)
    artifactStoreAssertLogicalKeyAvailable(store, logical_key)
    artifact_id <- artifactOpaqueId("art")
    operation_id <- artifactOpaqueId("op")
    managed_paths <- artifactStoreManagedPaths(store, logical_key, artifact_id)
    temp_token <- artifactOpaqueId("tmp")
    payload_parent <- artifactNormalizeRelativePath(
        dirname(managed_paths$payload_directory)
    )
    temporary_paths <- list(
        payload_directory = artifactNormalizeRelativePath(file.path(
            payload_parent,
            paste0(".", artifact_id, ".", temp_token, ".tmp")
        )),
        payload = NULL,
        sidecar = artifactNormalizeRelativePath(paste0(
            managed_paths$sidecar,
            ".",
            temp_token,
            ".tmp"
        )),
        intent = artifactNormalizeRelativePath(paste0(
            managed_paths$intent,
            ".",
            temp_token,
            ".tmp"
        ))
    )
    temporary_paths$payload <- artifactNormalizeRelativePath(file.path(
        temporary_paths$payload_directory,
        "payload.parquet"
    ))
    shape <- list(
        kind = metadata$kind,
        rows = as.integer(metadata$dimensions$rows),
        columns = as.integer(metadata$dimensions$columns),
        payloads = 1L,
        bytes = NA_real_
    )
    staged_ref <- newArtifactRef(
        logical_key = logical_key,
        relative_path = managed_paths$payload,
        codec_id = metadata$codec$id,
        payload_schema_id = metadata$payload_schema$id,
        shape = shape,
        semantic_digest = metadata$semantic_digest,
        provenance_ids = provenance_ids,
        status = "staged",
        artifact_id = artifact_id
    )
    intent <- artifactStoreNewIntent(
        staged_ref,
        metadata,
        managed_paths,
        temporary_paths,
        operation_id
    )
    artifactStoreWriteJson(store, intent, temporary_paths$intent)
    artifactStorePublishFile(store, temporary_paths$intent, managed_paths$intent)
    artifactStoreInvokeFailure(failure_injector, "before_write", intent)

    artifactStoreEnsureDirectory(store, temporary_paths$payload_directory)
    payload_path <- artifactStoreResolveFile(store, temporary_paths$payload)
    tryCatch(
        do.call(
            write_parquet_fn,
            c(
                list(x = encoded$payload, sink = payload_path),
                artifactParquetWriteArgs(encoded)
            )
        ),
        error = \(error) artifactStoreAbort(
            sprintf("artifact payload write failed: %s", conditionMessage(error)),
            "multischolar_artifact_write_failed",
            parent = error
        )
    )
    artifactStoreInvokeFailure(failure_injector, "after_temp_write", intent)

    actual_shape <- artifactStorePayloadShape(payload_path, metadata)
    before <- decodeArtifactRectangular(encoded$payload, metadata)
    after <- decodeArtifactRectangular(
        arrow::read_parquet(payload_path, as_data_frame = FALSE),
        metadata
    )
    if (!identical(after, before)) {
        artifactStoreAbort(
            "artifact payload changed during Parquet serialization",
            "multischolar_artifact_validation_failed"
        )
    }
    committed_ref <- newArtifactRef(
        logical_key = logical_key,
        relative_path = managed_paths$payload,
        codec_id = metadata$codec$id,
        payload_schema_id = metadata$payload_schema$id,
        shape = actual_shape,
        semantic_digest = metadata$semantic_digest,
        byte_digest = artifactByteDigest(payload_path),
        provenance_ids = provenance_ids,
        status = "committed",
        artifact_id = artifact_id,
        created_at = staged_ref$created_at,
        updated_at = artifactRefUtcNow()
    )
    sidecar <- artifactStoreNewSidecar(committed_ref, metadata, operation_id)
    artifactStoreWriteJson(store, sidecar, temporary_paths$sidecar)
    artifactStoreReadSidecar(store, temporary_paths$sidecar, validate_payload = FALSE)
    artifactStoreInvokeFailure(failure_injector, "after_validation", sidecar)

    artifactStorePublishDirectory(
        store,
        temporary_paths$payload_directory,
        managed_paths$payload_directory
    )
    validateArtifactRefPayload(committed_ref, store$project_root, actual_shape)
    artifactStoreInvokeFailure(failure_injector, "after_payload_rename", sidecar)
    artifactStorePublishFile(store, temporary_paths$sidecar, managed_paths$sidecar)
    artifactStoreInvokeFailure(failure_injector, "after_sidecar_rename", sidecar)
    unlink(artifactStoreResolveFile(store, managed_paths$intent), force = FALSE)
    committed_ref
}
