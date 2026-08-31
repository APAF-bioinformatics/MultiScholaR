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

artifactStoreValidateArrowShape <- function(payload, metadata) {
    fields <- payload$schema$fields
    actual_schema <- lapply(fields, \(field) {
        list(
            name = field$name,
            type = field$type$ToString(),
            nullable = field$nullable
        )
    })
    expected_rows <- as.integer(metadata$dimensions$rows)
    expected_columns <- length(metadata$physical_schema)
    valid <- identical(actual_schema, metadata$physical_schema) &&
        identical(as.integer(payload$num_rows), expected_rows) &&
        identical(as.integer(payload$num_columns), as.integer(expected_columns))
    if (!isTRUE(valid)) {
        artifactStoreAbort(
            "artifact Parquet shape differs from its codec metadata",
            "multischolar_artifact_shape_mismatch"
        )
    }
    row_order_index <- match(
        .artifactRowOrderColumn,
        vapply(fields, `[[`, character(1), "name")
    )
    row_order <- as.numeric(payload$column(row_order_index - 1L)$as_vector())
    if (!identical(row_order, as.numeric(seq_len(expected_rows)))) {
        artifactStoreAbort(
            "artifact Parquet row order is missing, duplicated, or reordered",
            "multischolar_artifact_order_mismatch"
        )
    }
    invisible(TRUE)
}

artifactStorePayloadShape <- function(path, metadata, validate_hydration = TRUE) {
    metadata <- artifactValidateRectangularMetadata(metadata)
    payload <- arrow::read_parquet(path, as_data_frame = FALSE)
    artifactStoreValidateArrowShape(payload, metadata)
    if (isTRUE(validate_hydration)) decodeArtifactRectangular(payload, metadata)
    list(
        kind = metadata$kind,
        rows = as.integer(metadata$dimensions$rows),
        columns = as.integer(metadata$dimensions$columns),
        payloads = 1L,
        bytes = unname(as.numeric(file.info(path)$size))
    )
}

artifactExactDigestValue <- function(value) {
    if (is.null(value)) return(NULL)
    if (is.atomic(value)) {
        value_attributes <- attributes(value)
        attributes(value) <- NULL
        materialized <- vector(typeof(value), length(value))
        materialized[] <- value
        return(list(
            storage_type = typeof(materialized),
            values = materialized,
            attributes = artifactExactDigestValue(value_attributes)
        ))
    }
    if (is.list(value)) {
        value_names <- names(value)
        return(list(
            names = if (is.null(value_names)) NULL else unname(value_names),
            values = unname(lapply(value, artifactExactDigestValue))
        ))
    }
    artifactStoreAbort(
        "exact hydration digest contains an unsupported value",
        "multischolar_invalid_artifact_hydration_digest"
    )
}

artifactExactHydrationDigest <- function(value) {
    if (!is.data.frame(value) && !is.matrix(value)) {
        artifactStoreAbort(
            "exact hydration digest requires a data frame or matrix",
            "multischolar_invalid_artifact_hydration_digest"
        )
    }
    descriptor <- if (is.data.frame(value)) {
        list(
            kind = "data.frame",
            class = unname(class(value)),
            dimensions = unname(dim(value)),
            names = unname(names(value)),
            row_names = artifactRowNamesMetadata(value),
            columns = unname(lapply(value, artifactExactDigestValue))
        )
    } else {
        list(
            kind = "matrix",
            class = unname(class(value)),
            storage_mode = storage.mode(value),
            dimensions = unname(dim(value)),
            dimnames = artifactExactDigestValue(dimnames(value)),
            values = artifactExactDigestValue(as.vector(value))
        )
    }
    artifactSemanticDigest(list(
        schema = "multischolar.exact_hydration_digest",
        schema_version = 1L,
        descriptor = descriptor
    ))
}

artifactStoreVerifyExactRef <- function(store, ref) {
    store <- validateArtifactStore(store)
    ref <- artifactStoreNormalizeRef(ref)
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(
        store,
        managed$sidecar,
        validate_payload = FALSE
    )
    if (!identical(artifactStoreNormalizeRef(sidecar$artifact_ref), ref)) {
        artifactStoreAbort(
            "artifact verification reference differs from its immutable sidecar",
            "multischolar_artifact_ref_mismatch"
        )
    }
    payload_path <- artifactStoreResolveFile(
        store,
        ref$relative_path,
        must_exist = TRUE
    )
    if (!identical(artifactByteDigest(payload_path), ref$hash_policy$byte$digest)) {
        artifactStoreAbort(
            "artifact verification byte digest differs from its reference",
            "multischolar_artifact_byte_digest_mismatch"
        )
    }
    payload <- arrow::read_parquet(payload_path, as_data_frame = FALSE)
    artifactStoreValidateArrowShape(payload, sidecar$codec_metadata)
    shape <- list(
        kind = sidecar$codec_metadata$kind,
        rows = as.integer(sidecar$codec_metadata$dimensions$rows),
        columns = as.integer(sidecar$codec_metadata$dimensions$columns),
        payloads = 1L,
        bytes = unname(as.numeric(file.info(payload_path)$size))
    )
    validateArtifactRefPayload(ref, store$project_root, shape)
    streaming <- identical(
        sidecar$codec_metadata$codec,
        list(
            id = .artifactStreamingCodec,
            version = .artifactStreamingCodecVersion
        )
    )
    if (isTRUE(streaming)) {
        verifyArtifactStreamingPayload(payload, sidecar$codec_metadata)
    }
    hydrated <- decodeArtifactRectangular(payload, sidecar$codec_metadata)
    payload <- NULL
    if (!isTRUE(streaming)) {
        stable_key <- sidecar$codec_metadata$stable_key$logical_columns
        if (length(stable_key) == 0L) stable_key <- NULL
        reencoded <- if (identical(sidecar$codec_metadata$kind, "matrix")) {
            encodeArtifactMatrix(hydrated, sidecar$codec_metadata$owner)
        } else {
            encodeArtifactTable(
                hydrated,
                stable_key = stable_key,
                owner = sidecar$codec_metadata$owner
            )
        }
        reencoded_metadata <- artifactStoreNormalizeCodecMetadata(reencoded$metadata)
        metadata_fields <- setdiff(names(reencoded_metadata), "semantic_digest")
        if (!identical(
            reencoded_metadata[metadata_fields],
            sidecar$codec_metadata[metadata_fields]
        )) {
            artifactStoreAbort(
                "artifact verification codec metadata is not an exact fixed point",
                "multischolar_artifact_validation_failed"
            )
        }
    }
    hydrated_shape <- c(rows = nrow(hydrated), columns = ncol(hydrated))
    expected_shape <- c(rows = ref$shape$rows, columns = ref$shape$columns)
    if (!identical(as.integer(hydrated_shape), as.integer(expected_shape))) {
        artifactStoreAbort(
            "artifact verification hydration shape differs from its reference",
            "multischolar_artifact_shape_mismatch"
        )
    }
    list(
        schema = "multischolar.artifact_exact_verification",
        schema_version = 1L,
        artifact_id = ref$artifact_id,
        semantic_digest = ref$hash_policy$semantic$digest,
        byte_digest = ref$hash_policy$byte$digest,
        hydration_digest = artifactExactHydrationDigest(hydrated),
        rows = as.integer(ref$shape$rows),
        columns = as.integer(ref$shape$columns),
        verifier_pid = as.integer(Sys.getpid())
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
                                      write_parquet_fn = arrow::write_parquet,
                                      verification = c("inline_exact", "deferred_exact")) {
    store <- validateArtifactStore(store)
    verification <- match.arg(verification)
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

    actual_shape <- artifactStorePayloadShape(
        payload_path,
        metadata,
        validate_hydration = FALSE
    )
    if (identical(verification, "inline_exact")) {
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
