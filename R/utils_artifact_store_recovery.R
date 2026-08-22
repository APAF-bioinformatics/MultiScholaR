# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactStoreReadJson <- function(store, relative_path) {
    path <- artifactStoreResolveFile(store, relative_path, must_exist = TRUE)
    tryCatch(
        jsonlite::fromJSON(
            path,
            simplifyVector = TRUE,
            simplifyDataFrame = FALSE,
            simplifyMatrix = FALSE
        ),
        error = \(error) artifactStoreAbort(
            sprintf("artifact metadata is not valid JSON: %s", relative_path),
            "multischolar_malformed_artifact_metadata",
            parent = error
        )
    )
}

artifactStoreCharacterVector <- function(value) {
    if (is.null(value) || (is.list(value) && length(value) == 0L)) {
        return(character())
    }
    unname(as.character(unlist(value, recursive = FALSE, use.names = FALSE)))
}

artifactStoreNormalizeCodecMetadata <- function(metadata) {
    metadata <- artifactValidateRectangularMetadata(metadata)
    metadata$stable_key$logical_columns <- artifactStoreCharacterVector(
        metadata$stable_key$logical_columns
    )
    metadata$stable_key$physical_columns <- artifactStoreCharacterVector(
        metadata$stable_key$physical_columns
    )
    metadata$schema_evolution$allowed_changes <- artifactStoreCharacterVector(
        metadata$schema_evolution$allowed_changes
    )
    metadata$writer_settings$dictionary_columns <- artifactStoreCharacterVector(
        metadata$writer_settings$dictionary_columns
    )
    numeric_settings <- c(
        "data_page_size", "target_file_size_bytes", "minimum_file_size_bytes",
        "maximum_file_size_bytes"
    )
    for (field in numeric_settings) {
        metadata$writer_settings[[field]] <- as.numeric(
            metadata$writer_settings[[field]]
        )
    }
    metadata
}

artifactStoreNormalizeRef <- function(ref) {
    if (!is.list(ref)) {
        artifactStoreAbort(
            "artifact metadata does not contain a reference",
            "multischolar_malformed_artifact_sidecar"
        )
    }
    ref$schema_version <- as.integer(ref$schema_version)
    ref$codec$version <- as.integer(ref$codec$version)
    ref$payload_schema$version <- as.integer(ref$payload_schema$version)
    ref$shape$rows <- as.integer(ref$shape$rows)
    ref$shape$columns <- as.integer(ref$shape$columns)
    ref$shape$payloads <- as.integer(ref$shape$payloads)
    ref$shape$bytes <- if (is.null(ref$shape$bytes)) {
        NA_real_
    } else {
        as.numeric(ref$shape$bytes)
    }
    if (is.null(ref$hash_policy$byte$digest)) {
        ref$hash_policy$byte$digest <- NA_character_
    }
    ref$hash_policy$semantic$version <- as.integer(
        ref$hash_policy$semantic$version
    )
    ref$provenance_ids <- artifactStoreCharacterVector(ref$provenance_ids)
    structure(
        validateArtifactRef(ref),
        class = c("MultiScholaRArtifactRef", "list")
    )
}

artifactStoreValidateGuarantees <- function(guarantees) {
    expected <- list(
        visibility = "same_filesystem_atomic_publication",
        durability = "power_loss_durability_not_claimed_without_fsync",
        overwrite = "forbidden"
    )
    if (!identical(guarantees, expected)) {
        artifactStoreAbort(
            "artifact metadata has unsupported storage guarantees",
            "multischolar_unsupported_artifact_store_protocol"
        )
    }
    invisible(guarantees)
}

artifactStoreValidateSidecar <- function(sidecar) {
    required <- c(
        "schema", "schema_version", "protocol_version", "artifact_ref",
        "codec_metadata", "registration", "operation_id", "guarantees",
        "created_at"
    )
    if (!is.list(sidecar) || !identical(names(sidecar), required) ||
        !identical(sidecar$schema, .artifactSidecarSchema)) {
        artifactStoreAbort(
            "artifact sidecar has unknown or missing fields",
            "multischolar_malformed_artifact_sidecar"
        )
    }
    schema_version <- as.integer(sidecar$schema_version)
    protocol_version <- as.integer(sidecar$protocol_version)
    if (is.na(schema_version) || is.na(protocol_version) ||
        schema_version > .artifactSidecarSchemaVersion ||
        protocol_version > .artifactStoreProtocolVersion) {
        artifactStoreAbort(
            "artifact sidecar uses a future schema or protocol",
            "multischolar_unsupported_artifact_sidecar_version"
        )
    }
    if (!identical(schema_version, .artifactSidecarSchemaVersion) ||
        !identical(protocol_version, .artifactStoreProtocolVersion)) {
        artifactStoreAbort(
            "artifact sidecar version is unsupported",
            "multischolar_unsupported_artifact_sidecar_version"
        )
    }
    sidecar$schema_version <- schema_version
    sidecar$protocol_version <- protocol_version
    sidecar$artifact_ref <- artifactStoreNormalizeRef(sidecar$artifact_ref)
    sidecar$codec_metadata <- artifactStoreNormalizeCodecMetadata(
        sidecar$codec_metadata
    )
    ref <- sidecar$artifact_ref
    registration <- sidecar$registration
    registration$provenance_ids <- artifactStoreCharacterVector(
        registration$provenance_ids
    )
    expected_registration <- list(
        artifact_id = ref$artifact_id,
        logical_key = ref$logical_key,
        provenance_ids = ref$provenance_ids,
        status = ref$status
    )
    if (!identical(registration, expected_registration) ||
        !identical(ref$status, "committed") ||
        !identical(sidecar$codec_metadata$codec, ref$codec) ||
        !identical(sidecar$codec_metadata$payload_schema, ref$payload_schema) ||
        !identical(
            sidecar$codec_metadata$semantic_digest,
            ref$hash_policy$semantic$digest
        )) {
        artifactStoreAbort(
            "artifact sidecar registration or codec metadata disagrees with its ref",
            "multischolar_malformed_artifact_sidecar"
        )
    }
    artifactRefValidateId(sidecar$operation_id, "operation_id", "op")
    artifactStoreValidateGuarantees(sidecar$guarantees)
    if (!artifactRefValidUtc(sidecar$created_at)) {
        artifactStoreAbort(
            "artifact sidecar timestamp is invalid",
            "multischolar_malformed_artifact_sidecar"
        )
    }
    sidecar$registration <- registration
    sidecar
}

artifactStoreReadSidecar <- function(store, relative_path,
                                     validate_payload = TRUE) {
    store <- validateArtifactStore(store)
    sidecar <- artifactStoreValidateSidecar(
        artifactStoreReadJson(store, relative_path)
    )
    artifactStoreValidateLogicalKey(store, sidecar$artifact_ref$logical_key)
    expected <- artifactStoreManagedPaths(
        store,
        sidecar$artifact_ref$logical_key,
        sidecar$artifact_ref$artifact_id
    )
    if (!identical(sidecar$artifact_ref$relative_path, expected$payload)) {
        artifactStoreAbort(
            "artifact sidecar payload path is not canonical for its logical key",
            "multischolar_malformed_artifact_sidecar"
        )
    }
    if (isTRUE(validate_payload)) {
        payload_path <- artifactStoreResolveFile(
            store,
            sidecar$artifact_ref$relative_path,
            must_exist = TRUE
        )
        shape <- artifactStorePayloadShape(payload_path, sidecar$codec_metadata)
        validateArtifactRefPayload(
            sidecar$artifact_ref,
            store$project_root,
            shape
        )
    }
    sidecar
}

artifactStoreValidateIntent <- function(store, intent) {
    required <- c(
        "schema", "schema_version", "protocol_version", "operation_id",
        "artifact_ref", "codec_metadata", "managed_paths", "temporary_paths",
        "guarantees", "created_at"
    )
    if (!is.list(intent) || !identical(names(intent), required) ||
        !identical(intent$schema, .artifactIntentSchema)) {
        artifactStoreAbort(
            "artifact write intent has unknown or missing fields",
            "multischolar_malformed_artifact_intent"
        )
    }
    schema_version <- as.integer(intent$schema_version)
    protocol_version <- as.integer(intent$protocol_version)
    if (is.na(schema_version) || is.na(protocol_version) ||
        schema_version != .artifactIntentSchemaVersion ||
        protocol_version != .artifactStoreProtocolVersion) {
        artifactStoreAbort(
            "artifact write intent version is unsupported",
            "multischolar_unsupported_artifact_intent_version"
        )
    }
    intent$schema_version <- schema_version
    intent$protocol_version <- protocol_version
    intent$artifact_ref <- artifactStoreNormalizeRef(intent$artifact_ref)
    if (!identical(intent$artifact_ref$status, "staged")) {
        artifactStoreAbort(
            "artifact write intent must contain a staged ref",
            "multischolar_malformed_artifact_intent"
        )
    }
    artifactStoreValidateLogicalKey(store, intent$artifact_ref$logical_key)
    intent$codec_metadata <- artifactStoreNormalizeCodecMetadata(
        intent$codec_metadata
    )
    if (!identical(intent$codec_metadata$codec, intent$artifact_ref$codec) ||
        !identical(
            intent$codec_metadata$payload_schema,
            intent$artifact_ref$payload_schema
        ) ||
        !identical(
            intent$codec_metadata$semantic_digest,
            intent$artifact_ref$hash_policy$semantic$digest
        )) {
        artifactStoreAbort(
            "artifact write intent codec metadata disagrees with its ref",
            "multischolar_malformed_artifact_intent"
        )
    }
    expected <- artifactStoreManagedPaths(
        store,
        intent$artifact_ref$logical_key,
        intent$artifact_ref$artifact_id
    )
    if (!identical(intent$managed_paths, expected) ||
        !identical(intent$artifact_ref$relative_path, expected$payload)) {
        artifactStoreAbort(
            "artifact write intent contains non-canonical managed paths",
            "multischolar_malformed_artifact_intent"
        )
    }
    temporary_names <- c(
        "payload_directory", "payload", "sidecar", "intent"
    )
    if (!identical(names(intent$temporary_paths), temporary_names)) {
        artifactStoreAbort(
            "artifact write intent temporary paths are incomplete",
            "multischolar_malformed_artifact_intent"
        )
    }
    intent$temporary_paths <- lapply(
        intent$temporary_paths,
        artifactNormalizeRelativePath
    )
    same_payload_parent <- identical(
        dirname(intent$temporary_paths$payload_directory),
        dirname(expected$payload_directory)
    )
    path_relationships <- same_payload_parent &&
        identical(
            intent$temporary_paths$payload,
            file.path(intent$temporary_paths$payload_directory, "payload.parquet")
        ) &&
        identical(dirname(intent$temporary_paths$sidecar), dirname(expected$sidecar)) &&
        identical(dirname(intent$temporary_paths$intent), dirname(expected$intent))
    if (!isTRUE(path_relationships)) {
        artifactStoreAbort(
            "artifact write intent temporary paths leave their final filesystem",
            "multischolar_malformed_artifact_intent"
        )
    }
    artifactRefValidateId(intent$operation_id, "operation_id", "op")
    artifactStoreValidateGuarantees(intent$guarantees)
    if (!artifactRefValidUtc(intent$created_at)) {
        artifactStoreAbort(
            "artifact write intent timestamp is invalid",
            "multischolar_malformed_artifact_intent"
        )
    }
    intent
}

artifactStoreReadIntent <- function(store, relative_path) {
    artifactStoreValidateIntent(
        store,
        artifactStoreReadJson(store, relative_path)
    )
}

artifactStoreListRelative <- function(store, root_relative, pattern) {
    root <- artifactStoreResolveFile(store, root_relative)
    if (!dir.exists(root)) return(character())
    candidates <- list.files(
        root,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    if (length(candidates) == 0L) return(character())
    normalized <- vapply(candidates, function(path) {
        normalizePath(path, winslash = "/", mustWork = TRUE)
    }, character(1))
    project_root <- store$project_root
    if (!all(vapply(
        normalized,
        artifactPathIsContained,
        logical(1),
        project_root = project_root
    ))) {
        artifactStoreAbort(
            "managed artifact inventory escapes through a symlink",
            "multischolar_artifact_path_escape"
        )
    }
    substring(normalized, nchar(project_root) + 2L)
}

artifactStoreSidecarPaths <- function(store) {
    artifactStoreListRelative(
        store,
        store$relative_paths$generations,
        "[.]artifact[.]json$"
    )
}

artifactStoreIntentPaths <- function(store) {
    artifactStoreListRelative(
        store,
        store$relative_paths$staging,
        "[.]intent[.]json$"
    )
}

artifactStoreAssertLogicalKeyAvailable <- function(store, logical_key) {
    artifactStoreValidateLogicalKey(store, logical_key)
    sidecar_paths <- artifactStoreSidecarPaths(store)
    for (path in sidecar_paths) {
        sidecar <- artifactStoreReadSidecar(store, path, validate_payload = FALSE)
        if (identical(sidecar$artifact_ref$logical_key, logical_key)) {
            artifactStoreAbort(
                "immutable artifact logical key already has a generation payload",
                "multischolar_duplicate_artifact_generation"
            )
        }
    }
    intent_paths <- artifactStoreIntentPaths(store)
    for (path in intent_paths) {
        intent <- artifactStoreReadIntent(store, path)
        if (identical(intent$artifact_ref$logical_key, logical_key)) {
            artifactStoreAbort(
                "artifact logical key already has an incomplete write",
                "multischolar_duplicate_artifact_generation"
            )
        }
    }
    invisible(TRUE)
}

artifactStoreRefForPayload <- function(intent, payload_path) {
    metadata <- intent$codec_metadata
    shape <- artifactStorePayloadShape(payload_path, metadata)
    staged <- intent$artifact_ref
    newArtifactRef(
        logical_key = staged$logical_key,
        relative_path = staged$relative_path,
        codec_id = staged$codec$id,
        payload_schema_id = staged$payload_schema$id,
        shape = shape,
        semantic_digest = staged$hash_policy$semantic$digest,
        byte_digest = artifactByteDigest(payload_path),
        provenance_ids = staged$provenance_ids,
        status = "committed",
        artifact_id = staged$artifact_id,
        created_at = staged$created_at,
        updated_at = artifactRefUtcNow()
    )
}

artifactStoreValidateRefAtPath <- function(store, ref, metadata, relative_path) {
    path <- artifactStoreResolveFile(store, relative_path, must_exist = TRUE)
    shape <- artifactStorePayloadShape(path, metadata)
    temporary_ref <- ref
    temporary_ref$relative_path <- relative_path
    validateArtifactRefPayload(temporary_ref, store$project_root, shape)
    invisible(shape)
}

artifactStoreRecoverIntent <- function(store, intent_path) {
    intent <- artifactStoreReadIntent(store, intent_path)
    managed <- intent$managed_paths
    temporary <- intent$temporary_paths
    exists_at <- function(path) {
        resolved <- artifactStoreResolveFile(store, path)
        file.exists(resolved) || dir.exists(resolved)
    }
    final_payload <- exists_at(managed$payload)
    final_sidecar <- exists_at(managed$sidecar)
    temporary_payload <- exists_at(temporary$payload)
    temporary_sidecar <- exists_at(temporary$sidecar)

    if (final_payload && final_sidecar) {
        artifactStoreReadSidecar(store, managed$sidecar, validate_payload = TRUE)
        unlink(artifactStoreResolveFile(store, managed$intent), force = FALSE)
        return(invisible("committed"))
    }
    if (!final_payload && !temporary_payload) {
        return(invisible("not_repairable"))
    }

    payload_relative <- if (final_payload) managed$payload else temporary$payload
    if (temporary_sidecar) {
        sidecar <- artifactStoreReadSidecar(
            store,
            temporary$sidecar,
            validate_payload = FALSE
        )
        ref <- sidecar$artifact_ref
    } else if (final_sidecar) {
        sidecar <- artifactStoreReadSidecar(
            store,
            managed$sidecar,
            validate_payload = FALSE
        )
        ref <- sidecar$artifact_ref
    } else {
        payload_path <- artifactStoreResolveFile(
            store,
            payload_relative,
            must_exist = TRUE
        )
        ref <- artifactStoreRefForPayload(intent, payload_path)
        sidecar <- artifactStoreNewSidecar(
            ref,
            intent$codec_metadata,
            intent$operation_id
        )
        artifactStoreWriteJson(store, sidecar, temporary$sidecar)
        artifactStoreReadSidecar(
            store,
            temporary$sidecar,
            validate_payload = FALSE
        )
        temporary_sidecar <- TRUE
    }
    artifactStoreValidateRefAtPath(
        store,
        ref,
        intent$codec_metadata,
        payload_relative
    )
    if (!final_payload) {
        artifactStorePublishDirectory(
            store,
            temporary$payload_directory,
            managed$payload_directory
        )
    }
    if (!final_sidecar) {
        artifactStorePublishFile(store, temporary$sidecar, managed$sidecar)
    }
    artifactStoreReadSidecar(store, managed$sidecar, validate_payload = TRUE)
    unlink(artifactStoreResolveFile(store, managed$intent), force = FALSE)
    invisible("recovered")
}
