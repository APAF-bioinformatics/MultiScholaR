# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactRefSchema <- "multischolar.artifact_ref"
.artifactRefSchemaVersion <- 1L
.artifactSemanticDigestVersion <- 1L
.artifactHashAlgorithm <- "sha256"
.artifactRefStatuses <- c("staged", "validated", "committed", "trashed")
.artifactLogicalKeyFields <- c(
    "project_id",
    "omic_type",
    "workflow_slug",
    "stage_id",
    "state_role",
    "generation_id"
)

artifactRefAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_ref_error"),
        ...
    )
}

artifactRefUtcNow <- function() {
    format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

artifactRefValidUtc <- function(value) {
    workflowCapabilityScalarString(value) &&
        grepl(
            "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]+)?Z$",
            value
        ) &&
        !is.na(as.POSIXct(
            value,
            format = "%Y-%m-%dT%H:%M:%OSZ",
            tz = "UTC"
        ))
}

artifactOpaqueId <- function(prefix = "art") {
    artifactValidatePathComponent(prefix, "artifact_id_prefix")
    entropy <- list(
        timestamp = artifactRefUtcNow(),
        pid = Sys.getpid(),
        elapsed = unname(proc.time()[[3L]]),
        temporary_name = tempfile(pattern = "multischolar-artifact-")
    )
    paste0(
        tolower(prefix),
        "_",
        digest::digest(entropy, algo = .artifactHashAlgorithm)
    )
}

artifactRefScalarVersion <- function(value, field) {
    valid <- length(value) == 1L && !is.na(value) &&
        is.numeric(value) && is.finite(value) &&
        identical(as.numeric(value), floor(as.numeric(value))) && value >= 1L
    if (!isTRUE(valid)) {
        artifactRefAbort(
            sprintf("artifact reference field '%s' must be a positive integer", field),
            "multischolar_invalid_artifact_ref_version"
        )
    }
    as.integer(value)
}

artifactRefValidateId <- function(value, field, prefix = NULL) {
    valid <- workflowCapabilityScalarString(value) &&
        grepl("^[a-z][a-z0-9_]*_[a-f0-9]{64}$", value)
    if (!is.null(prefix)) {
        valid <- valid && startsWith(value, paste0(prefix, "_"))
    }
    if (!isTRUE(valid)) {
        artifactRefAbort(
            sprintf("artifact reference field '%s' is not an opaque ID", field),
            "multischolar_invalid_artifact_id"
        )
    }
    value
}

artifactRefValidateLogicalKey <- function(logical_key) {
    if (!is.list(logical_key) ||
        !identical(names(logical_key), .artifactLogicalKeyFields)) {
        artifactRefAbort(
            "artifact logical key must contain the exact scoped key fields",
            "multischolar_invalid_artifact_logical_key"
        )
    }
    valid <- vapply(logical_key, workflowCapabilityScalarString, logical(1))
    if (!all(valid)) {
        artifactRefAbort(
            "artifact logical key fields must be non-empty scalar strings",
            "multischolar_invalid_artifact_logical_key"
        )
    }
    invisible(logical_key)
}

artifactRefValidateShape <- function(shape, field = "shape") {
    required <- c("kind", "rows", "columns", "payloads", "bytes")
    if (!is.list(shape) || !identical(names(shape), required) ||
        !workflowCapabilityScalarString(shape$kind)) {
        artifactRefAbort(
            sprintf("artifact reference field '%s' is malformed", field),
            "multischolar_invalid_artifact_shape"
        )
    }
    counts <- shape[c("rows", "columns", "payloads", "bytes")]
    valid_count <- function(value) {
        length(value) == 1L && is.numeric(value) &&
            (is.na(value) || (is.finite(value) && value >= 0 && value == floor(value)))
    }
    if (!all(vapply(counts, valid_count, logical(1)))) {
        artifactRefAbort(
            sprintf("artifact reference field '%s' has invalid dimensions", field),
            "multischolar_invalid_artifact_shape"
        )
    }
    invisible(shape)
}

artifactCanonicalizeValue <- function(value) {
    if (is.null(value)) {
        return(NULL)
    }
    if (is.list(value)) {
        value_names <- names(value)
        if (!is.null(value_names)) {
            if (any(!nzchar(value_names)) || anyDuplicated(value_names) > 0L) {
                artifactRefAbort(
                    "semantic digest input has ambiguous list names",
                    "multischolar_invalid_semantic_digest_input"
                )
            }
            value <- value[order(value_names, method = "radix")]
        }
        return(lapply(value, artifactCanonicalizeValue))
    }
    if (is.atomic(value)) {
        value_attributes <- attributes(value)
        unclassed <- value
        attributes(unclassed) <- NULL
        if (is.null(value_attributes)) {
            return(unclassed)
        }
        return(list(
            artifact_atomic_value = unclassed,
            artifact_atomic_attributes = artifactCanonicalizeValue(value_attributes)
        ))
    }
    artifactRefAbort(
        sprintf(
            "semantic digest input has unsupported type '%s' with class '%s'",
            typeof(value),
            paste(class(value), collapse = "/")
        ),
        "multischolar_invalid_semantic_digest_input"
    )
}

artifactSemanticDigest <- function(value, version = .artifactSemanticDigestVersion) {
    version <- artifactRefScalarVersion(version, "semantic_digest_version")
    if (!identical(version, .artifactSemanticDigestVersion)) {
        artifactRefAbort(
            "unsupported semantic digest version",
            "multischolar_unsupported_semantic_digest_version"
        )
    }
    canonical <- artifactCanonicalizeValue(list(
        schema = "multischolar.semantic_lineage",
        version = version,
        value = value
    ))
    encoded <- jsonlite::toJSON(
        canonical,
        auto_unbox = TRUE,
        null = "null",
        na = "string",
        digits = NA,
        POSIXt = "ISO8601",
        UTC = TRUE
    )
    digest::digest(
        as.character(encoded),
        algo = .artifactHashAlgorithm,
        serialize = FALSE
    )
}

artifactByteDigest <- function(path) {
    if (!workflowCapabilityScalarString(path) || !file.exists(path) || dir.exists(path)) {
        artifactRefAbort(
            "artifact byte digest requires one existing file",
            "multischolar_missing_artifact_payload"
        )
    }
    digest::digest(
        file = path,
        algo = .artifactHashAlgorithm,
        serialize = FALSE
    )
}

artifactRefValidateDigest <- function(value, field, allow_missing = FALSE) {
    if (isTRUE(allow_missing) && length(value) == 1L && is.na(value)) {
        return(value)
    }
    if (!workflowCapabilityScalarString(value) || !grepl("^[a-f0-9]{64}$", value)) {
        artifactRefAbort(
            sprintf("artifact reference field '%s' is not a SHA-256 digest", field),
            "multischolar_invalid_artifact_digest"
        )
    }
    value
}

newArtifactRef <- function(
    logical_key,
    relative_path,
    codec_id,
    payload_schema_id,
    shape,
    semantic_digest,
    codec_version = 1L,
    payload_schema_version = 1L,
    byte_digest = NA_character_,
    provenance_ids = character(),
    status = "staged",
    artifact_id = artifactOpaqueId(),
    created_at = artifactRefUtcNow(),
    updated_at = created_at
) {
    artifactRefValidateLogicalKey(logical_key)
    relative_path <- artifactNormalizeRelativePath(relative_path)
    if (!workflowCapabilityScalarString(codec_id) ||
        !workflowCapabilityScalarString(payload_schema_id)) {
        artifactRefAbort(
            "artifact codec and payload schema IDs must be scalar strings",
            "multischolar_invalid_artifact_codec"
        )
    }
    artifactRefValidateShape(shape)
    artifactRefValidateDigest(semantic_digest, "semantic_digest")
    if (length(provenance_ids) > 0L) {
        vapply(
            provenance_ids,
            artifactRefValidateId,
            character(1),
            field = "provenance_ids"
        )
    }
    ref <- structure(
        list(
            schema = .artifactRefSchema,
            schema_version = .artifactRefSchemaVersion,
            artifact_id = artifactRefValidateId(artifact_id, "artifact_id", "art"),
            logical_key = logical_key,
            relative_path = relative_path,
            codec = list(
                id = codec_id,
                version = artifactRefScalarVersion(codec_version, "codec_version")
            ),
            payload_schema = list(
                id = payload_schema_id,
                version = artifactRefScalarVersion(
                    payload_schema_version,
                    "payload_schema_version"
                )
            ),
            shape = shape,
            hash_policy = list(
                byte = list(
                    algorithm = .artifactHashAlgorithm,
                    digest = byte_digest
                ),
                semantic = list(
                    algorithm = .artifactHashAlgorithm,
                    version = .artifactSemanticDigestVersion,
                    digest = semantic_digest
                )
            ),
            provenance_ids = unname(as.character(provenance_ids)),
            status = status,
            created_at = created_at,
            updated_at = updated_at
        ),
        class = c("MultiScholaRArtifactRef", "list")
    )
    validateArtifactRef(ref)
}

validateArtifactRef <- function(ref) {
    required <- c(
        "schema", "schema_version", "artifact_id", "logical_key",
        "relative_path", "codec", "payload_schema", "shape", "hash_policy",
        "provenance_ids", "status", "created_at", "updated_at"
    )
    if (!is.list(ref) || !identical(names(ref), required)) {
        artifactRefAbort(
            "artifact reference has an unknown or incomplete structure",
            "multischolar_malformed_artifact_ref"
        )
    }
    if (!identical(ref$schema, .artifactRefSchema) ||
        !identical(ref$schema_version, .artifactRefSchemaVersion)) {
        artifactRefAbort(
            "unsupported artifact reference schema version",
            "multischolar_unsupported_artifact_ref_version"
        )
    }
    artifactRefValidateId(ref$artifact_id, "artifact_id", "art")
    artifactRefValidateLogicalKey(ref$logical_key)
    artifactNormalizeRelativePath(ref$relative_path)
    for (field in c("codec", "payload_schema")) {
        value <- ref[[field]]
        if (!is.list(value) || !identical(names(value), c("id", "version")) ||
            !workflowCapabilityScalarString(value$id)) {
            artifactRefAbort(
                sprintf("artifact reference field '%s' is malformed", field),
                "multischolar_malformed_artifact_ref"
            )
        }
        artifactRefScalarVersion(value$version, paste0(field, "_version"))
        if (!identical(value$version, 1L)) {
            artifactRefAbort(
                sprintf("unsupported artifact %s version", field),
                "multischolar_unsupported_artifact_codec_version"
            )
        }
    }
    artifactRefValidateShape(ref$shape)
    hash_policy <- ref$hash_policy
    valid_hash_policy <- is.list(hash_policy) &&
        identical(names(hash_policy), c("byte", "semantic")) &&
        identical(hash_policy$byte$algorithm, .artifactHashAlgorithm) &&
        identical(hash_policy$semantic$algorithm, .artifactHashAlgorithm) &&
        identical(
            hash_policy$semantic$version,
            .artifactSemanticDigestVersion
        )
    if (!isTRUE(valid_hash_policy)) {
        artifactRefAbort(
            "artifact reference hash policy is unsupported",
            "multischolar_unsupported_artifact_hash_policy"
        )
    }
    allow_missing_byte <- identical(ref$status, "staged")
    artifactRefValidateDigest(
        hash_policy$byte$digest,
        "byte_digest",
        allow_missing = allow_missing_byte
    )
    artifactRefValidateDigest(hash_policy$semantic$digest, "semantic_digest")
    if (!is.character(ref$provenance_ids) || anyDuplicated(ref$provenance_ids) > 0L) {
        artifactRefAbort(
            "artifact provenance IDs must be a unique character vector",
            "multischolar_invalid_artifact_provenance"
        )
    }
    if (length(ref$provenance_ids) > 0L) {
        vapply(
            ref$provenance_ids,
            artifactRefValidateId,
            character(1),
            field = "provenance_ids"
        )
    }
    if (!workflowCapabilityScalarString(ref$status) ||
        !ref$status %in% .artifactRefStatuses) {
        artifactRefAbort(
            "artifact reference status is unsupported",
            "multischolar_invalid_artifact_status"
        )
    }
    if (!artifactRefValidUtc(ref$created_at) || !artifactRefValidUtc(ref$updated_at)) {
        artifactRefAbort(
            "artifact reference timestamps must be UTC ISO-8601 values",
            "multischolar_invalid_artifact_timestamp"
        )
    }
    ref
}

validateArtifactRefPayload <- function(ref, project_root, actual_shape) {
    ref <- validateArtifactRef(ref)
    artifactRefValidateShape(actual_shape, "actual_shape")
    compared <- c("kind", "rows", "columns", "payloads", "bytes")
    mismatched <- compared[vapply(compared, function(field) {
        expected <- ref$shape[[field]]
        actual <- actual_shape[[field]]
        if (is.na(expected)) {
            return(FALSE)
        }
        if (identical(field, "kind")) {
            return(!identical(expected, actual))
        }
        !isTRUE(as.numeric(expected) == as.numeric(actual))
    }, logical(1))]
    if (length(mismatched) > 0L) {
        artifactRefAbort(
            sprintf(
                "artifact payload shape mismatch for: %s",
                paste(mismatched, collapse = ", ")
            ),
            "multischolar_artifact_shape_mismatch"
        )
    }
    path <- artifactResolveContainedPath(
        project_root,
        ref$relative_path,
        must_exist = TRUE
    )
    actual_digest <- artifactByteDigest(path)
    if (!identical(actual_digest, ref$hash_policy$byte$digest)) {
        artifactRefAbort(
            "artifact payload byte digest does not match its reference",
            "multischolar_artifact_hash_mismatch"
        )
    }
    invisible(ref)
}
