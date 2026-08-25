# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactLipidBindingSchema <- "multischolar.lipidomics_s4_binding"
.artifactLipidBindingVersion <- 1L
.artifactLipidRoleManifestSchema <- "multischolar.lipidomics_s4_roles"
.artifactLipidRoleManifestVersion <- 1L
.artifactLipidNestedS4Codec <- "multischolar.s4.nested_exact_role"
.artifactLipidNestedS4CodecVersion <- 1L
.artifactLipidSerializedRoleCodec <- "multischolar.r.serialized_role"
.artifactLipidSerializedRoleCodecVersion <- 1L
.artifactLipidSerializationVersion <- 3L
.artifactLipidSerializedChunkBytes <- 512L * 1024L

.artifactLipidSerializedSlots <- c(
    "fit.eb" = "lipidomics.da.fit_model",
    "num_sig_diff_exp_bar_plot" = "lipidomics.da.significant_count_plots",
    "volcano_plot" = "lipidomics.da.volcano_plots",
    "interactive_volcano_plot" = "lipidomics.da.interactive_volcano_plots",
    "p_value_dist_plot" = "lipidomics.da.p_value_distribution_plots"
)

artifactLipidCodecDescriptor <- function(class_name) {
    declarations <- artifactLipidomicsCodecDeclarations()
    matches <- Filter(\(codec) identical(codec$class_name, class_name), declarations)
    if (length(matches) != 1L) {
        artifactCodecAbort(
            sprintf("S4 class '%s' has no certified lipidomics codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    list(id = matches[[1L]]$codec_id, version = matches[[1L]]$codec_version)
}

artifactLipidLogicalKey <- function(logical_key) {
    artifactRefValidateLogicalKey(logical_key)
    if (!identical(logical_key$omic_type, "lipidomics")) {
        artifactCodecAbort(
            "lipidomics S4 codecs require a lipidomics logical key",
            "multischolar_invalid_lipidomics_codec_identity"
        )
    }
    logical_key
}

artifactLipidAssayIdentity <- function(logical_key, owner, assay_name, position) {
    paste0("lipassay_", artifactSemanticDigest(list(
        logical_key = logical_key,
        owner = owner,
        assay_name = assay_name,
        position = as.integer(position)
    )))
}

artifactLipidIdentityBinding <- function(logical_key, owner, assay_names) {
    logical_key <- artifactLipidLogicalKey(logical_key)
    if (!workflowCapabilityScalarString(owner) || anyNA(assay_names) ||
        any(!nzchar(assay_names)) || anyDuplicated(assay_names) > 0L) {
        artifactCodecAbort(
            "lipidomics assay identity requires unique named assays",
            "multischolar_invalid_lipidomics_codec_identity"
        )
    }
    assays <- lapply(seq_along(assay_names), \(index) {
        list(
            name = assay_names[[index]],
            position = as.integer(index),
            assay_id = artifactLipidAssayIdentity(
                logical_key,
                owner,
                assay_names[[index]],
                index
            )
        )
    })
    binding <- list(
        schema = .artifactLipidBindingSchema,
        schema_version = .artifactLipidBindingVersion,
        logical_key = logical_key,
        owner = owner,
        assays = assays,
        binding_digest = NULL
    )
    binding$binding_digest <- artifactSemanticDigest(
        binding[names(binding) != "binding_digest"]
    )
    binding
}

artifactLipidValidateIdentityBinding <- function(binding, owner, assay_names) {
    required <- c(
        "schema", "schema_version", "logical_key", "owner", "assays",
        "binding_digest"
    )
    valid <- is.list(binding) && identical(names(binding), required) &&
        identical(binding$schema, .artifactLipidBindingSchema) &&
        identical(binding$schema_version, .artifactLipidBindingVersion) &&
        identical(binding$owner, owner)
    if (!isTRUE(valid)) {
        artifactCodecAbort(
            "lipidomics S4 identity binding is malformed or unsupported",
            "multischolar_invalid_lipidomics_codec_identity"
        )
    }
    expected <- artifactLipidIdentityBinding(binding$logical_key, owner, assay_names)
    if (!identical(binding, expected)) {
        artifactCodecAbort(
            "lipidomics S4 assay identity binding does not match its logical key",
            "multischolar_lipidomics_codec_identity_mismatch"
        )
    }
    binding
}

artifactLipidAssayNames <- function(node, owner = "lipid_data") {
    if (!is.list(node) || !identical(node$node_type, "list") ||
        is.null(node$names)) {
        artifactCodecAbort(
            sprintf("%s must be an ordered named assay list", owner),
            "multischolar_lipidomics_codec_shape_mismatch"
        )
    }
    assay_names <- artifactDecodeInlineCharacter(node$names)
    if (length(assay_names) != length(node$values) || anyNA(assay_names) ||
        any(!nzchar(assay_names)) || anyDuplicated(assay_names) > 0L) {
        artifactCodecAbort(
            sprintf("%s contains missing or duplicate assay names", owner),
            "multischolar_lipidomics_codec_shape_mismatch"
        )
    }
    assay_names
}

artifactLipidNodePayloadKeys <- function(node) {
    if (!is.list(node) || !workflowCapabilityScalarString(node$node_type)) {
        return(character())
    }
    direct <- if (workflowCapabilityScalarString(node$payload_key)) {
        node$payload_key
    } else {
        character()
    }
    if (identical(node$node_type, "list")) {
        nested <- unlist(lapply(
            c(node$values, list(node$attributes)),
            artifactLipidNodePayloadKeys
        ), use.names = FALSE)
        return(c(direct, nested))
    }
    if (identical(node$node_type, "s4_exact_role")) {
        return(unlist(node$payload_keys, use.names = FALSE))
    }
    direct
}

artifactLipidNodeStorage <- function(node) {
    switch(node$node_type,
        s4_exact_role = "nested_exact_s4",
        immutable_r_role = "immutable_r_serialized",
        rectangular = "rectangular",
        atomic_rectangular = "rectangular",
        nested_rectangular = "rectangular_nested_metadata",
        list = "structured_list",
        "inline"
    )
}

artifactLipidRoleManifest <- function(slot_values, binding) {
    slot_roles <- lapply(names(slot_values), \(slot_name) {
        node <- slot_values[[slot_name]]
        list(
            slot_name = slot_name,
            storage = artifactLipidNodeStorage(node),
            role_id = if (workflowCapabilityScalarString(node$role_id)) {
                node$role_id
            } else {
                paste0("lipidomics.s4.slot.", slot_name)
            },
            payload_keys = unname(artifactLipidNodePayloadKeys(node))
        )
    })
    manifest <- list(
        schema = .artifactLipidRoleManifestSchema,
        schema_version = .artifactLipidRoleManifestVersion,
        slots = slot_roles,
        assays = binding$assays
    )
    artifactWorkflowStateAssertSafeMetadata(manifest, "lipidomics role manifest")
    manifest
}

artifactLipidAttachContract <- function(bundle, logical_key, owner, assay_names) {
    binding <- artifactLipidIdentityBinding(logical_key, owner, assay_names)
    bundle$metadata$identity_binding <- binding
    bundle$metadata$role_manifest <- artifactLipidRoleManifest(
        bundle$metadata$slot_values,
        binding
    )
    bundle$metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(bundle$metadata)
    )
    bundle
}

artifactLipidRoleClassSignature <- function(role_value) {
    member_classes <- if (is.list(role_value)) {
        lapply(role_value, \(item) unname(as.character(class(item))))
    } else {
        list()
    }
    list(
        class = unname(as.character(class(role_value))),
        names = if (is.null(names(role_value))) NULL else unname(names(role_value)),
        member_classes = member_classes
    )
}

artifactLipidRawChunks <- function(value) {
    if (length(value) == 0L) return(list(raw()))
    starts <- seq.int(1L, length(value), by = .artifactLipidSerializedChunkBytes)
    lapply(starts, \(start) {
        end <- min(length(value), start + .artifactLipidSerializedChunkBytes - 1L)
        value[start:end]
    })
}

artifactLipidEncodeSerializedRole <- function(
    role_value,
    state,
    slot_name,
    role_id
) {
    serialized <- serialize(
        role_value,
        connection = NULL,
        ascii = FALSE,
        version = .artifactLipidSerializationVersion
    )
    chunks <- artifactLipidRawChunks(serialized)
    table <- data.frame(
        chunk_index = seq_along(chunks),
        serialized_chunk = vapply(chunks, jsonlite::base64_enc, character(1)),
        stringsAsFactors = FALSE
    )
    encoded <- encodeArtifactTable(
        table,
        stable_key = "chunk_index",
        owner = paste("lipidomics immutable role", role_id)
    )
    list(
        node_type = "immutable_r_role",
        codec = list(
            id = .artifactLipidSerializedRoleCodec,
            version = .artifactLipidSerializedRoleCodecVersion
        ),
        role_id = role_id,
        slot_name = slot_name,
        payload_key = artifactAddRectangularPayload(state, encoded),
        serialization_version = .artifactLipidSerializationVersion,
        byte_count = as.numeric(length(serialized)),
        byte_digest = digest::digest(
            serialized,
            algo = .artifactHashAlgorithm,
            serialize = FALSE
        ),
        class_signature = artifactLipidRoleClassSignature(role_value)
    )
}

artifactLipidRemapPayloadKeys <- function(value, key_map) {
    if (!is.list(value)) return(value)
    remapped <- lapply(value, artifactLipidRemapPayloadKeys, key_map = key_map)
    names(remapped) <- names(value)
    payload_key <- remapped[["payload_key", exact = TRUE]]
    if (workflowCapabilityScalarString(payload_key)) {
        remapped$payload_key <- unname(key_map[[payload_key]])
    }
    payload_keys <- remapped[["payload_keys", exact = TRUE]]
    if (!is.null(payload_keys)) {
        remapped$payload_keys <- unname(key_map[unlist(payload_keys)])
    }
    remapped
}

artifactLipidMergeNestedBundle <- function(state, bundle) {
    old_keys <- names(bundle$payloads)
    new_keys <- vapply(old_keys, \(old_key) {
        artifactAddRectangularPayload(state, structure(
            list(
                payload = bundle$payloads[[old_key]],
                metadata = bundle$metadata$payloads[[old_key]]
            ),
            class = c("MultiScholaRArtifactRectangular", "list")
        ))
    }, character(1))
    key_map <- stats::setNames(new_keys, old_keys)
    metadata <- artifactLipidRemapPayloadKeys(bundle$metadata, key_map)
    metadata$payloads <- state$payload_metadata[new_keys]
    metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(metadata)
    )
    list(metadata = metadata, payload_keys = unname(new_keys))
}

artifactLipidEncodeNestedS4 <- function(value, state, logical_key, inline_limit_bytes) {
    nested <- dehydrateLipidomicsS4Artifact(
        value,
        logical_key,
        inline_limit_bytes = inline_limit_bytes
    )
    merged <- artifactLipidMergeNestedBundle(state, nested)
    list(
        node_type = "s4_exact_role",
        codec = list(
            id = .artifactLipidNestedS4Codec,
            version = .artifactLipidNestedS4CodecVersion
        ),
        role_id = "lipidomics.da.embedded_assay",
        slot_name = "theObject",
        metadata = merged$metadata,
        payload_keys = merged$payload_keys
    )
}

artifactLipidDaBundle <- function(value, codec, logical_key, inline_limit_bytes) {
    state <- new.env(parent = emptyenv())
    state$payloads <- list()
    state$payload_metadata <- list()
    slot_names <- methods::slotNames(value)
    slot_values <- lapply(slot_names, \(slot_name) {
        slot_value <- methods::slot(value, slot_name)
        if (identical(slot_name, "theObject")) {
            return(artifactLipidEncodeNestedS4(
                slot_value,
                state,
                logical_key,
                inline_limit_bytes
            ))
        }
        role_id <- if (slot_name %in% names(.artifactLipidSerializedSlots)) {
            .artifactLipidSerializedSlots[[slot_name]]
        } else {
            NULL
        }
        if (!is.null(role_id)) {
            return(artifactLipidEncodeSerializedRole(
                slot_value,
                state,
                slot_name,
                role_id
            ))
        }
        artifactEncodeValueNode(
            slot_value,
            state,
            paste0(class(value), "@", slot_name),
            inline_limit_bytes
        )
    })
    names(slot_values) <- slot_names
    metadata <- list(
        schema = .artifactS4BundleSchema,
        schema_version = .artifactS4BundleSchemaVersion,
        codec = list(id = codec$codec_id, version = codec$codec_version),
        class_name = class(value)[[1L]],
        slot_names = slot_names,
        slot_values = slot_values,
        payloads = state$payload_metadata,
        semantic_digest = NULL,
        created_at = artifactRefUtcNow()
    )
    bundle <- structure(
        list(metadata = metadata, payloads = state$payloads),
        class = c("MultiScholaRArtifactS4Bundle", "list")
    )
    assay_node <- slot_values$theObject$metadata$slot_values$lipid_data
    artifactLipidAttachContract(
        bundle,
        logical_key,
        class(value)[[1L]],
        artifactLipidAssayNames(assay_node, "embedded lipid_data")
    )
}

dehydrateLipidomicsS4Artifact <- function(
    value,
    logical_key,
    inline_limit_bytes = .artifactInlineLimitBytes
) {
    if (!isS4(value) || length(class(value)) != 1L) {
        artifactCodecAbort(
            "lipidomics codec requires one supported S4 object",
            "multischolar_missing_exact_s4_codec"
        )
    }
    logical_key <- artifactLipidLogicalKey(logical_key)
    class_name <- class(value)[[1L]]
    descriptor <- artifactLipidCodecDescriptor(class_name)
    codec <- artifactLipidomicsCodecDeclarations()[[descriptor$id]]
    validity <- methods::validObject(value, test = TRUE)
    if (!identical(validity, TRUE)) {
        artifactCodecAbort(
            sprintf("lipidomics S4 object is invalid: %s", paste(validity, collapse = "; ")),
            "multischolar_invalid_lipidomics_s4_object",
            owner = class_name
        )
    }
    if (identical(class_name, "LipidomicsAssayData")) {
        bundle <- dehydrateExactS4Artifact(value, codec, inline_limit_bytes)
        assay_names <- names(methods::slot(value, "lipid_data"))
        return(artifactLipidAttachContract(
            bundle,
            logical_key,
            class_name,
            assay_names
        ))
    }
    artifactLipidDaBundle(value, codec, logical_key, inline_limit_bytes)
}

artifactLipidPayloadExists <- function(key, payloads, payload_metadata) {
    workflowCapabilityScalarString(key) && key %in% names(payloads) &&
        key %in% names(payload_metadata)
}

artifactLipidPreflightNode <- function(node, payloads, payload_metadata) {
    if (!is.list(node) || !workflowCapabilityScalarString(node$node_type)) {
        artifactCodecAbort(
            "lipidomics S4 bundle contains a malformed value node",
            "multischolar_lipidomics_codec_shape_mismatch"
        )
    }
    payload_types <- c(
        "rectangular", "atomic_rectangular", "nested_rectangular",
        "immutable_r_role"
    )
    if (node$node_type %in% payload_types &&
        !artifactLipidPayloadExists(node$payload_key, payloads, payload_metadata)) {
        artifactCodecAbort(
            "lipidomics S4 value node references a missing payload",
            "multischolar_invalid_artifact_payload"
        )
    }
    if (identical(node$node_type, "nested_rectangular") && !identical(
        node$codec,
        list(id = .artifactNestedNodeCodec, version = .artifactNestedNodeCodecVersion)
    )) {
        artifactCodecAbort(
            "nested lipidomics metadata codec version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    if (identical(node$node_type, "immutable_r_role") && !identical(
        node$codec,
        list(
            id = .artifactLipidSerializedRoleCodec,
            version = .artifactLipidSerializedRoleCodecVersion
        )
    )) {
        artifactCodecAbort(
            "lipidomics immutable role codec version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    if (identical(node$node_type, "immutable_r_role")) {
        required <- c(
            "node_type", "codec", "role_id", "slot_name", "payload_key",
            "serialization_version", "byte_count", "byte_digest",
            "class_signature"
        )
        valid <- identical(names(node), required) &&
            workflowCapabilityScalarString(node$role_id) &&
            workflowCapabilityScalarString(node$slot_name) &&
            identical(
                node$serialization_version,
                .artifactLipidSerializationVersion
            ) && length(node$byte_count) == 1L && is.numeric(node$byte_count) &&
            !is.na(node$byte_count) && is.finite(node$byte_count) &&
            node$byte_count >= 0 && node$byte_count == floor(node$byte_count)
        if (!isTRUE(valid)) {
            artifactCodecAbort(
                "lipidomics immutable role descriptor is malformed",
                "multischolar_lipidomics_codec_role_mismatch"
            )
        }
        artifactRefValidateDigest(node$byte_digest, "serialized role byte digest")
        artifactWorkflowStateAssertSafeMetadata(
            node$class_signature,
            "serialized role class signature"
        )
    }
    if (identical(node$node_type, "s4_exact_role")) {
        expected <- list(
            id = .artifactLipidNestedS4Codec,
            version = .artifactLipidNestedS4CodecVersion
        )
        keys <- unlist(node$payload_keys, use.names = FALSE)
        required <- c(
            "node_type", "codec", "role_id", "slot_name", "metadata",
            "payload_keys"
        )
        nested_metadata <- node$metadata
        nested_version <- is.list(nested_metadata) &&
            identical(nested_metadata$schema, .artifactS4BundleSchema) &&
            identical(
                nested_metadata$schema_version,
                .artifactS4BundleSchemaVersion
            ) && identical(
                nested_metadata$class_name,
                "LipidomicsAssayData"
            ) && identical(
                nested_metadata$codec,
                artifactLipidCodecDescriptor("LipidomicsAssayData")
            ) && identical(names(nested_metadata$payloads), keys)
        if (!identical(names(node), required) ||
            !workflowCapabilityScalarString(node$role_id) ||
            !workflowCapabilityScalarString(node$slot_name) ||
            !identical(node$codec, expected) || !isTRUE(nested_version) ||
            !all(vapply(
            keys,
            artifactLipidPayloadExists,
            logical(1),
            payloads = payloads,
            payload_metadata = payload_metadata
        ))) {
            artifactCodecAbort(
                "nested lipidomics S4 role is malformed or unsupported",
                "multischolar_unsupported_artifact_codec_version"
            )
        }
        lapply(
            node$metadata$slot_values,
            artifactLipidPreflightNode,
            payloads = payloads,
            payload_metadata = payload_metadata
        )
    }
    if (identical(node$node_type, "list")) {
        children <- node$values
        if (!is.null(node$attributes)) {
            children <- c(children, list(node$attributes))
        }
        lapply(
            children,
            artifactLipidPreflightNode,
            payloads = payloads,
            payload_metadata = payload_metadata
        )
    }
    allowed <- c(
        "null", "raw_inline", "atomic_inline", "rectangular",
        "atomic_rectangular", "nested_rectangular", "list",
        "s4_exact_role", "immutable_r_role"
    )
    if (!node$node_type %in% allowed) {
        artifactCodecAbort(
            "lipidomics S4 value node type is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    invisible(node)
}

artifactLipidSerializedRaw <- function(node, payloads, payload_metadata) {
    table <- decodeArtifactRectangular(
        payloads[[node$payload_key]],
        payload_metadata[[node$payload_key]]
    )
    expected_names <- c("chunk_index", "serialized_chunk")
    if (!identical(names(table), expected_names) || nrow(table) < 1L ||
        !identical(table$chunk_index, seq_len(nrow(table))) ||
        anyNA(table$serialized_chunk)) {
        artifactCodecAbort(
            "lipidomics immutable role payload is malformed",
            "multischolar_invalid_artifact_payload",
            owner = node$role_id
        )
    }
    chunks <- lapply(table$serialized_chunk, jsonlite::base64_dec)
    serialized <- do.call(c, chunks)
    digest <- digest::digest(
        serialized,
        algo = .artifactHashAlgorithm,
        serialize = FALSE
    )
    if (!identical(as.numeric(length(serialized)), as.numeric(node$byte_count)) ||
        !identical(digest, node$byte_digest)) {
        artifactCodecAbort(
            sprintf(
                paste(
                    "lipidomics immutable role '%s' payload digest does not match",
                    "(%s/%s bytes; %s/%s digest)"
                ),
                node$role_id,
                node$byte_count,
                length(serialized),
                node$byte_digest,
                digest
            ),
            "multischolar_artifact_hash_mismatch",
            owner = node$role_id,
            expected_byte_count = node$byte_count,
            observed_byte_count = length(serialized),
            expected_digest = node$byte_digest,
            observed_digest = digest
        )
    }
    serialized
}

artifactLipidDecodeSerializedRole <- function(node, payloads, payload_metadata) {
    serialized <- artifactLipidSerializedRaw(node, payloads, payload_metadata)
    value <- tryCatch(
        unserialize(serialized),
        error = \(error) artifactCodecAbort(
            "lipidomics immutable R role could not be hydrated",
            "multischolar_invalid_artifact_payload",
            parent = error,
            owner = node$role_id
        )
    )
    if (!identical(artifactLipidRoleClassSignature(value), node$class_signature)) {
        artifactCodecAbort(
            "lipidomics immutable R role class signature changed",
            "multischolar_inexact_s4_hydration",
            owner = node$role_id
        )
    }
    value
}

artifactLipidNestedBundle <- function(node, payloads) {
    keys <- unlist(node$payload_keys, use.names = FALSE)
    structure(
        list(metadata = node$metadata, payloads = payloads[keys]),
        class = c("MultiScholaRArtifactS4Bundle", "list")
    )
}

artifactLipidValidateNodePayload <- function(node, payloads, payload_metadata) {
    if (identical(node$node_type, "immutable_r_role")) {
        artifactLipidDecodeSerializedRole(node, payloads, payload_metadata)
    } else if (identical(node$node_type, "s4_exact_role")) {
        validateLipidomicsS4Bundle(artifactLipidNestedBundle(node, payloads))
    } else if (identical(node$node_type, "list")) {
        children <- node$values
        if (!is.null(node$attributes)) {
            children <- c(children, list(node$attributes))
        }
        lapply(
            children,
            artifactLipidValidateNodePayload,
            payloads = payloads,
            payload_metadata = payload_metadata
        )
    } else {
        artifactDecodeValueNode(node, payloads, payload_metadata)
    }
    invisible(node)
}

artifactLipidValidateSpecialRoles <- function(metadata) {
    if (!identical(metadata$class_name, "LipidomicsDifferentialAbundanceResults")) {
        return(invisible(metadata))
    }
    embedded <- metadata$slot_values$theObject
    valid_embedded <- identical(embedded$node_type, "s4_exact_role") &&
        identical(embedded$role_id, "lipidomics.da.embedded_assay") &&
        identical(embedded$slot_name, "theObject")
    role_slots <- names(.artifactLipidSerializedSlots)
    valid_roles <- vapply(role_slots, \(slot_name) {
        node <- metadata$slot_values[[slot_name]]
        identical(node$node_type, "immutable_r_role") &&
            identical(node$role_id, .artifactLipidSerializedSlots[[slot_name]]) &&
            identical(node$slot_name, slot_name) &&
            identical(node$serialization_version, .artifactLipidSerializationVersion)
    }, logical(1))
    if (!isTRUE(valid_embedded) || !all(valid_roles)) {
        artifactCodecAbort(
            "lipidomics DA model, plot, or embedded assay role is inconsistent",
            "multischolar_lipidomics_codec_role_mismatch"
        )
    }
    invisible(metadata)
}

artifactLipidBundleAssayNames <- function(metadata) {
    node <- if (identical(metadata$class_name, "LipidomicsAssayData")) {
        metadata$slot_values$lipid_data
    } else {
        metadata$slot_values$theObject$metadata$slot_values$lipid_data
    }
    artifactLipidAssayNames(node)
}

validateLipidomicsS4Bundle <- function(bundle, validate_payloads = TRUE) {
    if (!is.logical(validate_payloads) || length(validate_payloads) != 1L ||
        is.na(validate_payloads)) {
        artifactCodecAbort(
            "lipidomics payload validation flag must be TRUE or FALSE",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    if (!is.list(bundle) || !is.list(bundle$metadata) ||
        !workflowCapabilityScalarString(bundle$metadata$class_name)) {
        artifactCodecAbort(
            "lipidomics S4 artifact bundle is malformed",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    bundle$metadata <- artifactLipidNormalizeMetadataTransport(bundle$metadata)
    descriptor <- artifactLipidCodecDescriptor(bundle$metadata$class_name)
    codec <- artifactLipidomicsCodecDeclarations()[[descriptor$id]]
    bundle <- validateExactS4Bundle(bundle, codec)
    metadata <- bundle$metadata
    assay_names <- artifactLipidBundleAssayNames(metadata)
    binding <- artifactLipidValidateIdentityBinding(
        metadata$identity_binding,
        metadata$class_name,
        assay_names
    )
    expected_manifest <- artifactLipidRoleManifest(metadata$slot_values, binding)
    if (!identical(metadata$role_manifest, expected_manifest)) {
        artifactCodecAbort(
            "lipidomics S4 role manifest does not match its value nodes",
            "multischolar_lipidomics_codec_role_mismatch"
        )
    }
    artifactLipidValidateSpecialRoles(metadata)
    lapply(
        metadata$slot_values,
        artifactLipidPreflightNode,
        payloads = bundle$payloads,
        payload_metadata = metadata$payloads
    )
    if (isTRUE(validate_payloads)) {
        lapply(
            metadata$slot_values,
            artifactLipidValidateNodePayload,
            payloads = bundle$payloads,
            payload_metadata = metadata$payloads
        )
    }
    bundle
}

artifactLipidDecodeNode <- function(node, payloads, payload_metadata) {
    if (identical(node$node_type, "immutable_r_role")) {
        return(artifactLipidDecodeSerializedRole(node, payloads, payload_metadata))
    }
    if (identical(node$node_type, "s4_exact_role")) {
        return(hydrateLipidomicsS4Artifact(
            artifactLipidNestedBundle(node, payloads)
        ))
    }
    artifactDecodeValueNode(node, payloads, payload_metadata)
}

hydrateLipidomicsS4Artifact <- function(bundle) {
    bundle <- validateLipidomicsS4Bundle(bundle, validate_payloads = FALSE)
    metadata <- bundle$metadata
    slot_values <- lapply(
        metadata$slot_values,
        artifactLipidDecodeNode,
        payloads = bundle$payloads,
        payload_metadata = metadata$payloads
    )
    object <- do.call(
        methods::new,
        c(list(Class = metadata$class_name), slot_values)
    )
    validity <- methods::validObject(object, test = TRUE)
    if (!identical(validity, TRUE)) {
        artifactCodecAbort(
            sprintf(
                "hydrated lipidomics object is invalid: %s",
                paste(validity, collapse = "; ")
            ),
            "multischolar_invalid_hydrated_s4_object",
            owner = metadata$class_name
        )
    }
    object
}

encodeLipidomicsS4Artifact <- dehydrateLipidomicsS4Artifact
decodeLipidomicsS4Artifact <- hydrateLipidomicsS4Artifact
