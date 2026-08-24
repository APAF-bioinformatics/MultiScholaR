# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactMetabBindingSchema <- "multischolar.metabolomics_s4_binding"
.artifactMetabBindingVersion <- 1L
.artifactMetabRoleManifestSchema <- "multischolar.metabolomics_s4_roles"
.artifactMetabRoleManifestVersion <- 1L
.artifactMetabNestedS4Codec <- "multischolar.s4.nested_exact_role"
.artifactMetabNestedS4CodecVersion <- 1L
.artifactMetabSerializedRoleCodec <- "multischolar.r.serialized_role"
.artifactMetabSerializedRoleCodecVersion <- 1L
.artifactMetabSerializationVersion <- 3L
.artifactMetabSerializedChunkBytes <- 512L * 1024L

.artifactMetabSerializedSlots <- c(
    "fit.eb" = "metabolomics.da.fit_model",
    "num_sig_diff_exp_bar_plot" = "metabolomics.da.significant_count_plots",
    "volcano_plot" = "metabolomics.da.volcano_plots",
    "interactive_volcano_plot" = "metabolomics.da.interactive_volcano_plots",
    "p_value_dist_plot" = "metabolomics.da.p_value_distribution_plots"
)

artifactMetabCodecDescriptor <- function(class_name) {
    declarations <- artifactMetabolomicsCodecDeclarations()
    matches <- Filter(\(codec) identical(codec$class_name, class_name), declarations)
    if (length(matches) != 1L) {
        artifactCodecAbort(
            sprintf("S4 class '%s' has no certified metabolomics codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    list(id = matches[[1L]]$codec_id, version = matches[[1L]]$codec_version)
}

artifactMetabLogicalKey <- function(logical_key) {
    artifactRefValidateLogicalKey(logical_key)
    if (!identical(logical_key$omic_type, "metabolomics")) {
        artifactCodecAbort(
            "metabolomics S4 codecs require a metabolomics logical key",
            "multischolar_invalid_metabolomics_codec_identity"
        )
    }
    logical_key
}

artifactMetabAssayIdentity <- function(logical_key, owner, assay_name, position) {
    paste0("metassay_", artifactSemanticDigest(list(
        logical_key = logical_key,
        owner = owner,
        assay_name = assay_name,
        position = as.integer(position)
    )))
}

artifactMetabIdentityBinding <- function(logical_key, owner, assay_names) {
    logical_key <- artifactMetabLogicalKey(logical_key)
    if (!workflowCapabilityScalarString(owner) || anyNA(assay_names) ||
        any(!nzchar(assay_names)) || anyDuplicated(assay_names) > 0L) {
        artifactCodecAbort(
            "metabolomics assay identity requires unique named assays",
            "multischolar_invalid_metabolomics_codec_identity"
        )
    }
    assays <- lapply(seq_along(assay_names), \(index) {
        list(
            name = assay_names[[index]],
            position = as.integer(index),
            assay_id = artifactMetabAssayIdentity(
                logical_key,
                owner,
                assay_names[[index]],
                index
            )
        )
    })
    binding <- list(
        schema = .artifactMetabBindingSchema,
        schema_version = .artifactMetabBindingVersion,
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

artifactMetabValidateIdentityBinding <- function(binding, owner, assay_names) {
    required <- c(
        "schema", "schema_version", "logical_key", "owner", "assays",
        "binding_digest"
    )
    valid <- is.list(binding) && identical(names(binding), required) &&
        identical(binding$schema, .artifactMetabBindingSchema) &&
        identical(binding$schema_version, .artifactMetabBindingVersion) &&
        identical(binding$owner, owner)
    if (!isTRUE(valid)) {
        artifactCodecAbort(
            "metabolomics S4 identity binding is malformed or unsupported",
            "multischolar_invalid_metabolomics_codec_identity"
        )
    }
    expected <- artifactMetabIdentityBinding(binding$logical_key, owner, assay_names)
    if (!identical(binding, expected)) {
        artifactCodecAbort(
            "metabolomics S4 assay identity binding does not match its logical key",
            "multischolar_metabolomics_codec_identity_mismatch"
        )
    }
    binding
}

artifactMetabAssayNames <- function(node, owner = "metabolite_data") {
    if (!is.list(node) || !identical(node$node_type, "list") ||
        is.null(node$names)) {
        artifactCodecAbort(
            sprintf("%s must be an ordered named assay list", owner),
            "multischolar_metabolomics_codec_shape_mismatch"
        )
    }
    assay_names <- artifactDecodeInlineCharacter(node$names)
    if (length(assay_names) != length(node$values) || anyNA(assay_names) ||
        any(!nzchar(assay_names)) || anyDuplicated(assay_names) > 0L) {
        artifactCodecAbort(
            sprintf("%s contains missing or duplicate assay names", owner),
            "multischolar_metabolomics_codec_shape_mismatch"
        )
    }
    assay_names
}

artifactMetabNodePayloadKeys <- function(node) {
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
            artifactMetabNodePayloadKeys
        ), use.names = FALSE)
        return(c(direct, nested))
    }
    if (identical(node$node_type, "s4_exact_role")) {
        return(unlist(node$payload_keys, use.names = FALSE))
    }
    direct
}

artifactMetabNodeStorage <- function(node) {
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

artifactMetabRoleManifest <- function(slot_values, binding) {
    slot_roles <- lapply(names(slot_values), \(slot_name) {
        node <- slot_values[[slot_name]]
        list(
            slot_name = slot_name,
            storage = artifactMetabNodeStorage(node),
            role_id = if (workflowCapabilityScalarString(node$role_id)) {
                node$role_id
            } else {
                paste0("metabolomics.s4.slot.", slot_name)
            },
            payload_keys = unname(artifactMetabNodePayloadKeys(node))
        )
    })
    manifest <- list(
        schema = .artifactMetabRoleManifestSchema,
        schema_version = .artifactMetabRoleManifestVersion,
        slots = slot_roles,
        assays = binding$assays
    )
    artifactWorkflowStateAssertSafeMetadata(manifest, "metabolomics role manifest")
    manifest
}

artifactMetabAttachContract <- function(bundle, logical_key, owner, assay_names) {
    binding <- artifactMetabIdentityBinding(logical_key, owner, assay_names)
    bundle$metadata$identity_binding <- binding
    bundle$metadata$role_manifest <- artifactMetabRoleManifest(
        bundle$metadata$slot_values,
        binding
    )
    bundle$metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(bundle$metadata)
    )
    bundle
}

artifactMetabRoleClassSignature <- function(role_value) {
    member_classes <- if (is.list(role_value)) {
        lapply(role_value, \(item) unname(class(item)))
    } else {
        list()
    }
    list(
        class = unname(class(role_value)),
        names = if (is.null(names(role_value))) NULL else unname(names(role_value)),
        member_classes = member_classes
    )
}

artifactMetabRawChunks <- function(value) {
    if (length(value) == 0L) return(list(raw()))
    starts <- seq.int(1L, length(value), by = .artifactMetabSerializedChunkBytes)
    lapply(starts, \(start) {
        end <- min(length(value), start + .artifactMetabSerializedChunkBytes - 1L)
        value[start:end]
    })
}

artifactMetabEncodeSerializedRole <- function(
    role_value,
    state,
    slot_name,
    role_id
) {
    serialized <- serialize(
        role_value,
        connection = NULL,
        ascii = FALSE,
        version = .artifactMetabSerializationVersion
    )
    chunks <- artifactMetabRawChunks(serialized)
    table <- data.frame(
        chunk_index = seq_along(chunks),
        serialized_chunk = vapply(chunks, jsonlite::base64_enc, character(1)),
        stringsAsFactors = FALSE
    )
    encoded <- encodeArtifactTable(
        table,
        stable_key = "chunk_index",
        owner = paste("metabolomics immutable role", role_id)
    )
    list(
        node_type = "immutable_r_role",
        codec = list(
            id = .artifactMetabSerializedRoleCodec,
            version = .artifactMetabSerializedRoleCodecVersion
        ),
        role_id = role_id,
        slot_name = slot_name,
        payload_key = artifactAddRectangularPayload(state, encoded),
        serialization_version = .artifactMetabSerializationVersion,
        byte_count = as.numeric(length(serialized)),
        byte_digest = digest::digest(
            serialized,
            algo = .artifactHashAlgorithm,
            serialize = FALSE
        ),
        class_signature = artifactMetabRoleClassSignature(role_value)
    )
}

artifactMetabRemapPayloadKeys <- function(value, key_map) {
    if (!is.list(value)) return(value)
    remapped <- lapply(value, artifactMetabRemapPayloadKeys, key_map = key_map)
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

artifactMetabMergeNestedBundle <- function(state, bundle) {
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
    metadata <- artifactMetabRemapPayloadKeys(bundle$metadata, key_map)
    metadata$payloads <- state$payload_metadata[new_keys]
    metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(metadata)
    )
    list(metadata = metadata, payload_keys = unname(new_keys))
}

artifactMetabEncodeNestedS4 <- function(value, state, logical_key, inline_limit_bytes) {
    nested <- dehydrateMetabolomicsS4Artifact(
        value,
        logical_key,
        inline_limit_bytes = inline_limit_bytes
    )
    merged <- artifactMetabMergeNestedBundle(state, nested)
    list(
        node_type = "s4_exact_role",
        codec = list(
            id = .artifactMetabNestedS4Codec,
            version = .artifactMetabNestedS4CodecVersion
        ),
        role_id = "metabolomics.da.embedded_assay",
        slot_name = "theObject",
        metadata = merged$metadata,
        payload_keys = merged$payload_keys
    )
}

artifactMetabDaBundle <- function(value, codec, logical_key, inline_limit_bytes) {
    state <- new.env(parent = emptyenv())
    state$payloads <- list()
    state$payload_metadata <- list()
    slot_names <- methods::slotNames(value)
    slot_values <- lapply(slot_names, \(slot_name) {
        slot_value <- methods::slot(value, slot_name)
        if (identical(slot_name, "theObject")) {
            return(artifactMetabEncodeNestedS4(
                slot_value,
                state,
                logical_key,
                inline_limit_bytes
            ))
        }
        role_id <- if (slot_name %in% names(.artifactMetabSerializedSlots)) {
            .artifactMetabSerializedSlots[[slot_name]]
        } else {
            NULL
        }
        if (!is.null(role_id)) {
            return(artifactMetabEncodeSerializedRole(
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
    assay_node <- slot_values$theObject$metadata$slot_values$metabolite_data
    artifactMetabAttachContract(
        bundle,
        logical_key,
        class(value)[[1L]],
        artifactMetabAssayNames(assay_node, "embedded metabolite_data")
    )
}

dehydrateMetabolomicsS4Artifact <- function(
    value,
    logical_key,
    inline_limit_bytes = .artifactInlineLimitBytes
) {
    if (!isS4(value) || length(class(value)) != 1L) {
        artifactCodecAbort(
            "metabolomics codec requires one supported S4 object",
            "multischolar_missing_exact_s4_codec"
        )
    }
    logical_key <- artifactMetabLogicalKey(logical_key)
    class_name <- class(value)[[1L]]
    descriptor <- artifactMetabCodecDescriptor(class_name)
    codec <- artifactMetabolomicsCodecDeclarations()[[descriptor$id]]
    validity <- methods::validObject(value, test = TRUE)
    if (!identical(validity, TRUE)) {
        artifactCodecAbort(
            sprintf("metabolomics S4 object is invalid: %s", paste(validity, collapse = "; ")),
            "multischolar_invalid_metabolomics_s4_object",
            owner = class_name
        )
    }
    if (identical(class_name, "MetaboliteAssayData")) {
        bundle <- dehydrateExactS4Artifact(value, codec, inline_limit_bytes)
        assay_names <- names(methods::slot(value, "metabolite_data"))
        return(artifactMetabAttachContract(
            bundle,
            logical_key,
            class_name,
            assay_names
        ))
    }
    artifactMetabDaBundle(value, codec, logical_key, inline_limit_bytes)
}

artifactMetabPayloadExists <- function(key, payloads, payload_metadata) {
    workflowCapabilityScalarString(key) && key %in% names(payloads) &&
        key %in% names(payload_metadata)
}

artifactMetabPreflightNode <- function(node, payloads, payload_metadata) {
    if (!is.list(node) || !workflowCapabilityScalarString(node$node_type)) {
        artifactCodecAbort(
            "metabolomics S4 bundle contains a malformed value node",
            "multischolar_metabolomics_codec_shape_mismatch"
        )
    }
    payload_types <- c(
        "rectangular", "atomic_rectangular", "nested_rectangular",
        "immutable_r_role"
    )
    if (node$node_type %in% payload_types &&
        !artifactMetabPayloadExists(node$payload_key, payloads, payload_metadata)) {
        artifactCodecAbort(
            "metabolomics S4 value node references a missing payload",
            "multischolar_invalid_artifact_payload"
        )
    }
    if (identical(node$node_type, "nested_rectangular") && !identical(
        node$codec,
        list(id = .artifactNestedNodeCodec, version = .artifactNestedNodeCodecVersion)
    )) {
        artifactCodecAbort(
            "nested metabolomics metadata codec version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    if (identical(node$node_type, "immutable_r_role") && !identical(
        node$codec,
        list(
            id = .artifactMetabSerializedRoleCodec,
            version = .artifactMetabSerializedRoleCodecVersion
        )
    )) {
        artifactCodecAbort(
            "metabolomics immutable role codec version is unsupported",
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
                .artifactMetabSerializationVersion
            ) && length(node$byte_count) == 1L && is.numeric(node$byte_count) &&
            !is.na(node$byte_count) && is.finite(node$byte_count) &&
            node$byte_count >= 0 && node$byte_count == floor(node$byte_count)
        if (!isTRUE(valid)) {
            artifactCodecAbort(
                "metabolomics immutable role descriptor is malformed",
                "multischolar_metabolomics_codec_role_mismatch"
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
            id = .artifactMetabNestedS4Codec,
            version = .artifactMetabNestedS4CodecVersion
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
                "MetaboliteAssayData"
            ) && identical(
                nested_metadata$codec,
                artifactMetabCodecDescriptor("MetaboliteAssayData")
            ) && identical(names(nested_metadata$payloads), keys)
        if (!identical(names(node), required) ||
            !workflowCapabilityScalarString(node$role_id) ||
            !workflowCapabilityScalarString(node$slot_name) ||
            !identical(node$codec, expected) || !isTRUE(nested_version) ||
            !all(vapply(
            keys,
            artifactMetabPayloadExists,
            logical(1),
            payloads = payloads,
            payload_metadata = payload_metadata
        ))) {
            artifactCodecAbort(
                "nested metabolomics S4 role is malformed or unsupported",
                "multischolar_unsupported_artifact_codec_version"
            )
        }
        lapply(
            node$metadata$slot_values,
            artifactMetabPreflightNode,
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
            artifactMetabPreflightNode,
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
            "metabolomics S4 value node type is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    invisible(node)
}

artifactMetabSerializedRaw <- function(node, payloads, payload_metadata) {
    table <- decodeArtifactRectangular(
        payloads[[node$payload_key]],
        payload_metadata[[node$payload_key]]
    )
    expected_names <- c("chunk_index", "serialized_chunk")
    if (!identical(names(table), expected_names) || nrow(table) < 1L ||
        !identical(table$chunk_index, seq_len(nrow(table))) ||
        anyNA(table$serialized_chunk)) {
        artifactCodecAbort(
            "metabolomics immutable role payload is malformed",
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
                    "metabolomics immutable role '%s' payload digest does not match",
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

artifactMetabDecodeSerializedRole <- function(node, payloads, payload_metadata) {
    serialized <- artifactMetabSerializedRaw(node, payloads, payload_metadata)
    value <- tryCatch(
        unserialize(serialized),
        error = \(error) artifactCodecAbort(
            "metabolomics immutable R role could not be hydrated",
            "multischolar_invalid_artifact_payload",
            parent = error,
            owner = node$role_id
        )
    )
    if (!identical(artifactMetabRoleClassSignature(value), node$class_signature)) {
        artifactCodecAbort(
            "metabolomics immutable R role class signature changed",
            "multischolar_inexact_s4_hydration",
            owner = node$role_id
        )
    }
    value
}

artifactMetabNestedBundle <- function(node, payloads) {
    keys <- unlist(node$payload_keys, use.names = FALSE)
    structure(
        list(metadata = node$metadata, payloads = payloads[keys]),
        class = c("MultiScholaRArtifactS4Bundle", "list")
    )
}

artifactMetabValidateNodePayload <- function(node, payloads, payload_metadata) {
    if (identical(node$node_type, "immutable_r_role")) {
        artifactMetabDecodeSerializedRole(node, payloads, payload_metadata)
    } else if (identical(node$node_type, "s4_exact_role")) {
        validateMetabolomicsS4Bundle(artifactMetabNestedBundle(node, payloads))
    } else if (identical(node$node_type, "list")) {
        children <- node$values
        if (!is.null(node$attributes)) {
            children <- c(children, list(node$attributes))
        }
        lapply(
            children,
            artifactMetabValidateNodePayload,
            payloads = payloads,
            payload_metadata = payload_metadata
        )
    } else {
        artifactDecodeValueNode(node, payloads, payload_metadata)
    }
    invisible(node)
}

artifactMetabValidateSpecialRoles <- function(metadata) {
    if (!identical(metadata$class_name, "MetabolomicsDifferentialAbundanceResults")) {
        return(invisible(metadata))
    }
    embedded <- metadata$slot_values$theObject
    valid_embedded <- identical(embedded$node_type, "s4_exact_role") &&
        identical(embedded$role_id, "metabolomics.da.embedded_assay") &&
        identical(embedded$slot_name, "theObject")
    role_slots <- names(.artifactMetabSerializedSlots)
    valid_roles <- vapply(role_slots, \(slot_name) {
        node <- metadata$slot_values[[slot_name]]
        identical(node$node_type, "immutable_r_role") &&
            identical(node$role_id, .artifactMetabSerializedSlots[[slot_name]]) &&
            identical(node$slot_name, slot_name) &&
            identical(node$serialization_version, .artifactMetabSerializationVersion)
    }, logical(1))
    if (!isTRUE(valid_embedded) || !all(valid_roles)) {
        artifactCodecAbort(
            "metabolomics DA model, plot, or embedded assay role is inconsistent",
            "multischolar_metabolomics_codec_role_mismatch"
        )
    }
    invisible(metadata)
}

artifactMetabBundleAssayNames <- function(metadata) {
    node <- if (identical(metadata$class_name, "MetaboliteAssayData")) {
        metadata$slot_values$metabolite_data
    } else {
        metadata$slot_values$theObject$metadata$slot_values$metabolite_data
    }
    artifactMetabAssayNames(node)
}

validateMetabolomicsS4Bundle <- function(bundle, validate_payloads = TRUE) {
    if (!is.logical(validate_payloads) || length(validate_payloads) != 1L ||
        is.na(validate_payloads)) {
        artifactCodecAbort(
            "metabolomics payload validation flag must be TRUE or FALSE",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    if (!is.list(bundle) || !is.list(bundle$metadata) ||
        !workflowCapabilityScalarString(bundle$metadata$class_name)) {
        artifactCodecAbort(
            "metabolomics S4 artifact bundle is malformed",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    descriptor <- artifactMetabCodecDescriptor(bundle$metadata$class_name)
    codec <- artifactMetabolomicsCodecDeclarations()[[descriptor$id]]
    bundle <- validateExactS4Bundle(bundle, codec)
    metadata <- bundle$metadata
    assay_names <- artifactMetabBundleAssayNames(metadata)
    binding <- artifactMetabValidateIdentityBinding(
        metadata$identity_binding,
        metadata$class_name,
        assay_names
    )
    expected_manifest <- artifactMetabRoleManifest(metadata$slot_values, binding)
    if (!identical(metadata$role_manifest, expected_manifest)) {
        artifactCodecAbort(
            "metabolomics S4 role manifest does not match its value nodes",
            "multischolar_metabolomics_codec_role_mismatch"
        )
    }
    artifactMetabValidateSpecialRoles(metadata)
    lapply(
        metadata$slot_values,
        artifactMetabPreflightNode,
        payloads = bundle$payloads,
        payload_metadata = metadata$payloads
    )
    if (isTRUE(validate_payloads)) {
        lapply(
            metadata$slot_values,
            artifactMetabValidateNodePayload,
            payloads = bundle$payloads,
            payload_metadata = metadata$payloads
        )
    }
    bundle
}

artifactMetabDecodeNode <- function(node, payloads, payload_metadata) {
    if (identical(node$node_type, "immutable_r_role")) {
        return(artifactMetabDecodeSerializedRole(node, payloads, payload_metadata))
    }
    if (identical(node$node_type, "s4_exact_role")) {
        return(hydrateMetabolomicsS4Artifact(
            artifactMetabNestedBundle(node, payloads)
        ))
    }
    artifactDecodeValueNode(node, payloads, payload_metadata)
}

hydrateMetabolomicsS4Artifact <- function(bundle) {
    bundle <- validateMetabolomicsS4Bundle(bundle, validate_payloads = FALSE)
    metadata <- bundle$metadata
    slot_values <- lapply(
        metadata$slot_values,
        artifactMetabDecodeNode,
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
                "hydrated metabolomics object is invalid: %s",
                paste(validity, collapse = "; ")
            ),
            "multischolar_invalid_hydrated_s4_object",
            owner = metadata$class_name
        )
    }
    object
}

encodeMetabolomicsS4Artifact <- dehydrateMetabolomicsS4Artifact
decodeMetabolomicsS4Artifact <- hydrateMetabolomicsS4Artifact
