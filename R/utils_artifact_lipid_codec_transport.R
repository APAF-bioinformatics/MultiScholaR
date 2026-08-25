# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactLipidTransportVector <- function(value, mode) {
    if (is.null(value)) return(NULL)
    flattened <- unlist(value, use.names = FALSE)
    switch(mode,
        character = as.character(flattened),
        integer = as.integer(flattened),
        numeric = as.numeric(flattened)
    )
}

artifactLipidNormalizeInlineCharacter <- function(value) {
    if (is.null(value)) return(NULL)
    value$values <- artifactLipidTransportVector(value$values, "character")
    value$missing <- artifactLipidTransportVector(value$missing, "integer")
    value
}

artifactLipidNormalizeNodeTransport <- function(node) {
    if (!is.list(node) || !workflowCapabilityScalarString(node$node_type)) {
        return(node)
    }
    if (is.list(node$codec) && !is.null(node$codec$version)) {
        node$codec$version <- as.integer(node$codec$version)
    }
    if (is.list(node$descriptor)) {
        node$descriptor["names"] <- list(
            artifactLipidNormalizeInlineCharacter(node$descriptor$names)
        )
    }
    if (!is.null(node$names)) {
        node$names <- artifactLipidNormalizeInlineCharacter(node$names)
    }
    if (identical(node$node_type, "atomic_inline") &&
        is.list(node$descriptor) &&
        workflowCapabilityScalarString(node$descriptor$r_type) &&
        node$descriptor$r_type %in% c(
            "character", "factor", "ordered", "integer64"
        )) {
        node$values <- artifactLipidNormalizeInlineCharacter(node$values)
    }
    if (identical(node$node_type, "list")) {
        node$values <- lapply(node$values, artifactLipidNormalizeNodeTransport)
        if (!is.null(node$attributes)) {
            node$attributes <- artifactLipidNormalizeNodeTransport(
                node$attributes
            )
        }
    }
    if (identical(node$node_type, "immutable_r_role")) {
        node$serialization_version <- as.integer(node$serialization_version)
        node$byte_count <- as.numeric(node$byte_count)
        node$class_signature$class <- artifactLipidTransportVector(
            node$class_signature$class,
            "character"
        )
        if (!is.null(node$class_signature$names)) {
            node$class_signature$names <- artifactLipidTransportVector(
                node$class_signature$names,
                "character"
            )
        }
        node$class_signature$member_classes <- lapply(
            node$class_signature$member_classes,
            artifactLipidTransportVector,
            mode = "character"
        )
    }
    if (identical(node$node_type, "s4_exact_role")) {
        node$payload_keys <- artifactLipidTransportVector(
            node$payload_keys,
            "character"
        )
        node$metadata <- artifactLipidNormalizeMetadataTransport(node$metadata)
    }
    node
}

artifactLipidNormalizeAssayBindings <- function(assays) {
    lapply(assays, \(assay) {
        assay$position <- as.integer(assay$position)
        assay
    })
}

artifactLipidNormalizeMetadataTransport <- function(metadata) {
    if (!is.list(metadata)) return(metadata)
    metadata$schema_version <- as.integer(metadata$schema_version)
    metadata$codec$version <- as.integer(metadata$codec$version)
    metadata$slot_names <- artifactLipidTransportVector(
        metadata$slot_names,
        "character"
    )
    metadata$slot_values <- lapply(
        metadata$slot_values,
        artifactLipidNormalizeNodeTransport
    )
    if (is.list(metadata$identity_binding)) {
        metadata$identity_binding$schema_version <- as.integer(
            metadata$identity_binding$schema_version
        )
        metadata$identity_binding$assays <- artifactLipidNormalizeAssayBindings(
            metadata$identity_binding$assays
        )
    }
    if (is.list(metadata$role_manifest)) {
        metadata$role_manifest$schema_version <- as.integer(
            metadata$role_manifest$schema_version
        )
        metadata$role_manifest$slots <- lapply(
            metadata$role_manifest$slots,
            \(slot) {
                slot$payload_keys <- artifactLipidTransportVector(
                    slot$payload_keys,
                    "character"
                )
                slot
            }
        )
        metadata$role_manifest$assays <- artifactLipidNormalizeAssayBindings(
            metadata$role_manifest$assays
        )
    }
    metadata
}
