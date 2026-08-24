# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactS4BundleSchema <- "multischolar.s4_bundle"
.artifactS4BundleSchemaVersion <- 1L
.artifactDiaCodecVersion <- 1L
.artifactCodecCatalogueSchema <- "multischolar.artifact_codec_catalogue"
.artifactCodecCatalogueVersion <- 1L
.artifactNestedNodeCodec <- "multischolar.r_nested_value"
.artifactNestedNodeCodecVersion <- 1L

artifactMetabolomicsCodecDeclarations <- function() {
    list(
        "multischolar.s4.metabolite_assay_data" = list(
            codec_id = "multischolar.s4.metabolite_assay_data",
            codec_version = 1L,
            class_name = "MetaboliteAssayData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        ),
        "multischolar.s4.metabolomics_da_results" = list(
            codec_id = "multischolar.s4.metabolomics_da_results",
            codec_version = 1L,
            class_name = "MetabolomicsDifferentialAbundanceResults",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    )
}

artifactDiaCodecDeclarations <- function() {
    artifactDiaWorkflowCodecDeclarations()
}

artifactDiaCodecDescriptor <- function(class_name) {
    declarations <- artifactDiaCodecDeclarations()
    matches <- Filter(function(value) {
        identical(value$class_name, class_name)
    }, declarations)
    if (length(matches) != 1L) {
        artifactCodecAbort(
            sprintf("S4 class '%s' has no certified DIA artifact codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    list(id = matches[[1L]]$codec_id, version = matches[[1L]]$codec_version)
}

artifactDiaBundleSemanticInput <- function(metadata) {
    input <- list(
        codec = metadata$codec,
        class_name = metadata$class_name,
        slot_names = metadata$slot_names,
        slot_values = unname(metadata$slot_values),
        payload_semantic_digests = unname(lapply(
            metadata$payloads,
            `[[`,
            "semantic_digest"
        ))
    )
    if (!is.null(metadata$identity_binding)) {
        input$identity_binding <- metadata$identity_binding
    }
    if (!is.null(metadata$role_manifest)) {
        input$role_manifest <- metadata$role_manifest
    }
    input
}

artifactInlineCharacter <- function(value) {
    missing <- is.na(value)
    value[missing] <- ""
    list(values = unname(value), missing = unname(which(missing)))
}

artifactDecodeInlineCharacter <- function(encoded) {
    value <- unlist(encoded$values, use.names = FALSE)
    value <- as.character(value)
    if (length(encoded$missing) > 0L) {
        value[as.integer(encoded$missing)] <- NA_character_
    }
    value
}

artifactInlineDouble <- function(value) {
    vapply(value, function(item) {
        if (is.na(item) && !is.nan(item)) {
            return("NA")
        }
        if (is.nan(item)) {
            return("NaN")
        }
        if (is.infinite(item) && item > 0) {
            return("+Inf")
        }
        if (is.infinite(item) && item < 0) {
            return("-Inf")
        }
        sprintf("%a", item)
    }, character(1))
}

artifactDecodeInlineDouble <- function(value) {
    vapply(value, function(item) {
        switch(item,
            `NA` = NA_real_,
            `NaN` = NaN,
            `+Inf` = Inf,
            `-Inf` = -Inf,
            as.numeric(item)
        )
    }, numeric(1))
}

artifactAtomicAttributes <- function(value, owner) {
    type <- artifactSupportedColumnType(value, owner)
    allowed <- switch(type,
        factor = c("names", "levels", "class", "units"),
        ordered = c("names", "levels", "class", "units"),
        Date = c("names", "class", "units"),
        POSIXct = c("names", "class", "tzone", "units"),
        integer64 = c("names", "class", "units"),
        c("names", "units")
    )
    unsupported <- setdiff(names(attributes(value)), allowed)
    if (length(unsupported) > 0L) {
        artifactCodecAbort(
            sprintf(
                "artifact value '%s' has unsupported attribute(s): %s",
                owner,
                paste(unsupported, collapse = ", ")
            ),
            "multischolar_artifact_externalization_required",
            owner = owner,
            attributes = unsupported
        )
    }
    list(
        r_type = type,
        names = if (is.null(names(value))) NULL else artifactInlineCharacter(names(value)),
        levels = if (type %in% c("factor", "ordered")) levels(value) else NULL,
        timezone = if (identical(type, "POSIXct")) {
            attr(value, "tzone", exact = TRUE)
        } else {
            NULL
        },
        units = artifactColumnUnits(value, owner)
    )
}

artifactEncodeInlineAtomic <- function(value, owner) {
    descriptor <- artifactAtomicAttributes(value, owner)
    type <- descriptor$r_type
    encoded_values <- switch(type,
        logical = ifelse(is.na(value), "NA", ifelse(value, "TRUE", "FALSE")),
        integer = ifelse(is.na(value), "NA", as.character(value)),
        integer64 = artifactInlineCharacter(as.character(value)),
        double = artifactInlineDouble(value),
        Date = artifactInlineDouble(unclass(value)),
        POSIXct = artifactInlineDouble(unclass(value)),
        character = artifactInlineCharacter(value),
        factor = artifactInlineCharacter(as.character(value)),
        ordered = artifactInlineCharacter(as.character(value))
    )
    if (!is.list(encoded_values)) encoded_values <- unname(encoded_values)
    list(
        node_type = "atomic_inline",
        descriptor = descriptor,
        values = encoded_values
    )
}

artifactDecodeInlineAtomic <- function(node) {
    descriptor <- node$descriptor
    type <- descriptor$r_type
    value <- switch(type,
        logical = c("NA" = NA, "TRUE" = TRUE, "FALSE" = FALSE)[unlist(node$values)],
        integer = {
            values <- unlist(node$values, use.names = FALSE)
            as.integer(ifelse(values == "NA", NA_character_, values))
        },
        integer64 = {
            if (!requireNamespace("bit64", quietly = TRUE)) {
                artifactCodecAbort(
                    "bit64 is required to hydrate an integer64 artifact value",
                    "multischolar_missing_artifact_codec_dependency"
                )
            }
            bit64::as.integer64(artifactDecodeInlineCharacter(node$values))
        },
        double = artifactDecodeInlineDouble(unlist(node$values, use.names = FALSE)),
        Date = structure(
            artifactDecodeInlineDouble(unlist(node$values, use.names = FALSE)),
            class = "Date"
        ),
        POSIXct = structure(
            artifactDecodeInlineDouble(unlist(node$values, use.names = FALSE)),
            class = c("POSIXct", "POSIXt"),
            tzone = descriptor$timezone
        ),
        character = artifactDecodeInlineCharacter(node$values),
        factor = factor(
            artifactDecodeInlineCharacter(node$values),
            levels = descriptor$levels
        ),
        ordered = ordered(
            artifactDecodeInlineCharacter(node$values),
            levels = descriptor$levels
        )
    )
    if (is.null(descriptor$names)) {
        value <- unname(value)
    } else {
        names(value) <- artifactDecodeInlineCharacter(descriptor$names)
    }
    if (!is.null(descriptor$units)) attr(value, "units") <- descriptor$units
    value
}

artifactAddRectangularPayload <- function(state, encoded) {
    key <- sprintf("payload_%08d", length(state$payloads) + 1L)
    state$payloads[[key]] <- encoded$payload
    state$payload_metadata[[key]] <- encoded$metadata
    key
}

artifactAddNestedNodePayload <- function(state, node, owner) {
    artifactWorkflowStateAssertSafeMetadata(node, owner)
    serialized <- as.character(jsonlite::serializeJSON(node, digits = NA))
    encoded <- encodeArtifactTable(
        data.frame(serialized_node = serialized, stringsAsFactors = FALSE),
        owner = paste(owner, "nested value")
    )
    artifactAddRectangularPayload(state, encoded)
}

artifactNodeOwner <- function(owner, names, index) {
    label <- if (!is.null(names) && nzchar(names[[index]])) names[[index]] else index
    paste0(owner, "[[", label, "]]")
}

artifactEncodeValueNode <- function(value, state, owner, inline_limit_bytes) {
    if (is.null(value)) {
        return(list(node_type = "null"))
    }
    if (isS4(value)) {
        artifactCodecAbort(
            sprintf("nested S4 value '%s' requires its own exact codec", owner),
            "multischolar_artifact_externalization_required",
            owner = owner
        )
    }
    if (is.data.frame(value)) {
        encoded <- encodeArtifactTable(value, owner = owner)
        return(list(
            node_type = "rectangular",
            payload_key = artifactAddRectangularPayload(state, encoded)
        ))
    }
    if (is.matrix(value)) {
        encoded <- encodeArtifactMatrix(value, owner = owner)
        return(list(
            node_type = "rectangular",
            payload_key = artifactAddRectangularPayload(state, encoded)
        ))
    }
    if (is.raw(value)) {
        return(list(
            node_type = "raw_inline",
            value = jsonlite::base64_enc(value),
            names = if (is.null(names(value))) NULL else artifactInlineCharacter(names(value))
        ))
    }
    if (is.atomic(value)) {
        if (as.numeric(object.size(value)) > inline_limit_bytes) {
            frame <- structure(
                list(value),
                names = "value",
                class = "data.frame",
                row.names = c(NA_integer_, -length(value))
            )
            encoded <- encodeArtifactTable(frame, owner = owner)
            return(list(
                node_type = "atomic_rectangular",
                payload_key = artifactAddRectangularPayload(state, encoded),
                names = if (is.null(names(value))) {
                    NULL
                } else {
                    artifactInlineCharacter(names(value))
                }
            ))
        }
        return(artifactEncodeInlineAtomic(value, owner))
    }
    if (is.list(value)) {
        value_names <- names(value)
        nodes <- lapply(seq_along(value), function(index) {
            artifactEncodeValueNode(
                value[[index]],
                state,
                artifactNodeOwner(owner, value_names, index),
                inline_limit_bytes
            )
        })
        extra_attributes <- attributes(value)
        extra_attributes[c("names", "class")] <- NULL
        attribute_node <- if (length(extra_attributes) == 0L) {
            NULL
        } else {
            artifactEncodeValueNode(
                extra_attributes,
                state,
                paste0(owner, "@attributes"),
                inline_limit_bytes
            )
        }
        node <- list(
            node_type = "list",
            names = if (is.null(value_names)) NULL else artifactInlineCharacter(value_names),
            class = unname(class(value)),
            values = nodes,
            attributes = attribute_node
        )
        if (as.numeric(object.size(node)) > inline_limit_bytes) {
            return(list(
                node_type = "nested_rectangular",
                codec = list(
                    id = .artifactNestedNodeCodec,
                    version = .artifactNestedNodeCodecVersion
                ),
                payload_key = artifactAddNestedNodePayload(state, node, owner)
            ))
        }
        return(node)
    }
    artifactCodecAbort(
        sprintf("artifact value '%s' requires a declared external codec", owner),
        "multischolar_artifact_externalization_required",
        owner = owner,
        value_type = typeof(value)
    )
}

artifactDecodeNestedNode <- function(
    node,
    payloads,
    payload_metadata,
    visited_payload_keys
) {
    expected_codec <- list(
        id = .artifactNestedNodeCodec,
        version = .artifactNestedNodeCodecVersion
    )
    if (!identical(node$codec, expected_codec) ||
        !workflowCapabilityScalarString(node$payload_key)) {
        artifactCodecAbort(
            "nested artifact value codec is malformed or unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    if (node$payload_key %in% visited_payload_keys) {
        artifactCodecAbort(
            "nested artifact value payload lineage contains a cycle",
            "multischolar_invalid_artifact_payload"
        )
    }
    frame <- decodeArtifactRectangular(
        payloads[[node$payload_key]],
        payload_metadata[[node$payload_key]]
    )
    if (!identical(names(frame), "serialized_node") || nrow(frame) != 1L ||
        !workflowCapabilityScalarString(frame$serialized_node[[1L]])) {
        artifactCodecAbort(
            "nested artifact value payload is malformed",
            "multischolar_invalid_artifact_payload"
        )
    }
    decoded <- tryCatch(
        jsonlite::unserializeJSON(frame$serialized_node[[1L]]),
        error = \(error) artifactCodecAbort(
            "nested artifact value payload could not be decoded",
            "multischolar_invalid_artifact_payload",
            parent = error
        )
    )
    if (!is.list(decoded) || !workflowCapabilityScalarString(decoded$node_type)) {
        artifactCodecAbort(
            "nested artifact value payload does not contain one value node",
            "multischolar_invalid_artifact_payload"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(decoded, "nested artifact value")
    artifactDecodeValueNode(
        decoded,
        payloads,
        payload_metadata,
        c(visited_payload_keys, node$payload_key)
    )
}

artifactDecodeValueNode <- function(
    node,
    payloads,
    payload_metadata,
    visited_payload_keys = character()
) {
    switch(node$node_type,
        null = NULL,
        raw_inline = {
            value <- jsonlite::base64_dec(node$value)
            if (!is.null(node$names)) names(value) <- artifactDecodeInlineCharacter(node$names)
            value
        },
        atomic_inline = artifactDecodeInlineAtomic(node),
        rectangular = decodeArtifactRectangular(
            payloads[[node$payload_key]],
            payload_metadata[[node$payload_key]]
        ),
        atomic_rectangular = {
            value <- decodeArtifactRectangular(
                payloads[[node$payload_key]],
                payload_metadata[[node$payload_key]]
            )[[1L]]
            if (!is.null(node$names)) names(value) <- artifactDecodeInlineCharacter(node$names)
            value
        },
        nested_rectangular = artifactDecodeNestedNode(
            node,
            payloads,
            payload_metadata,
            visited_payload_keys
        ),
        list = {
            value <- lapply(
                node$values,
                artifactDecodeValueNode,
                payloads = payloads,
                payload_metadata = payload_metadata,
                visited_payload_keys = visited_payload_keys
            )
            if (!is.null(node$names)) names(value) <- artifactDecodeInlineCharacter(node$names)
            if (!is.null(node$attributes)) {
                extra <- artifactDecodeValueNode(
                    node$attributes,
                    payloads,
                    payload_metadata,
                    visited_payload_keys
                )
                attributes(value) <- c(attributes(value), extra)
            }
            if (length(node$class) > 0L && !identical(node$class, "list")) {
                class(value) <- node$class
            }
            value
        },
        artifactCodecAbort(
            "artifact value node has an unsupported version or type",
            "multischolar_unsupported_artifact_codec_version"
        )
    )
}

artifactExactS4CodecDeclaration <- function(codec) {
    artifactDescriptorAssertDataOnly(codec, "codec declaration")
    if (!is.list(codec) || !workflowCapabilityScalarString(codec$codec_id)) {
        artifactCodecAbort(
            "exact S4 codec declaration is malformed",
            "multischolar_invalid_artifact_codec_catalogue"
        )
    }
    codec <- artifactDescriptorValidateCodec(codec, codec$codec_id)
    supported <- identical(codec$codec_version, 1L) &&
        identical(codec$payload_schema_id, "multischolar.rectangular") &&
        identical(codec$payload_schema_version, 1L)
    if (!isTRUE(supported)) {
        artifactCodecAbort(
            "exact S4 codec or payload schema version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    codec
}

dehydrateExactS4Artifact <- function(
    value,
    codec,
    inline_limit_bytes = .artifactInlineLimitBytes
) {
    codec <- artifactExactS4CodecDeclaration(codec)
    if (!isS4(value) || length(class(value)) != 1L) {
        artifactCodecAbort(
            "exact S4 codec requires one supported S4 object",
            "multischolar_missing_exact_s4_codec"
        )
    }
    class_name <- class(value)[[1L]]
    if (!identical(class_name, codec$class_name)) {
        artifactCodecAbort(
            sprintf("S4 class '%s' does not match its exact artifact codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    descriptor <- list(id = codec$codec_id, version = codec$codec_version)
    slot_names <- methods::slotNames(value)
    state <- new.env(parent = emptyenv())
    state$payloads <- list()
    state$payload_metadata <- list()
    slot_values <- lapply(slot_names, function(slot_name) {
        artifactEncodeValueNode(
            methods::slot(value, slot_name),
            state,
            paste0(class_name, "@", slot_name),
            inline_limit_bytes
        )
    })
    names(slot_values) <- slot_names
    metadata <- list(
        schema = .artifactS4BundleSchema,
        schema_version = .artifactS4BundleSchemaVersion,
        codec = descriptor,
        class_name = class_name,
        slot_names = slot_names,
        slot_values = slot_values,
        payloads = state$payload_metadata,
        semantic_digest = NULL,
        created_at = artifactRefUtcNow()
    )
    metadata$semantic_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(metadata)
    )
    structure(
        list(metadata = metadata, payloads = state$payloads),
        class = c("MultiScholaRArtifactS4Bundle", "list")
    )
}

validateExactS4Bundle <- function(bundle, codec) {
    codec <- artifactExactS4CodecDeclaration(codec)
    if (!is.list(bundle) || !is.list(bundle$metadata) || !is.list(bundle$payloads)) {
        artifactCodecAbort(
            "exact S4 artifact bundle is malformed",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    metadata <- bundle$metadata
    if (!identical(metadata$schema, .artifactS4BundleSchema) ||
        !identical(metadata$schema_version, .artifactS4BundleSchemaVersion)) {
        artifactCodecAbort(
            "exact S4 artifact bundle schema version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    descriptor <- list(id = codec$codec_id, version = codec$codec_version)
    if (!identical(metadata$codec, descriptor) ||
        !identical(metadata$class_name, codec$class_name)) {
        artifactCodecAbort(
            "exact S4 artifact codec version or class is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    expected_slots <- methods::slotNames(metadata$class_name)
    if (!identical(metadata$slot_names, expected_slots) ||
        !identical(names(metadata$slot_values), expected_slots) ||
        !identical(names(metadata$payloads), names(bundle$payloads))) {
        artifactCodecAbort(
            "exact S4 artifact bundle shape does not match the class contract",
            "multischolar_artifact_shape_mismatch"
        )
    }
    expected_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(metadata)
    )
    if (!identical(metadata$semantic_digest, expected_digest)) {
        artifactCodecAbort(
            "exact S4 artifact semantic lineage digest does not match its metadata",
            "multischolar_artifact_semantic_digest_mismatch"
        )
    }
    bundle
}

hydrateExactS4Artifact <- function(bundle, codec) {
    bundle <- validateExactS4Bundle(bundle, codec)
    metadata <- bundle$metadata
    slot_values <- lapply(metadata$slot_values, artifactDecodeValueNode,
        payloads = bundle$payloads,
        payload_metadata = metadata$payloads
    )
    object <- do.call(
        methods::new,
        c(list(Class = metadata$class_name), slot_values)
    )
    exact_slots <- vapply(metadata$slot_names, function(slot_name) {
        identical(methods::slot(object, slot_name), slot_values[[slot_name]])
    }, logical(1))
    if (!all(exact_slots)) {
        artifactCodecAbort(
            sprintf(
                "exact S4 hydration changed slot(s): %s",
                paste(metadata$slot_names[!exact_slots], collapse = ", ")
            ),
            "multischolar_inexact_s4_hydration",
            owner = metadata$class_name
        )
    }
    validity <- methods::validObject(object, test = TRUE)
    if (!identical(validity, TRUE)) {
        artifactCodecAbort(
            sprintf("hydrated S4 object is invalid: %s", paste(validity, collapse = "; ")),
            "multischolar_invalid_hydrated_s4_object",
            owner = metadata$class_name
        )
    }
    object
}

newArtifactCodecCatalogue <- function(codecs) {
    if (!is.list(codecs) || length(codecs) == 0L) {
        artifactCodecAbort(
            "artifact codec catalogue requires one or more declarations",
            "multischolar_invalid_artifact_codec_catalogue"
        )
    }
    codecs <- lapply(codecs, artifactExactS4CodecDeclaration)
    ids <- vapply(codecs, `[[`, character(1), "codec_id")
    if (anyDuplicated(ids) > 0L) {
        artifactCodecAbort(
            "artifact codec catalogue contains duplicate codec IDs",
            "multischolar_duplicate_artifact_codec"
        )
    }
    catalogue <- new.env(parent = emptyenv())
    catalogue$schema <- .artifactCodecCatalogueSchema
    catalogue$schema_version <- .artifactCodecCatalogueVersion
    catalogue$codecs <- stats::setNames(codecs, ids)
    class(catalogue) <- c("MultiScholaRArtifactCodecCatalogue", "environment")
    lockEnvironment(catalogue, bindings = TRUE)
    catalogue
}

validateArtifactCodecCatalogue <- function(catalogue) {
    valid <- inherits(catalogue, "MultiScholaRArtifactCodecCatalogue") &&
        is.environment(catalogue) && environmentIsLocked(catalogue) &&
        identical(catalogue$schema, .artifactCodecCatalogueSchema) &&
        identical(catalogue$schema_version, .artifactCodecCatalogueVersion)
    if (!isTRUE(valid)) {
        artifactCodecAbort(
            "artifact codec catalogue is invalid or mutable",
            "multischolar_invalid_artifact_codec_catalogue"
        )
    }
    catalogue
}

artifactCodecForClass <- function(catalogue, codec_ids, class_name) {
    catalogue <- validateArtifactCodecCatalogue(catalogue)
    unknown <- setdiff(codec_ids, names(catalogue$codecs))
    if (length(unknown) > 0L) {
        artifactCodecAbort(
            "workflow descriptor references an unavailable exact S4 codec",
            "multischolar_missing_exact_s4_codec",
            codec_ids = unknown
        )
    }
    matches <- Filter(function(codec) {
        identical(codec$class_name, class_name)
    }, catalogue$codecs[codec_ids])
    if (length(matches) != 1L) {
        artifactCodecAbort(
            sprintf("S4 class '%s' has no unambiguous exact codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    matches[[1L]]
}

artifactCodecForBundle <- function(catalogue, codec_ids, bundle) {
    if (!is.list(bundle) || !is.list(bundle$metadata) ||
        !workflowCapabilityScalarString(bundle$metadata$class_name)) {
        artifactCodecAbort(
            "exact S4 bundle cannot be resolved to a codec",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    artifactCodecForClass(catalogue, codec_ids, bundle$metadata$class_name)
}

artifactCodecAdapter <- function(descriptor, catalogue) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    validateArtifactCodecCatalogue(catalogue)
    codec_ids <- names(descriptor$codecs)
    for (codec_id in codec_ids) {
        declared <- descriptor$codecs[[codec_id]]
        available <- catalogue$codecs[[codec_id]]
        if (is.null(available) || !identical(declared, available)) {
            artifactCodecAbort(
                "workflow and codec catalogues disagree on an exact codec",
                "multischolar_incompatible_artifact_codec_catalogue",
                codec_id = codec_id
            )
        }
    }
    list(
        dehydrate = function(value) {
            codec <- artifactCodecForClass(catalogue, codec_ids, class(value)[[1L]])
            dehydrateExactS4Artifact(value, codec)
        },
        validate = function(bundle) {
            codec <- artifactCodecForBundle(catalogue, codec_ids, bundle)
            validateExactS4Bundle(bundle, codec)
        },
        hydrate = function(bundle) {
            codec <- artifactCodecForBundle(catalogue, codec_ids, bundle)
            hydrateExactS4Artifact(bundle, codec)
        }
    )
}

dehydrateDiaS4Artifact <- function(value, inline_limit_bytes = .artifactInlineLimitBytes) {
    declarations <- artifactDiaCodecDeclarations()
    descriptor <- artifactDiaCodecDescriptor(class(value)[[1L]])
    dehydrateExactS4Artifact(
        value,
        declarations[[descriptor$id]],
        inline_limit_bytes = inline_limit_bytes
    )
}

validateDiaS4Bundle <- function(bundle) {
    if (!is.list(bundle) || !is.list(bundle$metadata) ||
        !workflowCapabilityScalarString(bundle$metadata$class_name)) {
        artifactCodecAbort(
            "DIA S4 artifact bundle is malformed",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    descriptor <- artifactDiaCodecDescriptor(bundle$metadata$class_name)
    validateExactS4Bundle(bundle, artifactDiaCodecDeclarations()[[descriptor$id]])
}

hydrateDiaS4Artifact <- function(bundle) {
    bundle <- validateDiaS4Bundle(bundle)
    descriptor <- artifactDiaCodecDescriptor(bundle$metadata$class_name)
    hydrateExactS4Artifact(bundle, artifactDiaCodecDeclarations()[[descriptor$id]])
}

.ARTIFACT_S4_CODEC_CATALOGUE <- newArtifactCodecCatalogue(
    c(
        artifactDiaCodecDeclarations(),
        artifactProteomicsNonDiaCodecDeclarations(),
        artifactMetabolomicsCodecDeclarations()
    )
)

artifactS4CodecCatalogue <- function() {
    .ARTIFACT_S4_CODEC_CATALOGUE
}

encodeDiaS4Artifact <- dehydrateDiaS4Artifact
decodeDiaS4Artifact <- hydrateDiaS4Artifact
