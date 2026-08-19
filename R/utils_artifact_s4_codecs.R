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

artifactDiaCodecDescriptor <- function(class_name) {
    descriptors <- list(
        PeptideQuantitativeData = list(
            id = "multischolar.s4.peptide_quantitative_data.diann",
            version = .artifactDiaCodecVersion
        ),
        ProteinQuantitativeData = list(
            id = "multischolar.s4.protein_quantitative_data.diann",
            version = .artifactDiaCodecVersion
        )
    )
    descriptor <- descriptors[[class_name]]
    if (is.null(descriptor)) {
        artifactCodecAbort(
            sprintf("S4 class '%s' has no certified DIA artifact codec", class_name),
            "multischolar_missing_exact_s4_codec",
            owner = class_name
        )
    }
    descriptor
}

artifactDiaBundleSemanticInput <- function(metadata) {
    list(
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
            artifactCodecAbort(
                sprintf(
                    "nested owner '%s' exceeds the inline limit and requires an external codec",
                    owner
                ),
                "multischolar_artifact_externalization_required",
                owner = owner,
                inline_limit_bytes = inline_limit_bytes
            )
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

artifactDecodeValueNode <- function(node, payloads, payload_metadata) {
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
        list = {
            value <- lapply(node$values, artifactDecodeValueNode,
                payloads = payloads,
                payload_metadata = payload_metadata
            )
            if (!is.null(node$names)) names(value) <- artifactDecodeInlineCharacter(node$names)
            if (!is.null(node$attributes)) {
                extra <- artifactDecodeValueNode(
                    node$attributes,
                    payloads,
                    payload_metadata
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

dehydrateDiaS4Artifact <- function(
    value,
    inline_limit_bytes = .artifactInlineLimitBytes
) {
    if (!isS4(value) || length(class(value)) != 1L) {
        artifactCodecAbort(
            "DIA S4 codec requires one supported S4 object",
            "multischolar_missing_exact_s4_codec"
        )
    }
    class_name <- class(value)[[1L]]
    descriptor <- artifactDiaCodecDescriptor(class_name)
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

validateDiaS4Bundle <- function(bundle) {
    if (!is.list(bundle) || !is.list(bundle$metadata) || !is.list(bundle$payloads)) {
        artifactCodecAbort(
            "DIA S4 artifact bundle is malformed",
            "multischolar_invalid_s4_artifact_bundle"
        )
    }
    metadata <- bundle$metadata
    if (!identical(metadata$schema, .artifactS4BundleSchema) ||
        !identical(metadata$schema_version, .artifactS4BundleSchemaVersion)) {
        artifactCodecAbort(
            "DIA S4 artifact bundle schema version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    descriptor <- artifactDiaCodecDescriptor(metadata$class_name)
    if (!identical(metadata$codec, descriptor) ||
        !identical(metadata$codec$version, .artifactDiaCodecVersion)) {
        artifactCodecAbort(
            "DIA S4 artifact codec version is unsupported",
            "multischolar_unsupported_artifact_codec_version"
        )
    }
    expected_slots <- methods::slotNames(metadata$class_name)
    if (!identical(metadata$slot_names, expected_slots) ||
        !identical(names(metadata$slot_values), expected_slots) ||
        !identical(names(metadata$payloads), names(bundle$payloads))) {
        artifactCodecAbort(
            "DIA S4 artifact bundle shape does not match the class contract",
            "multischolar_artifact_shape_mismatch"
        )
    }
    expected_digest <- artifactSemanticDigest(
        artifactDiaBundleSemanticInput(metadata)
    )
    if (!identical(metadata$semantic_digest, expected_digest)) {
        artifactCodecAbort(
            "DIA S4 artifact semantic lineage digest does not match its metadata",
            "multischolar_artifact_semantic_digest_mismatch"
        )
    }
    bundle
}

hydrateDiaS4Artifact <- function(bundle) {
    bundle <- validateDiaS4Bundle(bundle)
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
                "DIA S4 hydration changed slot(s): %s",
                paste(metadata$slot_names[!exact_slots], collapse = ", ")
            ),
            "multischolar_inexact_s4_hydration",
            owner = metadata$class_name
        )
    }
    validity <- methods::validObject(object, test = TRUE)
    if (!identical(validity, TRUE)) {
        artifactCodecAbort(
            sprintf("hydrated DIA S4 object is invalid: %s", paste(validity, collapse = "; ")),
            "multischolar_invalid_hydrated_s4_object",
            owner = metadata$class_name
        )
    }
    object
}

encodeDiaS4Artifact <- dehydrateDiaS4Artifact
decodeDiaS4Artifact <- hydrateDiaS4Artifact
