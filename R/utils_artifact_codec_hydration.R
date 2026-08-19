# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactMetadataAbort <- function(message, owner = NULL) {
    artifactCodecAbort(
        message,
        "multischolar_invalid_artifact_metadata",
        owner = owner
    )
}

artifactMetadataCharacter <- function(value) {
    if (is.list(value) && length(value) == 0L) return(character())
    value
}

artifactExpectedPhysicalSchema <- function(metadata) {
    schema <- list(list(
        name = .artifactRowOrderColumn,
        type = "int64",
        nullable = FALSE
    ))
    for (descriptor in metadata$columns) {
        schema[[length(schema) + 1L]] <- list(
            name = descriptor$physical_name,
            type = artifactArrowType(descriptor$r_type)$ToString(),
            nullable = descriptor$nullable
        )
        if (!is.null(descriptor$status_name)) {
            schema[[length(schema) + 1L]] <- list(
                name = descriptor$status_name,
                type = "int32",
                nullable = FALSE
            )
        }
    }
    schema
}

artifactValidateRectangularMetadata <- function(metadata) {
    if (!is.list(metadata) || !workflowCapabilityScalarString(metadata$kind) ||
        !metadata$kind %in% c("data.frame", "matrix")) {
        artifactMetadataAbort("artifact rectangular metadata is malformed")
    }
    expected_names <- if (identical(metadata$kind, "data.frame")) {
        c(
            "codec", "payload_schema", "kind", "owner", "dimensions",
            "data_frame_class", "logical_names", "row_names", "columns",
            "physical_schema", "stable_key", "schema_evolution",
            "writer_settings", "internal_columns", "semantic_digest"
        )
    } else {
        c(
            "codec", "payload_schema", "kind", "owner", "dimensions",
            "columns", "physical_schema", "stable_key", "schema_evolution",
            "writer_settings", "internal_columns", "matrix", "semantic_digest"
        )
    }
    if (!identical(names(metadata), expected_names)) {
        artifactMetadataAbort(
            "artifact rectangular metadata has unknown or missing fields",
            metadata$owner
        )
    }
    valid_version <- identical(
        metadata$codec,
        list(id = .artifactRectangularCodec, version = .artifactRectangularCodecVersion)
    ) && identical(
        metadata$payload_schema,
        list(id = "multischolar.parquet_table", version = .artifactParquetSchemaVersion)
    )
    if (!isTRUE(valid_version)) {
        artifactCodecAbort(
            "artifact rectangular metadata version is unsupported",
            "multischolar_unsupported_artifact_codec_version",
            owner = metadata$owner
        )
    }
    dimensions <- metadata$dimensions
    valid_dimensions <- is.list(dimensions) &&
        identical(names(dimensions), c("rows", "columns")) &&
        all(vapply(dimensions, function(value) {
            length(value) == 1L && is.numeric(value) && !is.na(value) &&
                is.finite(value) && value >= 0 && value == floor(value)
        }, logical(1)))
    if (!isTRUE(valid_dimensions) ||
        !is.list(metadata$columns) ||
        length(metadata$columns) != dimensions$columns) {
        artifactMetadataAbort("artifact rectangular dimensions are invalid", metadata$owner)
    }
    descriptor_names <- c(
        "physical_name", "status_name", "r_type", "nullable", "levels",
        "timezone", "units"
    )
    descriptors_valid <- vapply(metadata$columns, function(descriptor) {
        is.list(descriptor) && identical(names(descriptor), descriptor_names) &&
            workflowCapabilityScalarString(descriptor$physical_name) &&
            workflowCapabilityScalarString(descriptor$r_type) &&
            length(descriptor$nullable) == 1L && is.logical(descriptor$nullable)
    }, logical(1))
    if (!all(descriptors_valid)) {
        artifactMetadataAbort("artifact column descriptors are malformed", metadata$owner)
    }
    expected_schema <- artifactExpectedPhysicalSchema(metadata)
    if (!identical(metadata$physical_schema, expected_schema)) {
        artifactMetadataAbort(
            "artifact physical schema does not match its logical descriptors",
            metadata$owner
        )
    }
    artifactRefValidateDigest(metadata$semantic_digest, "semantic_digest")
    metadata$logical_names <- artifactMetadataCharacter(metadata$logical_names)
    metadata
}

artifactPayloadDataFrame <- function(payload, metadata) {
    if (inherits(payload, "ArrowTabular")) {
        fields <- payload$schema$fields
        actual_schema <- lapply(fields, function(field) {
            list(
                name = field$name,
                type = field$type$ToString(),
                nullable = field$nullable
            )
        })
        if (!identical(actual_schema, metadata$physical_schema)) {
            artifactCodecAbort(
                "Arrow payload schema does not match the declared physical schema",
                "multischolar_artifact_physical_schema_mismatch",
                owner = metadata$owner
            )
        }
        return(as.data.frame(payload))
    }
    if (is.data.frame(payload)) return(payload)
    artifactCodecAbort(
        "artifact payload must be an Arrow table or data.frame",
        "multischolar_invalid_artifact_payload"
    )
}

artifactValidatePhysicalPayload <- function(payload, metadata) {
    expected_names <- vapply(metadata$physical_schema, `[[`, character(1), "name")
    if (!identical(names(payload), expected_names) ||
        nrow(payload) != metadata$dimensions$rows) {
        artifactCodecAbort(
            sprintf("artifact payload shape does not match '%s' metadata", metadata$owner),
            "multischolar_artifact_shape_mismatch",
            owner = metadata$owner
        )
    }
    for (field in metadata$physical_schema) {
        if (!isTRUE(field$nullable) && anyNA(payload[[field$name]])) {
            artifactCodecAbort(
                sprintf("artifact payload field '%s' violates nullability", field$name),
                "multischolar_artifact_nullability_mismatch",
                owner = metadata$owner
            )
        }
    }
    row_order <- as.numeric(payload[[.artifactRowOrderColumn]])
    if (!identical(row_order, as.numeric(seq_len(nrow(payload))))) {
        artifactCodecAbort(
            "artifact payload row order is missing, duplicated, or reordered",
            "multischolar_artifact_order_mismatch",
            owner = metadata$owner
        )
    }
    invisible(payload)
}

artifactDecodeColumn <- function(payload, descriptor) {
    value <- payload[[descriptor$physical_name]]
    if (!is.null(descriptor$status_name)) {
        status <- as.integer(payload[[descriptor$status_name]])
        if (anyNA(status) || any(!status %in% 0:4)) {
            artifactCodecAbort(
                "artifact non-finite status column is invalid",
                "multischolar_invalid_artifact_payload"
            )
        }
        value <- as.numeric(value)
        value[status == 1L] <- NA_real_
        value[status == 2L] <- NaN
        value[status == 3L] <- Inf
        value[status == 4L] <- -Inf
    }
    value <- switch(descriptor$r_type,
        logical = as.logical(value),
        integer = as.integer(value),
        integer64 = {
            if (!requireNamespace("bit64", quietly = TRUE)) {
                artifactCodecAbort(
                    "bit64 is required to hydrate an integer64 artifact column",
                    "multischolar_missing_artifact_codec_dependency"
                )
            }
            bit64::as.integer64(value)
        },
        double = as.numeric(value),
        character = as.character(value),
        factor = factor(value, levels = descriptor$levels),
        ordered = ordered(value, levels = descriptor$levels),
        Date = structure(as.numeric(value), class = "Date"),
        POSIXct = structure(
            as.numeric(value),
            class = c("POSIXct", "POSIXt"),
            tzone = descriptor$timezone
        )
    )
    if (!is.null(descriptor$units)) attr(value, "units") <- descriptor$units
    value
}

decodeArtifactRectangular <- function(payload, metadata) {
    metadata <- artifactValidateRectangularMetadata(metadata)
    payload <- artifactPayloadDataFrame(payload, metadata)
    artifactValidatePhysicalPayload(payload, metadata)
    columns <- lapply(metadata$columns, function(descriptor) {
        artifactDecodeColumn(payload, descriptor)
    })
    value <- structure(columns, names = metadata$logical_names)
    value <- as.data.frame(value, optional = TRUE, stringsAsFactors = FALSE)
    if (identical(metadata$kind, "matrix")) {
        matrix_value <- as.matrix(value)
        storage.mode(matrix_value) <- metadata$matrix$storage_mode
        dim(matrix_value) <- as.integer(metadata$matrix$dimensions)
        dimnames(matrix_value) <- metadata$matrix$dimnames
        return(matrix_value)
    }
    if (identical(metadata$row_names$kind, "automatic")) {
        attr(value, "row.names") <- c(NA_integer_, -nrow(value))
    } else {
        attr(value, "row.names") <- metadata$row_names$value
    }
    class(value) <- metadata$data_frame_class
    value
}
