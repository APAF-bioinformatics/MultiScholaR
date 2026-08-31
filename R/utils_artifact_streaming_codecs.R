# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactStreamingCodec <- "multischolar.rectangular_streaming"
.artifactStreamingCodecVersion <- 1L
.artifactStreamingSemanticVersion <- 1L
.artifactStreamingBlockRows <- 4096L

artifactStreamingRType <- function(field) {
    type <- field$type$ToString()
    if (identical(type, "bool")) return("logical")
    if (type %in% c("int8", "int16", "int32")) return("integer")
    if (type %in% c("int64", "uint8", "uint16", "uint32", "uint64")) {
        return("double")
    }
    if (type %in% c("float", "double")) return("double")
    if (type %in% c("string", "large_string")) return("character")
    if (startsWith(type, "date32") || startsWith(type, "date64")) return("Date")
    if (startsWith(type, "timestamp")) return("POSIXct")
    artifactCodecAbort(
        sprintf("streaming artifact field '%s' has unsupported type '%s'", field$name, type),
        "multischolar_unsupported_artifact_column",
        owner = field$name
    )
}

artifactStreamingDescriptor <- function(field) {
    list(
        physical_name = field$name,
        status_name = NULL,
        r_type = artifactStreamingRType(field),
        nullable = isTRUE(field$nullable),
        levels = NULL,
        timezone = NULL,
        units = NULL
    )
}

artifactStreamingCanonicalAtomic <- function(value, r_type) {
    switch(r_type,
        logical = ifelse(is.na(value), "NA", ifelse(value, "TRUE", "FALSE")),
        integer = ifelse(is.na(value), "NA", as.character(value)),
        double = artifactInlineDouble(as.numeric(value)),
        character = artifactInlineCharacter(as.character(value)),
        Date = artifactInlineDouble(as.numeric(value)),
        POSIXct = artifactInlineDouble(as.numeric(value)),
        artifactCodecAbort(
            "streaming artifact canonicalization has an unsupported type",
            "multischolar_unsupported_artifact_column"
        )
    )
}

artifactStreamingBlockDigest <- function(payload, metadata, offset, length) {
    columns <- lapply(seq_along(metadata$columns), function(index) {
        descriptor <- metadata$columns[[index]]
        field_index <- match(descriptor$physical_name, payload$ColumnNames())
        value <- payload$column(field_index - 1L)$Slice(offset, length)$as_vector()
        artifactStreamingCanonicalAtomic(value, descriptor$r_type)
    })
    artifactSemanticDigest(list(
        schema = "multischolar.streaming_semantic_block",
        schema_version = .artifactStreamingSemanticVersion,
        row_start = as.numeric(offset + 1L),
        row_count = as.integer(length),
        columns = unname(columns)
    ))
}

artifactStreamingBlockDigests <- function(payload, metadata) {
    rows <- as.integer(payload$num_rows)
    if (rows == 0L) return(list())
    starts <- seq.int(0L, rows - 1L, by = .artifactStreamingBlockRows)
    unname(lapply(starts, function(offset) {
        artifactStreamingBlockDigest(
            payload,
            metadata,
            offset,
            min(.artifactStreamingBlockRows, rows - offset)
        )
    }))
}

artifactStreamingSemanticInput <- function(metadata, block_digests) {
    list(
        schema = "multischolar.streaming_rectangular_semantic",
        schema_version = .artifactStreamingSemanticVersion,
        codec = metadata$codec,
        kind = metadata$kind,
        dimensions = metadata$dimensions,
        data_frame_class = metadata$data_frame_class,
        logical_names = metadata$logical_names,
        row_names = metadata$row_names,
        columns = metadata$columns,
        stable_key = metadata$stable_key,
        block_rows = .artifactStreamingBlockRows,
        block_digests = block_digests
    )
}

artifactStreamingWriterSettings <- function(rows) {
    compression <- if (arrow::codec_is_available("zstd")) "zstd" else "snappy"
    list(
        parquet_version = "2.6",
        compression = compression,
        compression_level = if (identical(compression, "zstd")) 3L else NULL,
        chunk_size = as.integer(max(1L, min(65536L, max(1L, rows)))),
        use_dictionary = TRUE,
        dictionary_columns = character(),
        write_statistics = TRUE,
        data_page_size = 1024 * 1024,
        partitioning = "none",
        target_file_size_bytes = 128 * 1024 * 1024,
        minimum_file_size_bytes = 32 * 1024 * 1024,
        maximum_file_size_bytes = 512 * 1024 * 1024,
        append_policy = "immutable_new_generation",
        small_file_policy = "one_payload_file_per_generation"
    )
}

encodeArtifactStreamingParquet <- function(path, owner = "streaming table") {
    payload <- arrow::read_parquet(path, as_data_frame = FALSE)
    fields <- payload$schema$fields
    field_names <- vapply(fields, `[[`, character(1), "name")
    if (!identical(field_names[[1L]], .artifactRowOrderColumn)) {
        artifactCodecAbort(
            "streaming Parquet lacks the canonical row-order field",
            "multischolar_artifact_order_mismatch",
            owner = owner
        )
    }
    logical_fields <- fields[-1L]
    descriptors <- lapply(logical_fields, artifactStreamingDescriptor)
    logical_names <- vapply(logical_fields, `[[`, character(1), "name")
    rows <- as.integer(payload$num_rows)
    metadata <- list(
        codec = list(
            id = .artifactStreamingCodec,
            version = .artifactStreamingCodecVersion
        ),
        payload_schema = list(
            id = "multischolar.parquet_table",
            version = .artifactParquetSchemaVersion
        ),
        kind = "data.frame",
        owner = owner,
        dimensions = list(rows = rows, columns = length(descriptors)),
        data_frame_class = c("tbl_df", "tbl", "data.frame"),
        logical_names = logical_names,
        row_names = list(kind = "automatic", value = NULL),
        columns = descriptors,
        physical_schema = artifactPhysicalSchemaMetadata(fields),
        stable_key = list(
            kind = "artifact_row_order",
            logical_columns = character(),
            physical_columns = .artifactRowOrderColumn
        ),
        schema_evolution = list(
            policy = "reject_unknown",
            allowed_changes = character()
        ),
        writer_settings = artifactStreamingWriterSettings(rows),
        internal_columns = .artifactRowOrderColumn,
        streaming_semantic = NULL,
        semantic_digest = NULL
    )
    blocks <- artifactStreamingBlockDigests(payload, metadata)
    metadata$streaming_semantic <- list(
        version = .artifactStreamingSemanticVersion,
        block_rows = .artifactStreamingBlockRows,
        block_digests = blocks
    )
    metadata$semantic_digest <- artifactSemanticDigest(
        artifactStreamingSemanticInput(metadata, blocks)
    )
    structure(
        list(payload = payload, metadata = metadata, source_path = path),
        class = c(
            "MultiScholaRArtifactStreamingRectangular",
            "MultiScholaRArtifactRectangular",
            "list"
        )
    )
}

artifactStorePublishStreamingParquet <- function(
    store,
    encoded,
    logical_key,
    provenance_ids = character(),
    failure_injector = NULL,
    verification = "deferred_exact"
) {
    if (!inherits(encoded, "MultiScholaRArtifactStreamingRectangular") ||
        !file.exists(encoded$source_path)) {
        artifactStoreAbort(
            "streaming publication requires one encoded Parquet file",
            "multischolar_invalid_artifact_store_payload"
        )
    }
    source_path <- encoded$source_path
    copy_parquet <- function(x, sink, ...) {
        copied <- file.copy(source_path, sink, overwrite = FALSE)
        if (!isTRUE(copied)) stop("streaming Parquet copy failed", call. = FALSE)
        invisible(sink)
    }
    artifactStoreWriteParquet(
        store,
        encoded,
        logical_key,
        provenance_ids = provenance_ids,
        failure_injector = failure_injector,
        write_parquet_fn = copy_parquet,
        verification = verification
    )
}

validateArtifactStreamingMetadata <- function(metadata) {
    if (!is.list(metadata)) {
        artifactCodecAbort(
            "streaming rectangular metadata is malformed",
            "multischolar_invalid_artifact_metadata"
        )
    }
    block_digests <- as.list(unlist(
        metadata$streaming_semantic$block_digests,
        use.names = FALSE
    ))
    metadata$streaming_semantic$block_digests <- block_digests
    expected <- c(
        "codec", "payload_schema", "kind", "owner", "dimensions",
        "data_frame_class", "logical_names", "row_names", "columns",
        "physical_schema", "stable_key", "schema_evolution",
        "writer_settings", "internal_columns", "streaming_semantic",
        "semantic_digest"
    )
    dimensions <- metadata$dimensions
    descriptor_names <- c(
        "physical_name", "status_name", "r_type", "nullable", "levels",
        "timezone", "units"
    )
    descriptors_valid <- is.list(metadata$columns) && all(vapply(
        metadata$columns,
        function(descriptor) {
            is.list(descriptor) && identical(names(descriptor), descriptor_names) &&
                workflowCapabilityScalarString(descriptor$physical_name) &&
                descriptor$r_type %in% c(
                    "logical", "integer", "double", "character", "Date", "POSIXct"
                ) && is.logical(descriptor$nullable) &&
                length(descriptor$nullable) == 1L
        },
        logical(1)
    ))
    expected_physical_names <- c(
        .artifactRowOrderColumn,
        vapply(metadata$columns, `[[`, character(1), "physical_name")
    )
    physical_names <- if (is.list(metadata$physical_schema)) {
        vapply(metadata$physical_schema, `[[`, character(1), "name")
    } else {
        character()
    }
    valid <- is.list(metadata) && identical(names(metadata), expected) &&
        identical(metadata$codec, list(
            id = .artifactStreamingCodec,
            version = .artifactStreamingCodecVersion
        )) && identical(metadata$payload_schema, list(
            id = "multischolar.parquet_table",
            version = .artifactParquetSchemaVersion
        )) && identical(metadata$kind, "data.frame") &&
        is.list(dimensions) &&
        identical(names(dimensions), c("rows", "columns")) &&
        all(vapply(dimensions, function(value) {
            is.numeric(value) && length(value) == 1L && is.finite(value) &&
                value >= 0 && value == floor(value)
        }, logical(1))) && descriptors_valid &&
        length(metadata$columns) == dimensions$columns &&
        identical(metadata$logical_names, expected_physical_names[-1L]) &&
        identical(physical_names, expected_physical_names) &&
        identical(metadata$streaming_semantic$version, 1L) &&
        identical(
            metadata$streaming_semantic$block_rows,
            .artifactStreamingBlockRows
        ) && all(vapply(
            block_digests,
            function(value) workflowPolicyDigestValid(value),
            logical(1)
        ))
    if (!isTRUE(valid)) {
        artifactCodecAbort(
            "streaming rectangular metadata is malformed",
            "multischolar_invalid_artifact_metadata"
        )
    }
    artifactRefValidateDigest(metadata$semantic_digest, "semantic_digest")
    metadata
}

verifyArtifactStreamingPayload <- function(payload, metadata) {
    metadata <- validateArtifactStreamingMetadata(metadata)
    blocks <- artifactStreamingBlockDigests(payload, metadata)
    observed <- artifactSemanticDigest(
        artifactStreamingSemanticInput(metadata, blocks)
    )
    valid <- identical(blocks, metadata$streaming_semantic$block_digests) &&
        identical(observed, metadata$semantic_digest)
    if (!isTRUE(valid)) {
        artifactCodecAbort(
            "streaming artifact semantic blocks differ",
            "multischolar_artifact_semantic_digest_mismatch"
        )
    }
    invisible(TRUE)
}
