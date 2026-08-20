# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
.artifactRectangularCodec <- "multischolar.rectangular"
.artifactRectangularCodecVersion <- 1L
.artifactParquetSchemaVersion <- 1L
.artifactRowOrderColumn <- ".multischolar_row_order"
.artifactValuePrefix <- ".multischolar_value_"
.artifactStatusPrefix <- ".multischolar_status_"
.artifactInlineLimitBytes <- 1024L * 1024L
artifactCodecAbort <- function(message, class, ..., owner = NULL) {
    details <- if (is.null(owner)) list(...) else c(list(owner = owner), list(...))
    do.call(
        rlang::abort,
        c(
            list(
                message = message,
                class = c(class, "multischolar_artifact_codec_error")
            ),
            details
        )
    )
}
artifactPhysicalName <- function(prefix, index) {
    sprintf("%s%08d", prefix, as.integer(index))
}
artifactSupportedColumnType <- function(value, owner) {
    if (inherits(value, "integer64")) {
        return("integer64")
    }
    if (is.factor(value)) {
        return(if (is.ordered(value)) "ordered" else "factor")
    }
    if (inherits(value, "Date")) {
        return("Date")
    }
    if (inherits(value, "POSIXct")) {
        return("POSIXct")
    }
    explicit_class <- attr(value, "class", exact = TRUE)
    if (is.logical(value) && is.null(explicit_class)) {
        return("logical")
    }
    if (is.integer(value) && is.null(explicit_class)) {
        return("integer")
    }
    if (is.double(value) && is.null(explicit_class)) {
        return("double")
    }
    if (is.character(value) && is.null(explicit_class)) {
        return("character")
    }
    artifactCodecAbort(
        sprintf("artifact column '%s' has unsupported R type", owner),
        "multischolar_unsupported_artifact_column",
        owner = owner
    )
}
artifactColumnUnits <- function(value, owner) {
    units <- attr(value, "units", exact = TRUE)
    if (is.null(units)) {
        return(NULL)
    }
    if (!is.character(units) || anyNA(units)) {
        artifactCodecAbort(
            sprintf("artifact column '%s' has unsupported units metadata", owner),
            "multischolar_unsupported_artifact_attribute",
            owner = owner
        )
    }
    unname(units)
}
artifactNonFiniteStatus <- function(value) {
    status <- integer(length(value))
    status[is.na(value)] <- 1L
    status[is.nan(value)] <- 2L
    status[is.infinite(value) & value > 0] <- 3L
    status[is.infinite(value) & value < 0] <- 4L
    status
}
artifactArrowType <- function(type) {
    switch(type,
        logical = arrow::boolean(),
        integer = arrow::int32(),
        integer64 = arrow::int64(),
        double = arrow::float64(),
        Date = arrow::float64(),
        POSIXct = arrow::float64(),
        factor = arrow::utf8(),
        ordered = arrow::utf8(),
        character = arrow::utf8(),
        artifactCodecAbort(
            sprintf("unsupported Arrow mapping for R type '%s'", type),
            "multischolar_unsupported_artifact_column"
        )
    )
}
artifactEncodeColumn <- function(value, index, owner) {
    type <- artifactSupportedColumnType(value, owner)
    physical_name <- artifactPhysicalName(.artifactValuePrefix, index)
    status_name <- NULL
    physical <- switch(type,
        factor = as.character(value),
        ordered = as.character(value),
        Date = unclass(value),
        POSIXct = unclass(value),
        value
    )
    status <- NULL
    if (type %in% c("double", "Date", "POSIXct")) {
        status <- artifactNonFiniteStatus(physical)
        if (any(status != 0L)) {
            status_name <- artifactPhysicalName(.artifactStatusPrefix, index)
            physical[status != 0L] <- NA_real_
        }
    }
    descriptor <- list(
        physical_name = physical_name,
        status_name = status_name,
        r_type = type,
        nullable = anyNA(physical),
        levels = if (type %in% c("factor", "ordered")) levels(value) else NULL,
        timezone = if (identical(type, "POSIXct")) {
            attr(value, "tzone", exact = TRUE)
        } else {
            NULL
        },
        units = artifactColumnUnits(value, owner)
    )
    payload <- list(physical)
    names(payload) <- physical_name
    fields <- list(arrow::field(
        physical_name,
        artifactArrowType(type),
        nullable = descriptor$nullable
    ))
    if (!is.null(status_name)) {
        payload[[status_name]] <- status
        fields[[length(fields) + 1L]] <- arrow::field(
            status_name,
            arrow::int32(),
            nullable = FALSE
        )
    }
    list(payload = payload, fields = fields, descriptor = descriptor)
}

artifactWriterSettings <- function(payload, descriptors, target_rows = 65536L) {
    row_count <- if (length(payload) == 0L) 0L else length(payload[[1L]])
    dictionary_columns <- character()
    character_columns <- character()
    for (descriptor in descriptors) {
        if (!descriptor$r_type %in% c("character", "factor", "ordered")) next
        character_columns <- c(character_columns, descriptor$physical_name)
        values <- payload[[descriptor$physical_name]]
        observed <- values[!is.na(values)]
        cardinality <- length(unique(observed))
        if (length(observed) > 0L && cardinality <= 65536L &&
            cardinality / length(observed) <= 0.5) {
            dictionary_columns <- c(dictionary_columns, descriptor$physical_name)
        }
    }
    compression <- if (arrow::codec_is_available("zstd")) "zstd" else "snappy"
    list(
        parquet_version = "2.6",
        compression = compression,
        compression_level = if (identical(compression, "zstd")) 3L else NULL,
        chunk_size = as.integer(max(1L, min(target_rows, max(1L, row_count)))),
        use_dictionary = length(character_columns) > 0L &&
            identical(dictionary_columns, character_columns),
        dictionary_columns = dictionary_columns,
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

artifactPhysicalSchemaMetadata <- function(fields) {
    lapply(fields, function(field) {
        list(
            name = field$name,
            type = field$type$ToString(),
            nullable = field$nullable
        )
    })
}

artifactRowNamesMetadata <- function(value) {
    row_names <- attr(value, "row.names", exact = TRUE)
    automatic <- is.integer(row_names) && (
        (length(row_names) == 0L && nrow(value) == 0L) ||
            (length(row_names) == 2L && is.na(row_names[[1L]]) &&
                row_names[[2L]] <= 0L)
    )
    list(
        kind = if (automatic) "automatic" else "explicit",
        value = if (automatic) NULL else unname(row_names)
    )
}

artifactRectangularSemanticDigest <- function(metadata, payload) {
    column_digests <- lapply(payload, function(column) {
        values <- if (inherits(column, "integer64")) as.character(column) else unclass(column)
        digest::digest(
            values,
            algo = .artifactHashAlgorithm,
            serialize = TRUE,
            ascii = FALSE,
            serializeVersion = 3L
        )
    })
    artifactSemanticDigest(list(
        metadata = metadata,
        physical_column_digests = unname(column_digests)
    ))
}

encodeArtifactTable <- function(value, stable_key = NULL, owner = "table") {
    if (!is.data.frame(value)) {
        artifactCodecAbort(
            "artifact table codec requires a data.frame",
            "multischolar_invalid_artifact_table",
            owner = owner
        )
    }
    payload <- list(.multischolar_row_order = seq_len(nrow(value)))
    fields <- list(arrow::field(
        .artifactRowOrderColumn,
        arrow::int64(),
        nullable = FALSE
    ))
    descriptors <- vector("list", ncol(value))
    for (index in seq_len(ncol(value))) {
        encoded <- artifactEncodeColumn(
            value[[index]],
            index,
            paste0(owner, "$", names(value)[[index]])
        )
        payload <- c(payload, encoded$payload)
        fields <- c(fields, encoded$fields)
        descriptors[[index]] <- encoded$descriptor
    }
    schema <- do.call(arrow::schema, fields)
    arrow_table <- do.call(
        arrow::Table$create,
        c(payload, list(schema = schema))
    )
    if (is.null(stable_key)) {
        stable_key_metadata <- list(
            kind = "artifact_row_order",
            logical_columns = character(),
            physical_columns = .artifactRowOrderColumn
        )
    } else {
        if (!is.character(stable_key) || length(stable_key) == 0L ||
            anyNA(stable_key) || anyDuplicated(stable_key) > 0L) {
            artifactCodecAbort(
                sprintf("artifact table '%s' has an invalid stable key", owner),
                "multischolar_invalid_artifact_stable_key",
                owner = owner
            )
        }
        key_occurrences <- vapply(stable_key, function(key) {
            sum(names(value) == key)
        }, integer(1))
        if (any(key_occurrences != 1L)) {
            artifactCodecAbort(
                sprintf("artifact table '%s' has an ambiguous stable key", owner),
                "multischolar_invalid_artifact_stable_key",
                owner = owner
            )
        }
        key_indexes <- match(stable_key, names(value))
        key_frame <- value[key_indexes]
        invalid_key <- any(!stats::complete.cases(key_frame)) ||
            anyDuplicated(key_frame) > 0L
        if (any(invalid_key)) {
            artifactCodecAbort(
                sprintf("artifact table '%s' stable key is missing or duplicated", owner),
                "multischolar_invalid_artifact_stable_key",
                owner = owner
            )
        }
        stable_key_metadata <- list(
            kind = "logical_columns",
            logical_columns = stable_key,
            physical_columns = vapply(
                descriptors[key_indexes],
                `[[`,
                character(1),
                "physical_name"
            )
        )
    }
    writer_settings <- artifactWriterSettings(payload, descriptors)
    metadata <- list(
        codec = list(
            id = .artifactRectangularCodec,
            version = .artifactRectangularCodecVersion
        ),
        payload_schema = list(
            id = "multischolar.parquet_table",
            version = .artifactParquetSchemaVersion
        ),
        kind = "data.frame",
        owner = owner,
        dimensions = list(rows = nrow(value), columns = ncol(value)),
        data_frame_class = unname(class(value)),
        logical_names = names(value),
        row_names = artifactRowNamesMetadata(value),
        columns = descriptors,
        physical_schema = artifactPhysicalSchemaMetadata(fields),
        stable_key = stable_key_metadata,
        schema_evolution = list(
            policy = "reject_unknown",
            allowed_changes = character()
        ),
        writer_settings = writer_settings,
        internal_columns = setdiff(names(payload), vapply(
            descriptors,
            `[[`,
            character(1),
            "physical_name"
        ))
    )
    metadata$semantic_digest <- artifactRectangularSemanticDigest(metadata, payload)
    structure(
        list(payload = arrow_table, metadata = metadata),
        class = c("MultiScholaRArtifactRectangular", "list")
    )
}

encodeArtifactMatrix <- function(value, owner = "matrix") {
    if (!is.matrix(value) || !typeof(value) %in% c("logical", "integer", "double", "character")) {
        artifactCodecAbort(
            "artifact matrix codec requires a supported atomic matrix",
            "multischolar_invalid_artifact_matrix",
            owner = owner
        )
    }
    table_value <- as.data.frame(value, optional = TRUE, stringsAsFactors = FALSE)
    names(table_value) <- rep("", ncol(table_value))
    encoded <- encodeArtifactTable(table_value, owner = owner)
    encoded$metadata$kind <- "matrix"
    encoded$metadata$data_frame_class <- NULL
    encoded$metadata$logical_names <- NULL
    encoded$metadata$row_names <- NULL
    encoded$metadata$matrix <- list(
        storage_mode = storage.mode(value),
        dimensions = unname(dim(value)),
        dimnames = lapply(dimnames(value), function(item) {
            if (is.null(item)) NULL else unname(item)
        })
    )
    encoded$metadata$semantic_digest <- artifactRectangularSemanticDigest(
        encoded$metadata[names(encoded$metadata) != "semantic_digest"],
        lapply(seq_len(encoded$payload$num_columns), function(index) {
            encoded$payload$column(index - 1L)$as_vector()
        })
    )
    encoded$metadata <- encoded$metadata[c(
        "codec", "payload_schema", "kind", "owner", "dimensions", "columns",
        "physical_schema", "stable_key", "schema_evolution", "writer_settings",
        "internal_columns", "matrix", "semantic_digest"
    )]
    encoded
}

artifactParquetWriteArgs <- function(encoded) {
    if (!inherits(encoded, "MultiScholaRArtifactRectangular")) {
        artifactCodecAbort(
            "Parquet writer settings require an encoded rectangular artifact",
            "multischolar_invalid_artifact_table"
        )
    }
    settings <- encoded$metadata$writer_settings
    list(
        version = settings$parquet_version,
        compression = settings$compression,
        compression_level = settings$compression_level,
        chunk_size = settings$chunk_size,
        use_dictionary = settings$use_dictionary,
        write_statistics = settings$write_statistics,
        data_page_size = settings$data_page_size
    )
}
