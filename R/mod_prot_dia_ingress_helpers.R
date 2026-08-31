# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaValidateImportSummary <- function(summary) {
    required <- c(
        "rows", "columns", "column_names", "run_values", "protein_count",
        "peptide_count"
    )
    scalar_count <- function(value) {
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value >= 0 && value == floor(value)
    }
    valid <- is.list(summary) && identical(names(summary), required) &&
        all(vapply(
            summary[c("rows", "columns", "protein_count", "peptide_count")],
            scalar_count,
            logical(1)
        )) && is.list(summary$column_names) &&
        all(vapply(summary$column_names, workflowCapabilityScalarString, logical(1))) &&
        is.list(summary$run_values) &&
        all(vapply(summary$run_values, workflowCapabilityScalarString, logical(1))) &&
        length(summary$column_names) == summary$columns &&
        !anyDuplicated(unlist(summary$column_names, use.names = FALSE))
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN deferred ingress summary is invalid",
            "multischolar_invalid_prot_dia_ingress_summary"
        )
    }
    summary
}

protDiaIngressTsvColumns <- function(source_path, max_header_bytes) {
    connection <- file(source_path, open = "rb")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header) != 1L || nchar(header, type = "bytes") > max_header_bytes) {
        protDiaArtifactAbort(
            "DIA-NN ingress header is missing or exceeds its bound",
            "multischolar_invalid_prot_dia_ingress_header"
        )
    }
    strsplit(header, "\t", fixed = TRUE)[[1L]]
}

protDiaIngressParquetColumns <- function(source_path) {
    dataset <- arrow::open_dataset(source_path, format = "parquet")
    dataset$schema$names
}

protDiaIngressPreflight <- function(
    source_path,
    use_precursor_norm = TRUE,
    max_header_bytes = 1024^2
) {
    if (!artifactResourceScalarString(source_path) || !file.exists(source_path) ||
        dir.exists(source_path)) {
        protDiaArtifactAbort(
            "DIA-NN ingress source must be one regular file",
            "multischolar_invalid_prot_dia_ingress_source"
        )
    }
    extension <- tolower(tools::file_ext(source_path))
    columns <- if (identical(extension, "parquet")) {
        protDiaIngressParquetColumns(source_path)
    } else {
        protDiaIngressTsvColumns(source_path, max_header_bytes)
    }
    required <- c("Protein.Group", "Stripped.Sequence", "Run")
    missing <- setdiff(required, columns)
    quantity <- if (isTRUE(use_precursor_norm) &&
        "Precursor.Normalised" %in% columns) {
        "Precursor.Normalised"
    } else if ("Precursor.Quantity" %in% columns) {
        "Precursor.Quantity"
    } else {
        NULL
    }
    if (length(missing) > 0L || is.null(quantity) || anyDuplicated(columns)) {
        protDiaArtifactAbort(
            "DIA-NN ingress schema is missing, duplicated, or unsupported",
            "multischolar_invalid_prot_dia_ingress_schema",
            missing_columns = missing
        )
    }
    info <- file.info(source_path)
    list(
        schema = "multischolar.prot_dia_ingress_preflight",
        schema_version = 1L,
        source_size_bytes = unname(as.numeric(info$size)),
        source_digest = artifactByteDigest(source_path),
        column_names = as.list(columns),
        quantity_column = quantity,
        complete_payload_materialized = FALSE
    )
}

validateProtDiaIngressBinding <- function(preflight, writer) {
    valid <- is.list(preflight) && is.list(writer$source) &&
        identical(
            preflight$source_size_bytes,
            writer$source$source_size_bytes
        ) && identical(
            preflight$source_digest,
            writer$source$source_digest
        ) && identical(
            unlist(preflight$column_names, use.names = FALSE),
            unlist(writer$import_summary$column_names, use.names = FALSE)
        ) && identical(
            preflight$quantity_column,
            writer$column_mapping$quantity_col
        ) && !isTRUE(preflight$complete_payload_materialized)
    if (!isTRUE(valid)) {
        protDiaArtifactAbort(
            "DIA-NN source changed between preflight and artifact staging",
            "multischolar_mutated_prot_dia_ingress_source"
        )
    }
    invisible(TRUE)
}

protDiaDuckdbCsvExpression <- function(connection, source_path) {
    source <- as.character(DBI::dbQuoteString(connection, source_path))
    paste0(
        "read_csv(", source, ", delim = chr(9), header = true, ",
        "auto_detect = true, sample_size = -1, ",
        "nullstr = ['', 'NA'], strict_mode = true)"
    )
}

protDiaStreamingSelect <- function(connection, source_path, use_precursor_norm) {
    source <- protDiaDuckdbCsvExpression(connection, source_path)
    description <- DBI::dbGetQuery(
        connection,
        paste("DESCRIBE SELECT * FROM", source)
    )
    columns <- as.character(description$column_name)
    quoted <- as.character(DBI::dbQuoteIdentifier(connection, columns))
    expressions <- quoted
    normalized <- match("Precursor.Normalised", columns)
    if (isTRUE(use_precursor_norm) && !is.na(normalized)) {
        expressions[[normalized]] <- paste0(
            "TRY_CAST(", quoted[[normalized]], " AS DOUBLE) AS ",
            quoted[[normalized]]
        )
    }
    list(
        source = source,
        columns = columns,
        expressions = expressions,
        sql = paste0(
            "SELECT row_number() OVER ()::BIGINT AS ",
            as.character(DBI::dbQuoteIdentifier(
                connection,
                .artifactRowOrderColumn
            )),
            ", ",
            paste(expressions, collapse = ", "),
            " FROM ", source
        )
    )
}

protDiaStreamingSanitizeRuns <- function(connection, select) {
    run_index <- match("Run", select$columns)
    if (is.na(run_index)) return(select)
    run_column <- as.character(DBI::dbQuoteIdentifier(connection, "Run"))
    observed <- DBI::dbGetQuery(
        connection,
        paste0(
            "SELECT DISTINCT ", run_column, " AS value FROM ",
            select$source, " ORDER BY value"
        )
    )$value
    cleaned <- janitor::make_clean_names(as.character(observed))
    clauses <- Map(function(before, after) {
        paste(
            "WHEN",
            as.character(DBI::dbQuoteString(connection, before)),
            "THEN",
            as.character(DBI::dbQuoteString(connection, after))
        )
    }, as.character(observed), cleaned)
    select$expressions[[run_index]] <- paste0(
        "CASE ", run_column, " ",
        paste(unlist(clauses, use.names = FALSE), collapse = " "),
        " ELSE ", run_column, " END AS ", run_column
    )
    select$sql <- paste0(
        "SELECT row_number() OVER ()::BIGINT AS ",
        as.character(DBI::dbQuoteIdentifier(
            connection,
            .artifactRowOrderColumn
        )),
        ", ", paste(select$expressions, collapse = ", "),
        " FROM ", select$source
    )
    select
}

protDiaStreamingImportSummary <- function(connection, select, mapping) {
    quote <- function(value) {
        as.character(DBI::dbQuoteIdentifier(connection, value))
    }
    query <- paste0(
        "SELECT count(*) AS rows, ",
        "list_sort(list_distinct(list(", quote(mapping$run_col),
        "))) AS runs, count(DISTINCT ", quote(mapping$protein_col),
        ") AS proteins, count(DISTINCT ", quote(mapping$peptide_col),
        ") AS peptides FROM (", select$sql, ")"
    )
    observed <- DBI::dbGetQuery(connection, query)
    list(
        rows = as.numeric(observed$rows[[1L]]),
        columns = as.integer(length(select$columns)),
        column_names = as.list(select$columns),
        run_values = as.list(as.character(observed$runs[[1L]])),
        protein_count = as.numeric(observed$proteins[[1L]]),
        peptide_count = as.numeric(observed$peptides[[1L]])
    )
}

protDiaStreamingColumnMapping <- function(columns, use_precursor_norm) {
    quantity <- if (isTRUE(use_precursor_norm) &&
        "Precursor.Normalised" %in% columns) {
        "Precursor.Normalised"
    } else {
        "Precursor.Quantity"
    }
    list(
        protein_col = "Protein.Group",
        peptide_col = "Stripped.Sequence",
        run_col = "Run",
        quantity_col = quantity,
        qvalue_col = "Q.Value",
        pg_qvalue_col = "PG.Q.Value"
    )
}

writeProtDiaStreamingParquet <- function(
    source_path,
    output_path,
    use_precursor_norm = TRUE,
    sanitize_names = FALSE,
    row_group_rows = 65536L,
    memory_limit_bytes = 512 * 1024^2
) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    database <- tempfile("prot-dia-streaming-", fileext = ".duckdb")
    temporary <- tempfile("prot-dia-spill-")
    dir.create(temporary)
    on.exit(unlink(c(database, temporary), recursive = TRUE, force = TRUE), add = TRUE)
    connection <- DBI::dbConnect(duckdb::duckdb(database, shared_home = FALSE))
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(connection, "SET threads = 1")
    DBI::dbExecute(connection, "SET preserve_insertion_order = true")
    DBI::dbExecute(
        connection,
        paste0("SET memory_limit = '", as.integer(memory_limit_bytes), "B'")
    )
    DBI::dbExecute(
        connection,
        paste0(
            "SET temp_directory = ",
            as.character(DBI::dbQuoteString(connection, temporary))
        )
    )
    select <- protDiaStreamingSelect(
        connection,
        source_path,
        use_precursor_norm
    )
    if (isTRUE(sanitize_names)) {
        select <- protDiaStreamingSanitizeRuns(connection, select)
    }
    mapping <- protDiaStreamingColumnMapping(
        select$columns,
        use_precursor_norm
    )
    summary <- protDiaStreamingImportSummary(connection, select, mapping)
    destination <- as.character(DBI::dbQuoteString(connection, output_path))
    compression <- if (arrow::codec_is_available("zstd")) "ZSTD" else "SNAPPY"
    statement <- paste0(
        "COPY (", select$sql, ") TO ", destination,
        " (FORMAT PARQUET, COMPRESSION ", compression, ", ROW_GROUP_SIZE ",
        as.integer(row_group_rows), ")"
    )
    DBI::dbExecute(connection, statement)
    list(
        path = output_path,
        column_mapping = mapping,
        import_summary = protDiaValidateImportSummary(summary)
    )
}

protDiaImportIsDeferred <- function(data_import_result) {
    is.list(data_import_result) &&
        isTRUE(data_import_result$parent_hydration_deferred) &&
        is.null(data_import_result$data)
}

protDiaImportSummary <- function(data_import_result) {
    if (!protDiaImportIsDeferred(data_import_result)) return(NULL)
    protDiaValidateImportSummary(data_import_result$import_summary)
}

protDiaImportRowCount <- function(data_import_result) {
    summary <- protDiaImportSummary(data_import_result)
    if (is.null(summary)) nrow(data_import_result$data) else summary$rows
}

protDiaImportColumnNames <- function(data_import_result) {
    summary <- protDiaImportSummary(data_import_result)
    if (is.null(summary)) names(data_import_result$data) else
        unlist(summary$column_names, use.names = FALSE)
}

protDiaHydrateDeferredImport <- function(artifact_import) {
    if (!is.list(artifact_import$pending_stage) ||
        !protDiaImportIsDeferred(artifact_import$result)) {
        return(artifact_import)
    }
    pending <- validateProtDiaArtifactPendingStage(
        artifact_import$pending_stage
    )
    tables <- readProtDiaArtifactStage(
        pending$stage$store,
        pending$stage
    )
    artifact_import$result$data <- tables$canonical_data
    artifact_import$result$parent_hydration_deferred <- FALSE
    artifact_import
}

protDiaDeferredWorkflowSummary <- function(workflow_data) {
    summary <- workflow_data$artifact_import_summary
    if (is.null(summary)) return(NULL)
    protDiaValidateImportSummary(summary)
}

protDiaDeferredDesignInput <- function(workflow_data) {
    summary <- protDiaDeferredWorkflowSummary(workflow_data)
    if (is.null(summary)) return(NULL)
    runs <- unlist(summary$run_values, use.names = FALSE)
    data.frame(
        Run = runs,
        .multischolar_original_run = runs,
        stringsAsFactors = FALSE
    )
}

resolveProtDesignBuilderPayload <- function(workflow_data) {
    deferred <- protDiaDeferredDesignInput(workflow_data)
    if (!is.null(deferred)) return(deferred)
    resolveProtDiaWorkflowTable(workflow_data, "data_tbl")
}

protDiaDeferredRunMapping <- function(data_cln) {
    original <- ".multischolar_original_run"
    if (!is.data.frame(data_cln) || !all(c("Run", original) %in% names(data_cln))) {
        return(NULL)
    }
    mapping <- unique(data_cln[c(original, "Run")])
    names(mapping)[[1L]] <- "original_run"
    if (anyNA(mapping) || anyDuplicated(mapping$original_run) > 0L ||
        anyDuplicated(mapping$Run) > 0L) {
        protDiaArtifactAbort(
            "DIA-NN deferred design run mapping is invalid",
            "multischolar_invalid_prot_dia_design_mapping"
        )
    }
    mapping
}

protDiaArtifactDistinctColumn <- function(workflow_data, logical_name) {
    ref <- protDiaWorkflowPayloadRef(workflow_data, "data_tbl")
    if (!is.list(ref)) return(NULL)
    identity <- workflow_data$workflow_context$getIdentity()
    store <- newArtifactStore(
        workflow_data$workflow_context$getPaths(),
        identity$project_id
    )
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(
        store,
        managed$sidecar,
        validate_payload = FALSE
    )
    metadata <- artifactValidateRectangularMetadata(sidecar$codec_metadata)
    index <- match(logical_name, metadata$logical_names)
    if (is.na(index)) {
        protDiaArtifactAbort(
            "DIA-NN bounded column is absent from the import artifact",
            "multischolar_invalid_prot_dia_payload_role"
        )
    }
    descriptor <- metadata$columns[[index]]
    selected <- c(descriptor$physical_name, descriptor$status_name)
    selected <- selected[!vapply(selected, is.null, logical(1))]
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, must_exist = TRUE),
        col_select = selected,
        as_data_frame = TRUE
    )
    value <- artifactDecodeColumn(payload, descriptor)
    unique(value)
}

protDiaDeferredAnnotationInput <- function(workflow_data) {
    mapping <- workflow_data$column_mapping
    proteins <- protDiaArtifactDistinctColumn(
        workflow_data,
        mapping$protein_col
    )
    if (is.null(proteins)) return(NULL)
    stats::setNames(
        data.frame(value = proteins, stringsAsFactors = FALSE),
        mapping$protein_col
    )
}
