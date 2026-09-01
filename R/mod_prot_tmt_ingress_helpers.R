# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Normalize one Proteome Discoverer reporter column
#'
#' @param column Exact input column name.
#'
#' @return The established reporter channel and sample name.
#' @noRd
protTmtIngressReporterName <- function(column) {
    if (grepl("^Abundance: ", column)) {
        return(gsub(
            "Abundance: F[0-9]+: ([0-9]+[A-Z]?), (.+)",
            "\\1_\\2",
            column
        ))
    }
    sub("^Abundance\\.", "", column)
}

#' Resolve one bounded PD-TMT TSV projection
#'
#' @param source_path Existing reviewed TSV source.
#'
#' @return Source-bound reporter, metadata, and mapping contract.
#' @noRd
protTmtIngressPreflight <- function(source_path) {
    columns <- protNonDiaIngressHeader(source_path)
    protein_candidates <- c("Accession", "Master.Protein.Accessions")
    protein <- protein_candidates[protein_candidates %in% columns][1L]
    quantity_columns <- grep(
        "^Abundance: |^Abundance\\.",
        columns,
        value = TRUE
    )
    reporter_names <- vapply(
        quantity_columns,
        protTmtIngressReporterName,
        character(1)
    )
    duplicates <- unique(reporter_names[duplicated(reporter_names)])
    if (length(duplicates)) {
        protNonDiaArtifactAbort(
            paste(
                "duplicate TMT reporter channel/sample names after",
                "normalization"
            ),
            "multischolar_invalid_prot_tmt_ingress_schema",
            duplicate_reporters = duplicates
        )
    }
    selected <- grepl("^[0-9]+[A-Z]?(_|$)", reporter_names)
    quantity_columns <- quantity_columns[selected]
    reporter_names <- reporter_names[selected]
    if (is.na(protein) || !length(quantity_columns)) {
        protNonDiaArtifactAbort(
            "PD-TMT ingress lacks an accession or reporter quantity column",
            "multischolar_invalid_prot_tmt_ingress_schema"
        )
    }
    metadata_columns <- setdiff(columns, quantity_columns)
    metadata_names <- metadata_columns
    metadata_names[metadata_columns == protein] <- "Protein.Ids"
    if (anyDuplicated(metadata_names)) {
        protNonDiaArtifactAbort(
            "PD-TMT ingress metadata names collide after accession mapping",
            "multischolar_invalid_prot_tmt_ingress_schema"
        )
    }
    list(
        source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
        source_digest = artifactByteDigest(source_path),
        columns = columns,
        protein_column = protein,
        metadata_columns = metadata_columns,
        metadata_names = metadata_names,
        quantity_columns = quantity_columns,
        reporter_names = reporter_names,
        complete_payload_materialized = FALSE
    )
}

#' Build ordered wide-to-long PD-TMT SQL
#'
#' @param connection Open DuckDB connection.
#' @param source_path Existing reviewed TSV source.
#' @param preflight Validated bounded TMT preflight.
#'
#' @return Canonical select SQL and exact public column mapping.
#' @noRd
protTmtIngressSelect <- function(connection, source_path, preflight) {
    source_literal <- as.character(DBI::dbQuoteString(connection, source_path))
    source_table <- protNonDiaIngressQuote(
        connection,
        ".multischolar_tmt_source"
    )
    source_order <- protNonDiaIngressQuote(connection, ".source_row")
    DBI::dbExecute(connection, paste0(
        "CREATE TEMP TABLE ", source_table, " AS SELECT ",
        "row_number() OVER ()::BIGINT AS ", source_order,
        ", * FROM read_csv(",
        source_literal, ", delim = chr(9), header = true, auto_detect = true, ",
        "sample_size = -1, nullstr = ['', 'NA'], strict_mode = true)"
    ))
    values <- vapply(seq_along(preflight$quantity_columns), \(index) {
        paste0(
            "(", index, ", ",
            as.character(DBI::dbQuoteString(
                connection,
                preflight$reporter_names[[index]]
            )), ", CAST(",
            protNonDiaIngressQuote(
                connection,
                preflight$quantity_columns[[index]]
            ), " AS DOUBLE))"
        )
    }, character(1))
    metadata <- vapply(seq_along(preflight$metadata_columns), \(index) {
        source <- protNonDiaIngressQuote(
            connection,
            preflight$metadata_columns[[index]]
        )
        target <- protNonDiaIngressQuote(
            connection,
            preflight$metadata_names[[index]]
        )
        if (identical(source, target)) source else paste(source, "AS", target)
    }, character(1))
    projections <- c(
        metadata,
        paste0("unpivoted.run AS ", protNonDiaIngressQuote(connection, "Run")),
        paste0(
            "unpivoted.abundance AS ",
            protNonDiaIngressQuote(connection, "Abundance")
        )
    )
    list(
        sql = paste0(
            "SELECT row_number() OVER (ORDER BY ", source_order,
            ", unpivoted.ordinal)::BIGINT AS ",
            protNonDiaIngressQuote(connection, .artifactRowOrderColumn),
            ", ", paste(projections, collapse = ", "), " FROM ",
            source_table, " CROSS JOIN LATERAL (VALUES ",
            paste(values, collapse = ", "),
            ") AS unpivoted(ordinal, run, abundance) ORDER BY ",
            source_order, ", unpivoted.ordinal"
        ),
        column_mapping = list(
            protein_col = "Protein.Ids",
            run_col = "Run",
            quantity_col = "Abundance",
            batch_col = NULL
        )
    )
}

#' Stream one reviewed PD-TMT TSV directly to canonical Parquet
#'
#' @param source_path Existing reviewed TSV source.
#' @param output_path Destination Parquet path.
#' @param row_group_rows Parquet row-group size.
#' @param memory_limit_bytes DuckDB memory limit.
#'
#' @return Source binding, mapping, and output path.
#' @noRd
writeProtTmtStreamingParquet <- function(
    source_path,
    output_path,
    row_group_rows = 65536L,
    memory_limit_bytes = 128 * 1024^2
) {
    preflight <- protTmtIngressPreflight(source_path)
    database <- tempfile("prot-tmt-streaming-", fileext = ".duckdb")
    temporary <- tempfile("prot-tmt-spill-")
    dir.create(temporary)
    on.exit(
        unlink(c(database, temporary), recursive = TRUE, force = TRUE),
        add = TRUE
    )
    connection <- DBI::dbConnect(duckdb::duckdb(database, shared_home = FALSE))
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    DBI::dbExecute(connection, "SET threads = 1")
    DBI::dbExecute(connection, "SET preserve_insertion_order = true")
    DBI::dbExecute(
        connection,
        paste0("SET memory_limit = '", as.integer(memory_limit_bytes), "B'")
    )
    DBI::dbExecute(connection, paste0(
        "SET temp_directory = ",
        as.character(DBI::dbQuoteString(connection, temporary))
    ))
    select <- protTmtIngressSelect(connection, source_path, preflight)
    destination <- as.character(DBI::dbQuoteString(connection, output_path))
    compression <- if (arrow::codec_is_available("zstd")) "ZSTD" else "SNAPPY"
    DBI::dbExecute(connection, paste0(
        "COPY (", select$sql, ") TO ", destination,
        " (FORMAT PARQUET, COMPRESSION ", compression,
        ", ROW_GROUP_SIZE ", as.integer(row_group_rows), ")"
    ))
    list(
        path = output_path,
        preflight = preflight,
        column_mapping = select$column_mapping
    )
}
