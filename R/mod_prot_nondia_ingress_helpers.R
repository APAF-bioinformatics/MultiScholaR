# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Read one bounded non-DIA TSV header
#'
#' @param source_path Existing source file.
#' @param maximum_bytes Maximum accepted header bytes.
#'
#' @return Exact input column names.
#' @noRd
protNonDiaIngressHeader <- function(
    source_path,
    maximum_bytes = 1024^2
) {
    connection <- file(source_path, open = "rb")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header) != 1L || nchar(header, type = "bytes") > maximum_bytes) {
        protNonDiaArtifactAbort(
            "non-DIA ingress header is missing or exceeds its bound",
            "multischolar_invalid_prot_nondia_ingress_header"
        )
    }
    columns <- strsplit(header, "\t", fixed = TRUE)[[1L]]
    if (anyDuplicated(columns)) {
        protNonDiaArtifactAbort(
            "non-DIA ingress columns must be unique",
            "multischolar_invalid_prot_nondia_ingress_schema"
        )
    }
    columns
}

#' Resolve one case-insensitive exact column candidate
#'
#' @param columns Input column names.
#' @param candidates Ordered accepted candidates.
#'
#' @return The first exact or case-insensitive match, or `NULL`.
#' @noRd
protNonDiaIngressColumn <- function(columns, candidates) {
    exact <- candidates[candidates %in% columns]
    if (length(exact)) return(exact[[1L]])
    positions <- match(tolower(candidates), tolower(columns), nomatch = 0L)
    positions <- positions[positions > 0L]
    if (!length(positions)) NULL else columns[[positions[[1L]]]]
}

#' Resolve one bounded non-DIA import preflight
#'
#' @param source_path Existing source file.
#' @param format Exact LFQ format.
#' @param parser_parameters Validated parser parameters.
#'
#' @return A source-bound parser and projection contract.
#' @noRd
protNonDiaIngressPreflight <- function(
    source_path,
    format,
    parser_parameters
) {
    columns <- protNonDiaIngressHeader(source_path)
    if (identical(format, "maxquant")) {
        protein <- protNonDiaIngressColumn(columns, c(
            "Majority.protein.IDs", "Protein.IDs", "Protein.Ids",
            "Protein IDs", "Protein ID"
        ))
        lfq <- grep("^LFQ\\.intensity\\.", columns, value = TRUE)
        intensity <- grep("^Intensity\\.", columns, value = TRUE)
        quantities <- if (isTRUE(parser_parameters$use_lfq) && length(lfq)) {
            lfq
        } else {
            intensity
        }
        metadata <- intersect(c("Gene.names", "Protein.names"), columns)
        filters <- intersect(c("Potential.contaminant", "Reverse"), columns)
    } else {
        protein <- protNonDiaIngressColumn(columns, c(
            "Protein ID", "Protein.ID", "Protein_ID", "Protein",
            "protein id", "protein.id", "protein_id"
        ))
        all_intensity <- grep(
            "Intensity$",
            columns,
            value = TRUE,
            ignore.case = TRUE
        )
        maxlfq <- grep(
            "MaxLFQ.*Intensity$",
            all_intensity,
            value = TRUE,
            ignore.case = TRUE
        )
        regular <- setdiff(all_intensity, maxlfq)
        quantities <- if (isTRUE(parser_parameters$use_maxlfq) &&
            length(maxlfq)) {
            maxlfq
        } else {
            regular
        }
        metadata <- character()
        filters <- character()
    }
    if (is.null(protein) || !length(quantities)) {
        protNonDiaArtifactAbort(
            "non-DIA ingress schema lacks protein or quantity columns",
            "multischolar_invalid_prot_nondia_ingress_schema"
        )
    }
    list(
        source_size_bytes = unname(as.numeric(file.info(source_path)$size)),
        source_digest = artifactByteDigest(source_path),
        columns = columns,
        protein_column = protein,
        metadata_columns = metadata,
        quantity_columns = quantities,
        filter_columns = filters,
        complete_payload_materialized = FALSE
    )
}

#' Quote one DuckDB identifier
#'
#' @param connection Open DuckDB connection.
#' @param value Identifier value.
#'
#' @return A quoted SQL identifier.
#' @noRd
protNonDiaIngressQuote <- function(connection, value) {
    as.character(DBI::dbQuoteIdentifier(connection, value))
}

#' Build ordered wide-to-long LFQ SQL
#'
#' @param connection Open DuckDB connection.
#' @param source_path Existing source file.
#' @param format Exact LFQ format.
#' @param parser_parameters Validated parser parameters.
#' @param preflight Bounded source preflight.
#'
#' @return Canonical select SQL and column mapping.
#' @noRd
protNonDiaIngressSelect <- function(
    connection,
    source_path,
    format,
    parser_parameters,
    preflight
) {
    source_literal <- as.character(DBI::dbQuoteString(connection, source_path))
    source_table <- protNonDiaIngressQuote(
        connection,
        ".multischolar_nondia_source"
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
        column <- preflight$quantity_columns[[index]]
        run <- if (identical(format, "maxquant")) {
            sub("^(LFQ\\.intensity\\.|Intensity\\.)", "", column)
        } else {
            gsub(
                "\\s+(MaxLFQ\\s+)?Intensity$",
                "",
                column,
                ignore.case = TRUE
            )
        }
        paste0(
            "(", index, ", ",
            as.character(DBI::dbQuoteString(connection, run)), ", ",
            protNonDiaIngressQuote(connection, column), ")"
        )
    }, character(1))
    protein <- protNonDiaIngressQuote(connection, preflight$protein_column)
    projections <- c(paste0(protein, " AS ",
        protNonDiaIngressQuote(connection, "Protein.Ids")))
    projections <- c(projections, vapply(
        preflight$metadata_columns,
        \(column) protNonDiaIngressQuote(connection, column),
        character(1)
    ))
    projections <- c(
        projections,
        paste0("unpivoted.run AS ", protNonDiaIngressQuote(connection, "Run")),
        paste0(
            "TRY_CAST(unpivoted.intensity AS DOUBLE) AS ",
            protNonDiaIngressQuote(connection, "Intensity")
        )
    )
    filter_sql <- ""
    if (identical(format, "maxquant") &&
        isTRUE(parser_parameters$filter_contaminants) &&
        length(preflight$filter_columns)) {
        clauses <- vapply(preflight$filter_columns, \(column) {
            quoted <- protNonDiaIngressQuote(connection, column)
            paste0("(", quoted, " IS NULL OR ", quoted, " <> '+')")
        }, character(1))
        filter_sql <- paste("WHERE", paste(clauses, collapse = " AND "))
    }
    mapping <- if (identical(format, "maxquant")) {
        list(
            protein_col = "Protein.Ids",
            peptide_col = NULL,
            run_col = "Run",
            quantity_col = "Intensity",
            qvalue_col = "Q.value"
        )
    } else {
        list(
            protein_col = "Protein.Ids",
            run_col = "Run",
            quantity_col = "Intensity",
            qvalue_col = NULL
        )
    }
    list(
        sql = paste0(
            "SELECT row_number() OVER (ORDER BY ", source_order,
            ", unpivoted.ordinal)",
            "::BIGINT AS ", protNonDiaIngressQuote(
                connection,
                .artifactRowOrderColumn
            ), ", ", paste(projections, collapse = ", "), " FROM ",
            source_table, " CROSS JOIN LATERAL (VALUES ",
            paste(values, collapse = ", "),
            ") AS unpivoted(ordinal, run, intensity) ", filter_sql,
            " ORDER BY ", source_order, ", unpivoted.ordinal"
        ),
        column_mapping = mapping
    )
}

#' Stream one non-DIA LFQ source directly to canonical Parquet
#'
#' @param source_path Existing source file.
#' @param output_path Destination Parquet path.
#' @param format Exact LFQ format.
#' @param parser_parameters Validated parser parameters.
#' @param row_group_rows Parquet row-group size.
#' @param memory_limit_bytes DuckDB memory limit.
#'
#' @return Source binding, mapping, and output path.
#' @noRd
writeProtNonDiaStreamingParquet <- function(
    source_path,
    output_path,
    format,
    parser_parameters,
    row_group_rows = 65536L,
    memory_limit_bytes = 512 * 1024^2
) {
    preflight <- protNonDiaIngressPreflight(
        source_path,
        format,
        parser_parameters
    )
    database <- tempfile("prot-nondia-streaming-", fileext = ".duckdb")
    temporary <- tempfile("prot-nondia-spill-")
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
    DBI::dbExecute(connection, paste0(
        "SET temp_directory = ",
        as.character(DBI::dbQuoteString(connection, temporary))
    ))
    select <- protNonDiaIngressSelect(
        connection,
        source_path,
        format,
        parser_parameters,
        preflight
    )
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
