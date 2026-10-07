#' Reduce a sorted precursor stream with one live scientific group.
#' @param scope,input Existing native query owner and declared scientific roles.
#' @param path Owned temporary aggregate file, removed by the parent publication.
#' @return Reducer counters; creates qc_aggregated over the completed stage.
#' @noRd
prepareDIAOrderedAggregate <- function(scope, input, path) {
    con <- scope$connection
    roles <- input$roles
    keys <- unlist(roles[c("sample", "protein", "peptide")], use.names = FALSE)
    q <- function(value) as.character(workflowQueryQuoteIdentifier(con, value))
    prototype <- data.frame(sample = character(), protein = character(), peptide = character(),
        Peptide.RawQuantity = double(), Peptide.Normalised = double(), peptidoform_count = integer(),
        raw_observed = integer(), normalised_observed = integer(), source_count = integer())
    names(prototype)[seq_len(3L)] <- keys
    writer <- fastduckdb::open_parquet_writer_fast(prototype, path, n_threads = 8L)
    on.exit(fastduckdb::abort_parquet_writer_fast(writer), add = TRUE)
    reducer <- .Call(C_workflow_dia_ordered_new)
    on.exit(.Call(C_workflow_dia_ordered_close, reducer), add = TRUE)
    query <- paste0("SELECT ", paste(paste0("CAST(", q(keys), " AS VARCHAR) AS ",
        c("sample", "protein", "peptide")), collapse = ","),
        ", ", q(roles$raw), " AS raw, ", q(roles$normalised), " AS normalised, CASE WHEN ",
        "regexp_full_match(CAST(", q(roles$modified), " AS VARCHAR), '[ \\t\\r\\n]*') ",
        "THEN NULL ELSE CAST(", q(roles$modified), " AS VARCHAR) END AS form, ",
        'CAST(".multischolar_order" AS DOUBLE) AS source_order FROM qc_aggregate_input ',
        "ORDER BY sample, protein, peptide, source_order")
    result <- fastduckdb::send_query_fast(con, query)
    on.exit(fastduckdb::clear_query_fast(result), add = TRUE)
    append <- function(part) {
        if (!nrow(part)) return(invisible(NULL))
        names(part)[seq_len(3L)] <- keys
        attr(part, "ordered_reducer_statistics") <- NULL
        fastduckdb::append_parquet_writer_fast(writer, part)
        invisible(NULL)
    }
    repeat {
        batch <- fastduckdb::fetch_query_fast(result)
        if (!nrow(batch)) break
        append(.Call(C_workflow_dia_ordered_push, reducer, batch, FALSE))
    }
    final <- .Call(C_workflow_dia_ordered_push, reducer, batch, TRUE)
    statistics <- attr(final, "ordered_reducer_statistics")
    append(final)
    fastduckdb::clear_query_fast(result)
    fastduckdb::close_parquet_writer_fast(writer)
    workflowQueryExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_aggregated AS SELECT * FROM read_parquet(",
        workflowQueryQuoteString(con, path), ")"))
    statistics
}

#' Refuse unsupported experimental ordered-group schemas before computation.
#' @param scope,input Existing physical read scope and scientific roles.
#' @return NULL invisibly for declared character/factor grouping keys.
#' @noRd
assertDIAOrderedSchema <- function(scope, input) {
    keys <- unlist(input$roles[c("sample", "protein", "peptide")], use.names = FALSE)
    allowed <- vapply(scope$metadata$column_metadata[keys], function(spec)
        identical(spec$typeof, "character") || !is.null(spec$levels), logical(1))
    if (length(allowed) != 3L || any(!allowed) || anyDuplicated(keys))
        cli::cli_abort("ordered-group experiment requires distinct character or factor keys",
            class = "workflow_execution_unsupported")
    invisible(NULL)
}
