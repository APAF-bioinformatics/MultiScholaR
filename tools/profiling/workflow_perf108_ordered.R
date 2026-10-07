# Fixed P14 follow-up probes. These are diagnostic phases, not full-operation claims.

#' Compare complete aggregation and staging with one ordered-group reducer.
#' @param stage RS0 for retained list aggregation or RS1 for ordered reduction.
#' @param run_dir,repo Existing guarded member and source snapshot.
#' @return Exact-output receipt covering all sort/reduction/staging work.
#' @noRd
runPerf108OrderedProbe <- function(stage, run_dir, repo) {
    stopifnot(stage %in% c("RS0", "RS1"))
    fixture <- readRDS(file.path(dirname(dirname(run_dir)), "structural", "structural-fixture.rds"))
    helpers <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(repo, "tools/profiling/workflow_perf108.R"), helpers)
    sys.source(file.path(repo, "tools/profiling/workflow_perf108_structural.R"), helpers)
    finish <- helpers$startPerf108Attribution(run_dir)
    on.exit(finish(), add = TRUE)
    saveRDS(list(epoch = as.numeric(Sys.time()), uptime = as.numeric(strsplit(
        readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])), file.path(run_dir, "clock-anchor.rds"))
    receipt <- list(status = "failed", stage = stage, library = find.package("MultiScholaR"),
        fixture_sha256 = digest::digest(fixture, algo = "sha256"), engine_threads = 8L)
    on.exit(saveRDS(receipt, file.path(run_dir, "component-result.rds")), add = TRUE)
    owner <- readSavedStoreIdentity(file.path(fixture$project, ".workflow-entities/store.json"))
    store <- newWorkflowEntityStore(fixture$project, owner$project_id, owner$workflow_id)
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    scope <- newDIARelationalScope(store, fixture$input, profile = "P14", native = TRUE)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
    assertDIAOrderedSchema(scope, fixture$input)
    stopifnot(validateDIARelationalInput(scope, fixture$input, fixture$samples) == fixture$cohort_input_rows)
    con <- scope$connection
    workflowQueryExecute(con, "CREATE TEMP VIEW qc_aggregate_input AS SELECT * FROM qc_input")
    keys <- unlist(fixture$input$roles[c("sample", "protein", "peptide")], use.names = FALSE)
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    order <- paste(q(keys), collapse = ",")
    counts <- c("peptidoform_count", "raw_observed", "normalised_observed", "source_count")
    # Normalize auxiliary count storage types, preserving exact integer values.
    query <- paste0("SELECT ", paste(q(c(keys, "Peptide.RawQuantity", "Peptide.Normalised")), collapse = ","),
        ",", paste(paste0("CAST(", q(counts), " AS BIGINT) AS ", q(counts)), collapse = ","),
        ", row_number() OVER (ORDER BY ", order, ') AS ".probe_order" FROM qc_aggregated ORDER BY ', order)
    path <- file.path(run_dir, "aggregate.parquet")
    intermediate <- file.path(run_dir, "ordered-group-stage.parquet")
    on.exit(unlink(intermediate), add = TRUE)
    phase <- getOption("multischolar.dia_rollup_phase")
    phase("aggregation_and_stage", function() {
        if (stage == "RS1") {
            receipt$reducer <<- prepareDIAOrderedAggregate(scope, fixture$input, intermediate)
            receipt$extra_stage_bytes <<- file.size(intermediate)
        } else prepareDIARelationalAggregate(scope, fixture$input, "F")
        workflowQueryExecute(con, paste0("COPY (", query, ") TO ", workflowQueryQuoteString(con, path),
            " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240)"))
        workflowQueryExecute(con, "DROP VIEW qc_aggregated")
        if (stage == "RS1") unlink(intermediate)
    })
    fastduckdb::close_query_fast(con)
    scope$connection <- scope$action$connection <- NULL
    receipt$fingerprints <- helpers$fingerprintPerf108ProbeFile(path, ".probe_order", fixture$cohort_output_rows)
    receipt$status <- "complete"
    receipt
}

#' Compare persistent decoding with matched range reads on the retained cohort.
#' @param stage RC0 for range reads or RC1 for the repaired persistent reader.
#' @param run_dir,repo Existing guarded member and source snapshot.
#' @return Exact-output receipt with phase boundaries and native mechanism counters.
#' @noRd
runPerf108ReaderProbe <- function(stage, run_dir, repo) {
    stopifnot(stage %in% c("RC0", "RC1"))
    fixture_path <- file.path(dirname(dirname(run_dir)), "structural", "structural-fixture.rds")
    fixture <- readRDS(fixture_path)
    stopifnot(identical(digest::digest(file = fixture$output_path, algo = "sha256"), fixture$output_sha256))
    helpers <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(repo, "tools/profiling/workflow_perf108.R"), helpers)
    sys.source(file.path(repo, "tools/profiling/workflow_perf108_structural.R"), helpers)
    finish <- helpers$startPerf108Attribution(run_dir)
    on.exit(finish(), add = TRUE)
    saveRDS(list(epoch = as.numeric(Sys.time()), uptime = as.numeric(strsplit(
        readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])), file.path(run_dir, "clock-anchor.rds"))
    receipt <- list(status = "failed", stage = stage, library = find.package("MultiScholaR"),
        fixture_sha256 = digest::digest(fixture, algo = "sha256"), threads = 8L)
    on.exit(saveRDS(receipt, file.path(run_dir, "component-result.rds")), add = TRUE)
    metadata <- readSavedWorkflowTableMetadata(fixture$output_path)
    prototype <- as.data.frame(lapply(metadata$column_metadata, function(x) vector(x$typeof, 0L)), check.names = FALSE)
    prototype <- restoreSavedWorkflowTable(prototype, metadata, metadata$columns)
    metadata["rows"] <- list(NULL)
    path <- file.path(run_dir, "output.parquet")
    writer <- openWorkflowTableWriter(path, metadata, backend = "fast_r")
    on.exit(abortWorkflowTableWriter(writer), add = TRUE)
    appendWorkflowTableWriter(writer, prototype)
    counts <- newDIAFilteringAccumulator(prototype, fixture$input$roles$protein)
    phase <- getOption("multischolar.dia_rollup_phase")
    phase("reader_summary_write", function() {
        fastduckdb::clear_parquet_cache()
        on.exit(fastduckdb::clear_parquet_cache(), add = TRUE)
        cursor <- if (stage == "RC1") fastduckdb::open_parquet_cursor_fast(fixture$output_path, n_threads = 8L) else NULL
        if (!is.null(cursor)) on.exit(fastduckdb::close_parquet_cursor_fast(cursor), add = TRUE)
        for (start in seq.int(1L, fixture$cohort_output_rows, by = 10000L)) {
            part <- if (is.null(cursor)) fastduckdb::read_parquet_fast(fixture$output_path,
                n_threads = 8L, lazy_strings = FALSE,
                filter = list(.multischolar_order = c(start, start + 9999L))) else
                fastduckdb::fetch_parquet_cursor_fast(cursor)
            stopifnot(identical(as.integer(part$.multischolar_order), start:(start + 9999L)))
            data <- restoreSavedWorkflowTable(part, metadata, metadata$columns)
            counts$push(data)
            appendWorkflowTableWriter(writer, data)
        }
        if (!is.null(cursor)) receipt$cursor <<- fastduckdb::parquet_cursor_info_fast(cursor)
        finishWorkflowTableWriter(writer)
    })
    receipt$counts <- counts$finish()
    receipt$fingerprints <- helpers$fingerprintPerf108ProbeFile(path, ".multischolar_order",
        fixture$cohort_output_rows, metadata = readSavedWorkflowTableMetadata(path))
    stopifnot(identical(receipt$fingerprints, fixture$output_fingerprints))
    receipt$status <- "complete"
    receipt
}
