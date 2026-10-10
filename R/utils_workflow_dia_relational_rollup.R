# DIA-025 qualified relational implementation; automatic selection remains gated.

#' Resolve the source-bound native provider and fastduckdb query capability.
#' @return Native provider DLL path and execution binding; not activation authority.
#' @noRd
diaRollupNativeInfo <- function() {
    fail <- function() cli::cli_abort("native DIA rollup requires the bound fastduckdb query runtime",
        class = "dia_rollup_unavailable")
    # The exact sum is defined by 80-bit long-double accumulation (moved here from the removed
    # extension check, CORE-015).
    if (!identical(Sys.info()[["sysname"]], "Linux") || !identical(Sys.info()[["machine"]], "x86_64") ||
        !isTRUE(.Machine$longdouble.digits == 64)) cli::cli_abort(
        "native DIA rollup requires Linux x86_64 with 80-bit long double", class = "dia_rollup_unavailable")
    writer <- workflowFastWriterInfo()
    path <- system.file("duckdb", "linux_amd64", "manifest.rds", package = "MultiScholaR")
    manifest <- readRDS(path)
    native <- manifest$native_query
    query_path <- system.file("query-manifest.rds", package = "fastduckdb")
    required <- c("open_query_fast", "send_query_fast", "fetch_query_fast", "clear_query_fast",
        "close_query_fast", "execute_query_fast", "query_data_fast", "register_query_native_fast")
    if (!isTRUE(manifest$supported) || is.null(native) || !nzchar(query_path) ||
        file.size(query_path) > 65536 || !all(required %in% getNamespaceExports("fastduckdb")) ||
        !identical(native$api, "duckdb-c-scope-v1") ||
        !identical(native$registration, "dia_ordered_sum_v1") ||
        !identical(native$unsigned_loading, FALSE) ||
        !identical(digest::digest(file = query_path, algo = "sha256"), native$query_manifest_sha256)) fail()
    query_binary <- system.file("libs", paste0("fastduckdb_query", .Platform$dynlib.ext), package = "fastduckdb")
    if (!nzchar(query_binary) || !identical(digest::digest(file = query_binary, algo = "sha256"),
        native$query_binary_sha256)) fail()
    capability <- readRDS(query_path)
    if (!identical(capability$api, native$api) || !identical(capability$registration_abi, 1L) ||
        !identical(as.character(utils::packageVersion("duckdb")), capability$duckdb_package) ||
        !identical(capability$threads, 8L) || !identical(capability$unsigned_extensions, FALSE)) fail()
    dependency <- Filter(\(value) identical(value$package, "duckdb"), manifest$dependencies)
    if (length(dependency) != 1L) fail()
    dependency <- dependency[[1L]]
    files <- names(dependency$files)
    if (!length(files) || any(grepl("(^/|[.][.])", files))) fail()
    actual <- vapply(file.path(system.file(package = "duckdb"), files),
        \(file) digest::digest(file = file, algo = "sha256"), character(1))
    if (!identical(unname(actual), unname(dependency$files))) fail()
    provider <- getLoadedDLLs()[["MultiScholaR"]][["path"]]
    binding <- digest::digest(list(provider = digest::digest(file = provider, algo = "sha256"),
        kernel = native$kernel_sha256, query = native$query_manifest_sha256,
        query_binary = native$query_binary_sha256,
        writer = writer$binary_sha256), algo = "sha256")
    list(path = normalizePath(provider), manifest = list(binary_sha256 = binding),
        kernel_sha256 = native$kernel_sha256, runtime_release_status = native$runtime_release_status)
}

#' Describe the P14 rollup profile, the only supported route (CORE-015).
#' @param profile Profile name; only "P14" is supported.
#' @return Minimal action settings, never serialized as a live engine object.
#' @noRd
diaRollupProfile <- function(profile = "P14") {
    if (!identical(profile, "P14")) cli::cli_abort("unsupported DIA rollup profile",
        class = "dia_rollup_unavailable")
    batches <- getOption("multischolar.dia_rollup_batches", 3L)
    if (!is.numeric(batches) || length(batches) != 1L || is.na(batches) ||
        !batches %in% c(3L, 6L)) cli::cli_abort("unsupported experimental P14 batch count",
            class = "dia_rollup_unavailable")
    list(name = profile, preserve_order = TRUE, streaming = FALSE, batches = as.integer(batches))
}

#' Open the action scope on the native fastduckdb query owner with the exact-sum provider.
#' @param store,input Existing validated owner and input specification.
#' @param profile Rollup profile; only "P14" is supported.
#' @return Scope whose normal release owns connection teardown.
#' @noRd
newDIARelationalScope <- function(store, input, profile = "P14") {
    settings <- diaRollupProfile(profile)
    artifact <- diaRollupNativeInfo()
    scope <- newWorkflowTableReadScope(store, input$reference, reader = "arrow")
    ready <- FALSE
    on.exit(if (!ready) releaseWorkflowTableReadScope(scope), add = TRUE)
    stopifnot(is.null(scope$action$connection), .Machine$longdouble.digits == 64L,
        as.character(utils::packageVersion("duckdb")) == "1.5.5")
    con <- fastduckdb::open_query_fast()
    class(con) <- "workflow_fast_query"
    scope$action$connection <- scope$connection <- con
    scope$dia_rollup_profile <- settings
    scope$dia_rollup_extension <- artifact$manifest$binary_sha256
    configureWorkflowBulkQueryScope(scope)
    workflowQueryExecute(con, "SET preserve_insertion_order = true")
    fastduckdb::register_query_native_fast(con, "MultiScholaR", "dia_ordered_sum_v1")
    ready <- TRUE
    scope
}

#' Validate the same declared-double scientific contract before aggregation.
#' @param scope,input Validated table owner and declared roles/evidence.
#' @param samples Complete sample selector or NULL.
#' @return Selected input row count; creates only the qc_input view.
#' @noRd
validateDIARelationalInput <- function(scope, input, samples) {
    assertWorkflowPhysicalTableScope(scope)
    con <- scope$connection
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    roles <- input$roles
    fields <- unique(c(unlist(roles), input$evidence, input$annotation))
    metadata <- scope$metadata
    if (any(!fields %in% metadata$columns) ||
        any(vapply(metadata$column_metadata[c(roles$raw, roles$normalised)],
            function(x) !identical(x$typeof, "double"), logical(1))))
        cli::cli_abort("saved rollup requires declared double-precision quantities")
    selector <- if (is.null(samples)) "" else paste0(" WHERE ", q(roles$sample), " IN (",
        paste(workflowQueryQuoteString(con, samples), collapse = ","), ")")
    workflowQueryExecute(con, paste0("CREATE TEMP VIEW qc_input AS SELECT ",
        paste(q(c(fields, ".multischolar_order")), collapse = ","), " FROM read_parquet(",
        workflowQueryQuoteString(con, scope$path), ")", selector))
    keys <- c(roles$sample, roles$protein, roles$peptide, roles$modified, roles$precursor)
    bad_keys <- vapply(keys, function(name) {
        blank <- if (identical(metadata$column_metadata[[name]]$typeof, "character"))
            paste0(" OR regexp_full_match(", q(name), ", '[ \\t\\r\\n]*')") else ""
        paste0("(", q(name), " IS NULL", blank, ")")
    }, character(1))
    if (nrow(workflowQueryCollect(con, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_keys, collapse = " OR "), " LIMIT 1")))) stop("missing or blank identity keys")
    quantities <- q(c(roles$raw, roles$normalised))
    bad_values <- paste0("(", quantities, " IS NOT NULL AND (NOT isfinite(", quantities,
        ") OR ", quantities, " < 0))")
    if (nrow(workflowQueryCollect(con, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_values, collapse = " OR "), " LIMIT 1")))) stop("negative or non-finite abundance")
    # The identity includes the sample, so complete-sample batches bound the check's grouping state.
    batches <- scope$dia_rollup_profile$batches %||% 1L
    if (hasWorkflowPartitionDuplicates(scope, "qc_input", roles$sample,
        c(roles$sample, roles$protein, roles$precursor), batches = batches))
        stop("duplicate precursor identities")
    if (length(input$evidence)) {
        fields <- paste(q(c(roles$protein, input$evidence)), collapse = ",")
        if (nrow(workflowQueryCollect(con, paste0("SELECT 1 FROM (SELECT DISTINCT ", fields,
            " FROM qc_input) GROUP BY ", q(roles$protein), " HAVING count(*) > 1 LIMIT 1"))))
            stop("conflicting frozen identification evidence")
    }
    as.integer(workflowQueryCollect(con, "SELECT count(*) n FROM qc_input")$n)
}

#' Stream a DuckDB query through Arrow rather than DBI's full R resultset.
#' @param con,query Owned connection and fixed quoted query.
#' @param consume One bounded-batch consumer.
#' @return NULL invisibly, after closing result and Arrow reader.
#' @noRd
streamDIARelationalQuery <- function(con, query, consume) {
    result <- fastduckdb::send_query_fast(con, query)
    on.exit(fastduckdb::clear_query_fast(result), add = TRUE)
    repeat {
        data <- fastduckdb::fetch_query_fast(result)
        if (!nrow(data)) break
        stopifnot(nrow(data) <= 10000L, as.numeric(object.size(data)) <= 67108864)
        consume(data)
    }
    invisible(NULL)
}

#' Collapse bounded distinct annotations with the established native token helper.
#' @param scope,input Current query scope and declared annotation roles.
#' @return Small annotation lookup only, never the assay.
#' @noRd
readDIARelationalAnnotations <- function(scope, input) {
    con <- scope$connection
    protein <- input$roles$protein
    key <- as.character(workflowQueryQuoteIdentifier(con, protein))
    out <- stats::setNames(data.frame(rep(NA_character_, 10000L), rep(NA_character_, 10000L)),
        c(protein, "Protein.Ids"))
    values <- rep(NA_character_, 10000L)
    count <- groups <- 0L
    current <- NULL
    finish <- function() {
        if (!count) return(invisible(NULL))
        groups <<- groups + 1L
        if (groups > 10000L) stop("annotation lookup exceeds key limit")
        out[groups, ] <<- list(current, collapseRollupProteinIds(values[seq_len(count)]))
        count <<- 0L
    }
    streamDIARelationalQuery(con, paste0("SELECT DISTINCT ", key,
        ', "Protein.Ids" FROM qc_input ORDER BY ', key, ', "Protein.Ids"'), function(part) {
        for (i in seq_len(nrow(part))) {
            next_key <- part[[protein]][[i]]
            if (count && !identical(current, next_key)) finish()
            current <<- next_key
            count <<- count + 1L
            if (count > 10000L) stop("annotation lookup exceeds values limit")
            values[[count]] <<- part$Protein.Ids[[i]]
        }
        stopifnot(as.numeric(object.size(out)) <= 67108864,
            as.numeric(object.size(values)) <= 67108864)
    })
    finish()
    out[seq_len(groups), , drop = FALSE]
}

#' Build the one lazy P14 aggregation on the native query owner.
#' @param scope,input Scientific input scope and roles.
#' @return NULL invisibly; qc_aggregated is a lazy view, not a temporary assay.
#' @noRd
prepareDIARelationalAggregate <- function(scope, input) {
    con <- scope$connection
    roles <- input$roles
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    group <- c(roles$sample, roles$protein, roles$peptide)
    # Blank factor labels are allowed by the native key check but not counted as forms.
    workflowQueryExecute(con, "DROP VIEW IF EXISTS qc_aggregated")
    workflowQueryExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_prepared AS SELECT *, CASE WHEN ",
        "regexp_full_match(CAST(", q(roles$modified), " AS VARCHAR), '[ \\t\\r\\n]*') ",
        "THEN NULL ELSE CAST(", q(roles$modified), " AS VARCHAR) END AS dia_modified FROM qc_aggregate_input"))
    query <- paste0("SELECT ", paste(q(group), collapse = ","),
        ", dia_r_sum(list(", q(roles$raw), ' ORDER BY ".multischolar_order")) AS "Peptide.RawQuantity",',
        " dia_r_sum(list(", q(roles$normalised), ' ORDER BY ".multischolar_order")) AS "Peptide.Normalised",',
        " list_unique(list(dia_modified)) AS peptidoform_count, count(", q(roles$raw),
        ") AS raw_observed, count(", q(roles$normalised), ") AS normalised_observed,",
        " count(*) AS source_count FROM qc_prepared GROUP BY ", paste(q(group), collapse = ","))
    workflowQueryExecute(con, paste0("CREATE TEMP VIEW qc_aggregated AS ", query))
    invisible(NULL)
}

#' Drain one ordered result file through bounded eight-thread fastduckdb reads.
#' @param path Private result Parquet with explicit consecutive output ordinals.
#' @param consume Existing codec/summary/writer callback.
#' @param rows Trusted row count recorded before closing the native query scope.
#' @return Number of bounded read calls after all readers have been released.
#' @noRd
streamDIARelationalResultFile <- function(path, consume, rows) {
    stopifnot(is.numeric(rows), length(rows) == 1L, is.finite(rows), rows >= 0, rows == floor(rows))
    # The worker owns this immutable staged file and consumes it serially.
    # Retain footer metadata across its pages, never across file lifetimes.
    fastduckdb::clear_parquet_cache()
    on.exit(fastduckdb::clear_parquet_cache(), add = TRUE)
    calls <- 0L
    if (!rows) return(calls)
    for (start in seq.int(1L, rows, by = 10000L)) {
        end <- min(rows, start + 9999L)
        part <- fastduckdb::read_parquet_fast(path, n_threads = 8L,
            lazy_strings = FALSE, filter = list(.dia_output_order = c(start, end)))
        if (is.unsorted(part$.dia_output_order)) part <- part[order(part$.dia_output_order), , drop = FALSE]
        stopifnot(identical(as.integer(part$.dia_output_order), seq.int(start, end)),
            nrow(part) <= 10000L, as.numeric(object.size(part)) <= 67108864)
        part$.dia_output_order <- NULL
        consume(part)
        calls <- calls + 1L
    }
    calls
}

#' Stream the ordered aggregation and joins into the existing staged writer.
#' @param scope,input,persist Existing scope, roles and native output schema/identity.
#' @param rows Selected precursor row count.
#' @return Existing result reference, exact summary and dependency/count structures.
#' @noRd
persistDIARelationalRollup <- function(scope, input, persist, rows) {
    con <- scope$connection
    roles <- input$roles
    persist$identity$software <- digest::digest(list(persist$identity$software,
        exact_sum = scope$dia_rollup_extension), algo = "sha256")
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    joins <- ""
    select <- paste0("a.", q(c(roles$sample, roles$protein, roles$peptide,
        "Peptide.RawQuantity", "Peptide.Normalised", "peptidoform_count")))
    if (length(input$evidence)) {
        fields <- paste(q(c(roles$protein, input$evidence)), collapse = ",")
        joins <- paste0(joins, " LEFT JOIN (SELECT DISTINCT ", fields, " FROM qc_input) e ON a.",
            q(roles$protein), " = e.", q(roles$protein))
        select <- c(select, paste0("e.", q(input$evidence)))
    }
    if (length(input$annotation)) {
        workflowQueryWriteLookup(con, "qc_annotations", readDIARelationalAnnotations(scope, input), temporary = TRUE)
        joins <- paste0(joins, " LEFT JOIN qc_annotations n ON a.", q(roles$protein), " = n.", q(roles$protein))
        select <- c(select, 'n."Protein.Ids"')
    }
    order <- vapply(c(roles$sample, roles$protein, roles$peptide), function(field) {
        column <- paste0("a.", q(field))
        levels <- scope$metadata$column_metadata[[field]]$levels
        if (is.null(levels)) return(column)
        paste0("CASE ", column, " ", paste(paste0("WHEN ", workflowQueryQuoteString(con, levels),
            " THEN ", seq_along(levels)), collapse = " "), " END")
    }, character(1))
    sample_order <- sub("a[.]", "", order[[1L]], fixed = FALSE)
    sample_values <- collectWorkflowBulkMetadata(scope, paste0("SELECT DISTINCT ",
        q(roles$sample), " FROM qc_input ORDER BY ", sample_order),
        max_rows = 10000L, max_bytes = 8388608)[[1L]]
    batch_sizes <- integer()
    max_group_rows <- 0L
    intermediate_bytes <- 0
    intermediate_bytes_written <- 0
    fastduckdb_calls <- 0L
    phase <- scope$dia_rollup_phase
    if (!is.function(phase)) phase <- \(label, fn) fn()
    # The native public-C connection owns computation and writes one private,
    # ordered stage per batch. fastduckdb file reads then preserve P14's 10,000-row
    # writer boundary without thousands of public result-chunk conversions.
    select <- c(select, paste0("row_number() OVER (ORDER BY ",
        paste(order, collapse = ","), ') AS ".dia_output_order"'))
    query <- paste0("SELECT ", paste(c(select, "a.raw_observed", "a.normalised_observed", "a.source_count"),
        collapse = ","), " FROM qc_aggregated a", joins, " ORDER BY ", paste(order, collapse = ","))
    metadata <- workflowSavedTableMetadata(persist$prototype, ".rownames", 3L)
    metadata$rows <- NULL
    counts <- newDIAFilteringAccumulator(persist$prototype, roles$protein)
    summary <- list(aggregation_rule = "linear_sum_modified_unmodified_and_charge_states",
        active_protein_key = roles$protein, precursor_identity_mode = "declared_precursor_id",
        input_precursor_rows = rows, unique_precursor_identities = rows,
        output_stripped_peptide_identities = 0L, raw_partial_missing_outputs = 0L,
        raw_all_missing_outputs = 0L, normalised_partial_missing_outputs = 0L,
        normalised_all_missing_outputs = 0L,
        frozen_identification_evidence = if (length(input$evidence)) "preserved" else "absent",
        protein_ids_provenance = if (length(input$annotation)) "preserved" else
            if (roles$protein == "Protein.Ids") "active_key" else "absent")
    reference <- saveWorkflowTable(scope$store, NULL, basename(tempfile("relational-rollup-")),
        persist$identity, metadata$key_column, metadata = metadata, write_fn = function(path, restoration) {
            writer <- openWorkflowTableWriter(path, restoration, scope$store$progress, backend = "fast_r",
                optimization = "none")
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            appendWorkflowTableWriter(writer, persist$prototype, integer())
            stage_paths <- character()
            on.exit(unlink(unique(stage_paths[file.exists(stage_paths)])), add = TRUE)
            record_stage <- function(file) {
                intermediate_bytes_written <<- intermediate_bytes_written + file.size(file)
                intermediate_bytes <<- max(intermediate_bytes,
                    sum(file.size(stage_paths[file.exists(stage_paths)])))
            }
            consume <- function(part) {
                summary$output_stripped_peptide_identities <<-
                    summary$output_stripped_peptide_identities + nrow(part)
                max_group_rows <<- max(max_group_rows, as.integer(part$source_count))
                for (quantity in c("raw", "normalised")) {
                    observed <- part[[paste0(quantity, "_observed")]]
                    partial <- paste0(quantity, "_partial_missing_outputs")
                    missing <- paste0(quantity, "_all_missing_outputs")
                    summary[[partial]] <<- summary[[partial]] +
                        sum(observed > 0 & observed < part$source_count)
                    summary[[missing]] <<- summary[[missing]] + sum(observed == 0)
                }
                output <- restoreSavedWorkflowTable(part, restoration, restoration$columns)
                counts$push(output)
                appendWorkflowTableWriter(writer, output)
            }
            # Complete-sample batches (KIT-006, the generic partitioner): each batch stages one ordered
            # aggregate file, the engine closes once, then every stage is read into the writer in order.
            partition <- partitionBySamples(sample_values, scope$dia_rollup_profile$batches,
                stage = function(sample_selection, batch_index) {
                    result_file <- file.path(dirname(path), paste0("ordered-result-", batch_index, ".parquet"))
                    stage_paths <<- unique(c(stage_paths, result_file))
                    filter <- if (!length(sample_selection)) "FALSE" else paste0(q(roles$sample), " IN (",
                        paste(workflowQueryQuoteString(con, sample_selection), collapse = ","), ")")
                    workflowQueryExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_aggregate_input AS ",
                        "SELECT * FROM qc_input WHERE ", filter))
                    prepareDIARelationalAggregate(scope, input)
                    phase("ordered_result_stage", \() workflowQueryExecute(con, paste0("COPY (", query,
                        ") TO ", workflowQueryQuoteString(con, result_file),
                        " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240)")))
                    record_stage(result_file)
                    result_rows <- as.numeric(workflowQueryCollect(con,
                        paste0("SELECT count(*) AS n FROM read_parquet(",
                            workflowQueryQuoteString(con, result_file), ")"))$n)
                    list(path = result_file, rows = result_rows)
                },
                between = function() {
                    phase("query_engine_close", \() fastduckdb::close_query_fast(con))
                    scope$connection <- scope$action$connection <- NULL
                },
                consume = function(result, batch_index) {
                    calls <- phase("staged_result_read_summary_write", \()
                        streamDIARelationalResultFile(result$path, consume, result$rows))
                    fastduckdb_calls <<- fastduckdb_calls + calls
                    unlink(result$path)
                })
            batch_sizes <<- partition$batch_sizes
            phase("final_writer_finish", \() finishWorkflowTableWriter(writer))
        })
    list(reference = reference, summary = summary, columns = metadata$columns,
        dependencies = list(input$reference), counts = counts$finish(),
        execution = list(profile = scope$dia_rollup_profile, samples = length(sample_values),
            query_owner = "fastduckdb_c", aggregation_frontend = "public_sql",
            output_transport = "staged_or_R",
            projected_row_bytes = workflowProjectedRowBytes(scope$metadata,
                unique(c(unlist(roles), input$evidence, input$annotation))),
            max_group_rows = max_group_rows, output_groups = summary$output_stripped_peptide_identities,
            staged_bytes = reference$payload$bytes, engine_threads = 8L,
            intermediate_bytes = intermediate_bytes,
            intermediate_bytes_written = intermediate_bytes_written,
            fastduckdb_calls = fastduckdb_calls, fastduckdb_threads = 8L,
            writer_backend = "fast_r", writer_optimization = "none", writer_threads = 8L,
            evidence_rows = 0, evidence_payload_estimate = 0,
            exact_sum_sha256 = scope$dia_rollup_extension,
            independence_unit = diaPrecursorRollupOperation()$unit,
            global_validation = "completed_before_batching",
            batch_sizes = batch_sizes))
}

#' Execute the P14 rollup through the common scientific/output boundary.
#' @param store,input,persist Existing production boundary values.
#' @param arm Must be "F", the native fastduckdb owner; arms C and D were removed in CORE-015.
#' @param samples Complete sample selection or NULL.
#' @param phase Natural-boundary timing callback.
#' @param extension Must be NULL; external extensions are not loaded.
#' @param profile Must be "P14".
#' @param worker_cap Optional admitted worker bytes for explicit guarded-adapter validation.
#' @return Existing accepted rollup descriptor and summary contract.
#' @noRd
computeDIARelationalRollup <- function(store, input, persist, arm = "F", samples = NULL,
    phase = function(label, fn) fn(), extension = NULL, profile = "P14", worker_cap = NULL) {
    if (!identical(arm, "F") || !is.null(extension)) cli::cli_abort(
        "only the native P14 rollup route is supported", class = "dia_rollup_unavailable")
    scope <- newDIARelationalScope(store, input, profile)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    scope$dia_rollup_phase <- phase
    rows <- phase("relational_validation", function() validateDIARelationalInput(scope, input, samples))
    admission <- if (!is.null(worker_cap)) phase("working_memory_admission", function()
        estimateDIARollupWorkingMemory(scope, input, worker_cap)) else NULL
    if (!is.null(admission)) scope$dia_rollup_profile$batches <- admission$batches
    result <- phase("relational_aggregate_order_summary_write", function()
        persistDIARelationalRollup(scope, input, persist, rows))
    result$execution$admission <- admission
    result
}

#' Bytes per input row for work over the whole input before batching (validation, admission,
#' layout). Calibrated 2026-10-08 (CORE-010) by a two-point fit of worker anonymous memory on
#' `full_p14` (24,000,000 rows): three batches from the CORE-017 traces (`core_005`, about 5.0
#' GiB after removing about 1.5 GiB of controller memory) and six batches from PERF-108's
#' unprofiled run (another build, quiet host, 5.04 GiB whole-slice). The fit gave about 2.1 GiB
#' outside the per-batch terms (61 bytes per row, rounded up to 64); its per-batch coefficient
#' (0.61) is not used, so the per-batch terms stay at their conservative weight of 1. The
#' kernel cgroup guard remains authoritative.
#' @noRd
diaRollupRowBytes <- 64

#' Estimate the bounded P14 working set before constructing ordered aggregate lists.
#' @param scope Validated relational scope with qc_input view.
#' @param input Scientific roles and bounded annotations.
#' @param worker_cap Actual admitted worker bytes; engine limits are separate.
#' @return Explicit conservative estimate and observed group/width statistics.
#' @noRd
estimateDIARollupWorkingMemory <- function(scope, input, worker_cap) {
    con <- scope$connection
    q <- function(x) as.character(workflowQueryQuoteIdentifier(con, x))
    roles <- input$roles
    keys <- c(roles$sample, roles$protein, roles$peptide)
    bytes <- function(fields) paste(vapply(fields, function(field) paste0(
        "coalesce(octet_length(encode(CAST(", q(field), " AS VARCHAR))),0)"), character(1)), collapse = " + ")
    groups <- function(predicate = "TRUE") paste0("SELECT ", paste(q(keys), collapse = ","),
        ", count(*) AS n, sum(32 + ", bytes(roles$modified), ") AS values_bytes, ",
        "max(256 + ", bytes(keys), ") AS key_bytes FROM qc_input WHERE ", predicate,
        " GROUP BY ", paste(q(keys), collapse = ","))
    statistics <- function(predicate = "TRUE", params = list())
        collectWorkflowBulkMetadata(scope, paste0("SELECT ", q(roles$sample),
        " AS sample, count(*) AS groups, ",
        "sum(n) AS rows, max(n) AS max_group_rows, max(values_bytes + key_bytes) AS max_group_bytes, ",
        "sum(values_bytes + key_bytes) AS state_bytes FROM (", groups(predicate), ") GROUP BY ",
        q(roles$sample)), params = params, max_rows = 10000L, max_bytes = 8388608)
    # Per-sample statistics are unchanged by grouping complete samples into the profile's batches.
    parts <- list()
    visitWorkflowCompletePartitions(scope, "qc_input", roles$sample, function(predicate, params)
        parts[[length(parts) + 1L]] <<- statistics(predicate, params),
        batches = scope$dia_rollup_profile$batches %||% 1L)
    stats <- do.call(rbind, parts)
    if (any(stats$max_group_bytes > 67108864)) cli::cli_abort(
        "indivisible rollup group exceeds 64 MiB qualified aggregate-state estimate",
        class = "workflow_execution_unsupported")
    annotation_width <- if (length(input$annotation)) as.numeric(collectWorkflowBulkMetadata(scope,
        paste0("SELECT coalesce(max(n),0) AS n FROM (SELECT ", q(roles$protein),
            ", sum(1 + ", bytes(input$annotation), ") AS n FROM (SELECT DISTINCT ",
            paste(q(c(roles$protein, input$annotation)), collapse = ","),
            " FROM qc_input) GROUP BY ", q(roles$protein), ")"),
        max_rows = 1L, max_bytes = 65536)$n) else 0
    # Upper bound any one of `b` complete-sample batches by the largest ceil(N/b) samples.
    # Include hash/list slack, materialized joins/order, R transfers and writer headroom, plus a
    # whole-input term for work done before batching (CORE-010 calibration, see diaRollupRowBytes).
    fixed <- 768 * 1024^2 + 2 * 67108864 + diaRollupRowBytes * sum(as.numeric(stats$rows))
    estimate_for <- function(b) fixed + 2 * completeGroupBatchBound(stats$state_bytes, b) +
        2 * completeGroupBatchBound(stats$groups, b) * (256 + annotation_width)
    # P14 takes the fewest batches (at least its profile's three) whose estimate fits the cap;
    # an explicit experimental count or another profile keeps its fixed count (KIT-006: the choice
    # is the generic complete-group admission; the estimate stays P14's CORE-010 calibration).
    minimum <- scope$dia_rollup_profile$batches
    adaptive <- identical(scope$dia_rollup_profile$name, "P14") &&
        is.null(getOption("multischolar.dia_rollup_batches"))
    counts <- if (adaptive) seq.int(minimum, max(minimum, nrow(stats))) else minimum
    admission <- tryCatch(chooseCompleteGroupBatches(counts, estimate_for, worker_cap,
        "P14-adaptive-complete-sample-batches-v2"), workflow_execution_unsupported = function(error)
        cli::cli_abort("rollup working-memory estimate exceeds admitted worker capacity",
            class = "workflow_execution_unsupported"))
    c(admission[c("estimated_bytes", "batches", "estimates", "worker_cap")],
        list(samples = nrow(stats), rows = sum(as.numeric(stats$rows)),
            output_groups = sum(as.numeric(stats$groups)), max_group_rows = max(c(0, stats$max_group_rows)),
            max_group_bytes = max(c(0, stats$max_group_bytes)), annotation_width = annotation_width),
        admission[c("policy", "limitation")])
}

#' The P14 precursor rollup as a declaration on the complete-sample partitioner (KIT-006). The
#' kernel is the ordered-sum SQL aggregate (`prepareDIARelationalAggregate`, exact against the R
#' kernel by the CORE P14 oracle), run through the unchanged reference entry point; phases and
#' `execution` fields are those the CORE verifier reads. Admission keeps the CORE-010 calibration.
#' @return Declaration from `newWorkflowOperation`.
#' @noRd
diaPrecursorRollupOperation <- function() {
    newWorkflowOperation(name = "precursor_rollup", input_class = "PeptideQuantitativeData",
        fields = "design_matrix", unit = "complete_sample",
        kernel = function(object, parameters, store) computeDIAReferencePrecursorRollup(store$manager, store$source,
            backend = "P14", worker_cap = parameters$worker_cap),
        output_class = "PeptideQuantitativeData", parameters = character(), product_type = "dia.precursor_rollup_result",
        global_statistics = c(protein_annotation_lookup = "Protein.Ids collapse spans every sample of a protein",
            frozen_evidence_lookup = "per-protein evidence join and conflict check over the whole input",
            global_validation = "duplicate, key and finite checks over the whole input before batching"),
        merge = c(design_matrix = "design", peptide_data = "concatenate", rollup_summary = "sum"))
}
