# DIA-025 qualified relational implementation; automatic selection remains gated.

#' Resolve the fixed, installed exact-sum artifact and its qualified runtime.
#' @return Path and reproducible build manifest, after capability/integrity checks.
#' @noRd
diaRollupExtensionInfo <- function() {
    fail <- function() cli::cli_abort("relational DIA rollup requires the qualified Linux exact-sum runtime",
        class = "dia_rollup_unavailable")
    if (!identical(Sys.info()[["sysname"]], "Linux") || !identical(Sys.info()[["machine"]], "x86_64") ||
        !isTRUE(.Machine$longdouble.digits == 64)) fail()
    root <- system.file("duckdb", "linux_amd64", package = "MultiScholaR")
    manifest_path <- file.path(root, "manifest.rds")
    path <- file.path(root, "dia_rsum.duckdb_extension")
    if (!nzchar(root) || !file.exists(manifest_path) || file.size(manifest_path) > 65536 || !file.exists(path)) fail()
    manifest <- readRDS(manifest_path)
    if (!identical(manifest$schema, 3L) || !isTRUE(manifest$supported) ||
        !identical(manifest$api, "v1.2.0") ||
        !identical(manifest$longdouble_digits, 64L) ||
        !identical(manifest$runtime_release_status, "held") ||
        !identical(as.character(utils::packageVersion("duckdb")), manifest$duckdb) ||
        !identical(as.character(utils::packageVersion("duckplyr")), manifest$duckplyr) ||
        !identical(digest::digest(file = path, algo = "sha256"), manifest$binary_sha256)) fail()
    if (length(manifest$dependencies) != 2L) fail()
    for (dependency in manifest$dependencies) {
        if (!dependency$package %in% c("duckdb", "duckplyr")) fail()
        dependency_root <- system.file(package = dependency$package)
        files <- names(dependency$files)
        if (!length(files) || any(grepl("(^/|[.][.])", files)) ||
            any(!file.exists(file.path(dependency_root, files)))) fail()
        actual <- vapply(files, function(file) digest::digest(
            file = file.path(dependency_root, file), algo = "sha256"), character(1))
        if (!identical(actual, dependency$files)) fail()
    }
    list(path = normalizePath(path), manifest = manifest)
}

#' Describe the frozen qualification profiles; automatic routing stays off.
#' @param profile Explicit diagnostic profile selected before an action.
#' @return Minimal action settings, never serialized as a live engine object.
#' @noRd
diaRollupProfile <- function(profile = "P0") {
    if (!is.character(profile) || length(profile) != 1L || is.na(profile) ||
        !profile %in% paste0("P", 0:16)) cli::cli_abort("unsupported DIA rollup profile",
            class = "dia_rollup_unavailable")
    list(name = profile, preserve_order = profile != "P1", streaming = profile == "P2",
        batches = if (profile %in% c("P3", "P7", paste0("P", 8:16))) 3L else
            if (profile == "P4") 6L else 1L)
}

#' Open the existing action scope with the isolated exact-sum extension.
#' @param store,input Existing validated owner and input specification.
#' @param extension Frozen local extension artifact.
#' @return Scope whose normal release owns connection teardown.
#' @noRd
newDIARelationalScope <- function(store, input, extension = NULL, profile = "P0") {
    artifact <- diaRollupExtensionInfo()
    if (!is.null(extension) && !identical(normalizePath(extension, mustWork = TRUE), artifact$path))
        cli::cli_abort("relational DIA cannot load an external extension", class = "dia_rollup_unavailable")
    extension <- artifact$path
    settings <- diaRollupProfile(profile)
    scope <- newWorkflowTableReadScope(store, input$reference, reader = "arrow")
    ready <- FALSE
    on.exit(if (!ready) releaseWorkflowTableReadScope(scope), add = TRUE)
    stopifnot(is.null(scope$action$connection), .Machine$longdouble.digits == 64L,
        as.character(utils::packageVersion("duckdb")) == "1.5.5", file.exists(extension))
    con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE,
        config = list(allow_unsigned_extensions = "true", temp_directory = "")))
    scope$action$connection <- scope$connection <- con
    scope$dia_rollup_profile <- settings
    scope$dia_rollup_extension <- artifact$manifest$binary_sha256
    configureWorkflowBulkQueryScope(scope)
    DBI::dbExecute(con, paste("SET preserve_insertion_order =", if (settings$preserve_order) "true" else "false"))
    DBI::dbExecute(con, paste("LOAD", DBI::dbQuoteString(con, extension)))
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
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    roles <- input$roles
    fields <- unique(c(unlist(roles), input$evidence, input$annotation))
    metadata <- scope$metadata
    if (any(!fields %in% metadata$columns) ||
        any(vapply(metadata$column_metadata[c(roles$raw, roles$normalised)],
            function(x) !identical(x$typeof, "double"), logical(1))))
        cli::cli_abort("saved rollup requires declared double-precision quantities")
    selector <- if (is.null(samples)) "" else paste0(" WHERE ", q(roles$sample), " IN (",
        paste(DBI::dbQuoteString(con, samples), collapse = ","), ")")
    DBI::dbExecute(con, paste0("CREATE TEMP VIEW qc_input AS SELECT ",
        paste(q(c(fields, ".multischolar_order")), collapse = ","), " FROM read_parquet(",
        DBI::dbQuoteString(con, scope$path), ")", selector))
    keys <- c(roles$sample, roles$protein, roles$peptide, roles$modified, roles$precursor)
    bad_keys <- vapply(keys, function(name) {
        blank <- if (identical(metadata$column_metadata[[name]]$typeof, "character"))
            paste0(" OR regexp_full_match(", q(name), ", '[ \\t\\r\\n]*')") else ""
        paste0("(", q(name), " IS NULL", blank, ")")
    }, character(1))
    if (nrow(DBI::dbGetQuery(con, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_keys, collapse = " OR "), " LIMIT 1")))) stop("missing or blank identity keys")
    quantities <- q(c(roles$raw, roles$normalised))
    bad_values <- paste0("(", quantities, " IS NOT NULL AND (NOT isfinite(", quantities,
        ") OR ", quantities, " < 0))")
    if (nrow(DBI::dbGetQuery(con, paste0("SELECT 1 FROM qc_input WHERE ",
        paste(bad_values, collapse = " OR "), " LIMIT 1")))) stop("negative or non-finite abundance")
    identity <- paste(q(c(roles$sample, roles$protein, roles$precursor)), collapse = ",")
    if (nrow(DBI::dbGetQuery(con, paste0("SELECT 1 FROM qc_input GROUP BY ", identity,
        " HAVING count(*) > 1 LIMIT 1")))) stop("duplicate precursor identities")
    if (length(input$evidence)) {
        fields <- paste(q(c(roles$protein, input$evidence)), collapse = ",")
        if (nrow(DBI::dbGetQuery(con, paste0("SELECT 1 FROM (SELECT DISTINCT ", fields,
            " FROM qc_input) GROUP BY ", q(roles$protein), " HAVING count(*) > 1 LIMIT 1"))))
            stop("conflicting frozen identification evidence")
    }
    as.integer(DBI::dbGetQuery(con, "SELECT count(*) n FROM qc_input")$n)
}

#' Stream a DuckDB query through Arrow rather than DBI's full R resultset.
#' @param con,query Owned connection and fixed quoted query.
#' @param consume One bounded-batch consumer.
#' @return NULL invisibly, after closing result and Arrow reader.
#' @noRd
streamDIARelationalQuery <- function(con, query, consume, streaming = FALSE) {
    if (streaming) {
        if (!requireNamespace("nanoarrow", quietly = TRUE)) cli::cli_abort(
            "streaming DIA rollup requires the qualified nanoarrow dependency", class = "dia_rollup_unavailable")
        result <- DBI::dbSendQueryArrow(con, query)
        on.exit(DBI::dbClearResult(result), add = TRUE)
        repeat {
            part <- DBI::dbFetchArrowChunk(result, chunk_size = 10000L)
            data <- as.data.frame(part)
            if (!nrow(data)) break
            stopifnot(nrow(data) <= 10000L, as.numeric(object.size(data)) <= 67108864)
            consume(data)
        }
        return(invisible(NULL))
    }
    result <- DBI::dbSendQuery(con, query, arrow = TRUE)
    on.exit(DBI::dbClearResult(result), add = TRUE)
    reader <- duckdb::duckdb_fetch_record_batch(result, chunk_size = 10000L)
    on.exit(reader$.unsafe_delete(), add = TRUE, after = FALSE)
    repeat {
        batch <- reader$read_next_batch()
        if (is.null(batch)) break
        data <- as.data.frame(batch)
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
    key <- as.character(DBI::dbQuoteIdentifier(con, protein))
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
    }, streaming = scope$dia_rollup_profile$streaming)
    finish()
    out[seq_len(groups), , drop = FALSE]
}

#' Build exactly one lazy aggregation using SQL or duckplyr's relational API.
#' @param scope,input Scientific input scope and roles.
#' @param arm C or D.
#' @return NULL invisibly; qc_aggregated is a lazy view, not a temporary assay.
#' @noRd
prepareDIARelationalAggregate <- function(scope, input, arm) {
    con <- scope$connection
    roles <- input$roles
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    group <- c(roles$sample, roles$protein, roles$peptide)
    # Blank factor labels are allowed by the native key check but not counted as forms.
    DBI::dbExecute(con, "DROP VIEW IF EXISTS qc_aggregated")
    DBI::dbExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_prepared AS SELECT *, CASE WHEN ",
        "regexp_full_match(CAST(", q(roles$modified), " AS VARCHAR), '[ \\t\\r\\n]*') ",
        "THEN NULL ELSE CAST(", q(roles$modified), " AS VARCHAR) END AS dia_modified FROM qc_aggregate_input"))
    if (arm == "C") {
        query <- paste0("SELECT ", paste(q(group), collapse = ","),
            ", dia_r_sum(list(", q(roles$raw), ' ORDER BY ".multischolar_order")) AS "Peptide.RawQuantity",',
            " dia_r_sum(list(", q(roles$normalised), ' ORDER BY ".multischolar_order")) AS "Peptide.Normalised",',
            " count(DISTINCT dia_modified) AS peptidoform_count, count(", q(roles$raw), ") AS raw_observed, count(",
            q(roles$normalised), ") AS normalised_observed, count(*) AS source_count FROM qc_prepared GROUP BY ",
            paste(q(group), collapse = ","))
        DBI::dbExecute(con, paste0("CREATE TEMP VIEW qc_aggregated AS ", query))
    } else {
        stopifnot(arm == "D", as.character(utils::packageVersion("duckplyr")) == "1.2.1")
        frame <- duckplyr::read_sql_duckdb("SELECT * FROM qc_prepared",
            con = scope$connection, prudence = "stingy")
        relation <- duckdb:::rel_from_altrep_df(frame, strict = FALSE, allow_materialized = FALSE)
        stopifnot(!is.null(relation))
        ref <- duckplyr::relexpr_reference
        fun <- duckplyr::relexpr_function
        sum_expr <- function(column, alias) fun("dia_r_sum", list(fun("list", list(ref(column)),
            order_bys = list(ref(".multischolar_order")))), alias = alias)
        aggregates <- list(sum_expr(roles$raw, "Peptide.RawQuantity"),
            sum_expr(roles$normalised, "Peptide.Normalised"),
            fun("list_unique", list(fun("list", list(ref("dia_modified")))), alias = "peptidoform_count"),
            fun("count", list(ref(roles$raw)), alias = "raw_observed"),
            fun("count", list(ref(roles$normalised)), alias = "normalised_observed"),
            fun("count_star", list(), alias = "source_count"))
        relation <- duckplyr::rel_aggregate(relation, groups = lapply(group, ref), aggregates = aggregates)
        # Pinned DuckDB bridge exports the lazy relation without collecting to R.
        duckdb:::rel_to_view(relation, "main", "qc_aggregated", TRUE)
    }
    invisible(NULL)
}

#' Prepare a projected sample-aligned file without collecting the assay in R.
#' @param scope,input Existing scientific owner and roles.
#' @param samples Canonically ordered complete sample identifiers.
#' @param path Private staging path owned by the existing writer request.
#' @return NULL after synchronous Parquet publication; immutable source ordinals are retained.
#' @noRd
prepareDIARelationalLayout <- function(scope, input, samples, path) {
    con <- scope$connection
    q <- \(x) as.character(DBI::dbQuoteIdentifier(con, x))
    fields <- unique(c(unlist(input$roles[c("sample", "protein", "peptide", "modified", "raw", "normalised")]),
        ".multischolar_order"))
    stopifnot(!".dia_sample_rank" %in% fields)
    ranks <- data.frame(value = samples, rank = seq_along(samples))
    DBI::dbWriteTable(con, "qc_sample_ranks", ranks, temporary = TRUE)
    query <- paste0("SELECT ", paste(paste0("s.", q(fields)), collapse = ","),
        ', r.rank AS ".dia_sample_rank" FROM qc_input s JOIN qc_sample_ranks r ON s.',
        q(input$roles$sample), ' = r.value ORDER BY r.rank, s.".multischolar_order"')
    DBI::dbExecute(con, paste0("COPY (", query, ") TO ", DBI::dbQuoteString(con, path),
        " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 122880)"))
    invisible(NULL)
}

#' Stage one complete sample cohort through the explicitly selected native engine.
#' @param con Owned DuckDB connection.
#' @param source,path Existing aligned file and private cohort output path.
#' @param ranks Inclusive complete-sample rank range.
#' @param fast Use fastduckdb's file-to-file streaming API at eight threads.
#' @return NULL after the cohort has been completely written.
#' @noRd
stageDIARelationalCohort <- function(con, source, path, ranks, fast) {
    if (fast) {
        if (!requireNamespace("fastduckdb", quietly = TRUE)) cli::cli_abort(
            "fast DIA staging requires the qualified fastduckdb package", class = "dia_rollup_unavailable")
        on.exit(fastduckdb::clear_parquet_cache(), add = TRUE)
        fastduckdb::stream_parquet_filter_fast(source, path, ".dia_sample_rank", ranks[[1L]],
            max_val = ranks[[2L]], compression = "snappy", n_threads = 8L)
    } else {
        DBI::dbExecute(con, paste0("COPY (SELECT * FROM read_parquet(", DBI::dbQuoteString(con, source),
            ') WHERE ".dia_sample_rank" BETWEEN ', ranks[[1L]], " AND ", ranks[[2L]], ") TO ",
            DBI::dbQuoteString(con, path), " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 262144)"))
    }
    invisible(NULL)
}

#' Drain one ordered result file through bounded, matched eight-thread readers.
#' @param con Owned DuckDB connection, also used by the matched control.
#' @param path Private result Parquet with explicit consecutive output ordinals.
#' @param consume Existing codec/summary/writer callback.
#' @param fast Select fastduckdb; both arms collect at most 10000 rows / 64 MiB.
#' @return Number of bounded read calls after all readers have been released.
#' @noRd
streamDIARelationalResultFile <- function(con, path, consume, fast) {
    if (fast && !requireNamespace("fastduckdb", quietly = TRUE)) cli::cli_abort(
        "fast DIA result reads require the qualified fastduckdb package", class = "dia_rollup_unavailable")
    quoted <- as.character(DBI::dbQuoteString(con, path))
    rows <- as.integer(DBI::dbGetQuery(con, paste0("SELECT count(*) n FROM read_parquet(", quoted, ")"))$n)
    calls <- 0L
    if (!rows) return(calls)
    for (start in seq.int(1L, rows, by = 10000L)) {
        end <- min(rows, start + 9999L)
        part <- if (fast) tryCatch(fastduckdb::read_parquet_fast(path, n_threads = 8L,
            lazy_strings = FALSE, filter = list(.dia_output_order = c(start, end))),
            finally = fastduckdb::clear_parquet_cache()) else
            DBI::dbGetQuery(con, paste0('SELECT * FROM read_parquet(', quoted,
                ') WHERE ".dia_output_order" BETWEEN ', start, " AND ", end,
                ' ORDER BY ".dia_output_order"'))
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
persistDIARelationalRollup <- function(scope, input, persist, rows, arm = "D") {
    con <- scope$connection
    roles <- input$roles
    persist$identity$software <- digest::digest(list(persist$identity$software,
        exact_sum = scope$dia_rollup_extension), algo = "sha256")
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    joins <- ""
    evidence_rows <- evidence_bytes <- 0
    select <- paste0("a.", q(c(roles$sample, roles$protein, roles$peptide,
        "Peptide.RawQuantity", "Peptide.Normalised", "peptidoform_count")))
    if (length(input$evidence)) {
        fields <- paste(q(c(roles$protein, input$evidence)), collapse = ",")
        evidence <- paste0("(SELECT DISTINCT ", fields, " FROM qc_input)")
        if (scope$dia_rollup_profile$name == "P7") {
            # Validation has established one evidence tuple per protein globally.
            # Only this bounded lookup is materialized, never an assay intermediate.
            DBI::dbExecute(con, paste0("CREATE TEMP TABLE qc_evidence AS SELECT DISTINCT ",
                fields, " FROM qc_input LIMIT 10001"))
            widths <- paste0("coalesce(octet_length(encode(CAST(",
                q(c(roles$protein, input$evidence)), " AS VARCHAR))), 0)")
            footprint <- DBI::dbGetQuery(con, paste0("SELECT count(*) AS rows, coalesce(sum(64 + ",
                paste(widths, collapse = " + "), "), 0) AS bytes FROM qc_evidence"))
            if (footprint$rows > 10000 || footprint$bytes > 67108864) cli::cli_abort(
                "DIA evidence lookup exceeds qualified limits", class = "dia_rollup_unavailable")
            evidence_rows <- as.numeric(footprint$rows)
            evidence_bytes <- as.numeric(footprint$bytes)
            evidence <- "qc_evidence"
        }
        joins <- paste0(joins, " LEFT JOIN ", evidence, " e ON a.",
            q(roles$protein), " = e.", q(roles$protein))
        select <- c(select, paste0("e.", q(input$evidence)))
    }
    if (length(input$annotation)) {
        DBI::dbWriteTable(con, "qc_annotations", readDIARelationalAnnotations(scope, input), temporary = TRUE)
        joins <- paste0(joins, " LEFT JOIN qc_annotations n ON a.", q(roles$protein), " = n.", q(roles$protein))
        select <- c(select, 'n."Protein.Ids"')
    }
    order <- vapply(c(roles$sample, roles$protein, roles$peptide), function(field) {
        column <- paste0("a.", q(field))
        levels <- scope$metadata$column_metadata[[field]]$levels
        if (is.null(levels)) return(column)
        paste0("CASE ", column, " ", paste(paste0("WHEN ", DBI::dbQuoteString(con, levels),
            " THEN ", seq_along(levels)), collapse = " "), " END")
    }, character(1))
    sample_order <- sub("a[.]", "", order[[1L]], fixed = FALSE)
    sample_values <- collectWorkflowBulkMetadata(scope, paste0("SELECT DISTINCT ",
        q(roles$sample), " FROM qc_input ORDER BY ", sample_order),
        max_rows = 10000L, max_bytes = 8388608)[[1L]]
    batches <- scope$dia_rollup_profile$batches
    selections <- partitionWorkflowCompleteGroups(sample_values, batches)
    max_group_rows <- 0L
    intermediate_bytes <- 0
    intermediate_bytes_written <- 0
    fastduckdb_calls <- 0L
    profile <- scope$dia_rollup_profile$name
    phase <- scope$dia_rollup_phase
    if (!is.function(phase)) phase <- \(label, fn) fn()
    staged_output <- profile %in% c("P11", "P12")
    if (staged_output) select <- c(select, paste0("row_number() OVER (ORDER BY ",
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
            writer_backend <- switch(profile, P13 = "fast_arrow", P14 = "fast_r", P15 = "fast_r", P16 = "fast_r", "arrow")
            writer <- openWorkflowTableWriter(path, restoration, scope$store$progress, backend = writer_backend,
                optimization = switch(profile, P15 = "metadata", P16 = "batch", "none"))
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            appendWorkflowTableWriter(writer, persist$prototype, integer())
            intermediate <- file.path(dirname(path), "aggregate-intermediate.parquet")
            aligned <- file.path(dirname(path), "sample-aligned.parquet")
            cohort <- file.path(dirname(path), "sample-cohort.parquet")
            result_file <- file.path(dirname(path), "ordered-result.parquet")
            stage_paths <- c(intermediate, aligned, cohort, result_file)
            on.exit(unlink(stage_paths[file.exists(stage_paths)]), add = TRUE)
            record_stage <- function(file) {
                intermediate_bytes_written <<- intermediate_bytes_written + file.size(file)
                intermediate_bytes <<- max(intermediate_bytes,
                    sum(file.size(stage_paths[file.exists(stage_paths)])))
            }
            aligned_input <- profile %in% c("P8", "P9", "P10")
            if (aligned_input && length(sample_values)) {
                phase("sample_layout_prepare", \() prepareDIARelationalLayout(scope, input, sample_values, aligned))
                record_stage(aligned)
            }
            for (batch_index in seq_along(selections)) {
                sample_selection <- selections[[batch_index]]
                cohort <- file.path(dirname(path), paste0("sample-cohort-", batch_index, ".parquet"))
                result_file <- file.path(dirname(path), paste0("ordered-result-", batch_index, ".parquet"))
                stage_paths <- c(intermediate, aligned, cohort, result_file)
                if (!length(sample_selection) && profile %in% paste0("P", 8:12)) next
                filter <- if (!length(sample_selection)) "FALSE" else paste0(q(roles$sample), " IN (",
                    paste(DBI::dbQuoteString(con, sample_selection), collapse = ","), ")")
                # qc_input already contains the selected population. P6 isolates
                # the cost of selecting it a second time in the single-batch arm.
                selector <- if (scope$dia_rollup_profile$name == "P6") "" else paste0(" WHERE ", filter)
                source_query <- paste0("SELECT * FROM qc_input", selector)
                if (aligned_input) {
                    ranks <- range(match(sample_selection, sample_values))
                    stopifnot(all(is.finite(ranks)))
                    source_query <- paste0("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, aligned),
                        ') WHERE ".dia_sample_rank" BETWEEN ', ranks[[1L]], " AND ", ranks[[2L]])
                    if (profile %in% c("P9", "P10")) {
                        if (file.exists(cohort)) unlink(cohort)
                        phase("sample_cohort_stage", \() stageDIARelationalCohort(con, aligned, cohort,
                            ranks, fast = profile == "P10"))
                        record_stage(cohort)
                        if (profile == "P10") fastduckdb_calls <<- fastduckdb_calls + 1L
                        source_query <- paste0("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, cohort), ")")
                    }
                }
                DBI::dbExecute(con, paste0("CREATE OR REPLACE TEMP VIEW qc_aggregate_input AS ", source_query))
                prepareDIARelationalAggregate(scope, input, arm)
                if (scope$dia_rollup_profile$name == "P5") {
                    DBI::dbExecute(con, paste0("COPY (SELECT * FROM qc_aggregated) TO ",
                        DBI::dbQuoteString(con, intermediate), " (FORMAT PARQUET, COMPRESSION SNAPPY)"))
                    intermediate_bytes <<- file.size(intermediate)
                    DBI::dbExecute(con, "DROP VIEW qc_aggregated")
                    DBI::dbExecute(con, paste0("CREATE TEMP VIEW qc_aggregated AS SELECT * FROM read_parquet(",
                        DBI::dbQuoteString(con, intermediate), ")"))
                }
                consume <- function(part) {
                summary$output_stripped_peptide_identities <<- summary$output_stripped_peptide_identities + nrow(part)
                max_group_rows <<- max(max_group_rows, as.integer(part$source_count))
                for (quantity in c("raw", "normalised")) {
                    observed <- part[[paste0(quantity, "_observed")]]
                    partial <- paste0(quantity, "_partial_missing_outputs")
                    missing <- paste0(quantity, "_all_missing_outputs")
                    summary[[partial]] <<- summary[[partial]] + sum(observed > 0 & observed < part$source_count)
                    summary[[missing]] <<- summary[[missing]] + sum(observed == 0)
                }
                output <- restoreSavedWorkflowTable(part, restoration, restoration$columns)
                counts$push(output)
                appendWorkflowTableWriter(writer, output)
                }
                if (staged_output) {
                    phase("ordered_result_stage", \() DBI::dbExecute(con, paste0("COPY (", query,
                        ") TO ", DBI::dbQuoteString(con, result_file),
                        " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240)")))
                    record_stage(result_file)
                    calls <- phase("staged_result_read_summary_write", \()
                        streamDIARelationalResultFile(con, result_file, consume, fast = profile == "P12"))
                    if (profile == "P12") fastduckdb_calls <<- fastduckdb_calls + calls
                    unlink(result_file)
                } else streamDIARelationalQuery(con, query, consume,
                    streaming = scope$dia_rollup_profile$streaming)
                if (file.exists(cohort)) unlink(cohort)
            }
            finishWorkflowTableWriter(writer)
        })
    list(reference = reference, summary = summary, columns = metadata$columns,
        dependencies = list(input$reference), counts = counts$finish(),
        execution = list(profile = scope$dia_rollup_profile, samples = length(sample_values),
            projected_row_bytes = workflowProjectedRowBytes(scope$metadata,
                unique(c(unlist(roles), input$evidence, input$annotation))),
            max_group_rows = max_group_rows, output_groups = summary$output_stripped_peptide_identities,
            staged_bytes = reference$payload$bytes, engine_threads = 8L,
            intermediate_bytes = intermediate_bytes,
            intermediate_bytes_written = intermediate_bytes_written,
            fastduckdb_calls = fastduckdb_calls,
            fastduckdb_threads = if (profile %in% c("P10", "P12")) 8L else 0L,
            writer_backend = switch(profile, P13 = "fast_arrow", P14 = "fast_r", P15 = "fast_r", P16 = "fast_r", "arrow"),
            writer_optimization = switch(profile, P15 = "metadata", P16 = "batch", "none"),
            writer_threads = if (profile %in% c("P13", "P14", "P15", "P16")) 8L else NA_integer_,
            evidence_rows = evidence_rows, evidence_payload_estimate = evidence_bytes,
            exact_sum_sha256 = scope$dia_rollup_extension,
            independence_unit = "complete_sample",
            global_validation = "completed_before_batching",
            batch_sizes = as.integer(lengths(selections))))
}

#' Execute one bulk C/D operation through the common scientific/output boundary.
#' @param store,input,persist Existing production boundary values.
#' @param arm C or D; samples are complete sample identifiers or NULL.
#' @param samples Complete sample selection or NULL.
#' @param phase Natural-boundary timing callback.
#' @param extension Frozen experimental sum binary.
#' @param worker_cap Optional admitted worker bytes for explicit guarded-adapter validation.
#' @return Existing accepted rollup descriptor and summary contract.
#' @noRd
computeDIARelationalRollup <- function(store, input, persist, arm = "D", samples = NULL,
    phase = function(label, fn) fn(), extension = NULL, profile = "P0", worker_cap = NULL) {
    scope <- newDIARelationalScope(store, input, extension, profile)
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    scope$dia_rollup_phase <- phase
    rows <- phase("relational_validation", function() validateDIARelationalInput(scope, input, samples))
    admission <- if (!is.null(worker_cap)) phase("working_memory_admission", function()
        estimateDIARollupWorkingMemory(scope, input, worker_cap)) else NULL
    result <- phase("relational_aggregate_order_summary_write", function()
        persistDIARelationalRollup(scope, input, persist, rows, arm))
    result$execution$admission <- admission
    result
}

#' Estimate the bounded P14 working set before constructing ordered aggregate lists.
#' @param scope Validated relational scope with qc_input view.
#' @param input Scientific roles and bounded annotations.
#' @param worker_cap Actual admitted worker bytes; engine limits are separate.
#' @return Explicit conservative estimate and observed group/width statistics.
#' @noRd
estimateDIARollupWorkingMemory <- function(scope, input, worker_cap) {
    con <- scope$connection
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    roles <- input$roles
    keys <- c(roles$sample, roles$protein, roles$peptide)
    bytes <- function(fields) paste(vapply(fields, function(field) paste0(
        "coalesce(octet_length(encode(CAST(", q(field), " AS VARCHAR))),0)"), character(1)), collapse = " + ")
    groups <- paste0("SELECT ", paste(q(keys), collapse = ","),
        ", count(*) AS n, sum(32 + ", bytes(roles$modified), ") AS values_bytes, ",
        "max(256 + ", bytes(keys), ") AS key_bytes FROM qc_input GROUP BY ", paste(q(keys), collapse = ","))
    stats <- collectWorkflowBulkMetadata(scope, paste0("SELECT ", q(roles$sample),
        " AS sample, count(*) AS groups, ",
        "sum(n) AS rows, max(n) AS max_group_rows, max(values_bytes + key_bytes) AS max_group_bytes, ",
        "sum(values_bytes + key_bytes) AS state_bytes FROM (", groups, ") GROUP BY ",
        q(roles$sample)), max_rows = 10000L, max_bytes = 8388608)
    if (any(stats$max_group_bytes > 67108864)) cli::cli_abort(
        "indivisible rollup group exceeds 64 MiB qualified aggregate-state estimate",
        class = "workflow_execution_unsupported")
    annotation_width <- if (length(input$annotation)) as.numeric(collectWorkflowBulkMetadata(scope,
        paste0("SELECT coalesce(max(n),0) AS n FROM (SELECT ", q(roles$protein),
            ", sum(1 + ", bytes(input$annotation), ") AS n FROM (SELECT DISTINCT ",
            paste(q(c(roles$protein, input$annotation)), collapse = ","),
            " FROM qc_input) GROUP BY ", q(roles$protein), ")"),
        max_rows = 1L, max_bytes = 65536)$n) else 0
    # Upper bound any one of the three complete-sample batches by the largest ceil(N/3)
    # samples. Include hash/list slack, materialized joins/order, R transfers and writer headroom.
    n <- ceiling(nrow(stats) / 3)
    batch <- function(values) sum(utils::head(sort(as.numeric(values), decreasing = TRUE), n))
    estimate <- 768 * 1024^2 + 2 * batch(stats$state_bytes) +
        2 * batch(stats$groups) * (256 + annotation_width) + 2 * 67108864
    if (!is.numeric(worker_cap) || length(worker_cap) != 1L || !is.finite(worker_cap) ||
        worker_cap <= 0 || estimate > worker_cap) cli::cli_abort(
        "rollup working-memory estimate exceeds admitted worker capacity",
        class = "workflow_execution_unsupported")
    list(estimated_bytes = estimate, worker_cap = worker_cap, samples = nrow(stats),
        output_groups = sum(as.numeric(stats$groups)), max_group_rows = max(c(0, stats$max_group_rows)),
        max_group_bytes = max(c(0, stats$max_group_bytes)), annotation_width = annotation_width,
        policy = "P14-three-batches-state-and-output-slack-v1",
        limitation = "admission estimate; kernel cgroup guard remains authoritative")
}
