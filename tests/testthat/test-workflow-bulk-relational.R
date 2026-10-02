test_that("exact relational rollups preserve native science and release resources", {
    materialize_only <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:materialize-preflight"
    adaptation_only <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:adaptation-preflight"
    fast_only <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:fast-p3-preflight"
    writer_r <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:writer-r-preflight"
    writer_only <- writer_r || Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:writer-preflight"
    skip_if(!materialize_only && !adaptation_only && !fast_only && !writer_only && Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") != "gui-bulk-rollup:relational-preflight",
        "explicit isolated extension preflight only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    withr::local_libpaths(c(file.path(proof, "bulk-rollup-library"), .libPaths()))
    withr::local_envvar(c(DUCKPLYR_FALLBACK_AUTOUPLOAD = "0", DUCKPLYR_FALLBACK_COLLECT = "0"))
    artifact <- diaRollupExtensionInfo()
    extension <- artifact$path
    wfx_runtime <- writer_r && artifact$manifest$schema >= 2L
    base <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    invoke <- function(data) rollUpPrecursorToPeptideHelper(data, sample_id_column = "Run",
        protein_id_column = "Protein.Group", peptide_sequence_column = "Stripped.Sequence",
        modified_peptide_sequence_column = "Modified.Sequence", precursor_quantity_column = "Precursor.Quantity",
        precursor_normalised_column = "Precursor.Normalised", precursor_id_column = "Precursor.Id",
        core_utilisation = NA, return_rollup_result = TRUE)
    table <- base@peptide_data
    table$identification_peptide_count <- 3L
    table$identification_peptidoform_count <- 6L
    table$peptides_for_protein_count <- 3L
    table$peptidoforms_for_protein_count <- 6L
    cases <- list(ordinary = table, empty = table[0, , drop = FALSE])
    modified <- table
    modified$Precursor.Quantity[1:2] <- NA_real_
    modified$Precursor.Normalised[3] <- NA_real_
    modified$Precursor.Quantity[4] <- 0
    modified$Run <- factor(modified$Run, levels = c(rev(unique(modified$Run)), "unused"), ordered = TRUE)
    attr(modified$Protein.Group, "identity_role") <- "quantitative protein group"
    modified$Protein.Ids[1:3] <- c("P1; P2", " P2;P3;", "P1; P2")
    cases$missing_factors_annotations <- modified
    threshold <- table[rep(1L, 4L), , drop = FALSE]
    threshold$Precursor.Id <- paste0("boundary-", 1:4)
    threshold$Precursor.Quantity <- threshold$Precursor.Normalised <- c(1e16, 1, 1e16, 2)
    cases$boundary <- threshold
    threshold <- threshold[1:3, , drop = FALSE]
    threshold$Precursor.Quantity <- threshold$Precursor.Normalised <- c(1e16, 1, 1)
    cases$threshold <- threshold
    receipt <- list(status = "running", cases = list(), lifecycle = list(), phases = list(), failures = list())
    path <- file.path(root, "relational-preflight.json")
    flush <- function() jsonlite::write_json(receipt, path, pretty = TRUE, auto_unbox = TRUE, null = "null")
    on.exit(flush(), add = TRUE)
    phase <- function(label, fn) {
        started <- proc.time()[["elapsed"]]
        value <- tryCatch(fn(), error = function(error) {
            receipt$failures[[length(receipt$failures) + 1L]] <<- list(phase = label, error = conditionMessage(error))
            flush()
            stop(error)
        })
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(phase = label,
            seconds = proc.time()[["elapsed"]] - started)
        flush()
        value
    }
    tryCatch({
        identity <- list(input = "fixture", design = "paired", parameters = "exact-native",
            selectors = "all", software = "extension-v1")
        if (writer_r) {
            for (optimization in if (wfx_runtime) "none" else c("none", "metadata", "batch")) for (partial in c(FALSE, TRUE)) {
                owned <- file.path(root, paste(optimization, partial, sep = "-"))
                dir.create(owned)
                store <- newWorkflowEntityStore(owned, basename(owned), "DIA")
                on.exit(closeWorkflowEntityStore(store), add = TRUE)
                table <- data.frame(value = c(NA_real_, seq_len(63L)),
                    group = factor(rep(c("a", "b"), 32L), levels = c("a", "b", "unused")),
                    annotation = rep(c("alpha", "beta"), 32L))
                names(table$value) <- paste0("row-", seq_len(64L))
                attr(table$value, "unit") <- "linear"
                attr(table, "provenance") <- "bounded non-rollup writer consumer"
                expected <- serialize(table, NULL, version = 3)
                metadata <- workflowSavedTableMetadata(table, ".rownames", 3L)
                reference <- saveWorkflowTable(store, NULL, "consumer", identity, ".rownames",
                    metadata = metadata, write_fn = function(path, restoration) {
                        writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r",
                            optimization = optimization)
                        on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                        if (partial) {
                            appendWorkflowTableWriter(writer, table, 1:32)
                            appendWorkflowTableWriter(writer, table, 33:64)
                        } else appendWorkflowTableWriter(writer, table)
                        table$value[] <- -10
                        table$annotation[] <- "changed after append"
                        table <- NULL
                        invisible(gc())
                        if (optimization == "metadata") {
                            encoded <- restoration
                            encoded["rows"] <- list(NULL)
                            expect_identical(writer$schema_metadata$multischolar_table,
                                jsonlite::base64_enc(serialize(encoded, NULL, version = 3)))
                        }
                        finishWorkflowTableWriter(writer)
                    })
                page <- readSavedWorkflowTable(store, reference, metadata$columns, limit = 64L)
                expect_identical(page$source_row_names, as.character(seq_len(64L)))
                expect_identical(page$source_order, seq_len(64L))
                # Projected pages preserve source row labels. The existing output
                # codec restores the original container's automatic row-name type.
                actual <- completeWorkflowTableBatch(page$data, metadata)
                expect_identical(actual, unserialize(expected))
                expect_identical(store$open_resources, 0L)
                closeWorkflowEntityStore(store)
            }
        }
        if (writer_only) {
            for (backend in if (wfx_runtime) "fast_r" else c("fast_arrow", "fast_r")) {
                owned <- file.path(root, paste0("writer-api-", backend))
                dir.create(owned)
                store <- newWorkflowEntityStore(owned, paste0("writer-api-", backend), "DIA")
                on.exit(closeWorkflowEntityStore(store), add = TRUE)
                table <- cases$missing_factors_annotations[1:64, , drop = FALSE]
                names(table$Precursor.Quantity) <- paste0("quantity-", seq_len(nrow(table)))
                attr(table, "cohort_note") <- "preserved by the existing codec"
                metadata <- workflowSavedTableMetadata(table, ".rownames", 3L)
                reference <- saveWorkflowTable(store, NULL, "writer-values", identity, ".rownames",
                    metadata = metadata, write_fn = function(path, restoration) {
                        writer <- openWorkflowTableWriter(path, restoration, backend = backend)
                        on.exit(abortWorkflowTableWriter(writer), add = TRUE)
                        appendWorkflowTableWriter(writer, table, integer())
                        appendWorkflowTableWriter(writer, table, 1:32)
                        invisible(gc())
                        appendWorkflowTableWriter(writer, table, 33:64)
                        finishWorkflowTableWriter(writer)
                    })
                actual <- readSavedWorkflowTable(store, reference, names(table), limit = nrow(table))$data
                expect_identical(actual, table)
                writer_path <- file.path(owned, "rejected.parquet")
                writer <- openWorkflowTableWriter(writer_path, metadata, backend = backend)
                appendWorkflowTableWriter(writer, table, integer())
                bad <- table
                bad$Precursor.Quantity <- as.character(bad$Precursor.Quantity)
                expect_error(appendWorkflowTableWriter(writer, bad), "schema differs")
                abortWorkflowTableWriter(writer)
                expect_false(file.exists(writer_path))
                expect_false(file.exists(writer$write_path))
                writer_path <- file.path(owned, "finalizer.parquet")
                writer <- openWorkflowTableWriter(writer_path, metadata, backend = backend)
                appendWorkflowTableWriter(writer, table, integer())
                appendWorkflowTableWriter(writer, table, 1:32)
                temporary <- writer$write_path
                expect_true(file.exists(temporary))
                writer <- NULL
                invisible(gc())
                expect_false(file.exists(temporary))
                expect_false(file.exists(writer_path))
                expect_identical(store$open_resources, 0L)
                closeWorkflowEntityStore(store)
                receipt$lifecycle[[paste0("writer-api-", backend)]] <- list(metadata = TRUE,
                    gc_while_pending = TRUE, schema_rejection = TRUE, finalizer_cleanup = TRUE)
            }
        }
        if (wfx_runtime) {
            info <- workflowFastWriterInfo()
            expect_identical(info$capability$api, "incremental-parquet-v2")
            expect_identical(info$runtime_release_status, "held")
            expect_true(all(c("dia_rsum.c", "duckdb.h", "duckdb_extension.h",
                "build_dia_rsum.R", "configure") %in% names(artifact$manifest$source_sha256)))
            expect_false(any(grepl("_rust$", getNamespaceExports("fastduckdb"))))
            expect_error(testthat::with_mocked_bindings(workflowFastWriterInfo(),
                parquet_writer_capabilities = function() list(api = "wrong"),
                .package = "fastduckdb"), class = "saved_entity_unsupported")
            original_require_namespace <- base::requireNamespace
            expect_error(testthat::with_mocked_bindings(workflowFastWriterInfo(),
                requireNamespace = function(package, ...) {
                    if (identical(package, "fastduckdb")) return(FALSE)
                    original_require_namespace(package, ...)
                }, .package = "base"), class = "saved_entity_unsupported")
            expect_error(testthat::with_mocked_bindings(diaRollupExtensionInfo(),
                packageVersion = function(...) numeric_version("0.0.0"),
                .package = "utils"), class = "dia_rollup_unavailable")
            owned <- file.path(root, "wfx-runtime")
            dir.create(owned)
            data <- data.frame(value = c(NA_real_, NaN, 0, 1, 2, 3),
                label = c("alpha", "β", "你好", NA, "", "end"),
                group = factor(c("a", "b", NA, "b", "a", "a"), levels = c("b", "a", "unused")))
            attr(data, "source") <- "non-rollup codec consumer"
            attr(data$value, "unit") <- "linear"
            metadata <- workflowSavedTableMetadata(data, ".rownames", 3L)
            store <- newWorkflowEntityStore(owned, "wfx-runtime", "DIA")
            reference <- saveWorkflowTable(store, NULL, "codec", identity, ".rownames",
                metadata = metadata, write_fn = function(path, restoration) {
                    writer <- openWorkflowTableWriter(path, restoration, backend = "fast_r")
                    on.exit(abortWorkflowTableWriter(writer))
                    appendWorkflowTableWriter(writer, data, integer())
                    appendWorkflowTableWriter(writer, data, 1:2)
                    invisible(gc())
                    appendWorkflowTableWriter(writer, data, 3:6)
                    finishWorkflowTableWriter(writer)
                })
            actual <- readSavedWorkflowTable(store, reference, names(data), limit = 6L)
            expect_identical(actual$source_row_names, as.character(1:6))
            expect_identical(completeWorkflowTableBatch(actual$data, metadata), data)
            protected <- digest::digest(file = file.path(store$directory,
                reference$payload$relative_path), algo = "sha256")
            invalid <- list(date = as.Date("2026-01-01"), nested = list(1),
                matrix = matrix(1), classed = structure(1, class = "unsupported"))
            for (name in names(invalid)) {
                expect_error(fastduckdb::open_parquet_writer_fast(list(x = invalid[[name]]),
                    file.path(owned, paste0(name, ".parquet"))), "unsupported|prototype|length")
            }
            # R object.size discounts repeated strings; native expansion must not.
            huge <- list(x = rep(strrep("x", 7000L), 10000L))
            writer_path <- file.path(owned, "expanded.parquet")
            handle <- fastduckdb::open_parquet_writer_fast(list(x = character()), writer_path)
            temporary <- fastduckdb::flush_parquet_writer_fast(handle)$temporary_path
            expect_error(fastduckdb::append_parquet_writer_fast(handle, huge), "expanded")
            expect_false(file.exists(temporary))
            rm(huge, handle)
            expect_error(fastduckdb::open_parquet_writer_fast(list(x = double()),
                file.path(owned, "absent", "file.parquet")), "create|writer")
            for (failure in c("flush", "finish")) {
                writer_path <- file.path(owned, paste0(failure, ".parquet"))
                handle <- fastduckdb::open_parquet_writer_fast(list(x = double()), writer_path)
                temporary <- fastduckdb::flush_parquet_writer_fast(handle)$temporary_path
                fastduckdb::append_parquet_writer_fast(handle, list(x = c(1, 2)))
                if (failure == "flush") {
                    expect_error(testthat::with_mocked_bindings(
                        fastduckdb::flush_parquet_writer_fast(handle),
                        flush_incremental_rust = function(...) stop("injected owned flush failure"),
                        .package = "fastduckdb"), "injected owned flush failure")
                } else {
                    dir.create(writer_path)
                    expect_error(fastduckdb::close_parquet_writer_fast(handle), "move|exists|finish|publish")
                    expect_true(dir.exists(writer_path))
                    unlink(writer_path, recursive = TRUE)
                }
                expect_false(file.exists(temporary))
                expect_true(fastduckdb::abort_parquet_writer_fast(handle)$closed)
            }
            expect_identical(digest::digest(file = file.path(store$directory,
                reference$payload$relative_path), algo = "sha256"), protected)
            expect_identical(store$open_resources, 0L)
            closeWorkflowEntityStore(store)
            # Publication must never replace a destination acquired by another producer.
            for (competitor in c("sentinel", "writer")) {
                target <- file.path(owned, paste0("collision-", competitor, ".parquet"))
                handle <- fastduckdb::open_parquet_writer_fast(list(x = double()), target)
                temporary <- fastduckdb::flush_parquet_writer_fast(handle)$temporary_path
                fastduckdb::append_parquet_writer_fast(handle, list(x = 1))
                if (competitor == "sentinel") writeLines("owned by another producer", target) else
                    callr::r(function(target) {
                        second <- fastduckdb::open_parquet_writer_fast(list(x = double()), target)
                        fastduckdb::append_parquet_writer_fast(second, list(x = 2))
                        fastduckdb::close_parquet_writer_fast(second)
                        TRUE
                    }, args = list(target))
                before <- digest::digest(file = target, algo = "sha256")
                expect_error(fastduckdb::close_parquet_writer_fast(handle), "publish")
                expect_identical(digest::digest(file = target, algo = "sha256"), before)
                expect_false(file.exists(temporary))
            }
            target <- file.path(owned, "lazy-collision.parquet")
            writer <- openWorkflowTableWriter(target, metadata, backend = "fast_r")
            writeLines("unowned destination", target)
            before <- digest::digest(file = target, algo = "sha256")
            expect_error(appendWorkflowTableWriter(writer, data), "destination exists")
            abortWorkflowTableWriter(writer)
            expect_identical(digest::digest(file = target, algo = "sha256"), before)
            # The legacy public writer retains its intentional replacement behavior.
            target <- file.path(owned, "legacy-replace.parquet")
            fastduckdb::write_parquet_fast(data.frame(x = 1), target, n_threads = 8L)
            fastduckdb::write_parquet_fast(data.frame(x = 2), target, n_threads = 8L)
            expect_equal(arrow::read_parquet(target)$x, 2)
            expect_identical(info$capability$input_modes, "r_columns")
            expect_error(workflowFastWriterInfo("fast_arrow"), class = "saved_entity_unsupported")
            for (stream_data in list(data.frame(x = 1), data.frame(x = as.Date("2026-01-01")),
                data.frame(x = I(list(1:2))))) {
                stream_table <- arrow::Table$create(stream_data)
                target <- file.path(owned, "unsupported-stream.parquet")
                expect_error(fastduckdb::open_parquet_writer_fast(
                    nanoarrow::as_nanoarrow_array_stream(stream_table$Slice(0, 0)), target),
                    "outside the qualified")
                expect_false(file.exists(target))
                handle <- fastduckdb::open_parquet_writer_fast(list(x = double()), target)
                temporary <- fastduckdb::flush_parquet_writer_fast(handle)$temporary_path
                expect_error(fastduckdb::append_parquet_writer_fast(handle,
                    nanoarrow::as_nanoarrow_array_stream(stream_table), owner = stream_table,
                    logical_bytes = 0), "outside the qualified")
                expect_false(file.exists(temporary))
            }
            # Bind installed R bytes as well as the unchanged DLL/capability.
            isolated <- file.path(owned, "altered-r-library")
            dir.create(isolated)
            expect_true(file.copy(system.file(package = "fastduckdb"), isolated, recursive = TRUE))
            mismatch <- callr::r(function(isolated) {
                original_path <- base::system.file
                target <- file.path(isolated, "fastduckdb", "R", "fastduckdb.rdb")
                con <- file(target, "ab")
                writeBin(as.raw(0), con)
                close(con)
                testthat::with_mocked_bindings(
                    tryCatch({ MultiScholaR:::workflowFastWriterInfo(); FALSE },
                        saved_entity_unsupported = function(error) TRUE),
                    system.file = function(..., package = "base", lib.loc = NULL, mustWork = FALSE) {
                        if (identical(package, "fastduckdb")) return(file.path(isolated, "fastduckdb", ...))
                        original_path(..., package = package, lib.loc = lib.loc, mustWork = mustWork)
                    }, .package = "base")
            }, args = list(isolated))
            expect_true(mismatch)
            receipt$runtime <- list(runtime_release_status = "held", capability = info,
                constructor_append_flush_finish_failures = TRUE, original_payload_unchanged = TRUE)
        }
        for (arm in if (materialize_only || adaptation_only || fast_only || writer_only) character() else c("C", "D")) {
            owned <- file.path(root, paste0("science-", arm))
            dir.create(owned)
            store <- newWorkflowEntityStore(owned, paste0("relational-", arm), "DIA")
            on.exit(closeWorkflowEntityStore(store), add = TRUE)
            identity <- list(input = "fixture", design = "paired", parameters = "exact-native",
                selectors = "all", software = "extension-v1")
            for (case in names(cases)) {
                object <- base
                object@peptide_data <- cases[[case]]
                expected <- invoke(object@peptide_data)
                input <- saveWorkflowPrecursorRollupInput(store, object, paste0("source-", case), identity)
                result <- computeDIARelationalRollup(store, input,
                    list(prototype = expected$data[0, , drop = FALSE], identity = identity), arm, NULL, phase, extension)
                expect_identical(result$summary, expected$summary, info = paste(arm, case))
                actual <- readSavedWorkflowTable(store, result$reference, columns = names(expected$data),
                    offset = 0L, limit = max(1L, nrow(expected$data)), reader = "duckdb")$data
                expect_identical(actual, expected$data, info = paste(arm, case))
                expect_identical(store$open_resources, 0L)
                receipt$cases[[paste(arm, case)]] <- list(rows = nrow(actual),
                    exact = identical(actual, expected$data) && identical(result$summary, expected$summary))
            }
            invalid_cases <- list(negative = -1, infinite = Inf, nan = NaN, integer = 1L)
            for (case in names(invalid_cases)) {
                object <- base
                if (case == "integer") object@peptide_data$Precursor.Quantity <-
                    as.integer(object@peptide_data$Precursor.Quantity) else
                    object@peptide_data$Precursor.Quantity[[1L]] <- invalid_cases[[case]]
                input <- saveWorkflowPrecursorRollupInput(store, object, paste0("invalid-", case), identity)
                expect_error(computeDIARelationalRollup(store, input,
                    list(prototype = invoke(base@peptide_data)$data[0, ], identity = identity), arm, NULL,
                    function(label, fn) fn(), extension), "negative or non-finite|double-precision")
                expect_identical(store$open_resources, 0L)
            }
            for (case in c("duplicate", "evidence")) {
                object <- base
                object@peptide_data <- if (case == "duplicate") rbind(table, table[1L, ]) else table
                if (case == "evidence") object@peptide_data$identification_peptide_count[[1L]] <- 4L
                input <- saveWorkflowPrecursorRollupInput(store, object, paste0("invalid-", case), identity)
                expect_error(computeDIARelationalRollup(store, input,
                    list(prototype = invoke(table)$data[0, ], identity = identity), arm, NULL,
                    function(label, fn) fn(), extension), "duplicate precursor|conflicting frozen")
                expect_identical(store$open_resources, 0L)
            }
            logged <- base
            logged@is_logged_data <- TRUE
            expect_error(saveWorkflowPrecursorRollupInput(store, logged, "logged", identity), "linear peptide")
            input <- saveWorkflowPrecursorRollupInput(store, base, "cancel-source", identity)
            expect_error(local({
                scope <- newDIARelationalScope(store, input, extension)
                on.exit(releaseWorkflowTableReadScope(scope))
                streamDIARelationalQuery(scope$connection, "SELECT 1 AS value",
                    function(...) stop("injected relational cancellation"))
            }), "injected relational cancellation")
            expect_identical(store$open_resources, 0L)
            closeWorkflowEntityStore(store)
            lifecycle_root <- file.path(root, paste0("lifecycle-", arm))
            dir.create(lifecycle_root)
            manager <- newDIAStateAdapter(newWorkflowEntityStore(lifecycle_root, paste0("lifecycle-", arm), "DIA"),
                WorkflowState$new(audit_enabled = FALSE))
            on.exit(closeDIAStateAdapter(manager), add = TRUE)
            manager$.context$reference_states <- TRUE
            manager$saveState("protein_peptide_filtered", base, list(), "relational fixture")
            value <- runDIABulkOperation(manager, lifecycle_root, proof, arm, NULL, phase, failure = TRUE)
            closeDIAStateAdapter(manager)
            receipt$lifecycle[[arm]] <- list(summary = value$summary, figures = length(value$figures))
        }
        profiles <- if (writer_r) if (wfx_runtime) "P14" else c("P14", "P15", "P16") else if (writer_only) c("P13", "P14") else if (fast_only) paste0("P", 8:12) else if (adaptation_only) c("P6", "P7") else
            if (materialize_only) "P5" else paste0("P", 0:5)
        for (profile in profiles) {
            owned <- file.path(root, paste0("package-", profile))
            dir.create(owned)
            store <- newWorkflowEntityStore(owned, profile, "DIA")
            on.exit(closeWorkflowEntityStore(store), add = TRUE)
            for (case in names(cases)) {
                object <- base
                object@peptide_data <- cases[[case]]
                expected <- invoke(object@peptide_data)
                input <- saveWorkflowPrecursorRollupInput(store, object, case, identity)
                result <- withr::with_options(list(multischolar.dia_rollup_profile = profile),
                    computeSavedWorkflowPrecursorRollup(store, input,
                        list(prototype = expected$data[0, , drop = FALSE], identity = identity)))
                actual <- readSavedWorkflowTable(store, result$reference, names(expected$data),
                    limit = max(1L, nrow(expected$data)))$data
                expect_identical(actual, expected$data, info = paste(profile, case))
                expect_identical(result$summary, expected$summary, info = paste(profile, case))
                expect_identical(result$execution$engine_threads, 8L)
                if (writer_only) {
                    expect_identical(result$execution$writer_backend,
                        if (profile == "P13") "fast_arrow" else "fast_r")
                    expect_identical(result$execution$writer_threads, 8L)
                    expect_equal(result$execution$intermediate_bytes_written, 0)
                }
                if (profile %in% c("P10", "P12")) {
                    expect_identical(result$execution$fastduckdb_threads, 8L)
                    if (nrow(cases[[case]])) expect_gt(result$execution$fastduckdb_calls, 0L)
                }
                expect_identical(store$open_resources, 0L)
                if (fast_only) expect_length(list.files(store$directory,
                    pattern = "^(sample-aligned|sample-cohort-[0-9]+|ordered-result-[0-9]+)[.]parquet$",
                    recursive = TRUE), 0L)
                receipt$cases[[paste(profile, case)]] <- list(rows = nrow(actual),
                    exact = identical(actual, expected$data) && identical(result$summary, expected$summary))
            }
            expect_error(local({
                scope <- newDIARelationalScope(store, input, profile = profile)
                on.exit(releaseWorkflowTableReadScope(scope))
                streamDIARelationalQuery(scope$connection, "SELECT 1 AS value",
                    function(...) stop("injected profile cancellation"), streaming = profile == "P2")
            }), "injected profile cancellation")
            expect_identical(store$open_resources, 0L)
            if (fast_only || writer_only) {
                input <- saveWorkflowPrecursorRollupInput(store, base, "staging-failure", identity)
                original_append <- appendWorkflowTableWriter
                expect_error(testthat::with_mocked_bindings(
                    computeDIARelationalRollup(store, input,
                        list(prototype = invoke(base@peptide_data)$data[0, ], identity = identity), profile = profile),
                    appendWorkflowTableWriter = function(writer, data, ...) {
                        if (nrow(data)) stop("injected staged output failure")
                        original_append(writer, data, ...)
                    }), "injected staged output failure")
                expect_identical(store$open_resources, 0L)
                expect_length(list.files(store$directory,
                    pattern = "^(sample-aligned|sample-cohort-[0-9]+|ordered-result-[0-9]+)[.]parquet$|[.]tmp-[0-9]+$",
                    recursive = TRUE), 0L)
            }
            closeWorkflowEntityStore(store)
        }
        for (profile in if (writer_r) if (wfx_runtime) "P14" else c("P14", "P15", "P16") else if (writer_only) c("P13", "P14") else if (fast_only) paste0("P", 8:12) else if (adaptation_only) c("P6", "P7") else
            if (materialize_only) "P5" else character()) {
            lifecycle_root <- file.path(root, paste0("lifecycle-", profile))
            dir.create(lifecycle_root)
            manager <- newDIAStateAdapter(newWorkflowEntityStore(lifecycle_root, paste0("lifecycle-", profile), "DIA"),
                WorkflowState$new(audit_enabled = FALSE))
            on.exit(closeDIAStateAdapter(manager), add = TRUE)
            manager$.context$reference_states <- TRUE
            manager$saveState("protein_peptide_filtered", base, list(), "materialization fixture")
            value <- runDIABulkOperation(manager, lifecycle_root, proof, profile, NULL, phase, failure = TRUE)
            closeDIAStateAdapter(manager)
            receipt$lifecycle[[profile]] <- list(summary = value$summary, figures = length(value$figures))
        }
        expect_error(diaRollupProfile("unqualified"), class = "dia_rollup_unavailable")
        expect_null(getOption("multischolar.dia_rollup_profile"))
        stopifnot(all(vapply(receipt$cases, function(x) x$exact, logical(1))))
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        flush()
        stop(error)
    })
    flush()
})

# WFX-002 readback-only verification: never repeat accepted computation.
test_that("the packaged writer supports predecessor readback and pending-output rollback", {
    completed <- Sys.getenv("MULTISCHOLAR_WFX002_VERIFY_MEMBER", unset = "")
    skip_if(!nzchar(completed), "explicit bounded installed readback only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    stopifnot(startsWith(normalizePath(completed), paste0(normalizePath(proof), "/")))
    # Preserve the wrongly located original receipt, and bind a recovery record.
    receipt_path <- file.path(completed, "relational-preflight.json")
    if (!file.exists(receipt_path)) receipt_path <- file.path(completed, "wfx-runtime", "finish.parquet")
    receipt <- jsonlite::read_json(receipt_path)
    expect_identical(receipt$status, "complete")
    expect_true(all(vapply(receipt$cases, function(case) isTRUE(case$exact), logical(1))))
    results <- readRDS(file.path(completed, "reader-tests.rds"))
    expect_true(all(results$failed == 0L & !results$error))
    source_store <- file.path(completed, "lifecycle-P14", ".workflow-entities")
    original_head <- digest::digest(file = file.path(source_store, "dia-state-head.rds"), algo = "sha256")
    clone <- file.path(root, "rollback-owned")
    dir.create(clone)
    expect_true(file.copy(source_store, clone, recursive = TRUE))
    owner <- readSavedStoreIdentity(file.path(clone, ".workflow-entities", "store.json"))
    store <- newWorkflowEntityStore(clone, owner$project_id, owner$workflow_id)
    manager <- newDIAStateAdapter(store, reopen = TRUE)
    record <- readDIAStateMetadata(manager)
    reference <- record$data$components$peptide_data
    scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
    expect_lte(scope$metadata$rows, 10000L)
    actual <- readSavedWorkflowTable(store, reference, scope$metadata$columns, limit = 10000L,
        scope = scope)$data
    expected_hash <- digest::digest(actual, algo = "sha256")
    releaseWorkflowTableReadScope(scope)
    history <- manager$state_history
    output <- file.path(clone, "pending-output")
    dir.create(output)
    manager$.context$experiment_paths <- list(source_dir = output)
    expect_error(publishDIAActionOutputs(manager, getDIAStateReference(manager),
        "wfx-rollback", list(output = output), function(...) stop("owned pending output")),
        "owned pending output")
    expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
    closeDIAStateAdapter(manager)
    baseline <- file.path(dirname(proof), "browser-smoke.XsaPvT", "library")
    rollback <- callr::r(function(library, clone, owner) {
        .libPaths(c(library, .libPaths()))
        ns <- asNamespace("MultiScholaR")
        api <- function(name) get(name, ns)
        stopifnot(identical(normalizePath(getNamespaceInfo(ns, "path")),
            normalizePath(file.path(library, "MultiScholaR"))))
        store <- api("newWorkflowEntityStore")(clone, owner$project_id, owner$workflow_id)
        manager <- api("newDIAStateAdapter")(store, reopen = TRUE)
        on.exit(api("closeDIAStateAdapter")(manager))
        record <- api("readDIAStateMetadata")(manager)
        reference <- record$data$components$peptide_data
        scope <- api("newWorkflowTableReadScope")(store, reference, reader = "duckdb")
        on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE, after = FALSE)
        stopifnot(scope$metadata$rows <= 10000L)
        data <- api("readSavedWorkflowTable")(store, reference, scope$metadata$columns,
            limit = 10000L, scope = scope)$data
        pending <- tryCatch({ api("assertDIAExportReady")(manager); FALSE },
            saved_entity_outputs_pending = function(error) TRUE)
        list(hash = digest::digest(data, algo = "sha256"), history = manager$state_history,
            pending_recognized = pending)
    }, args = list(baseline, clone, owner))
    expect_identical(rollback$hash, expected_hash)
    expect_identical(rollback$history, history)
    expect_true(rollback$pending_recognized)
    expect_identical(digest::digest(file = file.path(source_store, "dia-state-head.rds"),
        algo = "sha256"), original_head)
    jsonlite::write_json(list(status = "complete", computation_repeated = FALSE,
        recovered_receipt = receipt_path, receipt_sha256 = digest::digest(file = receipt_path, algo = "sha256"),
        original_head_sha256 = original_head, rollback = rollback),
        file.path(root, "wfx-runtime-readback.json"), auto_unbox = TRUE, pretty = TRUE)
})
