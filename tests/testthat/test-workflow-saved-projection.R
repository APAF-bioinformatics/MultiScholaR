testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"),
    "durable saved-table backend is Linux-qualified only")

readerTrialTestEnvironment <- function() {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$recordedCheckpointApi <- function(library_path) {
        namespace <- asNamespace("MultiScholaR")
        \(name) get(name, envir = namespace, inherits = FALSE)
    }
    environment
}

test_that("installed adapter members preserve exact consumers and close on late failure", {
    skip_if_not_installed("fastduckdb")
    driver <- readerTrialTestEnvironment()
    root <- withr::local_tempdir()
    api <- driver$recordedCheckpointApi("unused")
    spec <- driver$prepareCheckpointReaderTable(driver$checkpointProjectionTable(300L, 4L),
        "small", root, api, adapter = TRUE)
    expect_true(spec$store_adapter)
    jsonlite::write_json(list(status = "prepared", schema = "dia013_store_adapter_v1",
        cases = list(small = spec)), file.path(root, "reader-inputs.json"), auto_unbox = TRUE)
    for (engine in c("adapter_rds", "adapter_duckdb", "adapter_fastduckdb")) {
        member <- withr::local_tempdir()
        result <- driver$executeCheckpointReader("unused", member, root, paste("small", engine, "1", sep = ":"))
        expect_identical(result$status, "complete")
        expect_true(result$store_adapter)
        expect_named(result$outputs, "projected")
        expect_equal(result$open_resources, 0)
        if (engine != "adapter_rds") {
            resources <- result$events[[1L]]$value$resources$store_metrics
            expect_identical(resources$open_resources, 0L)
            expect_identical(resources$table_reads, 1L)
            expect_identical(resources$table_rows, 300L)
        }
        expect_identical(driver$verifyCheckpointReader("unused", withr::local_tempdir(),
            member)$status, "passed_reader_exactness")
    }
    owner <- NULL
    original <- driver$openCheckpointReader
    driver$openCheckpointReader <- function(...) { owner <<- original(...); owner }
    driver$consumeCheckpointReader <- function(...) stop("adapter consumer failed")
    expect_error(driver$executeCheckpointReader("unused", withr::local_tempdir(), root,
        "small:adapter_fastduckdb:1"), "adapter consumer failed")
    expect_true(owner$closed)
    expect_null(owner$table_scope)
    expect_null(owner$store)
    expect_identical(owner$read_metrics$open_resources, 0L)
})


test_that("reader trial preserves exact values and closes all four small readers", {
    skip_if_not_installed("fastduckdb")
    driver <- readerTrialTestEnvironment()
    expect_no_warning(result <- driver$checkCheckpointReader("unused", withr::local_tempdir()))
    expect_identical(result$status, "passed_reader_checks")
    expect_named(result$engines, c("rds", "arrow", "duckdb", "fastduckdb"))
    expect_true(all(vapply(result$engines, \(engine) engine$exact && engine$closed, logical(1))))
})

test_that("reader trial rejects changed and foreign payloads without retaining resources", {
    skip_if_not_installed("fastduckdb")
    driver <- readerTrialTestEnvironment()
    api <- driver$recordedCheckpointApi("unused")
    root <- withr::local_tempdir()
    table <- driver$checkpointProjectionTable(100L, 2L)
    spec <- driver$prepareCheckpointReaderTable(table, "test", root, api)
    expect_error(driver$checkpointReaderPaths(spec, withr::local_tempdir(), "rds"))
    for (engine in c("rds", "arrow", "duckdb", "fastduckdb")) {
        scope <- driver$openCheckpointReader(spec, root, engine, api)
        expect_error(driver$readCheckpointReader(scope, "missing"), "unsupported")
        expect_error(driver$readCheckpointReader(scope, "logFC", c(3, 1)), "unsupported")
        expect_error(driver$readCheckpointReader(scope, "logFC", c("PG_01", "PG_02")), "unsupported")
        expect_error(driver$readCheckpointReader(scope, c("logFC", "logFC")), "unsupported")
        path <- scope$paths[["payload"]]
        bytes <- readBin(path, "raw", n = file.size(path))
        writeBin(as.raw(1), path)
        expect_error(driver$readCheckpointReader(scope, "logFC"), "changed")
        driver$closeCheckpointReader(scope)
        expect_true(scope$closed)
        expect_null(scope$connection)
        expect_null(scope$table)
        expect_error(driver$checkpointReaderPaths(spec, root, engine), "integrity")
        writeBin(bytes, path)
        driver$closeCheckpointReader(scope)
        expect_error(driver$readCheckpointReader(scope, "logFC"), "closed")
    }
})

test_that("reader members verify complete small outputs and close after consumer failure", {
    skip_if_not_installed("fastduckdb")
    driver <- readerTrialTestEnvironment()
    root <- withr::local_tempdir()
    spec <- driver$prepareCheckpointReaderTable(driver$checkpointProjectionTable(300L, 4L),
        "small", root, driver$recordedCheckpointApi("unused"))
    jsonlite::write_json(list(status = "prepared", cases = list(small = spec)),
        file.path(root, "reader-inputs.json"), auto_unbox = TRUE)
    for (engine in c("rds", "fastduckdb")) {
        member <- withr::local_tempdir()
        result <- driver$executeCheckpointReader("unused", member, root,
            paste("small", engine, "1", sep = ":"))
        expect_identical(result$status, "complete")
        expect_equal(result$open_resources, 0)
        verified <- driver$verifyCheckpointReader("unused", withr::local_tempdir(), member)
        expect_identical(verified$status, "passed_reader_exactness")
    }
    owner <- NULL
    original_open <- driver$openCheckpointReader
    driver$openCheckpointReader <- function(...) {
        owner <<- original_open(...)
        owner
    }
    driver$consumeCheckpointReader <- function(...) stop("late consumer failure")
    failed <- withr::local_tempdir()
    expect_error(driver$executeCheckpointReader("unused", failed, root, "small:fastduckdb:1"),
        "late consumer failure")
    expect_true(owner$closed)
    expect_null(owner$table)
    expect_false(file.exists(file.path(failed, "reader-member.json")))
})

makeSavedProjectionDaTable <- function() {
    set.seed(13013L)
    samples <- sprintf("S%02d", seq_len(8L))
    features <- sprintf("P%03d", seq_len(300L))
    design <- data.frame(group = factor(rep(c("control", "treated"), each = 4L)),
        row.names = samples)
    data <- matrix(stats::rnorm(length(features) * length(samples)),
        nrow = length(features), dimnames = list(features, samples))
    data[seq.int(10L, 300L, by = 10L), 5:8] <-
        data[seq.int(10L, 300L, by = 10L), 5:8] + 2
    result <- suppressMessages(runTestsContrasts(data,
        "treated_vs_control=grouptreated-groupcontrol", design, "~0+group",
        weights = NA, treat_lfc_cutoff = NA))$results[[1L]]
    result$Protein.Group <- rownames(result)
    result$comparison <- factor("treated_vs_control",
        levels = c("treated_vs_control", "unused"))
    result$direction <- ordered(ifelse(result$logFC >= 0, "up", "down"),
        levels = c("down", "up", "unchanged"))
    result$edge <- seq_len(nrow(result)) / 10
    result$edge[1:4] <- c(NA_real_, NaN, Inf, -Inf)
    attr(result$logFC, "units") <- "log2"
    attr(result, "provenance") <- "real runTestsContrasts fixture"
    result
}

projectionIdentity <- function(table, selector = "da-result-table") {
    newSavedWorkflowScientificIdentity(
        digest::digest(table, algo = "sha256"),
        "design-projection-v1",
        "qvalue-bh-v1",
        selector,
        "MultiScholaR-projection-test-v1")
}

test_that("version 2 table metadata keeps row-sized names in Parquet for both readers", {
    skip_if_not_installed("fastduckdb")
    table <- data.frame(id = paste0("P", 1:40),
        value = rep(c(NA_real_, NaN, Inf, -Inf, 0, 1.5, 3, 8), 5),
        group = ordered(rep(c("b", "a"), 20), levels = c("a", "unused", "b")),
        row.names = paste0("row_", 40:1))
    attr(table$value, "units") <- "unchanged units"
    attr(table, "provenance") <- "adapter exactness fixture"
    for (container in c("data.frame", "tibble")) {
        original <- if (container == "tibble") tibble::as_tibble(table) else table
        plain <- unclass(original)
        names(plain$id) <- paste0("element_", seq_len(nrow(original)))
        names(plain$value) <- rep(c("", NA_character_, "named", "duplicate"), 10)
        class(plain) <- class(original)
        original <- plain
        store <- newWorkflowEntityStore(withr::local_tempdir(), paste0("v2-", container), "DIA")
        ref <- saveWorkflowTable(store, original, "named-table", projectionIdentity(original),
            ".rownames", metadata_version = 2L)
        metadata <- readSavedWorkflowTableMetadata(file.path(store$directory, ref$payload$relative_path))
        expect_identical(metadata$schema, 2L)
        expect_null(metadata$row_names)
        expect_named(metadata$name_columns, c("id", "value"))
        expect_lt(as.numeric(utils::object.size(metadata)), 16000)
        for (reader in c("duckdb", "fastduckdb")) {
            scope <- newWorkflowTableReadScope(store, ref, reader)
            if (reader == "fastduckdb") expect_null(scope$connection)
            for (offset in c(0L, 7L, 39L, 40L)) {
                indices <- seq_len(nrow(original))
                indices <- indices[indices > offset & indices <= offset + 8L]
                expected <- original[indices, , drop = FALSE]
                expected <- unclass(expected)
                names(expected$id) <- names(original$id)[indices]
                names(expected$value) <- names(original$value)[indices]
                attr(expected$value, "units") <- attr(original$value, "units")
                attr(expected, "provenance") <- attr(original, "provenance")
                class(expected) <- class(original)
                actual <- hydrateSavedWorkflowProjection(store, ref, names(original),
                    list(offset = offset, limit = 8L), scope = scope)
                expect_identical(actual$data, expected)
                expect_identical(actual$source_row_names, rownames(original)[indices])
                expect_lte(actual$collected_bytes, 67108864)
            }
            releaseWorkflowTableReadScope(scope)
            expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
        }
        reopened <- newWorkflowEntityStore(dirname(store$directory), paste0("v2-", container), "DIA")
        expect_identical(recoverSavedWorkflowRequest(reopened, "named-table"), ref)
        closeWorkflowEntityStore(reopened)
        closeWorkflowEntityStore(store)
    }
})

test_that("fast adapter reads legacy payloads without DuckDB and validates once per action", {
    skip_if_not_installed("fastduckdb")
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "fast-scope", "DIA")
    ref <- saveWorkflowTable(store, table, "legacy", projectionIdentity(table), ".rownames")
    capabilities <- savedWorkflowEntityCapabilities(ref, reader = "fastduckdb")
    expect_true(capabilities$row_ranges)
    expect_false(capabilities$predicate_pushdown)
    expect_false(capabilities$key_lookup)
    expect_false(capabilities$ranked_selection)
    expect_false(capabilities$full_hydration)
    count <- 0L
    metadata_reads <- 0L
    validate <- validateSavedWorkflowReference
    read_metadata <- readSavedWorkflowTableMetadata
    testthat::with_mocked_bindings({
        scope <- newWorkflowTableReadScope(store, ref, "fastduckdb")
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        expect_null(scope$connection)
        for (offset in c(0L, 20L, 60L)) {
            actual <- readSavedWorkflowTable(store, ref, c("Protein.Group", "edge"),
                offset = offset, limit = 10L, scope = scope)
            expect_identical(actual$data$edge, table$edge[offset + seq_len(10L)])
        }
        expect_identical(count, 1L)
        expect_identical(metadata_reads, 1L)
        expect_error(readSavedWorkflowTable(store, ref, "edge", keys = "P1", scope = scope),
            class = "saved_entity_unsupported")
        expect_error(readSavedWorkflowDiaHeatmapTable(store, ref, "edge", 2L, scope = scope),
            class = "saved_entity_unsupported")
        expect_error(readSavedWorkflowTable(store, ref, "edge", limit = 10001L, scope = scope),
            class = "saved_entity_limit")
        expect_error(readSavedWorkflowTable(store, ref, "edge", max_bytes = 1, scope = scope),
            class = "saved_entity_limit")
        releaseWorkflowTableReadScope(scope)
        expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    }, validateSavedWorkflowReference = function(...) { count <<- count + 1L; validate(...) },
        readSavedWorkflowTableMetadata = function(...) {
            metadata_reads <<- metadata_reads + 1L
            read_metadata(...)
        }, .package = "MultiScholaR")
    scope <- newWorkflowTableReadScope(store, ref, "fastduckdb")
    path <- file.path(store$directory, ref$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.size(path))
    writeBin(as.raw(1), path)
    expect_error(readSavedWorkflowTable(store, ref, "edge", scope = scope),
        class = "saved_entity_corrupt")
    expect_error(releaseWorkflowTableReadScope(scope), class = "saved_entity_corrupt")
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    writeBin(bytes, path)
    expect_error(readSavedWorkflowTable(store, ref, "bad", reader = "fastduckdb"),
        class = "saved_entity_limit")
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    closeWorkflowEntityStore(store)
})

test_that("version 2 metadata stays bounded without relaxing the legacy ceiling", {
    table <- data.frame(id = seq_len(1000L), value = seq_len(1000L) / 10,
        row.names = paste0(strrep("r", 9000), seq_len(1000L)))
    expect_error(workflowSavedTableMetadata(table, ".rownames"),
        class = "saved_entity_limit")
    metadata <- workflowSavedTableMetadata(table, ".rownames", 2L)
    expect_null(metadata$row_names)
    expect_identical(metadata$row_name_column, ".multischolar_row_key")
    expect_lt(as.numeric(utils::object.size(metadata)), 16000)
    expect_error(workflowSavedTableMetadata(table, ".rownames", 4L),
        class = "saved_entity_unsupported")
    class(table$value) <- "integer64"
    expect_error(workflowSavedTableMetadata(table, ".rownames", 2L),
        class = "saved_entity_unsupported")
})

test_that("Parquet DuckDB pages preserve order types attributes and nonfinite values", {
    table <- makeSavedProjectionDaTable()
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "projection-project", "DIA")
    reference <- saveWorkflowTable(store, table, "da-table",
        projectionIdentity(table), ".rownames",
        factors = list(route = "iq", contrast = "treated_vs_control"))
    expect_identical(reference$payload$codec, "parquet_duckdb")
    capabilities <- savedWorkflowEntityCapabilities(reference)
    expect_true(capabilities$column_projection)
    expect_true(capabilities$predicate_pushdown)
    expect_false(capabilities$full_hydration)
    columns <- c("Protein.Group", "logFC", "fdr_qvalue", "comparison",
        "direction", "edge")
    page <- readSavedWorkflowTable(store, reference, columns,
        offset = 10L, limit = 30L)
    expected <- table[11:40, columns, drop = FALSE]
    attr(expected$logFC, "units") <- attr(table$logFC, "units")
    attr(expected, "provenance") <- attr(table, "provenance")
    expect_identical(page$data, expected)
    expect_identical(page$source_row_names, rownames(expected))
    expect_identical(page$source_order, 11:40)
    expect_true(is.factor(page$data$comparison))
    expect_identical(levels(page$data$comparison),
        c("treated_vs_control", "unused"))
    expect_true(is.ordered(page$data$direction))
    expect_identical(levels(page$data$direction), c("down", "up", "unchanged"))
    expect_identical(attr(page$data$logFC, "units"), "log2")
    expect_identical(attr(page$data, "provenance"),
        "real runTestsContrasts fixture")
    expect_true(page$projected)
    info <- workflowEntityCacheInfo(store)
    expect_identical(info$table_reads, 1L)
    expect_identical(info$table_rows, 30L)
    expect_identical(info$table_columns, length(columns))
    expect_identical(info$open_resources, 0L)
    expect_identical(info$payload_reads, 0L)
})

test_that("ordered duplicate selectors and missing SQL-like keys are safe", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(),
        "selector-project", "DIA")
    reference <- saveWorkflowTable(store, table, "selector-table",
        projectionIdentity(table), ".rownames")
    ids <- rownames(table)[c(5L, 2L, 5L)]
    requested <- c(ids, "P001' OR 1=1 --")
    result <- hydrateSavedWorkflowProjection(store, reference,
        c("Protein.Group", "logFC", "edge"), filters = list(keys = requested,
            max_rows = 10L, limit = 10L))
    expect_identical(result$requested_keys, requested)
    expect_identical(result$found, c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(result$source_row_names[1:3], ids)
    expect_true(is.na(result$source_row_names[[4L]]))
    expect_identical(result$data$logFC[1:3], table[ids, "logFC"])
    expect_identical(result$data$logFC[[1L]], result$data$logFC[[3L]])
    expect_true(is.na(result$data$logFC[[4L]]))
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    expect_error(hydrateSavedWorkflowProjection(store, reference,
        "logFC; DROP TABLE x", filters = list()), class = "saved_entity_limit")
    expect_error(hydrateSavedWorkflowProjection(store, reference, "logFC",
        filters = list(sql = "1=1")), class = "saved_entity_invalid")
    expect_error(readSavedWorkflowTable(store, reference, "logFC",
        limit = 10001L, max_rows = 10000L), class = "saved_entity_limit")
})

test_that("edge nonfinite values round-trip and corrupt or foreign tables fail", {
    table <- makeSavedProjectionDaTable()
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "integrity-project", "DIA")
    reference <- saveWorkflowTable(store, table, "edge-table",
        projectionIdentity(table), ".rownames")
    edge <- readSavedWorkflowTable(store, reference, "edge", offset = 0L, limit = 4L)
    expect_identical(edge$data$edge, c(NA_real_, NaN, Inf, -Inf))
    reopened <- newWorkflowEntityStore(root, "integrity-project", "DIA")
    recovered <- recoverSavedWorkflowRequest(reopened, "edge-table")
    expect_identical(recovered, reference)
    expect_identical(readSavedWorkflowTable(reopened, recovered, "edge",
        offset = 0L, limit = 4L)$data$edge, c(NA_real_, NaN, Inf, -Inf))

    foreign <- newWorkflowEntityStore(withr::local_tempdir(), "foreign", "DIA")
    expect_error(readSavedWorkflowTable(foreign, reference, "edge"),
        class = "saved_entity_scope")
    path <- file.path(store$directory, reference$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    changed <- bytes
    changed[[length(changed)]] <- as.raw(bitwXor(as.integer(changed[[length(changed)]]), 1L))
    writeBin(changed, path)
    expect_error(readSavedWorkflowTable(store, reference, "edge"),
        class = "saved_entity_corrupt")
    writeBin(bytes, path)
})

test_that("whole RDS remains full hydration and Limpa abundance projection is rejected", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "comparators", "DIA")
    envelope <- newSavedWorkflowResultTable(table, ".rownames",
        "treated_vs_control", "fdr_qvalue")
    rds_reference <- saveWorkflowConsumerEntity(store, envelope, "rds-table",
        projectionIdentity(table, "rds"), "result_table")
    expect_error(hydrateSavedWorkflowProjection(store, rds_reference,
        c("logFC", "fdr_qvalue")), class = "saved_entity_unsupported")
    expect_true(savedWorkflowEntityCapabilities(rds_reference)$full_hydration)

    parquet_reference <- saveWorkflowTable(store, table, "parquet-table",
        projectionIdentity(table, "parquet"), ".rownames")
    expect_error(readSavedWorkflowLimpaModelInput(store, parquet_reference),
        "complete aligned EList", class = "saved_entity_unsupported")
    expression <- matrix(seq_len(12), 3, 4,
        dimnames = list(paste0("P", 1:3), paste0("S", 1:4)))
    support <- expression / 10
    observations <- expression * 0 + 3L
    original <- methods::new("EList", list(E = expression,
        genes = data.frame(id = rownames(expression), row.names = rownames(expression)),
        other = list(standard.error = support, n.observations = observations)))
    changed_support <- original
    changed_support$other$standard.error[[1L, 1L]] <- 99
    changed_order <- original
    changed_order$other <- lapply(changed_order$other,
        \(value) value[rev(seq_len(nrow(value))), , drop = FALSE])
    design <- data.frame(group = rep(c("A", "B"), each = 2L),
        row.names = colnames(expression))
    ids <- lapply(list(original, changed_support, changed_order),
        \(elist) newWorkflowLimpaFitIdentity(elist, design, "~0+group", "limpa-v1"))
    expect_identical(original$E, changed_support$E)
    expect_identical(dim(original$other$standard.error),
        dim(changed_order$other$standard.error))
    expect_length(unique(vapply(ids, digest::digest, character(1),
        algo = "sha256")), 3L)
})

test_that("one table action reuses validation metadata and connection", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "scoped-project", "DIA")
    reference <- saveWorkflowTable(store, table, "scoped-table",
        projectionIdentity(table), ".rownames")
    validations <- 0L
    metadata_reads <- 0L
    original_validate <- validateSavedWorkflowReference
    original_metadata <- readSavedWorkflowTableMetadata
    testthat::with_mocked_bindings({
        scope <- newWorkflowTableReadScope(store, reference)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
        expect_identical(workflowEntityCacheInfo(store)$open_resources, 1L)
        page <- readSavedWorkflowTable(store, reference, "logFC",
            limit = 25L, scope = scope)
        volcano <- readSavedWorkflowTable(store, reference,
            c("logFC", "fdr_qvalue"), limit = nrow(table), scope = scope)
        keys <- rownames(table)[c(8L, 2L, 8L)]
        selected <- readSavedWorkflowTable(store, reference, "edge",
            keys = keys, limit = 3L, max_rows = 3L, scope = scope)
        expect_identical(as.numeric(page$data$logFC), as.numeric(table$logFC[1:25]))
        expect_identical(volcano$data$fdr_qvalue, table$fdr_qvalue)
        expect_identical(selected$source_row_names, keys)
        expect_identical(validations, 1L)
        expect_identical(metadata_reads, 1L)
        expect_error(closeWorkflowEntityStore(store), class = "saved_entity_scope")
        releaseWorkflowTableReadScope(scope)
        expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
        expect_true(releaseWorkflowTableReadScope(scope))
        expect_error(readSavedWorkflowTable(store, reference, "logFC", scope = scope),
            class = "saved_entity_scope")
    }, validateSavedWorkflowReference = function(...) {
        validations <<- validations + 1L
        original_validate(...)
    }, readSavedWorkflowTableMetadata = function(...) {
        metadata_reads <<- metadata_reads + 1L
        original_metadata(...)
    }, .package = "MultiScholaR")
    closeWorkflowEntityStore(store)
})

test_that("scoped tables fail closed after mutation or project replacement", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "scope-integrity", "DIA")
    reference <- saveWorkflowTable(store, table, "scope-integrity-table",
        projectionIdentity(table), ".rownames")
    scope <- newWorkflowTableReadScope(store, reference)
    expect_identical(as.numeric(readSavedWorkflowTable(store, reference, "logFC",
        limit = 2L, scope = scope)$data$logFC), as.numeric(table$logFC[1:2]))
    path <- file.path(store$directory, reference$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.size(path))
    changed <- bytes
    changed[[length(changed)]] <- as.raw(bitwXor(as.integer(changed[[length(changed)]]), 1L))
    writeBin(changed, path)
    expect_error(readSavedWorkflowTable(store, reference, "logFC", scope = scope),
        class = "saved_entity_corrupt")
    expect_error(releaseWorkflowTableReadScope(scope), class = "saved_entity_corrupt")
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    writeBin(bytes, path)
    replacement <- newWorkflowTableReadScope(store, reference)
    marker <- file.path(store$directory, "store.json")
    identity <- jsonlite::read_json(marker)
    identity$incarnation <- strrep("a", 64L)
    jsonlite::write_json(identity, marker, auto_unbox = TRUE)
    expect_error(readSavedWorkflowTable(store, reference, "logFC", scope = replacement),
        class = "saved_entity_scope")
    expect_error(releaseWorkflowTableReadScope(replacement), class = "saved_entity_scope")
    expect_identical(store$open_resources, 0L)
})

test_that("one-query DIA heatmap preserves the R full-population ranking", {
    table <- makeSavedProjectionDaTable()
    table$fdr_qvalue[1:5] <- c(NA_real_, NaN, Inf, -Inf, 0.5)
    table$logFC[6:9] <- c(NA_real_, NaN, Inf, -Inf)
    table$fdr_qvalue[10:12] <- 0.25
    table$logFC[10:12] <- 1
    store <- newWorkflowEntityStore(withr::local_tempdir(), "heatmap-project", "DIA")
    reference <- saveWorkflowTable(store, table, "heatmap-table",
        projectionIdentity(table), ".rownames")
    columns <- c("Protein.Group", "logFC", "fdr_qvalue", "edge")
    ordered <- order(table$fdr_qvalue, -abs(table$logFC), seq_len(nrow(table)),
        na.last = TRUE)
    expected <- table[ordered[1:50], columns, drop = FALSE]
    attr(expected$logFC, "units") <- attr(table$logFC, "units")
    attr(expected, "provenance") <- attr(table, "provenance")
    scope <- newWorkflowTableReadScope(store, reference)
    result <- readSavedWorkflowDiaHeatmapTable(store, reference, columns,
        50L, scope = scope)
    expect_identical(result$data, expected)
    expect_identical(result$source_row_names, rownames(expected))
    expect_identical(result$source_order, ordered[1:50])
    expect_identical(result$full_population_rows, nrow(table))
    expect_false(result$statistics_recomputed)
    expect_error(readSavedWorkflowDiaHeatmapTable(store, reference, columns,
        10001L, scope = scope), class = "saved_entity_limit")
    releaseWorkflowTableReadScope(scope)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    expect_error(readSavedWorkflowTable(store, reference, "not_a_column"),
        class = "saved_entity_limit")
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("query table production preserves the codec and commits without collecting output", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "query-writer", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    table <- data.frame(id = 1:4, value = c(NA_real_, NaN, Inf, -Inf),
        group = ordered(c("b", "a", "b", "a"), levels = c("b", "a", "unused")))
    attr(table, "source") <- "query-codec"
    reference <- saveWorkflowTable(store, table, "query-source", projectionIdentity(table),
        ".rownames", metadata_version = 2L)
    scope <- newWorkflowTableReadScope(store, reference)
    withr::defer(releaseWorkflowTableReadScope(scope))
    sql <- paste0("SELECT * FROM read_parquet(", DBI::dbQuoteString(scope$connection, scope$path),
        ") ORDER BY \".multischolar_order\"")
    saved <- saveWorkflowTableQuery(scope, sql, scope$metadata, "query-result", projectionIdentity(table))
    expect_identical(readSavedWorkflowTable(store, saved, names(table), limit = 4L)$data,
        readSavedWorkflowTable(store, reference, names(table), limit = 4L)$data)
    expect_error(testthat::with_mocked_bindings(
        saveWorkflowTableQuery(scope, sql, scope$metadata, "query-failed", projectionIdentity(table)),
        savedWorkflowStorageBoundary = function(phase) if (phase == "after_payload") stop("interrupted")),
        "interrupted")
    expect_null(recoverSavedWorkflowRequest(store, "query-failed"))
    parser <- prepareDIAParserTable(readr::read_tsv(I("Run\tx\na\t1\nb\t2"), show_col_types = FALSE))
    parsed <- saveWorkflowTable(store, parser, "parser", projectionIdentity(parser),
        ".rownames", metadata_version = 2L)
    expect_identical(readSavedWorkflowTable(store, parsed, names(parser), limit = 2L)$data, parser)
})

test_that("oversized queries fail without spill and preserve accepted data", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "spill", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    table <- data.frame(id = 1L, payload = "template")
    source <- saveWorkflowTable(store, table, "spill-source", projectionIdentity(table),
        ".rownames", metadata_version = 2L)
    before <- readSavedWorkflowTable(store, source, names(table))$data
    scope <- newWorkflowTableReadScope(store, source)
    withr::defer(releaseWorkflowTableReadScope(scope))
    spill <- scope$spill_directory
    connection <- scope$connection
    expect_identical(DBI::dbGetQuery(connection,
        "SELECT current_setting('max_temp_directory_size') AS cap")$cap, "0 bytes")
    DBI::dbExecute(connection, "SET threads = 1")
    DBI::dbExecute(connection, "SET memory_limit = '32MB'")
    expect_error(DBI::dbExecute(connection, paste0("CREATE TEMP TABLE spill_input AS SELECT ",
        "i::INTEGER AS id, repeat(md5(i::VARCHAR), 4) AS payload FROM range(1, 500001) AS t(i)")),
        "memory|Memory|temporary directory")
    expect_equal(nrow(DBI::dbGetQuery(connection, "SELECT * FROM duckdb_temporary_files()")), 0L)
    expect_length(list.files(spill, all.files = TRUE, no.. = TRUE), 0L)
    DBI::dbExecute(connection, "SET memory_limit = '1GB'")
    query <- paste0("SELECT * FROM read_parquet(", DBI::dbQuoteString(connection, scope$path), ")")
    expect_error(testthat::with_mocked_bindings(
        saveWorkflowTableQuery(scope, query, scope$metadata, "failed-output", projectionIdentity(table)),
        savedWorkflowStorageBoundary = function(phase) {
            if (phase == "after_payload") stop("writer interrupted")
        }), "writer interrupted")
    expect_null(recoverSavedWorkflowRequest(store, "failed-output"))
    releaseWorkflowTableReadScope(scope)
    expect_false(dir.exists(spill))
    expect_identical(store$open_resources, 0L)
    expect_identical(readSavedWorkflowTable(store, source, names(table))$data, before)
})

test_that("complete table streams respect batch limits and release failed cursors", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "stream", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    table <- data.frame(id = seq_len(20001L), quantity = rep(c(NA_real_, NaN, Inf, 0), length.out = 20001L))
    reference <- saveWorkflowTable(store, table, "stream-source", projectionIdentity(table),
        ".rownames", metadata_version = 2L)
    scope <- newWorkflowTableReadScope(store, reference)
    withr::defer(releaseWorkflowTableReadScope(scope))
    expect_identical(readSavedWorkflowTablePrototype(scope), table[0, , drop = FALSE])
    batches <- integer()
    expect_identical(streamSavedWorkflowTable(scope, names(table), function(data, offset, row_names) {
        batches <<- c(batches, nrow(data))
        at <- seq.int(offset + 1L, offset + nrow(data))
        expect_identical(data$id, table$id[at])
        expect_identical(data$quantity, table$quantity[at])
    }), 20001L)
    expect_identical(batches, c(10000L, 10000L, 1L))
    expect_error(streamSavedWorkflowTable(scope, names(table), function(...) stop("consumer failed")),
        "consumer failed")
    expect_identical(DBI::dbGetQuery(scope$connection, "SELECT 1::INTEGER AS n")$n, 1L)
})

test_that("streaming v3 tables share validation without recovery I/O", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "shared-io", "matrix")
    reference <- saveWorkflowTable(store, table, "matrix", projectionIdentity(table),
        ".rownames", metadata_version = 3L)
    path <- file.path(store$directory, reference$payload$relative_path)
    expect_identical(readSavedWorkflowTableMetadata(path)$schema, 3L)
    expect_identical(readSavedWorkflowTable(store, reference, names(table), limit = nrow(table))$data, table)
    before <- store$io_counts
    metadata <- newWorkflowTableReadScope(store, reference, metadata_only = TRUE)
    prototype <- readSavedWorkflowTablePrototype(metadata)
    expect_equal(nrow(prototype), 0L)
    expect_identical(names(prototype), names(table))
    expect_identical(attr(prototype$logFC, "units"), attr(table$logFC, "units"))
    expect_identical(levels(prototype$comparison), levels(table$comparison))
    expect_null(metadata$connection)
    releaseWorkflowTableReadScope(metadata)
    expect_identical(store$io_counts$payload_hash_bytes, before$payload_hash_bytes)
    expect_identical(store$io_counts$recovery_flushes, before$recovery_flushes)
    action <- beginWorkflowTableAction(store)
    first <- newWorkflowTableReadScope(store, reference)
    second <- newWorkflowTableReadScope(store, reference)
    expect_identical(first$connection, second$connection)
    expect_equal(store$io_counts$payload_hash_bytes - before$payload_hash_bytes, reference$payload$bytes)
    releaseWorkflowTableReadScope(first)
    expect_true(DBI::dbIsValid(second$connection))
    releaseWorkflowTableReadScope(second)
    endWorkflowTableAction(store, action)
    expect_null(store$table_action)
    expect_identical(store$open_resources, 0L)
    expect_identical(store$io_counts$recovery_flushes, before$recovery_flushes)
})

test_that("shared group cursor and change writer preserve interleaved source order", {
    table <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))@peptide_data
    # Reorder the existing fixture so complete groups cross source and fetch boundaries.
    table <- table[order(table$Run, table$Stripped.Sequence), , drop = FALSE]
    store <- newWorkflowEntityStore(withr::local_tempdir(), "group-io", "long-table")
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table),
        ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference)
    withr::defer(releaseWorkflowTableReadScope(scope))
    cursor <- testthat::with_mocked_bindings(
        openWorkflowTableGroups(scope, c("Protein.Group", "Run", "Precursor.Quantity"), "Protein.Group"),
        dbSendQuery = function(...) stop("group production ran a global SQL query"), .package = "DBI")
    withr::defer(cursor$close())
    cached <- testthat::with_mocked_bindings(saveWorkflowGroupIndex(scope, "Protein.Group"),
        streamSavedWorkflowTable = function(...) stop("valid group index was rebuilt"))
    expect_identical(cached$reference, cursor$dependencies[[1L]])
    offset <- 0L
    observed <- integer()
    repeat {
        group <- cursor$next_group()
        if (is.null(group)) break
        positions <- attr(group, "workflow_source_order")
        expected <- table[positions, c("Protein.Group", "Run", "Precursor.Quantity"), drop = FALSE]
        attr(group, "workflow_source_order") <- NULL
        expect_equal(group, expected)
        for (key in unique(group$Protein.Group)) {
            expect_equal(sum(group$Protein.Group == key), sum(table$Protein.Group == key))
        }
        observed <- c(observed, positions)
        map <- data.frame(source_row = positions)
        for (start in seq.int(1L, nrow(map), by = 10000L)) {
            at <- seq.int(start, min(nrow(map), start + 9999L))
            offset <- appendWorkflowTableBatch(scope, "selected_rows", map[at, , drop = FALSE], offset)
        }
    }
    cursor$close()
    expect_identical(sort(observed), seq_len(nrow(table)))
    expect_gte(store$table_rows, nrow(table))
    expect_gt(store$table_bytes, 0)
    expect_length(DBI::dbListTables(scope$connection), 0L)
    map_path <- finishWorkflowTableBatch(scope, "selected_rows", projectionIdentity(table, "selection"))
    metadata <- workflowSavedTableMetadata(table, ".rownames", 3L)
    expect_true(workflowTableChangesAreIdentity(scope, map_path, metadata, "source_row"))
    changed_metadata <- metadata
    changed_metadata$attributes$changed <- TRUE
    expect_false(workflowTableChangesAreIdentity(scope, map_path, changed_metadata, "source_row"))
    result <- saveWorkflowTable(store, NULL, "result", projectionIdentity(table, "identity-map"),
        ".rownames", metadata = metadata, write_fn = function(path, restoration)
            writeWorkflowTableChanges(scope, map_path, NULL, restoration, path, source_column = "source_row"))
    releaseWorkflowTableReadScope(scope)
    expect_identical(readSavedWorkflowTable(store, result, names(table), limit = min(10000L, nrow(table)))$data,
        table[seq_len(min(10000L, nrow(table))), , drop = FALSE])
    expect_identical(store$io_counts$query_streams, 0L)
    expect_lt(store$io_counts$staged_bytes, reference$payload$bytes)
    expect_identical(store$open_resources, 0L)
})

test_that("incremental writers abort without replacing an accepted payload", {
    table <- makeSavedProjectionDaTable()
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "abort-io", "matrix")
    reference <- saveWorkflowTable(store, table, "accepted", projectionIdentity(table), ".rownames", metadata_version = 3L)
    accepted <- file.path(store$directory, reference$payload$relative_path)
    checksum <- digest::digest(file = accepted, algo = "sha256")
    path <- file.path(root, "incomplete.parquet")
    writer <- openWorkflowTableWriter(path, workflowSavedTableMetadata(table, ".rownames", 3L))
    appendWorkflowTableWriter(writer, table, 1:10)
    abortWorkflowTableWriter(writer)
    expect_false(file.exists(path))
    expect_identical(digest::digest(file = accepted, algo = "sha256"), checksum)
    expect_identical(recoverSavedWorkflowRequest(store, "accepted"), reference)
})


test_that("bounded read batches share larger Parquet row groups", {
    root <- withr::local_tempdir()
    table <- data.frame(id = seq_len(25001L), value = rep(c(NA_real_, NaN, Inf, 0),
        length.out = 25001L))
    path <- file.path(root, "rows.parquet")
    writeSavedWorkflowParquet(table, ".rownames", path, metadata_version = 3L)
    reader <- arrow::ParquetFileReader$create(path)
    withr::defer(reader$.unsafe_delete())
    expect_equal(reader$num_row_groups, 1L)
    expect_equal(readSavedWorkflowTableMetadata(path)$rows, nrow(table))
    physical <- as.data.frame(reader$ReadTable())
    expect_identical(physical$id, table$id)
    expect_identical(physical$value, table$value)
})


test_that("cancelled Arrow streams destroy the reader without draining its tail", {
    table <- data.frame(id = seq_len(10L))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "cancel-stream", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table),
        ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference)
    withr::defer(releaseWorkflowTableReadScope(scope))
    fetched <- closed <- released <- 0L
    release <- function() released <<- released + 1L
    reader <- list(read_next_batch = function() {
        fetched <<- fetched + 1L
        arrow::record_batch(id = 1L, .multischolar_order = 1L, .multischolar_row_key = "1")
    }, Close = function() closed <<- closed + 1L, .unsafe_delete = release)
    scanner <- list(ToRecordBatchReader = function() reader, .unsafe_delete = release)
    testthat::local_mocked_bindings(
        open_dataset = function(...) list(.unsafe_delete = release),
        Scanner = list(create = function(...) scanner), .package = "arrow")
    expect_error(streamSavedWorkflowTable(scope, "id", function(...) stop("cancelled")), "cancelled")
    expect_identical(fetched, 1L)
    expect_identical(closed, 0L)
    expect_identical(released, 3L)
})

test_that("cohort reads materialize only selected index and assay rows", {
    table <- makeWorkflowPrecursorChunk(readWorkflowCheckpointRecipe("representative"), seq_len(256L))
    table <- table[order(table$Run, table$Stripped.Sequence), , drop = FALSE]
    store <- newWorkflowEntityStore(withr::local_tempdir(), "cohort-selection", "long-table")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table),
        ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    index <- saveWorkflowGroupIndex(scope, "Protein.Group")
    directory <- readWorkflowIndexDirectory(store, index$value$groups)
    columns <- c("Protein.Group", "Run", "Precursor.Quantity")
    before <- store$table_rows
    streams <- store$io_counts$sequential_streams
    groups <- split(seq_len(nrow(directory)), seq_len(nrow(directory)) %% 4L)
    original_stream <- streamSavedWorkflowTable
    control_rows <- store$table_rows
    control_time <- system.time(testthat::with_mocked_bindings({
        for (selected in groups) readWorkflowGroupCohort(scope, columns, index$value$positions,
            directory[rev(selected), , drop = FALSE])
    }, streamSavedWorkflowTable = function(scope, columns, consume, ..., source_rows = NULL) {
        if (identical(scope$reference, reference)) {
            delivered <- 0L
            original_stream(scope, columns, function(data, offset, row_names) {
                physical <- offset + seq_len(nrow(data))
                keep <- physical %in% source_rows
                part <- data[keep, , drop = FALSE]
                attr(part, "dia_source_order") <- physical[keep]
                consume(part, delivered, row_names[keep])
                delivered <<- delivered + nrow(part)
            }, ...)
        } else original_stream(scope, columns, consume, ..., source_rows = source_rows)
    }))[["elapsed"]]
    expect_equal(store$table_rows - control_rows, (length(groups) + 1L) * nrow(table))
    before <- store$table_rows
    streams <- store$io_counts$sequential_streams
    elapsed <- system.time(for (selected in groups) {
        # Reversed groups exercise index physical order versus scientific order.
        group_rows <- directory[rev(selected), , drop = FALSE]
        value <- readWorkflowGroupCohort(scope, columns, index$value$positions, group_rows)
        expected_rows <- unlist(lapply(group_rows$value, function(key)
            which(table$Protein.Group == key)), use.names = FALSE)
        expect_identical(attr(value, "workflow_source_order"), expected_rows)
        attr(value, "workflow_source_order") <- NULL
        expect_equal(value, table[expected_rows, columns, drop = FALSE])
    })[["elapsed"]]
    # One index row and one assay row per selected row, independent of cohort count.
    expect_equal(store$table_rows - before, 2 * nrow(table))
    expect_equal(store$io_counts$sequential_streams - streams, 2L * length(groups))
    message(sprintf(paste0("cohort selection: %d rows, %d cohorts, ",
        "full-assay materialization %.3f s; pushed selection %.3f s"),
        nrow(table), length(groups), control_time, elapsed))
})

test_that("ordinal filter streams preserve sparse, empty and nonfinite values", {
    table <- makeSavedProjectionDaTable()
    store <- newWorkflowEntityStore(withr::local_tempdir(), "ordinal-filter", "table")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table),
        ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    for (rows in list(integer(), seq.int(1L, nrow(table), by = 2L), seq_len(nrow(table)))) {
        observed <- integer()
        streamSavedWorkflowTable(scope, names(table), function(data, offset, row_names) {
            actual <- attr(data, "dia_source_order")
            if (length(actual)) {
                attr(data, "dia_source_order") <- NULL
                expected <- table[actual, , drop = FALSE]
                attr(expected$logFC, "units") <- attr(table$logFC, "units")
                expect_equal(data, expected)
                expect_lte(nrow(data), 10000L)
                observed <<- c(observed, actual)
            } else expect_equal(nrow(data), 0L)
        }, source_rows = rows)
        expect_identical(observed, rows)
    }
    expect_error(streamSavedWorkflowTable(scope, names(table), identity,
        source_rows = c(2L, 1L)), class = "saved_entity_limit")
    # Exercise the bounded membership-expression branch without a large assay.
    selected <- seq.int(1L, 4097L, by = 2L)
    scan <- arrow::Scanner$create(arrow::Table$create(data.frame(.multischolar_order = 1:4097)),
        filter = workflowTableRowFilter(selected), batch_size = 10000L)
    expect_identical(as.data.frame(scan$ToTable())$.multischolar_order, selected)
})

test_that("sparse Arrow selections coalesce before restoration and ordered reads own one cursor", {
    table <- data.frame(id = seq_len(40000L), group = rep(1:400, 100L), value = seq_len(40000L) / 10)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "cursor-batches", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table), ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    selected <- which(table$group == 1L)
    batches <- 0L
    observed <- integer()
    streamSavedWorkflowTable(scope, "id", function(data, offset, row_names) {
        batches <<- batches + 1L
        observed <<- c(observed, data$id)
    }, source_rows = selected)
    expect_identical(observed, selected)
    expect_identical(batches, 1L)
    expect_length(as.list(scope$readers), 0L)
    before <- store$io_counts$sequential_streams
    cursor <- openWorkflowOrderedRows(store, reference, names(table), ".multischolar_order")
    withr::defer(cursor$close())
    pages <- testthat::with_mocked_bindings({
        out <- list()
        repeat {
            part <- cursor$next_batch()
            if (is.null(part)) break
            out[[length(out) + 1L]] <- part
        }
        out
    }, readSavedWorkflowTable = function(...) stop("ordered stream reopened a page reader"))
    expected <- table
    rownames(expected) <- as.character(seq_len(nrow(expected)))
    expect_identical(do.call(rbind, pages), expected)
    expect_equal(store$io_counts$sequential_streams - before, 1L)
    cursor <- openWorkflowArrowCursor(scope, "id")
    releaseWorkflowTableReadScope(scope)
    expect_true(cursor$closed)
    expect_null(cursor$fetch())
})

test_that("contiguous groups use one cursor and preserve complete large groups", {
    table <- data.frame(id = seq_len(32000L), group = rep(c("B", "A", "C"), c(12001L, 9999L, 10000L)))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "contiguous-groups", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, table, "input", projectionIdentity(table), ".rownames", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    saveWorkflowGroupIndex(scope, "group")
    before <- store$io_counts$sequential_streams
    cursor <- openWorkflowTableGroups(scope, names(table), "group")
    withr::defer(cursor$close())
    observed <- testthat::with_mocked_bindings({
        parts <- list()
        repeat {
            data <- cursor$next_group()
            if (is.null(data)) break
            for (key in unique(data$group)) expect_equal(sum(data$group == key), sum(table$group == key))
            parts[[length(parts) + 1L]] <- data$id
        }
        unlist(parts, use.names = FALSE)
    }, readWorkflowGroupCohort = function(...) stop("contiguous groups used scattered gathering"))
    expect_identical(observed, table$id)
    expect_equal(store$io_counts$sequential_streams - before, 2L) # directory and one assay cursor
    expect_length(as.list(scope$readers), 0L)
})
test_that("private table buffers preserve destinations and restoration through scattered writes", {
    rows <- 20000L
    source <- data.frame(quantity = rep(c(NA_real_, NaN, Inf, -Inf, 0, 2.5), length.out = rows),
        named = seq_len(rows), group = ordered(rep(c("b", "a"), length.out = rows), levels = c("b", "a")),
        flag = rep(c(TRUE, FALSE, NA), length.out = rows), text = rep(c("x", NA_character_), length.out = rows))
    names(source$named) <- paste0("n", seq_len(rows))
    rownames(source) <- paste0("r", seq_len(rows))
    attr(source, "provenance") <- list(origin = "buffer restoration fixture")
    metadata <- workflowSavedTableMetadata(source, ".rownames", 3L)
    buffer <- newWorkflowTableBuffer(metadata, names(source), rows)
    quantity_id <- tracemem(buffer$columns$quantity)
    position_id <- tracemem(buffer$positions)
    traces <- utils::capture.output(for (start in c(1L, 10001L)) {
        from <- seq.int(start, length.out = 10000L)
        to <- if (start == 1L) from + 10000L else from - 10000L
        assignWorkflowTableBuffer(buffer, source[from, , drop = FALSE], to, from, rownames(source)[from])
        expect_identical(tracemem(buffer$columns$quantity), quantity_id)
        expect_identical(tracemem(buffer$positions), position_id)
    })
    expect_false(any(startsWith(traces, "tracemem[")))
    observed <- finishWorkflowTableBuffer(buffer)
    order <- c(10001:20000, 1:10000)
    expected <- source[order, , drop = FALSE]
    attr(expected, "workflow_source_order") <- order
    expect_identical(observed, expected)
})
