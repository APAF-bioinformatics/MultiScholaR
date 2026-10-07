testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"),
    "durable saved tables require Linux")

test_that("selected tables preserve bounded logical rows and typed columns", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "perf103-selection", "DIA")
    on.exit({
        for (name in c("scope", "selected", "child_scope", "export_scope",
            "unicode_scope")) {
            if (exists(name, inherits = FALSE)) {
                owned <- get(name, inherits = FALSE)
                if (!isTRUE(owned$released)) releaseWorkflowTableReadScope(owned)
            }
        }
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    source <- data.frame(id = c("a", "b", "c", "d"),
        value = c(1, NA_real_, 3, 4),
        group = factor(c("x", "y", "x", "y"), levels = c("x", "unused", "y")))
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-source", "design", "selection", selector, "perf103-test")
    physical <- saveWorkflowTable(store, source, "source", identity("source"),
        "id", metadata_version = 3L)
    map <- data.frame(logical_ordinal = as.double(1:2),
        source_ordinal = c(3, 1))
    row_map <- saveWorkflowTable(store, map, "map", identity("map"),
        "logical_ordinal", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    expected <- data.frame(renamed = c(3, 1), group = c("x", "x"))
    metadata <- workflowSavedTableMetadata(expected, "renamed", 3L)
    mapping <- data.frame(output = c("renamed", "group"),
        source = c("value", "group"),
        cast = c("identity", "factor_to_character"))
    view <- saveWorkflowTableSelection(scope, row_map, mapping, metadata,
        "view", identity("view"))
    malformed <- metadata
    malformed$row_name_column <- ".multischolar_row_names"
    expect_error(saveWorkflowTableSelection(scope, row_map, mapping, malformed,
        "invalid-row-names", identity("invalid-row-names")),
        class = "saved_entity_corrupt")
    malformed <- metadata
    malformed$attributes <- new.env(parent = emptyenv())
    expect_error(saveWorkflowTableSelection(scope, row_map, mapping, malformed,
        "invalid-attributes", identity("invalid-attributes")),
        class = "saved_entity_corrupt")
    unicode_expected <- expected
    names(unicode_expected)[[1L]] <- "β_quantity"
    unicode_mapping <- mapping
    unicode_mapping$output[[1L]] <- "β_quantity"
    unicode_view <- saveWorkflowTableSelection(scope, row_map, unicode_mapping,
        workflowSavedTableMetadata(unicode_expected, "β_quantity", 3L),
        "unicode-view", identity("unicode-view"))
    releaseWorkflowTableReadScope(scope)
    unicode_scope <- newWorkflowTableReadScope(store, unicode_view, reader = "arrow")
    expect_identical(readSavedWorkflowTable(store, unicode_view,
        names(unicode_expected), limit = 2L, scope = unicode_scope)$data,
        unicode_expected)
    releaseWorkflowTableReadScope(unicode_scope)
    selected <- newWorkflowTableReadScope(store, view, reader = "arrow")
    page <- readSavedWorkflowTable(store, view, c("renamed", "group"),
        limit = 2L, scope = selected)
    expect_identical(page$data, expected)
    expect_identical(as.double(page$source_order), c(3, 1))
    streamed <- list()
    streamSavedWorkflowTable(selected, c("renamed", "group"),
        function(data, offset, row_names) streamed[[length(streamed) + 1L]] <<- data)
    expect_identical(streamed[[1L]], expected)
    expect_identical(workflowTablePopulationId(selected), selected$view_id)
    expect_identical(materializeDIATable(list(store = store), view, scope = selected),
        expected)
    reads_before <- store$table_reads
    testthat::with_mocked_bindings({
        expect_error(materializeDIATable(list(store = store), view, scope = selected),
            class = "saved_entity_limit")
    }, workflowProjectedRowBytes = function(...) 268435456)
    expect_identical(store$table_reads, reads_before)
    expect_error(materializeDIATable(list(store = store), view, scope = selected,
        relation = "physical_input", filter_column = "renamed", filter_value = 1),
        class = "saved_entity_unsupported")
    sql_scope <- newWorkflowTableReadScope(store, view, reader = "duckdb")
    expect_error(saveWorkflowTableQuery(sql_scope, "SELECT 1", metadata,
        "forbidden-query", identity("forbidden-query")),
        class = "saved_entity_unsupported")
    releaseWorkflowTableReadScope(sql_scope)
    expect_error(readWorkflowSelectedOrdinals(store, view, "renamed", 1),
        class = "saved_entity_unsupported")
    child_map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(2, 1)),
        "child-map", identity("child-map"), "logical_ordinal", metadata_version = 3L)
    child_expected <- expected[2:1, , drop = FALSE]
    row.names(child_expected) <- NULL
    child_metadata <- workflowSavedTableMetadata(child_expected, "renamed", 3L)
    child <- saveWorkflowTableSelection(selected, child_map,
        data.frame(output = c("renamed", "group"),
            source = c("renamed", "group"), cast = c("identity", "identity")),
        child_metadata, "child-view", identity("child-view"))
    releaseWorkflowTableReadScope(selected)
    child_scope <- newWorkflowTableReadScope(store, child, reader = "arrow")
    expect_identical(readSavedWorkflowTable(store, child,
        c("renamed", "group"), limit = 2L, scope = child_scope)$data,
        child_expected)
    expect_identical(as.double(readSavedWorkflowTable(store, child,
        "renamed", limit = 2L, scope = child_scope)$source_order), c(1, 3))
    invalid_map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(1, 1)),
        "invalid-map", identity("invalid-map"), "logical_ordinal",
        metadata_version = 3L)
    expect_error(saveWorkflowTableSelection(child_scope, invalid_map,
        data.frame(output = c("renamed", "group"),
            source = c("renamed", "group"), cast = c("identity", "identity")),
        child_metadata, "invalid-view", identity("invalid-view")),
        class = "saved_entity_corrupt")
    invalid_range <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(5, 1)),
        "invalid-range", identity("invalid-range"), "logical_ordinal",
        metadata_version = 3L)
    expect_error(saveWorkflowTableSelection(child_scope, invalid_range,
        data.frame(output = c("renamed", "group"),
            source = c("renamed", "group"), cast = c("identity", "identity")),
        child_metadata, "invalid-range-view", identity("invalid-range-view")),
        class = "saved_entity_corrupt")
    qualified <- tryCatch(workflowFastWriterInfo(),
        saved_entity_unsupported = function(error) NULL)
    if (is.null(qualified)) {
        expect_error(saveWorkflowTablePhysicalExport(child_scope, "export",
            identity("export")), class = "saved_entity_unsupported")
        expect_null(recoverSavedWorkflowRequest(store, "export"))
    } else {
        exported <- saveWorkflowTablePhysicalExport(child_scope, "export",
            identity("export"))
    }
    releaseWorkflowTableReadScope(child_scope)
    if (!is.null(qualified)) {
        export_scope <- newWorkflowTableReadScope(store, exported, reader = "duckdb")
        exported_page <- readSavedWorkflowTable(store, exported,
            c("renamed", "group"), limit = 2L, scope = export_scope)$data
        expect_identical(unname(exported_page$renamed), child_expected$renamed)
        expect_identical(unname(exported_page$group), child_expected$group)
        expect_identical(row.names(exported_page), row.names(child_expected))
        releaseWorkflowTableReadScope(export_scope)
        old_library <- Sys.getenv("PERF103_OLD_LIB")
        if (nzchar(old_library)) {
            old_result <- callr::r(function(root, exported, view, old_library) {
                library(MultiScholaR, lib.loc = old_library)
                api <- function(name) get(name, envir = asNamespace("MultiScholaR"))
                old_store <- api("newWorkflowEntityStore")(root, "perf103-selection", "DIA")
                on.exit(api("closeWorkflowEntityStore")(old_store), add = TRUE)
                old_scope <- api("newWorkflowTableReadScope")(old_store, exported,
                    reader = "duckdb")
                on.exit(api("releaseWorkflowTableReadScope")(old_scope), add = TRUE)
                values <- api("readSavedWorkflowTable")(old_store, exported,
                    c("renamed", "group"), limit = 2L, scope = old_scope)$data
                api("releaseWorkflowTableReadScope")(old_scope)
                refusal <- tryCatch({
                    api("newWorkflowTableReadScope")(old_store, view, reader = "duckdb")
                    FALSE
                }, error = function(error)
                    inherits(error, "saved_entity_unsupported") ||
                        inherits(error, "saved_entity_corrupt"))
                list(values = values, refusal = refusal)
            }, args = list(root = root, exported = exported, view = view,
                old_library = old_library), libpath = c(old_library, .libPaths()))
            expect_identical(unname(old_result$values$renamed), child_expected$renamed)
            expect_identical(unname(old_result$values$group), child_expected$group)
            expect_identical(row.names(old_result$values), row.names(child_expected))
            expect_true(old_result$refusal)
        }
    }
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("empty selected tables keep their schema and produce a physical export", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-empty", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-empty", "design", "selection", selector, "perf103-test")
    source <- data.frame(id = character(), q = double(),
        group = factor(character(), levels = c("x", "unused")))
    physical <- saveWorkflowTable(store, source, "source", identity("source"),
        "id", metadata_version = 3L)
    map <- saveWorkflowTable(store, data.frame(logical_ordinal = double(),
        source_ordinal = double()), "map", identity("map"),
        "logical_ordinal", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    expected <- data.frame(id = character(), q = double(), group = character())
    view <- saveWorkflowTableSelection(scopes[[1L]], map,
        data.frame(output = c("id", "q", "group"), source = c("id", "q", "group"),
            cast = c("identity", "identity", "factor_to_character")),
        workflowSavedTableMetadata(expected, "id", 3L), "view", identity("view"))
    releaseWorkflowTableReadScope(scopes[[1L]])
    scopes[[2L]] <- newWorkflowTableReadScope(store, view, reader = "arrow")
    expect_identical(readSavedWorkflowTable(store, view, names(expected),
        limit = 1L, scope = scopes[[2L]])$data, expected)
    batches <- list()
    streamSavedWorkflowTable(scopes[[2L]], names(expected),
        function(data, offset, row_names) batches[[length(batches) + 1L]] <<- data)
    expect_identical(batches, list(expected))
    exported <- saveWorkflowTablePhysicalExport(scopes[[2L]], "export", identity("export"))
    releaseWorkflowTableReadScope(scopes[[2L]])
    scopes[[3L]] <- newWorkflowTableReadScope(store, exported, reader = "duckdb")
    exported_data <- readSavedWorkflowTable(store, exported, names(expected),
        limit = 1L, scope = scopes[[3L]])$data
    expect_identical(unclass(exported_data)[names(expected)],
        unclass(expected)[names(expected)])
    expect_identical(class(exported_data), class(expected))
    expect_identical(row.names(exported_data), row.names(expected))
    releaseWorkflowTableReadScope(scopes[[3L]])
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("wide selected tables stream ordered columns without changing their source", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-wide", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-wide", "design", "selection", selector, "perf103-test")
    source <- data.frame(id = paste0("P", 1:4))
    for (i in seq_len(130L)) source[[paste0("sample_", i)]] <- as.double(seq_len(4L) + i)
    physical <- saveWorkflowTable(store, source, "source", identity("source"),
        "id", metadata_version = 3L)
    physical_path <- file.path(store$directory, physical$payload$relative_path)
    source_hash <- digest::digest(file = physical_path, algo = "sha256")
    map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(4, 2)),
        "map", identity("map"), "logical_ordinal", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    expected <- source[c(4L, 2L), , drop = FALSE]
    row.names(expected) <- NULL
    fields <- names(source)
    view <- saveWorkflowTableSelection(scopes[[1L]], map,
        data.frame(output = fields, source = fields, cast = rep("identity", length(fields))),
        workflowSavedTableMetadata(expected, "id", 3L), "view", identity("view"))
    releaseWorkflowTableReadScope(scopes[[1L]])
    scopes[[2L]] <- newWorkflowTableReadScope(store, view, reader = "arrow")
    batches <- list()
    streamWorkflowTableRows(scopes[[2L]], fields,
        function(data, offset, row_names) batches[[length(batches) + 1L]] <<- data)
    expect_identical(batches, list(expected))
    releaseWorkflowTableReadScope(scopes[[2L]])
    expect_identical(digest::digest(file = physical_path, algo = "sha256"), source_hash)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("large filtering counts cross the worker boundary by reference", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-counts", "DIA")
    scope <- NULL
    on.exit({
        if (!is.null(scope) && !isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    n <- 3200L
    table <- data.frame(Run = rep("S1", n),
        Protein.Ids = paste0(sprintf("P%04d", seq_len(n)), strrep("A", 1600L)),
        Stripped.Sequence = sprintf("PEP%04d", seq_len(n)))
    identity <- newSavedWorkflowScientificIdentity("perf103-counts", "design",
        "counts", "peptide-table", "perf103-test")
    reference <- saveWorkflowTable(store, table, "peptides", identity,
        "Protein.Ids", metadata_version = 3L)
    accumulator <- newDIAFilteringAccumulator(table[0, , drop = FALSE], "Protein.Ids")
    accumulator$push(table)
    counted <- accumulator$finish()
    expect_gt(as.numeric(object.size(counted)), 4 * 1024^2)
    saved <- saveDIALargeFilteringCountsReference(store,
        list(reference = reference, counts = counted), "Protein.Ids")
    expect_true(isSavedWorkflowReference(saved))
    expect_lt(as.numeric(object.size(saved)), 65536)
    scope <- newWorkflowTableReadScope(store, reference, reader = "arrow")
    reads_before <- store$table_reads
    restored <- summarizeSavedDIAFiltering(scope, NULL, "Protein.Ids", saved)
    expect_identical(restored$summary, counted$summary)
    expect_identical(restored$rows, counted$rows)
    expect_identical(restored$proteins, counted$proteins)
    expect_identical(store$table_reads, reads_before)
    expect_identical(restored$dependencies[[1L]], saved)
    releaseWorkflowTableReadScope(scope)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("a real wide-annotation rollup returns its counts by reference", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(),
        "perf103-wide-rollup", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    output_scope <- NULL
    on.exit({
        if (!is.null(output_scope) && !isTRUE(output_scope$released))
            releaseWorkflowTableReadScope(output_scope)
        closeDIAStateAdapter(manager)
    }, add = TRUE)
    manager$setWorkflowType("DIA")
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    n <- nrow(object@peptide_data)
    object@peptide_data$Protein.Group <- sprintf("PG%05d", seq_len(n))
    object@peptide_data$Protein.Ids <- paste0(sprintf("ACC%05d", seq_len(n)),
        strrep("A", 1600L))
    manager$saveState("protein_peptide_filtered", object, list(), "wide annotations")
    before <- getDIAStateReference(manager)
    native <- rollUpPrecursorToPeptide(object, core_utilisation = NA)
    result <- computeDIAReferencePrecursorRollup(manager, before)
    expect_true(isSavedWorkflowReference(result$counts))
    expect_lt(as.numeric(object.size(result$counts)), 65536)
    expect_identical(result$summary,
        native@args$rollUpPrecursorToPeptide$rollup_summary)
    output <- materializeDIATable(manager$.context, result$reference)
    expect_identical(output, native@peptide_data)
    accumulator <- newDIAFilteringAccumulator(native@peptide_data[0, , drop = FALSE],
        object@protein_id_column)
    accumulator$push(native@peptide_data)
    expected_counts <- accumulator$finish()
    output_scope <- newWorkflowTableReadScope(manager$.context$store,
        result$reference, reader = "arrow")
    reads_before <- manager$.context$store$table_reads
    counted <- summarizeSavedDIAFiltering(output_scope, NULL,
        object@protein_id_column, result$counts)
    expect_identical(counted$summary, expected_counts$summary)
    expect_identical(counted$proteins, expected_counts$proteins)
    expect_identical(manager$.context$store$table_reads, reads_before)
    expect_identical(getDIAStateReference(manager), before)
    releaseWorkflowTableReadScope(output_scope)
    expect_identical(manager$.context$store$open_resources, 0L)
})

test_that("distinct logical populations cannot share filtering counts", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-populations", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-populations", "design", "counts", selector, "perf103-test")
    table <- data.frame(Run = rep("S1", 4L),
        Protein.Ids = c("P1", "P1", "P2", "P2"),
        Stripped.Sequence = c("A", "B", "C", "D"))
    physical <- saveWorkflowTable(store, table, "source", identity("source"),
        "Stripped.Sequence", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    mapping <- data.frame(output = names(table), source = names(table),
        cast = rep("identity", ncol(table)))
    views <- vector("list", 2L)
    maps <- vector("list", 2L)
    for (i in seq_len(2L)) {
        at <- if (i == 1L) c(1, 2) else c(1, 3)
        maps[[i]] <- saveWorkflowTable(store,
            data.frame(logical_ordinal = as.double(1:2), source_ordinal = as.double(at)),
            paste0("map-", i), identity(paste0("map-", i)),
            "logical_ordinal", metadata_version = 3L)
        expected <- table[at, , drop = FALSE]
        row.names(expected) <- NULL
        views[[i]] <- saveWorkflowTableSelection(scopes[[1L]], maps[[i]], mapping,
            workflowSavedTableMetadata(expected, "Stripped.Sequence", 3L),
            paste0("view-", i), identity(paste0("view-", i)))
    }
    releaseWorkflowTableReadScope(scopes[[1L]])
    counted <- vector("list", 2L)
    for (i in seq_len(2L)) {
        scopes[[i + 1L]] <- newWorkflowTableReadScope(store, views[[i]], reader = "arrow")
        counted[[i]] <- summarizeSavedDIAFiltering(scopes[[i + 1L]], NULL, "Protein.Ids")
        repeated <- summarizeSavedDIAFiltering(scopes[[i + 1L]], NULL, "Protein.Ids")
        expect_identical(repeated$dependencies[[1L]], counted[[i]]$dependencies[[1L]])
    }
    expect_false(identical(scopes[[2L]]$view_id, scopes[[3L]]$view_id))
    expect_false(identical(counted[[1L]]$dependencies[[1L]]$entity_id,
        counted[[2L]]$dependencies[[1L]]$entity_id))
    expect_identical(counted[[1L]]$summary$proteins, 1L)
    expect_identical(counted[[2L]]$summary$proteins, 2L)
    closure <- expandWorkflowTableReferences(store,
        stats::setNames(list(views[[1L]], views[[2L]]),
            vapply(views, `[[`, character(1), "entity_id")))
    expect_true(all(vapply(c(list(physical), maps), `[[`, character(1), "entity_id") %in%
        names(closure)))
    for (scope in scopes[-1L]) releaseWorkflowTableReadScope(scope)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("interrupted physical export leaves no accepted result or open owner", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-failure", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-failure", "design", "export", selector, "perf103-test")
    table <- data.frame(id = paste0("P", 1:4), quantity = as.double(1:4))
    physical <- saveWorkflowTable(store, table, "source", identity("source"),
        "id", metadata_version = 3L)
    physical_path <- file.path(store$directory, physical$payload$relative_path)
    input_hash <- digest::digest(file = physical_path, algo = "sha256")
    map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(4, 1)),
        "map", identity("map"), "logical_ordinal", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    expected <- table[c(4, 1), , drop = FALSE]
    row.names(expected) <- NULL
    fields <- names(table)
    view <- saveWorkflowTableSelection(scopes[[1L]], map,
        data.frame(output = fields, source = fields, cast = rep("identity", length(fields))),
        workflowSavedTableMetadata(expected, "id", 3L), "view", identity("view"))
    releaseWorkflowTableReadScope(scopes[[1L]])
    scopes[[2L]] <- newWorkflowTableReadScope(store, view, reader = "arrow")
    original <- streamWorkflowTableRows
    testthat::with_mocked_bindings({
        expect_error(saveWorkflowTablePhysicalExport(scopes[[2L]], "failed-export",
            identity("failed-export")), "interrupted export")
    }, streamWorkflowTableRows = function(scope, columns, consume) {
        original(scope, columns, function(data, offset, row_names) {
            consume(data, offset, row_names)
            stop("interrupted export")
        })
    })
    expect_null(recoverSavedWorkflowRequest(store, "failed-export"))
    expect_identical(digest::digest(file = physical_path, algo = "sha256"), input_hash)
    releaseWorkflowTableReadScope(scopes[[2L]])
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("a selected table and its map survive store relocation", {
    source_root <- withr::local_tempdir()
    destination_root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(source_root, "perf103-relocation", "DIA")
    on.exit(if (!isTRUE(store$closed)) closeWorkflowEntityStore(store), add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-relocation", "design", "selection", selector, "perf103-test")
    table <- data.frame(id = c("A", "B", "C"), quantity = c(1, NA_real_, 3))
    physical <- saveWorkflowTable(store, table, "source", identity("source"),
        "id", metadata_version = 3L)
    map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(3, 1)),
        "map", identity("map"), "logical_ordinal", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    on.exit(if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope), add = TRUE)
    expected <- table[c(3, 1), , drop = FALSE]
    row.names(expected) <- NULL
    fields <- names(table)
    view <- saveWorkflowTableSelection(scope, map,
        data.frame(output = fields, source = fields, cast = c("identity", "identity")),
        workflowSavedTableMetadata(expected, "id", 3L), "view", identity("view"))
    releaseWorkflowTableReadScope(scope)
    closeWorkflowEntityStore(store)
    expect_true(file.copy(file.path(source_root, ".workflow-entities"),
        destination_root, recursive = TRUE, copy.date = TRUE))
    relocated <- newWorkflowEntityStore(destination_root, "perf103-relocation", "DIA")
    on.exit(closeWorkflowEntityStore(relocated), add = TRUE)
    references <- lapply(c("source", "map", "view"),
        function(request) recoverSavedWorkflowRequest(relocated, request))
    expect_identical(references, list(physical, map, view))
    reopened <- newWorkflowTableReadScope(relocated, view, reader = "arrow")
    on.exit(if (!isTRUE(reopened$released)) releaseWorkflowTableReadScope(reopened), add = TRUE)
    expect_identical(readSavedWorkflowTable(relocated, view, fields,
        limit = 2L, scope = reopened)$data, expected)
    releaseWorkflowTableReadScope(reopened)
    expect_identical(workflowEntityCacheInfo(relocated)$open_resources, 0L)
})

test_that("group and order indexes follow each logical population", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-indexes", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-indexes", "design", "selection", selector, "perf103-test")
    table <- data.frame(id = letters[1:4],
        Run = factor(c("S2", "S1", "S2", "S1"),
            levels = c("S1", "unused", "S2")), score = c(4, 1, 3, 2))
    physical <- saveWorkflowTable(store, table, "source", identity("source"),
        "id", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    positions <- list(c(4, 1, 3), c(2, 3, 1))
    views <- vector("list", 2L)
    for (i in seq_len(2L)) {
        map <- saveWorkflowTable(store,
            data.frame(logical_ordinal = as.double(1:3),
                source_ordinal = as.double(positions[[i]])),
            paste0("map-", i), identity(paste0("map-", i)),
            "logical_ordinal", metadata_version = 3L)
        expected <- table[positions[[i]], , drop = FALSE]
        row.names(expected) <- NULL
        views[[i]] <- saveWorkflowTableSelection(scopes[[1L]], map,
            data.frame(output = names(table), source = names(table),
                cast = rep("identity", ncol(table))),
            workflowSavedTableMetadata(expected, "id", 3L),
            paste0("view-", i), identity(paste0("view-", i)))
    }
    releaseWorkflowTableReadScope(scopes[[1L]])
    group_ids <- order_ids <- character(2L)
    for (i in seq_len(2L)) {
        scopes[[i + 1L]] <- newWorkflowTableReadScope(store, views[[i]], reader = "arrow")
        group <- saveWorkflowGroupIndex(scopes[[i + 1L]], "Run")
        again <- saveWorkflowGroupIndex(scopes[[i + 1L]], "Run")
        expect_identical(again$reference, group$reference)
        group_ids[[i]] <- group$reference$entity_id
        order <- saveWorkflowOrderIndex(scopes[[i + 1L]], "score")
        again <- saveWorkflowOrderIndex(scopes[[i + 1L]], "score")
        expect_identical(again$reference, order$reference)
        order_ids[[i]] <- order$reference$entity_id
        grouped <- openWorkflowTableGroups(scopes[[i + 1L]], names(table), "Run")
        parts <- list()
        repeat {
            part <- grouped$next_group()
            if (is.null(part)) break
            parts[[length(parts) + 1L]] <- part
        }
        grouped$close()
        expect_identical(unlist(lapply(parts, `[[`, "id"), use.names = FALSE),
            table$id[positions[[i]]][order(table$Run[positions[[i]]],
                method = "radix")])
        ordered <- openWorkflowOrderedRows(store, views[[i]], c("id", "score"), "score")
        rows <- list()
        repeat {
            part <- ordered$next_batch()
            if (is.null(part)) break
            rows[[length(rows) + 1L]] <- part
        }
        ordered$close()
        expect_identical(unlist(lapply(rows, `[[`, "id"), use.names = FALSE),
            table$id[positions[[i]]][order(table$score[positions[[i]]],
                method = "radix")])
    }
    expect_false(identical(group_ids[[1L]], group_ids[[2L]]))
    expect_false(identical(order_ids[[1L]], order_ids[[2L]]))
    for (scope in scopes[-1L]) releaseWorkflowTableReadScope(scope)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("selected reads reject source mutation and release their cursor", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf103-corruption", "DIA")
    scopes <- list()
    on.exit({
        for (scope in scopes) if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }, add = TRUE)
    identity <- function(selector) newSavedWorkflowScientificIdentity(
        "perf103-corruption", "design", "selection", selector, "perf103-test")
    table <- data.frame(id = c("A", "B", "C"), quantity = c(1, 2, 3))
    physical <- saveWorkflowTable(store, table, "source", identity("source"),
        "id", metadata_version = 3L)
    path <- file.path(store$directory, physical$payload$relative_path)
    original <- readBin(path, "raw", n = file.size(path))
    map <- saveWorkflowTable(store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = c(3, 1)),
        "map", identity("map"), "logical_ordinal", metadata_version = 3L)
    scopes[[1L]] <- newWorkflowTableReadScope(store, physical, reader = "arrow")
    expected <- table[c(3, 1), , drop = FALSE]
    row.names(expected) <- NULL
    view <- saveWorkflowTableSelection(scopes[[1L]], map,
        data.frame(output = names(table), source = names(table),
            cast = c("identity", "identity")),
        workflowSavedTableMetadata(expected, "id", 3L),
        "view", identity("view"))
    releaseWorkflowTableReadScope(scopes[[1L]])
    scopes[[2L]] <- newWorkflowTableReadScope(store, view, reader = "arrow")
    writeBin(as.raw(1L), path)
    expect_error(readSavedWorkflowTable(store, view, names(table),
        limit = 2L, scope = scopes[[2L]]), class = "saved_entity_corrupt")
    expect_error(releaseWorkflowTableReadScope(scopes[[2L]]),
        class = "saved_entity_corrupt")
    expect_true(scopes[[2L]]$released)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
    writeBin(original, path)
    scopes[[3L]] <- newWorkflowTableReadScope(store, view, reader = "arrow")
    expect_identical(readSavedWorkflowTable(store, view, names(table),
        limit = 2L, scope = scopes[[3L]])$data, expected)
    releaseWorkflowTableReadScope(scopes[[3L]])
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("native DIA session archive retains a selected checkpoint table", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "perf103-session", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$setWorkflowType("DIA")
    original <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw_data_s4", original, list(), "raw")
    convertDIANativeLegacyCheckpoints(manager)
    checkpoint <- readDIAStateMetadata(manager)$data
    physical <- checkpoint$components$peptide_data
    scope <- newWorkflowTableReadScope(manager$.context$store, physical, reader = "arrow")
    on.exit(if (!isTRUE(scope$released)) releaseWorkflowTableReadScope(scope), add = TRUE)
    at <- c(3L, 1L)
    expected <- original@peptide_data[at, , drop = FALSE]
    map <- saveWorkflowTable(manager$.context$store,
        data.frame(logical_ordinal = as.double(1:2), source_ordinal = as.double(at)),
        "session-map", newSavedWorkflowScientificIdentity(physical$entity_id,
            "session", "row-map", "perf103", "perf103-test"),
        "logical_ordinal", metadata_version = 3L)
    fields <- scope$metadata$columns
    view <- saveWorkflowTableSelection(scope, map,
        data.frame(output = fields, source = fields, cast = rep("identity", length(fields))),
        workflowSavedTableMetadata(expected, scope$metadata$key_column, 3L),
        "session-view", newSavedWorkflowScientificIdentity(physical$entity_id,
            "session", "view", "perf103", "perf103-test"))
    releaseWorkflowTableReadScope(scope)
    checkpoint$components$peptide_data <- view
    checkpoint$dimensions <- list(rows = 2L, columns = fields)
    state <- commitDIAReferenceState(manager, "selected", checkpoint, list(),
        "selected checkpoint")
    acknowledgeDIAOutputs(manager, state, "correlation")
    published <- runDIANativeSessionExport(manager, list(workflow_type = "DIA"),
        withr::local_tempdir())
    imported <- openDIANativeSession(published$publication$archive_path,
        withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(imported$manager))
    restored <- readDIAStateMetadata(imported$manager)$data$components$peptide_data
    expect_identical(restored, view)
    saved <- imported$manager$.context$store
    expect_true(all(c(physical$entity_id, map$entity_id) %in%
        ls(saved$records, all.names = TRUE)))
    reopened <- newWorkflowTableReadScope(saved, restored, reader = "arrow")
    on.exit(if (!isTRUE(reopened$released)) releaseWorkflowTableReadScope(reopened), add = TRUE)
    expect_identical(readSavedWorkflowTable(saved, restored, fields,
        limit = 2L, scope = reopened)$data, expected)
    releaseWorkflowTableReadScope(reopened)
    expect_identical(workflowEntityCacheInfo(saved)$open_resources, 0L)
})
