testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))

#' Own the small saved tables and leases used by the six audit repairs.
#' @param table Native oracle table.
#' @param envir Test cleanup owner.
#' @return Fixture with a store, scoped reader and selection constructor.
#' @noRd
perf103RepairFixture <- function(table = data.frame(id = letters[1:4],
    quantity = as.double(1:4), group = factor(c("a", "b", "a", "b"),
        levels = c("a", "unused", "b"))), envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = envir)
    store <- newWorkflowEntityStore(root, "perf103-repairs", "DIA")
    owners <- list()
    withr::defer({
        for (owner in rev(owners)) if (!isTRUE(owner$released)) releaseWorkflowTableReadScope(owner)
        closeWorkflowEntityStore(store)
    }, envir = envir)
    identity <- function(label) newSavedWorkflowScientificIdentity(
        "audit-repairs", "design", label, label, "perf103-repairs")
    source <- saveWorkflowTable(store, table, "source", identity("source"),
        "id", metadata_version = 3L)
    open <- function(reference = source, reader = "arrow") {
        owner <- newWorkflowTableReadScope(store, reference, reader)
        owners[[length(owners) + 1L]] <<- owner
        owner
    }
    select <- function(owner, positions, expected, label, mapping = NULL) {
        map <- if (is.null(positions)) NULL else saveWorkflowTable(store,
            data.frame(logical_ordinal = as.double(seq_along(positions)),
                source_ordinal = as.double(positions)), paste0(label, "-map"),
            identity(paste0(label, "-map")), "logical_ordinal", metadata_version = 3L)
        if (is.null(mapping)) mapping <- data.frame(output = names(expected),
            source = names(expected), cast = rep("identity", ncol(expected)))
        saveWorkflowTableSelection(owner, map, mapping,
            workflowSavedTableMetadata(expected, if ("id" %in% names(expected)) "id" else ".rownames", 3L),
            label, identity(label))
    }
    list(root = root, store = store, source = source, table = table, open = open,
        select = select, identity = identity)
}

test_that("selection pages match physical row-name storage at every offset", {
    for (kind in c("automatic", "integer", "character")) {
        table <- data.frame(id = letters[1:4], quantity = as.double(1:4))
        if (kind == "integer") attr(table, "row.names") <- c(11L, 22L, 33L, 44L)
        if (kind == "character") rownames(table) <- paste0("row-", letters[1:4])
        f <- perf103RepairFixture(table)
        expected <- table[c(4L, 2L, 1L), , drop = FALSE]
        if (kind == "automatic") rownames(expected) <- NULL
        view <- f$select(f$open(), c(4L, 2L, 1L), expected, "page-view")
        physical <- saveWorkflowTable(f$store, expected, "page-physical",
            f$identity("page-physical"), "id", metadata_version = 3L)
        for (reader in c("arrow", "duckdb", "fastduckdb")) {
            selected_scope <- f$open(view, reader)
            physical_scope <- f$open(physical,
                if (reader == "arrow") "fastduckdb" else reader)
            for (offset in c(0L, 1L, 2L, 3L)) {
                selected <- readSavedWorkflowTable(f$store, view, names(table),
                    offset = offset, limit = 2L, scope = selected_scope)
                control <- readSavedWorkflowTable(f$store, physical, names(table),
                    offset = offset, limit = 2L, scope = physical_scope)
                expect_identical(selected$data, control$data)
                expect_identical(digest::digest(selected$data), digest::digest(control$data))
            }
            expect_identical(materializeDIATable(list(store = f$store), view,
                scope = selected_scope), expected)
            releaseWorkflowTableReadScope(selected_scope)
            releaseWorkflowTableReadScope(physical_scope)
        }
    }
})

test_that("v2 metadata edits preserve reordered, filtered and projected populations", {
    f <- perf103RepairFixture()
    source <- f$open()
    for (i in seq_len(2L)) {
        positions <- if (i == 1L) c(4L, 2L, 1L, 3L) else c(4L, 2L)
        expected <- f$table[positions, c("id", "quantity"), drop = FALSE]
        row.names(expected) <- NULL
        view <- f$select(source, positions, expected, paste0("selection", i))
        scope <- f$open(view)
        metadata <- scope$metadata
        metadata$attributes$note <- "logical edit"
        metadata$attribute_order <- c(metadata$attribute_order, "note")
        changed <- saveWorkflowTableComponent(scope, metadata, paste0("edited", i),
            f$identity(paste0("edited", i)))
        actual <- materializeDIATable(list(store = f$store), changed)
        attr(expected, "note") <- "logical edit"
        expect_identical(actual, expected)
        component <- readWorkflowTableComponent(f$store, changed)
        expect_identical(component$view$rows, scope$view$rows)
        expect_identical(component$source, f$source)
        releaseWorkflowTableReadScope(scope)
    }
    releaseWorkflowTableReadScope(source)
    expect_identical(f$store$open_resources, 0L)
})

test_that("the DIA design metadata caller retains the v2 selection", {
    table <- prepareDIAParserTable(readr::read_tsv(I("id\tquantity\na\t1\nb\t2\nc\t3\nd\t4\n"),
        show_col_types = FALSE))
    f <- perf103RepairFixture(table)
    source <- f$open()
    expected <- table[c(4L, 2L, 1L, 3L), , drop = FALSE]
    row.names(expected) <- NULL
    view <- f$select(source, c(4L, 2L, 1L, 3L), expected, "design-selection")
    result <- transformDIADesignTable(list(.context = list(store = f$store,
        software = "audit-repair")), view, strip_parser = TRUE)
    class(expected) <- c("tbl_df", "tbl", "data.frame")
    attr(expected, "spec") <- attr(expected, "problems") <- NULL
    expect_identical(materializeDIATable(list(store = f$store), result), expected)
    expect_identical(readWorkflowTableComponent(f$store, result)$view$rows,
        readWorkflowTableComponent(f$store, view)$view$rows)
})

test_that("v1 physical export preserves projection, factors and table attributes", {
    f <- perf103RepairFixture()
    scope <- f$open()
    metadata <- scope$metadata
    metadata$columns <- c("quantity", "group")
    metadata$column_metadata <- metadata$column_metadata[metadata$columns]
    metadata$attributes$note <- "legacy logical output"
    metadata$attribute_order <- c(metadata$attribute_order, "note")
    view <- saveWorkflowTableComponent(scope, metadata, "legacy-overlay", f$identity("overlay"))
    logical <- f$open(view)
    exported <- saveWorkflowTablePhysicalExport(logical, "legacy-export", f$identity("export"))
    expected <- f$table[c("quantity", "group")]
    attr(expected, "note") <- "legacy logical output"
    expect_identical(materializeDIATable(list(store = f$store), exported), expected)
    expect_identical(exported$payload$codec, "parquet_duckdb")
    expect_false(identical(exported, f$source))
    metadata$attributes$note <- "different overlay"
    second <- saveWorkflowTableComponent(scope, metadata, "legacy-overlay-two", f$identity("overlay-two"))
    second_export <- saveWorkflowTablePhysicalExport(f$open(second), "legacy-export-two", f$identity("export"))
    expect_false(identical(exported$entity_id, second_export$entity_id))
    expect_identical(attr(materializeDIATable(list(store = f$store), second_export), "note"),
        "different overlay")
    old_library <- Sys.getenv("PERF103_OLD_LIB")
    if (nzchar(old_library)) {
        legacy <- callr::r(function(root, reference, old_library) {
            ns <- loadNamespace("MultiScholaR", lib.loc = old_library)
            store <- get("newWorkflowEntityStore", ns)(root, "perf103-repairs", "DIA")
            on.exit(get("closeWorkflowEntityStore", ns)(store))
            get("materializeDIATable", ns)(list(store = store), reference)
        }, args = list(f$root, exported, old_library), libpath = c(old_library, .libPaths()))
        expect_identical(legacy, expected)
    }
})

test_that("unrepresentable parent row names refuse and physical-parent export restores native behavior", {
    f <- perf103RepairFixture()
    source <- f$open()
    parent <- f$table[c(4L, 2L, 1L, 3L), ]
    row.names(parent) <- NULL
    view <- f$select(source, c(4L, 2L, 1L, 3L), parent, "parent")
    parent_scope <- f$open(view)
    expected <- parent[c(2L, 1L), ]
    expect_error(f$select(parent_scope, c(2L, 1L), expected, "unrepresentable"),
        class = "saved_entity_unsupported")
    exported <- saveWorkflowTablePhysicalExport(parent_scope, "parent-export", f$identity("export"))
    child <- f$select(f$open(exported), c(2L, 1L), expected, "supported-child")
    expect_identical(materializeDIATable(list(store = f$store), child), expected)
    for (kind in c("integer", "character")) {
        table <- f$table
        attr(table, "row.names") <- if (kind == "integer") c(11L, 22L, 33L, 44L) else
            c("row-a", "row-b", "row-c", "row-d")
        nested <- perf103RepairFixture(table)
        parent <- table[c(4L, 2L, 1L, 3L), ]
        parent_ref <- nested$select(nested$open(), c(4L, 2L, 1L, 3L), parent, "explicit-parent")
        expected <- parent[c(2L, 1L), ]
        child <- nested$select(nested$open(parent_ref), c(2L, 1L), expected, "explicit-child")
        expect_identical(materializeDIATable(list(store = nested$store), child), expected)
    }
})

test_that("v2 fast reads use fastduckdb with truthful counters and reusable scopes", {
    f <- perf103RepairFixture()
    source <- f$open()
    original <- fastduckdb::read_parquet_fast
    for (i in seq_len(3L)) {
        positions <- list(NULL, c(2L, 3L), c(4L, 2L))[[i]]
        expected <- if (is.null(positions)) f$table[c("id", "quantity")] else
            f$table[positions, c("id", "quantity")]
        row.names(expected) <- NULL
        view <- f$select(source, positions, expected, paste0("fast", i))
        scope <- f$open(view, "fastduckdb")
        calls <- 0L
        threads <- integer()
        before <- f$store$io_counts$fast_reads
        testthat::with_mocked_bindings({
            testthat::with_mocked_bindings({
                page <- readSavedWorkflowTable(f$store, view, names(expected),
                    limit = 4L, scope = scope)
                expect_identical(page$data, workflowPhysicalPageOracle(expected))
                child <- scope$view_scopes$source
                readSavedWorkflowTable(f$store, view, names(expected), limit = 4L, scope = scope)
                expect_identical(scope$view_scopes$source, child)
            }, dbGetQuery = function(...) stop("unexpected DuckDB query"), .package = "DBI")
        }, read_parquet_fast = function(..., n_threads) {
            calls <<- calls + 1L
            threads <<- c(threads, n_threads)
            original(..., n_threads = n_threads)
        }, .package = "fastduckdb")
        expect_identical(calls, c(2L, 4L, 6L)[[i]])
        expect_identical(f$store$io_counts$fast_reads - before, calls)
        expect_true(all(threads == 8L))
        releaseWorkflowTableReadScope(scope)
    }
    releaseWorkflowTableReadScope(source)
    expect_identical(f$store$open_resources, 0L)
})

test_that("wide v2 output counts actual fast backend calls", {
    table <- data.frame(id = letters[1:4])
    for (i in seq_len(65L)) table[[paste0("v", i)]] <- as.double(1:4 + i)
    f <- perf103RepairFixture(table)
    source <- f$open()
    expected <- table[c(4L, 2L), ]
    row.names(expected) <- NULL
    view <- f$select(source, c(4L, 2L), expected, "wide")
    scope <- f$open(view)
    before <- f$store$io_counts$fast_reads
    rows <- NULL
    streamWorkflowTableRows(scope, names(expected), function(data, offset, row_names) rows <<- data)
    expect_identical(rows, expected)
    expect_identical(f$store$io_counts$fast_reads - before, 6L)
})

test_that("failed fast v2 reads release their owned physical leases", {
    f <- perf103RepairFixture()
    source <- f$open()
    expected <- f$table[c(4L, 2L), ]
    row.names(expected) <- NULL
    view <- f$select(source, c(4L, 2L), expected, "fast-failure")
    releaseWorkflowTableReadScope(source)
    testthat::with_mocked_bindings({
        expect_error(readSavedWorkflowTable(f$store, view, "id", reader = "fastduckdb"),
            "injected fast read failure")
    }, read_parquet_fast = function(...) stop("injected fast read failure"),
        .package = "fastduckdb")
    expect_identical(f$store$open_resources, 0L)
})

test_that("S4 string accumulation refuses before retaining all batches", {
    n <- 30000L
    table <- data.frame(id = seq_len(n), text = paste0(sprintf("row%05d", seq_len(n)),
        strrep("x", 80L)))
    f <- perf103RepairFixture(table)
    expected <- table
    names(expected)[[2L]] <- "label"
    view <- f$select(f$open(), NULL, expected, "strings",
        data.frame(output = names(expected), source = names(table), cast = c("identity", "identity")))
    original <- readWorkflowViewPage
    offsets <- integer()
    testthat::with_mocked_bindings({
        expect_error(materializeDIATable(list(store = f$store), view, max_bytes = 2500000),
            class = "saved_entity_limit")
    }, readWorkflowViewPage = function(scope, columns, offset, limit, max_bytes) {
        offsets <<- c(offsets, offset)
        original(scope, columns, offset, limit, max_bytes)
    })
    expect_identical(as.integer(offsets), c(0L, 10000L))
    expect_identical(f$store$open_resources, 1L)
    validateSavedWorkflowReference(f$store, f$source, payload = TRUE)
})
