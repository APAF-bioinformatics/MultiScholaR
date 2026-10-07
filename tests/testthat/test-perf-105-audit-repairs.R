test_that("metadata admission rejects wide payloads before R transfer and cleans up", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf105-admission", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, data.frame(id = "a"), "input",
        newSavedWorkflowScientificIdentity("input", "design", "metadata", "all", "perf105"),
        "id", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
    withr::defer(releaseWorkflowTableReadScope(scope))
    configureWorkflowBulkQueryScope(scope)
    before <- DBI::dbListTables(scope$connection)
    expect_identical(collectWorkflowBulkMetadata(scope, "SELECT 1 AS id")$id, 1L)
    expect_error(collectWorkflowBulkMetadata(scope,
        "SELECT repeat('x', 1048576) AS payload", max_rows = 1L, max_bytes = 8192),
        class = "saved_entity_limit")
    expect_identical(DBI::dbListTables(scope$connection), before)
    if (isTRUE(capabilities("profmem"))) {
        path <- withr::local_tempfile()
        Rprofmem(path, threshold = 524288)
        tryCatch(collectWorkflowBulkMetadata(scope,
            "SELECT repeat('x', 1048576) AS payload", max_rows = 1L, max_bytes = 64),
            saved_entity_limit = \(error) NULL, finally = Rprofmem(NULL))
        allocations <- suppressWarnings(as.numeric(sub(" .*", "", readLines(path))))
        expect_lt(max(c(0, allocations), na.rm = TRUE), 1048576)
    }
    expect_error(collectWorkflowBulkMetadata(scope, "SELECT [1,2] AS nested"),
        class = "saved_entity_unsupported")
    expect_error(collectWorkflowBulkMetadata(scope,
        "SELECT range AS id FROM range(5)", max_rows = 4L), class = "saved_entity_limit")
    expect_error(collectWorkflowBulkMetadata(scope, "SELECT absent_column"))
    expect_identical(DBI::dbListTables(scope$connection), before)
    releaseWorkflowTableReadScope(scope)
    expect_identical(store$open_resources, 0L)
})

test_that("bounded metadata executes once and preserves order types and parameters", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf105-order", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    reference <- saveWorkflowTable(store, data.frame(id = "a"), "input",
        newSavedWorkflowScientificIdentity("input", "design", "order", "all", "perf105"),
        "id", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
    withr::defer(releaseWorkflowTableReadScope(scope))
    configureWorkflowBulkQueryScope(scope)
    con <- scope$connection
    DBI::dbExecute(con, "CREATE TEMP SEQUENCE metadata_once START 1")
    withr::defer(DBI::dbExecute(con, "DROP SEQUENCE metadata_once"))
    value <- collectWorkflowBulkMetadata(scope, paste(
        "SELECT CAST(i AS INTEGER) AS id, nextval('metadata_once') AS call,",
        "CAST(NULL AS VARCHAR) AS missing, ? AS label",
        "FROM range(4) t(i) ORDER BY id DESC"), params = list("\u03b2"))
    expect_identical(value$id, 3:0)
    expect_identical(value$missing, rep(NA_character_, 4L))
    expect_identical(value$label, rep("\u03b2", 4L))
    expect_equal(as.numeric(DBI::dbGetQuery(con,
        "SELECT nextval('metadata_once') AS n")$n), 5)
    empty <- collectWorkflowBulkMetadata(scope,
        "SELECT CAST(NULL AS INTEGER) AS id WHERE FALSE", max_rows = 0L)
    expect_identical(empty, data.frame(id = integer()))
    expect_false(any(startsWith(DBI::dbListTables(con), "workflow_metadata_")))
})

test_that("complete-group partition rejects overflow without losing valid groups", {
    groups <- c("A", "B", "C")
    expect_error(partitionWorkflowCompleteGroups(groups, .Machine$integer.max + 1),
        class = "saved_entity_limit")
    expect_error(partitionWorkflowCompleteGroups(matrix(groups), 3),
        class = "saved_entity_limit")
    expect_identical(unlist(partitionWorkflowCompleteGroups(groups,
        .Machine$integer.max), use.names = FALSE), groups)
    expect_identical(unlist(partitionWorkflowCompleteGroups(groups, 3L),
        use.names = FALSE), groups)
})
