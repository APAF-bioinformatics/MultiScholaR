test_that("persistent scope transitions remove closed cursors before teardown", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "C")
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "cursor-transitions", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    data <- data.frame(value = as.double(seq_len(20005L)),
        label = factor(rep(c("a", "\u03b2"), length.out = 20005L), levels = c("a", "\u03b2", "unused")))
    identity <- newSavedWorkflowScientificIdentity("cursor", "fixture", "unchanged", "all", "reader-repair")
    reference <- saveWorkflowTable(store, data, "input", identity, ".rownames", metadata_version = 3L)
    for (transition in c("projection", "random", "restart", "disable")) {
        scope <- newWorkflowTableReadScope(store, reference, reader = "fastduckdb")
        scope$read_threads <- 8L
        first <- readSavedWorkflowTable(store, reference, names(data), limit = 10000L, scope = scope)
        expect_identical(first$data$value, data$value[seq_len(10000L)])
        old <- scope$readers$persistent_native$handle
        columns <- if (transition %in% c("projection", "restart")) "value" else names(data)
        offset <- switch(transition, random = 7L, restart = 0L, 10000L)
        page <- withr::with_options(list(multischolar.perf108.persistent_reader = transition != "disable"),
            readSavedWorkflowTable(store, reference, columns, offset = offset, limit = 10000L, scope = scope))
        expect_identical(page$data$value, data$value[offset + seq_len(10000L)])
        expect_true(fastduckdb::parquet_cursor_info_fast(old)$closed)
        expect_false(fastduckdb::parquet_cursor_info_fast(old)$retained_file)
        expect_identical(exists("persistent_native", scope$readers, inherits = FALSE), transition == "restart")
        expect_no_error(releaseWorkflowTableReadScope(scope))
        expect_identical(store$open_resources, 0L)
    }
})

test_that("persistent staged readers close when their consumer fails", {
    skip_if(Sys.getenv("MULTISCHOLAR_PERF108_STRUCTURAL") != "C")
    path <- tempfile(fileext = ".parquet")
    withr::defer(unlink(path))
    arrow::write_parquet(data.frame(.dia_output_order = seq_len(20005L), value = seq_len(20005L)), path)
    closed <- NULL
    original <- fastduckdb::close_parquet_cursor_fast
    testthat::local_mocked_bindings(close_parquet_cursor_fast = function(cursor) {
        original(cursor)
        closed <<- fastduckdb::parquet_cursor_info_fast(cursor)
    }, .package = "fastduckdb")
    expect_error(streamDIARelationalResultFile(NULL, path, function(part) stop("owned consumer failure"),
        fast = TRUE, rows = 20005L), "owned consumer failure")
    expect_true(closed$closed)
    expect_false(closed$retained_file)
    expect_false(closed$retained_metadata)
})
