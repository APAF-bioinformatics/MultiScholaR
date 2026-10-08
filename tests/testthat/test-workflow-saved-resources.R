testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"),
    "durable saved-entity backend is Linux-qualified only")

resourceEntityFixture <- function(store, value, request_id) {
    saveWorkflowCheckpoint(store, value, request_id,
        identity = list(input = "resource-fixture-v1", design = "paired-design-v1",
            parameters = "none", selectors = request_id, software = "test-build-v1"),
        metadata = list(class = class(value), shape = as.integer(dim(value)), factors = list()))
}

test_that("completed output prefixes release cache without changing file bytes", {
    path <- withr::local_tempfile()
    connection <- file(path, "w+b")
    seek(connection, 67108864 + 4096, origin = "start", rw = "write")
    writeBin(as.raw(17L), connection)
    close(connection)
    before <- digest::digest(file = path, algo = "sha256")
    owner <- list2env(list(path = path, advised_bytes = 0), parent = emptyenv())
    releaseWorkflowWrittenPrefix(owner)
    expect_equal(owner$advised_bytes, 67108864)
    releaseWorkflowWrittenPrefix(owner, final = TRUE)
    expect_equal(owner$advised_bytes, file.size(path))
    expect_identical(digest::digest(file = path, algo = "sha256"), before)
    expect_error(releaseWorkflowFilePages(path, length = file.size(path) + 1),
        class = "saved_entity_write")
    expect_error(releaseWorkflowFilePages(dirname(path)), class = "saved_entity_write")
})

test_that("import cache ownership and reclamation wait for durable idle handoff", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "import-cache", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$.context$reclaim_enabled <- TRUE
    source <- withr::local_tempfile()
    writeLines("immutable report", source)
    reference <- saveDIAOwnedValue(data.frame(Run = "A", value = 1), "data_tbl", manager)
    registerDIAInputFiles(manager, source)
    completeDIAInputFiles(manager)
    expect_identical(manager$.context$reclaim_pending$label, "import")
    expect_identical(manager$.context$reclaim_pending$state, reference$entity_id)
    releases <- character()
    testthat::local_mocked_bindings(releaseWorkflowFilePages = function(path, ...) {
        releases <<- c(releases, path)
        TRUE
    }, reclaimWorkflowMemory = function(...) list(status = "released"))
    manager$.context$scope_depth <- 1L
    releaseDIAInputFiles(manager)
    expect_length(releases, 0L)
    expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    manager$.context$scope_depth <- 0L
    expect_identical(reclaimDIACompletedMemory(manager)$status, "released")
    releaseDIAInputFiles(manager)
    releaseDIAInputFiles(manager)
    expect_identical(releases, normalizePath(source))
    expect_identical(readLines(source), "immutable report")
    registerDIAInputFiles(manager, source)
    writeLines("changed original report", source)
    expect_error(completeDIAInputFiles(manager), class = "saved_entity_corrupt")
})

test_that("allocator release reports capability and never changes live scientific data", {
    live <- list(matrix = matrix(c(NA_real_, NaN, Inf, -Inf, 0, 1), 2L),
        labels = factor(c("b", "a"), levels = c("unused", "a", "b")))
    before <- serialize(live, NULL)
    support <- workflowAllocatorRelease(FALSE)
    expect_true(identical(support, TRUE) || identical(support, NA))
    expect_true(reclaimWorkflowMemory()$status %in%
        c("released", "no_pages_released", "unsupported"))
    expect_identical(serialize(live, NULL), before)
    expect_error(workflowAllocatorRelease(NA), "non-missing logical")
    expect_error(workflowAllocatorRelease(1), "non-missing logical")
    expect_identical(reclaimWorkflowMemory(function(release) NA,
        function(...) stop("unsupported platform must not collect"))$status, "unsupported")
    calls <- character()
    native <- function(release) { calls <<- c(calls, if (release) "release" else "probe"); !release }
    collect <- function(full) { expect_true(full); calls <<- c(calls, "gc") }
    expect_identical(reclaimWorkflowMemory(native, collect)$status, "no_pages_released")
    expect_identical(calls, c("probe", "gc", "release"))
    expect_error(reclaimWorkflowMemory(function(...) stop("native failure")),
        class = "workflow_reclaim_error")
    dll <- getLoadedDLLs()[["MultiScholaR"]]
    expect_false(dll[["dynamicLookup"]])
    expect_identical(getDLLRegisteredRoutines(dll)$.Call$workflow_reclaim_pages$numParameters, 1L)
})

test_that("native library identity detects changed bytes independently of R code", {
    root <- withr::local_tempdir()
    expect_length(workflowNativeLibraryIdentity(root), 0L)
    dir.create(file.path(root, "libs"))
    path <- file.path(root, "libs", "MultiScholaR.so")
    writeBin(as.raw(1:3), path)
    first <- workflowNativeLibraryIdentity(root)
    writeBin(as.raw(3:1), path)
    expect_named(first, "MultiScholaR.so")
    expect_false(identical(first, workflowNativeLibraryIdentity(root)))
    expect_identical(first[[1L]], digest::digest(as.raw(1:3), serialize = FALSE, algo = "sha256"))
})

test_that("zero and oversized caches decode without retaining payloads", {
    root <- withr::local_tempdir()
    value <- list(signal = rep(1L, 10000L), missing = c(NA, NaN, Inf, -Inf))
    zero <- newWorkflowEntityStore(root, "project", "DIA")
    reference <- resourceEntityFixture(zero, value, "zero")
    expect_identical(hydrateSavedWorkflowEntity(zero, reference), value)
    expect_identical(hydrateSavedWorkflowEntity(zero, reference), value)
    info <- workflowEntityCacheInfo(zero)
    expect_identical(info$entries, 0L)
    expect_identical(info$payload_reads, 2L)
    expect_identical(info$payload_decodes, 2L)
    expect_identical(info$cache_hits, 0L)
    expect_identical(info$open_resources, 0L)
    closeWorkflowEntityStore(zero)

    small <- newWorkflowEntityStore(root, "project", "DIA", max_cache_bytes = 1024)
    expect_identical(hydrateSavedWorkflowEntity(small, reference), value)
    expect_identical(hydrateSavedWorkflowEntity(small, reference), value)
    expect_identical(workflowEntityCacheInfo(small)$entries, 0L)
    expect_identical(workflowEntityCacheInfo(small)$payload_reads, 2L)
    expect_error(newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA",
        max_cache_bytes = 67108865), class = "saved_entity_invalid")
})

test_that("bounded resident entries hit, evict and reload without removing records", {
    value_one <- list(signal = rep(1L, 1024L))
    value_two <- list(signal = rep(2L, 1024L))
    one_entry <- as.numeric(object.size(value_one)) + 128
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "project", "DIA", max_cache_bytes = one_entry)
    first <- resourceEntityFixture(store, value_one, "first")
    second <- resourceEntityFixture(store, value_two, "second")
    expect_identical(hydrateSavedWorkflowEntity(store, first), value_one)
    expect_identical(hydrateSavedWorkflowEntity(store, first), value_one)
    expect_identical(workflowEntityCacheInfo(store)$cache_hits, 1L)
    expect_identical(hydrateSavedWorkflowEntity(store, second), value_two)
    expect_identical(workflowEntityCacheInfo(store)$entries, 1L)
    expect_identical(hydrateSavedWorkflowEntity(store, first), value_one)
    info <- workflowEntityCacheInfo(store)
    expect_lte(info$bytes, one_entry)
    expect_lte(info$bytes, 67108864)
    expect_identical(info$payload_reads, 3L)
    expect_identical(info$payload_decodes, 3L)
    expect_identical(info$cache_hits, 1L)
    expect_identical(info$top_level_aliases, 0L)
    expect_identical(info$unique_top_level_bytes, info$bytes)
    alias <- hydrateSavedWorkflowEntity(store, first)
    held_elsewhere <- alias
    expect_identical(rlang::obj_address(alias), rlang::obj_address(held_elsewhere))
    releaseWorkflowEntityCache(store)
    expect_identical(workflowEntityCacheInfo(store)$entries, 0L)
    expect_identical(held_elsewhere, value_one)
    expect_identical(hydrateSavedWorkflowEntity(store, first), value_one)
    expect_identical(workflowEntityCacheInfo(store)$payload_reads, 4L)
    expect_identical(describeSavedWorkflowEntity(store, second$entity_id), second)
    path <- file.path(store$directory, first$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
    writeBin(bytes, path)
    expect_error(hydrateSavedWorkflowEntity(store, first), class = "saved_entity_corrupt")
})

test_that("mutable environment owners bypass resident reuse", {
    context <- new.env(parent = emptyenv())
    context$value <- 1L
    value <- list(context = context)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA",
        max_cache_bytes = 67108864)
    reference <- resourceEntityFixture(store, value, "mutable")
    first <- hydrateSavedWorkflowEntity(store, reference)
    first$context$value <- 2L
    second <- hydrateSavedWorkflowEntity(store, reference)
    expect_identical(second$context$value, 1L)
    info <- workflowEntityCacheInfo(store)
    expect_identical(info$entries, 0L)
    expect_identical(info$payload_reads, 2L)
    expect_identical(info$cache_hits, 0L)
})

test_that("accounting deduplicates identical top-level owners but keeps a safe bound", {
    value <- list(signal = rep(1L, 1024L))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA",
        max_cache_bytes = 67108864)
    first <- resourceEntityFixture(store, value, "alias-one")
    second <- resourceEntityFixture(store, value, "alias-two")
    cacheWorkflowEntityPayload(store, first$entity_id, value)
    cacheWorkflowEntityPayload(store, second$entity_id, value)
    info <- workflowEntityCacheInfo(store)
    expect_identical(info$entries, 2L)
    expect_identical(info$top_level_aliases, 1L)
    expect_identical(info$unique_top_level_bytes * 2, info$bytes)
    expect_lte(info$bytes, info$max_bytes)
    releaseWorkflowEntityCache(store)
    expect_identical(value, list(signal = rep(1L, 1024L)))
})

test_that("ten close/reopen rounds retain committed and pending files under both ceilings", {
    for (cache_limit in c(0, 67108864)) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "project", "DIA", max_cache_bytes = cache_limit)
        pending <- file.path(store$directory, ".pending")
        dir.create(pending, showWarnings = FALSE)
        sentinel <- file.path(pending, "unpublished-input")
        writeLines("recoverable", sentinel)
        for (round in seq_len(10L)) {
            request_id <- sprintf("round%02d", round)
            value <- list(round = round, group = ordered(c("control", "treated")),
                intensity = c(NA_real_, NaN, Inf, -Inf, round))
            reference <- resourceEntityFixture(store, value, request_id)
            expect_identical(hydrateSavedWorkflowEntity(store, reference), value)
            info <- workflowEntityCacheInfo(store)
            expect_lte(info$bytes, cache_limit)
            expect_identical(info$open_resources, 0L)
            releaseWorkflowEntityCache(store)
            expect_identical(workflowEntityCacheInfo(store)$entries, 0L)
            expect_true(file.exists(sentinel))
            closeWorkflowEntityStore(store)
            closeWorkflowEntityStore(store)
            expect_error(hydrateSavedWorkflowEntity(store, reference), class = "saved_entity_scope")
            store <- newWorkflowEntityStore(root, "project", "DIA", max_cache_bytes = cache_limit)
            recovered <- recoverSavedWorkflowRequest(store, request_id)
            expect_identical(recovered, reference)
            expect_identical(hydrateSavedWorkflowEntity(store, recovered), value)
            expect_lte(workflowEntityCacheInfo(store)$bytes, cache_limit)
            expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
        }
        expect_true(file.exists(sentinel))
        expect_identical(readLines(sentinel), "recoverable")
        closeWorkflowEntityStore(store)
    }
})

test_that("separate project stores cannot share cache or references", {
    value <- list(group = factor(c("A", "B")), signal = c(1, NA_real_))
    first <- newWorkflowEntityStore(withr::local_tempdir(), "project-a", "DIA",
        max_cache_bytes = 67108864)
    second <- newWorkflowEntityStore(withr::local_tempdir(), "project-b", "DIA",
        max_cache_bytes = 67108864)
    reference <- resourceEntityFixture(first, value, "isolated")
    expect_identical(hydrateSavedWorkflowEntity(first, reference), value)
    expect_identical(workflowEntityCacheInfo(first)$entries, 1L)
    expect_identical(workflowEntityCacheInfo(second)$entries, 0L)
    expect_error(hydrateSavedWorkflowEntity(second, reference), class = "saved_entity_scope")
    closeWorkflowEntityStore(first)
    expect_identical(workflowEntityCacheInfo(second)$open_resources, 0L)
    closeWorkflowEntityStore(second)
})

test_that("progress receipts stay bounded, throttle, and do not imply scientific completion", {
    root <- withr::local_tempdir()
    observer <- newWorkflowProgressObserver(root, function() "raw_data_s4", interval = 3600)
    path <- file.path(root, "gui-progress-latest.json")
    observer("table_read", 10, 100, "source")
    first <- readLines(path)
    observer("table_read", 20, 100, "source")
    expect_identical(readLines(path), first)
    observer("table_read", 100, 100, "source")
    value <- jsonlite::read_json(path)
    expect_true(value$phase_complete)
    expect_identical(value$current_state, "raw_data_s4")
    expect_equal(value$completed, 100)
    observer("input_groups", 1, 300, "source", "groups:Run")
    value <- jsonlite::read_json(path)
    expect_false(value$phase_complete)
    expect_identical(value$unit, "groups:Run")
    expect_lt(file.size(path), 8192)
    expect_length(list.files(root, all.files = TRUE, pattern = "^[.]progress-"), 0L)
    expect_error(observer("bad", data.frame(x = 1), 2, "source"))
})

test_that("streamed writer progress does not treat prototype rows as output totals", {
    skip_if_not_installed("arrow")
    root <- withr::local_tempdir()
    observer <- newWorkflowProgressObserver(root, function() "raw_data_s4", interval = 0)
    metadata <- workflowSavedTableMetadata(data.frame(id = character()), "id", 3L)
    writer <- openWorkflowTableWriter(file.path(root, "payload.parquet"), metadata, observer)
    on.exit(abortWorkflowTableWriter(writer), add = TRUE)
    appendWorkflowTableWriter(writer, data.frame(id = c("a", "b", "c")))
    result <- finishWorkflowTableWriter(writer)
    expect_identical(result$rows, 3L)
    progress <- jsonlite::read_json(file.path(root, "gui-progress-latest.json"))
    expect_equal(progress$completed, 3)
    expect_null(progress$total)
    expect_false(progress$phase_complete)
    expect_null(writer$progress)
})

test_that("a failing cursor close still releases the scope, its action and connection", {
    table <- data.frame(id = 1:3, value = c(1, 2, 3))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "release-fault", "DIA")
    reference <- saveWorkflowTable(store, table, "release-fault",
        newSavedWorkflowScientificIdentity(digest::digest(table, algo = "sha256"),
            "release-fault-design-v1", "none", "release-fault", "MultiScholaR-release-test-v1"),
        ".rownames")
    scope <- newWorkflowTableReadScope(store, reference, "duckdb")
    action <- scope$action
    expect_false(is.null(action$connection))
    closed <- character()
    assign("a_failing", list(close = \() stop("cursor close failed")), envir = scope$readers)
    assign("b_healthy", list(close = \() closed <<- c(closed, "b_healthy")), envir = scope$readers)
    expect_error(releaseWorkflowTableReadScope(scope), "cursor close failed")
    expect_true(scope$released)
    expect_identical(closed, "b_healthy")
    expect_identical(store$open_resources, 0L)
    expect_null(store$table_action)
    expect_null(action$connection)
    expect_identical(action$scopes, 0L)
    expect_invisible(releaseWorkflowTableReadScope(scope))
})
