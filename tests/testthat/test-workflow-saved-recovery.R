testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"), "durable saved-entity backend is Linux-qualified only")

checkpointStorageFixture <- function(store, object, request_id = "checkpoint") {
    saveWorkflowEntity(store, object, request_id,
        identity = list(input = "complete-input-v1", design = "design-v1", parameters = "parameters-v1",
            selectors = "ordered-samples-v1", software = "test-build-v1"),
        metadata = list(class = class(object), shape = as.integer(dim(object)), factors = list()))
}

test_that("qs2 checkpoints use the same publication and recovery boundaries", {
    skip_if_not_installed("qs2", "0.2.2")
    value <- list(values = c(NA_real_, NaN, Inf, -Inf, 0), group = factor(c("B", "A")))
    identity <- list(input = "qs2-recovery", design = "paired", parameters = "unchanged",
        selectors = "all", software = "test-v1")
    metadata <- list(class = "list", shape = integer(), factors = list())
    for (phase in c("after_payload", "after_record", "after_commit")) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "qs2-recovery", "DIA")
        testthat::with_mocked_bindings({
            expect_error(saveWorkflowCheckpoint(store, value, "interrupted", identity,
                metadata, codec = "qs2"), "injected qs2 interruption")
        }, savedWorkflowStorageBoundary = function(current) {
            if (current == phase) stop("injected qs2 interruption")
        }, .package = "MultiScholaR")
        reopened <- newWorkflowEntityStore(root, "qs2-recovery", "DIA")
        recovered <- recoverSavedWorkflowRequest(reopened, "interrupted")
        if (phase == "after_commit") {
            expect_identical(hydrateSavedWorkflowCheckpoint(reopened, recovered), value)
        } else expect_null(recovered)
        accepted <- saveWorkflowCheckpoint(reopened, value, "interrupted", identity,
            metadata, codec = "qs2")
        expect_identical(hydrateSavedWorkflowCheckpoint(reopened, accepted), value)
    }
})

test_that("incomplete stages stay invisible and committed requests recover after lost acknowledgement", {
    for (phase in c("before_write", "after_payload", "after_record", "after_commit")) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "project", "DIA")
        object <- list(values = c(1, NA_real_, Inf), group = factor(c("A", "B")))
        testthat::with_mocked_bindings({
            expect_error(checkpointStorageFixture(store, object), "injected storage failure")
        }, savedWorkflowStorageBoundary = function(current) {
            if (identical(current, phase)) stop("injected storage failure")
        }, .package = "MultiScholaR")
        reopened <- newWorkflowEntityStore(root, "project", "DIA")
        recovered <- recoverSavedWorkflowRequest(reopened, "checkpoint")
        if (phase == "after_commit") {
            expect_identical(hydrateSavedWorkflowEntity(reopened, recovered), object)
            expect_identical(reopened$count, 1L)
        } else {
            expect_null(recovered)
            expect_identical(reopened$count, 0L)
        }
        accepted <- checkpointStorageFixture(reopened, object)
        expect_identical(hydrateSavedWorkflowEntity(reopened, accepted), object)
        expect_identical(reopened$count, 1L)
    }
})

test_that("partial serialization preserves prior accepted data", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    accepted <- checkpointStorageFixture(store, list(value = 1), "previous")
    testthat::with_mocked_bindings({
        expect_error(checkpointStorageFixture(store, list(value = 2), "new"), class = "saved_entity_write")
    }, writeSavedWorkflowPayload = function(object, path) {
        writeBin(as.raw(c(1, 2, 3)), path)
        stop("simulated serialization/disk failure")
    }, .package = "MultiScholaR")
    expect_null(recoverSavedWorkflowRequest(store, "new"))
    expect_identical(hydrateSavedWorkflowEntity(store, accepted), list(value = 1))
})

test_that("permission and disk-full write failures preserve the accepted checkpoint", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    accepted <- checkpointStorageFixture(store, list(value = 1), "previous")
    for (reason in c("Permission denied", "No space left on device")) {
        testthat::with_mocked_bindings({
            expect_error(checkpointStorageFixture(store, list(value = 2), "new"),
                class = "saved_entity_write")
        }, writeSavedWorkflowPayload = function(object, path) stop(reason),
            .package = "MultiScholaR")
        expect_null(recoverSavedWorkflowRequest(store, "new"))
        expect_identical(hydrateSavedWorkflowEntity(store, accepted), list(value = 1))
    }
})

test_that("flush failure cannot acknowledge success and published requests remain recoverable", {
    for (failure_call in c(1L, 2L)) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "project", "DIA")
        prior <- checkpointStorageFixture(store, list(value = 1), "previous")
        calls <- 0L
        real_sync <- syncSavedWorkflowPaths
        testthat::with_mocked_bindings({
            expect_error(checkpointStorageFixture(store, list(value = 2), "flush-failure"), class = "saved_entity_write")
        }, syncSavedWorkflowPaths = function(store, paths) {
            calls <<- calls + 1L
            if (calls == failure_call) cli::cli_abort("injected flush failure", class = "saved_entity_write")
            real_sync(store, paths)
        }, .package = "MultiScholaR")
        expect_identical(hydrateSavedWorkflowEntity(store, prior), list(value = 1))
        reopened <- newWorkflowEntityStore(root, "project", "DIA")
        recovered <- recoverSavedWorkflowRequest(reopened, "flush-failure")
        if (failure_call == 1L) {
            expect_null(recovered)
        } else {
            expect_identical(hydrateSavedWorkflowEntity(reopened, recovered), list(value = 2))
            expect_identical(reopened$count, 1L)
        }
    }
})

test_that("flush targets can never mean the whole filesystem or unrelated files", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    outside <- withr::local_tempfile()
    writeLines("unrelated", outside)
    expect_error(syncSavedWorkflowPaths(store, character()), class = "saved_entity_path")
    expect_error(syncSavedWorkflowPaths(store, outside), class = "saved_entity_path")
    expect_error(syncSavedWorkflowPaths(store, "/"), class = "saved_entity_path")
    expect_true(syncSavedWorkflowPaths(store, file.path(store$directory, "store.json")))
})

test_that("cancellation is distinct and cannot activate a staged checkpoint", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    prior <- checkpointStorageFixture(store, list(value = 1), "previous")
    testthat::with_mocked_bindings({
        expect_error(checkpointStorageFixture(store, list(value = 2), "cancelled"), class = "saved_entity_cancelled")
    }, savedWorkflowStorageBoundary = function(phase) {
        if (phase == "after_payload") signalCondition(structure(list(message = "cancel"), class = c("interrupt", "condition")))
    }, .package = "MultiScholaR")
    expect_null(recoverSavedWorkflowRequest(store, "cancelled"))
    expect_identical(hydrateSavedWorkflowEntity(store, prior), list(value = 1))
})

test_that("retry cleanup cannot remove a store root or unexpected files", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    accepted <- checkpointStorageFixture(store, list(value = 1), "previous")
    expect_error(discardRedundantWorkflowStage(store, store$directory, accepted), class = "saved_entity_path")
    staging <- file.path(store$directory, ".pending", "previous-owned-test")
    dir.create(staging)
    expect_true(file.copy(file.path(store$directory, accepted$payload$relative_path),
        file.path(staging, "payload.rds")))
    writeLines("preserve", file.path(staging, "unexpected.txt"))
    expect_error(discardRedundantWorkflowStage(store, staging, accepted), class = "saved_entity_path")
    expect_true(file.exists(file.path(staging, "unexpected.txt")))
    expect_identical(hydrateSavedWorkflowEntity(store, accepted), list(value = 1))
})

test_that("SIGKILL at storage boundaries is recovered by a fresh R process", {
    skip_on_os("windows")
    library <- dirname(getNamespaceInfo(asNamespace("MultiScholaR"), "path"))
    writer <- function(root, library, phase, ready) {
        .libPaths(c(library, .libPaths()))
        api <- asNamespace("MultiScholaR")
        stopifnot(identical(normalizePath(getNamespaceInfo(api, "path")),
            normalizePath(file.path(library, "MultiScholaR"))))
        pause <- function() {
            writeLines("ready", ready)
            repeat Sys.sleep(0.05)
        }
        testthat::local_mocked_bindings(savedWorkflowStorageBoundary = function(current) {
            if (identical(current, phase)) pause()
        }, .package = "MultiScholaR")
        if (phase == "during_serialization") {
            testthat::local_mocked_bindings(writeSavedWorkflowPayload = function(object, path) {
                bytes <- serialize(object, NULL, version = 3)
                connection <- file(path, "wb")
                writeBin(bytes[seq_len(10L)], connection)
                flush(connection)
                pause()
            }, .package = "MultiScholaR")
        }
        store <- api$newWorkflowEntityStore(root, "project", "DIA")
        object <- list(values = c(1, NA_real_, Inf), group = ordered(c("A", "B")))
        api$saveWorkflowEntity(store, object, "checkpoint",
            list(input = "input", design = "design", parameters = "parameters", selectors = "order", software = "test"),
            list(class = "list", shape = integer(), factors = list()))
    }
    reader <- function(root, library) {
        .libPaths(c(library, .libPaths()))
        api <- asNamespace("MultiScholaR")
        stopifnot(identical(normalizePath(getNamespaceInfo(api, "path")),
            normalizePath(file.path(library, "MultiScholaR"))))
        store <- api$newWorkflowEntityStore(root, "project", "DIA")
        record <- api$recoverSavedWorkflowRequest(store, "checkpoint")
        list(found = !is.null(record), count = store$count,
            value = if (is.null(record)) NULL else api$hydrateSavedWorkflowEntity(store, record))
    }
    environment(writer) <- baseenv()
    environment(reader) <- baseenv()
    for (phase in c("before_write", "during_serialization", "after_payload", "after_record", "after_rename", "after_commit")) {
        root <- withr::local_tempdir()
        ready <- file.path(root, "writer-ready")
        process <- callr::r_bg(writer, list(root, library, phase, ready), stdout = "|", stderr = "|")
        withr::defer(if (process$is_alive()) process$kill())
        deadline <- Sys.time() + 15
        while (!file.exists(ready) && process$is_alive() && Sys.time() < deadline) Sys.sleep(0.02)
        expect_true(file.exists(ready), info = phase)
        expect_true(process$is_alive(), info = phase)
        process$signal(9L)
        process$wait(timeout = 5000)
        expect_false(process$is_alive())
        expect_equal(process$get_exit_status(), -9L)
        recovered <- callr::r(reader, list(root, library))
        published <- phase %in% c("after_rename", "after_commit")
        expect_identical(recovered$found, published)
        expect_identical(recovered$count, if (published) 1L else 0L)
        if (published) {
            expect_identical(recovered$value,
                list(values = c(1, NA_real_, Inf), group = ordered(c("A", "B"))))
        }
    }
})

test_that("concurrent writers cannot overwrite an accepted request", {
    library <- dirname(getNamespaceInfo(asNamespace("MultiScholaR"), "path"))
    writer <- function(root, library, label, value) {
        .libPaths(c(library, .libPaths()))
        api <- asNamespace("MultiScholaR")
        stopifnot(identical(normalizePath(getNamespaceInfo(api, "path")),
            normalizePath(file.path(library, "MultiScholaR"))))
        testthat::local_mocked_bindings(savedWorkflowStorageBoundary = function(phase) {
            if (phase == "after_record") {
                writeLines("ready", file.path(root, paste0(label, "-ready")))
                while (!file.exists(file.path(root, "release"))) Sys.sleep(0.02)
            }
        }, .package = "MultiScholaR")
        store <- api$newWorkflowEntityStore(root, "project", "DIA")
        tryCatch({
            record <- api$saveWorkflowEntity(store, list(value = value), "race",
                list(input = "input", design = "design", parameters = "parameters", selectors = "order", software = "test"),
                list(class = "list", shape = integer(), factors = list()))
            list(status = "accepted", id = record$entity_id, value = value)
        }, error = function(error) list(status = class(error)[[1L]], message = conditionMessage(error)))
    }
    environment(writer) <- baseenv()
    for (same_value in c(TRUE, FALSE)) {
        root <- withr::local_tempdir()
        store <- newWorkflowEntityStore(root, "project", "DIA")
        first <- callr::r_bg(writer, list(root, library, "first", 1L), stdout = "|", stderr = "|")
        second <- callr::r_bg(writer, list(root, library, "second", if (same_value) 1L else 2L),
            stdout = "|", stderr = "|")
        withr::defer({if (first$is_alive()) first$kill(); if (second$is_alive()) second$kill()})
        ready <- file.path(root, c("first-ready", "second-ready"))
        deadline <- Sys.time() + 15
        while (!all(file.exists(ready)) && first$is_alive() && second$is_alive() && Sys.time() < deadline) {
            Sys.sleep(0.02)
        }
        expect_true(all(file.exists(ready)))
        writeLines("release", file.path(root, "release"))
        first$wait(timeout = 10000)
        second$wait(timeout = 10000)
        results <- list(first$get_result(), second$get_result())
        statuses <- vapply(results, `[[`, character(1), "status")
        if (same_value) {
            expect_identical(statuses, c("accepted", "accepted"))
            expect_identical(results[[1L]]$id, results[[2L]]$id)
        } else {
            expect_setequal(statuses, c("accepted", "saved_entity_conflict"))
        }
        accepted <- recoverSavedWorkflowRequest(store, "race")
        winner <- results[[which(statuses == "accepted")[[1L]]]]
        expect_identical(hydrateSavedWorkflowEntity(store, accepted), list(value = winner$value))
        expect_identical(store$count, 1L)
        expect_length(list.dirs(file.path(store$directory, "entities"), recursive = FALSE), 1L)
        expect_length(list.dirs(file.path(store$directory, ".pending"), recursive = FALSE),
            if (same_value) 0L else 1L)
    }
})
