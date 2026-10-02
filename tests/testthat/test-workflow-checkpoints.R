testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"), "durable saved-entity backend is Linux-qualified only")

checkpointStorageFixture <- function(store, object, request_id = "checkpoint") {
    saveWorkflowCheckpoint(store, object, request_id,
        identity = list(input = "complete-input-v1", design = "design-v1", parameters = "parameters-v1",
            selectors = "ordered-samples-v1", software = "test-build-v1"),
        metadata = list(class = class(object), shape = as.integer(dim(object)), factors = list()))
}

test_that("explicit qs2 and default RDS checkpoints reopen with exact S4 contents", {
    skip_if_not_installed("qs2", "0.2.2")
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = c("P3", "P1", "P2"),
            S1 = c(NA_real_, NaN, Inf), S2 = c(-Inf, 0, -0)),
        design_matrix = data.frame(Run = c("S1", "S2"),
            group = ordered(c("B", "A"), levels = c("B", "A", "unused"))),
        sample_id = "Run", args = list(nested = list(values = c(x = 1L, y = 2L))))
    attr(object@protein_quant_table$S1, "units") <- "log2"
    value <- list(data = object, config = list(flag = TRUE),
        timestamp = as.POSIXct("2026-09-16", tz = "UTC"))
    identity <- list(input = "qs2-fixture", design = "paired", parameters = "unchanged",
        selectors = "exact-order", software = "test-v1")
    metadata <- list(class = "list", shape = integer(), factors = list())
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "mixed-codecs", "DIA")
    native <- saveWorkflowCheckpoint(store, value, "native", identity, metadata)
    compressed <- saveWorkflowCheckpoint(store, value, "compressed", identity, metadata,
        codec = "qs2")
    expect_identical(native$payload$codec, "rds")
    expect_identical(compressed$payload$codec, "qs2")
    expect_match(compressed$payload$relative_path, "payload.qs2", fixed = TRUE)
    expect_identical(saveWorkflowCheckpoint(store, value, "compressed", identity, metadata,
        codec = "qs2"), compressed)
    closeWorkflowEntityStore(store)
    reopened <- newWorkflowEntityStore(root, "mixed-codecs", "DIA")
    for (reference in list(native, compressed)) {
        expect_identical(hydrateSavedWorkflowCheckpoint(reopened, reference), value)
        expect_true(savedWorkflowEntityCapabilities(reference)$full_hydration)
        expect_false(savedWorkflowEntityCapabilities(reference)$column_projection)
    }
    expect_identical(recoverSavedWorkflowRequest(reopened, "compressed"), compressed)
    path <- file.path(reopened$directory, compressed$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.size(path))
    bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
    writeBin(bytes, path)
    expect_error(hydrateSavedWorkflowCheckpoint(reopened, compressed), class = "saved_entity_corrupt")
    expect_identical(hydrateSavedWorkflowCheckpoint(reopened, native), value)
})

test_that("DIA worker configuration splits and reconnects without changing scientific settings", {
    worker <- new.env(parent = emptyenv())
    worker$is_alive <- function() TRUE
    cluster <- structure(list(worker, worker), class = "multidplyr_cluster")
    config <- list(globalParameters = list(number_of_cpus = 2L),
        peptideIntensityFiltering = list(core_utilisation = cluster, cutoff = 0.5),
        rollUpPrecursorToPeptide = list(core_utilisation = cluster, method = "iq"),
        untouched = list(core_utilisation = NA, groups = ordered(c("B", "A"))))
    split <- splitDIAWorkerConfig(config)
    expect_identical(split$runtime, cluster)
    expect_identical(config$peptideIntensityFiltering$core_utilisation, cluster)
    expect_identical(split$config$untouched, config$untouched)
    expect_identical(split$config$globalParameters, config$globalParameters)
    expect_identical(split$config$peptideIntensityFiltering$cutoff, 0.5)
    expect_identical(bindDIAWorkerConfig(split$config, cluster), config)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, split$config)
    expect_identical(hydrateSavedWorkflowEntity(store, reference), split$config)
    expect_identical(bindDIAWorkerConfig(hydrateSavedWorkflowEntity(store, reference), cluster), config)
    expect_error(bindDIAWorkerConfig(split$config, NULL), class = "saved_entity_resource")
    wrong_count <- config
    wrong_count$globalParameters$number_of_cpus <- 3L
    expect_error(splitDIAWorkerConfig(wrong_count), class = "saved_entity_conflict")
    worker$is_alive <- function() FALSE
    expect_error(bindDIAWorkerConfig(split$config, cluster), class = "saved_entity_resource")
})

test_that("an explicitly owned DIA worker can be recreated and released", {
    caller_search <- search()
    spec <- newDIAWorkerSpec(1L)
    cluster <- startDIAWorkerCluster(spec)
    on.exit(stopDIAWorkerCluster(cluster), add = TRUE)
    expect_identical(search(), caller_search)
    config <- list(globalParameters = list(number_of_cpus = 1L),
        rollUpPrecursorToPeptide = list(core_utilisation = spec, method = "iq"))
    bound <- bindDIAWorkerConfig(config, cluster)
    expect_identical(bound$rollUpPrecursorToPeptide$core_utilisation, cluster)
    result <- cluster[[1L]]$run(function() {
        list(value = sum(c(1, 2, NA_real_), na.rm = TRUE),
            packages = all(paste0("package:", c("tidyverse", "glue", "rlang", "lazyeval")) %in% search()))
    })
    expect_identical(result, list(value = 3, packages = TRUE))
    expect_identical(splitDIAWorkerConfig(bound)$config, config)
    stopDIAWorkerCluster(cluster)
    expect_false(cluster[[1L]]$is_alive())
    expect_error(bindDIAWorkerConfig(config, cluster), class = "saved_entity_resource")
})

test_that("failed DIA worker initialization releases the newly owned cluster", {
    killed <- 0L
    worker <- new.env(parent = emptyenv())
    worker$kill <- function() killed <<- killed + 1L
    cluster <- structure(list(worker), class = "multidplyr_cluster")
    testthat::with_mocked_bindings({
        expect_error(startDIAWorkerCluster(newDIAWorkerSpec(1L)), "injected package initialization failure")
    }, new_cluster = function(n) cluster,
        cluster_send = function(...) stop("injected package initialization failure"),
        .package = "multidplyr")
    expect_identical(killed, 1L)
})

test_that("invalid DIA backend settings fail before startup or rebinding", {
    for (workers in list(0, -1, NA_real_, Inf, 1.5, "2", c(1, 2))) {
        expect_error(newDIAWorkerSpec(workers), class = "saved_entity_invalid")
    }
    expect_error(splitDIAWorkerConfig(NULL), class = "saved_entity_invalid")
    expect_error(bindDIAWorkerConfig(NULL, NULL), class = "saved_entity_invalid")
    expect_error(startDIAWorkerCluster(list(workers = 1L)), class = "saved_entity_invalid")
    spec <- newDIAWorkerSpec(1L)
    spec$packages <- "unapproved-package"
    expect_error(startDIAWorkerCluster(spec), class = "saved_entity_invalid")
    expect_error(bindDIAWorkerConfig(list(step = list(core_utilisation = spec)), NULL),
        class = "saved_entity_invalid")
    first <- structure(list(new.env(parent = emptyenv())), class = "multidplyr_cluster")
    second <- structure(list(new.env(parent = emptyenv())), class = "multidplyr_cluster")
    expect_error(splitDIAWorkerConfig(list(first = list(core_utilisation = first),
        second = list(core_utilisation = second))), class = "saved_entity_unsupported")
})

test_that("DIA history preparation preserves science and only rebinds selected state settings", {
    worker <- new.env(parent = emptyenv())
    worker$is_alive <- function() TRUE
    cluster <- structure(list(worker), class = "multidplyr_cluster")
    config <- list(globalParameters = list(number_of_cpus = 1L),
        peptideIntensityFiltering = list(core_utilisation = cluster, cutoff = 0.5))
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = "P1", S1 = NA_real_),
        design_matrix = data.frame(Run = "S1", group = factor("A")), sample_id = "Run", args = config)
    history <- list(raw = list(data = object, config = config, description = "original"))
    prepared <- prepareDIAWorkflowHistory(history)
    expect_identical(prepared$raw$data@protein_quant_table, object@protein_quant_table)
    expect_identical(prepared$raw$data@design_matrix, object@design_matrix)
    expect_identical(prepared$raw$description, "original")
    expect_identical(history$raw$data@args, config)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, prepared)
    restored <- hydrateSavedWorkflowEntity(store, reference)
    expect_identical(restored, prepared)
    expect_identical(bindDIAWorkflowState(restored$raw, cluster), history$raw)
})

test_that("native RDS readr empty parser markers are preserved without admitting live pointers", {
    table <- readr::read_tsv(I("id\tvalue\nP1\t1\n"), show_col_types = FALSE, progress = FALSE)
    expect_error(checkSavedWorkflowPayload(table), class = "saved_entity_unsupported")
    restored <- unserialize(serialize(table, NULL))
    expect_identical(attr(restored, "problems"), methods::new("externalptr"))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, restored)
    expect_identical(hydrateSavedWorkflowEntity(store, reference), restored)
    expect_error(checkSavedWorkflowPayload(list(pipe = methods::new("externalptr"))),
        class = "saved_entity_unsupported")
})

test_that("DIA parser preparation preserves diagnostics, values and all other attributes", {
    for (bad_input in c(FALSE, TRUE)) {
        input <- if (bad_input) "id\tvalue\nP1\t1\nP2\tinvalid\n" else "id\tvalue\nP1\t1\nP2\tNA\n"
        # The malformed numeric field intentionally emits a parser warning.
        table <- suppressWarnings(readr::read_tsv(I(input),
            col_types = readr::cols(id = readr::col_character(), value = readr::col_double()),
            progress = FALSE))
        diagnostics <- readr::problems(table)
        marker <- attr(table, "problems")
        other_attributes <- attributes(table)
        other_attributes$problems <- NULL
        prepared <- prepareDIAParserTable(table)
        expect_identical(attr(table, "problems"), marker)
        expect_identical(readr::problems(prepared), diagnostics)
        expect_identical(nrow(diagnostics), if (bad_input) 1L else 0L)
        expect_identical(prepared$id, table$id)
        expect_identical(prepared$value, table$value)
        prepared_attributes <- attributes(prepared)
        prepared_attributes$problems <- NULL
        expect_identical(prepared_attributes, other_attributes)
        expect_identical(prepareDIAParserTable(prepared), prepared)
        store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
        reference <- checkpointStorageFixture(store, prepared)
        restored <- hydrateSavedWorkflowCheckpoint(store, reference)
        expect_identical(restored, prepared)
        expect_identical(readr::problems(restored), diagnostics)
    }
})

test_that("DIA history preparation handles parser pointers in S4 and raw table records", {
    table <- readr::read_tsv(I("Protein.Ids\tS1\nP1\t1\n"), show_col_types = FALSE, progress = FALSE)
    design <- readr::read_tsv(I("Run\tgroup\nS1\tA\n"), show_col_types = FALSE, progress = FALSE)
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = table,
        design_matrix = design, sample_id = "Run")
    history <- list(raw = list(data = table, config = list()),
        protein = list(data = object, config = list(), description = "parser-bearing S4"))
    prepared <- prepareDIAWorkflowHistory(history)
    expect_type(attr(history$raw$data, "problems"), "externalptr")
    expect_identical(prepared$raw$data, prepareDIAParserTable(table))
    expect_identical(prepared$protein$data@design_matrix, prepareDIAParserTable(design))
    expect_true(methods::validObject(prepared$protein$data))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, prepared)
    expect_identical(hydrateSavedWorkflowCheckpoint(store, reference), prepared)
    expect_identical(prepareDIAWorkflowHistory(prepared), prepared)
})

test_that("complete S4 checkpoints hydrate exactly after reopening the store", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "project", "DIA")
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2"), S1 = c(NA, Inf), S2 = c(NaN, -Inf)),
        design_matrix = data.frame(Run = c("S1", "S2"), group = ordered(c("A", "B"))), sample_id = "Run")
    before <- serialize(object, NULL)
    reference <- checkpointStorageFixture(store, object)
    expect_identical(serialize(object, NULL), before)
    expect_identical(hydrateSavedWorkflowEntity(store, reference), object)
    expect_identical(hydrateSavedWorkflowCheckpoint(store, reference), object)
    wrong_kind <- reference
    wrong_kind$kind <- "da_fit"
    expect_error(hydrateSavedWorkflowCheckpoint(store, wrong_kind), class = "saved_entity_invalid")
    reopened <- newWorkflowEntityStore(root, "project", "DIA")
    recovered <- recoverSavedWorkflowRequest(reopened, "checkpoint")
    expect_identical(recovered, reference)
    expect_identical(hydrateSavedWorkflowEntity(reopened, recovered), object)
    expect_identical(checkpointStorageFixture(reopened, object), reference)
    expect_identical(reopened$count, 1L)
    expect_length(list.dirs(file.path(reopened$directory, ".pending"), recursive = FALSE), 0L)
    changed <- object
    changed@args$changed <- TRUE
    expect_error(checkpointStorageFixture(reopened, changed), class = "saved_entity_conflict")
    expect_length(list.dirs(file.path(reopened$directory, ".pending"), recursive = FALSE), 1L)
    expect_identical(hydrateSavedWorkflowEntity(reopened, reference), object)
})

test_that("corrupt payloads and cross-project references fail before hydration", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, list(values = c(1, NA, NaN, Inf)))
    other <- newWorkflowEntityStore(withr::local_tempdir(), "other", "DIA")
    expect_error(hydrateSavedWorkflowEntity(other, reference), class = "saved_entity_scope")
    path <- file.path(store$directory, reference$payload$relative_path)
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
    writeBin(bytes, path)
    expect_error(hydrateSavedWorkflowEntity(store, reference), class = "saved_entity_corrupt")
    expect_null(recoverSavedWorkflowRequest(store, "absent"))
})

test_that("unsupported resources cannot become accepted checkpoints", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    connection <- file(withr::local_tempfile(), "wb")
    on.exit(close(connection), add = TRUE)
    expect_error(checkpointStorageFixture(store, list(connection = connection)), class = "saved_entity_unsupported")
    expect_null(recoverSavedWorkflowRequest(store, "checkpoint"))
    expect_identical(store$count, 0L)
})

test_that("worker handles embedded in scientific arguments are rejected without mutation", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    pipe <- methods::new("externalptr")
    worker <- new.env(parent = emptyenv())
    worker$pipe <- pipe
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = "P1", S1 = 1),
        design_matrix = data.frame(Run = "S1", group = "A"), sample_id = "Run",
        args = list(peptideIntensityFiltering = list(
            core_utilisation = structure(list(worker), class = "multidplyr_cluster"))))
    expect_error(checkpointStorageFixture(store, object), class = "saved_entity_unsupported")
    expect_identical(worker$pipe, pipe)
    expect_identical(object@args$peptideIntensityFiltering$core_utilisation[[1L]], worker)
    expect_identical(store$count, 0L)
    expect_null(recoverSavedWorkflowRequest(store, "checkpoint"))
})

test_that("metadata corruption and missing requests do not damage other accepted data", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    good <- checkpointStorageFixture(store, list(value = 1), "good")
    broken <- checkpointStorageFixture(store, list(value = 2), "broken")
    directory <- dirname(file.path(store$directory, broken$payload$relative_path))
    metadata_path <- file.path(directory, "descriptor.rds")
    bytes <- readBin(metadata_path, "raw", n = file.info(metadata_path)$size)
    bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
    writeBin(bytes, metadata_path)
    expect_error(hydrateSavedWorkflowEntity(store, broken), class = "saved_entity_corrupt")
    expect_identical(hydrateSavedWorkflowEntity(store, good), list(value = 1))
    absent <- checkpointStorageFixture(store, list(value = 3), "absent")
    directory <- dirname(file.path(store$directory, absent$payload$relative_path))
    expect_true(file.rename(directory, file.path(store$directory, "preserved-request")))
    expect_error(hydrateSavedWorkflowEntity(store, absent), class = "saved_entity_missing")
    malformed <- checkpointStorageFixture(store, list(value = 4), "malformed")
    writeLines("{broken", file.path(dirname(file.path(store$directory, malformed$payload$relative_path)), "commit.json"))
    expect_error(hydrateSavedWorkflowEntity(store, malformed), class = "saved_entity_corrupt")
})

test_that("a complete resident R6 history can be saved and reverted after hydration", {
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = "P1", S1 = 1, S2 = NA_real_),
        design_matrix = data.frame(Run = c("S1", "S2"), group = factor(c("A", "B"))), sample_id = "Run")
    manager <- WorkflowState$new()
    manager$saveState("raw", object, list(), "raw")
    changed <- object
    changed@args$revision <- 1L
    manager$saveState("changed", changed, list(), "changed")
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    reference <- checkpointStorageFixture(store, manager, "history")
    restored <- hydrateSavedWorkflowEntity(store, reference)
    expect_s3_class(restored, "WorkflowState")
    expect_identical(restored$getState(), changed)
    expect_identical(restored$revertToState("raw"), object)
    expect_identical(manager$getState(), changed)
})

test_that("hydrated history returns to the public resident manager without a store dependency", {
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = "P1", S1 = 1, S2 = 2),
        design_matrix = data.frame(Run = c("S1", "S2"), group = factor(c("A", "B"))), sample_id = "Run")
    manager <- WorkflowState$new()
    manager$saveState("raw", object, list(), "raw")
    changed <- object
    changed@args$revision <- 1L
    manager$saveState("changed", changed, list(), "changed")
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    prepared <- prepareDIAWorkflowHistory(manager$states)
    reference <- checkpointStorageFixture(store, prepared)
    restored <- WorkflowState$new()
    restored$states <- hydrateSavedWorkflowCheckpoint(store, reference)
    restored$state_history <- manager$state_history
    restored$current_state <- manager$current_state
    expect_identical(names(restored$states), names(manager$states))
    expect_identical(lapply(restored$states, \(state) state$data),
        lapply(manager$states, \(state) state$data))
    expect_identical(getDataMatrix(restored$getState()), getDataMatrix(changed))
    restored$states$raw$description <- "external replacement"
    expect_identical(restored$states$raw$description, "external replacement")
    expect_identical(restored$revertToState("raw"), object)
    restored$saveState("branch", object, list(), "branch")
    expect_true("changed" %in% names(restored$states))
    expect_false("changed" %in% restored$state_history)
    legacy <- withr::local_tempfile()
    saveRDS(restored, legacy)
    independent <- readRDS(legacy)
    expect_identical(independent$getState(), object)
    expect_identical(independent$states, restored$states)
})
