savedIdentityFixture <- function(store, value = "one", factors = list()) {
    fingerprint <- digest::digest(value, algo = "sha256")
    newSavedWorkflowDescriptor(store, "checkpoint",
        identity = list(input = fingerprint, design = "design-v1", parameters = "parameters-v1",
            selectors = "ordered-selection-v1", software = "MultiScholaR-test"),
        payload = list(relative_path = paste0(fingerprint, ".rds"), bytes = 123,
            sha256 = fingerprint, codec = "rds"),
        metadata = list(class = "ProteinQuantitativeData", shape = c(100L, 12L), factors = factors))
}

test_that("metadata lookup distinguishes omitted, NA and exact typed factors", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    expect_s3_class(store, "R6")
    expect_null(store$clone)
    missing <- savedIdentityFixture(store, "missing", list(batch = NA_character_))
    typed <- savedIdentityFixture(store, "factor", list(batch = factor("A", levels = c("A", "unused"))))
    for (record in list(missing, typed)) registerSavedWorkflowMetadata(store, record)
    expect_identical(findSavedWorkflowEntities(store)$status, "multiple")
    expect_identical(findSavedWorkflowEntities(store, list(factors = missing$metadata$factors))$status, "hit")
    expect_identical(findSavedWorkflowEntities(store, list(factors = list(batch = "A")))$status, "miss")
    expect_identical(describeSavedWorkflowEntity(store, typed$entity_id), typed)
    expect_error(findSavedWorkflowEntities(store, list(latest = TRUE)), class = "saved_entity_invalid")
    expect_error(findSavedWorkflowEntities(store, list(revision = "0")), class = "saved_entity_invalid")
    expect_error(findSavedWorkflowEntities(store, list(factors = list(batch = new.env()))), class = "saved_entity_invalid")
    expect_error(findSavedWorkflowEntities(store, list(identity = list(unknown = "value"))), class = "saved_entity_invalid")
    expect_error(findSavedWorkflowEntities(store, max_matches = 1), class = "saved_entity_limit")
    expect_identical(registerSavedWorkflowMetadata(store, typed), typed$entity_id)
    expect_identical(store$count, 2L)
})

test_that("identities bind complete effective inputs without opening payloads", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    first <- savedIdentityFixture(store, list(E = matrix(1, 2, 2), se = rep(1, 2), order = 1:2))
    second <- savedIdentityFixture(store, list(E = matrix(1, 2, 2), se = rep(2, 2), order = 1:2))
    third <- savedIdentityFixture(store, list(E = matrix(1, 2, 2), se = rep(1, 2), order = 2:1))
    expect_length(unique(c(first$entity_id, second$entity_id, third$entity_id)), 3L)
    registerSavedWorkflowMetadata(store, first)
    # No payload file exists. Queries/descriptions cannot rely on decoding one.
    expect_false(file.exists(file.path(store$directory, first$payload$relative_path)))
    expect_identical(findSavedWorkflowEntities(store, list(identity = second$identity))$status, "miss")
    expect_identical(findSavedWorkflowEntities(store, list(identity = first$identity))$matches[[1]], first)
    changed <- first
    changed$payload$bytes <- 124
    expect_error(registerSavedWorkflowMetadata(store, changed), class = "saved_entity_conflict")
    for (field in names(first$identity)) {
        identity <- first$identity
        identity[[field]] <- paste0(identity[[field]], "-changed")
        descriptor <- newSavedWorkflowDescriptor(store, first$kind, identity, first$payload, first$metadata)
        expect_false(identical(descriptor$entity_id, first$entity_id))
        expect_identical(findSavedWorkflowEntities(store, list(identity = identity))$status, "miss")
    }
})

test_that("EList uncertainty and aligned ordering participate in input identities", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    values <- matrix(1, 2, 2, dimnames = list(c("P1", "P2"), c("S1", "S2")))
    # Identity fixture, not a claim that dpcQuant was executed.
    original <- methods::new("EList", list(E = values,
        other = list(standard.error = values * 0.2, n.observations = values * 3)))
    se <- original
    se$other$standard.error[1, 1] <- 0.4
    counts <- original
    counts$other$n.observations[1, 1] <- 4
    reordered <- original
    reordered$E <- reordered$E[2:1, , drop = FALSE]
    reordered$other <- lapply(reordered$other, \(x) x[2:1, , drop = FALSE])
    records <- lapply(list(original, se, counts, reordered), \(x) savedIdentityFixture(store, x))
    expect_identical(original$E, se$E)
    expect_identical(original$E, counts$E)
    expect_length(unique(vapply(records, `[[`, character(1), "entity_id")), 4L)
    registerSavedWorkflowMetadata(store, records[[1L]])
    for (record in records[-1L]) {
        expect_identical(findSavedWorkflowEntities(store, list(identity = record$identity))$status, "miss")
    }
})

test_that("replacing the project directory at the same path invalidates old handles", {
    directory <- withr::local_tempdir()
    project <- file.path(directory, "project")
    dir.create(project)
    old <- newWorkflowEntityStore(project, "project", "DIA")
    record <- savedIdentityFixture(old)
    registerSavedWorkflowMetadata(old, record)
    expect_true(file.rename(project, file.path(directory, "preserved-project")))
    dir.create(project)
    replacement <- newWorkflowEntityStore(project, "project", "DIA")
    expect_false(identical(old$identity$incarnation, replacement$identity$incarnation))
    expect_error(describeSavedWorkflowEntity(old, record$entity_id), class = "saved_entity_scope")
    expect_error(registerSavedWorkflowMetadata(replacement, record), class = "saved_entity_scope")
})

test_that("ownership, path traversal and index ceilings are enforced", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "project", "DIA", max_entries = 1L)
    other <- newWorkflowEntityStore(withr::local_tempdir(), "other", "DIA")
    first <- savedIdentityFixture(store)
    expect_error(registerSavedWorkflowMetadata(other, first), class = "saved_entity_scope")
    for (path in c("../outside", "/absolute", "a/../b", "a\\b", "a//b")) {
        expect_error(checkSavedEntityPath(store, path), class = "saved_entity_path")
    }
    registerSavedWorkflowMetadata(store, first)
    expect_error(registerSavedWorkflowMetadata(store, savedIdentityFixture(store, "two")), class = "saved_entity_limit")
    expect_error(newWorkflowEntityStore(root, "wrong-project", "DIA"), class = "saved_entity_scope")
    marker <- file.path(store$directory, "store.json")
    replacement <- store$identity
    replacement$incarnation <- digest::digest("new-project-incarnation", algo = "sha256")
    jsonlite::write_json(replacement, marker, auto_unbox = TRUE)
    expect_error(describeSavedWorkflowEntity(store, first$entity_id), class = "saved_entity_scope")
})

test_that("malformed markers and hidden reference metadata fail with classed errors", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "project", "DIA")
    record <- savedIdentityFixture(store)
    metadata <- record$metadata
    attr(metadata$class, "hidden") <- new.env()
    expect_error(newSavedWorkflowDescriptor(store, record$kind, record$identity, record$payload, metadata),
        class = "saved_entity_invalid")
    expect_error(registerSavedWorkflowMetadata(store, 1), class = "saved_entity_invalid")
    writeLines("{broken", file.path(store$directory, "store.json"))
    expect_error(describeSavedWorkflowEntity(store, record$entity_id), class = "saved_entity_invalid")
    expect_error(newWorkflowEntityStore(root, "project", "DIA"), class = "saved_entity_invalid")
})

test_that("symlinks cannot redirect payloads outside the store", {
    skip_on_os("windows")
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    outside <- withr::local_tempdir()
    expect_true(file.symlink(outside, file.path(store$directory, "escape")))
    expect_error(checkSavedEntityPath(store, "escape/payload.rds"), class = "saved_entity_path")
    expect_true(file.symlink(file.path(outside, "absent"), file.path(store$directory, "dangling")))
    expect_error(checkSavedEntityPath(store, "dangling/payload.rds"), class = "saved_entity_path")
})

test_that("a thousand descriptors remain bounded and lookup never calls the payload reader", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA", max_entries = 1000L)
    for (index in seq_len(1000L)) {
        record <- savedIdentityFixture(store, as.character(index), list(batch = as.integer(index)))
        registerSavedWorkflowMetadata(store, record)
    }
    expect_equal(store$count, 1000L)
    expect_lte(store$metadata_bytes, store$max_metadata_bytes)
    reads <- 0L
    testthat::local_mocked_bindings(readRDS = function(...) {
        reads <<- reads + 1L
        stop("payload decoding is forbidden in metadata lookup")
    }, .package = "base")
    for (index in c(1L, 500L, 1000L)) {
        result <- findSavedWorkflowEntities(store, list(factors = list(batch = index)))
        expect_identical(result$status, "hit")
        expect_identical(describeSavedWorkflowEntity(store, result$matches[[1L]]$entity_id), result$matches[[1L]])
    }
    expect_identical(reads, 0L)
    expect_error(findSavedWorkflowEntities(store), class = "saved_entity_limit")
})

test_that("metadata byte ceilings reject admission without changing accepted records", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    record <- savedIdentityFixture(store)
    budget <- as.numeric(object.size(record)) + 512
    bounded <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA", max_metadata_bytes = budget)
    accepted <- savedIdentityFixture(bounded)
    registerSavedWorkflowMetadata(bounded, accepted)
    expect_error(registerSavedWorkflowMetadata(bounded, savedIdentityFixture(bounded, "second")),
        class = "saved_entity_limit")
    expect_identical(describeSavedWorkflowEntity(bounded, accepted$entity_id), accepted)
    expect_identical(bounded$count, 1L)
    expect_error(newWorkflowEntityStore(withr::local_tempdir(), "p", "DIA", max_metadata_bytes = 67108865),
        class = "saved_entity_invalid")
})

test_that("standalone metadata does not change the public resident S4 history path", {
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2"), S1 = c(1, NA), S2 = c(2, Inf)),
        design_matrix = data.frame(Run = c("S1", "S2"), group = factor(c("A", "B"))), sample_id = "Run")
    manager <- WorkflowState$new()
    manager$saveState("raw", object, list(), "public resident state")
    changed <- object
    changed@args$revision <- 1L
    manager$saveState("changed", changed, list(), "public changed state")
    store <- newWorkflowEntityStore(withr::local_tempdir(), "project", "DIA")
    registerSavedWorkflowMetadata(store, savedIdentityFixture(store))
    expect_identical(manager$getState(), changed)
    path <- withr::local_tempfile()
    saveRDS(manager, path)
    restored <- readRDS(path)
    expect_identical(restored$revertToState("raw"), object)
    expect_true(methods::validObject(restored$getState()))
    expect_identical(levels(restored$getState()@design_matrix$group), c("A", "B"))
})

test_that("a commit-time input descriptor matches the archive descriptor and is hashed once", {
    rm(list = ls(.workflowFileHashes, all.names = TRUE), envir = .workflowFileHashes)
    path <- tempfile(fileext = ".tsv")
    on.exit(unlink(path), add = TRUE)
    writeLines(c("Run\tPrecursor.Id", "a\tp1"), path)
    expect_identical(describeDIAInputFile(path), describeDIANativeFile(path))
    real <- digest::digest
    reads <- 0L
    testthat::local_mocked_bindings(digest = function(...) {
        if (!is.null(list(...)$file)) reads <<- reads + 1L
        real(...)
    }, .package = "digest")
    describeDIAInputFile(path)
    expect_identical(reads, 0L)
    writeLines(c("Run\tPrecursor.Id", "b\tp2"), path)
    expect_identical(describeDIAInputFile(path), describeDIANativeFile(path))
    expect_error(describeDIAInputFile(tempdir()))
})
