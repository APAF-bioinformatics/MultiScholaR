testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"),
    "durable saved export is Linux-qualified only")

makeSavedExportProtein <- function(revision = 0L) {
    object <- methods::new("ProteinQuantitativeData",
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2"),
            S1 = c(1, NA_real_), S2 = c(2, Inf)),
        design_matrix = data.frame(Run = c("S1", "S2"),
            group = factor(c("control", "treated"))),
        sample_id = "Run", group_id = "group")
    methods::slot(object, "args")$revision <- revision
    object
}

saveExportCheckpoint <- function(store, object, request_id, selector) {
    saveWorkflowCheckpoint(store, object, request_id,
        identity = list(input = digest::digest(object, algo = "sha256"),
            design = "two-sample-design", parameters = "none",
            selectors = selector, software = "export-test-v1"),
        metadata = list(class = class(object),
            shape = as.integer(dim(methods::slot(object, "protein_quant_table"))),
            factors = list(route = "iq")))
}

makeExportSessionTemplate <- function(first, second) {
    record <- function(data, previous) list(timestamp = as.POSIXct(1, origin = "1970-01-01"),
        data = data, config = list(route = "iq"), description = "fixture",
        previous_state = previous, s4_class = "ProteinQuantitativeData",
        audit_metadata = list(provenance_status = "fixture"))
    list(r6_current_state_name = "changed",
        r6_complete_states = list(raw = record(first, NA_character_),
            raw_alias = record(first, "raw"), changed = record(second, "raw_alias")),
        r6_state_history = as.list(c("raw", "raw_alias", "changed")),
        current_s4_object = second, correlation_filtered_s4 = second,
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        audit_records = list(), final_protein_count = 2L, final_sample_count = 2L)
}

test_that("export hydration decodes identical references once and releases ownership", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "export-project", "DIA")
    first <- saveExportCheckpoint(store, makeSavedExportProtein(0L), "raw", "raw")
    second <- saveExportCheckpoint(store, makeSavedExportProtein(1L), "changed", "changed")
    template <- makeExportSessionTemplate(first, second)
    scope <- newWorkflowExportHydrationScope(store)
    session <- hydrateWorkflowExportSession(scope, template)
    info <- workflowExportHydrationInfo(scope)
    expect_identical(info$unique_decodes, 2L)
    expect_identical(info$reference_uses, 5L)
    expect_identical(info$resident_entries, 2L)
    expect_s4_class(session$r6_complete_states$raw$data,
        "ProteinQuantitativeData")
    expect_identical(session$r6_complete_states$raw$data,
        session$r6_complete_states$raw_alias$data)
    changed_raw <- session$r6_complete_states$raw$data
    methods::slot(changed_raw, "args")$mutated <- TRUE
    expect_null(methods::slot(session$r6_complete_states$raw_alias$data,
        "args")$mutated)
    expect_identical(methods::slot(session$current_s4_object, "args")$revision, 1L)
    releaseWorkflowExportHydrationScope(scope)
    releaseWorkflowExportHydrationScope(scope)
    expect_identical(workflowExportHydrationInfo(scope)$resident_entries, 0L)
    expect_identical(workflowEntityCacheInfo(store)$open_resources, 0L)
})

test_that("one serialization publishes two complete sessions and required copies", {
    output <- withr::local_tempdir()
    session <- makeExportSessionTemplate(makeSavedExportProtein(0L),
        makeSavedExportProtein(1L))
    session$r6_complete_states <- lapply(session$r6_complete_states, \(state) {
        if (isSavedWorkflowReference(state$data)) stop("unexpected reference")
        state
    })
    latest <- file.path(output, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    previous_hash <- digest::digest(file = latest, algo = "sha256")
    tsv <- file.path(output, "source.tsv")
    figure <- file.path(output, "source.png")
    writeLines(c("id\tvalue", "P1\t1"), tsv)
    writeBin(as.raw(c(137L, 80L, 78L, 71L, 1L, 2L, 3L)), figure)
    serializations <- 0L
    serialize_once <- function(object, path) {
        serializations <<- serializations + 1L
        saveRDS(object, path)
    }
    result <- runWorkflowLegacyExport(session, output, "request-one",
        "20260916_160000",
        additional_files = list(table = list(source = tsv,
            destination = "Publication_tables/normalised_results.tsv"),
            figure = list(source = figure,
                destination = "Publication_figures/QC/figure.png")),
        serializeFn = serialize_once)
    expect_identical(serializations, 1L)
    expect_identical(result$status, "complete")
    expect_identical(result$serialization_count, 1L)
    expect_identical(result$scientific_checkpoint_count, 0L)
    expect_false(identical(digest::digest(file = latest, algo = "sha256"),
        previous_hash))
    expect_identical(digest::digest(file = result$timestamped_path, algo = "sha256"),
        digest::digest(file = result$latest_path, algo = "sha256"))
    timestamped <- validateExternalWorkflowLegacySession(result$timestamped_path,
        names(session$r6_complete_states))
    latest_session <- validateExternalWorkflowLegacySession(result$latest_path,
        names(session$r6_complete_states))
    expect_identical(timestamped, session)
    expect_identical(latest_session, session)
    expect_identical(readLines(file.path(output,
        "Publication_tables/normalised_results.tsv")), readLines(tsv))
    expect_identical(digest::digest(file = file.path(output,
        "Publication_figures/QC/figure.png"), algo = "sha256"),
        digest::digest(file = figure, algo = "sha256"))
    expect_true(file.exists(file.path(output, "export_complete_request-one.json")))
})

test_that("publication retry does output work only after timestamp failure", {
    output <- withr::local_tempdir()
    session <- makeExportSessionTemplate(makeSavedExportProtein(0L),
        makeSavedExportProtein(1L))
    serializations <- 0L
    serialize_once <- function(object, path) {
        serializations <<- serializations + 1L
        saveRDS(object, path)
    }
    fail_once <- local({
        failed <- FALSE
        function(phase) {
            if (!failed && identical(phase, "after_timestamped_publication")) {
                failed <<- TRUE
                stop("injected between publications")
            }
        }
    })
    saveRDS(list(previous = TRUE),
        file.path(output, "filtered_session_data_latest.rds"))
    old_hash <- digest::digest(file = file.path(output,
        "filtered_session_data_latest.rds"), algo = "sha256")
    expect_error(runWorkflowLegacyExport(session, output, "retry-output",
        "20260916_160100", serializeFn = serialize_once,
        boundaryFn = fail_once), "injected between publications")
    expect_identical(serializations, 1L)
    expect_identical(digest::digest(file = file.path(output,
        "filtered_session_data_latest.rds"), algo = "sha256"), old_hash)
    expect_true(file.exists(file.path(output,
        "filtered_session_data_20260916_160100.rds")))
    completed <- runWorkflowLegacyExport(NULL, output, "retry-output",
        "20260916_160100", serializeFn = function(...) stop("must not serialize"))
    expect_identical(completed$status, "complete")
    expect_identical(completed$serialization_count, 1L)
    expect_identical(completed$scientific_checkpoint_count, 0L)
})

test_that("pre-ack retry validates published latest without serialization", {
    output <- withr::local_tempdir()
    session <- makeExportSessionTemplate(makeSavedExportProtein(0L),
        makeSavedExportProtein(1L))
    fail_once <- local({
        failed <- FALSE
        function(phase) {
            if (!failed && identical(phase, "before_completion_acknowledgement")) {
                failed <<- TRUE
                stop("injected before acknowledgement")
            }
        }
    })
    expect_error(runWorkflowLegacyExport(session, output, "retry-ack",
        "20260916_160200", boundaryFn = fail_once), "before acknowledgement")
    latest <- file.path(output, "filtered_session_data_latest.rds")
    expect_true(file.exists(latest))
    before <- digest::digest(file = latest, algo = "sha256")
    completed <- runWorkflowLegacyExport(NULL, output, "retry-ack",
        "20260916_160200", serializeFn = function(...) stop("must not serialize"))
    expect_identical(digest::digest(file = latest, algo = "sha256"), before)
    expect_identical(completed$status, "complete")
})

test_that("serialization failure retains previous latest and can restart staging", {
    output <- withr::local_tempdir()
    session <- makeExportSessionTemplate(makeSavedExportProtein(0L),
        makeSavedExportProtein(1L))
    latest <- file.path(output, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    before <- digest::digest(file = latest, algo = "sha256")
    expect_error(runWorkflowLegacyExport(session, output, "retry-serialize",
        "20260916_160300", serializeFn = function(...) stop("serialization boom")),
        "serialization boom")
    expect_identical(digest::digest(file = latest, algo = "sha256"), before)
    expect_false(file.exists(file.path(output,
        "filtered_session_data_20260916_160300.rds")))
    completed <- runWorkflowLegacyExport(session, output, "retry-serialize",
        "20260916_160300")
    expect_identical(completed$status, "complete")
    expect_identical(completed$serialization_count, 1L)
})

test_that("DIA GUI exports serialize once and reuse only their verified owned generation", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "owned-export", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("raw", object, list(), "source")
    manager$saveState("correlation_filtered", object, list(skipped = TRUE), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$design_matrix <- object@design_matrix
    workflow$contrasts_tbl <- data.frame(contrasts = "grouptreated-groupcontrol")
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    serializations <- 0L
    writer <- function(value, path) {
        if (is.list(value) && "r6_complete_states" %in% names(value))
            serializations <<- serializations + 1L
        saveRDS(value, path)
    }
    save <- function(...) saveProtNormExportArtifacts(..., saveRdsFn = writer)
    export <- function(collect = collectProtNormExportSessionData) runProtNormExportSessionWorkflow(
        workflow, norm, input, directory,
        withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
        collectSessionDataFn = collect, saveExportArtifactsFn = save,
        serializeSessionFn = writer)
    head <- manager$.context$head
    before <- diaStateOwnership(manager)
    result <- export()
    expect_identical(serializations, 1L)
    expect_identical(manager$.context$head, head)
    expect_null(result$sessionData)
    expect_identical(diaStateOwnership(manager)$full_reads - before$full_reads, 1L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_gt(manager$.context$export_shared_tables, 0L)
    session <- validateExternalWorkflowLegacySession(result$exportArtifacts$latestFilepath,
        c("initial", "raw", "correlation_filtered"))
    expect_identical(session$current_s4_object, object)
    expect_identical(session$r6_complete_states, manager$states)
    before <- diaStateOwnership(manager)
    reused <- withDIAStateScope(manager, function()
        resolveDIAOwnedExport(manager, result$exportArtifacts$latestFilepath))
    expect_null(reused$r6_complete_states)
    expect_identical(reused$current_s4_object, session$current_s4_object)
    expect_identical(reused$correlation_filtered_s4, session$correlation_filtered_s4)
    expect_identical(reused$design_matrix, session$design_matrix)
    expect_identical(reused$contrasts_tbl, session$contrasts_tbl)
    expect_identical(diaStateOwnership(manager)$full_reads, before$full_reads)
    expect_identical(manager$.context$head, head)
    expect_null(manager$.context$export_pending)
    external <- file.path(directory, "external.rds")
    file.copy(result$exportArtifacts$latestFilepath, external)
    expect_null(resolveDIAOwnedExport(manager, external))
    session$current_s4_object@protein_quant_table$S1[[1L]] <- 99
    saveRDS(session, result$exportArtifacts$latestFilepath)
    expect_null(resolveDIAOwnedExport(manager, result$exportArtifacts$latestFilepath))
    expect_identical(manager$getState(), object)
    manager$revertToState("raw")
    expect_null(resolveDIAOwnedExport(manager, result$exportArtifacts$latestFilepath))
})

test_that("export table sharing preserves attributes and copy-on-modify history isolation", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "sharing", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("first", object, list(), "first")
    object@args$revision <- 1L
    manager$saveState("second", object, list(), "same table, different parameters")
    attr(object@protein_quant_table, "units") <- "different attribute"
    manager$saveState("third", object, list(), "different attribute")
    native <- manager$states
    records <- withDIAStateScope(manager, function()
        materializeDIAStates(manager$.context, share_immutable = TRUE))
    expect_identical(records, native)
    expect_identical(lobstr::obj_addr(records$first$data@protein_quant_table),
        lobstr::obj_addr(records$second$data@protein_quant_table))
    expect_false(identical(lobstr::obj_addr(records$second$data@protein_quant_table),
        lobstr::obj_addr(records$third$data@protein_quant_table)))
    records$first$data@protein_quant_table$S1[[1L]] <- 100
    expect_identical(records$second$data@protein_quant_table$S1[[1L]], 1)
    expect_identical(manager$getState("first")@protein_quant_table$S1[[1L]], 1)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
})

test_that("DIA publication failures retain pending references and retry output work only", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "pending", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("correlation_filtered", object, list(), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    latest <- file.path(directory, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    old_hash <- digest::digest(file = latest, algo = "sha256")
    export <- function(collect = collectProtNormExportSessionData) runProtNormExportSessionWorkflow(
        workflow, norm, input, directory,
        withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
        collectSessionDataFn = collect)
    before <- manager$.context$head
    testthat::with_mocked_bindings({
        expect_error(export(), "injected publication interruption")
    }, workflowExportBoundary = function(phase) {
        if (identical(phase, "after_timestamped_publication")) stop("injected publication interruption")
    })
    expect_identical(digest::digest(file = latest, algo = "sha256"), old_hash)
    expect_identical(manager$.context$head, before)
    expect_error(resolveDIAOwnedExport(manager, latest), class = "saved_entity_outputs_pending")
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    completed <- export(function(...) stop("must not reconstruct history on retry"))
    expect_identical(completed$exportArtifacts$publication$serialization_count, 1L)
    expect_identical(manager$.context$head, before)
    expect_null(manager$.context$export_pending)
    expect_identical(readRDS(latest)$current_s4_object, object)
})

test_that("reference action publication resumes after process loss without repeating science", {
    root <- withr::local_tempdir()
    source <- file.path(root, "source")
    output <- file.path(root, "outputs")
    dir.create(source)
    dir.create(output)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "reference-publication", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- list(source_dir = source)
    manager$saveState("ready", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(), "ready")
    manager$.context$reclaim_enabled <- TRUE
    reference <- getDIAStateReference(manager)
    writeLines("previous a", file.path(output, "a.tsv"))
    writeLines("previous b", file.path(output, "b.tsv"))
    producer <- function(paths) {
        writeLines("accepted a", file.path(paths$protein_qc_dir, "a.tsv"))
        writeLines("accepted b", file.path(paths$protein_qc_dir, "b.tsv"))
        list(status = "ready")
    }
    expect_error(testthat::with_mocked_bindings(
        publishDIAActionOutputs(manager, reference, "normalization", list(protein_qc_dir = output), producer),
        workflowExportBoundary = function(phase) if (phase == "between_action_publications") stop("lost process")),
        "lost process")
    expect_identical(getDIAStateReference(manager), reference)
    expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
    expect_null(manager$.context$reclaim_pending)
    request <- manager$.context$snapshot$action_outputs[[reference$entity_id]]
    backup <- file.path(source, ".workflow-export", request$request_id, "previous", "protein_qc_dir", "a.tsv")
    expect_identical(readLines(backup), "previous a")
    closeDIAStateAdapter(manager)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "reference-publication", "DIA"), reopen = TRUE)
    manager$.context$experiment_paths <- list(source_dir = source)
    manager$.context$reclaim_enabled <- TRUE
    result <- testthat::with_mocked_bindings(
        publishDIAActionOutputs(manager, reference, "normalization", list(protein_qc_dir = output),
            function(...) stop("output rendering repeated")),
        commitDIAReferenceState = function(...) stop("science repeated"),
        materializeDIAStateRecord = function(...) stop("assay hydrated"))
    expect_identical(result, list(status = "ready"))
    expect_identical(getDIAStateReference(manager), reference)
    expect_identical(readLines(file.path(output, "a.tsv")), "accepted a")
    expect_identical(readLines(file.path(output, "b.tsv")), "accepted b")
    expect_silent(assertDIAExportReady(manager))
    expect_identical(manager$.context$store$open_resources, 0L)
    expect_identical(manager$.context$reclaim_pending$request, request$request_id)
    testthat::with_mocked_bindings(reclaimDIACompletedMemory(manager),
        reclaimWorkflowMemory = function() list(status = "released"))
    expect_null(manager$.context$reclaim_pending)
    publishDIAActionOutputs(manager, reference, "normalization", list(protein_qc_dir = output),
        function(...) stop("completed output was rendered again"))
    expect_null(manager$.context$reclaim_pending)
    expect_identical(getDIAStateReference(manager), reference)
})

test_that("export shares unchanged columns and support tables without changing history", {
    skip_if_not(exists("sameDIAExportValue", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 sharing trial; retained for its installed candidates")
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "columns", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    support <- data.frame(protein = c("P1", "P2"), value = c(NA_real_, Inf),
        group = ordered(c("b", "a"), levels = c("b", "a", "unused")))
    attr(support, "meaning") <- "unchanged scientific support"
    sections <- c("peptideIntensityFiltering", "removePeptidesWithOnlyOneReplicate",
        "peptideMissingValueImputation")
    for (section in sections) object@args[[section]] <- list(support_table = support)
    manager$saveState("first", object, list(), "first")
    object@protein_quant_table$S1[[1L]] <- 9
    manager$saveState("second", object, list(), "different intensity")
    attr(object@args$peptideIntensityFiltering$support_table, "meaning") <- "different meaning"
    manager$saveState("third", object, list(), "different support attribute")
    expected <- manager$states
    for (iteration in 1:3) {
        reads <- manager$.context$hydrates
        records <- withDIAStateScope(manager, function() {
            records <- materializeDIAStates(manager$.context, share_immutable = TRUE)
            for (name in names(records)) {
                reference <- manager$.context$snapshot$states[[name]]
                expect_identical(rlang::obj_address(get(reference$entity_id,
                    manager$.context$working, inherits = FALSE)), rlang::obj_address(records[[name]]))
            }
            records
        })
        expect_true(sameDIAExportValue(records, expected))
        expect_identical(manager$.context$hydrates - reads, as.integer(length(records)))
        expect_identical(rlang::obj_address(records$first$data@protein_quant_table$S2),
            rlang::obj_address(records$second$data@protein_quant_table$S2))
        expect_false(identical(rlang::obj_address(records$first$data@protein_quant_table$S1),
            rlang::obj_address(records$second$data@protein_quant_table$S1)))
        for (section in sections) {
            expect_identical(rlang::obj_address(records$first$data@args[[section]]$support_table),
                rlang::obj_address(records$second$data@args[[section]]$support_table))
        }
        expect_gt(manager$.context$export_shared_columns, 0L)
        expect_gte(manager$.context$export_shared_support_tables, 5L)
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        records$first$data@args$peptideIntensityFiltering$support_table$value[[1L]] <- 88
        records$first$data@protein_quant_table$S2[[1L]] <- 88
        expect_true(is.na(records$second$data@args$peptideIntensityFiltering$support_table$value[[1L]]))
        expect_identical(records$second$data@protein_quant_table$S2[[1L]], 2)
        path <- tempfile(fileext = ".rds")
        saveRDS(records, path)
        expect_true(sameDIAExportValue(readRDS(path), records))
    }
    expect_identical(manager$revertToState("first"), expected$first$data)
    expect_identical(diaStateOwnership(manager)$store$open_resources, 0L)
})

test_that("export equality rejects signed-zero, missingness, encoding and attribute differences", {
    skip_if_not(exists("sameDIAExportValue", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 sharing trial; retained for its installed candidates")
    nested <- list(named = structure(list(value = 1), units = "signal", source = list(step = "QC")))
    expect_true(sameDIAExportValue(nested, nested))
    expect_true(sameDIAExportValue(data.frame(x = 1), data.frame(x = 1)))
    expect_false(sameDIAExportValue(c(0, 1), c(-0, 1)))
    expect_false(sameDIAExportValue(NA_real_, NaN))
    expect_false(sameDIAExportValue(1L, 1))
    expect_false(sameDIAExportValue(factor("a", levels = c("a", "b")), factor("a")))
    utf8 <- enc2utf8("\u00e9")
    latin1 <- iconv(utf8, from = "UTF-8", to = "latin1")
    expect_identical(utf8, latin1)
    expect_false(sameDIAExportValue(utf8, latin1))
    expect_false(sameDIAExportValue(structure(1, units = utf8), structure(1, units = latin1)))
    expect_false(sameDIAExportValue(structure(1, first = TRUE, second = 2),
        structure(1, second = 2, first = TRUE)))
    expect_true(sameDIAExportValue(rep(utf8, 65537L), rep(utf8, 65537L)))
    expect_true(sameDIAExportValue(c(NA, NaN, Inf, -Inf, -0), c(NA, NaN, Inf, -Inf, -0)))
})

test_that("export sharing skips mutable and unsupported tables and preserves containers", {
    skip_if_not(exists("sameDIAExportValue", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 sharing trial; retained for its installed candidates")
    pool <- new.env(parent = emptyenv())
    pool$tables <- pool$columns <- list()
    pool$shared_tables <- pool$shared_columns <- pool$shared_support_tables <- 0L
    original <- tibble::tibble(id = c("a", "b"), value = c(-0, NA_real_))
    attr(original, "units") <- "signal"
    first <- shareDIAExportTable(original, "assay", pool)
    second <- original
    second$added <- 1:2
    shared <- shareDIAExportTable(second, "assay", pool)
    expect_true(sameDIAExportValue(shared, second))
    expect_identical(rlang::obj_address(first$value), rlang::obj_address(shared$value))
    changed <- second
    changed$value[[1L]] <- 0
    changed <- shareDIAExportTable(changed, "assay", pool)
    expect_identical(1 / changed$value[[1L]], Inf)
    expect_identical(1 / first$value[[1L]], -Inf)
    owner <- new.env(parent = emptyenv())
    mutable <- original
    attr(mutable, "owner") <- owner
    before <- pool$shared_tables
    expect_identical(shareDIAExportTable(mutable, "assay", pool), mutable)
    unsupported <- original
    class(unsupported) <- c("custom_table", class(unsupported))
    expect_identical(shareDIAExportTable(unsupported, "assay", pool), unsupported)
    expect_identical(pool$shared_tables, before)
})

test_that("DIA export collects once before fresh hydration and never on output-only retry", {
    skip_if(!"gcFn" %in% names(formals(runDIAOwnedSessionExport)),
        "rejected pre-export collection is retained only for its installed trial build")
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "export-gc", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("raw", object, list(), "raw")
    manager$saveState("correlation_filtered", object, list(), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    events <- character()
    head <- manager$.context$head
    export <- function() withDIAStateScope(manager, function() runDIAOwnedSessionExport(
        manager, workflow, norm, input, directory,
        collectFn = function(...) {
            events <<- c(events, "hydrate")
            collectProtNormExportSessionData(...)
        },
        saveFn = function(...) {
            events <<- c(events, "save")
            saveProtNormExportArtifacts(...)
        },
        gcFn = function(full) {
            expect_true(full)
            expect_identical(manager$.context$full_reads, 0L)
            expect_identical(diaStateOwnership(manager)$working_records, 0L)
            events <<- c(events, "gc")
        }))
    testthat::with_mocked_bindings({
        expect_error(export(), "publication interrupted")
    }, workflowExportBoundary = function(phase) {
        if (identical(phase, "after_timestamped_publication")) stop("publication interrupted")
    })
    expect_identical(events, c("gc", "hydrate", "save"))
    expect_identical(manager$.context$head, head)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    result <- export()
    expect_identical(events, c("gc", "hydrate", "save", "save"))
    expect_identical(manager$.context$full_reads, 1L)
    expect_identical(manager$.context$head, head)
    expect_identical(result$exportArtifacts$publication$serialization_count, 1L)
    session <- validateExternalWorkflowLegacySession(result$exportArtifacts$latestFilepath,
        c("initial", "raw", "correlation_filtered"))
    expect_identical(session$r6_complete_states, manager$states)
    expect_identical(session$current_s4_object, object)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$store$open_resources, 0L)
})

test_that("DIA history and working-map roots expire before the first legacy copy", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "release-export", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("raw", object, list(), "raw")
    manager$saveState("correlation_filtered", object, list(), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    history_ref <- cache_ref <- NULL
    serialized <- copied <- 0L
    head <- manager$.context$head
    collector <- function(...) {
        session <- collectProtNormExportSessionData(...)
        history_owner <- new.env(parent = emptyenv())
        history_owner$payload <- numeric(100000L)
        session$r6_complete_states$raw$data@args$export_probe <- history_owner
        history_ref <<- rlang::new_weakref(history_owner)
        cache_owner <- new.env(parent = emptyenv())
        cache_owner$payload <- numeric(100000L)
        key <- manager$.context$snapshot$states$raw$entity_id
        record <- get(key, manager$.context$working, inherits = FALSE)
        record$data@args$cache_probe <- cache_owner
        assign(key, record, manager$.context$working)
        cache_ref <<- rlang::new_weakref(cache_owner)
        session
    }
    writer <- function(value, path) {
        serialized <<- serialized + 1L
        expect_true(is.environment(rlang::wref_key(history_ref)))
        expect_true(is.environment(rlang::wref_key(cache_ref)))
        saveRDS(value, path)
    }
    copy_file <- function(from, to, overwrite) {
        copied <<- copied + 1L
        expect_null(rlang::wref_key(history_ref))
        expect_null(rlang::wref_key(cache_ref))
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        expect_identical(manager$.context$head, head)
        file.copy(from, to, overwrite = overwrite)
    }
    original_export <- runWorkflowLegacyExport
    testthat::local_mocked_bindings(runWorkflowLegacyExport = function(...) {
        original_export(..., copyFn = copy_file)
    })
    result <- runProtNormExportSessionWorkflow(workflow, norm, input, directory,
        withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
        collectSessionDataFn = collector, serializeSessionFn = writer)
    expect_identical(serialized, 1L)
    expect_identical(copied, 1L)
    expect_null(result$sessionData)
    expect_identical(manager$.context$full_reads, 1L)
    session <- validateExternalWorkflowLegacySession(result$exportArtifacts$latestFilepath,
        c("initial", "raw", "correlation_filtered"))
    expect_identical(session$current_s4_object, object)
    expect_identical(session$r6_complete_states$raw$data@args$export_probe$payload, numeric(100000L))
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
})

test_that("DIA partial latest-copy retry never reconstructs or serializes history", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "copy-retry", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("correlation_filtered", object, list(), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    latest <- file.path(directory, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    previous <- digest::digest(file = latest, algo = "sha256")
    head <- manager$.context$head
    serializations <- 0L
    writer <- function(value, path) {
        serializations <<- serializations + 1L
        saveRDS(value, path)
    }
    export <- function(collect = collectProtNormExportSessionData, serialize = writer) {
        runProtNormExportSessionWorkflow(workflow, norm, input, directory,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
            collectSessionDataFn = collect, serializeSessionFn = serialize)
    }
    original_export <- runWorkflowLegacyExport
    testthat::with_mocked_bindings({
        expect_error(export(), class = "saved_entity_write")
    }, runWorkflowLegacyExport = function(...) original_export(..., copyFn = function(from, to, overwrite) {
        writeBin(as.raw(c(1L, 2L)), to)
        FALSE
    }))
    request <- manager$.context$export_pending
    paths <- workflowExportPaths(directory, request$request_id)
    journal <- jsonlite::read_json(file.path(paths$request_dir, "journal.json"))
    expect_null(journal$latest_sha256)
    expect_identical(journal$status, "outputs_pending")
    expect_gt(journal$session_bytes, file.size(journal$staged_latest))
    expect_identical(serializations, 1L)
    expect_identical(manager$.context$full_reads, 1L)
    expect_identical(manager$.context$head, head)
    expect_identical(digest::digest(file = latest, algo = "sha256"), previous)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_error(publishWorkflowLegacyExport(paths), class = "saved_entity_outputs_pending")
    expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
    result <- export(function(...) stop("unexpected hydration"), function(...) stop("unexpected serialization"))
    expect_identical(result$exportArtifacts$publication$serialization_count, 1L)
    expect_identical(result$exportArtifacts$publication$scientific_checkpoint_count, 0L)
    expect_identical(manager$.context$full_reads, 1L)
    expect_identical(manager$.context$head, head)
    expect_null(manager$.context$export_pending)
    expect_identical(readRDS(latest)$current_s4_object, object)
})

test_that("serialized-stage recovery validates its owned source before publication", {
    directory <- withr::local_tempdir()
    paths <- workflowExportPaths(directory, "serialized-recovery")
    session <- makeExportSessionTemplate(makeSavedExportProtein(), makeSavedExportProtein(1L))
    latest <- file.path(directory, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    previous <- digest::digest(file = latest, algo = "sha256")
    expect_error(stageWorkflowLegacyExport(session, paths, "serialized-recovery", "20260922_180000",
        serialize_only = TRUE, boundaryFn = function(phase) {
            if (identical(phase, "after_serialized_checkpoint")) stop("interrupted after serialization")
        }), "interrupted after serialization")
    expect_true(file.exists(file.path(paths$request_dir, "journal.json")))
    expect_false(file.exists(file.path(directory, "export_complete_serialized-recovery.json")))
    journal <- jsonlite::read_json(file.path(paths$request_dir, "journal.json"))
    connection <- file(journal$staged_session, "r+b")
    first <- readBin(connection, "raw", 1L)
    seek(connection, 0, rw = "write")
    writeBin(as.raw(bitwXor(as.integer(first), 1L)), connection)
    close(connection)
    expect_error(runWorkflowLegacyExport(NULL, directory, "serialized-recovery", "20260922_180000"),
        class = "saved_entity_corrupt")
    expect_identical(digest::digest(file = latest, algo = "sha256"), previous)
    connection <- file(journal$staged_session, "r+b")
    writeBin(first, connection)
    close(connection)
    result <- runWorkflowLegacyExport(NULL, directory, "serialized-recovery", "20260922_180000",
        serializeFn = function(...) stop("must not serialize"))
    expect_identical(result$status, "complete")
    expect_identical(readRDS(result$latest_path), session)
})

test_that("DIA serialization failures cannot certify or retain a partial history", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "serialize-failure", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeSavedExportProtein()
    manager$saveState("correlation_filtered", object, list(), "final")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "none", ruv_mode = "skip")
    directory <- withr::local_tempdir()
    latest <- file.path(directory, "filtered_session_data_latest.rds")
    saveRDS(list(previous = TRUE), latest)
    previous <- digest::digest(file = latest, algo = "sha256")
    head <- manager$.context$head
    export <- function(collect = collectProtNormExportSessionData, writer = saveRDS) {
        runProtNormExportSessionWorkflow(workflow, norm, input, directory,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
            collectSessionDataFn = collect, serializeSessionFn = writer)
    }
    expect_error(export(function(...) list(r6_complete_states = list())), class = "saved_entity_corrupt")
    expect_error(export(writer = function(...) stop("serialization failed")), "serialization failed")
    request <- manager$.context$export_pending
    paths <- workflowExportPaths(directory, request$request_id)
    expect_false(file.exists(file.path(paths$request_dir, "journal.json")))
    expect_identical(manager$.context$head, head)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(digest::digest(file = latest, algo = "sha256"), previous)
    result <- export()
    expect_identical(result$exportArtifacts$publication$status, "complete")
    expect_identical(readRDS(latest)$current_s4_object, object)
    expect_identical(manager$.context$head, head)
})

test_that("ten small normalization revert style export cycles release action resources", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "cycles", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    expected <- normaliseBetweenSamples(original, normalisation_method = "cyclicloess")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$design_matrix <- original@design_matrix
    workflow$contrasts_tbl <- data.frame(contrasts = "grouptreated-groupcontrol")
    output <- withr::local_tempdir()
    for (cycle in seq_len(10L)) {
        manager$revertToState("protein_replicate_filtered")
        normalized <- normaliseBetweenSamples(manager$getState(), normalisation_method = "cyclicloess")
        expect_identical(normalized@protein_quant_table, expected@protein_quant_table)
        manager$saveState("correlation_filtered", normalized, list(skipped = TRUE), "final")
        product <- buildProtNormQcPlotProduct(normalized,
            list(color_var = "group", shape_var = NULL), "correlation")
        before <- product$data
        styled <- styleWorkflowDiaQcPlot(product, list(title = paste("Cycle", cycle), font_size = 10))
        expect_identical(product$data, before)
        expect_warning(ggplot2::ggsave(file.path(output, paste0("cycle-", cycle, ".png")), styled,
            width = 4, height = 3, dpi = 100), "Removed 1 row")
        norm <- list(correlation_filtered_obj = normalized, best_k = NA_integer_)
        result <- runProtNormExportSessionWorkflow(workflow, norm,
            list(norm_method = "cyclicloess", ruv_mode = "skip"), output,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
        expect_identical(readRDS(result$exportArtifacts$latestFilepath)$current_s4_object, normalized)
        info <- diaStateOwnership(manager)
        expect_identical(info$working_records, 0L)
        expect_identical(info$native_records, 0L)
        expect_identical(info$store$open_resources, 0L)
        expect_equal(info$store$bytes, 0)
        expect_lte(info$store$metadata_bytes, 64 * 1024^2)
        expect_length(manager$.context$owned_workers, 0L)
        expect_null(manager$.context$export_pending)
    }
})

test_that("native artifact reuse skips copying and failed replacement preserves prior bytes", {
    root <- withr::local_tempdir()
    source <- file.path(root, "source")
    dir.create(source)
    name <- "filtered_session_20260925_000000_native-export-copy-fixture.msdia.tar"
    archive <- file.path(source, name)
    # This fixture exercises the completed-file copy boundary, not archive import.
    writeBin(serialize(list(generation = "current", values = 1:100), NULL), archive)
    pointer <- list(schema = "multischolar.dia_native_latest.v1", status = "complete",
        archive = name, bytes = file.size(archive), sha256 = digest::digest(file = archive, algo = "sha256"))
    jsonlite::write_json(pointer, file.path(source, "filtered_session_latest.json"), auto_unbox = TRUE)
    destination <- file.path(root, "integration.msdia.tar")
    expect_identical(copyDIANativeSessionArtifact(source, destination), destination)
    expect_identical(digest::digest(file = destination, algo = "sha256"), pointer$sha256)
    testthat::with_mocked_bindings(
        expect_identical(copyDIANativeSessionArtifact(source, destination), destination),
        runDIANativeIO = function(...) stop("unchanged copy was repeated"))
    prior <- serialize(list(generation = "previous", values = 1:50), NULL)
    writeBin(prior, destination)
    testthat::with_mocked_bindings(
        expect_error(copyDIANativeSessionArtifact(source, destination), "copy failed"),
        runDIANativeIO = function(...) stop("copy failed"))
    expect_identical(readBin(destination, "raw", n = file.size(destination)), prior)
    expect_identical(copyDIANativeSessionArtifact(source, destination), destination)
    expect_identical(digest::digest(file = destination, algo = "sha256"), pointer$sha256)
})


test_that("owned snapshots remain independent of later source writes", {
    root <- withr::local_tempdir()
    source <- file.path(root, "source")
    destination <- file.path(root, "copy")
    writeBin(charToRaw("original"), source)
    snapshotWorkflowFile(source, destination, file.path(root, "copy.log"))
    writeBin(charToRaw("modified"), source)
    expect_identical(readBin(destination, "raw", n = 8L), charToRaw("original"))
    expect_error(snapshotWorkflowFile(source, destination, file.path(root, "copy.log")),
        class = "saved_entity_conflict")
    dangling <- file.path(root, "dangling")
    target <- file.path(root, "missing-target")
    expect_true(file.symlink(target, dangling))
    expect_error(snapshotWorkflowFile(source, dangling, file.path(root, "copy.log")),
        class = "saved_entity_conflict")
    expect_identical(Sys.readlink(dangling), target)
    expect_false(file.exists(target))
})
