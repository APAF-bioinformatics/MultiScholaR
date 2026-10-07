# Installed WFX-004 preflight; explicit runner supplies fixture helper and bound root.
makeWFXAdapterManager <- function(object, root) {
    dir.create(root)
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "adapter fixture")
    manager
}

readWFXAdapterTable <- function(manager, reference) {
    scope <- newWorkflowTableReadScope(manager$.context$store, reference)
    on.exit(releaseWorkflowTableReadScope(scope))
    if (!scope$metadata$rows) return(readSavedWorkflowTablePrototype(scope))
    stopifnot(scope$metadata$rows <= 10000L)
    materializeDIATable(manager$.context, reference, scope = scope)
}

test_that("compute-only sample filtering matches frozen native science", {
    skip_if(Sys.getenv("MULTISCHOLAR_WFX004_PREFLIGHT") != "science", "explicit installed preflight")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    fixture <- makeWorkflowExecutionFilterFixture()
    base <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
    base@peptide_data <- fixture$input
    base@peptide_data$Q.Value <- rep(0, nrow(fixture$input))
    base@q_value_column <- "Q.Value"
    base@raw_quantity_column <- base@norm_quantity_column <- "Quantity"
    base@sample_id <- "Run"
    base@protein_id_column <- "Protein.Group"
    base@peptide_sequence_column <- "Stripped.Sequence"
    base@design_matrix <- rbind(fixture$design[4:1, ], data.frame(Run = "extra", group = "case"))
    cases <- list(mixed = list(object = base, cutoff = 2L, inclusion = fixture$inclusion),
        all = list(object = base, cutoff = 1L, inclusion = character()),
        none = list(object = base, cutoff = 5L, inclusion = character()))
    empty <- base
    empty@peptide_data <- base@peptide_data[0, ]
    cases$empty <- list(object = empty, cutoff = 2L, inclusion = fixture$inclusion)
    wide <- base
    wide@peptide_data$Protein.Ids <- rep(strrep("A", 2048L), nrow(wide@peptide_data))
    cases$wide <- list(object = wide, cutoff = 2L, inclusion = fixture$inclusion)
    for (label in names(cases)) {
        item <- cases[[label]]
        native <- filterMinNumPeptidesPerSample(item$object,
            peptides_per_sample_cutoff = item$cutoff, inclusion_list = item$inclusion, core_utilisation = NA)
        manager <- makeWFXAdapterManager(item$object, file.path(root, paste0("filter-", label)))
        before <- workflowExecutionHead(manager$.context$store)
        value <- computeDIAReferenceSampleFilter(manager, getDIAStateReference(manager),
            list(peptides_per_sample_cutoff = item$cutoff, inclusion_list = item$inclusion))
        expect_identical(readWFXAdapterTable(manager, value$reference), native@peptide_data, info = label)
        args <- native@args$filterMinNumPeptidesPerSample
        expect_identical(value$summary, args$filter_summary, info = label)
        for (field in c("support_table", "removal_ledger")) expect_identical(
            readWFXAdapterTable(manager, value$checkpoint$support$filterMinNumPeptidesPerSample[[field]]),
            args[[field]], info = paste(label, field))
        expect_identical(workflowExecutionHead(manager$.context$store), before)
        expect_identical(manager$.context$store$open_resources, 0L)
        closeDIAStateAdapter(manager)
    }
    for (kind in c("missing-key", "blank-key", "duplicate-design", "blank-design", "unmapped-design")) {
        invalid <- base
        if (kind == "missing-key") invalid@peptide_data$Protein.Group[1] <- NA_character_
        if (kind == "blank-key") invalid@peptide_data$Stripped.Sequence[1] <- " "
        if (kind == "duplicate-design") invalid@design_matrix <- rbind(base@design_matrix, base@design_matrix[1, ])
        if (kind == "blank-design") invalid@design_matrix$Run[1] <- " "
        if (kind == "unmapped-design") invalid@design_matrix <- base@design_matrix[-1, ]
        manager <- makeWFXAdapterManager(invalid, file.path(root, paste0("invalid-", kind)))
        before <- workflowExecutionHead(manager$.context$store)
        expect_error(computeDIAReferenceSampleFilter(manager, getDIAStateReference(manager),
            list(peptides_per_sample_cutoff = 2L, inclusion_list = fixture$inclusion)))
        expect_identical(workflowExecutionHead(manager$.context$store), before)
        expect_identical(manager$.context$store$open_resources, 0L)
        closeDIAStateAdapter(manager)
    }
    manager <- makeWFXAdapterManager(base, file.path(root, "invalid-threshold"))
    for (k in list(0, -1, 1.5, NA_real_, Inf)) expect_error(computeDIAReferenceSampleFilter(
        manager, getDIAStateReference(manager), list(peptides_per_sample_cutoff = k,
            inclusion_list = fixture$inclusion)), "positive integer")
    closeDIAStateAdapter(manager)
})

test_that("compute-only P14 preserves native quantities and head ownership", {
    skip_if(Sys.getenv("MULTISCHOLAR_WFX004_PREFLIGHT") != "science", "explicit installed preflight")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    base <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
    manager <- makeWFXAdapterManager(base, file.path(root, "rollup-small"))
    on.exit(closeDIAStateAdapter(manager))
    before <- workflowExecutionHead(manager$.context$store)
    native <- rollUpPrecursorToPeptide(base, core_utilisation = NA)
    value <- computeDIAReferencePrecursorRollup(manager, getDIAStateReference(manager))
    expect_identical(readWFXAdapterTable(manager, value$reference), native@peptide_data)
    expect_identical(value$summary, native@args$rollUpPrecursorToPeptide$rollup_summary)
    expect_identical(workflowExecutionHead(manager$.context$store), before)
    expect_identical(manager$.context$store$open_resources, 0L)
})

test_that("both scientific adapters use guarded controller acceptance and publication", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "adapters", "installed managed preflight")
    app <- workflow_test_application
    skip_if(!isTRUE(app$config$policy$diagnostic), "timed member has a separate fixed entry")
    withr::local_options(multischolar.workflow_execution_boundary = NULL)
    await <- function(owner) {
        deadline <- proc.time()[["elapsed"]] + 120
        repeat {
            value <- pollWorkflowAction(owner)
            if (!value$status %in% c("admitting", "running")) return(value)
            if (proc.time()[["elapsed"]] > deadline) stop("adapter lifecycle exceeded bounded wait")
            Sys.sleep(0.02)
        }
    }
    for (operation in c("precursor_rollup", "minimum_peptides_per_sample")) {
        object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
        parameters <- list()
        if (operation == "minimum_peptides_per_sample") {
            fixture <- makeWorkflowExecutionFilterFixture()
            object@peptide_data <- fixture$input
            object@peptide_data$Q.Value <- rep(0, nrow(fixture$input))
            object@q_value_column <- "Q.Value"
            object@raw_quantity_column <- object@norm_quantity_column <- "Quantity"
            object@sample_id <- "Run"
            object@protein_id_column <- "Protein.Group"
            object@peptide_sequence_column <- "Stripped.Sequence"
            object@design_matrix <- fixture$design
            parameters <- list(peptides_per_sample_cutoff = 2L, inclusion_list = fixture$inclusion)
        }
        manager <- makeWFXAdapterManager(object, file.path(app$config$root, operation))
        owner <- newWorkflowExecutionController(manager, app)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        native <- if (operation == "precursor_rollup") rollUpPrecursorToPeptide(object, core_utilisation = NA) else
            filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 2L,
                inclusion_list = parameters$inclusion_list, core_utilisation = NA)
        id <- runWorkflowAction(operation, list(state = getDIAStateReference(manager)), parameters, owner)
        paths <- owner$job$paths
        if (operation == "precursor_rollup") options(multischolar.workflow_execution_boundary = function(stage) {
            if (stage == "publication") stop("injected adapter publication failure")
        })
        done <- await(owner)
        options(multischolar.workflow_execution_boundary = NULL)
        expect_identical(done$status, if (operation == "precursor_rollup") "publication_pending" else "complete")
        worker_file <- file.path(paths$request_dir, "worker-result.rds")
        worker_hash <- digest::digest(file = worker_file, algo = "sha256")
        accepted <- getDIAStateReference(manager)
        result <- resumeWorkflowAction(owner, id)
        expect_identical(getDIAStateReference(manager), accepted)
        expect_identical(digest::digest(file = worker_file, algo = "sha256"), worker_hash)
        expect_identical(acceptWorkflowActionResult(owner, readRDS(worker_file))$status, "duplicate")
        expect_identical(readWFXAdapterTable(manager, result$reference), native@peptide_data)
        design <- readDIAExecutionComponents(manager, fields = "design_matrix")@design_matrix
        expect_identical(design, native@design_matrix)
        expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
        expect_null(app$active)
        expect_identical(manager$.context$store$open_resources, 0L)
        files <- list.files(paths$request_dir, recursive = TRUE, full.names = TRUE)
        figures <- files[grepl("[.](png|pdf)$", files) & !grepl("/figures/", files)]
        expect_length(figures, 12L)
        expect_true(all(file.size(figures) > 0))
        closeDIAStateAdapter(manager)
    }
})


test_that("admission estimates and reader boundaries preserve scientific membership", {
    skip_if(Sys.getenv("MULTISCHOLAR_WFX004_PREFLIGHT") != "domain", "explicit installed domain preflight")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    fixture <- makeWorkflowExecutionFilterFixture()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
    object@peptide_data <- fixture$input
    object@peptide_data$Q.Value <- rep(0, nrow(fixture$input))
    object@q_value_column <- "Q.Value"
    object@raw_quantity_column <- object@norm_quantity_column <- "Quantity"
    object@sample_id <- "Run"
    object@protein_id_column <- "Protein.Group"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@design_matrix <- fixture$design
    manager <- makeWFXAdapterManager(object, file.path(root, "domain-filter"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    reference <- getDIAStateReference(manager)
    params <- list(peptides_per_sample_cutoff = 2L, inclusion_list = fixture$inclusion)
    expected <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 2L,
        inclusion_list = fixture$inclusion, core_utilisation = NA)
    sys.source(file.path(workflow_test_root, "tests/testthat/helper-workflow-bulk-rollup.R"), environment())
    control <- computeWFXNativeFilterControl(manager, reference, params)
    expect_identical(readWFXAdapterTable(manager, control$reference), expected@peptide_data)
    expect_identical(control$summary, expected@args$filterMinNumPeptidesPerSample$filter_summary)
    expect_error(computeDIAReferenceSampleFilter(manager, reference, params, worker_cap = 1),
        class = "workflow_execution_unsupported")
    stream <- streamSavedWorkflowTable
    for (size in c(1L, 2L, 3L, 7L)) {
        value <- testthat::with_mocked_bindings(computeDIAReferenceSampleFilter(manager, reference, params),
            streamSavedWorkflowTable = function(scope, columns, consume, ...) stream(scope, columns,
                function(data, offset, row_names) {
                    for (at in split(seq_len(nrow(data)), ceiling(seq_len(nrow(data)) / size)))
                        consume(data[at, , drop = FALSE], offset + min(at) - 1L, row_names[at])
                }, ...), .package = "MultiScholaR")
        actual <- readWFXAdapterTable(manager, value$reference)
        expect_identical(actual, expected@peptide_data)
        expect_identical(actual$source_ordinal, fixture$expected$retained_source_ordinal)
        expect_gt(value$execution$admission$estimated_bytes, 0)
    }
    dropped <- computeDIAReferenceSampleFilter(manager, reference,
        list(peptides_per_sample_cutoff = 2L, inclusion_list = character()))
    expect_false(identical(readWFXAdapterTable(manager, dropped$reference)$source_ordinal,
        fixture$expected$retained_source_ordinal))
    # Mutate the count only, preserving the original input rows and inclusion list.
    query <- DBI::dbGetQuery
    wrong <- testthat::with_mocked_bindings(computeDIAReferenceSampleFilter(manager, reference, params),
        dbGetQuery = function(conn, statement, ...) {
            value <- query(conn, statement, ...)
            if (grepl("AS identities", statement, fixed = TRUE))
                value$identities[startsWith(value$sample, "A_")] <- 0L
            value
        }, .package = "DBI")
    expect_false(identical(readWFXAdapterTable(manager, wrong$reference)$source_ordinal,
        fixture$expected$retained_source_ordinal))
    expect_identical(getDIAStateReference(manager), reference)
    expect_identical(manager$.context$store$open_resources, 0L)

    base <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small", workflow_test_root))
    base@peptide_data$Protein.Ids[1:2] <- c("P1;P2", "P3;P4")
    rollup <- makeWFXAdapterManager(base, file.path(root, "domain-rollup"))
    on.exit(closeDIAStateAdapter(rollup), add = TRUE)
    input <- getDIAStateReference(rollup)
    expect_error(computeDIAReferencePrecursorRollup(rollup, input, worker_cap = 1),
        class = "workflow_execution_unsupported")
    value <- computeDIAReferencePrecursorRollup(rollup, input)
    native <- rollUpPrecursorToPeptide(base, core_utilisation = NA)
    expect_identical(readWFXAdapterTable(rollup, value$reference), native@peptide_data)
    expect_gt(value$execution$admission$output_groups, 0)
    expect_gte(value$execution$admission$annotation_width,
        max(nchar(native@peptide_data$Protein.Ids, type = "bytes")))
    expect_identical(getDIAStateReference(rollup), input)

    oversized <- base
    oversized@peptide_data <- base@peptide_data[rep(1L, 2000L), , drop = FALSE]
    oversized@peptide_data$Precursor.Id <- paste0("large-group-", seq_len(2000L))
    oversized@peptide_data$Modified.Sequence <- rep(strrep("A", 40000L), 2000L)
    large <- makeWFXAdapterManager(oversized, file.path(root, "oversized-group"))
    on.exit(closeDIAStateAdapter(large), add = TRUE)
    input <- getDIAStateReference(large)
    expect_error(computeDIAReferencePrecursorRollup(large, input), "indivisible rollup group")
    expect_identical(getDIAStateReference(large), input)
    expect_identical(large$.context$store$open_resources, 0L)
})


test_that("frozen WFX-004 comparison member emits science and complete-operation evidence", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "adapters", "explicit installed member")
    app <- workflow_test_application
    skip_if(isTRUE(app$config$policy$diagnostic), "production guard required for timing")
    request <- readRDS(file.path(app$config$root, "timed-member.rds"))
    skip_if(identical(request$programme, "perf106-repair-v2"), "PERF-106 v2 has its own bound member")
    stopifnot(identical(request$schema, 1L), request$cell %in% c("primary", "filter"),
        request$arm %in% c("control", "candidate"), request$rep %in% 1:2)
    sys.source(file.path(workflow_test_root, "tests/testthat/helper-workflow-bulk-rollup.R"), environment())
    manifest_path <- file.path(workflow_test_root, "tests/testdata/workflow-execution/qualification.json")
    stopifnot(identical(digest::digest(file = manifest_path, algo = "sha256"), request$manifest_sha256))
    manifest <- jsonlite::read_json(manifest_path)
    owned <- file.path(app$config$root, "project")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(request$source_store, owned)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    if (manager$current_state != "protein_peptide_filtered")
        selectDIAReplayCheckpoint(manager, "protein_peptide_filtered")
    record <- readDIAStateMetadata(manager)
    stopifnot(identical(record$data$components$peptide_data$entity_id, request$input_entity_id))
    owner <- newWorkflowExecutionController(manager, app)
    original <- digest::digest(file = file.path(request$source_store, "dia-state-head.rds"), algo = "sha256")
    stopifnot(identical(original, request$source_head_sha256))
    operation <- if (request$cell == "primary") "precursor_rollup" else "minimum_peptides_per_sample"
    receipt <- list(status = "running", request = request, output = NULL, science = FALSE,
        build = workflowExecutionBuildIdentity(), promotion_authorized = FALSE)
    path <- file.path(app$config$root, "member-result.rds")
    on.exit(saveRDS(receipt, path), add = TRUE)
    id <- tryCatch(runWorkflowAction(operation, list(state = getDIAStateReference(manager)), request$parameters, owner),
        error = function(error) {
            receipt$status <<- if (inherits(error, "workflow_execution_refused")) "admission_refused" else "dispatch_failed"
            receipt$error <<- conditionMessage(error)
            stop(error)
        })
    repeat {
        done <- pollWorkflowAction(owner)
        if (!done$status %in% c("admitting", "running")) break
        Sys.sleep(0.02)
    }
    receipt$execution <- done
    receipt$status <- done$status
    expect_identical(done$status, "complete")
    if (done$status != "complete") return(invisible(NULL))
    result <- resumeWorkflowAction(owner, id)
    receipt$output <- result
    # Independent bounded oracle work starts only after the action window closes.
    if (request$cell == "primary") {
        fingerprint <- fingerprintDIAWriterOutput(owned, result$reference)
        expected <- unlist(manifest$receipts$primary_native$output_hashes, use.names = FALSE)
        stopifnot(length(expected) == 1200L)
        expect_identical(fingerprint, expected)
        expect_equal(result$summary, manifest$receipts$primary_native$summary)
        receipt$fingerprints <- fingerprint
        receipt$science <- identical(fingerprint, expected) &&
            isTRUE(all.equal(result$summary, manifest$receipts$primary_native$summary))
    } else {
        fixture <- makeWorkflowExecutionFilterFixture(1000L)
        fixture$input$Q.Value <- rep(0, nrow(fixture$input))
        expected <- filterMinNumPeptidesPerSampleHelper(fixture$input,
            peptides_per_sample_cutoff = 2L, inclusion_list = fixture$inclusion,
            sample_id_column = "Run", protein_id_column = "Protein.Group",
            peptide_sequence_column = "Stripped.Sequence", design_matrix = fixture$design,
            design_sample_id_column = "Run", core_utilisation = NA, return_filter_result = TRUE)
        data <- readWFXAdapterTable(manager, result$reference)
        expect_identical(data, expected$data)
        expect_identical(result$summary, expected$summary)
        receipt$science <- identical(data, expected$data) && identical(result$summary, expected$summary)
        receipt$output_fingerprint <- digest::digest(data, algo = "sha256")
    }
    files <- list.files(file.path(owner$manager$.context$store$directory, ".workflow-export", id),
        recursive = TRUE, full.names = TRUE)
    figures <- files[grepl("[.](png|pdf)$", files) & !grepl("/figures/", files)]
    expect_length(figures, 12L)
    receipt$required_outputs <- length(figures) == 12L && all(file.size(figures) > 0)
    expect_identical(digest::digest(file = file.path(request$source_store, "dia-state-head.rds"),
        algo = "sha256"), original)
    receipt$input_unchanged <- identical(digest::digest(file = file.path(request$source_store,
        "dia-state-head.rds"), algo = "sha256"), original)
    saveRDS(receipt, path)
})

test_that("PERF-106 repair uses actual installed action dispatch and durable measurements", {
    skip_if(!exists("workflow_test_mode") || workflow_test_mode != "adapters",
        "explicit installed application only")
    app <- workflow_test_application
    request_path <- file.path(app$config$root, "timed-member.rds")
    skip_if(!file.exists(request_path), "no frozen member")
    request <- readRDS(request_path)
    skip_if(!identical(request$programme, "perf106-repair-v2"), "different programme")
    stopifnot(identical(request$schema, 1L),
        request$operation %in% c("minimum_peptides_per_sample", "precursor_rollup"),
        request$variant %in% c("native", "physical", "selection", "p14"),
        is.character(request$fault), length(request$fault) == 1L,
        request$fault %in% c("none", "publication", "worker_failure", "cancel", "worker_loss", "monitor_loss"))
    sys.source(file.path(workflow_test_root, "tests/testthat/helper-workflow-bulk-rollup.R"),
        envir = environment())
    owned <- file.path(app$config$root, "project")
    dir.create(owned)
    original_head <- NULL
    if (isTRUE(request$preflight)) {
        recipe <- readWorkflowCheckpointRecipe("small", workflow_test_root)
        if (isTRUE(request$summary_heavy)) recipe$n_proteins <- 2000L
        object <- makeWorkflowCheckpointPeptide(recipe)
        if (request$operation == "minimum_peptides_per_sample") {
            object@peptide_data$Protein.Ids <- paste0(object@peptide_data$Protein.Ids, strrep("x", 2048L))
            if (identical(request$fault, "worker_failure")) object@peptide_data$Stripped.Sequence[[1L]] <- NA_character_
        }
        manager <- newDIAStateAdapter(newWorkflowEntityStore(owned, basename(app$config$root), "DIA"),
            WorkflowState$new(audit_enabled = FALSE))
        manager$.context$reference_states <- TRUE
        manager$saveState("protein_peptide_filtered", object, list(), "installed preflight")
    } else {
        original_head <- digest::digest(file = file.path(request$source_store, "dia-state-head.rds"),
            algo = "sha256")
        stopifnot(identical(original_head, request$source_head_sha256))
        manager <- copyDIAReplayCheckpoint(request$source_store, owned)
        manager$.context$reference_states <- TRUE
        selectDIAStateReference(manager, "protein_peptide_filtered")
        stopifnot(identical(readDIAStateMetadata(manager)$data$components$peptide_data$entity_id,
            request$input_entity_id))
    }
    withr::defer(closeDIAStateAdapter(manager))
    owner <- newWorkflowExecutionController(manager, app)
    withr::defer(closeWorkflowExecutionController(owner))
    if (identical(request$fault, "publication")) withr::local_options(
        multischolar.workflow_execution_boundary = function(stage) {
            if (stage == "publication") stop("owned PERF-106 publication failure")
        })
    before <- getDIAStateReference(manager)
    if (isTRUE(request$profile_worker)) {
        sys.source(file.path(workflow_test_root, "tools/profiling/workflow_perf106.R"),
            envir = environment())
        testthat::local_mocked_bindings(
            workflowLaunchService = perf106ProfileLauncher(workflowLaunchService),
            .package = "MultiScholaR")
    }
    id <- runWorkflowAction(request$operation, list(state = before), request$parameters, owner)
    paths <- owner$job$paths
    injected <- FALSE
    repeat {
        done <- pollWorkflowAction(owner)
        if (!done$status %in% c("admitting", "running")) break
        if (!injected && identical(done$status, "running") &&
            request$fault %in% c("cancel", "worker_loss", "monitor_loss")) {
            injected <- TRUE
            stale_token <- owner$job$request$token
            if (request$fault == "cancel") {
                cancelWorkflowAction(owner)
                done <- readRDS(file.path(paths$request_dir, "completion.rds"))
                expect_error(acceptWorkflowActionResult(owner, list(token = stale_token)),
                    class = "workflow_execution_stale")
                break
            }
            workflowSystemd(c("stop", if (request$fault == "worker_loss")
                owner$job$worker_unit else owner$job$monitor_unit))
        }
        Sys.sleep(0.02)
    }
    record <- list(status = done$status, member = request, execution = done,
        source_head_sha256 = original_head, output = NULL)
    if (isTRUE(request$summary_heavy) && done$status == "complete") {
        worker <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
        expect_true(isSavedWorkflowReference(worker$result$counts))
        record$large_summary_by_reference <- isSavedWorkflowReference(worker$result$counts)
    }
    output_path <- file.path(app$config$root, "member-result.rds")
    # Persist the measured terminal outcome before any scientific verification.
    saveRDS(record, output_path)
    failed_fixture <- request$fault %in% c("worker_failure", "worker_loss", "monitor_loss", "cancel")
    if (failed_fixture) {
        expect_identical(done$status, switch(request$fault, cancel = "cancelled",
            monitor_loss = "monitor_lost", "worker_failed"))
        expect_identical(getDIAStateReference(manager), before)
    } else if (identical(request$fault, "publication")) {
        expect_identical(done$status, "publication_pending")
        worker_file <- file.path(paths$request_dir, "worker-result.rds")
        worker_hash <- digest::digest(file = worker_file, algo = "sha256")
        accepted <- getDIAStateReference(manager)
        options(multischolar.workflow_execution_boundary = NULL)
        result <- resumeWorkflowAction(owner, id)
        expect_identical(getDIAStateReference(manager), accepted)
        expect_identical(digest::digest(file = worker_file, algo = "sha256"), worker_hash)
        expect_identical(acceptWorkflowActionResult(owner, readRDS(worker_file))$status, "duplicate")
        record$status <- "complete"
        record$output <- result
        record$publication_retry_without_science <- TRUE
    } else if (done$status == "complete") record$output <- resumeWorkflowAction(owner, id)
    if (isTRUE(request$preflight) && !failed_fixture) {
        expect_identical(record$status, "complete")
        if (!is.null(record$output)) {
            expected <- if (request$operation == "precursor_rollup")
                rollUpPrecursorToPeptide(object, core_utilisation = NA) else
                filterMinNumPeptidesPerSample(object,
                    peptides_per_sample_cutoff = request$parameters$peptides_per_sample_cutoff,
                    inclusion_list = request$parameters$inclusion_list, core_utilisation = NA)
            actual <- readDIAExecutionComponents(manager, fields = c("peptide_data", "design_matrix"))
            expect_identical(actual@peptide_data, expected@peptide_data)
            expect_identical(actual@design_matrix, expected@design_matrix)
            record$small_science_exact <- identical(actual@peptide_data, expected@peptide_data) &&
                identical(actual@design_matrix, expected@design_matrix)
        }
    }
    if (!is.null(original_head)) expect_identical(digest::digest(
        file = file.path(request$source_store, "dia-state-head.rds"), algo = "sha256"), original_head)
    expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
    expect_identical(manager$.context$store$open_resources, 0L)
    record$project_root <- owned
    saveRDS(record, output_path)
})
