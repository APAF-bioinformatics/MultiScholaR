test_that("unchanged Design changes parser metadata without reading or rewriting the assay", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "metadata-design", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    table <- prepareDIAParserTable(readr::read_tsv(I("Run\tPrecursor.Normalised\nA\t1\nB\t2\nA\t3\n"),
        show_col_types = FALSE))
    source <- saveDIATableComponent(manager$.context, table)
    recordDIAImportMetadata(manager, table, list(data = table, data_type = "peptide",
        column_mapping = list(run_col = "Run", quantity_col = "Precursor.Normalised")), source)
    before <- manager$.context$store$io_counts
    path <- file.path(manager$.context$store$directory, source$payload$relative_path)
    identity <- workflowTableFileIdentity(path)
    validate <- validateWorkflowPayloadOnce
    binding <- testthat::with_mocked_bindings(
        transformDIADesignTable(manager, source, samples = c("B", "A"), strip_parser = TRUE),
        saveWorkflowTableQuery = function(...) stop("metadata edit ran SQL production"),
        saveWorkflowTable = function(...) stop("metadata edit rewrote an assay"),
        streamSavedWorkflowTable = function(...) stop("metadata edit scanned assay rows"),
        validateWorkflowPayloadOnce = function(store, record) {
            if (identical(record$payload$codec, "parquet_duckdb")) stop("metadata edit hashed the assay")
            validate(store, record)
        })
    expect_true(isWorkflowTableComponent(binding))
    expect_identical(workflowTableFileIdentity(path), identity)
    expect_equal(manager$.context$store$io_counts$payload_hash_bytes - before$payload_hash_bytes, binding$payload$bytes)
    component <- readWorkflowTableComponent(manager$.context$store, binding)
    expect_identical(component$source, source)
    expected <- buildProtDesignInitialState(table, list(), list(quantity_col = "Precursor.Normalised"))$data_cln
    expect_identical(materializeDIATable(manager$.context, binding), expected)
    expect_identical(readSavedWorkflowTable(manager$.context$store, binding, "Run", limit = 3L)$data$Run,
        table$Run)
    expect_identical(transformDIADesignTable(manager, binding, samples = c("A", "B"), strip_parser = TRUE), binding)
    references <- expandWorkflowTableReferences(manager$.context$store,
        stats::setNames(list(binding), binding$entity_id))
    expect_setequal(names(references), c(binding$entity_id, source$entity_id))
    expect_identical(manager$.context$store$open_resources, 0L)
})

test_that("DIA adapter scopes complete records and preserves native defaults", {
    skip_if_not_installed("qs2", "0.2.2")
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "test", "DIA"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager$setWorkflowType("DIA")
    manager$saveState("raw", object, list(flag = TRUE), "source")
    expect_identical(manager$getState(), object)
    expect_null(manager$getState("missing"))
    expect_identical(manager$getHistory(), c("initial", "raw"))
    before <- diaStateOwnership(manager)
    expect_true(hasDIAState(manager, "raw"))
    expect_false(hasDIAState(manager, "missing"))
    expect_identical(diaStateOwnership(manager)$hydrates, before$hydrates)
    withDIAStateScope(manager, function() {
        expect_identical(manager$getState(), object)
        expect_identical(manager$getState(), object)
        expect_identical(getDIAStateRecord(manager, "raw")$config, list(flag = TRUE))
        expect_identical(diaStateOwnership(manager)$working_records, 1L)
    })
    after <- diaStateOwnership(manager)
    expect_identical(after$hydrates - before$hydrates, 1L)
    expect_identical(after$working_records, 0L)
    expect_identical(after$native_records, 0L)
    expect_identical(after$native_audits, 0L)
    expect_equal(after$store$bytes, 0)
    expect_error(withDIAStateScope(manager, function() {
        manager$getState()
        stop("action failed")
    }), "action failed")
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    ordinary <- WorkflowState$new()
    ordinary$saveState("raw", object, list(flag = TRUE), "source")
    path <- withr::local_tempfile()
    saveRDS(ordinary, path)
    expect_identical(readRDS(path)$getState(), object)
    expect_false(inherits(ordinary, "DIAStateAdapter"))
    expect_error(manager$setWorkflowType("LFQ"), class = "saved_entity_scope")
})

test_that("DIA local helpers and SQL rollup do not materialize unused workers", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "local", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    sections <- c("peptideIntensityFiltering", "filterMinNumPeptidesPerProtein",
        "filterMinNumPeptidesPerSample", "removePeptidesWithOnlyOneReplicate",
        "peptideMissingValueImputation", "rollUpPrecursorToPeptide")
    for (section in sections) object@args[[section]]$core_utilisation <- newDIAWorkerSpec(8L)
    manager$saveState("raw", object, list(), "raw")
    manager$.context$sql_rollup <- TRUE
    testthat::with_mocked_bindings({
        for (section in sections) {
            executable <- getDIAExecutionState(manager, section)
            expect_identical(executable@args[[section]]$core_utilisation, NA)
            expect_identical(executable@peptide_data, object@peptide_data)
        }
    }, startDIAWorkerCluster = function(...) stop("unused worker pool"))
    expect_length(manager$.context$owned_workers, 0L)
    expect_identical(manager$getState(), object)
})

test_that("DIA protein count is local with identical science and truthful execution metadata", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "count", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(), "raw")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    paths <- list(source_dir = withr::local_tempdir(), protein_qc_dir = withr::local_tempdir())
    result <- runProteinReplicateFilterApplyStep(workflow, paths, "group", 4L,
        newClusterFn = function(...) stop("unexpected pool"), isTestModeFn = function() FALSE)
    expected <- removeProteinsWithOnlyOneReplicate(object, NA, "group")
    expect_identical(result$filteredS4@protein_quant_table, expected@protein_quant_table)
    expect_identical(result$filteredS4@design_matrix, expected@design_matrix)
    expect_identical(workflow$qc_params$protein_qc$replicate_filter$parallel_cores, 4L)
    expect_identical(workflow$qc_params$protein_qc$replicate_filter$effective_parallel_cores, 1L)
    expect_match(result$resultText, "Parallel cores used: 1", fixed = TRUE)
    expect_identical(manager$revertToState("raw"), object)
})

test_that("DIA public list replacement and branch audit history remain complete", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "test", "DIA"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(x = 1L), "raw", list(record_id = "r1", value = 1L))
    manager$saveState("filtered", object, list(x = 2L), "filtered", list(record_id = "r2", value = 2L))
    first <- manager$states
    expect_type(first, "list")
    expect_identical(lapply(first, \(x) x$data), list(initial = NULL, raw = object, filtered = object))
    expect_identical(manager$revertToState("raw"), object)
    expect_identical(manager$getHistory(), c("initial", "raw"))
    expect_true(hasDIAState(manager, "filtered"))
    expect_identical(names(manager$getAuditRecords()), c("r1", "r2"))
    manager$saveState("branch", object, list(), "branch", list(record_id = "r3"))
    manager$saveState("branch", object, list(), "replacement", list(record_id = "r4"))
    expect_identical(manager$getHistory(), c("initial", "raw", "branch"))
    expect_identical(names(manager$getAuditRecords()), c("r1", "r2", "r3", "r4"))
    expect_identical(manager$getStateAudit("raw"), list(record_id = "r1", value = 1L))
    manager$states$raw$description <- "replaced through public list"
    expect_identical(getDIAStateRecord(manager, "raw")$description, "replaced through public list")
    path <- withr::local_tempfile()
    saveRDS(manager$states, path)
    expect_identical(readRDS(path), manager$states)
    manager$states <- first
    expect_identical(manager$states, first)
    manager$state_history <- as.list(names(first))
    manager$current_state <- "filtered"
    closeDIAStateAdapter(manager)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "test", "DIA"), reopen = TRUE)
    expect_identical(manager$states, first)
    expect_identical(manager$getHistory(), names(first))
    expect_identical(manager$getState(), object)
    expect_identical(names(manager$getAuditRecords()), c("r1", "r2", "r3", "r4"))
})

test_that("failed checkpoint or head publication cannot replace accepted state", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "test", "DIA"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(), "raw")
    before <- manager$states
    before_head <- manager$.context$head
    testthat::with_mocked_bindings({
        expect_error(manager$saveState("next", object, list(), "next"), "write rejected")
        expect_error(manager$states <- list(next_state = before$raw), "write rejected")
    }, writeSavedWorkflowQsPayload = function(...) stop("write rejected"))
    expect_identical(manager$.context$head, before_head)
    expect_identical(manager$states, before)
    testthat::with_mocked_bindings({
        expect_error(manager$saveState("next", object, list(), "next"), "head rejected")
        expect_error(manager$states <- list(next_state = before$raw), "head rejected")
    }, diaStateBoundary = function(phase) stop("head rejected"))
    expect_identical(manager$.context$head, before_head)
    expect_identical(manager$states, before)
    expect_identical(manager$current_state, "raw")
    expect_identical(diaStateOwnership(manager)$native_records, 0L)
    expect_error(manager$states <- list(bad = 2), class = "saved_entity_invalid")
    expect_identical(manager$states, before)
})

test_that("corruption and project replacement refuse hydration without changing topology", {
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "test", "DIA")
    manager <- newDIAStateAdapter(store)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(), "raw")
    manager$saveState("final", object, list(), "final")
    reference <- manager$.context$snapshot$states$raw
    payload <- file.path(store$directory, reference$payload$relative_path)
    writeBin(as.raw(1L), payload)
    expect_error(manager$revertToState("raw"), class = "saved_entity_corrupt")
    expect_identical(manager$current_state, "final")
    expect_identical(manager$getHistory(), c("initial", "raw", "final"))
    expect_identical(manager$getState(), object)
    identity <- store$identity
    identity$incarnation <- strrep("f", 64)
    jsonlite::write_json(identity, file.path(store$directory, "store.json"), auto_unbox = TRUE)
    expect_error(manager$getState(), class = "saved_entity_scope")
    expect_error(manager$saveState("new", object, list(), "new"), class = "saved_entity_scope")
})

test_that("portable records keep section-specific worker counts without starting workers", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    worker <- new.env(parent = emptyenv())
    worker$is_alive <- function() TRUE
    cluster <- structure(list(worker, worker), class = "multidplyr_cluster")
    object@args <- list(globalParameters = list(number_of_cpus = 8L),
        first = list(core_utilisation = cluster), second = list(core_utilisation = 1L))
    record <- list(data = object, config = object@args, timestamp = Sys.time())
    portable <- prepareDIAStateRecord(record)
    expect_identical(portable$data@args$first$core_utilisation, newDIAWorkerSpec(2L))
    expect_identical(portable$data@args$globalParameters, object@args$globalParameters)
    expect_identical(portable$data@args$second, object@args$second)
    expect_identical(record$data@args$first$core_utilisation, cluster)
    expect_identical(portable$config, portable$data@args)
})

test_that("session-owned GUI holders stay cold and resolve inside their original reactive contract", {
    root <- withr::local_tempdir()
    source <- file.path(root, "project", "scripts")
    dir.create(source, recursive = TRUE)
    withr::local_options(multischolar.dia_capacity_trial = list(
        scope = "full-journey-capacity-v1", route = "iq", treatment = "qs2_duckdb", receipt_dir = root))
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(state_manager = WorkflowState$new(audit_enabled = FALSE))
        workflow$state_manager <- initializeDIACapacitySession(workflow$state_manager,
            list(source_dir = source), session, workflow)
        manager <- workflow$state_manager
        expect_false(manager$audit_enabled)
        manager$setWorkflowType("DIA")
        saved <- newDIASavedReactiveVal(NULL, "test-design")
        shiny::observeEvent(input$apply, {
          withDIAStateScope(manager, function() {
            decodes <- manager$.context$hydrates
            manager$saveState("raw", object, list(), "raw")
            workflow <- setDIAOwnedPayload(workflow, "data_cln", object@protein_quant_table)
            saved(object)
            expect_identical(saved(), object)
            expect_identical(manager$.context$hydrates, decodes)
          })
        })
        output$rows <- shiny::renderText({
            shiny::req(saved())
            nrow(saved()@protein_quant_table)
        })
    }, {
        session$setInputs(apply = 1L)
        expect_identical(getDIAActiveManager(), manager)
        expect_true(isSavedWorkflowReference(workflow$data_cln))
        expect_identical(getDIAOwnedPayload(workflow, "data_cln"), object@protein_quant_table)
        expect_identical(saved(), object)
        expect_identical(output$rows, "64")
        session$flushReact()
        expect_identical(diaStateOwnership(manager)$native_records, 0L)
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        expect_true(file.exists(file.path(root, "gui-state-latest.json")))
        expect_true(file.exists(file.path(root, "gui-state-phases.jsonl")))
        previous <- manager
        replacement <- initializeDIACapacitySession(WorkflowState$new(),
            list(source_dir = source), session, workflow)
        expect_true(previous$.context$store$closed)
        expect_null(previous$.context$on_commit)
        expect_length(previous$.context$callbacks, 0L)
        expect_null(environment(replacement$.context$on_commit)$previous)
        expect_type(replacement$.context$on_commit, "closure")
        expect_false(identical(previous$.context$store$identity, replacement$.context$store$identity))
        expect_true(dir.exists(previous$.context$store$directory))
        expect_error(previous$getState(), class = "saved_entity_scope")
        expect_identical(replacement$getHistory(), "initial")
        expect_identical(getDIAActiveManager(), replacement)
    })
})

test_that("candidate three-step production Apply helpers preserve native science and audit", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    # Same scientific input contract as makeRecordedPilotTestSession, not a
    # column-less default or mocked scientific transformation.
    object@args$srlQvalueProteotypicPeptideClean <- list(
        input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
            "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
            "Precursor.Quantity", "Precursor.Normalised"),
        exclude_contaminants = FALSE, qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
        choose_only_proteotypic_peptide = 1, confidence_provenance_mode = "diann_main_report")
    object@args$filterMinNumPeptidesPerProtein <- list(
        peptides_per_protein_cutoff = 1, peptidoforms_per_protein_cutoff = 2)
    object@args$rollUpPrecursorToPeptide <- list(core_utilisation = NA)
    candidate <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "steps", "DIA"))
    on.exit(closeDIAStateAdapter(candidate), add = TRUE)
    candidate$.context$sql_rollup <- TRUE
    native <- WorkflowState$new()
    now <- function() as.POSIXct("2026-09-21 01:00:00", tz = "UTC")
    managers <- list(native, candidate)
    outputs <- lapply(managers, function(manager) {
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        config_env <- new.env(parent = emptyenv())
        config_env$config_list <- object@args
        update <- function(...) updateConfigParameter(..., env = config_env)
        manager$saveState("raw_data_s4", object, list(), "raw")
        withDIAStateScope(manager, function() runPeptideQvalueApplyStep(workflow,
            0.01, 0.01, TRUE, updateConfigParameterFn = update,
            globalPGQvalueThreshold = 0.01, nowFn = now))
        withDIAStateScope(manager, function() runProteinPeptideApplyStep(workflow,
            1, 2, updateConfigParameterFn = update, nowFn = now))
        withDIAStateScope(manager, function() runPeptideRollupApplyStep(workflow, nowFn = now))
        manager$getState()
    })
    expect_identical(outputs[[2L]], outputs[[1L]])
    expect_identical(candidate$getHistory(), native$getHistory())
    expect_identical(candidate$getAuditRecords(), native$getAuditRecords())
    expect_identical(diaStateOwnership(candidate)$native_records, 0L)
    expect_identical(diaStateOwnership(candidate)$working_records, 0L)
    expect_identical(candidate$.context$store$open_resources, 0L)
})

test_that("external session adoption is atomic and native export remains self-contained", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "restore", "DIA"))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("correlation_filtered", object, list(), "final")
    before <- manager$.context$head
    records <- manager$states
    testthat::with_mocked_bindings({
        expect_error(restoreDIAStateSession(manager, records, as.list(names(records)),
            "correlation_filtered"), "restore interrupted")
    }, diaStateBoundary = function(...) stop("restore interrupted"))
    expect_identical(manager$.context$head, before)
    restoreDIAStateSession(manager, records, as.list(names(records)), "correlation_filtered")
    output <- withr::local_tempdir()
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    input <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    result <- runProtNormExportSessionWorkflow(workflow, norm, input, output,
        withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) NULL, messageFn = function(...) NULL)
    expect_null(result$sessionData)
    actual <- readRDS(result$exportArtifacts$latestFilepath)
    expect_identical(actual$r6_complete_states, records)
    expect_identical(actual$current_s4_object, object)
    expect_identical(actual$correlation_filtered_s4, object)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$scope_depth, 0L)
    restored <- WorkflowState$new()
    restored$states <- actual$r6_complete_states
    restored$state_history <- actual$r6_state_history
    restored$current_state <- actual$r6_current_state_name
    expect_identical(restored$getState(), object)
})

test_that("sampler distinguishes a vanished terminal cgroup from a live missing observation", {
    shell <- readLines(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.sh"))
    prefix <- shell[seq_len(which(grepl("^proof_script_dir=", shell)) - 1L)]
    counters <- withr::local_tempfile()
    writeLines(c("anon 50", "file 20", "kernel 5", "file_dirty 0", "file_writeback 0"), counters)
    run <- function(command) processx::run("bash", c("-c", paste(c(prefix, command), collapse = "\n")),
        error_on_status = FALSE)
    sample <- run(sprintf("checkpointMemorySample %s 1 75 80", shQuote(counters)))
    expect_identical(sample$status, 0L)
    expect_identical(trimws(sample$stdout), "1\t75\t80\t50\t20\t5\t0\t0")
    expect_true(run(sprintf("checkpointMemorySample %s 1 75 80", shQuote(paste0(counters, "-missing"))))$status != 0L)
    expect_identical(run("checkpointServiceTerminal 'SubState=exited'")$status, 0L)
    expect_identical(run("checkpointServiceTerminal 'SubState=failed'")$status, 0L)
    expect_true(run("checkpointServiceTerminal 'SubState=running'")$status != 0L)
    expect_true(run("checkpointServiceTerminal 'SubState=stop-sigterm'")$status != 0L)
})

test_that("all peptide reference adapters preserve native science and configuration", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$srlQvalueProteotypicPeptideClean <- list(input_matrix_column_ids = names(object@peptide_data))
    object@args$rollUpPrecursorToPeptide <- list(core_utilisation = NA)
    # A genuinely absent run/feature observation exercises generated-row ordering.
    table <- object@peptide_data
    removed <- table$Run == unique(table$Run)[[2L]] & table$Protein.Group == table$Protein.Group[[1L]] &
        table$Stripped.Sequence == table$Stripped.Sequence[[1L]]
    object@peptide_data <- table[!removed, , drop = FALSE]
    rlang::local_bindings(config_list = object@args, .env = .GlobalEnv)
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("raw_data_s4", object, object@args, "raw")
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "reference-journey", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    table_reference <- saveDIATableComponent(manager$.context, object@peptide_data)
    checkpoint <- newDIAPeptideCheckpoint(manager, table_reference, object@design_matrix, object@args)
    commitDIAReferenceState(manager, "raw_data_s4", checkpoint, object@args, "raw")
    workflow <- list2env(list(state_manager = manager), parent = emptyenv())
    expected_workflow <- list2env(list(state_manager = native), parent = emptyenv())
    actions <- list(
        qvalue = function(workflow) runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE),
        evidence = function(workflow) runProteinPeptideApplyStep(workflow, 1, 2),
        rollup = function(workflow) runPeptideRollupApplyStep(workflow),
        intensity = function(workflow) runPeptideIntensityApplyStep(workflow, FALSE, 1, 1, 0),
        sample = function(workflow) runPeptideSampleApplyStep(workflow, 1),
        replicate = function(workflow) runPeptideReplicateApplyStep(workflow, "tech_rep_group"),
        imputation = function(workflow) runPeptideImputationStep(workflow, 1)
    )
    summary_reader <- summarizeSavedDIAFiltering
    for (name in names(actions)) {
        expected <- actions[[name]](expected_workflow)
        result <- testthat::with_mocked_bindings(actions[[name]](workflow),
            materializeDIAStateRecord = function(...) stop("GUI science hydrated a complete state"),
            summarizeSavedDIAFiltering = function(scope, prototype, protein_column, counted = NULL) {
                testthat::with_mocked_bindings(summary_reader(scope, prototype, protein_column, counted),
                    streamSavedWorkflowTable = function(...) stop("post-production summary scanned assay"))
            })
        expect_true(isSavedWorkflowReference(result$stateReference), info = name)
        expect_false(any(vapply(result, isS4, logical(1))), info = name)
        expect_equal(manager$getState(), native$getState(), tolerance = 1e-12, info = name)
        expect_equal(readDIAStateMetadata(manager)$config, native$states[[native$current_state]]$config, info = name)
        expect_identical(result$resultText, expected$resultText, info = name)
        expect_identical(manager$.context$store$open_resources, 0L, info = name)
        # This fixture isolates science. The real output/retry fixtures exercise publication.
        acknowledgeDIAOutputs(manager, result$stateReference, "qc")
    }
    before <- getDIAStateReference(manager)
    expect_true(isSavedWorkflowReference(runPeptideImputationRevertStep(workflow)$revertedS4))
    expect_false(identical(getDIAStateReference(manager), before))
})

test_that("wide components retain types and metadata-only checkpoints reuse assay files", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "wide-reference", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    table <- data.frame(id = c("p1", "p2", "p3"), check.names = FALSE)
    for (index in seq_len(130L)) table[[paste0("sample_", index)]] <- c(NA_real_, NaN, Inf)
    reference <- saveDIATableComponent(manager$.context, table)
    expect_identical(materializeDIATable(manager$.context, reference), table)
    parsed <- prepareDIAParserTable(readr::read_tsv(I("id\tvalue\nA\t1\nB\t2"), show_col_types = FALSE))
    parsed$id <- stats::setNames(parsed$id, c("first", "second"))
    attrs <- attributes(parsed)
    attributes(parsed) <- attrs[c("spec", "problems", "names", "row.names", "class")]
    saved <- saveDIATableComponent(manager$.context, parsed)
    restored <- materializeDIATable(manager$.context, saved)
    expect_identical(attributes(restored), attributes(parsed))
    expect_identical(attributes(restored$id), attributes(parsed$id))
    manager$saveState("protein_replicate_filtered", object, list(), "input")
    before <- readDIAStateMetadata(manager)$data$components
    view <- readDIAExecutionComponents(manager)
    view@args$gui_note <- list(value = "metadata only")
    testthat::with_mocked_bindings(manager$saveState("annotated", view, list(), "metadata"),
        writeSavedWorkflowParquet = function(...) stop("unchanged assay was rewritten"))
    expect_identical(readDIAStateMetadata(manager)$data$components, before)
})

test_that("reference import design Q-value and Revert preserve the full scientific state", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@peptide_data$Decoy <- FALSE
    object@peptide_data$Contaminant <- FALSE
    object@peptide_data$Decoy[5L] <- TRUE
    object@peptide_data$Protein.Ids[6L] <- "REV__P00001"
    object@peptide_data$Protein.Group[7:8] <- c("", NA_character_)
    object@peptide_data$Precursor.Quantity[1:4] <- c(NA_real_, NaN, Inf, 0)
    object@peptide_data$Q.Value[2:4] <- c(NA_real_, Inf, NaN)
    names(object@peptide_data$Run) <- paste0("import_", object@peptide_data$Run)
    object@design_matrix$group <- factor(object@design_matrix$group,
        levels = c("B", "A", "unused"), ordered = TRUE)
    object@design_matrix <- rbind(object@design_matrix, object@design_matrix[1L, ])
    object@design_matrix$Run[nrow(object@design_matrix)] <- "absent_design_observation"
    attr(object@peptide_data, "scientific_provenance") <- "reference fixture"
    object@args$srlQvalueProteotypicPeptideClean <- list(
        input_matrix_column_ids = names(object@peptide_data),
        contaminant_manifest = c("P00002", "P00003"))
    rlang::local_bindings(config_list = object@args, .env = .GlobalEnv)
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("raw_data_s4", object, object@args, "raw")
    original <- new.env(parent = emptyenv())
    original$state_manager <- native
    clock <- function() as.POSIXct("2026-09-24", tz = "UTC")
    expected <- runPeptideQvalueApplyStep(original, 0.01, 0.01, TRUE, nowFn = clock)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "reference", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    source <- saveDIATableComponent(manager$.context, object@peptide_data)
    expect_identical(materializeDIATable(manager$.context, source), object@peptide_data)
    empty_table <- object@peptide_data[0, , drop = FALSE]
    empty_reference <- saveDIATableComponent(manager$.context, empty_table)
    expect_identical(materializeDIATable(manager$.context, empty_reference), empty_table)
    streamed <- withr::local_tempfile()
    original_file <- withr::local_tempfile()
    writeDIADesignTable(manager, source, streamed)
    utils::write.table(object@peptide_data, original_file, sep = "\t", quote = FALSE, row.names = FALSE)
    expect_identical(digest::digest(file = streamed, algo = "sha256"),
        digest::digest(file = original_file, algo = "sha256"))
    writeLines("previous completed table", streamed)
    expect_error(testthat::with_mocked_bindings(writeDIADesignTable(manager, source, streamed),
        streamSavedWorkflowTable = function(scope, columns, consume) {
            consume(object@peptide_data[1:4, ], 0L, as.character(1:4))
            stop("table stream interrupted")
        }), "table stream interrupted")
    expect_identical(readLines(streamed), "previous completed table")
    writes <- manager$.context$writes
    expect_identical(transformDIADesignTable(manager, source,
        samples = unique(object@peptide_data$Run)), source)
    expect_identical(manager$.context$writes, writes)
    selected <- transformDIADesignTable(manager, source,
        samples = unique(object@peptide_data$Run)[1:2])
    writes <- manager$.context$writes
    expect_identical(transformDIADesignTable(manager, source,
        samples = unique(object@peptide_data$Run)[1:2]), selected)
    expect_identical(manager$.context$writes, writes)
    checkpoint <- newDIAPeptideCheckpoint(manager, source, object@design_matrix, object@args)
    raw <- commitDIAReferenceState(manager, "raw_data_s4", checkpoint, object@args, "raw")
    expect_identical(manager$getState(), object)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    reader <- materializeDIATable
    result <- testthat::with_mocked_bindings(
        runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE, nowFn = clock),
        materializeDIAStateRecord = function(...) stop("full S4 access is forbidden"),
        materializeDIATable = function(context, reference) {
            if (!identical(reference, checkpoint$components$design_matrix)) stop("assay hydration is forbidden")
            reader(context, reference)
        })
    expect_null(result$filteredS4)
    expect_identical(result$resultText, expected$resultText)
    expect_identical(manager$getState(), expected$filteredS4)
    expect_identical(workflow$qc_params, original$qc_params)
    expected_table <- expected$filteredS4@peptide_data
    expect_identical(result$filteringSummary$proteins, countUniqueProteins(expected_table))
    expect_equal(result$filteringSummary$proteins_per_run,
        as.data.frame(countProteinsPerRun(expected_table)))
    plot_root <- withr::local_tempdir()
    graph_root <- withr::local_tempdir()
    manager$.context$experiment_paths <- list(time_dir = plot_root)
    arguments <- list(data = result$filteringSummary, step_name = "2_qval_filtered",
        omic_type = "proteomics", experiment_label = "fixture", return_grid = TRUE, overwrite = TRUE,
        project_dirs = list(proteomics_fixture = list(time_dir = plot_root,
            publication_graphs_dir = graph_root)))
    completed_file <- file.path(plot_root, "2_qval_filtered_proteins_total.png")
    writeLines("previous completed figure", completed_file)
    head_before_output <- manager$.context$head
    expect_error(testthat::with_mocked_bindings(
        publishDIAQvalueOutput(manager, result, arguments, updateProteinFiltering),
        workflowExportBoundary = function(phase) {
            if (phase == "after_qvalue_figure_publication") stop("publication interrupted")
        }), "publication interrupted")
    expect_identical(manager$.context$head, head_before_output)
    expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
    request <- file.path(plot_root, ".workflow-export",
        paste0("qvalue-", substr(result$stateReference$entity_id, 1L, 40L)))
    expect_identical(readLines(file.path(request, "previous", basename(completed_file))),
        "previous completed figure")
    store_root <- dirname(manager$.context$store$directory)
    closeDIAStateAdapter(manager)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root, "reference", "DIA"), reopen = TRUE)
    manager$.context$experiment_paths <- list(time_dir = plot_root)
    workflow$state_manager <- manager
    recovered <- testthat::with_mocked_bindings({
        resumed <- runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE, nowFn = clock)
        publishDIAQvalueOutput(manager, resumed, arguments, function(...) stop("render repeated"))
    }, readSavedWorkflowTable = function(...) stop("assay read during publication retry"),
        commitDIAReferenceState = function(...) stop("science repeated"))
    expect_s3_class(recovered, "gtable")
    expect_identical(getDIAStateReference(manager), result$stateReference)
    expect_true(manager$.context$snapshot$completed_outputs[[result$stateReference$entity_id]]$qc)
    expect_silent(assertDIAExportReady(manager))
    head <- manager$.context$head
    testthat::with_mocked_bindings({
        expect_identical(selectDIAStateReference(manager, "raw_data_s4"), raw)
        expect_identical(readDIAStateMetadata(manager)$data$components$peptide_data, source)
    }, materializeDIAStateRecord = function(...) stop("Revert hydrated a state"))
    expect_false(identical(head, manager$.context$head))
    page <- readSavedWorkflowTable(manager$.context$store, source,
        c("Run", "Precursor.Quantity"), limit = 4L)$data
    expected_page <- object@peptide_data[1:4, c("Run", "Precursor.Quantity")]
    rownames(expected_page) <- as.character(1:4)
    attr(expected_page, "scientific_provenance") <- attr(object@peptide_data, "scientific_provenance")
    expect_identical(page, expected_page)
    expect_identical(diaStateOwnership(manager)$store$open_resources, 0L)
})

test_that("component metadata commits do not consume historical assays", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "commit-io", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", object, list(), "science")
    record <- readDIAStateMetadata(manager)
    checkpoint <- record$data
    checkpoint$slots$args$gui_note <- "metadata only"
    original <- validateWorkflowPayloadOnce
    testthat::with_mocked_bindings({
        commitDIAReferenceState(manager, "annotated", checkpoint, list(), "metadata")
        selectDIAStateReference(manager, "protein_replicate_filtered")
        expect_identical(readDIAStateMetadata(manager)$data$components, record$data$components)
    }, validateWorkflowPayloadOnce = function(store, reference) {
        if (reference$payload$codec == "parquet_duckdb") stop("metadata consumed an assay")
        original(store, reference)
    })
    expect_identical(manager$getState(), object)
})

test_that("streamed progress preserves exact native populations and caches without reads", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data
    table$Run[1L] <- NA_character_
    table$Protein.Ids[2:3] <- c(NA_character_, "")
    table$Stripped.Sequence[4:5] <- c(NA_character_, "")
    table$Protein.Group[6L] <- NA_character_
    table$Run <- factor(table$Run, levels = c(rev(unique(stats::na.omit(table$Run))), "absent"))
    table$Precursor.Quantity[1:4] <- c(NA_real_, NaN, Inf, 0)
    for (rows in list(seq_len(nrow(table)), integer())) {
        data <- table[rows, , drop = FALSE]
        accumulator <- newDIAFilteringAccumulator(data[0, , drop = FALSE], "Protein.Group")
        for (start in seq.int(1L, max(1L, nrow(data)), by = 127L)) {
            at <- if (nrow(data)) seq.int(start, min(start + 126L, nrow(data))) else integer()
            accumulator$push(data[at, , drop = FALSE])
        }
        counted <- accumulator$finish()
        object@peptide_data <- data
        expect_identical(counted$summary$proteins, countUniqueProteins(object))
        expect_identical(counted$summary$total_peptides, calcTotalPeptides(object))
        expect_identical(counted$proteins, length(unique(data$Protein.Group)))
        expect_equal(as.character(counted$runs), as.character(unique(data$Run)))
        for (field in c("peptides_per_protein", "proteins_per_run", "peptides_per_run")) {
            oracle <- switch(field, peptides_per_protein = calcPeptidesPerProtein(object),
                proteins_per_run = countProteinsPerRun(object), peptides_per_run = countPeptidesPerRun(object))
            oracle <- as.data.frame(oracle)
            key <- names(oracle)[[1L]]
            oracle[[key]] <- as.character(oracle[[key]])
            oracle <- oracle[order(oracle[[key]], method = "radix", na.last = TRUE), , drop = FALSE]
            rownames(oracle) <- NULL
            expect_equal(counted$summary[[field]], oracle)
        }
        store <- newWorkflowEntityStore(withr::local_tempdir(), "summary", "DIA")
        ref <- saveWorkflowTable(store, data, "source",
            newSavedWorkflowScientificIdentity(digest::digest(data, algo = "sha256"),
                "observed", "summary-test", "all-rows", "test"), ".rownames", metadata_version = 3L)
        scope <- newWorkflowTableReadScope(store, ref, reader = "arrow")
        saved <- testthat::with_mocked_bindings(
            summarizeSavedDIAFiltering(scope, NULL, "Protein.Group", counted),
            streamSavedWorkflowTable = function(...) stop("producer summary reread assay"))
        hit <- testthat::with_mocked_bindings(summarizeSavedDIAFiltering(scope, NULL, "Protein.Group"),
            streamSavedWorkflowTable = function(...) stop("summary hit reread assay"))
        expect_identical(hit, saved)
        releaseWorkflowTableReadScope(scope)
        closeWorkflowEntityStore(store)
    }
})

test_that("protein progress is cached by declared population and excludes only missing quantities", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    samples <- as.character(object@design_matrix[[object@sample_id]])
    object@protein_quant_table[[samples[[1L]]]][1:4] <- c(NA_real_, NaN, Inf, 0)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "protein-counts", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein", object, object@args, "protein")
    reference <- getDIAStateReference(manager)
    value <- readDIAProteinProgress(manager, reference)
    expect_equal(value$proteins_per_run, countProteinsPerRun(object))
    expect_equal(value$proteins, nrow(object@protein_quant_table))
    hit <- testthat::with_mocked_bindings(
        testthat::with_mocked_bindings(readDIAProteinProgress(manager, reference),
            dbGetQuery = function(...) stop("count hit queried assay"), .package = "DBI"),
        readDIAExecutionComponents = function(...) stop("count hit read design"))
    expect_identical(hit, value)
})
