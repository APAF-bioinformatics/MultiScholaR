#' Build the existing small DIA input with its actual confidence/evidence contract.
#' @return Peptide S4 fixture for the two unchanged scientific Apply helpers.
#' @noRd
makeDIALifetimePeptide <- function() {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$srlQvalueProteotypicPeptideClean <- list(
        input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
            "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
            "Precursor.Quantity", "Precursor.Normalised"),
        exclude_contaminants = FALSE, qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
        choose_only_proteotypic_peptide = 1, confidence_provenance_mode = "diann_main_report")
    object@args$filterMinNumPeptidesPerProtein <- list(
        peptides_per_protein_cutoff = 1, peptidoforms_per_protein_cutoff = 2)
    object
}

test_that("owned operations release replaced records before their enclosing GUI flush", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "owner-action", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    owner <- new.env(parent = emptyenv())
    callbacks <- new.env(parent = emptyenv())
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- function(fn, once) { callbacks$before <- fn; function() invisible(NULL) }
    session$onFlushed <- function(fn, once) { callbacks$after <- fn; function() invisible(NULL) }
    session$onSessionEnded <- function(fn) { callbacks$end <- fn; function() invisible(NULL) }
    attachDIAStateSession(manager, session, owner)
    testthat::local_mocked_bindings(getDIAActiveManager = function() manager)
    first <- data.frame(Run = c("a", "b"), value = c(1, NA_real_))
    second <- first
    second$value <- c(2, Inf)
    callbacks$before()
    expect_identical(manager$.context$scope_depth, 0L)
    withDIAStateScope(manager, function() {
        setDIAOwnedPayload(owner, "data_tbl", first)
        previous <- owner$data_tbl
        alias <- previous
        setDIAOwnedPayload(owner, "data_tbl", second)
        expect_false(exists(previous$entity_id, manager$.context$working, inherits = FALSE))
        expect_identical(diaStateOwnership(manager)$working_records, 1L)
        reads <- manager$.context$hydrates
        withDIAStateScope(manager, function() {
            expect_identical(getDIAOwnedPayload(owner, "data_tbl"), second)
        })
        expect_identical(manager$.context$hydrates, reads)
        expect_identical(resolveDIAOwnedValue(alias), first)
        expect_true(file.exists(file.path(manager$.context$store$directory, alias$payload$relative_path)))
    })
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$scope_depth, 0L)
    callbacks$after()
    previous <- owner$data_tbl
    expect_error(withDIAStateScope(manager, function() {
        testthat::with_mocked_bindings(setDIAOwnedPayload(owner, "data_tbl", first),
            saveDIAStateValue = function(...) stop("injected owner save failure"))
    }), "injected owner save failure")
    expect_identical(owner$data_tbl, previous)
    expect_identical(getDIAOwnedPayload(owner, "data_tbl"), second)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$scope_depth, 0L)
})

test_that("completed publication release respects real owners and adapter replacement", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "reclaim", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    manager$saveState("raw", data.frame(value = c(1, NA_real_)), list(), "fixture")
    context <- manager$.context
    reference <- context$snapshot$states$raw
    head <- readBin(context$head_path, "raw", n = file.size(context$head_path))
    callbacks <- new.env(parent = emptyenv())
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- function(fn, once) { callbacks$before <- fn; function() NULL }
    session$onFlushed <- function(fn, once) { callbacks$after <- fn; function() NULL }
    session$onSessionEnded <- function(fn) { callbacks$end <- fn; function() NULL }
    attachDIAStateSession(manager, session, new.env(parent = emptyenv()))
    calls <- 0L
    testthat::local_mocked_bindings(reclaimWorkflowMemory = function() {
        calls <<- calls + 1L
        list(status = "released")
    })
    expect_false(requestDIAReclamation(manager, reference, "one", "qc"))
    context$reclaim_enabled <- TRUE
    expect_false(requestDIAReclamation(manager, reference, "one", "unmeasured"))
    expect_true(requestDIAReclamation(manager, reference, "one", "qc"))
    withDIAStateScope(manager, function() {
        callbacks$after()
        expect_identical(context$scope_depth, 1L)
        expect_identical(calls, 0L)
        expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    })
    context$store$open_resources <- 1L
    expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    context$store$open_resources <- 0L
    context$working$held <- 1L
    expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    rm("held", envir = context$working)
    context$borrowed_workers <- list("active")
    expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    context$borrowed_workers <- list()
    context$owned_workers <- list("active")
    expect_identical(reclaimDIACompletedMemory(manager)$status, "busy")
    context$owned_workers <- list()
    callbacks$after()
    callbacks$after()
    expect_identical(calls, 1L)
    expect_false(requestDIAReclamation(manager, reference, "one", "qc"))
    expect_error(withDIAStateScope(manager, function() {
        requestDIAReclamation(manager, reference, "failed", "normalization")
        stop("cancelled action")
    }), "cancelled action")
    expect_null(context$reclaim_pending)
    requestDIAReclamation(manager, reference, "superseded", "qc")
    context$reclaim_pending$incarnation <- "old"
    expect_identical(reclaimDIACompletedMemory(manager)$status, "superseded")
    requestDIAReclamation(manager, reference, "replacement", "qc")
    session$userData$dia_capacity_manager <- NULL
    callbacks$after()
    expect_identical(calls, 1L)
    session$userData$dia_capacity_manager <- manager
    testthat::with_mocked_bindings({
        expect_identical(reclaimDIACompletedMemory(manager)$status, "failed")
        expect_identical(reclaimDIACompletedMemory(manager)$status, "not_requested")
    }, reclaimWorkflowMemory = function() cli::cli_abort("injected", class = "workflow_reclaim_error"))
    expect_identical(readBin(context$head_path, "raw", n = file.size(context$head_path)), head)
    expect_identical(manager$getState(), data.frame(value = c(1, NA_real_)))
    requestDIAReclamation(manager, reference, "closing", "qc")
    closeDIAStateAdapter(manager)
    callbacks$after()
    expect_identical(calls, 1L)
    expect_null(context$reclaim_pending)
})

test_that("active-owner diagnostics record cleanup and failure without payload traversal", {
    root <- withr::local_tempdir()
    source <- file.path(root, "project", "scripts")
    dir.create(source, recursive = TRUE)
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(state_manager = WorkflowState$new(audit_enabled = FALSE))
        # KIT-011: the on-disk owner attached as the import does it, with the harness receipt directory.
        state_root <- tempfile(".dia-native-", tmpdir = source)
        dir.create(state_root, mode = "0700")
        manager <- newDIAStateAdapter(newWorkflowEntityStore(state_root, "DIA-IQ-native", "DIA", max_cache_bytes = 0),
            workflow$state_manager)
        # The retired capacity trial ran the adapter in the legacy session format without reference states;
        # these owner diagnostics are exercised in that context (the native-format path has its own tests).
        manager$.context$sql_rollup <- TRUE
        manager$.context$session_format <- "legacy"
        manager$.context$reference_states <- FALSE
        manager$.context$experiment_paths <- list(source_dir = source)
        manager <- attachDIAStateSession(manager, session, workflow, root)
        workflow$state_manager <- manager
        saved <- newDIASavedReactiveVal(NULL, "late-owner")
        shiny::observeEvent(input$apply, {
            withDIAStateScope(manager, function() {
                saved(data.frame(value = 1:3))
                saved(data.frame(value = 4:6))
                stopifnot(diaStateOwnership(manager)$working_records == 1L)
            })
            stopifnot(diaStateOwnership(manager)$working_records == 0L)
        })
    }, {
        session$setInputs(apply = 1L)
        expect_error(withDIAStateScope(manager, function() {
            getDIAOwnedPayload(list(value = attr(saved, "dia_reference")()), "value")
            stop("injected action failure")
        }), "injected action failure")
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        records <- lapply(readLines(file.path(root, "gui-owner-phases.jsonl")), jsonlite::fromJSON)
        expect_true(any(vapply(records, function(x) identical(x$stage, "operation_released") &&
            x$working_records == 0, logical(1))))
        expect_true(any(vapply(records, function(x) identical(x$stage, "consumed_reference_released"), logical(1))))
        expect_true(all(vapply(records, function(x) !any(c("data", "state_sizes", "payload", "object") %in% names(x)), logical(1))))
        expect_lt(file.size(file.path(root, "gui-owner-phases.jsonl")), 65536)
    })
})

test_that("DIA import persists one sanitized table and releases it before FASTA", {
    root <- withr::local_tempdir()
    original <- data.frame(Run = c("A B", "a-b", "A B"), Protein.Group = c("p1", NA, "p1"),
        Stripped.Sequence = factor(c("pep1", "pep2", "pep1"), levels = c("unused", "pep1", "pep2")),
        Precursor.Normalised = c(Inf, NA_real_, 0),
        Batch = factor(c("b2", "b1", "b2"), levels = c("b3", "b1", "b2")))
    attr(original, "scientific_annotation") <- "retained"
    imported <- list(data = original, data_type = "peptide", column_mapping = list(
        protein_col = "Protein.Group", peptide_col = "Stripped.Sequence", run_col = "Run",
        quantity_col = "Precursor.Normalised", qvalue_col = "Q.Value", pg_qvalue_col = "PG.Q.Value"))
    expected <- sanitizeProtImportTable(original, "Run")$data
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(state_manager = newDIAStateAdapter(
            newWorkflowEntityStore(root, "import-life", "DIA"), WorkflowState$new(audit_enabled = FALSE)),
            tab_status = list(), processing_log = list())
        manager <- workflow$state_manager
        attachDIAStateSession(manager, session, workflow)
        local_data <- createProtImportLocalData()
        observed <- new.env(parent = emptyenv())
        observed$fasta_calls <- 0L
        registerProtImportObservers(input, output, session, workflow, local_data,
            list(source_dir = root), FALSE, function() "fixture.tsv", function() "fixture.fasta",
            function() "diann", startProcessingFn = function(...) NULL,
            runFormatDetectionFn = function(...) NULL,
            readDataWithStatusFn = function(...) imported,
            loadConfigurationResourcesFn = function(...) NULL,
            processFastaDataFn = function(...) {
                frame <- parent.frame()
                expect_null(get("data_import_result", frame)$data)
                expect_null(get("apply_arguments", frame))
                expect_identical(diaStateOwnership(manager)$working_records, 0L)
                observed$fasta_calls <- observed$fasta_calls + 1L
            })
    }, {
        session$setInputs(sanitize_names = TRUE, organism_name = "fixture", taxon_id = 9606,
            process_data = 1L)
        expect_identical(observed$fasta_calls, 1L)
        expect_identical(getDIAOwnedPayload(workflow, "data_tbl"), expected)
        expect_identical(workflow$data_cln, workflow$data_tbl)
        descriptors <- as.list(manager$.context$store$records, all.names = TRUE)
        expect_equal(sum(vapply(descriptors, function(d)
            identical(d$identity$parameters, "owner:data_tbl"), logical(1))), 1)
        metadata <- getDIAImportMetadata(manager)
        expect_identical(metadata$import_summary, describeProtImportData(imported))
        expect_identical(metadata$reference_id, workflow$data_tbl$entity_id)
        expect_identical(metadata$protein_ids, unique(original$Protein.Group))
        expect_equal(nrow(metadata$samples), 2)
        expect_length(metadata$needs_numeric, 0L)
        expect_identical(workflow$processing_log$setup_import$n_rows, 3L)
        expect_identical(manager$getHistory(), "initial")
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        expect_lt(as.numeric(object.size(metadata)), 65536)
    })
})

test_that("DIA design metadata stays cold through late activation, edits and reset", {
    rlang::local_bindings(config_list = NULL, contrasts_tbl = NULL, .env = .GlobalEnv)
    root <- withr::local_tempdir()
    table <- data.frame(Run = c("s2", "s1", "s2"), Quantity = c(NA_real_, Inf, 0),
        Batch = factor(c("b2", "b1", "b2"), levels = c("unused", "b1", "b2")))
    attr(table, "annotation") <- "keep"
    mapping <- list(run_col = "Run", quantity_col = "Quantity")
    config <- list(deAnalysisParameters = list(formula_string = "~0+group"))
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(data_tbl = NULL, config_list = config,
            column_mapping = mapping, state_manager = WorkflowState$new(audit_enabled = FALSE))
        # The real parent provider and saved reactive shells exist before activation.
        builder <- registerProtDesignBuilderModule(workflow, builderServerExists = TRUE,
            builderServerFn = function(id, data_tbl, config_list, column_mapping) {
                initializeProtDesignBuilderServerState(input, session, data_tbl, config_list, column_mapping)
            })
        registerProtDesignAssignMetadataObserver(input, builder$designMatrix, builder$groups)
        registerProtDesignAssignTechRepsObserver(input, builder$designMatrix)
        registerProtDesignRenameSampleObserver(input, builder$designMatrix, builder$dataClnReactive, session)
        registerProtDesignBulkRenameObserver(input, builder$designMatrix, builder$dataClnReactive)
        registerProtDesignSaveResultsObserver(input, builder$designMatrix, builder$removedSamples,
            builder$dataClnReactive, builder$contrasts, function() config, builder$resultRv)
    }, {
        session$flushReact()
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "design-life", "DIA"))
        attachDIAStateSession(manager, session, workflow)
        workflow$state_manager <- manager
        setDIAOwnedPayload(workflow, "data_tbl", table)
        recordDIAImportMetadata(manager, table,
            list(data = table, data_type = "peptide", column_mapping = mapping), workflow$data_tbl)
        reads <- manager$.context$hydrates
        session$flushReact()
        initial <- builder$initialState()
        expect_null(initial$data_cln)
        expect_identical(initial$data_reference, workflow$data_tbl)
        expect_identical(getDIAReactiveReference(builder$dataClnReactive), workflow$data_tbl)
        expect_identical(builder$designMatrix()$batch, c("b1", "b2"))
        expect_identical(manager$.context$hydrates, reads)
        displays <- new.env(parent = emptyenv())
        deferred <- function(expr) {
            expr <- substitute(expr)
            owner <- parent.frame()
            function() eval(expr, owner)
        }
        workflow$processing_log <- list(setup_import = c(
            getDIAImportMetadata(manager)$import_summary,
            list(detected_format = "diann", data_type = "peptide", organism = "human", taxon_id = 9606)))
        registerProtImportOutputs(displays, workflow, list(processing = FALSE), session,
            function() "diann", renderUiFn = deferred, renderDtFn = deferred)
        expect_match(as.character(displays$import_status()), "Data imported successfully", fixed = TRUE)
        expect_match(as.character(displays$data_summary()), "Data Summary", fixed = TRUE)
        expect_identical(manager$.context$hydrates, reads)
        payload_count <- function() sum(vapply(as.list(manager$.context$store$records),
            function(x) identical(x$payload$codec, "qs2"), logical(1)))
        before_count <- payload_count()
        session$setInputs(selected_runs = c("s1", "s2"), factor1_select = "A",
            factor2_select = "", replicate_start = 1L, assign_metadata = 1L)
        session$setInputs(tech_rep_samples = c("s1", "s2"), tech_rep_assignment_mode = "lowest",
            assign_tech_reps = 1L)
        expect_identical(builder$designMatrix()$group, c("A", "A"))
        expect_identical(builder$designMatrix()$replicates, c(1L, 1L))
        expect_identical(manager$.context$hydrates, reads)
        expect_equal(payload_count(), before_count)

        session$setInputs(sample_to_rename = "s2", new_sample_name = "renamed", rename_sample = 1L)
        expect_identical(manager$.context$hydrates, reads + 1L)
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        renamed <- applyProtDesignSingleRenameUpdate(builder$designMatrix(), table, "s2", "renamed")$dataCln
        expect_identical(builder$dataClnReactive(), renamed)
        expect_identical(manager$.context$snapshot$design_current,
            getDIAReactiveReference(builder$dataClnReactive))
        reads <- manager$.context$hydrates
        applyProtDesignResetState("sample_names", initial,
            builder$designMatrix, builder$designMatrix, builder$dataClnReactive,
            builder$removedSamples(), builder$removedSamples, builder$factors,
            builder$groups, builder$contrasts)
        expect_identical(getDIAReactiveReference(builder$dataClnReactive), workflow$data_tbl)
        expect_identical(manager$.context$hydrates, reads)
        expect_identical(resolveDIAOwnedValue(initial$data_reference), table)

        session$setInputs(selected_runs = c("s1", "s2"), assign_metadata = 2L)
        builder$removedSamples("s2")
        expected <- buildProtDesignSaveResultsPayload(builder$designMatrix(), "s2", table,
            builder$contrasts(), config, "~0+group")
        reads <- manager$.context$hydrates
        session$setInputs(formula_string = "~0+group", save_results = 1L)
        expect_identical(manager$.context$hydrates, reads + 1L)
        result <- builder$resultRv()
        expect_null(result$data_cln)
        expect_true(isSavedWorkflowReference(result$data_reference))
        expect_identical(resolveDIAOwnedValue(result$data_reference), expected$data_cln)
        expect_identical(result$design_matrix, expected$design_matrix)
        expect_identical(result$config_list, expected$config_list)
        reads <- manager$.context$hydrates
        hydrateProtDesignBuilderResults(result, workflow)
        expect_identical(workflow$data_cln, result$data_reference)
        expect_identical(manager$.context$hydrates, reads)
        persistProtDesignBuilderArtifacts(result, workflow, root)
        expected_path <- file.path(root, "expected.tab")
        utils::write.table(expected$data_cln, expected_path, sep = "\t", row.names = FALSE, quote = FALSE)
        expect_identical(readLines(file.path(root, "data_cln.tab")), readLines(expected_path))
        expect_true(all(file.exists(file.path(root, c("design_matrix.tab", "manifest.json", "config.ini")))))
        expect_identical(manager$getHistory(), "initial")
        expect_identical(diaStateOwnership(manager)$working_records, 0L)

        replacement <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "replacement", "DIA"))
        closeDIAStateAdapter(manager)
        attachDIAStateSession(replacement, session, workflow)
        workflow$state_manager <- replacement
        setDIAOwnedPayload(workflow, "data_tbl", table[1, , drop = FALSE])
        session$flushReact()
        expect_identical(getDIAReactiveReference(builder$dataClnReactive)$owner,
            replacement$.context$store$identity)
        expect_identical(builder$initialState()$design_matrix$Run, "s2")
        expect_identical(diaStateOwnership(replacement)$working_records, 0L)
    })
})

test_that("DIA design converts only necessary successors and restores their roots", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "design-convert", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    mapping <- list(run_col = "Run", quantity_col = "Quantity")
    config <- list(deAnalysisParameters = list(formula_string = "~0+group"))
    for (quantity in list(c("1", NA, "Inf"), c(1L, NA_integer_, 0L),
        factor(c("2", NA, "1"), levels = c("unused", "1", "2")))) {
        table <- data.frame(Run = c("s1", "s2", "s1"), Quantity = quantity)
        reference <- saveDIAOwnedValue(table, "input", manager)
        expected <- buildProtDesignInitialState(table, config, mapping)
        reads <- manager$.context$hydrates
        initial <- withDIAStateScope(manager, function()
            buildDIAReferenceDesignState(manager, reference, config, mapping))
        expect_identical(manager$.context$hydrates, reads + 1L)
        expect_identical(resolveDIAOwnedValue(initial$data_reference, manager), expected$data_cln)
        expect_identical(initial$design_matrix, expected$design_matrix)
        expect_false(identical(initial$data_reference, reference))
        reads <- manager$.context$hydrates
        repeated <- withDIAStateScope(manager, function()
            buildDIAReferenceDesignState(manager, reference, config, mapping))
        expect_identical(repeated, initial)
        expect_identical(manager$.context$hydrates, reads)
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
    }
    saved <- runDIANativeSessionExport(manager, list(workflow_type = "DIA"), withr::local_tempdir())
    restored <- openDIANativeSession(saved$publication$archive_path, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(restored$manager))
    closeDIAStateAdapter(manager)
    restored_initial <- restored$manager$.context$snapshot$design_initial
    reads <- restored$manager$.context$hydrates
    rebuilt <- withDIAStateScope(restored$manager, function() buildDIAReferenceDesignState(
        restored$manager, restored_initial$source_reference, config, mapping))
    expect_identical(restored$manager$.context$hydrates, reads + 1L)
    expect_identical(resolveDIAOwnedValue(rebuilt$data_reference, restored$manager), expected$data_cln)
    expect_identical(restored$manager$getHistory(), "initial")
})

test_that("consumed-reference release removes only the old owned binding", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "release", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(), "raw")
    reference <- manager$.context$snapshot$states$raw
    path <- file.path(manager$.context$store$directory, reference$payload$relative_path)
    withDIAStateScope(manager, function() {
        manager$getState()
        expect_false(releaseDIAWorkingReference(manager, reference))
        manager$saveState("next_state", object, list(), "next")
        current <- manager$.context$snapshot$states$next_state
        foreign <- reference
        foreign$owner$incarnation <- "different-owner"
        expect_error(releaseDIAWorkingReference(manager, foreign), class = "saved_entity_scope")
        reads <- manager$.context$hydrates
        withDIAStateScope(manager, function() {
            expect_true(releaseDIAWorkingReference(manager, reference))
            expect_false(exists(reference$entity_id, manager$.context$working, inherits = FALSE))
            expect_true(exists(current$entity_id, manager$.context$working, inherits = FALSE))
            expect_false(releaseDIAWorkingReference(manager, reference))
        })
        expect_identical(manager$.context$scope_depth, 1L)
        expect_identical(manager$getState(), object)
        expect_identical(manager$.context$hydrates, reads)
        expect_true(file.exists(path))
    })
    expect_identical(manager$revertToState("raw"), object)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_false(releaseDIAWorkingReference(WorkflowState$new(), reference))
})

test_that("peptide Apply releases consumed inputs before counts and Parquet preparation", {
    run_case <- function(audit_enabled) {
        manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "apply", "DIA"),
            WorkflowState$new(audit_enabled = audit_enabled))
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        manager$.context$sql_rollup <- TRUE
        initial <- getDIAStateRecord(manager, "initial")
        object <- makeDIALifetimePeptide()
        manager$saveState("raw_data_s4", object, list(), "raw")
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        config <- new.env(parent = emptyenv())
        config$config_list <- object@args
        update <- function(...) updateConfigParameter(..., env = config)
        count_original <- .countPeptideProteinGroups
        observations <- character()
        observe <- function(frame, label) {
            predecessor <- get("predecessor", frame, inherits = FALSE)
            expect_null(get("currentS4", frame, inherits = FALSE))
            expect_false(exists(predecessor$entity_id, manager$.context$working, inherits = FALSE))
            current <- manager$.context$snapshot$states[[manager$current_state]]
            expect_true(exists(current$entity_id, manager$.context$working, inherits = FALSE))
            expect_identical(readRDS(manager$.context$head_path), manager$.context$head)
            observations <<- c(observations, label)
        }
        testthat::local_mocked_bindings(
            .countPeptideProteinGroups = function(theObject) {
                if (exists("predecessor", parent.frame(), inherits = FALSE)) observe(parent.frame(), "count")
                count_original(theObject)
            })
        now <- function() as.POSIXct("2026-09-22 00:00:00", tz = "UTC")
        withDIAStateScope(manager, function() {
            result <- runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE,
                updateConfigParameterFn = update, nowFn = now)
            reads <- manager$.context$hydrates
            expect_identical(manager$getState(), result$filteredS4)
            expect_identical(manager$.context$hydrates, reads)
        })
        for (iteration in 1:2) withDIAStateScope(manager, function() {
            result <- runProteinPeptideApplyStep(workflow, 1, 2,
                updateConfigParameterFn = update, nowFn = now)
            reads <- manager$.context$hydrates
            expect_identical(manager$getState(), result$filteredS4)
            expect_identical(manager$.context$hydrates, reads)
            expect_identical(diaStateOwnership(manager)$working_records, 1L)
        })
        # KIT-011: the capacity arm's rollup-input preparation after the protein-evidence filter is gone.
        expect_identical(observations, c("count", "count", "count"))
        expect_identical(manager$revertToState("raw_data_s4"), object)
        info <- diaStateOwnership(manager)
        expect_identical(info$working_records, 0L)
        expect_identical(info$store$open_resources, 0L)
        expect_equal(info$store$bytes, 0)
        manager$states <- list(initial = initial)
        manager$current_state <- "initial"
        manager$state_history <- list("initial")
        expect_null(manager$getState())
        expect_identical(manager$getHistory(), "initial")
    }
    invisible(lapply(c(FALSE, TRUE), run_case))
})

test_that("failed peptide saves retain the active predecessor until scope cleanup", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "failure", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeDIALifetimePeptide()
    manager$saveState("raw_data_s4", object, list(), "raw")
    reference <- manager$.context$snapshot$states$raw_data_s4
    head <- manager$.context$head
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    config <- new.env(parent = emptyenv())
    config$config_list <- object@args
    update <- function(...) updateConfigParameter(..., env = config)
    withDIAStateScope(manager, function() {
        manager$getState()
        apply <- function() runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE,
            updateConfigParameterFn = update)
        testthat::with_mocked_bindings({
            expect_error(apply(), "payload failure")
        }, writeSavedWorkflowQsPayload = function(...) stop("payload failure"))
        testthat::with_mocked_bindings({
            expect_error(apply(), "head failure")
        }, diaStateBoundary = function(...) stop("head failure"))
        expect_identical(manager$.context$head, head)
        expect_identical(readRDS(manager$.context$head_path), head)
        expect_true(exists(reference$entity_id, manager$.context$working, inherits = FALSE))
        expect_false(releaseDIAWorkingReference(manager, reference))
        expect_identical(manager$getState(), object)
    })
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
})

test_that("DIA Save Design preserves legacy parser attribute cleanup without resaving initial data", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "parser-design", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(getDIAActiveManager = function() manager)
    table <- readr::read_tsv(I("Run\tQuantity\ns1\t10\ns2\t20\n"),
        col_types = readr::cols(Run = readr::col_character(), Quantity = readr::col_double()))
    attr(table, "scientific_annotation") <- "preserve as in the native builder"
    portable <- prepareDIAParserTable(table)
    mapping <- list(run_col = "Run", quantity_col = "Quantity")
    config <- list(deAnalysisParameters = list(formula_string = "~0+group"))
    reference <- saveDIAOwnedValue(portable, "import", manager)
    recordDIAImportMetadata(manager, portable,
        list(data = table, data_type = "peptide", column_mapping = mapping), reference)
    native <- buildProtDesignInitialState(portable, config, mapping)
    native$design_matrix$group <- c("control", "treated")
    expected <- buildProtDesignSaveResultsPayload(native$design_matrix, character(),
        native$data_cln, native$contrasts, config, "~0+group")
    reads <- manager$.context$hydrates
    initial <- withDIAStateScope(manager, function()
        buildDIAReferenceDesignState(manager, reference, config, mapping))
    expect_identical(initial$data_reference, reference)
    expect_identical(manager$.context$hydrates, reads)
    result <- withDIAStateScope(manager, function() runProtDesignSaveResultsObserverShell(
        native$design_matrix, character(), resolveDIAOwnedValue(initial$data_reference),
        native$contrasts, config, "~0+group", resultSetter = function(...) NULL,
        showNotification = function(...) NULL))
    actual <- resolveDIAOwnedValue(result$data_reference)
    expect_identical(actual, expected$data_cln)
    expect_identical(attributes(actual), attributes(expected$data_cln))
    expect_identical(resolveDIAOwnedValue(reference), portable)
    expect_identical(attr(resolveDIAOwnedValue(reference), "problems"), attr(portable, "problems"))
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
})

test_that("native checkpoint I/O drops owned file pages without changing payloads", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "guarded disk-backed proof root required")
    root <- withr::local_tempdir(tmpdir = dirname(receipt))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "checkpoint-pages", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$session_format <- "native_v1"
    set.seed(902006)
    payload <- as.raw(sample.int(256L, 8L * 1024L^2L, replace = TRUE) - 1L)
    reference <- saveDIAOwnedValue(payload, "page-fixture", manager)
    path <- file.path(manager$.context$store$directory, reference$payload$relative_path)
    resident_bytes <- function() {
        output <- system2("fincore", c("--json", "--bytes", "--output", "RES", "--", shQuote(path)),
            stdout = TRUE, stderr = TRUE)
        expect_null(attr(output, "status"))
        jsonlite::fromJSON(paste(output, collapse = "\n"))$fincore$res[[1L]]
    }
    expect_gt(reference$payload$bytes, 4194304)
    after_save <- resident_bytes()
    expect_lte(after_save, 4194304)
    expect_identical(digest::digest(file = path, algo = "sha256"), reference$payload$sha256)
    restored <- resolveDIAOwnedValue(reference, manager)
    after_read <- resident_bytes()
    expect_identical(restored, payload)
    expect_lte(after_read, 4194304)
    expect_true(file.exists(path))
    expect_identical(manager$getHistory(), "initial")
    head <- manager$.context$head
    testthat::with_mocked_bindings({
        expect_error(withDIAStateScope(manager, function() {
            manager$saveState("failed", NULL, list(), "injected I/O failure")
        }), "injected cache I/O failure")
    }, releaseDIANativeFileCache = function(...) stop("injected cache I/O failure"))
    expect_identical(manager$.context$head, head)
    expect_identical(manager$getHistory(), "initial")
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    jsonlite::write_json(list(payload_bytes = reference$payload$bytes, after_save = after_save,
        after_read = after_read, scientific_performance_acceptance = FALSE),
        file.path(dirname(receipt), "checkpoint-file-residency.json"), auto_unbox = TRUE)
})

test_that("DIA QC observer owns its Apply, display and error paths", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "observer-life", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(getDIAActiveManager = function() manager)
    object <- makeDIALifetimePeptide()
    manager$saveState("raw_data_s4", object, list(), "raw")
    predecessor <- getDIACurrentReference(manager)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    config <- new.env(parent = emptyenv())
    config$config_list <- object@args
    update <- function(...) updateConfigParameter(..., env = config)
    reads <- manager$.context$hydrates
    result <- runPeptideQvalueApplyObserver(workflow, 0.01, 0.01, TRUE,
        new.env(parent = emptyenv()), function(...) NULL, "proteomics", "fixture",
        runApplyStepFn = function(...) {
            expect_identical(manager$.context$scope_depth, 1L)
            runPeptideQvalueApplyStep(..., updateConfigParameterFn = update)
        }, updateOutputsFn = function(...) {
            expect_identical(manager$.context$scope_depth, 1L)
            expect_false(exists(predecessor$entity_id, manager$.context$working, inherits = FALSE))
            before <- manager$.context$hydrates
            expect_identical(manager$getState(), manager$getState())
            expect_identical(manager$.context$hydrates, before)
            "rendered"
        }, showNotificationFn = function(...) NULL, removeNotificationFn = function(...) NULL)
    expect_identical(result$status, "success")
    expect_identical(manager$.context$hydrates, reads + 1L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    head <- manager$.context$head
    failed <- runPeptideQvalueApplyObserver(workflow, 0.01, 0.01, TRUE,
        new.env(parent = emptyenv()), function(...) NULL, "proteomics", "fixture",
        runApplyStepFn = function(...) { manager$getState(); stop("injected observer failure") },
        showNotificationFn = function(...) NULL, removeNotificationFn = function(...) NULL,
        logErrorFn = function(...) NULL)
    expect_identical(failed$status, "error")
    expect_identical(manager$.context$head, head)
    expect_identical(manager$.context$scope_depth, 0L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
})

test_that("DIA style regeneration reuses each saved input within one operation", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "style-life", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(getDIAActiveManager = function() manager)
    norm <- new.env(parent = emptyenv())
    norm$pre_norm_qc_generated <- FALSE
    norm$normalization_complete <- norm$ruv_complete <- TRUE
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    setDIAOwnedPayload(norm, "normalized_protein_obj", object)
    norm$ruv_normalized_obj <- norm$normalized_protein_obj
    reads <- manager$.context$hydrates
    callback <- function(value) {
        expect_identical(manager$.context$scope_depth, 1L)
        expect_identical(value, object)
        expect_identical(getDIAOwnedPayload(norm, "normalized_protein_obj"), value)
    }
    regenerateProtNormQcForAestheticChange(norm, function(...) stop("pre-QC not requested"),
        callback, callback, messageFn = function(...) NULL)
    expect_identical(manager$.context$hydrates, reads + 1L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$scope_depth, 0L)
})

test_that("DIA normalization releases obsolete inputs before QC delivery", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "normalization-life", "DIA"))
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(getDIAActiveManager = function() manager)
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    rlang::local_bindings(config_list = object@args, .env = .GlobalEnv)
    manager$saveState("protein_replicate_filtered", object, list(), "input")
    norm <- new.env(parent = emptyenv())
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    paths <- list(source_dir = withr::local_tempdir(), protein_qc_dir = withr::local_tempdir())
    action_frame <- function() {
        frames <- Filter(function(frame) exists("input_reference", frame, inherits = FALSE) &&
            exists("run_context", frame, inherits = FALSE), sys.frames())
        expect_length(frames, 1L)
        frames[[1L]]
    }
    result <- runProtNormNormalizationWorkflow(list(norm_method = "cyclicloess", ruv_mode = "skip"),
        workflow, norm, paths, "proteomics", "fixture", checkMemoryUsageFn = function(...) NULL,
        generatePostNormalizationQcFn = function(...) {
            frame <- action_frame()
            expect_null(get("run_context", frame, inherits = FALSE))
            expect_null(get("current_s4", frame, inherits = FALSE))
        }, generateRuvCorrectedQcFn = function(...) NULL, getRuvGroupingVariableFn = function() "group",
        generateCompositeQcFigureFn = function(...) {
            frame <- action_frame()
            for (name in c("normalized_s4", "ruv_corrected_s4", "ruv_corrected_s4_clean", "step6_object"))
                expect_null(get(name, frame, inherits = FALSE))
            expect_identical(diaStateOwnership(manager)$working_records, 1L)
        }, withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) NULL, showNotificationFn = function(...) NULL,
        handleErrorFn = function(error) stop(error), messageFn = function(...) NULL)
    expect_true(result)
    expect_identical(manager$current_state, "normalized")
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$getState("protein_replicate_filtered"), object)
})

test_that("current IQ path releases staging data before fitting and source after commit", {
    for (audited in c(FALSE, TRUE)) {
        manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "iq-current", "DIA"),
            WorkflowState$new(audit_enabled = audited))
        withr::defer(closeDIAStateAdapter(manager))
        peptide <- module_ci_prot_peptide_rollup_object()
        manager$saveState("imputed", peptide, list(), "input")
        predecessor <- getDIACurrentReference(manager)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        paths <- list(peptide_qc_dir = withr::local_tempdir(), protein_qc_dir = withr::local_tempdir())
        processor <- function(...) {
            frame <- parent.frame()
            expect_null(get("peptideDataForIq", frame, inherits = FALSE))
            expect_null(get("sampleIndex", frame, inherits = FALSE))
            expect_true(exists(predecessor$entity_id, manager$.context$working, inherits = FALSE))
            module_ci_prot_iq_process_stub(...)
        }
        result <- withDIAStateScope(manager, function() runProteinIqRollupApplyStep(workflow, paths,
            processLongFormatFn = processor,
            readTsvFn = function(file, .name_repair) readr::read_tsv(file, show_col_types = FALSE),
            captureCheckpointFn = function(...) {
                expect_null(get("peptideS4", parent.frame(), inherits = FALSE))
                expect_null(get("proteinLog2Quant", parent.frame(), inherits = FALSE))
                expect_false(exists(predecessor$entity_id, manager$.context$working, inherits = FALSE))
                expect_identical(manager$current_state, "protein_s4_created")
            }, showNotificationFn = function(...) NULL))
        # The established codec strips parser pointers, not scientific attributes.
        expect_identical(manager$getState(), prepareDIAStateRecord(list(data = result$proteinObj))$data)
        expect_identical(manager$getState("imputed"), peptide)
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        closeDIAStateAdapter(manager)
    }
})

test_that("IQ input preparation has unchanged TSV contents and no surviving preparation frame", {
    skip_if_not(exists("writeProteinIqInput", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 combined candidate; retained for its installed trial")
    peptide <- data.frame(protein = c("P1", "P2"), sample = c("two", "one"),
        quantity = c(NA_real_, Inf), annotation = c("alpha", "beta"), Q.Value = c(1, 1))
    mapping <- data.frame(Original = c("one", "two"), Alias = c("S_001", "S_002"))
    expected <- peptide
    expected$Q.Value <- 0.0009
    expected$PG.Q.Value <- 0.009
    expected$quantity[is.na(expected$quantity)] <- 0
    expected$sample <- mapping$Alias[match(as.character(expected$sample), mapping$Original)]
    root <- withr::local_tempdir()
    expected_path <- file.path(root, "expected.tsv")
    actual_path <- file.path(root, "actual.tsv")
    vroom::vroom_write(expected, expected_path)
    frame_ref <- NULL
    writer <- function(value, path) {
        frame_ref <<- rlang::new_weakref(parent.frame())
        expect_identical(value, expected)
        vroom::vroom_write(value, path)
    }
    expect_null(writeProteinIqInput(peptide, mapping, "sample", "quantity", actual_path, writer))
    gc(full = TRUE) # Ownership assertion only; production adds no collection.
    expect_null(rlang::wref_key(frame_ref))
    expect_identical(readBin(actual_path, "raw", file.size(actual_path)),
        readBin(expected_path, "raw", file.size(expected_path)))
    expect_identical(peptide$sample, c("two", "one"))
    expect_true(is.na(peptide$quantity[[1L]]))
    expect_error(writeProteinIqInput(peptide, mapping[1L, ], "sample", "quantity", actual_path, writer),
        "runs absent from the design matrix")
})

test_that("IQ keeps its source through commit and releases it before checkpoint output", {
    skip_if_not(exists("writeProteinIqInput", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 combined candidate; retained for its installed trial")
    run_case <- function(audit_enabled) {
        manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "iq-life", "DIA"),
            WorkflowState$new(audit_enabled = audit_enabled))
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        peptide <- module_ci_prot_peptide_rollup_object()
        manager$saveState("imputed", peptide, list(), "input")
        reference <- manager$.context$snapshot$states$imputed
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        paths <- list(peptide_qc_dir = withr::local_tempdir(), protein_qc_dir = withr::local_tempdir())
        frame_ref <- NULL
        writer <- function(value, path) {
            frame_ref <<- rlang::new_weakref(parent.frame())
            vroom::vroom_write(value, path)
        }
        processor <- function(...) {
            gc(full = TRUE) # Fixture checks reachability, not performance.
            expect_null(rlang::wref_key(frame_ref))
            expect_true(exists(reference$entity_id, manager$.context$working, inherits = FALSE))
            module_ci_prot_iq_process_stub(...)
        }
        apply <- function() runProteinIqRollupApplyStep(workflow, paths,
            writeTsvFn = writer, processLongFormatFn = processor,
            readTsvFn = function(file, .name_repair) readr::read_tsv(file, show_col_types = FALSE),
            captureCheckpointFn = function(...) {
                expect_identical(manager$current_state, "protein_s4_created")
                expect_null(get("peptideS4", parent.frame(), inherits = FALSE))
                expect_false(exists(reference$entity_id, manager$.context$working, inherits = FALSE))
            }, showNotificationFn = function(...) NULL, logInfoFn = function(...) NULL,
            logWarnFn = function(...) NULL)
        withDIAStateScope(manager, function() {
            head <- manager$.context$head
            testthat::with_mocked_bindings({
                expect_error(apply(), "injected IQ commit failure")
            }, diaStateBoundary = function(...) stop("injected IQ commit failure"))
            expect_identical(manager$.context$head, head)
            expect_identical(manager$current_state, "imputed")
            expect_identical(manager$getState(), peptide)
            result <- apply()
            expect_identical(manager$getState(), result$proteinObj)
        })
        expect_identical(diaStateOwnership(manager)$working_records, 0L)
        expect_identical(manager$revertToState("imputed"), peptide)
    }
    invisible(lapply(c(FALSE, TRUE), run_case))
})
