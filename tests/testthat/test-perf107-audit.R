test_that("DIA owner replacement retains the existing QC layout", {
    renders <- 0L
    testthat::local_mocked_bindings(
        mod_prot_qc_peptide_ui = function(...) {
            renders <<- renders + 1L
            shiny::div("same scientific workflow")
        },
        mod_prot_qc_protein_ui = function(...) shiny::div(),
        .package = "MultiScholaR"
    )
    workflow <- shiny::reactiveValues(state_manager = list(workflow_type = "DIA"))
    shiny::testServer(mod_prot_qc_server,
        args = list(workflow_data = workflow, experiment_paths = list(),
            omic_type = "proteomics", experiment_label = "audit",
            qc_trigger = function() FALSE), {
        session$flushReact()
        first <- output$dynamic_qc_tabs
        expect_identical(renders, 1L)
        workflow$state_manager <- list(workflow_type = "DIA", adopted = TRUE)
        session$flushReact()
        expect_identical(output$dynamic_qc_tabs, first)
        expect_identical(renders, 1L)
    })
    native <- WorkflowState$new(audit_enabled = FALSE)
    workflow <- shiny::reactiveValues(state_manager = native)
    trigger <- shiny::reactiveVal(FALSE)
    testthat::local_mocked_bindings(
        mod_prot_qc_peptide_server = function(...) NULL,
        mod_prot_qc_protein_server = function(...) NULL,
        .package = "MultiScholaR"
    )
    shiny::testServer(mod_prot_qc_server,
        args = list(workflow_data = workflow, experiment_paths = list(),
            omic_type = "proteomics", experiment_label = "audit", qc_trigger = trigger), {
        session$flushReact()
        first <- output$dynamic_qc_tabs
        expect_identical(qc_workflow_type(), "LFQ")
        native$setWorkflowType("DIA")
        trigger(TRUE)
        session$flushReact()
        expect_identical(qc_workflow_type(), "DIA")
        expect_false(identical(output$dynamic_qc_tabs, first))
    })
})

test_that("browser acknowledgement checks the accepted visible content", {
    node <- Sys.getenv("MULTISCHOLAR_PERF107_NODE", Sys.which("node"))
    skip_if(!nzchar(node), "Node is required for the JavaScript contract")
    script <- tempfile(fileext = ".js")
    withr::defer(unlink(script))
    writeLines(workflowPeptideSampleAckScript(), script)
    result <- system2(node, shQuote(c(test_path("fixtures", "perf107-ack-contract.cjs"),
        script)), stdout = TRUE, stderr = TRUE)
    expect_identical(attr(result, "status") %||% 0L, 0L, info = paste(result, collapse = "\n"))
    expect_true(any(grepl("ack lifecycle contract passed", result, fixed = TRUE)))
})

test_that("the supervisor cleans only its validated request-owned slice", {
    runner <- readLines(test_path("..", "..", "tools", "profiling",
        "run_workflow_checkpoint_proof.sh"))
    start <- which(runner == "checkpointStopOwnedApplication() {")
    end <- which(seq_along(runner) > start & runner == "}")[[1L]]
    root <- withr::local_tempdir()
    dir.create(file.path(root, "application"))
    marker <- file.path(root, "application", "application-slice.txt")
    unit <- "multischolarwfx0123456789abcdef.slice"
    shell <- tempfile(fileext = ".sh")
    withr::defer(unlink(shell))
    # A metadata-only fake systemctl exercises target validation and stop order.
    writeLines(c("set -eu", "root=$1", "log=$2",
        'systemctl() { printf "%s\\n" "$*" >> "$log"; if [[ "$*" == *--value* ]]; then echo inactive; fi; }',
        runner[start:end], 'checkpointStopOwnedApplication "$root"'), shell)
    log <- file.path(root, "commands.txt")
    writeLines(unit, marker)
    expect_identical(system2("bash", shQuote(c(shell, root, log))), 0L)
    commands <- readLines(log)
    expect_true(any(commands == paste("--user stop", unit)))
    expect_length(commands, 3L)
    writeLines("user.slice", marker)
    before <- readLines(log)
    expect_identical(system2("bash", shQuote(c(shell, root, log))), 2L)
    expect_identical(readLines(log), before)
})

test_that("receipt write failure releases the GUI reservation", {
    root <- withr::local_tempdir()
    application <- WorkflowExecutionApplication$new()
    application$config <- list(root = root, stem = "multischolarwfxaaaaaaaaaaaaaaaa")
    application$cgroup <- root
    application$pid <- Sys.getpid()
    peak_path <- file.path(root, "peak")
    writeLines("0", peak_path)
    testthat::local_mocked_bindings(
        checkWorkflowApplication = function(...) invisible(TRUE),
        openWorkflowPeakWindow = function(...) file(peak_path, "r+"),
        closeWorkflowPeakWindow = function(con) { close(con); 123 },
        workflowCgroupAccounting = function(...) list(cpu_seconds = 0),
        workflowMemoryCounter = function(...) 123,
        writeWorkflowExecutionMessage = function(...) stop("injected receipt failure"),
        .package = "MultiScholaR"
    )
    action <- beginWorkflowGUIAction(application, WorkflowState$new(audit_enabled = FALSE),
        "minimum_peptides_per_sample", list(experiment_label = "test",
            experiment_paths = list(source_dir = root), omic_type = "proteomics"))
    expect_error(finishWorkflowGUIAction(action, "cancelled"), "injected receipt failure")
    expect_null(application$active)
    expect_null(action$admission)
    expect_true(action$closed)
    lock <- filelock::lock(file.path(root, ".heavy-action.lock"), timeout = 0)
    expect_false(is.null(lock))
    filelock::unlock(lock)
})

test_that("Apply retries pending publication without dispatching science again", {
    manager <- new.env()
    class(manager) <- c("DIAStateAdapter", "environment")
    owner <- new.env()
    owner$manager <- manager
    owner$job <- NULL
    manager$.context <- list(store = list(execution_owner = owner), reference_states = TRUE)
    manager$current_state <- "sample_filtered"
    workflow <- shiny::reactiveValues(state_manager = manager, qc_params = list())
    dispatches <- retries <- 0L
    completed <- character()
    reference <- list(entity_id = "accepted")
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = withr::local_tempdir())
    testthat::local_mocked_bindings(
        # KIT-008: the dispatch rule is the managed application on the session, not a qualification.
        workflowSessionApplication = function(...) application,
        usesDIAReferenceStates = function(...) TRUE,
        beginWorkflowGUIAction = function(...) list2env(list(id = "gui-one", closed = FALSE)),
        finishWorkflowGUIAction = function(reservation, status, ...) {
            completed <<- c(completed, status)
            reservation$closed <- TRUE
            list()
        },
        getDIAStateReference = function(...) reference,
        runWorkflowAction = function(...) { dispatches <<- dispatches + 1L; "worker-one" },
        pollWorkflowAction = function(...) list(status = "publication_pending"),
        resumeWorkflowAction = function(...) { retries <<- retries + 1L; list(state_reference = reference) },
        readManagedPeptideResult = function(...) list(resultText = "accepted\n", stateReference = reference),
        withDIAStateScope = function(manager, operation) operation(),
        updatePeptideSampleOutputs = function(output, samplePlot, ...) {
            output$sample_results <- shiny::renderText("accepted")
            samplePlot(grid::textGrob("accepted"))
        },
        .package = "MultiScholaR"
    )
    shiny::testServer(mod_prot_qc_peptide_sample_server,
        args = list(workflow_data = workflow, omic_type = "proteomics",
            experiment_label = "audit", experiment_paths = list()), {
        session$setInputs(min_peptides_per_sample = 100, apply_sample_filter = 1)
        session$flushReact()
        expect_identical(dispatches, 1L)
        expect_identical(retries, 0L)
        session$setInputs(apply_sample_filter = 2)
        session$flushReact()
        expect_identical(dispatches, 1L)
        expect_identical(retries, 1L)
        action <- managed_action()
        session$setInputs(workflow_ack = "malformed")
        expect_identical(managed_action()$phase, "awaiting_browser")
        expect_identical(managed_action()$expected_text, "accepted")
        session$setInputs(workflow_ack = list(request_id = "wrong", content_id = action$content_id,
            plot_ready = TRUE, stable_ms = 200))
        expect_identical(managed_action()$phase, "awaiting_browser")
        for (invalid in list(NA_real_, Inf, numeric(), c(200, 201))) {
            session$setInputs(workflow_ack = list(request_id = action$reservation$id,
                content_id = action$content_id, plot_ready = TRUE, stable_ms = invalid))
            expect_identical(managed_action()$phase, "awaiting_browser")
        }
        workflow$state_manager <- new.env()
        session$flushReact()
        expect_null(managed_action())
        expect_identical(completed, c("publication_pending", "stale"))
    })
})

test_that("adoption failure preserves live assays and successful adoption releases holders", {
    root <- withr::local_tempdir()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("h0", object, list(), "h0")
    object@peptide_data$Precursor.Quantity[[1L]] <- 7
    native$saveState("h1", object, list(), "h1")
    original <- digest::digest(native$states)
    workflow <- shiny::reactiveValues(state_manager = native,
        data_tbl = object@peptide_data, data_cln = object@peptide_data,
        peptide_data = NULL, data_type = "peptide",
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
    shiny::testServer(function(input, output, session) {}, {
        cleaned <- newDIASavedReactiveVal(object@peptide_data, "design-cleaned-data",
            snapshot_field = "design_current")
        staged <- stageDIAGUIReferenceOwner(native, list(source_dir = root))
        testthat::with_mocked_bindings({
            expect_error(adoptDIAGUIReferenceOwner(staged, workflow, session), "injected swap")
        }, workflowExecutionBoundary = function(stage) stop("injected swap"),
            .package = "MultiScholaR")
        expect_identical(workflow$state_manager, native)
        expect_identical(workflow$data_tbl, object@peptide_data)
        expect_true(staged$.context$store$closed)
        expect_identical(digest::digest(native$states), original)
        staged <- stageDIAGUIReferenceOwner(native, list(source_dir = root))
        adoptDIAGUIReferenceOwner(staged, workflow, session)
        expect_true(isSavedWorkflowReference(workflow$data_tbl))
        expect_true(isSavedWorkflowReference(workflow$data_cln))
        metadata <- getDIAImportMetadata(staged)
        expect_identical(metadata$reference_id, workflow$data_tbl$entity_id)
        expect_identical(metadata$column_mapping, workflow$column_mapping)
        adopted_head <- staged$.context$head
        expect_true(is.list(buildDIAReferenceDesignState(staged, workflow$data_tbl,
            list(globalParameters = list(workflow_type = "DIA")), workflow$column_mapping)))
        preserve_design <- function(...) stop("adoption must not reset the design")
        observer <- registerProtDesignInitialStateSyncObserver(input, session,
            initialState = function() buildDIAReferenceDesignState(staged, workflow$data_tbl,
                list(globalParameters = list(workflow_type = "DIA")), workflow$column_mapping),
            dataTbl = function() workflow$data_tbl,
            designMatrix = preserve_design,
            dataClnReactive = cleaned,
            groups = preserve_design, factors = preserve_design, contrasts = preserve_design,
            removedSamples = preserve_design)
        session$flushReact()
        expect_identical(getDIAReactiveReference(cleaned), workflow$data_cln)
        expect_identical(staged$.context$head, adopted_head)
        expect_null(staged$.context$design_adoption_pending)
        observer$destroy()
        expect_identical(resolveDIAOwnedValue(workflow$data_tbl, staged), object@peptide_data)
        expect_identical(staged$getState("h1")@peptide_data, object@peptide_data)
        expect_identical(staged$getState("h0")@peptide_data$Precursor.Quantity[[1L]],
            native$getState("h0")@peptide_data$Precursor.Quantity[[1L]])
        expect_identical(digest::digest(native$states), original)
        export_dir <- file.path(root, "export"); dir.create(export_dir)
        saved <- runDIANativeSessionExport(staged, list(workflow_type = "DIA"),
            export_dir, extra_references = list(data_tbl = workflow$data_tbl,
                data_cln = workflow$data_cln))
        relocated <- file.path(root, "relocated"); dir.create(relocated)
        restored <- openDIANativeSession(saved$publication$archive_path, relocated)
        expect_identical(restored$manager$getState("h1")@peptide_data, object@peptide_data)
        expect_identical(restored$manager$getState("h0")@peptide_data,
            native$getState("h0")@peptide_data)
        selectDIAStateReference(restored$manager, "h0")
        expect_identical(restored$manager$current_state, "h0")
        expect_identical(resolveDIAOwnedValue(restored$extra_references$data_tbl,
            restored$manager), object@peptide_data)
        closeDIAStateAdapter(restored$manager)
        closeDIAStateAdapter(staged)
    })
})

test_that("the precursor rollup takes the managed worker under a managed application and the in-process declaration otherwise (KIT-008)", {
    manager <- new.env()
    class(manager) <- c("DIAStateAdapter", "environment")
    owner <- new.env()
    owner$manager <- manager
    owner$job <- NULL
    manager$.context <- list(store = list(execution_owner = owner), reference_states = TRUE)
    manager$current_state <- "protein_peptide_filtered"
    workflow <- shiny::reactiveValues(state_manager = manager, qc_params = list())
    reference <- list(entity_id = "accepted")
    dispatched <- list()
    in_process <- 0L
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = withr::local_tempdir())
    managed <- application
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(...) managed,
        usesDIAReferenceStates = function(...) TRUE,
        beginWorkflowGUIAction = function(application, manager, operation, publication) {
            expect_identical(operation, "precursor_rollup")
            expect_identical(publication$omic_type, "proteomics")
            list2env(list(id = "gui-rollup", closed = FALSE))
        },
        finishWorkflowGUIAction = function(reservation, status, ...) { reservation$closed <- TRUE; list() },
        getDIAStateReference = function(...) reference,
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            dispatched[[length(dispatched) + 1L]] <<- list(operation = operation, parameters = parameters,
                context = completion_context)
            "worker-rollup"
        },
        pollWorkflowAction = function(...) list(status = "complete"),
        resumeWorkflowAction = function(...) list(state_reference = reference),
        readManagedPeptideResult = function(...) list(resultText = "rolled\n", stateReference = reference),
        withDIAStateScope = function(manager, operation) operation(),
        updatePeptideRollupOutputs = function(output, rollupPlot, ...) {
            output$rollup_results <- shiny::renderText("rolled")
            rollupPlot(grid::textGrob("rolled"))
        },
        runDIAReferencePrecursorRollup = function(workflow) { in_process <<- in_process + 1L; list(resultText = "in process\n") },
        .package = "MultiScholaR"
    )
    shiny::testServer(mod_prot_qc_peptide_rollup_server,
        args = list(workflow_data = workflow, omic_type = "proteomics",
            experiment_label = "audit", experiment_paths = list()), {
        session$setInputs(apply_rollup = 1)
        session$flushReact()
        expect_length(dispatched, 1L)
        expect_identical(dispatched[[1L]]$operation, "precursor_rollup")
        expect_identical(dispatched[[1L]]$parameters, list())
        expect_identical(dispatched[[1L]]$context$peptide_qc, list())
        expect_identical(in_process, 0L)
        action <- managed_step$active()
        expect_identical(action$phase, "awaiting_browser")
        expect_identical(action$expected_text, "rolled")
        session$setInputs(workflow_ack = list(request_id = action$reservation$id,
            content_id = action$content_id, plot_ready = TRUE, stable_ms = 250))
        expect_null(managed_step$active())
        # No managed application on the session: the same declaration runs in process.
        managed <<- NULL
        session$setInputs(apply_rollup = 2)
        session$flushReact()
        expect_length(dispatched, 1L)
        expect_identical(in_process, 1L)
    })
})

test_that("the managed completion carries the bus's QC parameters and prior progress, and the result reader merges them (KIT-008 audit PE-02)", {
    progress <- methods::new("FilteringProgress")
    progress@steps <- "3_protein_peptide_filtered"
    progress@proteins <- 5
    workflow <- shiny::reactiveValues(qc_params = list(peptide_qc = list(qvalue = list(threshold = 0.01),
        protein_peptide = list(min = 2))), filtering_progress = progress)
    context <- workflowCompletionContext(workflow)
    expect_identical(names(context$peptide_qc), c("qvalue", "protein_peptide"))
    expect_identical(context$prior_progress$steps, "3_protein_peptide_filtered")
    expect_identical(validateWorkflowCompletionContext(context), context)
    expect_error(validateWorkflowCompletionContext(list(peptide_qc = "x")), class = "workflow_execution_invalid")
    manager <- structure(list(), class = "DIAStateAdapter")
    completion <- newWorkflowCompletionWorkflow(manager, context)
    expect_identical(completion$qc_params$peptide_qc, context$peptide_qc)
    expect_s4_class(completion$filtering_progress, "FilteringProgress")
    expect_identical(completion$filtering_progress@steps, "3_protein_peptide_filtered")
    bare <- newWorkflowCompletionWorkflow(manager, NULL)
    expect_identical(bare$qc_params, list(peptide_qc = list()))
    expect_null(bare$filtering_progress)
    # The result reader merges the step's entry into the bus instead of replacing the bus.
    testthat::local_mocked_bindings(
        readDIAStateMetadata = function(...) list(data = list(products = list(qc_progress = "product"))),
        hydrateSavedWorkflowEntity = function(...) list(value = list(resultText = "ok", qcParams = list(peptide_qc = list(
            qvalue = list(threshold = 0.01), protein_peptide = list(min = 2), sample_filter = list(cutoff = 100))))),
        checkSavedWorkflowConsumerPayload = function(...) invisible(TRUE),
        .package = "MultiScholaR")
    reference <- list(entity_id = "accepted")
    value <- readManagedPeptideResult(list(.context = list(store = NULL)), list(state_reference = reference), workflow)
    expect_identical(value$stateReference, reference)
    expect_identical(names(shiny::isolate(workflow$qc_params)$peptide_qc), c("qvalue", "protein_peptide", "sample_filter"))
    expect_identical(shiny::isolate(workflow$qc_params)$peptide_qc$qvalue$threshold, 0.01)
})

test_that("a managed launch with an in-memory owner runs the declaration in process (KIT-008 audit PE-04)", {
    native <- WorkflowState$new(audit_enabled = FALSE)
    workflow <- shiny::reactiveValues(state_manager = native, qc_params = list())
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = withr::local_tempdir())
    dispatched <- 0L
    in_process <- 0L
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(...) application,
        beginWorkflowGUIAction = function(...) { dispatched <<- dispatched + 1L; stop("must not reserve") },
        runPeptideRollupApplyObserver = function(...) { in_process <<- in_process + 1L; list(status = "success") },
        .package = "MultiScholaR")
    shiny::testServer(mod_prot_qc_peptide_rollup_server,
        args = list(workflow_data = workflow, omic_type = "proteomics", experiment_label = "audit", experiment_paths = list()), {
        session$setInputs(apply_rollup = 1)
        session$flushReact()
        expect_identical(dispatched, 0L)
        expect_identical(in_process, 1L)
    })
})

# KIT-014: every peptide QC module dispatches its registered operation to the managed worker under a managed
# application and runs the same declaration in process otherwise (the KIT-008 rollup pattern above).
checkManagedPeptideModule <- function(spec) {
    manager <- new.env()
    class(manager) <- c("DIAStateAdapter", "environment")
    owner <- new.env()
    owner$manager <- manager
    owner$job <- NULL
    manager$.context <- list(store = list(execution_owner = owner), reference_states = TRUE)
    manager$current_state <- "raw_data_s4"
    workflow <- shiny::reactiveValues(state_manager = manager, qc_params = list())
    reference <- list(entity_id = "accepted")
    dispatched <- list()
    in_process <- 0L
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = withr::local_tempdir())
    managed <- application
    mocks <- list(
        workflowSessionApplication = function(...) managed,
        usesDIAReferenceStates = function(...) TRUE,
        beginWorkflowGUIAction = function(application, manager, operation, publication) {
            expect_identical(operation, spec$operation)
            expect_identical(publication$omic_type, "proteomics")
            list2env(list(id = paste0("gui-", spec$operation), closed = FALSE))
        },
        finishWorkflowGUIAction = function(reservation, status, ...) { reservation$closed <- TRUE; list() },
        getDIAStateReference = function(...) reference,
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            dispatched[[length(dispatched) + 1L]] <<- list(operation = operation, parameters = parameters, context = completion_context)
            paste0("worker-", operation)
        },
        pollWorkflowAction = function(...) list(status = "complete"),
        resumeWorkflowAction = function(...) list(state_reference = reference),
        readManagedPeptideResult = function(...) list(resultText = "done\n", stateReference = reference),
        withDIAStateScope = function(manager, operation) operation())
    mocks[[spec$outputs]] <- function(output, ...) {
        arguments <- list(...)
        output[[spec$text]] <- shiny::renderText("done")
        arguments[[spec$plot]](grid::textGrob("done"))
    }
    mocks[[spec$in_process]] <- function(...) { in_process <<- in_process + 1L; list(resultText = "in process\n") }
    do.call(testthat::local_mocked_bindings, c(mocks, list(.package = "MultiScholaR", .env = environment())))
    shiny::testServer(spec$server,
        args = list(workflow_data = workflow, omic_type = "proteomics", experiment_label = "audit", experiment_paths = list()), {
        do.call(session$setInputs, c(spec$inputs, stats::setNames(list(1), spec$apply)))
        session$flushReact()
        expect_length(dispatched, 1L)
        expect_identical(dispatched[[1L]]$operation, spec$operation)
        expect_identical(dispatched[[1L]]$parameters, spec$expected)
        expect_identical(in_process, 0L)
        action <- managed_step$active()
        expect_identical(action$phase, "awaiting_browser")
        expect_identical(action$expected_text, "done")
        session$setInputs(workflow_ack = list(request_id = action$reservation$id,
            content_id = action$content_id, plot_ready = TRUE, stable_ms = 250))
        expect_null(managed_step$active())
        managed <<- NULL
        do.call(session$setInputs, stats::setNames(list(2), spec$apply))
        session$flushReact()
        expect_length(dispatched, 1L)
        expect_identical(in_process, 1L)
    })
    invisible(TRUE)
}

test_that("the KIT-014 peptide QC modules take the managed worker under a managed application and the in-process declaration otherwise", {
    cutoff <- MultiScholaR:::.diaNnIdentificationQvalueCutoff
    specs <- list(
        list(server = mod_prot_qc_peptide_qvalue_server, apply = "apply_qvalue_filter", inputs = list(),
            operation = "peptide_qvalue_filter", in_process = "runPeptideQvalueApplyObserver", outputs = "updatePeptideQvalueOutputs",
            text = "qvalue_results", plot = "qvaluePlot",
            expected = list(qvalue_threshold = cutoff, global_qvalue_threshold = cutoff, global_pg_qvalue_threshold = cutoff,
                proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report")),
        list(server = mod_prot_qc_peptide_protein_server, apply = "apply_protein_peptide_filter",
            inputs = list(min_peptides_per_protein = 2, min_peptidoforms_per_protein = 1),
            operation = "peptide_protein_evidence_filter", in_process = "runProteinPeptideApplyObserver", outputs = "updateProteinPeptideOutputs",
            text = "protein_peptida_results", plot = "proteinPeptidePlot",
            expected = list(num_peptides_per_protein_thresh = 2, num_peptidoforms_per_protein_thresh = 1)),
        list(server = mod_prot_qc_peptide_intensity_server, apply = "apply_intensity_filter",
            inputs = list(use_strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1, intensity_cutoff_percentile = 0.1),
            operation = "peptide_intensity_filter", in_process = "runPeptideIntensityApplyObserver", outputs = "updatePeptideIntensityOutputs",
            text = "intensity_results", plot = "intensityPlot",
            expected = list(strict_mode = FALSE, min_reps_per_group = 1, min_groups = 1, peptides_intensity_cutoff_percentile = 0.1,
                grouping_variable = "group")),
        list(server = mod_prot_qc_peptide_replicate_server, apply = "apply_replicate_filter",
            inputs = list(replicate_group_column = "replicates"),
            operation = "peptide_replicate_filter", in_process = "runPeptideReplicateApplyObserver", outputs = "updatePeptideReplicateOutputs",
            text = "replicate_results", plot = "replicatePlot", expected = list(replicate_group_column = "replicates")),
        list(server = mod_prot_qc_peptide_impute_server, apply = "apply_imputation",
            inputs = list(proportion_missing_values = 0.5, exclusion_column = ""),
            operation = "peptide_imputation", in_process = "runPeptideImputationApplyObserver", outputs = "updatePeptideImputationOutputs",
            text = "imputation_results", plot = "imputationPlot", expected = list(proportion_missing_values = 0.5)))
    for (spec in specs) checkManagedPeptideModule(spec)
})
