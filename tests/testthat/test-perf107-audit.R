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

test_that("lifecycle classification consumes exact receipts and service outcomes", {
    harness <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), harness)
    decide <- harness$perf107LifecycleSucceeded
    for (scenario in c("diagnostic", "science")) {
        receipt <- list(schema = 1L, scenario = scenario,
            status = if (scenario == "science") "science_complete" else "topology_complete")
        expect_true(decide(scenario, 0L, receipt, NULL))
        expect_false(decide(scenario, 1L, receipt, NULL))
        expect_false(decide(scenario, 0L, NULL, NULL))
        expect_false(decide(scenario, 0L, receipt, list(reason = "memory_guard")))
        receipt$status <- "failed"
        expect_false(decide(scenario, 0L, receipt, NULL))
    }
    expect_false(decide("topology", 0L, list(status = "topology_complete"), NULL))
    expect_false(decide("guard-stop", 1L, NULL, list(schema = 1L, reason = "memory_guard")))
    expect_true(decide("guard-stop", 1L, NULL,
        list(schema = 1L, reason = "injected_guard_stop")))
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

test_that("reference native control really calls the native scientific method", {
    rlang::local_bindings(config_list = list(filterMinNumPeptidesPerSample =
        list(peptides_per_sample_cutoff = 1)),
        .env = .GlobalEnv)
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("input", object, list(), "input")
    root <- withr::local_tempdir()
    manager <- stageDIAGUIReferenceOwner(native, list(source_dir = root))
    withr::defer(closeDIAStateAdapter(manager))
    workflow <- list2env(list(state_manager = manager, qc_params = list()))
    expected <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = 1, core_utilisation = NA)
    testthat::local_mocked_bindings(
        runDIAReferencePeptideGroups = function(...) stop("incorrect grouped control"),
        .package = "MultiScholaR"
    )
    result <- runWorkflowNativeSampleControl(workflow, 1)
    expect_identical(result$filteredS4@peptide_data, expected@peptide_data)
    expect_identical(manager$current_state, "sample_filtered")
    expect_gt(diaStateOwnership(manager)$full_s4_reads, 0L)
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

test_that("small scientific verification rejects changed values and evidence", {
    harness <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), harness)
    expected <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    expect_true(harness$verifyPerf107ScientificState(expected, expected, "filter"))
    changed <- expected
    changed@peptide_data$Precursor.Quantity[[1L]] <- -1
    expect_error(harness$verifyPerf107ScientificState(changed, expected, "filter"),
        "peptide_data")
    changed <- expected
    changed@args$filter <- list(support_table = data.frame(support = 1L))
    expect_error(harness$verifyPerf107ScientificState(changed, expected, "filter"),
        "support_table")
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
    application$config$qualification <- list(programme = "perf107-phase-a")
    action <- beginWorkflowGUIAction(application, WorkflowState$new(audit_enabled = FALSE),
        "minimum_peptides_per_sample", list(experiment_label = "test",
            experiment_paths = list(source_dir = root), omic_type = "proteomics"))
    expect_error(finishWorkflowGUIAction(action, "complete", list(request_id = action$id,
        content_id = "expected", plot_ready = TRUE, stable_seconds = 0.2)),
        class = "workflow_execution_unsupported")
    expect_null(application$active)
    expect_true(action$closed)
})

test_that("Apply retries pending publication without dispatching science again", {
    manager <- new.env()
    class(manager) <- c("DIAStateAdapter", "environment")
    owner <- new.env()
    owner$manager <- manager
    owner$job <- NULL
    manager$.context <- list(store = list(execution_owner = owner))
    manager$current_state <- "sample_filtered"
    workflow <- shiny::reactiveValues(state_manager = manager, qc_params = list())
    dispatches <- retries <- 0L
    completed <- character()
    reference <- list(entity_id = "accepted")
    testthat::local_mocked_bindings(
        workflowPeptideSampleQualification = function(...) list(entry = "reference", variant = "physical"),
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
        readManagedPeptideSampleResult = function(...) list(resultText = "accepted", stateReference = reference),
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
