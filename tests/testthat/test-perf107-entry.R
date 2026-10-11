test_that("small entry fixtures refuse ordinary, active and large qualification calls", {
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = withr::local_tempdir(),
        entry_fixture_test = TRUE, qualification = list(programme = "perf107-phase-a",
            variant = "native", entry = "resident", member_sha256 = strrep("a", 64), deadline = 900))
    testthat::local_mocked_bindings(checkWorkflowApplication = function(...) invisible(TRUE),
        .package = "MultiScholaR")
    native <- WorkflowState$new(audit_enabled = FALSE)
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L
    recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L
    recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    native$saveState("sample_filtered", object, list(), "fixture")
    workflow <- shiny::reactiveValues(state_manager = native)
    shiny::testServer(function(input, output, session) {}, {
        session$userData$workflow_execution <- application
        application$config$entry_fixture_test <- FALSE
        expect_error(runWorkflowGUIEntryFixture("capture_sample_filtered", workflow, session,
            list(source_dir = application$config$root)), "not authorized")
        application$config$entry_fixture_test <- TRUE
        application$active <- "reserved"
        expect_error(runWorkflowGUIEntryFixture("capture_sample_filtered", workflow, session,
            list(source_dir = application$config$root)), "overlap")
        application$active <- NULL
        expect_error(runWorkflowGUIEntryFixture("../unbound", workflow, session,
            list(source_dir = application$config$root)), "not authorized")
        expect_error(runWorkflowGUIEntryFixture("capture_imputed", workflow, session,
            list(source_dir = application$config$root)), "state differs")
        application$config$qualification$deadline <- 1800
        expect_error(runWorkflowGUIEntryFixture("capture_sample_filtered", workflow, session,
            list(source_dir = application$config$root)), "not authorized")
        application$config$qualification$deadline <- 900
        large <- object
        large@peptide_data <- large@peptide_data[rep(seq_len(nrow(large@peptide_data)),
            length.out = 2001L), , drop = FALSE]
        native$saveState("sample_filtered", large, list(), "oversized")
        expect_error(runWorkflowGUIEntryFixture("capture_sample_filtered", workflow, session,
            list(source_dir = application$config$root)), "small verification boundary")
    })
})

test_that("native inspection follows completion and reference setup preserves native history", {
    root <- withr::local_tempdir()
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = root, entry_fixture_test = TRUE,
        qualification = list(programme = "perf107-phase-a", variant = "native",
            entry = "resident", member_sha256 = strrep("b", 64), deadline = 900))
    testthat::local_mocked_bindings(checkWorkflowApplication = function(...) invisible(TRUE),
        .package = "MultiScholaR")
    native <- WorkflowState$new(audit_enabled = FALSE)
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L
    recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L
    recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    native$saveState("intensity_filtered", object, list(), "input")
    native$saveState("sample_filtered", object, list(), "accepted")
    workflow <- shiny::reactiveValues(state_manager = native)
    original <- digest::digest(native$states)
    shiny::testServer(function(input, output, session) {}, {
        session$userData$workflow_execution <- application
        captured <- runWorkflowGUIEntryFixture("capture_sample_filtered", workflow, session,
            list(source_dir = root))
        expect_identical(workflow$state_manager, native)
        expect_identical(digest::digest(native$states), original)
        expect_false(captured$measured_action_active)
        expect_true(captured$verification_only)
        expect_lt(as.numeric(object.size(captured)), 65536)
        reopened <- newDIAStateAdapter(newWorkflowEntityStore(captured$root,
            "DIA-GUI-pilot", "DIA", max_cache_bytes = 0), reopen = TRUE)
        expect_identical(reopened$getState("sample_filtered")@peptide_data, object@peptide_data)
        closeDIAStateAdapter(reopened)
        native$current_state <- "intensity_filtered"
        application$config$qualification$entry <- "reference"
        prepared <- runWorkflowGUIEntryFixture("prepare_reference", workflow, session,
            list(source_dir = root))
        expect_true(usesDIAReferenceStates(workflow$state_manager))
        expect_identical(workflow$state_manager$current_state, "intensity_filtered")
        expect_identical(digest::digest(native$states), original)
        expect_identical(prepared$current_state, "intensity_filtered")
        expect_error(runWorkflowGUIEntryFixture("prepare_reference", workflow, session,
            list(source_dir = root)), "already prepared")
    })
})

test_that("entry suite recognizes only the three explicit successful case receipts", {
    code <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools/profiling/workflow_checkpoint_dia_driver.R"), code)
    for (case in c("entry-native-resident", "entry-native-reference", "entry-physical-reference")) {
        receipt <- list(schema = 1L, scenario = case, status = "science_complete")
        expect_true(code$perf107LifecycleSucceeded(case, 0L, receipt, NULL))
        expect_false(code$perf107LifecycleSucceeded(case, 1L, receipt, NULL))
        expect_false(code$perf107LifecycleSucceeded(case, 0L, receipt, list(reason = "guard")))
    }
})

test_that("the sample module runs only the hook the application carries and writes its receipt (KIT-013)", {
    root <- withr::local_tempdir()
    application <- WorkflowExecutionApplication$new()
    application$config <- list(mode = "gui", root = root)
    manager <- new.env()
    class(manager) <- c("DIAStateAdapter", "environment")
    manager$current_state <- "sample_filtered"
    manager$.context <- list(store = list())
    workflow <- shiny::reactiveValues(state_manager = manager, qc_params = list())
    calls <- list()
    shiny::testServer(mod_prot_qc_peptide_sample_server,
        args = list(workflow_data = workflow, omic_type = "proteomics",
            experiment_label = "hook", experiment_paths = list(source_dir = root)), {
        session$userData$workflow_execution <- application
        # No hook installed: the input is ignored and nothing is written.
        session$setInputs(entry_hook = "probe_one")
        session$flushReact()
        expect_false(file.exists(file.path(root, "entry-probe_one.rds")))
        application$entry_hook <- function(command, workflow_data, session, experiment_paths) {
            calls[[length(calls) + 1L]] <<- list(command = command, root = experiment_paths$source_dir)
            if (command == "probe_fail") stop("hook refused")
            list(status = "complete", command = command)
        }
        session$setInputs(entry_hook = "probe_two")
        session$flushReact()
        receipt <- readRDS(file.path(root, "entry-probe_two.rds"))
        expect_identical(receipt$status, "complete")
        expect_identical(calls[[1L]], list(command = "probe_two", root = root))
        # A receipt is written once; a repeated command does not rerun the hook.
        session$setInputs(entry_hook = "probe_two")
        session$flushReact()
        expect_length(calls, 1L)
        # A hook failure is recorded, not raised.
        session$setInputs(entry_hook = "probe_fail")
        session$flushReact()
        failed <- readRDS(file.path(root, "entry-probe_fail.rds"))
        expect_identical(failed$status, "failed")
        expect_match(failed$message, "hook refused")
        # A command that is not a plain identifier never reaches the hook or the file system.
        session$setInputs(entry_hook = "../escape")
        session$flushReact()
        expect_length(calls, 2L)
        expect_identical(list.files(root, "^entry-"), c("entry-probe_fail.rds", "entry-probe_two.rds"))
    })
})
