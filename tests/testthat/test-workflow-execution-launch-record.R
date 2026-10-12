# GUI launch records: the PERF-107 qualification record and the KIT journey record (KIT-001).

test_that("a KIT journey record switches on native sessions with a bounded deadline", {
    record <- workflowGuiLaunchRecord(list(schema = 1L, programme = "kit-journey-v1", deadline = 5400,
        guard_stop_test = TRUE), strrep("a", 64))
    expect_null(record$qualification)
    expect_identical(record$journey$programme, "kit-journey-v1")
    expect_true(record$journey$native_sessions)
    expect_equal(record$journey$deadline, 5400)
    expect_identical(record$journey$member_sha256, strrep("a", 64))
    expect_true(record$allow_guard_injection)
    expect_error(workflowGuiLaunchRecord(list(schema = 1L, programme = "kit-journey-v1", deadline = 14401), "x"),
        "KIT journey record is invalid")
    expect_error(workflowGuiLaunchRecord(list(programme = "kit-journey-v1", deadline = 900), "x"),
        "GUI launch record is invalid")
})

test_that("the PERF-107 qualification record is retired (KIT-008): a GUI launch carries no qualification", {
    expect_error(workflowGuiLaunchRecord(list(schema = 1L, programme = "perf107-phase-a", variant = "physical",
        entry = "resident", deadline = 900, entry_fixture_test = TRUE, guard_stop_test = FALSE), strrep("b", 64)),
        class = "workflow_execution_unsupported")
    # Without a qualification the registry's production variant applies, for the GUI and for the adapters vehicle.
    expect_identical(workflowQualificationVariant(list(mode = "gui"), "minimum_peptides_per_sample"), "physical")
    expect_identical(workflowQualificationVariant(list(mode = "gui", journey = list(native_sessions = TRUE)),
        "minimum_peptides_per_sample"), "physical")
    expect_identical(workflowQualificationVariant(list(mode = "adapters"), "precursor_rollup"), "p14")
})

test_that("the application carries an optional entry hook and installs none by itself (KIT-008)", {
    expect_null(WorkflowExecutionApplication$new()$entry_hook)
    expect_null(workflowSessionApplication(NULL))
    expect_null(workflowSessionApplication(list(userData = list())))
    application <- WorkflowExecutionApplication$new()
    expect_identical(workflowSessionApplication(list(userData = list(workflow_execution = application))), application)
    application$closed <- TRUE
    expect_null(workflowSessionApplication(list(userData = list(workflow_execution = application))))
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
