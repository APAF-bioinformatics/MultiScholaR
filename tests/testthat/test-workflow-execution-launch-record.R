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
    expect_false(record$entry_fixture_test)
    expect_error(workflowGuiLaunchRecord(list(schema = 1L, programme = "kit-journey-v1", deadline = 14401), "x"),
        "KIT journey record is invalid")
    expect_error(workflowGuiLaunchRecord(list(schema = 1L, programme = "kit-journey-v1", deadline = 900,
        entry_fixture_test = TRUE), "x"), "KIT journey record is invalid")
    expect_error(workflowGuiLaunchRecord(list(programme = "kit-journey-v1", deadline = 900), "x"),
        "GUI qualification member is invalid")
})

test_that("the PERF-107 qualification record is interpreted as before", {
    record <- workflowGuiLaunchRecord(list(schema = 1L, programme = "perf107-phase-a", variant = "physical",
        entry = "resident", deadline = 900, entry_fixture_test = TRUE, guard_stop_test = FALSE), strrep("b", 64))
    expect_null(record$journey)
    expect_identical(record$qualification, list(programme = "perf107-phase-a", variant = "physical",
        entry = "resident", member_sha256 = strrep("b", 64), deadline = 900))
    expect_true(record$entry_fixture_test)
    expect_false(record$allow_guard_injection)
    config <- c(list(mode = "gui"), record)
    expect_identical(workflowQualificationVariant(config, "minimum_peptides_per_sample"), "physical")
    expect_identical(workflowQualificationVariant(list(mode = "gui", journey = list(native_sessions = TRUE)),
        "minimum_peptides_per_sample"), "legacy")
})

test_that("the launch record's fixture flag installs the entry hook on the application (KIT-013)", {
    expect_identical(workflowEntryHook(list(mode = "gui", entry_fixture_test = TRUE)), runWorkflowGUIEntryFixture)
    expect_null(workflowEntryHook(list(mode = "gui", entry_fixture_test = FALSE)))
    expect_null(workflowEntryHook(list(mode = "gui")))
    expect_null(workflowEntryHook(NULL))
    record <- workflowGuiLaunchRecord(list(schema = 1L, programme = "kit-journey-v1", deadline = 900), strrep("c", 64))
    expect_null(workflowEntryHook(c(list(mode = "gui"), record)))
    expect_null(WorkflowExecutionApplication$new()$entry_hook)
})
