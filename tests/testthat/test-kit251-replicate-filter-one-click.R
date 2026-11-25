# KIT-251 (KIT-233 A2): one click on the proteomics Protein Replicate Filter Apply button is handled once. The module
# observed the button twice, through its own input and through the root session's copy of the same input, so one click
# ran the apply path twice in one flush: under a managed application the second call met the "already running" guard
# (a false warning), a refusal before the request was shown twice, and in process the second call was skipped after a
# success and ran the step again after a failure. Rows E1, E2, E3 (failing first: twice on the base commit) and guard
# G1 of the ticket (section 10).
# `shiny::testServer` registers both observers (the module session's rootScope() is the mock session, not the module's
# proxy), as the base commit's test-perf107-audit.R in-process count of 2 shows.

library(testthat)

k251Args <- function(workflow) list(workflow_data = workflow, experiment_paths = list(source_dir = tempdir(),
    protein_qc_dir = tempdir()), omic_type = "proteomics", experiment_label = "kit251")

k251Workflow <- function(state) {
    manager <- new.env()
    manager$current_state <- state
    shiny::reactiveValues(state_manager = manager, qc_params = list(), tab_status = list())
}

# The notifications shown, by text and type (shiny::showNotification replaced in its namespace).
k251Notifications <- function(env = parent.frame()) {
    shown <- new.env(parent = emptyenv())
    shown$rows <- list()
    testthat::local_mocked_bindings(
        showNotification = function(ui, ..., type = "default", id = NULL) {
            shown$rows[[length(shown$rows) + 1L]] <- list(text = paste(as.character(ui), collapse = " "), type = type)
            invisible(id)
        },
        removeNotification = function(...) invisible(NULL),
        .package = "shiny", .env = env)
    shown
}

k251Count <- function(shown, pattern) {
    sum(vapply(shown$rows, \(r) grepl(pattern, r$text, fixed = TRUE), logical(1)))
}

# A managed application on the session (KIT-015): the managed step is replaced by a recorder of its apply calls.
k251Managed <- function(calls, env = parent.frame()) {
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(...) structure(list(), class = "kit251_application"),
        installWorkflowManagedStep = function(...) list(
            apply = function(parameters, ...) {
                calls$apply[[length(calls$apply) + 1L]] <- parameters
                invisible(NULL)
            },
            cancel = function(...) invisible(NULL),
            active = function() NULL),
        withDIAStateScope = function(manager, operation) operation(),
        .package = "MultiScholaR", .env = env)
}

k251Click <- function(workflow) {
    shiny::testServer(mod_prot_qc_protein_replicate_server, args = k251Args(workflow), {
        session$setInputs(protein_grouping_variable = "group", parallel_cores = 4)
        session$flushReact()
        session$setInputs(apply_protein_replicate_filter = 1)
        session$flushReact()
    })
}

test_that("E1: under a managed application one click starts the managed step once", {
    calls <- new.env(parent = emptyenv())
    calls$apply <- list()
    k251Notifications()
    k251Managed(calls)
    testthat::local_mocked_bindings(protQcManagedApplyAllowed = function(...) TRUE, .package = "MultiScholaR")
    k251Click(k251Workflow("duplicates_removed"))
    expect_identical(length(calls$apply), 1L)
    expect_identical(calls$apply[[1L]], list(grouping_variable = "group", parallel_cores = 4))
})

test_that("E2: under a managed application one click on a normalised table shows the refusal once", {
    calls <- new.env(parent = emptyenv())
    calls$apply <- list()
    shown <- k251Notifications()
    k251Managed(calls)
    # KIT-182's gate (protQcManagedApplyAllowed, R/utils_workflow_dia_products.R) runs as it is; its check refuses.
    testthat::local_mocked_bindings(
        assertProtQcHeadNotNormalised = function(...) rlang::abort("kit251 normalised head",
            class = "prot_qc_input_normalised"),
        .package = "MultiScholaR")
    k251Click(k251Workflow("normalized"))
    expect_identical(k251Count(shown, "Error applying protein replicate filter: kit251 normalised head"), 1L)
    expect_identical(length(calls$apply), 0L)
})

test_that("E3: in process one click runs the in-process apply path once", {
    calls <- new.env(parent = emptyenv())
    calls$in_process <- 0L
    k251Notifications()
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(...) NULL,
        runProteinReplicateFilterApplyObserver = function(...) {
            calls$in_process <- calls$in_process + 1L
            invisible(NULL)
        },
        .package = "MultiScholaR")
    k251Click(k251Workflow("duplicates_removed"))
    expect_identical(calls$in_process, 1L)
})

test_that("G1: in process one click completes the step and the head is protein_replicate_filtered", {
    calls <- new.env(parent = emptyenv())
    calls$step <- 0L
    calls$completed <- 0L
    shown <- k251Notifications()
    # The in-process observer (runProteinReplicateFilterApplyObserver) runs as it is; its step commits the head as the
    # real step does on success (the computation itself is the R2 files' subject), its outputs are not drawn.
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(...) NULL,
        getDIAActiveManager = function(...) NULL,
        withDIAStateScope = function(manager, operation) operation(),
        runProteinReplicateFilterApplyStep = function(workflowData, ...) {
            calls$step <- calls$step + 1L
            workflowData$state_manager$current_state <- "protein_replicate_filtered"
            list(resultText = "kit251 filtered\n", filteredS4 = NULL)
        },
        updateProteinReplicateFilterOutputs = function(...) invisible(NULL),
        completeProtQcWorkflowStatus = function(workflowData, ...) {
            calls$completed <- calls$completed + 1L
            invisible(NULL)
        },
        .package = "MultiScholaR")
    workflow <- k251Workflow("duplicates_removed")
    k251Click(workflow)
    expect_identical(calls$step, 1L)
    expect_identical(shiny::isolate(workflow$state_manager)$current_state, "protein_replicate_filtered")
    expect_gte(calls$completed, 1L)
    expect_identical(k251Count(shown, "Protein replicate filter applied successfully"), 1L)
})
