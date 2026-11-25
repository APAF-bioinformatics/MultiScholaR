# KIT-261: on a managed session the normalisation tab's "Run Normalization Pipeline" button follows the omic's runner:
# disabled, with the label "Please wait: a step is still running..." (a DRAFT for the owner's wording pass, VAL-41),
# while any step of the omic runs; enabled with its own label when none does. KIT-194's refusal stays the server-side
# guard. Lipidomics, metabolomics and transcriptomics (the frame's omics).
#
# The harness (ticket section 10.1): the real runner (`installWorkflowOmicsActions`), the real step state machine and the
# real tab server (`mod_omic_norm_server`) with its real entry observer, on the KIT-078 chain mocks
# (`omicsManagedChainMocks`: the worker's entry point in this process). The poll is mocked to answer "running" while the
# test holds the last reserved operation, "failed" when the test fails it, else "complete"; `session$elapse(100)` fires it.
# `shinyjs::toggleState`, `shiny::updateActionButton` and `shiny::showNotification` are recorded.
#
# Predictions, written before R1 (the tests-only commit = the base code; ticket section 10.2):
#   L1, L2, L3, G1 (lipid, metab), G2/G3 (lipid, metab), G4, G5: pass on base and head.
#   E1/E2 (lipid, metab): 4 failed expectations each on base (no toggleState or updateActionButton call).
#   E3: 3 failed on base. E5: 2 failed on base. E6 (trx): 4 failed on base. E7: 3 failed on base. In all 20 failed
#   expectations on base, each a toggleState or label expectation of an E row, and no error.
#   Every precondition expectation ("k261 precondition") passes on both; one that fails is stop rule 1.
#   Head (the fix): every expectation passes.

k261Busy <- "Please wait: a step is still running..."
k261Default <- "Run Normalization Pipeline"
k261Refusal <- "A step is still running. Press \"Run Normalization Pipeline\" again when it has finished."

k261Omics <- list(
    lipidomics = list(prefix = "lipid", owner = function(root) lipidDeclarationsOwner(root),
        intensity = function(...) runLipidReferenceIntensityFilter(...),
        duplicates = function(...) runLipidReferenceDuplicateResolution(...)),
    metabolomics = list(prefix = "metab", owner = function(root) metabDeclarationsOwner(root),
        intensity = function(...) runMetabReferenceIntensityFilter(...),
        duplicates = function(...) runMetabReferenceDuplicateResolution(...)))

#' A declaration fixture after quality control (intensity filter, duplicates, Finalize QC), in process: the head is
#' `<prefix>_qc_complete`, where the user opens the normalisation tab.
k261Fixture <- function(omic, root) {
    spec <- k261Omics[[omic]]
    fixture <- spec$owner(file.path(root, spec$prefix))
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    suppressWarnings(suppressMessages(utils::capture.output({
        spec$intensity(owner, 1, 0.5, config, hydrate = FALSE)
        spec$duplicates(owner, config, hydrate = FALSE)
        commitWorkflowOmicsStoreStep(owner, paste0(spec$prefix, "_qc_finalize"), paste0(spec$prefix, "_qc_complete"), config,
            "QC processing complete - ready for normalization")
    })))
    fixture
}

k261QuietProgress <- function(message) {
    list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL), close = function(...) invisible(NULL))
}

#' The tab under test: the omic's runner when managed (installed before the tab server, as R/mod_omic_workflow.R does),
#' then the normalisation frame with `selected_tab` the test's reactive value.
k261Server <- function(fixture, omic, tab, managed) {
    force(fixture)
    force(omic)
    force(tab)
    force(managed)
    function(input, output, session) {
        if (managed) {
            runner <- installWorkflowOmicsActions(session, fixture$workflow, omic, "kit261", fixture$paths)
            runner$progress <- k261QuietProgress
        }
        session$userData$k261_norm <- mod_omic_norm_server("norm", fixture$workflow, fixture$paths, omic, "kit261",
            selected_tab = function() tab(), config = omicGuiConfig(omic)$norm)
    }
}

#' Record the Run button's state and label calls, and the notifications, in `rec`.
k261Capture <- function(rec, env) {
    rec$toggles <- logical()
    rec$labels <- character()
    rec$other <- character()
    rec$notes <- character()
    testthat::local_mocked_bindings(toggleState = function(id = NULL, condition = NULL, selector = NULL, asis = FALSE) {
        if (identical(id, "run_normalization")) rec$toggles <- c(rec$toggles, isTRUE(condition)) else rec$other <- c(rec$other, id)
        invisible(NULL)
    }, .package = "shinyjs", .env = env)
    testthat::local_mocked_bindings(
        updateActionButton = function(session = NULL, inputId, label = NULL, icon = NULL, ...) {
            if (identical(inputId, "run_normalization")) rec$labels <- c(rec$labels, label %||% "<no label>") else rec$other <- c(rec$other, inputId)
            invisible(NULL)
        },
        showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
            type = c("default", "message", "warning", "error"), session = NULL) {
            rec$notes <- c(rec$notes, paste0("[", type[[1L]], "] ", paste(as.character(ui), collapse = " ")))
            invisible(id %||% "note")
        },
        .package = "shiny", .env = env)
    invisible(rec)
}

#' The poll of the chain mocks, held: "running" while the last reserved operation is in `hold$ops`, "failed" while it is
#' in `hold$fail`, else "complete".
k261HeldPoll <- function(probe, hold, env) {
    testthat::local_mocked_bindings(pollWorkflowAction = function(owner) {
        operation <- utils::tail(probe$operations, 1L)
        list(status = if (isTRUE(operation %in% hold$fail)) "failed" else if (isTRUE(operation %in% hold$ops)) "running" else "complete",
            error = NULL)
    }, .package = "MultiScholaR", .env = env)
}

#' One tab session of `omic` (managed or in process) on a fresh fixture; `drive(ctx)` runs inside testServer, where
#' `ctx` holds the session, the fixture, the probe, the hold, the records, the tab value and the omic.
k261Run <- function(omic, drive, managed = TRUE) {
    root <- withr::local_tempdir()
    ctx <- new.env(parent = emptyenv())
    ctx$omic <- omic
    ctx$prefix <- k261Omics[[omic]]$prefix
    ctx$fixture <- k261Fixture(omic, root)
    withr::defer(closeDIAStateAdapter(ctx$fixture$owner))
    ctx$probe <- new.env()
    ctx$hold <- new.env()
    ctx$hold$ops <- character()
    ctx$hold$fail <- character()
    ctx$rec <- new.env()
    here <- environment()
    if (managed) {
        omicsManagedChainMocks(ctx$fixture, ctx$probe, env = here)
        k261HeldPoll(ctx$probe, ctx$hold, here)
    }
    k261Capture(ctx$rec, here)
    ctx$tab <- shiny::reactiveVal("qc")
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(k261Server(ctx$fixture, omic, ctx$tab, managed), {
        ctx$session <- session
        drive(ctx)
    }))))
    invisible(ctx)
}

k261Runner <- function(ctx) ctx$session$userData$workflow_omics_actions[[ctx$omic]]
k261Norm <- function(ctx) ctx$session$userData$k261_norm
k261Active <- function(ctx) shiny::isolate(k261Runner(ctx)$step$active())

k261Tick <- function(ctx) {
    ctx$session$elapse(100)
    ctx$session$flushReact()
}

k261Enter <- function(ctx, value = "norm") {
    ctx$tab(value)
    ctx$session$flushReact()
}

#' The precondition of G2 and E1 to E7 (ticket 10.1): the last reserved operation is the one held, and the step is active.
k261Held <- function(ctx, operation, count = NULL) {
    expect_identical(utils::tail(ctx$probe$operations, 1L), operation, info = "k261 precondition: the last reserved operation")
    if (!is.null(count)) expect_identical(length(ctx$probe$operations), count, info = "k261 precondition: one new reservation")
    expect_false(is.null(k261Active(ctx)), info = "k261 precondition: the runner's step is active")
}

#' The lane's run inputs (ITSD off, no between-sample method, RUV skipped), the others from `omicsTabInputs()`.
k261SetRunInputs <- function(ctx) {
    inputs <- utils::modifyList(omicsTabInputs(), list(apply_itsd = FALSE, norm_method = "none", ruv_mode = "skip"))
    do.call(ctx$session$setInputs, stats::setNames(inputs, paste0("norm-", names(inputs))))
    ctx$session$flushReact()
}

k261Press <- function(ctx, value) {
    ctx$session$setInputs(`norm-run_normalization` = value)
    ctx$session$flushReact()
}

k261Log <- function(ctx) sub("^\\[[0-9:]+\\] ", "", shiny::isolate(k261Norm(ctx)$normalization_log))

k261Snapshot <- function(ctx) {
    owner <- ctx$fixture$owner
    list(head = owner$current_state, history = as.character(unlist(owner$getHistory())),
        complete = shiny::isolate(k261Norm(ctx)$normalization_complete))
}

k261Since <- function(x, n) utils::tail(x, length(x) - n)

k261Errors <- function(notes) sub("^\\[error\\] ", "", grep("^\\[error\\] ", notes, value = TRUE))

# ------------------------------------------------------------------------------------------------
# L rows: the lifecycle of the instrument.
# ------------------------------------------------------------------------------------------------

#' The runner alone (no tab), managed, on the lipid fixture.
k261RunnerOnly <- function(drive) {
    root <- withr::local_tempdir()
    ctx <- new.env(parent = emptyenv())
    ctx$omic <- "lipidomics"
    ctx$fixture <- k261Fixture("lipidomics", root)
    withr::defer(closeDIAStateAdapter(ctx$fixture$owner))
    ctx$probe <- new.env()
    ctx$hold <- new.env()
    ctx$hold$ops <- character()
    ctx$hold$fail <- character()
    ctx$rec <- new.env()
    here <- environment()
    omicsManagedChainMocks(ctx$fixture, ctx$probe, env = here)
    k261HeldPoll(ctx$probe, ctx$hold, here)
    k261Capture(ctx$rec, here)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        installWorkflowOmicsActions(session, ctx$fixture$workflow, "lipidomics", "kit261", ctx$fixture$paths)
    }, {
        ctx$session <- session
        drive(ctx)
    }))))
    invisible(ctx)
}

test_that("KIT-261 L1: a held runner action is active until it is released, then delivered", {
    skip_if_not(workflowOwnerCapable())
    ctx <- k261RunnerOnly(function(ctx) {
        ctx$hold$ops <- "lipid_filtering_metrics"
        ctx$delivered <- FALSE
        k261Runner(ctx)$run(list(omicsManagedProduct("lipid_filtering_metrics", list(), function(result) ctx$delivered <- TRUE)))
        ctx$session$flushReact()
        k261Tick(ctx)
        ctx$held_operation <- utils::tail(ctx$probe$operations, 1L)
        ctx$held_active <- !is.null(k261Active(ctx))
        ctx$held_delivered <- ctx$delivered
        ctx$hold$ops <- character()
        k261Tick(ctx)
        ctx$released_active <- !is.null(k261Active(ctx))
    })
    expect_identical(ctx$held_operation, "lipid_filtering_metrics")
    expect_true(ctx$held_active)
    expect_false(ctx$held_delivered)
    expect_false(ctx$released_active)
    expect_true(ctx$delivered)
})

test_that("KIT-261 L2: a held runner action whose poll answers 'failed' is told, and its request's failed is called", {
    skip_if_not(workflowOwnerCapable())
    ctx <- k261RunnerOnly(function(ctx) {
        ctx$hold$ops <- "lipid_filtering_metrics"
        ctx$failed <- NULL
        k261Runner(ctx)$run(list(omicsManagedProduct("lipid_filtering_metrics", list(), function(result) NULL)),
            failed = function(error) ctx$failed <- conditionMessage(error))
        ctx$session$flushReact()
        ctx$held_active <- !is.null(k261Active(ctx))
        ctx$hold$fail <- "lipid_filtering_metrics"
        k261Tick(ctx)
        ctx$released_active <- !is.null(k261Active(ctx))
    })
    expect_true(ctx$held_active)
    expect_false(ctx$released_active)
    expect_false(is.null(ctx$failed))
    expect_true(any(startsWith(k261Errors(ctx$rec$notes), "Error applying lipidomics step: ")), info = paste(ctx$rec$notes, collapse = " | "))
})

test_that("KIT-261 L3: the managed harness: opening the tab dispatches the entry figures, held then released; one press completes a run", {
    skip_if_not(workflowOwnerCapable())
    ctx <- k261Run("lipidomics", function(ctx) {
        k261SetRunInputs(ctx)
        ctx$hold$ops <- "lipid_qc_figures"
        before <- length(ctx$probe$operations)
        k261Enter(ctx)
        expect_identical(k261Since(ctx$probe$operations, before), "lipid_qc_figures")
        expect_false(is.null(k261Active(ctx)))
        trigger <- shiny::isolate(k261Norm(ctx)$plot_refresh_trigger)
        ctx$hold$ops <- character()
        k261Tick(ctx)
        expect_null(k261Active(ctx))
        expect_true(shiny::isolate(k261Norm(ctx)$pre_norm_qc_generated))
        expect_identical(shiny::isolate(k261Norm(ctx)$plot_refresh_trigger), trigger + 1)
        k261Press(ctx, 1L)
        k261Tick(ctx)
        expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete), info = paste(k261Log(ctx), collapse = " | "))
        expect_identical(utils::tail(k261Log(ctx), 1L), "Normalization pipeline complete!")
    })
})

# ------------------------------------------------------------------------------------------------
# E rows: fail on base (no call for the button), pass with the fix.
# ------------------------------------------------------------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
    local({
        omic <- omic
        test_that(paste("KIT-261 E1, E2:", omic, "the button is disabled with the busy label while the entry figures run, then enabled"), {
            skip_if_not(workflowOwnerCapable())
            k261Run(omic, function(ctx) {
                qc_figures <- paste0(ctx$prefix, "_qc_figures")
                ctx$hold$ops <- qc_figures
                k261Enter(ctx)
                k261Held(ctx, qc_figures)
                expect_identical(utils::tail(ctx$rec$toggles, 1L), FALSE, info = "E1 toggleState")
                expect_identical(utils::tail(ctx$rec$labels, 1L), k261Busy, info = "E1 label")
                ctx$hold$ops <- character()
                k261Tick(ctx)
                expect_null(k261Active(ctx), info = "E2: the entry step ended")
                expect_identical(utils::tail(ctx$rec$toggles, 1L), TRUE, info = "E2 toggleState")
                expect_identical(utils::tail(ctx$rec$labels, 1L), k261Default, info = "E2 label")
            })
        })
    })
}

test_that("KIT-261 E3: lipidomics, Finalize QC's figure then the entry figures (queued behind it): busy throughout, enabled at the end", {
    skip_if_not(workflowOwnerCapable())
    k261Run("lipidomics", function(ctx) {
        ctx$hold$ops <- c("lipid_filtering_metrics", "lipid_qc_figures")
        k261Runner(ctx)$run(list(omicsManagedProduct("lipid_filtering_metrics", list(), function(result) NULL)))
        ctx$session$flushReact()
        k261Held(ctx, "lipid_filtering_metrics")
        k261Enter(ctx)
        # The entry figures wait behind the metrics figure.
        expect_identical(utils::tail(ctx$probe$operations, 1L), "lipid_filtering_metrics")
        ctx$hold$ops <- "lipid_qc_figures"
        k261Tick(ctx)
        k261Held(ctx, "lipid_qc_figures")
        expect_identical(ctx$rec$toggles, FALSE, info = "E3: disabled once, while the metrics and then the entry figures ran")
        ctx$hold$ops <- character()
        k261Tick(ctx)
        expect_null(k261Active(ctx))
        expect_identical(ctx$rec$toggles, c(FALSE, TRUE), info = "E3 toggleState")
        expect_identical(ctx$rec$labels, c(k261Busy, k261Default), info = "E3 labels")
    })
})

test_that("KIT-261 E5: lipidomics, the entry figures fail: the button is enabled again with its own label", {
    skip_if_not(workflowOwnerCapable())
    k261Run("lipidomics", function(ctx) {
        ctx$hold$ops <- "lipid_qc_figures"
        k261Enter(ctx)
        k261Held(ctx, "lipid_qc_figures")
        ctx$hold$fail <- "lipid_qc_figures"
        k261Tick(ctx)
        ctx$hold$fail <- character()
        ctx$hold$ops <- character()
        expect_null(k261Active(ctx))
        expect_true(any(startsWith(k261Log(ctx), "Error generating Pre-QC: ")), info = paste(k261Log(ctx), collapse = " | "))
        expect_identical(ctx$rec$toggles, c(FALSE, TRUE), info = "E5 toggleState")
        expect_identical(utils::tail(ctx$rec$labels, 1L), k261Default, info = "E5 label")
    })
})

test_that("KIT-261 E6: transcriptomics, the button is disabled with the busy label while a transcriptomics step runs, then enabled", {
    skipWithoutTrxDeclReference()
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    ctx <- new.env(parent = emptyenv())
    ctx$omic <- "transcriptomics"
    ctx$probe <- new.env()
    ctx$hold <- new.env()
    ctx$hold$ops <- "trx_qc_figures"
    ctx$hold$fail <- character()
    ctx$rec <- new.env()
    here <- environment()
    omicsManagedChainMocks(fixture, ctx$probe, env = here)
    k261HeldPoll(ctx$probe, ctx$hold, here)
    k261Capture(ctx$rec, here)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(module_ci_trx_tabs_server(fixture, managed = TRUE), {
        ctx$session <- session
        ctx$session$flushReact()
        k261Runner(ctx)$run(list(omicsManagedProduct("trx_qc_figures", list(stage = "raw"), function(result) NULL)))
        ctx$session$flushReact()
        k261Held(ctx, "trx_qc_figures")
        expect_identical(utils::tail(ctx$rec$toggles, 1L), FALSE, info = "E6 toggleState")
        expect_identical(utils::tail(ctx$rec$labels, 1L), k261Busy, info = "E6 label")
        ctx$hold$ops <- character()
        k261Tick(ctx)
        expect_null(k261Active(ctx), info = "E6: the step ended")
        expect_identical(utils::tail(ctx$rec$toggles, 1L), TRUE, info = "E6 toggleState after")
        expect_identical(utils::tail(ctx$rec$labels, 1L), k261Default, info = "E6 label after")
    }))))
})

test_that("KIT-261 E7: lipidomics, during the run the button stays disabled when the tab is left and opened again", {
    skip_if_not(workflowOwnerCapable())
    k261Run("lipidomics", function(ctx) {
        k261SetRunInputs(ctx)
        k261Enter(ctx)
        expect_null(k261Active(ctx))
        ctx$hold$ops <- "lipid_qc_figures"
        count <- length(ctx$probe$operations)
        toggles <- length(ctx$rec$toggles)
        k261Press(ctx, 1L)
        k261Held(ctx, "lipid_qc_figures", count + 1L)
        k261Enter(ctx, "qc")
        k261Enter(ctx, "norm")
        k261Held(ctx, "lipid_qc_figures", count + 1L)
        expect_identical(k261Since(ctx$rec$toggles, toggles), FALSE, info = "E7: disabled once, kept disabled across the re-entry")
        expect_identical(utils::tail(ctx$rec$labels, 1L), k261Busy, info = "E7 label")
        ctx$hold$ops <- character()
        k261Tick(ctx)
        expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete))
        expect_identical(utils::tail(ctx$rec$toggles, 1L), TRUE, info = "E7: enabled after the run")
    })
})

# ------------------------------------------------------------------------------------------------
# G rows: pass on base and head.
# ------------------------------------------------------------------------------------------------

for (omic in c("lipidomics", "metabolomics")) {
    local({
        omic <- omic
        test_that(paste("KIT-261 G1:", omic, "in process (no runner): no call for the button; the run's log as before"), {
            skip_if_not(workflowOwnerCapable())
            ctx <- k261Run(omic, function(ctx) {
                k261SetRunInputs(ctx)
                k261Enter(ctx)
                lines <- length(k261Log(ctx))
                k261Press(ctx, 1L)
                log <- k261Since(k261Log(ctx), lines)
                expect_identical(utils::head(log, 2L), c("Starting normalization pipeline...", "Post-filtering state captured"))
                expect_identical(utils::tail(log, 1L), "Normalization pipeline complete!", info = paste(log, collapse = " | "))
                expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete))
            }, managed = FALSE)
            expect_identical(ctx$rec$toggles, logical())
            expect_identical(ctx$rec$labels, character())
        })

        test_that(paste("KIT-261 G2, G3:", omic, "a press while the entry figures run is refused by KIT-194's guard; after them one press runs"), {
            skip_if_not(workflowOwnerCapable())
            k261Run(omic, function(ctx) {
                qc_figures <- paste0(ctx$prefix, "_qc_figures")
                k261SetRunInputs(ctx)
                ctx$hold$ops <- qc_figures
                k261Enter(ctx)
                k261Held(ctx, qc_figures)
                operations <- ctx$probe$operations
                before <- k261Snapshot(ctx)
                notes <- length(ctx$rec$notes)
                # G2: a click already in flight (setInputs bypasses a disabled button).
                k261Press(ctx, 1L)
                pressed <- k261Since(ctx$rec$notes, notes)
                expect_identical(k261Errors(pressed), k261Refusal, info = paste(pressed, collapse = " | "))
                expect_identical(ctx$probe$operations, operations)
                expect_identical(k261Snapshot(ctx), before)
                ctx$hold$ops <- character()
                k261Tick(ctx)
                expect_null(k261Active(ctx))
                # G3: one press after the entry figures.
                notes <- length(ctx$rec$notes)
                k261Press(ctx, 2L)
                k261Tick(ctx)
                run <- k261Since(ctx$probe$operations, length(operations))
                expect_identical(run[[1L]], qc_figures, info = paste(run, collapse = " "))
                expect_false(any(k261Errors(k261Since(ctx$rec$notes, notes)) == k261Refusal))
                expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete), info = paste(k261Log(ctx), collapse = " | "))
                expect_true(length(k261Snapshot(ctx)$history) > length(before$history))
                expect_true("Normalization pipeline complete!" %in% k261Log(ctx))
            })
        })
    })
}

test_that("KIT-261 G4: lipidomics, re-entering the tab during the run dispatches no second entry figures", {
    skip_if_not(workflowOwnerCapable())
    k261Run("lipidomics", function(ctx) {
        k261SetRunInputs(ctx)
        k261Enter(ctx)
        ctx$hold$ops <- "lipid_qc_figures"
        count <- length(ctx$probe$operations)
        k261Press(ctx, 1L)
        k261Held(ctx, "lipid_qc_figures", count + 1L)
        operations <- ctx$probe$operations
        k261Enter(ctx, "qc")
        k261Enter(ctx, "norm")
        expect_identical(ctx$probe$operations, operations)
        expect_identical(sum(ctx$probe$operations == "lipid_qc_figures"), 2L)
        ctx$hold$ops <- character()
        k261Tick(ctx)
        expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete))
    })
})

test_that("KIT-261 G5: lipidomics, after the entry figures failed one press runs", {
    skip_if_not(workflowOwnerCapable())
    k261Run("lipidomics", function(ctx) {
        k261SetRunInputs(ctx)
        ctx$hold$fail <- "lipid_qc_figures"
        k261Enter(ctx)
        ctx$hold$fail <- character()
        expect_null(k261Active(ctx))
        expect_true(any(startsWith(k261Log(ctx), "Error generating Pre-QC: ")))
        notes <- length(ctx$rec$notes)
        k261Press(ctx, 1L)
        k261Tick(ctx)
        expect_false(any(k261Errors(k261Since(ctx$rec$notes, notes)) == k261Refusal))
        expect_true(shiny::isolate(k261Norm(ctx)$normalization_complete), info = paste(k261Log(ctx), collapse = " | "))
    })
})
