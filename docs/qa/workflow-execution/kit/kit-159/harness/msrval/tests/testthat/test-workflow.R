# The default-user driver under shiny::testServer (KIT-159 7.3: valRunWorkflow; step outcomes A21,
# server updates A01, modals A18, re-rendered uiOutputs A33, the ledger, forced settings N9), on the
# toy module of helper-driver.R. The routes themselves are T-TAB's (test-tab.R).

toySteps <- function() {
    list(
        valStep("load", tab = NULL, declared = list(head = "loaded"), press = c("load", "confirm"),
                run = function(api) {
                    api$defaults(toyUi, NULL)
                    api$press("load")
                    api$modal_confirm("Confirm Selection")
                }),
        valStep("run", tab = NULL, declared = list(head = "normalised"), press = "run",
                run = function(api) api$press("run")),
        valStep("diagnostic", tab = NULL, declared = "diagnostic", press = "refuse",
                run = function(api) api$press("refuse")),
        valStep("widen", tab = NULL, declared = "diagnostic", structural = "^mode$",
                run = function(api) api$set(list(mode = "wide"))),
        valStep("da", tab = NULL, declared = list(tab = c(da = "complete")), press = "finish",
                run = function(api) api$press("finish"))
    )
}

runToy <- function(steps, ...) {
    valRunWorkflow(toyScript(steps), files = list(), design = NULL, spec = toySpec(),
                   observe = toyObserve, project = FALSE, ...)
}

test_that("a toy cell runs to its end: outcomes per step, the diagnostic refusal recorded", {
    res <- runToy(toySteps())
    expect_identical(res$exit, "completed")
    expect_identical(res$steps$step, c("load", "run", "diagnostic", "widen", "da"))
    expect_identical(res$steps$outcome, c("ok", "ok", "refused", "ok", "ok"))
    diag <- res$steps[res$steps$step == "diagnostic", ]
    expect_identical(diag$class, "unclassed")
    expect_match(diag$message_sha256, "^[0-9a-f]{64}$")
})

test_that("the ledger: rendered defaults, a choices-only server update, modal inputs, a re-render", {
    res <- runToy(toySteps())
    led <- res$ledger
    expect_named(led, c("step", "action", "input_id", "value", "source", "location"),
                 ignore.order = FALSE)
    expect_true(all(led$source %in% c("ui_default", "server_update", "structural", "forced",
                                      "override")))
    at <- function(step, id) led[led$step == step & led$input_id == id, , drop = FALSE]
    # The UI's first option is median; the server's choices-only update selects its first choice.
    expect_identical(at("load", "method")$value[[1L]], "median")
    expect_identical(at("load", "method")$source[[1L]], "ui_default")
    run_method <- at("run", "method")
    expect_identical(run_method$value[[1L]], "cyclicloess")
    expect_identical(run_method$source[[1L]], "server_update")
    # The modal's inputs hold their rendered values while it is open, and go with removeModal.
    modal <- led[led$step == "load" & led$action == "confirm" & led$input_id == "filter_org", ]
    expect_identical(modal$location[[1L]], "modal")
    expect_identical(modal$value[[1L]], "FALSE")
    expect_identical(nrow(at("run", "filter_org")), 0L)
    # The uiOutput's input appears at its rendered value; after `mode` changes it re-renders and the
    # input is reset to the new rendered value (A33).
    expect_identical(at("run", "k")$value[[1L]], "3")
    expect_identical(at("run", "k")$location[[1L]], "uiOutput:dynamic")
    expect_identical(at("da", "k")$value[[1L]], "9")
    expect_identical(at("da", "k")$source[[1L]], "ui_default")
    expect_identical(at("da", "mode")$source[[1L]], "structural")
    expect_match(res$ledger_sha256, "^[0-9a-f]{64}$")
    expect_true(nrow(res$actions) >= 6L)
    expect_named(res$actions, c("step", "kind", "input_id", "value"))
})

test_that("a refused step ends the cell; later steps do not run", {
    steps <- toySteps()
    steps[[1L]] <- valStep("run_early", tab = NULL, declared = list(head = "normalised"),
                           press = "run", run = function(api) {
                               api$defaults(toyUi, NULL)
                               api$press("run")
                           })
    res <- runToy(steps[1:2])
    expect_identical(res$exit, "refused")
    expect_identical(res$steps$step, "run_early")
    expect_identical(res$steps$class, "toy_method_refused")
})

test_that("an R error is did_not_complete (error); a wrong state is (unexpected_state)", {
    boom <- list(valStep("boom", tab = NULL, declared = list(head = "loaded"), press = "boom",
                         run = function(api) api$press("boom")))
    res <- runToy(boom)
    expect_identical(res$exit, "did_not_complete")
    expect_identical(res$reason, "error")
    wrong <- list(valStep("wrong", tab = NULL, declared = list(head = "elsewhere"),
                          press = "finish", run = function(api) api$press("finish")))
    res <- runToy(wrong)
    expect_identical(res$exit, "did_not_complete")
    expect_identical(res$steps$reason, "unexpected_state")
})

test_that("a script setting a non-structural input is refused; an enrichment press is refused", {
    bad <- list(valStep("bad", tab = NULL, declared = "diagnostic",
                        run = function(api) {
                            api$defaults(toyUi, NULL)
                            api$set(list(method = "median"))
                        }))
    expect_error(runToy(bad), class = "val_spec_invalid")
    enrich <- list(valStep("enrich", tab = NULL, declared = "diagnostic",
                           press = "run_enrichment_analysis", run = function(api) NULL))
    expect_error(runToy(enrich), class = "val_enrichment_test_mode")
})

test_that("forced options and inputs are scoped to the cell and recorded (KIT-163 N9)", {
    withr::local_options(list(msrval.toy_forced = NULL))
    res <- runToy(toySteps(), forced = list(msrval.toy_forced = TRUE),
                  forced_inputs = list(label = "forced-label"))
    expect_identical(res$options_at_start[["msrval.toy_forced"]], "unset")
    expect_true("MultiScholaR.ruv_design_panel" %in% names(res$options_at_start))
    expect_null(getOption("msrval.toy_forced"))
    lab <- res$ledger[res$ledger$step == "run" & res$ledger$input_id == "label", ]
    expect_identical(lab$value[[1L]], "forced-label")
    expect_identical(lab$source[[1L]], "forced")
    withr::local_options(list(msrval.toy_forced = FALSE))
    expect_error(runToy(toySteps(), forced = list(msrval.toy_forced = TRUE)),
                 class = "val_spec_invalid")
})

test_that("an arm's overrides are set with source override; one not rendered is refused", {
    arm <- list(id = "toy_arm", kind = "workflow", options = list(),
                overrides = list(mode = "wide"))
    steps <- list(valStep("set", tab = NULL, declared = "diagnostic", run = function(api) {
        api$defaults(toyUi, NULL)
        api$overrides()
        api$press("run")
    }))
    res <- runToy(steps, arm = arm)
    mode <- res$ledger[res$ledger$input_id == "mode", ]
    expect_true(any(mode$source == "override"))
    absent <- list(id = "toy_arm", kind = "workflow", options = list(),
                   overrides = list(not_rendered = 1))
    expect_error(runToy(steps, arm = absent), class = "val_spec_invalid")
})

test_that("the hook-side callback runs inside the session after the steps, with the cell's context", {
    seen <- new.env()
    res <- runToy(toySteps(), on_complete = function(ctx) {
        seen$head <- ctx$observe()$head
        seen$has_inputs <- is.environment(ctx$inputs)
        "returned"
    })
    expect_identical(seen$head, "normalised")
    expect_true(seen$has_inputs)
    expect_identical(res$value, "returned")
})

test_that("a not_permitted ledger masks the listed inputs' identifiers by SHA-256 (5.2)", {
    led <- data.frame(step = "norm", action = "run", input_id = c("ruv_control_list", "k"),
                      value = c("P1\nP2", "3"), source = c("override", "ui_default"),
                      location = "ui", stringsAsFactors = FALSE)
    masked <- valMaskLedger(led, "ruv_control_list")
    expect_identical(masked$value[[2L]], "3")
    lines <- strsplit(masked$value[[1L]], "\n", fixed = TRUE)[[1L]]
    expect_identical(lines, vapply(c("P1", "P2"), function(x) {
        digest::digest(x, algo = "sha256", serialize = FALSE)
    }, character(1), USE.NAMES = FALSE))
})
