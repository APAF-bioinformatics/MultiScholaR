# The default-user driver's parts (KIT-159 7.3; section 10.3: T-DRV-1 to T-DRV-9) and the matrix
# stage rule (7.3, revision 1.8 item 4).

test_that("T-DRV-1: the defaults reader reads every tag type as a browser sends it", {
    ui <- htmltools::tagList(
        shiny::textInput("t", NULL, "abc"),
        shiny::numericInput("n", NULL, 5),
        shiny::checkboxInput("c", NULL, TRUE),
        shiny::checkboxInput("c0", NULL, FALSE),
        shiny::radioButtons("r", NULL, c("a", "b"), selected = "b"),
        shiny::sliderInput("s", NULL, 0, 10, 3),
        shiny::selectInput("sel", NULL, c("x", "y")),
        shiny::selectInput("m", NULL, c("x", "y"), multiple = TRUE),
        shiny::selectInput("m2", NULL, c("x", "y", "z"), selected = c("x", "z"), multiple = TRUE),
        shiny::textAreaInput("ta", NULL, "line1\nline2"),
        shiny::checkboxGroupInput("g", NULL, c("p", "q"), selected = "q"),
        shiny::checkboxGroupInput("g0", NULL, c("p", "q"))
    )
    d <- valUiDefaults(as.character(ui))
    expect_identical(d[["t"]], "abc")
    expect_identical(d[["n"]], 5)
    expect_true(d[["c"]])
    expect_false(d[["c0"]])
    expect_identical(d[["r"]], "b")
    expect_identical(d[["s"]], 3)
    expect_identical(d[["sel"]], "x")
    expect_true("m" %in% names(d))
    expect_null(d[["m"]])
    expect_identical(d[["m2"]], c("x", "z"))
    expect_identical(d[["ta"]], "line1\nline2")
    expect_identical(d[["g"]], "q")
    expect_true("g0" %in% names(d))
    expect_null(d[["g0"]])
    expect_false(any(c("r1", "") %in% names(d)))
})

test_that("T-DRV-2 and T-DRV-3: a choices-only update selects the first choice; a server value is not overwritten", {
    st <- valInputState()
    valApplyRendered(st, list(norm = "median"))
    valApplyServerUpdate(st, "norm", choices = c("cyclicloess", "median"))
    expect_identical(valInputValue(st, "norm"), "cyclicloess")
    expect_identical(valInputSource(st, "norm"), "server_update")
    valApplyServerUpdate(st, "ruv_mode", selected = "covariates")
    valApplyRendered(st, list(ruv_mode = "automatic"))
    expect_identical(valInputValue(st, "ruv_mode"), "covariates")
    expect_error(valSetInput(st, "ruv_mode", "automatic", source = "ui_default"),
                 class = "val_spec_invalid")
})

test_that("T-DRV-4: a script may set only structural inputs, forced settings and overrides", {
    st <- valInputState()
    valApplyRendered(st, list(file = NULL, q_threshold = 0.05))
    valScriptSet(st, "file", "data.tsv", structural = "file")
    expect_identical(valInputSource(st, "file"), "structural")
    expect_error(valScriptSet(st, "q_threshold", 0.01, structural = "file"),
                 class = "val_spec_invalid")
    valScriptSet(st, "q_threshold", 0.01, structural = "file", forced = "q_threshold")
    expect_identical(valInputSource(st, "q_threshold"), "forced")
    ledger <- valLedgerSnapshot(st, step = "import", action = "import_btn")
    expect_named(ledger, c("step", "action", "input_id", "value", "source", "location"))
    expect_true(all(ledger$source %in% LEDGER_SOURCES))
})

test_that("T-DRV-5: a modal's inputs at their rendered values; Cancel is refused", {
    modal <- as.character(shiny::modalDialog(
        shiny::selectInput("selected_organism", NULL, c("9606", "559292"), selected = "559292"),
        shiny::checkboxInput("filter_to_organism", NULL, FALSE),
        footer = htmltools::tagList(shiny::modalButton("Cancel"),
                                    shiny::actionButton("confirm", "Confirm Selection",
                                                        class = "btn-primary"))
    ))
    st <- valInputState()
    valApplyModal(st, modal)
    expect_identical(valInputValue(st, "selected_organism"), "559292")
    led <- valLedgerSnapshot(st, "import", "confirm")
    expect_identical(led$location[led$input_id == "selected_organism"], "modal")
    expect_identical(valModalPress(modal, "confirm"), "confirm")
    expect_error(valModalPress(modal, "Cancel"), class = "val_spec_invalid")
})

test_that("T-DRV-6: a re-rendered uiOutput resets its inputs to the new rendered values", {
    st <- valInputState()
    valApplyRendered(st, list(k_rule = "auto"), location = "uiOutput:ruv_panel")
    valScriptSet(st, "k_rule", "user", structural = character(), override = "k_rule")
    valApplyRendered(st, list(k_rule = "auto"), location = "uiOutput:ruv_panel", rerender = TRUE)
    expect_identical(valInputValue(st, "k_rule"), "auto")
    expect_identical(valInputSource(st, "k_rule"), "ui_default")
})

test_that("T-DRV-7 and T-DRV-8: step outcomes; a refused diagnostic step lets the cell continue", {
    o1 <- valStepOutcome(declared = "normalised", before = "qc", after = "normalised",
                         history_changed = TRUE, error_notification = NULL)
    expect_identical(o1$outcome, "ok")
    o2 <- valStepOutcome(declared = "normalised", before = "qc", after = "qc",
                         history_changed = FALSE, error_notification = "too few features")
    expect_identical(o2$outcome, "refused")
    expect_identical(o2$class, "unclassed")
    expect_match(o2$message_sha256, "^[0-9a-f]{64}$")
    o3 <- valStepOutcome(declared = "normalised", before = "qc", after = "qc",
                         history_changed = FALSE, error_notification = NULL)
    expect_identical(o3$outcome, "did_not_complete")
    expect_identical(o3$reason, "unexpected_state")
    ran <- character()
    steps <- list(
        list(id = "s1", declared = "a", run = function() list(after = "a", history_changed = TRUE)),
        list(id = "s2", declared = "diagnostic", run = function() {
            list(after = "a", history_changed = FALSE, error_notification = "itsd failed")
        }),
        list(id = "s3", declared = "b", run = function() list(after = "b", history_changed = TRUE))
    )
    res <- valRunSteps(steps, start = "start", on_step = function(id) ran <<- c(ran, id))
    expect_identical(ran, c("s1", "s2", "s3"))
    expect_identical(res$outcome, c("ok", "refused", "ok"))
    expect_identical(attr(res, "exit"), "completed")
    stop_steps <- steps
    stop_steps[[2]]$declared <- "b"
    res2 <- valRunSteps(stop_steps, start = "start")
    expect_identical(attr(res2, "exit"), "refused")
    expect_identical(nrow(res2), 2L)
})

test_that("T-DRV-9: forced options are scoped to the cell, recorded, and refused when already set", {
    withr::local_options(list(msrval.test.forced = NULL, MultiScholaR.ruv_design_panel = NULL))
    seen <- valWithForcedOptions(list(msrval.test.forced = TRUE), getOption("msrval.test.forced"))
    expect_true(seen$value)
    expect_identical(seen$options_at_start$msrval.test.forced, "unset")
    expect_identical(seen$options_at_start$MultiScholaR.ruv_design_panel, "unset")
    expect_null(getOption("msrval.test.forced"))
    withr::local_options(list(msrval.test.forced = FALSE))
    expect_error(valWithForcedOptions(list(msrval.test.forced = TRUE), 1), class = "val_spec_invalid")
})

test_that("enrichment buttons are refused in an action script", {
    expect_error(valCheckActionScript(list(list(tab = "enrich", press = "run_enrichment"))),
                 class = "val_enrichment_test_mode")
    expect_true(valCheckActionScript(list(list(tab = "da", press = "run_da"))))
})

test_that("the matrix stage rule: factors in the model test the uncorrected matrix; display is labelled", {
    lip <- valMatrixStages("lipidomics", ruv = "factors_in_model")
    expect_identical(lip$matrix, "uncorrected")
    expect_identical(lip$display, "corrected")
    expect_identical(lip$display_label, "display; not what is tested")
    leg <- valMatrixStages("metabolomics", ruv = "correction")
    expect_identical(leg$matrix, "corrected")
    expect_null(leg$display)
    trx <- valMatrixStages("transcriptomics", ruv = "factors_in_model")
    expect_identical(trx$matrix, "uncorrected")
    expect_identical(trx$display, "corrected")
    off <- valMatrixStages("proteomics", ruv = "off")
    expect_identical(off$matrix, "normalised")
})
