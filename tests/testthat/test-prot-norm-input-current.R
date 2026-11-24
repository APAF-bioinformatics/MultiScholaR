# KIT-182 (VAL-48, VAL-91, VAL-99): proteomics normalisation and the correlation filter use the current quality-control
# result, and quality control is refused on a normalised table. On the native DIA owner of helper-prot-norm-core.R (the
# "small" recipe: 64 proteins x 12 samples, committed as `protein_replicate_filtered`, PRF1 below) and its DIA_LIMPA
# owner. "QC change" is the Protein Intensity Filter through its kernel and the step's display publication (the
# scripted path, which carries no QC guard), so it builds the same states on the tests-only commit and on the fix.
# Calls to functions whose formals the fix extends pass only the arguments the function has, so the tests-only commit
# runs the same tests and fails at a value. Lines starting `KIT-182 effect` and `KIT-182 F1` are the evidence of A6 and
# F1 (docs/qa/workflow-execution/kit/kit-182/).

k182QcParameters <- list(strict_mode = FALSE, min_reps_per_group = 2, min_groups = 2,
    proteins_intensity_cutoff_percentile = 10)

k182Settings <- function(mode = "skip", ...) utils::modifyList(list(norm_method = "quantile", ruv_mode = mode,
    auto_percentage_min = NULL, auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
    adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL), list(...))

k182Aesthetics <- function() list(color_var = "group", shape_var = NULL)

#' A call with only the arguments `fn` has (the fix adds trailing formals).
k182Call <- function(fn, args) do.call(fn, args[intersect(names(args), names(formals(fn)))])

k182Scope <- function(fixture, code) withDIAStateScope(fixture$manager, code)

k182Rows <- function(fixture, reference) k182Scope(fixture, function()
    readDIAStateMetadata(fixture$manager, reference)$data$dimensions$rows)

k182Head <- function(fixture) getDIAStateReference(fixture$manager)

k182History <- function(fixture) diaStateHistory(fixture$manager)

k182Revert <- function(fixture, state) k182Scope(fixture, function() selectDIAStateReference(fixture$manager, state))

#' The tab's in-process normalisation click.
k182Click <- function(fixture, norm, mode = "skip", ..., pre_qc = NULL) {
    protNormCoreQuiet(k182Call(runDIAReferenceNormalization, list(settings = k182Settings(mode, ...),
        workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL),
        ruv_qc = function(reference) invisible(NULL), grouping = function() "group", pre_qc = pre_qc)))
}

#' The tab's managed normalisation click on a runner that records the actions it is given and runs none.
k182ManagedClick <- function(fixture, norm, mode = "skip") {
    seen <- new.env(parent = emptyenv())
    seen$actions <- list()
    runner <- list(run = function(actions, done = NULL, failed = NULL) {
        if (is.function(actions)) actions <- actions()
        seen$actions <- c(seen$actions, actions)
        invisible(TRUE)
    })
    protNormCoreQuiet(k182Call(dispatchProtNormNormalization, list(runner = runner, settings = k182Settings(mode),
        workflowData = fixture$workflow, normData = norm, experimentPaths = fixture$paths, omicType = "proteomics",
        experimentLabel = "core", postQc = function(reference) invisible(NULL), ruvQc = function(reference) invisible(NULL),
        grouping = function() "group", finish = function(completed) completed, failed = function(error) stop(error))))
    seen$actions
}

k182ActionInput <- function(actions, operation) {
    hits <- Filter(function(action) identical(action$operation, operation), actions)
    if (length(hits)) hits[[1L]]$input$entity_id else NA_character_
}

#' The correlation filter in process (Apply, or Skip with `skipped`).
k182Correlation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    protNormCoreQuiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped, "proteomics",
        "core", if (skipped) function() NULL else function() "group", k182Aesthetics))
}

#' The tab's normalisation run (`runProtNormNormalizationWorkflow`), its notices and errors recorded.
k182Tab <- function(fixture, norm, mode = "skip", generatePreNormalizationQcFn = NULL) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$errors <- list()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    protNormCoreQuiet(shiny::isolate(k182Call(runProtNormNormalizationWorkflow, list(input = k182Settings(mode),
        workflowData = fixture$workflow, normData = norm, experimentPaths = fixture$paths, omicType = "proteomics",
        experimentLabel = "core", generatePostNormalizationQcFn = function(reference) invisible(NULL),
        generateRuvCorrectedQcFn = function(reference) invisible(NULL), getRuvGroupingVariableFn = function() "group",
        showNotificationFn = note, handleErrorFn = function(error) {
            seen$errors <- c(seen$errors, list(error))
            invisible(FALSE)
        }, messageFn = function(...) NULL, generatePreNormalizationQcFn = generatePreNormalizationQcFn))))
    as.list(seen)
}

#' The QC change (the third re-check's C01): the Protein Intensity Filter on the head through its kernel, committed as
#' the step commits it, and the step's display published as the sub-tab publishes it after an Apply. Should the filter
#' remove no protein (F2), the head's object less its first protein row is committed instead. The new head.
k182QcChange <- function(fixture) {
    manager <- fixture$manager
    rlang::local_bindings(filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    before <- k182Rows(fixture, k182Head(fixture))
    protNormCoreQuiet(k182Scope(fixture, function() {
        value <- computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager), k182QcParameters)
        if (identical(as.integer(value$protein_count), as.integer(before))) {
            object <- manager$getState()
            object@protein_quant_table <- object@protein_quant_table[-1L, , drop = FALSE]
            manager$saveState("protein_intensity_filtered", object, list(), "qc change")
            return(invisible(NULL))
        }
        result <- finishDIAReferenceProteinStep(fixture$workflow, value, publish = FALSE)
        publishDIAProteinProgress(result, "11_protein_intensity_filtered", "proteomics", "core", updateProteinFiltering)
        invisible(NULL)
    }))
    expect_null(assertDIAExportReady(manager))
    head <- k182Head(fixture)
    cat(sprintf("KIT-182 F2 before=%d after=%d removed=%d\n", as.integer(before), as.integer(k182Rows(fixture, head)),
        as.integer(before - k182Rows(fixture, head))))
    head
}

#' A protein QC kernel on the head, committed as its step commits it, its display published.
k182QcKernel <- function(fixture, kernel = "intensity") {
    if (identical(kernel, "intensity")) return(k182QcChange(fixture))
    manager <- fixture$manager
    rlang::local_bindings(filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    protNormCoreQuiet(k182Scope(fixture, function() {
        value <- computeDIAReferenceProteinReplicateFilter(manager, getDIAStateReference(manager),
            list(grouping_variable = "group", parallel_cores = 1))
        result <- finishDIAReferenceProteinStep(fixture$workflow, value, publish = FALSE)
        publishDIAProteinProgress(result, "13_protein_replicate_filtered", "proteomics", "core", updateProteinFiltering)
        invisible(NULL)
    }))
    expect_null(assertDIAExportReady(manager))
    k182Head(fixture)
}

#' The DIA_LIMPA Peptide Count Filter through its kernel (as the owner's fixture commits it); no publication needed.
k182LimpaPeptideCount <- function(fixture, threshold = 2) {
    manager <- fixture$manager
    protNormCoreQuiet(k182Scope(fixture, function() {
        value <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = threshold, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = value$config[1:2])
        finishDIAReferenceProteinStep(fixture$workflow, value, publish = FALSE)
    }))
    k182Head(fixture)
}

k182PassThrough <- function(fixture, norm) protNormCoreQuiet(runDIAReferenceLimpaNormalisationPassThrough(
    k182Settings("covariates"), fixture$workflow, norm, function(reference) invisible(NULL)))

k182RuvCovariates <- function(fixture, norm, mode) protNormCoreQuiet(runDIAReferenceRuvCovariates(list(ruv_mode = mode,
    color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
    function(reference) invisible(NULL), function() "group"))

#' The texts proposed in the ticket's 7.4 (VAL-41), for the expectations.
k182RunText <- function(step) paste0("this normalisation would start from a table that has already been normalised ",
    "(quality-control steps were applied after normalisation, or the project recorded a normalised table as its ",
    "starting table), so it cannot be normalised again. In the Quality Control tab press Revert until the current ",
    "step is \"", step, "\", the last quality-control result before normalisation, apply again from there any ",
    "quality-control steps you want, then press \"Run Normalization & RUV\".")
k182QcText <- function(step) paste0("quality-control steps cannot be applied to a table that has already been ",
    "normalised. To change quality control, press Revert in this tab until the current step is \"", step, "\", the last ",
    "quality-control result before normalisation, then apply the steps again and press \"Run Normalization & RUV\" again.")
k182CorrelationText <- paste0("the table in the Quality Control tab has changed since normalisation was last run (a ",
    "step was applied or reverted), so the correlation filter would use a table that does not match it. Press \"Run ",
    "Normalization & RUV\" again, then apply or skip the correlation filter.")
k182NoticeText <- paste("This normalisation uses your current quality-control result. The previous normalisation used",
    "an earlier table, from before your latest quality-control changes.")

k182Error <- function(code) tryCatch({
    force(code)
    NULL
}, error = function(e) e)

#' The F1 predicate as the ticket's 7.1 defines it, clause by clause (the tests-only commit has no such function); on
#' the fix each case also compares the package's own predicate.
k182Clauses <- function(manager, spec, head) {
    labels <- c("normalized", "ruv_corrected", "correlation_filtered", "owner:normalized-science")
    history <- diaStateHistory(manager)
    c(i = !is.null(spec$products$normalization_science),
        ii = any(vapply(spec$dependencies, function(d) is.list(d) && isTRUE(d$identity$parameters %in% labels), NA)),
        iii = head && isTRUE(manager$.context$snapshot$step_provenance[[manager$current_state]]$adopted) &&
            any(match(c("normalized", "ruv_corrected"), history) < match(manager$current_state, history), na.rm = TRUE))
}

k182F1 <- function(fixture, case, expected, reference = NULL) {
    manager <- fixture$manager
    head <- is.null(reference)
    spec <- k182Scope(fixture, function() readDIAStateMetadata(manager, reference %||% getDIAStateReference(manager))$data)
    clauses <- k182Clauses(manager, spec, head)
    cat(sprintf("KIT-182 F1 %s value=%s expected=%s clauses=%s\n", case, any(clauses), expected,
        paste(names(clauses)[clauses], collapse = ",")))
    expect_identical(any(clauses), expected, info = case)
    ns <- asNamespace("MultiScholaR")
    if (exists("diaNormalisedLineage", envir = ns, inherits = FALSE)) {
        shipped <- k182Scope(fixture, function() get("diaNormalisedLineage", envir = ns)(manager, spec, head = head))
        expect_identical(shipped, any(clauses), info = paste("shipped", case))
    }
    clauses
}

# ------------------------------------------------------------------------------------------------
# F1: the lineage predicate's clauses on every path (feasibility gate; passes on both commits)
# ------------------------------------------------------------------------------------------------

test_that("KIT-182 F1 (a) to (f), (i): QC computed downstream of a normalisation has normalisation in its lineage; QC from PRF1 has not", {
    # (a) and (i): from PRF1, before and after a normalisation and a Revert.
    fixture <- protNormCoreOwner()
    k182QcChange(fixture)
    k182F1(fixture, "a", FALSE)
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    k182Revert(fixture, "protein_replicate_filtered")
    k182QcChange(fixture)
    k182F1(fixture, "i", FALSE)
    # (b) the filter on head `normalized`; (c) the replicate filter's `rows` variant there.
    for (kernel in c("intensity", "replicate")) {
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, "skip")
        k182QcKernel(fixture, kernel)
        k182F1(fixture, if (kernel == "intensity") "b" else "c", TRUE)
    }
    # (d) after the correlation applied, (e) skipped; (f) after the DA re-commits `correlation_filtered`.
    for (case in c("d", "e", "f")) {
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, "manual")
        k182Correlation(fixture, norm, 0.5, skipped = identical(case, "e"))
        if (identical(case, "f")) protNormCoreQuiet(k182Scope(fixture, function() commitDIADaArguments(fixture$manager,
            list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, formula_string = "~ 0 + group"), fixture$workflow,
            data.frame(contrasts = "groupB-groupA", full_format = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
            "Protein.Ids")))
        expect_identical(fixture$manager$current_state, "correlation_filtered", info = case)
        k182QcChange(fixture)
        clauses <- k182F1(fixture, case, TRUE)
        # Clause (i) alone misses everything downstream of the correlation filter.
        expect_false(clauses[["i"]], info = case)
    }
})

test_that("KIT-182 F1 (g): the QC step on the worker path carries the same lineage", {
    skip_if_not(workflowOwnerCapable())
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR")
    for (case in c("b", "d")) {
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, if (case == "b") "skip" else "manual")
        if (case == "d") k182Correlation(fixture, norm, 0.5)
        controller <- omicsWorkerController(list(owner = fixture$manager))
        # The worker computes and the controller's acceptance commits; what the helper reads afterwards is a display
        # handoff, which this case does not need (an error there is logged, the committed head is what is checked).
        # R1 instrument repair (VAL-54): a committing operation's completion names its state, as W2's dispatch does.
        error <- k182Error(protNormCoreQuiet(protNormW2Dispatch(fixture, controller, "protein_intensity_filter",
            k182QcParameters, k182Head(fixture), list(state_name = "protein_intensity_filtered", config = list(),
                description = "dispatch"))))
        if (!is.null(error)) cat(sprintf("KIT-182 F1 g%s dispatch tail: %s\n", case, conditionMessage(error)))
        expect_identical(fixture$manager$current_state, "protein_intensity_filtered", info = case)
        k182F1(fixture, paste0("g", case), TRUE)
    }
})

test_that("KIT-182 F1 (h): after the tab's restore, the state keeps its lineage", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    k182Correlation(fixture, norm, 0.5)
    # R1 instrument repairs (VAL-54, runs 45294678011 and 45302634386): the session export stopped on this fixture on
    # both commits, before and after the QC change ("saved table changed during read" from describeDIANativeFile: a dependency source that is not a
    # regular file), so no export runs here; the case checks the lineage after the tab's restore only, and
    # the export and reopen are recorded as not checked.
    k182QcChange(fixture)
    restored <- new.env(parent = emptyenv())
    protNormCoreQuiet(restoreProtNormTabState(NULL, fixture$workflow, restored))
    k182F1(fixture, "h", TRUE)
})

test_that("KIT-182 F1 (j), (k): an adopted legacy head after RUV has normalisation in its lineage; a re-saved design has not", {
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    k182Scope(fixture, function() {
        manager$saveState("protein_intensity_filtered", fixture$original, list(), "legacy adopted")
        snapshot <- manager$.context$snapshot
        snapshot$step_provenance$protein_intensity_filtered$adopted <- TRUE
        publishDIAStateHead(manager$.context, snapshot)
    })
    clauses <- k182F1(fixture, "j", TRUE)
    expect_identical(names(clauses)[clauses], "iii")
    # (k) the design saved again after RUV: the import table under its own name, then a QC change on it.
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    expect_true("ruv_corrected" %in% k182History(fixture))
    k182Scope(fixture, function() {
        setDIAStepProvenance(manager, "design_builder")
        manager$saveState("protein_replicate_filtered", fixture$original, list(), "design re-saved")
    })
    k182F1(fixture, "k", FALSE)
    k182QcChange(fixture)
    k182F1(fixture, "k-qc", FALSE)
})

# ------------------------------------------------------------------------------------------------
# T1 to T8: the normalisation's source (VAL-48)
# ------------------------------------------------------------------------------------------------

test_that("KIT-182 T1 (S2): a re-run after a Revert and a QC change reads the QC result, on DIA, LFQ and TMT, in process and managed (fails before)", {
    for (route in c("DIA", "LFQ", "TMT")) {
        fixture <- protNormCoreOwner()
        if (!identical(route, "DIA")) fixture$manager$setWorkflowType(route)
        prf1 <- k182Head(fixture)
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, "skip")
        k182Revert(fixture, "protein_replicate_filtered")
        q <- k182QcChange(fixture)
        k <- 64L - as.integer(k182Rows(fixture, q))
        expect_gte(k, 1L)
        # Managed: the action's input, on a copy of the tab's state (the runner runs nothing).
        managed <- k182ManagedClick(fixture, list2env(as.list(norm), parent = emptyenv()), "skip")
        input <- k182ActionInput(managed, "normalization-ruv")
        cat(sprintf("KIT-182 effect %s managed source=%s\n", route,
            if (identical(input, q$entity_id)) "current" else if (identical(input, prf1$entity_id)) "PRF1" else "other"))
        expect_identical(input, q$entity_id, info = route)
        # In process: the source recorded and the rows normalised.
        error <- k182Error(k182Click(fixture, norm, "skip"))
        if (!is.null(error)) {
            # F3: a failure unrelated to the source on a route the DIA fixture was not built for is recorded.
            cat(sprintf("KIT-182 F3 %s in_process error=%s\n", route, conditionMessage(error)))
            expect_identical(route %in% c("LFQ", "TMT"), TRUE, info = conditionMessage(error))
            next
        }
        expect_identical(norm$normalization_input_reference$entity_id, q$entity_id, info = route)
        rows <- as.integer(k182Rows(fixture, getDIAStateReference(fixture$manager, "normalized")))
        cat(sprintf("KIT-182 effect %s in_process rows=%d k=%d source=%s\n", route, rows, k,
            if (identical(norm$normalization_input_reference$entity_id, q$entity_id)) "current" else "PRF1"))
        expect_identical(rows, 64L - k, info = route)
    }
})

test_that("KIT-182 T1 (e) (S9): after the design is saved again, a QC change and a click read the QC result (fails before)", {
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    k182Scope(fixture, function() {
        setDIAStepProvenance(manager, "design_builder")
        manager$saveState("protein_replicate_filtered", fixture$original, list(), "design re-saved")
    })
    q <- k182QcChange(fixture)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, q$entity_id)
})

test_that("KIT-182 T2: settings-only re-runs read PRF1 every time and the repeat commits nothing (as before)", {
    fixture <- protNormCoreOwner()
    prf1 <- k182Head(fixture)$entity_id
    norm <- new.env(parent = emptyenv())
    source <- function() norm$normalization_input_reference$entity_id
    expect_true(k182Click(fixture, norm, "skip"))
    expect_identical(source(), prf1)
    head <- k182Head(fixture)
    history <- k182History(fixture)
    expect_true(k182Click(fixture, norm, "skip"))
    expect_identical(source(), prf1)
    expect_identical(k182Head(fixture), head)
    expect_identical(k182History(fixture), history)
    expect_true(k182Click(fixture, norm, "manual"))
    expect_identical(source(), prf1)
    expect_true(isSavedWorkflowReference(k182Correlation(fixture, norm, 0.5)))
    expect_true(k182Click(fixture, norm, "manual"))
    expect_identical(source(), prf1)
    expect_true(isSavedWorkflowReference(k182Correlation(fixture, norm, 0, skipped = TRUE)))
    expect_true(k182Click(fixture, norm, "manual", ruv_k = 1))
    expect_identical(source(), prf1)
})

test_that("KIT-182 T3 (S4, S5): a first run after a QC change reads the QC result, with and without tab entry (fails before)", {
    for (entry in c(TRUE, FALSE)) {
        fixture <- protNormCoreOwner()
        if (entry) protNormCoreQuiet(k182Scope(fixture, function() resolveProtNormQcStateObject(fixture$manager)))
        q <- k182QcChange(fixture)
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, "skip")
        expect_identical(norm$normalization_input_reference$entity_id, q$entity_id, info = paste("tab entry", entry))
    }
})

test_that("KIT-182 T4 (a) to (d) (S1): QC applied to a normalised table is refused at the next normalisation, in process, managed and on the tab (fails before)", {
    cases <- list(a = list(mode = "skip", kernel = "intensity"), b = list(mode = "manual", kernel = "intensity"),
        c = list(mode = "skip", kernel = "replicate"), d_apply = list(mode = "manual", correlation = FALSE),
        d_skip = list(mode = "manual", correlation = TRUE))
    for (case in names(cases)) {
        spec <- cases[[case]]
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        k182Click(fixture, norm, spec$mode)
        if (!is.null(spec$correlation)) k182Correlation(fixture, norm, 0.5, skipped = spec$correlation)
        k182QcKernel(fixture, spec$kernel %||% "intensity")
        head <- k182Head(fixture)
        history <- k182History(fixture)
        states <- names(fixture$manager$.context$snapshot$states)
        # In process.
        error <- k182Error(k182Click(fixture, norm, "skip"))
        expect_s3_class(error, "prot_norm_input_after_normalisation")
        # (c) re-commits `protein_replicate_filtered` itself, so the last table before normalisation by that name is the
        # QC-on-normalised one; the walk names the step before it.
        if (!identical(case, "c") && inherits(error, "condition")) expect_identical(conditionMessage(error),
            k182RunText("Protein Replicate Filter"), info = case)
        expect_identical(k182Head(fixture), head, info = case)
        expect_identical(k182History(fixture), history, info = case)
        expect_identical(names(fixture$manager$.context$snapshot$states), states, info = case)
        # Managed: refused in the controller; nothing dispatched.
        actions <- NULL
        error <- k182Error(actions <- k182ManagedClick(fixture, norm, "skip"))
        expect_s3_class(error, "prot_norm_input_after_normalisation")
        expect_length(actions, 0L)
        # The tab: the error handler once, no completion notice.
        tab <- k182Tab(fixture, norm, "skip")
        expect_length(tab$errors, 1L)
        expect_true(all(vapply(tab$errors, inherits, NA, "prot_norm_input_after_normalisation")), info = case)
        expect_false(any(grepl(buildProtNormCompletionNotification("skip"), tab$notes, fixed = TRUE)), info = case)
        expect_identical(k182Head(fixture), head, info = case)
    }
})

test_that("KIT-182 T4 (e): a remembered or QC-source table that is itself normalised is refused (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    norm$normalization_input_reference <- NULL
    protNormCoreQuiet(k182Scope(fixture, function() resolveProtNormQcStateObject(fixture$manager)))
    expect_identical(fixture$manager$.context$snapshot$qc_sources$pre_normalization$entity_id, k182Head(fixture)$entity_id)
    head <- k182Head(fixture)
    error <- k182Error(k182Click(fixture, norm, "skip"))
    expect_s3_class(error, "prot_norm_input_after_normalisation")
    expect_identical(k182Head(fixture), head)
})

test_that("KIT-182 T4 (f): on DIA_LIMPA, a table to which QC was applied after the RUV correction is not designated (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k182PassThrough(fixture, norm)
    k182Correlation(fixture, norm, 0.5)
    k182RuvCovariates(fixture, norm, "correction")
    expect_identical(fixture$manager$current_state, "ruv_corrected")
    protNormCoreQuiet(k182Scope(fixture, function() {
        value <- computeDIAReferenceProteinNormalisation(fixture$manager, getDIAStateReference(fixture$manager),
            list(norm_method = "none"))
        finishDIAReferenceProteinStep(fixture$workflow, value, publish = FALSE)
    }))
    designated <- norm$normalization_input_reference
    tab <- k182Tab(fixture, norm, "covariates")
    expect_length(tab$errors, 1L)
    expect_true(all(vapply(tab$errors, inherits, NA, "prot_norm_input_after_normalisation")))
    if (length(tab$errors)) expect_identical(conditionMessage(tab$errors[[1L]]), k182RunText("Peptide Count Filter"))
    expect_identical(norm$normalization_input_reference, designated)
})

test_that("KIT-182 T4 (g): an adopted legacy head after RUV is refused (fails before)", {
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    prf1 <- k182Head(fixture)$entity_id
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    k182Scope(fixture, function() {
        manager$saveState("protein_intensity_filtered", fixture$original, list(), "legacy adopted")
        snapshot <- manager$.context$snapshot
        snapshot$step_provenance$protein_intensity_filtered$adopted <- TRUE
        publishDIAStateHead(manager$.context, snapshot)
    })
    error <- k182Error(k182Click(fixture, norm, "skip"))
    expect_s3_class(error, "prot_norm_input_after_normalisation")
    if (is.null(error)) expect_identical(norm$normalization_input_reference$entity_id, prf1)
})

test_that("KIT-182 T5: Reset then run reads PRF1 (as before); after Reset, a QC change and a run read the QC result (fails before)", {
    reset <- function(fixture, norm) protNormCoreQuiet(k182Scope(fixture, function() runProtNormResetWorkflow(fixture$workflow,
        norm, output = new.env(), resetOutputsFn = function(...) NULL, showNotificationFn = function(...) NULL,
        messageFn = function(...) NULL)))
    # (a)
    fixture <- protNormCoreOwner()
    prf1 <- k182Head(fixture)$entity_id
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    reset(fixture, norm)
    expect_identical(k182Head(fixture)$entity_id, prf1)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, prf1)
    # (b) S6
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    protNormCoreQuiet(k182Scope(fixture, function() resolveProtNormQcStateObject(fixture$manager)))
    k182Click(fixture, norm, "skip")
    reset(fixture, norm)
    q <- k182QcChange(fixture)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, q$entity_id)
})

test_that("KIT-182 T6: a restored session re-runs from its own state as before, and from a QC change reads the QC result (fails before)", {
    restore <- function(fixture, prf1) {
        fixture$manager$.context$native_extra_references <- list(normalization_input_reference = prf1)
        restored <- new.env(parent = emptyenv())
        protNormCoreQuiet(restoreProtNormTabState(NULL, fixture$workflow, restored))
        restored
    }
    fixture <- protNormCoreOwner()
    prf1 <- k182Head(fixture)
    k182Click(fixture, new.env(parent = emptyenv()), "skip")
    # (a) at `normalized`, then at `correlation_filtered`.
    norm <- restore(fixture, prf1)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, prf1$entity_id)
    k182Correlation(fixture, norm, 0.5)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    norm <- restore(fixture, prf1)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, prf1$entity_id)
    # (b) a Revert and a QC change after the restore.
    norm <- restore(fixture, prf1)
    k182Revert(fixture, "protein_replicate_filtered")
    q <- k182QcChange(fixture)
    k182Click(fixture, norm, "skip")
    expect_identical(norm$normalization_input_reference$entity_id, q$entity_id)
})

test_that("KIT-182 T7 (a): the tab says when the run reads another table than the previous run, and only then (fails before)", {
    completion <- paste0("[message 10] ", buildProtNormCompletionNotification("skip"))
    # S2 through the tab.
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    expect_identical(k182Tab(fixture, norm, "skip")$notes, completion)
    k182Revert(fixture, "protein_replicate_filtered")
    k182QcChange(fixture)
    expect_identical(k182Tab(fixture, norm, "skip")$notes, c(paste0("[message 10] ", k182NoticeText), completion))
    # T2's sequence through the tab: the completion notice only.
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    for (mode in c("skip", "skip", "manual")) {
        expect_identical(k182Tab(fixture, norm, mode)$notes, paste0("[message 10] ", buildProtNormCompletionNotification(mode)),
            info = mode)
    }
    # A restored session with no normalisation, a QC change and a run: no notice.
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    protNormCoreQuiet(restoreProtNormTabState(NULL, fixture$workflow, norm))
    expect_false(isTRUE(norm$normalization_complete))
    k182QcChange(fixture)
    expect_identical(k182Tab(fixture, norm, "skip")$notes, completion)
})

test_that("KIT-182 T7 (b): the DIA_LIMPA pass-through designates the current QC result, as before", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k182PassThrough(fixture, norm)
    k182Revert(fixture, "protein_normalised")
    p2 <- k182LimpaPeptideCount(fixture, 2)
    k182PassThrough(fixture, norm)
    expect_identical(norm$normalization_input_reference$entity_id, p2$entity_id)
})

test_that("KIT-182 T8 (in process): the pre-normalisation figures are redrawn from the table normalised, before the publication (fails before)", {
    log <- new.env(parent = emptyenv())
    publish <- publishDIAActionOutputs
    testthat::local_mocked_bindings(publishDIAActionOutputs = function(...) {
        log$calls <- c(log$calls, "publish")
        publish(...)
    }, .package = "MultiScholaR")
    for (sequence in c("S2", "S4")) {
        fixture <- protNormCoreOwner()
        manager <- fixture$manager
        norm <- new.env(parent = emptyenv())
        generators <- protNormQcGenerators(fixture$workflow, norm, fixture$paths, k182Aesthetics)
        pre_qc <- function() {
            log$calls <- c(log$calls, "pre_qc")
            generators$pre()
        }
        image <- file.path(fixture$paths$protein_qc_dir, "pre_norm_pca.png")
        # Tab entry with head PRF1: the QC source is set and the figures drawn.
        protNormCoreQuiet(k182Scope(fixture, function() {
            resolveProtNormQcStateObject(manager)
            generators$pre()
        }))
        drawn <- unname(tools::md5sum(image))
        if (identical(sequence, "S2")) {
            k182Click(fixture, norm, "skip", pre_qc = pre_qc)
            k182Revert(fixture, "protein_replicate_filtered")
        }
        log$calls <- character()
        k182QcChange(fixture)
        k182Click(fixture, norm, "skip", pre_qc = pre_qc)
        calls <- log$calls
        source <- norm$normalization_input_reference$entity_id
        # Consistent on both commits: the figures' source is the run's source.
        expect_identical(manager$.context$snapshot$qc_sources$pre_normalization$entity_id, source, info = sequence)
        expect_identical(sum(calls == "pre_qc"), 1L, info = sequence)
        expect_true(isTRUE(match("pre_qc", calls) < match("publish", calls)), info = sequence)
        expect_false(identical(unname(tools::md5sum(image)), drawn), info = sequence)
    }
})

# ------------------------------------------------------------------------------------------------
# T9: a protein QC Apply on a normalised table is refused (VAL-91 D1)
# ------------------------------------------------------------------------------------------------

k182QcSteps <- function(fixture) list(
    cleanup = function() runProteinAccessionCleanupStep(fixture$workflow, ";", "mean"),
    intensity = function() runProteinIntensityFilterApplyStep(fixture$workflow, FALSE, 2, 2, 10),
    duplicates = function() runProteinDuplicateRemovalStep(fixture$workflow, "mean"),
    replicate = function() runProteinReplicateFilterApplyStep(fixture$workflow, fixture$paths, "group", 1))

k182ExpectQcRefused <- function(fixture, step, name, case, title) {
    head <- k182Head(fixture)
    history <- k182History(fixture)
    error <- k182Error(protNormCoreQuiet(k182Scope(fixture, step)))
    expect_s3_class(error, "prot_qc_input_normalised")
    if (inherits(error, "prot_qc_input_normalised")) expect_identical(conditionMessage(error), k182QcText(title),
        info = paste(case, name))
    expect_identical(k182Head(fixture), head, info = paste(case, name))
    expect_identical(k182History(fixture), history, info = paste(case, name))
}

test_that("KIT-182 T9 (a): the four DIA protein QC steps refuse a normalised head and run on PRF1 (fails before)", {
    heads <- list(
        normalized = function(fixture, norm) k182Click(fixture, norm, "skip"),
        correlation_filtered = function(fixture, norm) {
            k182Click(fixture, norm, "skip")
            k182Correlation(fixture, norm, 0.5)
        },
        qc_on_normalised = function(fixture, norm) {
            k182Click(fixture, norm, "skip")
            k182QcChange(fixture)
        },
        enrichment_completed = function(fixture, norm) {
            k182Click(fixture, norm, "skip")
            k182Scope(fixture, function() {
                record <- readDIAStateMetadata(fixture$manager)
                commitDIAReferenceState(fixture$manager, "enrichment_completed", record$data, record$config, "fixture")
            })
        })
    # One owner per head; on the fix each refused step leaves the head as it was for the next.
    for (case in names(heads)) {
        fixture <- protNormCoreOwner()
        norm <- new.env(parent = emptyenv())
        heads[[case]](fixture, norm)
        steps <- k182QcSteps(fixture)
        for (name in names(steps)) k182ExpectQcRefused(fixture, steps[[name]], name, case, "Protein Replicate Filter")
    }
    # On PRF1 the step runs.
    fixture <- protNormCoreOwner()
    error <- k182Error(protNormCoreQuiet(k182Scope(fixture, k182QcSteps(fixture)$intensity)))
    expect_null(error)
    expect_identical(fixture$manager$current_state, "protein_intensity_filtered")
})

test_that("KIT-182 T9 (b): DIA_LIMPA's between-sample normalisation and peptide count filter refuse a normalised head (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k182PassThrough(fixture, norm)
    k182Correlation(fixture, norm, 0.5)
    k182RuvCovariates(fixture, norm, "covariates")
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    steps <- list(normalisation = function() runProteinNormalisationStep(fixture$workflow, "none"),
        peptide_count = function() runProteinPeptideCountStep(fixture$workflow, 2, 2))
    for (name in names(steps)) k182ExpectQcRefused(fixture, steps[[name]], name, "ruv_covariates", "Peptide Count Filter")
    k182Revert(fixture, "correlation_filtered")
    for (name in names(steps)) k182ExpectQcRefused(fixture, steps[[name]], name, "correlation_filtered", "Peptide Count Filter")
    k182Revert(fixture, "protein_normalised")
    error <- k182Error(protNormCoreQuiet(k182Scope(fixture, steps$peptide_count)))
    expect_null(error)
    expect_identical(fixture$manager$current_state, "peptide_count_filtered")
})

#' A managed protein QC sub-tab under testServer: `apply` recorded, error notices recorded.
k182ManagedServer <- function(fixture, server, button) {
    seen <- new.env(parent = emptyenv())
    seen$applied <- 0L
    seen$errors <- character()
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(session) "managed",
        installWorkflowManagedStep = function(...) {
            step <- new.env(parent = emptyenv())
            step$apply <- function(...) seen$applied <- seen$applied + 1L
            step$cancel <- function(...) invisible(NULL)
            step$active <- function() NULL
            step
        }, .package = "MultiScholaR")
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        if (identical(type, "error")) seen$errors <- c(seen$errors, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }, removeNotification = function(...) invisible(NULL), .package = "shiny")
    protNormCoreQuiet(shiny::testServer(server, args = list(workflow_data = fixture$workflow,
        experiment_paths = fixture$paths, omic_type = "proteomics", experiment_label = "core"), {
        do.call(session$setInputs, stats::setNames(list(1), button))
        session$flushReact()
    }))
    as.list(seen)
}

test_that("KIT-182 T9 (c): managed, the protein intensity filter refuses before apply with the step's prefix, and applies on an unnormalised head (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    refused <- k182ManagedServer(fixture, mod_prot_qc_protein_intensity_server, "apply_protein_intensity_filter")
    expect_identical(refused$applied, 0L)
    expect_identical(refused$errors, paste0("Error applying protein intensity filter: ", k182QcText("Protein Replicate Filter")))
    k182Revert(fixture, "protein_replicate_filtered")
    allowed <- k182ManagedServer(fixture, mod_prot_qc_protein_intensity_server, "apply_protein_intensity_filter")
    expect_identical(allowed$applied, 1L)
    expect_length(allowed$errors, 0L)
})

test_that("KIT-182 T9 (c): managed, DIA_LIMPA's peptide count filter refuses before apply and applies after a Revert (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k182PassThrough(fixture, norm)
    k182Correlation(fixture, norm, 0.5)
    k182RuvCovariates(fixture, norm, "covariates")
    refused <- k182ManagedServer(fixture, mod_prot_qc_protein_peptide_count_server, "apply_peptide_count")
    expect_identical(refused$applied, 0L)
    expect_identical(refused$errors, paste0("Error applying peptide count filter: ", k182QcText("Peptide Count Filter")))
    k182Revert(fixture, "protein_normalised")
    allowed <- k182ManagedServer(fixture, mod_prot_qc_protein_peptide_count_server, "apply_peptide_count")
    expect_identical(allowed$applied, 1L)
})

test_that("KIT-182 T9 (d), (e): a re-saved design and a non-reference manager are not refused (as before)", {
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "manual")
    k182Scope(fixture, function() {
        setDIAStepProvenance(manager, "design_builder")
        manager$saveState("protein_replicate_filtered", fixture$original, list(), "design re-saved")
    })
    managed <- k182ManagedServer(fixture, mod_prot_qc_protein_intensity_server, "apply_protein_intensity_filter")
    expect_identical(managed$applied, 1L)
    error <- k182Error(protNormCoreQuiet(k182Scope(fixture, k182QcSteps(fixture)$intensity)))
    expect_null(error)
    expect_identical(manager$current_state, "protein_intensity_filtered")
    # (e) the module contract tests' harness: an environment, not the reference owner.
    harness <- new.env(parent = emptyenv())
    harness$current_state <- "normalized"
    harness$getHistory <- function() list("protein_replicate_filtered", "normalized")
    workflow <- list2env(list(state_manager = harness, qc_params = list()))
    error <- k182Error(protNormCoreQuiet(runProteinIntensityFilterApplyStep(workflow, FALSE, 2, 2, 10)))
    expect_false(inherits(error, "prot_qc_input_normalised"))
})

# ------------------------------------------------------------------------------------------------
# T10: the correlation filter after a QC change (VAL-91 D2) and after a Revert between two runs (VAL-99 (4))
# ------------------------------------------------------------------------------------------------

#' Normalised (RUV skipped), reverted to PRF1 and the QC changed: `final` is the run's result the tab holds.
k182Stale <- function(envir = parent.frame()) {
    fixture <- protNormCoreOwner(envir)
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    final <- norm$ruv_normalized_obj
    k182Revert(fixture, "protein_replicate_filtered")
    q <- k182QcChange(fixture)
    list(fixture = fixture, norm = norm, final = final, q = q)
}

k182ExpectCorrelationRefused <- function(fixture, run, case) {
    head <- k182Head(fixture)
    history <- k182History(fixture)
    error <- k182Error(run())
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    if (inherits(error, "prot_norm_correlation_input_stale")) expect_identical(conditionMessage(error), k182CorrelationText,
        info = case)
    expect_identical(k182Head(fixture), head, info = case)
    expect_identical(k182History(fixture), history, info = case)
    expect_null(getDIAStateReference(fixture$manager, "correlation_filtered"), info = case)
}

test_that("KIT-182 T10 (a) to (c): Apply, Skip and the managed Apply refuse after a QC change (fails before)", {
    for (skipped in c(FALSE, TRUE)) {
        stale <- k182Stale()
        expect_identical(as.integer(k182Rows(stale$fixture, stale$final)), 64L)
        k182ExpectCorrelationRefused(stale$fixture, function() k182Correlation(stale$fixture, stale$norm, 0.5, skipped),
            if (skipped) "skip" else "apply")
    }
    stale <- k182Stale()
    seen <- new.env(parent = emptyenv())
    seen$actions <- list()
    seen$failed <- list()
    runner <- list(run = function(actions, done = NULL, failed = NULL) {
        if (is.function(actions)) actions <- tryCatch(actions(), error = function(e) {
            failed(e)
            NULL
        })
        seen$actions <- c(seen$actions, actions)
        invisible(TRUE)
    })
    protNormCoreQuiet(dispatchProtNormCorrelation(runner, stale$fixture$workflow, stale$norm, stale$fixture$paths, 0.5,
        "proteomics", "core", function() "group", k182Aesthetics, function(final) final,
        function(error) seen$failed <- c(seen$failed, list(error))))
    expect_length(seen$failed, 1L)
    expect_true(all(vapply(seen$failed, inherits, NA, "prot_norm_correlation_input_stale")))
    expect_length(seen$actions, 0L)
})

test_that("KIT-182 T10 (d): the tab's Apply and Skip show the refusal through the correlation error handler (fails before)", {
    stale <- k182Stale()
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$completed <- 0L
    handler <- function(error, logPrefix, notificationPrefix) handleProtNormCorrelationError(error, logPrefix,
        notificationPrefix, showNotificationFn = function(ui, ...) seen$notes <- c(seen$notes, ui),
        messageFn = function(...) NULL)
    complete <- function(...) seen$completed <- seen$completed + 1L
    protNormCoreQuiet(runProtNormApplyCorrelationObserver(list(min_pearson_correlation_threshold = 0.5), new.env(),
        stale$fixture$workflow, stale$norm, stale$fixture$paths, "proteomics", "core", function() "group", k182Aesthetics,
        function(...) NULL, completeCorrelationWorkflowFn = complete, handleCorrelationErrorFn = handler))
    protNormCoreQuiet(runProtNormSkipCorrelationObserver(new.env(), stale$fixture$workflow, stale$norm, stale$fixture$paths,
        "proteomics", "core", k182Aesthetics, function(...) NULL, completeCorrelationWorkflowFn = complete,
        handleCorrelationErrorFn = handler))
    expect_identical(seen$completed, 0L)
    expect_identical(seen$notes, rep(paste("Error in correlation filtering:", k182CorrelationText), 2L))
})

test_that("KIT-182 T10 (e): the correlation at the tab's own states commits as before; at a QC-on-normalised head it is refused", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    for (step in list(list(0.5, FALSE), list(0.7, FALSE), list(0, TRUE))) {
        final <- k182Correlation(fixture, norm, step[[1L]], step[[2L]])
        expect_true(isSavedWorkflowReference(final))
        expect_identical(fixture$manager$current_state, "correlation_filtered")
    }
    # A T4 (d) head: the QC change applied at `correlation_filtered` (built by the kernel).
    k182QcChange(fixture)
    head <- k182Head(fixture)
    error <- k182Error(k182Correlation(fixture, norm, 0.5))
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    expect_identical(k182Head(fixture), head)
})

test_that("KIT-182 T10 (f): on DIA_LIMPA the correlation follows the designated table; after a QC change it waits for the pass-through (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- protNormW1LimpaOwner()
    norm <- new.env(parent = emptyenv())
    k182PassThrough(fixture, norm)
    p <- norm$normalization_input_reference
    expect_identical(p$entity_id, k182Head(fixture)$entity_id)
    expect_null(k182Error(k182Correlation(fixture, norm, 0.5)))
    k182RuvCovariates(fixture, norm, "covariates")
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    expect_null(k182Error(k182Correlation(fixture, norm, 0.5)))
    k182Revert(fixture, "protein_normalised")
    p2 <- k182LimpaPeptideCount(fixture, 2)
    error <- k182Error(k182Correlation(fixture, norm, 0.5))
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    expect_identical(k182Head(fixture), p2)
    # The pass-through again from the QC result (where the refused press left it).
    if (!identical(fixture$manager$current_state, "peptide_count_filtered")) k182Revert(fixture, "peptide_count_filtered")
    k182PassThrough(fixture, norm)
    expect_identical(norm$normalization_input_reference$entity_id, k182Head(fixture)$entity_id)
    designated <- norm$normalization_input_reference$entity_id
    expect_null(k182Error(k182Correlation(fixture, norm, 0.5)))
    expect_identical(norm$normalization_input_reference$entity_id, designated)
})

test_that("KIT-182 T10 (g) (VAL-99 (4)): after two runs and a Revert to the first, the correlation is refused (fails before)", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    k182Click(fixture, norm, "skip")
    k182Click(fixture, norm, "manual")
    manual_final <- norm$ruv_normalized_obj
    protNormCoreQuiet(k182Scope(fixture, function() runProteinIntensityFilterRevertStep(fixture$workflow)))
    expect_identical(fixture$manager$current_state, "normalized")
    head <- k182Head(fixture)
    error <- k182Error(k182Correlation(fixture, norm, 0.5))
    expect_s3_class(error, "prot_norm_correlation_input_stale")
    expect_identical(k182Head(fixture), head)
    expect_identical(norm$ruv_normalized_obj, manual_final)
})
