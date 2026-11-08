# KIT-104 W2: the QC core for lipid, metab and trx, written before it exists. Each shared function is compared with
# today's lipid, metab and trx function on the MCI fixtures (helper-module-ci-{lipid,metab}-qc.R) and the trx_screen
# corpus: the counts table, the duplicate detection and resolution statistics, the ITSD analysis and its figures, the
# finalise validation and its alignment, and the QC action's two branches (D6). Back-reference #9
# (`resolveMetabDuplicateAssayData`) and trx's QC figure functions are checked where they moved.

#' Run `code` with a namespace function wrapped so that its result is kept in `captured[[name]]`.
qcCoreSpy <- function(name, captured, env = parent.frame()) {
    real <- get(name, envir = asNamespace("MultiScholaR"))
    spy <- function(...) {
        result <- real(...)
        captured[[name]] <- result
        result
    }
    do.call(testthat::local_mocked_bindings, c(stats::setNames(list(spy), name), list(.package = "MultiScholaR", .env = env)))
}

#' Capture shiny's notifications as "[type] text" (test-omic-qc-module.R's construction).
qcCoreNotes <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        seen$notes <- c(seen$notes, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }, .package = "shiny", .env = env)
    seen
}

#' The ggplot of a figure with its built layer data, for comparing two builders of the same figure.
qcCorePlotParts <- function(plot) {
    list(data = plot$data, labels = plot$labels, layers = ggplot2::ggplot_build(plot)$data)
}

#' The omic's qc configuration with the finalise's alignment (the lead wires it at the phase B merge).
qcCoreConfig <- function(omic_type) {
    qc <- omicGuiConfig(omic_type)$qc
    qc$finalize_alignment <- switch(omic_type,
        transcriptomics = trxDesignAlignment,
        function(shell) omicsQcFinalizeAlignment(shell, omic_type))
    qc
}

#' An object whose first design run has a blank group: valid for the class, refused by the alignment.
qcCoreBlankGroup <- function(object) {
    object@design_matrix[[object@group_id]][[1L]] <- NA_character_
    object
}

#' The message of an error, or NULL when `code` succeeds.
qcCoreError <- function(code) tryCatch({
    force(code)
    NULL
}, error = function(error) conditionMessage(error))

# ------------------------------------------------------------------------------------------------
# The intensity filter's counts table.
# ------------------------------------------------------------------------------------------------

test_that("omicsQcCountsTable is today's lipid (inline) and metab (buildMetabIntensityFilterStats) table", {
    workflow <- module_ci_lipid_qc_workflow()
    captured <- new.env(parent = emptyenv())
    qcCoreSpy("runLipidReferenceIntensityFilter", captured)
    stats <- shiny::reactiveVal(NULL)
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
        registerLipidIntensityApplyFilterObserver(input, workflow, output, stats, shiny::reactiveVal(NULL), "lipidomics",
            updateLipidFilteringFn = function(...) NULL)
    }, {
        session$setInputs(intensity_cutoff_percentile = 40L, proportion_below_cutoff = 0.5, apply_filter = 1L)
        captured$stats <- stats()
    }))))
    counts <- captured$runLipidReferenceIntensityFilter$value$counts
    expect_true(is.data.frame(captured$stats))
    expect_identical(omicsQcCountsTable(counts), captured$stats)
    expect_identical(omicsQcCountsTable(counts), buildMetabIntensityFilterStats(originalCounts = counts$before,
        filteredCounts = counts$after))
    # Rounding to one decimal, as both tabs round.
    odd <- list(before = c(A = 3L, B = 7L), after = c(A = 1L, B = 7L))
    expect_identical(omicsQcCountsTable(odd)$Percent_Retained, c(33.3, 100))
    expect_identical(omicsQcCountsTable(odd), buildMetabIntensityFilterStats(originalCounts = odd$before,
        filteredCounts = odd$after))
})

# ------------------------------------------------------------------------------------------------
# Duplicates.
# ------------------------------------------------------------------------------------------------

test_that("detectWorkflowOmicsDuplicates finds today's duplicates (lipid and metab detection)", {
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    today <- handleLipidDuplicateDetection(lipid, logInfoFn = function(...) invisible(NULL))
    core <- detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, lipid$state_manager)
    expect_identical(core$duplicates, today$duplicatesList)
    expect_equal(core$total, today$totalDuplicates)
    expect_gt(core$total, 0)

    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    today <- detectMetabDuplicateFeatures(metab$state_manager)
    core <- detectWorkflowOmicsDuplicates(metabolomicsDeclarationBinding, metab$state_manager)
    expect_identical(core$duplicates, today$duplicatesList)
    expect_equal(core$total, today$totalDuplicates)

    none <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(include_duplicates = FALSE))
    expect_equal(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, none$state_manager)$total,
        handleLipidDuplicateDetection(none, logInfoFn = function(...) invisible(NULL))$totalDuplicates)
    # Another omic's state is refused as today (the id reader fails on the binding's id slot).
    today <- qcCoreError(handleLipidDuplicateDetection(metab, logInfoFn = function(...) invisible(NULL)))
    expect_false(is.null(today))
    expect_identical(qcCoreError(detectWorkflowOmicsDuplicates(lipidomicsDeclarationBinding, metab$state_manager)), today)
})

test_that("omicsDuplicateResolutionStats reads both products as today's lipid and metab resolutions do", {
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    captured <- new.env(parent = emptyenv())
    qcCoreSpy("runLipidReferenceDuplicateResolution", captured)
    names_before <- readLipidWorkflowMetadata(lipid$state_manager)$assay_names
    today <- handleLipidDuplicateResolution(lipid, "lipidomics", updateFilteringFn = function(...) NULL,
        logWarnFn = function(...) invisible(NULL))
    value <- captured$runLipidReferenceDuplicateResolution$value
    expect_true("counts" %in% names(value) && !"stats" %in% names(value))
    expect_identical(omicsDuplicateResolutionStats(value, names_before), today$statsList)

    metab <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    qcCoreSpy("runMetabReferenceDuplicateResolution", captured)
    shown <- new.env(parent = emptyenv())
    runMetabDuplicateResolutionWorkflow(metab, "metabolomics", setResolutionStatsFn = function(stats) shown$stats <- stats,
        setFilterPlotFn = function(plot) invisible(NULL),
        applyResolutionStateFn = function(currentS4, statsList, setResolutionStatsFn, ...) {
            setResolutionStatsFn(statsList)
            list(stateName = "metab_duplicates_resolved")
        })
    value <- captured$runMetabReferenceDuplicateResolution$value
    expect_true("stats" %in% names(value))
    expect_identical(omicsDuplicateResolutionStats(value, readMetabWorkflowMetadata(metab$state_manager)$assay_names),
        shown$stats)
    # An empty statistics list stays empty (metab's "No assay statistics available").
    expect_identical(omicsDuplicateResolutionStats(list(stats = list()), character()), list())
})

test_that("#9: metabolomicsDeclarationBinding$resolve_duplicates works through the moved resolveMetabDuplicateAssayData", {
    expect_identical(basename(utils::getSrcFilename(resolveMetabDuplicateAssayData, full.names = TRUE)),
        "func_metab_s4_qc_methods.R")
    object <- module_ci_metab_qc_object(layout = "combined")
    resolved <- metabolomicsDeclarationBinding$resolve_duplicates(object)
    direct <- resolveMetabDuplicateAssayData(assayList = object@metabolite_data,
        metaboliteIdCol = object@metabolite_id_column, logWarnFn = function(...) NULL)
    expect_identical(resolved$value, list(stats = direct$statsList))
    expect_identical(resolved$object@metabolite_data, direct$resolvedAssayList)
    expect_true(all(vapply(direct$statsList, `[[`, numeric(1), "removed") >= 2))
    expect_false(any(vapply(resolved$object@metabolite_data, function(assay) anyDuplicated(assay$database_identifier) > 0L,
        logical(1))))
})

# ------------------------------------------------------------------------------------------------
# ITSD.
# ------------------------------------------------------------------------------------------------

#' Today's lipid analysis: the exported ITSD server's inline analysis, its values and its text.
qcCoreLipidItsd <- function(workflow, pattern) {
    notes <- qcCoreNotes()
    out <- new.env(parent = emptyenv())
    suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(mod_lipid_qc_itsd_server,
        args = list(workflow_data = workflow, omic_type = "lipidomics", experiment_label = "kit104"), {
            session$setInputs(is_pattern = pattern, analyze_is = 1L)
            out$metrics <- is_metrics()
            out$long <- is_data_long()
            out$text <- tryCatch(output$is_results, error = function(error) NULL)
        }))))
    out$notes <- notes$notes
    as.list(out)
}

test_that("analyzeWorkflowOmicsItsd gives lipid's inline analysis: values, text and refusals", {
    itsd <- omicGuiConfig("lipidomics")$qc$itsd
    for (pattern in c("", "^IS_By_Name$", "^NO_MATCH$")) {
        workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
        today <- qcCoreLipidItsd(workflow, pattern)
        core <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, pattern, itsd)
        expect_identical(core$metrics, today$metrics, info = pattern)
        expect_identical(core$long_data, today$long, info = pattern)
        expect_identical(core$result_text, today$text, info = pattern)
        expect_identical(sprintf("Found %d internal standards", core$total), sub("^\\[message\\] ", "", today$notes[[2L]]))
    }
    # The refusals, with their texts: no standard found (the tip names lipid's column), no pattern at all.
    cases <- list(list(object = module_ci_lipid_qc_object(layout = "gc", include_itsd = FALSE), pattern = "^IS_"),
        list(object = module_ci_lipid_qc_object(layout = "gc", internal_standard_regex = NA_character_), pattern = ""))
    for (case in cases) {
        workflow <- module_ci_lipid_qc_workflow(case$object)
        today <- qcCoreLipidItsd(workflow, case$pattern)
        message <- qcCoreError(analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, case$pattern, itsd))
        expect_false(is.null(message))
        expect_identical(paste("[error] Error analyzing internal standards:", message), today$notes[[2L]])
    }
    expect_match(message, "No internal standard pattern provided", fixed = TRUE)
})

test_that("F104 (W2): a NULL pattern is read as an empty one, as metab reads it (lipid's inline analysis errors)", {
    itsd <- omicGuiConfig("lipidomics")$qc$itsd
    workflow <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    # Today's lipid server fails on a NULL input (nzchar(NULL)); only a test session leaves the input NULL.
    today <- qcCoreLipidItsd(workflow, NULL)
    expect_identical(today$notes[[2L]], "[error] Error analyzing internal standards: argument is of length zero")
    expect_identical(analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, NULL, itsd),
        analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, workflow$state_manager, "", itsd))
})

test_that("analyzeWorkflowOmicsItsd gives metab's analysis (analyzeMetabQcItsdData on readMetabQcItsdState)", {
    itsd <- omicGuiConfig("metabolomics")$qc$itsd
    for (pattern in list("", "^IS_", "^NO_MATCH$", NULL)) {
        workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
        manager <- workflow$state_manager
        today <- analyzeMetabQcItsdData(readMetabQcItsdState(manager, pattern), inputPattern = pattern,
            logInfoFn = function(...) invisible(NULL))
        core <- analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, pattern, itsd)
        info <- paste("pattern", format(pattern))
        expect_identical(core$metrics, today$metrics, info = info)
        expect_identical(core$long_data, today$longData, info = info)
        expect_identical(core$result_text, today$resultText, info = info)
        expect_identical(core$total, today$nIsTotal, info = info)
        expect_identical(core$pattern, today$pattern, info = info)
        expect_identical(core$metrics_by_assay, today$metricsByAssay, info = info)
    }
    for (case in list(list(object = module_ci_metab_qc_object(layout = "gc", include_itsd = FALSE), pattern = "^IS_"),
        list(object = module_ci_metab_qc_object(layout = "gc", internal_standard_regex = NA_character_), pattern = ""))) {
        manager <- module_ci_metab_qc_workflow(case$object)$state_manager
        today <- qcCoreError(analyzeMetabQcItsdData(readMetabQcItsdState(manager, case$pattern), inputPattern = case$pattern,
            logInfoFn = function(...) invisible(NULL)))
        expect_false(is.null(today))
        expect_identical(qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, case$pattern, itsd)), today)
    }
    # Another omic's state is refused as today (the shell reader fails on the binding's id slot).
    manager <- module_ci_lipid_qc_workflow()$state_manager
    today <- qcCoreError(analyzeMetabQcItsdData(readMetabQcItsdState(manager, ""), inputPattern = ""))
    expect_false(is.null(today))
    expect_identical(qcCoreError(analyzeWorkflowOmicsItsd(metabolomicsDeclarationBinding, manager, "", itsd)), today)
})

test_that("omicsItsdCvPlot and omicsItsdIntensityPlot are today's ITSD figures", {
    workflow <- module_ci_metab_qc_workflow(module_ci_metab_qc_object(layout = "combined"))
    analysis <- analyzeMetabQcItsdData(readMetabQcItsdState(workflow$state_manager, ""), inputPattern = "",
        logInfoFn = function(...) invisible(NULL))
    expect_identical(qcCorePlotParts(omicsItsdCvPlot(analysis$metrics)), qcCorePlotParts(buildMetabQcItsdCvPlot(analysis$metrics)))
    expect_identical(qcCorePlotParts(omicsItsdIntensityPlot(analysis$longData)),
        qcCorePlotParts(buildMetabQcItsdIntensityPlot(analysis$longData)))
    # Lipid's inline CV figure: the data its output draws.
    drawn <- new.env(parent = emptyenv())
    ggplot <- ggplot2::ggplot
    testthat::local_mocked_bindings(ggplot = function(data = NULL, ...) {
        drawn$data <- c(drawn$data, list(data))
        ggplot(data, ...)
    }, .package = "ggplot2")
    lipid <- module_ci_lipid_qc_workflow(module_ci_lipid_qc_object(layout = "combined"))
    metrics <- analyzeWorkflowOmicsItsd(lipidomicsDeclarationBinding, lipid$state_manager, "",
        omicGuiConfig("lipidomics")$qc$itsd)$metrics
    suppressWarnings(utils::capture.output(shiny::testServer(function(input, output, session) {
        registerLipidItsdCvPlotOutput(output, shiny::reactiveVal(metrics))
    }, {
        output$cv_plot
    })))
    expect_length(drawn$data, 1L)
    expect_identical(omicsItsdCvPlot(metrics)$data, drawn$data[[1L]])
})

# ------------------------------------------------------------------------------------------------
# The finalise.
# ------------------------------------------------------------------------------------------------

test_that("validateWorkflowOmicsQcFinalize is lipid's finalise validation: passes, refusals and their texts", {
    qc <- qcCoreConfig("lipidomics")
    aligned <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager, measures = "features")
    expect_identical(validateWorkflowOmicsQcFinalize(aligned, qc, lipidomicsDeclarationBinding),
        validateLipidQcS4FinalizeState(aligned))
    misaligned <- readLipidWorkflowShell(module_ci_lipid_qc_workflow(qcCoreBlankGroup(module_ci_lipid_qc_object()))$state_manager)
    zero_rows <- aligned
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    no_design <- aligned
    no_design@design_matrix <- no_design@design_matrix[0L, , drop = FALSE]
    no_assays <- aligned
    no_assays@lipid_data <- list()
    attr(no_assays, "assay_rows") <- integer()
    other_class <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager)
    for (case in list(misaligned, zero_rows, no_design, no_assays, other_class)) {
        today <- qcCoreError(validateLipidQcS4FinalizeState(case))
        expect_false(is.null(today))
        expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(case, qc, lipidomicsDeclarationBinding)), today)
    }
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(misaligned, qc, lipidomicsDeclarationBinding)),
        "Current LipidomicsAssayData object has assay/design mismatch: Design matrix group column 'group' contains blank assignments.")
})

test_that("validateWorkflowOmicsQcFinalize for metab: today's passes and refusals, and the alignment (parity)", {
    qc <- qcCoreConfig("metabolomics")
    aligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager, measures = "features")
    expect_identical(validateWorkflowOmicsQcFinalize(aligned, qc, metabolomicsDeclarationBinding),
        validateMetabQcS4FinalizeState(aligned))
    # D2: metab refuses only an empty assay list; a zero-row assay passes, as today.
    zero_rows <- aligned
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    expect_identical(validateWorkflowOmicsQcFinalize(zero_rows, qc, metabolomicsDeclarationBinding),
        validateMetabQcS4FinalizeState(zero_rows))
    no_design <- aligned
    no_design@design_matrix <- no_design@design_matrix[0L, , drop = FALSE]
    no_assays <- aligned
    no_assays@metabolite_data <- list()
    attr(no_assays, "assay_rows") <- integer()
    other_class <- readLipidWorkflowShell(module_ci_lipid_qc_workflow()$state_manager)
    for (case in list(no_design, no_assays, other_class)) {
        today <- qcCoreError(validateMetabQcS4FinalizeState(case))
        expect_false(is.null(today))
        expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(case, qc, metabolomicsDeclarationBinding)), today)
    }
    # Parity (D11): a blank group passes today's metab finalise and is refused by the core with lipid's text.
    misaligned <- readMetabWorkflowShell(module_ci_metab_qc_workflow(qcCoreBlankGroup(module_ci_metab_qc_object()))$state_manager)
    expect_s4_class(validateMetabQcS4FinalizeState(misaligned), "MetaboliteAssayData")
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(misaligned, qc, metabolomicsDeclarationBinding)),
        "Current MetaboliteAssayData object has assay/design mismatch: Design matrix group column 'group' contains blank assignments.")
    # Every omic's alignment is called: a configuration without one is refused, not skipped.
    qc$finalize_alignment <- NULL
    expect_error(validateWorkflowOmicsQcFinalize(aligned, qc, metabolomicsDeclarationBinding), class = "omic_gui_config_invalid")
})

test_that("omicsQcFinalizeAlignment is lipid's alignment call (validateLipidDesignAlignment with its arguments)", {
    for (object in list(module_ci_lipid_qc_object(), qcCoreBlankGroup(module_ci_lipid_qc_object()))) {
        shell <- readLipidWorkflowShell(module_ci_lipid_qc_workflow(object)$state_manager)
        runs <- as.character(shell@design_matrix[[shell@sample_id]])
        today <- validateLipidDesignAlignment(designMatrix = shell@design_matrix, assayList = shell@lipid_data,
            columnMapping = list(lipid_id_col = shell@lipid_id_column, annotation_col = shell@annotation_id_column,
                sample_columns = runs), sampleIdCol = shell@sample_id, groupCol = shell@group_id, requireReplicates = FALSE)
        expect_identical(omicsQcFinalizeAlignment(shell, "lipidomics"), today[c("valid", "errors")])
    }
    shell <- readMetabWorkflowShell(module_ci_metab_qc_workflow()$state_manager)
    expect_identical(omicsQcFinalizeAlignment(shell, "metabolomics"), list(valid = TRUE, errors = character()))
})

test_that("trxDesignAlignment passes on the trx_screen shell and refuses what validateTrxDesign refuses", {
    skip_if_not(workflowOwnerCapable())
    fixture <- trxDeclarationsOwner(withr::local_tempdir(), module_ci_trx_corpus())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    # The trx class's validity refuses a zero-row count table (W2 finding B: readWorkflowOmicsShell cannot read a trx
    # state), so the shell is the object with its zero-row prototype and stored row count, built without validity.
    shell <- fixture$object
    shell@count_data <- lapply(shell@count_data, function(table) table[0L, , drop = FALSE])
    attr(shell, "assay_rows") <- vapply(fixture$object@count_data, nrow, integer(1))
    expect_identical(trxDesignAlignment(shell), list(valid = TRUE, errors = character()))
    # The same check as every trx step (trxAssayInputs) on the whole object.
    expect_type(trxAssayInputs(fixture$object), "list")
    blank <- qcCoreBlankGroup(shell)
    expect_false(trxDesignAlignment(blank)$valid)
    expect_match(trxDesignAlignment(blank)$errors, "outside the declared levels", fixed = TRUE)
    expect_identical(trxDesignAlignment(blank)$errors, trxConditionText(tryCatch(trxAssayInputs(qcCoreBlankGroup(fixture$object)),
        trx_error = identity)))
    absent <- shell
    extra <- absent@design_matrix[1L, , drop = FALSE]
    extra[[absent@sample_id]] <- "S_absent"
    absent@design_matrix <- rbind(absent@design_matrix, extra)
    expect_match(trxDesignAlignment(absent)$errors, "S_absent", fixed = TRUE)

    # The finalise: today's trx checks (on the metadata) give the same passes and texts; the alignment is the parity.
    qc <- qcCoreConfig("transcriptomics")
    today <- trxQcFinalizeCollaborators()$validateFinalizeStateFn
    metadata <- readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, fixture$owner)
    expect_identical(today(metadata), metadata)
    expect_identical(validateWorkflowOmicsQcFinalize(shell, qc, transcriptomicsDeclarationBinding), shell)
    zero_rows <- shell
    attr(zero_rows, "assay_rows")[[1L]] <- 0L
    metadata_zero <- metadata
    metadata_zero$assay_rows[[1L]] <- 0L
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(zero_rows, qc, transcriptomicsDeclarationBinding)),
        qcCoreError(today(metadata_zero)))
    expect_identical(qcCoreError(validateWorkflowOmicsQcFinalize(blank, qc, transcriptomicsDeclarationBinding)),
        paste("Current TranscriptomicsAssayData object has assay/design mismatch:", trxDesignAlignment(blank)$errors))
})

# ------------------------------------------------------------------------------------------------
# The QC action (D6).
# ------------------------------------------------------------------------------------------------

test_that("runWorkflowOmicsQcAction in process: today's reference step, then the figure from the state read for it", {
    workflow <- module_ci_lipid_qc_workflow()
    twin <- module_ci_lipid_qc_workflow()
    step <- omicsIntensityFilterStep(lipidomicsDeclarationBinding, 40L, 0.5)
    seen <- new.env(parent = emptyenv())
    runWorkflowOmicsQcAction(lipidomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$result <- result
        seen$history_at_delivery <- workflow$state_manager$getHistory()
        seen$figure <- figure()
    }, draw = function(arguments) arguments, failed = function(error) stop("not called in process"))
    today <- runLipidReferenceIntensityFilter(twin$state_manager, 40L, 0.5, twin$config_list)
    # The value's `final` is the owner's own product reference; the counts are the step's.
    expect_identical(names(seen$result$value), names(today$value))
    expect_identical(seen$result$value$counts, today$value$counts)
    expect_identical(workflow$state_manager$getHistory(), twin$state_manager$getHistory())
    expect_identical(seen$history_at_delivery, twin$state_manager$getHistory())
    expect_identical(seen$result$reference$entity_id, getWorkflowStateReference(workflow$state_manager)$entity_id)
    expect_named(seen$figure, "theObject")
    expect_identical(seen$figure$theObject@lipid_data, readLipidWorkflowState(twin$state_manager, today$reference)@lipid_data)
})

test_that("runWorkflowOmicsQcAction on a runner: the managed commit, then the figure from the worker's metrics", {
    seen <- new.env(parent = emptyenv())
    metrics <- new.env(parent = emptyenv())
    metrics$value <- list(metric = 1)
    runner <- list(run = function(actions, done = NULL, failed = NULL, dropped = NULL) {
        seen$actions <- actions
        seen$failed <- failed
        actions[[1L]]$deliver(list(value = list(counts = "counts"), stateReference = "reference"))
    })
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) {
        seen$omic_type <- omic_type
        runner
    }, readWorkflowOmicsFilteringMetrics = function(manager, result) metrics$value, .package = "MultiScholaR")
    workflow <- list2env(list(state_manager = "owner", config_list = list(config = TRUE)))
    step <- omicsDuplicateResolutionStep(metabolomicsDeclarationBinding)
    failed <- function(error) "failed"
    runWorkflowOmicsQcAction(metabolomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$result <- result
        seen$figure <- tryCatch(figure(), error = function(error) conditionMessage(error))
    }, draw = function(arguments) arguments, failed = failed)
    expect_identical(seen$omic_type, "metabolomics")
    expect_length(seen$actions, 1L)
    expected <- omicsManagedCommit(step, workflow$config_list, NULL)
    expect_identical(seen$actions[[1L]][c("operation", "parameters", "before", "context")], expected[c("operation", "parameters", "before", "context")])
    expect_identical(seen$failed, failed)
    expect_identical(seen$result, list(value = list(counts = "counts"), reference = "reference"))
    expect_identical(seen$figure, list(theObject = NULL, metrics = list(metric = 1)))
    # No metrics from the worker: the figure fails with today's message (the tab logs it and shows none).
    metrics$value <- NULL
    runWorkflowOmicsQcAction(metabolomicsDeclarationBinding, workflow, step, deliver = function(result, figure) {
        seen$figure <- tryCatch(figure(), error = function(error) conditionMessage(error))
    }, draw = function(arguments) arguments, failed = failed)
    expect_identical(seen$figure, "the worker computed no filtering metrics")
})

# ------------------------------------------------------------------------------------------------
# Trx's QC figure functions, moved with their names.
# ------------------------------------------------------------------------------------------------

test_that("trx's QC figure functions are in func_trx_qc_plots.R, text-identical to any module copy left", {
    definitions <- function(file) {
        exprs <- parse(testthat::test_path("..", "..", "R", file), keep.source = FALSE)
        found <- list()
        for (e in exprs) if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.name(e[[2L]])) found[[as.character(e[[2L]])]] <- e[[3L]]
        found
    }
    moved <- c("trxQcRawStats", "trxQcFilterText", "trxQcMetricsTable", "trxQcFilteringFigure")
    core <- definitions("func_trx_qc_plots.R")
    expect_true(all(moved %in% names(core)))
    module_file <- testthat::test_path("..", "..", "R", "mod_trx_qc_filter.R")
    module <- if (file.exists(module_file)) definitions("mod_trx_qc_filter.R") else list()
    for (name in intersect(moved, names(module))) expect_identical(core[[name]], module[[name]], info = name)
    counts <- list(before = 10L, after_filter_by_expr = 8L, after_zero_variance = 7L)
    expect_match(trxQcFilterText(counts), "After the Zero-variance step: 7", fixed = TRUE)
    metrics <- list(library_size = c(R1 = 100, R2 = 200), detected_genes = c(R1 = 5L, R2 = 6L))
    expect_identical(names(trxQcMetricsTable(metrics)), c("Run", "Library size", "Genes detected"))
    expect_s3_class(trxQcFilteringFigure(counts, metrics), "ggplot")
})
