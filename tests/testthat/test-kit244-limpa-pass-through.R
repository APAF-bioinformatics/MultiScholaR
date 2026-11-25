# KIT-244 (VAL-220 (1); KIT-225 implementation audit I01, I02): DIA with limpa, "Run Normalization & RUV".
# (A) Pressed again while the correlation result is current (after the filter, or after RUV set to Skip, which commits
#     no state), the pass-through keeps the starting table it recorded before, as it already did at the RUV states
#     (GUI audit P-35). A re-applied filter then reads that table, a Skip commits it, and Reset returns it.
# (B) Pressed on a peptide table (before DPC Quantification) the press refuses in plain words, before anything is
#     recorded, drawn or queued (class prot_norm_limpa_peptide_input; the text is a draft, VAL-41).
# On protNormW1LimpaOwner's recipe (helper-prot-norm-core.R:53-86; the DPC result computed once per file and committed
# into each fresh owner), in process and through the worker (the correlation's pairs and selection and the RUV
# covariates through protNormW2Dispatch; the pass-through, Skip, RUV Skip and Reset in process on both paths, as
# test-prot-norm-core.R's W2 runs them). The builders are the G0 script's
# (docs/qa/workflow-execution/kit/kit-244/g0-base-values.R; run 45583972146). Tables are named by the state that
# committed them; "<label>-state" is the state a correlation committed, "<label>" its result.

kit244Quiet <- function(code) suppressWarnings(suppressMessages(code))
kit244Noop <- function(...) invisible(NULL)
kit244Aesthetics <- function() list(color_var = "group", shape_var = NULL)
kit244Settings <- function(mode = "covariates") list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
kit244Input <- function(mode = "covariates") list(norm_method = "quantile", ruv_mode = mode, ruv_percentage = 50,
    ruv_k = 1, color_variable = "group", shape_variable = NULL)
kit244Paths <- function() c("in process", if (workflowOwnerCapable()) "worker")

#' The refusal's text (ticket section 7; the owner sets the wording, VAL-41).
kit244PeptideText <- paste0("on DIA with limpa, \"Run Normalization & RUV\" needs the protein table that DPC ",
    "Quantification makes, and the current table is still a peptide table. In the Quality Control tab apply the ",
    "remaining peptide steps and DPC Quantification (normally followed by Between-Sample Normalisation and the Peptide ",
    "Count Filter), then press \"Run Normalization & RUV\" again.")
kit244ResetText <- paste0("[warning 5] Normalisation reset: the current table is again the one the normalisation ",
    "started from (Peptide Count Filter).")

# ---- fixtures ----

kit244Cache <- new.env(parent = emptyenv())
kit244Objects <- function() {
    if (is.null(kit244Cache$rolled)) {
        kit244Cache$precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        kit244Cache$rolled <- rollUpPrecursorToPeptide(kit244Cache$precursors, core_utilisation = NA)
    }
    kit244Cache
}
kit244Quantified <- function() {
    cache <- kit244Objects()
    if (is.null(cache$quantified)) cache$quantified <- kit244Quiet(proteinMissingValueImputationLimpa(
        log2TransformPeptideMatrix(calcPeptideMatrix(cache$rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000,
        engine = "limpa"))
    cache$quantified
}

kit244HeadRef <- function(fixture) withDIAStateScope(fixture$manager, function() getDIAStateReference(fixture$manager))
kit244Name <- function(fixture, label, reference) {
    assign(label, reference, envir = fixture$refs)
    invisible(reference)
}
kit244Label <- function(fixture, entity_id) {
    if (is.null(entity_id)) return("<null>")
    labels <- sort(ls(fixture$refs))
    hit <- labels[vapply(labels, function(label) identical(get(label, envir = fixture$refs)$entity_id, entity_id), NA)]
    if (length(hit)) paste(hit, collapse = "|") else paste0("other:", substr(entity_id, 1L, 8L))
}
kit244Known <- function(fixture, entity_id) !startsWith(kit244Label(fixture, entity_id), "other:")
kit244HeadLabel <- function(fixture) kit244Label(fixture, kit244HeadRef(fixture)$entity_id)
kit244Recorded <- function(fixture, norm) kit244Label(fixture, norm$normalization_input_reference$entity_id)

kit244Commit <- function(fixture, state, object) {
    kit244Quiet(fixture$manager$saveState(state, object, list(), "kit244"))
    kit244Name(fixture, state, kit244HeadRef(fixture))
}

#' A DIA with limpa owner. `peptide = FALSE`: protNormW1LimpaOwner's commits (protein_dpc_quant, protein_normalised,
#' peptide_count_filtered). `peptide = TRUE`: kit225LimpaOwner's peptide states first. Stops at `upto`.
kit244LimpaOwner <- function(peptide = FALSE, upto = "peptide_count_filtered", envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "kit244-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager,
        refs = new.env(parent = emptyenv()), seen = new.env(parent = emptyenv()))
    cache <- kit244Objects()
    states <- if (peptide) c("raw_data_s4", "qvalue_filtered", "precursor_rollup", "log2_matrix",
        "peptide_accession_cleaned") else character()
    for (state in states) {
        object <- if (state %in% c("raw_data_s4", "qvalue_filtered")) cache$precursors else cache$rolled
        kit244Commit(fixture, state, object)
        if (identical(state, upto)) return(fixture)
    }
    setDIAStepProvenance(manager, "design_builder")
    kit244Commit(fixture, "protein_dpc_quant", kit244Quantified())
    if (identical(upto, "protein_dpc_quant")) return(fixture)
    kit244Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        kit244Name(fixture, manager$current_state, getDIAStateReference(manager))
        if (identical(upto, "protein_normalised")) return(invisible(NULL))
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        kit244Name(fixture, manager$current_state, getDIAStateReference(manager))
    }))
    fixture
}

#' A native proteomics owner of another route (kit225Owner, test-kit225-prot-norm-reset-source.R:20-46), copied.
kit244Owner <- function(type = "LFQ", envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$design_matrix <- original@design_matrix
    workflow$config_list <- original@args
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    manager <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$setWorkflowType(type)
    list(root = root, paths = paths, workflow = workflow, manager = manager, original = original,
        refs = new.env(parent = emptyenv()), seen = new.env(parent = emptyenv()))
}

#' protNormW2Mocks (test-prot-norm-core.R:1217-1221, file-local there), copied.
kit244W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The (A) owner on a path; on the worker path the W2 mocks in the test's frame and a controller.
kit244Start <- function(path, envir = parent.frame()) {
    fixture <- kit244LimpaOwner(envir = envir)
    if (identical(path, "worker")) {
        kit244W2Mocks(env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
    }
    fixture
}

# ---- reading a table (protNormW1State, test-prot-norm-core.R:361-371) ----

kit244Value <- function(fixture, reference) {
    manager <- fixture$manager
    state <- withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        list(slots = lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot)),
            args = object@args[intersect(c("limpa_dpc_quant_results", "ruvCovariates"), names(object@args))],
            samples = setdiff(colnames(object@protein_quant_table), object@protein_id_column))
    })
    list(samples = state$samples, sha = digest::digest(state[c("slots", "args")], algo = "sha256"))
}
kit244Components <- function(fixture, reference) withDIAStateScope(fixture$manager, function()
    readDIAStateMetadata(fixture$manager, reference)$data$components)

#' A failure as "class: message", NA when the code completed.
kit244Error <- function(code) tryCatch({
    kit244Quiet(code)
    NA_character_
}, error = function(e) paste0(class(e)[[1L]], ": ", conditionMessage(e)))

# ---- the tab's actions ----

#' The pass-through (in process on both paths, as protNormW2Clicks$pass_through); the tables its figures were drawn from.
kit244Press <- function(fixture, norm) {
    drawn <- character()
    post <- function(reference) {
        drawn <<- c(drawn, kit244Label(fixture, reference$entity_id))
        invisible(NULL)
    }
    error <- kit244Error(runDIAReferenceLimpaNormalisationPassThrough(kit244Settings(), fixture$workflow, norm, post))
    list(error = error, figures = drawn)
}

#' The managed figure source (dispatchProtNormPassThrough's).
kit244ManagedSource <- function(fixture, norm) kit244Label(fixture, withDIAStateScope(fixture$manager, function()
    diaLimpaPassThroughSource(fixture$manager, norm))$entity_id)

kit244Dispatch <- function(fixture) function(operation, parameters, input, state_name = NULL) {
    context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
    fixture$seen$inputs <- c(fixture$seen$inputs, paste0(operation, ":", kit244Label(fixture, input$entity_id)))
    protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
}

#' The managed correlation click (protNormW2Correlation, test-prot-norm-core.R:1135-1149).
kit244WorkerCorrelation <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm)
        readDIAStateMetadata(manager, source)
        source
    })
    dispatch <- kit244Dispatch(fixture)
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        kit244Aesthetics))
}

#' Correlation Apply at `threshold`: its source, the worker's inputs, whether it committed, its result.
kit244Correlate <- function(fixture, norm, path, threshold, label) {
    manager <- fixture$manager
    before <- kit244HeadRef(fixture)$entity_id
    history <- diaStateHistory(manager)
    source <- tryCatch(withDIAStateScope(manager, function() diaCorrelationSource(manager, norm)), error = function(e) NULL)
    fixture$seen$inputs <- character()
    final <- NULL
    error <- kit244Error(final <- if (identical(path, "worker")) kit244WorkerCorrelation(fixture, norm, threshold) else
        runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
            function() "group", kit244Aesthetics))
    head <- kit244HeadRef(fixture)
    out <- list(error = error, source = kit244Label(fixture, source$entity_id), inputs = fixture$seen$inputs,
        committed = !identical(before, head$entity_id), history_changed = !identical(history, diaStateHistory(manager)))
    if (is.null(final)) return(out)
    if (!kit244Known(fixture, final$entity_id)) kit244Name(fixture, label, final)
    if (!kit244Known(fixture, head$entity_id)) kit244Name(fixture, paste0(label, "-state"), head)
    c(out, list(head = kit244Label(fixture, head$entity_id)), kit244Value(fixture, final))
}

#' Correlation Skip (in process on both paths).
kit244Skip <- function(fixture, norm) {
    error <- kit244Error(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0, TRUE, "proteomics", "core",
        function() NULL, kit244Aesthetics))
    head <- kit244HeadRef(fixture)
    config <- withDIAStateScope(fixture$manager, function() readDIAStateMetadata(fixture$manager)$config)
    c(list(error = error, state = fixture$manager$current_state, head = head, config = config,
        components = kit244Components(fixture, head)), kit244Value(fixture, head))
}

#' The RUV covariates (in process, or through the worker as protNormW2RuvCovariates, test-prot-norm-core.R:1153-1166).
kit244RuvCovariates <- function(fixture, norm, path, mode) {
    manager <- fixture$manager
    settings <- list(ruv_mode = mode, color_variable = "group", shape_variable = NULL)
    if (!identical(path, "worker")) return(kit244Error(runDIAReferenceRuvCovariates(settings, fixture$workflow, norm,
        fixture$paths, "proteomics", "core", kit244Noop, function() "group")))
    kit244Error({
        request <- withDIAStateScope(manager, function() {
            assertDIAExportReady(manager, allow_action = "ruv-covariates")
            list(source = diaRuvCovariatesSource(manager, norm),
                parameters = resolveDIARuvCovariatesParameters(fixture$workflow$config_list, "group", mode))
        })
        result <- kit244Dispatch(fixture)("protein_ruv_covariates", request$parameters, request$source,
            if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
        withDIAStateScope(manager, function() publishDIAReferenceRuvCovariates(manager, settings,
            result[c("value", "reference")], fixture$workflow, norm, fixture$paths, kit244Noop))
    })
}

#' RUV covariates set to Skip through the tab's observer, with the flag the tab sets after the correlation (R1-03).
kit244RuvSkip <- function(fixture, norm) {
    notes <- character()
    note <- function(ui, ...) {
        notes <<- c(notes, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }
    norm$correlation_filtering_complete <- TRUE
    kit244Quiet(runProtNormRuvCovariatesObserver(list(ruv_mode = "skip", color_variable = "group", shape_variable = NULL),
        new.env(), fixture$workflow, norm, fixture$paths, "proteomics", "core", kit244Noop, function() "group",
        renderTextFn = function(x) x, showNotificationFn = note, messageFn = kit244Noop))
    notes
}

#' The tab's Reset (kit225Reset, test-kit225-prot-norm-reset-source.R:77-87), its notices recorded.
kit244Reset <- function(fixture, norm) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    error <- kit244Error(withDIAStateScope(fixture$manager, function() runProtNormResetWorkflow(fixture$workflow, norm,
        output = new.env(), resetOutputsFn = function(...) NULL, showNotificationFn = note,
        messageFn = function(...) NULL)))
    list(error = error, notes = seen$notes, state = fixture$manager$current_state, head = kit244HeadLabel(fixture),
        history = diaStateHistory(fixture$manager))
}

#' The press through runProtNormNormalizationWorkflow: in process, or managed with a recording runner.
kit244FullPress <- function(fixture, norm, managed = FALSE, mode = "covariates") {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$classes <- character()
    seen$runs <- 0L
    seen$figures <- character()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    handle <- function(error, ...) {
        seen$classes <- c(seen$classes, class(error)[[1L]])
        handleProtNormNormalizationError(error, showNotificationFn = note, messageFn = kit244Noop)
    }
    runner <- list(run = function(fn, failed = NULL) {
        seen$runs <- seen$runs + 1L
        invisible(TRUE)
    })
    figure <- function(reference) {
        seen$figures <- c(seen$figures, kit244Label(fixture, reference$entity_id))
        invisible(NULL)
    }
    press <- function() shiny::isolate(runProtNormNormalizationWorkflow(kit244Input(mode), fixture$workflow, norm,
        fixture$paths, "proteomics", "core", generatePostNormalizationQcFn = figure, generateRuvCorrectedQcFn = kit244Noop,
        getRuvGroupingVariableFn = function() "group", showNotificationFn = note, handleErrorFn = handle,
        messageFn = kit244Noop))
    kit244Quiet(if (managed) testthat::with_mocked_bindings(press(),
        workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR") else press())
    as.list(seen)
}

#' t_hi by the rule fixed in G0 (ticket section 10): from the pairs product of correlation 0.5 on the Peptide Count
#' Filter's table, b(s) is a sample's best within-group partner correlation (missing pairs ignored); m is the smallest
#' b(s) among samples with b(s) >= 0.5, m2 the next larger value; t_hi = (m + m2) / 2. NULL when there is no m2.
#' Also kept: the samples correlation 0.5 keeps after the first press alone, and the samples t_hi removes.
kit244ThresholdHigh <- function() {
    if (!is.null(kit244Cache$t_hi)) return(kit244Cache$t_hi)
    fixture <- kit244Start("in process", envir = environment())
    norm <- new.env(parent = emptyenv())
    kit244Press(fixture, norm)
    first <- kit244Correlate(fixture, norm, "in process", 0.5, "C1")
    hydrated <- hydrateSavedWorkflowEntity(fixture$manager$.context$store, norm$correlation_vector)
    pairs <- as.data.frame(if (is.data.frame(hydrated)) hydrated else hydrated$value)
    x <- grep("[.]x$", names(pairs), value = TRUE)[[1L]]
    y <- sub("[.]x$", ".y", x)
    in_pairs <- unique(c(as.character(pairs[[x]]), as.character(pairs[[y]])))
    best <- vapply(in_pairs, function(s) {
        r <- pairs$pearson_correlation[pairs[[x]] == s | pairs[[y]] == s]
        r <- r[!is.na(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
    candidates <- best[!is.na(best) & best >= 0.5]
    if (!length(candidates)) return(NULL)
    m <- min(candidates)
    above <- best[!is.na(best) & best > m]
    if (!length(above)) return(NULL)
    kit244Cache$samples_05 <- first$samples
    kit244Cache$removed_hi <- names(best)[!is.na(best) & best == m]
    kit244Cache$t_hi <- (m + min(above)) / 2
    kit244Cache$t_hi
}

kit244NeedLimpa <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
}

# ---- (A) ----

test_that("KIT-244 E1, G1: pressed again at the correlation result, the recorded start stays the Peptide Count Filter's table (fails before)", {
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        expect_identical(kit244Press(fixture, norm)$error, NA_character_, info = path)
        # G1: the first press designates the QC result.
        expect_identical(kit244Recorded(fixture, norm), "peptide_count_filtered", info = path)
        first <- kit244Correlate(fixture, norm, path, 0.5, "C1")
        expect_identical(first$error, NA_character_, info = path)
        expect_identical(fixture$manager$current_state, "correlation_filtered", info = path)
        press <- kit244Press(fixture, norm)
        expect_identical(press$error, NA_character_, info = path)
        expect_identical(kit244Recorded(fixture, norm), "peptide_count_filtered", info = path)
        expect_identical(press$figures, "peptide_count_filtered", info = path)
        # The managed figure source agrees with the pass-through.
        expect_identical(kit244ManagedSource(fixture, norm), "peptide_count_filtered", info = path)
    }
})

test_that("KIT-244 E1b: after RUV set to Skip (no state committed) the second press keeps the recorded start (fails before)", {
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        kit244Correlate(fixture, norm, path, 0.5, "C1")
        head <- kit244HeadLabel(fixture)
        history <- diaStateHistory(fixture$manager)
        notes <- kit244RuvSkip(fixture, norm)
        expect_identical(notes, "RUV-III covariates skipped. Ready for differential expression analysis.", info = path)
        expect_identical(fixture$manager$current_state, "correlation_filtered", info = path)
        expect_identical(kit244HeadLabel(fixture), head, info = path)
        expect_identical(diaStateHistory(fixture$manager), history, info = path)
        expect_identical(kit244Press(fixture, norm)$error, NA_character_, info = path)
        expect_identical(kit244Recorded(fixture, norm), "peptide_count_filtered", info = path)
        expect_identical(kit244ManagedSource(fixture, norm), "peptide_count_filtered", info = path)
    }
})

test_that("KIT-244 E2: a filter at t_hi, the second press, then 0.5 reads the recorded start and keeps the samples t_hi removed (fails before)", {
    kit244NeedLimpa()
    t_hi <- kit244ThresholdHigh()
    skip_if(is.null(t_hi), "no t_hi on this fixture (stop rule 6)")
    expect_length(kit244Cache$removed_hi, 1L)
    expect_true(all(kit244Cache$removed_hi %in% kit244Cache$samples_05))
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        high <- kit244Correlate(fixture, norm, path, t_hi, "Chi")
        expect_false(any(kit244Cache$removed_hi %in% high$samples), info = path)
        kit244Press(fixture, norm)
        again <- kit244Correlate(fixture, norm, path, 0.5, "R2")
        expect_identical(again$error, NA_character_, info = path)
        expect_identical(again$source, "peptide_count_filtered", info = path)
        if (identical(path, "worker")) expect_identical(again$inputs,
            c("pearson-pairs:peptide_count_filtered", "pearson-selection:peptide_count_filtered"), info = path)
        expect_identical(again$samples, kit244Cache$samples_05, info = path)
    }
})

test_that("KIT-244 G row (E2b): the filter at 0.5 again after the second press gives the first value and commits nothing (passes before: G0)", {
    # Kept as a G row: G0 (run 45583972146) observed it on base, where the re-apply read the committed correlation state.
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        first <- kit244Correlate(fixture, norm, path, 0.5, "C1")
        kit244Press(fixture, norm)
        head <- kit244HeadLabel(fixture)
        again <- kit244Correlate(fixture, norm, path, 0.5, "R2")
        expect_identical(again$error, NA_character_, info = path)
        expect_identical(again$sha, first$sha, info = path)
        expect_identical(again$samples, first$samples, info = path)
        expect_false(again$committed, info = path)
        expect_false(again$history_changed, info = path)
        expect_identical(kit244HeadLabel(fixture), head, info = path)
    }
})

test_that("KIT-244 E3, E3b: a filter at t_hi, the second press, then Skip commits the recorded start; Reset returns it (fails before)", {
    kit244NeedLimpa()
    t_hi <- kit244ThresholdHigh()
    skip_if(is.null(t_hi), "no t_hi on this fixture (stop rule 6)")
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        l7 <- fixture$refs$peptide_count_filtered
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        kit244Correlate(fixture, norm, path, t_hi, "Chi")
        kit244Press(fixture, norm)
        skip <- kit244Skip(fixture, norm)
        expect_identical(skip$error, NA_character_, info = path)
        expect_identical(skip$state, "correlation_filtered", info = path)
        expect_identical(skip$components, kit244Components(fixture, l7), info = path)
        expect_identical(skip$samples, kit244Value(fixture, l7)$samples, info = path)
        expect_true(isTRUE(skip$config$skipped), info = path)
        # E3b: Reset after the Skip.
        reset <- kit244Reset(fixture, norm)
        expect_identical(reset$error, NA_character_, info = path)
        expect_identical(reset$state, "peptide_count_filtered", info = path)
        expect_identical(reset$head, "peptide_count_filtered", info = path)
        expect_false("correlation_filtered" %in% reset$history, info = path)
    }
})

test_that("KIT-244 E4: after the second press, Reset returns the Peptide Count Filter's table and names it (fails before)", {
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        kit244Correlate(fixture, norm, path, 0.5, "C1")
        kit244Press(fixture, norm)
        reset <- kit244Reset(fixture, norm)
        expect_identical(reset$error, NA_character_, info = path)
        expect_identical(reset$state, "peptide_count_filtered", info = path)
        expect_identical(reset$head, "peptide_count_filtered", info = path)
        expect_false("correlation_filtered" %in% reset$history, info = path)
        expect_identical(reset$notes, kit244ResetText, info = path)
    }
})

test_that("KIT-244 E4b: after the second press and a filter at t_hi that commits again, Reset returns the recorded start (fails before)", {
    kit244NeedLimpa()
    t_hi <- kit244ThresholdHigh()
    skip_if(is.null(t_hi), "no t_hi on this fixture (stop rule 6)")
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        kit244Correlate(fixture, norm, path, 0.5, "C1")
        kit244Press(fixture, norm)
        high <- kit244Correlate(fixture, norm, path, t_hi, "C2")
        expect_true(high$committed, info = path)
        reset <- kit244Reset(fixture, norm)
        expect_identical(reset$error, NA_character_, info = path)
        expect_identical(reset$state, "peptide_count_filtered", info = path)
        expect_identical(reset$head, "peptide_count_filtered", info = path)
        expect_false("correlation_filtered" %in% reset$history, info = path)
        expect_identical(reset$notes, kit244ResetText, info = path)
    }
})

test_that("KIT-244 G2b: P-35's sequence (second press at the RUV state, filter again) keeps the recorded start and the first value", {
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        for (mode in c("covariates", "correction")) {
            info <- paste(path, mode)
            fixture <- kit244Start(path)
            norm <- new.env(parent = emptyenv())
            kit244Press(fixture, norm)
            first <- kit244Correlate(fixture, norm, path, 0.5, "C1")
            expect_identical(kit244RuvCovariates(fixture, norm, path, mode), NA_character_, info = info)
            expect_identical(kit244Press(fixture, norm)$error, NA_character_, info = info)
            expect_identical(kit244Recorded(fixture, norm), "peptide_count_filtered", info = info)
            again <- kit244Correlate(fixture, norm, path, 0.5, "R2")
            expect_identical(again$error, NA_character_, info = info)
            expect_identical(again$sha, first$sha, info = info)
            if (identical(path, "worker")) {
                pairs <- unique(grep("^pearson-pairs:", c(first$inputs, again$inputs), value = TRUE))
                expect_identical(pairs, "pearson-pairs:peptide_count_filtered", info = info)
            }
        }
    }
})

test_that("KIT-244 G3, G7: a second press at the QC result designates it; without a recorded start the head is designated", {
    kit244NeedLimpa()
    for (path in kit244Paths()) {
        fixture <- kit244Start(path)
        norm <- new.env(parent = emptyenv())
        kit244Press(fixture, norm)
        kit244Press(fixture, norm)
        expect_identical(kit244Recorded(fixture, norm), "peptide_count_filtered", info = path)
        # G7 (F3): a restored session without the recorded start designates the head, as before.
        kit244Correlate(fixture, norm, path, 0.5, "C1")
        norm$normalization_input_reference <- NULL
        expect_identical(kit244Press(fixture, norm)$error, NA_character_, info = path)
        expect_identical(kit244Recorded(fixture, norm), "C1-state", info = path)
    }
})

# ---- (B) ----

test_that("KIT-244 E5, E5b: pressed on a peptide table, the press refuses in plain words before anything is recorded or queued (fails before)", {
    for (head in c("raw_data_s4", "qvalue_filtered", "precursor_rollup", "log2_matrix", "peptide_accession_cleaned")) {
        for (managed in c(FALSE, TRUE)) {
            info <- paste(head, if (managed) "managed" else "in process")
            fixture <- kit244LimpaOwner(peptide = TRUE, upto = head)
            history <- diaStateHistory(fixture$manager)
            norm <- new.env(parent = emptyenv())
            out <- kit244FullPress(fixture, norm, managed = managed)
            expect_identical(out$classes, "prot_norm_limpa_peptide_input", info = info)
            expect_identical(out$notes, paste0("[error 10] Error in normalization: ", kit244PeptideText), info = info)
            expect_null(norm$normalization_input_reference, info = info)
            expect_identical(out$figures, character(), info = info)
            expect_identical(out$runs, 0L, info = info)
            expect_identical(fixture$manager$current_state, head, info = info)
            expect_identical(diaStateHistory(fixture$manager), history, info = info)
        }
    }
})

test_that("KIT-244 G4: pressed on the protein tables before the Peptide Count Filter, the press designates them", {
    kit244NeedLimpa()
    for (head in c("protein_dpc_quant", "protein_normalised", "peptide_count_filtered")) {
        fixture <- kit244LimpaOwner(upto = head)
        norm <- new.env(parent = emptyenv())
        out <- kit244FullPress(fixture, norm)
        expect_identical(out$classes, character(), info = head)
        expect_identical(kit244Recorded(fixture, norm), head, info = head)
        expect_identical(out$figures, head, info = head)
    }
})

test_that("KIT-244 G5: on DIA (iq) and label-free the limpa guard is not reached", {
    calls <- new.env(parent = emptyenv())
    calls$n <- 0L
    local_mocked_bindings(assertDIALimpaDesignationNotNormalised = function(manager, normData) {
        calls$n <- calls$n + 1L
        invisible(TRUE)
    })
    # DIA (iq) at its protein head (protNormCoreOwner) and at a peptide head.
    core <- protNormCoreOwner()
    core$refs <- new.env(parent = emptyenv())
    out <- kit244FullPress(core, new.env(parent = emptyenv()), mode = "skip")
    expect_false("prot_norm_limpa_peptide_input" %in% out$classes)
    peptide <- kit244Owner("DIA")
    kit244Commit(peptide, "raw_data_s4", kit244Objects()$precursors)
    out <- kit244FullPress(peptide, new.env(parent = emptyenv()), mode = "skip")
    expect_false("prot_norm_limpa_peptide_input" %in% out$classes)
    # Label-free at its QC result (KIT-225's chain without the replicate filter).
    lfq <- kit244Owner("LFQ")
    for (state in c("protein_s4_initial", "protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed"))
        kit244Commit(lfq, state, lfq$original)
    out <- kit244FullPress(lfq, new.env(parent = emptyenv()), mode = "skip")
    expect_identical(out$classes, character())
    expect_identical(calls$n, 0L)
})
