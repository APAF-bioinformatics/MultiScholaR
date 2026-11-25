# KIT-244 gate G0 (r1.2, section 10): base values, recorded with no assertion, of the DIA with limpa
# "Run Normalization & RUV" pass-through at origin/kit-225 600dcc6e (coordinator's ruling). Nothing under R/ is changed.
#
# (A) Sequence S at 0.5 and, with t_hi, the sequences A-hi, A-skip, A-C2; each in process and through the worker (the
#     correlation's pairs and selection through protNormW2Dispatch; the pass-through, Skip, RUV Skip and Reset in
#     process on both paths, as test-prot-norm-core.R's W2 runs them).
# t_hi: fixed here before any value is read (ticket 10, G0): from the 0.5 run's pairs product on L7, b(s) = a sample's
#     best within-group partner correlation (NA pairs ignored; a sample whose pairs are all NA has no b(s));
#     m = the smallest b(s) among samples with b(s) >= 0.5; m2 = the smallest b(s) strictly greater than m;
#     t_hi = (m + m2) / 2; no m2 (or no m): "no t_hi" (stop rule 6).
# (B1) The press through runProtNormNormalizationWorkflow at each of the five peptide heads (fixture states, saveState),
#     in process and managed (a recording runner); the class of raw_data_s4.
# (B2) The peptide states log2_matrix and peptide_accession_cleaned made by the route's own operations (the Quality
#     Control tab's observers) on the owner's rolled object; press, then correlation Apply; fresh copy: press, Skip,
#     RUV covariates.
#
# Usage (from the package root, on a runner): Rscript docs/qa/workflow-execution/kit/kit-244/g0-base-values.R
#   writes docs/qa/workflow-execution/kit/kit-244/g0-out/{g0-base-values,g0-best-partner,g0-pairs}.tsv

kit244Out <- "docs/qa/workflow-execution/kit/kit-244/g0-out"
dir.create(kit244Out, recursive = TRUE, showWarnings = FALSE)
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

kit244Rows <- new.env(parent = emptyenv())
kit244Rows$rows <- list()
kit244Record <- function(case, path, step, value) {
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = ","))
    if (!nzchar(value)) value <- "<none>"
    kit244Rows$rows[[length(kit244Rows$rows) + 1L]] <- data.frame(case = case, path = path, step = step, value = value)
    cat(sprintf("KIT-244 G0 %s %s %s %s\n", case, path, step, value))
}

kit244Quiet <- function(code) suppressWarnings(suppressMessages(code))
kit244Settings <- function(mode = "covariates") list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
kit244Input <- list(norm_method = "quantile", ruv_mode = "covariates", ruv_percentage = 50, ruv_k = 1,
    color_variable = "group", shape_variable = NULL)
kit244Aesthetics <- function() list(color_var = "group", shape_var = NULL)
kit244Noop <- function(...) invisible(NULL)

kit244Outcome <- function(case, path, step, code) {
    outcome <- tryCatch({
        kit244Quiet(code)
        "completed"
    }, error = function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e)))
    kit244Record(case, path, step, outcome)
    outcome
}

kit244Case <- function(case, path, fn) {
    tryCatch(fn(path), error = function(e) kit244Record(case, path, "INSTRUMENT_ERROR",
        paste0("[", paste(class(e), collapse = "/"), "] ", conditionMessage(e))))
    invisible(NULL)
}

# ---- fixtures (protNormW1LimpaOwner, helper-prot-norm-core.R:53-86, and kit225LimpaOwner's peptide states,
# test-kit225-prot-norm-reset-source.R:98-140; the objects computed once and committed into each fresh owner) ----

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

#' Name a reference (labels are kept with their reference; one entity may carry several labels).
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

kit244Commit <- function(fixture, state, object) {
    kit244Quiet(fixture$manager$saveState(state, object, list(), "kit244"))
    kit244Name(fixture, sprintf("%02d-%s", length(ls(fixture$refs)) + 1L, state), kit244HeadRef(fixture))
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
        kit244Name(fixture, sprintf("%02d-%s", length(ls(fixture$refs)) + 1L, manager$current_state),
            getDIAStateReference(manager))
        if (identical(upto, "protein_normalised")) return(invisible(NULL))
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        kit244Name(fixture, sprintf("%02d-%s", length(ls(fixture$refs)) + 1L, manager$current_state),
            getDIAStateReference(manager))
    }))
    fixture
}

#' protNormW2Mocks (test-prot-norm-core.R:1217-1221, file-local there), copied.
kit244W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The (A) owner on a path: L7 named; on the worker path the W2 mocks in the case's frame and a controller.
kit244Start <- function(path, envir = parent.frame()) {
    fixture <- kit244LimpaOwner(envir = envir)
    kit244Name(fixture, "L7", kit244HeadRef(fixture))
    if (identical(path, "worker")) {
        kit244W2Mocks(env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
    }
    fixture
}

# ---- reading a table (protNormW1State, test-prot-norm-core.R:361-371) ----

kit244State <- function(manager, reference) {
    withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        list(slots = lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot)),
            args = object@args[intersect(c("limpa_dpc_quant_results", "ruvCovariates"), names(object@args))],
            samples = if (methods::.hasSlot(object, "protein_quant_table"))
                setdiff(colnames(object@protein_quant_table), object@protein_id_column) else
                as.character(object@design_matrix[[object@sample_id]]))
    })
}
kit244Value <- function(fixture, reference) {
    state <- kit244State(fixture$manager, reference)
    list(samples = state$samples,
        sha = substr(digest::digest(state[c("slots", "args")], algo = "sha256"), 1L, 16L))
}
kit244Components <- function(fixture, reference) withDIAStateScope(fixture$manager, function()
    readDIAStateMetadata(fixture$manager, reference)$data$components)

kit244RecordTable <- function(fixture, case, path, step, reference) {
    value <- tryCatch(kit244Value(fixture, reference), error = function(e) {
        kit244Record(case, path, paste0(step, "_table_read"), paste("error", conditionMessage(e)))
        NULL
    })
    if (is.null(value)) return(invisible(NULL))
    kit244Record(case, path, paste0(step, "_samples"), value$samples)
    kit244Record(case, path, paste0(step, "_n_samples"), length(value$samples))
    kit244Record(case, path, paste0(step, "_value_sha"), value$sha)
    invisible(value)
}

# ---- the tab's actions ----

#' The pass-through (in process on both paths, as protNormW2Clicks$pass_through).
kit244Press <- function(fixture, norm, case, path, step) {
    drawn <- character()
    post <- function(reference) {
        drawn <<- c(drawn, kit244Label(fixture, reference$entity_id))
        invisible(NULL)
    }
    kit244Outcome(case, path, step, runDIAReferenceLimpaNormalisationPassThrough(kit244Settings("covariates"),
        fixture$workflow, norm, post))
    kit244Record(case, path, paste0(step, "_head"), paste(fixture$manager$current_state, kit244HeadLabel(fixture)))
    kit244Record(case, path, paste0(step, "_recorded_start"),
        kit244Label(fixture, norm$normalization_input_reference$entity_id))
    kit244Record(case, path, paste0(step, "_figures_from"), drawn)
    managed <- tryCatch(kit244Label(fixture, withDIAStateScope(fixture$manager, function()
        diaLimpaPassThroughSource(fixture$manager, norm))$entity_id), error = function(e) paste("error", conditionMessage(e)))
    kit244Record(case, path, paste0(step, "_managed_figure_source"), managed)
}

#' The managed correlation click (protNormW2Correlation, test-prot-norm-core.R:1135-1149), the inputs recorded.
kit244WorkerCorrelation <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- diaCorrelationSource(manager, norm)
        readDIAStateMetadata(manager, source)
        source
    })
    dispatch <- function(operation, parameters, input, state_name = NULL) {
        context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
        fixture$seen$inputs <- c(fixture$seen$inputs, paste0(operation, ":", kit244Label(fixture, input$entity_id)))
        protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
    }
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        kit244Aesthetics))
}

#' Correlation Apply at `threshold`: its source, whether it committed, its result (named `label` when new).
kit244Correlate <- function(fixture, norm, case, path, step, threshold, label) {
    manager <- fixture$manager
    before_head <- kit244HeadRef(fixture)$entity_id
    before_history <- diaStateHistory(manager)
    source <- tryCatch(withDIAStateScope(manager, function() diaCorrelationSource(manager, norm)), error = function(e) NULL)
    kit244Record(case, path, paste0(step, "_threshold"), format(threshold, digits = 17))
    kit244Record(case, path, paste0(step, "_source"), kit244Label(fixture, source$entity_id))
    fixture$seen$inputs <- character()
    final <- NULL
    kit244Outcome(case, path, step, final <- if (identical(path, "worker")) kit244WorkerCorrelation(fixture, norm,
        threshold) else runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics",
        "core", function() "group", kit244Aesthetics))
    if (identical(path, "worker")) kit244Record(case, path, paste0(step, "_worker_inputs"), fixture$seen$inputs)
    head <- kit244HeadRef(fixture)
    kit244Record(case, path, paste0(step, "_committed"), !identical(before_head, head$entity_id))
    kit244Record(case, path, paste0(step, "_history_changed"), !identical(before_history, diaStateHistory(manager)))
    if (is.null(final)) return(invisible(NULL))
    if (!kit244Known(fixture, final$entity_id)) kit244Name(fixture, label, final)
    if (!kit244Known(fixture, head$entity_id)) kit244Name(fixture, paste0(label, "-state"), head)
    kit244Record(case, path, paste0(step, "_result"), kit244Label(fixture, final$entity_id))
    kit244Record(case, path, paste0(step, "_head_after"), paste(manager$current_state, kit244Label(fixture, head$entity_id)))
    kit244RecordTable(fixture, case, path, step, final)
    invisible(final)
}

#' Correlation Skip (in process on both paths): what it committed, its config and which named table's components.
kit244Skip <- function(fixture, norm, case, path, step) {
    manager <- fixture$manager
    before <- kit244HeadRef(fixture)$entity_id
    source <- tryCatch(withDIAStateScope(manager, function() diaCorrelationSource(manager, norm)), error = function(e) NULL)
    kit244Record(case, path, paste0(step, "_source"), kit244Label(fixture, source$entity_id))
    kit244Outcome(case, path, step, runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0, TRUE,
        "proteomics", "core", function() NULL, kit244Aesthetics))
    head <- kit244HeadRef(fixture)
    kit244Record(case, path, paste0(step, "_committed"), !identical(before, head$entity_id))
    if (!kit244Known(fixture, head$entity_id)) kit244Name(fixture, "S-state", head)
    kit244Record(case, path, paste0(step, "_head_after"), paste(manager$current_state, kit244Label(fixture, head$entity_id)))
    config <- withDIAStateScope(manager, function() readDIAStateMetadata(manager)$config)
    kit244Record(case, path, paste0(step, "_config"), paste(names(config), unlist(config), sep = "=", collapse = ";"))
    components <- kit244Components(fixture, head)
    same <- Filter(function(label) identical(kit244Components(fixture, get(label, envir = fixture$refs)), components),
        setdiff(sort(ls(fixture$refs)), "S-state"))
    kit244Record(case, path, paste0(step, "_components_as"), same)
    kit244RecordTable(fixture, case, path, step, head)
}

#' RUV covariates set to Skip through the tab's observer (R1-03), with the flag the tab sets after the correlation.
kit244RuvSkip <- function(fixture, norm, case, path, step) {
    manager <- fixture$manager
    notes <- character()
    note <- function(ui, ...) {
        notes <<- c(notes, paste(as.character(ui), collapse = " "))
        invisible(NULL)
    }
    norm$correlation_filtering_complete <- TRUE
    before <- kit244HeadRef(fixture)$entity_id
    before_history <- diaStateHistory(manager)
    kit244Outcome(case, path, step, runProtNormRuvCovariatesObserver(list(ruv_mode = "skip", color_variable = "group",
        shape_variable = NULL), new.env(), fixture$workflow, norm, fixture$paths, "proteomics", "core",
        kit244Noop, function() "group", renderTextFn = function(x) x, showNotificationFn = note, messageFn = kit244Noop))
    kit244Record(case, path, paste0(step, "_notices"), notes)
    kit244Record(case, path, paste0(step, "_committed"), !identical(before, kit244HeadRef(fixture)$entity_id))
    kit244Record(case, path, paste0(step, "_history_changed"), !identical(before_history, diaStateHistory(manager)))
    kit244Record(case, path, paste0(step, "_head_after"), paste(manager$current_state, kit244HeadLabel(fixture)))
    kit244Record(case, path, paste0(step, "_ruv_covariates_complete"), isTRUE(norm$ruv_covariates_complete))
}

#' The tab's Reset (kit225Reset, test-kit225-prot-norm-reset-source.R:77-87), its notices recorded.
kit244Reset <- function(fixture, norm, case, path, step) {
    manager <- fixture$manager
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    kit244Outcome(case, path, step, seen$result <- withDIAStateScope(manager, function() runProtNormResetWorkflow(
        fixture$workflow, norm, output = new.env(), resetOutputsFn = function(...) NULL, showNotificationFn = note,
        messageFn = function(...) NULL)))
    kit244Record(case, path, paste0(step, "_current_state"), manager$current_state)
    kit244Record(case, path, paste0(step, "_head"), kit244HeadLabel(fixture))
    kit244Record(case, path, paste0(step, "_history"), diaStateHistory(manager))
    kit244Record(case, path, paste0(step, "_correlation_filtered_in_history"),
        "correlation_filtered" %in% diaStateHistory(manager))
    previous <- if (is.list(seen$result)) seen$result$previousState else NULL
    kit244Record(case, path, paste0(step, "_previous_state"), previous %||% "<null>")
    kit244Record(case, path, paste0(step, "_notices"), seen$notes)
}

#' The press through runProtNormNormalizationWorkflow (B, E5): in process, or managed with a recording runner.
kit244FullPress <- function(fixture, norm, case, path, step, managed = FALSE) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$errors <- character()
    seen$runs <- 0L
    seen$figures <- character()
    note <- function(ui, type = "default", ...) {
        seen$notes <- c(seen$notes, sprintf("[%s] %s", type, paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    handle <- function(error, ...) {
        seen$errors <- c(seen$errors, paste0("[", paste(class(error), collapse = "/"), "] ", conditionMessage(error)))
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
    press <- function() runProtNormNormalizationWorkflow(kit244Input, fixture$workflow, norm, fixture$paths,
        "proteomics", "core", generatePostNormalizationQcFn = figure, generateRuvCorrectedQcFn = kit244Noop,
        getRuvGroupingVariableFn = function() "group", showNotificationFn = note, handleErrorFn = handle,
        messageFn = kit244Noop)
    kit244Outcome(case, path, step, if (managed) testthat::with_mocked_bindings(press(),
        workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR") else press())
    kit244Record(case, path, paste0(step, "_head"), paste(fixture$manager$current_state, kit244HeadLabel(fixture)))
    kit244Record(case, path, paste0(step, "_handled_errors"), seen$errors)
    kit244Record(case, path, paste0(step, "_notices"), seen$notes)
    kit244Record(case, path, paste0(step, "_recorded_start"),
        kit244Label(fixture, norm$normalization_input_reference$entity_id))
    kit244Record(case, path, paste0(step, "_figures_from"), seen$figures)
    if (managed) kit244Record(case, path, paste0(step, "_runner_actions"), seen$runs)
}

# ---- (A) ----

kit244Prefix <- function(fixture, norm, case, path, threshold, label) {
    kit244Press(fixture, norm, case, path, "press1")
    kit244Correlate(fixture, norm, case, path, "corr1", threshold, label)
    kit244Press(fixture, norm, case, path, "press2")
}

kit244T <- function(path) {
    case <- "T"
    fixture <- kit244Start(path)
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    kit244RecordTable(fixture, case, path, "L7", fixture$refs$L7)
    kit244Press(fixture, norm, case, path, "press1")
    kit244Correlate(fixture, norm, case, path, "corr1", 0.5, "C1")
    kit244Record(case, path, "pairs_reference", kit244Label(fixture, norm$correlation_vector$entity_id))
    hydrated <- hydrateSavedWorkflowEntity(manager$.context$store, norm$correlation_vector)
    kit244Record(case, path, "pairs_hydrated_class", class(hydrated))
    kit244Record(case, path, "pairs_hydrated_names", names(hydrated))
    pairs <- if (is.data.frame(hydrated)) hydrated else if (is.data.frame(hydrated$value)) hydrated$value else
        Filter(is.data.frame, hydrated)[[1L]]
    pairs <- as.data.frame(pairs)
    x <- grep("[.]x$", names(pairs), value = TRUE)[[1L]]
    y <- sub("[.]x$", ".y", x)
    kit244Record(case, path, "pairs_columns", names(pairs))
    kit244Record(case, path, "pairs_rows", nrow(pairs))
    utils::write.table(data.frame(sample_x = pairs[[x]], sample_y = pairs[[y]],
        pearson_correlation = format(pairs$pearson_correlation, digits = 17)),
        file.path(kit244Out, "g0-pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    l7 <- kit244Value(fixture, fixture$refs$L7)$samples
    in_pairs <- unique(c(as.character(pairs[[x]]), as.character(pairs[[y]])))
    best <- vapply(in_pairs, function(s) {
        r <- pairs$pearson_correlation[pairs[[x]] == s | pairs[[y]] == s]
        r <- r[!is.na(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
    utils::write.table(data.frame(sample = in_pairs, best_partner = format(best, digits = 17),
        pairs = vapply(in_pairs, function(s) sum(pairs[[x]] == s | pairs[[y]] == s), integer(1))),
        file.path(kit244Out, "g0-best-partner.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    kit244Record(case, path, "samples_in_no_pair", setdiff(l7, in_pairs))
    kit244Record(case, path, "samples_without_b", in_pairs[is.na(best)])
    kept <- function(t) l7[!(l7 %in% in_pairs) | l7 %in% in_pairs[!is.na(best) & best >= t]]
    observed <- kit244Value(fixture, fixture$refs$C1)$samples
    kit244Record(case, path, "predicted_kept_0.5", kept(0.5))
    kit244Record(case, path, "observed_equals_predicted_0.5", setequal(observed, kept(0.5)))
    candidates <- best[!is.na(best) & best >= 0.5]
    m <- if (length(candidates)) min(candidates) else NA_real_
    above <- best[!is.na(best) & best > m]
    m2 <- if (!is.na(m) && length(above)) min(above) else NA_real_
    kit244Record(case, path, "m", format(m, digits = 17))
    kit244Record(case, path, "m_samples", names(best)[!is.na(best) & best == m])
    kit244Record(case, path, "m2", format(m2, digits = 17))
    if (is.na(m2)) {
        kit244Record(case, path, "t_hi", "no t_hi (stop rule 6)")
        return(invisible(NULL))
    }
    kit244Cache$t_hi <- (m + m2) / 2
    kit244Record(case, path, "t_hi", format(kit244Cache$t_hi, digits = 17))
    kit244Record(case, path, "predicted_kept_t_hi", kept(kit244Cache$t_hi))
    kit244Record(case, path, "predicted_removed_t_hi_not_0.5", setdiff(kept(0.5), kept(kit244Cache$t_hi)))
}

kit244A05 <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A05", path, 0.5, "C1")
    kit244Correlate(fixture, norm, "A05", path, "corr2", 0.5, "R2")
    kit244Record("A05", path, "history_end", diaStateHistory(fixture$manager))
}

kit244A05Skip <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A05-skip", path, 0.5, "C1")
    kit244Skip(fixture, norm, "A05-skip", path, "skip")
    kit244Reset(fixture, norm, "A05-skip", path, "reset")
}

kit244A05Reset <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A05-reset", path, 0.5, "C1")
    kit244Reset(fixture, norm, "A05-reset", path, "reset")
}

kit244A05RuvSkip <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Press(fixture, norm, "A05-ruvskip", path, "press1")
    kit244Correlate(fixture, norm, "A05-ruvskip", path, "corr1", 0.5, "C1")
    kit244RuvSkip(fixture, norm, "A05-ruvskip", path, "ruv_skip")
    kit244Press(fixture, norm, "A05-ruvskip", path, "press2")
}

kit244AHi <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A-hi", path, kit244Cache$t_hi, "Chi")
    kit244Correlate(fixture, norm, "A-hi", path, "corr2", 0.5, "R2")
}

kit244ASkip <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A-skip", path, kit244Cache$t_hi, "Chi")
    kit244Skip(fixture, norm, "A-skip", path, "skip")
    kit244Reset(fixture, norm, "A-skip", path, "reset")
}

kit244AC2 <- function(path) {
    fixture <- kit244Start(path)
    norm <- new.env(parent = emptyenv())
    kit244Prefix(fixture, norm, "A-C2", path, 0.5, "C1")
    kit244Correlate(fixture, norm, "A-C2", path, "corr2", kit244Cache$t_hi, "C2")
    kit244Reset(fixture, norm, "A-C2", path, "reset")
}

# ---- (B) ----

kit244B1 <- function(path) {
    managed <- identical(path, "managed")
    for (head in c("raw_data_s4", "qvalue_filtered", "precursor_rollup", "log2_matrix", "peptide_accession_cleaned")) {
        case <- paste0("B1-", head)
        fixture <- kit244LimpaOwner(peptide = TRUE, upto = head)
        head_class <- withWorkflowStateScope(fixture$manager, function()
            readWorkflowStateMetadata(fixture$manager, kit244HeadRef(fixture))$data$class)
        kit244Record(case, path, "head_class", head_class)
        kit244FullPress(fixture, new.env(parent = emptyenv()), case, path, "press", managed = managed)
    }
}

#' The peptide states made by the Quality Control tab's own observers on the rolled object (B2).
kit244RouteStates <- function(fixture, case, path) {
    manager <- fixture$manager
    log2 <- withDIAStateScope(manager, function() runPeptideLog2ApplyObserver(fixture$workflow, new.env(), kit244Noop,
        "proteomics", "core", showNotificationFn = kit244Noop, removeNotificationFn = kit244Noop, logInfoFn = kit244Noop,
        logErrorFn = kit244Noop, experimentPaths = fixture$paths))
    kit244Record(case, path, "log2_status", paste(log2$status, log2$errorMessage %||% ""))
    kit244Record(case, path, "log2_head", paste(manager$current_state, kit244HeadLabel(fixture)))
    if (!identical(manager$current_state, "log2_matrix")) return(FALSE)
    kit244Name(fixture, "route-log2_matrix", kit244HeadRef(fixture))
    cleaned <- withDIAStateScope(manager, function() runPeptideAccessionApplyObserver(fixture$workflow, ";", "mean",
        new.env(), kit244Noop, "proteomics", "core", showNotificationFn = kit244Noop, removeNotificationFn = kit244Noop,
        logInfoFn = kit244Noop, logErrorFn = kit244Noop, experimentPaths = fixture$paths))
    kit244Record(case, path, "accession_status", paste(cleaned$status, cleaned$errorMessage %||% ""))
    kit244Record(case, path, "accession_head", paste(manager$current_state, kit244HeadLabel(fixture)))
    if (!identical(manager$current_state, "peptide_accession_cleaned")) return(FALSE)
    kit244Name(fixture, "route-peptide_accession_cleaned", kit244HeadRef(fixture))
    kit244Record(case, path, "head_class", withWorkflowStateScope(manager, function()
        readWorkflowStateMetadata(manager, kit244HeadRef(fixture))$data$class))
    TRUE
}

kit244B2 <- function(path) {
    case <- "B2-apply"
    fixture <- kit244LimpaOwner(peptide = TRUE, upto = "precursor_rollup")
    if (!kit244RouteStates(fixture, case, path)) return(kit244Record(case, path, "stop_rule_5", "fixture states only"))
    norm <- new.env(parent = emptyenv())
    kit244FullPress(fixture, norm, case, path, "press")
    kit244Correlate(fixture, norm, case, path, "corr1", 0.5, "P1")
    case <- "B2-skip"
    fixture <- kit244LimpaOwner(peptide = TRUE, upto = "precursor_rollup")
    if (!kit244RouteStates(fixture, case, path)) return(kit244Record(case, path, "stop_rule_5", "fixture states only"))
    norm <- new.env(parent = emptyenv())
    kit244FullPress(fixture, norm, case, path, "press")
    kit244Skip(fixture, norm, case, path, "skip")
    kit244Outcome(case, path, "ruv_covariates", runDIAReferenceRuvCovariates(list(ruv_mode = "covariates",
        color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
        kit244Noop, function() "group"))
    kit244Record(case, path, "ruv_covariates_head", paste(fixture$manager$current_state, kit244HeadLabel(fixture)))
}

# ---- run ----

kit244Record("env", "-", "package_commit", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    error = function(e) "unknown"))
kit244Record("env", "-", "R", R.version.string)
for (package in c("limpa", "mixOmics", "testthat")) kit244Record("env", "-", package,
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "missing")
kit244Record("env", "-", "workflowOwnerCapable", workflowOwnerCapable())

paths_a <- c("in_process", if (workflowOwnerCapable()) "worker")
kit244Case("T", "in_process", kit244T)
for (path in paths_a) {
    kit244Case("A05", path, kit244A05)
    kit244Case("A05-skip", path, kit244A05Skip)
    kit244Case("A05-reset", path, kit244A05Reset)
    kit244Case("A05-ruvskip", path, kit244A05RuvSkip)
    if (!is.null(kit244Cache$t_hi)) {
        kit244Case("A-hi", path, kit244AHi)
        kit244Case("A-skip", path, kit244ASkip)
        kit244Case("A-C2", path, kit244AC2)
    }
}
for (path in c("in_process", "managed")) kit244Case("B1", path, kit244B1)
kit244Case("B2", "in_process", kit244B2)

rows <- do.call(rbind, kit244Rows$rows)
utils::write.table(rows, file.path(kit244Out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("KIT-244 G0: %d rows; instrument errors: %d\n", nrow(rows), sum(rows$step == "INSTRUMENT_ERROR")))
