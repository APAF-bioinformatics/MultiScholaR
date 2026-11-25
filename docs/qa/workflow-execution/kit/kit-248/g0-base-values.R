# KIT-248 gate G0 (ticket r1.2, section 10): base values, recorded with no assertion, of the proteomics correlation
# filter's Skip in the analysis record, the study parameters and a replay, at v0.9-staging 16112917. Observation only:
# branch kit-248-g0 carries this script and nothing under R/; it is never merged.
#
# Owners. Classic routes (workflow types LFQ, TMT, DIA (iq)): the KIT-182 owner's object (the "small" recipe, 64 proteins
# x 12 samples, helper-prot-norm-core.R protNormCoreOwner), committed as protein_replicate_filtered with design_builder
# provenance on a native owner of that workflow type (the record test pattern of test-workflow-analysis-record.R, so the
# record's preconditions are import and design_builder and replay can run); normalisation = the tab's in-process
# "Run Normalization & RUV" (runDIAReferenceNormalization, quantile, RUV skip). DIA with limpa: the limpa owner
# (protNormW1LimpaOwner's commits: protein_dpc_quant, protein_normalised, peptide_count_filtered); normalisation = the
# pass-through (runDIAReferenceLimpaNormalisationPassThrough); after the correlation step, RUV covariates (the tab's
# runner, mode covariates), as the route requires before its DA.
# Apply / Skip: runDIAReferenceCorrelation in process (Apply at t, grouping "group"; Skip = threshold 0, skipped TRUE);
# the managed Apply through the worker's entry point in this process (protNormW2Dispatch with protNormW2Mocks, as
# test-prot-norm-core.R W2 and KIT-244's G0), the Skip always in process (the Skip is controller-only).
# DA: runDIAReferenceDaAnalysis, "~ 0 + group", treated.vs.control, q 0.05, treat_lfc_cutoff 0, publish FALSE.
#
# Per route: SCAN (fresh owner): normalise, then Apply at 0.5, 0.9, 0.95, 0.99 (each Apply reads the same normalised
#   source); the samples each keeps. t = the smallest listed threshold that removes some but not all samples. If none
#   does (stop rule 2), the (A) cases run at t = 0.5 and their numbers rows are record-only. Supplementary (not a
#   ticket observable; chooses nothing): each sample's best within-group partner correlation from the 0.5 Apply.
# Cases, each in a fresh owner: B (normalise, Skip, DA); A (normalise, Apply at t, Skip, DA), in process and with the
#   Apply through the worker; A2 (normalise, Apply at t, normalise again, Skip, DA). (Limpa: RUV covariates before DA.)
# Per case: step_provenance$correlation_filtered; the record's correlation step (operation, parameters, outcomes,
#   preceded_by); validateWorkflowAnalysisRecord's incomplete; the study parameters' correlation lines; the live
#   correlation_filtered table's samples; then a replay into a prepared copy (kitPrepareReplayProject): its outcome
#   (completed, or class and message), the replayed correlation_filtered samples and provenance, and the first
#   contrast's DA table digest against the live one (limpa: the ruv_covariates design digest too).
#
# PREDICTIONS (written before any value was read; from the ticket's section 2 and 10 by reading the code; the observed
# values are compared with these in outcome.md, not here). They are also written to g0-predictions.tsv first.
#
# Usage (run-tests scripts lane, from the package root): Rscript docs/qa/workflow-execution/kit/kit-248/g0-base-values.R <out>

args <- commandArgs(trailingOnly = TRUE)
kit248Out <- if (length(args)) args[[1L]] else "kit248-g0"
dir.create(kit248Out, recursive = TRUE, showWarnings = FALSE)

kit248Predictions <- data.frame(stringsAsFactors = FALSE, rbind(
    c("SCAN", "kept at 0.5 and 0.9", "all 12 samples on the classic routes (fixture: baseline sd 2, noise sd 0.35, so between-sample r near 0.97); limpa not predicted numerically"),
    c("SCAN", "kept at 0.99", "classic: most or all samples removed; limpa not predicted"),
    c("SCAN", "kept at 0.95; t", "not predicted (stop rule 2 may apply)"),
    c("B", "prov_operation", "legacy (no provenance pending at the Skip's commit; takeDIAStepProvenance)"),
    c("B", "record_operation", "legacy; parameters = the state's args with provenance unknown"),
    c("B", "record_incomplete", "contains 'legacy step: correlation_filtered'"),
    c("B", "record_outcomes", "skipped = TRUE"),
    c("B", "study_params_correlation", "<none> (a legacy step is skipped by the writer)"),
    c("B", "live_config", "min_pearson_correlation_threshold=0;skipped=TRUE"),
    c("B", "live_samples", "all 12 (the Skip commits the normalised source)"),
    c("B", "replay_outcome", "error [workflow_replay_unsupported] 'step correlation_filtered cannot be replayed: operation legacy'"),
    c("A", "prov_operation", "pearson-selection (the Apply's entry kept by the keep rule, R/utils_workflow_dia_state.R:838); values threshold = t; preceded_by pearson-pairs"),
    c("A", "record_operation", "pearson-selection; parameters include threshold = t; preceded_by pearson-pairs"),
    c("A", "record_outcomes", "skipped = TRUE"),
    c("A", "record_incomplete", "no 'legacy step: correlation_filtered'"),
    c("A", "study_params_correlation", "a '[n. correlation_filtered: pearson-selection]' block with the pairs parameters and threshold = t; nothing about the skip"),
    c("A", "live_config", "min_pearson_correlation_threshold=0;skipped=TRUE"),
    c("A", "live_samples", "all 12"),
    c("A", "replay_outcome", "completed"),
    c("A", "replay_prov_operation", "pearson-selection"),
    c("A", "replay_samples", "the samples kept at t (fewer than live when t removes a sample)"),
    c("A", "da_equal", "FALSE when t removes a sample; TRUE when it removes none (stop rule 2)"),
    c("A (limpa)", "ruv_design_equal", "FALSE when t removes a sample"),
    c("A worker", "all of A's", "as A (INFERRED: the managed completion writes the request's operation pearson-selection)"),
    c("A2", "all of A's", "as A (the second normalisation commits normalized/ruv_corrected without truncating; the Skip then keeps the pearson-selection entry)")))
names(kit248Predictions) <- c("case", "observable", "predicted")
utils::write.table(kit248Predictions, file.path(kit248Out, "g0-predictions.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
kit248Kit <- new.env(parent = globalenv())
sys.source(file.path("tools", "profiling", "workflow_kit.R"), kit248Kit)
kit248Api <- function(name) get(name, asNamespace("MultiScholaR"))

# ---- recording ----

kit248Rows <- new.env(parent = emptyenv())
kit248Rows$rows <- list()
kit248Text <- function(x) gsub("[\t\r\n]+", " ", gsub("[0-9a-f]{64}", "<id>", x))
kit248Record <- function(case, route, path, step, value) {
    value <- kit248Text(paste(as.character(unlist(value)), collapse = ","))
    if (!length(value) || !nzchar(value)) value <- "<none>"
    kit248Rows$rows[[length(kit248Rows$rows) + 1L]] <- data.frame(case = case, route = route, path = path, step = step,
        value = value)
    cat(sprintf("KIT-248 G0 %s %s %s %s %s\n", case, route, path, step, value))
}
kit248Json <- function(x) {
    if (is.null(x)) return("<null>")
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, force = TRUE))
}
kit248Quiet <- function(code) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(code)))
    value
}
kit248Error <- function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e))
kit248Outcome <- function(case, route, path, step, code) {
    outcome <- tryCatch({
        kit248Quiet(code)
        "completed"
    }, error = kit248Error)
    kit248Record(case, route, path, step, outcome)
    outcome
}
kit248Case <- function(case, route, path, fn) {
    tryCatch(fn(case, route, path), error = function(e) kit248Record(case, route, path, "INSTRUMENT_ERROR", kit248Error(e)))
    invisible(NULL)
}

kit248Noop <- function(...) invisible(NULL)
kit248Aesthetics <- function() list(color_var = "group", shape_var = NULL)
kit248Settings <- function(mode) list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
kit248Contrasts <- data.frame(contrasts = "treated.vs.control", full_format = "treated.vs.control=grouptreated-groupcontrol",
    friendly_names = "treated.vs.control", stringsAsFactors = FALSE)

# ---- owners ----

kit248Cache <- new.env(parent = emptyenv())
kit248Cache$t <- list()
kit248Protein <- function() {
    if (is.null(kit248Cache$protein)) kit248Cache$protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    kit248Cache$protein
}
kit248Quantified <- function() {
    if (is.null(kit248Cache$quantified)) {
        rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")),
            core_utilisation = NA)
        kit248Cache$quantified <- kit248Quiet(proteinMissingValueImputationLimpa(
            log2TransformPeptideMatrix(calcPeptideMatrix(rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000,
            engine = "limpa"))
    }
    kit248Cache$quantified
}

kit248Paths <- function(root) {
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"),
        da_output_dir = file.path(root, "results", "proteomics", "da"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    paths
}

#' A native owner of the route's workflow type, its store at paths$source_dir, the route's start committed.
kit248Owner <- function(route, envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- kit248Paths(root)
    limpa <- identical(route, "DIA_LIMPA")
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, paste0("kit248-", tolower(route)),
        if (limpa) "DIA" else route), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType(route)
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager, route = route,
        seen = new.env(parent = emptyenv()))
    setDIAStepProvenance(manager, "design_builder")
    if (!limpa) {
        original <- kit248Protein()
        workflow$design_matrix <- original@design_matrix
        workflow$config_list <- original@args
        kit248Quiet(manager$saveState("protein_replicate_filtered", original, list(), "fixture"))
        return(fixture)
    }
    kit248Quiet(manager$saveState("protein_dpc_quant", kit248Quantified(), list(), "fixture"))
    kit248Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
    }))
    fixture
}

#' protNormW2Mocks (test-prot-norm-core.R, file-local there), copied.
kit248W2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

# ---- reading ----

kit248SamplesRef <- function(manager, reference) withWorkflowStateScope(manager, function() {
    object <- resolveDIAStateValue(manager$.context, reference)$data
    if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
    if (methods::.hasSlot(object, "protein_quant_table"))
        setdiff(colnames(object@protein_quant_table), object@protein_id_column) else
        as.character(object@design_matrix[[object@sample_id]])
})
kit248Samples <- function(manager, state) kit248SamplesRef(manager, manager$.context$snapshot$states[[state]])
kit248DesignSha <- function(manager, state) withWorkflowStateScope(manager, function() {
    object <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$states[[state]])$data
    if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
    design <- as.data.frame(object@design_matrix)
    rownames(design) <- NULL
    substr(digest::digest(design, algo = "sha256"), 1L, 16L)
})
kit248DaTable <- function(manager) withDIAStateScope(manager, function() {
    manifest <- resolveWorkflowStateValue(manager$.context, manager$.context$snapshot$da_results)
    table <- as.data.frame(readDIAContrastOutput(manifest$results[[1L]], "da_proteins_long")$da_proteins_long)
    key <- intersect(c("comparison", "Protein.Group", "Protein.Ids"), names(table))
    if (length(key)) table <- table[do.call(order, unname(as.list(table[key]))), , drop = FALSE]
    rownames(table) <- NULL
    table[sort(names(table))]
})
kit248Sha <- function(x) substr(digest::digest(x, algo = "sha256"), 1L, 16L)

# ---- the tab's actions ----

kit248Normalise <- function(fixture, norm) {
    if (identical(fixture$route, "DIA_LIMPA")) return(runDIAReferenceLimpaNormalisationPassThrough(
        kit248Settings("covariates"), fixture$workflow, norm, kit248Noop))
    runDIAReferenceNormalization(kit248Settings("skip"), fixture$workflow, norm, fixture$paths, "proteomics", "core",
        kit248Noop, kit248Noop, function() "group")
}

#' The managed correlation click (KIT-244 G0's kit244WorkerCorrelation, from test-prot-norm-core.R protNormW2Correlation).
kit248WorkerApply <- function(fixture, norm, threshold) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        diaCorrelationSource(manager, norm)
    })
    dispatch <- function(operation, parameters, input, state_name = NULL) {
        context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
        protNormW2Dispatch(fixture, fixture$seen$controller, operation, parameters, input, context)
    }
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        kit248Aesthetics))
}

kit248Apply <- function(fixture, norm, threshold, worker = FALSE) {
    if (worker) return(kit248WorkerApply(fixture, norm, threshold))
    runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core",
        function() "group", kit248Aesthetics)
}

kit248Skip <- function(fixture, norm) runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0, TRUE,
    "proteomics", "core", function() NULL, kit248Aesthetics)

kit248RuvCovariates <- function(fixture, norm) runDIAReferenceRuvCovariates(list(ruv_mode = "covariates",
    color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
    kit248Noop, function() "group")

kit248Da <- function(fixture) withDIAStateScope(fixture$manager, function() {
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(fixture$manager)
    da$contrasts_available <- kit248Contrasts$contrasts
    fixture$workflow$contrasts_tbl <- kit248Contrasts
    runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da,
        fixture$workflow, fixture$paths, contrasts = kit248Contrasts, publish = FALSE)
})

# ---- SCAN ----

kit248Scan <- function(case, route, path) {
    fixture <- kit248Owner(route)
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    if (!identical(kit248Outcome(case, route, path, "normalise", kit248Normalise(fixture, norm)), "completed")) return()
    source <- withDIAStateScope(manager, function() diaCorrelationSource(manager, norm))
    all <- kit248SamplesRef(manager, source)
    kit248Record(case, route, path, "source_samples", all)
    kit248Record(case, route, path, "source_n", length(all))
    t <- NA_real_
    for (threshold in c(0.5, 0.9, 0.95, 0.99)) {
        step <- paste0("apply_", threshold)
        final <- NULL
        outcome <- kit248Outcome(case, route, path, step, final <- kit248Apply(fixture, norm, threshold))
        if (!identical(outcome, "completed") || is.null(final)) next
        kept <- kit248SamplesRef(manager, final)
        kit248Record(case, route, path, paste0(step, "_kept_n"), length(kept))
        kit248Record(case, route, path, paste0(step, "_removed"), setdiff(all, kept))
        if (is.na(t) && length(kept) > 0L && length(kept) < length(all)) t <- threshold
        if (identical(threshold, 0.5)) kit248BestPartner(fixture, norm, route)
    }
    kit248Record(case, route, path, "t", if (is.na(t)) "none (stop rule 2): the (A) cases run at 0.5" else t)
    kit248Cache$t[[route]] <- if (is.na(t)) 0.5 else t
}

#' Supplementary: each sample's best within-group partner correlation from the 0.5 Apply's pairs (KIT-244 G0's form).
kit248BestPartner <- function(fixture, norm, route) tryCatch({
    hydrated <- hydrateSavedWorkflowEntity(fixture$manager$.context$store, norm$correlation_vector)
    pairs <- if (is.data.frame(hydrated)) hydrated else if (is.data.frame(hydrated$value)) hydrated$value else
        Filter(is.data.frame, hydrated)[[1L]]
    pairs <- as.data.frame(pairs)
    x <- grep("[.]x$", names(pairs), value = TRUE)[[1L]]
    y <- sub("[.]x$", ".y", x)
    samples <- unique(c(as.character(pairs[[x]]), as.character(pairs[[y]])))
    best <- vapply(samples, function(s) {
        r <- pairs$pearson_correlation[pairs[[x]] == s | pairs[[y]] == s]
        r <- r[!is.na(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
    utils::write.table(data.frame(route = route, sample = samples, best_partner = format(best, digits = 17)),
        file.path(kit248Out, paste0("g0-best-partner-", tolower(route), ".tsv")), sep = "\t", quote = FALSE,
        row.names = FALSE)
}, error = function(e) kit248Record("SCAN", route, "in_process", "best_partner", kit248Error(e)))

# ---- a case's observations ----

kit248Observe <- function(fixture, case, route, path) {
    manager <- fixture$manager
    record <- function(step, value) kit248Record(case, route, path, step, value)
    snapshot <- manager$.context$snapshot
    history <- diaStateHistory(manager)
    record("history", history)
    record("head", manager$current_state)
    record("pending_after", manager$.context$pending_provenance$operation %||% "<null>")
    record("prov_all", vapply(history, function(state) paste0(state, ":",
        snapshot$step_provenance[[state]]$operation %||% "<none>"), ""))
    entry <- snapshot$step_provenance$correlation_filtered
    record("prov_operation", entry$operation %||% "<null>")
    record("prov_values", kit248Json(entry$values))
    record("prov_fixed", kit248Json(entry$fixed))
    record("prov_preceded_by", vapply(entry$preceded_by %||% list(), function(p) p$operation, ""))
    record("prov_committed_at", entry$committed_at %||% "<null>")
    config <- withDIAStateScope(manager, function() readDIAStateMetadata(manager, snapshot$states$correlation_filtered)$config)
    record("live_config", paste(names(config), vapply(config, function(v) paste(format(v), collapse = "|"), ""),
        sep = "=", collapse = ";"))
    live <- kit248Samples(manager, "correlation_filtered")
    record("live_samples", live)
    record("live_n", length(live))
    # The analysis record, written and read back as replay reads it.
    json <- file.path(fixture$root, "analysis_record.json")
    kit248Quiet(writeWorkflowAnalysisRecord(manager, json))
    rec <- jsonlite::read_json(json, simplifyVector = FALSE)
    check <- validateWorkflowAnalysisRecord(rec)
    record("record_valid", check$valid)
    record("record_incomplete", check$incomplete)
    record("record_steps", vapply(rec$steps, function(s) paste0(s$state, ":", s$operation), ""))
    step <- Filter(function(s) identical(s$state, "correlation_filtered"), rec$steps)[[1L]]
    record("record_operation", step$operation)
    record("record_parameters", vapply(step$parameters, function(p) sprintf("%s[%s/%s]=%s", p$name, p$source,
        p$provenance, substr(workflowRecordValueText(p$value), 1L, 60L)), ""))
    record("record_outcomes", kit248Json(step$outcomes))
    record("record_preceded_by", vapply(step$preceded_by %||% list(), function(p) paste0(p$operation, "(",
        paste(vapply(p$parameters, function(q) q$name, ""), collapse = "|"), ")"), ""))
    record("record_committed_at", step$committed_at %||% "<null>")
    # The study parameters' correlation lines.
    study <- file.path(fixture$root, "study-parameters")
    dir.create(study, showWarnings = FALSE)
    kit248Quiet(createWorkflowArgsFromConfig("kit248", source_dir_path = study, final_s4_object = getDIAStateReference(manager),
        workflow_data = list(state_manager = manager, config_list = list())))
    text <- readLines(file.path(study, "study_parameters.txt"))
    record("study_params_has_record", "Analysis Record:" %in% text)
    start <- grep("^\\[[0-9]+\\. correlation_filtered:", text)
    block <- if (length(start)) {
        end <- which(!nzchar(text) & seq_along(text) > start[[1L]])
        text[start[[1L]]:(if (length(end)) end[[1L]] - 1L else length(text))]
    } else "<none>"
    record("study_params_correlation", paste(trimws(block), collapse = " | "))
    record("study_params_threshold_lines", grep("threshold", text, value = TRUE))
    # The live DA.
    da_live <- tryCatch(kit248DaTable(manager), error = function(e) {
        record("live_da_read", kit248Error(e))
        NULL
    })
    if (!is.null(da_live)) {
        record("live_da_rows", nrow(da_live))
        record("live_da_sha", kit248Sha(da_live))
    }
    limpa <- identical(route, "DIA_LIMPA")
    design_live <- if (limpa) tryCatch(kit248DesignSha(manager, "ruv_covariates"), error = kit248Error)
    if (limpa) record("live_ruv_design_sha", design_live)
    # Replay into a prepared copy holding only the record's preconditions.
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(fixture$paths$source_dir, "."), copy)) == 0L)
    states <- vapply(rec$steps, function(s) s$state, "")
    operations <- vapply(rec$steps, function(s) s$operation, "")
    preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
    record("replay_preconditions", preconditions)
    prepared <- kit248Kit$kitPrepareReplayProject(copy, preconditions, kit248Api)
    replay_paths <- kit248Paths(withr::local_tempdir())
    scratch <- withr::local_tempdir()
    replay <- NULL
    outcome <- tryCatch({
        kit248Quiet(replay <- replayWorkflowAnalysisRecord(rec, prepared, character(), replay_paths, scratch))
        "completed"
    }, error = kit248Error)
    record("replay_outcome", outcome)
    if (!identical(outcome, "completed")) return(invisible(NULL))
    record("replay_states", replay$states)
    identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
    reopened <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeDIAStateAdapter(reopened), add = TRUE)
    record("replay_prov_operation", reopened$.context$snapshot$step_provenance$correlation_filtered$operation %||% "<null>")
    replayed <- tryCatch(kit248Samples(reopened, "correlation_filtered"), error = kit248Error)
    record("replay_samples", replayed)
    record("replay_n", length(replayed))
    record("replay_samples_equal_live", setequal(replayed, live) && length(replayed) == length(live))
    record("replay_removed_vs_live", setdiff(live, replayed))
    replay_config <- tryCatch(withDIAStateScope(reopened, function() readDIAStateMetadata(reopened,
        reopened$.context$snapshot$states$correlation_filtered)$config), error = function(e) list(error = kit248Error(e)))
    record("replay_config", paste(names(replay_config), vapply(replay_config, function(v) paste(format(v), collapse = "|"), ""),
        sep = "=", collapse = ";"))
    da_replay <- tryCatch(kit248DaTable(reopened), error = function(e) {
        record("replay_da_read", kit248Error(e))
        NULL
    })
    if (!is.null(da_replay)) {
        record("replay_da_rows", nrow(da_replay))
        record("replay_da_sha", kit248Sha(da_replay))
        if (!is.null(da_live)) {
            same_columns <- setequal(names(da_live), names(da_replay))
            record("da_equal", same_columns && nrow(da_live) == nrow(da_replay) &&
                isTRUE(all.equal(da_live, da_replay[names(da_live)], tolerance = 0, check.attributes = FALSE)))
        }
    }
    if (limpa) {
        design_replay <- tryCatch(kit248DesignSha(reopened, "ruv_covariates"), error = kit248Error)
        record("replay_ruv_design_sha", design_replay)
        record("ruv_design_equal", identical(design_live, design_replay))
    }
    invisible(NULL)
}

# ---- cases ----

kit248Finish <- function(fixture, norm, case, route, path) {
    if (identical(route, "DIA_LIMPA"))
        kit248Outcome(case, route, path, "ruv_covariates", kit248RuvCovariates(fixture, norm))
    kit248Outcome(case, route, path, "da", kit248Da(fixture))
}

kit248B <- function(case, route, path) {
    fixture <- kit248Owner(route)
    norm <- new.env(parent = emptyenv())
    kit248Outcome(case, route, path, "normalise", kit248Normalise(fixture, norm))
    kit248Outcome(case, route, path, "skip", kit248Skip(fixture, norm))
    kit248Finish(fixture, norm, case, route, path)
    kit248Observe(fixture, case, route, path)
}

kit248A <- function(case, route, path) {
    fixture <- kit248Owner(route)
    norm <- new.env(parent = emptyenv())
    t <- kit248Cache$t[[route]]
    kit248Record(case, route, path, "t_used", t)
    live <- function() {
        worker <- identical(path, "worker")
        if (worker) {
            kit248W2Mocks(env = environment())
            fixture$seen$controller <- omicsWorkerController(list(owner = fixture$manager))
        }
        kit248Outcome(case, route, path, "normalise", kit248Normalise(fixture, norm))
        kit248Outcome(case, route, path, "apply", kit248Apply(fixture, norm, t, worker))
        kit248Record(case, route, path, "after_apply_prov",
            fixture$manager$.context$snapshot$step_provenance$correlation_filtered$operation %||% "<null>")
        kit248Record(case, route, path, "after_apply_n", tryCatch(length(kit248Samples(fixture$manager,
            "correlation_filtered")), error = kit248Error))
        if (identical(case, "A2")) {
            kit248Outcome(case, route, path, "normalise_again", kit248Normalise(fixture, norm))
            kit248Record(case, route, path, "after_normalise_again_head", fixture$manager$current_state)
            kit248Record(case, route, path, "after_normalise_again_history", diaStateHistory(fixture$manager))
        }
        kit248Outcome(case, route, path, "skip", kit248Skip(fixture, norm))
        kit248Finish(fixture, norm, case, route, path)
    }
    live()
    kit248Observe(fixture, case, route, path)
}

# ---- run ----

kit248Record("env", "-", "-", "package_commit", tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    error = function(e) "unknown"))
kit248Record("env", "-", "-", "R", R.version.string)
for (package in c("limpa", "mixOmics", "testthat")) kit248Record("env", "-", "-", package,
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "missing")
kit248Record("env", "-", "-", "workflowOwnerCapable", workflowOwnerCapable())

routes <- c("LFQ", "TMT", "DIA", if (requireNamespace("limpa", quietly = TRUE) && requireNamespace("mixOmics", quietly = TRUE))
    "DIA_LIMPA")
for (route in routes) {
    started <- Sys.time()
    kit248Case("SCAN", route, "in_process", kit248Scan)
    if (is.null(kit248Cache$t[[route]])) kit248Cache$t[[route]] <- 0.5
    kit248Case("B", route, "in_process", kit248B)
    kit248Case("A", route, "in_process", kit248A)
    if (workflowOwnerCapable()) kit248Case("A", route, "worker", kit248A)
    kit248Case("A2", route, "in_process", kit248A)
    kit248Record("timing", route, "-", "seconds", round(as.numeric(difftime(Sys.time(), started, units = "secs"))))
}

rows <- do.call(rbind, kit248Rows$rows)
utils::write.table(rows, file.path(kit248Out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("KIT-248 G0: %d rows; instrument errors: %d\n", nrow(rows), sum(rows$step == "INSTRUMENT_ERROR")))
