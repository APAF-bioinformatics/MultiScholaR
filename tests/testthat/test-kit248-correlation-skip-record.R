# KIT-248 (VAL-238; AUDIT-KIT-244-plan P14, U1): the proteomics correlation filter's Skip records itself. When it commits,
# the Skip sets its own step provenance, `pearson-skip` (registered in process, no parameters, the fixed constant
# skipped = TRUE), so the analysis record and the study parameters state the skip, and a replay commits the live source
# for it. Records made before the fix (OD-1 (a), (a')): a `pearson-selection` correlation step or a `legacy` one whose
# outcome says skipped replays as the skip it was.
# Owners as the G0 script builds them (docs/qa/workflow-execution/kit/kit-248/g0-base-values.R; runs 45592410449 and
# 45597013061): the KIT-182 object ("small" recipe, 64 proteins x 12 samples) committed as protein_replicate_filtered with
# design_builder provenance on a native owner of the workflow type (label-free, TMT, DIA (iq)); the DIA with limpa owner
# of protNormW1LimpaOwner (helper-prot-norm-core.R). t_hi is the G0b rule (ticket r1.3, section 10): from the 0.5 Apply's
# pairs, the midpoint of the smallest and second-smallest best within-group partner correlation, which removes exactly
# one sample. The managed Apply runs through the worker's entry point in this process (protNormW2Dispatch); the Skip is
# in process on both paths, as the tab runs it.

kit248Quiet <- function(code) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(code)))
    value
}
kit248Noop <- function(...) invisible(NULL)
kit248Aesthetics <- function() list(color_var = "group", shape_var = NULL)
kit248Settings <- function(mode) list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
kit248Contrasts <- data.frame(contrasts = "treated.vs.control", full_format = "treated.vs.control=grouptreated-groupcontrol",
    friendly_names = "treated.vs.control", stringsAsFactors = FALSE)
#' A failure as "class: message", NA when there was none.
kit248ErrorText <- function(error) if (is.null(error)) NA_character_ else
    paste0(class(error)[[1L]], ": ", conditionMessage(error))
kit248NeedLimpa <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
}

# ---- owners (the G0 script's) ----

kit248Cache <- new.env(parent = emptyenv())
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
kit248Kit <- function() {
    if (is.null(kit248Cache$kit)) {
        kit <- new.env(parent = globalenv())
        sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), kit)
        kit248Cache$kit <- kit
    }
    kit248Cache$kit
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

#' A native owner of the route's workflow type, its store at paths$source_dir, the route's start committed; on the
#' worker path the W2 mocks in `envir` and a controller.
kit248Owner <- function(route, path = "in process", envir = parent.frame()) {
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
    if (identical(path, "worker")) {
        testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
            finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
            .package = "MultiScholaR", .env = envir)
        fixture$seen$controller <- omicsWorkerController(list(owner = manager))
    }
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

# ---- the tab's actions ----

#' "Run Normalization & RUV" (classic: in process, quantile, RUV skipped; DIA with limpa: the pass-through).
kit248Normalise <- function(fixture, norm) kit248Quiet(if (identical(fixture$route, "DIA_LIMPA"))
    runDIAReferenceLimpaNormalisationPassThrough(kit248Settings("covariates"), fixture$workflow, norm, kit248Noop) else
    runDIAReferenceNormalization(kit248Settings("skip"), fixture$workflow, norm, fixture$paths, "proteomics", "core",
        kit248Noop, kit248Noop, function() "group"))

#' The managed correlation click (protNormW2Correlation, test-prot-norm-core.R; KIT-244's kit244WorkerCorrelation).
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

kit248Apply <- function(fixture, norm, threshold, path = "in process") kit248Quiet(if (identical(path, "worker"))
    kit248WorkerApply(fixture, norm, threshold) else runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths,
    threshold, FALSE, "proteomics", "core", function() "group", kit248Aesthetics))

kit248Skip <- function(fixture, norm) kit248Quiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, 0,
    TRUE, "proteomics", "core", function() NULL, kit248Aesthetics))

#' The Skip as the code before KIT-248 committed it (`diaCorrelationCommit` and the commit, no provenance set): builds
#' the records made before the fix (OD-1 (a), (a')) on any tree.
kit248PreFixSkip <- function(fixture, norm) {
    manager <- fixture$manager
    kit248Quiet(withDIAStateScope(manager, function() {
        commit <- diaCorrelationCommit(manager, diaCorrelationSource(manager, norm), 0, TRUE)
        commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description)
    }))
}

kit248RuvCovariates <- function(fixture, norm) kit248Quiet(runDIAReferenceRuvCovariates(list(ruv_mode = "covariates",
    color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core",
    kit248Noop, function() "group"))

kit248Da <- function(fixture) kit248Quiet(withDIAStateScope(fixture$manager, function() {
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(fixture$manager)
    da$contrasts_available <- kit248Contrasts$contrasts
    fixture$workflow$contrasts_tbl <- kit248Contrasts
    runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), da,
        fixture$workflow, fixture$paths, contrasts = kit248Contrasts, publish = FALSE)
}))

# ---- reading ----

#' A committed state's table: its samples, quantified table and design.
kit248Table <- function(manager, state) withWorkflowStateScope(manager, function() {
    object <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$states[[state]])$data
    if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
    design <- as.data.frame(object@design_matrix)
    rownames(design) <- NULL
    list(samples = setdiff(colnames(object@protein_quant_table), object@protein_id_column),
        table = object@protein_quant_table, design = design)
})
kit248Entry <- function(manager, state = "correlation_filtered") manager$.context$snapshot$step_provenance[[state]]
kit248Config <- function(manager) withDIAStateScope(manager, function() readDIAStateMetadata(manager,
    manager$.context$snapshot$states$correlation_filtered)$config)
kit248Components <- function(manager, reference) withDIAStateScope(manager, function()
    readDIAStateMetadata(manager, reference)$data$components)
kit248DaTable <- function(manager) withDIAStateScope(manager, function() {
    manifest <- resolveWorkflowStateValue(manager$.context, manager$.context$snapshot$da_results)
    table <- as.data.frame(readDIAContrastOutput(manifest$results[[1L]], "da_proteins_long")$da_proteins_long)
    key <- intersect(c("comparison", "Protein.Group", "Protein.Ids"), names(table))
    if (length(key)) table <- table[do.call(order, unname(as.list(table[key]))), , drop = FALSE]
    rownames(table) <- NULL
    table[sort(names(table))]
})

#' The analysis record, written and read back as replay reads it.
kit248Record <- function(fixture) {
    json <- file.path(fixture$root, paste0("analysis_record_", length(list.files(fixture$root, "^analysis_record_")), ".json"))
    kit248Quiet(writeWorkflowAnalysisRecord(fixture$manager, json))
    jsonlite::read_json(json, simplifyVector = FALSE)
}
kit248Step <- function(record, state = "correlation_filtered")
    Filter(function(step) identical(step$state, state), record$steps)[[1L]]
kit248ParameterNames <- function(step) vapply(step$parameters, function(p) p$name, "")

#' The study parameters' block for the correlation step (its header line to the blank line), or character().
kit248StudyBlock <- function(fixture) {
    study <- withr::local_tempdir()
    kit248Quiet(createWorkflowArgsFromConfig("kit248", source_dir_path = study,
        final_s4_object = getDIAStateReference(fixture$manager),
        workflow_data = list(state_manager = fixture$manager, config_list = list())))
    text <- readLines(file.path(study, "study_parameters.txt"))
    start <- grep("^\\[[0-9]+\\. correlation_filtered:", text)
    if (!length(start)) return(character())
    end <- which(!nzchar(text) & seq_along(text) > start[[1L]])
    text[start[[1L]]:(if (length(end)) end[[1L]] - 1L else length(text))]
}

#' Replay `record` into a copy of the fixture's store holding only its preconditions; the reopened replay owner.
kit248Replay <- function(fixture, record, envir = parent.frame()) {
    copy <- withr::local_tempdir(.local_envir = envir)
    stopifnot(system2("cp", c("-a", file.path(fixture$paths$source_dir, "."), copy)) == 0L)
    states <- vapply(record$steps, function(step) step$state, "")
    operations <- vapply(record$steps, function(step) step$operation, "")
    preconditions <- states[seq_len(sum(cumprod(operations %in% c("import", "design_builder"))))]
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    prepared <- kit248Kit()$kitPrepareReplayProject(copy, preconditions, api)
    paths <- kit248Paths(withr::local_tempdir(.local_envir = envir))
    scratch <- withr::local_tempdir(.local_envir = envir)
    error <- tryCatch({
        kit248Quiet(replayWorkflowAnalysisRecord(record, prepared, character(), paths, scratch))
        NULL
    }, error = function(e) e)
    if (!is.null(error)) return(list(error = error))
    identity <- readSavedStoreIdentity(file.path(prepared, ".workflow-entities", "store.json"))
    reopened <- newDIAStateAdapter(newWorkflowEntityStore(prepared, identity$project_id, identity$workflow_id), reopen = TRUE)
    withr::defer(closeDIAStateAdapter(reopened), envir = envir)
    list(error = NULL, manager = reopened)
}

#' t_hi by the G0b rule on the route's data (classic routes share the KIT-182 object), with the sample it removes.
kit248THi <- function(route) {
    key <- if (identical(route, "DIA_LIMPA")) "limpa" else "classic"
    if (!is.null(kit248Cache[[key]])) return(kit248Cache[[key]])
    fixture <- kit248Owner(route, envir = environment())
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, 0.5)
    all <- kit248Table(fixture$manager, "correlation_filtered")$samples
    hydrated <- hydrateSavedWorkflowEntity(fixture$manager$.context$store, norm$correlation_vector)
    pairs <- as.data.frame(if (is.data.frame(hydrated)) hydrated else hydrated$value)
    x <- grep("[.]x$", names(pairs), value = TRUE)[[1L]]
    y <- sub("[.]x$", ".y", x)
    samples <- unique(c(as.character(pairs[[x]]), as.character(pairs[[y]])))
    best <- vapply(samples, function(s) {
        r <- pairs$pearson_correlation[pairs[[x]] == s | pairs[[y]] == s]
        r <- r[!is.na(r)]
        if (length(r)) max(r) else NA_real_
    }, numeric(1))
    best <- sort(best[!is.na(best)])
    value <- if (length(best) >= 2L && best[[1L]] < best[[2L]])
        list(t_hi = (best[[1L]] + best[[2L]]) / 2, removed = names(best)[[1L]], all = all) else NULL
    kit248Cache[[key]] <- value
    value
}

#' What E1 asserts on the record and the study parameters of a Skip (fails before: the step is the Apply's or legacy).
kit248ExpectSkipRecorded <- function(fixture, info) {
    expect_identical(kit248Entry(fixture$manager)$operation, "pearson-skip", info = info)
    expect_identical(unlist(kit248Entry(fixture$manager)$fixed), c(skipped = TRUE), info = info)
    expect_length(kit248Entry(fixture$manager)$preceded_by, 0L)
    record <- kit248Record(fixture)
    step <- kit248Step(record)
    expect_identical(step$operation, "pearson-skip", info = info)
    expect_false("threshold" %in% kit248ParameterNames(step), info = info)
    skipped <- Filter(function(p) identical(p$name, "skipped"), step$parameters)
    expect_length(skipped, 1L)
    if (length(skipped)) {
        expect_identical(skipped[[1L]]$source, "fixed", info = info)
        expect_identical(workflowRecordValueText(skipped[[1L]]$value), "TRUE", info = info)
    }
    expect_true(isTRUE(unlist(step$outcomes$skipped)), info = info)
    expect_false("legacy step: correlation_filtered" %in% validateWorkflowAnalysisRecord(record)$incomplete, info = info)
    block <- kit248StudyBlock(fixture)
    expect_match(if (length(block)) block[[1L]] else "<no correlation block>", "^\\[[0-9]+\\. correlation_filtered: pearson-skip\\]$", info = info)
    expect_true("  skipped [fixed] = TRUE" %in% block, info = info)
    expect_false(any(grepl("threshold", block)), info = info)
    record
}

# ---- E rows ----

test_that("KIT-248 E1, E2, E3, E7, G2, G3: on label-free, TMT and DIA (iq), Apply at t_hi then Skip is recorded as the skip and replays on every sample (fails before)", {
    thi <- kit248THi("LFQ")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    for (route in c("LFQ", "TMT", "DIA")) {
        fixture <- kit248Owner(route)
        manager <- fixture$manager
        norm <- new.env(parent = emptyenv())
        kit248Normalise(fixture, norm)
        source <- withDIAStateScope(manager, function() diaCorrelationSource(manager, norm))
        kit248Apply(fixture, norm, thi$t_hi)
        # Instrument: t_hi removes exactly the one sample (G0b).
        expect_identical(setdiff(thi$all, kit248Table(manager, "correlation_filtered")$samples), thi$removed, info = route)
        kit248Skip(fixture, norm)
        # G2: the Skip's table and config are unchanged by the fix.
        expect_equal(kit248Config(manager), list(min_pearson_correlation_threshold = 0, skipped = TRUE), info = route)
        expect_identical(kit248Components(manager, getDIAStateReference(manager)), kit248Components(manager, source), info = route)
        live <- kit248Table(manager, "correlation_filtered")
        expect_identical(live$samples, thi$all, info = route)
        # E1, E7.
        kit248ExpectSkipRecorded(fixture, route)
        # G3 / E3: the DA's argument commit keeps the Skip's entry.
        kit248Da(fixture)
        expect_identical(kit248Entry(manager)$operation, "pearson-skip", info = paste(route, "after DA"))
        da_live <- kit248DaTable(manager)
        # E2, E3: the replay of the record with its DA.
        replay <- kit248Replay(fixture, kit248Record(fixture))
        expect_identical(kit248ErrorText(replay$error), NA_character_, info = route)
        if (is.null(replay$manager)) next
        expect_identical(kit248Entry(replay$manager)$operation, "pearson-skip", info = route)
        replayed <- kit248Table(replay$manager, "correlation_filtered")
        expect_identical(replayed$samples, live$samples, info = route)
        expect_identical(replayed$table, live$table, info = route)
        expect_identical(kit248DaTable(replay$manager), da_live, info = route)
    }
})

test_that("KIT-248 E4, E5: label-free Skip with no Apply before it is recorded as the skip and the record replays (fails before: legacy, refused)", {
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Skip(fixture, norm)
    live <- kit248Table(fixture$manager, "correlation_filtered")
    record <- kit248ExpectSkipRecorded(fixture, "E4")
    expect_false("legacy step: correlation_filtered" %in% validateWorkflowAnalysisRecord(record)$incomplete)
    replay <- kit248Replay(fixture, record)
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    skip_if(is.null(replay$manager), "replay refused")
    expect_identical(kit248Entry(replay$manager)$operation, "pearson-skip")
    replayed <- kit248Table(replay$manager, "correlation_filtered")
    expect_identical(replayed$samples, live$samples)
    expect_identical(replayed$table, live$table)
})

test_that("KIT-248 E8: the Apply through the worker, then Skip, is recorded as the skip (fails before)", {
    skip_if_not(workflowOwnerCapable())
    thi <- kit248THi("LFQ")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("LFQ", "worker")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi, "worker")
    expect_identical(kit248Entry(fixture$manager)$operation, "pearson-selection")
    kit248Skip(fixture, norm)
    kit248ExpectSkipRecorded(fixture, "E8")
})

test_that("KIT-248 E12: Apply at t_hi, Run Normalization & RUV again, Skip is recorded as the skip (fails before)", {
    thi <- kit248THi("LFQ")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi)
    kit248Normalise(fixture, norm)
    expect_identical(fixture$manager$current_state, "normalized")
    kit248Skip(fixture, norm)
    expect_identical(kit248Table(fixture$manager, "correlation_filtered")$samples, thi$all)
    kit248ExpectSkipRecorded(fixture, "E12")
})

test_that("KIT-248 E11: a second Skip commits nothing and leaves no provenance pending to be chained (fails before at its first assertion)", {
    noSkipElsewhere <- function(manager, info) {
        entries <- manager$.context$snapshot$step_provenance
        for (state in setdiff(names(entries), "correlation_filtered")) {
            expect_false(identical(entries[[state]]$operation, "pearson-skip"), info = paste(info, state))
            expect_false("pearson-skip" %in% vapply(entries[[state]]$preceded_by %||% list(), function(p) p$operation, ""),
                info = paste(info, state))
        }
        expect_length(entries$correlation_filtered$preceded_by, 0L)
    }
    # Label-free: Skip, Skip, then a normalisation re-run.
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Skip(fixture, norm)
    expect_identical(kit248Entry(fixture$manager)$operation, "pearson-skip", info = "LFQ")
    history <- diaStateHistory(fixture$manager)
    kit248Skip(fixture, norm)
    expect_identical(diaStateHistory(fixture$manager), history, info = "LFQ second Skip commits nothing")
    expect_null(fixture$manager$.context$pending_provenance)
    kit248Normalise(fixture, norm)
    noSkipElsewhere(fixture$manager, "LFQ")
    # DIA with limpa: Skip, Skip, then the RUV covariates.
    kit248NeedLimpa()
    fixture <- kit248Owner("DIA_LIMPA")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Skip(fixture, norm)
    expect_identical(kit248Entry(fixture$manager)$operation, "pearson-skip", info = "DIA_LIMPA")
    kit248Skip(fixture, norm)
    expect_null(fixture$manager$.context$pending_provenance)
    kit248RuvCovariates(fixture, norm)
    expect_identical(fixture$manager$current_state, "ruv_covariates")
    noSkipElsewhere(fixture$manager, "DIA_LIMPA")
})

test_that("KIT-248 E6: DIA with limpa, Apply at t_hi, Skip, RUV covariates, DA: recorded as the skip; the replay's table, covariates and DA equal the live ones (fails before)", {
    kit248NeedLimpa()
    thi <- kit248THi("DIA_LIMPA")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("DIA_LIMPA")
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi)
    expect_identical(setdiff(thi$all, kit248Table(manager, "correlation_filtered")$samples), thi$removed)
    kit248Skip(fixture, norm)
    live <- kit248Table(manager, "correlation_filtered")
    expect_identical(live$samples, thi$all)
    kit248ExpectSkipRecorded(fixture, "E6")
    kit248RuvCovariates(fixture, norm)
    kit248Da(fixture)
    covariates_live <- kit248Table(manager, "ruv_covariates")$design
    da_live <- kit248DaTable(manager)
    replay <- kit248Replay(fixture, kit248Record(fixture))
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    skip_if(is.null(replay$manager), "replay failed")
    expect_identical(kit248Entry(replay$manager)$operation, "pearson-skip")
    replayed <- kit248Table(replay$manager, "correlation_filtered")
    expect_identical(replayed$samples, live$samples)
    expect_identical(replayed$table, live$table)
    expect_identical(kit248Table(replay$manager, "ruv_covariates")$design, covariates_live)
    expect_identical(kit248DaTable(replay$manager), da_live)
})

test_that("KIT-248 E10: DIA with limpa, press, Apply at t_hi, press again, Skip: live and replay agree on the recorded start (fails before)", {
    kit248NeedLimpa()
    thi <- kit248THi("DIA_LIMPA")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    for (path in c("in process", if (workflowOwnerCapable()) "worker")) {
        fixture <- kit248Owner("DIA_LIMPA", path)
        manager <- fixture$manager
        l7 <- getDIAStateReference(manager)
        norm <- new.env(parent = emptyenv())
        kit248Normalise(fixture, norm)
        kit248Apply(fixture, norm, thi$t_hi, path)
        kit248Normalise(fixture, norm)
        kit248Skip(fixture, norm)
        # KIT-244 (merged): the Skip commits the recorded start, L7.
        expect_identical(kit248Components(manager, getDIAStateReference(manager)), kit248Components(manager, l7), info = path)
        live <- kit248Table(manager, "correlation_filtered")
        expect_identical(live$samples, thi$all, info = path)
        record <- kit248Record(fixture)
        expect_identical(kit248Step(record)$operation, "pearson-skip", info = path)
        replay <- kit248Replay(fixture, record)
        expect_identical(kit248ErrorText(replay$error), NA_character_, info = path)
        if (is.null(replay$manager)) next
        replayed <- kit248Table(replay$manager, "correlation_filtered")
        expect_identical(replayed$samples, live$samples, info = path)
        expect_identical(replayed$table, live$table, info = path)
    }
})

test_that("KIT-248 E9 (OD-1 (a)): a record made before the fix by Apply then Skip replays as the skip, on every sample (fails before)", {
    thi <- kit248THi("LFQ")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi)
    kit248PreFixSkip(fixture, norm)
    live <- kit248Table(fixture$manager, "correlation_filtered")
    expect_identical(live$samples, thi$all)
    record <- kit248Record(fixture)
    # The record's shape before the fix (G0b): the Apply's step, its threshold, the outcome skipped.
    step <- kit248Step(record)
    expect_identical(step$operation, "pearson-selection")
    expect_true("threshold" %in% kit248ParameterNames(step))
    expect_true(isTRUE(unlist(step$outcomes$skipped)))
    replay <- kit248Replay(fixture, record)
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    skip_if(is.null(replay$manager), "replay failed")
    expect_identical(kit248Entry(replay$manager)$operation, "pearson-skip")
    expect_equal(kit248Config(replay$manager), list(min_pearson_correlation_threshold = 0, skipped = TRUE))
    replayed <- kit248Table(replay$manager, "correlation_filtered")
    expect_identical(replayed$samples, live$samples)
    expect_identical(replayed$table, live$table)
})

test_that("KIT-248 E13 (OD-1 (a')): a legacy correlation step whose outcome says skipped replays as the skip; an earlier legacy step is still refused (fails before)", {
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248PreFixSkip(fixture, norm)
    live <- kit248Table(fixture$manager, "correlation_filtered")
    record <- kit248Record(fixture)
    step <- kit248Step(record)
    expect_identical(step$operation, "legacy")
    expect_true(isTRUE(unlist(step$outcomes$skipped)))
    replay <- kit248Replay(fixture, record)
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    if (!is.null(replay$manager)) {
        expect_identical(kit248Entry(replay$manager)$operation, "pearson-skip")
        replayed <- kit248Table(replay$manager, "correlation_filtered")
        expect_identical(replayed$samples, live$samples)
        expect_identical(replayed$table, live$table)
    }
    # An earlier legacy step (the normalisation's) is refused at that step, before the skip is reached.
    earlier <- record
    position <- which(vapply(earlier$steps, function(s) identical(s$state, "normalized"), logical(1)))
    expect_length(position, 1L)
    earlier$steps[[position]]$operation <- "legacy"
    refused <- kit248Replay(fixture, earlier)
    expect_s3_class(refused$error, "workflow_replay_unsupported")
    expect_match(kit248ErrorText(refused$error), "normalized")
})

# ---- G rows ----

test_that("KIT-248 G1: Apply only is recorded as the Apply and replays at its threshold (passes before and after)", {
    thi <- kit248THi("LFQ")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("LFQ")
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi)
    live <- kit248Table(fixture$manager, "correlation_filtered")
    expect_identical(setdiff(thi$all, live$samples), thi$removed)
    record <- kit248Record(fixture)
    step <- kit248Step(record)
    expect_identical(step$operation, "pearson-selection")
    threshold <- Filter(function(p) identical(p$name, "threshold"), step$parameters)
    expect_length(threshold, 1L)
    expect_equal(decodeWorkflowRecordValue(threshold[[1L]]$value), thi$t_hi, tolerance = 1e-12)
    expect_identical(vapply(step$preceded_by, function(p) p$operation, ""), "pearson-pairs")
    replay <- kit248Replay(fixture, record)
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    skip_if(is.null(replay$manager), "replay failed")
    expect_identical(kit248Entry(replay$manager)$operation, "pearson-selection")
    replayed <- kit248Table(replay$manager, "correlation_filtered")
    expect_identical(replayed$samples, live$samples)
    expect_identical(replayed$table, live$table)
})

test_that("KIT-248 E10b (G row; AUDIT-KIT-244-impl I05): DIA with limpa, press, Apply at t_hi, press again, Apply at 0.5: replay equals live (passes before and after)", {
    kit248NeedLimpa()
    thi <- kit248THi("DIA_LIMPA")
    skip_if(is.null(thi), "no t_hi on this fixture (stop rule 2)")
    fixture <- kit248Owner("DIA_LIMPA")
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, thi$t_hi)
    kit248Normalise(fixture, norm)
    kit248Apply(fixture, norm, 0.5)
    live <- kit248Table(manager, "correlation_filtered")
    expect_identical(live$samples, thi$all)
    record <- kit248Record(fixture)
    step <- kit248Step(record)
    expect_identical(step$operation, "pearson-selection")
    threshold <- Filter(function(p) identical(p$name, "threshold"), step$parameters)
    expect_equal(decodeWorkflowRecordValue(threshold[[1L]]$value), 0.5)
    replay <- kit248Replay(fixture, record)
    expect_identical(kit248ErrorText(replay$error), NA_character_)
    skip_if(is.null(replay$manager), "replay failed")
    replayed <- kit248Table(replay$manager, "correlation_filtered")
    expect_identical(replayed$samples, live$samples)
    expect_identical(replayed$table, live$table)
})
