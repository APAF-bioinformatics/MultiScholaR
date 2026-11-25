# KIT-225 (KIT-182 F-b; triage H-5): the proteomics "Reset to Pre-Normalization" returns to the table the normalisation
# started from (the recorded source, matched by its saved identity in the active history), changes no state when the
# active history holds no normalisation-tab state, and refuses (nothing changed) when the recorded source is in no state
# of the active history. On native owners as protNormCoreOwner makes them (helper-prot-norm-core.R), with each case's
# states committed by `manager$saveState` (the fixture's own form) and the workflow type of the case; E3 on the DIA with
# limpa recipe of protNormW1LimpaOwner (peptide-class peptide states). The builders are the G0 script's
# (docs/qa/workflow-execution/kit/kit-225/g0-base-values.R). The texts are drafts (VAL-41; kit-225/wording.md).

kit225Quiet <- function(code) suppressWarnings(suppressMessages(code))
kit225Call <- function(fn, args) do.call(fn, args[intersect(names(args), names(formals(fn)))])
kit225Scope <- function(fixture, code) withDIAStateScope(fixture$manager, code)
kit225Head <- function(fixture) getDIAStateReference(fixture$manager)
kit225History <- function(fixture) diaStateHistory(fixture$manager)
kit225Settings <- function(mode = "skip") list(norm_method = "quantile", ruv_mode = mode, auto_percentage_min = NULL,
    auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)

#' A native proteomics owner as protNormCoreOwner makes it (helper-prot-norm-core.R:10-31), with no state committed after
#' the import and the workflow type of the case: each case commits its own states (`kit225Commit`).
kit225Owner <- function(type = "LFQ", envir = parent.frame()) {
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
        ids = new.env(parent = emptyenv()))
}

#' The head recorded under the label `NN-<state>` (NN: the commit's order in the case).
kit225Remember <- function(fixture) {
    reference <- kit225Head(fixture)
    label <- sprintf("%02d-%s", length(ls(fixture$ids)) + 1L, fixture$manager$current_state)
    assign(label, reference$entity_id, envir = fixture$ids)
    reference
}

#' A state committed as the fixture commits it (`manager$saveState`), recorded.
kit225Commit <- function(fixture, state, object = fixture$original) {
    kit225Quiet(fixture$manager$saveState(state, object, list(), "kit225"))
    kit225Remember(fixture)
}

#' The tab's in-process normalisation click (KIT-182's k182Click).
kit225Click <- function(fixture, norm, mode = "skip") {
    kit225Quiet(kit225Call(runDIAReferenceNormalization, list(settings = kit225Settings(mode),
        workflow = fixture$workflow, norm_data = norm, paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", post_qc = function(reference) invisible(NULL),
        ruv_qc = function(reference) invisible(NULL), grouping = function() "group")))
}

kit225Correlation <- function(fixture, norm, threshold = 0.5, skipped = FALSE) {
    kit225Quiet(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped, "proteomics",
        "core", if (skipped) function() NULL else function() "group", function() list(color_var = "group", shape_var = NULL)))
}

kit225PassThrough <- function(fixture, norm) kit225Quiet(runDIAReferenceLimpaNormalisationPassThrough(
    kit225Settings("covariates"), fixture$workflow, norm, function(reference) invisible(NULL)))

#' The tab's Reset (`runProtNormResetWorkflow`, as KIT-182's T5 presses it), its notices recorded.
kit225Reset <- function(fixture, norm) {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    note <- function(ui, type = "default", duration = 5, ...) {
        seen$notes <- c(seen$notes, sprintf("[%s %s] %s", type, format(duration), paste(as.character(ui), collapse = " ")))
        invisible(NULL)
    }
    seen$result <- kit225Quiet(kit225Scope(fixture, function() runProtNormResetWorkflow(fixture$workflow, norm,
        output = new.env(), resetOutputsFn = function(...) NULL, showNotificationFn = note, messageFn = function(...) NULL)))
    as.list(seen)
}

#' The label-free QC chain without the Protein Replicate Filter: P0 (design), accession, intensity, P3 (duplicates).
kit225LabelFreeChain <- function(fixture) {
    for (state in c("protein_s4_initial", "protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed"))
        kit225Commit(fixture, state)
    kit225Head(fixture)
}

#' The DIA with limpa owner: protNormW1LimpaOwner's recipe, objects and kernels (helper-prot-norm-core.R:63-85), with the
#' peptide states committed first (the helper itself commits protein_dpc_quant first, so it is not called).
kit225LimpaOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "kit225-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- rollUpPrecursorToPeptide(precursors, core_utilisation = NA)
    quantified <- kit225Quiet(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    fixture <- list(root = root, paths = paths, workflow = workflow, manager = manager, ids = new.env(parent = emptyenv()))
    kit225Commit(fixture, "raw_data_s4", precursors)
    kit225Commit(fixture, "qvalue_filtered", precursors)
    for (state in c("precursor_rollup", "log2_matrix", "peptide_accession_cleaned")) kit225Commit(fixture, state, rolled)
    setDIAStepProvenance(manager, "design_builder")
    kit225Commit(fixture, "protein_dpc_quant", quantified)
    kit225Quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        kit225Remember(fixture)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
        kit225Remember(fixture)
    }))
    fixture
}

#' The drafts of the ticket's section 7 (VAL-41).
kit225NoticeText <- function(title) paste0("[warning 5] Normalisation reset: the current table is again the one the ",
    "normalisation started from (", title, ").")
kit225NothingText <- paste("[warning 5] Normalisation reset: no normalisation result was current, so the Quality",
    "Control steps were left as they are.")

kit225Unchanged <- function(fixture, history, head) {
    expect_identical(kit225History(fixture), history)
    expect_identical(kit225Head(fixture)$entity_id, head$entity_id)
}

test_that("KIT-225 E1, E2, E8: label-free and TMT normalised without the replicate filter, Reset returns its source (fails before)", {
    for (type in c("LFQ", "TMT")) {
        fixture <- kit225Owner(type)
        p3 <- kit225LabelFreeChain(fixture)
        norm <- new.env(parent = emptyenv())
        kit225Click(fixture, norm)
        expect_identical(norm$normalization_input_reference$entity_id, p3$entity_id, info = type)
        reset <- kit225Reset(fixture, norm)
        expect_identical(fixture$manager$current_state, "duplicates_removed", info = type)
        expect_identical(kit225Head(fixture)$entity_id, p3$entity_id, info = type)
        expect_identical(reset$result$previousState, "duplicates_removed", info = type)
        # E8: the notice names the step as the Quality Control tab does.
        expect_identical(reset$notes, kit225NoticeText("Duplicate Removal"), info = type)
        expect_no_match(reset$notes, "duplicates_removed", fixed = TRUE)
        kit225Click(fixture, norm)
        expect_identical(norm$normalization_input_reference$entity_id, p3$entity_id, info = type)
    }
})

test_that("KIT-225 E3: DIA with limpa at its normal ending, Reset returns the Peptide Count Filter's table (fails before)", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- kit225LimpaOwner()
    l7 <- kit225Head(fixture)
    expect_identical(fixture$manager$current_state, "peptide_count_filtered")
    norm <- new.env(parent = emptyenv())
    kit225PassThrough(fixture, norm)
    kit225Correlation(fixture, norm, skipped = TRUE)
    expect_identical(fixture$manager$current_state, "correlation_filtered")
    reset <- kit225Reset(fixture, norm)
    expect_identical(fixture$manager$current_state, "peptide_count_filtered")
    expect_identical(kit225Head(fixture)$entity_id, l7$entity_id)
    expect_false("correlation_filtered" %in% kit225History(fixture))
    expect_identical(reset$notes, kit225NoticeText("Peptide Count Filter"))
})

test_that("KIT-225 E4: DIA (iq) normalised before the replicate filter, Reset returns the protein table (fails before)", {
    fixture <- kit225Owner("DIA")
    precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    for (state in c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
        "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed")) kit225Commit(fixture, state, precursors)
    for (state in c("protein_s4_created", "protein_accession_cleaned", "protein_intensity_filtered"))
        d3 <- kit225Commit(fixture, state)
    norm <- new.env(parent = emptyenv())
    kit225Click(fixture, norm)
    kit225Reset(fixture, norm)
    expect_identical(fixture$manager$current_state, "protein_intensity_filtered")
    expect_identical(kit225Head(fixture)$entity_id, d3$entity_id)
})

test_that("KIT-225 E5: after a design saved again, Reset returns the second run's source by identity, not by name (fails before)", {
    fixture <- kit225Owner("LFQ")
    kit225LabelFreeChain(fixture)
    kit225Commit(fixture, "protein_replicate_filtered")
    norm <- new.env(parent = emptyenv())
    kit225Click(fixture, norm)
    p3b <- kit225LabelFreeChain(fixture)
    kit225Click(fixture, norm)
    expect_identical(norm$normalization_input_reference$entity_id, p3b$entity_id)
    kit225Reset(fixture, norm)
    expect_identical(kit225Head(fixture)$entity_id, p3b$entity_id)
    expect_identical(fixture$manager$current_state, "duplicates_removed")
})

test_that("KIT-225 E6, E6b, E6c: with no normalisation-tab state in the history Reset leaves Quality Control as it is (fails before)", {
    # E6: no run.
    fixture <- kit225Owner("LFQ")
    p3 <- kit225LabelFreeChain(fixture)
    history <- kit225History(fixture)
    reset <- kit225Reset(fixture, new.env(parent = emptyenv()))
    kit225Unchanged(fixture, history, p3)
    expect_null(reset$result$previousState)
    expect_identical(reset$notes, kit225NothingText)
    # E6b: the pre-normalisation figures' source is an earlier QC table (as entering the tab earlier leaves it).
    fixture <- kit225Owner("LFQ")
    p3 <- kit225LabelFreeChain(fixture)
    kit225Scope(fixture, function() {
        snapshot <- fixture$manager$.context$snapshot
        snapshot$qc_sources$pre_normalization <- getDIAStateReference(fixture$manager, "protein_intensity_filtered")
        publishDIAStateHead(fixture$manager$.context, snapshot)
    })
    history <- kit225History(fixture)
    reset <- kit225Reset(fixture, new.env(parent = emptyenv()))
    kit225Unchanged(fixture, history, p3)
    expect_identical(reset$notes, kit225NothingText)
    # E6c: the recorded starting table is a reference no active state holds (the restore's fill, then a Revert).
    fixture <- kit225Owner("LFQ")
    kit225LabelFreeChain(fixture)
    prf <- kit225Commit(fixture, "protein_replicate_filtered")
    p3 <- kit225Scope(fixture, function() selectDIAStateReference(fixture$manager, "duplicates_removed"))
    history <- kit225History(fixture)
    norm <- new.env(parent = emptyenv())
    norm$normalization_input_reference <- prf
    reset <- kit225Reset(fixture, norm)
    kit225Unchanged(fixture, history, p3)
    expect_identical(reset$notes, kit225NothingText)
})

test_that("KIT-225 E7: a recorded source in no state of the active history is refused and nothing changes (fails before)", {
    fixture <- kit225Owner("DIA")
    prf1 <- kit225Commit(fixture, "protein_replicate_filtered")
    norm <- new.env(parent = emptyenv())
    kit225Click(fixture, norm)
    history <- kit225History(fixture)
    head <- kit225Head(fixture)
    unknown <- prf1
    unknown$entity_id <- strrep("0", 64L)
    norm$normalization_input_reference <- unknown
    expect_error(kit225Reset(fixture, norm), class = "prot_norm_reset_source_unknown")
    kit225Unchanged(fixture, history, head)
    expect_identical(fixture$manager$current_state, "normalized")
    expect_identical(norm$normalization_input_reference$entity_id, unknown$entity_id)
})

test_that("KIT-225 G1, G2: QC complete, run, correlation Apply, Reset gives back the state before the run (as before)", {
    for (type in c("DIA", "LFQ")) {
        fixture <- kit225Owner(type)
        if (identical(type, "LFQ")) kit225LabelFreeChain(fixture)
        prf1 <- kit225Commit(fixture, "protein_replicate_filtered")
        history <- kit225History(fixture)
        provenance <- names(fixture$manager$.context$snapshot$step_provenance)
        norm <- new.env(parent = emptyenv())
        kit225Click(fixture, norm)
        kit225Correlation(fixture, norm, 0.5)
        expect_identical(fixture$manager$current_state, "correlation_filtered", info = type)
        kit225Reset(fixture, norm)
        expect_identical(kit225History(fixture), history, info = type)
        expect_identical(fixture$manager$current_state, "protein_replicate_filtered", info = type)
        expect_identical(kit225Head(fixture)$entity_id, prf1$entity_id, info = type)
        expect_identical(names(fixture$manager$.context$snapshot$step_provenance), provenance, info = type)
        expect_null(norm$normalization_input_reference)
        expect_identical(fixture$workflow$tab_status$normalization, "pending", info = type)
        expect_identical(fixture$workflow$tab_status$differential_expression, "disabled", info = type)
    }
})
