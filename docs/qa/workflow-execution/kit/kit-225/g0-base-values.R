# KIT-225 gate G0 (r1.2, section 10): the base values of the proteomics "Reset to Pre-Normalization" for the cases
# E1 to E7 and the guards G1, G2, recorded with no assertion. For each case: the state Reset selects, its entity (as
# the label of the commit that made it, `NN-state`), the active history after it, the notice; for E1 and E2 the next
# run's source; for E3 (DIA with limpa, peptide-class states, A1) the outcome of the next pass-through, then of the
# correlation Apply, then in a fresh copy of the correlation Skip followed by the RUV covariates. The package is the
# base: this commit changes nothing under R/. The builders are those of the new test file
# tests/testthat/test-kit225-prot-norm-reset-source.R (copied; the test file is written after this gate).
#
# Usage (from the package root, on a runner): Rscript docs/qa/workflow-execution/kit/kit-225/g0-base-values.R
#   writes docs/qa/workflow-execution/kit/kit-225/g0-out/g0-base-values.tsv (case, step, value).

kit225Out <- "docs/qa/workflow-execution/kit/kit-225/g0-out"
dir.create(kit225Out, recursive = TRUE, showWarnings = FALSE)
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)

kit225Rows <- new.env(parent = emptyenv())
kit225Rows$rows <- list()
kit225Record <- function(case, step, value) {
    value <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(value)), collapse = ","))
    if (!nzchar(value)) value <- "<none>"
    kit225Rows$rows[[length(kit225Rows$rows) + 1L]] <- data.frame(case = case, step = step, value = value)
    cat(sprintf("KIT-225 G0 %s %s %s\n", case, step, value))
}

# ---- builders (as in the test file) ----

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

kit225Label <- function(fixture, entity_id) {
    if (is.null(entity_id)) return("<null>")
    labels <- ls(fixture$ids)
    hit <- labels[vapply(labels, function(label) identical(get(label, envir = fixture$ids), entity_id), NA)]
    if (length(hit)) paste(hit, collapse = "|") else paste0("other:", substr(entity_id, 1L, 8L))
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

kit225RuvCovariates <- function(fixture, norm, mode = "covariates") kit225Quiet(runDIAReferenceRuvCovariates(
    list(ruv_mode = mode, color_variable = "group", shape_variable = NULL), fixture$workflow, norm, fixture$paths,
    "proteomics", "core", function(reference) invisible(NULL), function() "group"))

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

# ---- G0 recording ----

kit225After <- function(case, fixture, reset) {
    force(reset) # the Reset runs before the state is read (G0 attempt 1 read it first: VAL-54 instrument defect)
    kit225Record(case, "current_state", fixture$manager$current_state)
    kit225Record(case, "head", kit225Label(fixture, kit225Head(fixture)$entity_id))
    kit225Record(case, "history", kit225History(fixture))
    kit225Record(case, "previous_state", reset$result$previousState %||% "<null>")
    kit225Record(case, "notice", reset$notes)
}

kit225Outcome <- function(case, step, code) {
    outcome <- tryCatch({
        force(code)
        "completed"
    }, error = function(e) paste0("error [", paste(class(e), collapse = "/"), "] ", conditionMessage(e)))
    kit225Record(case, step, outcome)
    outcome
}

kit225Case <- function(case, code) {
    tryCatch(force(code), error = function(e) kit225Record(case, "INSTRUMENT_ERROR",
        paste0("[", paste(class(e), collapse = "/"), "] ", conditionMessage(e))))
}

kit225RunCases <- function() {
    # E1, E2: label-free and TMT, normalised without the Protein Replicate Filter, Reset, run again.
    for (type in c("LFQ", "TMT")) {
        case <- if (identical(type, "LFQ")) "E1" else "E2"
        kit225Case(case, {
            fixture <- kit225Owner(type)
            kit225LabelFreeChain(fixture)
            norm <- new.env(parent = emptyenv())
            kit225Click(fixture, norm)
            kit225Record(case, "first_run_source", kit225Label(fixture, norm$normalization_input_reference$entity_id))
            kit225Record(case, "history_before_reset", kit225History(fixture))
            reset <- kit225Reset(fixture, norm)
            kit225After(case, fixture, reset)
            kit225Click(fixture, norm)
            kit225Record(case, "second_run_source", kit225Label(fixture, norm$normalization_input_reference$entity_id))
        })
    }
    # E3: DIA with limpa at its normal ending, pass-through, correlation Skip, Reset; then A1's outcomes.
    kit225Case("E3", {
        fixture <- kit225LimpaOwner()
        norm <- new.env(parent = emptyenv())
        kit225PassThrough(fixture, norm)
        kit225Record("E3", "designated", kit225Label(fixture, norm$normalization_input_reference$entity_id))
        kit225Correlation(fixture, norm, skipped = TRUE)
        kit225Record("E3", "history_before_reset", kit225History(fixture))
        reset <- kit225Reset(fixture, norm)
        kit225After("E3", fixture, reset)
        kit225Outcome("E3", "A1_pass_through", kit225PassThrough(fixture, norm))
        kit225Record("E3", "A1_designated", kit225Label(fixture, norm$normalization_input_reference$entity_id))
        kit225Outcome("E3", "A1_correlation_apply", kit225Correlation(fixture, norm, 0.5))
        kit225Record("E3", "A1_after_apply_state", fixture$manager$current_state)
    })
    kit225Case("E3-fresh", {
        fixture <- kit225LimpaOwner()
        norm <- new.env(parent = emptyenv())
        kit225PassThrough(fixture, norm)
        kit225Correlation(fixture, norm, skipped = TRUE)
        reset <- kit225Reset(fixture, norm)
        kit225Record("E3-fresh", "current_state_after_reset", fixture$manager$current_state)
        kit225Outcome("E3-fresh", "A1_pass_through", kit225PassThrough(fixture, norm))
        kit225Outcome("E3-fresh", "A1_correlation_skip", kit225Correlation(fixture, norm, skipped = TRUE))
        kit225Record("E3-fresh", "A1_after_skip_state", fixture$manager$current_state)
        kit225Outcome("E3-fresh", "A1_ruv_covariates", kit225RuvCovariates(fixture, norm, "covariates"))
        kit225Record("E3-fresh", "A1_after_covariates_state", fixture$manager$current_state)
    })
    # E4: DIA (iq), peptide states through imputed, protein states through the intensity filter (D3), run, Reset.
    kit225Case("E4", {
        fixture <- kit225Owner("DIA")
        precursors <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        for (state in c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
            "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed")) kit225Commit(fixture, state, precursors)
        for (state in c("protein_s4_created", "protein_accession_cleaned", "protein_intensity_filtered"))
            kit225Commit(fixture, state)
        norm <- new.env(parent = emptyenv())
        kit225Click(fixture, norm)
        kit225Record("E4", "first_run_source", kit225Label(fixture, norm$normalization_input_reference$entity_id))
        kit225After("E4", fixture, kit225Reset(fixture, norm))
    })
    # E5: label-free, design saved again after a run (KIT-182 S9).
    kit225Case("E5", {
        fixture <- kit225Owner("LFQ")
        kit225LabelFreeChain(fixture)
        kit225Commit(fixture, "protein_replicate_filtered")
        norm <- new.env(parent = emptyenv())
        kit225Click(fixture, norm)
        kit225LabelFreeChain(fixture)
        kit225Record("E5", "history_after_second_chain", kit225History(fixture))
        kit225Click(fixture, norm)
        kit225Record("E5", "second_run_source", kit225Label(fixture, norm$normalization_input_reference$entity_id))
        kit225After("E5", fixture, kit225Reset(fixture, norm))
    })
    # E6: label-free at duplicates_removed, no run.
    kit225Case("E6", {
        fixture <- kit225Owner("LFQ")
        kit225LabelFreeChain(fixture)
        kit225After("E6", fixture, kit225Reset(fixture, new.env(parent = emptyenv())))
    })
    # E6b: as E6, the pre-normalisation figures' source set to the intensity filter's table.
    kit225Case("E6b", {
        fixture <- kit225Owner("LFQ")
        kit225LabelFreeChain(fixture)
        kit225Scope(fixture, function() {
            snapshot <- fixture$manager$.context$snapshot
            snapshot$qc_sources$pre_normalization <- getDIAStateReference(fixture$manager, "protein_intensity_filtered")
            publishDIAStateHead(fixture$manager$.context, snapshot)
        })
        kit225After("E6b", fixture, kit225Reset(fixture, new.env(parent = emptyenv())))
    })
    # E6c: as E6, the recorded starting table a reference no active state holds (a Revert past the replicate filter).
    kit225Case("E6c", {
        fixture <- kit225Owner("LFQ")
        kit225LabelFreeChain(fixture)
        prf <- kit225Commit(fixture, "protein_replicate_filtered")
        kit225Scope(fixture, function() selectDIAStateReference(fixture$manager, "duplicates_removed"))
        norm <- new.env(parent = emptyenv())
        norm$normalization_input_reference <- prf
        kit225After("E6c", fixture, kit225Reset(fixture, norm))
    })
    # E7: recorded starting table of no state, head at normalized.
    kit225Case("E7", {
        fixture <- kit225Owner("DIA")
        prf1 <- kit225Commit(fixture, "protein_replicate_filtered")
        norm <- new.env(parent = emptyenv())
        kit225Click(fixture, norm)
        kit225Record("E7", "history_before_reset", kit225History(fixture))
        unknown <- prf1
        unknown$entity_id <- strrep("0", 64L)
        norm$normalization_input_reference <- unknown
        outcome <- kit225Outcome("E7", "reset", reset <- kit225Reset(fixture, norm))
        kit225Record("E7", "current_state", fixture$manager$current_state)
        kit225Record("E7", "history", kit225History(fixture))
    })
    # G1, G2: QC complete, run, correlation Apply, Reset; the whole state against the state before the run.
    for (type in c("DIA", "LFQ")) {
        case <- if (identical(type, "DIA")) "G1" else "G2"
        kit225Case(case, {
            fixture <- kit225Owner(type)
            if (identical(type, "LFQ")) kit225LabelFreeChain(fixture)
            prf1 <- kit225Commit(fixture, "protein_replicate_filtered")
            history <- kit225History(fixture)
            provenance <- names(fixture$manager$.context$snapshot$step_provenance)
            norm <- new.env(parent = emptyenv())
            kit225Click(fixture, norm)
            kit225Correlation(fixture, norm, 0.5)
            kit225Record(case, "history_before_reset", kit225History(fixture))
            reset <- kit225Reset(fixture, norm)
            kit225After(case, fixture, reset)
            kit225Record(case, "history_equal", identical(kit225History(fixture), history))
            kit225Record(case, "head_is_prf1", identical(kit225Head(fixture)$entity_id, prf1$entity_id))
            kit225Record(case, "provenance_equal",
                identical(names(fixture$manager$.context$snapshot$step_provenance), provenance))
            kit225Record(case, "input_reference_null", is.null(norm$normalization_input_reference))
            kit225Record(case, "tab_status", unlist(fixture$workflow$tab_status[c("normalization", "differential_expression")]))
        })
    }
}

kit225RunCases()
rows <- do.call(rbind, kit225Rows$rows)
utils::write.table(rows, file.path(kit225Out, "g0-base-values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("KIT-225 G0: %d rows; instrument errors: %d\n", nrow(rows), sum(rows$step == "INSTRUMENT_ERROR")))
