# KIT-174 part (a): the proteomics RUV grouping help says what the column does (T1); a
# demonstration of both halves of what it says (T2); the route it points to for a known batch
# runs (T9); and the run-time warning when the chosen column's values span more than one group
# (T10; the owner's decision VAL-29, at the three run actions).
library(testthat)

# The approved help (ticket section 5.1; implemented as proposed, the owner edits before
# release, VAL-41).
kit174HelpText <- paste(
    "Choose the column that holds your experimental groups (normally 'group'). Samples that",
    "share a value in this column are treated as replicates of each other: RUV-III takes what",
    "differs between them to be unwanted variation and builds its correction from it. In the",
    "Automatic and Manual modes, and in both RUV modes of DIA with limpa, the same column is also",
    "used to pick the control proteins (those that differ least between its values), and, except",
    "in Manual mode, the number of factors; the sample-correlation filter on this tab compares",
    "samples within it. Do not choose 'batch', or any column in which a value contains samples of",
    "more than one group: the correction would then be built from the differences between your",
    "groups, and the differences between batches would stay in the data. If you know each",
    "sample's batch, add it to the Model Formula on the Differential Expression tab instead, for",
    "example ~ 0 + group + batch. If the design declares technical replicates, those are the",
    "replicates RUV-III uses; this column still picks the control proteins."
)

# The approved run-time warning (5.1): the column, the values that span groups, all values.
kit174WarningText <- function(column, n, m) {
    sprintf(
        paste(
            "The RUV grouping column '%s' puts samples of different groups together (%d of its %d",
            "values contain samples of more than one group). RUV-III treats samples with the same",
            "value as replicates, so the correction would be built from the differences between",
            "your groups. Choose 'group' unless each value of the column lies inside one group."
        ),
        column, n, m
    )
}

kit174Quiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressMessages(suppressWarnings(expr)))
    value
}

kit174HasCheck <- function() {
    exists("omicNormRuvGroupingSpansGroups", envir = asNamespace("MultiScholaR"), inherits = FALSE)
}

#' 2 groups x 2 batches x 3 runs. `declared`: each group's six runs are three units measured
#' twice (technical replicates declared through `tech_rep_group`, as the design step writes
#' it). Other columns: `nested` (each value inside one group), `partly` (value `x` inside
#' group A, value `y` across both groups), `empty` (all NA).
kit174Design <- function(declared = FALSE) {
    group <- rep(c("A", "B"), each = 6L)
    batch <- rep(rep(c("b1", "b2"), each = 3L), 2L)
    replicates <- if (declared) rep(c(1L, 1L, 2L, 2L, 3L, 3L), 2L) else rep(1:6, 2L)
    design <- data.frame(
        Run = sprintf("S%02d", 1:12), group = group, batch = batch, replicates = replicates,
        tech_rep_group = paste(group, replicates, sep = "_"),
        nested = paste(group, batch, sep = "_"),
        partly = c(rep("x", 3L), rep("y", 9L)), empty = NA_character_, stringsAsFactors = FALSE
    )
    rownames(design) <- design$Run
    design
}

#' What a run action shows and logs: notifications "[type] text", logger warnings, the tab's
#' log, the error handler's messages, whether the run ran.
kit174Seen <- function() {
    seen <- new.env(parent = emptyenv())
    seen$notes <- character()
    seen$warn <- character()
    seen$log <- character()
    seen$errors <- character()
    seen$ran <- FALSE
    seen
}

kit174Count <- function(seen, text) {
    c(
        notifications = sum(seen$notes == paste0("[warning] ", text)),
        log_warn = sum(seen$warn == text),
        warnings = sum(startsWith(seen$notes, "[warning]"))
    )
}

kit174Once <- c(notifications = 1L, log_warn = 1L, warnings = 1L)

kit174Note <- function(seen) {
    function(ui, ..., type = "default") {
        text <- paste(as.character(ui), collapse = " ")
        seen$notes <- c(seen$notes, paste0("[", type, "] ", text))
        invisible(NULL)
    }
}

kit174MockLogger <- function(seen, env = parent.frame()) {
    testthat::local_mocked_bindings(
        log_warn = function(...) {
            seen$warn <- c(seen$warn, paste0(...))
            invisible(NULL)
        },
        .package = "logger", .env = env
    )
}

#' Classic proteomics: the normalisation action in process, the runner stubbed.
kit174ClassicRun <- function(column, design, declared = FALSE, input = list()) {
    seen <- kit174Seen()
    testthat::local_mocked_bindings(
        runDIAReferenceNormalization = function(...) {
            seen$ran <- TRUE
            TRUE
        },
        recordProtNormMethod = \(...) invisible(NULL),
        updateProtNormRuvAuditTrail = \(...) invisible(NULL),
        workflowOmicsManagedActions = \(...) NULL,
        omicNormRuvHeadTechnicalReplicates = \(manager) declared
    )
    kit174MockLogger(seen)
    input <- utils::modifyList(list(
        norm_method = "cyclicloess", ruv_mode = "automatic", auto_percentage_min = 1,
        auto_percentage_max = 20, separation_metric = "max_difference", k_penalty_weight = 0.5,
        adaptive_k_penalty = TRUE, ruv_percentage = 5, ruv_k = 1, color_variable = "group",
        shape_variable = "group", ruv_grouping_variable = column
    ), input)
    workflow <- list2env(list(state_manager = list(workflow_type = "LFQ"), design_matrix = design))
    kit174Quiet(runProtNormNormalizationWorkflow(
        input, workflow, list2env(list()), list(), "proteomics", "k174", \(...) NULL, \(...) NULL,
        \() column,
        showNotificationFn = kit174Note(seen),
        handleErrorFn = function(error) seen$errors <- c(seen$errors, conditionMessage(error)),
        messageFn = \(...) NULL
    ))
    seen
}

#' DIA with limpa: the RUV-covariates action, the runner stubbed.
kit174DiaRun <- function(column, design, mode = "covariates", declared = FALSE) {
    seen <- kit174Seen()
    testthat::local_mocked_bindings(
        runDIAReferenceRuvCovariates = function(settings, workflowData, normData, ...) {
            seen$ran <- TRUE
            normData$best_k <- 1L
            normData$control_genes_index <- c(TRUE, FALSE, TRUE)
            normData$ruv_optimization_result <- list(
                best_percentage = 5, replicate_structure = "group"
            )
            "accepted"
        },
        updateProtNormRuvAuditTrail = \(...) invisible(NULL),
        workflowOmicsManagedActions = \(...) NULL,
        omicNormRuvHeadTechnicalReplicates = \(manager) declared
    )
    kit174MockLogger(seen)
    workflow <- list2env(list(
        state_manager = list(workflow_type = "DIA_LIMPA"), design_matrix = design
    ))
    norm_data <- list2env(list(correlation_filtering_complete = TRUE))
    input <- list(
        ruv_mode = mode, color_variable = NULL, shape_variable = NULL,
        ruv_grouping_variable = column
    )
    kit174Quiet(runProtNormRuvCovariatesObserver(
        input, list2env(list()), workflow, norm_data, list(), "proteomics", "k174",
        \(reference) NULL, \() column,
        renderTextFn = \(text) text, showNotificationFn = kit174Note(seen), messageFn = \(...) NULL
    ))
    seen
}

#' The shared frame (lipidomics): "Run Normalization Pipeline" in process, the omic's runner
#' stubbed.
kit174FrameRun <- function(column, design, input = list()) {
    seen <- kit174Seen()
    testthat::local_mocked_bindings(
        showNotification = kit174Note(seen), withProgress = \(expr, ...) expr,
        incProgress = \(...) invisible(NULL), .package = "shiny"
    )
    testthat::local_mocked_bindings(
        omicNormApplyRunStart = \(...) NULL, workflowOmicsManagedActions = \(...) NULL,
        omicNormRuvHeadTechnicalReplicates = \(manager) FALSE
    )
    kit174MockLogger(seen)
    in_process <- function(...) {
        seen$ran <- TRUE
        invisible(TRUE)
    }
    config <- list(
        binding = list(binding = list(
            omic_type = "lipidomics", prefix = "lipid", qc_dir_role = "qc"
        )),
        norm = list(run = list(in_process = in_process), panels = list(log2 = list(offset = 0)))
    )
    input <- utils::modifyList(list(
        apply_itsd = FALSE, itsd_aggregation = "median", log_offset = 0, norm_method = "quantile",
        ruv_mode = "automatic", auto_percentage_min = 1, auto_percentage_max = 20,
        ruv_grouping_variable = column, separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 1, ruv_percentage = 5
    ), input)
    workflow <- list2env(list(
        state_manager = new.env(), design_matrix = design, config_list = list()
    ))
    norm_data <- list2env(list(assay_names = "assay1", plot_refresh_trigger = 0))
    kit174Quiet(runOmicNormPipeline(
        input, workflow, list(qc = tempdir()), norm_data,
        \(text) seen$log <- c(seen$log, text),
        \() list(color_var = "group", shape_var = NULL), config
    ))
    seen
}

# Gate-on inputs of the RUV design section (KIT-129): the design and its k rule.
kit174GateInput <- function(design) {
    list(
        ruv_design = design, ruv_k_rule = "num_sv_be", ruv_k_user = NA, ruv_k_seed = 7,
        ruv_control_list = ""
    )
}

#' A ProteinQuantitativeData object of 2 groups x 2 batches x 3 runs from a protein x sample
#' matrix; no technical replicates declared.
kit174ProteinObject <- function(values) {
    samples <- colnames(values)
    group <- rep(c("A", "B"), each = 6L)
    ids <- rownames(values)
    table <- data.frame(Protein.Ids = ids, values, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(table) <- NULL # the kernels move the id column to the row names
    design <- data.frame(
        Run = samples, group = group, batch = rep(rep(c("b1", "b2"), each = 3L), 2L),
        replicates = rep(1:6, 2L), tech_rep_group = paste(group, rep(1:6, 2L), sep = "_"),
        stringsAsFactors = FALSE
    )
    ProteinQuantitativeData(
        protein_quant_table = table, protein_id_column = "Protein.Ids",
        protein_id_table = data.frame(Protein.Ids = ids, stringsAsFactors = FALSE),
        design_matrix = design, sample_id = "Run", group_id = "group",
        technical_replicate_id = "tech_rep_group", args = list()
    )
}

test_that("T1 (KIT-174 (a)): the proteomics RUV grouping help is the approved text", {
    help <- omicGuiConfig("proteomics")$norm$panels$ruv_iii$help[["grouping"]]
    expect_identical(help, kit174HelpText)
    expect_false(grepl("select 'batch'", help, fixed = TRUE))
    expect_true(grepl("Model Formula on the Differential Expression tab", help, fixed = TRUE))
})

test_that("T2 (KIT-174 (a), shown): group as replicate sets removes batch, batch keeps it", {
    set.seed(174L)
    samples <- sprintf("S%02d", 1:12)
    group <- rep(c("A", "B"), each = 6L)
    batch <- rep(rep(c("b1", "b2"), each = 3L), 2L)
    proteins <- 300L
    ids <- sprintf("P%03d", seq_len(proteins))
    loading <- stats::rnorm(proteins)
    effect <- sample(c(-1, 1), proteins, replace = TRUE)
    noise <- matrix(stats::rnorm(proteins * 12L, 0, 0.2), proteins, 12L)
    values <- 20 + outer(loading, ifelse(batch == "b2", 0.5, -0.5)) +
        outer(effect, as.numeric(group == "B")) + noise
    dimnames(values) <- list(ids, samples)
    object <- kit174ProteinObject(values)
    X <- cbind(1, as.numeric(group == "B"), as.numeric(batch == "b2"))
    coefficients <- function(m) {
        B <- solve(crossprod(X), crossprod(X, t(m[ids, samples, drop = FALSE])))
        c(group = stats::median(abs(B[2L, ])), batch = stats::median(abs(B[3L, ])))
    }
    corrected <- function(column) {
        controls <- kit174Quiet(getNegCtrlProtAnova(
            object, ruv_grouping_variable = column, percentage_as_neg_ctrl = 10,
            ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH"
        ))
        out <- kit174Quiet(ruvIII_C_Varying(
            object, ruv_grouping_variable = column, ruv_number_k = 1L, ctrl = controls
        ))
        m <- as.matrix(out@protein_quant_table[, samples])
        rownames(m) <- out@protein_quant_table$Protein.Ids
        list(matrix = m, controls = sum(controls))
    }
    before <- coefficients(values)
    by_group <- corrected("group")
    by_batch <- corrected("batch")
    expect_setequal(rownames(by_group$matrix), ids)
    expect_setequal(rownames(by_batch$matrix), ids)
    after_group <- coefficients(by_group$matrix)
    after_batch <- coefficients(by_batch$matrix)
    observed <- sprintf(
        paste(
            "KIT-174 T2 observed: before group %.4f batch %.4f; with 'group' (%d controls) group",
            "%.4f batch %.4f; with 'batch' (%d controls) group %.4f batch %.4f"
        ),
        before[["group"]], before[["batch"]], by_group$controls, after_group[["group"]],
        after_group[["batch"]], by_batch$controls, after_batch[["group"]], after_batch[["batch"]]
    )
    cat(observed, "\n")
    # With 'group': the batch component removed, the group effect kept.
    expect_lt(after_group[["batch"]], before[["batch"]] / 4, label = observed)
    expect_gt(after_group[["group"]], 0.8, label = observed)
    # With 'batch': the batch component kept, the group effect eroded.
    expect_gt(after_batch[["batch"]], 0.75 * before[["batch"]], label = observed)
    expect_lt(after_batch[["group"]], 0.5, label = observed)
})

test_that("T9 (KIT-174 (a)): the classic route's differential analysis runs with a batch term", {
    set.seed(1740L)
    samples <- sprintf("S%02d", 1:12)
    group <- rep(c("A", "B"), each = 6L)
    batch <- rep(rep(c("b1", "b2"), each = 3L), 2L)
    ids <- sprintf("P%02d", 1:40)
    values <- 20 + outer(stats::rnorm(40L), as.numeric(batch == "b2")) +
        outer(stats::rnorm(40L), as.numeric(group == "B")) +
        matrix(stats::rnorm(40L * 12L, 0, 0.3), 40L, 12L)
    dimnames(values) <- list(ids, samples)
    result <- kit174Quiet(differentialAbundanceAnalysis(
        kit174ProteinObject(values),
        contrasts_tbl = data.frame(contrasts = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
        formula_string = "~ 0 + group + batch", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = FALSE, eBayes_robust = FALSE, args_row_id = "Protein.Ids",
        qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"
    ))
    tested <- result$contrasts_results$results[[1L]]
    expect_setequal(rownames(tested), ids)
    expect_true(all(is.finite(tested$logFC)))
    expect_true(all(is.finite(tested$raw_pvalue)))
    expect_true("batchb2" %in% colnames(result$contrasts_results$fit.eb$design))
    expect_identical(result$da_model$formula_used, "~ 0 + group + batch")
})

test_that("T10 (KIT-174 (a), VAL-29): the check names the column's values that span groups", {
    expect_true(kit174HasCheck())
    if (kit174HasCheck()) {
        spans <- get("omicNormRuvGroupingSpansGroups", envir = asNamespace("MultiScholaR"))
        design <- kit174Design()
        expect_length(spans(design, "group"), 0L)
        expect_setequal(spans(design, "batch"), c("b1", "b2"))
        expect_length(spans(design, "nested"), 0L)
        expect_identical(as.character(spans(design, "empty")), NA_character_)
        expect_identical(as.character(spans(design, "partly")), "y")
        expect_length(spans(design, NULL), 0L)
    }
})

test_that("T10 (KIT-174 (a), VAL-29): classic proteomics warns once, only for replicate sets", {
    expect_true(kit174HasCheck())
    batch <- kit174WarningText("batch", 2L, 2L)
    crossed <- kit174ClassicRun("batch", kit174Design())
    expect_true(crossed$ran)
    expect_identical(crossed$errors, character())
    expect_identical(kit174Count(crossed, batch), kit174Once)
    partly <- kit174ClassicRun("partly", kit174Design())
    expect_identical(kit174Count(partly, kit174WarningText("partly", 1L, 2L)), kit174Once)
    grouped <- kit174ClassicRun("group", kit174Design())
    expect_true(grouped$ran)
    expect_identical(kit174Count(grouped, batch)[["warnings"]], 0L)
    skipped <- kit174ClassicRun("batch", kit174Design(), input = list(ruv_mode = "skip"))
    expect_true(skipped$ran)
    expect_identical(kit174Count(skipped, batch)[["warnings"]], 0L)
    declared <- kit174ClassicRun("batch", kit174Design(declared = TRUE), declared = TRUE)
    expect_true(declared$ran)
    expect_identical(kit174Count(declared, batch)[["warnings"]], 0L)
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    technical <- kit174ClassicRun("batch", kit174Design(declared = TRUE), declared = TRUE,
        input = kit174GateInput("technical_replicates"))
    expect_true(technical$ran)
    expect_identical(kit174Count(technical, batch)[["warnings"]], 0L)
    # C01: under groups_all_features the classic routes withdraw the declaration for the RUV
    # step, so the column is the replicate sets whatever is declared.
    recommended <- kit174ClassicRun("batch", kit174Design(declared = TRUE), declared = TRUE,
        input = kit174GateInput("groups_all_features"))
    expect_true(recommended$ran)
    expect_identical(kit174Count(recommended, batch), kit174Once)
})

test_that("T10 (KIT-174 (a), VAL-29): DIA with limpa warns once at the RUV-covariates action", {
    expect_true(kit174HasCheck())
    batch <- kit174WarningText("batch", 2L, 2L)
    crossed <- kit174DiaRun("batch", kit174Design())
    expect_true(crossed$ran)
    expect_identical(kit174Count(crossed, batch), kit174Once)
    grouped <- kit174DiaRun("group", kit174Design())
    expect_true(grouped$ran)
    expect_identical(kit174Count(grouped, batch)[["warnings"]], 0L)
    skipped <- kit174DiaRun("batch", kit174Design(), mode = "skip")
    expect_identical(kit174Count(skipped, batch)[["warnings"]], 0L)
    declared <- kit174DiaRun("batch", kit174Design(declared = TRUE), declared = TRUE)
    expect_true(declared$ran)
    expect_identical(kit174Count(declared, batch)[["warnings"]], 0L)
})

test_that("T10 (KIT-174 (a), VAL-29): the shared frame warns once (lipidomics); trx never", {
    expect_true(kit174HasCheck())
    batch <- kit174WarningText("batch", 2L, 2L)
    crossed <- kit174FrameRun("batch", kit174Design())
    expect_true(crossed$ran)
    expect_identical(kit174Count(crossed, batch), kit174Once)
    expect_identical(sum(grepl(batch, crossed$log, fixed = TRUE)), 1L)
    grouped <- kit174FrameRun("group", kit174Design())
    expect_true(grouped$ran)
    expect_identical(kit174Count(grouped, batch)[["warnings"]], 0L)
    skipped <- kit174FrameRun("batch", kit174Design(), input = list(ruv_mode = "skip"))
    expect_true(skipped$ran)
    expect_identical(kit174Count(skipped, batch)[["warnings"]], 0L)
    withr::with_options(list(MultiScholaR.ruv_design_panel = TRUE), {
        technical <- kit174FrameRun("batch", kit174Design(),
            input = kit174GateInput("technical_replicates"))
    })
    expect_true(technical$ran)
    expect_identical(kit174Count(technical, batch)[["warnings"]], 0L)
    # Transcriptomics has no selector: its frame's run (runOmicNormHandler) raises no warning.
    seen <- kit174Seen()
    testthat::local_mocked_bindings(showNotification = kit174Note(seen), .package = "shiny")
    kit174MockLogger(seen)
    handler <- function(...) {
        seen$ran <- TRUE
        invisible(TRUE)
    }
    trx <- list(
        binding = list(binding = list(omic_type = "transcriptomics")),
        norm = list(run = list(settings = \(input) list(), handler = handler))
    )
    workflow <- list2env(list(
        state_manager = new.env(), design_matrix = kit174Design(), config_list = list()
    ))
    kit174Quiet(runOmicNormHandler(
        list(ruv_mode = "automatic"), workflow, list(), list2env(list()),
        \(text) seen$log <- c(seen$log, text), trx
    ))
    expect_true(seen$ran)
    expect_identical(kit174Count(seen, batch)[["warnings"]], 0L)
})
