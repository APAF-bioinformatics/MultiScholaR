# KIT-126: classic proteomics RUV with the factors in the model (KIT-124's gate, arm 7 of the simulator S1,
# `docs/qa/workflow-execution/kit/kit-124/sim/s1_classic.R`): groups as replicate sets, every protein a control, W
# computed on the proteins observed in every sample, the test on the normalised matrix with W as covariates. Today's
# design (the corrected matrix tested) stays the default and is not changed.

# The simulator's centre cell at a smaller size: 20% effects of 1 log2, one factor of strength 3 orthogonal to the
# groups, intensity-dependent missingness, proteins kept when half of each group is observed.
ruvModelFixture <- function(seed = 20261118L, proteins = 300L, per_group = 6L, missing = 0.25) {
    set.seed(seed)
    groups <- rep(c("control", "treated"), each = per_group)
    n_samples <- length(groups)
    mu <- stats::rnorm(proteins, 22, 2)
    s <- sqrt(4 * 0.25^2 / stats::rchisq(proteins, 4))
    effect <- numeric(proteins)
    where <- sample.int(proteins, round(0.2 * proteins))
    effect[where] <- sample(rep(c(1, -1), length.out = length(where)))
    scale1 <- function(v) { v <- v - mean(v); v / stats::sd(v) }
    w <- scale1(stats::lm.fit(stats::model.matrix(~ factor(groups)), stats::rnorm(n_samples))$residuals)
    tau <- sqrt(3 * 4 * 0.25^2 / stats::qchisq(0.5, 4))
    x <- matrix(mu, proteins, n_samples) + outer(effect, as.numeric(groups == "treated")) +
        outer(stats::rnorm(proteins, 0, tau), w) + matrix(stats::rnorm(proteins * n_samples), proteins, n_samples) * s
    if (missing > 0) {
        lod <- stats::uniroot(function(l) mean(stats::plogis(l - x)) - missing, c(0, 60))$root
        x[stats::runif(length(x)) < stats::plogis(lod - x)] <- NA_real_
    }
    dimnames(x) <- list(sprintf("PG_%05d", seq_len(proteins)), sprintf("run_%03d", seq_len(n_samples)))
    keep <- rowSums(is.finite(x[, groups == "control"])) >= per_group / 2 &
        rowSums(is.finite(x[, groups == "treated"])) >= per_group / 2
    design <- data.frame(Run = colnames(x), group = groups, replicates = seq_len(n_samples),
        tech_rep_group = paste0("u", seq_len(n_samples)), stringsAsFactors = FALSE)
    list(x = x[keep, , drop = FALSE], design = design)
}

ruvModelObject <- function(x, design, technical = "tech_rep_group") {
    quant <- as.data.frame(c(list(Protein.Group = rownames(x)), as.data.frame(x)), check.names = FALSE, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = quant, protein_id_column = "Protein.Group",
        protein_id_table = data.frame(Protein.Group = rownames(x), stringsAsFactors = FALSE), design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = technical,
        args = list(globalParameters = list(workflow_type = "LFQ")))
}

ruvModelContrast <- "treated_vs_control=grouptreated-groupcontrol"

# The package's path at the kernel: the factors, the corrected matrix (the same controls, k and replicate sets), the
# corrected object carrying the factors, then the DA entry point.
ruvModelPackageState <- function(fixture, k, k_source = "user", technical = "tech_rep_group") {
    normalised <- ruvModelObject(fixture$x, fixture$design, technical)
    factors <- computeProtRuvModelFactors(normalised, k = k, k_source = k_source, ruv_grouping_variable = "group")
    corrected <- suppressMessages(ruvIII_C_Varying(protRuvGroupReplicateObject(normalised), ruv_grouping_variable = "group",
        ruv_number_k = factors$k, ctrl = factors$controls))
    list(state = attachProtRuvModelFactors(corrected, normalised, factors), factors = factors, corrected = corrected)
}

ruvModelDa <- function(object) {
    result <- NULL
    utils::capture.output(result <- suppressWarnings(suppressMessages(differentialAbundanceAnalysis(object,
        contrasts_tbl = data.frame(contrasts = ruvModelContrast, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))))
    result
}

# The limma kernel as the default path calls it (the simulator's `da`): `covariates` join the design and the formula.
ruvModelKernel <- function(matrix, design, covariates = NULL) {
    rownames(design) <- design$Run
    formula <- "~ 0 + group"
    if (!is.null(covariates)) {
        design <- cbind(design, as.data.frame(covariates[rownames(design), , drop = FALSE]))
        formula <- paste(c(formula, colnames(covariates)), collapse = " + ")
    }
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(runTestsContrasts(matrix, contrast_strings = ruvModelContrast,
        design_matrix = design, formula_string = formula, weights = NA, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
        eBayes_robust = TRUE, group_id = "group", technical_replicate_id = "tech_rep_group"))))
    value
}

# Arm 7 as the simulator computes it (`in_model(all_controls, k, "no_tech")`): W of RUV-III on the complete proteins,
# every one a control, the groups as replicate sets; then the kernel on the uncorrected matrix with W.
ruvModelArm7 <- function(fixture, k, rows = rownames(fixture$x)) {
    x <- fixture$x
    design <- fixture$design
    design$no_tech <- paste0("i", seq_len(nrow(design)))
    complete <- rowSums(!is.finite(x)) == 0
    value <- suppressMessages(computeRuvCovariates(ruvModelObject(x[complete, , drop = FALSE], design, "no_tech"), k = k,
        controls = rep(TRUE, sum(complete)), ruv_grouping_variable = "group"))
    list(W = value$W, test = ruvModelKernel(x[rows, , drop = FALSE], fixture$design, value$W))
}

test_that("the factors-in-the-model path reproduces the simulator's arm 7 on a fixture with missing values", {
    fixture <- ruvModelFixture()
    expect_gt(mean(!is.finite(fixture$x)), 0.1)
    for (k in 1:2) {
        package <- ruvModelPackageState(fixture, k)
        state <- package$state
        ids <- state@protein_quant_table$Protein.Group
        # every protein of the normalised matrix is in the corrected state here, so the test is arm 7's own
        expect_identical(ids, rownames(fixture$x))
        result <- ruvModelDa(state)
        reference <- ruvModelArm7(fixture, k, rows = ids)
        expect_equal(unname(package$factors$W), unname(reference$W), tolerance = 1e-12)
        table <- result$contrasts_results$results[[1L]]
        expected <- reference$test$results[[1L]]
        expect_identical(rownames(table), rownames(expected))
        for (column in c("logFC", "raw_pvalue", "fdr_qvalue")) {
            expect_equal(table[[column]], expected[[column]], tolerance = 1e-12, info = paste(column, "k =", k))
        }
        # the model that ran is named on the result
        w_terms <- paste0("W", seq_len(k))
        expect_identical(result$da_model$model, "limma")
        expect_identical(result$da_model$formula_used, paste(c("~ 0 + group", w_terms), collapse = " + "))
        expect_identical(result$da_model$w_terms, w_terms)
        expect_identical(result$da_model$k, k)
        expect_identical(result$da_model$ruv_mode, "covariates")
        expect_identical(result$da_model$k_source, "user")
        expect_identical(result$da_model$ruv_controls, "all")
        expect_identical(result$da_model$tested_matrix, "uncorrected")
        # it is not the corrected matrix that was tested, and not the uncorrected matrix without W
        on_corrected <- ruvModelKernel(as.matrix(tibble::column_to_rownames(package$corrected@protein_quant_table, "Protein.Group")),
            fixture$design)$results[[1L]]
        without_w <- ruvModelKernel(fixture$x, fixture$design)$results[[1L]]
        expect_gt(max(abs(table$raw_pvalue - on_corrected$raw_pvalue), na.rm = TRUE), 1e-6)
        expect_gt(max(abs(table$raw_pvalue - without_w$raw_pvalue), na.rm = TRUE), 1e-6)
        # proteins with missing values are tested, on the samples they have
        partial <- rownames(fixture$x)[rowSums(!is.finite(fixture$x)) > 0]
        expect_gt(length(partial), 10L)
        expect_true(all(is.finite(table[partial, "raw_pvalue"])))
    }
})

test_that("the corrected state carries the corrected table, the factors and the normalised matrix", {
    fixture <- ruvModelFixture(proteins = 120L)
    package <- ruvModelPackageState(fixture, 2L)
    state <- package$state
    samples <- fixture$design$Run
    expect_identical(state@protein_quant_table, package$corrected@protein_quant_table)
    expect_identical(state@args$ruv_uncorrected_quant, fixture$x)
    expect_gt(max(abs(as.matrix(state@protein_quant_table[, samples]) - fixture$x[state@protein_quant_table$Protein.Group, samples]),
        na.rm = TRUE), 1e-3)
    expect_identical(unname(as.matrix(state@design_matrix[, c("W1", "W2")])), unname(package$factors$W[samples, ]))
    expect_identical(state@technical_replicate_id, "tech_rep_group")
    declared <- state@args$ruvCovariates
    expect_identical(declared[c("ruv_mode", "ruv_design", "k", "k_source", "k_estimated", "k_bounded", "ruv_controls",
        "replicate_structure", "ruv_grouping_variable", "tested_matrix")],
        list(ruv_mode = "covariates", ruv_design = "factors_in_model", k = 2L, k_source = "user", k_estimated = 2L,
            k_bounded = FALSE, ruv_controls = "all", replicate_structure = "group", ruv_grouping_variable = "group",
            tested_matrix = "uncorrected"))
    complete <- sum(rowSums(!is.finite(fixture$x)) == 0)
    expect_identical(declared$control_count, complete)
    expect_identical(declared$protein_count, nrow(fixture$x))
})

test_that("W comes from the proteins observed in every sample, with the groups as replicate sets", {
    fixture <- ruvModelFixture(proteins = 120L)
    complete <- rowSums(!is.finite(fixture$x)) == 0
    expect_true(any(!complete))
    factors <- computeProtRuvModelFactors(ruvModelObject(fixture$x, fixture$design), k = 2L, k_source = "user")
    reference <- suppressMessages(computeRuvCovariates(ruvModelObject(fixture$x[complete, ], fixture$design), k = 2L,
        controls = rep(TRUE, sum(complete)), ruv_grouping_variable = "group"))
    expect_equal(factors$W, reference$W, tolerance = 1e-12)
    expect_identical(dimnames(factors$W), list(fixture$design$Run, c("W1", "W2")))
    expect_identical(factors$controls, stats::setNames(rep(TRUE, nrow(fixture$x)), rownames(fixture$x)))
    expect_identical(factors$replicate_structure, "group")
    # declared technical replicates do not change the replicate sets of this design (arm 7; arms 15/16 are another option)
    paired <- fixture$design
    paired$tech_rep_group <- paste0(paired$group, "_", rep(1:3, each = 2L, times = 2L))
    declared <- ruvModelObject(fixture$x, paired)
    expect_identical(resolveProtRuvReplicateBlocks(paired, "Run", "group", "tech_rep_group", "group")$structure, "technical_replicates")
    expect_equal(computeProtRuvModelFactors(declared, k = 2L, k_source = "user")$W, factors$W, tolerance = 1e-12)
})

test_that("k is a parameter: the sweep's estimate is bounded and recorded, a user's k is refused above the bound", {
    fixture <- ruvModelFixture(proteins = 120L)
    object <- ruvModelObject(fixture$x, fixture$design)
    k_max <- 12L - 2L - 1L
    expect_gt(sum(rowSums(!is.finite(fixture$x)) == 0), k_max + 1L)
    bounded <- computeProtRuvModelFactors(object, k = 40L, k_source = "sweep")
    expect_identical(bounded[c("k", "k_source", "k_estimated", "k_bounded", "k_max")],
        list(k = k_max, k_source = "sweep", k_estimated = 40L, k_bounded = TRUE, k_max = k_max))
    expect_identical(ncol(bounded$W), k_max)
    kept <- computeProtRuvModelFactors(object, k = 3L, k_source = "sweep")
    expect_identical(kept[c("k", "k_estimated", "k_bounded")], list(k = 3L, k_estimated = 3L, k_bounded = FALSE))
    expect_identical(ncol(kept$W), 3L)
    expect_error(computeProtRuvModelFactors(object, k = 40L, k_source = "user"), class = "prot_ruv_factors_unavailable")
    expect_error(computeProtRuvModelFactors(object, k = 0L, k_source = "user"), "invalid")
    expect_error(computeProtRuvModelFactors(object, k = 1L, k_source = "num.sv"))
    # the bound follows the complete proteins when they are the scarce side
    scarce <- fixture$x
    complete <- which(rowSums(!is.finite(scarce)) == 0)
    scarce[complete[-(1:4)], 1L] <- NA_real_
    few <- computeProtRuvModelFactors(ruvModelObject(scarce, fixture$design), k = 8L, k_source = "sweep")
    expect_identical(few[c("k", "k_bounded", "k_max", "control_count")], list(k = 3L, k_bounded = TRUE, k_max = 3L, control_count = 4L))
    # no factor can be computed: refused with the reason, never a silent fallback
    scarce[complete[-1L], 1L] <- NA_real_
    expect_error(computeProtRuvModelFactors(ruvModelObject(scarce, fixture$design), k = 1L, k_source = "sweep"),
        "observed in every sample", class = "prot_ruv_factors_unavailable")
    tiny <- ruvModelObject(fixture$x[, c(1L, 2L, 7L)], fixture$design[c(1L, 2L, 7L), ])
    expect_error(computeProtRuvModelFactors(tiny, k = 1L, k_source = "user"), class = "prot_ruv_factors_unavailable")
})

test_that("an object without the declaration is tested exactly as today", {
    fixture <- ruvModelFixture(proteins = 120L)
    normalised <- ruvModelObject(fixture$x, fixture$design)
    ctrl <- stats::setNames(rep(TRUE, nrow(fixture$x)), rownames(fixture$x))
    legacy <- suppressMessages(ruvIII_C_Varying(normalised, ruv_grouping_variable = "group", ruv_number_k = 1L, ctrl = ctrl))
    expect_null(legacy@args$ruv_uncorrected_quant)
    expect_null(legacy@args$ruvCovariates)
    result <- ruvModelDa(legacy)
    today <- ruvModelKernel(as.matrix(tibble::column_to_rownames(legacy@protein_quant_table, "Protein.Group")), legacy@design_matrix)
    expect_identical(result$contrasts_results, today)
    expect_identical(result$da_model, list(model = "limma", formula_used = "~ 0 + group", w_terms = character(), k = 0L,
        ruv_mode = "none", treat_applied = FALSE, treat_lfc_cutoff_requested = 0))
    # W columns alone select nothing: only the stored normalised matrix does (a DIA_LIMPA state on the limma path)
    with_w <- legacy
    with_w@design_matrix$W1 <- seq_len(nrow(with_w@design_matrix))
    with_w@args$ruvCovariates <- list(k = 1L, ruv_mode = "covariates")
    w_result <- ruvModelDa(with_w)
    expect_identical(w_result$contrasts_results, today)
    expect_identical(w_result$da_model$formula_used, "~ 0 + group")
})

test_that("an inconsistent factors-in-the-model object is refused, not tested", {
    fixture <- ruvModelFixture(proteins = 120L)
    state <- ruvModelPackageState(fixture, 2L)$state
    refused <- function(object) expect_error(resolveProtRuvModelInput(object, "~ 0 + group"), class = "prot_ruv_model_invalid")
    expect_null(resolveProtRuvModelInput(ruvModelObject(fixture$x, fixture$design), "~ 0 + group"))
    accepted <- resolveProtRuvModelInput(state, "~ 0 + group")
    expect_identical(accepted$formula_used, "~ 0 + group + W1 + W2")
    expect_identical(accepted$matrix, fixture$x[state@protein_quant_table$Protein.Group, fixture$design$Run])
    no_w <- state; no_w@design_matrix$W2 <- NULL
    refused(no_w)
    other_k <- state; other_k@args$ruvCovariates$k <- 3L
    refused(other_k)
    corrected_mode <- state; corrected_mode@args$ruvCovariates$ruv_mode <- "correction"
    refused(corrected_mode)
    undeclared <- state; undeclared@args$ruvCovariates <- NULL
    refused(undeclared)
    uncovered <- state; uncovered@args$ruv_uncorrected_quant <- fixture$x[-1L, ]
    refused(uncovered)
    missing_w <- state; missing_w@design_matrix$W1[[2L]] <- NA_real_
    refused(missing_w)
    expect_error(ruvModelDa(no_w), class = "prot_ruv_model_invalid")
    # samples removed after the RUV step (the correlation filter) are dropped from the test with their W rows
    fewer <- state
    fewer@design_matrix <- fewer@design_matrix[-1L, ]
    fewer@protein_quant_table$run_001 <- NULL
    expect_identical(colnames(resolveProtRuvModelInput(fewer, "~ 0 + group")$matrix), fixture$design$Run[-1L])
})

test_that("the normalisation declaration is today's by default and declares the design only when it is chosen", {
    today <- c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max", "separation_metric",
        "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k", "grouping")
    default <- diaNormalizationOperation()
    expect_identical(default$parameters, today)
    expect_identical(default[c("name", "input_class", "fields", "unit", "output_class", "product_type", "factors")],
        list(name = "normalization-ruv", input_class = "ProteinQuantitativeData",
            fields = c("protein_quant_table", "design_matrix", "protein_id_table"), unit = "whole_matrix",
            output_class = "ProteinQuantitativeData", product_type = "dia-normalization-ruv", factors = list()))
    expect_identical(diaNormalizationOperation("legacy")$parameters, today)
    chosen <- diaNormalizationOperation("factors_in_model")
    expect_identical(chosen$parameters, c(today, "ruv_design"))
    expect_identical(chosen$name, "normalization-ruv")
    expect_error(diaNormalizationOperation("subtract"))
    expect_identical(workflowClassDeclaration("ProteinQuantitativeData")$args_as_owned,
        c("limpa_dpc_quant_results", "ruv_uncorrected_quant"))
})

test_that("through the operations: the RUV state carries the factors and the DA reads the normalised matrix", {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "kit-126", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("protein_replicate_filtered", object, list(), "source")
    source <- getDIAStateReference(manager)
    contrasts <- normaliseProtDaContrastsTable(data.frame(contrasts = "grouptreated-groupcontrol"))
    parameters <- list(norm_method = "quantile", ruv_mode = "manual", auto_percentage_min = NULL,
        auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL,
        adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 2, grouping = "group")
    outcome <- NULL
    utils::capture.output(outcome <- suppressWarnings(suppressMessages(withDIAStateScope(manager, function() {
        contrast_on <- function(final) {
            spec <- readDIAStateMetadata(manager, final)$data
            runSavedOperation(manager, diaContrastOperation(), final, list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0,
                weights = NA, contrast = contrasts$full_format[[1L]], annotations = spec$components$protein_id_table$entity_id,
                name = contrasts$contrasts[[1L]], comparison = contrasts$friendly_names[[1L]]))
        }
        legacy <- runSavedOperation(manager, diaNormalizationOperation(), source, parameters)
        expect_error(runSavedOperation(manager, diaNormalizationOperation("factors_in_model"), source, parameters),
            class = "workflow_operation_parameters")
        chosen <- runSavedOperation(manager, diaNormalizationOperation("factors_in_model"), source,
            c(parameters, list(ruv_design = "factors_in_model")))
        again <- runSavedOperation(manager, diaNormalizationOperation(), source, parameters)
        list(legacy = legacy, chosen = chosen, again = again,
            legacy_state = readDIAExecutionComponents(manager, legacy$value$final),
            chosen_state = readDIAExecutionComponents(manager, chosen$value$final),
            chosen_record = readDIAStateMetadata(manager, chosen$value$final)$data,
            normalised = readDIAExecutionComponents(manager, chosen$value$normalized),
            legacy_contrast = contrast_on(legacy$value$final), chosen_contrast = contrast_on(chosen$value$final))
    }))))
    # today's declaration: computed under its own identity, hit again after the new design ran, no declaration on its state
    expect_identical(c(outcome$legacy$status, outcome$chosen$status, outcome$again$status), c("computed", "computed", "hit"))
    expect_identical(outcome$again$reference, outcome$legacy$reference)
    expect_false(identical(outcome$chosen$reference, outcome$legacy$reference))
    expect_identical(outcome$chosen$value$normalized, outcome$legacy$value$normalized)
    expect_null(outcome$legacy_state@args$ruvCovariates)
    expect_null(outcome$legacy_state@args$ruv_uncorrected_quant)
    expect_false(any(grepl("^W[0-9]+$", names(outcome$legacy_state@design_matrix))))
    expect_identical(outcome$legacy_contrast$value$model$formula_used, "~ 0 + group")
    expect_identical(outcome$legacy_contrast$value$model$ruv_mode, "none")
    # the new design: k and its source, every protein a control, the normalised matrix stored out of line and reattached
    state <- outcome$chosen_state
    expect_identical(outcome$chosen$value$best_k, 2L)
    expect_true(all(outcome$chosen$value$controls))
    expect_identical(outcome$chosen$value$ruv_model[c("ruv_design", "k", "k_source", "ruv_controls", "replicate_structure")],
        list(ruv_design = "factors_in_model", k = 2L, k_source = "user", ruv_controls = "all", replicate_structure = "group"))
    expect_true(isSavedWorkflowReference(outcome$chosen_record$slots$args$ruv_uncorrected_quant))
    samples <- as.character(object@design_matrix$Run)
    normalised <- as.matrix(tibble::column_to_rownames(outcome$normalised@protein_quant_table, "Protein.Group"))
    expect_equal(state@args$ruv_uncorrected_quant[, samples], normalised[, samples], tolerance = 1e-12)
    expect_identical(grep("^W[0-9]+$", names(state@design_matrix), value = TRUE), c("W1", "W2"))
    expect_identical(state@technical_replicate_id, "tech_rep_group")
    expect_gt(max(abs(as.matrix(state@protein_quant_table[, samples]) -
        normalised[state@protein_quant_table$Protein.Group, samples]), na.rm = TRUE), 1e-6)
    model <- outcome$chosen_contrast$value$model
    expect_identical(model[c("model", "formula_used", "w_terms", "k", "ruv_mode", "k_source", "ruv_controls", "tested_matrix")],
        list(model = "limma", formula_used = "~ 0 + group + W1 + W2", w_terms = c("W1", "W2"), k = 2L, ruv_mode = "covariates",
            k_source = "user", ruv_controls = "all", tested_matrix = "uncorrected"))
    expect_identical(outcome$chosen_contrast$status, "computed")
    expect_gt(outcome$chosen_contrast$value$rows, 0L)
})
