# KIT-129 W3: the control rule and the k rule of the DIA_LIMPA RUV covariates step (`protein_ruv_covariates`).
# Today's step (the limpa script's sweep gives the controls and k) is the default and is not changed. The declared
# rules are measured arms of KIT-124's simulator S2 (`docs/qa/workflow-execution/kit/kit-124/sim/s2_limpa.R`):
#   arm 7  all proteins as controls, k = E (`sva::num.sv`), the factors as covariates of `limpa::dpcDE`;
#   arm 16 the same on a design with declared technical replicates (they are RUV-III's replicate sets on this route);
#   arm 8  a supplied control list, k = E;
#   arm 5  the proteins outside the 25% most significant of a first test, the sweep's k.
# Each is reproduced here on seeded data from the simulator's own lines.

# The data of simulator S2 at a smaller size (`simulateData` and `quantify` of s2_limpa.R, the centre cell's terms):
# limpa's generator, group effects of +-0.5 log2 per group, unwanted factors with loadings N(0, strength * 0.09) on the
# complete precursor matrix, missingness redrawn with the generator's DPC coefficients, rows permuted, then the
# route's quantification.
limpaRuvSimulate <- function(seed, proteins = 200L, per_group = 6L, groups = 2L, factors = 1L, strength = 3, units = 0L) {
    set.seed(seed)
    sim <- limpa::simProteinDataSet(n.peptides = proteins * 4L, n.groups = groups, samples.per.group = per_group,
        peptides.per.protein = 4L, mu.range = c(2, 10), prop.de = 0, dpc.slope = 0.7, prop.missing = 0.25)
    y <- sim$other$E.complete
    samples <- ncol(y)
    ids <- unique(sim$genes$Protein)
    protein_of <- match(sim$genes$Protein, ids)
    group_index <- sim$targets$Group
    contrast <- c(-0.5, 0.5, 0)[group_index]
    n_effects <- round(0.2 * proteins)
    delta <- numeric(proteins)
    delta[sample.int(proteins, n_effects)] <- rep(c(1, -1), c(ceiling(0.5 * n_effects), n_effects - ceiling(0.5 * n_effects)))
    term <- outer(delta, contrast)
    if (factors > 0L) {
        scores <- matrix(stats::rnorm(samples * factors), samples)
        loadings <- sqrt(strength * 0.09) * matrix(stats::rnorm(proteins * factors), proteins)
        term <- term + loadings %*% t(scores)
    }
    within_group <- stats::ave(group_index, group_index, FUN = seq_along)
    replicates <- within_group
    if (units > 0L) {
        replicates <- ceiling(within_group / (per_group / units))
        unit <- (group_index - 1L) * units + replicates
        term <- term + matrix(stats::rnorm(proteins * max(unit), 0, 0.3), proteins)[, unit]
    }
    y <- y + term[protein_of, , drop = FALSE]
    probability <- stats::plogis(sim$dpc[[1L]] + y * sim$dpc[[2L]], lower.tail = FALSE)
    y[stats::rbinom(length(y), prob = probability, size = 1L) == 1L] <- NA
    order <- sample.int(nrow(y))
    runs <- sprintf("run_%03d", seq_len(samples))
    design <- data.frame(Run = runs, group = c("A", "B", "C")[group_index], replicates = replicates, stringsAsFactors = FALSE)
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    list(y = y[order, , drop = FALSE], protein = sim$genes$Protein[order], precursor = rownames(y)[order], design = design)
}

limpaRuvQuantify <- function(sim, engine) {
    y <- sim$y
    dimnames(y) <- list(paste(sim$protein, sim$precursor, sep = "%"), sim$design$Run)
    object <- methods::new("PeptideQuantitativeData")
    object@peptide_data <- data.frame(Protein.Ids = sim$protein, Stripped.Sequence = sim$precursor, stringsAsFactors = FALSE)
    object@peptide_matrix <- y
    object@protein_id_column <- "Protein.Ids"
    object@peptide_sequence_column <- "Stripped.Sequence"
    object@is_logged_data <- TRUE
    object@design_matrix <- sim$design
    object@sample_id <- "Run"
    object@group_id <- "group"
    object@technical_replicate_id <- "tech_rep_group"
    object@args <- list(globalParameters = list(workflow_type = "DIA_LIMPA"))
    counts <- tibble::tibble(Protein.Ids = unique(sim$protein), peptide_count = 4L, peptidoform_count = 4L)
    quantified <- NULL
    utils::capture.output(quantified <- suppressWarnings(suppressMessages(proteinMissingValueImputationLimpa(object,
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = engine, peptide_summary = counts))))
    quantified
}

limpaRuvEngines <- function() Filter(function(engine) requireNamespace(engine, quietly = TRUE), c("limpa", "fastlimpadpc"))

limpaRuvContrast <- "B_vs_A=groupB-groupA"

limpaRuvQuiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

# The route's test as the simulator calls it (`routeTest` of s2_limpa.R on `sim-k124-exact`): `covariates` join the
# design and the formula; `exact` is the exactness round's arm 25 (decision M29): the comparison as the first
# coefficient of the reparameterised model (`limma::contrastAsCoef`), the route's dpcDE and converter called with it.
limpaRuvRouteTest <- function(elist, design, covariates = NULL, exact = FALSE) {
    formula_used <- "~ 0 + group"
    if (!is.null(covariates)) {
        design <- cbind(design[setdiff(names(design), colnames(covariates))],
            as.data.frame(covariates)[rownames(design), , drop = FALSE])
        formula_used <- paste("~ 0 + group", "+", paste(colnames(covariates), collapse = " + "))
    }
    model <- stats::model.matrix(stats::as.formula(formula_used), design)
    samples <- rownames(model)
    replicate <- resolveProtTechnicalReplicateBlock(design[samples, , drop = FALSE], "group", "tech_rep_group")
    block <- if (isTRUE(replicate$technical)) replicate$block else NULL
    contrast_string <- limpaRuvContrast
    if (exact) {
        contrast <- limma::makeContrasts(contrasts = sub("^.*=", "", limpaRuvContrast), levels = colnames(model))
        model <- limma::contrastAsCoef(model, contrast)$design
        dimnames(model) <- list(samples, c("comparison", paste0("Q", seq_len(ncol(model))[-1L])))
        contrast_string <- sub("=.*$", "=comparison", limpaRuvContrast)
    }
    fit <- limpaRuvQuiet(limpa::dpcDE(elist[, samples], model, plot = FALSE, block = block))
    result <- limpaRuvQuiet(convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = contrast_string,
        design_matrix = model, eBayes_trend = TRUE, eBayes_robust = TRUE))
    table <- result$results[[1L]]
    data.frame(id = rownames(table), logFC = table$logFC, p = table$raw_pvalue, q = table$fdr_qvalue, stringsAsFactors = FALSE)
}

# One quantified data set with what the simulator's arms read from it.
limpaRuvCase <- function(engine, ...) {
    sim <- limpaRuvSimulate(...)
    quantified <- limpaRuvQuantify(sim, engine)
    design <- sim$design
    rownames(design) <- design$Run
    elist <- quantified@args$limpa_dpc_quant_results$quantified_elist
    stopifnot(!anyNA(elist$E), identical(colnames(elist$E), design$Run),
        identical(rownames(elist$E), quantified@protein_quant_table[["Protein.Ids"]]))
    list(quantified = quantified, elist = elist, design = design, ids = rownames(elist$E))
}

# E as the exactness round computes it (arms 31 and 32: 200 permutations, decision M30).
limpaRuvSimulatorE <- function(case, seed) {
    as.integer(sva::num.sv(case$elist$E, stats::model.matrix(~ group, case$design), method = "be", B = 200, seed = seed))
}

# The simulator's `inModel`: the factors of `computeRuvCovariates`, then the route's test with them, with the exact
# standard error of the declared rules' test (decision M29, arm 25).
limpaRuvInModel <- function(case, k, controls) {
    factors <- limpaRuvQuiet(computeRuvCovariates(case$quantified, k, controls, "group"))
    list(W = factors$W, replicate_structure = factors$replicate_structure,
        table = limpaRuvRouteTest(case$elist, case$design, factors$W, exact = TRUE))
}

# The package: the declared kernel on the quantified object (the store keeps what it checkpoints), and the DA entry
# point on the state the kernel returns as final.
limpaRuvKernel <- function(quantified, parameters, ruv_design = "declared") {
    saved <- list()
    store <- list(checkpoint = function(value, label, from = NULL) {
        saved[[label]] <<- value
        label
    })
    value <- diaRuvCovariatesOperation(ruv_design)$kernel(quantified, parameters, store)
    list(value = value, state = saved[[value$final]], saved = saved)
}

limpaRuvDa <- function(state) {
    limpaRuvQuiet(differentialAbundanceAnalysis(state,
        contrasts_tbl = data.frame(contrasts = limpaRuvContrast, stringsAsFactors = FALSE), formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
}

expectLimpaRuvTable <- function(result, reference, info) {
    table <- result$contrasts_results$results[[1L]]
    expect_identical(rownames(table), reference$id, info = info)
    # 1e-9 (decision M29): the package places the reparameterised fit's coefficient into the original dpcDE fit, the
    # simulator tests the reparameterised fit itself; dpcDE's iterated weights agree between the two to rounding only
    # (measured: 1e-11 relative on p with k = 2).
    expect_equal(table$logFC, reference$logFC, tolerance = 1e-9, info = info)
    expect_equal(table$raw_pvalue, reference$p, tolerance = 1e-9, info = info)
    expect_equal(table$fdr_qvalue, reference$q, tolerance = 1e-9, info = info)
}

limpaRuvSweep <- function(quantified) {
    limpaRuvQuiet(findBestNegCtrlPercentageCapped(quantified, percentage_range = seq(1, 30, by = 1),
        num_components_to_impute = 5, ruv_grouping_variable = "group", ruv_qval_cutoff = 0.05, ruv_fdr_method = "BH",
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, verbose = FALSE))
}

limpaRuvSweepParameters <- c("percentage_min", "percentage_max", "num_components_to_impute", "ruv_grouping_variable",
    "ruv_qval_cutoff", "ruv_fdr_method", "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_mode")
limpaRuvRuleParameters <- c("ruv_control_rule", "ruv_k_rule", "ruv_k", "ruv_control_ids", "ruv_k_seed")

test_that("the defaults are today's declaration and today's parameters; the rules are declared only when chosen", {
    expect_identical(names(resolveDIARuvCovariatesParameters(list(), "group")), limpaRuvSweepParameters)
    expect_identical(names(resolveDIARuvCovariatesParameters(list(), "group", "correction")), limpaRuvSweepParameters)
    expect_identical(resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "sweep_least_significant"),
        resolveDIARuvCovariatesParameters(list(), "group"))
    legacy <- diaRuvCovariatesOperation()
    expect_identical(legacy$parameters, limpaRuvSweepParameters)
    named <- diaRuvCovariatesOperation("legacy")
    expect_identical(named[setdiff(names(named), "kernel")], legacy[setdiff(names(legacy), "kernel")])
    declared <- diaRuvCovariatesOperation("declared")
    expect_identical(declared$parameters, c(limpaRuvSweepParameters, limpaRuvRuleParameters))
    expect_identical(declared[setdiff(names(declared), c("kernel", "parameters"))], legacy[setdiff(names(legacy), c("kernel", "parameters"))])
    # each design's pre-selected k rule is the one its figure was measured with
    recommended <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features")
    expect_identical(names(recommended), c(limpaRuvSweepParameters, limpaRuvRuleParameters))
    expect_identical(recommended[limpaRuvSweepParameters], resolveDIARuvCovariatesParameters(list(), "group"))
    expect_identical(recommended[limpaRuvRuleParameters], list(ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be",
        ruv_k = NULL, ruv_control_ids = NULL, ruv_k_seed = RUV_FACTOR_COUNT_SEED))
    expect_identical(resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "first_pass_drop_top_25")$ruv_k_rule, "sweep_bounded")
    supplied <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "supplied_list", ruv_control_ids = c("b", "a", "b"),
        ruv_k_rule = "user", ruv_k = 2, ruv_k_seed = 7)
    expect_identical(supplied[limpaRuvRuleParameters], list(ruv_control_rule = "supplied_list", ruv_k_rule = "user", ruv_k = 2L,
        ruv_control_ids = c("a", "b"), ruv_k_seed = 7L))
    # what is not offered is refused
    refused <- list(
        list(ruv_mode = "correction", ruv_control_rule = "all_features"),
        list(ruv_control_rule = "some_rule"),
        list(ruv_control_rule = "all_features", ruv_k_rule = "elbow"),
        list(ruv_control_rule = "all_features", ruv_k_rule = "user"),
        list(ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = 0),
        list(ruv_control_rule = "all_features", ruv_k = 2),
        list(ruv_control_rule = "supplied_list"),
        list(ruv_control_rule = "all_features", ruv_control_ids = "a"),
        list(ruv_k_rule = "num_sv_be"),
        list(ruv_control_rule = "sweep_least_significant", ruv_k = 2)
    )
    for (arguments in refused) {
        expect_error(do.call(resolveDIARuvCovariatesParameters, c(list(list(), "group"), arguments)), class = "workflow_execution_invalid",
            info = paste(names(arguments), collapse = ", "))
    }
})

test_that("all proteins as controls with the estimated count reproduces the simulator's arm 7 on each DPC engine", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    # the sweep is skipped on this path
    testthat::local_mocked_bindings(findBestNegCtrlPercentageCapped = function(...) stop("the sweep is not to run"))
    for (engine in limpaRuvEngines()) {
        for (fixture in list(list(seed = 101L), list(seed = 106L, factors = 2L, strength = 1.5), list(seed = 107L, groups = 3L))) {
            info <- paste(engine, "seed", fixture$seed)
            case <- do.call(limpaRuvCase, c(list(engine), fixture))
            E <- limpaRuvSimulatorE(case, 4242L)
            expect_identical(E, if (identical(fixture$factors, 2L)) 2L else 1L, info = info)
            reference <- limpaRuvInModel(case, E, rep(TRUE, length(case$ids)))
            parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L)
            package <- limpaRuvKernel(case$quantified, parameters)
            value <- package$value
            expect_identical(value$best_k, E, info = info)
            expect_identical(unname(value$controls), rep(TRUE, length(case$ids)), info = info)
            expect_identical(value$control_count, length(case$ids), info = info)
            expect_identical(value[c("ruv_control_rule", "ruv_k_rule", "ruv_k_estimate", "ruv_applied", "replicate_structure")],
                list(ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be", ruv_k_estimate = E, ruv_applied = TRUE,
                    replicate_structure = "group"), info = info)
            expect_equal(value$W, reference$W, tolerance = 1e-12, info = info)
            # the state: the quantified matrix unchanged, W in the design, the declaration on the object
            state <- package$state
            w_terms <- paste0("W", seq_len(E))
            expect_identical(state@protein_quant_table, case$quantified@protein_quant_table, info = info)
            expect_identical(setdiff(names(state@design_matrix), names(case$quantified@design_matrix)), w_terms, info = info)
            expect_identical(state@args$ruvCovariates[c("k", "ruv_mode", "ruv_control_rule", "ruv_k_rule", "ruv_applied")],
                list(k = E, ruv_mode = "covariates", ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be", ruv_applied = TRUE), info = info)
            # the test through the package's DA entry point is arm 7's table
            result <- limpaRuvDa(state)
            expectLimpaRuvTable(result, reference$table, info)
            expect_identical(result$da_model[c("model", "formula_used", "w_terms", "k", "ruv_mode")],
                list(model = "dpcDE", formula_used = paste("~ 0 + group +", paste(w_terms, collapse = " + ")), w_terms = w_terms, k = E,
                    ruv_mode = "covariates"), info = info)
            # ... and not the table without RUV
            expect_gt(max(abs(reference$table$p - limpaRuvRouteTest(case$elist, case$design)$p)), 1e-3)
            # the count's record travels with the product, for the diagnostic plot
            expect_identical(value$ruv_result$k_estimate$k, E, info = info)
            expect_identical(value$ruv_result$k_estimate$seed, 4242L, info = info)
            expect_s3_class(plotRuvFactorCount(value$ruv_result$k_estimate), "ggplot")
            expect_null(value$ruv_result$optimization_results)
        }
    }
})

test_that("declared technical replicates stay RUV-III's replicate sets on this route (arm 16; the count as arm 29, M28)", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 105L, per_group = 8L, units = 4L)
    # Decision M28: on this design the count removes the biological units (the exactness round's arm 29), not the groups
    E <- as.integer(sva::num.sv(case$elist$E, stats::model.matrix(~ tech_rep_group, case$design), method = "be", B = 200,
        seed = 4242L))
    expect_gt(E, 0L)
    reference <- limpaRuvInModel(case, E, rep(TRUE, length(case$ids)))
    expect_identical(reference$replicate_structure, "technical_replicates")
    package <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L))
    expect_identical(package$value$replicate_structure, "technical_replicates")
    expect_identical(package$value$blocks, 8L)
    expect_identical(package$value$best_k, E)
    expect_equal(package$value$W, reference$W, tolerance = 1e-12)
    result <- limpaRuvDa(package$state)
    expectLimpaRuvTable(result, reference$table, engine)
    expect_true(result$da_model$uses_duplicate_correlation)
    # 16 injections of 8 units in 2 groups: the units carry at most 8 - 2 - 1 = 5 factors (M28; 8 by the injections)
    for (k in c(6, 9)) {
        expect_error(computeDiaRuvCovariatesByRule(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
            ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = k)), "technical replicates", class = "ruv_k_above_limit")
    }
    within <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features",
        ruv_k_rule = "user", ruv_k = 5, ruv_k_seed = 4242L))
    expect_identical(within$value$best_k, 5L)
})

test_that("an estimated count of 0 is a run without RUV, with a notice and the reason on the state", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 103L, factors = 0L)
    expect_identical(limpaRuvSimulatorE(case, 4242L), 0L)
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L)
    expect_message(computeDiaRuvCovariatesByRule(case$quantified, parameters), class = "ruv_k_zero")
    package <- limpaRuvKernel(case$quantified, parameters)
    value <- package$value
    expect_identical(value[c("best_k", "ruv_applied", "ruv_k_rule", "ruv_k_estimate")],
        list(best_k = 0L, ruv_applied = FALSE, ruv_k_rule = "num_sv_be", ruv_k_estimate = 0L))
    expect_match(value$ruv_not_applied_reason, "no unwanted factor", fixed = TRUE)
    expect_null(value$W)
    expect_identical(value$final, value$covariates)
    expect_identical(value$adjusted, value$covariates)
    state <- package$state
    expect_identical(state@design_matrix, case$quantified@design_matrix)
    expect_identical(state@protein_quant_table, case$quantified@protein_quant_table)
    expect_identical(state@args$ruvCovariates[c("k", "ruv_mode", "ruv_applied")], list(k = 0L, ruv_mode = "none", ruv_applied = FALSE))
    # the simulator scores such a run as arm 1
    result <- limpaRuvDa(state)
    expectLimpaRuvTable(result, limpaRuvRouteTest(case$elist, case$design), engine)
    expect_identical(result$da_model[c("formula_used", "k", "ruv_mode")], list(formula_used = "~ 0 + group", k = 0L, ruv_mode = "none"))
})

test_that("a supplied control list with the estimated count is the simulator's arm 8 computation", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    set.seed(8L)
    listed <- sample(case$ids, 50L)
    E <- limpaRuvSimulatorE(case, 4242L)
    reference <- limpaRuvInModel(case, E, case$ids %in% listed)
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "supplied_list",
        ruv_control_ids = c(listed, "NOT_IN_THE_TABLE"), ruv_k_seed = 4242L)
    package <- NULL
    expect_warning(package <- limpaRuvKernel(case$quantified, parameters), class = "dia_ruv_control_ids_unmatched")
    value <- package$value
    expect_identical(unname(value$controls), case$ids %in% listed)
    expect_identical(names(value$controls), case$ids)
    expect_identical(value[c("best_k", "control_count", "ruv_control_rule", "ruv_k_rule")],
        list(best_k = E, control_count = 50L, ruv_control_rule = "supplied_list", ruv_k_rule = "num_sv_be"))
    expect_identical(value$ruv_result$controls_unmatched, "NOT_IN_THE_TABLE")
    expect_equal(value$W, reference$W, tolerance = 1e-12)
    expectLimpaRuvTable(limpaRuvDa(package$state), reference$table, engine)
    # the list's controls give other factors than all proteins do
    all_proteins <- limpaRuvInModel(case, E, rep(TRUE, length(case$ids)))
    expect_gt(max(abs(all_proteins$W - reference$W)), 1e-6)
    # the separation rule is defined with a list: it is today's estimate, the sweep's k
    sweep <- limpaRuvSweep(case$quantified)
    separation <- limpaRuvQuiet(computeDiaRuvCovariatesByRule(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
        ruv_control_rule = "supplied_list", ruv_control_ids = listed, ruv_k_rule = "sweep_bounded")))
    expect_identical(separation$k[c("k", "k_rule")], list(k = sweep$best_k, k_rule = "sweep_bounded"))
    expect_identical(unname(separation$controls), case$ids %in% listed)
    # a list that names no protein of the table is refused
    expect_error(computeDiaRuvCovariatesByRule(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
        ruv_control_rule = "supplied_list", ruv_control_ids = c("X1", "X2"))), class = "workflow_execution_invalid")
})

test_that("first-pass nulls with the sweep's k is the simulator's arm 5", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    first <- limpaRuvRouteTest(case$elist, case$design)
    controls <- rank(first$p, ties.method = "first") > 0.25 * nrow(first)
    sweep <- limpaRuvSweep(case$quantified)
    reference <- limpaRuvInModel(case, sweep$best_k, controls)
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "first_pass_drop_top_25")
    package <- limpaRuvQuiet(limpaRuvKernel(case$quantified, parameters))
    value <- package$value
    expect_identical(unname(value$controls), controls)
    expect_identical(value$control_count, sum(controls))
    expect_identical(value[c("best_k", "ruv_control_rule", "ruv_k_rule", "ruv_applied")],
        list(best_k = sweep$best_k, ruv_control_rule = "first_pass_drop_top_25", ruv_k_rule = "sweep_bounded", ruv_applied = TRUE))
    expect_equal(value$W, reference$W, tolerance = 1e-12)
    expectLimpaRuvTable(limpaRuvDa(package$state), reference$table, engine)
    # the sweep gave k only: its own record (its percentage, controls, scores and table) is kept under a name of its
    # own, and is not presented as what this run used (audit F2)
    k_sweep <- value$ruv_result$k_sweep
    expect_identical(k_sweep$optimization_results, sweep$optimization_results)
    expect_identical(k_sweep[c("best_k", "best_percentage", "control_count", "max_acceptable_k")],
        list(best_k = sweep$best_k, best_percentage = sweep$best_percentage,
            control_count = sum(sweep$best_control_genes_index, na.rm = TRUE), max_acceptable_k = sweep$max_acceptable_k))
    expect_null(value$ruv_result$optimization_results)
    expect_null(value$ruv_result$best_cancor_plot)
    expect_true(is.na(value$percentage) && is.na(value$ruv_result$best_percentage))
    expect_true(is.na(value$ruv_result$best_separation_score) && is.na(value$ruv_result$max_acceptable_k))
    expect_identical(value$ruv_controls_selected, sum(controls))
    expect_null(value$ruv_result$k_estimate)
    # the sweep's controls are not this design's
    expect_false(identical(controls, unname(sweep$best_control_genes_index)))
})

test_that("a user's k is used as given; what the data or the installation cannot give is refused or leaves the run without RUV", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    with_k <- function(k) resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = k)
    automatic <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features")
    # the separation rule is not defined with all proteins as controls
    expect_error(computeDiaRuvCovariatesByRule(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
        ruv_control_rule = "all_features", ruv_k_rule = "sweep_bounded")), class = "ruv_k_separation_undefined")
    # 12 samples in 2 groups carry at most 9 factors: a user's k above that is an error, not a reduced k
    expect_error(computeDiaRuvCovariatesByRule(case$quantified, with_k(10)), class = "ruv_k_above_limit")
    # an estimate above it leaves the run without RUV, with a warning that says why
    local({
        testthat::local_mocked_bindings(estimateRuvFactorCount = function(...) list(k = 10L, seed = 1L))
        refused <- NULL
        expect_warning(refused <- limpaRuvKernel(case$quantified, automatic), class = "dia_ruv_not_applied")
        expect_identical(refused$value[c("best_k", "ruv_applied")], list(best_k = 0L, ruv_applied = FALSE))
        expect_match(refused$value$ruv_not_applied_reason, "at most 9", fixed = TRUE)
        expect_identical(refused$state@design_matrix, case$quantified@design_matrix)
        expect_identical(refused$state@args$ruvCovariates$ruv_mode, "none")
    })
    # without sva the automatic rule stops with the reason; a user's k still runs, the missing estimate recorded
    testthat::local_mocked_bindings(ruvFactorCountSvaAvailable = function() FALSE)
    expect_error(computeDiaRuvCovariatesByRule(case$quantified, automatic), class = "ruv_factor_count_needs_sva")
    chosen <- limpaRuvKernel(case$quantified, with_k(2))
    expect_identical(chosen$value[c("best_k", "ruv_k_rule", "ruv_k_estimate", "ruv_applied")],
        list(best_k = 2L, ruv_k_rule = "user", ruv_k_estimate = NA_integer_, ruv_applied = TRUE))
    expect_match(chosen$value$ruv_result$ruv_k_estimate_reason, "sva", fixed = TRUE)
    reference <- limpaRuvInModel(case, 2L, rep(TRUE, length(case$ids)))
    expect_equal(chosen$value$W, reference$W, tolerance = 1e-12)
    expectLimpaRuvTable(limpaRuvDa(chosen$state), reference$table, engine)
    # the declared kernel does not run today's sweep by default, and the declaration and its parameters must agree
    expect_error(limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group")), class = "workflow_execution_invalid")
})

test_that("the declared rules run through the saved operation; today's declaration on the same source is another product", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("correlation_filtered", case$quantified, list(), "fixture")
    manager$setWorkflowType("DIA_LIMPA")
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L)
    in_memory <- limpaRuvKernel(case$quantified, parameters)
    withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager)
        run <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation("declared"), source, parameters))
        expect_identical(run$status, "computed")
        value <- run$value
        expect_identical(value$best_k, in_memory$value$best_k)
        expect_equal(value$W, in_memory$value$W, tolerance = 1e-12)
        expect_identical(names(value$parameters), c(limpaRuvSweepParameters, limpaRuvRuleParameters))
        final <- readDIAExecutionComponents(manager, value$final)
        samples <- as.character(final@design_matrix[[final@sample_id]])
        expect_equal(final@design_matrix$W1, unname(value$W[samples, "W1"]), tolerance = 1e-12)
        expect_identical(final@args$ruvCovariates$ruv_control_rule, "all_features")
        expect_identical(final@args$ruvCovariates$ruv_k_rule, "num_sv_be")
        # the hydrated state gives the same test as the kernel's state in memory
        hydrated <- limpaRuvDa(final)$contrasts_results$results[[1L]]
        expected <- limpaRuvDa(in_memory$state)$contrasts_results$results[[1L]]
        expect_identical(rownames(hydrated), rownames(expected))
        expect_equal(hydrated$raw_pvalue, expected$raw_pvalue, tolerance = 1e-12)
        # the adjusted matrix is still a product of the step
        adjusted <- readDIAExecutionComponents(manager, value$adjusted)
        expect_false(isTRUE(all.equal(as.matrix(adjusted@protein_quant_table[, samples]), as.matrix(final@protein_quant_table[, samples]))))
        # the same request is reused; the commit description is the existing one
        again <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation("declared"), source, parameters))
        expect_identical(again$status, "hit")
        commit <- diaRuvCovariatesCommit(manager, parameters, run)
        expect_identical(commit$state_name, "ruv_covariates")
        # the commit describes what this run used (audit F2): its control rule and count, k and the rule that gave it
        E <- in_memory$value$best_k
        expect_identical(commit$config[c("ruv_mode", "ruv_k", "control_count", "ruv_applied", "ruv_control_rule", "ruv_k_rule",
            "ruv_k_estimate")], list(ruv_mode = "covariates", ruv_k = E, control_count = length(case$ids), ruv_applied = TRUE,
            ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be", ruv_k_estimate = E))
        expect_true(is.na(commit$config$percentage_as_neg_ctrl) && is.na(commit$config$separation_score))
        expect_identical(commit$config$outcomes[c("k", "control_count", "ruv_mode", "ruv_applied", "control_rule", "k_rule", "k_estimate")],
            list(k = E, control_count = length(case$ids), ruv_mode = "covariates", ruv_applied = TRUE, control_rule = "all_features",
                k_rule = "num_sv_be", k_estimate = E))
        expect_match(commit$description, "all_features", fixed = TRUE)
        expect_match(commit$description, "num_sv_be", fixed = TRUE)
        # today's commit carries today's fields only
        legacy_fields <- c("ruv_mode", "ruv_k", "percentage_as_neg_ctrl", "control_count", "separation_score", "composite_score",
            "max_acceptable_k", "replicate_structure", "block_reason", "replicate_blocks", "ruv_grouping_variable", "ruv_applied", "outcomes")
        expect_identical(names(commit$config)[seq_along(legacy_fields)], legacy_fields)
        # today's declaration on the same source: its own identity, its ten parameters, the sweep's controls
        legacy <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation(), source, resolveDIARuvCovariatesParameters(list(), "group")))
        expect_identical(legacy$status, "computed")
        expect_false(identical(legacy$reference$entity_id, run$reference$entity_id))
        expect_identical(names(legacy$value$parameters), limpaRuvSweepParameters)
        expect_lt(legacy$value$control_count, length(case$ids))
        expect_null(legacy$value$ruv_control_rule)
        legacy_commit <- diaRuvCovariatesCommit(manager, resolveDIARuvCovariatesParameters(list(), "group"), legacy)
        expect_identical(names(legacy_commit$config), legacy_fields)
        expect_true(legacy_commit$config$ruv_applied)
        expect_identical(legacy_commit$config$percentage_as_neg_ctrl, legacy$value$percentage)
    })
})

# ---- Audit of 2026-11-18 (KIT-129, register F2, F3, F4, F6, F10) ----

limpaRuvManager <- function(quantified, env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("correlation_filtered", quantified, list(), "fixture")
    manager$setWorkflowType("DIA_LIMPA")
    manager
}

test_that("a run without RUV is recorded as one: product, state, commit, run description and integration input (audit F2)", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 103L, factors = 0L)
    expect_identical(limpaRuvSimulatorE(case, 4242L), 0L)
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L)
    manager <- limpaRuvManager(case$quantified)
    withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager)
        run <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation("declared"), source, parameters))
        value <- run$value
        # the product: no mode, no k, and nothing presented as used
        expect_identical(value[c("ruv_mode", "best_k", "ruv_applied", "ruv_control_rule", "ruv_k_rule", "ruv_k_estimate")],
            list(ruv_mode = "none", best_k = 0L, ruv_applied = FALSE, ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be",
                ruv_k_estimate = 0L))
        expect_identical(value$ruv_result$ruv_mode, "none")
        expect_identical(value$ruv_result$best_k, 0L)
        expect_true(is.na(value$percentage) && is.na(value$control_count) && is.na(value$ruv_result$best_percentage))
        expect_null(value$controls)
        expect_null(value$ruv_result$best_control_genes_index)
        # what the rule selected is a diagnostic under its own name
        expect_identical(value$ruv_controls_selected, length(case$ids))
        # the state's declaration agrees with the product
        final <- readDIAExecutionComponents(manager, value$final)
        expect_identical(final@args$ruvCovariates[c("k", "ruv_mode", "ruv_applied")], list(k = 0L, ruv_mode = "none", ruv_applied = FALSE))
        expect_true(is.na(final@args$ruvCovariates$control_count) && is.na(final@args$ruvCovariates$percentage_as_neg_ctrl))
        # the commit takes them from the product
        commit <- diaRuvCovariatesCommit(manager, parameters, run)
        config <- commit$config
        expect_identical(config[c("ruv_mode", "ruv_k", "ruv_applied", "ruv_control_rule", "ruv_k_rule", "ruv_k_estimate")],
            list(ruv_mode = "none", ruv_k = 0L, ruv_applied = FALSE, ruv_control_rule = "all_features", ruv_k_rule = "num_sv_be",
                ruv_k_estimate = 0L))
        expect_true(is.na(config$percentage_as_neg_ctrl) && is.na(config$control_count))
        expect_match(config$ruv_not_applied_reason, "no unwanted factor", fixed = TRUE)
        expect_identical(config$outcomes[c("k", "ruv_mode", "ruv_applied", "control_rule", "k_rule", "k_estimate")],
            list(k = 0L, ruv_mode = "none", ruv_applied = FALSE, control_rule = "all_features", k_rule = "num_sv_be", k_estimate = 0L))
        expect_true(is.na(config$outcomes$percentage_as_neg_ctrl) && is.na(config$outcomes$control_count))
        expect_match(config$outcomes$not_applied_reason, "no unwanted factor", fixed = TRUE)
        expect_match(commit$description, "RUV not applied", fixed = TRUE)
        expect_match(commit$description, "no unwanted factor", fixed = TRUE)
        expect_no_match(commit$description, "W1", fixed = TRUE)
        limpaRuvQuiet(commitDIAReferenceState(manager, commit$state_name, commit$checkpoint, commit$config, commit$description))
        # the tab's description of the run
        expect_identical(resolveDIARuvRunDescription(manager)[c("ruv_mode", "grouping")], list(ruv_mode = "none", grouping = "group"))
        # the integration input: the quantified matrix, normalised, no RUV method and no adjusted product
        input <- protIntegrationInputSource(manager)
        expect_identical(input$manifest$matrix_role, "normalised")
        expect_null(input$manifest$correction)
        expect_null(input$manifest$source_state$product)
        expect_identical(input$manifest$source_state$name, "ruv_covariates")
        samples <- case$design$Run
        expect_identical(input$matrix$feature_id, case$ids)
        expect_equal(as.matrix(input$matrix[samples]), as.matrix(case$quantified@protein_quant_table[samples]), tolerance = 1e-12)
    })
})

test_that("a refused automatic k is a run without RUV in the product and the commit, as an estimate of 0 is (audit F2)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    testthat::local_mocked_bindings(estimateRuvFactorCount = \(...) list(k = 10L, seed = 1L))
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features")
    refused <- NULL
    expect_warning(refused <- limpaRuvKernel(case$quantified, parameters), class = "dia_ruv_not_applied")
    expect_identical(refused$value[c("ruv_mode", "best_k", "ruv_applied")], list(ruv_mode = "none", best_k = 0L, ruv_applied = FALSE))
    expect_identical(refused$value$ruv_result$ruv_mode, "none")
    expect_true(is.na(refused$value$control_count))
    # the commit's fields come from the product (the state metadata stands in: no store on this path)
    testthat::local_mocked_bindings(readDIAStateMetadata = \(...) list(data = list(products = list())))
    commit <- diaRuvCovariatesCommit(list(current_state = "correlation_filtered"), parameters,
        list(value = refused$value, reference = "product"))
    expect_false(commit$config$ruv_applied)
    expect_identical(commit$config$ruv_mode, "none")
    expect_match(commit$config$ruv_not_applied_reason, "at most 9", fixed = TRUE)
    expect_match(commit$description, "RUV not applied", fixed = TRUE)
})

test_that("today's declaration refuses the rule parameters, on a computed run and on a saved product (audit F3)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    todays <- resolveDIARuvCovariatesParameters(list(), "group")
    declared <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = 1)
    # the kernel of today's declaration: any one rule parameter is refused, never run as the sweep
    local({
        testthat::local_mocked_bindings(findBestNegCtrlPercentageCapped = \(...) stop("the sweep is not to run"))
        refusal <- expect_error(limpaRuvKernel(case$quantified, declared, "legacy"), class = "workflow_execution_invalid")
        expect_s3_class(refusal, "workflow_execution_invalid")
        expect_match(conditionMessage(refusal), "ruv_control_rule", fixed = TRUE)
        for (extra in list(list(ruv_control_rule = "all_features"), list(ruv_k_rule = "num_sv_be"), list(ruv_k = 2L),
            list(ruv_control_ids = "P1"), list(ruv_k_seed = 7L))) {
            expect_error(limpaRuvKernel(case$quantified, c(todays, extra), "legacy"), class = "workflow_execution_invalid",
                info = names(extra))
        }
    })
    # the auditor's case: today's product is saved, then the same source is asked with rule parameters. The identity
    # of today's declaration reads its ten parameters only, so the saved product would be returned (a hit) before
    # any kernel runs; the route's compute half refuses first.
    manager <- limpaRuvManager(case$quantified)
    withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager)
        saved <- limpaRuvQuiet(computeDIAReferenceRuvCovariates(manager, source, todays))
        expect_identical(saved$status, "computed")
        expect_error(computeDIAReferenceRuvCovariates(manager, source, declared), class = "workflow_execution_invalid")
        expect_error(computeDIAReferenceRuvCovariates(manager, source, c(todays, list(ruv_k_seed = 7L))),
            class = "workflow_execution_invalid")
        # today's request is still served from the saved product
        expect_identical(limpaRuvQuiet(computeDIAReferenceRuvCovariates(manager, source, todays))$status, "hit")
    })
})

limpaRuvDaWithFormula <- function(state, formula_string) {
    differentialAbundanceAnalysis(state, contrasts_tbl = data.frame(contrasts = limpaRuvContrast, stringsAsFactors = FALSE),
        formula_string = formula_string, da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE,
        args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
}

test_that("the DA refuses, naming k, a model that leaves no residual degree of freedom (audit F4)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    with_k <- function(k) {
        state <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features",
            ruv_k_rule = "user", ruv_k = k))$state
        state@design_matrix$batch <- rep(c("x", "y", "z"), 4L)
        state
    }
    # k at the limit the groups leave (12 samples, 2 groups, 9 factors: one residual degree of freedom): the DA runs
    at_limit <- with_k(9L)
    result <- limpaRuvQuiet(limpaRuvDaWithFormula(at_limit, "~ 0 + group"))
    expect_identical(result$da_model$k, 9L)
    expect_false(anyNA(result$contrasts_results$results[[1L]]$raw_pvalue))
    # a batch covariate of the user's formula beside the factors: 13 columns, then 12, on 12 samples
    for (k in c(9L, 8L)) {
        state <- if (k == 9L) at_limit else with_k(k)
        refusal <- NULL
        limpaRuvQuiet(refusal <- tryCatch(limpaRuvDaWithFormula(state, "~ 0 + group + batch"), error = \(e) e))
        expect_s3_class(refusal, "omics_ruv_no_residual_df")
        expect_match(conditionMessage(refusal), sprintf("k = %d", k), fixed = TRUE)
        expect_match(conditionMessage(refusal), "12 samples", fixed = TRUE)
    }
    # with room left the same formula runs
    result <- limpaRuvQuiet(limpaRuvDaWithFormula(with_k(2L), "~ 0 + group + batch"))
    expect_identical(result$da_model$formula_used, "~ 0 + group + batch + W1 + W2")
})

test_that("a supplied list matched in part runs with a warning that gives the matched and listed counts (audit F6)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    listed <- c(case$ids[1:2], sprintf("UNKNOWN%02d", 1:48))
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "supplied_list", ruv_control_ids = listed,
        ruv_k_rule = "user", ruv_k = 1)
    package <- NULL
    warned <- expect_warning(package <- limpaRuvKernel(case$quantified, parameters), class = "dia_ruv_control_ids_unmatched")
    expect_identical(warned$matched, 2L)
    expect_identical(warned$listed, 50L)
    expect_match(conditionMessage(warned), "2 of 50", fixed = TRUE)
    expect_identical(package$value[c("best_k", "control_count", "ruv_applied")], list(best_k = 1L, control_count = 2L, ruv_applied = TRUE))
    expect_length(package$value$ruv_result$controls_unmatched, 48L)
    # a list matched in full does not warn
    full <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "supplied_list", ruv_control_ids = case$ids[1:30],
        ruv_k_rule = "user", ruv_k = 1)
    expect_no_warning(limpaRuvKernel(case$quantified, full), class = "dia_ruv_control_ids_unmatched")
})

test_that("with three groups the first pass is the F test of the group differences, not one contrast (audit F8, F10)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 107L, groups = 3L)
    model <- stats::model.matrix(~ 0 + group, case$design)
    fit <- limpaRuvQuiet(limpa::dpcDE(case$elist, model, plot = FALSE))
    pairs <- limma::makeContrasts(groupB - groupA, groupC - groupA, levels = model)
    f_test <- limma::eBayes(limma::contrasts.fit(fit, pairs), trend = TRUE, robust = TRUE)$F.p.value
    first <- limpaRuvQuiet(diaRuvFirstPassPValues(case$quantified, "group"))
    expect_identical(names(first), case$ids)
    expect_equal(unname(first), f_test, tolerance = 1e-12)
    expected <- rank(f_test, ties.method = "first") > 0.25 * length(case$ids)
    package <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
        ruv_control_rule = "first_pass_drop_top_25", ruv_k_rule = "user", ruv_k = 1))
    expect_identical(unname(package$value$controls), expected)
    expect_identical(package$value$control_count, sum(expected))
    # the simulator's arm ranks B against A only: with three groups that is another control set
    single <- limpaRuvRouteTest(case$elist, case$design)
    expect_false(identical(expected, rank(single$p, ties.method = "first") > 0.25 * nrow(single)))
})

test_that("an error of the estimate does not stop a run whose k the user gave (audit F11)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    testthat::local_mocked_bindings(estimateRuvFactorCount = \(...) stop("a plain error of the estimate"))
    chosen <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features",
        ruv_k_rule = "user", ruv_k = 2))
    expect_identical(chosen$value[c("best_k", "ruv_k_rule", "ruv_k_estimate", "ruv_applied")],
        list(best_k = 2L, ruv_k_rule = "user", ruv_k_estimate = NA_integer_, ruv_applied = TRUE))
    expect_match(chosen$value$ruv_result$ruv_k_estimate_reason, "a plain error", fixed = TRUE)
})

test_that("the count's record says what it removed: the biological units with technical replicates, else the groups (audit F1; M28)", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    parameters <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_seed = 4242L)
    technical <- limpaRuvCase(engine, seed = 105L, per_group = 8L, units = 4L)
    estimate <- limpaRuvKernel(technical$quantified, parameters)$value$ruv_result$k_estimate
    expect_identical(estimate[c("removed", "technical_replicates_ignored")], list(removed = "biological units", technical_replicates_ignored = FALSE))
    plain <- limpaRuvCase(engine, seed = 101L)
    estimate <- limpaRuvKernel(plain$quantified, parameters)$value$ruv_result$k_estimate
    expect_identical(estimate[c("removed", "technical_replicates_ignored")], list(removed = "groups", technical_replicates_ignored = FALSE))
    # without technical replicates the count is the simulator's E on the groups, as before
    expect_identical(estimate$k, limpaRuvSimulatorE(plain, 4242L))
})

test_that("the shared runner refuses today's declaration the rule parameters before a saved product is returned (KIT-129 B1)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 101L)
    todays <- resolveDIARuvCovariatesParameters(list(), "group")
    declared <- resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = 1)
    manager <- limpaRuvManager(case$quantified)
    withDIAStateScope(manager, function() {
        source <- getDIAStateReference(manager)
        saved <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation(), source, todays))
        expect_identical(saved$status, "computed")
        # the auditor's case through the runner itself (the replay's path), with nothing before it: each rule parameter
        # alone and the whole declared request are refused, never served today's saved sweep
        expect_error(runSavedOperation(manager, diaRuvCovariatesOperation(), source, declared), class = "workflow_execution_invalid")
        for (extra in list(list(ruv_control_rule = "all_features"), list(ruv_k_rule = "num_sv_be"), list(ruv_k = 2L),
            list(ruv_control_ids = "P1"), list(ruv_k_seed = 7L))) {
            expect_error(runSavedOperation(manager, diaRuvCovariatesOperation(), source, c(todays, extra)),
                class = "workflow_execution_invalid", info = names(extra))
        }
        # a plain request of today's design is still a hit; the declared rules compute their own product
        again <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation(), source, todays))
        expect_identical(again$status, "hit")
        expect_identical(again$reference$entity_id, saved$reference$entity_id)
        chosen <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation("declared"), source, declared))
        expect_identical(chosen$status, "computed")
        expect_false(identical(chosen$reference$entity_id, saved$reference$entity_id))
        expect_identical(chosen$value$ruv_control_rule, "all_features")
    })
    expect_null(workflowLegacyDesignRefusal(diaRuvCovariatesOperation("declared")))
    expect_identical(workflowLegacyDesignRefusal(diaRuvCovariatesOperation()), assertDiaRuvSweepParameters)
})

test_that("with technical replicates a declared state needs a residual degree of freedom between the units at the DA (M28)", {
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 105L, per_group = 8L, units = 4L)
    state <- limpaRuvKernel(case$quantified, resolveDIARuvCovariatesParameters(list(), "group", ruv_control_rule = "all_features",
        ruv_k_rule = "user", ruv_k = 5))$state
    expect_gt(nrow(limpaRuvDa(state)$contrasts_results$results[[1L]]), 0L)
    # a covariate between the units (one per unit) leaves the 8 units no degree of freedom beside 2 groups and 5 factors
    state@design_matrix$site <- c("s1", "s2")[(match(state@design_matrix$tech_rep_group, unique(state@design_matrix$tech_rep_group)) %% 2L) + 1L]
    refused <- tryCatch(limpaRuvDaWithFormula(state, "~ 0 + group + site"), error = function(e) e)
    expect_s3_class(refused, "omics_ruv_no_residual_df")
    expect_match(conditionMessage(refused), "between the biological units", fixed = TRUE)
})

# The legacy DA of this route on fixed inputs (decision M29: legacy unchanged). The digests were recorded on the tree
# before M29 was implemented (`kit-129-integration` d4b5db27), with the limpa engine: the quantified state without RUV,
# and today's covariates state (W in the design, no declared rule) with one factor.
limpaRuvLegacyDigests <- function() {
    case <- limpaRuvCase("limpa", seed = 101L, proteins = 120L)
    factors <- limpaRuvQuiet(computeRuvCovariates(case$quantified, 1L, rep(TRUE, length(case$ids)), "group"))
    todays <- case$quantified
    todays@design_matrix$W1 <- unname(factors$W[todays@design_matrix$Run, 1L])
    todays@args$ruvCovariates <- list(k = 1L, ruv_mode = "covariates")
    vapply(list(no_ruv = case$quantified, todays_covariates = todays), function(state) {
        result <- limpaRuvDa(state)
        digest::digest(list(result$contrasts_results$results, result$contrasts_results$fit.eb, result$da_model))
    }, character(1L))
}

test_that("the legacy DA (no RUV, today's covariates) is unchanged by the exact standard errors (M29)", {
    skip_if_not_installed("limpa")
    expect_identical(limpaRuvLegacyDigests(), c(no_ruv = "b1316e237d8b23073dba131c0daba919",
        todays_covariates = "a610463e3967a448bafeb71dfbf2e068"))
})
