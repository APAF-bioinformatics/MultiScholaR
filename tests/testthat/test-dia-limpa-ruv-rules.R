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

# The route's test as the simulator calls it (`routeTest` of s2_limpa.R): `covariates` join the design and the formula.
limpaRuvRouteTest <- function(elist, design, covariates = NULL) {
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
    fit <- limpaRuvQuiet(limpa::dpcDE(elist[, samples], model, plot = FALSE, block = block))
    result <- limpaRuvQuiet(convertDpcDAToStandardFormat(dpc_fit = fit, contrast_strings = limpaRuvContrast,
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

# E as the simulator computes it.
limpaRuvSimulatorE <- function(case, seed) {
    as.integer(sva::num.sv(case$elist$E, stats::model.matrix(~ group, case$design), method = "be", B = 20, seed = seed))
}

# The simulator's `inModel`: the factors of `computeRuvCovariates`, then the route's test with them.
limpaRuvInModel <- function(case, k, controls) {
    factors <- limpaRuvQuiet(computeRuvCovariates(case$quantified, k, controls, "group"))
    list(W = factors$W, replicate_structure = factors$replicate_structure,
        table = limpaRuvRouteTest(case$elist, case$design, factors$W))
}

# The package: the declared kernel on the quantified object (the store keeps what it checkpoints), and the DA entry
# point on the state the kernel returns as final.
limpaRuvKernel <- function(quantified, parameters, ruv_rules = "declared") {
    saved <- list()
    store <- list(checkpoint = function(value, label, from = NULL) {
        saved[[label]] <<- value
        label
    })
    value <- diaRuvCovariatesOperation(ruv_rules)$kernel(quantified, parameters, store)
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
    expect_equal(table$logFC, reference$logFC, tolerance = 1e-12, info = info)
    expect_equal(table$raw_pvalue, reference$p, tolerance = 1e-12, info = info)
    expect_equal(table$fdr_qvalue, reference$q, tolerance = 1e-12, info = info)
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

test_that("declared technical replicates stay RUV-III's replicate sets on this route (arm 16)", {
    skip_if_not_installed("sva")
    skip_if_not_installed("limpa")
    engine <- utils::tail(limpaRuvEngines(), 1L)
    case <- limpaRuvCase(engine, seed = 105L, per_group = 8L, units = 4L)
    E <- limpaRuvSimulatorE(case, 4242L)
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
    # 16 samples in 8 replicate sets carry at most 8 factors
    expect_error(computeDiaRuvCovariatesByRule(case$quantified, resolveDIARuvCovariatesParameters(list(), "group",
        ruv_control_rule = "all_features", ruv_k_rule = "user", ruv_k = 9)), class = "ruv_k_above_limit")
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
    package <- limpaRuvKernel(case$quantified, parameters)
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
    # the sweep's own record is kept for this design, as today
    expect_identical(value$ruv_result$optimization_results, sweep$optimization_results)
    expect_identical(value$ruv_result$max_acceptable_k, sweep$max_acceptable_k)
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
        expect_identical(diaRuvCovariatesCommit(manager, parameters, run)$state_name, "ruv_covariates")
        # today's declaration on the same source: its own identity, its ten parameters, the sweep's controls
        legacy <- limpaRuvQuiet(runSavedOperation(manager, diaRuvCovariatesOperation(), source, resolveDIARuvCovariatesParameters(list(), "group")))
        expect_identical(legacy$status, "computed")
        expect_false(identical(legacy$reference$entity_id, run$reference$entity_id))
        expect_identical(names(legacy$value$parameters), limpaRuvSweepParameters)
        expect_lt(legacy$value$control_count, length(case$ids))
        expect_null(legacy$value$ruv_control_rule)
    })
})
