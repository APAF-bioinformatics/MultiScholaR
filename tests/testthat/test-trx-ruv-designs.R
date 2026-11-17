# KIT-129 W4: the transcriptomics RUV designs beside today's grid, behind a declared parameter nothing selects yet.
# Each design is compared with the lines of KIT-124's simulator S3 that ran its arm
# (docs/qa/workflow-execution/kit/kit-124/sim/s3_trx.R), written out here with direct EDASeq, RUVSeq and edgeR calls:
#   arm 7  RUVs on upper-quartile counts, every gene a control, the groups as replicate sets, W in the edgeR model;
#   arm 8  RUVg on upper-quartile counts with a supplied control list, W in the model;
#   arm 5  RUVg with the genes outside the most significant quarter of the package's first pass, W in the model;
#   arm 1  no RUV, ~ 0 + group.
# k is an argument of every design (the shared rule for k is KIT-129 W6's). The data are the simulator's centre cell at
# a smaller size, generated here from a fixed seed.

trxDesignTestPackages <- function() {
    for (package in c("edgeR", "limma", "EDASeq", "RUVSeq")) testthat::skip_if_not_installed(package)
}

# The simulator's data model (s3_trx.R generate()): log-normal baseline counts, the voom dispersion trend, effects of
# 1 log2 in group B for a share of the genes, one unwanted factor orthogonal to the groups with strength ratio 3.
trxDesignTestData <- function(seed = 129L, genes = 1500L, per_group = 6L, groups = 2L, share = 0.2) {
    withr::local_seed(seed)
    levels <- LETTERS[seq_len(groups)]
    group <- factor(rep(levels, each = per_group), levels = levels)
    n <- length(group)
    log2_lambda <- stats::rnorm(genes, 6.5, 2)
    lambda <- 2^log2_lambda
    phi <- (0.2 + 1 / sqrt(lambda))^2 * 40 / stats::rchisq(genes, 40)
    n_effects <- round(share * genes)
    beta <- matrix(0, genes, groups, dimnames = list(NULL, levels))
    beta[seq_len(n_effects), "B"] <- rep(c(1, -1), length.out = n_effects)
    depth <- stats::rlnorm(n, -0.25^2 / 2, 0.25)
    z <- qr.resid(qr(stats::model.matrix(~ group)), matrix(stats::rnorm(n), n, 1L))
    w <- qr.Q(qr(z)) * sqrt(n - 1)
    alpha <- stats::rnorm(genes) * sqrt(3 * stats::median((phi + 1 / lambda) / log(2)^2))
    perm <- sample.int(genes)
    lambda <- lambda[perm]
    phi <- phi[perm]
    beta <- beta[perm, , drop = FALSE]
    alpha <- alpha[perm]
    mu <- outer(lambda, depth) * 2^(beta[, as.integer(group), drop = FALSE] + outer(alpha, w[, 1L]))
    counts <- matrix(as.double(stats::rnbinom(genes * n, mu = mu, size = 1 / phi)), genes, n)
    runs <- sprintf("%s_%02d", group, stats::ave(seq_len(n), group, FUN = seq_along))
    dimnames(counts) <- list(sprintf("g%05d", seq_len(genes)), runs)
    design <- data.frame(Run = runs, group = group, stringsAsFactors = FALSE)
    f <- trxExpressionFilter(counts, design$group)
    y <- f$filtered_counts
    list(y = y, design = design, log_cpm_tmm = trxTmmLogCpm(y, design$group)$log_cpm_tmm,
        null = rownames(y)[rowSums(beta[match(rownames(y), rownames(counts)), , drop = FALSE] != 0) == 0])
}

# The simulator's upper-quartile counts and first-pass p-values (s3_trx.R, the shared preparation of every arm).
trxDesignTestUq <- function(y) EDASeq::betweenLaneNormalization(y, which = "upper", round = TRUE)

trxDesignTestFirstPass <- function(y, design) {
    model <- stats::model.matrix(~ group, data = design)
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = y, group = design$group), method = "TMM")
    first <- edgeR::glmQLFit(edgeR::estimateDisp(dge, model), model)
    edgeR::topTags(edgeR::glmQLFTest(first, coef = 2:ncol(model)), n = nrow(y), sort.by = "none")$table$PValue
}

# The package's test by direct edgeR calls (TMM, the robust dispersion and QL fit, the QL F-test of B - A), with an arm's
# factors as covariates; the table in the order of the counts' rows.
trxDesignTestEdgeR <- function(y, design, W = NULL) {
    pdata <- data.frame(design, row.names = design$Run)
    formula <- "~ 0 + group"
    if (!is.null(W)) {
        W <- as.matrix(W)
        colnames(W) <- paste0("W_", seq_len(ncol(W)))
        pdata <- cbind(pdata, W)
        formula <- paste("~ 0 + group +", paste(colnames(W), collapse = " + "))
    }
    model <- stats::model.matrix(stats::as.formula(formula), data = pdata)
    colnames(model) <- make.names(colnames(model))
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = y, group = pdata$group), method = "TMM")
    fit <- edgeR::glmQLFit(edgeR::estimateDisp(dge, model, robust = TRUE), model, robust = TRUE)
    contrast <- limma::makeContrasts(contrasts = "groupB-groupA", levels = model)
    table <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast = contrast[, 1L]), n = Inf, sort.by = "PValue")$table
    table[rownames(y), c("logFC", "PValue", "FDR")]
}

# The package's entry: the design's choice (controls and k), the factors, then the session's fit and test halves on the
# factors' sample table with the formula realised from it.
trxDesignTestPackage <- function(d, ruv_design, k = NA_integer_, k_rule = NA_character_, control_ids = character()) {
    x <- list(counts = d$y, design = d$design)
    choice <- trxChooseRuvDesign(x, trxRuvDesignParameters(ruv_design, k, k_rule, control_ids))
    s <- trxRuvFactors(d$y, d$log_cpm_tmm, d$design, choice$chosen$control_genes, choice$chosen$k, choice$ruv_design$method)
    table <- trxDaApply(trxDaFit(s$ruvs_counts, s$pdata, trxDaFormula(s$pdata)), "B_vs_A=groupB-groupA")$results[[1L]]$top_table
    list(choice = choice, s = s, table = table[rownames(d$y), c("logFC", "PValue", "FDR")])
}

test_that("(KIT-129 W4) the recommended design is the simulator's arm 7: RUVs on upper-quartile counts, all genes, groups as replicate sets, W in the model", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData()))
    genes <- rownames(d$y)
    uq <- trxDesignTestUq(d$y)
    expect_identical(trxUpperQuartileCounts(d$y), uq)
    for (k in 1:2) {
        arm <- RUVSeq::RUVs(uq, genes, k, RUVSeq::makeGroups(d$design$group))
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "groups_all_features", k, "user")))
        expect_identical(p$choice$chosen$control_genes, genes)
        expect_identical(p$choice$chosen$k, as.integer(k))
        expect_identical(p$choice$ruv_design[c("id", "method", "control_rule", "k_rule", "applied", "test_form")],
            list(id = "groups_all_features", method = "RUVs", control_rule = "all_features", k_rule = "user", applied = TRUE,
                test_form = "factors_in_model"))
        expect_identical(unname(as.matrix(p$s$W)), unname(arm$W))
        expect_identical(names(p$s$W), paste0("W_", seq_len(k)))
        expect_identical(rownames(p$s$W), d$design$Run)
        # the test reads the filtered counts, not the upper-quartile or the corrected ones
        expect_identical(p$s$ruvs_counts, d$y)
        expect_identical(names(p$s$pdata), c("Run", "group", paste0("W_", seq_len(k))))
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))))
        # the display matrix, for figures and export: today's removeBatchEffect call with these factors
        expect_identical(p$s$log_cpm_ruvs, limma::removeBatchEffect(d$log_cpm_tmm, covariates = arm$W,
            design = stats::model.matrix(~ group, data = d$design)))
    }
})

test_that("(KIT-129 W4) the recommended design with the estimated count is arm 7 with E", {
    trxDesignTestPackages()
    testthat::skip_if_not_installed("sva")
    d <- suppressWarnings(suppressMessages(trxDesignTestData()))
    E <- as.integer(sva::num.sv(d$log_cpm_tmm, stats::model.matrix(~ group, data = d$design), method = "be", B = 20, seed = 129L))
    expect_gte(E, 1L)
    arm <- RUVSeq::RUVs(trxDesignTestUq(d$y), rownames(d$y), E, RUVSeq::makeGroups(d$design$group))
    p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "groups_all_features", E, "num_sv_be")))
    expect_identical(p$choice$ruv_design$k_rule, "num_sv_be")
    expect_identical(unname(as.matrix(p$s$W)), unname(arm$W))
    expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))))
})

test_that("(KIT-129 W4) the supplied-control design is the simulator's arm 8: RUVg on upper-quartile counts with the list, W in the model", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData()))
    uq <- trxDesignTestUq(d$y)
    supplied <- withr::with_seed(8L, sample(d$null, 50L))
    absent <- c("not_a_gene_1", "not_a_gene_2")
    matched <- rownames(d$y)[rownames(d$y) %in% supplied]
    for (k in 1:2) {
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "supplied_controls", k, "user", c(supplied, absent, supplied[1:3]))))
        # the controls are the listed genes found in the matrix, in its row order; the ids not found are reported
        expect_identical(p$choice$chosen$control_genes, matched)
        expect_identical(p$choice$chosen$n_controls, 50L)
        expect_identical(p$choice$ruv_design$unmatched, absent)
        expect_identical(p$choice$ruv_design[c("method", "control_rule")], list(method = "RUVg", control_rule = "supplied_list"))
        arm <- RUVSeq::RUVg(uq, matched, k)
        expect_identical(unname(as.matrix(p$s$W)), unname(arm$W))
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))))
        # the simulator passes the list in the order it was drawn: the same factors and table to rounding
        drawn <- RUVSeq::RUVg(uq, supplied, k)
        expect_equal(unname(as.matrix(p$s$W)), unname(drawn$W), tolerance = 1e-10)
        expect_equal(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, drawn$W))), tolerance = 1e-8)
    }
    # the same list in another order is the same declared parameter
    expect_identical(trxRuvDesignParameters("supplied_controls", 1L, "user", rev(supplied)),
        trxRuvDesignParameters("supplied_controls", 1L, "user", c(supplied, supplied[1:2])))
})

test_that("(KIT-129 W4) the first-pass design is the simulator's arm 5: RUVg with the genes outside the most significant quarter, k given", {
    trxDesignTestPackages()
    testthat::skip_if_not_installed("ruv")
    d <- suppressWarnings(suppressMessages(trxDesignTestData()))
    genes <- rownames(d$y)
    uq <- trxDesignTestUq(d$y)
    first_p <- suppressWarnings(suppressMessages(trxDesignTestFirstPass(d$y, d$design)))
    controls <- genes[order(first_p)[-seq_len(ceiling(0.25 * length(genes)))]]
    # the arm's k is today's grid's
    grid <- suppressWarnings(suppressMessages(trxRuvGrid(d$y, d$log_cpm_tmm, d$design)))
    for (k in unique(c(grid$best$best_k, 1L))) {
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "first_pass_nulls", k, "sweep_bounded")))
        expect_identical(p$choice$chosen$control_genes, controls)
        expect_identical(p$choice$ruv_design[c("method", "control_rule", "k_rule")],
            list(method = "RUVg", control_rule = "first_pass_drop_top_25", k_rule = "sweep_bounded"))
        arm <- RUVSeq::RUVg(uq, controls, k)
        expect_identical(unname(as.matrix(p$s$W)), unname(arm$W))
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))))
    }
    # the first pass is the package's own, over every group coefficient; today's ranking reads the same p-values
    three <- suppressWarnings(suppressMessages(trxDesignTestData(seed = 7L, genes = 600L, groups = 3L)))
    model <- stats::model.matrix(~ group, data = three$design)
    p3 <- suppressWarnings(suppressMessages(trxPrelimPValues(trxFilteredDge(three$y, three$design$group), model)))
    expect_identical(p3, suppressWarnings(suppressMessages(trxDesignTestFirstPass(three$y, three$design))))
    expect_identical(suppressWarnings(suppressMessages(trxPrelimOrder(trxFilteredDge(three$y, three$design$group), model))),
        order(p3, decreasing = TRUE))
})

test_that("(KIT-129 W4) the none design and a k of 0 are the simulator's arm 1: no factor, ~ 0 + group, the display matrix uncorrected", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData()))
    arm <- suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design)))
    none <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "none")))
    zero <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "groups_all_features", 0L, "num_sv_be")))
    for (p in list(none, zero)) {
        expect_identical(p$choice$chosen$k, 0L)
        expect_false(p$choice$ruv_design$applied)
        expect_identical(p$choice$ruv_design$test_form, "none")
        expect_identical(ncol(p$s$W), 0L)
        expect_identical(names(p$s$pdata), c("Run", "group"))
        expect_identical(trxDaFormula(p$s$pdata), "~ 0 + group")
        expect_identical(p$s$log_cpm_ruvs, d$log_cpm_tmm)
        expect_identical(p$s$ruvs_counts, d$y)
        expect_identical(p$table, arm)
    }
    expect_identical(none$choice$ruv_design[c("id", "method", "control_rule")], list(id = "none", method = NA_character_, control_rule = NA_character_))
    expect_match(none$choice$ruv_design$reason, "none")
    expect_match(zero$choice$ruv_design$reason, "no unwanted factor")
    # the outcomes the RUVs step records say so
    expect_identical(trxRuvDesignOutcomes(zero$choice)[c("ruv_design", "ruv_k", "ruv_k_rule", "ruv_test_form", "ruv_applied", "k")],
        list(ruv_design = "groups_all_features", ruv_k = 0L, ruv_k_rule = "num_sv_be", ruv_test_form = "none", ruv_applied = FALSE, k = 0L))
    # a design with factors still realises today's formula
    expect_identical(trxDaFormula(data.frame(Run = "a", group = "A", W_1 = 0, W_2 = 0)), "~ 0 + group + W_1 + W_2")
})

test_that("(KIT-129 W4) what a design cannot compute is refused by class, never reduced or replaced", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    x <- list(counts = d$y, design = d$design)
    genes <- rownames(d$y)
    choose <- function(...) trxChooseRuvDesign(x, trxRuvDesignParameters(...))
    expect_error(choose("legacy", 1L, "user"), class = "trx_ruv_design_unsupported")
    expect_error(choose("technical_replicates", 1L, "user"), class = "trx_ruv_design_unsupported")
    expect_error(choose("groups_all_features"), class = "trx_ruv_design_invalid")
    expect_error(choose("groups_all_features", 1L, "elbow"), class = "trx_ruv_design_invalid")
    expect_error(choose("groups_all_features", 0L, "user"), class = "trx_ruv_design_invalid")
    expect_error(choose("groups_all_features", 1.5, "user"), class = "trx_ruv_design_invalid")
    expect_error(choose("groups_all_features", 1L, "user", genes[1:20]), class = "trx_ruv_design_invalid")
    expect_error(choose("none", 1L, "user"), class = "trx_ruv_design_invalid")
    expect_error(choose("supplied_controls", 1L, "user"), class = "trx_ruv_design_invalid")
    expect_error(choose("supplied_controls", 1L, "user", c("not_a_gene_1", "not_a_gene_2")), class = "trx_ruv_controls_invalid")
    # RUVSeq returns fewer factors than asked without saying so; the kernel refuses instead
    samples <- ncol(d$y)
    expect_lt(ncol(RUVSeq::RUVs(trxDesignTestUq(d$y), genes, samples - 1L, RUVSeq::makeGroups(d$design$group))$W), samples - 1L)
    factors <- function(controls, k, method) suppressWarnings(suppressMessages(trxRuvFactors(d$y, d$log_cpm_tmm, d$design, controls, k, method)))
    # 12 samples in 2 groups: at most 9 factors leave a residual degree of freedom beside the groups
    expect_error(factors(genes, samples - 2L, "RUVs"), class = "trx_ruv_k_unavailable")
    expect_identical(ncol(factors(genes, samples - 3L, "RUVs")$W), samples - 3L)
    # fewer factors than controls
    expect_error(factors(genes[1:3], 3L, "RUVg"), class = "trx_ruv_k_unavailable")
    expect_error(factors(c(genes[1:20], "not_a_gene"), 1L, "RUVg"), class = "trx_ruv_controls_invalid")
    expect_error(factors(genes, 1L, "RUVr"))
    for (class in c("trx_ruv_design_unsupported", "trx_ruv_design_invalid", "trx_ruv_controls_invalid", "trx_ruv_k_unavailable")) {
        condition <- tryCatch(switch(class,
            trx_ruv_design_unsupported = choose("legacy", 1L, "user"),
            trx_ruv_design_invalid = choose("none", 1L, "user"),
            trx_ruv_controls_invalid = choose("supplied_controls", 1L, "user", "not_a_gene"),
            trx_ruv_k_unavailable = factors(genes[1:3], 3L, "RUVg")), error = function(e) e)
        expect_s3_class(condition, c(class, "trx_error"), exact = FALSE)
    }
})

test_that("(KIT-129 W4) today's declarations are the default: the grid's and the RUVs step's names and parameters are unchanged", {
    binding <- transcriptomicsDeclarationBinding
    todays_grid <- c("logcpm", "ruv_mode", "percentage_min", "percentage_max", "percentage_step", "k_penalty_weight", "manual_k",
        "manual_percentage")
    for (declaration in list(trxRuvGridOperation(binding), trxRuvGridOperation(binding, "legacy"), workflowOperationEntry("trx_ruv_grid")$declaration)) {
        expect_identical(declaration$name, "trx_ruv_grid")
        expect_identical(declaration$parameters, todays_grid)
        expect_identical(declaration$product_type, "trx.ruv_grid")
        expect_identical(names(declaration$fixed), c("max_k", "controls"))
    }
    declared <- trxRuvGridOperation(binding, "declared")
    expect_identical(declared$name, "trx_ruv_grid")
    expect_identical(declared$parameters, c("logcpm", "ruv_design", "ruv_k", "ruv_k_rule", "ruv_control_ids"))
    expect_identical(declared$product_type, "trx.ruv_grid")
    for (declaration in list(trxRuvsOperation(binding), workflowOperationEntry("trx_ruvs")$declaration)) {
        expect_identical(declaration$parameters, c("logcpm", "grid"))
    }
    expect_identical(TRX_RUV_DESIGNS, c("groups_all_features", "supplied_controls", "first_pass_nulls", "none"))
    expect_error(trxRuvGridOperation(binding, "groups_all_features"))
})

# The session's steps on a kit corpus (the helper of test-workflow-trx-declarations.R): the declared design's choice, the
# RUVs state and the DA, in process on a session owner.
trxDesignTestOwner <- function(env = parent.frame()) {
    skipWithoutTrxDeclReference()
    reference <- trxDeclReference("trx_edge_12", env = env)
    fixture <- trxDeclarationsOwner(withr::local_tempdir(.local_envir = env), reference$dir)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    suppressWarnings(suppressMessages(runTrxReferenceExpressionFilter(fixture$owner, list())))
    commitWorkflowOmicsStoreStep(fixture$owner, "trx_qc_finalize", "trx_qc_complete", list(), "QC processing complete - ready for normalization")
    fixture
}

test_that("(KIT-129 W4) through the declarations: the recommended design commits W to the RUVs state and the DA tests the counts with W in the model", {
    fixture <- trxDesignTestOwner()
    owner <- fixture$owner
    read <- function(name) trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, name)))
    filtered <- read("trx_qc_complete")
    # today's grid on the same source, run first: another product, with no declared design
    legacy <- suppressWarnings(suppressMessages(withWorkflowStateScope(owner, function() {
        source <- getWorkflowStateReference(owner)
        logcpm <- trxRunProduct(owner, "trx_tmm_logcpm", source, list())
        trxRunProduct(owner, "trx_ruv_grid", source, c(list(logcpm = logcpm$reference$entity_id), trxRuvGridParameters()), list(logcpm$reference))
    })))
    expect_null(legacy$value[["ruv_design", exact = TRUE]])
    expect_identical(legacy$value$mode, "automatic")
    run <- suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, NULL, list(),
        design_parameters = trxRuvDesignParameters("groups_all_features", 2L, "user"))))
    expect_identical(owner$current_state, "trx_ruvs")
    expect_false(identical(legacy$reference$entity_id, run$grid$reference$entity_id))
    choice <- run$grid$value
    expect_identical(choice$ruv_design$id, "groups_all_features")
    expect_identical(choice$chosen$control_genes, rownames(filtered$counts))
    arm <- RUVSeq::RUVs(trxDesignTestUq(filtered$counts), rownames(filtered$counts), 2L, RUVSeq::makeGroups(filtered$design$group))
    expect_identical(unname(as.matrix(run$value$W)), unname(arm$W))
    expect_identical(run$value$k, 2L)
    expect_identical(run$value$outcomes[c("ruv_design", "ruv_method", "ruv_control_rule", "ruv_control_count", "ruv_k", "ruv_k_rule",
        "ruv_test_form", "ruv_applied")], list(ruv_design = "groups_all_features", ruv_method = "RUVs", ruv_control_rule = "all_features",
        ruv_control_count = nrow(filtered$counts), ruv_k = 2L, ruv_k_rule = "user", ruv_test_form = "factors_in_model", ruv_applied = TRUE))
    state <- read("trx_ruvs")
    expect_identical(state$counts, filtered$counts)
    expect_identical(unname(as.matrix(state$design[, c("W_1", "W_2")])), unname(arm$W))
    # the display matrix is saved for the figures and the export
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$logcpm$value$log_cpm))
    corrected <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$value$corrected))
    expect_equal(corrected, limma::removeBatchEffect(log_cpm, covariates = arm$W, design = stats::model.matrix(~ group, data = filtered$design)),
        tolerance = 1e-12)
    da <- suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1), fixture$workflow$contrasts_tbl)))
    expect_identical(da$manifest$parameters$formula, "~ 0 + group + W_1 + W_2")
    expect_identical(da$manifest$model$w_columns, c("W_1", "W_2"))
    top <- trxDeclValue(owner, da$manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
    pdata <- cbind(data.frame(filtered$design, row.names = filtered$design$Run), arm$W)
    entry <- names(da$manifest$results)[[1L]]
    direct <- suppressWarnings(suppressMessages(trxDaApply(trxDaFit(filtered$counts, pdata, "~ 0 + group + W_1 + W_2"), entry)))$results[[1L]]$top_table
    expect_identical(top, direct)
})

test_that("(KIT-129 W4) through the declarations: the none design commits a RUVs state without factors and the DA fits ~ 0 + group", {
    fixture <- trxDesignTestOwner()
    owner <- fixture$owner
    read <- function(name) trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, name)))
    filtered <- read("trx_qc_complete")
    run <- suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters("none"))))
    expect_identical(owner$current_state, "trx_ruvs")
    expect_identical(run$value$k, 0L)
    expect_identical(ncol(run$value$W), 0L)
    expect_identical(run$value$outcomes[c("ruv_design", "ruv_test_form", "ruv_applied")],
        list(ruv_design = "none", ruv_test_form = "none", ruv_applied = FALSE))
    state <- read("trx_ruvs")
    expect_identical(state$counts, filtered$counts)
    expect_identical(names(state$design), names(filtered$design))
    log_cpm <- withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$logcpm$value$log_cpm))
    expect_identical(withWorkflowStateScope(owner, function() trxSavedMatrix(owner, run$value$corrected)), log_cpm)
    # the tab's final PCA reads the state alone: without factors its matrix is the TMM logCPM
    expect_equal(trxStateCorrectedLogCpm(state), log_cpm, tolerance = 1e-12)
    da <- suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1), fixture$workflow$contrasts_tbl)))
    expect_identical(da$manifest$parameters$formula, "~ 0 + group")
    expect_identical(da$manifest$model$w_columns, character())
    top <- trxDeclValue(owner, da$manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
    entry <- names(da$manifest$results)[[1L]]
    pdata <- data.frame(filtered$design, row.names = filtered$design$Run)
    direct <- suppressWarnings(suppressMessages(trxDaApply(trxDaFit(filtered$counts, pdata, "~ 0 + group"), entry)))$results[[1L]]$top_table
    expect_identical(top, direct)
})
