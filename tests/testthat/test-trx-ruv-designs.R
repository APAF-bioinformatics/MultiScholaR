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
    group <- factor(rep(levels, times = rep_len(per_group, groups)), levels = levels)
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

# The classed warnings an expression raises, with its value (audit F-06).
trxDesignTestWarnings <- function(expr) {
    caught <- list()
    value <- withCallingHandlers(expr, warning = function(w) {
        caught[[length(caught) + 1L]] <<- w
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = caught)
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
        # audit F-06: the list's hash is of the sorted ids, each once, whatever was matched
        expect_identical(p$choice$ruv_design$list_sha256, digest::digest(paste(sort(c(supplied, absent), method = "radix"), collapse = "\n"),
            algo = "sha256", serialize = FALSE))
        expect_identical(trxRuvDesignOutcomes(p$choice)$ruv_control_list_sha256, p$choice$ruv_design$list_sha256)
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
    # the arm passes its controls in the first pass's ranking
    ranked <- genes[order(first_p)[-seq_len(ceiling(0.25 * length(genes)))]]
    # audit F-06: the package passes the same genes in the matrix's row order, as it passes a supplied list
    controls <- genes[genes %in% ranked]
    expect_setequal(controls, ranked)
    expect_false(identical(controls, ranked))
    # the arm's k is today's grid's
    grid <- suppressWarnings(suppressMessages(trxRuvGrid(d$y, d$log_cpm_tmm, d$design)))
    for (k in unique(c(grid$best$best_k, 1L))) {
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "first_pass_nulls", k, "sweep_bounded")))
        expect_identical(p$choice$chosen$control_genes, controls)
        expect_identical(p$choice$ruv_design[c("method", "control_rule", "k_rule")],
            list(method = "RUVg", control_rule = "first_pass_drop_top_25", k_rule = "sweep_bounded"))
        row_order <- RUVSeq::RUVg(uq, controls, k)
        expect_identical(unname(as.matrix(p$s$W)), unname(row_order$W))
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, row_order$W))))
        # the arm's own call, with the ranking's order: RUVg does not depend on the order of its controls beyond
        # rounding, so the arm is reproduced at the tolerances of the supplied list (W 1e-10, the table 1e-8)
        arm <- RUVSeq::RUVg(uq, ranked, k)
        expect_equal(unname(as.matrix(p$s$W)), unname(arm$W), tolerance = 1e-10)
        expect_equal(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))), tolerance = 1e-8)
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
    # audit F-09: the "returned fewer factors" refusal is reached (above it is the residual rule that refuses): five
    # control genes of rank 1 carry one factor, and RUVg returns one where two were asked
    deficient <- d$y
    deficient[1:5, ] <- matrix(rep(deficient[1L, ], each = 5L), 5L)
    expect_identical(ncol(RUVSeq::RUVg(trxDesignTestUq(deficient), genes[1:5], 2L)$W), 1L)
    fewer <- tryCatch(suppressWarnings(suppressMessages(trxRuvFactors(deficient, d$log_cpm_tmm, d$design, genes[1:5], 2L, "RUVg"))), error = function(e) e)
    expect_s3_class(fewer, c("trx_ruv_k_unavailable", "trx_error"), exact = FALSE)
    expect_match(conditionMessage(fewer), "RUVSeq::RUVg returned 1 factor", fixed = TRUE)
    # audit F-09: a run whose upper quartile is 0 cannot be upper-quartile normalised
    sparse <- d$y
    sparse[-(1:50), 3L] <- 0
    expect_error(trxUpperQuartileCounts(sparse), class = "trx_ruv_upper_quartile_invalid")
    expect_error(suppressWarnings(suppressMessages(trxRuvFactors(sparse, d$log_cpm_tmm, d$design, genes, 1L, "RUVs"))), class = "trx_ruv_upper_quartile_invalid")
    for (class in c("trx_ruv_design_unsupported", "trx_ruv_design_invalid", "trx_ruv_controls_invalid", "trx_ruv_k_unavailable")) {
        condition <- tryCatch(switch(class,
            trx_ruv_design_unsupported = choose("legacy", 1L, "user"),
            trx_ruv_design_invalid = choose("none", 1L, "user"),
            trx_ruv_controls_invalid = choose("supplied_controls", 1L, "user", "not_a_gene"),
            trx_ruv_k_unavailable = factors(genes[1:3], 3L, "RUVg")), error = function(e) e)
        expect_s3_class(condition, c(class, "trx_error"), exact = FALSE)
    }
})

test_that("(KIT-129 audit F-05) the factors' kernel refuses a design, counts and logCPM that are not one set of runs and genes", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    genes <- rownames(d$y)
    shuffled <- d$design[c(7:12, 1:6), ]
    rownames(shuffled) <- NULL
    refused <- function(expr) expect_error(suppressWarnings(suppressMessages(expr)), class = "trx_ruv_inputs_misaligned")
    for (method in c("RUVs", "RUVg")) refused(trxRuvFactors(d$y, d$log_cpm_tmm, shuffled, genes, 1L, method))
    # without factors too: the returned sample table would carry the wrong runs
    refused(trxRuvFactors(d$y, d$log_cpm_tmm, shuffled, character(), 0L))
    refused(trxRuvFactors(d$y, d$log_cpm_tmm[1:50, ], d$design, genes, 1L, "RUVs"))
    refused(trxRuvFactors(d$y, d$log_cpm_tmm[rev(genes), ], d$design, genes, 1L, "RUVs"))
    refused(trxRuvFactors(d$y, d$log_cpm_tmm[, c(7:12, 1:6)], d$design, genes, 1L, "RUVs"))
    refused(trxRuvFactors(d$y, d$log_cpm_tmm, d$design[-1L, ], genes, 1L, "RUVs"))
    # the controls' kernel reads the groups by position for its first pass
    for (ruv_design in c("groups_all_features", "first_pass_nulls")) refused(trxRuvDesignControls(d$y, shuffled, ruv_design))
    refused(trxChooseRuvDesign(list(counts = d$y, design = shuffled), trxRuvDesignParameters("first_pass_nulls", 1L, "sweep_bounded")))
    condition <- tryCatch(trxRuvFactors(d$y, d$log_cpm_tmm, shuffled, genes, 1L, "RUVs"), error = function(e) e)
    expect_s3_class(condition, c("trx_ruv_inputs_misaligned", "trx_error"), exact = FALSE)
    # aligned inputs pass as before
    expect_identical(ncol(suppressWarnings(suppressMessages(trxRuvFactors(d$y, d$log_cpm_tmm, d$design, genes, 1L, "RUVs")))$W), 1L)
})

test_that("(KIT-129 audit F-01) the factors' kernel drops the W columns a design already carries", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    genes <- rownames(d$y)
    carried <- cbind(d$design, W_1 = seq_len(nrow(d$design)), W_2 = rev(seq_len(nrow(d$design))))
    clean <- suppressWarnings(suppressMessages(trxRuvFactors(d$y, d$log_cpm_tmm, d$design, genes, 1L, "RUVs")))
    again <- suppressWarnings(suppressMessages(trxRuvFactors(d$y, d$log_cpm_tmm, carried, genes, 1L, "RUVs")))
    expect_identical(names(again$pdata), c("Run", "group", "W_1"))
    expect_identical(again, clean)
    expect_identical(trxDaFormula(again$pdata), "~ 0 + group + W_1")
    none <- trxRuvFactors(d$y, d$log_cpm_tmm, carried, character(), 0L)
    expect_identical(names(none$pdata), c("Run", "group"))
    expect_identical(trxDaFormula(none$pdata), "~ 0 + group")
    expect_identical(none, trxRuvFactors(d$y, d$log_cpm_tmm, d$design, character(), 0L))
})

test_that("(KIT-129 audit F-06) a supplied control list: each id once, NA and empty ids refused, unmatched ids warned of, fewer than two matches refused, the list's hash", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    genes <- rownames(d$y)
    ids <- genes[c(40L, 3L, 17L, 250L)]
    once <- trxRuvDesignControls(d$y, d$design, "supplied_controls", ids)
    expect_identical(once$control_genes, genes[c(3L, 17L, 40L, 250L)])
    # duplicates are removed in the kernel, without a warning
    twice <- trxDesignTestWarnings(trxRuvDesignControls(d$y, d$design, "supplied_controls", c(ids, rev(ids), ids[1:2])))
    expect_length(twice$warnings, 0L)
    expect_identical(twice$value, once)
    expect_identical(once$unmatched, character())
    # NA and empty ids are refused, with their number
    for (bad in list(c(ids, NA), c(ids, ""), c(ids, "  "), c(NA_character_, NA_character_))) {
        condition <- tryCatch(trxRuvDesignControls(d$y, d$design, "supplied_controls", bad), error = function(e) e)
        expect_s3_class(condition, c("trx_ruv_controls_invalid", "trx_error"), exact = FALSE)
        expect_match(conditionMessage(condition), "NA or empty", fixed = TRUE)
    }
    # the constructor does not lose an NA on the way to that refusal
    expect_true(anyNA(trxRuvDesignParameters("supplied_controls", 1L, "user", c(ids, NA))$ruv_control_ids))
    expect_error(trxChooseRuvDesign(list(counts = d$y, design = d$design), trxRuvDesignParameters("supplied_controls", 1L, "user", c(ids, NA))),
        class = "trx_ruv_controls_invalid")
    # every unmatched id is warned of: the count, and the ids
    absent <- c("not_a_gene_2", "not_a_gene_1")
    partly <- trxDesignTestWarnings(trxRuvDesignControls(d$y, d$design, "supplied_controls", c(ids, absent, absent)))
    expect_length(partly$warnings, 1L)
    expect_s3_class(partly$warnings[[1L]], c("trx_ruv_controls_unmatched", "trx_warning"), exact = FALSE)
    expect_match(conditionMessage(partly$warnings[[1L]]), "2 of the 6 supplied control ids", fixed = TRUE)
    expect_match(conditionMessage(partly$warnings[[1L]]), "not_a_gene_1, not_a_gene_2", fixed = TRUE)
    expect_identical(partly$value$control_genes, once$control_genes)
    expect_identical(partly$value$unmatched, sort(absent))
    # the auditor's case, one match of 500: refused, and so is no match; two matches are the fewest taken
    many <- sprintf("absent_%03d", 1:499)
    expect_error(suppressWarnings(trxRuvDesignControls(d$y, d$design, "supplied_controls", c(ids[1L], many))), class = "trx_ruv_controls_invalid")
    expect_error(suppressWarnings(trxRuvDesignControls(d$y, d$design, "supplied_controls", many)), class = "trx_ruv_controls_invalid")
    expect_identical(TRX_MIN_SUPPLIED_CONTROLS, 2L)
    two <- trxDesignTestWarnings(trxRuvDesignControls(d$y, d$design, "supplied_controls", c(ids[1:2], many)))
    expect_identical(two$value$control_genes, genes[c(3L, 40L)])
    expect_match(conditionMessage(two$warnings[[1L]]), "499 of the 501 supplied control ids", fixed = TRUE)
    # the hash: sha256 of the sorted ids, each once, joined by a line feed (a value computed outside R: printf 'a\nb' | sha256sum)
    expect_identical(trxControlListSha256(c("b", "a", "b")), "7e18f737311b2dc3b2f269dd78396b0351f14fb66efa879f768cb23181883c78")
    expect_identical(once$list_sha256, trxControlListSha256(rev(ids)))
    expect_identical(twice$value$list_sha256, once$list_sha256)
    expect_false(identical(partly$value$list_sha256, once$list_sha256))
    for (ruv_design in c("groups_all_features", "none")) expect_identical(trxRuvDesignControls(d$y, d$design, ruv_design)$list_sha256, NA_character_)
    recommended <- trxChooseRuvDesign(list(counts = d$y, design = d$design), trxRuvDesignParameters("groups_all_features", 1L, "user"))
    expect_identical(trxRuvDesignOutcomes(recommended)$ruv_control_list_sha256, NA_character_)
})

test_that("(KIT-129 audit F-07) the choice does not depend on how its parameters were passed; a k that is not whole is refused", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    x <- list(counts = d$y, design = d$design)
    ids <- rownames(d$y)[c(40L, 3L, 17L, 250L)]
    expect_identical(trxRuvDesignParameters("groups_all_features", 2, "user"), trxRuvDesignParameters("groups_all_features", 2L, "user"))
    expect_identical(trxCanonicalRuvDesignParameters(list(ruv_design = "groups_all_features", ruv_k = 2, ruv_k_rule = "user")),
        list(ruv_design = "groups_all_features", ruv_k = 2L, ruv_k_rule = "user", ruv_control_ids = character()))
    # the design none takes no k: none given, NA and 0 are one declaration
    for (k in list(NULL, NA, NA_real_, 0, 0L)) {
        expect_identical(trxCanonicalRuvDesignParameters(list(ruv_design = "none", ruv_k = k)), trxRuvDesignParameters("none"))
    }
    # a plain list, as a request carries it, gives the constructor's choice
    built <- trxChooseRuvDesign(x, trxRuvDesignParameters("supplied_controls", 1L, "user", ids))
    plain <- trxChooseRuvDesign(x, list(ruv_design = "supplied_controls", ruv_k = 1, ruv_k_rule = "user", ruv_control_ids = c(rev(ids), ids[1:2])))
    expect_identical(plain, built)
    expect_identical(built$settings$ruv_control_ids, sort(ids, method = "radix"))
    expect_identical(built$settings$ruv_k, 1L)
    expect_identical(built$k_record$k, 1L)
    # the canonical order is the C locale's, the same on every machine
    expect_identical(trxCanonicalControlIds(c("b", "A", "a", "B", "a")), c("A", "B", "a", "b"))
    # a k that is not one whole number of 0 or more is refused, never rounded
    for (k in list(2.5, -1L, "2", c(1L, 2L), TRUE)) {
        expect_error(trxCanonicalRuvDesignParameters(list(ruv_design = "groups_all_features", ruv_k = k, ruv_k_rule = "user")), class = "trx_ruv_design_invalid")
        expect_error(trxChooseRuvDesign(x, list(ruv_design = "groups_all_features", ruv_k = k, ruv_k_rule = "user", ruv_control_ids = character())),
            class = "trx_ruv_design_invalid")
    }
})

test_that("(KIT-129 audit F-08) degenerate inputs are refused by class: a k of the wrong type at the kernel, identical runs within every group under RUVs", {
    trxDesignTestPackages()
    d <- suppressWarnings(suppressMessages(trxDesignTestData(genes = 600L)))
    genes <- rownames(d$y)
    for (k in list(2.5, NA_integer_, NA_real_, "2", -1L, c(1L, 2L), NULL, TRUE, 1e10)) {
        condition <- tryCatch(trxRuvFactors(d$y, d$log_cpm_tmm, d$design, genes, k, "RUVs"), error = function(e) e)
        expect_s3_class(condition, c("trx_ruv_k_invalid", "trx_error"), exact = FALSE)
    }
    # a whole k held as a double is a k
    expect_identical(ncol(suppressWarnings(suppressMessages(trxRuvFactors(d$y, d$log_cpm_tmm, d$design, genes, 1, "RUVs")))$W), 1L)
    # every group's runs identical: RUVs has no replicate difference to estimate from (RUVSeq fails on an empty index)
    same <- d$y
    same[, 2:6] <- same[, 1L]
    same[, 8:12] <- same[, 7L]
    log_cpm <- suppressWarnings(suppressMessages(trxTmmLogCpm(same, d$design$group)))$log_cpm_tmm
    expect_error(suppressWarnings(RUVSeq::RUVs(trxDesignTestUq(same), genes, 1L, RUVSeq::makeGroups(d$design$group))))
    condition <- tryCatch(suppressWarnings(suppressMessages(trxRuvFactors(same, log_cpm, d$design, genes, 1L, "RUVs"))), error = function(e) e)
    expect_s3_class(condition, c("trx_ruv_k_unavailable", "trx_error"), exact = FALSE)
    expect_match(conditionMessage(condition), "identical", fixed = TRUE)
    # one group with a difference is enough for RUVs to run
    partly <- d$y
    partly[, 2:6] <- partly[, 1L]
    log_cpm <- suppressWarnings(suppressMessages(trxTmmLogCpm(partly, d$design$group)))$log_cpm_tmm
    expect_identical(ncol(suppressWarnings(suppressMessages(trxRuvFactors(partly, log_cpm, d$design, genes, 1L, "RUVs")))$W), 1L)
})

test_that("(KIT-129 audit F-09) unequal groups, three groups and a group of one: the recommended and the supplied designs are the direct RUVSeq and edgeR calls", {
    trxDesignTestPackages()
    layouts <- list(unequal = list(per_group = c(6L, 3L), groups = 2L), three = list(per_group = 4L, groups = 3L),
        one = list(per_group = c(5L, 3L, 1L), groups = 3L), six_one = list(per_group = c(6L, 1L), groups = 2L))
    for (name in names(layouts)) {
        d <- suppressWarnings(suppressMessages(trxDesignTestData(seed = 131L, genes = 900L, per_group = layouts[[name]]$per_group,
            groups = layouts[[name]]$groups)))
        genes <- rownames(d$y)
        uq <- trxDesignTestUq(d$y)
        expect_identical(as.integer(table(d$design$group)), rep_len(layouts[[name]]$per_group, layouts[[name]]$groups), info = name)
        # the recommended design: RUVs with the groups as replicate sets (a group of one has no set)
        sets <- RUVSeq::makeGroups(d$design$group)
        arm <- RUVSeq::RUVs(uq, genes, 1L, sets)
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "groups_all_features", 1L, "user")))
        expect_identical(unname(as.matrix(p$s$W)), unname(arm$W), info = name)
        expect_identical(rownames(p$s$W), d$design$Run, info = name)
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))), info = name)
        # the supplied design: RUVg with the list in the matrix's row order
        supplied <- withr::with_seed(9L, sample(d$null, 40L))
        matched <- genes[genes %in% supplied]
        arm <- RUVSeq::RUVg(uq, matched, 1L)
        p <- suppressWarnings(suppressMessages(trxDesignTestPackage(d, "supplied_controls", 1L, "user", supplied)))
        expect_identical(p$choice$chosen$control_genes, matched, info = name)
        expect_identical(unname(as.matrix(p$s$W)), unname(arm$W), info = name)
        expect_identical(p$table, suppressWarnings(suppressMessages(trxDesignTestEdgeR(d$y, d$design, arm$W))), info = name)
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

test_that("(KIT-129 audit F-01) a second declared design on a RUVs head leaves the first one's factors out of the state and the test", {
    fixture <- trxDesignTestOwner()
    owner <- fixture$owner
    read <- function(name) trxAssayInputs(withWorkflowStateScope(owner, function() readTrxWorkflowState(owner, name)))
    formula <- function() withWorkflowStateScope(owner, function() trxRealisedFormula(owner, getWorkflowStateReference(owner)))
    declared <- function(...) suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, NULL, list(), design_parameters = trxRuvDesignParameters(...))))
    filtered <- read("trx_qc_complete")
    first <- declared("groups_all_features", 2L, "user")
    expect_identical(formula(), "~ 0 + group + W_1 + W_2")
    # none on that head: no factor in the state's design, the formula the DA realises is the record's
    none <- declared("none")
    expect_false(none$value$outcomes$ruv_applied)
    state <- read(owner$current_state)
    expect_identical(names(state$design), names(filtered$design))
    expect_identical(formula(), "~ 0 + group")
    expect_identical(state$counts, filtered$counts)
    da <- suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 1), fixture$workflow$contrasts_tbl)))
    expect_identical(da$manifest$parameters$formula, "~ 0 + group")
    expect_identical(da$manifest$model$w_columns, character())
    # k = 1 after k = 2, on a head that carries two factors again: exactly one W, this run's
    again <- declared("groups_all_features", 2L, "user")
    expect_identical(formula(), "~ 0 + group + W_1 + W_2")
    one <- declared("groups_all_features", 1L, "user")
    expect_identical(one$value$outcomes[c("ruv_k", "ruv_applied")], list(ruv_k = 1L, ruv_applied = TRUE))
    state <- read(owner$current_state)
    expect_identical(names(state$design), c(names(filtered$design), "W_1"))
    arm <- RUVSeq::RUVs(trxDesignTestUq(filtered$counts), rownames(filtered$counts), 1L, RUVSeq::makeGroups(filtered$design$group))
    expect_identical(unname(as.matrix(state$design[, "W_1", drop = FALSE])), unname(arm$W))
    expect_identical(unname(as.matrix(one$value$W)), unname(arm$W))
    expect_identical(formula(), "~ 0 + group + W_1")
})

test_that("(KIT-129 audit F-07) on a session owner, one declared design is one product however its parameters were passed", {
    fixture <- trxDesignTestOwner()
    owner <- fixture$owner
    runs <- suppressWarnings(suppressMessages(withWorkflowStateScope(owner, function() {
        source <- getWorkflowStateReference(owner)
        logcpm <- trxRunProduct(owner, "trx_tmm_logcpm", source, list())
        one <- function(parameters) trxRunRuvDesign(owner, source, logcpm$reference, parameters)
        genes <- rownames(trxAssayInputs(readTrxWorkflowState(owner, owner$current_state))$counts)[1:40]
        list(integer = one(trxRuvDesignParameters("groups_all_features", 2L, "user")),
            double = one(list(ruv_design = "groups_all_features", ruv_k = 2, ruv_k_rule = "user", ruv_control_ids = character())),
            listed = one(trxRuvDesignParameters("supplied_controls", 1L, "user", genes)),
            plain = one(list(ruv_design = "supplied_controls", ruv_k = 1, ruv_k_rule = "user", ruv_control_ids = c(rev(genes), genes[1:3]))),
            fractional = tryCatch(one(list(ruv_design = "groups_all_features", ruv_k = 1.5, ruv_k_rule = "user", ruv_control_ids = character())),
                error = function(e) e))
    })))
    expect_identical(runs$integer$status, "computed")
    expect_identical(runs$double$status, "hit")
    expect_identical(runs$double$reference$entity_id, runs$integer$reference$entity_id)
    expect_identical(runs$listed$status, "computed")
    expect_identical(runs$plain$status, "hit")
    expect_identical(runs$plain$reference$entity_id, runs$listed$reference$entity_id)
    expect_false(identical(runs$listed$reference$entity_id, runs$integer$reference$entity_id))
    expect_s3_class(runs$fractional, c("trx_ruv_design_invalid", "trx_error"), exact = FALSE)
})

