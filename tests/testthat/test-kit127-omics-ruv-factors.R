# KIT-127: lipidomics and metabolomics RUV with the factors in the model (KIT-124's gate, 2026-11-18), and the
# corrected matrix with each feature's level kept. The reference is KIT-124's simulator S4
# (docs/qa/workflow-execution/kit/kit-124/sim/s4_lipid_metab.R at 0b6fcb6d): its arm 7 (all features as controls) and
# arm 8 (internal standards as controls), W in the limma model, the test on the uncorrected values. The `sim*`
# functions below are that script's `ruvW`, `testWithW`, `daKernel` and `longTable`, copied; they call
# `ruv::RUVIII` and the DA kernel, never the functions under test.

kit127Quiet <- function(expr) {
    con <- file(nullfile(), "w")
    sink(con)
    sink(con, type = "message")
    on.exit({
        sink(type = "message")
        sink()
        close(con)
    })
    suppressWarnings(suppressMessages(expr))
}

# Two assays, two groups of six, one sample-wise factor uncorrelated with the group, values under a detection limit
# missing, every feature seen at least once per group (the simulator's presence rule), five standards per assay.
kit127Fixture <- function(seed = 127L) {
    set.seed(seed)
    design <- data.frame(Sample_ID = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L), stringsAsFactors = FALSE)
    n <- nrow(design)
    p <- 60L
    f <- stats::residuals(stats::lm(stats::rnorm(n) ~ factor(design$group)))
    f <- f / sqrt(mean(f^2))
    assays <- lapply(c(pos = "pos", neg = "neg"), function(assay) {
        level <- stats::rnorm(p, 20, 2)
        effect <- c(rep(1, 6L), rep(-1, 6L), rep(0, p - 12L))
        y <- level + outer(effect, as.numeric(design$group == "B")) + outer(stats::rnorm(p, sd = 0.6), f) +
            matrix(stats::rnorm(p * n, sd = 0.5), p)
        seen <- y >= stats::quantile(y, 0.15, names = FALSE)
        kept <- rowSums(seen[, design$group == "A"]) >= 1L & rowSums(seen[, design$group == "B"]) >= 1L
        y[!seen] <- NA_real_
        y <- y[kept, , drop = FALSE]
        dimnames(y) <- list(sprintf("%s_f%03d", assay, seq_len(nrow(y))), design$Sample_ID)
        is <- 22 + outer(stats::rnorm(5L, sd = 0.6), f) + matrix(stats::rnorm(5L * n, sd = 0.25), 5L)
        dimnames(is) <- list(sprintf("%s_ITSD_%d", assay, 1:5), design$Sample_ID)
        list(y = y, is = is)
    })
    list(design = design, assays = assays, contrasts = "groupB-groupA",
        contrasts_tbl = data.frame(contrasts = "groupB-groupA", friendly_names = "B vs A", stringsAsFactors = FALSE))
}

kit127AssayTable <- function(m) {
    table <- tibble::as_tibble(as.data.frame(m, check.names = FALSE))
    table <- tibble::add_column(table, database_identifier = rownames(m), .before = 1L)
    tibble::add_column(table, name = rownames(m), .after = 1L)
}
kit127Lipid <- function(mats, design) {
    methods::new("LipidomicsAssayData", lipid_data = lapply(mats, kit127AssayTable), lipid_id_column = "database_identifier",
        annotation_id_column = "name", design_matrix = design, sample_id = "Sample_ID", group_id = "group", args = list())
}
kit127Metab <- function(mats, design) {
    methods::new("MetaboliteAssayData", metabolite_data = lapply(mats, kit127AssayTable),
        metabolite_id_column = "database_identifier", annotation_id_column = "name", design_matrix = design,
        sample_id = "Sample_ID", group_id = "group", args = list())
}
kit127Values <- function(table, samples) {
    m <- as.matrix(as.data.frame(table[, samples], check.names = FALSE))
    rownames(m) <- table[["database_identifier"]]
    m
}

# ---- the simulator's computation (s4_lipid_metab.R) --------------------------------------------------------------
simLongTable <- function(results) {
    do.call(rbind, lapply(names(results), function(contrast) {
        t <- results[[contrast]]
        data.frame(feature = rownames(t), contrast = contrast, est = t$logFC, p = t$P.Value, q = t$fdr_qvalue,
            bh = t$fdr_value_bh, stringsAsFactors = FALSE)
    }))
}
# `exact`: the exactness round's arm 25 (decision M29; `daKernel` of s4_lipid_metab.R on `sim-k129-exact`): per
# comparison, the model reparameterised so that the comparison is its first coefficient (`limma::contrastAsCoef`).
simDaKernel <- function(m, design, formula, contrasts, exact = FALSE) {
    if (anyNA(m)) {
        x <- stats::model.matrix(stats::as.formula(formula), design)
        seen <- is.finite(m)
        pattern <- drop(seen %*% 2^(seq_len(ncol(m)) - 1))
        estimable <- vapply(split(seq_len(nrow(m)), pattern), function(i) qr(x[seen[i[[1L]], ], , drop = FALSE])$rank == ncol(x), logical(1L))
        m <- m[estimable[as.character(pattern)], , drop = FALSE]
    }
    if (exact) {
        x <- stats::model.matrix(stats::as.formula(formula), design)
        return(do.call(rbind, lapply(contrasts, function(contrast) {
            q <- limma::contrastAsCoef(x, limma::makeContrasts(contrasts = contrast, levels = colnames(x)))$design
            colnames(q) <- c("comparison", paste0("Q", seq_len(ncol(q))[-1L]))
            table <- simDaKernel(m, cbind(design[setdiff(names(design), colnames(q))], q),
                paste("~ 0 +", paste(colnames(q), collapse = " + ")), "comparison")
            table$contrast <- contrast
            table
        })))
    }
    simLongTable(kit127Quiet(runTestsContrastsLipidDA(data = m, contrast_strings = contrasts, design_matrix = design,
        formula_string = formula, sample_id_col = "Sample_ID", treat_lfc_cutoff = NA, eBayes_trend = TRUE,
        eBayes_robust = TRUE))$results)
}
simRuvW <- function(y, m, ctl, k, keep = rownames(y)) {
    ac <- ruv::RUVIII(Y = y, M = m[rownames(y), , drop = FALSE], ctl = ctl, k = k, return.info = TRUE)$fullalpha[seq_len(k), ctl, drop = FALSE]
    w <- sweep(y[, ctl, drop = FALSE], 2L, colMeans(y[, ctl, drop = FALSE])) %*% t(ac) %*% solve(ac %*% t(ac))
    rownames(w) <- rownames(y)
    w <- scale(w[keep, , drop = FALSE], center = TRUE, scale = FALSE)
    colnames(w) <- paste0("W", seq_len(k))
    w
}
simTestWithW <- function(fixture, y, w) {
    design <- cbind(fixture$design, as.data.frame(w[fixture$design$Sample_ID, , drop = FALSE]))
    simDaKernel(y, design, paste("~ 0 + group +", paste(colnames(w), collapse = " + ")), fixture$contrasts, exact = TRUE)
}
simGroupM <- function(design) {
    m <- ruv::replicate.matrix(stats::setNames(design$group, design$Sample_ID))
    rownames(m) <- design$Sample_ID
    m
}
# One assay of an arm: W and the long table. `standards` NULL is arm 7, the standards' matrix arm 8.
simArm <- function(fixture, assay, k, standards = NULL) {
    a <- fixture$assays[[assay]]
    complete <- rowSums(!is.finite(a$y)) == 0L
    y <- t(rbind(a$y[complete, , drop = FALSE], standards))
    ctl <- if (is.null(standards)) rep(TRUE, ncol(y)) else grepl("_ITSD_", colnames(y))
    w <- simRuvW(y, simGroupM(fixture$design), ctl, k, keep = fixture$design$Sample_ID)
    list(W = w, table = simTestWithW(fixture, a$y, w))
}

# The first pass of S4's arms 5 and 23 (`first_pass`, `daKernel`): arm 1's test, per feature the smallest p over the
# comparisons, every feature outside the 25% most significant a control.
simFirstPassControls <- function(fixture, y) {
    base <- simDaKernel(y, fixture$design, "~ 0 + group", fixture$contrasts)
    fp <- stats::aggregate(p ~ feature, base, min)
    fp$feature[rank(fp$p, na.last = "keep", ties.method = "first") > 0.25 * sum(is.finite(fp$p)) & is.finite(fp$p)]
}

expectKit127Table <- function(long, id_column, assay, reference) {
    got <- long[long$assay == assay, ]
    expect_setequal(got[[id_column]], reference$feature)
    at <- match(reference$feature, got[[id_column]])
    expect_equal(got$logFC[at], reference$est, tolerance = 1e-10)
    expect_equal(got$P.Value[at], reference$p, tolerance = 1e-10)
    expect_equal(got$fdr_qvalue[at], reference$q, tolerance = 1e-10)
    expect_equal(got$fdr_value_bh[at], reference$bh, tolerance = 1e-10)
}
expectKit127W <- function(w, reference) {
    expect_identical(dimnames(w), list(rownames(reference), colnames(reference)))
    expect_equal(unname(w), matrix(as.numeric(reference), nrow(reference)), tolerance = 1e-10)
}

K127_K <- list(pos = 2L, neg = 1L)

test_that("arm 7 (all features as controls, W in the model) is the simulator's computation, in lipidomics", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, ruv_design = "groups_all_features", ruv_k = K127_K, ruv_k_rule = "user"))
    da <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    not_estimable <- 0L
    for (assay in names(fixture$assays)) {
        reference <- simArm(fixture, assay, K127_K[[assay]])
        record <- run$factors[[assay]]
        expect_identical(record$status, "ok")
        expectKit127W(record$W, reference$W)
        expectKit127Table(da$da_lipids_long, "lipid_id", assay, reference$table)
        not_estimable <- not_estimable + nrow(fixture$assays[[assay]]$y) - nrow(reference$table)
        expect_identical(record$k, K127_K[[assay]])
        expect_identical(record$ruv_k_rule, "user")
        expect_identical(record$control_count, sum(rowSums(is.na(fixture$assays[[assay]]$y)) == 0L))
        expect_identical(record$not_tested, character())
    }
    # the fixture has features the model with W cannot estimate: they are not tested, as in the simulator
    expect_gt(not_estimable, 0L)
    expect_identical(run$parameters$ruv_design, "groups_all_features")
    expect_identical(run$parameters$ruv_grouping_variable, "group")
})

test_that("arm 8 (the standards as controls, W in the model, the standards not tested) is the simulator's computation", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    ids <- lapply(fixture$assays, function(a) rownames(a$is))
    run <- kit127Quiet(runOmicsRuvFactors(object, ruv_design = "supplied_controls", ruv_k = K127_K,
        ruv_k_rule = "user", ruv_control_ids = ids))
    da <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    for (assay in names(fixture$assays)) {
        reference <- simArm(fixture, assay, K127_K[[assay]], standards = fixture$assays[[assay]]$is)
        record <- run$factors[[assay]]
        expectKit127W(record$W, reference$W)
        expectKit127Table(da$da_lipids_long, "lipid_id", assay, reference$table)
        expect_identical(record$control_count, 5L)
        expect_setequal(record$control_ids, ids[[assay]])
        expect_setequal(record$not_tested, ids[[assay]])
        expect_false(any(ids[[assay]] %in% da$da_lipids_long$lipid_id))
        # the standards stay in the corrected matrix
        expect_true(all(ids[[assay]] %in% run$corrected@lipid_data[[assay]]$database_identifier))
    }
})

test_that("metabolomics gives the factors and the tables lipidomics gives", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    lipid <- kit127Lipid(mats, fixture$design)
    metab <- kit127Metab(mats, fixture$design)
    run_l <- kit127Quiet(runOmicsRuvFactors(lipid, "groups_all_features", K127_K, "user"))
    run_m <- kit127Quiet(runOmicsRuvFactors(metab, "groups_all_features", K127_K, "user"))
    expect_identical(run_m$factors, run_l$factors)
    expect_identical(run_m$corrected@metabolite_data, run_l$corrected@lipid_data)
    da_l <- kit127Quiet(runLipidsDA(lipid, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run_l$factors))$da_lipids_long
    da_m <- kit127Quiet(runMetabolitesDA(metab, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run_m$factors))$da_metabolites_long
    cols <- c("logFC", "P.Value", "fdr_qvalue", "fdr_value_bh", "assay", "comparison", "significant")
    expect_identical(da_m$metabolite_id, da_l$lipid_id)
    expect_identical(da_m[cols], da_l[cols])
    for (assay in names(fixture$assays)) {
        expectKit127Table(da_m, "metabolite_id", assay, simArm(fixture, assay, K127_K[[assay]])$table)
    }
})

test_that("the corrected matrix of the new path keeps each feature's level; today's correction keeps it too (KIT-224)", {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    samples <- fixture$design$Sample_ID
    for (make in list(kit127Lipid, kit127Metab)) {
        object <- make(mats, fixture$design)
        slot_name <- if (methods::is(object, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"
        run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", K127_K, "user"))
        legacy <- kit127Quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = K127_K,
            ctrl = lapply(mats, rownames)))
        # KIT-224: the earlier rule (records made before) keeps the shift
        legacy_shifted <- withr::with_options(list(MultiScholaR.ruv_feature_levels = "shifted"),
            kit127Quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = K127_K,
                ctrl = lapply(mats, rownames))))
        for (assay in names(mats)) {
            before <- mats[[assay]]
            after <- kit127Values(methods::slot(run$corrected, slot_name)[[assay]], samples)
            today <- kit127Values(methods::slot(legacy, slot_name)[[assay]], samples)
            earlier <- kit127Values(methods::slot(legacy_shifted, slot_name)[[assay]], samples)
            expect_identical(rownames(after), rownames(today))
            kept <- rownames(after)
            both <- is.finite(before[kept, ]) & is.finite(after)
            expect_identical(is.finite(after), is.finite(today))
            mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
            expect_equal(mean_in(after), mean_in(before[kept, ]), tolerance = 1e-9)
            # the correction itself is today's: the two differ by a constant per feature
            expect_equal((after - mean_in(after))[both], (today - mean_in(today))[both], tolerance = 1e-9)
            expect_gt(max(abs(after - before[kept, ])[both]), 0.01)
            # today's legacy path restores too (KIT-224); under the earlier rule its level shift is still there
            expect_equal(mean_in(today), mean_in(before[kept, ]), tolerance = 1e-9)
            expect_gt(max(abs(mean_in(earlier) - mean_in(before[kept, ]))), 0.01)
        }
        # the object the DA reads is not touched
        expect_identical(methods::slot(object, slot_name), lapply(mats, kit127AssayTable))
    }
})

test_that("without ruv_factors the DA entry points are today's: the unchanged kernel on the object's own matrix", {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    cases <- list(list(object = kit127Lipid(mats, fixture$design), run = runLipidsDA, long = "da_lipids_long", id = "lipid_id",
        kernel = runTestsContrastsLipidDA),
        list(object = kit127Metab(mats, fixture$design), run = runMetabolitesDA, long = "da_metabolites_long", id = "metabolite_id",
            kernel = runTestsContrastsMetabDA))
    for (case in cases) {
        default <- kit127Quiet(case$run(case$object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))
        explicit <- kit127Quiet(case$run(case$object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = NULL))
        expect_identical(default[[case$long]], explicit[[case$long]])
        expect_identical(names(default), c("theObject", "contrasts_results", case$long, "per_assay_results", "significant_counts",
            "qvalue_warnings", "da_q_val_thresh", "treat_lfc_cutoff"))
        expect_gt(nrow(default[[case$long]]), 0L)
        for (assay in names(mats)) {
            direct <- kit127Quiet(case$kernel(data = mats[[assay]], contrast_strings = fixture$contrasts, design_matrix = fixture$design,
                formula_string = "~ 0 + group", sample_id_col = "Sample_ID", treat_lfc_cutoff = NA, eBayes_trend = TRUE,
                eBayes_robust = TRUE))$results[[fixture$contrasts]]
            got <- default[[case$long]][default[[case$long]]$assay == assay, ]
            # the entry point's table has its own row order: compared by feature
            expect_setequal(got[[case$id]], rownames(direct))
            at <- match(rownames(direct), got[[case$id]])
            expect_identical(got$P.Value[at], direct$P.Value)
            expect_identical(got$logFC[at], direct$logFC)
            expect_identical(got$fdr_qvalue[at], direct$fdr_qvalue)
        }
    }
})

test_that("k has no default, 0 is no RUV, and an assay whose factors cannot be computed gets no RUV with a warning", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    object <- kit127Lipid(mats, fixture$design)
    expect_error(runOmicsRuvFactors(object, "groups_all_features"), class = "omics_ruv_factors_invalid")
    expect_error(runOmicsRuvFactors(object, "groups_all_features", 2L), class = "omics_ruv_factors_invalid")
    expect_error(runOmicsRuvFactors(object, ruv_k = 2L, ruv_k_rule = "user"), class = "omics_ruv_factors_invalid")
    plain <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_lipids_long

    # k = 0: no factors, the corrected object is the input, the test is the plain test
    none <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 0L, "user"))
    expect_identical(vapply(none$factors, `[[`, "", "status"), c(pos = "none", neg = "none"))
    expect_null(none$factors$pos$W)
    expect_identical(none$corrected@lipid_data, object@lipid_data)
    expect_identical(kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = none$factors))$da_lipids_long, plain)

    # 12 samples in 2 groups leave 9 factors at most. Decision M31: a typed k of 10 for pos stops the run, as in classic
    # proteomics and transcriptomics (before: pos analysed without RUV, with a warning)
    typed <- tryCatch(suppressMessages(runOmicsRuvFactors(object, "groups_all_features", list(pos = 10L, neg = 1L), "user")),
        error = function(e) e)
    expect_s3_class(typed, c("ruv_k_above_limit", "ruv_k_refused"))
    expect_match(conditionMessage(typed), "pos", fixed = TRUE)
    expect_match(conditionMessage(typed), "residual", fixed = TRUE)
    # under an automatic rule (the separation rule, first-pass controls) pos is refused and analysed without RUV, neg
    # is corrected
    expect_warning(refused <- suppressMessages(runOmicsRuvFactors(object, "first_pass_nulls", list(pos = 10L, neg = 1L),
        "sweep_bounded")), class = "omics_ruv_factors_refused")
    expect_identical(refused$factors$pos$status, "refused")
    expect_match(refused$factors$pos$reason, "residual")
    # KIT-129 item 2: the refused assay names its refusal by the class the other omics raise
    expect_identical(refused$factors$pos$reason_class, "ruv_k_above_limit")
    expect_identical(refused$factors$neg$reason_class, NA_character_)
    warned <- tryCatch(suppressMessages(runOmicsRuvFactors(object, "first_pass_nulls", list(pos = 10L, neg = 1L), "sweep_bounded")),
        omics_ruv_factors_refused = function(w) w)
    expect_match(conditionMessage(warned), "(ruv_k_above_limit)", fixed = TRUE)
    expect_identical(refused$factors$pos$k, 0L)
    expect_identical(refused$factors$pos$k_requested, 10L)
    expect_null(refused$factors$pos$W)
    expect_identical(refused$factors$neg$status, "ok")
    expect_identical(refused$corrected@lipid_data$pos, object@lipid_data$pos)
    expect_false(identical(refused$corrected@lipid_data$neg, object@lipid_data$neg))
    mixed <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = refused$factors))$da_lipids_long
    expect_identical(mixed[mixed$assay == "pos", ], plain[plain$assay == "pos", ])
    neg <- fixture$assays$neg$y
    y <- t(neg[rowSums(!is.finite(neg)) == 0L, , drop = FALSE])
    w <- simRuvW(y, simGroupM(fixture$design), colnames(y) %in% simFirstPassControls(fixture, neg), 1L, keep = fixture$design$Sample_ID)
    expectKit127Table(mixed, "lipid_id", "neg", simTestWithW(fixture, neg, w))

    # more factors than controls observed in every sample: typed, the run stops (M31); automatic, the assay is refused
    with_standards <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    expect_error(suppressMessages(runOmicsRuvFactors(with_standards, "supplied_controls", 6L, "user",
        ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is)))), "controls observed in every sample", class = "ruv_k_above_limit")
    expect_warning(few <- suppressMessages(runOmicsRuvFactors(with_standards, "supplied_controls", 6L, "sweep_bounded",
        ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is)))), class = "omics_ruv_factors_refused")
    expect_identical(unname(vapply(few$factors, `[[`, "", "status")), c("refused", "refused"))
    expect_match(few$factors$pos$reason, "controls observed in every sample")
    expect_identical(few$factors$pos$reason_class, "ruv_k_above_limit")

    # a k that could not be estimated is a refusal, not a guess
    expect_warning(unknown <- suppressMessages(runOmicsRuvFactors(with_standards, "supplied_controls",
        list(pos = NA_integer_, neg = 1L), "sweep_bounded", ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is)))),
        class = "omics_ruv_factors_refused")
    expect_identical(unknown$factors$pos$status, "refused")
    expect_identical(unknown$factors$pos$reason_class, "ruv_k_not_estimated")
    expect_identical(unknown$factors$neg$ruv_k_rule, "sweep_bounded")
})

test_that("wrong arguments stop", {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    object <- kit127Lipid(mats, fixture$design)
    ids <- lapply(mats, function(m) rownames(m)[1:5])
    invalid <- function(...) expect_error(kit127Quiet(runOmicsRuvFactors(object, ...)), class = "omics_ruv_factors_invalid")
    invalid("legacy", 1L, "user")
    invalid("groups_all_features", 1.5, "user")
    invalid("groups_all_features", -1L, "user")
    invalid("groups_all_features", list(pos = 1L), "user")
    invalid("groups_all_features", 1L, "guess")
    invalid("groups_all_features", 1L, "user", ruv_control_ids = ids)
    invalid("supplied_controls", 1L, "user")
    invalid("supplied_controls", 1L, "user", ruv_control_ids = list(pos = ids$pos))
    invalid("supplied_controls", 1L, "user", ruv_control_ids = list(pos = "absent", neg = ids$neg))
    invalid("groups_all_features", 1L, "user", ruv_grouping_variable = "absent")
    expect_error(runOmicsRuvFactors(fixture$design, "groups_all_features", 1L, "user"), class = "omics_ruv_factors_invalid")
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user"))
    # the DA refuses factors that do not cover an assay, or a design that already has a W column
    expect_error(kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = run$factors["pos"])), class = "omics_ruv_factors_invalid")
    clash <- object
    clash@design_matrix$W1 <- 0
    expect_error(kit127Quiet(runLipidsDA(clash, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = run$factors)), class = "omics_ruv_factors_invalid")
})

test_that("today's estimate of k is available and bounded by the sample-size limit", {
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    sweep <- list(pos = list(success = TRUE, best_k = 5L), neg = list(success = TRUE, best_k = 1L))
    bounded <- boundOmicsRuvSweepK(object, sweep)
    # 12 samples: `calculateAdaptiveMaxK` gives 2
    expect_identical(bounded$k, list(pos = 2L, neg = 1L))
    expect_identical(bounded$sweep_k, list(pos = 5L, neg = 1L))
    expect_identical(bounded$limit, list(pos = 2L, neg = 2L))
    expect_identical(bounded$ruv_k_rule, "sweep_bounded")
    failed <- boundOmicsRuvSweepK(object, list(pos = list(success = FALSE, error = "no valid percentage"), neg = sweep$neg))
    expect_identical(failed$k, list(pos = NA_integer_, neg = 1L))
    expect_error(boundOmicsRuvSweepK(object, sweep["pos"]), class = "omics_ruv_factors_invalid")

    skip_if_not_installed("ruv")
    params <- list(ruv_grouping_variable = "group", percentage_min = 10, percentage_max = 12, separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE)
    real <- boundOmicsRuvSweepK(object, kit127Quiet(runLipidPerAssayRuvOptimization(object, "automatic", params)))
    known <- !is.na(unlist(real$k))
    expect_true(all(unlist(real$k)[known] <= 2L))
    # the bounded sweep's k goes with supplied controls only (audit F1)
    with_standards <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(with_standards, "supplied_controls", real$k, real$ruv_k_rule,
        ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is))))
    expect_identical(unname(vapply(run$factors, `[[`, "", "ruv_k_rule")), c("sweep_bounded", "sweep_bounded"))
    expect_true(all(vapply(run$factors, `[[`, "", "status") %in% c("ok", "refused")))
})

# ---- the independent audit of 2026-11-18 (KIT-127, findings F2 to F11) ---------------------------------------------

# As `kit127Quiet`, but a warning of class `keep` reaches the caller (the other warnings are limma's own notices).
kit127Hush <- function(expr, keep) {
    con <- file(nullfile(), "w")
    sink(con)
    sink(con, type = "message")
    on.exit({
        sink(type = "message")
        sink()
        close(con)
    })
    withCallingHandlers(suppressMessages(expr), warning = function(w) if (!inherits(w, keep)) invokeRestart("muffleWarning"))
}

# One assay without missing values but the ones set here: `groups` of `n_per` samples, a batch column, one sample-wise
# factor. `no_c`: unobserved in the whole last group, complete elsewhere. `no_c_gaps`: the same with two more values
# missing. `sparse`: one value per group, too few to estimate a W coefficient.
kit127DesignFixture <- function(groups, n_per, p = 40L, seed = 1127L, patterns = FALSE) {
    set.seed(seed)
    n <- length(groups) * n_per
    design <- data.frame(Sample_ID = sprintf("S%02d", seq_len(n)), group = rep(groups, each = n_per),
        batch = rep(c("b1", "b2"), length.out = n), stringsAsFactors = FALSE)
    f <- stats::residuals(stats::lm(stats::rnorm(n) ~ factor(design$group)))
    y <- stats::rnorm(p, 20, 2) + outer(c(rep(1, 8L), rep(0, p - 8L)), match(design$group, groups) - 1) +
        outer(stats::rnorm(p, sd = 0.6), f / sqrt(mean(f^2))) + matrix(stats::rnorm(p * n, sd = 0.4), p)
    dimnames(y) <- list(sprintf("f%03d", seq_len(p)), design$Sample_ID)
    sets <- list()
    if (patterns) {
        last <- design$group == groups[[length(groups)]]
        sets <- list(no_c = rownames(y)[c(1:3, 21:23)], no_c_gaps = rownames(y)[c(4:5, 24:25)], sparse = rownames(y)[26:27])
        y[c(sets$no_c, sets$no_c_gaps), last] <- NA_real_
        y[sets$no_c_gaps, c(1L, n_per + 2L)] <- NA_real_
        y[sets$sparse, -(n_per * (seq_along(groups) - 1L) + 1L)] <- NA_real_
    }
    c(list(design = design, y = y), sets)
}
kit127Contrasts <- function(x) data.frame(contrasts = x, friendly_names = x, stringsAsFactors = FALSE)

test_that("F2: a model with no residual degree of freedom once W is added stops the DA, and never returns an empty table", {
    skip_if_not_installed("ruv")
    fx <- kit127DesignFixture(c("A", "B"), 4L)
    object <- kit127Lipid(list(pos = fx$y), fx$design)
    # 8 samples in 2 replicate sets: k = 5 is the kernel's bound; with the batch covariate the model has 8 columns
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 5L, "user"))
    expect_identical(run$factors$pos$status, "ok")
    plain <- kit127Quiet(runLipidsDA(object, kit127Contrasts("groupB-groupA"), "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    expect_identical(nrow(plain$da_lipids_long), nrow(fx$y))
    expect_error(kit127Quiet(runLipidsDA(object, kit127Contrasts("groupB-groupA"), "~ 0 + group + batch", 0.05, 0, TRUE, TRUE,
        ruv_factors = run$factors)), class = "omics_ruv_no_residual_df")
    expect_error(kit127Quiet(runMetabolitesDA(kit127Metab(list(pos = fx$y), fx$design), kit127Contrasts("groupB-groupA"),
        "~ 0 + group + batch", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors)), class = "omics_ruv_no_residual_df")
    # one factor fewer leaves a degree of freedom: the covariate formula is fitted
    run4 <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 4L, "user"))
    covariate <- kit127Quiet(runLipidsDA(object, kit127Contrasts("groupB-groupA"), "~ 0 + group + batch", 0.05, 0, TRUE, TRUE,
        ruv_factors = run4$factors))
    expect_identical(nrow(covariate$da_lipids_long), nrow(fx$y))
    expect_identical(colnames(covariate$contrasts_results$pos$fit.eb$design), c("groupA", "groupB", "batchb2", paste0("W", 1:4)))
})

test_that("F3: the corrected object is stamped, and a stamped object is refused together with factors", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    for (case in list(list(make = kit127Lipid, run = runLipidsDA), list(make = kit127Metab, run = runMetabolitesDA))) {
        object <- case$make(mats, fixture$design)
        run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", K127_K, "user"))
        expect_identical(run$corrected@args$ruv_design, "groups_all_features")
        expect_identical(run$corrected@args$ruv_number_k, K127_K)
        expect_identical(run$corrected@args$ruv_grouping_variable, "group")
        expect_identical(run$corrected@args$ctrl, lapply(mats, rownames))
        expect_null(object@args$ruv_design)
        expect_error(kit127Quiet(case$run(run$corrected, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
            ruv_factors = run$factors)), class = "omics_ruv_object_already_corrected")
        # without factors the corrected object is tested as any object is (today's design)
        expect_gt(nrow(kit127Quiet(case$run(run$corrected, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))[[3L]]), 0L)
    }
    # no assay corrected: nothing to stamp
    none <- kit127Quiet(runOmicsRuvFactors(kit127Lipid(mats, fixture$design), "groups_all_features", 0L, "user"))
    expect_null(none$corrected@args$ruv_design)
})

test_that("F5: a feature unobserved in one whole group is tested between the groups where it is observed, as a two-group fit", {
    skip_if_not_installed("ruv")
    fx <- kit127DesignFixture(c("A", "B", "C"), 6L, patterns = TRUE)
    object <- kit127Lipid(list(pos = fx$y), fx$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 2L, "user"))
    expect_warning(da <- kit127Hush(runLipidsDA(object, kit127Contrasts(c("groupB-groupA", "groupC-groupA")), "~ 0 + group", 0.05, 0,
        TRUE, TRUE, ruv_factors = run$factors), "omics_ruv_features_not_tested"), class = "omics_ruv_features_not_tested")
    long <- da$da_lipids_long
    no_c <- c(fx$no_c, fx$no_c_gaps)
    # only the features whose W coefficients cannot be estimated are left out, and they are named (F6)
    expect_identical(da$ruv_features_not_tested$pos$not_estimable, fx$sparse)
    expect_identical(da$ruv_features_not_tested$pos$supplied_controls, character())
    expect_setequal(long$lipid_id, setdiff(rownames(fx$y), fx$sparse))
    # flagged as not observed in group C
    flagged <- da$ruv_features_not_observed$pos
    expect_identical(flagged[c("feature", "group")], data.frame(feature = no_c[order(no_c)], group = "C", stringsAsFactors = FALSE))
    at <- function(contrast) match(paste(no_c, contrast), paste(long$lipid_id, long$comparison))
    expect_true(all(is.finite(long$logFC[at("groupB-groupA")]) & is.finite(long$P.Value[at("groupB-groupA")])))
    expect_true(all(is.na(long$logFC[at("groupC-groupA")]) & is.na(long$P.Value[at("groupC-groupA")])))

    # the guard: per feature, least squares on the samples of A and B only, with the same W
    eb <- da$contrasts_results$pos$fit.eb
    w <- run$factors$pos$W[fx$design$Sample_ID, , drop = FALSE]
    x <- cbind(groupA = as.numeric(fx$design$group == "A"), groupB = as.numeric(fx$design$group == "B"), w)
    reference <- t(vapply(no_c, function(feature) {
        s <- is.finite(fx$y[feature, ]) & fx$design$group %in% c("A", "B")
        fit <- stats::lm.fit(x[s, , drop = FALSE], fx$y[feature, s])
        contrast <- c(-1, 1, 0, 0)
        df <- sum(s) - ncol(x)
        unscaled <- sqrt(drop(contrast %*% solve(crossprod(x[s, , drop = FALSE])) %*% contrast))
        row <- match(feature, rownames(eb$coefficients))
        t_ref <- sum(contrast * fit$coefficients) / (sqrt(eb$s2.post[[row]]) * unscaled)
        c(logFC = sum(contrast * fit$coefficients), sigma = sqrt(sum(fit$residuals^2) / df), df = df, unscaled = unscaled,
            p = 2 * stats::pt(-abs(t_ref), eb$df.total[[row]]))
    }, numeric(5L)))
    rows <- match(no_c, rownames(eb$coefficients))
    # exact: the estimate, the residual standard deviation and its degrees of freedom
    expect_equal(long$logFC[at("groupB-groupA")], unname(reference[, "logFC"]), tolerance = 1e-10)
    expect_equal(unname(eb$sigma[rows]), unname(reference[, "sigma"]), tolerance = 1e-10)
    expect_equal(unname(eb$df.residual[rows]), unname(reference[, "df"]))
    # exact too (decision M29): the contrast's standard error, per feature on its own observed samples. limma's
    # contrast step alone (`limma::contrasts.fit`) formed it from the full design's correlations, 0.15 per cent off
    # on this fixture for the features complete in A and B.
    ratio <- unname(eb$stdev.unscaled[rows, "groupB-groupA"] / reference[, "unscaled"])
    expect_equal(ratio, rep(1, length(no_c)), tolerance = 1e-10)
    expect_equal(long$P.Value[at("groupB-groupA")], unname(reference[, "p"]), tolerance = 1e-10)
})

test_that("F6: the features not tested are returned by assay and announced by a classed warning", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, function(a) rbind(a$y, a$is))
    ids <- lapply(fixture$assays, function(a) rownames(a$is))
    object <- kit127Metab(mats, fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, "supplied_controls", K127_K, "user", ruv_control_ids = ids))
    expect_warning(da <- kit127Hush(runMetabolitesDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = run$factors), "omics_ruv_features_not_tested"), class = "omics_ruv_features_not_tested")
    expect_named(da$ruv_features_not_tested, names(mats))
    total <- 0L
    for (assay in names(mats)) {
        tested <- da$da_metabolites_long$metabolite_id[da$da_metabolites_long$assay == assay]
        record <- da$ruv_features_not_tested[[assay]]
        expect_setequal(record$supplied_controls, ids[[assay]])
        expect_setequal(record$not_estimable, setdiff(rownames(fixture$assays[[assay]]$y), tested))
        total <- total + length(record$not_estimable)
    }
    expect_gt(total, 0L)
    # today's path returns no such field
    expect_null(kit127Quiet(runMetabolitesDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))$ruv_features_not_tested)
})

test_that("F9: an incomplete supplied control is left out of W and of the correction alike, with a classed warning", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, function(a) rbind(a$y, a$is))
    mats$pos["pos_ITSD_2", "S03"] <- NA_real_
    ids <- lapply(fixture$assays, function(a) rownames(a$is))
    complete <- list(pos = setdiff(ids$pos, "pos_ITSD_2"), neg = ids$neg)
    object <- kit127Lipid(mats, fixture$design)
    expect_warning(run <- kit127Hush(runOmicsRuvFactors(object, "supplied_controls", K127_K, "user", ruv_control_ids = ids),
        "omics_ruv_controls_incomplete"), class = "omics_ruv_controls_incomplete", regexp = "pos_ITSD_2")
    expect_identical(run$factors$pos$control_ids, complete$pos)
    expect_identical(run$factors$pos$controls_left_out, "pos_ITSD_2")
    expect_identical(run$factors$neg$controls_left_out, character())
    expect_setequal(run$factors$pos$not_tested, ids$pos)
    expect_identical(run$corrected@args$ctrl, complete)
    # the correction is the one the complete controls give
    reference <- kit127Quiet(runOmicsRuvFactors(object, "supplied_controls", K127_K, "user", ruv_control_ids = complete))
    expect_identical(run$corrected@lipid_data, reference$corrected@lipid_data)
    expect_identical(run$factors$pos$W, reference$factors$pos$W)
})

test_that("F10: a non-finite k, a missing group and the factors of another assay are refused", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    object <- kit127Lipid(mats, fixture$design)
    expect_error(kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", Inf, "user")), class = "omics_ruv_factors_invalid")
    expect_error(kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", list(pos = NaN, neg = 1L), "user")),
        class = "omics_ruv_factors_invalid")
    no_group <- fixture$design
    no_group$group[[2L]] <- NA_character_
    expect_error(kit127Quiet(runOmicsRuvFactors(kit127Lipid(mats, no_group), "groups_all_features", 1L, "user")),
        class = "omics_ruv_factors_invalid")
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user"))
    expect_identical(vapply(run$factors, `[[`, "", "assay"), c(pos = "pos", neg = "neg"))
    expect_error(kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = list(pos = run$factors$pos, neg = run$factors$pos))), class = "omics_ruv_factors_invalid")
})

test_that("F1: the bounded sweep's k is refused with all features as controls", {
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    expect_error(kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "sweep_bounded")), class = "omics_ruv_factors_invalid")
})

test_that("F11: the test with W in the model is limma's own, computed here without the package's DA kernel", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", K127_K, "user"))
    da <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    for (assay in names(fixture$assays)) {
        y <- fixture$assays[[assay]]$y
        x <- cbind(groupA = as.numeric(fixture$design$group == "A"), groupB = as.numeric(fixture$design$group == "B"),
            run$factors[[assay]]$W[fixture$design$Sample_ID, , drop = FALSE])
        estimable <- apply(is.finite(y), 1L, function(seen) qr(x[seen, , drop = FALSE])$rank == ncol(x))
        # the comparison as the first coefficient of the reparameterised model: its exact standard error (M29)
        exact <- limma::contrastAsCoef(x, limma::makeContrasts(contrasts = "groupB-groupA", levels = colnames(x)))$design
        fit <- suppressWarnings(limma::lmFit(y[estimable, ], exact))[, 1L]
        top <- limma::topTable(limma::eBayes(fit, trend = TRUE, robust = TRUE), coef = 1L, number = Inf, sort.by = "none")
        expectKit127Table(da$da_lipids_long, "lipid_id", assay, data.frame(feature = rownames(top), est = top$logFC, p = top$P.Value,
            q = qvalue::qvalue(top$P.Value)$qvalues, bh = stats::p.adjust(top$P.Value, "BH")))
    }
})

test_that("F11: metabolomics with supplied standards is the simulator's arm 8", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Metab(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    ids <- lapply(fixture$assays, function(a) rownames(a$is))
    run <- kit127Quiet(runOmicsRuvFactors(object, "supplied_controls", K127_K, "user", ruv_control_ids = ids))
    da <- kit127Quiet(runMetabolitesDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    for (assay in names(fixture$assays)) {
        reference <- simArm(fixture, assay, K127_K[[assay]], standards = fixture$assays[[assay]]$is)
        expectKit127W(run$factors[[assay]]$W, reference$W)
        expectKit127Table(da$da_metabolites_long, "metabolite_id", assay, reference$table)
        expect_false(any(ids[[assay]] %in% da$da_metabolites_long$metabolite_id))
        expect_true(all(ids[[assay]] %in% run$corrected@metabolite_data[[assay]]$database_identifier))
    }
})

test_that("F11: runMetabolitesDA refuses factors that are not the kernel's per-assay records", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Metab(lapply(fixture$assays, `[[`, "y"), fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user"))
    refused <- function(factors) {
        expect_error(kit127Quiet(runMetabolitesDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = factors)),
            class = "omics_ruv_factors_invalid")
    }
    refused(run)
    refused(run$factors["pos"])
    refused(list(pos = run$factors$pos, neg = run$factors$pos))
    short <- run$factors
    short$neg$W <- short$neg$W[-1L, , drop = FALSE]
    refused(short)
})

test_that("F11: a correction that returns an assay unchanged stops the run", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    local_mocked_bindings(ruvIII_C_Varying = function(theObject, ...) theObject)
    expect_error(kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", 1L, "user")), class = "omics_ruv_correction_failed")
})

test_that("F11: the factors record survives a save and a restore", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    run <- kit127Quiet(runOmicsRuvFactors(object, "supplied_controls", K127_K, "user",
        ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is))))
    path <- withr::local_tempfile(fileext = ".rds")
    saveRDS(run[c("factors", "parameters")], path)
    restored <- readRDS(path)
    expect_identical(restored, run[c("factors", "parameters")])
    da <- function(factors) kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = factors))
    expect_identical(da(restored$factors)[c("da_lipids_long", "ruv_features_not_tested")],
        da(run$factors)[c("da_lipids_long", "ruv_features_not_tested")])
})

# KIT-129 B1 (the rehearsal's probe a): 2 groups of 4 biological units, each injected twice, the units declared as
# technical replicates; two assays.
techRepRuvObject <- function(omic, seed = 11L, features = 250L) {
    groups <- rep(c("A", "B"), each = 8L)
    unit <- rep(rep(1:4, each = 2L), times = 2L)
    d <- data.frame(Run = sprintf("%s%d_%d", groups, unit, rep(1:2, times = 8L)), group = groups, replicates = unit,
        tech_rep_group = paste(groups, unit, sep = "_"), stringsAsFactors = FALSE)
    id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
    assay <- function(seed) {
        u <- match(d$tech_rep_group, unique(d$tech_rep_group))
        values <- withr::with_seed(seed, {
            f <- stats::residuals(stats::lm(stats::rnorm(16L) ~ factor(d$group)))
            stats::rnorm(features, 20, 2) + matrix(stats::rnorm(features * 8L, 0, sqrt(0.1875)), nrow = features)[, u] +
                outer(stats::rnorm(features, sd = 0.5), f / sqrt(mean(f^2))) + matrix(stats::rnorm(features * 16L, 0, 0.25), nrow = features)
        })
        colnames(values) <- d$Run
        table <- data.frame(sprintf("F%d_%03d", seed, seq_len(features)), paste("Name", seq_len(features)), stringsAsFactors = FALSE)
        names(table) <- c(id, "name")
        tibble::as_tibble(cbind(table, as.data.frame(values, check.names = FALSE)))
    }
    assays <- list(pos = assay(seed), neg = assay(seed + 1000L))
    if (identical(omic, "lipid")) {
        methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = d, sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", lipid_id_column = "lipid_id", annotation_id_column = "name")
    } else {
        methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = d, sample_id = "Run", group_id = "group",
            technical_replicate_id = "tech_rep_group", metabolite_id_column = "metabolite_id", annotation_id_column = "name")
    }
}

test_that("with RUV factors, a block correlation that cannot be estimated is named with k in the failure (KIT-129 B1)", {
    # The rehearsal reached it at k = 12 of the injection bound; since M28 such a k is refused before the fit (below), so
    # the wrapper is shown on the condition itself.
    factors <- list(W = matrix(0, 16L, 3L))
    refused <- tryCatch(withOmicsRuvKOnBlockFailure(rlang::abort("the correlation could not be estimated",
        class = "omics_da_block_correlation_failed"), factors), error = function(e) e)
    expect_s3_class(refused, c("omics_ruv_k_block_correlation_failed", "omics_da_block_correlation_failed"))
    expect_match(conditionMessage(refused), "with k = 3 RUV factors in the model, the correlation could not be estimated", fixed = TRUE)
    # without factors, and for any other failure, the condition passes unchanged
    plain <- tryCatch(withOmicsRuvKOnBlockFailure(rlang::abort("x", class = "omics_da_block_correlation_failed"), NULL), error = function(e) e)
    expect_identical(conditionMessage(plain), "x")
    other <- tryCatch(withOmicsRuvKOnBlockFailure(stop("y"), factors), error = function(e) e)
    expect_identical(conditionMessage(other), "y")
    expect_identical(withOmicsRuvKOnBlockFailure(1 + 1, factors), 2)
})

test_that("with technical replicates the biological units limit k, in the kernel and at the test (M28)", {
    contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    for (omic in c("lipid", "metab")) {
        run <- if (identical(omic, "lipid")) runLipidsDA else runMetabolitesDA
        object <- techRepRuvObject(omic)
        kernel <- function(k, obj = object) suppressMessages(runOmicsRuvFactors(obj, "groups_all_features", k, "user"))
        # 8 units in 2 groups: at most 5 factors (the injections' bound was 13); a typed k = 6 stops the run (M31), never
        # reduced
        refused <- tryCatch(kernel(6L), error = function(e) e)
        expect_s3_class(refused, c("ruv_k_above_limit", "ruv_k_refused"))
        expect_match(conditionMessage(refused), "limit of 5: the design has technical replicates", fixed = TRUE, info = omic)
        within <- kernel(5L)
        expect_identical(vapply(within$factors, `[[`, "", "status"), c(pos = "ok", neg = "ok"), info = omic)
        da <- function(factors, formula = "~ 0 + group") {
            utils::capture.output(out <- tryCatch(suppressWarnings(suppressMessages(run(object, contrasts, formula, 0.05, 0,
                TRUE, TRUE, ruv_factors = factors))), error = function(e) e))
            out
        }
        fitted <- da(within$factors)
        expect_false(inherits(fitted, "error"), info = omic)
        expect_true(isTRUE(fitted$contrasts_results$pos$fit.eb$uses_duplicate_correlation), info = omic)
        # at the test: factors made without the declaration (the injections' bound) are refused between the units
        undeclared <- object
        undeclared@technical_replicate_id <- NA_character_
        many <- kernel(12L, undeclared)
        expect_identical(vapply(many$factors, `[[`, 0L, "k"), c(pos = 12L, neg = 12L), info = omic)
        stopped <- da(many$factors)
        expect_s3_class(stopped, "omics_ruv_no_residual_df")
        expect_match(conditionMessage(stopped), "between the biological units", fixed = TRUE, info = omic)
        # without the declaration the limit is the injections' one, as before
        expect_identical(vapply(kernel(12L, undeclared)$factors, `[[`, "", "status"), c(pos = "ok", neg = "ok"), info = omic)
    }
})

# The legacy DA (no RUV factors) on fixed inputs (decision M29: legacy unchanged). The digests were recorded on the tree
# before M29 was implemented (`kit-129-integration` d4b5db27).
kit127LegacyDigests <- function() {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    lipid <- kit127Quiet(runLipidsDA(kit127Lipid(mats, fixture$design), fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))
    metab <- kit127Quiet(runMetabolitesDA(kit127Metab(mats, fixture$design), fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))
    c(lipid = digest::digest(list(lipid$contrasts_results, lipid$da_lipids_long)),
        metab = digest::digest(list(metab$contrasts_results, metab$da_metabolites_long)))
}

test_that("the legacy DA without RUV factors is unchanged by the exact standard errors (M29)", {
    expect_identical(kit127LegacyDigests(), c(lipid = "2d923db5e9fc4ff2d193b200014e73d9", metab = "c6fa5d5cc4881d65a71296422749962c"))
})

test_that("(KIT-129 item 3) the estimated count is the kernel's own k, per assay; with technical replicates it removes the units", {
    skip_if_not_installed("ruv")
    skip_if_not_installed("sva")
    fixture <- kit127Fixture()
    for (omic in c("lipid", "metab")) {
        mats <- lapply(fixture$assays, function(a) rbind(a$y, a$is))
        object <- if (omic == "lipid") kit127Lipid(mats, fixture$design) else kit127Metab(mats, fixture$design)
        run <- suppressWarnings(kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 11L)))
        expect_identical(run$parameters[c("ruv_k_rule", "ruv_k_seed")], list(ruv_k_rule = "num_sv_be", ruv_k_seed = 11L), info = omic)
        for (assay in names(mats)) {
            record <- run$factors[[assay]]
            expected <- tryCatch(estimateRuvFactorCount(mats[[assay]], fixture$design$group, seed = 11L), ruv_factor_count_unavailable = \(e) e)
            expect_identical(record$ruv_k_rule, "num_sv_be", info = omic)
            if (inherits(expected, "condition")) {
                expect_identical(record[c("status", "reason", "reason_class")], list(status = "refused",
                    reason = conditionMessage(expected), reason_class = class(expected)[[1L]]), info = omic)
            } else {
                expect_identical(record$k, expected$k, info = omic)
                expect_identical(record$ruv_k_estimate, expected$k, info = omic)
                expect_identical(record$k_estimate$components, expected$components, info = omic)
                expect_s3_class(plotRuvFactorCount(record$k_estimate), "ggplot")
            }
        }
    }
    object <- kit127Lipid(lapply(fixture$assays, `[[`, "y"), fixture$design)
    expect_error(runOmicsRuvFactors(object, "groups_all_features", 1L, "num_sv_be"), class = "omics_ruv_factors_invalid")
    expect_error(runOmicsRuvFactors(object, "groups_all_features", ruv_k_rule = "num_sv_be", ruv_k_seed = 1.5), class = "omics_ruv_factors_invalid")
    # technical replicates: the count removes the biological units (M28)
    for (omic in c("lipid", "metab")) {
        tech <- techRepRuvObject(omic)
        run <- suppressWarnings(kit127Quiet(runOmicsRuvFactors(tech, "groups_all_features", ruv_k_rule = "num_sv_be")))
        slots <- omicsRuvAssaySlots(tech)
        for (assay in names(methods::slot(tech, slots$data))) {
            values <- omicsRuvAssayValues(methods::slot(tech, slots$data)[[assay]], methods::slot(tech, slots$id),
                tech@design_matrix$Run, assay)
            expected <- estimateRuvFactorCount(values, tech@design_matrix$group, units = tech@design_matrix$tech_rep_group)
            expect_identical(run$factors[[assay]]$k_estimate$removed, "biological units", info = omic)
            expect_identical(run$factors[[assay]]$ruv_k_estimate, expected$k, info = omic)
        }
    }
})

# The replicate matrix of the biological units (arms 15 and 16 of simulator S4: `extractM`).
simUnitM <- function(design, unit_column, sample_column) {
    m <- ruv::replicate.matrix(stats::setNames(design[[unit_column]], design[[sample_column]]))
    rownames(m) <- design[[sample_column]]
    m
}

test_that("the technical-replicate design: the biological units as replicate sets, k limited by the units and the groups across them (KIT-129 integration)", {
    skip_if_not_installed("ruv")
    contrasts <- data.frame(contrasts = "groupB-groupA", friendly_names = "B_vs_A", stringsAsFactors = FALSE)
    for (omic in c("lipid", "metab")) {
        object <- techRepRuvObject(omic)
        slots <- omicsRuvAssaySlots(object)
        kernel <- function(...) suppressMessages(runOmicsRuvFactors(object, "groups_all_features", ...,
            ruv_grouping_variable = "tech_rep_group"))
        # 8 units in 2 groups: k = 2 and k = 5 are computed (the 8 units counted as 8 groups left none)
        run <- kernel(2L, "user")
        expect_identical(vapply(run$factors, `[[`, "", "status"), c(pos = "ok", neg = "ok"), info = omic)
        expect_identical(vapply(kernel(5L, "user")$factors, `[[`, 0L, "k"), c(pos = 5L, neg = 5L), info = omic)
        for (assay in c("pos", "neg")) {
            # W is arm 16's: RUV-III on the complete features, every one a control, the units as replicate sets
            values <- omicsRuvAssayValues(methods::slot(object, slots$data)[[assay]], methods::slot(object, slots$id),
                object@design_matrix$Run, assay)
            y <- t(values[rowSums(is.na(values)) == 0L, , drop = FALSE])
            expectKit127W(run$factors[[assay]]$W, simRuvW(y, simUnitM(object@design_matrix, "tech_rep_group", "Run"),
                rep(TRUE, ncol(y)), 2L))
            expect_identical(run$factors[[assay]]$ruv_grouping_variable, "tech_rep_group", info = omic)
        }
        # the DA runs with W in the model and the units as block
        da <- NULL
        utils::capture.output(da <- suppressWarnings(suppressMessages((if (identical(omic, "lipid")) runLipidsDA else runMetabolitesDA)(
            object, contrasts, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))))
        expect_null(da$failed_assays)
        expect_true(isTRUE(da$contrasts_results$pos$fit.eb$uses_duplicate_correlation), info = omic)
        # the estimated count removes the units (arm 29's model, with M30's permutations)
        skip_if_not_installed("sva")
        counted <- suppressWarnings(kernel(ruv_k_rule = "num_sv_be"))
        for (assay in c("pos", "neg")) {
            values <- omicsRuvAssayValues(methods::slot(object, slots$data)[[assay]], methods::slot(object, slots$id),
                object@design_matrix$Run, assay)
            expected <- estimateRuvFactorCount(values, object@design_matrix$tech_rep_group, units = object@design_matrix$tech_rep_group)
            expect_identical(counted$factors[[assay]]$k_estimate$removed, "biological units", info = omic)
            expect_identical(counted$factors[[assay]]$k, expected$k, info = omic)
        }
    }
})

test_that("first-pass null controls: the features outside the 25% most significant of the no-RUV test, W in the model (M7; S4 arms 5 and 23)", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    for (omic in c("lipid", "metab")) {
        object <- if (identical(omic, "lipid")) kit127Lipid(mats, fixture$design) else kit127Metab(mats, fixture$design)
        id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
        run <- kit127Quiet(runOmicsRuvFactors(object, "first_pass_nulls", K127_K, "user"))
        da <- kit127Quiet((if (identical(omic, "lipid")) runLipidsDA else runMetabolitesDA)(object, fixture$contrasts_tbl,
            "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
        long <- if (identical(omic, "lipid")) da$da_lipids_long else da$da_metabolites_long
        for (assay in names(mats)) {
            record <- run$factors[[assay]]
            controls <- simFirstPassControls(fixture, mats[[assay]])
            complete <- rownames(mats[[assay]])[rowSums(!is.finite(mats[[assay]])) == 0L]
            expect_gt(length(controls), 0L)
            expect_lt(length(controls), nrow(mats[[assay]]))
            expect_setequal(record$control_ids, intersect(controls, complete))
            expect_identical(record[c("status", "k", "ruv_design", "not_tested", "controls_left_out")],
                list(status = "ok", k = K127_K[[assay]], ruv_design = "first_pass_nulls", not_tested = character(),
                    controls_left_out = character()), info = omic)
            # W from the complete controls, the groups as replicate sets; the test reads every feature with W
            y <- t(mats[[assay]][complete, , drop = FALSE])
            w <- simRuvW(y, simGroupM(fixture$design), colnames(y) %in% controls, K127_K[[assay]], keep = fixture$design$Sample_ID)
            expectKit127W(record$W, w)
            expectKit127Table(long, id, assay, simTestWithW(fixture, mats[[assay]], w))
        }
        # the controls are a subset, so the separation rule is defined here (its k given, as for supplied controls)
        separation <- kit127Quiet(runOmicsRuvFactors(object, "first_pass_nulls", 1L, "sweep_bounded"))
        expect_identical(unname(vapply(separation$factors, `[[`, "", "status")), c("ok", "ok"), info = omic)
        expect_identical(unname(vapply(separation$factors, `[[`, "", "ruv_k_rule")), c("sweep_bounded", "sweep_bounded"))
        expect_error(runOmicsRuvFactors(object, "first_pass_nulls", 1L, "user", ruv_control_ids = lapply(mats, rownames)),
            class = "omics_ruv_factors_invalid")
    }
    # the estimated count
    skip_if_not_installed("sva")
    object <- kit127Lipid(mats, fixture$design)
    counted <- suppressWarnings(kit127Quiet(runOmicsRuvFactors(object, "first_pass_nulls", ruv_k_rule = "num_sv_be", ruv_k_seed = 5L)))
    for (assay in names(mats)) {
        expected <- tryCatch(estimateRuvFactorCount(mats[[assay]], fixture$design$group, seed = 5L), ruv_factor_count_unavailable = \(e) e)
        if (!inherits(expected, "condition")) expect_identical(counted$factors[[assay]]$ruv_k_estimate, expected$k)
    }
})

test_that("k must be below the number of controls observed in every sample (M32)", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    with_standards <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    standards <- lapply(fixture$assays, function(a) rownames(a$is))
    run <- function(k, rule) suppressMessages(runOmicsRuvFactors(with_standards, "supplied_controls", k, rule, ruv_control_ids = standards))
    # five standards per assay: k = 5 is not below them; typed, the run stops (M31); automatic, the assay is refused
    typed <- tryCatch(run(5L, "user"), error = function(e) e)
    expect_s3_class(typed, c("ruv_k_above_limit", "ruv_k_refused"))
    expect_match(conditionMessage(typed), "k = 5 is not below the 5 controls observed in every sample", fixed = TRUE)
    expect_warning(automatic <- run(5L, "sweep_bounded"), class = "omics_ruv_factors_refused")
    expect_identical(unname(vapply(automatic$factors, `[[`, "", "status")), c("refused", "refused"))
    expect_identical(unname(vapply(automatic$factors, `[[`, "", "reason_class")), c("ruv_k_above_limit", "ruv_k_above_limit"))
    # k = 4 is computed
    expect_identical(unname(vapply(run(4L, "user")$factors, `[[`, 0L, "k")), c(4L, 4L))
})

test_that("a feature with no value in any sample is not tested and is reported, the assay is analysed (KIT-127 found and left)", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    for (omic in c("lipid", "metab")) {
        mats <- lapply(fixture$assays, `[[`, "y")
        empty <- matrix(NA_real_, 1L, ncol(mats$pos), dimnames = list("pos_empty", colnames(mats$pos)))
        mats$pos <- rbind(mats$pos, empty)
        object <- if (identical(omic, "lipid")) kit127Lipid(mats, fixture$design) else kit127Metab(mats, fixture$design)
        run <- kit127Quiet(runOmicsRuvFactors(object, "groups_all_features", K127_K, "user"))
        warned <- NULL
        da <- withCallingHandlers(kit127Hush((if (identical(omic, "lipid")) runLipidsDA else runMetabolitesDA)(object,
            fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors), "omics_ruv_features_not_tested"),
            omics_ruv_features_not_tested = function(w) {
                warned <<- c(warned, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        long <- if (identical(omic, "lipid")) da$da_lipids_long else da$da_metabolites_long
        id <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
        # the assay is analysed (before: it failed with "prior.weights contain NA values")
        expect_null(da$failed_assays)
        expect_true("pos_empty" %in% da$ruv_features_not_tested$pos$not_estimable, info = omic)
        expect_false("pos_empty" %in% long[[id]], info = omic)
        expect_true(any(grepl("assay 'pos'.*1 with no value in any sample", warned)), info = omic)
        # the other features are tested as without it (arm 7's computation)
        expectKit127Table(long, id, "pos", simArm(fixture, "pos", K127_K$pos)$table)
    }
})
