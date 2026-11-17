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
simDaKernel <- function(m, design, formula, contrasts) {
    if (anyNA(m)) {
        x <- stats::model.matrix(stats::as.formula(formula), design)
        seen <- is.finite(m)
        pattern <- drop(seen %*% 2^(seq_len(ncol(m)) - 1))
        estimable <- vapply(split(seq_len(nrow(m)), pattern), function(i) qr(x[seen[i[[1L]], ], , drop = FALSE])$rank == ncol(x), logical(1L))
        m <- m[estimable[as.character(pattern)], , drop = FALSE]
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
    simDaKernel(y, design, paste("~ 0 + group +", paste(colnames(w), collapse = " + ")), fixture$contrasts)
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
    run <- kit127Quiet(runOmicsRuvFactors(object, ruv_design = "factors_all_features", ruv_number_k = K127_K, ruv_k_source = "user"))
    da <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE, ruv_factors = run$factors))
    for (assay in names(fixture$assays)) {
        reference <- simArm(fixture, assay, K127_K[[assay]])
        record <- run$factors[[assay]]
        expect_identical(record$status, "ok")
        expectKit127W(record$W, reference$W)
        expectKit127Table(da$da_lipids_long, "lipid_id", assay, reference$table)
        # the fixture has features the model with W cannot estimate: they are not tested, as in the simulator
        expect_lt(nrow(reference$table), nrow(fixture$assays[[assay]]$y))
        expect_identical(record$k, K127_K[[assay]])
        expect_identical(record$k_source, "user")
        expect_identical(record$control_count, sum(rowSums(is.na(fixture$assays[[assay]]$y)) == 0L))
        expect_identical(record$not_tested, character())
    }
    expect_identical(run$parameters$ruv_design, "factors_all_features")
    expect_identical(run$parameters$ruv_grouping_variable, "group")
})

test_that("arm 8 (the standards as controls, W in the model, the standards not tested) is the simulator's computation", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    object <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    ids <- lapply(fixture$assays, function(a) rownames(a$is))
    run <- kit127Quiet(runOmicsRuvFactors(object, ruv_design = "factors_supplied_controls", ruv_number_k = K127_K,
        ruv_k_source = "user", ruv_control_ids = ids))
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
    run_l <- kit127Quiet(runOmicsRuvFactors(lipid, "factors_all_features", K127_K, "user"))
    run_m <- kit127Quiet(runOmicsRuvFactors(metab, "factors_all_features", K127_K, "user"))
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

test_that("the corrected matrix of the new path keeps each feature's level; today's correction still shifts it", {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    samples <- fixture$design$Sample_ID
    for (make in list(kit127Lipid, kit127Metab)) {
        object <- make(mats, fixture$design)
        slot_name <- if (methods::is(object, "LipidomicsAssayData")) "lipid_data" else "metabolite_data"
        run <- kit127Quiet(runOmicsRuvFactors(object, "factors_all_features", K127_K, "user"))
        legacy <- kit127Quiet(ruvIII_C_Varying(object, ruv_grouping_variable = "group", ruv_number_k = K127_K,
            ctrl = lapply(mats, rownames)))
        for (assay in names(mats)) {
            before <- mats[[assay]]
            after <- kit127Values(methods::slot(run$corrected, slot_name)[[assay]], samples)
            today <- kit127Values(methods::slot(legacy, slot_name)[[assay]], samples)
            expect_identical(rownames(after), rownames(today))
            kept <- rownames(after)
            both <- is.finite(before[kept, ]) & is.finite(after)
            expect_identical(is.finite(after), is.finite(today))
            mean_in <- function(m) rowSums(replace(m, !both, 0)) / rowSums(both)
            expect_equal(mean_in(after), mean_in(before[kept, ]), tolerance = 1e-9)
            # the correction itself is today's: the two differ by a constant per feature
            expect_equal((after - mean_in(after))[both], (today - mean_in(today))[both], tolerance = 1e-9)
            expect_gt(max(abs(after - before[kept, ])[both]), 0.01)
            # today's path is not changed by this ticket: its level shift is still there
            expect_gt(max(abs(mean_in(today) - mean_in(before[kept, ]))), 0.01)
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
            expect_identical(got[[case$id]], rownames(direct))
            expect_identical(got$P.Value, direct$P.Value)
            expect_identical(got$logFC, direct$logFC)
            expect_identical(got$fdr_qvalue, direct$fdr_qvalue)
        }
    }
})

test_that("k has no default, 0 is no RUV, and an assay whose factors cannot be computed gets no RUV with a warning", {
    skip_if_not_installed("ruv")
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    object <- kit127Lipid(mats, fixture$design)
    expect_error(runOmicsRuvFactors(object, "factors_all_features"), class = "omics_ruv_factors_invalid")
    expect_error(runOmicsRuvFactors(object, "factors_all_features", 2L), class = "omics_ruv_factors_invalid")
    expect_error(runOmicsRuvFactors(object, ruv_number_k = 2L, ruv_k_source = "user"), class = "omics_ruv_factors_invalid")
    plain <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE))$da_lipids_long

    # k = 0: no factors, the corrected object is the input, the test is the plain test
    none <- kit127Quiet(runOmicsRuvFactors(object, "factors_all_features", 0L, "user"))
    expect_identical(vapply(none$factors, `[[`, "", "status"), c(pos = "none", neg = "none"))
    expect_null(none$factors$pos$W)
    expect_identical(none$corrected@lipid_data, object@lipid_data)
    expect_identical(kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = none$factors))$da_lipids_long, plain)

    # 12 samples in 2 replicate sets leave 9 factors at most; pos asks for 10 and is refused, neg is corrected
    expect_warning(refused <- suppressMessages(runOmicsRuvFactors(object, "factors_all_features", list(pos = 10L, neg = 1L), "user")),
        class = "omics_ruv_factors_refused")
    expect_identical(refused$factors$pos$status, "refused")
    expect_match(refused$factors$pos$reason, "residual")
    expect_identical(refused$factors$pos$k, 0L)
    expect_identical(refused$factors$pos$k_requested, 10L)
    expect_null(refused$factors$pos$W)
    expect_identical(refused$factors$neg$status, "ok")
    expect_identical(refused$corrected@lipid_data$pos, object@lipid_data$pos)
    expect_false(identical(refused$corrected@lipid_data$neg, object@lipid_data$neg))
    mixed <- kit127Quiet(runLipidsDA(object, fixture$contrasts_tbl, "~ 0 + group", 0.05, 0, TRUE, TRUE,
        ruv_factors = refused$factors))$da_lipids_long
    expect_identical(mixed[mixed$assay == "pos", ], plain[plain$assay == "pos", ])
    expectKit127Table(mixed, "lipid_id", "neg", simArm(fixture, "neg", 1L)$table)

    # more factors than controls observed in every sample
    with_standards <- kit127Lipid(lapply(fixture$assays, function(a) rbind(a$y, a$is)), fixture$design)
    expect_warning(few <- suppressMessages(runOmicsRuvFactors(with_standards, "factors_supplied_controls", 6L, "user",
        ruv_control_ids = lapply(fixture$assays, function(a) rownames(a$is)))), class = "omics_ruv_factors_refused")
    expect_identical(unname(vapply(few$factors, `[[`, "", "status")), c("refused", "refused"))
    expect_match(few$factors$pos$reason, "controls observed in every sample")

    # a k that could not be estimated is a refusal, not a guess
    expect_warning(unknown <- suppressMessages(runOmicsRuvFactors(object, "factors_all_features", list(pos = NA_integer_, neg = 1L),
        "sweep_bounded")), class = "omics_ruv_factors_refused")
    expect_identical(unknown$factors$pos$status, "refused")
    expect_identical(unknown$factors$neg$k_source, "sweep_bounded")
})

test_that("wrong arguments stop", {
    fixture <- kit127Fixture()
    mats <- lapply(fixture$assays, `[[`, "y")
    object <- kit127Lipid(mats, fixture$design)
    ids <- lapply(mats, function(m) rownames(m)[1:5])
    invalid <- function(...) expect_error(kit127Quiet(runOmicsRuvFactors(object, ...)), class = "omics_ruv_factors_invalid")
    invalid("legacy", 1L, "user")
    invalid("factors_all_features", 1.5, "user")
    invalid("factors_all_features", -1L, "user")
    invalid("factors_all_features", list(pos = 1L), "user")
    invalid("factors_all_features", 1L, "guess")
    invalid("factors_all_features", 1L, "user", ruv_control_ids = ids)
    invalid("factors_supplied_controls", 1L, "user")
    invalid("factors_supplied_controls", 1L, "user", ruv_control_ids = list(pos = ids$pos))
    invalid("factors_supplied_controls", 1L, "user", ruv_control_ids = list(pos = "absent", neg = ids$neg))
    invalid("factors_all_features", 1L, "user", ruv_grouping_variable = "absent")
    expect_error(runOmicsRuvFactors(fixture$design, "factors_all_features", 1L, "user"), class = "omics_ruv_factors_invalid")
    run <- kit127Quiet(runOmicsRuvFactors(object, "factors_all_features", 1L, "user"))
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
    expect_identical(bounded$source, "sweep_bounded")
    failed <- boundOmicsRuvSweepK(object, list(pos = list(success = FALSE, error = "no valid percentage"), neg = sweep$neg))
    expect_identical(failed$k, list(pos = NA_integer_, neg = 1L))
    expect_error(boundOmicsRuvSweepK(object, sweep["pos"]), class = "omics_ruv_factors_invalid")

    skip_if_not_installed("ruv")
    params <- list(ruv_grouping_variable = "group", percentage_min = 10, percentage_max = 12, separation_metric = "max_difference",
        k_penalty_weight = 0.5, adaptive_k_penalty = TRUE)
    real <- boundOmicsRuvSweepK(object, kit127Quiet(runLipidPerAssayRuvOptimization(object, "automatic", params)))
    known <- !is.na(unlist(real$k))
    expect_true(all(unlist(real$k)[known] <= 2L))
    run <- kit127Quiet(runOmicsRuvFactors(object, "factors_all_features", real$k, real$source))
    expect_identical(unname(vapply(run$factors, `[[`, "", "k_source")), c("sweep_bounded", "sweep_bounded"))
    expect_true(all(vapply(run$factors, `[[`, "", "status") %in% c("ok", "refused")))
})
