# KIT-193 group R (ticket sections 7.2, 7.4 "Paired R tests for group R", 5.5): the RUV algebra the Lean files prove
# (Spec/RuvAlgebra.lean, RuvDownstream.lean, RuvIIIC.lean, RuvControlRules.lean), checked on the shipped functions and
# the locked libraries they call. Fixtures: integer entries, 12 runs in 3 groups of 4, 30 features, 10 controls
# (features 1 to 10, no group effect), k = 1, 2; a variant with 3 missing values (non-control features 16, 18, 20,
# runs 2, 7, 11); for R-K1 the same with 14 controls. Each check to relative 1e-10 of the largest absolute value
# unless stated. A tolerance check is evidence that the code implements the model on these fixtures, not a proof
# about floating point (ticket 9.1). Kinds and planted defects: models.md section 2. Not here, with the reason:
# R-I1 (KIT-184's writer, not merged), the lipid method's cases (L1) and R-N1's sweep-helper cases (KIT-190, not
# merged), R-J3's clean-up counts (the pattern is checked; the counts need the clean-up step's own fixture).

kit193RuvData <- function(controls = 10L, missing = FALSE) {
    runs <- sprintf("R%02d", 1:12)
    group <- rep(c("A", "B", "C"), each = 4L)
    ids <- sprintf("P%02d", 1:30)
    Y <- withr::with_seed(20261193L, {
        level <- sample(15:30, 30L, replace = TRUE)
        u1 <- sample(-3:3, 12L, replace = TRUE)
        u2 <- sample(-2:2, 12L, replace = TRUE)
        b1 <- sample(-2:2, 30L, replace = TRUE)
        b2 <- sample(-2:2, 30L, replace = TRUE)
        noise <- matrix(sample(-2:2, 360L, replace = TRUE), 12L, 30L)
        effect <- matrix(0, 12L, 30L)
        effect[group == "B", 21:30] <- 3
        effect[group == "C", 26:30] <- -2
        outer(rep(1, 12L), level) + outer(u1, b1) + outer(u2, b2) + noise + effect
    })
    dimnames(Y) <- list(runs, ids)
    if (missing) Y[cbind(c(2L, 7L, 11L), c(16L, 18L, 20L))] <- NA
    M <- outer(group, unique(group), "==") * 1
    dimnames(M) <- list(runs, unique(group))
    X <- M
    colnames(X) <- paste0("group", unique(group))
    list(Y = Y, runs = runs, group = group, ids = ids, ctl = seq_along(ids) <= controls, M = M, X = X)
}

kit193RuvObject <- function(d, Y = d$Y) {
    table <- data.frame(Protein.Ids = colnames(Y), t(Y), check.names = FALSE, row.names = NULL)
    ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids",
        protein_id_table = table["Protein.Ids"], design_matrix = data.frame(Run = rownames(Y), group = d$group),
        sample_id = "Run", group_id = "group", args = list())
}

#' One expectation: max |a - b| <= rel * max |b|
kit193Close <- function(a, b, rel = 1e-10, info = NULL) {
    a <- unname(as.matrix(a))
    b <- unname(as.matrix(b))
    expect_lte(max(abs(a - b)), rel * max(abs(b), 1), label = paste("max |difference|", info))
}

#' RUV-III-C as the package calls it (complete Y or not), runs x features
kit193RuvIIIC <- function(d, k, Y = d$Y, ctl = d$ctl) {
    out <- RUVIIIC::RUVIII_C_Varying(k = k, Y = Y, M = d$M, toCorrect = colnames(Y), potentialControls = d$ids[ctl],
        progress = FALSE)
    out[rownames(Y), colnames(Y), drop = FALSE]
}

#' The corrected table of the package's ruvIII_C_Varying (restored levels), features x runs
kit193RuvMethod <- function(d, k, Y = d$Y) {
    out <- suppressMessages(ruvIII_C_Varying(kit193RuvObject(d, Y), ruv_grouping_variable = "group", ruv_number_k = k,
        ctrl = stats::setNames(d$ctl, d$ids)))
    m <- as.matrix(out@protein_quant_table[, d$runs])
    rownames(m) <- out@protein_quant_table$Protein.Ids
    m[d$ids, , drop = FALSE]
}

#' A = ac' (ac ac')^-1 alpha (controls x features), from computeRuvCovariates()$alpha
kit193RuvA <- function(alpha, ctl) {
    ac <- alpha[, ctl, drop = FALSE]
    t(ac) %*% solve(ac %*% t(ac)) %*% alpha
}

kit193RowShift <- function(Y, c) Y + matrix(c, nrow(Y), ncol(Y), byrow = TRUE)

test_that("KIT-193 R-A1, R-A2 pair: ruv::RUVIII's Y - newY = W alpha, rank <= k, alpha = U'Y, U'newY = 0, the controls' form", {
    expect_identical(format(utils::packageVersion("ruv")), "0.9.7.1")
    d <- kit193RuvData()
    Y0 <- d$Y - d$M %*% solve(t(d$M) %*% d$M) %*% t(d$M) %*% d$Y
    for (k in 1:2) {
        r <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = k, return.info = TRUE)
        alpha <- r$fullalpha[seq_len(k), , drop = FALSE]
        ac <- alpha[, d$ctl, drop = FALSE]
        W <- d$Y[, d$ctl] %*% t(ac) %*% solve(ac %*% t(ac))
        kit193Close(d$Y - r$newY, W %*% alpha, info = paste("A1 k", k))
        expect_lte(qr(d$Y - r$newY, tol = 1e-10)$rank, k)
        U <- svd(Y0 %*% t(Y0))$u[, seq_len(k), drop = FALSE]
        projected <- t(U) %*% d$Y
        sign <- sign(rowSums(projected * alpha))
        kit193Close(sign * projected, alpha, info = paste("alpha = U'Y k", k))
        expect_lte(max(abs(t(U) %*% r$newY)), 1e-10 * max(abs(d$Y)))
        kit193Close(r$newY[, d$ctl], d$Y[, d$ctl] %*% (diag(sum(d$ctl)) - t(ac) %*% solve(ac %*% t(ac)) %*% ac),
            info = paste("A2 controls k", k))
    }
})

test_that("KIT-193 R-A1 pair, planted defect: ruv::RUVIII returning newY = Y fails Y - newY = W alpha", {
    d <- kit193RuvData()
    original <- ruv::RUVIII
    testthat::local_mocked_bindings(RUVIII = function(Y, M, ctl, k = NULL, ...) {
        out <- original(Y = Y, M = M, ctl = ctl, k = k, return.info = TRUE)
        out$newY <- Y
        out
    }, .package = "ruv")
    r <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = 1L, return.info = TRUE)
    alpha <- r$fullalpha[1L, , drop = FALSE]
    ac <- alpha[, d$ctl, drop = FALSE]
    W <- d$Y[, d$ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    testthat::expect_failure(kit193Close(d$Y - r$newY, W %*% alpha))
})

test_that("KIT-193 R-A3 pair: k = 0 writes no factor (computeOmicsRuvAssayFactors status none; computeProtRuvModelFactors no W)", {
    d <- kit193RuvData()
    omics <- computeOmicsRuvAssayFactors(t(d$Y), d$M, k = 0L)
    expect_identical(omics$status, "none")
    expect_null(omics$W)
    # The user rule refuses k = 0 (validateProtRuvModelK); the k = 0 branch is reached through the estimated count.
    testthat::local_mocked_bindings(estimateRuvFactorCount = function(...) list(k = 0L))
    prot <- suppressMessages(computeProtRuvModelFactors(kit193RuvObject(d), ruv_k_rule = "num_sv_be"))
    expect_identical(prot$k, 0L)
    expect_null(prot$W)
})

test_that("KIT-193 R-A3 pair, planted defect: a k = 0 that returns a zero W with status ok fails", {
    d <- kit193RuvData()
    testthat::local_mocked_bindings(computeOmicsRuvAssayFactors = function(values, replicate_matrix, k, ...) {
        list(status = "ok", W = matrix(0, ncol(values), 1L), k = 0L)
    })
    omics <- computeOmicsRuvAssayFactors(t(d$Y), d$M, k = 0L)
    testthat::expect_failure(expect_identical(omics$status, "none"))
})

test_that("KIT-193 R-B1, R-B2, R-E1 pair: level shifts move RUV-III-C's output by c_j - v'alpha_j; restored, by exactly c", {
    expect_identical(format(utils::packageVersion("RUVIIIC")), "1.0.19")
    d <- kit193RuvData()
    shift <- withr::with_seed(5L, sample(-8:8, 30L, replace = TRUE)) / 4
    for (k in 1:2) {
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        A <- kit193RuvA(cov$alpha, d$ctl)
        r0 <- kit193RuvIIIC(d, k)
        r1 <- kit193RuvIIIC(d, k, kit193RowShift(d$Y, shift))
        r2 <- kit193RuvIIIC(d, k, d$Y + 2)
        amount1 <- shift - colSums(shift[d$ctl] * A)
        amount2 <- 2 * (1 - colSums(A))
        kit193Close(r1 - r0, matrix(amount1, 12L, 30L, byrow = TRUE), info = paste("B1 k", k))
        kit193Close(r2 - r0, matrix(amount2, 12L, 30L, byrow = TRUE), info = paste("E1 k", k))
        expect_gt(diff(range(amount2)), 1e-6)
        restore <- function(r, Y) t(restoreRuvFeatureLevels(t(r), t(Y)))
        kit193Close(restore(r1, kit193RowShift(d$Y, shift)) - restore(r0, d$Y), matrix(shift, 12L, 30L, byrow = TRUE),
            info = paste("B2 k", k))
        kit193Close(restore(r2, d$Y + 2) - restore(r0, d$Y), matrix(2, 12L, 30L), info = paste("E1 restored k", k))
        cov1 <- computeRuvCovariates(kit193RuvObject(d, kit193RowShift(d$Y, shift)), k, d$ctl)
        kit193Close(cov1$adjusted - cov$adjusted, matrix(shift, 30L, 12L), info = paste("adjusted shift k", k))
        kit193Close(cov1$W, cov$W, info = paste("W under a shift k", k))
    }
})

test_that("KIT-193 R-B2 pair, planted defect: restoreRuvFeatureLevels returning its input fails the exact shift", {
    d <- kit193RuvData()
    shift <- withr::with_seed(5L, sample(-8:8, 30L, replace = TRUE)) / 4
    testthat::local_mocked_bindings(restoreRuvFeatureLevels = function(corrected, original) corrected)
    restore <- function(r, Y) t(restoreRuvFeatureLevels(t(r), t(Y)))
    moved <- restore(kit193RuvIIIC(d, 1L, kit193RowShift(d$Y, shift)), kit193RowShift(d$Y, shift)) -
        restore(kit193RuvIIIC(d, 1L), d$Y)
    testthat::expect_failure(kit193Close(moved, matrix(shift, 12L, 30L, byrow = TRUE)))
})

test_that("KIT-193 R-C1 pair: restoreRuvFeatureLevels(newY, Y) is computeRuvCovariates' adjusted; the second centring is the identity", {
    d <- kit193RuvData()
    for (k in 1:2) {
        r <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = k)
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        kit193Close(restoreRuvFeatureLevels(t(r), t(d$Y)), cov$adjusted, info = paste("C1 k", k))
        omics <- computeOmicsRuvAssayFactors(t(d$Y), d$M, k = k, control_ids = d$ids[d$ctl])
        expect_identical(omics$status, "ok")
        kit193Close(omics$W, cov$W, info = paste("omics W k", k))
        expect_lte(max(abs(colSums(cov$W))), 1e-10 * max(abs(cov$W)))
    }
})

test_that("KIT-193 R-C1 pair, planted defect (ticket): restoreRuvFeatureLevels returning its input fails", {
    d <- kit193RuvData()
    r <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = 1L)
    cov <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)
    testthat::local_mocked_bindings(restoreRuvFeatureLevels = function(corrected, original) corrected)
    testthat::expect_failure(kit193Close(restoreRuvFeatureLevels(t(r), t(d$Y)), cov$adjusted))
})

test_that("KIT-193 R-D1, R-D2 pair: a rotated basis gives the same newY; 2Y gives the same W and twice the adjusted values", {
    d <- kit193RuvData()
    for (k in 1:2) {
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        cov2 <- computeRuvCovariates(kit193RuvObject(d, 2 * d$Y), k, d$ctl)
        kit193Close(cov2$W, cov$W, info = paste("D2 W k", k))
        kit193Close(cov2$adjusted, 2 * cov$adjusted, info = paste("D2 adjusted k", k))
        r <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = k, return.info = TRUE)
        Q <- if (k == 1L) matrix(-3) else matrix(c(2, 1, 1, 1), 2L)
        rotated <- ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = k, fullalpha = Q %*% r$fullalpha[seq_len(k), , drop = FALSE])
        kit193Close(rotated, r$newY, info = paste("D1 k", k))
    }
})

test_that("KIT-193 R-D2 pair, planted defect: a computeRuvCovariates that doubles W on the scaled input fails 'W unchanged'", {
    d <- kit193RuvData()
    cov <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)
    original <- computeRuvCovariates
    testthat::local_mocked_bindings(computeRuvCovariates = function(...) {
        out <- original(...)
        out$W <- 2 * out$W
        out
    })
    cov2 <- computeRuvCovariates(kit193RuvObject(d, 2 * d$Y), 1L, d$ctl)
    testthat::expect_failure(kit193Close(cov2$W, cov$W))
})

#' Residual sum of squares per feature of an MArrayLM fit
kit193Rss <- function(fit) fit$sigma^2 * fit$df.residual

test_that("KIT-193 R-F1, R-F2, R-F3, R-G2 pair: limma on [groups, W]: the same group coefficients, residuals and df for y and y*", {
    d <- kit193RuvData()
    y <- t(d$Y)
    weights <- matrix(rep(c(1, 2, 0.5, 4), 3L), 30L, 12L, byrow = TRUE)
    for (k in 1:2) {
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        prot <- computeProtRuvModelFactors(kit193RuvObject(d), ruv_k = k, ruv_k_rule = "user")
        # (W, y*): computeRuvCovariates' W and adjusted values (P5, P6); computeProtRuvModelFactors' W and the display
        # matrix of the same design, RUV-III-C with every protein a control, levels restored (P3)
        all_ctl <- d
        all_ctl$ctl <- rep(TRUE, 30L)
        pairs <- list(covariates = list(W = cov$W, ystar = cov$adjusted),
            model_factors = list(W = prot$W[d$runs, , drop = FALSE], ystar = kit193RuvMethod(all_ctl, k)))
        for (name in names(pairs)) {
            info <- paste(name, "k", k)
            full <- cbind(d$X, pairs[[name]]$W)
            ystar <- pairs[[name]]$ystar
            f_y <- limma::lmFit(y, full)
            f_star <- limma::lmFit(ystar, full)
            kit193Close(f_star$coefficients[, 1:3], f_y$coefficients[, 1:3], info = info)
            kit193Close(f_star$sigma, f_y$sigma, info = info)
            expect_identical(f_star$df.residual, f_y$df.residual)
            expect_equal(unique(f_y$df.residual), 12 - 3 - k)
            kit193Close(stats::residuals(f_star, ystar), stats::residuals(f_y, y), info = info)
            w_y <- limma::lmFit(y, full, weights = weights)
            w_star <- limma::lmFit(ystar, full, weights = weights)
            kit193Close(w_star$coefficients[, 1:3], w_y$coefficients[, 1:3], info = paste(info, "weighted"))
            kit193Close(w_star$sigma, w_y$sigma, info = paste(info, "weighted"))
            # R-G2: RSS(y* on groups) >= RSS(y on [groups, W]) = RSS(y* on [groups, W])
            f_groups <- limma::lmFit(ystar, d$X)
            expect_true(all(kit193Rss(f_groups) >= kit193Rss(f_y) * (1 - 1e-10)), info = info)
            kit193Close(kit193Rss(f_star), kit193Rss(f_y), info = info)
        }
    }
})

test_that("KIT-193 R-G1 pair (an illustration): with C orthogonal to the groups, the naive t is sqrt((n - g)/(n - g - q)) times the full model's", {
    d <- kit193RuvData()
    y <- t(d$Y)
    z <- c(3, -1, 2, -4, 1, 0, -2, 1, 5, -3, -1, 2)
    C <- matrix(z - d$X %*% solve(crossprod(d$X), crossprod(d$X, z)), ncol = 1L, dimnames = list(d$runs, "C1"))
    expect_lte(max(abs(crossprod(d$X, C))), 1e-12)
    ystar <- y - t(C %*% solve(crossprod(C), t(C) %*% t(y)))
    contrast <- matrix(c(-1, 1, 0), ncol = 1L, dimnames = list(colnames(d$X), "BvA"))
    t_of <- function(fit) as.numeric(fit$coefficients / (fit$stdev.unscaled * fit$sigma))
    naive <- t_of(limma::contrasts.fit(limma::lmFit(ystar, d$X), contrast))
    full <- t_of(limma::contrasts.fit(limma::lmFit(y, cbind(d$X, C)), rbind(contrast, C1 = 0)))
    kit193Close(naive, sqrt(9 / 8) * full)
})

test_that("KIT-193 R-F1 pair, planted defect: adjusted values Y - 2 W alpha (W from the uncentred controls) change the group coefficients", {
    d <- kit193RuvData()
    y <- t(d$Y)
    cov <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)
    original <- computeRuvCovariates
    testthat::local_mocked_bindings(computeRuvCovariates = function(theObject, k, controls, ...) {
        out <- original(theObject, k, controls, ...)
        Y <- t(out$adjusted) + out$W %*% out$alpha
        ac <- out$alpha[, controls, drop = FALSE]
        W <- Y[, controls, drop = FALSE] %*% t(ac) %*% solve(ac %*% t(ac))
        out$adjusted <- t(Y - 2 * W %*% out$alpha)
        out
    })
    planted <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)
    full <- cbind(d$X, cov$W)
    testthat::expect_failure(kit193Close(limma::lmFit(planted$adjusted, full)$coefficients[, 1:3],
        limma::lmFit(y, full)$coefficients[, 1:3]))
})

#' Rank of the residuals of a features x runs matrix on the group design (runs x features), qr at tol 1e-8
kit193ResidRank <- function(m, X) {
    r <- t(m) - X %*% solve(crossprod(X), crossprod(X, t(m)))
    qr(r, tol = 1e-8)$rank
}

test_that("KIT-193 R-G0 pair: the corrected values' residuals on the groups have rank n - g - k (this fixture; at most, proved)", {
    # The fixture is in RUV-III-C's general branch (k below the 9 residual dimensions, 10 complete controls) with the
    # default formula; "equal" is a property of the fixture, the theorem gives "at most" (proofs audit P05 item 4).
    d <- kit193RuvData()
    expect_identical(kit193ResidRank(t(d$Y), d$X), 9L)
    for (k in 1:2) {
        expect_identical(kit193ResidRank(kit193RuvMethod(d, k), d$X), 9L - k, info = paste("ruvIII_C_Varying k", k))
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        expect_identical(kit193ResidRank(cov$adjusted, d$X), 9L - k, info = paste("adjusted k", k))
    }
})

test_that("KIT-193 R-G0 pair, planted defect: the uncorrected values in place of the corrected fail rank n - g - k", {
    d <- kit193RuvData()
    testthat::expect_failure(expect_identical(kit193ResidRank(t(d$Y), d$X), 9L - 1L))
})

#' The exact contrast and its unscaled SD per feature on its observed runs, with run weights w
kit193ExactContrasts <- function(y, design, contrasts, w = rep(1, ncol(y))) {
    rows <- lapply(seq_len(nrow(y)), function(i) {
        o <- is.finite(y[i, ])
        Xo <- design[o, , drop = FALSE]
        V <- solve(t(Xo) %*% (w[o] * Xo))
        b <- V %*% t(Xo) %*% (w[o] * y[i, o])
        list(coef = as.numeric(t(contrasts) %*% b), su = sqrt(diag(t(contrasts) %*% V %*% contrasts)))
    })
    list(coef = do.call(rbind, lapply(rows, `[[`, "coef")), su = do.call(rbind, lapply(rows, `[[`, "su")))
}

test_that("KIT-193 R-H1 pair: ruvExactContrastsFit gives c'beta and sqrt(c'(X_o'X_o)^-1 c) on features with missing runs, weighted too", {
    d <- kit193RuvData()
    dm <- kit193RuvData(missing = TRUE)
    W <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)$W
    design <- cbind(d$X, W)
    colnames(design) <- c("A", "B", "C", "W1")
    contrasts <- cbind(BvA = c(-1, 1, 0, 0), CvA = c(-1, 0, 1, 0))
    rownames(contrasts) <- colnames(design)
    y <- t(dm$Y)
    expect_identical(sum(!is.finite(y)), 3L)
    fit <- limma::lmFit(y, design)
    exact <- ruvExactContrastsFit(fit, contrasts, design, function(x) limma::lmFit(y, x))
    expected <- kit193ExactContrasts(y, design, contrasts)
    kit193Close(exact$coefficients, expected$coef, info = "coefficients")
    kit193Close(exact$stdev.unscaled, expected$su, info = "stdev.unscaled")
    w <- rep(c(1, 2, 0.5, 4), 3L)
    weights <- matrix(w, nrow(y), ncol(y), byrow = TRUE)
    fit_w <- limma::lmFit(y, design, weights = weights)
    exact_w <- ruvExactContrastsFit(fit_w, contrasts, design, function(x) limma::lmFit(y, x, weights = weights))
    expected_w <- kit193ExactContrasts(y, design, contrasts, w)
    kit193Close(exact_w$coefficients, expected_w$coef, info = "weighted coefficients")
    kit193Close(exact_w$stdev.unscaled, expected_w$su, info = "weighted stdev.unscaled")
})

test_that("KIT-193 R-H1 pair, planted defect: contrasts.fit's approximation in place of the refit fails on the missing runs", {
    d <- kit193RuvData()
    dm <- kit193RuvData(missing = TRUE)
    W <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)$W
    design <- cbind(d$X, W)
    colnames(design) <- c("A", "B", "C", "W1")
    contrasts <- cbind(BvA = c(-1, 1, 0, 0), CvA = c(-1, 0, 1, 0))
    rownames(contrasts) <- colnames(design)
    y <- t(dm$Y)
    testthat::local_mocked_bindings(ruvExactContrastsFit = function(fit, contrasts, design, refit) {
        limma::contrasts.fit(fit, contrasts)
    })
    approx <- ruvExactContrastsFit(limma::lmFit(y, design), contrasts, design, function(x) limma::lmFit(y, x))
    testthat::expect_failure(kit193Close(approx$stdev.unscaled, kit193ExactContrasts(y, design, contrasts)$su))
})

test_that("KIT-193 R-J1, R-J2, R-J3 pair: RUV-III-C with missing values: one W for complete features, upstream's zero-filled equality, the pattern", {
    dm <- kit193RuvData(missing = TRUE)
    Y0 <- dm$Y
    Y0[is.na(Y0)] <- 0
    for (k in 1:2) {
        r <- kit193RuvIIIC(dm, k)
        complete <- colSums(is.na(dm$Y)) == 0L
        s <- svd((r - dm$Y)[, complete])$d
        expect_lte(s[[k + 1L]], 1e-10 * s[[1L]])
        expect_identical(unname(is.finite(r)), unname(is.finite(dm$Y)))
        for (j in seq_along(dm$ids)) {
            rows <- is.finite(dm$Y[, j])
            upstream <- ruv::RUVIII(Y = Y0[rows, , drop = FALSE], M = dm$M[rows, , drop = FALSE], ctl = dm$ctl, k = k)
            kit193Close(r[rows, j], upstream[, j], rel = 1e-8, info = paste("J2 feature", j, "k", k))
        }
    }
})

test_that("KIT-193 R-J2 pair, planted defect: mean-filled in place of zero-filled fails upstream's equality", {
    dm <- kit193RuvData(missing = TRUE)
    filled <- dm$Y
    for (j in which(colSums(is.na(filled)) > 0L)) filled[is.na(filled[, j]), j] <- mean(filled[, j], na.rm = TRUE)
    r <- kit193RuvIIIC(dm, 1L)
    upstream <- ruv::RUVIII(Y = filled, M = dm$M, ctl = dm$ctl, k = 1L)
    testthat::expect_failure(kit193Close(r[, 1L], upstream[, 1L], rel = 1e-8))
})

test_that("KIT-193 R-K1 pair: k = 9 = n - sets with 14 controls: RUV-III-C equals ruv's k = 9 and its closed form; constant within sets", {
    d <- kit193RuvData(controls = 14L)
    rc <- kit193RuvIIIC(d, 9L)
    kit193Close(rc, ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = 9L), rel = 1e-8, info = "k = 9")
    kit193Close(rc, ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = NULL), rel = 1e-8, info = "closed form")
    spread <- vapply(unique(d$group), function(g) max(apply(rc[d$group == g, , drop = FALSE], 2L, function(x) diff(range(x)))),
        numeric(1L))
    expect_lte(max(spread), 1e-8)
})

test_that("KIT-193 R-K1 pair, planted defect: ruv's newY at k = 8 in place of the closed form fails", {
    d <- kit193RuvData(controls = 14L)
    rc <- kit193RuvIIIC(d, 9L)
    testthat::expect_failure(kit193Close(rc, ruv::RUVIII(Y = d$Y, M = d$M, ctl = d$ctl, k = 8L), rel = 1e-8))
})

test_that("KIT-193 R-L1, R-L2 pair: a group contrast of ruvIII_C_Varying's values changes by e_j = -d_C'A a_j; |e_j| within the bound", {
    d <- kit193RuvData()
    y <- t(d$Y)
    lambda <- rep(c(-1 / 4, 1 / 4, 0), each = 4L)
    for (k in 1:2) {
        cov <- computeRuvCovariates(kit193RuvObject(d), k, d$ctl)
        A <- kit193RuvA(cov$alpha, d$ctl)
        dC <- as.numeric(t(d$Y[, d$ctl]) %*% lambda)
        e <- -as.numeric(dC %*% A)
        change <- as.numeric((kit193RuvMethod(d, k) - y) %*% lambda)
        kit193Close(change, e, info = paste("L1 k", k))
        ac <- cov$alpha[, d$ctl, drop = FALSE]
        kappa <- min(eigen(ac %*% t(ac), symmetric = TRUE, only.values = TRUE)$values)
        bound <- sqrt(sum(dC^2)) * sqrt(colSums(cov$alpha^2)) / sqrt(kappa)
        expect_true(all(abs(e) <= bound * (1 + 1e-10)), info = paste("L2 k", k))
    }
})

test_that("KIT-193 R-L1 pair, planted defect: a contrast whose weights do not sum to 0 fails on the restored values", {
    d <- kit193RuvData()
    lambda <- rep(c(0, 1 / 4, 0), each = 4L)
    cov <- computeRuvCovariates(kit193RuvObject(d), 1L, d$ctl)
    e <- -as.numeric(as.numeric(t(d$Y[, d$ctl]) %*% lambda) %*% kit193RuvA(cov$alpha, d$ctl))
    change <- as.numeric((kit193RuvMethod(d, 1L) - t(d$Y)) %*% lambda)
    testthat::expect_failure(kit193Close(change, e))
})

#' edgeR's quasi-likelihood fit and F test of B - A, as the transcriptomics DA calls it
kit193EdgeR <- function(counts, design, contrast) {
    dge <- edgeR::estimateDisp(edgeR::DGEList(counts), design)
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    list(fitted = fit$fitted.values, F = edgeR::glmQLFTest(fit, contrast = contrast)$table$F)
}

kit193Counts <- function(d) {
    counts <- matrix(withr::with_seed(7L, stats::rnbinom(360L, mu = 200, size = 10)), 30L, 12L,
        dimnames = list(d$ids, d$runs))
    counts[21:30, d$group == "B"] <- 2L * counts[21:30, d$group == "B"]
    counts
}

test_that("KIT-193 R-M1 pair: edgeR's fit with [X, W] and [X, W Q] (Q invertible): equal fitted values and F statistics", {
    d <- kit193RuvData()
    W <- computeRuvCovariates(kit193RuvObject(d), 2L, d$ctl)$W
    counts <- kit193Counts(d)
    Q <- matrix(c(1, 0.5, 0, 1), 2L)
    a <- kit193EdgeR(counts, cbind(d$X, W), c(-1, 1, 0, 0, 0))
    b <- kit193EdgeR(counts, cbind(d$X, W %*% Q), c(-1, 1, 0, 0, 0))
    kit193Close(b$fitted, a$fitted, rel = 1e-8, info = "fitted")
    kit193Close(b$F, a$F, rel = 1e-8, info = "F")
})

test_that("KIT-193 R-M1 pair, planted defect: W q of rank 1 (one column) in place of W Q changes the fitted values", {
    # models.md proposes a non-invertible Q; [X, W Q] with Q of rank 1 and two columns is not of full rank, which
    # edgeR refuses, so the planted span is W q, its one column (the same rank-1 map).
    d <- kit193RuvData()
    W <- computeRuvCovariates(kit193RuvObject(d), 2L, d$ctl)$W
    counts <- kit193Counts(d)
    a <- kit193EdgeR(counts, cbind(d$X, W), c(-1, 1, 0, 0, 0))
    b <- kit193EdgeR(counts, cbind(d$X, W %*% c(1, 0.5)), c(-1, 1, 0, 0))
    testthat::expect_failure(kit193Close(b$fitted, a$fitted, rel = 1e-8))
})

kit193FirstPassTable <- function() {
    table <- utils::read.delim(testthat::test_path("..", "testdata", "kit193", "first-pass.tsv"), sep = "\t", quote = "",
        colClasses = "character", na.strings = character())
    table$p <- ifelse(table$p8 == "NA", NA_real_, suppressWarnings(as.numeric(table$p8)) / 8)
    table
}

#' The three shipped first-pass rules driven with a case's p-values (the p-value computations replaced by the case's)
kit193FirstPassShipped <- function(table, state) {
    out <- lapply(split(table, factor(table$case, levels = unique(table$case))), function(rows) {
        n <- nrow(rows)
        ids <- sprintf("f%02d", seq_len(n))
        state$pv <- rows$p
        values <- matrix(withr::with_seed(n, stats::rnorm(n * 6L)), n, 6L, dimnames = list(ids, sprintf("s%d", 1:6)))
        groups <- rep(c("a", "b"), each = 3L)
        prot <- protRuvFirstPassControls(values, groups)
        omics <- omicsRuvFirstPassControls(values, data.frame(Run = colnames(values), group = groups), "Run", "group",
            NA_character_, fitFn = function(values, ...) list(ids = rownames(values)),
            applyFn = function(fit, contrasts, ...) list(results = list(data.frame(P.Value = state$pv, row.names = fit$ids))))
        trx <- trxRuvDesignControls(values, data.frame(Run = colnames(values), group = groups), "first_pass_nulls")
        data.frame(prot = as.integer(prot), omics = as.integer(ids %in% omics), trx = as.integer(ids %in% trx$control_genes))
    })
    do.call(rbind, out)
}

test_that("KIT-193 R-N1 pair (first-pass rules): protRuvFirstPassControls, omicsRuvFirstPassControls and trxRuvDesignControls equal the Lean model", {
    table <- kit193FirstPassTable()
    state <- new.env()
    testthat::local_mocked_bindings(eBayes = function(fit, ...) fit,
        topTable = function(fit, ...) data.frame(P.Value = state$pv), .package = "limma")
    testthat::local_mocked_bindings(trxFilteredDge = function(...) NULL, trxPrelimPValues = function(...) state$pv)
    got <- kit193FirstPassShipped(table, state)
    expect_identical(got$prot, as.integer(table$first_pass))
    expect_identical(got$omics, as.integer(table$first_pass))
    expect_identical(got$trx, as.integer(table$trx))
})

test_that("KIT-193 R-N1 pair, planted defect: ties ranked last in place of first fail the first-pass table", {
    table <- kit193FirstPassTable()
    rows <- table[table$case == "tied_blocks", , drop = FALSE]
    testthat::local_mocked_bindings(protRuvFirstPassControls = function(values, groups) {
        wide <- rank(rows$p, ties.method = "last", na.last = "keep") > 0.25 * sum(is.finite(rows$p))
        wide[is.na(wide)] <- FALSE
        wide
    })
    got <- protRuvFirstPassControls(matrix(0, nrow(rows), 6L), rep(c("a", "b"), each = 3L))
    testthat::expect_failure(expect_identical(as.integer(got), as.integer(rows$first_pass)))
})
