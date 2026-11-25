# Hand-computed metric cases (KIT-159 section 10.1: T-M1 to T-M10).

test_that("T-M1/M2: bias and median absolute error with one missing estimate", {
    s <- scoredRows(
        key = paste0("e", 1:5), expected = rep(-2, 5), class = "ecoli",
        estimate = c(-1.5, -2.0, -2.5, -1.0, NA), p = rep(0.01, 5)
    )
    b <- valBias(s, class = "ecoli", n_min = 1)
    m <- valMedianAbsError(s, class = "ecoli", n_min = 1)
    expect_equal(b$value, 0.25, tolerance = 1e-12)
    expect_equal(m$value, 0.5, tolerance = 1e-12)
    expect_identical(b$n, 4L)
    expect_identical(b$n_missing, 1L)
    expect_identical(m$n, 4L)
    expect_identical(m$n_missing, 1L)
    b5 <- valBias(s, class = "ecoli", n_min = 5)
    m5 <- valMedianAbsError(s, class = "ecoli", n_min = 5)
    expect_true(is.na(b5$value))
    expect_true(is.na(m5$value))
    expect_identical(b5$status, "too_few")
    expect_identical(m5$status, "too_few")
})

test_that("T-M2-pooled: pooled median and the mean of per-class medians", {
    s <- scoredRows(
        key = paste0("k", 1:8), expected = rep(0, 8),
        estimate = c(0.1, 0.2, 0.3, 1, 1, 1, 1, 1), p = rep(0.5, 8),
        class = c(rep("A", 3), rep("B", 5))
    )
    m <- valMedianAbsError(s, class = NULL, n_min = 1)
    expect_equal(m$value, 1, tolerance = 1e-12)
    expect_equal(m$mean_of_classes, 0.6, tolerance = 1e-12)
    b <- valBias(s, class = NULL, n_min = 1)
    expect_equal(b$value, 1, tolerance = 1e-12)
    expect_equal(b$mean_of_classes, 0.6, tolerance = 1e-12)
})

test_that("T-M3: linear and geometric CV; a group with one finite value is skipped", {
    x <- matrix(c(10, 11, 12), nrow = 1, dimnames = list("f1", c("s1", "s2", "s3")))
    g <- c("a", "a", "a")
    expect_equal(valCv(x, g, scale = "linear", n_min = 1)$value, 0.6546536707, tolerance = 1e-9)
    expect_equal(valCv(x, g, scale = "geometric", n_min = 1)$value, 0.7853704045, tolerance = 1e-9)
    x2 <- matrix(c(10, 11, 12, 9, NA), nrow = 1,
                 dimnames = list("f1", c("s1", "s2", "s3", "s4", "s5")))
    g2 <- c("a", "a", "a", "b", "b")
    expect_equal(valCv(x2, g2, scale = "linear", n_min = 1)$value, 0.6546536707, tolerance = 1e-9)
})

test_that("T-M4/M5/M6: confusion, FDP and sensitivity in the own and union views", {
    conf <- valConfusion(tm4Rows(), alpha = 0.05, call_rule = "own", positive_rule = "expected")
    own <- conf[conf$view == "own", ]
    uni <- conf[conf$view == "union", ]
    expect_identical(own$TP, 2L)
    expect_identical(own$FP, 1L)
    expect_identical(own$TN, 5L)
    expect_identical(own$FN, 1L)
    expect_identical(uni$FN, 2L)
    expect_identical(uni$no_statistic, 0L)
    expect_identical(uni$sign_errors, 0L)
    expect_equal(valFdp(conf)$value, 1 / 3, tolerance = 1e-12)
    expect_equal(valSensitivity(conf, view = "own")$value, 2 / 3, tolerance = 1e-12)
    expect_equal(valSensitivity(conf, view = "union")$value, 2 / 4, tolerance = 1e-12)
})

test_that("T-M4-sign: a wrong-sign call is a false discovery and its key stays a miss", {
    conf <- valConfusion(tm4Rows(estimate_f2 = -0.8), alpha = 0.05, call_rule = "own")
    own <- conf[conf$view == "own", ]
    uni <- conf[conf$view == "union", ]
    expect_identical(own$TP, 1L)
    expect_identical(own$FP, 2L)
    expect_identical(own$sign_errors, 1L)
    expect_identical(own$FN, 2L)
    expect_identical(uni$FN, 3L)
    expect_identical(uni$TP + uni$FN, 4L)
    expect_identical(uni$FP + uni$TN, 6L + uni$sign_errors)
    expect_equal(valFdp(conf)$value, 2 / 3, tolerance = 1e-12)
})

test_that("T-M4-nostat: a finite estimate without p is quantified, not testable, and a union miss", {
    s <- tm4Rows()
    s$p[3] <- NA
    s$s_own[3] <- NA
    s$testable <- is.finite(s$p)
    conf <- valConfusion(s, alpha = 0.05, call_rule = "own")
    uni <- conf[conf$view == "union", ]
    own <- conf[conf$view == "own", ]
    expect_identical(uni$no_statistic, 1L)
    expect_true(s$quantified[3])
    expect_false(s$testable[3])
    expect_identical(uni$FN, 2L)
    expect_identical(own$FN, 0L)
})

test_that("T-M5-zero: no call gives FDP 0; an all-null truth with one call gives FDP 1", {
    s <- scoredRows(key = paste0("k", 1:4), expected = c(1, 1, 0, 0),
                    estimate = c(1, 1, 0, 0), p = rep(0.5, 4))
    f <- valFdp(valConfusion(s, alpha = 0.05))
    expect_identical(f$value, 0)
    expect_true(f$no_calls)
    s0 <- scoredRows(key = paste0("k", 1:4), expected = rep(0, 4),
                     estimate = c(1, 0, 0, 0), p = c(0.001, 0.5, 0.5, 0.5))
    f0 <- valFdp(valConfusion(s0, alpha = 0.05))
    expect_identical(f0$value, 1)
    expect_false(f0$no_calls)
})

test_that("T-M7: partial AUC from the origin, normalised by fpr_max, and the null floor", {
    s <- rankRows(c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE),
                  c(0.002, 0.005, 0.02, 0.03, 0.05, 0.08))
    a <- valPartialAuc(s, fpr_max = 0.5, min_nulls_in_region = 0, min_events = 1)
    expect_equal(a$raw, 0.375, tolerance = 1e-12)
    expect_equal(a$value, 0.75, tolerance = 1e-12)
    expect_equal(a$full_auc, 0.875, tolerance = 1e-12)
    a25 <- valPartialAuc(s, fpr_max = 0.25, min_nulls_in_region = 0, min_events = 1)
    expect_equal(a25$value, 0.5, tolerance = 1e-12)
    expect_error(valPartialAuc(s, fpr_max = 0.5, min_events = 1), class = "val_spec_invalid")
})

test_that("T-M7-ties: a tied block is one straight segment", {
    s <- rankRows(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
                  c(0.002, 0.02, 0.01, 0.03, 0.05, 0.08))
    a <- valPartialAuc(s, fpr_max = 0.5, min_nulls_in_region = 0, min_events = 1)
    expect_equal(a$raw, 0.4375, tolerance = 1e-12)
    expect_equal(a$value, 0.875, tolerance = 1e-12)
})

test_that("T-M7-sign: a wrong-sign positive ranks by its one-sided p near the end", {
    s <- scoredRows(
        key = c("a", "b", "c", "d", "e", "f"),
        expected = c(1, 1, 0, 0, 0, 0),
        estimate = c(0.5, -0.5, 0.1, -0.1, 0.2, -0.2),
        p = c(0.01, 0.02, 0.04, 0.3, 0.5, 0.8)
    )
    a <- valPartialAuc(s, fpr_max = 0.5, min_nulls_in_region = 0, min_events = 1)
    expect_equal(a$ranking_p, c(a = 0.005, b = 0.99, c = 0.04, d = 0.3, e = 0.5, f = 0.8),
                 tolerance = 1e-12)
    expect_equal(a$full_auc, 0.5, tolerance = 1e-12)
    expect_equal(a$value, 0.5, tolerance = 1e-12)
    expect_equal(a$roc$fpr, c(0, 0, 0.25, 0.5, 0.75, 1, 1), tolerance = 1e-12)
    expect_equal(a$roc$tpr, c(0, 0.5, 0.5, 0.5, 0.5, 0.5, 1), tolerance = 1e-12)
})

test_that("T-M7-prev: 5 nulls in the region are refused at the default floor", {
    s <- rankRows(c(rep(TRUE, 10), rep(FALSE, 50)), seq(0.001, 0.6, length.out = 60))
    expect_error(valPartialAuc(s, fpr_max = 0.1, min_events = 1), class = "val_spec_invalid")
})

test_that("T-M8: Pearson and Spearman correlations", {
    s <- scoredRows(key = paste0("k", 1:4), expected = c(1, 3, 2, 4),
                    estimate = c(1, 2, 3, 4), p = rep(0.01, 4))
    expect_equal(valCorrelation(s, method = "pearson", n_min = 1)$value, 0.8, tolerance = 1e-12)
    expect_equal(valCorrelation(s, method = "spearman", n_min = 1)$value, 0.8, tolerance = 1e-12)
    expect_true(is.na(valCorrelation(s, method = "pearson", n_min = 5)$value))
})

test_that("T-M9: Quartet SNR from scores", {
    scores <- cbind(PC1 = c(-2, -1, 1, 2), PC2 = c(0.5, -0.5, 0.5, -0.5))
    snr <- valSnrFromScores(scores, weights = c(0.8, 0.2), groups = c("A", "A", "B", "B"))
    expect_equal(snr, 8.864907252, tolerance = 1e-9)
})

test_that("T-M9-pca: the matrix form equals the score form with prcomp's weights", {
    withr::local_seed(159009)
    x <- matrix(rnorm(50 * 12), nrow = 50, dimnames = list(paste0("f", 1:50), paste0("s", 1:12)))
    groups <- rep(c("D5", "D6", "F7", "M8"), each = 3)
    z <- t(scale(t(x)))
    pc <- stats::prcomp(t(z), scale. = FALSE)
    w5 <- summary(pc)$importance[2, 1:2]
    w_exact <- (pc$sdev^2 / sum(pc$sdev^2))[1:2]
    expect_equal(valQuartetSnr(x, groups, weights = "rounded5")$value,
                 valSnrFromScores(pc$x[, 1:2], w5, groups), tolerance = 1e-12)
    expect_equal(valQuartetSnr(x, groups, weights = "exact")$value,
                 valSnrFromScores(pc$x[, 1:2], w_exact, groups), tolerance = 1e-12)
})

test_that("T-M10: Matthews correlation from the confusion table", {
    s <- scoredRows(key = paste0("k", 1:6), expected = c(1, 0, 0, 0, 1, 1),
                    estimate = rep(1, 6), p = c(0.01, 0.01, 0.5, 0.5, 0.01, 0.5))
    conf <- valConfusion(s, alpha = 0.05)
    own <- conf[conf$view == "own", ]
    expect_identical(c(own$TP, own$FP, own$FN, own$TN), c(2L, 1L, 1L, 2L))
    expect_equal(valMcc(conf, view = "own")$value, 1 / 3, tolerance = 1e-9)
})

test_that("T-M10-sign: a wrong-sign key counts in FP and FN", {
    conf <- valConfusion(tm4Rows(estimate_f2 = -0.8), alpha = 0.05)
    own <- conf[conf$view == "own", ]
    expect_identical(c(own$TP, own$FP, own$FN, own$TN), c(1L, 2L, 2L, 5L))
    expect_identical(own$TP + own$FP + own$FN + own$TN, 10L)
    expect_equal(valMcc(conf, view = "own")$value, 1 / 21, tolerance = 1e-9)
})

test_that("M11: coverage counts quantified, testable and called keys per class", {
    s <- tm4Rows()
    s$class <- c(rep("changed", 4), rep("null", 6))
    x <- matrix(c(1, NA, 3, 4), nrow = 2)
    cov <- valCoverage(s, x, alpha = 0.05)
    ch <- cov[cov$class == "changed", ]
    expect_identical(ch$quantified, 3L)
    expect_identical(ch$testable, 3L)
    expect_identical(ch$called_own, 2L)
    expect_equal(unique(cov$completeness), 0.75)
})
