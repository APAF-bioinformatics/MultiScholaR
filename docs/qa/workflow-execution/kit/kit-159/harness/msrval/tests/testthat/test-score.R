# The scored table, the matched BH rule, views, key status and blinding (KIT-159 4.1, 4.3, 7.8;
# section 10.1: T-BH, T-BH-assay, T-VIEW, T-KEY, T-BLIND-1).

armTable <- function(key, p, assay = "all", estimate = 1, feature_id = key, key_status = "ok",
                     contrast = "B_vs_A") {
    data.frame(feature_id = as.character(feature_id), key = as.character(key),
               key_status = key_status, assay = assay, contrast = contrast,
               estimate = estimate, p = p, s_own = p, stringsAsFactors = FALSE)
}

test_that("T-BH: the matched rule is Benjamini-Hochberg of the raw p", {
    arms <- list(a1 = armTable(paste0("k", 1:4), c(0.01, 0.04, 0.03, 0.2)))
    universe <- data.frame(key = paste0("k", 1:4), assay = "all")
    truth <- data.frame(key = paste0("k", 1:4), class = "c", expected = 0, truth_status = "scored")
    s <- valScoreTable(arms, truth, universe, alpha = 0.05)
    s <- s[order(s$key), ]
    expect_equal(s$s_bh, c(0.04, 0.0533333333, 0.0533333333, 0.2), tolerance = 1e-9)
    expect_identical(s$called_bh, c(TRUE, FALSE, FALSE, FALSE))
    expect_named(s, c("truth_id", "unit", "arm", "assay", "contrast", "key", "key_status", "class",
                      "truth_status", "expected", "estimate", "p", "s_own", "s_bh", "called_own",
                      "called_bh", "quantified", "testable"))
})

test_that("T-BH-assay: one family per assay; the same key in two assays is two rows", {
    arms <- list(a1 = armTable(c("x", "y", "x", "z"), c(0.01, 0.04, 0.03, 0.2),
                               assay = c("lc", "lc", "gc", "gc")))
    universe <- data.frame(key = c("x", "y", "x", "z"), assay = c("lc", "lc", "gc", "gc"))
    truth <- data.frame(key = c("x", "y", "z"), class = "c", expected = 0, truth_status = "scored")
    s <- valScoreTable(arms, truth, universe, alpha = 0.05)
    lc <- s[s$assay == "lc", ]
    gc <- s[s$assay == "gc", ]
    expect_equal(lc$s_bh[order(lc$key)], c(0.02, 0.04), tolerance = 1e-12)
    expect_equal(gc$s_bh[order(gc$key)], c(0.06, 0.2), tolerance = 1e-12)
    expect_identical(sum(s$key == "x"), 2L)
    expect_true(all(s$key_status == "ok"))
})

test_that("T-VIEW: intersection, union, and union metrics unchanged by a third arm", {
    universe <- data.frame(key = paste0("k", 1:5), assay = "all")
    truth <- data.frame(key = paste0("k", 1:5), class = "c", expected = c(1, 1, 0, 0, 1),
                        truth_status = "scored")
    a1 <- armTable(paste0("k", 1:3), c(0.01, 0.5, 0.5))
    a2 <- armTable(paste0("k", 2:4), c(0.01, 0.5, 0.01))
    a3 <- armTable("k5", 0.01)
    s12 <- valScoreTable(list(a1 = a1, a2 = a2), truth, universe)
    inter <- valViews(s12, arms = c("a1", "a2"), view = "intersection", metric_type = "detection")
    expect_setequal(unique(inter$key), c("k2", "k3"))
    uni <- valViews(s12, arms = c("a1", "a2"), view = "union", metric_type = "detection")
    expect_setequal(unique(uni$key), paste0("k", 1:5))
    sens <- function(s) {
        u <- valViews(s, arms = c("a1", "a2"), view = "union", metric_type = "detection")
        vapply(c("a1", "a2"), function(a) {
            valSensitivity(valConfusion(u[u$arm == a, ], alpha = 0.05), view = "union")$value
        }, numeric(1))
    }
    s123 <- valScoreTable(list(a1 = a1, a2 = a2, a3 = a3), truth, universe)
    expect_identical(sens(s12), sens(s123))
})

test_that("T-KEY: a duplicate in one arm and an unmapped key in another leave the comparison", {
    universe <- data.frame(key = c("6", "7", "8"), assay = "all")
    truth <- data.frame(key = c("6", "7", "8"), class = "c", expected = 1, truth_status = "scored")
    a1 <- armTable(c("6", "7", "7", "8"), c(0.01, 0.01, 0.02, 0.01),
                   feature_id = c("p6", "p7a", "p7b", "p8"))
    a2 <- armTable(c("6", "7", "8"), c(0.01, 0.01, NA), feature_id = c("q6", "q7", "q8"),
                   key_status = c("ok", "ok", "unmapped"))
    s <- valScoreTable(list(a1 = a1, a2 = a2), truth, universe)
    expect_identical(s$key_status[s$arm == "a1" & s$key == "7"], "duplicate")
    v <- valViews(s, arms = c("a1", "a2"), view = "union", metric_type = "detection")
    expect_setequal(unique(v$key), "6")
    ex <- attr(v, "excluded")
    expect_identical(ex$n[ex$reason == "duplicate"], 1L)
    expect_identical(ex$n[ex$reason == "unmapped"], 1L)
})

test_that("T-VIEW: an own-view comparison between arms is refused by the classifier", {
    iv <- valNewInterval(ci_low = -0.005, ci_high = 0.008, method = "paired_t", n_units = 8,
                         conf = 0.95)
    expect_error(valClassify(iv, margin = 0.01, higher_is_better = TRUE, view = "own",
                             regimes = regimeRows("paired_t", 8)),
                 class = "val_spec_invalid")
})

test_that("T-BLIND-1: two truths with one key set and label multiset blind identically", {
    keys <- paste0("k", 1:8)
    t1 <- data.frame(key = keys, class = rep(c("changed", "null"), each = 4),
                     expected = rep(c(1, 0), each = 4), truth_status = "scored")
    t2 <- t1
    t2$class <- rep(c("null", "changed"), 4)
    t2$expected <- rep(c(0, 1), 4)
    b1 <- valBlindTruth(t1, seed = 159101)
    b2 <- valBlindTruth(t2, seed = 159101)
    expect_identical(b1, b2)
    expect_true(all(b1$blinded))
    expect_identical(sort(paste(b1$class, b1$expected)), sort(paste(t1$class, t1$expected)))
    arm <- list(a1 = armTable(keys, c(0.001, 0.2, 0.01, 0.5, 0.03, 0.6, 0.04, 0.9)))
    universe <- data.frame(key = keys, assay = "all")
    c1 <- valConfusion(valScoreTable(arm, b1, universe), alpha = 0.05)
    c2 <- valConfusion(valScoreTable(arm, b2, universe), alpha = 0.05)
    expect_identical(c1, c2)
})
