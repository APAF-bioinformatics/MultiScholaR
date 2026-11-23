# Interval methods (KIT-159 4.4; section 10.1: T-ORD, T-PT, T-JK-*, T-CP, T-ONE, T-WIL, T-NEW,
# T-BOOT-1 to T-BOOT-4, T-REF-FILL).

d8 <- c(0.02, -0.01, 0.03, 0, 0.01, 0.04, -0.02, 0.01)

test_that("T-ORD: order-statistic intervals and their achieved levels", {
    a <- valInterval(c(-0.02, -0.01, 0, 0.01, 0.01, 0.02, 0.03, 0.04), method = "exact_median",
                     conf = 0.95)
    expect_equal(c(a$ci_low, a$ci_high), c(-0.02, 0.04), tolerance = 1e-12)
    expect_equal(a$achieved_level, 0.9921875, tolerance = 1e-12)
    b <- valInterval(1:13, method = "exact_median", conf = 0.95)
    expect_equal(c(b$ci_low, b$ci_high), c(3, 11))
    expect_equal(b$achieved_level, 0.9775390625, tolerance = 1e-12)
    c5 <- valInterval(1:5, method = "exact_median", conf = 0.95)
    expect_identical(c5$status, "too_few_units")
    expect_true(is.na(c5$ci_low))
    expect_identical(a$scope, "run_level")
})

test_that("T-PT: paired t on eight per-unit differences", {
    a <- valInterval(d8, method = "paired_t", conf = 0.95)
    expect_equal(a$estimate, 0.01, tolerance = 1e-12)
    expectClose(a$se, 0.0070710678, 10)
    expectClose(a$critical, 2.364624252, 9)
    expectClose(c(a$ci_low, a$ci_high), c(-0.0067204, 0.0267204), 7)
    p <- valInterval(d8 + 1, y = rep(1, 8), method = "paired_t", conf = 0.95)
    expect_equal(c(p$ci_low, p$ci_high), c(a$ci_low, a$ci_high), tolerance = 1e-12)
})

test_that("T-JK-1: the delete-one jackknife of the mean has SE sd / sqrt(n)", {
    j <- valJackknife(values = d8, statistic = mean, form = "unit", conf = 0.95)
    expect_equal(j$se, stats::sd(d8) / sqrt(8), tolerance = 1e-12)
    expectClose(j$se, 0.0070710678, 10)
})

test_that("T-JK-2: jackknife from full and deletion values", {
    j <- valJackknife(full = 0.031, deletions = c(0.028, 0.035, 0.030, 0.033), form = "replicate",
                      conf = 0.95)
    expectClose(j$se, 0.0046636895, 10)
    expectClose(j$critical, 3.182446305, 9)
    expectClose(c(j$ci_low, j$ci_high), c(0.0161581, 0.0458419), 7)
})

test_that("T-JK-2-bias: the bias estimate and its flag against half the margin", {
    j2 <- valJackknife(full = 0.031, deletions = c(0.028, 0.035, 0.030, 0.033), form = "replicate",
                       conf = 0.95, margin = 0.002)
    expect_equal(j2$bias, 0.0015, tolerance = 1e-12)
    expect_true("depends_on_replicate_count" %in% j2$flags)
    j4 <- valJackknife(full = 0.031, deletions = c(0.028, 0.035, 0.030, 0.033), form = "replicate",
                       conf = 0.95, margin = 0.004)
    expect_false("depends_on_replicate_count" %in% j4$flags)
})

test_that("T-JK-log: the log transform and an undefined deletion", {
    j <- valJackknife(full = 0.05, deletions = c(0.04, 0.06, 0.05, 0.05), form = "replicate",
                      conf = 0.95, transform = "log")
    expectClose(j$se, 0.2489241, 7)
    expectClose(c(j$ci_low, j$ci_high), c(0.0226427, 0.1104111), 7)
    z <- valJackknife(full = 0.05, deletions = c(0.04, 0.06, 0, 0.05), form = "replicate",
                      conf = 0.95, transform = "log")
    expect_identical(z$status, "undefined_deletion")
    expect_true(is.na(z$ci_low))
})

test_that("T-JK-crit: a declared critical value may be larger than t, never smaller", {
    del <- c(0.028, 0.035, 0.030, 0.033)
    expect_error(valJackknife(full = 0.031, deletions = del, form = "replicate", conf = 0.95,
                              critical = 3.0),
                 class = "val_spec_invalid")
    j <- valJackknife(full = 0.031, deletions = del, form = "replicate", conf = 0.95,
                      critical = 4.0)
    expectClose(c(j$ci_low, j$ci_high), c(0.0123452, 0.0496548), 7)
    expect_true("calibrated_critical" %in% j$flags)
})

test_that("T-JK-3: the jackknife of a median is refused", {
    expect_error(valJackknife(values = d8, statistic = stats::median, form = "unit"),
                 class = "val_spec_invalid")
    expect_error(valJackknife(values = d8, statistic = "median", form = "unit"),
                 class = "val_spec_invalid")
})

test_that("the jackknife needs four deletions; dependent units are refused", {
    j <- valJackknife(full = 0.031, deletions = c(0.028, 0.035, 0.030), form = "replicate")
    expect_identical(j$status, "too_few_units")
    expect_error(valInterval(d8, method = "paired_t", unit_role = "dependent"),
                 class = "val_spec_invalid")
})

test_that("T-CP: Clopper-Pearson intervals", {
    a <- valInterval(20, n = 21, method = "clopper_pearson", conf = 0.95)
    expectClose(c(a$ci_low, a$ci_high), c(0.7618401, 0.9987951), 7)
    b <- valInterval(0, n = 10, method = "clopper_pearson", conf = 0.95)
    expectClose(c(b$ci_low, b$ci_high), c(0, 0.3084971), 7)
    c3 <- valInterval(3, n = 250, method = "clopper_pearson", conf = 0.95)
    expectClose(c(c3$ci_low, c3$ci_high), c(0.0024816, 0.0346666), 7)
})

test_that("T-ONE: the one-sample paired t against a fixed value", {
    a <- valInterval(d8, method = "paired_t", conf = 0.95, fixed = 0)
    b <- valInterval(d8, method = "paired_t", conf = 0.95)
    expect_equal(c(a$ci_low, a$ci_high), c(b$ci_low, b$ci_high), tolerance = 1e-12)
    f <- valInterval(d8, method = "paired_t", conf = 0.95, fixed = 0.01)
    expectClose(c(f$ci_low, f$ci_high), c(-0.0167204, 0.0167204), 7)
})

test_that("T-WIL: Wilson intervals without an event floor", {
    a <- valInterval(21, n = 21, method = "wilson", conf = 0.95)
    expectClose(c(a$ci_low, a$ci_high), c(0.8453610, 1), 7)
    b <- valInterval(20, n = 21, method = "wilson", conf = 0.95)
    expectClose(c(b$ci_low, b$ci_high), c(0.7733064, 0.9915440), 7)
    expect_identical(a$scope, "within_run")
})

test_that("T-NEW: Newcombe's paired interval, method 10", {
    a <- valInterval(c(both = 18, first_only = 2, second_only = 0, neither = 1),
                     method = "newcombe_paired", conf = 0.95)
    expectClose(a$estimate, 0.0952381, 7)
    expect_equal(a$phi, 7.5 / sqrt(1080), tolerance = 1e-9)
    expectClose(c(a$ci_low, a$ci_high), c(-0.0867513, 0.2935048), 7)
})

wideErrors <- function(n = 200) {
    withr::with_seed(159013, {
        e <- abs(stats::rnorm(n)) + 0.01
    })
    data.frame(key = paste0("k", seq_len(n)), class = rep(c("x", "y"), length.out = n),
               cluster = rep(seq_len(n / 4), each = 4), err_a = e, err_b = e + 0.1)
}

test_that("T-BOOT-1: the same seed and B give identical intervals", {
    w <- wideErrors()
    stat <- function(d) stats::median(abs(d$err_b)) - stats::median(abs(d$err_a)) +
        stats::sd(d$err_a)
    a <- valBootstrap(w, stat, B = 2000, conf = 0.95, seed = 159001)
    b <- valBootstrap(w, stat, B = 2000, conf = 0.95, seed = 159001)
    expect_identical(c(a$ci_low, a$ci_high), c(b$ci_low, b$ci_high))
    expect_identical(a$scope, "within_run")
})

test_that("T-BOOT-2: a constant paired shift gives a zero-width interval at the shift", {
    w <- wideErrors()
    stat <- function(d) stats::median(abs(d$err_b)) - stats::median(abs(d$err_a))
    a <- valBootstrap(w, stat, B = 2000, conf = 0.95, seed = 159002)
    expect_equal(c(a$ci_low, a$ci_high), c(0.1, 0.1), tolerance = 1e-12)
})

test_that("T-BOOT-3: strata keep their sizes; clusters are resampled whole", {
    w <- wideErrors()
    sizes <- new.env()
    sizes$ok <- TRUE
    stat <- function(d) {
        tab <- table(d$class)
        if (!identical(as.integer(tab[c("x", "y")]), c(100L, 100L))) sizes$ok <- FALSE
        mean(d$err_a)
    }
    valBootstrap(w, stat, B = 2000, conf = 0.95, seed = 159003, strata = "class")
    expect_true(sizes$ok)
    whole <- new.env()
    whole$ok <- TRUE
    stat_c <- function(d) {
        per <- tapply(d$key, d$cluster, function(k) {
            length(unique(k)) == 4L && length(unique(as.integer(table(k)))) == 1L
        })
        counts <- table(d$cluster) %% 4L
        if (!all(per) || any(counts != 0L)) whole$ok <- FALSE
        mean(d$err_a)
    }
    valBootstrap(w, stat_c, B = 2000, conf = 0.95, seed = 159004, cluster = "cluster")
    expect_true(whole$ok)
})

test_that("T-BOOT-4: B below 50 / tail probability is refused", {
    w <- wideErrors()
    expect_error(valBootstrap(w, function(d) mean(d$err_a), B = 1000, conf = 0.95, seed = 1),
                 class = "val_spec_invalid")
    expect_identical(valBootstrapMinB(0.975), 4000)
    expect_identical(valBootstrapMinB(1 - 0.05 / 6), 12000)
})

test_that("T-REF-FILL: the least favourable fill decides whether the verdict stands", {
    d <- c(seq(0.001, 0.029, by = 0.001), NA)
    reg <- rbind(regimeRows("paired_t", 29), regimeRows("paired_t", 30))
    r <- valRefusalFill(d, tolerance = 1, higher_is_better = TRUE, margin = 0.01,
                        method = "paired_t", conf = 0.95, regimes = reg)
    expectClose(c(r$interval_dropped$ci_low, r$interval_dropped$ci_high),
                c(0.0117612, 0.0182388), 7)
    expectClose(c(r$interval_filled$ci_low, r$interval_filled$ci_high),
                c(0.0112667, 0.0178000), 7)
    expect_identical(r$outcome, "better")
    expect_identical(r$status, "ok")
    r2 <- valRefusalFill(d, tolerance = 1, higher_is_better = TRUE, margin = 0.0115,
                         method = "paired_t", conf = 0.95, regimes = reg)
    expect_identical(r2$outcome_dropped, "better")
    expect_identical(r2$outcome_filled, "inconclusive")
    expect_identical(r2$status, "not_evaluable")
    d2 <- c(seq(0.001, 0.028, by = 0.001), NA, NA)
    r3 <- valRefusalFill(d2, tolerance = 1, higher_is_better = TRUE, margin = 0.01,
                         method = "paired_t", conf = 0.95, regimes = reg)
    expect_identical(r3$status, "not_evaluable")
})

test_that("refusals under 30 units: one is dropped with a substitute required, two are not evaluable", {
    reg <- regimeRows("paired_t", 7)
    one <- valRefusalFill(c(d8[1:7], NA), tolerance = 0, higher_is_better = TRUE, margin = 0.01,
                          method = "paired_t", conf = 0.95, regimes = reg)
    expect_identical(one$n_refused, 1L)
    expect_true(one$refusal_substitute_required)
    two <- valRefusalFill(c(d8[1:6], NA, NA), tolerance = 0, higher_is_better = TRUE,
                          margin = 0.01, method = "paired_t", conf = 0.95, regimes = reg)
    expect_identical(two$status, "not_evaluable")
})
