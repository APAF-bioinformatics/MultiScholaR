# The decision classifier (KIT-159 4.5, 4.6; section 10.1: T-CLS, T-CLS-regime, T-CLS-compound).

iv <- function(lo, hi, method = "paired_t", n = 8, scope = "run_level", class = "accuracy",
               status = "ok") {
    valNewInterval(ci_low = lo, ci_high = hi, method = method, n_units = n, conf = 0.95,
                   scope = scope, class = class, status = status)
}
reg8 <- regimeRows("paired_t", 8)
cls <- function(lo, hi, ...) valClassify(iv(lo, hi), margin = 0.01, higher_is_better = TRUE,
                                          regimes = reg8, ...)

test_that("T-CLS: the four outcomes, no_worse and inclusive ends", {
    expect_identical(cls(-0.005, 0.008)$outcome, "match")
    expect_identical(cls(0.012, 0.03)$outcome, "better")
    expect_identical(cls(-0.03, -0.012)$outcome, "worse")
    a <- cls(-0.02, 0.005)
    expect_identical(a$outcome, "inconclusive")
    expect_false(a$no_worse)
    b <- cls(-0.008, 0.05)
    expect_identical(b$outcome, "inconclusive")
    expect_true(b$no_worse)
    expect_identical(cls(-0.01, 0.01)$outcome, "match")
    e <- cls(0.01, 0.03)
    expect_identical(e$outcome, "inconclusive")
    expect_true(e$no_worse)
})

test_that("T-CLS: refusals and a one-arm comparison against a fixed value", {
    expect_error(valClassify(iv(-0.005, 0.008), margin = NULL, higher_is_better = TRUE,
                             regimes = reg8),
                 class = "val_spec_invalid")
    wil <- iv(0.8, 0.95, method = "wilson", scope = "within_run")
    expect_error(valClassify(wil, margin = 0.01, higher_is_better = TRUE,
                             regimes = regimeRows("wilson", 8)),
                 class = "val_spec_invalid")
    expect_error(valClassify(iv(-0.005, 0.008), margin = 0.01, higher_is_better = TRUE,
                             view = "own", regimes = reg8),
                 class = "val_spec_invalid")
    one <- valClassify(iv(0.03, 0.045), margin = 0, higher_is_better = FALSE, fixed = 0.05,
                       view = "own", regimes = reg8)
    expect_identical(one$outcome, "better")
    zero <- iv(0.01, 0.01, status = "too_few_units")
    expect_error(valClassify(zero, margin = 0.01, higher_is_better = TRUE, regimes = reg8),
                 class = "val_spec_invalid")
})

test_that("T-CLS: a within-run interval decides only when declared on that comparison", {
    wil <- iv(0.02, 0.05, method = "wilson", scope = "within_run")
    r <- valClassify(wil, margin = 0.01, higher_is_better = TRUE,
                     regimes = regimeRows("wilson", 8), feature_interval_decides = TRUE)
    expect_identical(r$outcome, "better")
    expect_match(r$label, "conditional on the samples")
})

test_that("T-CLS-regime: no regime, no harness regime for the class, flagged without consent", {
    tk <- iv(-0.005, 0.008, method = "ticket:k16xInterval")
    expect_error(valClassify(tk, margin = 0.01, higher_is_better = TRUE, regimes = reg8),
                 class = "val_spec_invalid")
    jk <- valNewInterval(ci_low = -0.005, ci_high = 0.008, method = "jackknife", n_units = 5,
                         conf = 0.95, class = "rank", form = "replicate")
    expect_error(valClassify(jk, margin = 0.01, higher_is_better = TRUE,
                             regimes = regimeRows("jackknife", 4, class = "rank",
                                                  form = "replicate")),
                 class = "val_spec_invalid")
    fl <- regimeRows("paired_t", 8, status = "flagged")
    expect_error(valClassify(iv(-0.005, 0.008), margin = 0.01, higher_is_better = TRUE,
                             regimes = fl),
                 class = "val_spec_invalid")
    ok <- valClassify(iv(-0.005, 0.008), margin = 0.01, higher_is_better = TRUE, regimes = fl,
                      accept_flagged = list(label = "skewed case not covered",
                                            ledger_row = "VAL-99"))
    expect_identical(ok$outcome, "match")
    expect_match(ok$label, "skewed case not covered")
    expect_identical(ok$coverage, 0.95)
})

test_that("T-CLS-compound: conjunctions, worse_when, controls and guards", {
    out <- list(
        c1 = cls(-0.005, 0.008), c2 = cls(0.012, 0.03), c3 = cls(-0.03, -0.012),
        ctl_ok = cls(-0.005, 0.008), ctl_worse = cls(-0.03, -0.012)
    )
    m1 <- valClassifyMember(list(all_of = c("c1", "c2")), out)
    expect_identical(m1$outcome, "inconclusive")
    expect_true(m1$no_worse)
    m2 <- valClassifyMember(list(all_of = c("c1", "c2", "c3")), out)
    expect_identical(m2$outcome, "worse")
    m3 <- valClassifyMember(list(all_of = c("c1", "c2", "c3"),
                                 worse_when = list(c("c1", "c2"), "c3")), out)
    expect_identical(m3$outcome, "worse")
    m3b <- valClassifyMember(list(all_of = c("c1", "c3"), worse_when = list(c("c1", "c3"))), out)
    expect_identical(m3b$outcome, "inconclusive")
    m4 <- valClassifyMember(list(all_of = c("c1", "c2"), controls = list(c2 = "ctl_ok")), out)
    expect_true(m4$components$not_sensitive[m4$components$component == "c2"])
    expect_identical(m4$outcome, "match")
    m4w <- valClassifyMember(list(all_of = c("c1", "c2"), controls = list(c2 = "ctl_worse")), out)
    expect_false(m4w$components$not_sensitive[m4w$components$component == "c2"])
    m5 <- valClassifyMember(list(all_of = c("c1", "c2"),
                                 guards = list(list(component = "c1")),
                                 controls = list(c1 = "ctl_ok")), out)
    expect_identical(m5$outcome, "inconclusive")
    expect_true("c1" %in% m5$components$component)
    m6 <- valClassifyMember(list(all_of = "c2", controls = list(c2 = "ctl_ok")), out)
    expect_identical(m6$outcome, "no_verdict")
})

test_that("T-CLS-compound: a guard requiring match, and a contradicted guard", {
    out <- list(c1 = cls(-0.005, 0.008), c2 = cls(0.012, 0.03), c3 = cls(-0.03, -0.012),
                g_in = cls(-0.005, 0.008), g_shift = cls(0.012, 0.03))
    pass <- valClassifyMember(list(all_of = c("c2", "g_in"),
                                   guards = list(list(component = "g_in", require = "match"))),
                              out)
    expect_identical(pass$outcome, "better")
    fail <- valClassifyMember(list(all_of = c("c2", "g_shift"),
                                   guards = list(list(component = "g_shift", require = "match"))),
                              out)
    expect_identical(fail$outcome, "inconclusive")
    con <- valClassifyMember(list(all_of = c("c1", "c2", "c3"),
                                  guards = list(list(component = "c3"))), out)
    expect_identical(con$outcome, "inconclusive")
    expect_true("guard c3: contradicted" %in% con$guard_text)
})
