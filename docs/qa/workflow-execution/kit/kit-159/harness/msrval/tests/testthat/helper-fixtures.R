# Fixtures shared by the scoring tests (KIT-159 section 10.1). Every table is built by hand here so
# that each expected value in the tests can be checked against the arithmetic in the ticket.

# A scored table (section 4.1) for one truth, unit, arm and contrast.
scoredRows <- function(key, expected, estimate, p, s_own = p, class = "c",
                       arm = "a1", assay = "all", contrast = "B_vs_A",
                       truth_status = "scored", key_status = "ok", s_bh = NA_real_) {
    n <- length(key)
    data.frame(
        truth_id = "primary",
        unit = "u01",
        arm = rep(arm, length.out = n),
        assay = rep(assay, length.out = n),
        contrast = contrast,
        key = as.character(key),
        key_status = rep(key_status, length.out = n),
        class = rep(class, length.out = n),
        truth_status = rep(truth_status, length.out = n),
        expected = expected,
        estimate = estimate,
        p = p,
        s_own = s_own,
        s_bh = rep(s_bh, length.out = n),
        called_own = NA,
        called_bh = NA,
        quantified = is.finite(estimate),
        testable = is.finite(p),
        stringsAsFactors = FALSE
    )
}

# The T-M4 table: P = f1..f4 (theta = +1), N = f5..f10 (theta = 0).
tm4Rows <- function(estimate_f2 = 0.9) {
    scoredRows(
        key = paste0("f", 1:10),
        expected = c(rep(1, 4), rep(0, 6)),
        estimate = c(1.1, estimate_f2, 0.4, NA, 0.2, 0.1, -0.1, 0.3, -0.2, 0.05),
        p = c(0.001, 0.02, 0.2, NA, 0.01, 0.5, 0.9, 0.3, 0.6, 0.7)
    )
}

# Ranking fixture for M7: keys in the order given, positives expected +1 with positive estimates
# (one-sided ranking p = p / 2), nulls expected 0 (two-sided p as given).
rankRows <- function(is_pos, p) {
    scoredRows(
        key = paste0("k", seq_along(is_pos)),
        expected = ifelse(is_pos, 1, 0),
        estimate = ifelse(is_pos, 1, 0.1),
        p = p
    )
}

# A regime table as `coverage-table.tsv` holds it (section 10.2), for the classifier tests.
regimeRows <- function(method, n, level = 0.95, class = NA_character_, form = NA_character_,
                       status = "pass", source = "harness", regime = "test") {
    data.frame(
        regime = regime, method = method, class = class, form = form, n = n, level = level,
        status = status, coverage = level, mc_se = 0.005, source = source,
        statistic = "test fixture", stringsAsFactors = FALSE
    )
}

# The ticket prints hand-computed values to a fixed number of decimals; a value matches when it is
# within one unit of the last printed decimal (the printed value is rounded, so the true one lies
# within half a unit).
expectClose <- function(actual, expected, digits) {
    testthat::expect_true(
        all(abs(actual - expected) <= 10^-digits),
        info = paste("actual", paste(format(actual, digits = 15), collapse = ", "))
    )
}
