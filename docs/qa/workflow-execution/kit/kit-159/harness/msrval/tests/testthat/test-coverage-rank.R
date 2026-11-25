# The rank cells re-pinned (KIT-159 r1.13, VAL-69): the generator's effect per g calibrated so that
# the plain arm's mean normalised pAUC at FPR 0.01 lies in [0.4, 0.8]; the `degenerate` rule;
# run-level `conservative` read on the miss rate.

test_that("T-RANK-1: the rank effect is calibrated into the declared pAUC band and recorded", {
    cal <- valRankEffect("replicate", 4L, nsim = 20L)
    expect_true(cal$effect %in% msrval:::RANK_EFFECT_GRID)
    expect_gte(cal$mean_pauc, 0.4)
    expect_lte(cal$mean_pauc, 0.8)
    again <- valRankEffect("replicate", 4L, nsim = 20L)
    expect_identical(again, cal)
})

test_that("T-RANK-2: a cell with zero-width intervals or a target at a bound is degenerate", {
    cell <- list(cell = "C4-144", regime = "C4", family = "c4:4840", n = 4L,
                 condition = "matrix4840")
    hits <- matrix(TRUE, 100L, length(msrval:::COVERAGE_LEVELS))
    zero <- matrix(rep(c(TRUE, FALSE), c(50L, 50L)), 100L, length(msrval:::COVERAGE_LEVELS))
    rows <- valCoverageRows(cell, "pauc_diff", "rank", hits, "approx_run", target = 0.3,
                            zero_width = zero, target_se = 0.001, range = c(-1, 1),
                            effect = 0.35)
    expect_true(all(rows$status == "degenerate"))
    expect_equal(rows$zero_width, rep(0.5, length(msrval:::COVERAGE_LEVELS)))
    expect_equal(unique(rows$effect), 0.35)
    none <- matrix(FALSE, 100L, length(msrval:::COVERAGE_LEVELS))
    bound <- valCoverageRows(cell, "pauc_diff", "rank", hits, "approx_run", target = 0.9995,
                             zero_width = none, target_se = 0.001, range = c(-1, 1))
    expect_true(all(bound$status == "degenerate"))
    fine <- valCoverageRows(cell, "pauc_diff", "rank", hits[, ] & FALSE | TRUE, "approx_run",
                            target = 0.3, zero_width = none, target_se = 0.001, range = c(-1, 1))
    expect_false(any(fine$status == "degenerate"))
})

test_that("T-RANK-3: run-level conservative is read on the miss rate", {
    expect_identical(valCoverageStatus("approx_run", 0.97, 0.95), "pass")
    expect_identical(valCoverageStatus("approx_run", 0.99, 0.95), "conservative")
    expect_identical(valCoverageStatus("approx_run", 0.999, 0.99375), "conservative")
    expect_identical(valCoverageStatus("approx_run", 0.995, 0.99375), "pass")
})

test_that("T-RANK-4: a degenerate condition makes the regime row degenerate (not offered)", {
    row <- function(st) {
        data.frame(cell = "C4-146", regime = "C4", family = "c4:4840", n = 6L,
                   condition = "matrix4840", statistic = "pauc_diff", class = "rank",
                   level = 0.95, nsim = 2000L, coverage = 1, mc_se = 0, achieved = NA_real_,
                   undefined = 0, target = 0, status = st, kind = "approx_run",
                   stringsAsFactors = FALSE)
    }
    tab <- valSummariseCoverage(row("degenerate"))
    expect_identical(tab$status, "degenerate")
})
