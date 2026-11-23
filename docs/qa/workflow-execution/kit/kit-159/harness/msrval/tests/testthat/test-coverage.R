# The T-COV machinery (KIT-159 4.4 pass bands, 10.2): bands, the regime summary and small runs.
# The regimes' own coverage results come from T2, not from these tests.

test_that("pass bands: approximate run-level, within-run, exact and Clopper-Pearson", {
    expect_identical(valCoverageStatus("approx_run", 0.92, 0.95), "fail")
    expect_identical(valCoverageStatus("approx_run", 0.93, 0.95), "pass")
    expect_identical(valCoverageStatus("approx_run", 0.995, 0.95), "conservative")
    expect_identical(valCoverageStatus("approx_within", 0.995, 0.95), "fail")
    expect_identical(valCoverageStatus("approx_run", 1, 1 - 0.05 / 15), "pass")
    expect_identical(valCoverageStatus("exact", 0.99, 0.95, achieved = 0.9921875), "pass")
    expect_identical(valCoverageStatus("exact", 0.98, 0.95, achieved = 0.9921875), "fail")
    expect_identical(valCoverageStatus("clopper", 0.94, 0.95), "conservative")
    expect_identical(valCoverageStatus("clopper", 0.93, 0.95), "fail")
})

test_that("the summary flags a paired t count passing on normal and failing on skewed data", {
    row <- function(condition, coverage) {
        data.frame(cell = "C2-x", regime = "C2", family = "paired_t", n = 8L,
                   condition = condition, statistic = "paired_t", class = NA_character_,
                   level = 0.95, nsim = 2000L, coverage = coverage, mc_se = 0.005,
                   achieved = NA_real_, undefined = 0, target = 0,
                   status = valCoverageStatus("approx_run", coverage, 0.95), kind = "approx_run",
                   stringsAsFactors = FALSE)
    }
    tab <- valSummariseCoverage(rbind(row("normal", 0.948), row("skewed", 0.90)))
    expect_identical(tab$status, "flagged")
    expect_equal(tab$coverage, 0.90)
    both <- valSummariseCoverage(rbind(row("normal", 0.90), row("skewed", 0.90)))
    expect_identical(both$status, "fail")
})

test_that("the catalogue has unique cells and every regime of 10.2", {
    cells <- valCoverageCells()
    expect_false(anyDuplicated(cells$cell) > 0L)
    expect_setequal(unique(cells$regime), c("C1", "C2", "C3", "C4", "C4e", "C5", "C6", "C7", "C8"))
})

test_that("small runs of a C1 and a C3 cell return one row per level", {
    cells <- valCoverageCells()
    c1 <- valCoverageRegime(cells[cells$regime == "C1", ][1L, ], nsim = 30L)
    expect_identical(nrow(c1), 8L)
    offered <- c1$status != "not_offered"
    expect_identical(offered, c(TRUE, rep(FALSE, 7)))
    expect_equal(c1$achieved[offered], 0.96875)
    c3 <- valCoverageRegime(cells[cells$regime == "C3", ][1L, ], nsim = 5L, target_sims = 50L)
    expect_identical(nrow(c3), 8L)
    expect_true(all(c3$coverage >= 0 & c3$coverage <= 1))
})
