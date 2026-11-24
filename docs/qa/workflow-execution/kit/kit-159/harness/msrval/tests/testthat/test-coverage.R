# The T-COV machinery (KIT-159 4.4 pass bands, 10.2): bands, the regime summary and small runs.
# The regimes' own coverage results come from T2, not from these tests.

test_that("pass bands: approximate run-level, within-run, exact and Clopper-Pearson", {
    expect_identical(valCoverageStatus("approx_run", 0.92, 0.95), "fail")
    expect_identical(valCoverageStatus("approx_run", 0.93, 0.95), "pass")
    expect_identical(valCoverageStatus("approx_run", 0.995, 0.95), "conservative")
    expect_identical(valCoverageStatus("approx_within", 0.995, 0.95), "fail")
    # r1.13 (VAL-69): run-level `conservative` is read on the miss rate (misses below a quarter of
    # the nominal miss rate), so it can be seen at levels where level + 0.04 exceeds 1
    expect_identical(valCoverageStatus("approx_run", 1, 1 - 0.05 / 15), "conservative")
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

# The range split of T2 (coordinator's ruling of 2026-11-23, point 3): one seed per simulation and
# per target simulation, fixed by the cell and the index alone, so a cell run as ranges gives the
# same simulated data sets, counts and coverage as the cell run whole.
splitCompare <- function(cell, nsim, ntarget) {
    drop_time <- function(x) x[, setdiff(names(x), c("seconds_target", "seconds_sims"))]
    whole_part <- valCoverageRange(cell, seq_len(nsim), seq_len(ntarget))
    whole <- valCombineCoverage(cell, list(whole_part), nsim, ntarget)
    half <- nsim %/% 2L
    th <- ntarget %/% 2L
    parts <- list(valCoverageRange(cell, (half + 1L):nsim, if (ntarget) (th + 1L):ntarget),
                  valCoverageRange(cell, seq_len(half), if (ntarget) seq_len(th)))
    split <- valCombineCoverage(cell, parts, nsim, ntarget)
    expect_identical(drop_time(split), drop_time(whole))
    expect_identical(attr(split, "data_digests"), attr(whole, "data_digests"))
    expect_identical(length(attr(whole, "data_digests")), as.integer(nsim))
    expect_error(valCombineCoverage(cell, parts[1], nsim, ntarget), class = "val_spec_invalid")
    expect_error(valCombineCoverage(cell, c(parts, parts[2]), nsim, ntarget),
                 class = "val_spec_invalid")
}

test_that("a coverage cell run as ranges equals the cell run whole; a dropped or repeated range is refused", {
    cells <- valCoverageCells()
    splitCompare(cells[cells$regime == "C3", ][1L, ], 6L, 6L)
    splitCompare(cells[cells$regime == "C1", ][1L, ], 6L, 0L)
    splitCompare(cells[cells$family == "c4:500", ][1L, ], 4L, 4L)
})

test_that("seeds depend on the cell and the index only", {
    cells <- valCoverageCells()
    a <- cells[cells$regime == "C3", ][1L, ]
    expect_identical(valSimSeed(a, 7L), valSimSeed(a, 7L))
    expect_false(identical(valSimSeed(a, 7L), valSimSeed(a, 8L)))
    expect_false(identical(valSimSeed(a, 7L), valTargetSeed(a, 7L)))
})

test_that("a coverage plan splits each cell into fixed consecutive ranges", {
    cells <- valCoverageCells()
    plan <- valCoveragePlan(cells[cells$regime == "C3", ][1:2, ], ranges = c(2L, 1L), nsim = 10L,
                            target_sims = 7L)
    expect_identical(nrow(plan), 3L)
    a <- plan[plan$cell == plan$cell[1L], ]
    expect_identical(c(a$sim_from, a$sim_to), c(1L, 6L, 5L, 10L))
    expect_identical(c(a$target_from, a$target_to), c(1L, 5L, 4L, 7L))
})
