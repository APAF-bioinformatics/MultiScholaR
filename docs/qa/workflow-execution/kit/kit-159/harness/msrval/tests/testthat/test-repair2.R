# Second pass of the pre-release repair (audit register AUDIT-KIT-159-impl-prerelease-20261124;
# owner rulings of methods ledger VAL-109): O8 at every unit count, O5 (both labels), O4 (zero-miss
# conservative rows; the half-width ratio) and the I10 diagnostic generator (fixed per-cell feature
# draw, changes split +1 and -1).

oneVerdict <- function(...) {
    data.frame(member = NA_character_, metric = "mae", arm = "good", comparator = "null",
               outcome = "match", conf = 0.95, regime = "C2", coverage = 0.95, mc_se = 0.005,
               statistic = "paired_t@normal", family = NA_character_, m = NA_integer_, ...,
               stringsAsFactors = FALSE)
}

oneAttribution <- function() {
    data.frame(role = "corpus", kind = "data", licence = "generated", redistribution = "permitted",
               citation = "synthetic", stringsAsFactors = FALSE)
}

test_that("R-O8: over independent units nothing is widened at any count; report only, stated", {
    decl <- repairDeclaration()
    decl$metrics[[1L]]$within_run <- list(B = 2000L, seed = 159103L)
    for (n in c(5L, 8L)) {
        u <- repairUnits(n)
        v <- 0.1 + 0.001 * seq_len(n)
        metrics <- rbind(repairMetrics("mae", "good", v),
                         repairMetrics("mae", "null", v - 0.002 + 0.0005 * (-1)^seq_len(n)))
        within <- data.frame(arm = "good", comparator = "null", metric = "mae",
                             scope = "within_run", view = "intersection(good+null)",
                             ci_low = -0.5, ci_high = 0.5,
                             flags = sprintf("unit=%s;n=20;status=ok", u$unit),
                             stringsAsFactors = FALSE)
        cmp <- valRunLevel(decl, metrics, u, within = within,
                           regimes = regimeRows("paired_t", n, class = "accuracy"))
        run <- cmp[cmp$scope %in% "run_level", ]
        expect_lt(run$ci_high[[1L]], 0.5)
        expect_identical(lengths(regmatches(run$flags[[1L]],
                                            gregexpr("within_run_reported=", run$flags[[1L]]))),
                         n)
        want <- if (n < 6L) "below_6" else "6_or_more_open"
        expect_match(run$flags[[1L]], sprintf("within_run_widening=report_only(O8:%s)", want),
                     fixed = TRUE)
        vr <- attr(cmp, "verdicts")
        expect_identical(vr$within_note[[1L]], sprintf("report_only(O8:%s)", want))
    }
    path <- withr::local_tempfile(fileext = ".md")
    valWriteResults(path, oneVerdict(within_note = "report_only(O8:6_or_more_open)"),
                    oneAttribution())
    txt <- readLines(path)
    expect_true(any(grepl(paste("within-run intervals reported beside the verdict, not widened:",
                                "at 6 or more independent units the widening rule is open (O8)"),
                          txt, fixed = TRUE)))
})

test_that("R-O5: T2 rows reading pass under T2's band and conservative under r1.13's; both printed", {
    reg <- msrval:::valHarnessRegimes()
    expect_true("status_r113" %in% names(reg))
    two <- reg$status == "pass" & reg$status_r113 == "conservative"
    # 117 in T2; the 19 C4 rows among them were re-pinned by the I10 diagnostic (VAL-109 (2)),
    # classified under r1.13's rule, so 98 remain
    expect_identical(sum(two %in% TRUE), 98L)
    tab <- table(paste(reg$regime[two %in% TRUE], reg$class[two %in% TRUE]))
    expect_identical(as.integer(tab[c("C3 accuracy", "C3 threshold_count")]), c(45L, 53L))
    same <- !(two %in% TRUE)
    expect_identical(reg$status_r113[same], reg$status[same])
    path <- withr::local_tempfile(fileext = ".md")
    valWriteResults(path, oneVerdict(regime_status = "pass", regime_status_r113 = "conservative"),
                    oneAttribution())
    expect_true(any(grepl("[pass under T2's band; conservative under r1.13's miss-rate rule]",
                          readLines(path), fixed = TRUE)))
})

test_that("R-O4: a conservative row with no observed miss licenses only the statistic it simulated", {
    reg <- data.frame(regime = "C4e", method = "jackknife", class = "rank", form = "block",
                      n = 14L, level = 0.95, status = "conservative", coverage = 1, mc_se = 0,
                      source = "harness", statistic = "pauc_diff@matrix4840_block",
                      coverage_min = 1, coverage_max = 1, stringsAsFactors = FALSE)
    iv <- valNewInterval(0.01, 0.03, "jackknife", n_units = 14L, conf = 0.95, form = "block",
                         class = "rank")
    expect_error(valClassify(iv, 0.05, TRUE, regimes = reg), "O4", class = "val_spec_invalid")
    iv$statistic <- "rejection_share"
    expect_error(valClassify(iv, 0.05, TRUE, regimes = reg), "O4", class = "val_spec_invalid")
    iv$statistic <- "pauc_diff"
    expect_identical(valClassify(iv, 0.05, TRUE, regimes = reg)$regime, "C4e")
    one_miss <- reg
    one_miss$coverage <- one_miss$coverage_min <- 0.9995
    iv$statistic <- NULL
    expect_identical(valClassify(iv, 0.05, TRUE, regimes = one_miss)$regime, "C4e")
    # The declared statistic reaches the interval in valRunLevel.
    decl <- repairDeclaration()
    decl$metrics[[1L]]$interval <- list(method = "paired_t", statistic = "mae")
    zr <- regimeRows("paired_t", 5L, class = "accuracy", status = "conservative")
    zr$coverage <- 1
    zr$coverage_min <- 1
    zr$statistic <- "mae@normal"
    v <- c(0.10, 0.11, 0.09, 0.10, 0.12)
    metrics <- rbind(repairMetrics("mae", "good", v), repairMetrics("mae", "null", v - 0.001))
    cmp <- valRunLevel(decl, metrics, repairUnits(5L), regimes = zr)
    expect_identical(cmp$regime[[1L]], "test")
})

test_that("R-O4: the half-width ratio is mean half-width over the statistic's SD times z", {
    cell <- list(cell = "C4-143", regime = "C4", family = "c4:500", n = 4L,
                 condition = "matrix500")
    lv <- msrval:::COVERAGE_LEVELS
    est <- rep(c(-1, 1), 1000L)
    z <- stats::qnorm(1 - (1 - lv) / 2)
    hw <- matrix(rep(z * stats::sd(est), each = 2000L), 2000L)
    hits <- matrix(TRUE, 2000L, length(lv))
    rows <- msrval:::valCoverageRows(cell, "mae", "accuracy", hits, "approx_run", target = 0,
                                     half_width = hw, estimate = est)
    expect_equal(rows$half_width_ratio, rep(1, length(lv)), tolerance = 1e-12)
    rows0 <- msrval:::valCoverageRows(cell, "mae", "accuracy", hits, "approx_run", target = 0)
    expect_true(all(is.na(rows0$half_width_ratio)))
    cells <- valCoverageCells()
    part <- valCoverageRange(cells[cells$cell == "C4-143", ], 1:2, 1:2)
    expect_identical(dim(part$estimate), c(2L, 6L))
    expect_true(all(is.finite(part$estimate)))
    det <- valCombineCoverage(cells[cells$cell == "C4-143", ], list(part), 2L, 2L)
    expect_true("half_width_ratio" %in% names(det))
    expect_true("half_width_ratio" %in% names(valSummariseCoverage(det)))
})

test_that("R-I10: the redraw generator is unchanged; the fixed generator draws features once per cell", {
    lay <- msrval:::valCovLayout("replicate", 4L)
    got <- msrval:::valWithSeed(7L, msrval:::valCovMatrixData(lay, 500L, "replicate"))
    old <- msrval:::valWithSeed(7L, {
        changed <- seq_len(500L) <= 100L
        gi <- match(lay$group, c("A", "B"))
        x <- matrix(stats::rnorm(500L, 20, 1), 500L, nrow(lay)) +
            outer(changed, c(0, 1)[gi]) +
            matrix(stats::rnorm(nrow(lay), 0, 0.1), 500L, nrow(lay), byrow = TRUE) +
            matrix(stats::rnorm(500L * nrow(lay), 0, 0.3), 500L, nrow(lay))
        list(x = x, truth = outer(changed, 1))
    })
    expect_identical(got, old)
    cells <- valCoverageCells()
    cell <- as.list(cells[cells$cell == "C4-143", ])
    cell$generator <- "fixed_symmetric"
    f1 <- msrval:::valCovFixedFeatures(cell, 500L)
    expect_identical(f1, msrval:::valCovFixedFeatures(cell, 500L))
    expect_identical(length(f1$baseline), 500L)
    expect_identical(as.integer(table(f1$sign[seq_len(100L)])), c(50L, 50L))
    expect_true(all(f1$sign[-seq_len(100L)] == 0))
    d1 <- msrval:::valWithSeed(1L, msrval:::valCovMatrixData(lay, 500L, "replicate", fixed = f1))
    d2 <- msrval:::valWithSeed(2L, msrval:::valCovMatrixData(lay, 500L, "replicate", fixed = f1))
    expect_identical(d1$truth, d2$truth)
    expect_identical(as.integer(table(d1$truth[, 1L])), c(50L, 400L, 50L))
    # the feature baselines are common: rows differ between draws only by noise of SD < 0.5
    expect_lt(stats::sd(rowMeans(d1$x) - rowMeans(d2$x)), 0.5)
    # split and unsplit ranges agree under the fixed generator
    whole <- valCoverageRange(cell, 1:2, 1:2)
    a <- valCoverageRange(cell, 1L, 1L)
    b <- valCoverageRange(cell, 2L, 2L)
    expect_identical(whole$lo[2L, , ], b$lo[1L, , ])
    expect_identical(whole$target_values, rbind(a$target_values, b$target_values))
    redraw <- valCoverageRange(cells[cells$cell == "C4-143", ], 1L, integer())
    expect_false(identical(redraw$lo, a$lo))
})

test_that("R-I10: under the fixed symmetric generator the plain arm's null rejection share is near 0.05", {
    cells <- valCoverageCells()
    cell <- as.list(cells[cells$cell == "C4e-147", ])
    expect_identical(cell$n, 5L)
    cell$generator <- "fixed_symmetric"
    spec <- msrval:::valCovSpec(cell)
    rej <- vapply(seq_len(100L), function(s) {
        msrval:::valWithSeed(1600000000L + s, spec$target_one())[[3L]]
    }, numeric(1))
    se <- stats::sd(rej) / sqrt(length(rej))
    expect_lte(abs(mean(rej) - 0.05), 3 * se + 1e-9)
})

test_that("R-I10: a plan's generator column reaches the cell in the chunk tasks and the detail", {
    dir <- withr::local_tempdir()
    cells <- valCoverageCells()
    plan <- valCoveragePlan(cells[cells$cell == "C4-143", ], 1L, nsim = 1L, target_sims = 1L)
    plan$generator <- "fixed_symmetric"
    pf <- file.path(dir, "plan.tsv")
    utils::write.table(plan, pf, sep = "\t", quote = FALSE, row.names = FALSE)
    valCoverageTasks(pf, plan$task, dir)
    got <- readRDS(file.path(dir, paste0(plan$task, ".rds")))
    cell <- as.list(cells[cells$cell == "C4-143", ])
    cell$generator <- "fixed_symmetric"
    expect_identical(got$lo, valCoverageRange(cell, 1L, 1L)$lo)
    # R3 showed the two lines above passing on the tests-only commit (both sides ignored the
    # generator), so the task must also differ from the redraw generator's and name it
    redraw <- valCoverageRange(cells[cells$cell == "C4-143", ], 1L, 1L)
    expect_false(identical(got$lo, redraw$lo))
    det <- valCombineCoverage(cell, list(got), 1L, 1L)
    expect_true("generator" %in% names(det))
    expect_true(all(det$generator == "fixed_symmetric"))
})
