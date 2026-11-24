# P1K (methods ledger VAL-176; coordinator's ruling of 2026-11-25): the confirmatory rank cell in
# the form of KIT-162's primary endpoint P1. One cell after DI10K's, with its own index (C4e-171, family
# `c4e:4840p1`): 4,840 features, 40 changed, all up (`fixed_onedir` from the plan), C4e's five-group
# block layout at g = 5, and `pauc_diff` (normalised pAUC to FPR 0.01, moderated minus ordinary) on
# ONE contrast (group C against A: 40 positives, 4,800 nulls). The effect scale is the rank cells'
# rule (`valRankEffect`: the grid value whose mean plain-arm pAUC lies in [0.4, 0.8], nearest 0.6)
# applied to that one contrast. T-COV's 169 cells, DI10K's C4e-170, their seeds and the pinned
# coverage table stay as they are.

p1Cell <- function(generator = "fixed_onedir") {
    conf <- valCoverageCells(confirmatory = TRUE)
    cell <- as.list(conf[conf$family == "c4e:4840p1", ])
    cell$generator <- generator
    cell
}

paucOf <- function(truth, est, p) {
    s <- data.frame(key = seq_along(truth), expected = as.vector(truth),
                    estimate = as.vector(est), p = as.vector(p), s_own = as.vector(p),
                    truth_status = "scored", key_status = "ok", stringsAsFactors = FALSE)
    valPartialAuc(s, fpr_max = 0.01)
}

test_that("P1K: the 170 earlier cells, their seeds and the pinned coverage table are unchanged", {
    base <- utils::read.delim(test_path("cov-cells-base.tsv"), stringsAsFactors = FALSE)
    base170 <- rbind(base, data.frame(regime = "C4e", family = "c4e:4840acc", n = 5L,
                                      condition = "matrix4840acc_block", cell = "C4e-170",
                                      stringsAsFactors = FALSE))
    cells <- valCoverageCells()
    expect_identical(as.list(cells[, names(base)]), as.list(base))
    conf <- valCoverageCells(confirmatory = TRUE)
    expect_identical(as.list(conf[seq_len(170L), names(base170)]), as.list(base170))
    idx <- as.integer(sub("^.*-", "", base170$cell))
    expect_identical(vapply(seq_len(170L), function(i) valSimSeed(conf[i, ], 3L), 0L),
                     1000000000L + idx * 100000L + 3L)
    expect_identical(vapply(seq_len(170L), function(i) valTargetSeed(conf[i, ], 3L), 0L),
                     1500000000L + idx * 100000L + 3L)
    path <- system.file("extdata", "coverage-table.tsv", package = "msrval")
    expect_identical(digest::digest(file = path, algo = "sha256"),
                     "894e42e5033386d4ed82a9292900aca78c3903e78342997efc464f0538c82db1")
})

test_that("P1K: the P1 cell has its own index after DI10K's; T-COV's catalogue does not hold it", {
    conf <- valCoverageCells(confirmatory = TRUE)
    expect_identical(nrow(conf), 171L)
    k <- conf[171L, ]
    expect_identical(k$cell, "C4e-171")
    expect_identical(k$regime, "C4e")
    expect_identical(k$family, "c4e:4840p1")
    expect_identical(k$n, 5L)
    expect_false(any(valCoverageCells()$family == "c4e:4840p1"))
    expect_false(msrval:::valCellIndex(k) %in% msrval:::valCellIndex(conf[1:170, ]))
})

test_that("P1K: the effect is the rank cells' rule on the one contrast; the pooled rule unchanged", {
    pooled <- valRankEffect("block", 5L)
    # the pinned C4e g = 5 rank rows carry effect 0.35 (T2R)
    expect_equal(pooled$effect, 0.35)
    one <- valRankEffect("block", 5L, contrast = 2L)
    expect_true(any(abs(one$effect - seq(0.05, 1, by = 0.05)) < 1e-9))
    expect_gte(one$mean_pauc, 0.4)
    expect_lte(one$mean_pauc, 0.8)
    cat(sprintf("[p1k-calibration] one-contrast effect %.2f mean pAUC %.4f; pooled %.2f %.4f\n",
                one$effect, one$mean_pauc, pooled$effect, pooled$mean_pauc))
})

test_that("P1K: 40 positives of 4,840 on one contrast at the one-contrast scale; pauc_diff on it", {
    calls <- list()
    local_mocked_bindings(valRankEffect = function(form, g, nsim = 50L, contrast = NULL) {
        calls[[length(calls) + 1L]] <<- list(form = form, g = g, contrast = contrast)  # nolint
        list(effect = if (is.null(contrast)) 0.35 else 0.4, mean_pauc = NA_real_)
    }, .package = "msrval")
    cell <- p1Cell()
    spec <- msrval:::valCovSpec(cell)
    expect_identical(calls[[length(calls)]]$contrast, 2L)
    expect_identical(spec$stats$statistic, "pauc_diff")
    expect_identical(spec$stats$class, "rank")
    expect_equal(spec$effect, 0.4)
    expect_identical(spec$range, c(-1, 1))
    f <- msrval:::valCovFixedFeatures(cell, 4840L)
    expect_identical(f$sign, c(rep(1, 40L), rep(0, 4800L)))
    lay <- msrval:::valCovLayout("block", 5L)
    sim <- msrval:::valWithSeed(11L, spec$sim(1L))
    dat <- msrval:::valWithSeed(11L, msrval:::valCovMatrixData(lay, 4840L, "block", 0.4,
                                                               fixed = f))
    expect_identical(sim$data, dat$x)
    expect_equal(dat$truth[1:40, 2L], rep(0.4, 40L))
    fit <- msrval:::valCovFit(dat$x, lay, "block")
    a <- paucOf(dat$truth[, 2L], fit$est[, 2L], fit$p[, 2L])
    b <- paucOf(dat$truth[, 2L], fit$est[, 2L], fit$p_ord[, 2L])
    expect_identical(c(a$n_pos, a$n_null), c(40L, 4800L))
    expect_equal(sim$estimate, a$value - b$value)
})

test_that("P1K: the existing rank and C4e cells dispatch as before", {
    local_mocked_bindings(valRankEffect = function(form, g, nsim = 50L, contrast = NULL) {
        list(effect = if (is.null(contrast)) 0.35 else 0.4, mean_pauc = NA_real_)
    }, .package = "msrval")
    conf <- valCoverageCells(confirmatory = TRUE)
    c148 <- as.list(conf[conf$cell == "C4e-148", ])
    s148 <- msrval:::valCovSpec(c148)
    expect_identical(s148$stats$statistic, "pauc_diff")
    expect_equal(s148$effect, 0.35)
    lay <- msrval:::valCovLayout("block", 5L)
    sim <- msrval:::valWithSeed(5L, s148$sim(1L))
    dat <- msrval:::valWithSeed(5L, msrval:::valCovMatrixData(lay, 4840L, "block", 0.35))
    fit <- msrval:::valCovFit(dat$x, lay, "block")
    a <- paucOf(dat$truth, fit$est, fit$p)
    expect_identical(a$n_pos, 160L)
    expect_equal(sim$estimate, a$value - paucOf(dat$truth, fit$est, fit$p_ord)$value)
    for (id in c("C4e-147", "C4e-170")) {
        s <- msrval:::valCovSpec(as.list(conf[conf$cell == id, ]))
        expect_identical(s$stats$statistic, c("mae", "bias", "rej"))
        expect_true(is.na(s$effect))
    }
})

test_that("P1K: a plan naming the cell runs it with the plan's generator in the chunk tasks", {
    local_mocked_bindings(valRankEffect = function(form, g, nsim = 50L, contrast = NULL) {
        list(effect = if (is.null(contrast)) 0.35 else 0.4, mean_pauc = NA_real_)
    }, .package = "msrval")
    dir <- withr::local_tempdir()
    conf <- valCoverageCells(confirmatory = TRUE)
    plan <- valCoveragePlan(conf[conf$family == "c4e:4840p1", ], 1L, nsim = 1L, target_sims = 1L,
                            request = "PK")
    plan$generator <- "fixed_onedir"
    pf <- file.path(dir, "plan.tsv")
    utils::write.table(plan, pf, sep = "\t", quote = FALSE, row.names = FALSE)
    valCoverageTasks(pf, plan$task, dir)
    got <- readRDS(file.path(dir, paste0(plan$task, ".rds")))
    expect_identical(got$lo, valCoverageRange(p1Cell(), 1L, 1L)$lo)
    expect_false(identical(got$lo, valCoverageRange(p1Cell(NULL), 1L, 1L)$lo))
    det <- valCombineCoverage(p1Cell(), list(got), 1L, 1L)
    expect_identical(unique(det$cell), "C4e-171")
    expect_identical(unique(det$statistic), "pauc_diff")
    expect_true(all(det$generator == "fixed_onedir"))
    expect_equal(unique(det$effect), 0.4)
})

# R-P1 (harness auditor, P1K check): KIT-162 6.7 computes its stated false rates from C4e's
# one-sided misses (KIT-159 C66), so the detail records, per statistic and level, the misses at the
# lower bound (`miss_low`: the target below the interval) and at the upper bound (`miss_high`: the
# target above it), as counts of simulations. The regime table does not carry them, so the pinned
# table's rows are unchanged.
missesOf <- function(part, target, k) {
    lo <- matrix(part$lo[, , k], nrow = length(part$sims))
    hi <- matrix(part$hi[, , k], nrow = length(part$sims))
    list(low = as.integer(colSums((lo > target) %in% TRUE & !is.na(lo))),
         high = as.integer(colSums((hi < target) %in% TRUE & !is.na(hi))))
}

test_that("R-P1: the detail records the misses at the lower and at the upper bound per level", {
    cells <- valCoverageCells()
    cell <- cells[cells$cell == "C4-143", ]
    part <- valCoverageRange(cell, 1:6, 1:3)
    det <- valCombineCoverage(cell, list(part), 6L, 3L)
    expect_true(all(c("miss_low", "miss_high") %in% names(det)))
    for (st in unique(det$statistic)) {
        r <- det[det$statistic == st, ]
        k <- match(st, unique(det$statistic))
        want <- missesOf(part, r$target[[1L]], k)
        expect_identical(r$miss_low, want$low)
        expect_identical(r$miss_high, want$high)
        # without undefined intervals, a miss is on exactly one side
        expect_equal(r$miss_low + r$miss_high, round((1 - r$coverage) * 6))
    }
})

test_that("R-P1: the P1 cell's detail carries one-sided misses; the regime table does not", {
    local_mocked_bindings(valRankEffect = function(form, g, nsim = 50L, contrast = NULL) {
        list(effect = if (is.null(contrast)) 0.35 else 0.4, mean_pauc = NA_real_)
    }, .package = "msrval")
    cell <- p1Cell()
    part <- valCoverageRange(cell, 1:2, 1:2)
    det <- valCombineCoverage(cell, list(part), 2L, 2L)
    want <- missesOf(part, det$target[[1L]], 1L)
    expect_identical(det$miss_low, want$low)
    expect_identical(det$miss_high, want$high)
    tab <- valSummariseCoverage(det)
    expect_false(any(c("miss_low", "miss_high") %in% names(tab)))
    pinned <- utils::read.delim(system.file("extdata", "coverage-table.tsv", package = "msrval"),
                                nrows = 1L)
    expect_false(any(c("miss_low", "miss_high") %in% names(pinned)))
})
