# DI10K (methods ledger VAL-163 (3); coordinator's ruling of 2026-11-25): the confirmatory cell in
# KIT-162's composition. One cell outside T-COV's catalogue, with its own index (so its own seeds):
# 4,840 features, 40 changed, all up at effect scale 1 (group effects 0, 0.5, 1, 1.5, 2), C4e's
# five-group block layout at g = 5, and C4e's statistics `mae`, `bias`, `rej` (not the rank
# statistic). T-COV's catalogue, its seeds and the pinned coverage table stay as they are.
# `cov-cells-base.tsv` is the catalogue at staging 8cb536cc (169 cells), written out from
# `valCoverageCells()`'s loops; the first test shows it equal to the code before and after.

dkCell <- function(generator = "fixed_onedir") {
    conf <- valCoverageCells(confirmatory = TRUE)
    cell <- as.list(conf[conf$family == "c4e:4840acc", ])
    cell$generator <- generator
    cell
}

test_that("DI10K: T-COV's catalogue, its seeds and the pinned coverage table are unchanged", {
    base <- utils::read.delim(test_path("cov-cells-base.tsv"), stringsAsFactors = FALSE)
    cells <- valCoverageCells()
    expect_identical(nrow(cells), 169L)
    expect_identical(as.list(cells[, names(base)]), as.list(base))
    idx <- as.integer(sub("^.*-", "", base$cell))
    expect_identical(vapply(seq_len(nrow(cells)), function(i) valSimSeed(cells[i, ], 3L), 0L),
                     1000000000L + idx * 100000L + 3L)
    expect_identical(vapply(seq_len(nrow(cells)), function(i) valTargetSeed(cells[i, ], 3L), 0L),
                     1500000000L + idx * 100000L + 3L)
    path <- system.file("extdata", "coverage-table.tsv", package = "msrval")
    expect_identical(digest::digest(file = path, algo = "sha256"),
                     "894e42e5033386d4ed82a9292900aca78c3903e78342997efc464f0538c82db1")
})

test_that("DI10K: the confirmatory cell has its own index; T-COV's catalogue does not hold it", {
    cells <- valCoverageCells()
    conf <- valCoverageCells(confirmatory = TRUE)
    expect_identical(nrow(conf), nrow(cells) + 1L)
    expect_identical(as.list(conf[seq_len(nrow(cells)), ]), as.list(cells))
    k <- conf[nrow(conf), ]
    expect_identical(k$cell, "C4e-170")
    expect_identical(k$regime, "C4e")
    expect_identical(k$family, "c4e:4840acc")
    expect_identical(k$n, 5L)
    expect_false(any(cells$family == "c4e:4840acc"))
    expect_false(msrval:::valCellIndex(k) %in% msrval:::valCellIndex(cells))
})

test_that("DI10K: 4,840 features, exactly 40 changed, all up at scale 1, C4e's g = 5 layout", {
    cell <- dkCell()
    f <- msrval:::valCovFixedFeatures(cell, 4840L)
    expect_identical(length(f$baseline), 4840L)
    expect_identical(f$sign, c(rep(1, 40L), rep(0, 4800L)))
    lay <- msrval:::valCovLayout("block", 5L)
    expect_identical(nrow(lay), 50L)
    expect_identical(as.vector(table(lay$group, lay$unit)), rep(2L, 25L))
    spec <- msrval:::valCovSpec(cell)
    sim <- msrval:::valWithSeed(11L, spec$sim(1L))
    want <- msrval:::valWithSeed(11L, msrval:::valCovMatrixData(lay, 4840L, "block", 1,
                                                                fixed = f))
    # the cell's simulated matrix is the generator's at scale 1 (not the rank-calibrated scale)
    expect_identical(sim$data, want$x)
    expect_identical(dim(sim$data), c(4840L, 50L))
    expect_identical(sum(rowSums(want$truth != 0) > 0), 40L)
    expect_true(all(want$truth[1:40, ] == matrix(c(0.5, 1, 1.5, 2), 40L, 4L, byrow = TRUE)))
    expect_true(all(want$truth[-(1:40), ] == 0))
})

test_that("DI10K: the cell computes C4e's mae, bias and rej on the simulated matrix", {
    cell <- dkCell()
    spec <- msrval:::valCovSpec(cell)
    expect_identical(spec$stats$statistic, c("mae", "bias", "rej"))
    expect_identical(spec$stats$class, c("accuracy", "accuracy", "threshold_count"))
    expect_null(spec$range)
    expect_true(is.na(spec$effect))
    part <- valCoverageRange(cell, 1L, 1L)
    expect_identical(dim(part$estimate), c(1L, 3L))
    expect_true(all(is.finite(part$target_values)))
    lay <- msrval:::valCovLayout("block", 5L)
    f <- msrval:::valCovFixedFeatures(cell, 4840L)
    dat <- msrval:::valWithSeed(valSimSeed(cell, 1L),
                                msrval:::valCovMatrixData(lay, 4840L, "block", 1, fixed = f))
    fit <- msrval:::valCovFit(dat$x, lay, "block")
    pos <- dat$truth != 0
    expect_identical(sum(!pos), 19200L)
    expect_equal(part$estimate[1L, ],
                 c(stats::median(abs(fit$est - dat$truth)), stats::median(fit$est - dat$truth),
                   mean(fit$p[!pos] < 0.05)))
})

test_that("DI10K: a plan naming the cell runs it with the plan's generator, in tasks and runs", {
    dir <- withr::local_tempdir()
    withr::local_envvar(RUNNER_TEMP = dir)
    conf <- valCoverageCells(confirmatory = TRUE)
    plan <- valCoveragePlan(conf[conf$family == "c4e:4840acc", ], 1L, nsim = 1L, target_sims = 1L,
                            request = "DK")
    plan$est_seconds <- 1
    plan$generator <- "fixed_onedir"
    pf <- file.path(dir, "plan.tsv")
    utils::write.table(plan, pf, sep = "\t", quote = FALSE, row.names = FALSE)
    valCoverageTasks(pf, plan$task, dir)
    got <- readRDS(file.path(dir, paste0(plan$task, ".rds")))
    expect_identical(got$lo, valCoverageRange(dkCell(), 1L, 1L)$lo)
    redraw <- valCoverageRange(dkCell(NULL), 1L, 1L)
    expect_false(identical(got$lo, redraw$lo))
    out <- file.path(dir, "out")
    dir.create(out)
    valRunCoverage(pf, request = "DK", nsim = 1L, target_sims = 1L, chunks = 1L, out_dir = out)
    det <- utils::read.delim(file.path(out, "coverage-detail.tsv"), stringsAsFactors = FALSE)
    expect_identical(unique(det$cell), "C4e-170")
    expect_identical(unique(det$statistic), c("mae", "bias", "rej"))
    expect_true(all(det$generator == "fixed_onedir"))
})
