# The repair after its re-audit (register AUDIT-KIT-159-impl-recheck-20261124; owner rulings of
# methods ledger VAL-137): the family bound's one convention (VAL-137 (2)), the jackknife widening
# only in forms (ii) and (iii) (VAL-137 (3)), the SE inflation beside the half-width ratio (E8),
# each regime's generator (E7), per-unit truth hashes in provenance (E4), and the one-directional
# fixed generator of the confirmatory cell DI10b.

r4Verdict <- function(...) {
    data.frame(member = NA_character_, metric = "mae", arm = "good", comparator = "null",
               outcome = "match", conf = 0.95, regime = "C4", coverage = 0.981, mc_se = 0.003,
               statistic = "mae@matrix500", family = NA_character_, m = NA_integer_, ...,
               stringsAsFactors = FALSE)
}

r4Attribution <- function() {
    data.frame(role = "corpus", kind = "data", licence = "generated", redistribution = "permitted",
               citation = "synthetic", stringsAsFactors = FALSE)
}

test_that("R-V137-2: the family bound sums each component's one-sided 95 % Clopper-Pearson upper bound", {
    cmp <- r4Verdict(family = "fam")[rep(1L, 3L), ]
    cmp$metric <- c("a", "b", "c")
    cmp$m <- 3L
    cmp$conf <- 0.975
    cmp$coverage <- c(0.972, 0.968, 1)
    path <- withr::local_tempfile(fileext = ".md")
    valWriteResults(path, cmp, r4Attribution())
    txt <- readLines(path)
    x <- c(56, 64, 0)
    want <- sum(stats::qbeta(0.95, x + 1, 2000 - x))
    line <- grep("^- fam: m = 3", txt, value = TRUE)
    expect_length(line, 1L)
    expect_match(line, paste("family error bound from measured misses", format(round(want, 4))),
                 fixed = TRUE)
    expect_match(line, "one-sided 95 % Clopper-Pearson upper bound", fixed = TRUE)
    expect_false(any(grepl("3/simulations", txt, fixed = TRUE)))
    # a regime's own simulation count is used where the row carries it
    cmp$nsim <- c(200, 2000, 2000)
    cmp$coverage <- c(0.97, 0.968, 1)
    valWriteResults(path, cmp, r4Attribution())
    x2 <- c(6, 64, 0)
    want2 <- sum(stats::qbeta(0.95, x2 + 1, c(200, 2000, 2000) - x2))
    expect_true(any(grepl(paste("family error bound from measured misses", format(round(want2, 4))),
                          readLines(path), fixed = TRUE)))
})

test_that("R-V137-3: the jackknife widens only in forms (ii) and (iii); form (i) reports only", {
    ju <- data.frame(unit = sprintf("u%02d", 1:5), unit_role = c("full", rep("deletion", 4)),
                     stringsAsFactors = FALSE)
    jm <- rbind(repairMetrics("mae", "good", c(0.10, 0.11, 0.12, 0.10, 0.11)),
                repairMetrics("mae", "null", rep(0.1, 5)))
    within <- data.frame(arm = "good", comparator = "null", metric = "mae", scope = "within_run",
                         view = "intersection(good+null)", ci_low = -0.5, ci_high = 0.5,
                         flags = "unit=u01;n=20;status=ok", stringsAsFactors = FALSE)
    for (form in c("unit", "replicate", "block")) {
        jk <- repairDeclaration()
        jk$metrics[[1L]]$interval <- list(method = "jackknife", form = form)
        cmp <- valRunLevel(jk, jm, ju, within = within,
                           regimes = regimeRows("jackknife", 4L, class = "accuracy", form = form))
        run <- cmp[cmp$scope %in% "run_level", ]
        expect_identical(nrow(run), 1L)
        if (form == "unit") {
            expect_lt(run$ci_high[[1L]], 0.5)
            expect_false(grepl("widened_to_within_run", run$flags[[1L]], fixed = TRUE))
            expect_match(run$flags[[1L]], "within_run_reported=u01:[-0.5,0.5]", fixed = TRUE)
        } else {
            expect_equal(c(run$ci_low[[1L]], run$ci_high[[1L]]), c(-0.5, 0.5))
            expect_match(run$flags[[1L]], "widened_to_within_run", fixed = TRUE)
        }
    }
    # the exported widening applies to forms (ii) and (iii) only
    iv <- valNewInterval(-0.01, 0.02, "jackknife", n_units = 4, conf = 0.975, form = "unit")
    w <- valWidenFromRequests(iv, within, "good", "null", "mae")
    expect_identical(c(w$ci_low, w$ci_high), c(-0.01, 0.02))
    expect_false("widened_to_within_run" %in% w$flags)
    expect_true(any(grepl("^within_run_reported", w$flags)))
})

test_that("R-E8: the SE inflation is the half-width ratio with the t factor removed", {
    cell <- list(cell = "C4-143", regime = "C4", family = "c4:500", n = 4L,
                 condition = "matrix500")
    lv <- msrval:::COVERAGE_LEVELS
    est <- rep(c(-1, 1), 1000L)
    tq <- stats::qt(1 - (1 - lv) / 2, 3)
    hw <- matrix(rep(tq * stats::sd(est), each = 2000L), 2000L)
    hits <- matrix(TRUE, 2000L, length(lv))
    rows <- msrval:::valCoverageRows(cell, "mae", "accuracy", hits, "approx_run", target = 0,
                                     half_width = hw, estimate = est)
    # an exactly estimated SE: the ratio is t(g - 1) / z, not 1; the SE inflation is 1
    expect_equal(rows$half_width_ratio, tq / stats::qnorm(1 - (1 - lv) / 2), tolerance = 1e-12)
    expect_equal(rows$se_inflation, rep(1, length(lv)), tolerance = 1e-12)
    rows0 <- msrval:::valCoverageRows(cell, "mae", "accuracy", hits, "approx_run", target = 0)
    expect_true(all(is.na(rows0$se_inflation)))
    cells <- valCoverageCells()
    part <- valCoverageRange(cells[cells$cell == "C4-143", ], 1:2, 1:2)
    det <- valCombineCoverage(cells[cells$cell == "C4-143", ], list(part), 2L, 2L)
    expect_true("se_inflation" %in% names(valSummariseCoverage(det)))
    # the pinned table: every measured row carries both, related by z / t(g - 1)
    reg <- msrval:::valHarnessRegimes()
    expect_true("se_inflation" %in% names(reg))
    ratio <- suppressWarnings(as.numeric(reg$half_width_ratio))
    m <- is.finite(ratio)
    expect_identical(sum(m), 64L)
    z <- stats::qnorm(1 - (1 - reg$level[m]) / 2)
    tt <- stats::qt(1 - (1 - reg$level[m]) / 2, reg$n[m] - 1)
    expect_equal(as.numeric(reg$se_inflation[m]), ratio[m] * z / tt, tolerance = 1e-9)
    expect_true(all(reg$se_inflation[!m] == "not_recorded"))
})

test_that("R-E7: the pinned table names each row's generator; the classifier and results.md carry it", {
    reg <- msrval:::valHarnessRegimes()
    expect_true("generator" %in% names(reg))
    mat <- reg$regime %in% c("C4", "C4e")
    expect_true(all(reg$generator[mat & reg$class != "rank"] == "fixed_symmetric"))
    expect_true(all(reg$generator[mat & reg$class == "rank"] == "redraw"))
    expect_true(all(reg$generator[!mat] == "redraw"))
    iv <- valNewInterval(-0.01, 0.01, "jackknife", n_units = 4L, conf = 0.95, form = "replicate",
                         class = "accuracy")
    expect_identical(valClassify(iv, 0.05, TRUE)$generator, "fixed_symmetric")
    decl <- repairDeclaration()
    v <- c(0.10, 0.11, 0.09, 0.10, 0.12)
    metrics <- rbind(repairMetrics("mae", "good", v), repairMetrics("mae", "null", v - 0.001))
    rg <- regimeRows("paired_t", 5L, class = "accuracy")
    rg$generator <- "redraw"
    cmp <- valRunLevel(decl, metrics, repairUnits(5L), regimes = rg)
    expect_identical(attr(cmp, "verdicts")$generator[[1L]], "redraw")
    path <- withr::local_tempfile(fileext = ".md")
    valWriteResults(path, rbind(r4Verdict(generator = "fixed_symmetric"),
                                r4Verdict(generator = "redraw")), r4Attribution())
    txt <- readLines(path)
    expect_true(any(grepl("generator fixed_symmetric", txt, fixed = TRUE)))
    expect_true(any(grepl("conditional on one feature draw", txt, fixed = TRUE)))
    expect_true(any(grepl("generator redraw", txt, fixed = TRUE)))
})

test_that("R-E4: per-unit truths are hashed into provenance", {
    dir <- withr::local_tempdir()
    code <- writeRepairTicket(dir)
    decl <- repairDeclaration()
    out <- valRunPhase(phase = "counted", decl = decl, dataset = "d1", units = repairUnits(),
                       label = "toy", code_files = code, data_dir = dir,
                       work_dir = file.path(dir, "work"), out_dir = file.path(dir, "val-out"),
                       chunks = 2L, regimes = regimeRows("paired_t", 5L, class = "accuracy"),
                       manifest = repairManifest())
    prov <- jsonlite::read_json(file.path(out, "provenance.json"))
    pu <- prov$truth$per_unit
    expect_setequal(names(pu), sprintf("u%02d", 1:5))
    for (u in names(pu)) {
        expect_match(pu[[u]]$primary$sha256, "^[0-9a-f]{64}$")
        expect_identical(pu[[u]]$primary$rows, 20L)
    }
    expect_true(valCheckProvenance(file.path(out, "provenance.json")))
})

test_that("R-DI10b: the one-directional fixed generator keeps the fixed draw and moves every changed feature up", {
    cells <- valCoverageCells()
    cell <- as.list(cells[cells$cell == "C4e-147", ])
    cell$generator <- "fixed_onedir"
    sym <- cell
    sym$generator <- "fixed_symmetric"
    f <- msrval:::valCovFixedFeatures(cell, 500L)
    fs <- msrval:::valCovFixedFeatures(sym, 500L)
    expect_identical(f$baseline, fs$baseline)
    expect_true(all(f$sign[seq_len(100L)] == 1))
    expect_true(all(f$sign[-seq_len(100L)] == 0))
    # the generator reaches the simulation: median normalisation is shifted as under T2's
    # one-directional changes (target bias -0.21 at g = 5), not centred as under the split
    spec <- msrval:::valCovSpec(cell)
    b <- vapply(seq_len(20L), function(s) {
        msrval:::valWithSeed(1600000000L + s, spec$target_one())[[2L]]
    }, numeric(1))
    expect_lt(mean(b), -0.1)
    lay <- msrval:::valCovLayout("block", 5L)
    d1 <- msrval:::valWithSeed(1L, msrval:::valCovMatrixData(lay, 500L, "block", fixed = f))
    expect_true(all(d1$truth[seq_len(100L), ] > 0))
})
