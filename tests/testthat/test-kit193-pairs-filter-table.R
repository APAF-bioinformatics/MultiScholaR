# KIT-193 target 1 (ticket sections 6.1, 5.5), second file: 6.1 (ii) and (v) on every row of the Lean oracle table
# tests/testdata/kit193/filter-thresholds.tsv, and gaps 2 to 5 of AUDIT-KIT-187-impl I10 (the full helper on a sample of
# rows, unequal groups included; two-group designs' stored P; the grouping column `group`; the group level through
# the helper). Header, bound code and the other parts: test-kit193-pairs-filter.R; helpers: helper-kit193-filter.R.

test_that("KIT-193 T1 pair (ii), (v): the decision and the peptide floor equal the model on every row (n* <= 60)", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    expect_identical(nrow(table), 73810L)
    expect_true(all(table$rounded_ok == 1L))
    # one comparison per (n*, m) over all its rows at once: every (n, k) with 0 <= k <= n, as 6.1 (ii) and (v) state
    for (key in split(seq_len(nrow(table)), list(table$nstar, table$m), drop = TRUE)) {
        rows <- table[key, , drop = FALSE]
        cutoff <- kit193Cutoffs(rows$nstar[[1L]], rows$m[[1L]], 1L)$groupwise_percentage_cutoff
        n <- rep(rows$n, rows$n + 1L)
        k <- unlist(lapply(rows$n, function(size) 0:size))
        kmax <- rep(rows$kmax, rows$n + 1L)
        info <- sprintf("n* = %d, m = %d", rows$nstar[[1L]], rows$m[[1L]])
        expect_identical(exceedsRoundedPercent(k, n, cutoff), k > kmax, info = info)
        expect_identical(floor(rows$n * roundedPercentNumerator(cutoff) / 200000), as.numeric(rows$kmax), info = info)
        expect_identical(allowedByRoundedPercent(rows$n, cutoff), as.numeric(rows$kmax), info = info)
    }
})

test_that("KIT-193 T1 pair (I10 gap 3): on two-group designs the stored cutoff's P is the larger group's", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    # every larger group b up to 60 and every m; the smaller group a in {1, b %/% 2, b - 1} (a sample of a < b)
    for (b in 2:60) {
        for (a in unique(c(1L, b %/% 2L, b - 1L))) {
            for (m in seq_len(b)) {
                P <- round(1000 * kit193Cutoffs(c(a, b), m, 1L)$groupwise_percentage_cutoff)
                expect_identical(P, as.numeric(kit193FilterRow(table, b, m, 1L)$P),
                    info = sprintf("a = %d, b = %d, m = %d", a, b, m))
            }
        }
    }
})

test_that("KIT-193 T1 pair (I10 gaps 2, 4): the full helper keeps exactly the rows the model keeps", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    # group G1 the larger (n* runs, it sets the cutoff), group G2 the one under test (n runs, k of them missing);
    # min_groups = 2, so a protein failing one group is removed and row kj is kept exactly when G2 passes
    for (nstar in c(3L, 7L, 20L, 33L, 60L)) {
        for (m in unique(c(1L, as.integer(ceiling(nstar / 2)), nstar))) {
            for (n in unique(c(1L, as.integer(ceiling(nstar / 3)), nstar - 1L, nstar))) {
                kmax <- kit193FilterRow(table, nstar, m, n)$kmax
                design <- kit193Design(c(nstar, n))
                cutoffs <- kit193Cutoffs(c(nstar, n), m, 2L)
                kept <- kit193Kept(kit193Table(design, kit193Missing("G2", n)), design, cutoffs)
                expect_setequal(kept, paste0("k", 0:kmax))
            }
        }
    }
})

test_that("KIT-193 T1 pair (I10 gap 5): the group level through the helper keeps exactly the model's rows", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    # G groups of one run each, one required per group (a group fails when its run is missing), g groups required;
    # row fj has the runs of groups 1 to j missing: the (G, g) instance of the same theorem, n* = n = G
    for (G in c(1:20, 40L)) {
        design <- kit193Design(rep(1L, G))
        missing <- stats::setNames(lapply(0:G, function(j) if (j) paste0("G", seq_len(j), "_1") else character()),
            paste0("f", 0:G))
        for (g in seq_len(G)) {
            kmax <- kit193FilterRow(table, G, g, G)$kmax
            cutoffs <- kit193Cutoffs(rep(1L, G), 1L, g)
            expect_identical(cutoffs$groupwise_percentage_cutoff, 0)
            kept <- kit193Kept(kit193Table(design, missing), design, cutoffs)
            expect_setequal(kept, paste0("f", 0:kmax))
        }
    }
})
