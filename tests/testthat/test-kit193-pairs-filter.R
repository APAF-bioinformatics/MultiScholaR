# KIT-193 target 1 (ticket sections 6.1, 5.5): the shipped missing-value rule against the Lean model of
# proofs/lean/MultiScholaRProofs/Spec/FilterCount.lean. tests/testdata/kit193/filter-thresholds.tsv is written by Lean's
# interpreter from the Spec definitions (Oracle/FilterCount.lean: P, RoundedOK decided, kmax by the proved floor form);
# the record job regenerates it and compares bytes. Kind: package (bindings.tsv rows T1). Bound code: KIT-187 as merged
# (6416dbd6; the lines the Spec cites are unchanged from e6c42d1d to d84b6dfd, ticket 8.2).
# This file holds 6.1 (i), (iii), (iv), the planted defect and gap 1; test-kit193-pairs-filter-table.R holds
# (ii), (v) and gaps 2 to 5 (split so that each file stays within the ticket's 15 minutes; helpers in
# helper-kit193-filter.R). The two files together close the five gaps of AUDIT-KIT-187-impl I10:
#   1. the earlier comparison's option is unset in every test, and asserted unset;
#   2. the full helper `removeRowsWithMissingValuesPercentHelper` is driven on a sample of oracle rows, unequal group
#      sizes included, not only `exceedsRoundedPercent`;
#   3. on two-group designs the stored cutoff's P is the model's P of the larger group (the `max` of rounded values);
#   4. the grouping column is `group` in the design read by `updateMissingValueParameters` and in the helper;
#   5. the group level (G groups, g required) runs through the helper.

test_that("KIT-193 T1 pair (I10 gap 1): the earlier comparison's option is unset", {
    option <- kit193NoLegacy()
    expect_identical(option, "MultiScholaR.missing_filter_rule")
    expect_null(getOption(option))
})

test_that("KIT-193 T1 pair (i): every shipped cutoff satisfies RoundedOK, groups of 1 to 316 and groups up to 40", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    for (nstar in 1:316) {
        for (m in 1:nstar) {
            cutoff <- kit193Cutoffs(nstar, m, 1L)$groupwise_percentage_cutoff
            P <- round(1000 * cutoff)
            expect_true(abs(2 * nstar * P - 200000 * (nstar - min(nstar, m))) <= nstar,
                info = sprintf("n* = %d, m = %d, P = %s", nstar, m, format(P)))
        }
    }
    for (G in 1:40) {
        for (g in 1:G) {
            cutoff <- kit193Cutoffs(rep(1L, G), 1L, g)$max_groups_percentage_cutoff
            P <- round(1000 * cutoff)
            expect_true(abs(2 * G * P - 200000 * (G - g)) <= G, info = sprintf("G = %d, g = %d", G, g))
        }
    }
})

test_that("KIT-193 T1 pair (iii), (iv): the counterexample rows and KIT-187 T3's rows", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    cut_320 <- kit193Cutoffs(320L, 137L, 1L)$groupwise_percentage_cutoff
    P <- round(1000 * cut_320)
    expect_true(P %in% c(57187, 57188))
    # the shipped decision at (313, 179) is the model's for that P (either P is a valid record)
    expect_identical(exceedsRoundedPercent(179, 313, cut_320), 200000 * 179 > 313 * (2 * P + 1))
    cut_321 <- kit193Cutoffs(321L, 320L, 1L)$groupwise_percentage_cutoff
    expect_identical(round(1000 * cut_321), 312)
    expect_false(exceedsRoundedPercent(1, 320, cut_321))
    table <- kit193FilterTable()
    for (key in list(c(20, 9), c(3, 2))) {
        rows <- table[table$nstar == key[[1L]] & table$m == key[[2L]], , drop = FALSE]
        cutoff <- kit193Cutoffs(key[[1L]], key[[2L]], 1L)$groupwise_percentage_cutoff
        for (i in seq_len(nrow(rows))) {
            expect_identical(exceedsRoundedPercent(0:rows$n[[i]], rows$n[[i]], cutoff), 0:rows$n[[i]] > rows$kmax[[i]])
        }
    }
})

test_that("KIT-193 T1 pair, planted defect (ticket 5.5): the half step dropped fails at n* = 3, m = 2, n = 3, k = 1", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    row <- kit193FilterRow(table, 3L, 2L, 3L)
    cutoff <- kit193Cutoffs(3L, 2L, 1L)$groupwise_percentage_cutoff
    design <- kit193Design(c(3L, 3L))
    cutoffs <- kit193Cutoffs(c(3L, 3L), 2L, 2L)
    helper_table <- kit193Table(design, kit193Missing("G2", 3L))
    # both comparisons pass on the shipped code first
    expect_identical(exceedsRoundedPercent(0:3, 3, cutoff), 0:3 > row$kmax)
    expect_setequal(kit193Kept(helper_table, design, cutoffs), paste0("k", 0:row$kmax))
    testthat::local_mocked_bindings(roundedPercentNumerator = function(cutoff_percent) 2 * round(1000 * cutoff_percent))
    testthat::expect_failure(expect_identical(exceedsRoundedPercent(0:3, 3, cutoff), 0:3 > row$kmax))
    testthat::expect_failure(expect_setequal(kit193Kept(helper_table, design, cutoffs), paste0("k", 0:row$kmax)))
})
