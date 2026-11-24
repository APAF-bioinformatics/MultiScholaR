# KIT-193 target 1 (ticket sections 6.1, 5.5): the shipped missing-value rule against the Lean model of
# proofs/lean/MultiScholaRProofs/Spec/FilterCount.lean. tests/testdata/kit193/filter-thresholds.tsv is written by Lean's
# interpreter from the Spec definitions (Oracle/FilterCount.lean: P, RoundedOK decided, kmax by the proved floor form);
# the record job regenerates it and compares bytes. Kind: package (bindings.tsv rows T1). Bound code: KIT-187 as merged
# (6416dbd6; the lines the Spec cites are unchanged from e6c42d1d to d84b6dfd, ticket 8.2).
# Besides 6.1 (i) to (v), this file closes the five gaps of AUDIT-KIT-187-impl I10:
#   1. the earlier comparison's option is unset in every test, and asserted unset;
#   2. the full helper `removeRowsWithMissingValuesPercentHelper` is driven on a sample of oracle rows, unequal group
#      sizes included, not only `exceedsRoundedPercent`;
#   3. on two-group designs the stored cutoff's P is the model's P of the larger group (the `max` of rounded values);
#   4. the grouping column is `group` in the design read by `updateMissingValueParameters` and in the helper;
#   5. the group level (G groups, g required) runs through the helper.

kit193FilterTable <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "filter-thresholds.tsv"), sep = "\t",
        stringsAsFactors = FALSE)
}

kit193FilterRow <- function(table, nstar, m, n) {
    row <- table[table$nstar == nstar & table$m == m & table$n == n, , drop = FALSE]
    if (nrow(row) != 1L) stop(sprintf("no unique oracle row for n* = %d, m = %d, n = %d", nstar, m, n))
    row
}

# I10 gap 1: the earlier comparison applies only when the marker's option holds its earlier token; the model describes
# the count rule, so the option is unset here.
kit193NoLegacy <- function(env = parent.frame()) {
    option <- workflowReplayRuleMarkers$missing_value_comparison$option
    withr::local_options(stats::setNames(list(NULL), option), .local_envir = env)
    invisible(option)
}

if (!methods::isClass("Kit193Design")) {
    methods::setClass("Kit193Design", slots = c(design_matrix = "data.frame", args = "list"))
}

kit193Design <- function(sizes) {
    groups <- rep(paste0("G", seq_along(sizes)), sizes)
    data.frame(Run = paste0(groups, "_", sequence(sizes)), group = groups, stringsAsFactors = FALSE)
}

#' The shipped cutoffs for a design: `updateMissingValueParameters()` on a minimal S4 object, grouping column `group`
kit193Cutoffs <- function(group_sizes, min_reps, min_groups) {
    object <- methods::new("Kit193Design", design_matrix = kit193Design(group_sizes), args = list())
    env <- new.env()
    assign("config_list", list(), envir = env)
    object <- suppressMessages(updateMissingValueParameters(object, min_reps_per_group = min_reps,
        min_groups = min_groups, grouping_variable = "group", env = env))
    object@args$removeRowsWithMissingValuesPercent
}

#' A table: `missing` maps a row id to the runs set to NA; every other value is 20. `P_floor`'s single value of 10 in
#' the first run sets the intensity threshold at percentile 0 (present means strictly above it).
kit193Table <- function(design, missing) {
    runs <- design$Run
    ids <- c(names(missing), "P_floor")
    values <- matrix(20, nrow = length(ids), ncol = length(runs), dimnames = list(ids, runs))
    for (id in names(missing)) values[id, missing[[id]]] <- NA_real_
    values["P_floor", 1L] <- 10
    data.frame(Protein.Ids = ids, values, check.names = FALSE, stringsAsFactors = FALSE)
}

#' The rows the shipped helper keeps (`P_floor` aside), grouping column `group` (I10 gap 4)
kit193Kept <- function(table, design, cutoffs) {
    kept <- NULL
    utils::capture.output(kept <- suppressMessages(removeRowsWithMissingValuesPercentHelper(table,
        cols = "Protein.Ids", design_matrix = design, sample_id = Run, row_id = Protein.Ids, grouping_variable = group,
        groupwise_percentage_cutoff = cutoffs$groupwise_percentage_cutoff,
        max_groups_percentage_cutoff = cutoffs$max_groups_percentage_cutoff, proteins_intensity_cutoff_percentile = 0)))
    setdiff(as.character(kept$Protein.Ids), "P_floor")
}

#' Rows `k0` to `kn`: row `kj` has the first j runs of group `group` missing
kit193Missing <- function(group, n) {
    stats::setNames(lapply(0:n, function(j) if (j) paste0(group, "_", seq_len(j)) else character()), paste0("k", 0:n))
}

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

test_that("KIT-193 T1 pair (ii), (v): the decision and the peptide floor equal the model on every row (n* <= 60)", {
    option <- kit193NoLegacy()
    expect_null(getOption(option))
    table <- kit193FilterTable()
    expect_identical(nrow(table), 73810L)
    expect_true(all(table$rounded_ok == 1L))
    cutoffs <- new.env()
    for (i in seq_len(nrow(table))) {
        key <- paste(table$nstar[[i]], table$m[[i]])
        if (is.null(cutoffs[[key]])) {
            cutoffs[[key]] <- kit193Cutoffs(table$nstar[[i]], table$m[[i]], 1L)$groupwise_percentage_cutoff
        }
        cutoff <- cutoffs[[key]]
        n <- table$n[[i]]
        k <- 0:n
        expect_identical(exceedsRoundedPercent(k, n, cutoff), k > table$kmax[[i]])
        expect_identical(floor(n * roundedPercentNumerator(cutoff) / 200000), as.numeric(table$kmax[[i]]))
        expect_identical(allowedByRoundedPercent(n, cutoff), as.numeric(table$kmax[[i]]))
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
