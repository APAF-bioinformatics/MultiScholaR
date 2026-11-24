# KIT-193 target 1 (ticket section 6.1, 5.5): the shipped missing-value rule against the Lean model.
# PENDING: moves to tests/testthat/ when KIT-187 merges (it calls exceedsRoundedPercent, roundedPercentNumerator and
# allowedByRoundedPercent, which exist only on kit-187-filter-by-count at e6c42d1d). Before it first runs, section 8.2:
# diff the merged code against the lines Spec/FilterCount.lean cites. Kind: package (bindings.tsv rows T1).
# Table: tests/testdata/kit193/filter-thresholds.tsv (Oracle/FilterCount.lean: P, RoundedOK decided, kmax by the
# proved floor form).

kit193FilterTable <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "filter-thresholds.tsv"), sep = "\t",
        stringsAsFactors = FALSE)
}

methods::setClass("Kit193Design", slots = c(design_matrix = "data.frame", args = "list"))

#' The shipped cutoffs for a design: `updateMissingValueParameters()` on a minimal S4 object
kit193Cutoffs <- function(group_sizes, min_reps, min_groups) {
    design <- data.frame(group = rep(paste0("G", seq_along(group_sizes)), group_sizes), stringsAsFactors = FALSE)
    object <- methods::new("Kit193Design", design_matrix = design, args = list())
    env <- new.env()
    assign("config_list", list(), envir = env)
    object <- suppressMessages(updateMissingValueParameters(object, min_reps_per_group = min_reps,
        min_groups = min_groups, env = env))
    object@args$removeRowsWithMissingValuesPercent
}

test_that("KIT-193 T1 pair (i): every shipped cutoff satisfies RoundedOK, groups of 1 to 316 and groups up to 40", {
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

test_that("KIT-193 T1 pair (iii), (iv): the counterexample rows and KIT-187 T3's rows", {
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
    table <- kit193FilterTable()
    row <- table[table$nstar == 3 & table$m == 2 & table$n == 3, , drop = FALSE]
    cutoff <- kit193Cutoffs(3L, 2L, 1L)$groupwise_percentage_cutoff
    testthat::local_mocked_bindings(roundedPercentNumerator = function(cutoff_percent) 2 * round(1000 * cutoff_percent))
    testthat::expect_failure(expect_identical(exceedsRoundedPercent(0:3, 3, cutoff), 0:3 > row$kmax[[1L]]))
})
