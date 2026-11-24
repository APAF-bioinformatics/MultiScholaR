# KIT-193 target 4 (ticket section 6.4, 5.5): q-value against Benjamini-Hochberg, the library-level parts.
# (i) qvalue::qvalue's q-values are pi0 times stats::p.adjust(p, "BH") (QvalueBh.qv, bhSorted; 4a's nesting);
# (iii) the finite counterexamples of 4d-i: tests/testdata/kit193/bh-counterexample.tsv is written by
# proofs/lean/Kit193Oracles.lean (Oracle/QvalueBh.lean, the Spec's callsKept and callsWhole evaluated) and lists the
# calls of each outcome of the N = 20 and N = 25 grids; stats::p.adjust is applied to the same outcomes here.
# Kind: library (bindings.tsv rows T4 qvalue and stats). Not here: (ii), the route's stored tables and the flag of
# processEnrichments (R/func_general_enrichment_orchestration.R), which KIT-203 changes; that part and its binding row
# wait for KIT-203's merge (ticket 8.1).

#' TRUE where q equals pi0 * BH to relative 1e-14 (exactly where both are 0)
kit193QvAgrees <- function(q, p, pi0) {
    bh <- pi0 * stats::p.adjust(p, "BH")
    q == bh | abs(q - bh) <= 1e-14 * abs(bh)
}

#' 20 seeded p-value vectors of length 1 to 200, with ties, zeros and ones
kit193QvVectors <- function() {
    lengths <- c(1L, 2L, 3L, 5L, 8L, 13L, 20L, 30L, 45L, 60L, 80L, 100L, 120L, 140L, 150L, 160L, 175L, 185L, 195L, 200L)
    lapply(seq_along(lengths), function(i) {
        m <- lengths[[i]]
        p <- withr::with_seed(1000L + i, {
            x <- stats::runif(m)^3
            if (m >= 4L) x[sample(m, max(1L, m %/% 5L))] <- sample(x, 1L)
            if (m >= 3L) {
                x[sample(m, 1L)] <- 0
                x[sample(m, 1L)] <- 1
            }
            x
        })
        p
    })
}

test_that("KIT-193 T4 pair (i): qvalue::qvalue with a given pi0 is pi0 * p.adjust(p, 'BH'); the q flag contains the BH flag", {
    expect_identical(format(utils::packageVersion("qvalue")), "2.44.0")
    vectors <- kit193QvVectors()
    expect_true(any(vapply(vectors, function(p) anyDuplicated(p) > 0L, logical(1L))))
    expect_true(any(vapply(vectors, function(p) any(p == 0) && any(p == 1), logical(1L))))
    for (i in seq_along(vectors)) for (pi0 in c(1, 0.5, 0.13)) {
        p <- vectors[[i]]
        q <- qvalue::qvalue(p, pi0 = pi0, lfdr.out = FALSE)$qvalues
        case <- sprintf("vector %d (length %d), pi0 %g", i, length(p), pi0)
        expect_true(all(kit193QvAgrees(q, p, pi0)), info = case)
        bh <- stats::p.adjust(p, "BH")
        flagged <- bh < 0.05
        expect_true(all(q[flagged] < 0.05), info = paste(case, "exceptions (BH, q):",
            paste(sprintf("(%.17g, %.17g)", bh[flagged & q >= 0.05], q[flagged & q >= 0.05]), collapse = " ")))
    }
})

test_that("KIT-193 T4 pair (i): with pi0 estimated by qvalue, q-values are pi0 * BH and nest the BH flag", {
    for (i in 1:5) {
        p <- withr::with_seed(2000L + i, {
            x <- c(stats::rbeta(60L, 0.3, 6), stats::runif(130L), stats::runif(10L, 0.951, 1))
            sample(x)
        })
        expect_gte(sum(p > 0.95), 10L)
        fit <- qvalue::qvalue(p, lfdr.out = FALSE)
        expect_true(fit$pi0 > 0 && fit$pi0 <= 1, info = i)
        expect_true(all(kit193QvAgrees(fit$qvalues, p, fit$pi0)), info = i)
        bh <- stats::p.adjust(p, "BH")
        expect_true(all(fit$qvalues[bh < 0.05] < 0.05), info = i)
    }
})

kit193BhTable <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "bh-counterexample.tsv"), sep = "\t", quote = "")
}

#' p-value of grid outcome j: 1 at j = 0, else j / (2N)
kit193BhValue <- function(N, j) ifelse(j == 0L, 1, j / (2 * N))

#' Calls of stats::p.adjust over the kept terms (p < 1) or the whole family, at "<= q" or "< q"
kit193BhCalls <- function(p, kept, strict, q = 0.05, kept_rule = kept) {
    family <- if (kept_rule) p[p < 1] else p
    if (!length(family)) return(0L)
    adjusted <- stats::p.adjust(family, "BH")
    sum(if (strict) adjusted < q else adjusted <= q)
}

#' The calls of every outcome of the table by p.adjust, and the false discovery rate under the complete null
kit193BhCompare <- function(table, kept_rule = NULL) {
    columns <- c(kept_le = "kept_le", kept_lt = "kept_lt", whole_le = "whole_le", whole_lt = "whole_lt")
    calls <- vapply(seq_len(nrow(table)), function(i) {
        p <- kit193BhValue(table$N[[i]], c(table$j1[[i]], table$j2[[i]]))
        c(kept_le = kit193BhCalls(p, TRUE, FALSE, kept_rule = kept_rule %||% TRUE),
            kept_lt = kit193BhCalls(p, TRUE, TRUE, kept_rule = kept_rule %||% TRUE),
            whole_le = kit193BhCalls(p, FALSE, FALSE), whole_lt = kit193BhCalls(p, FALSE, TRUE))
    }, integer(4L))
    weight <- ifelse(table$j1 == 0L, 1 / 2, 1 / (2 * table$N)) * ifelse(table$j2 == 0L, 1 / 2, 1 / (2 * table$N))
    fdr <- vapply(columns, function(column) sum(weight * (calls[column, ] > 0L)), numeric(1L))
    list(calls = calls, fdr = fdr)
}

test_that("KIT-193 T4 pair (iii): p.adjust on the 441 and 676 grid outcomes calls what the Lean model calls; FDRs 3/40, 1/20, 3/50, 1/25", {
    table <- kit193BhTable()
    expect_identical(as.integer(table(table$N)), c(441L, 676L))
    for (N in c(20L, 25L)) {
        rows <- table[table$N == N, , drop = FALSE]
        got <- kit193BhCompare(rows)
        for (column in rownames(got$calls)) {
            expect_identical(unname(got$calls[column, ]), as.integer(rows[[column]]), info = paste(N, column))
        }
        expected <- if (N == 20L) c(kept_le = 3 / 40, kept_lt = 41 / 1600, whole_le = 1 / 20) else
            c(kept_le = 3 / 50, kept_lt = 3 / 50, whole_le = 1 / 25, whole_lt = 1 / 25)
        for (column in names(expected)) {
            expect_lt(abs(got$fdr[[column]] - expected[[column]]), 1e-12, label = paste(N, column))
        }
    }
})

test_that("KIT-193 T4 pair (iii), planted defect: BH over the whole family where the kept family is asked fails 3/40", {
    table <- kit193BhTable()
    rows <- table[table$N == 20L, , drop = FALSE]
    planted <- kit193BhCompare(rows, kept_rule = FALSE)
    # one expectation, so expect_failure() sees exactly one failure
    testthat::expect_failure(expect_lt(abs(planted$fdr[["kept_le"]] - 3 / 40), 1e-12))
})

test_that("KIT-193 T4 pair (i), planted defect: qvalue::qvalue ignoring pi0 fails the comparison", {
    p <- kit193QvVectors()[[12L]]
    original <- qvalue::qvalue
    testthat::local_mocked_bindings(qvalue = function(p, ...) {
        out <- original(p, ...)
        out$qvalues <- stats::p.adjust(p, "BH")
        out
    }, .package = "qvalue")
    q <- qvalue::qvalue(p, pi0 = 0.5, lfdr.out = FALSE)$qvalues
    # one expectation, so expect_failure() sees exactly one failure
    testthat::expect_failure(expect_true(all(kit193QvAgrees(q, p, 0.5))))
})
