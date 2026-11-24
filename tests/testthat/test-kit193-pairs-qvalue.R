# KIT-193 target 4 (ticket section 6.4, 5.5): q-value against Benjamini-Hochberg.
# (i) qvalue::qvalue's q-values are pi0 times stats::p.adjust(p, "BH") (QvalueBh.qv, bhSorted; 4a's nesting);
# (ii) the route: KIT-173's effect fixture re-run through processEnrichments with significance_rule = "qvalue" (KIT-203
# made the q-value an option; QvalueBh.flag models that option's flag, `significant = qvalue < q_cutoff`,
# R/func_general_enrichment_orchestration.R:669): on enrichit's stored columns q <= BH and the BH flag nested in the
# q-value flag (4a), and the written flag equal to the model's flag (a missing q flags nothing);
# (iii) the finite counterexamples of 4d-i: tests/testdata/kit193/bh-counterexample.tsv is written by
# proofs/lean/Kit193Oracles.lean (Oracle/QvalueBh.lean, the Spec's callsKept and callsWhole evaluated) and lists the
# calls of each outcome of the N = 20 and N = 25 grids; stats::p.adjust is applied to the same outcomes here.
# Kind: library (bindings.tsv rows T4 qvalue and stats) and package (row T4 processEnrichments).

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

#' KIT-173's T7 project (helper-go-propagation.R) run through processEnrichments with significance_rule = "qvalue";
#' the stored up and down results and the written tables
kit193QvRouteRun <- function() {
    work <- tempfile("k193-route-")
    dir.create(work)
    annotations <- goPropAnnotate(goPropT7Raw(), file.path(work, "cache"))
    process <- function(...) {
        args <- list(...)
        args$significance_rule <- "qvalue"
        do.call(processEnrichments, args)
    }
    result <- NULL
    withr::with_envvar(c(MULTISCHOLAR_TEST_MODE = "false"), utils::capture.output(result <- suppressMessages(
        suppressWarnings(goPropRunProject(file.path(work, "project"), annotations, process_fn = process)))))
    tables <- unlist(lapply(result$state@enrichment_data, function(contrast) lapply(contrast, function(x)
        if (!is.null(x)) x@result)), recursive = FALSE)
    files <- list.files(result$paths$pathway_dir, pattern = "_enrichment_results[.]tsv$", full.names = TRUE)
    list(tables = Filter(Negate(is.null), tables), files = files)
}

#' q <= BH on every term with both values (4a: q = pi0 BH with pi0 <= 1)
kit193QvBelowBh <- function(table) {
    both <- !is.na(table$qvalue) & !is.na(table$p.adjust)
    all(table$qvalue[both] <= table$p.adjust[both] * (1 + 1e-12))
}

#' The Lean `flag` on an optional value: present and below alpha
kit193QvFlag <- function(q, alpha) !is.na(q) & q < alpha

test_that("KIT-193 T4 pair (ii): the route under significance_rule = 'qvalue': q <= BH, the BH flag nested, the written flag is the model's", {
    testthat::skip_if_not_installed("GO.db")
    testthat::skip_if_not_installed("clusterProfiler")
    route <- kit193QvRouteRun()
    expect_gt(length(route$tables), 0L)
    expect_gt(length(route$files), 0L)
    # not vacuous: the fixture's tables hold BH calls
    expect_gt(sum(vapply(route$tables, function(x) sum(x$p.adjust < 0.05, na.rm = TRUE), integer(1L))), 0L)
    for (name in names(route$tables)) {
        table <- route$tables[[name]]
        expect_true(kit193QvBelowBh(table), info = name)
        called <- !is.na(table$p.adjust) & table$p.adjust < 0.05
        expect_true(all(kit193QvFlag(table$qvalue[called], 0.05)), info = name)
    }
    for (file in route$files) {
        written <- utils::read.delim(file, sep = "\t", quote = "", stringsAsFactors = FALSE)
        expect_false("p_adjust_family" %in% names(written), label = basename(file))
        expect_identical(written$significant %in% TRUE, kit193QvFlag(written$qvalue, 0.05), info = basename(file))
    }
})

test_that("KIT-193 T4 pair (ii), planted defect: stored q-values pmin(1, 2 * p.adjust) (pi0 = 2) fail q <= BH", {
    testthat::skip_if_not_installed("GO.db")
    testthat::skip_if_not_installed("clusterProfiler")
    # addEnrichmentFamilyAdjustment() is applied to each up and down result before it is stored and written
    # (R/func_general_enrichment_orchestration.R:594-595); the planted version stores q = min(1, 2 BH).
    testthat::local_mocked_bindings(addEnrichmentFamilyAdjustment = function(result, significance_rule) {
        if (!is.null(result)) result@result$qvalue <- pmin(1, 2 * result@result$p.adjust)
        result
    })
    tables <- kit193QvRouteRun()$tables
    # one expectation, so expect_failure() sees exactly one failure
    testthat::expect_failure(expect_true(all(vapply(tables, kit193QvBelowBh, logical(1L)))))
})
