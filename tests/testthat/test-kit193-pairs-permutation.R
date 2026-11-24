# KIT-193 target 5 (ticket section 6.5, 5.5): the shipped permutation test of the integration enrichment against the
# Lean model. tests/testdata/kit193/perm-floor.tsv is written by proofs/lean/Kit193Oracles.lean
# (Oracle/Permutation.lean): for fixtures A, B and C and each of the 13,824 cases (tau, pi1, pi2) of four samples and
# two draws, d (draws that leave the factor g[tau] unchanged), r (draws that reflect it about its mean) and (1 + d)/3.
# The shipped function draws its reorderings itself, so the draws are supplied through the seed: the function's own
# expression is evaluated for seeds 1, 2, ... until each of the 576 pairs (pi1, pi2) has a seed. Calls go through
# integrationEnrichmentPathwayTest(), the entry point the KEGG and Reactome tests and the kernel call (the exported
# integrationEnrichFactor() needs a trained project). Kind: package (bindings.tsv rows T5); stats::p.adjust (5c').
# The statistic and p are not in the oracle (the model's T is abstract): the pair checks d, the floor, and the
# finite-core bound #{p <= k/3} <= 4608 k (perm_finite_core_stmt) on the shipped p-values. The tolerance's behaviour at
# a large |t| (proofs audit P05 item 1) is a recorded note, not tested.

kit193PermTable <- function() {
    utils::read.delim(testthat::test_path("..", "testdata", "kit193", "perm-floor.tsv"), sep = "\t", quote = "",
        colClasses = "character", na.strings = character())
}

#' One-based orderings written as digit strings, as a matrix (one row per string)
kit193PermDigits <- function(x) matrix(as.integer(unlist(strsplit(x, ""))), ncol = 4L, byrow = TRUE)

#' The shipped function's reorderings for a seed (`integrationEnrichmentPermutation`, n = 4, B = 2)
kit193PermDraws <- function(seed) {
    withr::with_seed(seed, vapply(seq_len(2L), function(i) sample(4L), integer(4L)), .rng_kind = "Mersenne-Twister",
        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

#' The first seed giving each pair (pi1, pi2), named "pi1|pi2"
kit193PermSeeds <- function(limit = 100000L) {
    seeds <- integer()
    seed <- 0L
    while (length(seeds) < 576L && seed < limit) {
        seed <- seed + 1L
        draws <- kit193PermDraws(seed)
        key <- paste(paste(draws[, 1L], collapse = ""), paste(draws[, 2L], collapse = ""), sep = "|")
        if (is.na(seeds[key])) seeds[[key]] <- seed
    }
    seeds
}

#' Fixture: 12 features by 4 samples, one set of five
kit193PermData <- function() {
    y <- matrix(withr::with_seed(20261124L, stats::rnorm(48L)), 12L, 4L)
    ids <- sprintf("C%05d", seq_len(12L))
    rownames(y) <- ids
    list(y = y, ids = ids, members = list(seq_len(5L)),
        term2gene = data.frame(term = "S1", gene = ids[1:5], stringsAsFactors = FALSE),
        term2name = data.frame(term = "S1", name = "set one", stringsAsFactors = FALSE),
        ranked = data.frame(identifier = ids, value = 0, stringsAsFactors = FALSE),
        settings = c(integrationGseaSettings(p_adjust_method = "BH", min_size = 3L, max_size = 500L, p_cutoff = 1,
            tie_jitter = 0), list(permutations = 2L, time_limit = 1e9)))
}

#' |t| of the observed-statistic route for every ordering of g (rows of `orders`)
kit193PermObservedT <- function(data, g, orders) {
    apply(orders, 1L, function(o) {
        abs(integrationEnrichmentSetT(integrationEnrichmentFisherZ(as.numeric(stats::cor(t(data$y), g[o]))),
            data$members)$t[1L, 1L])
    })
}

#' The shipped test on every case of one fixture: exceed count e = 3p - 1, identity_draws, min_p
kit193PermShipped <- function(data, g, rows, seeds) {
    tau <- kit193PermDigits(rows$tau)
    keys <- paste(rows$pi1, rows$pi2, sep = "|")
    out <- vapply(seq_len(nrow(rows)), function(i) {
        values <- list(y = data$y, f = g[tau[i, ]], counts = list(universe = 12L, samples = 4L))
        res <- integrationEnrichmentPathwayTest(data$ranked, data$term2gene, data$term2name, seeds[[keys[[i]]]],
            data$settings, test = "permutation", values = values)
        c(res$terms$pvalue, res$permutation$identity_draws, res$permutation$min_p)
    }, numeric(3L))
    list(p = out[1L, ], e = round(3 * out[1L, ]) - 1, d = as.integer(out[2L, ]), min_p = out[3L, ])
}

#' The table's (1 + d)/3, written by Lean as "1/3", "2/3" or "1"
kit193PermFraction <- function(x) {
    vapply(strsplit(x, "/", fixed = TRUE), function(s) if (length(s) == 1L) as.numeric(s) else
        as.numeric(s[[1L]]) / as.numeric(s[[2L]]), numeric(1L))
}

#' #{p <= k/3} for k = 1, 2, 3, from the exceed counts (integers, no tolerance)
kit193PermCounts <- function(e) vapply(1:3, function(k) sum(1 + e <= k), integer(1L))

kit193PermFixtures <- list(A = c(1, 2, 4, 8) / 8, B = c(1, 1, 4, 8) / 8, C = c(1, 2, 3, 4) / 4)

test_that("KIT-193 T5 pair: the draws of every pair (pi1, pi2) are reached through the seed", {
    seeds <- kit193PermSeeds()
    expect_identical(length(seeds), 576L)
    table <- kit193PermTable()
    expect_identical(nrow(table), 3L * 13824L)
    expect_true(all(paste(table$pi1, table$pi2, sep = "|") %in% names(seeds)))
    # the seed's draws are the pair named
    first <- names(seeds)[[1L]]
    expect_identical(kit193PermDraws(seeds[[first]]), unname(t(kit193PermDigits(strsplit(first, "|", fixed = TRUE)[[1L]]))))
})

test_that("KIT-193 T5 pair: fixtures A and B: identity_draws = d and min_p = (1 + d)/3 on every case; p >= (1 + d)/3; on A #{p <= k/3} <= 4608 k", {
    seeds <- kit193PermSeeds()
    table <- kit193PermTable()
    data <- kit193PermData()
    orders <- kit193PermDigits(unique(table$tau))
    # Condition of fixture A (not judged): the 24 values |t(g o sigma)| are pairwise more than 1e-8 apart, so the
    # tolerance 1e-9 does not act and the shipped count is the exact count.
    t_a <- kit193PermObservedT(data, kit193PermFixtures$A, orders)
    expect_gt(min(diff(sort(t_a))), 1e-8)
    for (name in c("A", "B")) {
        rows <- table[table$fixture == name, , drop = FALSE]
        got <- kit193PermShipped(data, kit193PermFixtures[[name]], rows, seeds)
        d <- as.integer(rows$d)
        expect_identical(got$d, d, info = name)
        expect_identical(got$min_p, kit193PermFraction(rows$min_p), info = name)
        expect_true(all(got$e >= d), info = name)
        if (identical(name, "A")) {
            counts <- kit193PermCounts(got$e)
            expect_true(all(counts <= 4608L * 1:3), info = paste(counts, collapse = ", "))
        }
    }
})

test_that("KIT-193 T5 pair: fixture C (symmetric): identity_draws = d; p >= (1 + d + r)/3; #{p <= k/3} <= 4608 k", {
    seeds <- kit193PermSeeds()
    table <- kit193PermTable()
    data <- kit193PermData()
    orders <- kit193PermDigits(unique(table$tau))
    # Condition of fixture C (not judged): the 12 orderings not related by reflection (g's first value 1/4 or 2/4)
    # have |t| pairwise more than 1e-8 apart; a reflection gives t* = -t exactly, which the tolerance must count.
    g <- kit193PermFixtures$C
    t_c <- kit193PermObservedT(data, g, orders[g[orders[, 1L]] <= 2 / 4, , drop = FALSE])
    expect_identical(length(t_c), 12L)
    expect_gt(min(diff(sort(t_c))), 1e-8)
    rows <- table[table$fixture == "C", , drop = FALSE]
    got <- kit193PermShipped(data, g, rows, seeds)
    d <- as.integer(rows$d)
    r <- as.integer(rows$r)
    expect_gt(sum(r), 0L)
    expect_identical(got$d, d)
    expect_true(all(got$e >= d + r))
    counts <- kit193PermCounts(got$e)
    expect_true(all(counts <= 4608L * 1:3), info = paste(counts, collapse = ", "))
})

test_that("KIT-193 T5 pair: the BH floor in whole numbers (5c): 60 at B = 1,000, d = 10, m = 273; ceiling(20 m (1 + d)/(1 + B)) on a grid", {
    expect_identical(integrationEnrichmentFloorMinimum(10, 1000, 273), 60L)
    texts <- integrationGuiTexts()
    expect_identical(integrationEnrichmentPermutationFloorText(list(permutations = 1000L, identity_draws = 10L, samples = 12L),
        273L, texts), sprintf(texts[["enrich_test_permutation_floor"]], 1000L, 12L, "11/1,001", 273L, 60L))
    # the exact ceiling of 20 m (1 + d) / (1 + B), computed as -floor(-a / b) on whole numbers (exact in doubles)
    grid <- expand.grid(m = 1:6000, d = 0:119, B = c(1000, 10000))
    exact <- -((-(20 * grid$m * (1 + grid$d))) %/% (1 + grid$B))
    expect_identical(integrationEnrichmentFloorMinimum(grid$d, grid$B, grid$m), as.integer(exact))
})

test_that("KIT-193 T5 pair: the floor sentence's m is the number of terms BH adjusted over (the kernel's family)", {
    n <- 12L
    g <- 40L
    f <- c(-1.5, -1.3, -1.0, -0.7, -0.4, -0.2, 0.1, 0.3, 0.6, 0.9, 1.2, 1.6)
    y <- 16 + withr::with_seed(20261125L, matrix(stats::rnorm(g * n), g, n))
    y[1:8, ] <- y[1:8, ] + 0.8 * matrix(f, 8L, n, byrow = TRUE)
    ids <- sprintf("C%05d", seq_len(g))
    rownames(y) <- ids
    sets <- list(1:8, 9:16, 17:24, 25:32, 33:37, 38:40, c(1:3, 20:23))
    pathway <- sprintf("path:hsa%05d", seq_along(sets))
    table <- data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", ids[unlist(sets)]),
        name = rep(paste("pathway", seq_along(sets)), lengths(sets)), stringsAsFactors = FALSE)
    testthat::local_mocked_bindings(integrationEnrichmentMappedList = function(ranking, p, features) {
        list(table = data.frame(identifier = ids, value = 0, feature_id = ids, stringsAsFactors = FALSE),
            counts = list(), rule = "kit193", mapping = data.frame())
    })
    p <- list(factor = "Factor1", view = "metab", omic = "metabolomics", database = "kegg", id_type = "kegg",
        id_column = "kegg", correlation_method = "pearson", min_samples = 5L, kegg_organism = "hsa", min_set_size = 5L,
        max_set_size = 500L, seed = 7L, test = "permutation", permutations = 1000L, time_limit = 1e9)
    out <- integrationEnrichmentKernel(NULL, list(p = p, ranking_value = list(direction = list(), counts = list()),
        mapping_table = table, mapping_meta = list(release = "kit193", sha256 = strrep("e", 64L), url = "kit193",
            downloaded_at = "2026-11-24T00:00:00+0000"),
        values = list(y = y, f = f, counts = list(universe = g, samples = n))))
    # six pathways of at least five members are tested; the seventh (three members) is not in the family
    expect_identical(nrow(out$terms), 6L)
    expect_identical(out$tests, nrow(out$terms))
    expect_identical(out$terms$p.adjust, stats::p.adjust(out$terms$pvalue, "BH"))
    d <- out$permutation$identity_draws
    texts <- integrationGuiTexts()
    expect_identical(integrationEnrichmentPermutationFloorText(out$permutation, as.integer(out$tests), texts),
        sprintf(texts[["enrich_test_permutation_floor"]], 1000L, n, integrationEnrichmentFloorFraction(d, 1000L), 6L,
            as.integer(-((-(20 * 6 * (1 + d))) %/% 1001))))
})

test_that("KIT-193 T5 pair: KIT-191's earlier floor sentence was false (5c'), by stats::p.adjust on its 100 values", {
    p <- c(rep(2 / 1001, 4L), rep(1, 96L))
    adjusted <- stats::p.adjust(p, "BH")
    expect_identical(sum(adjusted < 0.05), 4L)
    expect_equal(adjusted[1:4], rep(50 / 1001, 4L), tolerance = 1e-14)
    expect_identical(integrationEnrichmentFloorMinimum(0, 1000, 100), 2L)
    expect_identical(sum(p <= 0.05 * 2 / 100), 0L)
})

test_that("KIT-193 T5 pair, planted defect: permuted t shrunk by 1e-6 (only unchanged draws counted) breaks fixture A's bound", {
    seeds <- kit193PermSeeds()
    table <- kit193PermTable()
    data <- kit193PermData()
    rows <- table[table$fixture == "A", , drop = FALSE]
    original <- integrationEnrichmentSetT
    testthat::local_mocked_bindings(integrationEnrichmentSetT = function(stat, members) {
        out <- original(stat, members)
        if (ncol(as.matrix(stat)) > 1L) out$t <- out$t * 1e-6
        out
    })
    got <- kit193PermShipped(data, kit193PermFixtures$A, rows, seeds)
    # one expectation, so expect_failure() sees exactly one failure
    testthat::expect_failure(expect_lte(kit193PermCounts(got$e)[[1L]], 4608L))
})
