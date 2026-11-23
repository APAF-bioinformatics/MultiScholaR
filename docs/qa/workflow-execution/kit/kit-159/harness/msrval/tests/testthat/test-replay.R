# Replay (KIT-159 7.3: valReplay; statuses as G5: identical, differs, unsupported), with the
# record builder and the replay as functions so that the comparison is tested without a store.

daRows <- function(shift = 0) {
    data.frame(feature_id = c("P1", "P2"), assay = "all", contrast = "KO_vs_WT",
               log2fc = c(1, -1) + shift, p = c(0.01, 0.2), s_own = c(0.02, 0.2),
               call_statistic = "qvalue", called_own = c(TRUE, FALSE), mean_expression = c(10, 11))
}

test_that("a replay giving the same DA table is identical", {
    res <- valReplay(manager = NULL, original_da = daRows(),
                     build = function(manager) list(steps = list()),
                     replay = function(record) daRows())
    expect_identical(res$status, "identical")
    expect_true(is.na(res$class))
})

test_that("a replay giving other values differs, with the largest difference", {
    res <- valReplay(manager = NULL, original_da = daRows(),
                     build = function(manager) list(), replay = function(record) daRows(1e-3))
    expect_identical(res$status, "differs")
    expect_equal(res$max_abs_diff, 1e-3, tolerance = 1e-9)
    expect_identical(res$rows_only_one_side, 0L)
})

test_that("a replay the package refuses is unsupported, with its class", {
    refusing <- function(record) {
        cli::cli_abort("cannot replay", class = c("workflow_replay_invalid", "rlang_error"))
    }
    res <- valReplay(manager = NULL, original_da = daRows(), build = function(manager) list(),
                     replay = refusing)
    expect_identical(res$status, "unsupported")
    expect_identical(res$class, "workflow_replay_invalid")
})
