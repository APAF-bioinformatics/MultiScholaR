# Provenance, determinism and phase score (KIT-159 5.6, 7.11, 7.8; section 10.3: T-PROV-1,
# T-DET-1, T-DET-2, T-SCORE-1).

provFixture <- function() {
    valProvenance(
        ticket = "KIT-999", stage = "lifecycle", label = "life1", phase = "lifecycle",
        status = "completed", blinded = TRUE, started = "2026-11-23T00:00:00Z",
        ended = "2026-11-23T01:00:00Z",
        code = list(package_commit = strrep("a", 40), release_commit = strrep("a", 40),
                    product_tree_identical = TRUE,
                    harness = list(commit = strrep("b", 40), tree = strrep("c", 40)),
                    declaration = list(stage = "lifecycle", commit = strrep("d", 40),
                                       sha256_md = strrep("e", 64), sha256_json = strrep("f", 64)),
                    code_files = list(), helpers_extracted = list()),
        inputs = list(universe = list(sha256 = strrep("1", 64), rows = 20L)),
        truth = list(primary = list(sha256 = strrep("2", 64), rows = 20L),
                     blind = list(seed = 159101L, applied = TRUE)),
        cells = list(list(unit = "u01", arm = "shipped", exit = "completed", seed = 1L,
                          replay = list(status = "not_declared"))),
        scoring = list(metrics_sha256 = strrep("3", 64), comparisons_sha256 = strrep("4", 64)),
        outputs = list(list(path = "metrics.tsv", sha256 = strrep("3", 64)))
    )
}

test_that("T-PROV-1: a complete provenance passes; a missing field fails", {
    p <- provFixture()
    expect_true(valCheckProvenance(p))
    expect_true(nzchar(p$run$runner$image_os))
    expect_true(p$run$runner$cores >= 1)
    q <- p
    q$code$harness$tree <- NULL
    expect_error(valCheckProvenance(q), class = "val_provenance_invalid")
    r <- p
    r$status <- 3
    expect_error(valCheckProvenance(r), class = "val_provenance_invalid")
    # I07: 5.6's nested fields of a cell phase are checked (R0 showed the fixture alone passing
    # on the tests-only commit, so each nested field gets its own case)
    for (drop in list(c("inputs", "universe"), c("truth", "blind"), c("scoring", "metrics_sha256"),
                      c("code", "product_tree_identical"))) {
        s <- p
        s[[drop[[1L]]]][[drop[[2L]]]] <- NULL
        expect_error(valCheckProvenance(s), class = "val_provenance_invalid",
                     info = paste(drop, collapse = "."))
    }
    s <- p
    s$cells[[1L]]$seed <- NULL
    expect_error(valCheckProvenance(s), class = "val_provenance_invalid")
    path <- tempfile(fileext = ".json")
    valWriteProvenance(p, path)
    expect_true(valCheckProvenance(path))
})

test_that("T-PROV-2: every cell carries a replay status of G5; per-unit truths are hashed (I07, E4, E5)", {
    p <- provFixture()
    for (st in c("not_built", "", NA_character_)) {
        s <- p
        s$cells[[1L]]$replay$status <- st
        expect_error(valCheckProvenance(s), class = "val_provenance_invalid", info = st)
    }
    s <- p
    s$cells[[1L]]$replay <- NULL
    expect_error(valCheckProvenance(s), class = "val_provenance_invalid")
    s <- p
    s$cells[[1L]]$replay <- list(status = "differs", differences = list())
    expect_error(valCheckProvenance(s), class = "val_provenance_invalid")
    s$cells[[1L]]$replay$differences <- list("da_results")
    expect_true(valCheckProvenance(s))
    s <- p
    s$truth$per_unit <- list(u01 = list(unit = list(sha256 = strrep("5", 64), rows = 20L)))
    expect_true(valCheckProvenance(s))
    s$truth$per_unit$u01$unit$sha256 <- NULL
    expect_error(valCheckProvenance(s), class = "val_provenance_invalid")
})

detTable <- function(values) data.frame(key = paste0("k", seq_along(values)), value = values)

test_that("T-DET-1: identical, within rounding, and different tables", {
    a <- withr::local_tempdir()
    b <- withr::local_tempdir()
    v <- c(1.234567890123, 2.5, 1e-3)
    put <- function(values, dir) {
        tab <- detTable(values)
        tab$value <- sprintf("%.17g", tab$value)
        utils::write.table(tab, file.path(dir, "da.tsv"), sep = "\t", quote = FALSE,
                           row.names = FALSE)
    }
    put(v, a)
    put(v, b)
    keys <- list(da.tsv = "key")
    expect_identical(valCompareRuns(a, b, keys = keys)$result, "identical")
    put(v * (1 + c(1e-15, 0, 0)), b)
    expect_identical(valCompareRuns(a, b, keys = keys)$result, "within_rounding")
    put(v * (1 + c(1e-6, 0, 0)), b)
    res <- valCompareRuns(a, b, keys = keys)
    expect_identical(res$result, "differs")
    expect_equal(res$max_rel_diff, 1e-6, tolerance = 1e-3)
})

test_that("T-DET-2: fingerprints tell a real change from a rounding one", {
    v <- c(1.234567890123, 2.5, 1e-3, 7)
    f0 <- valFingerprint(detTable(v))
    f1 <- valFingerprint(detTable(v * c(1 + 1e-6, 1, 1, 1)))
    f2 <- valFingerprint(detTable(v * c(1 + 1e-14, 1, 1, 1)))
    expect_false(identical(f0$sha256_10, f1$sha256_10))
    expect_identical(f0$sha256_10, f2$sha256_10)
    expect_identical(f0$rows, 4L)
    expect_true(all(c("sum", "min", "max", "q01", "q50", "q99") %in% names(f0$numeric)))
})

scoreFolder <- function(root, label, commit, units, arms = c("a1", "a2")) {
    d <- file.path(root, "runs", label)
    dir.create(d, recursive = TRUE)
    grid <- expand.grid(unit = units, arm = arms, stringsAsFactors = FALSE)
    met <- data.frame(unit = grid$unit, arm = grid$arm, metric = "mae", value = seq_len(nrow(grid)))
    utils::write.table(met, file.path(d, "metrics.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    cells <- data.frame(unit = grid$unit, arm = grid$arm, exit = "completed")
    utils::write.table(cells, file.path(d, "cells.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    jsonlite::write_json(list(code = list(package_commit = commit,
                                          harness = list(tree = strrep("c", 40)),
                                          declaration = list(commit = strrep("d", 40)))),
                         file.path(d, "provenance.json"), auto_unbox = TRUE)
    # 5.7: every run folder holds comparisons.tsv (phase score reads its within-run rows, I06);
    # R1 showed the fixture without it (test mechanics, VAL-54)
    cmp <- as.data.frame(stats::setNames(rep(list(character()), length(
        msrval:::EVIDENCE_COLUMNS[["comparisons.tsv"]])), msrval:::EVIDENCE_COLUMNS[["comparisons.tsv"]]))
    utils::write.table(cmp, file.path(d, "comparisons.tsv"), sep = "\t", quote = FALSE,
                       row.names = FALSE)
    valWriteOutputsSha(d)
    d
}

test_that("T-SCORE-1: phase score combines requests; a missing cell or another commit is refused", {
    root <- withr::local_tempdir()
    commit <- strrep("a", 40)
    scoreFolder(root, "r1", commit, c("u01", "u02"))
    scoreFolder(root, "r2", commit, c("u03", "u04"))
    scoreFolder(root, "r3", commit, c("u05"))
    expected <- expand.grid(unit = sprintf("u%02d", 1:5), arm = c("a1", "a2"),
                            stringsAsFactors = FALSE)
    res <- valScoreRuns(root, c("r1", "r2", "r3"), expected_cells = expected)
    expect_identical(nrow(res$metrics), 10L)
    expect_setequal(unique(res$metrics$unit), sprintf("u%02d", 1:5))
    expected6 <- rbind(expected, data.frame(unit = "u06", arm = "a1"))
    expect_error(valScoreRuns(root, c("r1", "r2", "r3"), expected_cells = expected6),
                 class = "val_spec_invalid")
    scoreFolder(root, "r4", strrep("9", 40), c("u06"))
    expect_error(valScoreRuns(root, c("r1", "r4"), expected_cells = expected[1:4, ]),
                 class = "val_declaration_changed")
    writeLines("tampered", file.path(root, "runs/r1/metrics.tsv"))
    expect_error(valScoreRuns(root, c("r1", "r2", "r3"), expected_cells = expected),
                 class = "val_declaration_changed")
})
