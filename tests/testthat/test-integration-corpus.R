# KIT-137: tier 1 of the synthetic multi-omics corpus for integration (four omics on shared biological
# samples, known factors, per-omic batches, one sample missing from one omic). The composer and its
# committed output live in tests/testdata/integration/; the composer calls the generators of
# tools/profiling/, which are not shipped, so these tests skip under R CMD check.
# Thresholds, declared before the first run: Z2 (shared by every view) correlates with the first
# principal component of the concatenated, group- and batch-adjusted, feature-scaled views at
# |r| >= 0.9; per view and factor, the share of features associated at p < 0.001 is >= 5% where the
# factor is active and <= 1% where it is not.

intCorpusRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)

intCorpusDir <- function() file.path(intCorpusRoot(), "tests", "testdata", "integration", "corpus")

skipWithoutIntegrationComposer <- function() {
    root <- intCorpusRoot()
    testthat::skip_if_not(
        file.exists(file.path(root, "tools", "profiling", "kit_omics_corpus.R")),
        "the corpus generators are not shipped (R CMD check)"
    )
}

intCorpusComposer <- function() {
    h <- new.env(parent = globalenv())
    sys.source(
        file.path(intCorpusRoot(), "tests", "testdata", "integration", "compose_integration_corpus.R"),
        h
    )
    h
}

intCorpusHashes <- function(manifest) {
    vapply(manifest$files, \(f) f$sha256, "")
}

test_that("the composer is deterministic and the committed corpus is its output", {
    skipWithoutIntegrationComposer()
    h <- intCorpusComposer()
    committed <- jsonlite::read_json(file.path(intCorpusDir(), "manifest.json"))
    on_disk <- vapply(
        names(committed$files),
        \(path) digest::digest(file = file.path(intCorpusDir(), path), algo = "sha256"),
        ""
    )
    expect_identical(on_disk, intCorpusHashes(committed))

    fresh <- h$composeIntegrationCorpus(withr::local_tempdir(), repo = intCorpusRoot())
    expect_identical(intCorpusHashes(fresh), intCorpusHashes(committed))

    other <- h$composeIntegrationCorpus(
        withr::local_tempdir(),
        seed = h$INTEGRATION_CORPUS_SEED + 1L,
        repo = intCorpusRoot()
    )
    matrices <- grep("matrix_", names(committed$files), value = TRUE)
    expect_false(any(intCorpusHashes(other)[matrices] == intCorpusHashes(committed)[matrices]))
})

test_that("the corpus has the declared design", {
    skipWithoutIntegrationComposer()
    dir <- intCorpusDir()
    truth <- jsonlite::read_json(file.path(dir, "truth.json"), simplifyVector = TRUE)
    scores <- utils::read.delim(file.path(dir, "truth_scores.tsv"), stringsAsFactors = FALSE)
    loadings <- utils::read.delim(file.path(dir, "truth_loadings.tsv"), stringsAsFactors = FALSE)
    omics <- c("proteomics", "transcriptomics", "lipidomics", "metabolomics")
    expect_setequal(names(truth$views), omics)

    expect_setequal(unique(scores$group), c("WT", "KO"))
    expect_gte(nrow(scores), 16L)
    expect_identical(scores$Z1, as.numeric(scores$group == "KO"))

    samples <- lapply(omics, \(o) utils::read.delim(file.path(dir, o, "samples.tsv"), stringsAsFactors = FALSE))
    names(samples) <- omics
    runs <- unlist(lapply(samples, `[[`, "Run"))
    expect_false(anyDuplicated(runs) > 0L)
    absent <- lapply(samples, \(s) setdiff(scores$biological_sample, s$biological_sample))
    expect_identical(lengths(absent)[lengths(absent) > 0L], c(metabolomics = 1L))
    expect_identical(absent$metabolomics, truth$dropped$biological_sample)

    for (omic in omics) {
        for (view in truth$views[[omic]]) {
            m <- utils::read.delim(file.path(dir, omic, sprintf("matrix_%s.tsv", view)), check.names = FALSE)
            f <- utils::read.delim(file.path(dir, omic, sprintf("features_%s.tsv", view)), stringsAsFactors = FALSE)
            expect_identical(names(m), c("feature_id", samples[[omic]]$Run))
            expect_identical(m$feature_id, f$feature_id)
            expect_identical(m$feature_id, loadings$feature_id[loadings$omic == omic & loadings$view == view])
            expect_identical(anyNA(m[-1L]), omic != "transcriptomics")
        }
    }
    active <- \(k) sort(unique(loadings$omic[loadings[[k]] != 0]))
    expect_identical(active("Z1"), sort(omics))
    expect_identical(active("Z2"), sort(omics))
    expect_identical(active("Z3"), c("proteomics", "transcriptomics"))
})

test_that("the known factors are recoverable without a factor model", {
    skipWithoutIntegrationComposer()
    h <- intCorpusComposer()
    check <- h$checkIntegrationCorpusTruth(intCorpusDir())
    expect_gte(check$z2_pc1_abs_cor, 0.9)
    a <- check$association
    expect_true(all(a$fraction[a$active] >= 0.05))
    expect_true(all(a$fraction[!a$active] <= 0.01))
    expect_setequal(a$view_key[a$factor == "Z3" & a$active], c("proteomics/proteomics", "transcriptomics/transcriptomics"))
})
