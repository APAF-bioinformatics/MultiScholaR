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

# What a corpus file holds: a Parquet file's table, a JSON file's values (an omic manifest without its file
# hashes, which identify the bytes on the machine that wrote them), a TSV file's table.
intCorpusContent <- function(path) {
    if (grepl("\\.parquet$", path)) {
        return(as.data.frame(arrow::read_parquet(path)))
    }
    if (grepl("\\.json$", path)) {
        content <- jsonlite::read_json(path, simplifyVector = FALSE)
        if (identical(basename(path), "integration_input.json")) content$files <- names(content$files)
        return(content)
    }
    utils::read.delim(path, colClasses = "character", stringsAsFactors = FALSE)
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

    # The same bytes on one machine: composed twice, the same hashes.
    fresh_dir <- withr::local_tempdir()
    fresh <- h$composeIntegrationCorpus(fresh_dir, repo = intCorpusRoot())
    again <- h$composeIntegrationCorpus(withr::local_tempdir(), repo = intCorpusRoot())
    expect_identical(intCorpusHashes(again), intCorpusHashes(fresh))

    # The same content as the committed corpus on every platform: another platform's Arrow build may compress the
    # same pages into other bytes, so the fresh files are compared with the committed ones by what they hold.
    expect_identical(names(fresh$files), names(committed$files))
    for (path in names(committed$files)) {
        expect_identical(
            intCorpusContent(file.path(fresh_dir, path)),
            intCorpusContent(file.path(intCorpusDir(), path)),
            info = path
        )
    }

    other <- h$composeIntegrationCorpus(
        withr::local_tempdir(),
        seed = h$INTEGRATION_CORPUS_SEED + 1L,
        repo = intCorpusRoot()
    )
    matrices <- grep("matrix_", names(committed$files), value = TRUE)
    expect_false(any(intCorpusHashes(other)[matrices] == intCorpusHashes(fresh)[matrices]))
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
    expect_identical(as.numeric(scores$Z1), as.numeric(scores$group == "KO"))

    # KIT-147: the omics' inputs are contract v2 Parquet files, read here with Arrow alone.
    samples <- lapply(omics, \(o) as.data.frame(arrow::read_parquet(file.path(dir, o, "samples.parquet"))))
    names(samples) <- omics
    runs <- unlist(lapply(samples, `[[`, "Run"))
    expect_false(anyDuplicated(runs) > 0L)
    # KIT-062: the cross-omic key is tech_rep_group in the omics' format (KO_5).
    absent <- lapply(samples, \(s) setdiff(scores$biological_sample, s$tech_rep_group))
    expect_false(any(vapply(samples, \(s) "biological_sample" %in% names(s), logical(1))))
    expect_true(all(grepl("^(WT|KO)_[1-9][0-9]*$", scores$biological_sample)))
    expect_identical(lengths(absent)[lengths(absent) > 0L], c(metabolomics = 1L))
    expect_identical(absent$metabolomics, truth$dropped$biological_sample)

    for (omic in omics) {
        for (view in truth$views[[omic]]) {
            m <- as.data.frame(arrow::read_parquet(file.path(dir, omic, sprintf("matrix_%s.parquet", view))))
            f <- as.data.frame(arrow::read_parquet(file.path(dir, omic, sprintf("features_%s.parquet", view))))
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

test_that("each omic's folder is a KIT-061 integration input the contract's reader accepts", {
    dir <- intCorpusDir()
    testthat::skip_if_not(dir.exists(dir), "the committed corpus is not present")
    truth <- jsonlite::read_json(file.path(dir, "truth.json"), simplifyVector = TRUE)
    workflows <- c(proteomics = "LFQ", transcriptomics = "TRANSCRIPTOMICS", lipidomics = "LIPIDOMICS",
        metabolomics = "METABOLOMICS")
    for (omic in names(workflows)) {
        input <- readIntegrationInput(file.path(dir, omic))
        manifest <- input$manifest
        # KIT-147: contract v2.
        expect_identical(manifest$schema, INTEGRATION_INPUT_SCHEMA)
        expect_identical(manifest$omic_type, omic)
        expect_identical(manifest$workflow, workflows[[omic]])
        expect_identical(manifest$matrix_role, "normalised")
        expect_null(manifest$correction)
        expect_identical(manifest$biological_sample_column, "tech_rep_group")
        # INT-4: views as the writers name them.
        expected_views <- if (omic %in% c("lipidomics", "metabolomics")) {
            startsWith(names(input$views), paste0(omic, "_"))
        } else {
            names(input$views) == omic
        }
        expect_true(all(expected_views))
        expect_identical(names(input$views), truth$views[[omic]])
        expect_identical(nrow(input$samples), as.integer(truth$runs[[omic]]))
        # Nothing in the manifest varies by machine, package version or time of writing.
        expect_null(manifest$software$version)
        expect_identical(manifest$written_at, "2026-11-17T00:00:00Z")
    }
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
