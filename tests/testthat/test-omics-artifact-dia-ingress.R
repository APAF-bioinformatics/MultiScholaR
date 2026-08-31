test_that("streaming DIA ingestion is exact and row-group invariant", {
    testthat::skip_if_not_installed("arrow")
    testthat::skip_if_not_installed("DBI")
    testthat::skip_if_not_installed("duckdb")
    sources <- c(
        testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv"),
        testthat::test_path(
            "..", "testdata", "e2e", "prot_dia_limpa", "seed_report.tsv"
        )
    )
    for (source in sources) {
        expected <- suppressMessages(importDIANNData(source))$data
        digests <- character()
        for (rows in c(2048L, 65536L)) {
            path <- tempfile(fileext = ".parquet")
            withr::defer(unlink(path, force = TRUE))
            streamed <- writeProtDiaStreamingParquet(
                source,
                path,
                row_group_rows = rows
            )
            encoded <- encodeArtifactStreamingParquet(path, "DIA ingress test")
            observed <- decodeArtifactRectangular(
                encoded$payload,
                encoded$metadata
            )
            expect_identical(observed, expected)
            expect_silent(verifyArtifactStreamingPayload(
                encoded$payload,
                encoded$metadata
            ))
            expect_identical(
                unlist(streamed$import_summary$run_values, use.names = FALSE),
                unique(expected$Run)
            )
            digests <- c(digests, encoded$metadata$semantic_digest)
        }
        expect_identical(digests[[1L]], digests[[2L]])
    }
})

test_that("streaming semantic evidence rejects block mutation", {
    testthat::skip_if_not_installed("arrow")
    testthat::skip_if_not_installed("DBI")
    testthat::skip_if_not_installed("duckdb")
    source <- testthat::test_path(
        "..", "testdata", "e2e", "prot_dia", "report.tsv"
    )
    path <- tempfile(fileext = ".parquet")
    withr::defer(unlink(path, force = TRUE))
    writeProtDiaStreamingParquet(source, path)
    encoded <- encodeArtifactStreamingParquet(path, "DIA mutation test")
    changed <- encoded$metadata
    changed$streaming_semantic$block_digests[[1L]] <- strrep("0", 64L)

    expect_error(
        verifyArtifactStreamingPayload(encoded$payload, changed),
        class = "multischolar_artifact_semantic_digest_mismatch"
    )
})

test_that("streaming DIA ingress preserves sample sanitization", {
    testthat::skip_if_not_installed("arrow")
    testthat::skip_if_not_installed("DBI")
    testthat::skip_if_not_installed("duckdb")
    source <- testthat::test_path(
        "..", "testdata", "e2e", "prot_dia", "report.tsv"
    )
    expected <- sanitizeProtDiaArtifactImport(
        suppressMessages(importDIANNData(source))
    )$data
    path <- tempfile(fileext = ".parquet")
    withr::defer(unlink(path, force = TRUE))
    writeProtDiaStreamingParquet(source, path, sanitize_names = TRUE)
    encoded <- encodeArtifactStreamingParquet(path, "DIA sanitize test")

    expect_identical(
        decodeArtifactRectangular(encoded$payload, encoded$metadata),
        expected
    )
})

test_that("DIA preflight and streaming code remain payload bounded", {
    sources <- unlist(lapply(c(
        testthat::test_path("..", "..", "R", "mod_prot_dia_ingress_helpers.R"),
        testthat::test_path("..", "..", "R", "utils_artifact_streaming_codecs.R")
    ), readLines, warn = FALSE), use.names = FALSE)

    expect_false(any(grepl("object.size", sources, fixed = TRUE)))
    expect_false(any(grepl("gc(", sources, fixed = TRUE)))
    expect_false(any(nchar(sources, type = "width") > 100L))
})

test_that("deferred design state retains only run identity operations", {
    workflow <- new.env(parent = emptyenv())
    workflow$artifact_import_summary <- list(
        rows = 100,
        columns = 5L,
        column_names = as.list(c("Run", "Protein.Group", "Stripped.Sequence",
            "Precursor.Quantity", "Q.Value")),
        run_values = as.list(c("WT 1", "KO 1")),
        protein_count = 10,
        peptide_count = 20
    )
    input <- protDiaDeferredDesignInput(workflow)
    renamed <- applyProtDesignSingleRenameUpdate(
        designMatrix = data.frame(Run = input$Run),
        dataCln = input,
        originalName = "WT 1",
        newName = "WT_1"
    )$dataCln
    mapping <- protDiaDeferredRunMapping(renamed)

    expect_identical(names(input), c("Run", ".multischolar_original_run"))
    expect_identical(mapping$original_run, c("WT 1", "KO 1"))
    expect_identical(mapping$Run, c("WT_1", "KO 1"))
})
