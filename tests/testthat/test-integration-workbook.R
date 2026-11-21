# KIT-143: the integration workbook rewritten over the module (plan item H), checked without MOFA2. Declared before the
# run:
#   - the shipped workbook is the thin one: its code parses, every function it takes from the package exists there, it
#     reads the model's, the results' and the enrichment's settings from their tables, it holds no credential and reads
#     the STRING key from STRING_API_KEY (never a parameter); the stub
#     MultiOmics_workflow_starter.rmd is gone (fails on the old workbook: it names functions the package does not have
#     and carries a key);
#   - KIT-154: the workbook calls only the package's exported integration functions (no `:::`, no internal field),
#     each of them exported in NAMESPACE, documented in man/ and called (fails on the KIT-143 workbook, which called
#     51 internal functions through `:::`, and on a tree without the exports);
#   - the reference workbook's sort-order sample matching (`x[, sort(colnames(x))]` per omic, paired by position)
#     misaligns proteomics and transcriptomics when the transcriptomics runs carry acquisition-order names (run_01, ...),
#     while the module's matching
#     by the saved biological sample puts every run with its own sample (a characterisation of the reference
#     workbook's step, not of new code: it passes on the old code too).
# The workbook rendered end to end and the equivalence with the reference workbook's MOFA2 call are
# test-integration-workbook-mofa.R, on the MOFA2 lane.

intWbPath <- function(name = "integration_workflow_starter.rmd") {
    normalizePath(testthat::test_path("..", "..", "inst", "workbooks", "integration", name), mustWork = FALSE)
}
intWbCorpus <- function() normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)

# The workbook's R code, chunk by chunk, as knitr would run it.
intWbCode <- function() {
    lines <- readLines(intWbPath(), warn = FALSE)
    opens <- grep("^```\\{r", lines)
    closes <- grep("^```\\s*$", lines)
    unlist(lapply(opens, function(start) {
        end <- closes[closes > start][[1L]]
        if (end > start + 1L) lines[(start + 1L):(end - 1L)] else character()
    }))
}

test_that("KIT-143: the shipped integration workbook is the thin one over the module, and the stub is gone", {
    testthat::skip_if_not(file.exists(intWbPath()), "the package source tree is not present")
    expect_false(file.exists(intWbPath("MultiOmics_workflow_starter.rmd")))
    code <- intWbCode()
    expect_gt(length(code), 0L)
    expect_silent(parse(text = code, keep.source = FALSE))
    # KIT-154: the module's entry points are reached through the exported workbook functions (the test below).
    # No credential: no long literal string of letters and digits (a key), and no key assigned from a literal.
    text <- readLines(intWbPath(), warn = FALSE)
    expect_false(any(grepl("[\"'][A-Za-z0-9]{20,}[\"']", text)))
    expect_false(any(grepl("(api_?key|caller_identity|token)\\s*(=|<-)\\s*[\"']", text, ignore.case = TRUE)))
    # KIT-141, INT-21: the STRING key is read from STRING_API_KEY (or asked for), never a parameter of the workbook.
    expect_true(any(grepl("Sys.getenv(\"STRING_API_KEY\"", code, fixed = TRUE)))
    expect_false(any(grepl("key", sub(":.*$", "", text[seq_len(match("---", text[-1L]) + 1L)]), ignore.case = TRUE)))
})

test_that("KIT-154: the workbook calls only the exported integration functions, each exported, documented and called", {
    testthat::skip_if_not(file.exists(intWbPath()), "the package source tree is not present")
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
    testthat::skip_if_not(file.exists(file.path(root, "NAMESPACE")), "the package source tree is not present")
    code <- intWbCode()
    expect_false(any(grepl(":::", code, fixed = TRUE)))
    expect_false(any(grepl(".context", code, fixed = TRUE)))
    api <- c("integrationSettings", "integrationFindInputs", "integrationOpenProject", "integrationMatchSamples",
        "integrationConfirmMatch", "integrationTrainModel", "integrationGetResults", "integrationExportResults",
        "integrationEnrichFactor", "integrationWriteRecord", "integrationCloseProject", "integrationReplayProject")
    exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", readLines(file.path(root, "NAMESPACE")), value = TRUE))
    # Every function of the package the workbook calls (other packages' calls, `pkg::f(`, set aside).
    plain <- gsub("[A-Za-z][A-Za-z0-9.]*::[A-Za-z_.][A-Za-z0-9_.]*", "", code)
    called <- unique(unlist(regmatches(plain, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*(?=\\()", plain, perl = TRUE))))
    from_package <- intersect(called, ls(asNamespace("MultiScholaR"), all.names = TRUE))
    expect_setequal(from_package, api)
    for (name in api) {
        expect_true(name %in% exports, info = name)
        expect_true(file.exists(file.path(root, "man", paste0(name, ".Rd"))), info = name)
        expect_true(exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE), info = name)
    }
})

test_that("KIT-143: the reference workbook's sort-order matching misaligns samples whose run names differ; the module's does not", {
    testthat::skip_if_not(dir.exists(intWbCorpus()), "the tier-1 corpus is not present")
    testthat::skip_if_not_installed("arrow")
    # One run per biological sample, as the reference workbook's inputs had: proteomics as saved, and for
    # transcriptomics the first of each sample's two runs.
    one_run <- function(omic) {
        dir <- file.path(intWbCorpus(), omic)
        manifest <- jsonlite::read_json(file.path(dir, "integration_input.json"))
        samples <- as.data.frame(arrow::read_parquet(file.path(dir, manifest$samples_file)))
        run <- as.character(samples[[manifest$sample_id_column]])
        key <- as.character(samples[[manifest$biological_sample_column]])
        sorted <- order(run, method = "radix")
        first <- sorted[!duplicated(key[sorted])]
        stats::setNames(key[first], run[first])
    }
    proteomics <- one_run("proteomics")
    transcriptomics <- one_run("transcriptomics")
    expect_length(proteomics, 20L)
    expect_length(transcriptomics, 20L)
    # The corpus's run names sort in sample order in both omics (run-tests 44862823716: no pair misaligned). Runs named
    # by the instrument in acquisition order (the order of the file's columns, which the corpus shuffles), as run_01,
    # run_02, ..., do not.
    dir <- file.path(intWbCorpus(), "transcriptomics")
    manifest <- jsonlite::read_json(file.path(dir, "integration_input.json"))
    acquired <- as.character(as.data.frame(arrow::read_parquet(file.path(dir,
        manifest$samples_file)))[[manifest$sample_id_column]])
    names(transcriptomics) <- sprintf("run_%02d", match(names(transcriptomics), acquired))
    # The reference workbook: each omic's columns sorted by name, then the n-th column of each taken as one sample.
    by_order <- data.frame(proteomics = proteomics[sort(names(proteomics))],
        transcriptomics = transcriptomics[sort(names(transcriptomics))], stringsAsFactors = FALSE)
    misaligned <- sum(by_order$proteomics != by_order$transcriptomics)
    message(sprintf("KIT-143 sort-order matching: %d of %d pairs join two different biological samples", misaligned,
        nrow(by_order)))
    expect_gt(misaligned, 0L)
    # The module: each run placed by the biological sample its omic saved.
    read <- integrationReadInputs(file.path(intWbCorpus(), c("proteomics", "transcriptomics")))
    design <- read$design
    rows <- design[design$omic == "proteomics", , drop = FALSE]
    expect_identical(unname(rows$biological_sample[match(names(proteomics), rows$run)]), unname(proteomics))
    rows <- design[design$omic == "transcriptomics", , drop = FALSE]
    expect_true(all(unname(transcriptomics) %in% rows$biological_sample))
    aligned <- integrationAlignViews(integrationViewsFromInputs(file.path(intWbCorpus(), c("proteomics",
        "transcriptomics")))$object, list(views = c("proteomics", "transcriptomics"), sample_match_table = NULL,
        internal_standards = "keep", missing_samples = "keep_missing", feature_selection = "none", scaling = "none"))
    expect_setequal(colnames(aligned$views$proteomics), unname(proteomics))
    expect_identical(colnames(aligned$views$proteomics), colnames(aligned$views$transcriptomics))
})
