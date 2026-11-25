# KIT-143: the integration workbook rewritten over the module (plan item H), checked without MOFA2 (the shipped
# workbook's static tests were removed by KIT-234: the workbooks are archived). Declared before the run:
#   - KIT-154 (kept by KIT-234): each of the twelve exported integration functions is exported in NAMESPACE, documented
#     in man/ and present (under `load_all(export_all = TRUE)` a dropped `export()` passes every other test; NAMESPACE
#     is edited by hand).
#   - the reference workbook's sort-order sample matching (`x[, sort(colnames(x))]` per omic, paired by position)
#     misaligns proteomics and transcriptomics when the transcriptomics runs carry acquisition-order names (run_01, ...),
#     while the module's matching
#     by the saved biological sample puts every run with its own sample (a characterisation of the reference
#     workbook's step, not of new code: it passes on the old code too).
# The equivalence with the reference workbook's MOFA2 call, and the twelve functions run end to end with MOFA2, are
# test-integration-workbook-mofa.R, on the MOFA2 lane.

intWbCorpus <- function() normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)

test_that("KIT-154 (KIT-234): the twelve integration functions are exported, documented and present", {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
    testthat::skip_if_not(file.exists(file.path(root, "NAMESPACE")), "the package source tree is not present")
    api <- c("integrationSettings", "integrationFindInputs", "integrationOpenProject", "integrationMatchSamples",
        "integrationConfirmMatch", "integrationTrainModel", "integrationGetResults", "integrationExportResults",
        "integrationEnrichFactor", "integrationWriteRecord", "integrationCloseProject", "integrationReplayProject")
    exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", readLines(file.path(root, "NAMESPACE")), value = TRUE))
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
