# KIT-062: MOFA2 trained end to end on the KIT-137 tier-1 corpus (four omics, six views, twenty biological samples,
# KO_7 absent from metabolomics) with the provisional parameters. Skipped where MOFA2 or a Python with mofapy2 is
# absent (the package's own CI); run on the MOFA2 lane (branch run-mofa-smoke-k062) with the bundle's Python.
# Thresholds, declared before the first run:
#   - each planted factor (Z1, Z2, Z3) correlates with some trained factor at |r| > 0.95;
#   - the factor that matches Z3 explains at least 1% of the variance of proteomics and of transcriptomics and less
#     than 1% of every lipidomics and metabolomics view (the smoke lane's rule for an active factor, >= 1%);
#   - the factors, weights and variance explained equal the stored result (tests/testdata/integration/mofa_reference,
#     written by the lane's first run) within 1e-5 relative (Linux and Windows agree to 5 to 6 digits, INT-5);
#   - a record replayed on a copy of the inputs state gives the same product (identical, one platform).
# MULTISCHOLAR_MOFA_REFERENCE_OUT set: the trained tables are written there (the lane's recording run) and not compared.

intMofaRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
intMofaCorpus <- function() file.path(intMofaRoot(), "tests", "testdata", "integration", "corpus")
intMofaReference <- function() file.path(intMofaRoot(), "tests", "testdata", "integration", "mofa_reference")
intMofaOmics <- c("proteomics", "transcriptomics", "lipidomics", "metabolomics")

intMofaSkip <- function() {
    testthat::skip_if_not(dir.exists(intMofaCorpus()), "the tier-1 corpus is not present")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intMofaProduct <- function(root, identity) {
    owner <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeWorkflowStateOwner(owner), add = TRUE)
    withWorkflowStateScope(owner, function() {
        hydrateSavedWorkflowEntity(owner$.context$store, owner$.context$snapshot$exports$integration_mofa_model$reference)$value
    })
}

# One training on the corpus, its record, and the record replayed on a copy of the inputs state; kept for the file.
intMofaRun <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) return(cache)
        base <- tempfile("kit062-mofa-")
        root <- file.path(base, "project")
        dirs <- file.path(intMofaCorpus(), intMofaOmics)
        manager <- newIntegrationStateOwner(root)
        importIntegrationInputs(manager, dirs)
        closeWorkflowStateOwner(manager)
        prepared <- file.path(base, "prepared")
        dir.create(prepared)
        file.copy(list.files(root, all.files = TRUE, full.names = TRUE, no.. = TRUE), prepared, recursive = TRUE)
        identity <- readSavedStoreIdentity(file.path(root, ".workflow-entities", "store.json"))
        manager <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
        started <- Sys.time()
        run <- runIntegrationModel(manager)
        seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        json <- file.path(base, "analysis_record.json")
        writeWorkflowAnalysisRecord(manager, json)
        record <- buildWorkflowAnalysisRecord(manager)
        closeWorkflowStateOwner(manager)
        inputs <- stats::setNames(file.path(dirs, "integration_input.json"), paste0("integration_input_", intMofaOmics))
        scratch <- file.path(base, "scratch")
        dir.create(scratch)
        replay <- replayWorkflowAnalysisRecord(jsonlite::read_json(json, simplifyVector = FALSE), prepared, inputs, list(),
            scratch)
        cache <<- list(run = run, record = record, seconds = seconds, replay = replay,
            replayed = intMofaProduct(prepared, identity))
        cache
    }
})

intMofaWriteTable <- function(table, path) {
    numeric <- vapply(table, is.numeric, logical(1))
    table[numeric] <- lapply(table[numeric], function(x) ifelse(is.na(x), NA_character_, sprintf("%.17g", x)))
    utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

intMofaReadTable <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

test_that("KIT-062: MOFA2 trains on the tier-1 corpus at the bundle's versions and recovers the planted factors", {
    intMofaSkip()
    out <- intMofaRun()
    value <- out$run$value
    expect_identical(out$run$status, "computed")
    # INT-5: the versions the paper states.
    expect_identical(value$versions$MOFA2, "1.22.0")
    expect_identical(value$versions$mofapy2, "0.7.3")
    # Twenty biological samples, nineteen in every view; KO_7 kept as missing in the metabolomics views (D4).
    expect_identical(nrow(value$samples), 20L)
    expect_identical(value$shared_samples, 19L)
    expect_identical(value$samples$metabolomics_GCMS_runs[value$samples$biological_sample == "KO_7"], "")
    expect_setequal(value$view_summary$view, c("lipidomics_LCMS_Neg", "lipidomics_LCMS_Pos", "metabolomics_GCMS",
        "metabolomics_LCMS_Pos", "proteomics", "transcriptomics"))
    expect_length(value$notes, 0L)
    truth <- utils::read.delim(file.path(intMofaCorpus(), "truth_scores.tsv"), stringsAsFactors = FALSE)
    factors <- value$factors
    planted <- as.matrix(truth[match(factors$biological_sample, truth$biological_sample), c("Z1", "Z2", "Z3")])
    trained <- as.matrix(factors[-1L])
    correlation <- abs(stats::cor(planted, trained))
    best <- apply(correlation, 1L, max)
    message(sprintf("KIT-062 corpus: |r| Z1 %.4f, Z2 %.4f, Z3 %.4f; %d iterations; trained in %.1f s",
        best[["Z1"]], best[["Z2"]], best[["Z3"]], value$elbo$iterations, out$seconds))
    expect_true(all(best > 0.95), info = paste(names(best), round(best, 4), collapse = ", "))
    z3 <- colnames(trained)[[which.max(correlation["Z3", ])]]
    r2 <- value$variance_explained[value$variance_explained$factor == z3, ]
    r2 <- stats::setNames(r2$r2_percent, r2$view)
    active <- r2[c("proteomics", "transcriptomics")]
    inactive <- r2[setdiff(names(r2), names(active))]
    expect_true(all(active >= 1), info = paste(names(r2), round(r2, 3), collapse = ", "))
    expect_true(all(inactive < 1), info = paste(names(r2), round(r2, 3), collapse = ", "))
})

test_that("KIT-062: the trained factors, weights and variance explained equal the stored result within 1e-5 relative", {
    intMofaSkip()
    value <- intMofaRun()$run$value
    tables <- list(factors = value$factors, weights = value$weights, variance_explained = value$variance_explained)
    out <- Sys.getenv("MULTISCHOLAR_MOFA_REFERENCE_OUT", unset = "")
    if (nzchar(out)) {
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        for (name in names(tables)) intMofaWriteTable(tables[[name]], file.path(out, paste0(name, ".tsv")))
        succeed("recording run: the tables are written, not compared")
        return(invisible())
    }
    reference <- intMofaReference()
    for (name in names(tables)) {
        path <- file.path(reference, paste0(name, ".tsv"))
        if (!file.exists(path)) fail(sprintf("no stored result %s (written by the MOFA2 lane's recording run)", path))
        else expect_equal(tables[[name]], intMofaReadTable(path), tolerance = 1e-5, info = name)
    }
})

test_that("KIT-062: the record replays on a copy of the inputs state to the same model", {
    intMofaSkip()
    out <- intMofaRun()
    expect_true(validateWorkflowAnalysisRecord(out$record)$valid)
    expect_identical(out$replay$states, c("initial", "integration_inputs"))
    analysis <- out$record$analyses$integration_mofa_model
    expect_identical(analysis$outcomes$mofapy2, "0.7.3")
    expect_identical(analysis$outcomes$factors, 5L)
    original <- out$run$value
    replayed <- out$replayed
    for (name in c("factors", "weights", "variance_explained", "variance_total", "samples", "view_summary")) {
        expect_identical(replayed[[name]], original[[name]], info = name)
    }
    expect_identical(replayed$elbo$trace, original$elbo$trace)
    expect_identical(replayed$versions, original$versions)
})
