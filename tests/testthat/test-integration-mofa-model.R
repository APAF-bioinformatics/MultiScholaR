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
# MULTISCHOLAR_MOFA_REFERENCE_OUT set: the trained tables are written there (the lane's recording run) and not compared;
# KIT-151: the run also reports the stored (old) reference against the new tables, in the same job.
# KIT-151: the defaults are the decided ones (selection by variance above the trend, slow, 10,000 rounds, five seeds);
# training stops by its convergence rule, the five seeds are the declared seed and the next four, and the maths library
# is pinned.

intMofaRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
intMofaCorpus <- function() file.path(intMofaRoot(), "tests", "testdata", "integration", "corpus")
intMofaReference <- function() file.path(intMofaRoot(), "tests", "testdata", "integration", "mofa_reference")
intMofaOmics <- c("proteomics", "transcriptomics", "lipidomics", "metabolomics")

intMofaSkip <- function() {
    testthat::skip_if_not(dir.exists(intMofaCorpus()), "the tier-1 corpus is not present")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

# KIT-148: the product as held outside its store (the weights read from it before it is closed).
intMofaProduct <- function(root, identity) {
    owner <- newWorkflowStateOwner(newWorkflowEntityStore(root, identity$project_id, identity$workflow_id), reopen = TRUE)
    on.exit(closeWorkflowStateOwner(owner), add = TRUE)
    integrationModelComparable(owner, withWorkflowStateScope(owner, function() {
        hydrateSavedWorkflowEntity(owner$.context$store, owner$.context$snapshot$exports$integration_mofa_model$reference)$value
    }))
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
        product <- run$value
        run$value <- integrationModelComparable(manager, product)
        json <- file.path(base, "analysis_record.json")
        writeWorkflowAnalysisRecord(manager, json)
        record <- buildWorkflowAnalysisRecord(manager)
        closeWorkflowStateOwner(manager)
        inputs <- stats::setNames(file.path(dirs, "integration_input.json"), paste0("integration_input_", intMofaOmics))
        scratch <- file.path(base, "scratch")
        dir.create(scratch)
        replay <- replayWorkflowAnalysisRecord(jsonlite::read_json(json, simplifyVector = FALSE), prepared, inputs, list(),
            scratch)
        cache <<- list(run = run, product = product, record = record, seconds = seconds, replay = replay,
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
    message(sprintf("KIT-062 corpus: |r| Z1 %.4f, Z2 %.4f, Z3 %.4f; %d iterations (%s); trained in %.1f s",
        best[["Z1"]], best[["Z2"]], best[["Z3"]], value$elbo$iterations, value$elbo$stopped_by, out$seconds))
    # KIT-151: the models trained (one per seed) and each factor's agreement across them; the maths library's pin.
    seeds <- value$seeds
    message(sprintf("KIT-151 seeds %s: rounds %s (%s); lowest agreement %s; median %s", paste(seeds$seeds, collapse = ", "),
        paste(seeds$models$iterations, collapse = ", "), paste(seeds$models$stopped_by, collapse = ", "),
        paste(sprintf("%s %.4f", seeds$stability$factor, seeds$stability$agreement_lowest), collapse = ", "),
        paste(sprintf("%.4f", seeds$stability$agreement_median), collapse = ", ")))
    message(sprintf("KIT-151 maths pin applied %s (%s); libraries %s", value$versions$maths_pin$applied,
        paste(names(value$versions$maths_pin$variables), unlist(value$versions$maths_pin$variables), sep = "=", collapse = " "),
        paste(vapply(value$versions$maths_libraries, function(l) sprintf("%s %s threads %s %s", l$api, l$version, l$threads,
            l$architecture), character(1)), collapse = "; ")))
    expect_identical(value$elbo$stopped_by, "convergence_rule")
    expect_identical(seeds$seeds, 20261119L + 0:4)
    expect_true(isTRUE(value$versions$maths_pin$applied))
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
        # KIT-151: the recorded reference (old) against the new one, in the same job: shapes, the largest differences
        # where the shapes agree, and each new factor's best absolute correlation with an old one.
        for (name in names(tables)) {
            path <- file.path(intMofaReference(), paste0(name, ".tsv"))
            if (!file.exists(path)) next
            old <- intMofaReadTable(path)
            new <- tables[[name]]
            numeric <- names(old)[vapply(old, is.numeric, logical(1))]
            same <- identical(dim(old), dim(new)) && identical(names(old), names(new))
            message(sprintf("KIT-151 reference old against new, %s: old %s, new %s%s", name, paste(dim(old), collapse = " x "),
                paste(dim(new), collapse = " x "), if (same) sprintf("; largest absolute difference %.3g, relative %.3g",
                    max(abs(as.matrix(old[numeric]) - as.matrix(new[numeric]))),
                    max(abs(as.matrix(old[numeric]) - as.matrix(new[numeric]))) / max(abs(as.matrix(old[numeric])))) else ""))
            if (identical(name, "factors")) {
                o <- as.matrix(old[match(new$biological_sample, old$biological_sample), -1L])
                r <- abs(stats::cor(as.matrix(new[-1L]), o, use = "pairwise.complete.obs"))
                message(sprintf("KIT-151 reference old against new, factors: best |r| per new factor %s",
                    paste(sprintf("%s %.4f (%s)", rownames(r), apply(r, 1L, max), colnames(r)[apply(r, 1L, which.max)]),
                        collapse = ", ")))
            }
            if (identical(name, "weights")) {
                message(sprintf("KIT-151 reference old against new, weights: features per data set old %s; new %s",
                    paste(names(table(old$view)), table(old$view), collapse = ", "),
                    paste(names(table(new$view)), table(new$view), collapse = ", ")))
            }
        }
        succeed("recording run: the tables are written, not compared")
        return(invisible())
    }
    reference <- intMofaReference()
    for (name in names(tables)) {
        path <- file.path(reference, paste0(name, ".tsv"))
        if (!file.exists(path)) {
            fail(sprintf("no stored result %s (written by the MOFA2 lane's recording run)", path))
            next
        }
        stored <- intMofaReadTable(path)
        numeric <- names(stored)[vapply(stored, is.numeric, logical(1))]
        if (identical(dim(stored), dim(tables[[name]]))) {
            a <- as.matrix(tables[[name]][numeric])
            b <- as.matrix(stored[numeric])
            message(sprintf("KIT-062 stored-result comparison, %s: largest relative difference %.3g (absolute %.3g)",
                name, max(abs(a - b)) / max(abs(b)), max(abs(a - b))))
        }
        expect_equal(tables[[name]], stored, tolerance = 1e-5, info = name)
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

# KIT-148 item 1: the trainer as KIT-147 left it (save_data = TRUE, the model loaded a second time, its warnings kept),
# verbatim but for the hand-off names, which the kernel now sets at the read. The new trainer must give the same model
# (identical factors, weights, variance explained and trace) and the same warnings in the same order, its overall-level
# warning computed by MultiScholaR from the views.
intMofaKit147Trainer <- function(views, options) {
    object <- MOFA2::create_mofa(views)
    data_options <- MOFA2::get_default_data_options(object)
    data_options$scale_views <- options$scale_views
    data_options$center_groups <- options$center_groups
    data_options$use_float32 <- identical(options$precision, "single")
    model_options <- MOFA2::get_default_model_options(object)
    model_options$num_factors <- options$num_factors
    training_options <- MOFA2::get_default_training_options(object)
    training_options$seed <- options$seed
    training_options$maxiter <- options$maxiter
    training_options$convergence_mode <- options$convergence_mode
    training_options$verbose <- FALSE
    object <- MOFA2::prepare_mofa(object, data_options = data_options, model_options = model_options,
        training_options = training_options)
    outfile <- tempfile("integration-mofa-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    messages <- character()
    withCallingHandlers(MOFA2::run_mofa(object, outfile = outfile, save_data = TRUE, use_basilisk = FALSE),
        message = function(m) {
            messages <<- c(messages, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
            invokeRestart("muffleWarning")
        })
    model <- withCallingHandlers(MOFA2::load_model(outfile), warning = function(w) {
        text <- paste("WARNING:", trimws(conditionMessage(w)))
        if (!text %in% messages) messages <<- c(messages, text)
        invokeRestart("muffleWarning")
    })
    variance <- MOFA2::get_variance_explained(model)
    r2 <- variance$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(views))) r2 <- t(r2)
    stats_slot <- methods::slot(model, "training_stats")
    list(factors = MOFA2::get_factors(model)[[1L]], weights = MOFA2::get_weights(model), r2_per_factor = r2,
        r2_total = variance$r2_total[[1L]], elbo = as.numeric(stats_slot$elbo), iterations = length(stats_slot$time),
        messages = messages, model_bytes = file.size(outfile))
}

# The corpus aligned once, trained by the KIT-147 trainer and by this one (MOFA2's run_mofa hand-off; KIT-151, INT-22:
# the direct hand-off is removed); kept for the file.
intMofaHandoffRun <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) return(cache)
        manager <- newIntegrationStateOwner(file.path(tempfile("kit148-trainer-"), "project"))
        on.exit(closeWorkflowStateOwner(manager), add = TRUE)
        importIntegrationInputs(manager, file.path(intMofaCorpus(), intMofaOmics))
        object <- integrationImportedObject(manager)
        resolved <- integrationResolveMofaParameters(object)
        aligned <- integrationAlignViews(object, resolved, prefix_rows = TRUE)
        options <- integrationMofaOptions(resolved)
        cache <<- list(before = intMofaKit147Trainer(aligned$views, options), after = integrationTrainMofa(aligned$views, options))
        cache
    }
})

intMofaWarnings <- function(messages) messages[startsWith(messages, "WARNING:")]

test_that("KIT-148: the trainer without the data in the model file gives KIT-147's model and warnings", {
    intMofaSkip()
    out <- intMofaHandoffRun()
    before <- out$before
    after <- out$after
    for (name in c("factors", "weights", "r2_per_factor", "r2_total", "elbo", "iterations")) {
        expect_identical(after[[name]], before[[name]], info = name)
    }
    message(sprintf("KIT-148 trainer: model file %d bytes without the data, %d bytes with it", as.integer(after$model_file$size),
        as.integer(before$model_bytes)))
    message(paste(c("KIT-148 warnings, KIT-147 trainer:", intMofaWarnings(before$messages), "KIT-148 warnings, this trainer:",
        intMofaWarnings(after$messages)), collapse = "\n"))
    expect_identical(intMofaWarnings(after$messages), intMofaWarnings(before$messages))
    expect_lt(after$model_file$size, before$model_bytes)
})
