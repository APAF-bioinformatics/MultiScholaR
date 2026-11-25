# KIT-143 (plan item H; audit ENG-23, SCI-14): the workbook-equivalence proof and, since KIT-234, the twelve exported
# integration functions run end to end, with MOFA2 (this run replaces the shipped workbook's render).
# Skipped where MOFA2 or a Python with mofapy2 is absent (the package's own CI); run on the MOFA2 lane with the bundle's
# Python. Declared before the first run:
#   - Equivalence. Both sides run in this one process, after the maths-library pin (one thread, one code path, as the
#     module trains), on the tier-1 corpus. The reference workbook's side (APAF-bioinformatics/e_shivani_mofa,
#     scripts/mofa_ruvIII.Rmd, chunks 6 to 14) is fed from the corpus's own files, read here with arrow and never through
#     the module: each data set's matrix with its feature ids as row names (chunk 6); brought to the form of the
#     reference workbook's inputs, which had one column per biological sample named by it (each run renamed to its
#     saved biological sample and the runs of one sample averaged over the observed values; the reference workbook has
#     no such step, its files were already in that form); the samples not measured in every data set removed from all
#     (chunk 7's rule); each data set's columns sorted by name (chunk 7); each feature z-scored with t(scale(t(x)))
#     (chunk 9); then create_mofa, get_default_data_options with center_groups = FALSE, get_default_model_options with
#     num_factors = 5, get_default_training_options unchanged, prepare_mofa, run_mofa(use_basilisk = FALSE) (chunks 11
#     to 14). The module side is runIntegrationModel on the same corpus imported, with the settings that state the
#     same: every parameter of integrationMofaParameterSpec() given, feature_selection none, scaling feature_zscore,
#     scale_views FALSE, center_groups FALSE, num_factors 5, convergence_mode fast, maxiter 1000, seed 42 (MOFA2's
#     default, which the reference workbook used without stating it), one training seed, precision single
#     (use_float32 = TRUE, MOFA2's default), internal standards kept, samples not in every data set dropped. MOFA2's
#     defaults are checked to be those values, not assumed.
#     Compared: the options handed to MOFA2 (equal, all.equal: 5 and 5L are one value), the values handed to MOFA2,
#     then the factors, the weights per data set and the variance explained (per factor and data set, and in all).
#     Claim: identical() where the values handed to MOFA2 are identical; the pass criterion is the largest difference
#     relative to the table's largest absolute value at most 1e-5 (the replay's rule), because the two sides compute
#     the z-score by two formulas (scale's root mean square against stats::sd) that may differ in the last bit before
#     MOFA2 rounds to single precision; both the identity and the largest difference are reported.
#     A negative control: the module trained with seed 43 differs from the reference workbook's side beyond 1e-5, so
#     the comparison can fail.
#     The claim: "the module reproduces the reference workbook's MOFA2 configuration exactly on synthetic data"
#     (implementation conformance; nothing about the module's defaults, the collaborator's data or results).
#   - The twelve exported integration functions (integrationFindInputs() to integrationReplayProject()) called in
#     the order and with the arguments of the archived workbook's chunks (KIT-234; it rendered them until then), on the
#     same project holding the corpus, with the default settings and the docx report: they write the results, the
#     record and report, and the repeat from the record equals the saved model within 1e-5 relative; the enrichment
#     runs STRING (stood in, no network) for the protein and gene data sets with the key from STRING_API_KEY, which
#     reaches the query and no project file, and the repeat gives the same ranking and terms from the recorded answer.
# MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT set: the comparison, both sides' settings, the maths libraries and the run's
# record and study parameters are written to its kit143 folder.

intWbmRoot <- function() normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
intWbmCorpus <- function() file.path(intWbmRoot(), "tests", "testdata", "integration", "corpus")
intWbmOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intWbmSkip <- function() {
    testthat::skip_if_not(dir.exists(intWbmCorpus()), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    # The probe pins the maths library before Python starts; both sides train under it.
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intWbmEvidence <- function() {
    out <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    if (!nzchar(out)) return(NULL)
    dir <- file.path(out, "kit143")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir
}

# The reference workbook's inputs from the corpus's files: per data set, in the module's order (omics by name, their
# data sets in manifest order), features x biological samples.
intWbmWorkbookInputs <- function() {
    views <- list()
    for (omic in intWbmOmics) {
        dir <- file.path(intWbmCorpus(), omic)
        manifest <- jsonlite::read_json(file.path(dir, "integration_input.json"))
        samples <- as.data.frame(arrow::read_parquet(file.path(dir, manifest$samples_file)))
        key <- stats::setNames(as.character(samples[[manifest$biological_sample_column]]),
            as.character(samples[[manifest$sample_id_column]]))
        for (view in manifest$views) {
            table <- as.data.frame(arrow::read_parquet(file.path(dir, view$matrix_file)))
            # Chunk 6: column_to_rownames on the feature id.
            values <- as.matrix(table[setdiff(names(table), manifest$feature_id_column)])
            rownames(values) <- as.character(table[[manifest$feature_id_column]])
            # The reference workbook's input form: one column per biological sample, named by it.
            sample_of <- unname(key[colnames(values)])
            levels <- unique(sample_of)
            values <- vapply(levels, function(level) {
                means <- rowMeans(values[, sample_of == level, drop = FALSE], na.rm = TRUE)
                means[is.nan(means)] <- NA_real_
                means
            }, numeric(nrow(values)))
            values <- matrix(values, ncol = length(levels), dimnames = list(as.character(table[[manifest$feature_id_column]]),
                levels))
            views[[view$name]] <- values
        }
    }
    # Chunk 7: a sample missing from one omic is removed from all; each data set's columns sorted by name.
    shared <- Reduce(intersect, lapply(views, colnames))
    views <- lapply(views, function(x) {
        x <- x[, colnames(x) %in% shared, drop = FALSE]
        x[, sort(colnames(x)), drop = FALSE]
    })
    # Chunk 7's lookup table: the n-th column of every data set is one sample.
    lookup <- do.call(cbind, lapply(views, colnames))
    stopifnot(all(lookup == lookup[, 1L]))
    # Chunk 9: each feature z-scored.
    lapply(views, function(x) t(scale(t(x))))
}

# The reference workbook's MOFA2 calls (chunks 11 to 14), as written there; the options it handed to MOFA2 kept.
intWbmWorkbookTrain <- function(mofa_matrix_input) {
    MOFAobject <- MOFA2::create_mofa(mofa_matrix_input)
    data_opts <- MOFA2::get_default_data_options(MOFAobject)
    data_opts$center_groups <- FALSE
    model_opts <- MOFA2::get_default_model_options(MOFAobject)
    model_opts$num_factors <- 5
    train_opts <- MOFA2::get_default_training_options(MOFAobject)
    MOFAobject <- MOFA2::prepare_mofa(object = MOFAobject, data_options = data_opts, model_options = model_opts,
        training_options = train_opts)
    outfile <- tempfile("kit143-workbook-", fileext = ".hdf5")
    on.exit(unlink(outfile), add = TRUE)
    messages <- character()
    model <- withCallingHandlers(MOFA2::run_mofa(MOFAobject, outfile, use_basilisk = FALSE),
        message = function(m) {
            messages <<- c(messages, trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            messages <<- c(messages, paste("WARNING:", trimws(conditionMessage(w))))
            invokeRestart("muffleWarning")
        })
    stopifnot(!any(grepl("Setting use_basilisk = TRUE", messages, fixed = TRUE)))
    variance <- MOFA2::get_variance_explained(model)
    r2 <- variance$r2_per_factor[[1L]]
    if (all(rownames(r2) %in% names(mofa_matrix_input))) r2 <- t(r2)
    list(options = list(data = unclass(data_opts), model = unclass(model_opts), training = unclass(train_opts)),
        factors = MOFA2::get_factors(model)[[1L]], weights = MOFA2::get_weights(model), r2_per_factor = r2,
        r2_total = variance$r2_total[[1L]])
}

# The module with the settings equivalent to the reference workbook's, every parameter of the table given.
intWbmModuleParameters <- function(seed = 42L) {
    spec <- integrationMofaParameterSpec()
    parameters <- lapply(spec, `[[`, "default")
    parameters[c("technical_replicates", "internal_standards", "missing_samples", "feature_selection", "features_per_view",
        "min_observed_fraction", "scaling", "scale_views", "center_groups", "num_factors", "convergence_mode", "maxiter",
        "seed", "training_seeds", "precision")] <- list("average", "keep", "drop_incomplete", "none", NA_integer_, NA_real_,
        "feature_zscore", FALSE, FALSE, 5L, "fast", 1000L, seed, 1L, "single")
    parameters
}

# The module's model read from its store: the options, the values handed to MOFA2, factors, weights, variance.
intWbmModuleTables <- function(manager, value) {
    views <- value$view_summary$view
    matrix_of <- function(table) {
        m <- as.matrix(table[-1L])
        rownames(m) <- as.character(table[[1L]])
        m
    }
    factors <- as.matrix(value$factors[-1L])
    rownames(factors) <- value$factors$biological_sample
    r2 <- value$variance_explained
    r2 <- matrix(r2$r2_percent, nrow = length(unique(r2$factor)), dimnames = list(unique(r2$factor), unique(r2$view)))
    list(options = value$options, view_summary = value$view_summary, samples = value$samples,
        inputs = stats::setNames(lapply(views, function(view) matrix_of(integrationModelViewTable(manager, value,
            "aligned_values", view))), views),
        factors = factors,
        weights = stats::setNames(lapply(views, function(view) matrix_of(integrationModelViewTable(manager, value, "weights",
            view))), views),
        r2_per_factor = r2, r2_total = stats::setNames(value$variance_total$r2_percent, value$variance_total$view))
}

# The replay's rule (integrationCompareTable): the largest absolute difference, relative to the table's largest
# absolute value; and identical() of the numbers. Dimension names are compared apart.
intWbmCompare <- function(name, workbook, module) {
    a <- unname(as.matrix(workbook))
    b <- unname(as.matrix(module))
    storage.mode(a) <- "double"
    storage.mode(b) <- "double"
    same_shape <- identical(dim(a), dim(b))
    same_missing <- same_shape && identical(is.na(a), is.na(b))
    absolute <- if (same_missing) max(c(0, abs(a - b)), na.rm = TRUE) else Inf
    scale <- if (same_missing) max(c(0, abs(b)), na.rm = TRUE) else 0
    relative <- if (!is.finite(absolute)) Inf else if (absolute == 0) 0 else absolute / max(scale, .Machine$double.xmin)
    data.frame(table = name, rows = nrow(b), columns = ncol(b), identical = same_shape && identical(a, b),
        largest_absolute = absolute, largest_relative = relative, within_1e5 = relative <= 1e-5, stringsAsFactors = FALSE)
}

intWbmEquivalence <- function(workbook, module) {
    rows <- list(intWbmCompare("factors", workbook$factors, module$factors))
    for (view in names(module$weights)) {
        rows[[length(rows) + 1L]] <- intWbmCompare(paste0("weights:", view), workbook$weights[[view]], module$weights[[view]])
    }
    rows[[length(rows) + 1L]] <- intWbmCompare("variance_explained", workbook$r2_per_factor[rownames(module$r2_per_factor),
        colnames(module$r2_per_factor), drop = FALSE], module$r2_per_factor)
    rows[[length(rows) + 1L]] <- intWbmCompare("variance_total", workbook$r2_total[names(module$r2_total)],
        module$r2_total)
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

# Both sides once for the file: the reference workbook's, then the module's with seed 42 and with seed 43.
intWbmRun <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) return(cache)
        pin <- integrationPinMathsLibrary()
        inputs <- intWbmWorkbookInputs()
        workbook <- intWbmWorkbookTrain(inputs)
        libraries <- integrationMathsLibraries()
        defaults <- local({
            object <- MOFA2::create_mofa(inputs)
            list(data = MOFA2::get_default_data_options(object), training = MOFA2::get_default_training_options(object))
        })
        root <- tempfile("kit143-module-")
        manager <- newIntegrationStateOwner(root)
        on.exit(closeWorkflowStateOwner(manager), add = TRUE)
        importIntegrationInputs(manager, file.path(intWbmCorpus(), intWbmOmics))
        run <- runIntegrationModel(manager, intWbmModuleParameters(42L))
        module <- intWbmModuleTables(manager, run$value)
        control_run <- runIntegrationModel(manager, intWbmModuleParameters(43L))
        control <- intWbmModuleTables(manager, control_run$value)
        cache <<- list(pin = pin, libraries = libraries, defaults = defaults, inputs = inputs, workbook = workbook,
            module = module, control = control, versions = run$value$versions)
        cache
    }
})

test_that("KIT-143: MOFA2's defaults are the settings the reference workbook trained with", {
    intWbmSkip()
    defaults <- intWbmRun()$defaults
    expect_identical(as.numeric(defaults$training$seed), 42)
    expect_identical(as.numeric(defaults$training$maxiter), 1000)
    expect_identical(defaults$training$convergence_mode, "fast")
    expect_true(isTRUE(defaults$data$use_float32))
    expect_false(isTRUE(defaults$data$scale_views))
})

test_that("KIT-143: the module reproduces the reference workbook's MOFA2 configuration exactly on synthetic data", {
    intWbmSkip()
    out <- intWbmRun()
    workbook <- out$workbook
    module <- out$module
    # The same samples, data sets and features, nothing dropped by the module's own rules on this corpus.
    expect_identical(names(out$inputs), module$view_summary$view)
    expect_true(all(module$view_summary$sparse_dropped == 0L))
    expect_true(all(module$view_summary$constant_dropped == 0L))
    expect_true(all(module$view_summary$internal_standards_dropped == 0L))
    expect_identical(colnames(out$inputs[[1L]]), module$samples$biological_sample)
    expect_identical(rownames(workbook$factors), rownames(module$factors))
    expect_identical(colnames(workbook$factors), colnames(module$factors))
    # The options handed to MOFA2: the same values.
    expect_true(isTRUE(all.equal(workbook$options, module$options, check.attributes = FALSE)),
        info = paste(all.equal(workbook$options, module$options, check.attributes = FALSE), collapse = "; "))
    # The values handed to MOFA2.
    inputs <- do.call(rbind, lapply(names(out$inputs), function(view) {
        expect_identical(rownames(out$inputs[[view]]), rownames(module$inputs[[view]]), info = view)
        expect_identical(colnames(out$inputs[[view]]), colnames(module$inputs[[view]]), info = view)
        intWbmCompare(paste0("input:", view), out$inputs[[view]], module$inputs[[view]])
    }))
    equivalence <- intWbmEquivalence(workbook, module)
    control <- intWbmEquivalence(workbook, out$control)
    for (view in names(module$weights)) {
        expect_identical(sub(paste0("_", view, "$"), "", rownames(workbook$weights[[view]])), rownames(module$weights[[view]]),
            info = view)
    }
    show <- function(table) paste(sprintf("%s: identical %s, largest %.3g (relative %.3g)", table$table, table$identical,
        table$largest_absolute, table$largest_relative), collapse = "\n")
    message("KIT-143 values handed to MOFA2, reference workbook against module:\n", show(inputs))
    message("KIT-143 trained model, reference workbook against module (seed 42):\n", show(equivalence))
    message("KIT-143 negative control, reference workbook against module with seed 43:\n", show(control))
    message(sprintf("KIT-143 maths pin applied %s; libraries %s; MOFA2 %s, mofapy2 %s", out$pin$applied,
        paste(vapply(out$libraries, function(l) sprintf("%s %s threads %s %s", l$api, l$version, l$threads, l$architecture),
            character(1)), collapse = "; "), out$versions$MOFA2, out$versions$mofapy2))
    evidence <- intWbmEvidence()
    if (!is.null(evidence)) {
        utils::write.table(rbind(cbind(comparison = "inputs", inputs), cbind(comparison = "seed 42", equivalence),
            cbind(comparison = "control seed 43", control)), file.path(evidence, "equivalence.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
        jsonlite::write_json(list(workbook = workbook$options, module = module$options, maths_pin = out$pin,
            maths_libraries = out$libraries, versions = out$versions), file.path(evidence, "settings.json"),
            auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", na = "null")
    }
    expect_true(all(inputs$within_1e5), info = show(inputs))
    expect_true(all(equivalence$within_1e5), info = show(equivalence))
    expect_false(control$within_1e5[control$table == "factors"], info = show(control))
})

test_that("KIT-143 (KIT-234): the twelve integration functions run end to end on the corpus with MOFA2 and repeat the model from the record", {
    intWbmSkip()
    # The record's integration_report.docx is rendered by pandoc.
    testthat::skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
    project <- withr::local_tempdir()
    for (omic in intWbmOmics) {
        dir <- file.path(project, "integration", omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intWbmCorpus(), omic), full.names = TRUE), dir)))
    }
    scripts <- file.path(project, "scripts", "integration")
    dir.create(scripts, recursive = TRUE)
    # KIT-141 (INT-21): STRING stood in, in the client's form (no network); the key a sentinel from STRING_API_KEY.
    key <- "kit143-lane-sentinel-not-a-real-key-Qz7w"
    withr::local_envvar(STRING_API_KEY = key)
    seen <- new.env()
    seen$keys <- character()
    terms <- paste0("category\tterm\tdescription\tenrichment_score\tfalse_discovery_rate\n",
        "Process\tGO:0008150\tbiological_process\t1.4\t0.002\n", "KEGG\thsa04110\tCell cycle\t-1.1\t0.03\n")
    testthat::local_mocked_bindings(integrationStringRankedEnrichment = function(ranked, api_key, species, options) {
        seen$keys <- c(seen$keys, api_key)
        text <- integrationStringIdentifiersText(ranked)
        list(service = "STRING", method = "valuesranks_enrichment",
            version = list(string_version = "12.5", string_stable_address = options$host),
            request = list(host = options$host, method = "valuesranks_enrichment", species = as.character(species),
                ge_fdr = options$ge_fdr, rank_direction = as.integer(options$rank_direction),
                caller_identity = options$caller_identity, identifiers = nrow(ranked),
                identifiers_sha256 = digest::digest(text, algo = "sha256", serialize = FALSE), polls = 1L,
                queried_at = "2026-11-21T00:00:00.000+0000"),
            response = list(submission = list(job_id = "lane"), status = list(status = "success"), results_tsv = terms),
            terms = integrationStringTermsTable(terms))
    }, .package = "MultiScholaR")
    # The archived workbook's chunks (integration_workflow_starter.rmd:58-286 at 67569959), with its render's
    # parameters (project_dir, sample_match_file "", report_format "docx").
    inputs <- integrationFindInputs(project)
    project_obj <- integrationOpenProject(inputs)
    match <- integrationMatchSamples(project_obj, "")
    confirmed <- integrationConfirmMatch(project_obj, match)
    settings <- integrationSettings()
    input_matrix <- settings$model$default[settings$model$parameter == "input_matrix"]
    run <- integrationTrainModel(project_obj, confirmed, input_matrix = input_matrix, overrides = list())
    results <- integrationGetResults(project_obj, run, overrides = list())
    saved <- integrationExportResults(project_obj, run, results)
    enriched <- integrationEnrichFactor(project_obj, results, factor = results$factor,
        api_key = Sys.getenv("STRING_API_KEY"), id_columns = list(), overrides = list())
    bundle <- integrationWriteRecord(project_obj, format = "docx")
    store_root <- integrationCloseProject(project_obj)
    replay <- integrationReplayProject(store_root)
    expect_identical(run$status, "computed")
    expect_gt(length(saved$tables), 0L)
    expect_gt(length(saved$figures), 0L)
    expect_true(all(file.exists(bundle$files)))
    expect_identical(basename(bundle$files), c("analysis_record.json", "analysis_record.xlsx", "study_parameters.txt",
        "integration_report.docx", "model.hdf5", "sample_match.tsv"))
    comparison <- replay$comparison
    message(sprintf("KIT-143 (KIT-234) direct run; its repeat from the record: equal %s, identical %s, largest relative %.3g",
        comparison$equal, comparison$identical, comparison$largest_relative))
    expect_true(isTRUE(comparison$equal))
    # The enrichment: STRING for the protein and gene data sets with the key from the environment, KEGG and Reactome
    # skipped (no identifier column named); recorded and repeated offline from the recorded answer.
    expect_setequal(names(enriched$enrichments), c("proteomics string", "transcriptomics string"))
    expect_identical(seen$keys, rep(key, 2L))
    replayed <- replay$enrichment_comparison
    message("KIT-143 enrichment repeated from the record:\n", paste(utils::capture.output(print(replayed)), collapse = "\n"))
    expect_identical(nrow(replayed), 2L)
    expect_true(all(replayed$ranking_equal) && all(replayed$terms_equal))
    record <- jsonlite::read_json(bundle$files[[1L]], simplifyVector = FALSE)
    expect_length(record$analyses$integration_mofa_model$integration$enrichment, 2L)
    files <- list.files(project, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    holding <- files[vapply(files, function(path) length(grepRaw(key, readBin(path, "raw", file.info(path)$size),
        fixed = TRUE)) > 0L, logical(1))]
    expect_identical(holding, character())
    evidence <- intWbmEvidence()
    if (!is.null(evidence)) file.copy(bundle$files[c(1L, 3L)], evidence, overwrite = TRUE)
})
