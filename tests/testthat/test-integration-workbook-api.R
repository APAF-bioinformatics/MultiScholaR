# KIT-154: the exported integration workbook functions, each a thin wrapper over the module's internals, on the tier-1
# corpus without MOFA2. The model is the stored double-precision reference (tests/testdata/integration/mofa_reference)
# returned by a stand-in trainer (as test-integration-enrichment.R); STRING is stood in, in the client's form (no
# network; the key a sentinel); the report is a stand-in file (no pandoc). Declared before the run:
#   - each function returns what its help page says: the elements named there, tables as data frames, figures as
#     ggplot objects, files that exist in the folders named there;
#   - the trained model's parameters are those the store accepted, and the same settings trained again reuse the model;
#   - the enrichment's key reaches the query and is in no value returned; without a key STRING is skipped in words;
#   - closed and repeated from its record, the analysis gives the same model;
#   - arguments are refused in plain words (class `integration_api_invalid`) before anything runs.
# Fails on the base tree: the functions are absent. With MOFA2, the same twelve functions run end to end in
# test-integration-workbook-mofa.R (MOFA2 lane; KIT-234 replaced the workbook's render by that direct run).

intApiData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intApiStored <- function(name) utils::read.delim(file.path(intApiData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)

intApiTrainer <- function() {
    function(views, options) {
        factors <- intApiStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intApiStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intApiStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

# A project holding the corpus under integration/<omic>/, as each omic's Summary step writes it.
intApiProject <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intApiData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(root, "integration", omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intApiData(), "corpus", omic), full.names = TRUE), dir)))
    }
    root
}

test_that("KIT-154: the workbook functions refuse their arguments in plain words before anything runs", {
    refused <- "integration_api_invalid"
    fake <- list(manager = NULL, store_root = "absent", paths = list(), views = "a")
    fake_model <- list(value = list(), reference = list(), parameters = list())
    expect_error(integrationFindInputs(file.path(tempdir(), "kit154-no-such-project")), "existing project folder",
        class = refused)
    expect_error(integrationFindInputs(c("a", "b")), class = refused)
    expect_error(integrationOpenProject(list()), "integrationFindInputs", class = refused)
    expect_error(integrationMatchSamples(list(manager = NULL)), "integrationOpenProject", class = refused)
    expect_error(integrationMatchSamples(fake, file.path(tempdir(), "kit154-no-such-match.tsv")), "sample match table",
        class = refused)
    expect_error(integrationConfirmMatch(fake, list()), "integrationMatchSamples", class = refused)
    confirmed <- list(parameters = list(), samples = data.frame())
    expect_error(integrationTrainModel(fake, confirmed, input_matrix = "bogus"), "input_matrix must be one of",
        class = refused)
    expect_error(integrationTrainModel(fake, confirmed, overrides = list(nonsense = 1)), "nonsense", class = refused)
    expect_error(integrationTrainModel(fake, confirmed, overrides = list(8L)), "named", class = refused)
    expect_error(integrationTrainModel(fake, confirmed, overrides = list(missing_samples = "drop_incomplete")),
        "fixed by the confirmed sample match", class = refused)
    expect_error(integrationGetResults(fake, list()), "integrationTrainModel", class = refused)
    expect_error(integrationGetResults(fake, fake_model, overrides = list(nope = TRUE)), "nope", class = refused)
    expect_error(integrationExportResults(fake, fake_model, list()), "integrationGetResults", class = refused)
    results <- list(results = list(), factor = "Factor1")
    expect_error(integrationEnrichFactor(fake, results, api_key = 1), "api_key", class = refused)
    expect_error(integrationEnrichFactor(fake, results, id_columns = list("kegg_id")), "id_columns", class = refused)
    expect_error(integrationEnrichFactor(fake, results, overrides = list(api_key = "x")), "api_key", class = refused)
    expect_error(integrationWriteRecord(fake, format = "pdf"), "docx", class = refused)
    expect_error(integrationCloseProject(list()), class = refused)
    expect_error(integrationReplayProject(file.path(tempdir(), "kit154-no-such-store")), "existing integration store",
        class = refused)
    # Every refusal is also the module's refusal class, which the app and the workbook already catch.
    expect_error(integrationCloseProject(list()), class = "integration_refused")
})

test_that("KIT-154: integrationSettings() returns the model, results and enrichment settings tables", {
    settings <- integrationSettings()
    expect_named(settings, c("model", "results", "enrichment", "input_matrix_help"))
    expect_identical(settings$model, integrationMofaParameterTable())
    expect_identical(settings$results, integrationResultsParameterTable())
    expect_named(settings$model, c("parameter", "default", "choices", "settled_by"))
    expect_named(settings$results, c("parameter", "default", "choices", "settled_by"))
    expect_s3_class(settings$enrichment, "data.frame")
    expect_named(settings$enrichment, c("setting", "default", "settled_by"))
    expect_identical(settings$enrichment$setting, names(integrationEnrichmentParameterSpec()))
    expect_true("input_matrix" %in% settings$model$parameter)
    expect_identical(settings$input_matrix_help, integrationGuiTexts()[["help_input_matrix"]])
    expect_length(settings$input_matrix_help, 1L)
})

test_that("KIT-154: the workbook functions, step by step on the corpus, return what their help pages say", {
    root <- intApiProject()
    stand_in <- intApiTrainer()
    real_run <- runIntegrationModel
    real_replay <- replayIntegrationProject
    key <- "kit154-sentinel-not-a-real-key-Hb4x"
    seen <- new.env()
    seen$keys <- character()
    terms <- paste0("category\tterm\tdescription\tenrichment_score\tfalse_discovery_rate\n",
        "Process\tGO:0008150\tbiological_process\t1.4\t0.002\n", "KEGG\thsa04110\tCell cycle\t-1.1\t0.03\n")
    testthat::local_mocked_bindings(
        runIntegrationModel = function(manager, parameters = list(), ..., trainer = NULL) {
            real_run(manager, parameters, ..., trainer = if (is.null(trainer)) stand_in else trainer)
        },
        replayIntegrationProject = function(project_root, scratch, ...) {
            real_replay(project_root, scratch, ..., trainer = stand_in)
        },
        renderIntegrationReport = function(record_path, output_dir, format) {
            path <- file.path(output_dir, paste0("integration_report.", format))
            writeLines("stand-in report", path)
            path
        },
        integrationStringRankedEnrichment = function(ranked, api_key, species, options) {
            seen$keys <- c(seen$keys, api_key)
            text <- integrationStringIdentifiersText(ranked)
            list(service = "STRING", method = "valuesranks_enrichment",
                version = list(string_version = "12.5", string_stable_address = options$host),
                request = list(host = options$host, method = "valuesranks_enrichment", species = as.character(species),
                    ge_fdr = options$ge_fdr, rank_direction = as.integer(options$rank_direction),
                    caller_identity = options$caller_identity, identifiers = nrow(ranked),
                    identifiers_sha256 = digest::digest(text, algo = "sha256", serialize = FALSE), polls = 1L,
                    queried_at = "2026-11-22T00:00:00.000+0000"),
                response = list(submission = list(job_id = "kit154"), status = list(status = "success"), results_tsv = terms),
                terms = integrationStringTermsTable(terms))
        },
        .package = "MultiScholaR")

    # Find the inputs.
    inputs <- integrationFindInputs(root)
    expect_named(inputs, c("project_dir", "paths", "omics", "views", "scan"))
    expect_identical(inputs$project_dir, normalizePath(root, winslash = "/"))
    expect_named(inputs$paths, c("base_dir", "integration_dir", "source_dir", "results_dir", "multiomic_inputs_dir",
        "mofa_tables_dir", "mofa_plots_dir", "study_report_dir"))
    expect_true(dir.exists(inputs$paths$source_dir))
    expect_s3_class(inputs$omics, "data.frame")
    expect_identical(nrow(inputs$omics), 4L)
    expect_s3_class(inputs$views, "data.frame")
    expect_gt(nrow(inputs$views), 0L)

    # Open the store and import.
    project <- integrationOpenProject(inputs)
    closed <- FALSE
    withr::defer(if (!closed) closeWorkflowStateOwner(project$manager))
    expect_named(project, c("manager", "store_root", "paths", "views"))
    expect_true(dir.exists(project$store_root))
    expect_identical(dirname(project$store_root), inputs$paths$source_dir)
    expect_identical(project$paths, inputs$paths)
    expect_type(project$views, "character")
    expect_setequal(project$views, inputs$views$view)

    # Match the samples (an empty sample match file: matched by the saved biological samples).
    match <- integrationMatchSamples(project, "")
    expect_named(match, c("facts", "data_sets", "samples", "preview", "object"))
    expect_type(match$facts, "character")
    expect_gt(length(match$facts), 0L)
    expect_s3_class(match$data_sets, "data.frame")
    expect_identical(nrow(match$data_sets), length(project$views))
    expect_s3_class(match$samples, "data.frame")
    expect_gt(nrow(match$samples), 0L)

    # Confirm the match.
    confirmed <- integrationConfirmMatch(project, match)
    expect_true(all(c("samples", "view_summary", "parameters", "manifests", "file", "confirmed_at", "template_file") %in%
        names(confirmed)))
    expect_null(confirmed$sample_match_file)
    expect_identical(confirmed$template_file, file.path(inputs$paths$multiomic_inputs_dir,
        "integration_sample_match_template.tsv"))
    expect_true(file.exists(confirmed$template_file))
    expect_true(file.exists(confirmed$file))
    expect_identical(dirname(confirmed$file), inputs$paths$multiomic_inputs_dir)
    expect_s3_class(confirmed$samples, "data.frame")
    expect_s3_class(confirmed$confirmed_at, "POSIXct")

    # Train (the stand-in trainer), then the same settings again: the model is reused.
    run <- integrationTrainModel(project, confirmed)
    expect_named(run, c("status", "value", "reference", "parameters", "summary", "variance"))
    expect_identical(run$status, "computed")
    expect_identical(run$parameters, project$manager$.context$snapshot$exports$integration_mofa_model$parameters)
    expect_identical(run$parameters$input_matrix, integrationMofaParameterSpec()$input_matrix$default)
    expect_type(run$summary, "character")
    expect_gt(length(run$summary), 0L)
    expect_s3_class(run$variance, "data.frame")
    expect_identical(integrationTrainModel(project, confirmed)$status, "hit")

    # The results.
    results <- integrationGetResults(project, run)
    expect_named(results, c("messages", "claim_limit", "maps_text", "factor", "tables", "figures", "results"))
    expect_type(results$messages, "character")
    expect_type(results$claim_limit, "character")
    expect_identical(results$maps_text, character())
    expect_identical(results$factor, results$results$shown[[1L]])
    expect_named(results$tables, c("variance", "stability", "association", "correlation"))
    for (table in results$tables) expect_s3_class(table, "data.frame")
    expect_named(results$figures, c("variance", "variance_total", "factors_by_group", "weights", "factor_correlation",
        "factor_pairs", "top_feature_scatter", "top_feature_heatmap", "maps"))
    for (name in c("variance", "variance_total", "factors_by_group", "factor_correlation", "factor_pairs")) {
        expect_s3_class(results$figures[[name]], "ggplot")
    }
    for (name in c("weights", "top_feature_scatter", "top_feature_heatmap")) {
        expect_named(results$figures[[name]], results$results$views)
        for (figure in results$figures[[name]]) expect_s3_class(figure, "ggplot")
    }
    expect_length(results$figures$maps, 0L)

    # Save the results.
    saved <- integrationExportResults(project, run, results)
    expect_named(saved, c("tables", "figures", "tables_dir", "plots_dir"))
    expect_identical(saved$tables_dir, inputs$paths$mofa_tables_dir)
    expect_identical(saved$plots_dir, inputs$paths$mofa_plots_dir)
    expect_gt(length(saved$tables), 0L)
    expect_gt(length(saved$figures), 0L)
    expect_true(all(file.exists(c(saved$tables, saved$figures))))
    expect_true(all(dirname(saved$tables) == saved$tables_dir))

    # Enrich the first factor shown: without a key STRING is skipped in words; with one it reaches the query only.
    none <- integrationEnrichFactor(project, results)
    expect_named(none, c("enrichments", "skipped"))
    expect_length(none$enrichments, 0L)
    expect_true(any(grepl("STRING skipped (no key)", none$skipped, fixed = TRUE)))
    expect_identical(seen$keys, character())
    enriched <- integrationEnrichFactor(project, results, api_key = key)
    expect_named(enriched, c("enrichments", "skipped"))
    expect_setequal(names(enriched$enrichments), c("proteomics string", "transcriptomics string"))
    for (enrichment in enriched$enrichments) {
        expect_named(enrichment, c("summary", "terms"))
        expect_type(enrichment$summary, "character")
        expect_s3_class(enrichment$terms, "data.frame")
        expect_gt(nrow(enrichment$terms), 0L)
    }
    expect_true(all(grepl("skipped (no identifier column named)", enriched$skipped, fixed = TRUE)))
    expect_gt(length(enriched$skipped), 0L)
    expect_identical(seen$keys, rep(key, 2L))
    expect_length(grepRaw(key, serialize(enriched, NULL), fixed = TRUE), 0L)

    # The record and report.
    bundle <- integrationWriteRecord(project, format = "html")
    expect_named(bundle, c("directory", "files", "methods"))
    expect_identical(bundle$directory, inputs$paths$study_report_dir)
    expect_true(all(file.exists(bundle$files)))
    expect_true(all(c("analysis_record.json", "integration_report.html") %in% basename(bundle$files)))
    expect_type(bundle$methods, "character")

    # Close, then repeat from the record.
    store_root <- integrationCloseProject(project)
    closed <- TRUE
    expect_identical(store_root, project$store_root)
    replay <- suppressWarnings(integrationReplayProject(store_root))
    expect_named(replay, c("comparison", "outputs_comparison", "enrichment_comparison"))
    expect_true(isTRUE(replay$comparison$equal))
    expect_s3_class(replay$comparison$tables, "data.frame")
    expect_s3_class(replay$outputs_comparison$tables, "data.frame")
    expect_s3_class(replay$enrichment_comparison, "data.frame")
    expect_identical(nrow(replay$enrichment_comparison), 2L)
    expect_true(any(grepl("^integration-replay-", list.files(dirname(store_root)))))
})
