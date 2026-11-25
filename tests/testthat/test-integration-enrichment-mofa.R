# KIT-141 end to end with MOFA2: import, sample matching with an uploaded match table (the transcriptomics samples
# renamed `S-<name>`, as test-integration-record-mofa.R), training, results saved, the enrichment run through the
# Enrichment tab with a stand-in STRING (no network; an invented key), the record and the report saved, the project
# closed, reopened and the analysis repeated from its record, the enrichment offline from the recorded answer. Skipped
# where MOFA2 or a Python with mofapy2 is absent; run on the MOFA2 lane. Declared before the run:
#   - every measured transcriptomics feature is ranked over the model's biological samples with the uploaded match
#     applied (at least 90% of the features ranked; without the match none would be);
#   - the key reaches the stand-in and nothing else: not the head, the record bundle or any file of the project;
#   - the record's enrichment sheet names STRING 12.5, the report holds the enrichment section and its terms;
#   - the repeated analysis gives the same ranking (within 1e-5 relative) and the same terms from the recorded answer.

intEnrMofaData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intEnrMofaKey <- "kit141-lane-sentinel-not-a-real-key-Hb4m"

intEnrMofaSkip <- function() {
    testthat::skip_if_not(dir.exists(file.path(intEnrMofaData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intEnrMofaProject <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        mofa_tables_dir = file.path(results_dir, "mofa", "tables"), mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intEnrMofaData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(paths$integration_dir, "transcriptomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit141-renamed"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    samples <- input$samples
    samples[[manifest$biological_sample_column]] <- paste0("S-", samples[[manifest$biological_sample_column]])
    writeIntegrationInput(dir, manifest, samples, lapply(input$views, `[[`, "matrix"), lapply(input$views, `[[`, "features"))
    paths
}

intEnrMofaTerms <- paste0("category\tterm\tdescription\tenrichment_score\tfalse_discovery_rate\n",
    "Process\tGO:0008150\tbiological_process\t1.4\t0.002\n", "KEGG\thsa04110\tCell cycle\t-1.1\t0.03\n")

# STRING stood in, in the client's form; it notes the keys it was given.
intEnrMofaQuery <- function(seen) {
    function(ranked, api_key, species, options) {
        seen$keys <- c(seen$keys, api_key)
        text <- integrationStringIdentifiersText(ranked)
        list(service = "STRING", method = "valuesranks_enrichment",
            version = list(string_version = "12.5", string_stable_address = options$host),
            request = list(host = options$host, method = "valuesranks_enrichment", species = as.character(species),
                ge_fdr = options$ge_fdr, rank_direction = as.integer(options$rank_direction),
                caller_identity = options$caller_identity, identifiers = nrow(ranked),
                identifiers_sha256 = digest::digest(text, algo = "sha256", serialize = FALSE), polls = 1L,
                queried_at = "2026-11-21T00:00:00.000+0000"),
            response = list(submission = list(job_id = "lane"), status = list(status = "success"), results_tsv = intEnrMofaTerms),
            terms = integrationStringTermsTable(intEnrMofaTerms))
    }
}

test_that("KIT-141: with MOFA2, a factor is enriched through the tab, recorded, reported and repeated offline", {
    intEnrMofaSkip()
    paths <- intEnrMofaProject()
    no_progress <- function(message, value, expr) force(expr)
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    edited <- integrationSampleMatchTemplate(integrationImportedObject(manager))
    rows <- edited$omic == "transcriptomics"
    edited$biological_sample[rows] <- sub("^S-", "", edited$saved_biological_sample[rows])
    upload <- file.path(paths$base_dir, "edited.tsv")
    utils::write.table(edited, upload, sep = "\t", quote = FALSE, row.names = FALSE)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload), type = "text/tab-separated-values",
            datapath = upload, stringsAsFactors = FALSE))
        session$setInputs(confirm = 1)
    })
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths, withProgressFn = no_progress), {
        session$setInputs(run = 1)
        expect_identical(bus$integration_model$status, "computed", info = output$run_status$html)
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(save = 1)
        expect_false(is.null(bus$integration_results), info = output$save_status$html)
    })
    seen <- new.env()
    seen$keys <- character()
    shiny::testServer(mod_integration_enrichment_server, args = list(workflow_data = bus, paths = paths,
        query = intEnrMofaQuery(seen), withProgressFn = no_progress), {
        session$flushReact()
        session$setInputs(factor = "Factor1", view = "transcriptomics", database = "string", correlation_method = "pearson",
            min_samples = 5L, seed = 20261121L, string_species = 9606L, string_fdr = 0.05, kegg_organism = "hsa",
            reactome_organism = "Homo sapiens", min_set_size = 5L, max_set_size = 500L, string_key = intEnrMofaKey)
        session$setInputs(run = 1)
        expect_false(is.null(bus$integration_enrichment), info = output$result$html)
        expect_identical(bus$tab_status$integration_enrichment, "complete")
    })
    expect_identical(seen$keys, intEnrMofaKey)
    export <- manager$.context$snapshot$exports$integration_mofa_model
    entry <- export$enrichment[["Factor1|transcriptomics|string"]]
    value <- withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store, entry$reference)$value)
    ranking <- withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store, entry$ranking)$value)
    object <- integrationImportedObject(manager)
    expect_identical(nrow(ranking$ranking), nrow(object@views$transcriptomics))
    expect_gte(ranking$counts$ranked, 0.9 * ranking$counts$features)
    expect_identical(value$terms$term, c("GO:0008150", "hsa04110"))
    expect_false(grepl(intEnrMofaKey, paste(deparse(manager$.context$snapshot), collapse = ""), fixed = TRUE))
    # The record and the report (the package's template, as html).
    bundle <- integrationWriteRecordBundle(manager, paths$study_report_dir, format = "html")
    record <- jsonlite::read_json(bundle$files[[1L]], simplifyVector = FALSE)
    texts <- integrationRecordTexts()
    tables <- integrationRecordTables(record)
    expect_match(tables$integration_enrichment[[texts[["col_enrich_release"]]]], "STRING 12.5", fixed = TRUE)
    report <- bundle$files[[4L]]
    html <- paste(readLines(report, warn = FALSE), collapse = "\n")
    expect_match(html, texts[["section_enrichment"]], fixed = TRUE)
    expect_match(html, "GO:0008150", fixed = TRUE)
    files <- list.files(paths$base_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    holding <- files[vapply(files, function(path) length(grepRaw(intEnrMofaKey, readBin(path, "raw", file.info(path)$size),
        fixed = TRUE)) > 0L, logical(1))]
    expect_identical(holding, character())
    # Reopened and repeated from the record with MOFA2; the enrichment offline from the recorded answer.
    root <- dirname(manager$.context$store$directory)
    closeWorkflowStateOwner(manager)
    bus$state_manager <- NULL
    replay <- suppressWarnings(replayIntegrationProject(root, file.path(paths$base_dir, "replay"),
        enrichment_query = function(...) stop("STRING was queried during an offline replay")))
    expect_true(replay$comparison$equal)
    comparison <- replay$enrichment_comparison
    expect_identical(comparison$enrichment, "Factor1|transcriptomics|string")
    expect_true(comparison$ranking_equal)
    expect_true(comparison$terms_identical)
    evidence <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    if (nzchar(evidence)) {
        dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
        utils::write.table(comparison, file.path(evidence, "kit141-enrichment-replay.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
        writeLines(integrationMethodsText(record), file.path(evidence, "kit141-methods.txt"))
        file.copy(report, file.path(evidence, "kit141-integration_report.html"), overwrite = TRUE)
    }
})

# KIT-183 T14: the test of record on a model trained with MOFA2, through runIntegrationEnrichment(), recorded and
# repeated offline from the record (stand-ins that fail on any network call). Declared before the run: the KEGG entry
# is camera's (`Factor1|metabolomics_GCMS|kegg|camera`), its terms carry the camera columns, the record's Methods name
# camera, and the repeated analysis replays camera with terms equal within the tolerance.
test_that("KIT-183 T14: with MOFA2, camera through the store, recorded and repeated offline", {
    intEnrMofaSkip()
    skip_if_not_installed("clusterProfiler")
    root <- withr::local_tempdir()
    integration_dir <- file.path(root, "integration")
    for (omic in c("lipidomics", "metabolomics", "proteomics", "transcriptomics")) {
        dir <- file.path(integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intEnrMofaData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(integration_dir, "metabolomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit183-kegg"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    features <- lapply(input$views, `[[`, "features")
    ids <- sprintf("C%05d", 10000L + seq_len(nrow(features$metabolomics_GCMS)))
    features$metabolomics_GCMS$kegg_id <- ids
    writeIntegrationInput(dir, manifest, input$samples, lapply(input$views, `[[`, "matrix"), features)
    manager <- newIntegrationStateOwner(file.path(root, "store"))
    opened <- TRUE
    withr::defer(if (opened) closeWorkflowStateOwner(manager))
    scan <- integrationScanProject(integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, list())
    expect_identical(run$status, "computed")
    sets <- split(ids, rep(c("first", "second", "third", "rest"), c(10, 10, 10, length(ids) - 30L)))
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    mapping <- list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-23T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
    out <- runIntegrationEnrichment(manager, list(factor = "Factor1", view = "metabolomics_GCMS", database = "kegg",
        id_column = "kegg_id"), mapping = mapping)
    expect_identical(out$run$value$test, "camera")
    expect_true(all(c("Direction", "correlation", "effect", "pvalue", "p.adjust", "test") %in% names(out$run$value$terms)))
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(names(export$enrichment), "Factor1|metabolomics_GCMS|kegg|camera")
    record <- buildWorkflowAnalysisRecord(manager)
    texts <- integrationRecordTexts()
    expect_true(grepl(texts[["methods_test_name_camera"]], integrationMethodsText(record), fixed = TRUE))
    store_root <- dirname(manager$.context$store$directory)
    closeWorkflowStateOwner(manager)
    opened <- FALSE
    offline <- function(...) stop("a service was called during an offline replay")
    replay <- suppressWarnings(replayIntegrationProject(store_root, file.path(root, "replay"), enrichment_query = offline,
        enrichment_download = offline))
    comparison <- replay$enrichment_comparison
    expect_identical(comparison$enrichment, "Factor1|metabolomics_GCMS|kegg|camera")
    expect_identical(comparison$note, "")
    expect_true(comparison$ranking_equal)
    expect_true(comparison$terms_equal)
    evidence <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    if (nzchar(evidence)) {
        dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
        utils::write.table(comparison, file.path(evidence, "kit183-camera-replay.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
    }
})

# KIT-184 K13: a model MOFA2 trained on RUV-corrected values; the correlation-adjusted test with the removed factors in
# its design through the store, recorded, and repeated offline from the record. The corrected matrices are made here as
# in ticket 6.6: the corpus's proteomics and metabolomics matrices less two known directions over the runs (the parity
# of the replicate number, and the replicate number) times known per-feature coefficients, every second feature's level
# shifted by 0.5; the correction declares k = 2.
test_that("KIT-184 K13: with MOFA2 trained on corrected values, camera with the removed factors, recorded and repeated offline", {
    intEnrMofaSkip()
    skip_if_not_installed("clusterProfiler")
    root <- withr::local_tempdir()
    integration_dir <- file.path(root, "integration")
    ids <- NULL
    for (omic in c("metabolomics", "proteomics")) {
        dir <- file.path(integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intEnrMofaData(), "corpus", omic), full.names = TRUE), dir)))
        input <- readIntegrationInputWithMatrices(dir)
        read <- input$manifest
        manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
            correction = list(state = "ruv_corrected", method = "kit184-test", k = 2L, replicate_structure = "group",
                grouping_variable = "group"), source_state = list(name = "kit184-mofa"),
            sample_id_column = read$sample_id_column, group_column = read$group_column,
            technical_replicate_column = read$technical_replicate_column,
            biological_sample_column = read$biological_sample_column)
        features <- lapply(input$views, `[[`, "features")
        if (omic == "metabolomics") {
            ids <- sprintf("C%05d", 10000L + seq_len(nrow(features$metabolomics_GCMS)))
            features$metabolomics_GCMS$kegg_id <- ids
        }
        replicate <- as.integer(sub("^.*_", "", as.character(input$samples$tech_rep_group)))
        w <- cbind(ifelse(replicate %% 2L == 0L, 1, -1), replicate - 5.5)
        matrices <- lapply(input$views, `[[`, "matrix")
        corrected <- lapply(matrices, function(table) {
            m <- as.matrix(table[-1L])
            i <- seq_len(nrow(m))
            a <- cbind(0.30 + 0.01 * (i %% 13L), 0.02 + 0.05 * ((i %% 5L) - 2))
            m <- m - a %*% t(w) - ifelse(i %% 2L == 0L, 0.5, 0)
            out <- data.frame(feature_id = table$feature_id, m, check.names = FALSE, stringsAsFactors = FALSE)
            names(out) <- names(table)
            out
        })
        writeIntegrationInput(dir, manifest, input$samples, matrices, features, corrected = corrected)
        expect_true(all(file.exists(file.path(dir, sprintf("removed_%s.parquet", names(matrices))))))
    }
    manager <- newIntegrationStateOwner(file.path(root, "store"))
    opened <- TRUE
    withr::defer(if (opened) closeWorkflowStateOwner(manager))
    scan <- integrationScanProject(integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, list(input_matrix = "ruv_corrected"))
    expect_identical(run$status, "computed")
    sets <- split(ids, rep(c("first", "second", "third", "rest"), c(10, 10, 10, length(ids) - 30L)))
    pathway <- paste0("path:hsa", sprintf("%05d", seq_along(sets)))
    mapping <- list(database = "KEGG", release = "KEGG Release fixture", sha256 = strrep("e", 64L), url = "fixture",
        downloaded_at = "2026-11-24T00:00:00+0000",
        table = data.frame(pathway = rep(pathway, lengths(sets)), compound = paste0("cpd:", unlist(sets, use.names = FALSE)),
            name = rep(names(sets), lengths(sets)), stringsAsFactors = FALSE))
    # The adjustment is off unless chosen (VAL-160): the request switches it on.
    out <- runIntegrationEnrichment(manager, list(factor = "Factor1", view = "metabolomics_GCMS", database = "kegg",
        id_column = "kegg_id", adjust_ruv_factors = TRUE), mapping = mapping)
    value <- out$run$value
    expect_identical(value$test, "camera")
    expect_identical(value$design_columns, c("intercept", "factor", "ruv_W1", "ruv_W2"))
    expect_identical(value$design_detail$status, "exact")
    expect_identical(value$design_detail$residual_df, value$test_counts$samples - 4L)
    expect_true(out$entry$settings$adjust_ruv_factors)
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(names(export$enrichment), "Factor1|metabolomics_GCMS|kegg|camera|adjusted")
    record <- buildWorkflowAnalysisRecord(manager)
    expect_true(grepl("unwanted factors that the RUV correction removed from this data set", integrationMethodsText(record),
        fixed = TRUE))
    store_root <- dirname(manager$.context$store$directory)
    closeWorkflowStateOwner(manager)
    opened <- FALSE
    offline <- function(...) stop("a service was called during an offline replay")
    replay <- suppressWarnings(replayIntegrationProject(store_root, file.path(root, "replay"), enrichment_query = offline,
        enrichment_download = offline))
    comparison <- replay$enrichment_comparison
    expect_identical(comparison$enrichment, "Factor1|metabolomics_GCMS|kegg|camera|adjusted")
    expect_identical(comparison$note, "")
    expect_true(comparison$ranking_equal)
    expect_true(comparison$terms_equal)
    evidence <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    if (nzchar(evidence)) {
        dir.create(evidence, recursive = TRUE, showWarnings = FALSE)
        utils::write.table(comparison, file.path(evidence, "kit184-covariate-replay.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
    }
})
