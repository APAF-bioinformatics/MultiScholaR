# KIT-063 end to end with MOFA2 (items A and B): import, sample matching with an uploaded match table, training,
# results, through the tabs' server logic. The tier-1 corpus is copied into a project with its transcriptomics input
# written again under other sample names (`S-<name>`), as an omic that names its samples differently saves them; the
# edited table maps them back. Skipped where MOFA2 or a Python with mofapy2 is absent (the package's own CI); run on
# the MOFA2 lane with the bundle's Python. Declared before the run:
#   - without the table the samples cannot be matched (refused in words); with it, the trained model's factors,
#     weights and variance explained equal the stored double-precision result (tests/testdata/integration/
#     mofa_reference) within 1e-5 relative, the tolerance of test-integration-mofa-model.R, and the model records the
#     match;
#   - the results show the planted truth: the factor matched to Z1 has at least 90% of its variance between the groups
#     with an adjusted p-value below 0.001; the factor matched to Z3 explains at least 1% of the variance in
#     proteomics and transcriptomics only; the top 10 features of each planted factor's match hold more of its planted
#     features than their data set as a whole, one-sided hypergeometric p below 0.05, in every data set where it is
#     planted;
#   - MOFA2's warning that a factor follows the samples' overall level (seen in the console on this lane at KIT-140) is
#     kept with the model and shown on the Results tab in plain words;
#   - the results are saved to results/integration/mofa, every file written.

intE2eData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intE2eOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intE2eSkip <- function() {
    testthat::skip_if_not(dir.exists(file.path(intE2eData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intE2eProject <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intE2eOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intE2eData(), "corpus", omic), full.names = TRUE), dir)))
    }
    dir <- file.path(paths$integration_dir, "transcriptomics")
    input <- readIntegrationInputWithMatrices(dir)
    read <- input$manifest
    manifest <- list(omic_type = read$omic_type, workflow = read$workflow, matrix_role = read$matrix_role, scale = read$scale,
        correction = NULL, source_state = list(name = "kit063-renamed"), sample_id_column = read$sample_id_column,
        group_column = read$group_column, technical_replicate_column = read$technical_replicate_column,
        biological_sample_column = read$biological_sample_column)
    samples <- input$samples
    samples[[manifest$biological_sample_column]] <- paste0("S-", samples[[manifest$biological_sample_column]])
    writeIntegrationInput(dir, manifest, samples, lapply(input$views, `[[`, "matrix"), lapply(input$views, `[[`, "features"))
    paths
}

intE2eStored <- function(name) utils::read.delim(file.path(intE2eData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)

test_that("KIT-063: an uploaded match, a MOFA2 model and its results on the corpus, with the planted truth", {
    intE2eSkip()
    paths <- intE2eProject()
    texts <- integrationGuiTexts()
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    object <- integrationImportedObject(manager)
    # The table as downloaded, edited: every transcriptomics run back to its biological sample.
    edited <- integrationSampleMatchTemplate(object)
    rows <- edited$omic == "transcriptomics"
    edited$biological_sample[rows] <- sub("^S-", "", edited$saved_biological_sample[rows])
    upload <- file.path(paths$base_dir, "edited.tsv")
    utils::write.table(edited, upload, sep = "\t", quote = FALSE, row.names = FALSE)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        expect_s3_class(session$returned(), "integration_alignment_refused")
        session$setInputs(match_upload = data.frame(name = "edited.tsv", size = file.size(upload), type = "text/tab-separated-values",
            datapath = upload, stringsAsFactors = FALSE))
        expect_identical(session$returned()$match$report$runs_changed, sum(rows))
        session$setInputs(confirm = 1)
    })
    confirmed <- shiny::isolate(bus$integration_alignment)
    expect_false(is.null(confirmed$parameters$sample_match))
    started <- Sys.time()
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(run = 1)
        expect_identical(bus$integration_model$status, "computed", info = output$run_status$html)
    })
    model <- shiny::isolate(bus$integration_model)
    expect_identical(model$parameters$sample_match, confirmed$parameters$sample_match)
    value <- withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store, model$reference)$value)
    message(sprintf("KIT-063 training with the uploaded match: %d iterations, %.1f s; model sha256 %s", value$elbo$iterations,
        as.numeric(difftime(Sys.time(), started, units = "secs")), value$model_file$sha256))
    # The match restores the corpus: the stored double-precision result within 1e-5 relative.
    tables <- list(factors = value$factors, weights = integrationModelWeights(manager, value), variance_explained = value$variance_explained)
    for (name in names(tables)) {
        stored <- intE2eStored(name)
        numeric <- names(stored)[vapply(stored, is.numeric, logical(1))]
        if (identical(dim(stored), dim(tables[[name]]))) {
            message(sprintf("KIT-063 stored-result comparison with the uploaded match, %s: largest relative difference %.3g", name,
                max(abs(as.matrix(tables[[name]][numeric]) - as.matrix(stored[numeric]))) / max(abs(as.matrix(stored[numeric])))))
        }
        expect_equal(tables[[name]], stored, tolerance = 1e-5, info = name)
    }
    warnings <- value$messages[startsWith(value$messages, "WARNING:")]
    message(paste(c("KIT-063 MOFA2 warnings kept with the model:", warnings), collapse = "\n"))
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        results <- session$returned()
        expect_false(inherits(results, "integration_results_refused"))
        # The planted truth.
        truth <- jsonlite::read_json(file.path(intE2eData(), "corpus", "truth.json"))
        scores <- utils::read.delim(file.path(intE2eData(), "corpus", "truth_scores.tsv"), stringsAsFactors = FALSE)
        loadings <- utils::read.delim(file.path(intE2eData(), "corpus", "truth_loadings.tsv"), stringsAsFactors = FALSE)
        values <- results$factor_values
        scores <- scores[match(values$biological_sample, scores$biological_sample), ]
        planted <- c("Z1", "Z2", "Z3")
        r <- sapply(planted, function(z) vapply(results$factors, function(f) stats::cor(values[[f]], scores[[z]]), numeric(1)))
        rownames(r) <- results$factors
        matched <- stats::setNames(results$factors[apply(abs(r), 2L, which.max)], planted)
        message(sprintf("KIT-063 planted truth: matched %s (abs r %s)", paste(planted, matched, sep = "=", collapse = ", "),
            paste(sprintf("%.3f", apply(abs(r), 2L, max)), collapse = ", ")))
        expect_identical(anyDuplicated(matched), 0L)
        z1 <- results$association[results$association$factor == matched[["Z1"]], ]
        message(sprintf("KIT-063 planted truth: %s against the group R2 %.2f%%, F %.1f, adjusted p %.3g", z1$factor, z1$r2_percent,
            z1$f_statistic, z1$p_adjusted))
        expect_identical(z1$status, "tested")
        expect_gte(z1$r2_percent, 90)
        expect_lt(z1$p_adjusted, 0.001)
        active <- lapply(truth$factors, function(z) unlist(lapply(unlist(z$active), function(omic) unlist(truth$views[[omic]]))))
        z3 <- results$variance[results$variance$factor == matched[["Z3"]], ]
        message(sprintf("KIT-063 planted truth: %s variance explained %s", matched[["Z3"]],
            paste(z3$view, sprintf("%.2f", z3$r2_percent), collapse = "; ")))
        expect_setequal(z3$view[z3$r2_percent >= 1], active$Z3)
        for (z in planted) for (view in active[[z]]) {
            all <- results$weights[results$weights$factor == matched[[z]] & results$weights$view == view, ]
            loaded <- loadings$feature_id[loadings$view == view & loadings[[z]] != 0]
            top <- all[all$rank <= results$parameters$top_features, ]
            k_all <- sum(all$feature_id %in% loaded)
            k_top <- sum(top$feature_id %in% loaded)
            p <- stats::phyper(k_top - 1L, k_all, nrow(all) - k_all, nrow(top), lower.tail = FALSE)
            message(sprintf("KIT-063 planted truth, top features: %s (%s) %s: %d of %d planted (data set %.2f), p %.2g", z,
                matched[[z]], view, k_top, nrow(top), k_all / nrow(all), p))
            expect_gt(k_top / nrow(top), k_all / nrow(all), label = paste(z, view))
            expect_lt(p, 0.05, label = paste(z, view))
        }
        # MOFA2's warning that a factor follows the samples' overall level: kept and shown in words.
        kinds <- vapply(results$advisories, `[[`, "", "kind")
        expect_true("overall_level" %in% kinds, info = paste(warnings, collapse = " | "))
        shown <- gsub("&#39;", "'", output$advisories$html, fixed = TRUE)
        for (advisory in results$advisories) expect_match(shown, advisory$sentence, fixed = TRUE)
        # Saved to the project.
        session$setInputs(save = 1)
        saved <- bus$integration_results
        expect_false(is.null(saved), info = output$save_status$html)
        files <- c(file.path(paths$mofa_tables_dir, saved$tables), file.path(paths$mofa_plots_dir, saved$figures))
        expect_true(all(file.size(files) > 0))
        message(sprintf("KIT-063 results saved: %d tables, %d figures", length(saved$tables), length(saved$figures)))
    })
})
