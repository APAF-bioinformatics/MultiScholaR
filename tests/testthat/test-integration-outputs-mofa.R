# KIT-152 with MOFA2 (MOFA2 lane; audit ENG-22): the further outputs on the model MOFA2 trains on the tier-1 corpus,
# the UMAP and t-SNE maps drawn with uwot and Rtsne (MOFA2's dependencies, present on that lane only), the record and
# the report with them, and the analysis repeated from its record after the project is reopened. Skipped where MOFA2,
# a Python with mofapy2, uwot or Rtsne is absent. Declared before the run:
#   - each map holds every biological sample once, two finite coordinates each, its settings recorded (software
#     version, neighbours or perplexity, seed, one thread); drawn again with the same seed in the same process it is
#     identical, with another seed it differs; at an invalid setting it is refused in words (no stored coordinates are
#     compared: the maps' starts depend on the computer's numerical libraries);
#   - the top features' values equal the stored aligned values;
#   - the saved results hold the maps' tables and figures; the Word report (the package's template) holds the section
#     "Further figures", the maps' rows and their caution; the methods paragraph names uwot and Rtsne with their settings;
#   - repeated from the record of the reopened project, the model is the same within 1e-5 relative (INT-6) and the top
#     features' values and the maps' coordinates are identical (same computer).

intOutMofaData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intOutMofaOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intOutMofaSkip <- function() {
    testthat::skip_if_not(dir.exists(file.path(intOutMofaData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    testthat::skip_if_not_installed("uwot")
    testthat::skip_if_not_installed("Rtsne")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
}

intOutMofaProject <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"),
        study_report_dir = file.path(root, "results_summary", "integration", "Study_report"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intOutMofaOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intOutMofaData(), "corpus", omic), full.names = TRUE), dir)))
    }
    paths
}

test_that("KIT-152: with MOFA2, the outputs and both maps on the corpus model, recorded, reported and repeated from the record", {
    intOutMofaSkip()
    paths <- intOutMofaProject()
    evidence <- Sys.getenv("MULTISCHOLAR_INTEGRATION_EVIDENCE_OUT", unset = "")
    texts <- integrationGuiTexts()
    record_texts <- integrationRecordTexts()
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    withr::defer(if (!is.null(shiny::isolate(bus$state_manager))) closeWorkflowStateOwner(shiny::isolate(bus$state_manager)))
    shiny::testServer(mod_integration_import_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(views = bus$integration_scan$views$view, import = 1)
    })
    manager <- shiny::isolate(bus$state_manager)
    shiny::testServer(mod_integration_alignment_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(confirm = 1)
    })
    shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(run = 1)
        expect_identical(bus$integration_model$status, "computed", info = output$run_status$html)
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(param_umap = TRUE, param_tsne = TRUE)
        results <- session$returned()
        expect_false(inherits(results, "integration_results_refused"), info = paste(unlist(results), collapse = " "))
        samples <- results$factor_values$biological_sample
        value <- withWorkflowStateScope(manager, function() hydrateSavedWorkflowEntity(manager$.context$store,
            bus$integration_model$reference)$value)
        # The top features' values: the stored values given to the model.
        for (view in results$views) {
            aligned <- integrationModelViewTable(manager, value, "aligned_values", view)
            rows <- results$feature_values[results$feature_values$view == view, , drop = FALSE]
            expect_identical(rows$value, unname(as.matrix(aligned[-1L])[cbind(match(rows$feature_id, aligned$feature_id),
                match(rows$biological_sample, names(aligned)[-1L]))]), info = view)
        }
        # Both maps: every sample once, finite coordinates, the settings recorded.
        for (method in c("umap", "tsne")) {
            map <- results$maps[[method]]
            expect_null(map$refusal, info = method)
            expect_identical(map$coordinates$biological_sample, samples, info = method)
            expect_true(all(is.finite(c(map$coordinates$map_1, map$coordinates$map_2))), info = method)
            expect_identical(map$settings$package, c(umap = "uwot", tsne = "Rtsne")[[method]])
            expect_identical(map$settings$version, as.character(utils::packageVersion(map$settings$package)))
            expect_identical(c(map$settings$samples, map$settings$factors, map$settings$seed, map$settings$threads),
                c(20L, 5L, 20261119L, 1L))
            message(sprintf("KIT-152 %s map: %s %s; settings %s; notes: %s", method, map$settings$package, map$settings$version,
                paste(names(map$settings), vapply(map$settings, function(x) paste(format(x), collapse = ","), ""), sep = "=",
                    collapse = "; "), paste(map$notes, collapse = " | ")))
        }
        expect_identical(results$maps$umap$settings$n_neighbors, 15L)
        expect_identical(results$maps$tsne$settings$perplexity, 5)
        expect_identical(results$maps$tsne$settings$theta, 0)
        # Determinism within this process: the same seed the same map. Another seed gives another UMAP map; t-SNE starts
        # from the principal components, so whether its seed changes anything is reported, not assumed.
        p <- results$parameters
        for (method in c("umap", "tsne")) {
            again <- integrationFactorMap(results$factor_values, results$factors, method, p)
            expect_identical(again$coordinates, results$maps[[method]]$coordinates, info = method)
            other <- integrationFactorMap(results$factor_values, results$factors, method, utils::modifyList(p, list(map_seed = 7L)))
            changed <- !identical(other$coordinates, results$maps[[method]]$coordinates)
            message(sprintf("KIT-152 %s map with seed 7 instead of %d: %s", method, p$map_seed,
                if (changed) "different" else "identical"))
            if (method == "umap") expect_true(changed)
        }
        # Refused in words at settings invalid for 20 samples.
        expect_identical(conditionMessage(tryCatch(integrationFactorMap(results$factor_values, results$factors, "tsne",
            utils::modifyList(p, list(tsne_perplexity = 7))), error = function(e) e)),
            sprintf(texts[["map_refused_perplexity"]], "7", 20L, "6.33"))
        expect_identical(conditionMessage(tryCatch(integrationFactorMap(results$factor_values, results$factors, "umap",
            utils::modifyList(p, list(umap_n_neighbors = 20L))), error = function(e) e)),
            sprintf(texts[["map_refused_neighbors"]], 20L, 20L, 19L))
        # The figures build; the caution is on the tab.
        for (plot in list(integrationPlotFactorPairs(results), integrationPlotFactorMap(results, "umap"),
            integrationPlotFactorMap(results, "tsne"), integrationPlotTopFeatureScatter(results, "Factor1", "proteomics"),
            integrationPlotTopFeatureHeatmap(results, "Factor1", "proteomics"))) {
            expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
        }
        expect_match(output$maps_status$html, sprintf(texts[["maps_caution"]], 20L), fixed = TRUE)
        # Saved: the maps' tables and figures with the others.
        session$setInputs(save = 1)
        saved <- bus$integration_results
        expect_false(is.null(saved), info = output$save_status$html)
        expect_true(all(c("top_feature_values.tsv", "map_umap.tsv", "map_tsne.tsv", "map_settings.tsv") %in% saved$tables))
        expect_true(all(c("map_umap.png", "map_tsne.png", "factor_pairs.png") %in% saved$figures))
        message(sprintf("KIT-152 results saved: %d tables, %d figures", length(saved$tables), length(saved$figures)))
    })
    shiny::testServer(mod_integration_record_server, args = list(workflow_data = bus, paths = paths,
        withProgressFn = function(message, value, expr) force(expr)), {
        session$setInputs(save = 1)
        expect_false(is.null(bus$integration_record), info = output$record_status$html)
    })
    saved <- shiny::isolate(bus$integration_record)
    record <- jsonlite::read_json(saved$files[[1L]], simplifyVector = FALSE)
    outputs <- integrationRecordTables(record)$integration_outputs
    message("KIT-152 record, further figures:\n", paste(utils::capture.output(print(outputs, right = FALSE)), collapse = "\n"))
    expect_true(all(unname(record_texts[c("out_map_umap", "out_map_tsne", "out_caution")]) %in% outputs[[record_texts[["col_output"]]]]))
    methods <- saved$methods
    message("KIT-152 methods paragraph:\n", methods)
    for (fragment in c(sprintf("drawn with uwot %s (15 neighbours", as.character(utils::packageVersion("uwot"))),
        sprintf("drawn with Rtsne %s (perplexity 5, exact computation", as.character(utils::packageVersion("Rtsne"))),
        "not evidence of groups", sprintf(record_texts[["methods_top_features"]], 10L))) {
        expect_true(grepl(fragment, methods, fixed = TRUE), info = fragment)
    }
    # The Word report: the section, the maps' rows, the caution.
    unpacked <- withr::local_tempdir()
    utils::unzip(saved$files[[4L]], files = "word/document.xml", exdir = unpacked)
    words <- gsub("<[^>]+>", "", paste(readLines(file.path(unpacked, "word", "document.xml"), warn = FALSE), collapse = ""))
    for (fragment in c(record_texts[["section_outputs"]], record_texts[["out_map_umap"]], record_texts[["out_map_tsne"]],
        "not evidence of groups", record_texts[["out_scatter"]])) {
        expect_true(grepl(fragment, words, fixed = TRUE), info = fragment)
    }
    # Closed, reopened, repeated from the record: the same model, the same outputs on this computer.
    root <- dirname(manager$.context$store$directory)
    closeWorkflowStateOwner(manager)
    replay <- replayIntegrationProject(root, file.path(paths$base_dir, "replay"))
    table <- replay$outputs_comparison$tables
    message(sprintf("KIT-152 replay: model largest relative %.3g (identical %s); outputs: %s", replay$comparison$largest_relative,
        replay$comparison$identical, paste(sprintf("%s %.3g (identical %s)", table$table, table$largest_relative, table$identical),
            collapse = "; ")))
    expect_true(replay$comparison$equal)
    expect_identical(table$table, c("top_feature_values", "map_umap", "map_tsne"))
    expect_true(table$identical[[1L]])
    # The maps are drawn from the repeated model's factors: identical when the model is (INT-6: expected on one computer).
    if (isTRUE(replay$comparison$identical)) expect_true(replay$outputs_comparison$identical)
    if (nzchar(evidence)) {
        out <- file.path(evidence, "kit152")
        dir.create(out, recursive = TRUE, showWarnings = FALSE)
        writeLines(methods, file.path(out, "methods-paragraph.txt"))
        file.copy(saved$files[1:4], out, overwrite = TRUE)
        utils::write.table(table, file.path(out, "replay-outputs-comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        utils::write.table(outputs, file.path(out, "record-further-figures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
})
