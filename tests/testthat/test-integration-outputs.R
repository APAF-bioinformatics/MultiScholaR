# KIT-152: the integration Results' further outputs, without MOFA2: the factor pair plot, the top features' scatter
# and heatmap, the UMAP and t-SNE maps' settings and refusals, the Results tab's settings built from their parameter
# table, and the outputs' place in the saved results, the record, the report's tables, the methods paragraph and the
# replay. The model is the stored double-precision reference (tests/testdata/integration/mofa_reference) returned by a
# stand-in trainer on the tier-1 corpus imported into a store, as test-integration-results.R; a second store holds 70
# biological samples and 10,050 features (audit ENG-13). Declared before the run:
#   - the top features' values are the stored aligned values (the values given to the model), equal to an independent
#     whole read of that table, read by row and in blocks of at most 64 samples and 10,000 rows; the heatmap centres
#     (and scales) them for display only; the y axis names the scale (log2, log2 counts per million, z-score);
#   - a map is refused in plain words, with the class `integration_map_invalid`, at a number of neighbours outside 2 to
#     the samples minus 1, at a perplexity whose triple is not below the samples minus 1, with fewer than two factors,
#     and without its software; a refused map leaves the other results shown and stops the save;
#   - the saved results carry the new tables and figures; the head keeps plain lists and table references only, so the
#     project reopens; the record, its tables, the study parameters and the methods paragraph describe the outputs; a
#     replay from the reopened project gives the same top features' values.
# The maps themselves (uwot, Rtsne: MOFA2's dependencies, absent on the run-tests runner) are drawn, repeated and
# reported in test-integration-outputs-mofa.R on the MOFA2 lane.

intOutData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intOutOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")

intOutStored <- function(name) utils::read.delim(file.path(intOutData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)

# The stored reference as `integrationTrainMofa` returns a model.
intOutTrainer <- function() {
    function(views, options) {
        factors <- intOutStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intOutStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intOutStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = character(),
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

intOutProject <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intOutData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intOutOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intOutData(), "corpus", omic), full.names = TRUE), dir)))
    }
    paths
}

# The corpus imported into a store and the stored model accepted on it, as the Model tab leaves the bus.
intOutModel <- function(paths, close = TRUE, .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("store-", basename(tempfile()))))
    if (close) withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    # The values as saved, so that their scale is the omics' (KIT-151: the z-score is the default since INT-24).
    run <- runIntegrationModel(manager, list(scaling = "none"), trainer = intOutTrainer())
    export <- manager$.context$snapshot$exports$integration_mofa_model
    list(manager = manager, run = run, object = integrationImportedObject(manager),
        model = list(reference = run$reference, parameters = export$parameters, status = run$status, trained_at = Sys.time()))
}

intOutHtml <- function(tag) as.character(htmltools::renderTags(tag)$html)

# A factor table of 20 biological samples in two groups (for the maps' refusals; no model needed).
intOutFactorValues <- function(k = 5L) {
    values <- matrix(withr::with_seed(1152L, stats::rnorm(20L * k)), 20L, k, dimnames = list(NULL, paste0("Factor", seq_len(k))))
    data.frame(biological_sample = sprintf("S%02d", 1:20), group = rep(c("A", "B"), each = 10L), values,
        stringsAsFactors = FALSE)
}

test_that("KIT-152: the Results parameters gain the outputs' settings, each with a type, and the tab's inputs are built from them", {
    texts <- integrationGuiTexts()
    spec <- integrationResultsParameterSpec()
    expect_identical(vapply(spec, integrationParameterType, ""), c(factors_shown = "factors", top_features = "whole",
        weight_scale = "choice", design_covariate = "choice", association_test = "choice", p_adjust_method = "choice",
        figure_format = "choice", figure_width = "number", figure_height = "number", figure_dpi = "whole",
        feature_values = "choice", heatmap_scaling = "choice", umap = "flag", umap_n_neighbors = "whole", tsne = "flag",
        tsne_perplexity = "number", map_seed = "whole"))
    expect_true(all(grepl("^VAL-221 \\(decided", integrationResultsParameterTable()$settled_by)))
    factors <- paste0("Factor", 1:5)
    p <- integrationResolveResultsParameters(list(), factors)
    expect_identical(p[c("feature_values", "heatmap_scaling", "umap", "umap_n_neighbors", "tsne", "tsne_perplexity", "map_seed")],
        list(feature_values = "model_input", heatmap_scaling = "centre", umap = FALSE, umap_n_neighbors = 15L, tsne = FALSE,
            tsne_perplexity = 5, map_seed = 20261119L))
    for (bad in list(list(umap = "yes"), list(tsne = NA), list(tsne_perplexity = 0), list(umap_n_neighbors = 2.5),
        list(heatmap_scaling = "zscore"), list(feature_values = "raw"))) {
        expect_error(integrationResolveResultsParameters(bad, factors), class = "integration_parameters_invalid", info = names(bad))
    }
    # The inputs' values in the parameters' forms; nothing drawn yet is the default.
    expect_identical(integrationModelSettings(list(), spec, exclude = character()), lapply(spec, `[[`, "default"))
    typed <- integrationModelSettings(list(factors_shown = c("Factor1", "Factor3"), figure_width = 7.5, top_features = 5,
        umap = TRUE, tsne_perplexity = 4), spec, exclude = character())
    expect_identical(typed[c("factors_shown", "figure_width", "top_features", "umap", "tsne_perplexity")],
        list(factors_shown = c("Factor1", "Factor3"), figure_width = 7.5, top_features = 5L, umap = TRUE, tsne_perplexity = 4))
    expect_identical(integrationModelSettings(list(factors_shown = c("all", "Factor2")), spec, exclude = character())$factors_shown,
        "all")
    refusal <- tryCatch(integrationModelSettings(list(tsne_perplexity = NA), spec, exclude = character()), error = function(e) e)
    expect_s3_class(refusal, "integration_parameters_invalid")
    expect_identical(conditionMessage(refusal), sprintf(texts[["setting_not_number"]], texts[["param_tsne_perplexity"]]))
    expect_error(integrationModelSettings(list(umap_n_neighbors = 4.5), spec, exclude = character()),
        class = "integration_parameters_invalid", regexp = texts[["param_umap_n_neighbors"]], fixed = TRUE)
    # The tab: one input per parameter, in words, with the table's default.
    ui <- intOutHtml(mod_integration_results_ui("r"))
    for (name in names(spec)) {
        expect_match(ui, sprintf('data-testid="integration-results-param-%s"', name), fixed = TRUE)
        expect_match(ui, texts[[paste0("param_", name)]], fixed = TRUE)
    }
    expect_match(ui, 'id="r-param_factors_shown"[^>]*multiple')
    expect_match(ui, 'id="r-param_umap_n_neighbors"[^>]*value="15"')
    expect_match(ui, 'id="r-param_figure_width"[^>]*value="8"')
    expect_match(ui, sprintf('name="r-param_heatmap_scaling" value="centre" checked="checked"'), fixed = TRUE)
    expect_match(ui, 'id="r-param_umap" type="checkbox"', fixed = TRUE)
    expect_false(grepl('id="r-param_umap" type="checkbox"[^>]*checked', ui))
    for (id in c("integration-results-pairs-plot", "integration-results-scatter-plot", "integration-results-heatmap-plot",
        "integration-results-maps-status", "integration-results-umap-plot", "integration-results-tsne-plot")) {
        expect_match(ui, sprintf('data-testid="%s"', id), fixed = TRUE)
    }
    # The Model tab's inputs are those of before (one more argument each, all defaulted).
    model_ui <- intOutHtml(mod_integration_model_ui("m"))
    expect_false(grepl("integration-results-param-", model_ui, fixed = TRUE))
    expect_match(model_ui, 'id="m-param_num_factors"[^>]*value="5"')
})

test_that("KIT-152: the top features' values are the stored values given to the model; the heatmap centres them for display only", {
    paths <- intOutProject()
    fixture <- intOutModel(paths)
    manager <- fixture$manager
    value <- fixture$run$value
    texts <- integrationGuiTexts()
    results <- integrationResults(value, fixture$object, manager = manager)
    fv <- results$feature_values
    expect_identical(names(fv), c("factor", "view", "rank", "feature_id", "display_name", "weight", "biological_sample", "group",
        "factor_value", "value"))
    samples <- value$factors$biological_sample
    expect_identical(nrow(fv), nrow(results$top_weights) * length(samples))
    # Against an independent whole read of each data set's aligned table (the values handed to MOFA2).
    for (view in results$views) {
        aligned <- integrationModelViewTable(manager, value, "aligned_values", view)
        rows <- fv[fv$view == view, , drop = FALSE]
        expected <- as.matrix(aligned[-1L])[cbind(match(rows$feature_id, aligned$feature_id), match(rows$biological_sample,
            names(aligned)[-1L]))]
        expect_identical(rows$value, unname(expected), info = view)
    }
    factor_matrix <- as.matrix(value$factors[-1L])
    expect_identical(fv$factor_value, unname(factor_matrix[cbind(match(fv$biological_sample, samples),
        match(fv$factor, colnames(factor_matrix)))]))
    top <- results$top_weights
    expect_identical(as.list(unique(fv[c("factor", "view", "rank", "feature_id")])),
        as.list(top[c("factor", "view", "rank", "feature_id")]))
    # The scale each data set's values are on: the omic's (log2; log2 counts per million for transcriptomics), or
    # z-scores when the model scaled each feature.
    scales <- stats::setNames(results$value_scales$scale, results$value_scales$view)
    expect_identical(unname(scales["transcriptomics"]), "log2_cpm")
    expect_true(all(scales[names(scales) != "transcriptomics"] == "log2"))
    zscored <- integrationResults(value, fixture$object, manager = manager, model_parameters = list(scaling = "feature_zscore"))
    expect_true(all(zscored$value_scales$scale == "feature_zscore"))
    expect_identical(zscored$feature_values, fv)
    # The heatmap: centred per feature, samples in order of the factor; scaled per feature only on request; the saved
    # values unchanged.
    h <- integrationHeatmapValues(results, "Factor1", "proteomics")
    expect_identical(dim(h), c(10L, length(samples)))
    expect_identical(rownames(h), top$feature_id[top$factor == "Factor1" & top$view == "proteomics"])
    expect_identical(colnames(h), samples[order(value$factors$Factor1)])
    expect_equal(unname(rowMeans(h, na.rm = TRUE)), rep(0, 10L), tolerance = 1e-12)
    scaled <- integrationResults(value, fixture$object, list(heatmap_scaling = "centre_scale"), manager = manager)
    expect_identical(scaled$feature_values, fv)
    hs <- integrationHeatmapValues(scaled, "Factor1", "proteomics")
    expect_equal(unname(apply(hs, 1L, stats::sd, na.rm = TRUE)), rep(1, 10L), tolerance = 1e-12)
    # The figures: ten pairs of five factors; one panel per top feature, the axis naming the values; the heatmap's legend.
    pairs <- ggplot2::ggplot_build(integrationPlotFactorPairs(results))
    expect_identical(length(unique(pairs$data[[1L]]$PANEL)), 10L)
    expect_identical(nrow(pairs$data[[1L]]), 10L * length(samples))
    scatter <- integrationPlotTopFeatureScatter(results, "Factor2", "transcriptomics")
    expect_identical(scatter$labels$y, texts[["plot_values_log2_cpm"]])
    expect_identical(length(unique(ggplot2::ggplot_build(scatter)$data[[1L]]$PANEL)), 10L)
    expect_identical(integrationPlotTopFeatureScatter(results, "Factor2", "proteomics")$labels$y, texts[["plot_values_log2"]])
    heatmap <- integrationPlotTopFeatureHeatmap(results, "Factor1", "proteomics")
    expect_s3_class(ggplot2::ggplot_build(heatmap), "ggplot_built")
    expect_identical(heatmap$labels$fill, texts[["plot_heatmap_centre"]])
    expect_identical(integrationPlotTopFeatureHeatmap(scaled, "Factor1", "proteomics")$labels$fill,
        texts[["plot_heatmap_centre_scale"]])
    one <- integrationResults(value, fixture$object, list(factors_shown = "Factor2"), manager = manager)
    expect_identical(ggplot2::ggplot_build(integrationPlotFactorPairs(one))$data[[1L]]$label, texts[["plot_pairs_one"]])
})

# Audit ENG-13: the store reads at most 10,000 rows and 64 columns at a time; the top features' values of a study with
# 70 biological samples and 10,050 features (the top features beyond the 10,000th row) are read by row and in blocks.
test_that("KIT-152: 70 samples and 10,050 features: the top features' values are read by row, in blocks of 64 samples", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir()
    input <- function(omic, values, samples) {
        manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
            source_state = list(name = "kit152-test"), sample_id_column = "Run", group_column = "group",
            technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
        matrices <- list(data.frame(feature_id = rownames(values), values, check.names = FALSE, stringsAsFactors = FALSE))
        names(matrices) <- omic
        features <- list(data.frame(feature_id = rownames(values), display_name = tolower(rownames(values)),
            stringsAsFactors = FALSE))
        names(features) <- omic
        dir <- file.path(root, "integration", omic)
        writeIntegrationInput(dir, manifest, samples, matrices, features)
        dir
    }
    design <- function(runs, biological) data.frame(Run = runs, group = sub("_.*$", "", biological),
        replicates = as.integer(sub("^.*_", "", biological)), tech_rep_group = biological, stringsAsFactors = FALSE)
    runs <- sprintf("big_%03d", seq_len(70L))
    biological <- sprintf("%s_%d", rep(c("A", "B"), each = 35L), rep(seq_len(35L), 2L))
    n <- 10050L
    values <- matrix(withr::with_seed(2152L, stats::rnorm(n * length(runs), 20, 2)), n, length(runs),
        dimnames = list(sprintf("BIG%05d", seq_len(n)), runs))
    values[n, 3L] <- NA
    big <- input("proteomics", values, design(runs, biological))
    small_runs <- sprintf("s%d", 1:4)
    small <- matrix(c(20, 21, 22, 23, 5, 6, 8, 7), 2L, 4L, byrow = TRUE, dimnames = list(c("L1", "L2"), small_runs))
    lipid <- input("lipidomics", small, design(small_runs, c("A_1", "A_2", "B_1", "B_2")))
    # Weights rising with the row: the ten largest are the last ten rows, beyond the 10,000th.
    trainer <- function(views, options) {
        samples <- colnames(views[[1L]])
        k <- options$num_factors
        labels <- paste0("Factor", seq_len(k))
        list(factors = matrix(withr::with_seed(1152L, stats::rnorm(length(samples) * k)), length(samples), k,
            dimnames = list(samples, labels)),
            weights = lapply(views, function(v) matrix(seq_len(nrow(v)) / nrow(v), nrow(v), k, dimnames = list(rownames(v), labels))),
            r2_per_factor = matrix(1, k, length(views), dimnames = list(labels, names(views))),
            r2_total = stats::setNames(rep(k, length(views)), names(views)), elbo = c(-20, -10), iterations = 2L,
            options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
            model_file = list(bytes = raw(0), sha256 = NA_character_, size = 0))
    }
    manager <- newIntegrationStateOwner(file.path(root, "project"))
    withr::defer(closeWorkflowStateOwner(manager))
    importIntegrationInputs(manager, c(big, lipid))
    # Every feature, as saved (KIT-151: the defaults select 1,000 features and z-score them).
    value <- suppressWarnings(runIntegrationModel(manager, list(feature_selection = "none", scaling = "none"),
        trainer = trainer))$value
    samples <- value$factors$biological_sample
    expect_identical(length(samples), 70L)
    original <- readWorkflowSelectedOrdinals
    reads <- list()
    testthat::local_mocked_bindings(readWorkflowSelectedOrdinals = function(store, reference, columns, ordinals, scope = NULL) {
        reads[[length(reads) + 1L]] <<- c(columns = length(columns), rows = length(ordinals), largest = max(ordinals))
        original(store, reference, columns, ordinals, scope = scope)
    })
    results <- integrationResults(value, integrationImportedObject(manager), manager = manager)
    reads <- do.call(rbind, reads)
    expect_lte(max(reads[, "columns"]), 64L)
    expect_lte(max(reads[, "rows"]), 10000L)
    expect_true(all(c(1L, 64L, 6L) %in% reads[, "columns"]))
    expect_equal(max(reads[, "largest"]), n)
    fv <- results$feature_values[results$feature_values$view == "proteomics" & results$feature_values$factor == "Factor1", ]
    expect_identical(unique(fv$feature_id), sprintf("BIG%05d", n:(n - 9L)))
    run_of <- runs[match(fv$biological_sample, biological)]
    expect_identical(fv$value, unname(values[cbind(match(fv$feature_id, rownames(values)), match(run_of, runs))]))
    expect_true(anyNA(fv$value))
    h <- integrationHeatmapValues(results, "Factor1", "proteomics")
    expect_identical(dim(h), c(10L, 70L))
    expect_s3_class(ggplot2::ggplot_build(integrationPlotTopFeatureScatter(results, "Factor1", "proteomics")), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(integrationPlotTopFeatureHeatmap(results, "Factor1", "proteomics")), "ggplot_built")
})

test_that("KIT-152: a map is refused in plain words at settings invalid for the samples, with fewer than two factors, or without its software", {
    texts <- integrationGuiTexts()
    values <- intOutFactorValues()
    factors <- paste0("Factor", 1:5)
    p <- function(...) integrationResolveResultsParameters(list(...), factors)
    refused <- function(expr) tryCatch(expr, error = function(e) e)
    # UMAP: 2 to the samples minus 1 neighbours (20 samples: at most 19).
    for (k in c(1L, 20L, 50L)) {
        e <- refused(integrationFactorMap(values, factors, "umap", p(umap = TRUE, umap_n_neighbors = k)))
        expect_s3_class(e, "integration_map_invalid")
        expect_s3_class(e, "integration_refused")
        expect_identical(conditionMessage(e), sprintf(texts[["map_refused_neighbors"]], k, 20L, 19L))
    }
    # A sample with a factor missing is not on the map: 19 samples, at most 18 neighbours.
    gap <- values
    gap$Factor3[[7L]] <- NA
    expect_identical(conditionMessage(refused(integrationFactorMap(gap, factors, "umap", p(umap_n_neighbors = 19L)))),
        sprintf(texts[["map_refused_neighbors"]], 19L, 19L, 18L))
    # t-SNE: three times the perplexity below the samples minus 1 (20 samples: below 6.33).
    e <- refused(integrationFactorMap(values, factors, "tsne", p(tsne = TRUE, tsne_perplexity = 6.4)))
    expect_s3_class(e, "integration_map_invalid")
    expect_identical(conditionMessage(e), sprintf(texts[["map_refused_perplexity"]], "6.4", 20L, "6.33"))
    expect_identical(conditionMessage(refused(integrationFactorMap(values[1:10, ], factors, "tsne", p(tsne_perplexity = 5)))),
        sprintf(texts[["map_refused_perplexity"]], "5", 10L, "3"))
    # The caution says, per map shown, what it depends on (KIT-152 correction): UMAP the seed, t-SNE not.
    none <- list(claim = list(group_corrected = character()))
    expect_identical(integrationMapCaution(none, 20L, "tsne"), unname(c(sprintf(texts[["maps_caution"]], 20L),
        texts[["maps_caution_tsne"]])))
    expect_identical(integrationMapCaution(none, 20L, c("umap", "tsne")), unname(c(sprintf(texts[["maps_caution"]], 20L),
        texts[["maps_caution_umap"]], texts[["maps_caution_tsne"]])))
    expect_identical(integrationMapCaution(list(claim = list(group_corrected = "proteomics")), 20L, "umap")[[3L]],
        sprintf(texts[["maps_caution_group_corrected"]], "proteomics"))
    # One factor: no map.
    expect_identical(conditionMessage(refused(integrationFactorMap(values, "Factor1", "umap", p()))),
        sprintf(texts[["map_refused_factors"]], 1L))
    # In the results: the refusal kept in place of the map, the other results shown; the save stops with its sentence.
    paths <- intOutProject()
    fixture <- intOutModel(paths)
    results <- integrationResults(fixture$run$value, fixture$object, list(umap = TRUE, umap_n_neighbors = 25L, tsne = TRUE, tsne_perplexity = 7),
        manager = fixture$manager)
    expect_identical(results$maps$umap$refusal, sprintf(texts[["map_refused_neighbors"]], 25L, 20L, 19L))
    expect_identical(results$maps$tsne$refusal, sprintf(texts[["map_refused_perplexity"]], "7", 20L, "6.33"))
    expect_identical(nrow(results$association), 5L)
    expect_false(is.null(results$feature_values))
    e <- refused(integrationSaveResults(fixture$manager, fixture$model$reference, results, paths))
    expect_s3_class(e, "integration_map_invalid")
    expect_identical(integrationRefusalText(e)$sentence, results$maps$umap$refusal)
    expect_false(dir.exists(paths$mofa_tables_dir) && length(list.files(paths$mofa_tables_dir)) > 0L)
    # Not chosen: no map, nothing refused.
    expect_length(integrationResults(fixture$run$value, fixture$object, manager = fixture$manager)$maps, 0L)
    # Valid settings without the software: refused in words, naming it.
    testthat::local_mocked_bindings(optionalPackageAvailable = function(package) FALSE)
    e <- refused(integrationFactorMap(values, factors, "umap", p(umap_n_neighbors = 19L)))
    expect_s3_class(e, "integration_map_invalid")
    expect_identical(conditionMessage(e), sprintf(texts[["map_unavailable"]], "UMAP", "uwot"))
    expect_identical(conditionMessage(refused(integrationFactorMap(values, factors, "tsne", p(tsne_perplexity = 6.3)))),
        sprintf(texts[["map_unavailable"]], "t-SNE", "Rtsne"))
})

test_that("KIT-152: saved results carry the outputs; the head keeps plain lists and table references; record, methods and replay describe them", {
    paths <- intOutProject()
    fixture <- intOutModel(paths, close = FALSE)
    manager <- fixture$manager
    texts <- integrationRecordTexts()
    results <- integrationResults(fixture$run$value, fixture$object, manager = manager)
    display <- integrationSaveResults(manager, fixture$model$reference, results, paths)
    expect_true("top_feature_values.tsv" %in% display$tables)
    expect_false(any(grepl("^map_", display$tables)))
    expect_true(all(c("factor_pairs.png", "top_features_Factor1_proteomics.png", "top_features_heatmap_Factor5_transcriptomics.png")
        %in% display$figures))
    expect_identical(length(display$figures), 5L + 15L * length(results$views))
    written <- utils::read.delim(file.path(display$tables_dir, "top_feature_values.tsv"), stringsAsFactors = FALSE)
    expect_identical(written$feature_id, results$feature_values$feature_id)
    expect_equal(written$value, results$feature_values$value, tolerance = 1e-12)
    # On the head: what the outputs show as plain lists; the values as a table of the store, its reference only.
    export <- manager$.context$snapshot$exports$integration_mofa_model
    expect_identical(names(export$results$output_tables), "feature_values")
    expect_true(isWorkflowTableReference(export$results$output_tables$feature_values))
    expect_identical(export$results$outputs$value_scales$scale, results$value_scales$scale)
    expect_identical(export$results$outputs$heatmap_scaling, "centre")
    expect_length(export$results$outputs$maps, 0L)
    expect_no_error(collectWorkflowNativeReferences(export$results))
    stored <- integrationRecordedOutputs(manager)
    expect_identical(as.list(stored$feature_values), as.list(results$feature_values))
    # The record: no reference in it; the outputs' sheet, the study parameters and the methods paragraph describe them.
    record <- buildWorkflowAnalysisRecord(manager)
    block <- record$analyses$integration_mofa_model$integration
    expect_null(block$results$output_tables)
    json <- file.path(paths$base_dir, "record.json")
    jsonlite::write_json(record, json, auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
    read <- jsonlite::read_json(json, simplifyVector = FALSE)
    tables <- integrationRecordTables(read)
    expect_true("integration_outputs" %in% names(tables))
    outputs <- tables$integration_outputs
    expect_identical(outputs[[texts[["col_output"]]]], unname(texts[c("out_pairs", "out_top", "out_scatter", "out_heatmap")]))
    expect_match(outputs[[texts[["col_shows"]]]][[3L]], "transcriptomics: log2 counts per million", fixed = TRUE)
    expect_identical(integrationRecordTables(record)$integration_outputs, outputs)
    methods <- integrationMethodsText(read)
    expect_identical(integrationMethodsText(record), methods)
    for (fragment in c(sprintf(texts[["methods_top_features"]], 10L), texts[["methods_heatmap_centre"]],
        "before MOFA2 centres and scales each data set; ", "proteomics: log2", "transcriptomics: log2 counts per million")) {
        expect_true(grepl(fragment, methods, fixed = TRUE), info = fragment)
    }
    # A map not chosen is absent from the methods paragraph and the record.
    expect_false(grepl("UMAP|t-SNE", methods))
    lines <- integrationStudyParameterLines(read)
    expect_true(texts[["section_outputs"]] %in% lines)
    # Reopened (the head holds no live payload) and repeated from the record: the same top features' values.
    root <- dirname(manager$.context$store$directory)
    closeWorkflowStateOwner(manager)
    replay <- suppressWarnings(replayIntegrationProject(root, file.path(paths$base_dir, "replay"), trainer = intOutTrainer()))
    expect_true(replay$comparison$identical)
    comparison <- replay$outputs_comparison
    expect_identical(comparison$tables$table, "top_feature_values")
    expect_true(comparison$identical)
    expect_identical(comparison$tables$largest_relative, 0)
    expect_false(comparison$tables$same_computer_only)
    # A different value is told apart.
    changed <- replay$results
    at <- which(!is.na(changed$feature_values$value))[[1L]]
    changed$feature_values$value[[at]] <- changed$feature_values$value[[at]] + 1
    expect_false(integrationCompareOutputs(list(feature_values = results$feature_values), changed)$equal)
})

test_that("KIT-152: the Results tab's settings change the results; the outputs and the maps' status are shown in words", {
    paths <- intOutProject()
    fixture <- intOutModel(paths)
    texts <- integrationGuiTexts()
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate({
        bus$state_manager <- fixture$manager
        bus$integration_model <- fixture$model
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        results <- session$returned()
        expect_identical(results$parameters$top_features, 10L)
        expect_match(output$maps_status$html, texts[["maps_off"]], fixed = TRUE)
        session$setInputs(weights_factor = "Factor2", weights_view = "transcriptomics")
        intro <- gsub("&#39;", "'", output$scatter_intro$html, fixed = TRUE)
        expect_match(intro, sprintf(texts[["scatter_values"]], texts[["values_log2_cpm"]]), fixed = TRUE)
        expect_match(gsub("&#39;", "'", output$heatmap_intro$html, fixed = TRUE), texts[["heatmap_intro_centre"]], fixed = TRUE)
        # A setting changed on the tab: five features listed.
        session$setInputs(param_top_features = 5, param_heatmap_scaling = "centre_scale")
        results <- session$returned()
        expect_identical(results$parameters$top_features, 5L)
        expect_identical(nrow(results$feature_values), 5L * length(results$shown) * length(results$views) * 20L)
        expect_match(gsub("&#39;", "'", output$heatmap_intro$html, fixed = TRUE), texts[["heatmap_intro_centre_scale"]],
            fixed = TRUE)
        # A map chosen with too many neighbours: refused in words in its place; the rest is shown; the save says why.
        session$setInputs(param_umap = TRUE, param_umap_n_neighbors = 30)
        expect_match(output$maps_status$html, sprintf(texts[["map_refused_neighbors"]], 30L, 20L, 19L), fixed = TRUE)
        expect_match(output$variance_table, "Factor1", fixed = TRUE)
        session$setInputs(save = 1)
        expect_match(output$save_status$html, sprintf(texts[["map_refused_neighbors"]], 30L, 20L, 19L), fixed = TRUE)
        expect_null(bus$integration_results)
        # A number box left empty: refused in words.
        session$setInputs(param_umap = FALSE, param_tsne_perplexity = NA)
        expect_match(output$results_status$html, sprintf(texts[["setting_not_number"]], texts[["param_tsne_perplexity"]]),
            fixed = TRUE)
    })
})

test_that("KIT-152: factors picked for an earlier model that a newer model lacks fall back to all factors, said in one line", {
    paths <- intOutProject()
    fixture <- intOutModel(paths)
    texts <- integrationGuiTexts()
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate({
        bus$state_manager <- fixture$manager
        bus$integration_model <- fixture$model
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$flushReact()
        session$setInputs(param_factors_shown = "Factor5")
        expect_identical(session$returned()$shown, "Factor5")
        expect_false(grepl(sprintf(texts[["factors_shown_reset"]], 5L), output$results_status$html, fixed = TRUE))
        # A newer model with four factors: Factor5 is not in it.
        other <- runIntegrationModel(fixture$manager, list(num_factors = 4L), trainer = function(views, options) {
            fit <- intOutTrainer()(views, options)
            fit$factors <- fit$factors[, 1:4, drop = FALSE]
            fit$weights <- lapply(fit$weights, function(w) w[, 1:4, drop = FALSE])
            fit$r2_per_factor <- fit$r2_per_factor[1:4, , drop = FALSE]
            fit
        })
        bus$integration_model <- list(reference = other$reference, parameters = list(), status = "computed", trained_at = Sys.time())
        session$flushReact()
        results <- session$returned()
        expect_false(inherits(results, "integration_results_refused"))
        expect_identical(results$shown, paste0("Factor", 1:4))
        expect_identical(results$parameters$factors_shown, paste0("Factor", 1:4))
        status <- gsub("&quot;", "\"", gsub("&#39;", "'", output$results_status$html, fixed = TRUE), fixed = TRUE)
        expect_match(status, sprintf(texts[["factors_shown_reset"]], 4L), fixed = TRUE)
        # Factors this model has stay as picked, with no note.
        session$setInputs(param_factors_shown = c("Factor2", "Factor4"))
        expect_identical(session$returned()$shown, c("Factor2", "Factor4"))
        expect_null(session$returned()$settings_note)
    })
})
