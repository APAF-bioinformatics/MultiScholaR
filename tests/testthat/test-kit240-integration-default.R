# KIT-240: the integration tabs on by default in the release (INT-003, VAL-209), the Model tab's words when MOFA2, reticulate
# or mofapy2 is missing, the reworded marks of the integration defaults (VAL-221 (4), the values unchanged) and each
# training seed's final ELBO kept with the reported model's rank shown (VAL-221 (1)). Rows of the ticket, section 5, with
# the value each shows on the base commit (3eeaf100), where it fails first:
#   T1 the option unset: the tile live, the home script edited, the six tabs (base: off, disabled tile, the panel);
#   T3 the Model tab names the missing piece and how to get it (base: the generic sentence for all three);
#   T4 the real probe names what is missing (base: no `missing`); skipped where MOFA2 is installed;
#   T6 the marks name their decision; exactly the eight waiting enrichment settings say so; every default and choice as
#      frozen from base (base: (a), (b), (d) fail; (c) passes, a guard);
#   T7 each seed's final ELBO kept, the reported model's rank shown on the Results tab and in the record, the reported
#      model unchanged (base: no `elbo_last`, no `elbo_rank`, no sentence);
#   T8 seeds without an ELBO (the stand-ins of test-integration-model.R): no rank, no sentence, no error (passes on base);
#   T9 (MOFA2 lane) five finite final ELBOs and a rank in 1 to 5 on the tier-1 corpus (skipped without MOFA2).
# The small inputs are those of test-integration-model.R (KIT-062), written by the contract's writer.

k240Input <- function(root, omic, views, samples, standards = list()) {
    matrices <- lapply(views, function(m) data.frame(feature_id = rownames(m), m, check.names = FALSE,
        stringsAsFactors = FALSE))
    features <- lapply(stats::setNames(names(views), names(views)), function(view) {
        table <- data.frame(feature_id = rownames(views[[view]]), display_name = tolower(rownames(views[[view]])),
            stringsAsFactors = FALSE)
        if (!is.null(standards[[view]])) table$is_internal_standard <- standards[[view]]
        table
    })
    manifest <- list(omic_type = omic, workflow = "TEST", matrix_role = "normalised", scale = "log2", correction = NULL,
        source_state = list(name = "kit240-test"), sample_id_column = "Run", group_column = "group",
        technical_replicate_column = "tech_rep_group", biological_sample_column = "tech_rep_group")
    dir <- file.path(root, "integration", omic)
    writeIntegrationInput(dir, manifest, samples, matrices, features)
    dir
}

k240Samples <- function(runs, biological) {
    data.frame(Run = runs, group = sub("_.*$", "", biological), replicates = as.integer(sub("^.*_", "", biological)),
        tech_rep_group = biological, stringsAsFactors = FALSE)
}

# Proteomics (A_1..A_3, B_1, B_2) and lipidomics (two runs per sample, A_3 absent; ITSD1 an internal standard), imported
# into a store of their own.
k240Owner <- function(.local_envir = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = .local_envir)
    biological <- c("A_1", "A_2", "A_3", "B_1", "B_2")
    runs <- paste0("p", seq_along(biological))
    values <- rbind(P1 = c(10, 11, 12, 13, 14), P2 = c(5, 7, 6, 9, 8), P3 = c(1, 1.5, 2, 2.5, 3.5), P4 = c(NA, NA, 4, NA, NA))
    colnames(values) <- runs
    proteomics <- k240Input(root, "proteomics", list(proteomics = values), k240Samples(runs, biological))
    biological <- c("A_1", "A_1", "A_2", "A_2", "B_1", "B_1", "B_2", "B_2")
    runs <- paste0("l", seq_along(biological))
    values <- rbind(L1 = c(20, 22, 21, 23, 25, 27, 26, 28), L2 = c(3, NA, 4, 5, 6, 7, NA, NA), ITSD1 = rep(15, 8L))
    colnames(values) <- runs
    lipidomics <- k240Input(root, "lipidomics", list(lipidomics_LCMS_Pos = values), k240Samples(runs, biological),
        standards = list(lipidomics_LCMS_Pos = c(FALSE, FALSE, TRUE)))
    manager <- newIntegrationStateOwner(file.path(root, "project"))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    importIntegrationInputs(manager, c(proteomics, lipidomics))
    list(root = root, manager = manager)
}

# A stand-in trainer in `integrationTrainMofa`'s form: three factors (those of test-integration-model.R's seed test), the
# declared seed's ELBO trace ending at -100 (`declared_elbo`); each further seed (`factors_only`) gives its final ELBO from
# `others_elbo` in seed order, or none when `others_elbo` is NULL (the form of the stand-ins before KIT-240).
k240Trainer <- function(declared_elbo = c(-130, -100, NaN), others_elbo = c(-99, -101, -100.5, -102)) {
    function(views, options) {
        samples <- colnames(views[[1L]])
        labels <- paste0("Factor", seq_len(3L))
        f <- matrix(c(1, 4, 2, 5, 3, 0, 1, 0, 1, 1, 9, 7, 8, 6, 5), 5L, 3L, dimnames = list(samples, labels))
        if (isTRUE(options$factors_only)) {
            other <- list(factors = f, iterations = 7L)
            if (!is.null(others_elbo)) other$elbo_last <- others_elbo[[options$seed - 20261119L]]
            return(other)
        }
        list(factors = f, weights = lapply(views, function(v) matrix(0.5, nrow(v), 3L, dimnames = list(rownames(v), labels))),
            r2_per_factor = matrix(1, 3L, length(views), dimnames = list(labels, names(views))),
            r2_total = stats::setNames(rep(3, length(views)), names(views)), elbo = declared_elbo, iterations = 2L,
            options = options, versions = list(MOFA2 = "stand-in"), messages = character(),
            model_file = list(bytes = raw(0), sha256 = NA_character_, size = 0))
    }
}

k240Html <- function(tag) as.character(htmltools::renderTags(tag)$html)

# The HTML of a rendered output, with the escapes htmltools may apply to quotes undone (the texts are matched as written).
k240Plain <- function(html) gsub("&quot;", "\"", gsub("&#39;", "'", html, fixed = TRUE), fixed = TRUE)

# A wording-table entry, NA where the table has none (the base commit lacks the new entries; the test then fails, not
# errors).
k240Text <- function(texts, key) if (key %in% names(texts)) texts[[key]] else NA_character_

# The model trained on the small inputs at three factors and the default five seeds, with the stand-in trainer.
k240Run <- function(trainer, .local_envir = parent.frame()) {
    owner <- k240Owner(.local_envir = .local_envir)
    run <- suppressWarnings(runIntegrationModel(owner$manager, list(num_factors = 3L), trainer = trainer))
    export <- owner$manager$.context$snapshot$exports$integration_mofa_model
    c(owner, list(run = run, model = list(reference = run$reference, parameters = export$parameters, status = run$status,
        trained_at = Sys.time())))
}

# The Results tab's stability sentences for a trained model (the tab's own server, as the app runs it).
k240StabilityIntro <- function(fixture) {
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate({
        bus$state_manager <- fixture$manager
        bus$integration_model <- fixture$model
    })
    html <- NULL
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = list()), {
        session$flushReact()
        html <<- output$stability_intro$html
    })
    k240Plain(html)
}

test_that("KIT-240 T1: with the option unset the Integration tile is live and the six tabs are shown", {
    withr::local_options(MultiScholaR.integration_gui = NULL)
    expect_true(integrationGuiEnabled())
    page <- k240Html(app_ui(NULL))
    expect_match(page, paste0('<div id="integration_box" class="omic-selection-box" ',
        'onclick="toggleOmicSelection(&#39;integration&#39;)" data-testid="tile-integration">'), fixed = TRUE)
    expect_match(page, "if (true) { // KIT-138: the Integration tile stays enabled", fixed = TRUE)
    expect_false(grepl('<div id="integration_box" class="omic-selection-box disabled"', page, fixed = TRUE))
    ui <- k240Html(mod_integration_ui("integration_workflow"))
    expect_match(ui, 'data-testid="integration-workflow-tabs"', fixed = TRUE)
    for (tab in c("import", "alignment", "model", "results", "enrichment", "record")) {
        expect_match(ui, sprintf('data-testid="integration-tab-%s"', tab), fixed = TRUE)
    }
})

test_that("KIT-240 T3: the Model tab names the missing MOFA2, reticulate or mofapy2 and how to get it", {
    texts <- integrationGuiTexts()
    probes <- list(
        MOFA2 = list(available = FALSE, reason = "the R package MOFA2 is not installed", missing = "MOFA2"),
        reticulate = list(available = FALSE, reason = "the R package reticulate is not installed", missing = "reticulate"),
        mofapy2 = list(available = FALSE, reason = paste("no Python with the module mofapy2 is available to reticulate",
            "(RETICULATE_PYTHON is not set)"), missing = "mofapy2"))
    readiness_of <- function(probe) {
        bus <- shiny::reactiveValues(integration_alignment = list(parameters = list()), state_manager = NULL,
            integration_model = NULL, tab_status = list())
        html <- NULL
        shiny::testServer(mod_integration_model_server, args = list(workflow_data = bus, paths = list(),
            probe = function() probe, withProgressFn = function(message, value, expr) force(expr)), {
            session$flushReact()
            html <<- output$readiness$html
        })
        k240Plain(html)
    }
    for (name in names(probes)) {
        readiness <- readiness_of(probes[[name]])
        sentence <- k240Text(texts, paste0("readiness_missing_", tolower(name)))
        expect_false(grepl(texts[["refusal_mofa_unavailable"]], readiness, fixed = TRUE), info = name)
        expect_true(!is.na(sentence) && grepl(sentence, readiness, fixed = TRUE), info = name)
        expect_true(!is.na(sentence) && grepl(name, sentence, fixed = TRUE), info = name)
        # The probe's reason stays under "Technical detail".
        expect_match(readiness, sprintf("<summary>%s</summary>", texts[["technical_detail"]]), fixed = TRUE, info = name)
        expect_match(readiness, probes[[name]]$reason, fixed = TRUE, info = name)
    }
    # A probe that does not say what is missing (a stand-in): the general sentence, as before.
    readiness <- readiness_of(list(available = FALSE, reason = "no Python here"))
    expect_match(readiness, texts[["refusal_mofa_unavailable"]], fixed = TRUE)
    expect_match(readiness, "no Python here", fixed = TRUE)
})

test_that("KIT-240 T4: the software check says which piece is missing", {
    # The piece this machine lacks, checked in the probe's order: MOFA2, then reticulate, then a Python with mofapy2.
    expected <- if (!requireNamespace("MOFA2", quietly = TRUE)) "MOFA2" else
        if (!requireNamespace("reticulate", quietly = TRUE)) "reticulate" else "mofapy2"
    capability <- integrationMofaCapability()
    testthat::skip_if(isTRUE(capability$available), "MOFA2, reticulate and a Python with mofapy2 are all present here")
    expect_false(capability$available)
    expect_identical(capability$missing, expected)
    if (identical(expected, "MOFA2")) expect_identical(capability$reason, "the R package MOFA2 is not installed")
})

test_that("KIT-240 T6: the marks name the decision behind each default; the eight waiting enrichment settings say so", {
    # (a) Every model and Results default names its decision and ledger row; none is marked provisional.
    tables <- list(model = integrationMofaParameterTable(), results = integrationResultsParameterTable())
    for (name in names(tables)) {
        table <- tables[[name]]
        expect_false(any(grepl("workbook|provisional|may override|awaiting", table$settled_by)),
            info = paste(name, paste(table$parameter[grepl("workbook|provisional|may override|awaiting", table$settled_by)],
                collapse = ", ")))
        expect_true(all(grepl("INT-|VAL-|KIT-060", table$settled_by)), info = name)
    }
    for (spec in list(integrationMofaParameterSpec(), integrationResultsParameterSpec())) {
        expect_false(any(vapply(spec, integrationParameterProvisional, logical(1))))
    }
    # (b) Exactly the eight enrichment settings still awaiting the owner say so; the STRING pair names its ledger row.
    enrichment <- integrationSettings()$enrichment
    waiting <- c("correlation_method", "min_samples", "string_species", "kegg_organism", "reactome_organism",
        "min_set_size", "max_set_size", "seed")
    expect_setequal(enrichment$setting[grepl("awaiting the owner", enrichment$settled_by, fixed = TRUE)], waiting)
    expect_true(all(grepl("^KIT-141 \\(decided 2026-11-20",
        enrichment$settled_by[enrichment$setting %in% c("string_fdr", "string_rank_direction")])))
    # (c) No value changes: every default and choice as at the base commit (frozen here from 3eeaf100).
    d <- function(default, choices = NULL) list(default = default, choices = choices)
    values <- function(spec) lapply(spec, function(entry) entry[c("default", "choices")])
    expect_identical(values(integrationMofaParameterSpec()), list(
        input_matrix = d("uncorrected", c("uncorrected", "ruv_corrected")),
        technical_replicates = d("average", "average"),
        internal_standards = d("drop", c("drop", "keep")),
        missing_samples = d("keep_missing", c("keep_missing", "drop_incomplete")),
        feature_selection = d("variance_above_trend", c("variance_above_trend", "top_variance", "none")),
        features_per_view = d(1000L), min_observed_fraction = d(0.5),
        scaling = d("feature_zscore", c("none", "feature_zscore")),
        scale_views = d(TRUE), center_groups = d(TRUE), num_factors = d(5L),
        convergence_mode = d("slow", c("fast", "medium", "slow")),
        maxiter = d(10000L), seed = d(20261119L), training_seeds = d(5L),
        precision = d("double", c("double", "single"))))
    expect_identical(values(integrationResultsParameterSpec()), list(
        factors_shown = d("all"), top_features = d(10L), weight_scale = d("scaled", c("scaled", "raw")),
        design_covariate = d("group", "group"), association_test = d("anova", c("anova", "kruskal_wallis")),
        p_adjust_method = d("BH", c("BH", "holm", "bonferroni", "none")), figure_format = d("png", c("png", "pdf")),
        figure_width = d(8), figure_height = d(6), figure_dpi = d(300L), feature_values = d("model_input", "model_input"),
        heatmap_scaling = d("centre", c("centre", "centre_scale")), umap = d(FALSE), umap_n_neighbors = d(15L),
        tsne = d(FALSE), tsne_perplexity = d(5), map_seed = d(20261119L)))
    expect_identical(values(integrationEnrichmentParameterSpec()), list(
        correlation_method = d("pearson", c("pearson", "spearman")), min_samples = d(5L), string_species = d(9606L),
        string_fdr = d(0.05), string_rank_direction = d(-1L), kegg_organism = d("hsa"),
        reactome_organism = d("Homo sapiens"), min_set_size = d(5L), max_set_size = d(500L), seed = d(20261121L),
        test = d(NULL, c("camera", "camera_floor", "permutation", "gsea")), permutations = d(1000L, c(1000L, 10000L)),
        adjust_ruv_factors = d(FALSE, c(TRUE, FALSE))))
    # (d) Nothing on the Model and Results tabs says "provisional".
    expect_false(grepl("provisional", k240Html(mod_integration_model_ui("m")), ignore.case = TRUE))
    expect_false(grepl("provisional", k240Html(mod_integration_results_ui("r")), ignore.case = TRUE))
})

test_that("KIT-240 T7: each seed's final ELBO is kept and the reported model's rank shown, the model unchanged", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    fixture <- k240Run(k240Trainer())
    seeds <- fixture$run$value$seeds
    expect_identical(seeds$seeds, 20261119L + 0:4)
    expect_equal(seeds$models$elbo_last, c(-100, -99, -101, -100.5, -102))
    expect_equal(seeds$elbo_rank, list(rank = 2L, of = 5L, below_best = 1, range = 3, relative_spread = 3 / 99))
    # r1.4 (AUDIT-KIT-240-impl I01): a distance below 0.005 never prints as "0.00" (rank 2 or lower); zero (rank 1) does.
    expect_identical(integrationElboText(0.004), "less than 0.01")
    expect_identical(integrationElboText(0), "0.00")
    expect_identical(integrationElboText(1), "1.00")
    # The reported model is the one trained without the further seeds' ELBOs.
    plain <- k240Run(k240Trainer(declared_elbo = c(-130, -100, NaN), others_elbo = NULL))
    for (field in c("factors", "variance_explained", "variance_total", "elbo")) {
        expect_identical(fixture$run$value[[field]], plain$run$value[[field]], info = field)
    }
    expect_identical(integrationModelWeights(fixture$manager, fixture$run$value),
        integrationModelWeights(plain$manager, plain$run$value))
    # The Results list carries the rank; the Results tab shows it and the stability table's limitation.
    object <- integrationImportedObject(fixture$manager)
    results <- integrationResults(fixture$run$value, object, list(), manager = fixture$manager)
    expect_equal(results$elbo_rank, seeds$elbo_rank)
    texts <- integrationGuiTexts()
    intro <- k240StabilityIntro(fixture)
    rank_sentence <- k240Text(texts, "stability_elbo_rank")
    limit <- k240Text(texts, "stability_limit")
    expect_true(!is.na(rank_sentence) && grepl(sprintf(rank_sentence, 2L, 5L, "1.00", "3.00"), intro, fixed = TRUE))
    expect_true(!is.na(limit) && grepl(limit, intro, fixed = TRUE))
    # The record keeps the rank, and its methods paragraph states it, as built and from JSON.
    record <- buildWorkflowAnalysisRecord(fixture$manager)
    expect_equal(record$analyses$integration_mofa_model$integration$model$elbo_rank, seeds$elbo_rank)
    json <- file.path(fixture$root, "record.json")
    writeWorkflowAnalysisRecord(fixture$manager, json)
    methods <- integrationMethodsText(record)
    expect_identical(integrationMethodsText(jsonlite::read_json(json, simplifyVector = FALSE)), methods)
    methods_sentence <- k240Text(integrationRecordTexts(), "methods_seeds_elbo")
    expect_true(!is.na(methods_sentence) && grepl(sprintf(methods_sentence, 2L, 5L, "1.00", "3.00"), methods, fixed = TRUE))
})

test_that("KIT-240 T8: seeds trained without an ELBO give no rank and no sentence, and nothing fails", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    # The form of test-integration-model.R's stand-ins: the declared seed's trace ends at -10, the others give none.
    fixture <- k240Run(k240Trainer(declared_elbo = c(-20, -10, NaN), others_elbo = NULL))
    seeds <- fixture$run$value$seeds
    expect_null(seeds$elbo_rank)
    # Before KIT-240 the column does not exist; with it, the declared seed's value and NA for the others.
    elbo_last <- seeds$models$elbo_last
    expect_true(is.null(elbo_last) || identical(elbo_last, c(-10, NA, NA, NA, NA)))
    object <- integrationImportedObject(fixture$manager)
    results <- integrationResults(fixture$run$value, object, list(), manager = fixture$manager)
    expect_null(results$elbo_rank)
    intro <- k240StabilityIntro(fixture)
    expect_match(intro, sprintf(integrationGuiTexts()[["stability_intro"]], 5L, paste(20261119L + 0:4, collapse = ", "),
        20261119L), fixed = TRUE)
    expect_false(grepl("ELBO", intro, fixed = TRUE))
    record <- buildWorkflowAnalysisRecord(fixture$manager)
    expect_null(record$analyses$integration_mofa_model$integration$model$elbo_rank)
    methods <- integrationMethodsText(record)
    expect_false(grepl("ELBO", methods, fixed = TRUE))
    # r1.4 (I09): the stability figures' limitation is stated whenever several seeds were trained, with or without a rank.
    limit <- k240Text(integrationRecordTexts(), "methods_seeds_limit")
    expect_true(!is.na(limit) && grepl(limit, methods, fixed = TRUE))
})

test_that("KIT-240 T9 (MOFA2 lane): every seed's final ELBO is finite and the reported model's rank is in 1 to 5", {
    corpus <- normalizePath(testthat::test_path("..", "testdata", "integration", "corpus"), mustWork = FALSE)
    testthat::skip_if_not(dir.exists(corpus), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    capability <- integrationMofaCapability()
    testthat::skip_if_not(isTRUE(capability$available), capability$reason)
    root <- withr::local_tempdir()
    manager <- newIntegrationStateOwner(file.path(root, "project"))
    withr::defer(closeWorkflowStateOwner(manager))
    importIntegrationInputs(manager, file.path(corpus, c("proteomics", "transcriptomics", "lipidomics", "metabolomics")))
    run <- runIntegrationModel(manager)
    seeds <- run$value$seeds
    rank <- seeds$elbo_rank
    message(sprintf("KIT-240 T9 final ELBO by seed %s: %s; rank %s of %s, %s below the highest, range %s (relative %s)",
        paste(seeds$models$seed, collapse = ", "), paste(sprintf("%.6f", seeds$models$elbo_last), collapse = ", "),
        rank$rank, rank$of, format(rank$below_best, digits = 10), format(rank$range, digits = 10),
        format(rank$relative_spread, digits = 6)))
    expect_identical(seeds$seeds, 20261119L + 0:4)
    expect_length(seeds$models$elbo_last, 5L)
    expect_true(all(is.finite(seeds$models$elbo_last)))
    expect_identical(seeds$models$elbo_last[[1L]], run$value$elbo$last)
    expect_true(rank$rank %in% 1:5)
    expect_identical(rank$of, 5L)
    expect_true(rank$below_best >= 0 && rank$range >= rank$below_best)
})

# KIT-240 r1.3 (coordinator's ruling on AUDIT-KIT-240-plan Re-check 1, R1-04): the Enrichment tab marks the eight
# defaults that are not yet confirmed (their `settled_by` says they await confirmation), beside each setting's label; the
# confirmed ones carry no mark. Base: no mark on the tab (the tab never read `settled_by`).
test_that("KIT-240 T10: the Enrichment tab marks its eight defaults not yet confirmed, and only those", {
    texts <- integrationGuiTexts()
    ui <- k240Plain(k240Html(mod_integration_enrichment_ui("e")))
    marked <- c(correlation_method = "enrich_correlation_method", min_samples = "enrich_min_samples",
        string_species = "enrich_string_species", kegg_organism = "enrich_kegg_organism",
        reactome_organism = "enrich_reactome_organism", min_set_size = "enrich_min_set", max_set_size = "enrich_max_set",
        seed = "enrich_seed")
    mark <- sprintf(" (%s)", texts[["mark_provisional"]])
    for (name in names(marked)) {
        expect_match(ui, paste0(texts[[marked[[name]]]], mark, "</label>"), fixed = TRUE, info = name)
    }
    expect_false(grepl(paste0(texts[["enrich_string_fdr"]], mark), ui, fixed = TRUE))
    # Eight marks, no more; the mark's words say plainly that the default is not yet confirmed.
    expect_identical(lengths(regmatches(ui, gregexpr(mark, ui, fixed = TRUE))), 8L)
    expect_match(texts[["mark_provisional"]], "not yet confirmed", fixed = TRUE)
    expect_false(grepl("owner", texts[["mark_provisional"]], fixed = TRUE))
})

# r1.4 (AUDIT-KIT-240-impl-20261125 I03, I05, I10): the alignment note is scoped to the three rules it cites; the report's
# status legend uses a user's words; the Methods paragraph defines the ELBO and names its units. Base (fac10b51): the note
# covers "These rules", the legend says "the owner", the Methods sentence starts "By its final ELBO".
test_that("KIT-240 T12: the alignment note, the report legend and the ELBO sentence say exactly what is true", {
    gui <- integrationGuiTexts()
    record <- integrationRecordTexts()
    # I03: the note names the three rules it covers.
    for (rule in c("averaging of runs", "internal standards", "missing samples")) {
        expect_match(gui[["align_provisional"]], rule, fixed = TRUE, info = rule)
    }
    # I05: the legend and the status word, as on the Enrichment tab; no "owner".
    expect_false(grepl("owner", record[["legend_provisional_default"]], fixed = TRUE))
    expect_identical(record[["status_provisional_default"]], gui[["mark_provisional"]])
    # I10: the ELBO defined and its units named in the Methods sentence.
    expect_match(record[["methods_seeds_elbo"]], "evidence lower bound (ELBO", fixed = TRUE)
    expect_match(record[["methods_seeds_elbo"]], "units", fixed = TRUE)
    # I09: the rank sentence no longer carries the limitation; it is a sentence of its own.
    expect_false(grepl("stability figures", record[["methods_seeds_elbo"]], fixed = TRUE))
})
