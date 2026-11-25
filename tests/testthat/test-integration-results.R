# KIT-063: the Results tab and its helpers, without MOFA2. The model is the stored double-precision reference
# (tests/testdata/integration/mofa_reference, MOFA2 1.22.0 / mofapy2 0.7.3 on the tier-1 corpus at the provisional
# defaults), returned by a stand-in trainer through the model step on the corpus imported into a store, so the
# product is the one MOFA2 gave. Declared before the run:
#   - the variance, factor values, weights and factor correlation shown are the product's, with the features' display
#     names; the factor-group statistic equals an independent computation (summary(lm)$r.squared, oneway.test);
#   - the corpus's planted truth (tests/testdata/integration/corpus/truth*.tsv, truth.json) is visible in the results:
#     the factor matched to Z1 (the group) has at least 90% of its variance between the groups with an adjusted
#     p-value below 0.001; the factor matched to Z3 explains at least 1% of the variance in proteomics and
#     transcriptomics only; in every data set where Z1, Z2 or Z3 is planted, the top features of its factor hold more
#     of its planted features than the data set as a whole, with a one-sided hypergeometric p-value below 0.05;
#   - under a correction with the groups as replicate sets (KIT-060 D10) no p-value is given and the tab says why;
#   - MOFA2's warning that a factor follows the samples' overall level is shown in plain words;
#   - the saved tables and figures reach results/integration/mofa through the action-output publication, every file
#     checked; a figure the device did not write stops the save by name.
# The same checks on a model MOFA2 trains in the job are test-integration-results-mofa.R (MOFA2 lane).

intResData <- function() normalizePath(testthat::test_path("..", "testdata", "integration"), mustWork = FALSE)
intResOmics <- c("lipidomics", "metabolomics", "proteomics", "transcriptomics")
# MOFA2 1.22.0's text (MOFA2 lane run 44616769378, on the tier-1 corpus).
intResOverallLevel <- paste("WARNING: Factor(s) 2 are strongly correlated with the average expression of features for",
    "at least one of your omics. Such factors appear when there are differences in the total 'levels' between your",
    "samples, *sometimes* because of poor normalisation in the preprocessing steps.")

intResStored <- function(name) utils::read.delim(file.path(intResData(), "mofa_reference", paste0(name, ".tsv")),
    check.names = FALSE, stringsAsFactors = FALSE)

# The stored reference as `integrationTrainMofa` returns a model (the total per data set as the sum over factors).
intResTrainer <- function(messages = character()) {
    function(views, options) {
        factors <- intResStored("factors")
        f <- as.matrix(factors[-1L])
        rownames(f) <- factors$biological_sample
        f <- f[colnames(views[[1L]]), , drop = FALSE]
        w <- intResStored("weights")
        weights <- lapply(stats::setNames(names(views), names(views)), function(view) {
            rows <- w[w$view == view, , drop = FALSE]
            m <- as.matrix(rows[-(1:2)])
            rownames(m) <- paste(view, rows$feature_id, sep = ":")
            m[rownames(views[[view]]), , drop = FALSE]
        })
        v <- intResStored("variance_explained")
        r2 <- matrix(v$r2_percent, nrow = length(unique(v$factor)), dimnames = list(unique(v$factor), unique(v$view)))
        r2 <- r2[, names(views), drop = FALSE]
        list(factors = f, weights = weights, r2_per_factor = r2, r2_total = colSums(r2), elbo = c(-2, -1), iterations = 186L,
            options = options, versions = list(MOFA2 = "stored reference"), messages = messages,
            model_file = list(bytes = as.raw(1:3), sha256 = strrep("d", 64L), size = 3))
    }
}

intResProject <- function(.local_envir = parent.frame()) {
    testthat::skip_if_not(dir.exists(file.path(intResData(), "corpus")), "the tier-1 corpus is not present")
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    root <- withr::local_tempdir(.local_envir = .local_envir)
    results_dir <- file.path(root, "results", "integration")
    paths <- list(base_dir = root, integration_dir = file.path(root, "integration"),
        source_dir = file.path(root, "scripts", "integration"), results_dir = results_dir,
        multiomic_inputs_dir = file.path(results_dir, "multiomic_inputs"), mofa_tables_dir = file.path(results_dir, "mofa", "tables"),
        mofa_plots_dir = file.path(results_dir, "mofa", "plots"))
    dir.create(paths$source_dir, recursive = TRUE)
    for (omic in intResOmics) {
        dir <- file.path(paths$integration_dir, omic)
        dir.create(dir, recursive = TRUE)
        stopifnot(all(file.copy(list.files(file.path(intResData(), "corpus", omic), full.names = TRUE), dir)))
    }
    paths
}

# The corpus imported into a store and the stored model accepted on it, as the Model tab leaves the bus.
intResModel <- function(paths, messages = intResOverallLevel, .local_envir = parent.frame()) {
    manager <- newIntegrationStateOwner(file.path(paths$source_dir, paste0("store-", basename(tempfile()))))
    withr::defer(closeWorkflowStateOwner(manager), envir = .local_envir)
    scan <- integrationScanProject(paths$integration_dir)
    integrationCommitImport(manager, integrationImportRequest(scan, scan$views$view))
    run <- runIntegrationModel(manager, trainer = intResTrainer(messages))
    export <- manager$.context$snapshot$exports$integration_mofa_model
    list(manager = manager, run = run, object = integrationImportedObject(manager),
        model = list(reference = run$reference, parameters = export$parameters, status = run$status, trained_at = Sys.time()))
}

# The planted truth against a model's results: each planted factor matched to the model factor most correlated with
# its scores; the matched factors' group statistic and variance; the planted features among the top features.
intResPlanted <- function(results) {
    truth <- jsonlite::read_json(file.path(intResData(), "corpus", "truth.json"))
    scores <- utils::read.delim(file.path(intResData(), "corpus", "truth_scores.tsv"), stringsAsFactors = FALSE)
    loadings <- utils::read.delim(file.path(intResData(), "corpus", "truth_loadings.tsv"), stringsAsFactors = FALSE)
    values <- results$factor_values
    scores <- scores[match(values$biological_sample, scores$biological_sample), ]
    planted <- c("Z1", "Z2", "Z3")
    r <- sapply(planted, function(z) vapply(results$factors, function(f) stats::cor(values[[f]], scores[[z]]), numeric(1)))
    rownames(r) <- results$factors
    matched <- stats::setNames(results$factors[apply(abs(r), 2L, which.max)], planted)
    active <- lapply(truth$factors, function(z) unlist(lapply(unlist(z$active), function(omic) unlist(truth$views[[omic]]))))
    enrichment <- do.call(rbind, lapply(planted, function(z) do.call(rbind, lapply(active[[z]], function(view) {
        all <- results$weights[results$weights$factor == matched[[z]] & results$weights$view == view, , drop = FALSE]
        loaded <- loadings$feature_id[loadings$view == view & loadings[[z]] != 0]
        top <- all[all$rank <= results$parameters$top_features, , drop = FALSE]
        n_all <- nrow(all)
        k_all <- sum(all$feature_id %in% loaded)
        k_top <- sum(top$feature_id %in% loaded)
        data.frame(planted = z, factor = matched[[z]], view = view, top = nrow(top), top_planted = k_top,
            share_top = k_top / nrow(top), share_all = k_all / n_all,
            p = stats::phyper(k_top - 1L, k_all, n_all - k_all, nrow(top), lower.tail = FALSE), stringsAsFactors = FALSE)
    }))))
    variance <- results$variance[results$variance$factor == matched[["Z3"]], , drop = FALSE]
    list(r = r, matched = matched, association = results$association[match(matched, results$association$factor), ],
        z3_variance = stats::setNames(variance$r2_percent, variance$view), enrichment = enrichment,
        z3_active = active$Z3)
}

intResExpectPlanted <- function(results) {
    planted <- intResPlanted(results)
    message(sprintf("KIT-063 planted truth: matched %s (abs r %s)", paste(names(planted$matched), planted$matched, sep = "=", collapse = ", "),
        paste(sprintf("%.3f", apply(abs(planted$r), 2L, max)), collapse = ", ")))
    expect_identical(anyDuplicated(planted$matched), 0L)
    z1 <- planted$association[planted$association$factor == planted$matched[["Z1"]], ]
    message(sprintf("KIT-063 planted truth: %s against the group R2 %.2f%%, F %.1f, adjusted p %.3g", z1$factor, z1$r2_percent,
        z1$f_statistic, z1$p_adjusted))
    expect_identical(z1$status, "tested")
    expect_gte(z1$r2_percent, 90)
    expect_lt(z1$p_adjusted, 0.001)
    message(sprintf("KIT-063 planted truth: %s variance explained %s", planted$matched[["Z3"]],
        paste(names(planted$z3_variance), sprintf("%.2f", planted$z3_variance), sep = " ", collapse = "; ")))
    above <- names(planted$z3_variance)[planted$z3_variance >= 1]
    expect_setequal(above, planted$z3_active)
    e <- planted$enrichment
    message(paste(c("KIT-063 planted truth, top features:", sprintf("%s (%s) %s: %d of %d planted (data set %.2f), p %.2g",
        e$planted, e$factor, e$view, e$top_planted, e$top, e$share_all, e$p)), collapse = "\n"))
    for (i in seq_len(nrow(e))) {
        expect_gt(e$share_top[[i]], e$share_all[[i]])
        expect_lt(e$p[[i]], 0.05)
    }
    invisible(planted)
}

test_that("KIT-063: the Results parameters are one table, resolved and checked", {
    spec <- integrationResultsParameterSpec()
    table <- integrationResultsParameterTable()
    expect_identical(table$parameter, names(spec))
    # KIT-240 (VAL-221 (4)): every default kept by the owner's decision, which each row names.
    expect_true(all(grepl("^VAL-221 \\(decided 2026-11-25", table$settled_by)))
    factors <- paste0("Factor", 1:5)
    p <- integrationResolveResultsParameters(list(), factors)
    expect_identical(p$factors_shown, factors)
    expect_identical(p$top_features, 10L)
    expect_identical(p$p_adjust_method, "BH")
    expect_identical(integrationResolveResultsParameters(list(factors_shown = c("Factor3", "Factor1"), top_features = 5), factors)$factors_shown,
        c("Factor1", "Factor3"))
    expect_error(integrationResolveResultsParameters(list(factors_shown = "Factor9"), factors), class = "integration_parameters_invalid")
    expect_error(integrationResolveResultsParameters(list(top_features = 0), factors), class = "integration_parameters_invalid")
    expect_error(integrationResolveResultsParameters(list(weight_scale = "log"), factors), class = "integration_parameters_invalid")
})

test_that("KIT-063: the results are the stored model's, with display names, an independent group statistic and the planted truth", {
    paths <- intResProject()
    fixture <- intResModel(paths)
    value <- fixture$run$value
    results <- integrationResults(value, fixture$object, manager = fixture$manager)
    factors <- paste0("Factor", 1:5)
    expect_identical(results$factors, factors)
    expect_identical(results$views, unique(intResStored("variance_explained")$view))
    # (1) Variance explained: the product's.
    expect_equal(results$variance, intResStored("variance_explained"), tolerance = 1e-12)
    expect_identical(results$variance_total$view, value$variance_total$view)
    # (2) Factor values with the group, and the statistic against an independent computation.
    expect_identical(results$factor_values$group, value$samples$group[match(results$factor_values$biological_sample,
        value$samples$biological_sample)])
    expect_false(anyNA(results$factor_values$group))
    for (f in factors) {
        y <- results$factor_values[[f]]
        g <- factor(results$factor_values$group)
        row <- results$association[results$association$factor == f, ]
        expect_equal(row$r2_percent, 100 * summary(stats::lm(y ~ g))$r.squared, tolerance = 1e-10, info = f)
        expect_equal(row$p_value, stats::oneway.test(y ~ g, var.equal = TRUE)$p.value, tolerance = 1e-10, info = f)
        expect_identical(row$status, "tested")
    }
    expect_equal(results$association$p_adjusted, stats::p.adjust(results$association$p_value, "BH"))
    # (3) Weights: every feature, display names from the inputs, scaled to the largest, ranked; the top 10 each.
    stored <- intResStored("weights")
    expect_identical(nrow(results$weights), nrow(stored) * length(factors))
    one <- results$weights[results$weights$factor == "Factor2" & results$weights$view == "proteomics", ]
    features <- fixture$object@features$proteomics
    expect_identical(one$display_name, as.character(features$display_name[match(one$feature_id, features$feature_id)]))
    expect_identical(one$rank, seq_len(nrow(one)))
    expect_true(all(diff(abs(one$weight)) <= 0))
    expect_equal(max(abs(one$weight_scaled)), 1)
    expect_equal(one$weight, stored$Factor2[stored$view == "proteomics"][match(one$feature_id, stored$feature_id[stored$view == "proteomics"])])
    expect_identical(as.integer(table(paste(results$top_weights$factor, results$top_weights$view))),
        rep(10L, length(factors) * length(results$views)))
    # (4) Factor correlation.
    expect_equal(results$correlation, stats::cor(as.matrix(value$factors[factors])))
    # The tables as the tab shows them.
    texts <- integrationGuiTexts()
    expect_identical(names(integrationResultsVarianceTable(results)), c(texts[["col_factor"]], results$views))
    expect_identical(nrow(integrationResultsWeightsTable(results, "Factor1", "proteomics")), 10L)
    expect_identical(nrow(integrationResultsCorrelationTable(results)), 5L)
    # The figures build.
    for (plot in list(integrationPlotVariance(results), integrationPlotVarianceTotal(results), integrationPlotFactorsByGroup(results),
        integrationPlotWeights(results, "Factor3", "transcriptomics"), integrationPlotFactorCorrelation(results))) {
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
    }
    # The planted truth.
    intResExpectPlanted(results)
})

test_that("KIT-063: MOFA2's warnings are shown in plain words; D10 withholds the p-values under a group-based correction", {
    paths <- intResProject()
    fixture <- intResModel(paths, messages = c("Checking data options...", intResOverallLevel, "WARNING: something else"))
    texts <- integrationGuiTexts()
    results <- integrationResults(fixture$run$value, fixture$object, manager = fixture$manager)
    kinds <- vapply(results$advisories, `[[`, "", "kind")
    expect_identical(kinds, c("overall_level", "other"))
    expect_identical(results$advisories[[1L]]$sentence, sprintf(texts[["advisory_overall_level"]], "Factor2"))
    expect_match(results$advisories[[1L]]$detail, "average expression of features", fixed = TRUE)
    # The wording of earlier MOFA2 releases, and two factors named.
    earlier <- integrationModelAdvisories(list(messages = paste("WARNING: Factor(s) 2, 4 are strongly correlated with the",
        "total number of expressed features for at least one of your omics.")), texts)
    expect_identical(earlier[[1L]]$kind, "overall_level")
    expect_identical(earlier[[1L]]$sentence, sprintf(texts[["advisory_overall_level"]], "Factor2, Factor4"))
    expect_identical(results$advisories[[2L]]$sentence, texts[["advisory_other"]])
    # The corpus is normalised, not corrected: the factors are descriptive, the p-values given. KIT-153: the model took
    # the uncorrected matrices (the default), and the caution says that unwanted variation was not removed.
    expect_identical(results$claim$group_corrected, character())
    expect_true(results$claim$uncorrected)
    lines <- integrationClaimLimitText(results)
    expect_identical(lines[[1L]], texts[["claim_descriptive"]])
    expect_identical(lines[[2L]], texts[["claim_uncorrected"]])
    expect_identical(length(lines), 2L + length(intResOmics))
    # Proteomics corrected with the groups as replicate sets: no p-value, the reason stated.
    object <- fixture$object
    inputs <- object@inputs
    inputs$proteomics$manifest$matrix_role <- "corrected"
    inputs$proteomics$manifest$correction <- list(method = "RUVIII_C_Varying", replicate_structure = "group", ruv_k = 2L)
    object@inputs <- inputs
    # KIT-153: the claim limit follows the matrix the model used, here the corrected one.
    used <- utils::modifyList(fixture$manager$.context$snapshot$exports$integration_mofa_model$parameters,
        list(input_matrix = "ruv_corrected"))
    corrected <- integrationResults(fixture$run$value, object, manager = fixture$manager, model_parameters = used)
    expect_false(corrected$claim$uncorrected)
    # The same input under the uncorrected matrix: no group-corrected data set, the p-values given.
    expect_identical(integrationGroupCorrectedViews(object, results$views, "uncorrected")[["proteomics"]], FALSE)
    expect_identical(corrected$claim$group_corrected, "proteomics")
    expect_true(all(corrected$association$status == "not_tested"))
    expect_true(all(is.na(corrected$association$p_value)) && all(is.na(corrected$association$p_adjusted)))
    expect_equal(corrected$association$r2_percent, results$association$r2_percent)
    expect_identical(integrationClaimLimitText(corrected)[[2L]], sprintf(texts[["claim_group_corrected"]], "proteomics"))
    expect_true(all(integrationResultsAssociationTable(corrected)[[texts[["col_p"]]]] == texts[["cell_not_tested"]]))
    # Corrected with the technical replicates as replicate sets: the p-values are given.
    inputs$proteomics$manifest$correction$replicate_structure <- "technical_replicates"
    object@inputs <- inputs
    expect_identical(integrationGroupCorrectedViews(object, results$views)[["proteomics"]], FALSE)
    # A per-assay correction that was not applied to a data set leaves that data set out.
    inputs$lipidomics$manifest$matrix_role <- "corrected"
    inputs$lipidomics$manifest$correction <- list(method = "RUVIII_C_Varying", views = list(list(assay = "LCMS_Pos", applied = FALSE, k = 0L),
        list(assay = "LCMS_Neg", applied = TRUE, k = 1L)))
    object@inputs <- inputs
    flags <- integrationGroupCorrectedViews(object, results$views)
    expect_identical(unname(flags[c("lipidomics_LCMS_Pos", "lipidomics_LCMS_Neg")]), c(FALSE, TRUE))
})

test_that("KIT-063: the Results tab shows the latest model's results, saves them to the project and resets on a new model", {
    paths <- intResProject()
    fixture <- intResModel(paths)
    texts <- integrationGuiTexts()
    ui <- as.character(htmltools::renderTags(mod_integration_results_ui("r"))$html)
    for (id in c("integration-results-variance-plot", "integration-results-factors-plot", "integration-results-weights-plot",
        "integration-results-correlation-plot", "integration-results-claim", "integration-results-save")) {
        expect_match(ui, sprintf('data-testid="%s"', id), fixed = TRUE)
    }
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate(bus$state_manager <- fixture$manager)
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        expect_null(session$returned())
        expect_match(output$results_status$html, texts[["results_not_trained"]], fixed = TRUE)
        bus$integration_model <- fixture$model
        session$flushReact()
        results <- session$returned()
        expect_identical(results$factors, paste0("Factor", 1:5))
        expect_match(output$results_status$html, "Results of the model trained at", fixed = TRUE)
        advisories <- gsub("&#39;", "'", output$advisories$html, fixed = TRUE)
        expect_match(advisories, texts[["advisories_heading"]], fixed = TRUE)
        expect_match(advisories, sprintf(texts[["advisory_overall_level"]], "Factor2"), fixed = TRUE)
        expect_match(output$claim$html, texts[["claim_descriptive"]], fixed = TRUE)
        expect_match(output$association_intro$html, texts[["padjust_BH"]], fixed = TRUE)
        expect_match(output$association_table, sprintf("%.1f", results$association$r2_percent[[2L]]), fixed = TRUE)
        expect_match(gsub("&#39;", "'", output$weights_intro$html, fixed = TRUE), sprintf(texts[["weights_intro"]], 10L), fixed = TRUE)
        # The weights of the chosen factor and data set.
        session$setInputs(weights_factor = "Factor3", weights_view = "transcriptomics")
        top <- results$top_weights[results$top_weights$factor == "Factor3" & results$top_weights$view == "transcriptomics", ]
        expect_match(output$weights_table, top$feature_id[[1L]], fixed = TRUE)
        expect_match(output$weights_table, top$display_name[[10L]], fixed = TRUE)
        # Saved to the project: the tables and figures, every one written.
        session$setInputs(save = 1)
        saved <- bus$integration_results
        expect_identical(saved$model, fixture$model$reference$entity_id)
        expect_identical(saved$tables_dir, paths$mofa_tables_dir)
        expect_identical(bus$tab_status$integration_results, "complete")
        tables <- file.path(paths$mofa_tables_dir, saved$tables)
        figures <- file.path(paths$mofa_plots_dir, saved$figures)
        expect_true(all(file.size(c(tables, figures)) > 0))
        expect_true(all(c("variance_explained.tsv", "factor_group_association.tsv", "top_weights.tsv", "weights.tsv",
            "factor_values.tsv", "factor_correlation.tsv", "model_warnings.tsv") %in% saved$tables))
        # KIT-152: with the pair plot and, per factor and data set, the top features' scatter and heatmap (no map chosen).
        expect_identical(length(saved$figures), 4L + 5L * length(results$views) + 1L + 2L * 5L * length(results$views))
        written <- utils::read.delim(file.path(paths$mofa_tables_dir, "top_weights.tsv"), stringsAsFactors = FALSE)
        expect_identical(written$feature_id, results$top_weights$feature_id)
        expect_match(output$save_status$html, sprintf("%d tables to %s", length(saved$tables), paths$mofa_tables_dir), fixed = TRUE)
        # Staged under the project's scripts folder, journaled complete.
        journals <- list.files(file.path(paths$source_dir, ".workflow-export"), "journal.json", recursive = TRUE, full.names = TRUE)
        expect_length(journals, 1L)
        expect_identical(jsonlite::read_json(journals[[1L]])$status, "complete")
        expect_false(dir.exists(file.path(paths$mofa_plots_dir, ".workflow-export")))
        # Saved again: the same publication.
        session$setInputs(save = 2)
        expect_identical(list.files(file.path(paths$source_dir, ".workflow-export"), "journal.json", recursive = TRUE, full.names = TRUE),
            journals)
        # A new model: its results, not saved yet.
        other <- runIntegrationModel(fixture$manager, list(num_factors = 4L), trainer = function(views, options) {
            fit <- intResTrainer()(views, options)
            fit$factors <- fit$factors[, 1:4, drop = FALSE]
            fit$weights <- lapply(fit$weights, function(w) w[, 1:4, drop = FALSE])
            fit$r2_per_factor <- fit$r2_per_factor[1:4, , drop = FALSE]
            fit
        })
        bus$integration_model <- list(reference = other$reference, parameters = list(), status = "computed", trained_at = Sys.time())
        session$flushReact()
        expect_identical(session$returned()$factors, paste0("Factor", 1:4))
        expect_null(bus$integration_results)
        html <- function(value) paste(tryCatch(value$html, error = function(error) NULL) %||% "", collapse = "")
        expect_false(grepl("Saved at", html(output$save_status), fixed = TRUE))
        expect_false(grepl(texts[["advisories_heading"]], html(output$advisories), fixed = TRUE))
        # A new import clears the model: no results.
        bus$integration_model <- NULL
        session$flushReact()
        expect_match(output$results_status$html, texts[["results_not_trained"]], fixed = TRUE)
    })
})

test_that("KIT-063: a figure the device did not write stops the save by name, in plain words on the tab", {
    paths <- intResProject()
    fixture <- intResModel(paths)
    texts <- integrationGuiTexts()
    results <- integrationResults(fixture$run$value, fixture$object, manager = fixture$manager)
    staged <- list(mofa_tables_dir = withr::local_tempdir(), mofa_plots_dir = withr::local_tempdir())
    testthat::local_mocked_bindings(ggsave = function(filename, ...) invisible(filename), .package = "ggplot2")
    expect_error(integrationWriteResults(results, staged), class = "workflow_figure_unwritten", regexp = "were not written")
    bus <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(integrationGuiConfig()$shell))
    shiny::isolate({
        bus$state_manager <- fixture$manager
        bus$integration_model <- fixture$model
    })
    shiny::testServer(mod_integration_results_server, args = list(workflow_data = bus, paths = paths), {
        session$setInputs(save = 1)
        status <- output$save_status$html
        expect_match(status, sprintf(texts[["save_refused"]], texts[["refusal_other"]]), fixed = TRUE)
        expect_match(status, "variance_explained", fixed = TRUE)
        expect_null(bus$integration_results)
        expect_false(file.exists(file.path(paths$mofa_tables_dir, "variance_explained.tsv")))
    })
})

# KIT-151 (INT-20, audit SCI-10): the factor-against-group p-values are adjusted over every factor of the model, so the
# factors chosen for display do not change them; the family's size is kept with each row.
test_that("KIT-151 (INT-20): the adjusted p-values do not change with the factors shown", {
    values <- data.frame(biological_sample = sprintf("S%02d", 1:12), group = rep(c("A", "B"), each = 6L),
        Factor1 = c(1, 2, 1.5, 2.5, 1.2, 1.8, 3, 3.5, 2.8, 3.9, 3.1, 2.9), Factor2 = c(5, 4, 6, 5, 4, 6, 5, 6, 4, 5, 6, 4),
        Factor3 = c(0.1, 0.3, 0.2, 0.4, 0.3, 0.2, 0.5, 0.4, 0.6, 0.5, 0.4, 0.7), stringsAsFactors = FALSE)
    factors <- paste0("Factor", 1:3)
    all <- integrationFactorAssociation(values, factors)
    one <- integrationFactorAssociation(values, factors, shown = "Factor1")
    expect_identical(one$factor, "Factor1")
    expect_identical(one$p_adjusted, all$p_adjusted[[1L]])
    expect_equal(all$p_adjusted, stats::p.adjust(all$p_value, "BH"))
    expect_identical(one$adjusted_over, 3L)
    expect_true(one$p_adjusted > one$p_value)
    # Through the results: Factor1 alone shown, its adjusted p-value as with every factor shown.
    paths <- intResProject()
    fixture <- intResModel(paths)
    every <- integrationResults(fixture$run$value, fixture$object, manager = fixture$manager)
    single <- integrationResults(fixture$run$value, fixture$object, list(factors_shown = "Factor1"), manager = fixture$manager)
    expect_identical(single$association$factor, "Factor1")
    expect_identical(single$association$p_adjusted, every$association$p_adjusted[every$association$factor == "Factor1"])
    expect_identical(single$association$adjusted_over, 5L)
})

# KIT-151 (plan item E): the Kruskal-Wallis rank test as a declared alternative to the analysis of variance (still the
# default), epsilon-squared (H / (n - 1)) as its effect size, adjusted over every factor, and no p-value under the
# claim limit, as for the analysis of variance.
test_that("KIT-151: the factors can be related to the groups by the Kruskal-Wallis test, epsilon-squared as the effect", {
    expect_identical(integrationResultsParameterSpec()$association_test$default, "anova")
    expect_identical(integrationResultsParameterSpec()$association_test$choices, c("anova", "kruskal_wallis"))
    values <- data.frame(biological_sample = sprintf("S%02d", 1:12), group = rep(c("A", "B", "C"), each = 4L),
        Factor1 = c(1, 2, 1.5, 2.5, 3, 3.5, 2.8, 3.9, 5, 6, 5.5, 4.9), Factor2 = c(5, 4, 6, 5, 4, 6, 5, 6, 4, 5, 6, 4),
        stringsAsFactors = FALSE)
    factors <- c("Factor1", "Factor2")
    kw <- integrationFactorAssociation(values, factors, test = "kruskal_wallis")
    expect_identical(kw$test, c("kruskal_wallis", "kruskal_wallis"))
    for (f in factors) {
        independent <- stats::kruskal.test(values[[f]] ~ factor(values$group))
        row <- kw[kw$factor == f, ]
        expect_equal(row$h_statistic, unname(independent$statistic), tolerance = 1e-12, info = f)
        expect_equal(row$p_value, independent$p.value, tolerance = 1e-12, info = f)
        expect_identical(row$df_group, 2L)
        expect_equal(row$epsilon_squared_percent, 100 * unname(independent$statistic) / (12 - 1), tolerance = 1e-12, info = f)
        expect_true(is.na(row$r2_percent) && is.na(row$f_statistic))
    }
    expect_equal(kw$p_adjusted, stats::p.adjust(kw$p_value, "BH"))
    # The claim limit applies alike: no p-value.
    limited <- integrationFactorAssociation(values, factors, tested = FALSE, test = "kruskal_wallis")
    expect_true(all(limited$status == "not_tested") && all(is.na(limited$p_value)))
    expect_equal(limited$epsilon_squared_percent, kw$epsilon_squared_percent)
    # The analysis of variance's columns are unchanged in meaning.
    anova <- integrationFactorAssociation(values, factors)
    expect_identical(anova$test, c("anova", "anova"))
    expect_true(all(is.na(anova$h_statistic)))
    # Through the results and the tab's table.
    paths <- intResProject()
    fixture <- intResModel(paths)
    results <- integrationResults(fixture$run$value, fixture$object, list(association_test = "kruskal_wallis"),
        manager = fixture$manager)
    expect_identical(unique(results$association$test), "kruskal_wallis")
    texts <- integrationGuiTexts()
    expect_identical(names(integrationResultsAssociationTable(results, texts)), unname(texts[c("col_factor", "col_samples_tested",
        "col_groups_n", "col_epsilon2", "col_h", "col_df", "col_p", "col_p_adjusted")]))
})
