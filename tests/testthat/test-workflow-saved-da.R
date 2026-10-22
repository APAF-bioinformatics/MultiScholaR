makeSavedDaInputs <- function(technical_replicates = TRUE) {
    set.seed(12012L)
    samples <- sprintf("S%02d", seq_len(8L))
    features <- sprintf("P%03d", seq_len(200L))
    group <- factor(rep(c("control", "treated"), each = 4L),
        levels = c("control", "treated"))
    design <- data.frame(group = group, row.names = samples)
    if (technical_replicates) design$replicates <- rep(c(1L, 1L, 2L, 2L), 2L)
    quantities <- matrix(stats::rnorm(length(features) * length(samples)),
        nrow = length(features), dimnames = list(features, samples))
    quantities[seq.int(10L, length(features), by = 10L), group == "treated"] <-
        quantities[seq.int(10L, length(features), by = 10L), group == "treated"] + 2
    list(data = quantities, design = design, samples = samples, features = features)
}

buildSavedDaArtifacts <- function(inputs, contrast, treat_lfc = NA_real_,
                                  counters = NULL) {
    formula <- stats::as.formula("~0+group")
    model_frame <- stats::model.frame(formula, inputs$design)
    model <- stats::model.matrix(formula, model_frame)
    data <- inputs$data[, rownames(model), drop = FALSE]
    count <- if (is.null(counters)) \(operation) invisible(NULL) else
        \(operation) counters[[operation]] <- (counters[[operation]] %||% 0L) + 1L
    base <- fitProtDaBaseModel(data, model, inputs$design, rownames(model),
        counterFn = count, messageFn = \(...) invisible(NULL))
    contrast_matrix <- limma::makeContrasts(contrasts = contrast,
        levels = colnames(model))
    contrasted <- applyProtDaContrastModel(base, contrast_matrix,
        treat_lfc_cutoff = treat_lfc, counterFn = count,
        messageFn = \(...) invisible(NULL))
    tables <- buildProtDaResultTables(contrasted, contrast,
        counterFn = count, messageFn = \(...) invisible(NULL))
    list(base = base, contrasted = contrasted, tables = tables)
}

test_that("the reference DA GUI publishes native-equivalent tables and serves projected views", {
    inputs <- makeSavedDaInputs()
    design <- inputs$design
    design$Run <- rownames(design)
    table <- data.frame(Protein.Group = inputs$features, inputs$data, check.names = FALSE)
    rownames(table) <- NULL
    annotation <- data.frame(Protein.Group = inputs$features, gene_names = inputs$features)
    contrasts <- normaliseProtDaContrastsTable(data.frame(contrasts = "grouptreated-groupcontrol"))
    object <- methods::new("ProteinQuantitativeData", protein_quant_table = table,
        protein_id_table = annotation, protein_id_column = "Protein.Group", design_matrix = design,
        sample_id = "Run", group_id = "group", technical_replicate_id = "replicates")
    parameters <- list(contrasts_tbl = data.frame(contrasts = contrasts$full_format),
        formula_string = "~0+group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0,
        eBayes_trend = TRUE, eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
        args_row_id = "Protein.Group", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
    object@args$differentialAbundanceAnalysis <- parameters
    object@args$differentialAbundanceAnalysisHelper <- parameters
    rlang::local_bindings(config_list = object@args, .env = .GlobalEnv)
    expected <- suppressMessages(differentialAbundanceAnalysis(object,
        contrasts_tbl = parameters$contrasts_tbl, formula_string = "~0+group",
        da_q_val_thresh = 0.05, treat_lfc_cutoff = 0))
    object@args <- list()
    prior_plot <- ggplot2::get_last_plot()
    withr::defer(ggplot2::set_last_plot(prior_plot))
    for (kind in c("volcano", "heatmap", "glimma")) {
        ggplot2::set_last_plot(NULL)
        lifetime <- new.env(parent = emptyenv())
        lifetime$collected <- FALSE
        finalize <- function(value) lifetime$collected <- TRUE
        environment(finalize) <- list2env(list(lifetime = lifetime), parent = baseenv())
        display <- local({
            sentinel <- new.env(parent = emptyenv())
            reg.finalizer(sentinel, finalize, onexit = FALSE)
            value <- list(da_proteins_long = expected$da_proteins_long,
                theObject = object, ownership_sentinel = sentinel)
            value$da_proteins_long$Protein.Ids <- value$da_proteins_long$Protein.Group
            switch(kind,
                volcano = generateProtDAVolcanoStatic(value, contrasts$friendly_names),
                heatmap = generateProtDAHeatmap(value, contrasts$friendly_names,
                    top_n_genes = 20L, da_q_val_thresh = 1),
                glimma = generateProtDAVolcanoPlotGlimma(value, contrasts$friendly_names,
                    args_row_id = "Protein.Group"))
        })
        gc(full = TRUE)
        gc(full = TRUE)
        expect_true(lifetime$collected, info = paste(kind, "retains source result owner"))
        expect_false(is.null(display))
        if (kind == "volcano") {
            layers <- ggplot2::ggplot_build(display)$data
            expect_equal(layers[[1L]]$x, expected$da_proteins_long$log2FC, tolerance = 1e-12)
            expect_equal(layers[[1L]]$y, -log10(expected$da_proteins_long$fdr_qvalue), tolerance = 1e-12)
            expect_null(ggplot2::get_last_plot())
        }
        if (kind == "heatmap") expect_identical(nrow(methods::slot(display$plot, "matrix")), 20L)
        if (kind == "glimma") expect_s3_class(display, "htmlwidget")
        display <- NULL
    }
    root <- withr::local_tempdir()
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path)
    shiny::testServer(function(input, output, session) {
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "reference-da-gui", "DIA"),
            WorkflowState$new(audit_enabled = FALSE))
        session$onSessionEnded(function() closeDIAStateAdapter(manager))
        manager$.context$reference_states <- TRUE
        manager$.context$experiment_paths <- paths
        session$userData$dia_capacity_manager <- manager
        manager$saveState("correlation_filtered", object, list(), "input")
        workflow <- shiny::reactiveValues(state_manager = manager, contrasts_tbl = contrasts)
        da <- shiny::reactiveValues(current_s4_object = getDIAStateReference(manager))
        shiny::observeEvent(input$go, {
            testthat::with_mocked_bindings(runDIAReferenceDaAnalysis(
                list(formula_string = "~0+group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
                da, workflow, paths),
                materializeDIAStateRecord = function(...) stop("DA hydrated a full state"))
        })
    }, {
        session$setInputs(go = 1L)
        session$flushReact()
        expect_true(isDIAResultManifest(da$da_results_list))
        # KIT-025 audit PH-01: the contrasts commit no state, so the DA step leaves no provenance pending.
        expect_null(manager$.context$pending_provenance)
        expect_identical(manager$.context$snapshot$step_provenance$correlation_filtered$operation, "legacy")
        interactive_files <- list.files(file.path(paths$publication_graphs_dir,
            "Interactive_Volcano_Plots"), recursive = TRUE, full.names = TRUE)
        expect_true(length(interactive_files) > 0L && all(file.size(interactive_files) > 0))
        expect_length(interactive_files[endsWith(interactive_files, ".html")], nrow(contrasts))
        expect_identical(tail(names(readDIAStateMetadata(manager)$data$slots$args), 2L),
            c("differentialAbundanceAnalysis", "daAnalysisUI"))
        selected <- selectDIAContrastReference(da$da_results_list, contrasts$friendly_names)
        expect_equal(readDIAContrastOutput(selected, "da_proteins_long")$da_proteins_long,
            expected$da_proteins_long, tolerance = 1e-12)
        page <- readDIAResultPage(da$da_results_list, contrasts$friendly_names, limit = 10L)
        expect_equal(nrow(page$data), 10L)
        expect_equal(page$counts$total, nrow(expected$da_proteins_long))
        native_heatmap <- generateProtDAHeatmap(list(da_proteins_long = expected$da_proteins_long,
            theObject = object), contrasts$friendly_names, top_n_genes = 20L,
            da_q_val_thresh = 1, tree_cut_method = "k_clusters", n_clusters = 3L)
        reader <- materializeDIATable
        source_table <- readDIAStateMetadata(manager, da$da_results_list$source)$data$components$protein_quant_table
        heatmap <- testthat::with_mocked_bindings(
            generateProtDAHeatmap(da$da_results_list, contrasts$friendly_names, top_n_genes = 20L,
                da_q_val_thresh = 1, tree_cut_method = "k_clusters", n_clusters = 3L),
            materializeDIATable = function(context, reference, ...) {
                if (identical(reference, source_table) && is.null(list(...)$relation)) stop("full heatmap assay read")
                reader(context, reference, ...)
            })
        expect_equal(methods::slot(heatmap$plot, "matrix"), methods::slot(native_heatmap$plot, "matrix"),
            tolerance = 1e-12)
        expect_identical(heatmap$row_clusters, native_heatmap$row_clusters)
        styled <- testthat::with_mocked_bindings(
            generateProtDAHeatmap(da$da_results_list, contrasts$friendly_names, top_n_genes = 20L,
                da_q_val_thresh = 1, tree_cut_method = "k_clusters", n_clusters = 3L,
                color_scheme = "viridis", show_gene_names = TRUE),
            materializeDIATable = function(...) stop("heatmap product hit read an assay"),
            computeProtDAHeatmapProduct = function(...) stop("heatmap reclustered for colour"))
        expect_identical(styled$row_clusters, heatmap$row_clusters)
        testthat::with_mocked_bindings({
            session$setInputs(go = 2L)
            session$flushReact()
        }, differentialAbundanceAnalysis = function(...) stop("DA science repeated on a valid result hit"))
        expect_true(da$analysis_complete)
        expect_identical(manager$.context$full_s4_reads, 0L)
        expect_identical(manager$.context$store$open_resources, 0L)
    })
})

test_that("the reference DA dispatch reuses base and contrast fits with weights NA", {
    inputs <- makeSavedDaInputs(technical_replicates = TRUE)
    contrast <- "treated_vs_control=grouptreated-groupcontrol"
    expected <- suppressMessages(runTestsContrasts(inputs$data, contrast, inputs$design,
        "~0+group", weights = NA, treat_lfc_cutoff = 0))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "reference-fit", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$.context$da_reference_execution <- TRUE
    testthat::with_mocked_bindings({
        result <- suppressMessages(runTestsContrasts(inputs$data, contrast, inputs$design,
            "~0+group", weights = NA, treat_lfc_cutoff = 0))
        expect_equal(result, expected, tolerance = 1e-12)
        repeated <- testthat::with_mocked_bindings(
            suppressMessages(runTestsContrasts(inputs$data, contrast, inputs$design,
                "~0+group", weights = NA, treat_lfc_cutoff = 0)),
            fitProtDaBaseModel = function(...) stop("base fit repeated"),
            applyProtDaContrastModel = function(...) stop("contrast fit repeated"))
        expect_identical(repeated, result)
        expect_identical(manager$.context$store$open_resources, 0L)
        expect_identical(manager$.context$full_s4_reads, 0L)
    }, getDIAActiveManager = function() manager)
})

test_that("reusable DA layers agree with public no-TREAT and TREAT results", {
    inputs <- makeSavedDaInputs(technical_replicates = TRUE)
    contrast <- "treated_vs_control=grouptreated-groupcontrol"
    for (treat_lfc in c(NA_real_, 0.5)) {
        counters <- new.env(parent = emptyenv())
        layered <- buildSavedDaArtifacts(inputs, contrast, treat_lfc, counters)
        public <- suppressMessages(runTestsContrasts(inputs$data, contrast,
            inputs$design, "~0+group", weights = NA,
            treat_lfc_cutoff = treat_lfc))
        expect_identical(layered$tables$results, public$results)
        expect_identical(layered$tables$qvalue_warnings, public$qvalue_warnings)
        expect_equal(layered$tables$fit.eb$coefficients,
            public$fit.eb$coefficients, tolerance = 0)
        expect_identical(layered$base$feature_ids, inputs$features)
        expect_identical(layered$base$sample_ids, inputs$samples)
        expect_true(layered$base$uses_duplicate_correlation)
        expect_identical(counters$duplicateCorrelation, 1L)
        expect_identical(counters$lmFit, 1L)
        expect_identical(counters$contrasts.fit, 1L)
        expect_identical(counters$eBayes, 1L)
        expect_identical(counters$treat %||% 0L, if (is.na(treat_lfc)) 0L else 1L)
        expect_identical(counters$topTable %||% 0L, if (is.na(treat_lfc)) 1L else 0L)
        expect_identical(counters$topTreat %||% 0L, if (is.na(treat_lfc)) 0L else 1L)
        expect_identical(counters[["p.adjust:fdr"]], 1L)
    }
})

test_that("base-fit reuse recomputes contrast layers without fitting again", {
    inputs <- makeSavedDaInputs(technical_replicates = FALSE)
    counters <- new.env(parent = emptyenv())
    first <- buildSavedDaArtifacts(inputs,
        "treated_vs_control=grouptreated-groupcontrol", NA_real_, counters)
    expect_identical(counters$lmFit, 1L)
    model <- first$base$fit$design
    reverse_matrix <- limma::makeContrasts(
        contrasts = "control_vs_treated=groupcontrol-grouptreated",
        levels = colnames(model))
    reverse <- applyProtDaContrastModel(first$base, reverse_matrix,
        counterFn = \(operation) counters[[operation]] <-
            (counters[[operation]] %||% 0L) + 1L,
        messageFn = \(...) invisible(NULL))
    reverse_tables <- buildProtDaResultTables(reverse,
        "control_vs_treated=groupcontrol-grouptreated",
        counterFn = \(operation) counters[[operation]] <-
            (counters[[operation]] %||% 0L) + 1L,
        messageFn = \(...) invisible(NULL))
    expect_identical(counters$lmFit, 1L)
    expect_identical(counters$contrasts.fit, 2L)
    expect_identical(counters$eBayes, 2L)
    expect_equal(reverse_tables$results[[1L]]$logFC,
        -first$tables$results[[1L]]$logFC, tolerance = 1e-12)
})

test_that("adjustment masks and qvalue fallback retain full populations", {
    table <- data.frame(logFC = c(1, 2, 3, 4),
        P.Value = c(0.01, NA_real_, Inf, 0.2), row.names = paste0("P", 1:4))
    counters <- new.env(parent = emptyenv())
    adjusted <- adjustProtDaResultTable(table, "raw_pvalue", "fdr_qvalue",
        "fdr_value_bh_adjustment", qvalueFn = \(values) stop("forced qvalue failure"),
        counterFn = \(operation) counters[[operation]] <-
            (counters[[operation]] %||% 0L) + 1L,
        messageFn = \(...) invisible(NULL))
    expect_true(adjusted$qvalue_failed)
    expect_identical(adjusted$valid_p_value_mask, c(TRUE, FALSE, FALSE, TRUE))
    expected <- rep(NA_real_, 4L)
    expected[c(1L, 4L)] <- stats::p.adjust(c(0.01, 0.2), method = "BH")
    expect_identical(adjusted$table$fdr_qvalue, expected)
    expect_identical(adjusted$table$fdr_value_bh_adjustment, expected)
    expect_identical(adjusted$table$raw_pvalue,
        c(0.01, NA_real_, Inf, 0.2))
    expect_identical(counters[["p.adjust:qvalue_fallback"]], 1L)
    expect_identical(counters[["p.adjust:fdr"]], 1L)

    invalid <- table
    invalid$P.Value <- c(NA_real_, NaN, Inf, -Inf)
    untouched <- adjustProtDaResultTable(invalid, "raw_pvalue", "fdr_qvalue",
        "fdr_value_bh_adjustment")
    expect_false(untouched$qvalue_failed)
    expect_true(all(is.na(untouched$table$fdr_qvalue)))
    expect_identical(nrow(untouched$table), nrow(invalid))
})

test_that("saved base fits and full results invalidate at the correct layer", {
    inputs <- makeSavedDaInputs(technical_replicates = TRUE)
    contrast <- "treated_vs_control=grouptreated-groupcontrol"
    layered <- buildSavedDaArtifacts(inputs, contrast)
    software <- "MultiScholaR-da-test-v1"
    base_identity <- newWorkflowDiaBaseFitIdentity(inputs$data, inputs$design,
        "~0+group", layered$base$block, NA, "cyclicloess", "ruv-skip", software)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "da-project", "DIA",
        max_cache_bytes = 67108864)
    base_ref <- saveWorkflowDiaBaseFit(store, layered$base, "base-fit", base_identity)
    base_hit <- resolveWorkflowDiaBaseFit(store, base_identity,
        layered$base$uses_duplicate_correlation,
        layered$base$consensus_correlation)
    expect_identical(base_hit$status, "hit")
    expect_equal(base_hit$value$fit$coefficients,
        layered$base$fit$coefficients, tolerance = 0)

    result_identity <- newWorkflowDiaResultIdentity(base_ref, contrast, NA_real_,
        FALSE, FALSE, "qvalue-then-BH-v1", software)
    result_ref <- saveWorkflowDiaResult(store, layered$tables, "full-result",
        result_identity)
    releaseWorkflowEntityCache(store)
    before <- workflowEntityCacheInfo(store)
    result_hit <- resolveWorkflowDiaResult(store, result_identity)
    after <- workflowEntityCacheInfo(store)
    expect_identical(result_hit$status, "hit")
    expect_identical(result_hit$value$value$results, layered$tables$results)
    expect_identical(after$payload_reads, before$payload_reads + 1L)

    changed_contrast <- newWorkflowDiaResultIdentity(base_ref,
        "control_vs_treated=groupcontrol-grouptreated", NA_real_, FALSE, FALSE,
        "qvalue-then-BH-v1", software)
    expect_identical(resolveWorkflowDiaResult(store, changed_contrast)$status, "miss")
    changed_data <- inputs$data
    changed_data[[1L, 1L]] <- changed_data[[1L, 1L]] + 1
    changed_base <- newWorkflowDiaBaseFitIdentity(changed_data, inputs$design,
        "~0+group", layered$base$block, NA, "cyclicloess", "ruv-skip", software)
    expect_identical(resolveWorkflowDiaBaseFit(store, changed_base,
        layered$base$uses_duplicate_correlation,
        layered$base$consensus_correlation)$status, "miss")
    weighted_base <- newWorkflowDiaBaseFitIdentity(inputs$data, inputs$design,
        "~0+group", layered$base$block, matrix(1, nrow = 8L, ncol = 1L),
        "cyclicloess", "ruv-skip", software)
    expect_false(identical(weighted_base, base_identity))
    expect_identical(base_ref$kind, "da_fit")
    expect_identical(result_ref$kind, "scientific_product")
})

test_that("table volcano and heatmap views retain full-population statistics", {
    inputs <- makeSavedDaInputs(technical_replicates = FALSE)
    contrast <- "treated_vs_control=grouptreated-groupcontrol"
    layered <- buildSavedDaArtifacts(inputs, contrast)
    envelope <- newSavedWorkflowResultTable(layered$tables$results[[1L]],
        ".rownames", contrast, "fdr_qvalue")
    page <- pageWorkflowDiaResult(envelope, page = 1L, page_size = 25L)
    expect_identical(nrow(page$data), 25L)
    expect_identical(page$full_population_rows, 200L)
    expect_false(page$statistics_recomputed)
    expect_identical(page$data$fdr_qvalue,
        sort(layered$tables$results[[1L]]$fdr_qvalue, na.last = TRUE)[1:25])

    volcano <- filterWorkflowDiaVolcano(envelope, 0.05, 1)
    expect_identical(length(volcano$selected), 200L)
    expect_identical(volcano$full_population_rows, 200L)
    expect_true(all(volcano$data$fdr_qvalue <= 0.05))
    expect_true(all(abs(volcano$data$logFC) >= 1))
    expect_false(volcano$statistics_recomputed)

    heatmap <- selectWorkflowDiaHeatmap(envelope, inputs$data, 50L,
        scaling = "row", clustering = "both")
    expect_identical(dim(heatmap$matrix), c(50L, 8L))
    expect_identical(heatmap$sample_ids, inputs$samples)
    expect_identical(rownames(heatmap$matrix), heatmap$feature_ids)
    expect_identical(heatmap$full_population_rows, 200L)
    expect_identical(heatmap$full_tested_features, 200L)
    expect_identical(heatmap$cluster_method, "ward.D2")
    expect_identical(heatmap$distance_method, "euclidean")
    expect_identical(heatmap$tree_cut_method, "k_clusters")
    expect_identical(heatmap$n_clusters, 4L)
    expect_false(heatmap$statistics_recomputed)
    expect_true(all(heatmap$feature_ids %in% rownames(inputs$data)))
    selected <- selectWorkflowDiaExpression(inputs$data, heatmap$feature_ids[[1L]])
    expect_identical(selected$sample_id, inputs$samples)
    expect_identical(selected$value,
        unname(inputs$data[heatmap$feature_ids[[1L]], ]))
    expect_true(all(selected$feature_id == heatmap$feature_ids[[1L]]))
})

test_that("Limpa identity includes uncertainty support and forbids projection", {
    inputs <- makeSavedDaInputs(technical_replicates = FALSE)
    support <- matrix(0.2, nrow(inputs$data), ncol(inputs$data),
        dimnames = dimnames(inputs$data))
    observations <- matrix(3L, nrow(inputs$data), ncol(inputs$data),
        dimnames = dimnames(inputs$data))
    elist <- methods::new("EList", list(E = inputs$data,
        genes = data.frame(protein.id = rownames(inputs$data),
            row.names = rownames(inputs$data)),
        other = list(standard.error = support, n.observations = observations)))
    changed <- elist
    changed$other$standard.error[[1L, 1L]] <- 0.4
    original_id <- newWorkflowLimpaFitIdentity(elist, inputs$design,
        "~0+group", "limpa-test-v1")
    changed_id <- newWorkflowLimpaFitIdentity(changed, inputs$design,
        "~0+group", "limpa-test-v1")
    expect_identical(elist$E, changed$E)
    expect_false(identical(original_id, changed_id))

    layered <- buildSavedDaArtifacts(inputs,
        "treated_vs_control=grouptreated-groupcontrol")
    base_identity <- newWorkflowDiaBaseFitIdentity(inputs$data, inputs$design,
        "~0+group", layered$base$block, NA, "none", "ruv-skip", "test-v1")
    store <- newWorkflowEntityStore(withr::local_tempdir(), "projection-project", "DIA")
    reference <- saveWorkflowDiaBaseFit(store, layered$base, "base", base_identity)
    before <- workflowEntityCacheInfo(store)
    expect_error(hydrateSavedWorkflowProjection(store, reference,
        c("coefficients")), class = "saved_entity_unsupported")
    after <- workflowEntityCacheInfo(store)
    expect_identical(after$payload_reads, before$payload_reads)
    expect_false(savedWorkflowEntityCapabilities(reference)$column_projection)
})

test_that("invalidateDIAVisibleResults keeps accepted results on every DA-ready state, the DIA_LIMPA route's ruv_covariates included (KIT-053)", {
  fake_manager <- function(state) structure(list2env(list(current_state = state,
    .context = list(reference_states = TRUE))), class = "DIAStateAdapter")
  outcome <- function(state) {
    da_data <- list2env(list(da_results_list = list(kept = TRUE), analysis_complete = TRUE,
      current_s4_object = list(object = TRUE)))
    workflow <- list2env(list(state_manager = fake_manager(state), da_analysis_results_list = list(kept = TRUE)))
    invalidateDIAVisibleResults(da_data, workflow)
    list(complete = da_data$analysis_complete, results = !is.null(da_data$da_results_list), object = !is.null(da_data$current_s4_object))
  }
  for (state in c("ruv_covariates", "correlation_filtered", "ruv_corrected", "normalized")) {
    expect_identical(outcome(state), list(complete = TRUE, results = TRUE, object = TRUE), info = state)
  }
  expect_identical(outcome("peptide_count_filtered"), list(complete = FALSE, results = FALSE, object = FALSE))
})
