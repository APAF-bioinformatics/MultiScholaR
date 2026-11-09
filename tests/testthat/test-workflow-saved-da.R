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

# ------------------------------------------------------------------------------------------------
# KIT-112 W1 (a): `runDIAReferenceDaAnalysis` in process, at core level, on D6's cases, captured on the start commit
# (64e7c2ae) before the runner is split into its request, compute and accept halves. Each case's capture: the owner's
# and the head's current state and history; each committed state's record (its science, its `@args` with the DA
# arguments' time dropped, config, description, previous state and product names); the completed outputs; the head's
# products by fingerprint; the store's entities by kind and layer without the head snapshots (`state-head`) and the
# fits' factors; the runner's manifest and the head's, the model's unset fields; `da_dependencies` by key and value; the
# step provenance and the pending provenance; the analysis record's `differential_abundance`; the published files;
# whether the head's publication request is keyed by `diaDaPublicationIdentity` of the manifest; the display fields;
# the value or refusal; and the call counts of the halves' callees. Store-local ids, times and the root are masked.
# MULTISCHOLAR_KIT112_W1_DUMP=<file> appends each case's capture as R code (`1` prints it): how the frozen values were made.
# ------------------------------------------------------------------------------------------------

savedDaW1Callees <- c("runSavedOperation", "commitDIAReferenceState", "saveDIAStateValue", "publishDIAActionOutputs")

savedDaW1Sha <- function(x) digest::digest(x, algo = "sha256")

#' Store-local ids (32 to 64 hexadecimal characters) and the fixture's root masked in every character value of plain lists.
savedDaW1Mask <- function(x, root) {
    if (is.character(x)) return(gsub("[a-f0-9]{32,64}", "<id>", gsub(root, "<root>", x, fixed = TRUE)))
    if (is.list(x) && !is.object(x)) {
        x[] <- lapply(x, savedDaW1Mask, root = root)
        if (!is.null(names(x))) names(x) <- savedDaW1Mask(names(x), root)
    }
    x
}

#' A saved checkpoint's science: its slots, and its `@args` canonical (so the DA arguments' commit is seen, its time not).
savedDaW1State <- function(manager, reference, root, memo) {
    withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        list(slots = lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot)),
            args = savedDaW1Canonical(object@args, manager, root, memo))
    })
}

#' A value without store-local ids: a checkpoint by its science, a table by its rows, any other saved entity by its
#' hydrated envelope; plots by class (a ggplot also by its data and title); bookkeeping fields dropped.
savedDaW1Canonical <- function(value, manager, root, memo, depth = 0L) {
    if (depth > 32L) return("<depth>")
    if (is.null(value)) return(NULL)
    if (is.function(value) || is.environment(value)) return("<opaque>")
    if (inherits(value, c("gg", "ggplot"))) return(list(.gg = class(value),
        data = savedDaW1Sha(as.data.frame(value$data)), title = value$labels$title))
    if (inherits(value, c("grob", "gtable", "htmlwidget", "recordedplot"))) return(paste0("<", class(value)[[1L]], ">"))
    if (isSavedWorkflowReference(value)) {
        key <- value$entity_id
        if (!is.null(memo[[key]])) return(memo[[key]])
        out <- if (identical(value$kind, "checkpoint")) {
            tryCatch(list(.state = savedDaW1Sha(savedDaW1State(manager, value, root, memo))), error = function(e)
                list(.checkpoint = savedDaW1Canonical(resolveWorkflowStateValue(manager$.context, value), manager, root, memo,
                    depth + 1L)))
        } else if (identical(value$payload$codec, "parquet_duckdb")) {
            table <- as.data.frame(materializeWorkflowTable(manager$.context, value))
            attr(table, "dia_source_order") <- NULL
            rownames(table) <- NULL
            list(.table = savedDaW1Sha(table))
        } else if (identical(value$kind, "rendered_asset")) {
            list(.rendered = value$kind)
        } else {
            list(.entity = value$kind, value = savedDaW1Canonical(hydrateSavedWorkflowEntity(manager$.context$store, value),
                manager, root, memo, depth + 1L))
        }
        memo[[key]] <- out
        return(out)
    }
    if (isS4(value)) return(list(.s4 = class(value)[[1L]]))
    if (is.data.frame(value)) {
        value <- as.data.frame(value)
        rownames(value) <- NULL
        return(list(.df = savedDaW1Sha(savedDaW1Mask(as.list(value), root)), dim = dim(value)))
    }
    if (is.list(value)) {
        keep <- !(names(value) %||% rep("", length(value))) %in% c("timestamp", "mtime", "dependencies",
            "product_dependencies", "owned_components", "component_dependencies", "committed_at")
        return(lapply(value[keep], savedDaW1Canonical, manager = manager, root = root, memo = memo, depth = depth + 1L))
    }
    savedDaW1Mask(value, root)
}

#' One saved reference's fingerprint: its product type (or kind) and 16 characters of its canonical value's digest.
savedDaW1Fingerprint <- function(reference, manager, root, memo) {
    if (identical(reference$kind, "rendered_asset")) return("rendered_asset")
    canonical <- savedDaW1Canonical(reference, manager, root, memo)
    type <- canonical$value$product_type
    if (!is.character(type)) type <- reference$kind
    paste0(type, ":", substr(savedDaW1Sha(canonical), 1L, 16L))
}

#' One committed state's record: its science, config, description, previous state and product names.
savedDaW1Record <- function(manager, reference, root, memo) {
    record <- readDIAStateMetadata(manager, reference)
    list(state = savedDaW1State(manager, reference, root, memo), config = savedDaW1Canonical(record$config, manager, root, memo),
        description = record$description, previous = record$previous_state, products = sort(names(record$data$products)))
}

#' The store's committed entities read from disk (a reopened owner indexes only its roots), without the head snapshots.
savedDaW1Descriptors <- function(manager) {
    files <- sort(list.files(file.path(manager$.context$store$directory, "entities"), pattern = "^descriptor[.]rds$",
        recursive = TRUE, full.names = TRUE))
    records <- lapply(files, readRDS)
    Filter(function(record) !identical(record$identity$parameters, "state-head"), records)
}

#' Entity counts by kind, with a checkpoint's label or a product's layer.
savedDaW1Entities <- function(records) {
    keys <- vapply(records, function(record) paste0(record$kind, if (identical(record$kind, "checkpoint"))
        paste0(":", record$identity$parameters) else if (!is.null(record$metadata$factors$layer))
        paste0(":", record$metadata$factors$layer) else ""), "")
    counts <- table(keys)
    stats::setNames(as.integer(counts), names(counts))
}

#' The base fits and contrasted fits by layer and 16 characters of their factors' digest.
savedDaW1Fits <- function(records) {
    fits <- Filter(function(record) identical(record$kind, "da_fit") ||
        identical(record$metadata$factors$layer, "full_results"), records)
    sort(vapply(fits, function(record) paste0(record$kind, ":", record$metadata$factors$layer, ":",
        substr(savedDaW1Sha(record$metadata$factors), 1L, 16L)), ""))
}

#' `da_dependencies` by key (the head's manifest, a head product's fingerprint, or another entity) and value.
savedDaW1Dependencies <- function(manager, snapshot, root, memo) {
    dependencies <- snapshot$da_dependencies
    if (!length(dependencies)) return(character())
    products <- snapshot$products
    ids <- vapply(products, function(reference) reference$entity_id, "")
    key <- function(id) {
        if (identical(id, snapshot$da_results$entity_id)) return("manifest")
        hit <- products[ids == id]
        if (length(hit)) savedDaW1Fingerprint(hit[[1L]], manager, root, memo) else "other"
    }
    sort(vapply(names(dependencies), function(id) paste0(key(id), " <- ", substr(savedDaW1Sha(sort(unname(vapply(
        dependencies[[id]], savedDaW1Fingerprint, "", manager = manager, root = root, memo = memo)))), 1L, 16L)), "",
        USE.NAMES = FALSE))
}

#' Files under `root` (relative names), with the md5 of each text table; images, PDFs, workbooks and pages by name only.
savedDaW1Files <- function(root) {
    files <- sort(list.files(root, recursive = TRUE))
    files <- files[!grepl("preflight", files, fixed = TRUE)]
    tables <- files[grepl("[.](tsv|txt|tab)$", files)]
    list(files = files, tables = stats::setNames(unname(tools::md5sum(file.path(root, tables))), tables))
}

#' Run one click with the halves' callees counted; a refusal is returned as its message.
savedDaW1Run <- function(code) {
    calls <- stats::setNames(integer(length(savedDaW1Callees)), savedDaW1Callees)
    counters <- lapply(savedDaW1Callees, function(name) {
        original <- get(name, envir = asNamespace("MultiScholaR"))
        function(...) {
            calls[[name]] <<- calls[[name]] + 1L
            original(...)
        }
    })
    names(counters) <- savedDaW1Callees
    outcome <- NULL
    utils::capture.output(outcome <- rlang::inject(testthat::with_mocked_bindings(
        tryCatch(list(value = suppressWarnings(suppressMessages(code))), error = function(e) list(error = conditionMessage(e))),
        !!!counters, .package = "MultiScholaR")))
    c(outcome, list(calls = calls))
}

#' The capture after one click (see the section's header).
savedDaW1Capture <- function(fixture, run) {
    force(run)
    manager <- fixture$workflow$state_manager
    root <- fixture$root
    memo <- new.env(parent = emptyenv())
    canonical <- function(value) savedDaW1Canonical(value, manager, root, memo)
    sha16 <- function(value) substr(savedDaW1Sha(canonical(value)), 1L, 16L)
    # The head as an action scope reads it (the context's copy outside a scope can lag the last commit).
    snapshot <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
    history <- as.character(unlist(snapshot$state_history))
    committed <- setdiff(history, "initial")
    descriptors <- savedDaW1Descriptors(manager)
    manifest <- if (isDIAResultManifest(run$value)) run$value
    head_manifest <- if (isSavedWorkflowReference(snapshot$da_results)) resolveWorkflowStateValue(manager$.context, snapshot$da_results)
    head <- snapshot$states[[snapshot$current_state]]
    request <- snapshot$action_outputs[[head$entity_id]]
    record <- tryCatch(buildWorkflowAnalysisRecord(manager)$analyses$differential_abundance,
        error = function(e) paste("error:", conditionMessage(e)))
    list(
        current = manager$current_state,
        history = as.character(unlist(manager$state_history)),
        head_current = snapshot$current_state,
        head_history = history,
        states = vapply(stats::setNames(nm = committed), function(state)
            savedDaW1Sha(savedDaW1Record(manager, snapshot$states[[state]], root, memo)), ""),
        outputs = lapply(stats::setNames(nm = committed), function(state)
            unlist(snapshot$completed_outputs[[snapshot$states[[state]]$entity_id]])),
        products = sort(vapply(unname(snapshot$products), savedDaW1Fingerprint, "", manager = manager, root = root, memo = memo)),
        entities = savedDaW1Entities(descriptors),
        fits = savedDaW1Fits(descriptors),
        manifest = if (is.null(manifest)) NA_character_ else savedDaW1Sha(canonical(manifest)),
        model_unset = if (is.null(manifest)) NA_character_ else names(manifest$model)[vapply(manifest$model, is.null, NA)],
        head_manifest = if (is.null(head_manifest)) NA_character_ else savedDaW1Sha(canonical(head_manifest)),
        head_model_unset = if (is.null(head_manifest)) NA_character_ else
            names(head_manifest$model)[vapply(head_manifest$model, is.null, NA)],
        dependencies = savedDaW1Dependencies(manager, snapshot, root, memo),
        provenance = savedDaW1Sha(savedDaW1Mask(lapply(snapshot$step_provenance, function(entry)
            entry[setdiff(names(entry), c("committed_at", "build", "source_commit"))]), root)),
        pending = savedDaW1Sha(savedDaW1Mask(manager$.context$pending_provenance, root)),
        record = savedDaW1Sha(savedDaW1Mask(record, root)),
        files = savedDaW1Sha(savedDaW1Files(root)),
        publication = if (is.null(request) || is.null(manifest)) NA else identical(request$request_id, paste0("dia-da-",
            substr(digest::digest(list(state = head$entity_id, action = "da", identity = diaDaPublicationIdentity(manifest)),
                algo = "sha256"), 1L, 32L))),
        display = c(da_results_list = sha16(fixture$da$da_results_list),
            analysis_complete = as.character(fixture$da$analysis_complete %||% NA),
            current_s4_object = sha16(fixture$da$current_s4_object),
            da_analysis_results_list = sha16(fixture$workflow$da_analysis_results_list),
            tab_status = sha16(fixture$workflow$tab_status), da_ui_params = sha16(fixture$workflow$da_ui_params)),
        error = if (is.null(run$error)) NA_character_ else savedDaW1Mask(run$error, root),
        calls = run$calls)
}

savedDaW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT112_W1_DUMP")
    if (!nzchar(path)) NULL else if (identical(path, "1")) "" else path
}

#' Compare one case with its frozen value; under the dump switch, append it as R code first.
savedDaW1Expect <- function(case, capture) {
    path <- savedDaW1DumpPath()
    if (!is.null(path)) {
        fields <- vapply(names(capture), function(field) paste0("        ", field, " = ", paste(deparse(capture[[field]],
            width.cutoff = 100L, control = c("keepNA", "keepInteger", "niceNames", "showAttributes")),
            collapse = "\n            ")), "")
        cat(sprintf("    %s = list(\n%s),\n", case, paste(fields, collapse = ",\n")), file = path, append = TRUE)
    }
    frozen <- savedDaW1Frozen()[[case]]
    expect_false(is.null(frozen), info = case)
    for (field in union(names(frozen), names(capture))) {
        expect_identical(capture[[field]], frozen[[field]], info = paste(case, field))
    }
}

#' A DIA (`iq` route) owner on the "small" recipe's protein state, accepted at `correlation_filtered`, with the paths
#' the DA publication writes to, the tab's `da_data` and a plain-environment workflow. The owner is closed at the end.
savedDaW1Owner <- function(project = "kit112-w1", envir = parent.frame()) {
    root <- normalizePath(withr::local_tempdir(.local_envir = envir))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da"),
        publication_graphs_dir = file.path(root, "graphs"))
    for (path in paths) dir.create(path)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, project, "DIA"), WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA")
    manager$.context$experiment_paths <- paths
    manager$saveState("correlation_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        list(), "fixture")
    workflow <- list2env(list(state_manager = manager, tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(manager)
    list(root = root, project = project, paths = paths, workflow = workflow, da = da)
}

#' A DIA_LIMPA owner: KIT-111 W1's builder (test-prot-norm-core.R, the "small" recipe quantified by the limpa DPC
#' engine, its between-sample normalisation and peptide count filter), then the tab's pass-through, correlation (0.5)
#' and RUV covariates in the covariates mode through the in-process runners, so the DA source is `ruv_covariates`.
savedDaW1LimpaOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    quiet <- function(code) suppressWarnings(suppressMessages(utils::capture.output(code)))
    root <- normalizePath(withr::local_tempdir(.local_envir = envir))
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"),
        da_output_dir = file.path(root, "results", "proteomics", "da_proteins"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "kit112-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list(),
        tab_status = list()))
    withr::defer(closeDIAStateAdapter(workflow$state_manager), envir = envir)
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    quantified <- suppressWarnings(suppressMessages(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(
        calcPeptideMatrix(rolled)), dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa")))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    quiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
    }))
    norm <- new.env(parent = emptyenv())
    settings <- list(norm_method = "quantile", ruv_mode = "covariates", auto_percentage_min = NULL,
        auto_percentage_max = NULL, separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL,
        ruv_percentage = 50, ruv_k = 1, color_variable = "group", shape_variable = NULL)
    none <- function(reference) invisible(NULL)
    quiet(runDIAReferenceLimpaNormalisationPassThrough(settings, workflow, norm, none))
    quiet(runDIAReferenceCorrelation(workflow, norm, paths, 0.5, FALSE, "proteomics", "core", function() "group",
        function() list(color_var = "group", shape_var = NULL)))
    quiet(runDIAReferenceRuvCovariates(list(ruv_mode = "covariates", color_variable = "group", shape_variable = NULL),
        workflow, norm, paths, "proteomics", "core", none, function() "group"))
    stopifnot(identical(manager$current_state, "ruv_covariates"))
    da <- new.env(parent = emptyenv())
    da$current_s4_object <- getDIAStateReference(manager)
    list(root = root, project = "kit112-limpa", paths = paths, workflow = workflow, da = da)
}

savedDaW1Contrasts <- function(n = 1L) {
    data.frame(contrasts = c("treated.vs.control", "control.vs.treated")[seq_len(n)],
        full_format = c("treated.vs.control=grouptreated-groupcontrol", "control.vs.treated=groupcontrol-grouptreated")[seq_len(n)],
        friendly_names = c("treated.vs.control", "control.vs.treated")[seq_len(n)], stringsAsFactors = FALSE)
}

#' The tab's DA click on one owner, as the run plugin calls the runner in process (`publish = TRUE`); `fail` forces the
#' kernel to fail in process.
savedDaW1Click <- function(fixture, contrasts, formula = "~ 0 + group", qvalue = 0.05, lfc = 0, fail = FALSE) {
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    settings <- list(formula_string = formula, da_q_val_thresh = qvalue, treat_lfc_cutoff = lfc)
    click <- function() runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts,
        fixture$workflow$uniprot_dat_cln)
    run <- if (fail) savedDaW1Run(testthat::with_mocked_bindings(click(),
        differentialAbundanceAnalysis = function(...) stop("forced DA kernel failure (KIT-112 W1)"), .package = "MultiScholaR")) else
        savedDaW1Run(click())
    savedDaW1Capture(fixture, run)
}

# The frozen values: two dumps on the start commit (64e7c2ae) in separate processes, byte-identical (cmp; kit-112/w1/test-first.txt),
# pasted from the first.
savedDaW1Frozen <- function() list(
    dia_one_contrast = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:cd10d5254dc9df4b", "dia_contrast_results:4e64d874e252e4e5"
            ),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 1L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 8L, scientific_product = 3L, "scientific_product:contrast_references" = 1L, 
            "scientific_product:full_results" = 1L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b"),
        manifest = "0b61131aa2fafaae2f66773198bda5e089d7f9b4d696b269e4ed35bb66077eef",
        model_unset = character(0),
        head_manifest = "0b61131aa2fafaae2f66773198bda5e089d7f9b4d696b269e4ed35bb66077eef",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", "manifest <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "830e51cbaa1a574f079f5de593d9cd74d20625dc2e4dffb7af58d8de2fdf76c6",
        files = "eb6592a02532ef8b8be372cb584ad6bb8daa96c89b27bcacb1e618f175b5faab",
        publication = TRUE,
        display = c(da_results_list = "0b61131aa2fafaae", analysis_complete = "TRUE", current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "0b61131aa2fafaae", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIAStateValue = 10L, publishDIAActionOutputs = 1L
            )),
    dia_two_contrasts = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:7a1252e3b85240ac", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:a813da88d8cd0a89"),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 2L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 12L, scientific_product = 6L, "scientific_product:contrast_references" = 2L, 
            "scientific_product:full_results" = 2L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "a46b71d99173a59365f68f2cd1f7a323650072bb6e684204e2a28be1d61d88fd",
        model_unset = character(0),
        head_manifest = "a46b71d99173a59365f68f2cd1f7a323650072bb6e684204e2a28be1d61d88fd",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", 
            "manifest <- a341788b6cc6b281", "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "0a3492d186cc1aafb17b836c6748ca23a6faf441c19fd9a90fb95150ee0488c0",
        files = "68ed9d605f1e3214d0b7ee1de1bae6a86b15c67de8492452b450b90758534ec6",
        publication = TRUE,
        display = c(da_results_list = "a46b71d99173a593", analysis_complete = "TRUE", current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "a46b71d99173a593", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 0L, saveDIAStateValue = 7L, publishDIAActionOutputs = 1L
            )),
    dia_repeat = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:7a1252e3b85240ac", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:a813da88d8cd0a89"),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 3L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 12L, scientific_product = 6L, "scientific_product:contrast_references" = 2L, 
            "scientific_product:full_results" = 2L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "a46b71d99173a59365f68f2cd1f7a323650072bb6e684204e2a28be1d61d88fd",
        model_unset = character(0),
        head_manifest = "a46b71d99173a59365f68f2cd1f7a323650072bb6e684204e2a28be1d61d88fd",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", 
            "manifest <- a341788b6cc6b281", "other <- a341788b6cc6b281", "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "0a3492d186cc1aafb17b836c6748ca23a6faf441c19fd9a90fb95150ee0488c0",
        files = "68ed9d605f1e3214d0b7ee1de1bae6a86b15c67de8492452b450b90758534ec6",
        publication = TRUE,
        display = c(da_results_list = "a46b71d99173a593", analysis_complete = "TRUE", current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "a46b71d99173a593", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 0L, saveDIAStateValue = 4L, publishDIAActionOutputs = 1L
            )),
    dia_reopened_rerun = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:7a1252e3b85240ac", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:a813da88d8cd0a89"),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 4L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 12L, scientific_product = 6L, "scientific_product:contrast_references" = 2L, 
            "scientific_product:full_results" = 2L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "71d301329ba75b068706d07c5c9014b9c67171fa4bd7a0b1b5976056c29da33f",
        model_unset = c("group_id", "technical_replicate_id", "block_reason", "uses_duplicate_correlation", "consensus_correlation"
            ),
        head_manifest = "71d301329ba75b068706d07c5c9014b9c67171fa4bd7a0b1b5976056c29da33f",
        head_model_unset = c("group_id", "technical_replicate_id", "block_reason", "uses_duplicate_correlation", "consensus_correlation"
            ),
        dependencies = c("dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", 
            "manifest <- a341788b6cc6b281", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", "other <- ebf0e7b372c931b9"
            ),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "af820274ca19065c235d79240275bee252bbca52f2e934ed90320dbb97391dee",
        files = "68ed9d605f1e3214d0b7ee1de1bae6a86b15c67de8492452b450b90758534ec6",
        publication = TRUE,
        display = c(da_results_list = "71d301329ba75b06", analysis_complete = "TRUE", current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "71d301329ba75b06", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 0L, saveDIAStateValue = 4L, publishDIAActionOutputs = 1L
            )),
    dia_non_head = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:7a1252e3b85240ac", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:a813da88d8cd0a89"),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 5L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 12L, scientific_product = 6L, "scientific_product:contrast_references" = 2L, 
            "scientific_product:full_results" = 2L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "e4491943192ae29450909abf2d0c198b7106550d24aaa609670e7ed0abf06e15",
        model_unset = c("group_id", "technical_replicate_id", "block_reason", "uses_duplicate_correlation", "consensus_correlation"
            ),
        head_manifest = "e4491943192ae29450909abf2d0c198b7106550d24aaa609670e7ed0abf06e15",
        head_model_unset = c("group_id", "technical_replicate_id", "block_reason", "uses_duplicate_correlation", "consensus_correlation"
            ),
        dependencies = c("dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", 
            "manifest <- 8d51ea440481d031", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", 
            "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "f53317e704c4ec41d7adfb19c969e4adff7555cb581320a473b9cbb666b4f6ba",
        files = "68ed9d605f1e3214d0b7ee1de1bae6a86b15c67de8492452b450b90758534ec6",
        publication = TRUE,
        display = c(da_results_list = "e4491943192ae294", analysis_complete = "TRUE", current_s4_object = "6d9fd5ca0fe42b72", 
            da_analysis_results_list = "e4491943192ae294", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 0L, saveDIAStateValue = 4L, publishDIAActionOutputs = 1L
            )),
    dia_changed_arguments = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "c198349b104bcdd3e2506f9457fd99b9edd161bd9db982506219495a70b57698"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:7a1252e3b85240ac", "dia-contrast-references:9b1b0675a159d519", 
            "dia-contrast-references:cd10d5254dc9df4b", "dia-contrast-references:d330e022a5dd9434", "dia_contrast_results:4e64d874e252e4e5", 
            "dia_contrast_results:56c76221d51a948e", "dia_contrast_results:a813da88d8cd0a89", "dia_contrast_results:c8fa68a5d5e3ac68", 
            "rendered_asset", "rendered_asset"),
        entities = c("checkpoint:correlation_filtered" = 3L, "checkpoint:da-result-references" = 6L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, rendered_asset = 2L, result_table = 20L, scientific_product = 10L, "scientific_product:contrast_references" = 4L, 
            "scientific_product:full_results" = 4L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b", 
            "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "51141afd59bb292e292fa99ce78a6dcef0f4a9faab00d35bbec89a8a5a77f768",
        model_unset = character(0),
        head_manifest = "51141afd59bb292e292fa99ce78a6dcef0f4a9faab00d35bbec89a8a5a77f768",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", "dia-contrast-references:9b1b0675a159d519 <- 64c4c3b85efb5dd0", 
            "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", "dia-contrast-references:d330e022a5dd9434 <- d564e3af382febac", 
            "manifest <- b3569d64b04c2781", "other <- 8d51ea440481d031", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", 
            "other <- a341788b6cc6b281", "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "f740ed077607f0c239dde97d94d46835cafbb9cefca664e4ed3a11380f570d4d",
        files = "f5079a9e29f77a2489e71865c757fc3442631e095fa3ef321f6b2e0624531c20",
        publication = TRUE,
        display = c(da_results_list = "51141afd59bb292e", analysis_complete = "TRUE", current_s4_object = "78960fb694e3791b", 
            da_analysis_results_list = "51141afd59bb292e", tab_status = "2bce5494a485e1c6", da_ui_params = "a2e1ddc241cb1d1c"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIAStateValue = 14L, publishDIAActionOutputs = 1L
            )),
    dia_qvalue_only = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "a89f68852af9d70828456e6f74bde276dce671d2d186af8b64f237089bf5c8b2"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:6f0fd77cbec27fab", "dia-contrast-references:7a1252e3b85240ac", 
            "dia-contrast-references:8b88692e19d49482", "dia-contrast-references:9b1b0675a159d519", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia-contrast-references:d330e022a5dd9434", "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:56c76221d51a948e", 
            "dia_contrast_results:a813da88d8cd0a89", "dia_contrast_results:c8fa68a5d5e3ac68", "rendered_asset", "rendered_asset", 
            "rendered_asset", "rendered_asset"),
        entities = c("checkpoint:correlation_filtered" = 4L, "checkpoint:da-result-references" = 7L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, rendered_asset = 4L, result_table = 22L, scientific_product = 14L, "scientific_product:contrast_references" = 6L, 
            "scientific_product:full_results" = 4L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b", 
            "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = "93645907ca87d8542877a1ffc33eaf2f8d101f2e16a2d9eab393b2c356393ad9",
        model_unset = character(0),
        head_manifest = "93645907ca87d8542877a1ffc33eaf2f8d101f2e16a2d9eab393b2c356393ad9",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:6f0fd77cbec27fab <- 01147f4a3b2e0a1b", "dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", 
            "dia-contrast-references:8b88692e19d49482 <- 24821eb3c0104091", "dia-contrast-references:9b1b0675a159d519 <- 64c4c3b85efb5dd0", 
            "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", "dia-contrast-references:d330e022a5dd9434 <- d564e3af382febac", 
            "manifest <- f488c61c4db42822", "other <- 8d51ea440481d031", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", 
            "other <- a341788b6cc6b281", "other <- b3569d64b04c2781", "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "a9ec539f90c1241cacb7c8c29e9e55441b7cf5d66bed968eb1e3bd8be34036f9",
        files = "f5079a9e29f77a2489e71865c757fc3442631e095fa3ef321f6b2e0624531c20",
        publication = TRUE,
        display = c(da_results_list = "93645907ca87d854", analysis_complete = "TRUE", current_s4_object = "718d0ef3e308dbcd", 
            da_analysis_results_list = "93645907ca87d854", tab_status = "2bce5494a485e1c6", da_ui_params = "55dee18b660fca71"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIAStateValue = 12L, publishDIAActionOutputs = 1L
            )),
    dia_refused = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "a89f68852af9d70828456e6f74bde276dce671d2d186af8b64f237089bf5c8b2"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:6f0fd77cbec27fab", "dia-contrast-references:7a1252e3b85240ac", 
            "dia-contrast-references:8b88692e19d49482", "dia-contrast-references:9b1b0675a159d519", "dia-contrast-references:cd10d5254dc9df4b", 
            "dia-contrast-references:d330e022a5dd9434", "dia_contrast_results:4e64d874e252e4e5", "dia_contrast_results:56c76221d51a948e", 
            "dia_contrast_results:a813da88d8cd0a89", "dia_contrast_results:c8fa68a5d5e3ac68", "rendered_asset", "rendered_asset", 
            "rendered_asset", "rendered_asset"),
        entities = c("checkpoint:correlation_filtered" = 4L, "checkpoint:da-result-references" = 7L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, rendered_asset = 4L, result_table = 22L, scientific_product = 14L, "scientific_product:contrast_references" = 6L, 
            "scientific_product:full_results" = 4L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b", 
            "scientific_product:full_results:e48e9537858ead4b", "scientific_product:full_results:e48e9537858ead4b"
            ),
        manifest = NA_character_,
        model_unset = NA_character_,
        head_manifest = "93645907ca87d8542877a1ffc33eaf2f8d101f2e16a2d9eab393b2c356393ad9",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:6f0fd77cbec27fab <- 01147f4a3b2e0a1b", "dia-contrast-references:7a1252e3b85240ac <- a3714efd5b0f4d82", 
            "dia-contrast-references:8b88692e19d49482 <- 24821eb3c0104091", "dia-contrast-references:9b1b0675a159d519 <- 64c4c3b85efb5dd0", 
            "dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", "dia-contrast-references:d330e022a5dd9434 <- d564e3af382febac", 
            "manifest <- f488c61c4db42822", "other <- 8d51ea440481d031", "other <- a341788b6cc6b281", "other <- a341788b6cc6b281", 
            "other <- a341788b6cc6b281", "other <- b3569d64b04c2781", "other <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "a9ec539f90c1241cacb7c8c29e9e55441b7cf5d66bed968eb1e3bd8be34036f9",
        files = "f5079a9e29f77a2489e71865c757fc3442631e095fa3ef321f6b2e0624531c20",
        publication = NA,
        display = c(da_results_list = "93645907ca87d854", analysis_complete = "TRUE", current_s4_object = "718d0ef3e308dbcd", 
            da_analysis_results_list = "93645907ca87d854", tab_status = "2bce5494a485e1c6", da_ui_params = "55dee18b660fca71"
            ),
        error = "Invalid DA formula: object 'absent_factor' not found",
        calls = c(runSavedOperation = 0L, commitDIAReferenceState = 0L, saveDIAStateValue = 0L, publishDIAActionOutputs = 0L
            )),
    dia_failure = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = NULL),
        products = character(0),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:initial" = 1L, result_table = 3L),
        fits = character(0),
        manifest = NA_character_,
        model_unset = NA_character_,
        head_manifest = NA_character_,
        head_model_unset = NA_character_,
        dependencies = character(0),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "8f9426728496b00bd37ddc156ba8bb201a18d60efb3bd5611d4b79c39e27ae8c",
        record = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        files = "8ab794fecc83e91d92f2fee89cd8fdc808bd44a3e8c693472e90be9cbf830156",
        publication = NA,
        display = c(da_results_list = "71557d1c8bac9bbe", analysis_complete = NA, current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "71557d1c8bac9bbe", tab_status = "28df520b992b4a4b", da_ui_params = "7fba41c2c7de2596"
            ),
        error = "forced DA kernel failure (KIT-112 W1)",
        calls = c(runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIAStateValue = 2L, publishDIAActionOutputs = 0L
            )),
    dia_retry = list(
        current = "correlation_filtered",
        history = c("initial", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "correlation_filtered"),
        states = c(correlation_filtered = "f8014ed48f2a15fd9fd697dff280fc327fb438466755a885d196eefaceb5d88e"),
        outputs = list(correlation_filtered = c(da = TRUE)),
        products = c("da_fit:5e1ff6c465b978f5", "dia-contrast-references:cd10d5254dc9df4b", "dia_contrast_results:4e64d874e252e4e5"
            ),
        entities = c("checkpoint:correlation_filtered" = 2L, "checkpoint:da-result-references" = 1L, "checkpoint:initial" = 1L, 
            "da_fit:base_fit" = 1L, result_table = 8L, scientific_product = 3L, "scientific_product:contrast_references" = 1L, 
            "scientific_product:full_results" = 1L),
        fits = c("da_fit:base_fit:37aff968eddb0fe1", "scientific_product:full_results:e48e9537858ead4b"),
        manifest = "0b61131aa2fafaae2f66773198bda5e089d7f9b4d696b269e4ed35bb66077eef",
        model_unset = character(0),
        head_manifest = "0b61131aa2fafaae2f66773198bda5e089d7f9b4d696b269e4ed35bb66077eef",
        head_model_unset = character(0),
        dependencies = c("dia-contrast-references:cd10d5254dc9df4b <- ebf0e7b372c931b9", "manifest <- ebf0e7b372c931b9"),
        provenance = "25fda8778949dd0cf7a79409fbc976e93c89159d1a72b0f1c16813b9135237d1",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "830e51cbaa1a574f079f5de593d9cd74d20625dc2e4dffb7af58d8de2fdf76c6",
        files = "eb6592a02532ef8b8be372cb584ad6bb8daa96c89b27bcacb1e618f175b5faab",
        publication = TRUE,
        display = c(da_results_list = "0b61131aa2fafaae", analysis_complete = "TRUE", current_s4_object = "f41d17b919c2ce82", 
            da_analysis_results_list = "0b61131aa2fafaae", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 1L, commitDIAReferenceState = 0L, saveDIAStateValue = 8L, publishDIAActionOutputs = 1L
            )),
    limpa_one_contrast = list(
        current = "ruv_covariates",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        head_current = "ruv_covariates",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        states = c(protein_dpc_quant = "22f553887dc4edb1623ab2d0b556c048cd7f216ba22d1f250abd9d86eb5db726", protein_normalised = "8f9c8ef12ddf5740a442f8b23bf258aeb7a02793932943a0f1901d84a9e9867b", 
            peptide_count_filtered = "32b145aaa12d358ffff2ff6d857c643e9305792be84d2aeec59c0bb9c9fab661", correlation_filtered = "69b20b13edbf2f85ae088f8fa1e296b1412d2a1d7ea8021a85307b983e4b45fb", 
            ruv_covariates = "7ba49b7725df7e803f6e4587bda299c59c8002911538d7f6e7e0674640b641fa"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE), 
                ruv_covariates = c(da = TRUE)),
        products = c("dia-contrast-references:97aa6e4c96980b17", "dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:cb46f62c60963c21", 
            "dia-ruv-covariates:0055cc0a9db60428", "dia_normalization_pca:d9b5720af8d8a298", "rendered_asset"),
        entities = c("checkpoint:correlation_filtered" = 1L, "checkpoint:da-result-references" = 1L, "checkpoint:initial" = 1L, 
            "checkpoint:owned_args:limpa_dpc_quant_results" = 2L, "checkpoint:owner:pearson-filtered-science" = 1L, 
            "checkpoint:owner:ruv-adjusted-science" = 1L, "checkpoint:owner:ruv-covariates-science" = 1L, "checkpoint:peptide_count_filtered" = 1L, 
            "checkpoint:protein_dpc_quant" = 1L, "checkpoint:protein_normalised" = 1L, "checkpoint:ruv_covariates" = 2L, 
            rendered_asset = 1L, result_table = 14L, scientific_product = 10L, "scientific_product:contrast_references" = 1L
            ),
        fits = character(0),
        manifest = "dd86ce37c0f1928eb0946caa82b8643f6344f5959d774f844b4dd68f762885cd",
        model_unset = c("group_id", "technical_replicate_id", "consensus_correlation"),
        head_manifest = "dd86ce37c0f1928eb0946caa82b8643f6344f5959d774f844b4dd68f762885cd",
        head_model_unset = c("group_id", "technical_replicate_id", "consensus_correlation"),
        dependencies = c("dia-contrast-references:97aa6e4c96980b17 <- 52bb3bd82853a6b0", "manifest <- 52bb3bd82853a6b0"),
        provenance = "b277eb8197ce4789839d312cf291fed25eae7f1c6b5450795b7b6706b9c33067",
        pending = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        record = "e2c78f55bc3748ad827e2ac8cb66d8414304e4c035940dcc2bb41f98a548ebe8",
        files = "c983c2965d1b14c7ab12169b6d2a50dad95db689680af6bc3a5d439e1562f0b1",
        publication = TRUE,
        display = c(da_results_list = "dd86ce37c0f1928e", analysis_complete = "TRUE", current_s4_object = "5d94fbb7979cd1eb", 
            da_analysis_results_list = "dd86ce37c0f1928e", tab_status = "2bce5494a485e1c6", da_ui_params = "7fba41c2c7de2596"
            ),
        error = NA_character_,
        calls = c(runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIAStateValue = 9L, publishDIAActionOutputs = 1L
            ))
)

test_that("KIT-112 W1 (a): the DIA DA clicks in process give the recorded states, products, fits, manifest, record, files and display", {
    skip_if_not_installed("Glimma")
    fixture <- savedDaW1Owner()
    first <- fixture$da$current_s4_object
    savedDaW1Expect("dia_one_contrast", savedDaW1Click(fixture, savedDaW1Contrasts(1L)))
    # The first contrast's product is reused; the second computes on the saved base fit.
    savedDaW1Expect("dia_two_contrasts", savedDaW1Click(fixture, savedDaW1Contrasts(2L)))
    # The same click again: every contrast a hit, no new scientific entity.
    savedDaW1Expect("dia_repeat", savedDaW1Click(fixture, savedDaW1Contrasts(2L)))
    # A fully cached rerun after the owner is reopened from its store: the base-fit reference is unset. (At q 0.05
    # without TREAT this fixture publishes no Glimma page, so no saved page is met on the rerun; W1 finding F112-03.)
    closeDIAStateAdapter(fixture$workflow$state_manager)
    reopened <- newDIAStateAdapter(newWorkflowEntityStore(fixture$root, fixture$project, "DIA"),
        WorkflowState$new(audit_enabled = FALSE), reopen = TRUE)
    reopened$.context$experiment_paths <- fixture$paths
    fixture$workflow$state_manager <- reopened
    fixture$da$current_s4_object <- getDIAStateReference(reopened)
    savedDaW1Expect("dia_reopened_rerun", savedDaW1Click(fixture, savedDaW1Contrasts(2L)))
    # A non-head source: the tab on the state before the arguments' commit, the arguments unchanged.
    head <- fixture$da$current_s4_object
    fixture$da$current_s4_object <- first
    savedDaW1Expect("dia_non_head", savedDaW1Click(fixture, savedDaW1Contrasts(2L)))
    fixture$da$current_s4_object <- head
    # Changed DA arguments (TREAT at 0.5): the arguments commit, then the run.
    savedDaW1Expect("dia_changed_arguments", savedDaW1Click(fixture, savedDaW1Contrasts(2L), lfc = 0.5))
    # A q-value change only: new contrast products on the same base fit.
    savedDaW1Expect("dia_qvalue_only", savedDaW1Click(fixture, savedDaW1Contrasts(2L), qvalue = 0.01, lfc = 0.5))
    # A refusal by the controller's validation: nothing committed, nothing run.
    savedDaW1Expect("dia_refused", savedDaW1Click(fixture, savedDaW1Contrasts(2L), formula = "~ 0 + absent_factor",
        qvalue = 0.01, lfc = 0.5))
})

test_that("KIT-112 W1 (a): a DA kernel failure in process, then a retry with the same arguments, give the recorded values", {
    skip_if_not_installed("Glimma")
    fixture <- savedDaW1Owner("kit112-w1-retry")
    savedDaW1Expect("dia_failure", savedDaW1Click(fixture, savedDaW1Contrasts(1L), fail = TRUE))
    savedDaW1Expect("dia_retry", savedDaW1Click(fixture, savedDaW1Contrasts(1L)))
})

test_that("KIT-112 W1 (a): the DIA_LIMPA DA click (dpcDE, no base fit) in process gives the recorded values", {
    skip_if_not_installed("Glimma")
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    fixture <- savedDaW1LimpaOwner()
    savedDaW1Expect("limpa_one_contrast", savedDaW1Click(fixture, savedDaW1Contrasts(1L)))
})

# ------------------------------------------------------------------------------------------------
# KIT-112 W2 (W1 (b)): the split runner on the managed worker, at core level. Each W1 (a) click runs here as the managed
# tab will run it: the request half in the controller, one `DA-contrast-output` action per contrast through the
# worker's entry point (`workflowWorkerMain`) on the request's named input, the base fit adopted when delivered, the
# accept half in the controller, then the `da_publication` step through the worker and the controller's publication,
# and the step's delivery (`readManagedDaResult`, `completeDIADaDisplay`). The DIA clicks are also run on a twin owner
# through today's managed path (the whole runner in the controller without its publication, then the same publication
# step), and the two managed paths must be equal in every captured field but the `saveDIAStateValue` count (it counts
# the head snapshots, which the worker's completions write; D6). Against W1 (a)'s frozen in-process values the worker
# path must be equal too, except, by design: that count; the pending step provenance, which on both managed paths is the
# `da_publication` completion's (KIT-019, unchanged; in process the accept half leaves none) and, after a worker
# failure, is held by the worker's owner (D6); the committed step provenance of a DA-arguments commit that follows a
# managed publication, which records that pending entry on both managed paths (draft finding, W2 notes); and in the
# non-head case the refused publication and what it leaves (R7: today's managed outcome, the publication refuses a
# non-head manifest).
# ------------------------------------------------------------------------------------------------

#' One managed action in this process: the request (with its named input) through the worker's entry point, the
#' acceptance and what the tab reads. A worker failure is raised with the worker's message. Each delivered value's size
#' and names are kept in `log`.
savedDaW2Dispatch <- function(fixture, controller, operation, parameters, input, log = NULL) {
    store <- fixture$workflow$state_manager$.context$store
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
        generation = controller$control$generation, controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
        list(experiment_label = "core", experiment_paths = fixture$paths, omic_type = "proteomics"))
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
        expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
        input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
        project_root = dirname(store$directory), directory = paths$request_dir, config = list(mode = "adapters"),
        worker_cap = 12 * 1024^3)
    path <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    if (is.environment(log)) log$values[[length(log$values) + 1L]] <- list(operation = operation,
        bytes = as.numeric(object.size(message$result)), names = names(message$result), status = message$result$status,
        base_fit = !is.null(message$result$display$base_fit), mutable = hasMutableWorkflowEntityOwner(message$result))
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, NULL))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(fixture$workflow$state_manager, accepted, NULL), list(acceptance = receipt$status, accepted = accepted))
}

#' The managed DA click's front: the request half in the controller, each contrast on the worker (the route in the
#' request), the base fit adopted only when delivered (R1), the accept half in the controller. `today = TRUE` is
#' today's managed front instead: the whole runner in the controller without its publication.
savedDaW2Front <- function(fixture, dispatch, settings, contrasts, today = FALSE) {
    manager <- fixture$workflow$state_manager
    annotations <- fixture$workflow$uniprot_dat_cln
    if (today) return(runDIAReferenceDaAnalysis(settings, fixture$da, fixture$workflow, fixture$paths, contrasts, annotations,
        publish = FALSE))
    request <- withDIAStateScope(manager, function() diaDaRequest(manager, settings, fixture$da, fixture$workflow, contrasts,
        annotations))
    route <- withDIAStateScope(manager, function() diaActiveRoute(manager))
    results <- list()
    for (step in request$steps) {
        delivered <- dispatch("DA-contrast-output", c(step, list(route = route)), request$source)
        results[[step$name]] <- delivered$value
        if (!is.null(delivered$display$base_fit)) manager$.context$da_base_fit_reference <- delivered$display$base_fit
    }
    withDIAStateScope(manager, function() acceptDIADaManifest(manager, request$manifest, results))
}

#' One managed DA click: the front, then the publication step on the worker and its delivery to the tab.
savedDaW2Managed <- function(fixture, dispatch, settings, contrasts, today = FALSE) {
    manager <- fixture$workflow$state_manager
    manifest <- savedDaW2Front(fixture, dispatch, settings, contrasts, today)
    published <- dispatch("da_publication", list(annotations = NULL), getDIAStateReference(manager))
    delivered <- withDIAStateScope(manager, function() readManagedDaResult(manager, published$accepted, NULL))
    completeDIADaDisplay(fixture$da, fixture$workflow, delivered$manifest)
    manifest
}

#' The tab's DA click on the managed path, captured as `savedDaW1Click` captures the in-process one.
savedDaW2Click <- function(fixture, dispatch, contrasts, formula = "~ 0 + group", qvalue = 0.05, lfc = 0, fail = FALSE,
    today = FALSE) {
    fixture$workflow$contrasts_tbl <- contrasts
    fixture$da$contrasts_available <- contrasts$contrasts
    settings <- list(formula_string = formula, da_q_val_thresh = qvalue, treat_lfc_cutoff = lfc)
    click <- function() savedDaW2Managed(fixture, dispatch, settings, contrasts, today)
    run <- if (fail) savedDaW1Run(testthat::with_mocked_bindings(click(),
        differentialAbundanceAnalysis = function(...) stop("forced DA kernel failure (KIT-112 W1)"), .package = "MultiScholaR")) else
        savedDaW1Run(click())
    savedDaW1Capture(fixture, run)
}

#' A dispatcher on the fixture's current owner, with a controller per owner (a reopen gets its own).
savedDaW2Dispatcher <- function(fixture, log) {
    owner <- NULL
    controller <- NULL
    function(operation, parameters, input) {
        manager <- fixture$workflow$state_manager
        if (!identical(manager, owner)) {
            owner <<- manager
            controller <<- omicsWorkerController(list(owner = manager))
        }
        savedDaW2Dispatch(fixture, controller, operation, parameters, input, log)
    }
}

#' Compare one managed case with W1 (a)'s frozen in-process value: every field, with `saveDIAStateValue`'s count and
#' the fields named in `masked` excepted (D6).
savedDaW2Expect <- function(case, capture, masked = character(), frozen = savedDaW1Frozen()[[case]],
    masked_calls = character()) {
    expect_false(is.null(frozen), info = case)
    for (field in setdiff(union(names(frozen), names(capture)), c("calls", masked))) {
        expect_identical(capture[[field]], frozen[[field]], info = paste("worker", case, field))
    }
    keep <- setdiff(names(frozen$calls), c("saveDIAStateValue", masked_calls))
    expect_identical(capture$calls[keep], frozen$calls[keep], info = paste("worker", case, "calls"))
}

#' The operation of the owner's pending step provenance (NA when none is pending).
savedDaW2Pending <- function(fixture) fixture$workflow$state_manager$.context$pending_provenance$operation %||% NA_character_

savedDaW2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' Reopen a fixture's owner from its store, as W1 (a)'s reopened case does.
savedDaW2Reopen <- function(fixture) {
    closeDIAStateAdapter(fixture$workflow$state_manager)
    reopened <- newDIAStateAdapter(newWorkflowEntityStore(fixture$root, fixture$project, "DIA"),
        WorkflowState$new(audit_enabled = FALSE), reopen = TRUE)
    reopened$.context$experiment_paths <- fixture$paths
    fixture$workflow$state_manager <- reopened
    fixture$da$current_s4_object <- getDIAStateReference(reopened)
    invisible(reopened)
}

test_that("KIT-112 W2 (W1 (b)): the DIA DA clicks through the worker equal today's managed path and W1 (a)'s values", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    savedDaW2Mocks()
    log <- new.env(parent = emptyenv())
    fixture <- savedDaW1Owner("kit112-w2")
    twin <- savedDaW1Owner("kit112-w2-today")
    dispatch <- savedDaW2Dispatcher(fixture, log)
    twin_dispatch <- savedDaW2Dispatcher(twin, NULL)
    # One click on both managed paths: the worker against today's, then against W1 (a)'s frozen in-process value.
    both <- function(case, n, ..., masked = character(), masked_calls = character()) {
        worker <- savedDaW2Click(fixture, dispatch, savedDaW1Contrasts(n), ...)
        today <- savedDaW2Click(twin, twin_dispatch, savedDaW1Contrasts(n), ..., today = TRUE)
        savedDaW2Expect(case, worker, frozen = today)
        savedDaW2Expect(case, worker, masked = c("pending", "provenance", masked), masked_calls = masked_calls)
        worker
    }
    first <- list(worker = fixture$da$current_s4_object, today = twin$da$current_s4_object)
    both("dia_one_contrast", 1L)
    # The pending provenance on a managed path is the publication step's completion's.
    expect_identical(c(savedDaW2Pending(fixture), savedDaW2Pending(twin)), c("da_publication", "da_publication"))
    both("dia_two_contrasts", 2L)
    both("dia_repeat", 2L)
    # R1: after the reopen every contrast hits, no base fit is delivered, and the model's base-fit fields stay NULL.
    savedDaW2Reopen(fixture)
    savedDaW2Reopen(twin)
    both("dia_reopened_rerun", 2L)
    # R7: a non-head source. Both managed paths compute and accept it, and the publication refuses it.
    head <- list(worker = fixture$da$current_s4_object, today = twin$da$current_s4_object)
    fixture$da$current_s4_object <- first$worker
    twin$da$current_s4_object <- first$today
    worker <- both("dia_non_head", 2L, masked = c("outputs", "files", "publication", "display", "error", "manifest",
        "model_unset"), masked_calls = "publishDIAActionOutputs")
    expect_match(worker$error, "DA manifest does not belong to the head state", fixed = TRUE)
    fixture$da$current_s4_object <- head$worker
    twin$da$current_s4_object <- head$today
    both("dia_changed_arguments", 2L, lfc = 0.5)
    both("dia_qvalue_only", 2L, qvalue = 0.01, lfc = 0.5)
    both("dia_refused", 2L, formula = "~ 0 + absent_factor", qvalue = 0.01, lfc = 0.5)
    # R5: each contrast's value holds references only, under the acceptance's 8 MiB; the base fit is delivered by the
    # contrast that ran the kernel on a base fit, never by a hit.
    contrasts <- Filter(function(entry) identical(entry$operation, "DA-contrast-output"), log$values)
    expect_length(contrasts, 13L)
    for (entry in contrasts) {
        expect_lt(entry$bytes, 8388608)
        expect_false(entry$mutable)
        expect_identical(entry$names, c("checkpoint", "reference", "native", "display", "operation", "roots", "owned"))
    }
    expect_identical(vapply(contrasts, `[[`, NA, "base_fit"),
        c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
})

test_that("KIT-112 W2 (W1 (b)): a worker failure, then a retry with the same arguments, give W1 (a)'s values", {
    skip_if_not_installed("Glimma")
    skip_if_not(workflowOwnerCapable())
    savedDaW2Mocks()
    fixture <- savedDaW1Owner("kit112-w2-retry")
    dispatch <- savedDaW2Dispatcher(fixture, NULL)
    # The worker's owner, not the controller's, held the failed contrast's pending provenance (D6).
    savedDaW2Expect("dia_failure", savedDaW2Click(fixture, dispatch, savedDaW1Contrasts(1L), fail = TRUE), masked = "pending")
    expect_identical(savedDaW2Pending(fixture), NA_character_)
    savedDaW2Expect("dia_retry", savedDaW2Click(fixture, dispatch, savedDaW1Contrasts(1L)), masked = "pending")
    expect_identical(savedDaW2Pending(fixture), "da_publication")
})

test_that("KIT-112 W2 (W1 (b)): the DIA_LIMPA DA click through the worker (dpcDE, no base fit) gives W1 (a)'s values", {
    skip_if_not_installed("Glimma")
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    skip_if_not(workflowOwnerCapable())
    savedDaW2Mocks()
    log <- new.env(parent = emptyenv())
    fixture <- savedDaW1LimpaOwner()
    savedDaW2Expect("limpa_one_contrast", savedDaW2Click(fixture, savedDaW2Dispatcher(fixture, log), savedDaW1Contrasts(1L)),
        masked = "pending")
    expect_identical(savedDaW2Pending(fixture), "da_publication")
    expect_false(log$values[[1L]]$base_fit)
})

test_that("KIT-112 W2: the worker handler strips the route from the request, so the compute half takes the in-process parameters", {
    step <- list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0, weights = NA, contrast = "a=groupB-groupA",
        annotations = "entity", name = "a", comparison = "a")
    seen <- new.env(parent = emptyenv())
    source <- list(entity_id = "source", kind = "checkpoint")
    testthat::local_mocked_bindings(
        computeDIAReferenceDaContrast = function(manager, source, parameters, route) {
            seen$parameters <- parameters
            seen$route <- route
            list(run = list(status = "computed", value = list(), reference = list(entity_id = "product")), base_fit = NULL)
        },
        readWorkflowStateMetadata = function(manager, reference) list(data = list(record = reference$entity_id)),
        .package = "MultiScholaR")
    for (route in c("iq", "limpa")) {
        value <- computeWorkflowDiaContrastAction(NULL, list(input_references = list(state = source),
            parameters = c(step, list(route = route))), "whole")
        expect_identical(seen$parameters, step, info = route)
        expect_identical(seen$route, route)
        expect_identical(value$checkpoint, list(record = "source"))
        expect_identical(value$display, list(base_fit = NULL))
    }
    expect_error(computeWorkflowDiaContrastAction(NULL, list(input_references = list(state = source), parameters = step), "whole"),
        class = "workflow_operation_parameters")
})

test_that("KIT-112 W2 (D4): the DIA adaptor's `da` entry names the registered omics operation, the publication step and existing functions", {
    da <- proteomicsDiaBinding$da
    expect_identical(names(da), c("operation", "in_process", "request", "accept", "publication"))
    entry <- workflowActionSpec(da$operation)
    expect_identical(entry[c("compute", "handler", "completion", "named_input")],
        list(compute = "omics", handler = "computeWorkflowDiaContrastAction", completion = "omics", named_input = TRUE))
    expect_true(is.function(get(entry$handler, envir = asNamespace("MultiScholaR"), mode = "function")))
    expect_identical(workflowActionSpec(da$publication)$compute, "protein")
    for (field in c("in_process", "request", "accept")) expect_true(is.function(da[[field]]), info = field)
    for (name in c("runDIAReferenceDaAnalysis", "diaDaRequest", "acceptDIADaManifest", "computeDIAReferenceDaContrast"))
        expect_true(exists(name, envir = asNamespace("MultiScholaR"), mode = "function", inherits = FALSE), info = name)
    expect_identical(names(formals(runDIAReferenceDaAnalysis)),
        c("settings", "da_data", "workflow", "paths", "contrasts", "annotations", "publish"))
})
