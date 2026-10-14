makeSavedQcProteinObject <- function(constant = FALSE) {
  set.seed(11011L)
  samples <- sprintf("S%02d", seq_len(8L))
  features <- sprintf("PG%03d", seq_len(60L))
  values <- if (constant) {
    matrix(5, nrow = length(features), ncol = length(samples))
  } else {
    matrix(stats::rnorm(length(features) * length(samples), 20, 2),
      nrow = length(features), ncol = length(samples))
  }
  colnames(values) <- samples
  table <- data.frame(Protein.Group = features, values, check.names = FALSE)
  rownames(table) <- NULL
  design <- data.frame(
    Run = samples,
    group = factor(rep(c("control", "treated"), each = 4L)),
    batch = factor(rep(c("one", "two"), times = 4L)),
    tech_rep_group = rep(c("control_1", "control_2", "treated_1", "treated_2"),
      each = 2L),
    row.names = samples
  )
  methods::new("ProteinQuantitativeData",
    protein_quant_table = table,
    protein_id_column = "Protein.Group",
    design_matrix = design,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "tech_rep_group")
}

test_that("automatic RUV failures preserve GUI holders and release transient plot owners", {
  original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
  prior_plot <- ggplot2::get_last_plot()
  withr::defer(ggplot2::set_last_plot(prior_plot))
  root <- withr::local_tempdir()
  shiny::testServer(function(input, output, session) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "ruv-failure-owner", "DIA"),
      WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    workflow <- shiny::reactiveValues(state_manager = manager,
      ruv_optimization_result = list(previous = TRUE))
    attachDIAStateSession(manager, session, workflow)
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    norm <- shiny::reactiveValues(ruv_optimization_result = list(previous = TRUE))
  }, {
    rlang::local_bindings(config_list = original@args, .env = .GlobalEnv)
    settings <- list(norm_method = "quantile", ruv_mode = "automatic",
      auto_percentage_min = 1, auto_percentage_max = 20,
      separation_metric = "max_difference", k_penalty_weight = 0.1,
      adaptive_k_penalty = TRUE, color_variable = "group", shape_variable = NULL)
    reference <- getDIAStateReference(manager)
    run <- function() runDIAReferenceNormalization(settings, workflow, norm,
      list(source_dir = root), "proteomics", "fixture", function(...) NULL,
      function(...) NULL, function() "group")
    for (stage in c("validation", "correction")) {
      ggplot2::set_last_plot(prior_plot)
      reached <- FALSE
      fail <- function(...) { reached <<- TRUE; stop("injected RUV failure") }
      if (stage == "validation") {
        expect_error(testthat::with_mocked_bindings(run(), validateProtNormRuvK = fail),
          "injected RUV failure")
      } else {
        expect_error(testthat::with_mocked_bindings(run(), ruvIII_C_Varying = fail),
          "injected RUV failure")
      }
      expect_true(reached, info = stage)
      expect_true(identical(norm$ruv_optimization_result, list(previous = TRUE)), info = stage)
      expect_true(identical(workflow$ruv_optimization_result, list(previous = TRUE)), info = stage)
      expect_true(identical(ggplot2::get_last_plot(), prior_plot), info = stage)
      expect_identical(getDIAStateReference(manager), reference)
      expect_identical(manager$.context$scope_depth, 0L)
      expect_identical(manager$.context$store$open_resources, 0L)
      expect_length(ls(manager$.context$working), 0L)
      expect_null(manager$.context$reclaim_pending)
    }
  })
})

test_that("reference QC colour changes use saved science without reading an assay", {
  object <- makeSavedQcProteinObject()
  manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "reference-qc", "DIA"),
    WorkflowState$new(audit_enabled = FALSE))
  withr::defer(closeDIAStateAdapter(manager))
  manager$.context$reference_states <- TRUE
  manager$saveState("protein_replicate_filtered", object, list(), "input")
  reference <- getDIAStateReference(manager)
  testthat::with_mocked_bindings({
    first <- resolveDIAReferenceQcPlot(reference, list(color_var = "group", shape_var = NULL), "pca",
      buildProtNormQcPlotProduct, plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot)
    changed <- testthat::with_mocked_bindings(
      resolveDIAReferenceQcPlot(reference, list(color_var = "batch", shape_var = NULL), "pca",
        function(...) stop("PCA was refitted"), plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot),
      materializeDIATable = function(...) stop("assay or design was hydrated on a product hit"))
    expect_identical(first$data[c("Run", "PC1", "PC2")], changed$data[c("Run", "PC1", "PC2")])
    expect_match(rlang::as_label(changed$mapping$colour), "batch", fixed = TRUE)
    expect_identical(attr(first, "dia_product_reference"), attr(changed, "dia_product_reference"))
    expect_identical(manager$.context$full_s4_reads, 0L)
    expect_identical(manager$.context$store$open_resources, 0L)
  }, getDIAActiveManager = function() manager)
})

test_that("QC products preserve populations and density runs its own PCA", {
  object <- makeSavedQcProteinObject()
  aesthetics <- list(color_var = "group", shape_var = "batch")
  pca_calls <- 0L
  pca <- function(...) {
    pca_calls <<- pca_calls + 1L
    plotPca(...)
  }
  density <- function(s4Object, aesthetics) {
    buildProtNormDensityPlot(s4Object, aesthetics, plotPcaFn = pca)
  }
  products <- list(
    pca = buildProtNormQcPlotProduct(object, aesthetics, "pca",
      plotPcaFn = pca, buildDensityFn = density),
    rle = buildProtNormQcPlotProduct(object, aesthetics, "rle"),
    density = buildProtNormQcPlotProduct(object, aesthetics, "density",
      plotPcaFn = pca, buildDensityFn = density),
    correlation = buildProtNormQcPlotProduct(object, aesthetics, "correlation")
  )
  expect_identical(pca_calls, 2L)
  expect_identical(products$pca$population_role, "pca_existing_producer")
  expect_identical(products$density$population_role,
    "density_independent_pca_call")
  expect_identical(products$pca$feature_ids,
    methods::slot(object, "protein_quant_table")$Protein.Group)
  expect_identical(products$pca$sample_ids,
    as.character(methods::slot(object, "design_matrix")$Run))
  expect_identical(rownames(products$pca$annotations), products$pca$sample_ids)
  expect_identical(products$pca$data$Run, products$pca$sample_ids)
  expect_true(all(c("PC1", "PC2", "group", "batch") %in%
    names(products$pca$data)))
  expect_true(all(c("min", "lower", "middle", "upper", "max") %in%
    names(products$rle$data)))
  expect_true("pearson_correlation" %in% names(products$correlation$data))
  expect_identical(products$pca$plot$plot_env, rlang::base_env())
  expect_identical(products$rle$plot$plot_env, rlang::base_env())
  expect_identical(products$correlation$plot$plot_env, rlang::base_env())
  expect_false(any(vapply(products, methods::is, logical(1), "ProteinQuantitativeData")))
})

test_that("saved QC products and exact assets avoid repeat production and rendering", {
  object <- makeSavedQcProteinObject()
  aesthetics <- list(color_var = "group", shape_var = "batch")
  producer_calls <- 0L
  produce <- function() {
    producer_calls <<- producer_calls + 1L
    buildProtNormQcPlotProduct(object, aesthetics, "pca")
  }
  product <- produce()
  store <- newWorkflowEntityStore(withr::local_tempdir(), "qc-project", "DIA",
    max_cache_bytes = 67108864)
  identity <- newWorkflowDiaQcProductIdentity(object, "post_normalization",
    "pca", aesthetics, "MultiScholaR-qc-test-v1")
  product_ref <- saveWorkflowDiaQcProduct(store, product, "post-norm-pca",
    identity, "post_normalization")
  hit <- resolveWorkflowDiaQcProduct(store, identity, "post_normalization", "pca")
  expect_identical(hit$status, "hit")
  expect_identical(hit$value$value$data, product$data)
  expect_identical(producer_calls, 1L)

  style <- list(title = "Post-normalization PCA", font_size = 9)
  render_calls <- 0L
  save_plot <- function(...) {
    render_calls <<- render_calls + 1L
    ggplot2::ggsave(...)
  }
  expect_null(hit$value$value$plot)
  asset <- renderWorkflowDiaQcProduct(hit$value$value, style, "png", c(800L, 600L),
    savePlotFn = save_plot)
  render_identity <- newSavedWorkflowRenderIdentity(product_ref, style, "png",
    c(800L, 600L), "ggplot2-test-v1")
  asset_ref <- saveWorkflowDiaQcAsset(store, asset, "post-norm-pca-png",
    render_identity)
  releaseWorkflowEntityCache(store)
  before <- workflowEntityCacheInfo(store)
  asset_hit <- resolveSavedWorkflowAsset(store, product_ref, render_identity)
  after <- workflowEntityCacheInfo(store)
  expect_identical(asset_hit$value, asset)
  expect_identical(render_calls, 1L)
  expect_identical(after$payload_reads, before$payload_reads + 1L)
  expect_identical(asset_ref$kind, "rendered_asset")
  expect_gt(length(asset_hit$value$bytes), 1000L)

  changed_style <- list(title = "Changed title only", font_size = 9)
  changed_render_identity <- newSavedWorkflowRenderIdentity(product_ref,
    changed_style, "png", c(800L, 600L), "ggplot2-test-v1")
  expect_identical(resolveSavedWorkflowAsset(store, product_ref,
    changed_render_identity)$status, "miss")
  expect_identical(resolveWorkflowDiaQcProduct(store, identity,
    "post_normalization", "pca")$status, "hit")
  expect_identical(producer_calls, 1L)
  expect_identical(render_calls, 1L)

  changed_object <- object
  changed_table <- methods::slot(changed_object, "protein_quant_table")
  changed_table[[2L]][[1L]] <- changed_table[[2L]][[1L]] + 1
  methods::slot(changed_object, "protein_quant_table") <- changed_table
  changed_identity <- newWorkflowDiaQcProductIdentity(changed_object,
    "post_normalization", "pca", aesthetics, "MultiScholaR-qc-test-v1")
  expect_identical(resolveWorkflowDiaQcProduct(store, changed_identity,
    "post_normalization", "pca")$status, "miss")
})

test_that("constant inputs and unavailable panels preserve original behavior", {
  object <- makeSavedQcProteinObject(constant = TRUE)
  aesthetics <- list(color_var = "group", shape_var = "batch")
  baseline <- plotPca(object, grouping_variable = "group", shape_variable = "batch",
    label_column = "", title = "", font_size = 8)
  product <- buildProtNormQcPlotProduct(object, aesthetics, "pca")
  expect_identical(product$data, as.data.frame(baseline$data))
  expect_true(all(product$data$PC1 == 0))
  expect_true(all(product$data$PC2 == 0))

  unavailable <- newSavedWorkflowScientificProduct(NULL,
    "dia_normalization_correlation", provenance = list(route = "iq"),
    status = "unavailable", reason = "no finite within-group correlations")
  store <- newWorkflowEntityStore(withr::local_tempdir(), "qc-unavailable", "DIA")
  identity <- newWorkflowDiaQcProductIdentity(object, "post_normalization",
    "correlation", aesthetics, "MultiScholaR-qc-test-v1")
  reference <- saveWorkflowConsumerEntity(store, unavailable,
    "unavailable-correlation", identity, "scientific_product",
    factors = list(route = "iq", stage = "post_normalization",
      plot_type = "correlation"))
  hit <- resolveWorkflowDiaQcProduct(store, identity, "post_normalization",
    "correlation")
  expect_identical(hit$reference, reference)
  expect_identical(hit$value$status, "unavailable")
  expect_null(hit$value$value)
  expect_match(hit$value$reason, "no finite")
})

test_that("all-missing PCA behavior matches the unchanged producer", {
  object <- makeSavedQcProteinObject()
  table <- methods::slot(object, "protein_quant_table")
  table[-1L] <- lapply(table[-1L], \(value) rep(NA_real_, length(value)))
  methods::slot(object, "protein_quant_table") <- table
  aesthetics <- list(color_var = "group", shape_var = "batch")
  baseline <- tryCatch(
    plotPca(object, grouping_variable = "group", shape_variable = "batch",
      label_column = "", title = "", font_size = 8),
    error = \(condition) condition
  )
  candidate <- tryCatch(
    buildProtNormQcPlotProduct(object, aesthetics, "pca"),
    error = \(condition) condition
  )
  expect_identical(inherits(candidate, "error"), inherits(baseline, "error"))
  if (inherits(baseline, "error")) {
    expect_identical(conditionMessage(candidate), conditionMessage(baseline))
  } else {
    expect_identical(candidate$data, as.data.frame(baseline$data))
  }
})


test_that("QC product audit records independent plotting owners", {
    original <- ggplot2::get_last_plot()
    withr::defer(ggplot2::set_last_plot(original))
    object <- makeSavedQcProteinObject()
    inspect <- function(frame) {
        fields <- if (is.environment(frame) && !identical(frame, baseenv()) &&
            !identical(frame, globalenv()) && !isNamespace(frame))
                intersect(c("theObject", "s4Object", "data", "data_filtered", "pca.res", "working_matrix"),
                    ls(frame, all.names = TRUE)) else character()
        lapply(fields, function(field) {
            value <- get(field, frame, inherits = FALSE)
            list(name = field, class = class(value), shape = dim(value),
                same_scientific_input = identical(value, object))
        })
    }
    inspect_plot <- function(plot) {
        if (is.null(plot)) return(NULL)
        mappings <- lapply(plot$mapping, function(value) if (rlang::is_quosure(value))
            inspect(rlang::quo_get_env(value)) else list())
        labels <- lapply(plot$scales$scales, function(scale) {
            value <- scale$labels
            if (inherits(value, "ggproto_method")) value <- environment(value)$f
            if (is.function(value)) inspect(environment(value)) else list()
        })
        list(plot_env = inspect(plot$plot_env), mappings = mappings, labels = labels,
            children = lapply(plot$patches$plots, inspect_plot))
    }
    observations <- list()
    for (kind in c("pca", "rle", "density", "correlation")) {
        ggplot2::set_last_plot(NULL)
        product <- buildProtNormQcPlotProduct(object,
            list(color_var = "group", shape_var = "batch"), kind)
        global_plot <- ggplot2::get_last_plot()
        expect_null(global_plot)
        returned <- inspect_plot(product$plot)
        retained <- function(value) {
            if (!is.list(value)) return(FALSE)
            if (!is.null(value$same_scientific_input)) return(isTRUE(value$same_scientific_input) ||
                any(value$class %in% c("ProteinQuantitativeData", "matrix", "pca")))
            any(vapply(value, retained, logical(1)))
        }
        expect_false(retained(returned))
        observations[[kind]] <- list(global_plot_present = !is.null(global_plot),
            returned_plot_environment_is_base = identical(product$plot$plot_env, baseenv()),
            global = inspect_plot(global_plot), returned = inspect_plot(product$plot))
        expect_identical(product$plot$plot_env, baseenv())
    }
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    if (nzchar(receipt)) jsonlite::write_json(observations,
        file.path(dirname(receipt), "qc-global-owner-audit.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null")
})


test_that("QC ownership detachment preserves native built data and external expressions", {
    original <- ggplot2::get_last_plot()
    withr::defer(ggplot2::set_last_plot(original))
    object <- makeSavedQcProteinObject()
    aesthetics <- list(color_var = "group", shape_var = "batch")
    prior <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, x)) + ggplot2::geom_point()
    for (kind in c("pca", "rle", "density", "correlation")) {
        baseline <- switch(kind,
            pca = plotPca(object, grouping_variable = "group", label_column = "",
                shape_variable = "batch", title = "", font_size = 8),
            rle = plotRle(object, group = "group", yaxis_limit = c(-6, 6)),
            density = buildProtNormDensityPlot(object, aesthetics),
            correlation = buildProtNormCorrelationPlot(object, "group"))
        ggplot2::set_last_plot(prior)
        product <- buildProtNormQcPlotProduct(object, aesthetics, kind)
        expect_identical(ggplot2::get_last_plot(), prior)
        compare <- function(left, right) {
            left_built <- ggplot2::ggplot_build(left)
            right_built <- ggplot2::ggplot_build(right)
            expect_equal(left_built$data, right_built$data, tolerance = 1e-12)
            ticks <- function(built) lapply(built$layout$panel_params, function(panel)
                list(x = panel$x$get_labels(), y = panel$y$get_labels()))
            expect_identical(ticks(left_built), ticks(right_built))
            expect_equal(length(left$patches$plots), length(right$patches$plots))
            for (i in seq_along(left$patches$plots)) compare(left$patches$plots[[i]], right$patches$plots[[i]])
        }
        suppressWarnings(compare(baseline, product$plot))
    }
    ggplot2::set_last_plot(prior)
    expect_error(buildProtNormQcPlotProduct(object, aesthetics, "pca", plotPcaFn = function(...) {
        ggplot2::ggplot()
        stop("producer failed")
    }), "producer failed")
    expect_identical(ggplot2::get_last_plot(), prior)
    offset <- 7
    external <- ggplot2::ggplot(data.frame(x = 1:3), ggplot2::aes(x, x + offset)) + ggplot2::geom_point()
    external_product <- buildProtNormQcPlotProduct(object, aesthetics, "pca",
        plotPcaFn = function(...) external)
    expect_equal(ggplot2::ggplot_build(external_product$plot)$data,
        ggplot2::ggplot_build(external)$data)
})


test_that("retained native QC products allow their scientific input to be collected", {
    original <- ggplot2::get_last_plot()
    withr::defer(ggplot2::set_last_plot(original))
    for (kind in c("pca", "rle", "density", "correlation")) {
        ggplot2::set_last_plot(NULL)
        lifetime <- new.env(parent = emptyenv())
        lifetime$collected <- FALSE
        finalize <- function(value) lifetime$collected <- TRUE
        environment(finalize) <- list2env(list(lifetime = lifetime), parent = baseenv())
        product <- local({
            object <- makeSavedQcProteinObject()
            sentinel <- new.env(parent = emptyenv())
            reg.finalizer(sentinel, finalize, onexit = FALSE)
            args <- methods::slot(object, "args")
            args$ownership_sentinel <- sentinel
            methods::slot(object, "args") <- args
            buildProtNormQcPlotProduct(object, list(color_var = "group", shape_var = "batch"), kind)
        })
        gc(full = TRUE)
        gc(full = TRUE)
        expect_true(lifetime$collected, info = kind)
        expect_true(nrow(product$data) > 0L)
        expect_null(ggplot2::get_last_plot())
    }
})

test_that("the normalisation module edge writes the config_list method and RUV audit trail around the declaration", {
  original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
  root <- withr::local_tempdir()
  paths <- list(source_dir = root, protein_qc_dir = file.path(root, "qc"), time_dir = file.path(root, "time"),
    publication_graphs_dir = file.path(root, "graphs"))
  for (path in paths[-1]) dir.create(path)
  shiny::testServer(function(input, output, session) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "module-edge-audit", "DIA"),
      WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$.context$experiment_paths <- paths
    # KIT-013: the configuration lives on the workflow bus, which the module edge writes.
    workflow <- shiny::reactiveValues(state_manager = manager, config_list = list(normaliseBetweenSamples = list(method = "cyclicloess")))
    attachDIAStateSession(manager, session, workflow)
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    norm <- shiny::reactiveValues()
  }, {
    settings <- list(norm_method = "quantile", ruv_mode = "manual", ruv_k = 1, ruv_percentage = 50,
      color_variable = "group", shape_variable = NULL)
    notified <- character()
    for (pass in 1:2) {
      completed <- shiny::isolate(runProtNormNormalizationWorkflow(settings, workflow, norm, paths, "proteomics",
        "audit", checkMemoryUsageFn = function(...) NULL, generatePostNormalizationQcFn = function(...) NULL,
        generateRuvCorrectedQcFn = function(...) NULL, getRuvGroupingVariableFn = function() "group",
        showNotificationFn = function(message, ...) notified <<- c(notified, message),
        handleErrorFn = function(error, ...) stop(error), messageFn = function(...) NULL))
      expect_true(completed, info = paste("pass", pass))
      config <- shiny::isolate(workflow$config_list)
      expect_identical(config$normaliseBetweenSamples$method, "quantile", info = paste("pass", pass))
      expect_identical(config$ruvParameters$best_k, 1L, info = paste("pass", pass))
      expect_identical(config$ruvParameters$percentage_as_neg_ctrl, 50, info = paste("pass", pass))
      expect_identical(config$ruvParameters$num_neg_ctrl, length(norm$control_genes_index), info = paste("pass", pass))
      workflow$config_list <- list(normaliseBetweenSamples = list(method = "cyclicloess"))
    }
    # The second pass is an identity hit: the audit trail is written on hits as on misses (KIT-004 PC-04).
    expect_length(notified, 2L)
    expect_identical(manager$current_state, "ruv_corrected")
    closeDIAStateAdapter(manager)
  })
})
