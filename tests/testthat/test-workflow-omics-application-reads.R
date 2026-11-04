# KIT-082: the omics application reads no whole assay table on a managed session. Each read is proved against the whole-object
# read it replaces; the probe (with no window) runs the tabs' real outputs and observers, with the in-process tabs as its control.
library(testthat)

omicsReadCase <- function(omic) {
  if (identical(omic, "lipidomics")) {
    list(binding = lipidomicsDeclarationBinding, owner_of = lipidDeclarationsOwner, prefix = "lipid",
      summary = buildLipidQcS4DataSummaryUi, stats = function(object) buildLipidQcS4AssayStatsTable(object, datatableFn = function(table, ...) table),
      text = buildLipidQcS4FinalizeResultsText, validate = function(object) validateLipidQcS4FinalizeState(object, reqFn = identity),
      finder = findLipidDuplicateFeatureIDs,
      itsd = function(page, object) buildLipidItsdSelectionTable(assay_data = page, lipid_id_col = object@lipid_id_column, annotation_cols = object@annotation_id_column),
      volcano = "generateLipidDAVolcanoPlotGlimma")
  } else {
    list(binding = metabolomicsDeclarationBinding, owner_of = metabDeclarationsOwner, prefix = "metab",
      summary = buildMetabQcS4DataSummaryUi, stats = function(object) buildMetabQcS4AssayStatsDatatable(object, datatableFn = function(table, ...) table),
      text = buildMetabQcS4FinalizeResultsText, validate = function(object) validateMetabQcS4FinalizeState(object, reqFn = identity),
      finder = findMetabDuplicateFeatureIDs,
      itsd = function(page, object) buildItsdSelectionTable(assay_data = page, metabolite_id_col = object@metabolite_id_column, annotation_cols = object@annotation_id_column),
      volcano = "generateMetabDAVolcanoPlotGlimma")
  }
}

#' A managed tab's runner in this process: each action through the real worker entry point and acceptance.
omicsSyncRunner <- function(fixture, controller, calls = new.env()) {
  list(run = function(actions, done = NULL, failed = NULL) {
    for (action in actions) {
      calls$operations <- c(calls$operations, action$operation)
      result <- omicsWorkerAction(fixture, controller, action$operation, action$parameters, action$context)
      action$deliver(result)
    }
    if (is.function(done)) done()
    invisible(NULL)
  }, retrying = function() FALSE, progress = function(message) list(close = function() invisible(NULL)))
}

omicsReadMocks <- function(env = parent.frame()) testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
  finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR", .env = env)

for (omic in c("lipidomics", "metabolomics")) local({
  omic <- omic
  case <- omicsReadCase(omic)
  binding <- case$binding

  test_that(sprintf("KIT-082: the %s QC tab's counts, text and validation from the store equal the whole object's", omic), {
    skip_if_not(workflowOwnerCapable())
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    whole <- readWorkflowOmicsState(binding, owner)
    shell <- readWorkflowOmicsShell(binding, owner, measures = c("features", "missing"))
    expect_true(all(vapply(methods::slot(shell, binding$data_slot), nrow, 1L) == 0L))
    expect_identical(attr(shell, "assay_rows"), vapply(methods::slot(whole, binding$data_slot), nrow, 1L))
    expect_equal(attr(shell, "feature_counts"), omicsFeatureCounts(binding, whole))
    expect_identical(as.character(case$summary(shell)), as.character(case$summary(whole)))
    expect_equal(case$stats(shell), case$stats(whole))
    expect_identical(case$text(shell, "a"), case$text(whole, "a"))
    # The zero-row shell validates exactly as the object does.
    expect_s4_class(case$validate(shell), binding$class)
    expect_s4_class(case$validate(whole), binding$class)
    # The design, the slots and the parameters come with it.
    expect_identical(shell@design_matrix, whole@design_matrix)
    expect_identical(shell@args, whole@args)
    expect_identical(names(methods::slot(shell, binding$data_slot)[[1L]]), names(methods::slot(whole, binding$data_slot)[[1L]]))
  })

  test_that(sprintf("KIT-082: %s duplicate detection on the id column equals the whole object's", omic), {
    skip_if_not(workflowOwnerCapable())
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    whole <- readWorkflowOmicsState(binding, fixture$owner)
    ids <- readWorkflowOmicsIds(binding, fixture$owner)
    expect_identical(lengths(methods::slot(ids, binding$data_slot)), stats::setNames(rep(1L, length(methods::slot(whole, binding$data_slot))), names(methods::slot(whole, binding$data_slot))))
    found <- suppressMessages(case$finder(whole))
    expect_true(any(!vapply(found, is.null, logical(1))))
    expect_identical(suppressMessages(case$finder(ids)), found)
  })

  test_that(sprintf("KIT-082: the %s ITSD tables, missingness and counts read by pages equal the whole table's past one page", omic), {
    skip_if_not(workflowOwnerCapable())
    set.seed(82L)
    samples <- sprintf("S%02d", 1:6)
    assay <- function(rows, prefix) {
      ids <- c(sprintf("%s%05d", prefix, seq_len(rows)), sprintf("IS_%s%d", prefix, 1:5), sprintf("%s00003", prefix))
      values <- matrix(stats::rlnorm(length(ids) * length(samples), 8, 1), nrow = length(ids), dimnames = list(NULL, samples))
      values[sample(length(values), length(values) %/% 20L)] <- NA_real_
      values[sample(length(values), length(values) %/% 50L)] <- 0
      # Ties in the sort: a tenth of the rows have one value everywhere, so their spread is zero or missing.
      values[seq(1L, nrow(values), by = 10L), ] <- 5
      cbind(data.frame(feature = ids, stringsAsFactors = FALSE), as.data.frame(values, check.names = FALSE))
    }
    design <- data.frame(Run = samples, group = rep_len(c("A", "B"), 6L), stringsAsFactors = FALSE)
    assays <- list(big = assay(12700L, "P"), small = assay(40L, "N"))
    object <- if (identical(omic, "lipidomics")) methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = design, sample_id = "Run", group_id = "group",
      lipid_id_column = "feature", annotation_id_column = "feature") else methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = design,
      sample_id = "Run", group_id = "group", metabolite_id_column = "feature", annotation_id_column = "feature")
    frame <- environment()
    owner <- omicsTestOwner(object, "raw", omic_type = omic, config = list(qc = list(mode = "kit082")), env = frame)$state_manager
    whole <- readWorkflowOmicsState(binding, owner)
    tables <- readWorkflowOmicsItsdTables(binding, owner, c("big", "small", "absent"), function(page) case$itsd(page, whole))
    expect_identical(names(tables), c("big", "small"))
    for (name in names(tables)) expect_identical(tables[[name]], case$itsd(methods::slot(whole, binding$data_slot)[[name]], whole), info = name)
    shell <- readWorkflowOmicsShell(binding, owner, measures = c("features", "missing"))
    expect_equal(case$stats(shell), case$stats(whole))
    expect_identical(as.character(case$summary(shell)), as.character(case$summary(whole)))
    # The rows a filter keeps are in source order, past the page boundary.
    keep <- readWorkflowOmicsRows(binding, owner, NULL, function(page, id_column) grepl("^IS_|00003$", page[[id_column]]))
    table <- methods::slot(whole, binding$data_slot)$big
    expect_equal(methods::slot(keep, binding$data_slot)$big, {
      rows <- table[grepl("^IS_|00003$", table$feature), , drop = FALSE]
      rownames(rows) <- NULL
      rows
    })
  })

  test_that(sprintf("KIT-082: the %s DA heatmap reading its chosen rows equals the whole object's", omic), {
    skip_if_not(workflowOwnerCapable())
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(binding, owner, list())))
    source <- getWorkflowStateReference(owner)
    results <- suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDaAnalysis(binding, owner, source,
      list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$workflow$contrasts_tbl)))$results
    expect_false(is.null(results))
    state <- owner$current_state
    generate <- if (identical(omic, "lipidomics")) generateLipidDAHeatmap else generateMetabDAHeatmap
    with_whole <- results
    with_whole$theObject <- readWorkflowOmicsState(binding, owner, state)
    with_rows <- omicsDaResultsWithRowReader(binding, results, list(current_state_name = state), fixture$workflow)
    run <- function(list) suppressWarnings(suppressMessages(generate(list, selected_contrast = fixture$workflow$contrasts_tbl$contrasts[[1L]], selected_assay = names(methods::slot(with_whole$theObject, binding$data_slot))[[1L]], top_n = 10L, da_q_val_thresh = 1)))
    from_whole <- run(with_whole)
    from_rows <- run(with_rows)
    expect_false(is.null(from_whole))
    expect_identical(from_rows$plot@matrix, from_whole$plot@matrix)
  })

  test_that(sprintf("KIT-082: the %s DA tab holds no per-sample column and shows what the whole results show; its download has the same bytes", omic), {
    skip_if_not(workflowOwnerCapable())
    omicsReadMocks()
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(binding, owner, list())))
    source <- getWorkflowStateReference(owner)
    whole <- suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDaAnalysis(binding, owner, source,
      list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$workflow$contrasts_tbl)))$results
    controller <- omicsWorkerController(fixture)
    done <- omicsWorkerAction(fixture, controller, paste0(case$prefix, "_da_tab_results"))
    expect_identical(done$acceptance, "complete")
    tab <- resolveWorkflowStateValue(owner$.context, done$display$tab)
    long <- whole[[binding$long_table]]
    expect_gt(sum(grepl("^intensity[.]", names(long))), 0L)
    expect_false(any(grepl("^intensity[.]", names(tab[[binding$long_table]]))))
    expect_null(tab$contrasts_results)
    expect_null(tab$per_assay_results)
    expect_equal(as.data.frame(tab[[binding$long_table]]), as.data.frame(long[, !grepl("^intensity[.]", names(long)), drop = FALSE]))
    expect_identical(tab$significant_counts, whole$significant_counts)
    # Every display the tab builds from these rows, on each selection.
    contrast <- unique(long$comparison)[[1L]]
    assay <- unique(long$assay)[[1L]]
    for (selected in c(assay, "Combined", "All")) for (filter in c("all", "significant", "up", "down")) {
      shown <- function(results) if (identical(omic, "lipidomics")) {
        buildLipidDaResultsTableRenderWidget(results, contrast, selected, filter, 1, 0, 1000L, reqFn = identity)$x$data
      } else buildMetabDaResultsTableRenderOutput(results, contrast, selected, filter, 1, 0, 1000L)$x$data
      expect_identical(shown(tab), shown(whole), info = paste(selected, filter))
    }
    stats <- function(results) if (identical(omic, "lipidomics")) buildLipidDaSummaryStatsRenderText(results, contrast, assay, 0.05, reqFn = identity) else
      buildMetabDaSummaryStatsRenderOutput(results, contrast, assay, 0.05)
    expect_identical(stats(tab), stats(whole))
    static <- if (identical(omic, "lipidomics")) generateLipidDAVolcanoStatic else generateMetabDAVolcanoStatic
    # A plot's data carries every column of its input; the layers draw none of the per-sample ones.
    plotted <- function(results) {
      data <- as.data.frame(suppressWarnings(static(results, selected_contrast = contrast, selected_assay = assay, da_q_val_thresh = 1))$data)
      data[, !grepl("^intensity[.]", names(data)), drop = FALSE]
    }
    expect_equal(plotted(tab), plotted(whole))
    # The heatmap on the tab's rows chooses and draws what the whole results do.
    generate <- if (identical(omic, "lipidomics")) generateLipidDAHeatmap else generateMetabDAHeatmap
    state <- owner$current_state
    with_rows <- function(results) omicsDaResultsWithRowReader(binding, results, list(current_state_name = state), fixture$workflow)
    drawn <- function(results) suppressWarnings(suppressMessages(generate(with_rows(results), selected_contrast = contrast, selected_assay = assay, top_n = 10L, da_q_val_thresh = 1)))$plot@matrix
    expect_identical(drawn(tab), drawn(whole))
    # The download.
    manifest <- resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results)
    streamed <- tempfile(fileext = ".csv")
    writeWorkflowOmicsDaDownload(binding, owner, manifest, tab$download_columns, streamed)
    written <- tempfile(fileext = ".csv")
    utils::write.csv(long, written, row.names = FALSE)
    expect_identical(digest::digest(file = streamed, algo = "sha256"), digest::digest(file = written, algo = "sha256"))
    # Through each omic's handler builder, the managed list writes that file.
    tab$write_download <- function(file) writeWorkflowOmicsDaDownload(binding, owner, manifest, tab$download_columns, file)
    through <- tempfile(fileext = ".csv")
    if (identical(omic, "lipidomics")) {
      handler <- buildLipidDaResultsDownloadOutputHandler(tab, buildResultsDownloadHandlerFn = function(daResultsReactive) buildLipidDaResultsDownloadHandler(daResultsReactive, downloadHandlerFactory = function(filename, content) content, reqFn = identity))
      handler(through)
    } else {
      handler <- buildMetabDaResultsDownloadHandler(list(da_results_list = tab), downloadHandler = function(filename, content) content)
      handler(through)
    }
    expect_identical(digest::digest(file = through, algo = "sha256"), digest::digest(file = written, algo = "sha256"))
  })

  test_that(sprintf("KIT-082: a managed %s QC finalise completes when its tracking figure fails, as in process, and forwards the worker's metrics", omic), {
    finalize <- if (identical(omic, "lipidomics")) runLipidQcS4FinalizeWorkflow else runMetabQcS4FinalizeWorkflow
    tracking <- if (identical(omic, "lipidomics")) updateLipidQcS4TrackingPlot else updateMetabQcS4TrackingPlot
    runner <- list(run = function(actions, done = NULL, failed = NULL) failed(simpleError("the worker failed")))
    testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) runner, .package = "MultiScholaR")
    seen <- new.env()
    finalize(list(state_manager = "manager"), omic, function(plot) seen$plot <- list(plot), "output",
      getFinalizeStateFn = function(stateManager) "state", validateFinalizeStateFn = function(currentS4) currentS4,
      saveCompletedStateFn = function(...) NULL, completeTabStatusFn = function(workflowData) seen$complete <- TRUE,
      getFinalizeHistoryFn = function(stateManager) "history", reportFinalizeSuccessFn = function(...) seen$success <- TRUE,
      reportFinalizeErrorFn = function(error) seen$error <- conditionMessage(error))
    expect_true(seen$complete)
    expect_true(seen$success)
    expect_null(seen$error)
    expect_identical(seen$plot, list(NULL))
    # A failure in the tab's own completion still reaches the tab's report, as the in-process path's `tryCatch` does.
    seen <- new.env()
    finalize(list(state_manager = "manager"), omic, function(plot) NULL, "output",
      getFinalizeStateFn = function(stateManager) "state", validateFinalizeStateFn = function(currentS4) currentS4,
      saveCompletedStateFn = function(...) NULL, completeTabStatusFn = function(workflowData) stop("no status"),
      getFinalizeHistoryFn = function(stateManager) "history", reportFinalizeSuccessFn = function(...) seen$success <- TRUE,
      reportFinalizeErrorFn = function(error) seen$error <- conditionMessage(error))
    expect_identical(seen$error, "no status")
    expect_null(seen$success)
    # The metrics argument reaches the kernel, with no object.
    arguments <- NULL
    updater <- function(...) { arguments <<- list(...); "plot" }
    if (identical(omic, "lipidomics")) tracking(NULL, omic, function(plot) NULL, updateLipidFilteringFn = updater, metrics = "metrics") else
      tracking(NULL, omic, function(plot) NULL, updateMetaboliteFilteringFn = updater, metrics = "metrics")
    expect_identical(arguments$metrics, "metrics")
    expect_null(arguments$theObject)
  })

  test_that(sprintf("KIT-082: a store read error building the %s ITSD tables leaves the session running", omic), {
    skip_if_not(workflowOwnerCapable())
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    assays <- names(methods::slot(fixture$object, binding$data_slot))
    testthat::local_mocked_bindings(readWorkflowOmicsItsdTables = function(...) stop("the store could not be read"), .package = "MultiScholaR")
    output <- new.env()
    norm <- list(assay_names = assays)
    if (identical(omic, "lipidomics")) {
      expect_no_error(suppressWarnings(registerLipidNormItsdTableOutputs(output, fixture$workflow, norm, observeFn = function(code) eval(substitute(code), parent.frame()), reqFn = identity)))
    } else {
      expect_no_error(suppressWarnings(runMetabNormItsdSelectionTableObserverShell(norm, fixture$workflow, output, reqFn = identity)))
    }
    expect_length(ls(output), 0L)
  })

  test_that(sprintf("KIT-082: the %s worker's filtering metrics, PCA, session file and integration object equal the whole object's", omic), {
    skip_if_not(workflowOwnerCapable())
    omicsReadMocks()
    fixture <- case$owner_of(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$owner))
    owner <- fixture$owner
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(binding, owner, list())))
    controller <- omicsWorkerController(fixture)
    whole <- owner$getState(owner$current_state)
    state <- owner$current_state
    history <- owner$getHistory()
    # The filtering metrics.
    metrics <- omicsWorkerAction(fixture, controller, paste0(case$prefix, "_filtering_metrics"))
    expect_equal(withWorkflowStateScope(owner, function() readWorkflowOmicsFilteringMetrics(owner, metrics)), binding$filtering_metrics(whole))
    # The PCA: the same points, with no copy of the state in the plot.
    pca <- omicsWorkerAction(fixture, controller, paste0(case$prefix, "_final_pca"), list(state = owner$current_state, grouping_variable = "group", shape_variable = NULL))
    expect_null(pca$display$message)
    plot <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, pca$display$plot))
    in_process <- suppressWarnings(suppressMessages(binding$final_pca(readWorkflowOmicsState(binding, owner), "group", NULL)))
    draw <- function(paint) {
      path <- tempfile(fileext = ".png")
      grDevices::png(path, width = 600, height = 800)
      paint()
      grDevices::dev.off()
      digest::digest(file = path, algo = "sha256")
    }
    expect_identical(draw(function() { grid::grid.newpage(); grid::grid.draw(plot) }), draw(function() print(in_process)))
    expect_length(grepRaw(charToRaw(binding$class), serialize(plot, NULL), all = TRUE), 0L)
    # The integration object.
    path <- file.path(fixture$paths$source_dir, "integration.RDS")
    done <- omicsWorkerAction(fixture, controller, paste0(case$prefix, "_integration_export"), list(path = path, state_name = owner$current_state))
    expect_identical(done$acceptance, "complete")
    expect_equal(readRDS(path), whole)
    # The session file: the worker's is the list the controller collected with the state as the in-process export writes it.
    source_dir <- fixture$paths$source_dir
    collect_with <- function(state) function(workflowData, normData, input) list(r6_current_state_name = owner$current_state, current_s4_object = state, config_list = list(a = 1), omic_type = omic)
    staged <- NULL
    runner <- omicsSyncRunner(fixture, controller)
    startWorkflowManagedSessionExport(owner, runner, fixture$workflow, NULL, NULL, source_dir, collect_with(NULL), paste0(case$prefix, "_"),
      publish = function() {
        request <- owner$.context$export_pending
        staged <<- readRDS(file.path(workflowExportPaths(source_dir, request$request_id)$request_dir, "session.rds"))
      }, failed = function(error) stop(error))
    expect_false(is.null(staged))
    expect_identical(names(staged), names(collect_with(whole)()))
    expect_equal(staged$current_s4_object, whole)
    # None of them committed a state.
    expect_identical(owner$current_state, state)
    expect_identical(owner$getHistory(), history)
  })
})

#' The tabs' journey with every site the application read the matrix at: the module's real outputs and observers.
omicsAppServer <- function(fixture, omic, probe, managed) {
  base <- omicsTabServer(fixture, omic, probe, managed)
  function(input, output, session) {
    base(input, output, session)
    workflow <- fixture$workflow
    paths <- fixture$paths
    norm <- session$userData$probe_norm
    da <- session$userData$probe_da
    filter_plot <- shiny::reactiveVal(NULL)
    values <- shiny::reactiveValues(workflow_args_saved = FALSE)
    dirs <- stats::setNames(list(list(source_dir = paths$source_dir, integration_dir = file.path(paths$source_dir, "integration"), base_dir = paths$source_dir)), omic)
    quiet <- function(...) invisible(NULL)
    if (identical(omic, "lipidomics")) {
      registerLipidQcS4DataSummaryOutput(output, workflow)
      registerLipidQcS4AssayStatsOutput(output, workflow)
      shiny::observeEvent(input$detect, probe$detected <- handleLipidDuplicateDetection(workflow))
      shiny::observeEvent(input$finalize, runLipidQcS4FinalizeWorkflow(workflow, omic, filter_plot, output))
      registerLipidNormItsdTableOutputs(output, workflow, norm)
      mod_lipid_qc_itsd_server("itsd", workflow, omic, "kit082")
      output$final_qc <- buildLipidNormFinalQcPlotRenderer(norm, function(...) list(color_var = "group", shape_var = NULL), workflow)()
      # The test's DA holder is not reactive, so its outputs are registered once it holds the results.
      session$userData$register_da_outputs <- function() {
        registerLipidDaVolcanoGlimmaOutput(output, input, da, workflow)
        registerLipidDaHeatmapPlotOutput(output, input, da, workflow)
      }
      shiny::observeEvent(input$export_session, handleLipidNormExportSession(input, workflow, paths, "kit082", norm, addLog = quiet))
      shiny::observeEvent(input$save_args, handleLipidSummarySaveWorkflowArgs(input, values, output, dirs, omic, workflow))
    } else {
      output$data_summary <- shiny::renderUI(buildMetabQcS4DataSummaryRenderOutput(workflow$state_manager))
      output$assay_stats_table <- DT::renderDT(buildMetabQcS4AssayStatsRenderOutput(workflow$state_manager))
      shiny::observeEvent(input$detect, probe$detected <- detectMetabDuplicateFeatures(workflow$state_manager))
      shiny::observeEvent(input$finalize, runMetabQcS4FinalizeWorkflow(workflow, omic, filter_plot, output))
      shiny::observe(runMetabNormItsdSelectionTableObserverShell(norm, workflow, output))
      mod_metab_qc_itsd_server("itsd", workflow, omic, "kit082")
      output$final_qc <- renderMetabNormFinalQcPlot(norm, workflow, colorVariableFn = function() "group", shapeVariableFn = function() NULL)
      session$userData$register_da_outputs <- function() registerMetabDaVisualizationOutputs(input, output, da, workflow)
      shiny::observeEvent(input$export_session, runMetabNormExportSessionObserverShell(workflow, norm, input, paths, "kit082"))
      shiny::observeEvent(input$save_args, runMetabSummarySaveWorkflowArgsObserverShell(input, dirs, omic, workflow, values, output))
    }
  }
}

omicsAppJourney <- function(session, output, omic, prefix, assays, contrast) {
  do.call(session$setInputs, omicsTabInputs())
  apply_intensity <- if (identical(omic, "metabolomics")) list(`intensity-intensity_cutoff_percentile` = 1, `intensity-proportion_below_cutoff` = 0.5, `intensity-apply_filter` = 1) else list(apply_filter = 1)
  do.call(session$setInputs, apply_intensity)
  session$flushReact()
  session$setInputs(detect = 1)
  session$flushReact()
  session$setInputs(resolve_duplicates = 1)
  session$flushReact()
  force(output$data_summary)
  force(output$assay_stats_table)
  session$setInputs(`itsd-is_pattern` = "^IS_")
  session$setInputs(`itsd-analyze_is` = 1)
  session$flushReact()
  session$setInputs(finalize = 1)
  session$flushReact()
  session$userData$probe_norm$assay_names <- assays
  session$flushReact()
  session$setInputs(run_normalization = 1)
  session$flushReact()
  session$setInputs(apply_correlation_filter = 1)
  session$flushReact()
  # The first read asks the worker for the figure; the second draws what it delivered.
  try(force(output$final_qc), silent = TRUE)
  session$flushReact()
  force(output$final_qc)
  session$setInputs(experiment_label = "kit082", description = "probe")
  session$setInputs(export_session = 1)
  session$flushReact()
  session$setInputs(run_da = 1)
  session$flushReact()
  session$userData$register_da_outputs()
  session$setInputs(da_q_val_thresh = 1, volcano_contrast = contrast, volcano_assay = assays[[1L]], heatmap_contrast = contrast, heatmap_assay = assays[[1L]],
    heatmap_top_n = 10L, heatmap_clustering = "both", heatmap_scaling = "row", heatmap_color_scheme = "RdBu", heatmap_cluster_method = "complete",
    heatmap_distance_method = "euclidean", heatmap_show_labels = TRUE, heatmap_tree_cut_method = "none")
  session$flushReact()
  force(output$heatmap_plot)
  session$flushReact()
  page <- output$volcano_glimma
  session$setInputs(save_args = 1)
  session$flushReact()
  page
}

for (omic in c("lipidomics", "metabolomics")) local({
  omic <- omic
  test_that(sprintf("KIT-082: on a managed session the %s application reads no assay table on the whole journey; the same tabs in process do", omic), {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("Glimma")
    case <- omicsReadCase(omic)
    binding <- case$binding
    for (managed in c(TRUE, FALSE)) {
      fixture <- case$owner_of(file.path(withr::local_tempdir(), if (managed) "managed" else "in-process"))
      withr::defer(closeDIAStateAdapter(fixture$owner))
      probe <- new.env()
      assays <- names(methods::slot(fixture$object, binding$data_slot))
      # The selectors offer the results' comparison names.
      contrast <- fixture$workflow$contrasts_tbl$contrasts[[1L]]
      if (managed) omicsManagedChainMocks(fixture, probe, window = FALSE) else {
        table <- materializeWorkflowTable
        components <- readWorkflowExecutionComponents
        testthat::local_mocked_bindings(materializeWorkflowTable = function(context, reference, ...) {
          probe$tables <- c(probe$tables, reference$entity_id)
          table(context, reference, ...)
        }, materializeDIATable = function(context, reference, ...) {
          probe$tables <- c(probe$tables, reference$entity_id)
          table(context, reference, ...)
        }, .package = "MultiScholaR")
      }
      suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(omicsAppServer(fixture, omic, probe, managed), {
        session$userData$probe_owner <- fixture$owner
        probe$page <- omicsAppJourney(session, output, omic, case$prefix, assays, contrast)
      }))))
      materialised <- intersect(probe$tables, omicsAssayTableIds(fixture$owner, binding))
      if (!managed) {
        expect_gt(length(materialised), 0L)
        next
      }
      expect_length(materialised, 0L)
      expect_true(all(paste0(case$prefix, c("_da_tab_results", "_filtering_metrics", "_final_pca", "_volcano_page", "_session_export", "_integration_export")) %in% probe$operations),
        info = paste(omic, paste(probe$operations, collapse = " ")))
      # The tab holds the DA rows without reading a contrast's product; only the download, on the user's click, reads them.
      manifest <- resolveWorkflowStateValue(fixture$owner$.context, fixture$owner$.context$snapshot$da_results)
      contrast_products <- unlist(lapply(manifest$results, function(entry) entry$products$da_contrast$entity_id))
      expect_gt(length(contrast_products), 0L)
      expect_length(intersect(probe$hydrated, contrast_products), 0L)
      columns <- c(id = "character")
      streamed <- tempfile(fileext = ".csv")
      long <- readWorkflowOmicsDaResults(binding, fixture$owner, manifest)[[binding$long_table]]
      probe$hydrated <- NULL
      writeWorkflowOmicsDaDownload(binding, fixture$owner, manifest, vapply(long, typeof, character(1)), streamed)
      expect_gt(length(intersect(probe$hydrated, contrast_products)), 0L)
      # A read of the store keeps the rows or the column it was asked for, never the table: the counter above sees only whole reads.
      whole <- max(vapply(methods::slot(fixture$object, binding$data_slot), function(table) nrow(table) * ncol(table), integer(1)))
      expect_gt(length(probe$kept_cells), 0L)
      expect_lt(max(probe$kept_cells), whole)
      expect_true(all(probe$statuses == "delivered"), info = omic)
      expect_true(paste0(case$prefix, "_qc_complete") %in% unlist(fixture$owner$getHistory()), info = omic)
      # What the tabs show and write.
      expect_gt(sum(probe$detected$totalDuplicates), 0L)
      expect_match(paste(unlist(probe$page), collapse = " "), "<iframe")
      latest <- readRDS(file.path(fixture$paths$source_dir, paste0(case$prefix, "_filtered_session_data_latest.rds")))
      expect_s4_class(latest$current_s4_object, binding$class)
      integration <- readRDS(file.path(fixture$paths$source_dir, "integration", sprintf("%s_kit082_final_s4.RDS", omic)))
      expect_s4_class(integration, binding$class)
      expect_identical(nrow(methods::slot(integration, binding$data_slot)[[1L]]) > 0L, TRUE)
      expect_true(file.exists(file.path(fixture$paths$source_dir, "study_parameters.txt")))
    }
  })
})

test_that("KIT-082: the runner runs what waits in order, and tells what waits when an action fails", {
  skip_if_not(workflowOwnerCapable())
  fixture <- lipidDeclarationsOwner(withr::local_tempdir())
  withr::defer(closeDIAStateAdapter(fixture$owner))
  probe <- new.env()
  omicsManagedChainMocks(fixture, probe)
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
    runner <- installWorkflowOmicsActions(session, fixture$workflow, "lipidomics", "kit082", fixture$paths)
    session$userData$runner <- runner
  }, {
    runner <- session$userData$runner
    seen <- character()
    action <- function(name) omicsManagedProduct("lipid_filtering_metrics", list(), function(result) seen <<- c(seen, name))
    # Each `run` inside a delivery finds the action still active, so it waits.
    runner$run(list(omicsManagedProduct("lipid_filtering_metrics", list(), function(result) {
      seen <<- c(seen, "first")
      runner$run(list(action("second")))
      runner$run(list(action("third")))
    })))
    session$flushReact()
    expect_identical(seen, c("first", "second", "third"))
    # A display asks once per key, whatever renders it, and reads what was delivered.
    asked <- 0L
    read <- function(result) {
      asked <<- asked + 1L
      list(value = "delivered")
    }
    expect_null(shiny::isolate(workflowOmicsDisplay(runner, "key", "lipid_filtering_metrics", list(), read)))
    session$flushReact()
    expect_identical(shiny::isolate(workflowOmicsDisplay(runner, "key", "lipid_filtering_metrics", list(), read)), list(value = "delivered"))
    expect_identical(asked, 1L)
  }))))
  # A failure drops what waits and tells it.
  told <- character()
  hook <- new.env()
  testthat::local_mocked_bindings(pollWorkflowAction = function(owner) list(status = if (isTRUE(hook$fail)) "failed" else "running"), .package = "MultiScholaR")
  suppressWarnings(suppressMessages(utils::capture.output(shiny::testServer(function(input, output, session) {
    session$userData$runner <- installWorkflowOmicsActions(session, fixture$workflow, "lipidomics", "kit082", fixture$paths)
  }, {
    runner <- session$userData$runner
    action <- function() omicsManagedProduct("lipid_filtering_metrics", list(), function(result) stop("not delivered"))
    runner$run(list(action()), failed = function(error) told <<- c(told, "running"))
    session$flushReact()
    runner$run(list(action()), failed = function(error) told <<- c(told, "waiting"))
    expect_null(shiny::isolate(workflowOmicsDisplay(runner, "broken", "lipid_filtering_metrics", list())))
    hook$fail <- TRUE
    session$elapse(100)
    session$flushReact()
    # The failure was not this display's: it asks again, and only its own failure is its message, shown once.
    expect_null(shiny::isolate(workflowOmicsDisplay(runner, "broken", "lipid_filtering_metrics", list())))
    expect_false(is.null(shiny::isolate(runner$displays[[".asked"]])))
    session$elapse(1000)
    session$flushReact()
    # Read once into a variable: a failure is shown once, and `expect_match` evaluates its argument twice.
    shown <- shiny::isolate(workflowOmicsDisplay(runner, "broken", "lipid_filtering_metrics", list()))
    expect_match(shown$message, "managed action ended")
    hook$fail <- FALSE
    expect_null(shiny::isolate(workflowOmicsDisplay(runner, "broken", "lipid_filtering_metrics", list())))
  }))))
  expect_identical(told, c("running", "waiting"))
})

test_that("KIT-082: a failed publication's retry resumes the stopped action, keeps the first request's callbacks and lets other work wait", {
  applied <- character()
  pending <- shiny::reactiveVal(NULL)
  active <- shiny::reactiveVal(NULL)
  callbacks <- NULL
  testthat::local_mocked_bindings(workflowSessionApplication = function(session) "application",
    installWorkflowManagedStep = function(session, a, workflow_data, b, publication, ui, deliver, read_result) {
      callbacks <<- list(ui = ui, deliver = deliver)
      list(apply = function(parameters, operation, context) {
        applied <<- c(applied, operation)
        active(operation)
      }, active = active, pending = pending)
    }, .package = "MultiScholaR")
  session <- shiny::MockShinySession$new()
  runner <- installWorkflowOmicsActions(session, list(), "lipidomics", "label", list())
  first_done <- 0L
  second_done <- 0L
  action <- function(operation) omicsManagedProduct(operation, list(), function(result) NULL)
  runner$run(list(action("first")), done = function() first_done <<- first_done + 1L)
  expect_identical(applied, "first")
  expect_false(runner$retrying())
  shiny::isolate(pending("the publication failed"))
  expect_true(runner$retrying())
  # The same actions are the retry. Nothing else begins, and its progress has nothing left open.
  shown <- runner$progress("Exporting")
  expect_null(shown$close())
  runner$run(list(action("first")), done = function() second_done <<- second_done + 1L)
  expect_identical(applied, c("first", "first"))
  # A display during the wait says so, and asks for nothing.
  displays <- shiny::isolate(workflowOmicsDisplay(runner, "key", "other", list()))
  expect_match(displays$message, "Press that step's button again")
  expect_identical(applied, c("first", "first"))
  # Other work takes its turn after the retry.
  other_done <- 0L
  runner$run(list(action("other")), done = function() other_done <<- other_done + 1L)
  expect_identical(applied, c("first", "first", "first"))
  shiny::isolate({
    pending(NULL)
    active(NULL)
  })
  callbacks$ui$delivered()
  expect_identical(first_done, 1L)
  expect_identical(second_done, 0L)
  expect_identical(applied, c("first", "first", "first", "other"))
  shiny::isolate(active(NULL))
  callbacks$ui$delivered()
  expect_identical(other_done, 1L)
})

test_that("KIT-082: the DA download's writer handles an empty piece, a type that differs between pieces, and no rows at all", {
  pieces <- list(list(LCMS_Pos = data.frame(id = "a", n = 100000L, stringsAsFactors = FALSE)),
    list(LCMS_Pos = data.frame(id = character(), n = integer(), stringsAsFactors = FALSE), LCMS_Neg = data.frame(id = "b", n = 100000, extra = 2.5, stringsAsFactors = FALSE)))
  manifest <- list(assays = c("LCMS_Pos", "LCMS_Neg"), results = list(list(products = list(da_contrast = 1L)), list(products = list(da_contrast = 2L)), list(products = list())))
  manager <- list(.context = list(store = "store"))
  testthat::local_mocked_bindings(hydrateSavedWorkflowEntity = function(store, reference, ...) {
    list(value = list(assays = lapply(pieces[[reference]], function(table) list(table = table))))
  }, .package = "MultiScholaR")
  whole <- dplyr::bind_rows(pieces[[1L]]$LCMS_Pos, pieces[[2L]]$LCMS_Pos, pieces[[2L]]$LCMS_Neg)
  streamed <- tempfile(fileext = ".csv")
  writeWorkflowOmicsDaDownload(NULL, manager, manifest, vapply(whole, typeof, character(1)), streamed)
  written <- tempfile(fileext = ".csv")
  utils::write.csv(whole, written, row.names = FALSE)
  expect_identical(readLines(streamed), readLines(written))
  # No rows anywhere: the header alone.
  pieces[[1L]]$LCMS_Pos <- pieces[[1L]]$LCMS_Pos[0L, ]
  pieces[[2L]]$LCMS_Neg <- pieces[[2L]]$LCMS_Neg[0L, ]
  empty <- dplyr::bind_rows(pieces[[1L]]$LCMS_Pos, pieces[[2L]]$LCMS_Neg)
  streamed <- tempfile(fileext = ".csv")
  writeWorkflowOmicsDaDownload(NULL, manager, manifest, vapply(empty, typeof, character(1)), streamed)
  written <- tempfile(fileext = ".csv")
  utils::write.csv(empty, written, row.names = FALSE)
  expect_identical(readLines(streamed), readLines(written))
})

test_that("KIT-082: the lipid download serves each run's results through the registered handler, as metabolomics does", {
  da <- new.env()
  da$da_results_list <- NULL
  output <- new.env()
  contents <- NULL
  registerLipidDaResultsDownloadOutput(output, da, buildResultsDownloadOutputHandlerFn = function(...) {
    buildLipidDaResultsDownloadOutputHandler(..., buildResultsDownloadHandlerFn = function(daResultsReactive) daResultsReactive)
  })
  # A download before any run does not fix the answer.
  expect_null(output$download_da_results())
  da$da_results_list <- list(da_lipids_long = data.frame(run = 1L))
  expect_identical(output$download_da_results()$run, 1L)
  da$da_results_list <- list(da_lipids_long = data.frame(run = 2L))
  expect_identical(output$download_da_results()$run, 2L)
  writer <- function(file) NULL
  da$da_results_list <- list(da_lipids_long = data.frame(run = 3L), write_download = writer)
  expect_identical(output$download_da_results(), writer)
})

test_that("KIT-082: a managed DA whose tab rows cannot be read clears the tab; once it has them a failed publication only warns", {
  manifest <- list(accepted = TRUE)
  sequence <- NULL
  runner <- list(run = function(actions, done = NULL, failed = NULL) {
    if (is.null(sequence)) {
      sequence <<- "contrasts"
      return(done())
    }
    sequence <<- list(actions = actions, done = done, failed = failed)
  })
  testthat::local_mocked_bindings(getWorkflowStateReference = function(manager, ...) "source", omicsDaContrastSteps = function(...) list(),
    withWorkflowStateScope = function(manager, code) code(), acceptWorkflowOmicsDaManifest = function(...) list(manifest = manifest),
    resolveWorkflowStateValue = function(context, reference) list(download_columns = c(id = "character"), rows = reference), .package = "MultiScholaR")
  managed <- function() {
    sequence <<- NULL
    told <- new.env()
    runWorkflowOmicsManagedDa(lipidomicsDeclarationBinding, runner, list(.context = list()), "source", list(), NULL,
      deliver = function(results, outcome) told$delivered <- list(results = results, outcome = outcome),
      failed = function(error) told$failed <- conditionMessage(error), invalidate = function() told$cleared <- TRUE)
    told
  }
  told <- managed()
  sequence$failed(simpleError("the tab's read failed"))
  expect_true(told$cleared)
  expect_identical(told$failed, "the tab's read failed")
  expect_null(told$delivered)
  told <- managed()
  sequence$actions[[1L]]$deliver(list(display = list(tab = "reference")))
  sequence$failed(simpleError("the files failed"))
  expect_null(told$cleared)
  expect_null(told$failed)
  expect_false(told$delivered$outcome$success)
  expect_identical(told$delivered$results$rows, "reference")
  expect_true(is.function(told$delivered$results$write_download))
})
