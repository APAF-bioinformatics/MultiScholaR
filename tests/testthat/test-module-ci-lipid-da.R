library(testthat)

test_that("MCI-022.1 reload matrix covers valid lipid sessions and failure branches", {
  for (layout in c("single_pos", "pos_neg", "pos_gcms")) {
    tmp_root <- tempfile(paste0("mci-022-reload-", layout, "-"))
    dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
    withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

    session_file <- file.path(tmp_root, "lipid_filtered_session_data_latest.rds")
    session_data <- module_ci_lipid_da_session_data(layout = layout)
    saveRDS(session_data, session_file)

    resolved <- resolveLipidDaSessionFile(
      experimentPaths = list(source_dir = tmp_root, export_dir = tempfile("unused-"))
    )
    expect_identical(resolved$sessionFile, session_file)

    workflow <- module_ci_lipid_da_workflow(session_data$current_s4_object)
    da_data <- module_ci_lipid_da_data_env()
    capture <- module_ci_lipid_da_capture_ui()

    load_state <- loadLipidDaSessionFromFile(
      sessionFile = session_file,
      workflowData = workflow,
      daData = da_data,
      session = "module-ci-lipid-da",
      showLoadingFn = function() list(notificationId = "loading_session"),
      finalizeSuccessFn = function(...) list(successMessage = "Session loaded successfully!"),
      finalizeErrorFn = function(errorMessage, ...) stop(errorMessage),
      restorePostReadFn = module_ci_lipid_da_restore_post_read(capture)
    )

    expect_identical(load_state$status, "success")
    expect_s4_class(da_data$current_s4_object, "LipidomicsAssayData")
    expect_identical(names(da_data$current_s4_object@lipid_data), session_data$assay_names)
    expect_identical(da_data$assays_available, session_data$assay_names)
    expect_identical(workflow$contrasts_tbl$contrasts, session_data$contrasts_tbl$contrasts)
    expect_identical(da_data$formula_from_s4, "~ 0 + group")
    expect_true(any(vapply(capture$select_updates, function(update) {
      identical(update$inputId, "volcano_assay") &&
        identical(update$choices, c("Combined", session_data$assay_names))
    }, logical(1))))
  }

  stale_export_dir <- tempfile("mci-022-stale-export-")
  dir.create(stale_export_dir, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(stale_export_dir, recursive = TRUE, force = TRUE))
  stale_file <- file.path(stale_export_dir, "lipid_filtered_session_data_latest.rds")
  saveRDS(module_ci_lipid_da_session_data("pos_gcms"), stale_file)
  stale_resolved <- resolveLipidDaSessionFile(
    experimentPaths = list(source_dir = file.path(stale_export_dir, "missing"), export_dir = stale_export_dir)
  )
  expect_identical(stale_resolved$sessionFile, stale_file)

  missing_capture <- module_ci_lipid_da_capture_ui()
  missing_resolved <- resolveLipidDaSessionFile(
    experimentPaths = list(source_dir = stale_export_dir, export_dir = stale_export_dir),
    fileExistsFn = function(path) FALSE,
    notifySessionFileMissingFn = function(sessionFile) {
      missing_capture$showNotification(sprintf("Session file not found: %s", sessionFile), type = "error")
      list(sessionFile = sessionFile)
    }
  )
  expect_null(missing_resolved)
  expect_match(missing_capture$notifications[[1L]][[1L]], "Session file not found", fixed = TRUE)

  malformed_file <- file.path(stale_export_dir, "malformed.rds")
  writeLines("not an rds", malformed_file)
  malformed <- loadLipidDaSessionFromFile(
    sessionFile = malformed_file,
    workflowData = module_ci_lipid_da_workflow(),
    daData = module_ci_lipid_da_data_env(),
    session = "module-ci-lipid-da",
    showLoadingFn = function() list(notificationId = "loading_session"),
    finalizeSuccessFn = function(...) stop("malformed session should not succeed"),
    finalizeErrorFn = function(errorMessage, ...) list(errorMessage = errorMessage)
  )
  expect_identical(malformed$status, "error")
  expect_true(nzchar(malformed$errorMessage))

  mismatch_file <- file.path(stale_export_dir, "mismatch.rds")
  mismatch_session <- module_ci_lipid_da_session_data(
    layout = "pos_gcms",
    assay_names = c("LCMS_Pos", "STALE_ASSAY")
  )
  saveRDS(mismatch_session, mismatch_file)
  mismatch_capture <- module_ci_lipid_da_capture_ui()
  mismatch <- loadLipidDaSessionFromFile(
    sessionFile = mismatch_file,
    workflowData = module_ci_lipid_da_workflow(mismatch_session$current_s4_object),
    daData = module_ci_lipid_da_data_env(),
    session = "module-ci-lipid-da",
    showLoadingFn = function() list(notificationId = "loading_session"),
    finalizeSuccessFn = function(...) stop("mismatched session should not succeed"),
    finalizeErrorFn = function(errorMessage, ...) list(errorMessage = errorMessage),
    restorePostReadFn = module_ci_lipid_da_restore_post_read(mismatch_capture)
  )
  expect_identical(mismatch$status, "error")
  expect_match(mismatch$errorMessage, "assay_names do not match", fixed = TRUE)
})

test_that("MCI-022.2 formula and contrast preflight rejects invalid DA inputs before analysis", {
  object <- module_ci_lipid_da_object(layout = "pos_gcms")
  valid_cases <- list(
    two_group = list(
      contrasts = module_ci_lipid_da_contrasts("two_group"),
      formula = "~ 0 + group"
    ),
    multi_group = list(
      contrasts = module_ci_lipid_da_contrasts("multi_group"),
      formula = "~ 0 + group"
    ),
    batch_aware = list(
      contrasts = module_ci_lipid_da_contrasts("batch_aware"),
      formula = "~ 0 + group + batch"
    ),
    reversed = list(
      contrasts = module_ci_lipid_da_contrasts("reversed"),
      formula = "~ 0 + group"
    )
  )

  for (case in valid_cases) {
    preflight <- validateLipidDaAnalysisPreflight(
      currentS4 = object,
      contrastsTbl = case$contrasts,
      formulaString = case$formula
    )
    expect_true(preflight$valid)
    expect_true(all(case$contrasts$contrasts %in% preflight$contrast_validation$contrasts_tbl$contrasts))
  }

  duplicate <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("duplicate"),
    formulaString = "~ 0 + group"
  )
  expect_false(duplicate$valid)
  expect_true(any(grepl("duplicate friendly contrast labels", duplicate$errors, fixed = TRUE)))

  invalid_term <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("invalid_term"),
    formulaString = "~ 0 + group"
  )
  expect_false(invalid_term$valid)
  expect_true(any(grepl("references terms absent", invalid_term$errors, fixed = TRUE)))

  invalid_formula <- validateLipidDaAnalysisPreflight(
    currentS4 = object,
    contrastsTbl = module_ci_lipid_da_contrasts("two_group"),
    formulaString = "~ 0 + missing_factor"
  )
  expect_false(invalid_formula$valid)
  expect_true(any(grepl("cannot produce a model matrix", invalid_formula$errors, fixed = TRUE)))

  no_contrasts_context <- prepareLipidDaRunAnalysisContext(
    daData = list2env(list(current_s4_object = object), parent = emptyenv()),
    workflowData = module_ci_lipid_da_workflow(object, contrasts = module_ci_lipid_da_contrasts("no_contrasts")),
    notify = function(...) invisible(NULL),
    logInfo = function(...) invisible(NULL)
  )
  expect_null(no_contrasts_context)

  invalid_workflow <- module_ci_lipid_da_workflow(
    object,
    contrasts = module_ci_lipid_da_contrasts("invalid_term")
  )
  capture <- module_ci_lipid_da_capture_ui()
  shell_called <- FALSE
  analysis_state <- handleLipidDaRunAnalysisPreflight(
    analysisContext = list(
      currentS4 = object,
      contrastsTbl = invalid_workflow$contrasts_tbl,
      runningNotification = list(id = "da_running")
    ),
    formulaString = "~ 0 + group",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    daData = module_ci_lipid_da_data_env(),
    workflowData = invalid_workflow,
    session = "module-ci-lipid-da",
    experimentPaths = list(),
    removeNotificationFn = capture$removeNotification,
    notify = capture$showNotification,
    executeRunAnalysisFn = function(...) {
      shell_called <<- TRUE
      stop("analysis shell must not run for invalid preflight")
    }
  )
  expect_identical(analysis_state$status, "error")
  expect_false(shell_called)
  expect_match(analysis_state$errorMessage, "references terms absent", fixed = TRUE)
  expect_true("da_running" %in% capture$removed)
  expect_match(capture$notifications[[1L]][[1L]], "references terms absent", fixed = TRUE)
})

test_that("MCI-022.3 DA orchestration preserves schema, provenance, edge cases, and intensity columns", {
  object <- module_ci_lipid_da_object(layout = "triple")
  contrasts <- module_ci_lipid_da_contrasts("multi_group")
  observed <- module_ci_lipid_da_bind_orchestration_stubs("mixed")

  suppressMessages(capture.output(
    result <- runLipidsDA(
      theObject = object,
      contrasts_tbl = contrasts,
      formula_string = "~ 0 + group",
      da_q_val_thresh = 0.05,
      treat_lfc_cutoff = 0
    )
  ))

  module_ci_lipid_da_assert_long_schema(result$da_lipids_long)
  expect_named(result$per_assay_results, c("LCMS_Pos", "LCMS_Neg", "GCMS"))
  expect_identical(sort(unique(result$da_lipids_long$assay)), c("GCMS", "LCMS_Neg", "LCMS_Pos"))
  expect_identical(sort(unique(result$da_lipids_long$friendly_name)), c("KO_vs_WT", "RES_vs_WT"))
  expect_identical(sort(unique(result$da_lipids_long$comparison)), sort(contrasts$contrasts))
  expect_true(all(module_ci_lipid_da_feature_ids() %in% result$da_lipids_long$lipid_id))
  expect_true(all(c("intensity.S1.WT", "intensity.S3.KO", "intensity.S5.RES") %in% names(result$da_lipids_long)))
  expect_gt(result$significant_counts$LCMS_Pos$up, 0)
  expect_gt(result$significant_counts$LCMS_Pos$down, 0)
  expect_gt(result$significant_counts$LCMS_Pos$ns, 0)
  expect_identical(observed$calls[[1L]]$sample_cols, module_ci_lipid_da_samples())

  cases <- list(
    no_significant = function(res) all(res$da_lipids_long$significant == "NS"),
    all_significant = function(res) all(res$da_lipids_long$significant != "NS"),
    no_variance = function(res) all(res$da_lipids_long$logFC == 0) &&
      all(res$da_lipids_long$fdr_qvalue == 1),
    tied_pvalues = function(res) all(res$da_lipids_long$fdr_qvalue == 0.02)
  )
  for (case_name in names(cases)) {
    module_ci_lipid_da_bind_orchestration_stubs(case_name)
    suppressMessages(capture.output(
      case_result <- runLipidsDA(
        theObject = object,
        contrasts_tbl = contrasts,
        formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0
      )
    ))
    expect_true(cases[[case_name]](case_result), info = case_name)
  }

  module_ci_lipid_da_bind_orchestration_stubs("mixed")
  missing_metadata_object <- module_ci_lipid_da_object(
    layout = "single_pos",
    missing_metadata = TRUE
  )
  suppressMessages(capture.output(
    missing_metadata <- runLipidsDA(
      theObject = missing_metadata_object,
      contrasts_tbl = module_ci_lipid_da_contrasts("two_group"),
      formula_string = "~ 0 + group"
    )
  ))
  expect_identical(
    missing_metadata$da_lipids_long$lipid_name,
    missing_metadata$da_lipids_long$lipid_id
  )
})

test_that("MCI-022.4 render matrix covers combined/per-assay filters, plots, tables, thresholds, and selectors", {
  results_list <- module_ci_lipid_da_results_list("mixed")
  long <- results_list$da_lipids_long

  combined <- long[long$comparison == "groupKO-groupWT", , drop = FALSE]
  expect_setequal(unique(combined$assay), c("LCMS_Pos", "LCMS_Neg", "GCMS"))

  per_assay <- buildLipidDaResultsTableWidget(
    daResults = long,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "LCMS_Pos",
    tableSignificance = "up",
    daQValThresh = 0.05,
    lfcThreshold = 1,
    tableMaxRows = 10,
    datatableFactory = function(data, ...) data,
    formatRoundFactory = function(widget, ...) widget,
    formatStyleFactory = function(widget, ...) widget,
    styleEqualFactory = function(...) NULL
  )
  expect_s3_class(per_assay, "data.frame")
  expect_true(nrow(per_assay) > 0)
  expect_identical(unique(per_assay$assay), "LCMS_Pos")
  expect_true(all(per_assay$logFC > 1))

  empty <- buildLipidDaResultsTableWidget(
    daResults = long,
    selectedContrast = "DOES_NOT_EXIST",
    selectedAssay = "LCMS_Pos",
    tableSignificance = "all",
    daQValThresh = 0.05,
    lfcThreshold = 1,
    tableMaxRows = 10
  )
  expect_null(empty)

  captured_volcano <- NULL
  volcano <- buildLipidDaVolcanoStaticPlot(
    daResultsList = results_list,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "Combined",
    daQValThresh = 0.05,
    lfcThreshold = 1,
    plotFactory = function(...) {
      captured_volcano <<- list(...)
      "volcano-rendered"
    }
  )
  expect_identical(volcano, "volcano-rendered")
  expect_identical(captured_volcano$selected_assay, "Combined")
  expect_identical(captured_volcano$lfc_threshold, 1)

  da_data <- module_ci_lipid_da_data_env()
  captured_heatmap <- NULL
  heatmap_state <- buildLipidDaHeatmapRenderState(
    daResultsList = results_list,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "LCMS_Pos",
    topN = 4,
    heatmapClustering = "row",
    scaleData = "row",
    clusteringMethod = "ward.D2",
    distanceMethod = "euclidean",
    colorScheme = "RdBu",
    showLipidNames = TRUE,
    daQValThresh = 0.05,
    treeCutMethod = "kmeans",
    nClusters = 2,
    cutHeight = 0.5,
    minClusterSize = 2,
    heatmapFactory = function(...) {
      captured_heatmap <<- list(...)
      list(plot = "heatmap-rendered", row_clusters = c(L_UP = 1, L_DOWN = 2), col_clusters = NULL)
    }
  )
  expect_identical(heatmap_state$plot, "heatmap-rendered")
  expect_false(captured_heatmap$cluster_cols)
  expect_true(captured_heatmap$cluster_rows)
  expect_identical(captured_heatmap$top_n, 4)
  expect_identical(captured_heatmap$n_clusters, 2)
  expect_identical(storeLipidDaHeatmapRenderState(heatmap_state, da_data), "heatmap-rendered")
  expect_identical(da_data$current_row_clusters, c(L_UP = 1, L_DOWN = 2))

  summary_text <- buildLipidDaSummaryStatsText(
    daResults = long,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "LCMS_Pos",
    daQValThresh = 0.05
  )
  expect_match(summary_text, "Total lipids:", fixed = TRUE)
  expect_match(summary_text, "Up-regulated", fixed = TRUE)

  capture <- module_ci_lipid_da_capture_ui()
  workflow <- module_ci_lipid_da_workflow(module_ci_lipid_da_object(layout = "triple"))
  da_env <- module_ci_lipid_da_data_env()
  selector_state <- finalizeLipidDaRunAnalysisSuccess(
    results = results_list,
    daData = da_env,
    workflowData = workflow,
    session = "module-ci-lipid-da",
    experimentPaths = list(da_output_dir = NULL, publication_graphs_dir = NULL),
    daQValThresh = 0.05,
    treatLfcCutoff = 1,
    removeNotificationFn = capture$removeNotification,
    notify = capture$showNotification,
    updateSelectInputFn = capture$updateSelectInput,
    logInfo = function(...) invisible(NULL),
    logWarn = function(...) invisible(NULL),
    writeResultsFn = function(...) TRUE
  )
  expect_identical(selector_state$contrastChoices, c("KO_vs_WT", "RES_vs_WT"))
  expect_identical(selector_state$assayChoices, c("Combined", "LCMS_Pos", "LCMS_Neg", "GCMS"))
  expect_identical(selector_state$tableAssayChoices, c("All", "LCMS_Pos", "LCMS_Neg", "GCMS"))
  expect_identical(workflow$tab_status$differential_analysis, "complete")
})

test_that("MCI-022.5 export matrix preserves TSV/XLSX schema, assay provenance, intensity columns, and prefixes", {
  tmp_root <- tempfile("mci-022-export-")
  da_output_dir <- file.path(tmp_root, "da_lipids")
  publication_graphs_dir <- file.path(tmp_root, "Publication")
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  module_ci_lipid_da_local_binding(
    package_ns,
    "generateLipidDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  module_ci_lipid_da_local_binding(
    package_ns,
    "generateLipidDAHeatmap",
    function(...) NULL
  )

  expect_true(outputLipidDaResultsAllContrasts(
    da_results_list = module_ci_lipid_da_results_list("mixed"),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  pos_tsv <- file.path(da_output_dir, "de_posmode_lipids_groupKO-groupWT_long_annot.tsv")
  neg_tsv <- file.path(da_output_dir, "de_negmode_lipids_groupKO-groupWT_long_annot.tsv")
  gc_tsv <- file.path(da_output_dir, "de_gcms_lipids_groupKO-groupWT_long_annot.tsv")
  expect_true(file.exists(pos_tsv))
  expect_true(file.exists(neg_tsv))
  expect_true(file.exists(gc_tsv))
  expect_true(file.exists(file.path(da_output_dir, "de_posmode_lipids_groupKO-groupWT_long_annot.xlsx")))
  expect_true(file.exists(file.path(da_output_dir, "de_negmode_lipids_groupKO-groupWT_long_annot.xlsx")))

  exported <- read.delim(pos_tsv, check.names = FALSE)
  expect_equal(
    names(exported)[seq_len(7)],
    c("lipid_id", "lipid_name", "assay", "logFC", "raw_pvalue", "fdr_qvalue", "significant")
  )
  expect_true(all(exported$assay == "LCMS_Pos"))
  expect_true(all(c("comparison", "friendly_name", "numerator", "denominator") %in% names(exported)))
  expect_true(any(grepl("^intensity\\.", names(exported))))
  expect_identical(unique(exported$comparison), "groupKO-groupWT")
  expect_identical(unique(exported$friendly_name), "KO_vs_WT")

  numsig <- read.delim(
    file.path(publication_graphs_dir, "NumSigDeMolecules", "lipids_num_sig_de_molecules.tab"),
    check.names = FALSE
  )
  expect_setequal(numsig$mode, c("posmode", "negmode", "gcms"))
  expect_setequal(numsig$assay, c("LCMS_Pos", "LCMS_Neg", "GCMS"))
})

test_that("MCI-022.6 summary/report consumer smoke accepts exported lipidomics DA schemas", {
  tmp_root <- tempfile("mci-022-summary-smoke-")
  da_output_dir <- file.path(tmp_root, "da_lipids")
  publication_graphs_dir <- file.path(tmp_root, "Publication")
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  module_ci_lipid_da_local_binding(
    package_ns,
    "generateLipidDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  module_ci_lipid_da_local_binding(
    package_ns,
    "generateLipidDAHeatmap",
    function(...) NULL
  )

  expect_true(outputLipidDaResultsAllContrasts(
    da_results_list = module_ci_lipid_da_results_list("mixed"),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  tsv_paths <- file.path(
    da_output_dir,
    c(
      "de_posmode_lipids_groupKO-groupWT_long_annot.tsv",
      "de_negmode_lipids_groupKO-groupWT_long_annot.tsv",
      "de_gcms_lipids_groupKO-groupWT_long_annot.tsv"
    )
  )
  smoke <- module_ci_lipid_da_report_schema_smoke(tsv_paths)
  expect_identical(smoke$table_count, 3L)
  expect_setequal(smoke$assays, c("GCMS", "LCMS_Neg", "LCMS_Pos"))
  expect_true(length(smoke$intensity_columns) > 0L)
})
