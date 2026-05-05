test_that("MCI-016.1 reload matrix covers valid LC/GC/combined sessions and failure branches", {
  for (layout in c("lc", "gc", "combined")) {
    tmp_root <- tempfile(paste0("mci-016-reload-", layout, "-"))
    dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
    withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

    session_file <- file.path(tmp_root, "metab_filtered_session_data_latest.rds")
    session_data <- module_ci_metab_da_session_data(layout = layout)
    saveRDS(session_data, session_file)

    resolved <- resolveMetabDaLoadSessionFile(
      experimentPaths = list(source_dir = tmp_root, export_dir = tempfile("unused-"))
    )
    expect_true(resolved$ok)
    expect_identical(resolved$sessionFile, session_file)

    workflow <- module_ci_metab_da_workflow(session_data$current_s4_object)
    da_data <- module_ci_metab_da_data_env()
    capture <- module_ci_metab_da_capture_ui()

    load_state <- runMetabDaLoadSessionObserverShell(
      sessionFile = session_file,
      workflowData = workflow,
      daData = da_data,
      session = "module-ci-metab-da",
      restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
        restoreMetabDaLoadedSessionState(
          sessionData = sessionData,
          sessionFile = sessionFile,
          workflowData = workflowData,
          daData = daData,
          session = session,
          updateSelectInput = capture$updateSelectInput,
          updateTextAreaInput = capture$updateTextAreaInput,
          logInfo = function(...) invisible(NULL),
          logWarn = function(...) invisible(NULL),
          debugLog = debugLog
        )
      },
      showNotification = capture$showNotification,
      removeNotification = capture$removeNotification,
      logInfo = function(...) invisible(NULL),
      logError = function(...) invisible(NULL)
    )

    expect_identical(load_state$status, "success")
    expect_s4_class(da_data$current_s4_object, "MetaboliteAssayData")
    expect_identical(names(da_data$current_s4_object@metabolite_data), session_data$assay_names)
    expect_identical(da_data$assays_available, session_data$assay_names)
    expect_identical(workflow$contrasts_tbl$contrasts, session_data$contrasts_tbl$contrasts)
    expect_identical(da_data$formula_from_s4, "~ 0 + group")
    expect_true("loading_session" %in% capture$removed)
    expect_true(any(vapply(capture$select_updates, function(update) {
      identical(update$inputId, "volcano_assay") &&
        identical(update$choices, c("Combined", session_data$assay_names))
    }, logical(1))))
  }

  stale_export_dir <- tempfile("mci-016-stale-export-")
  dir.create(stale_export_dir, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(stale_export_dir, recursive = TRUE, force = TRUE))
  stale_file <- file.path(stale_export_dir, "metab_filtered_session_data_latest.rds")
  saveRDS(module_ci_metab_da_session_data("combined"), stale_file)
  stale_resolved <- resolveMetabDaLoadSessionFile(
    experimentPaths = list(source_dir = file.path(stale_export_dir, "missing"), export_dir = stale_export_dir)
  )
  expect_true(stale_resolved$ok)
  expect_identical(stale_resolved$sessionFile, stale_file)

  missing_resolved <- resolveMetabDaLoadSessionFile(
    experimentPaths = list(source_dir = stale_export_dir, export_dir = stale_export_dir),
    sessionFilename = "does-not-exist.rds"
  )
  expect_false(missing_resolved$ok)
  expect_match(missing_resolved$errorMessage, "Session file not found", fixed = TRUE)

  malformed_file <- file.path(stale_export_dir, "malformed.rds")
  writeLines("not an rds", malformed_file)
  malformed_capture <- module_ci_metab_da_capture_ui()
  malformed <- runMetabDaLoadSessionObserverShell(
    sessionFile = malformed_file,
    workflowData = module_ci_metab_da_workflow(),
    daData = module_ci_metab_da_data_env(),
    session = "module-ci-metab-da",
    showNotification = malformed_capture$showNotification,
    removeNotification = malformed_capture$removeNotification,
    logInfo = function(...) invisible(NULL),
    logError = function(...) invisible(NULL)
  )
  expect_identical(malformed$status, "error")
  expect_true(nzchar(malformed$errorMessage))

  mismatch_file <- file.path(stale_export_dir, "mismatch.rds")
  mismatch_session <- module_ci_metab_da_session_data(
    layout = "combined",
    assay_names = c("LCMS_Pos", "STALE_ASSAY")
  )
  saveRDS(mismatch_session, mismatch_file)
  mismatch_capture <- module_ci_metab_da_capture_ui()
  mismatch <- runMetabDaLoadSessionObserverShell(
    sessionFile = mismatch_file,
    workflowData = module_ci_metab_da_workflow(mismatch_session$current_s4_object),
    daData = module_ci_metab_da_data_env(),
    session = "module-ci-metab-da",
    restoreState = function(sessionData, sessionFile, workflowData, daData, session, debugLog) {
      restoreMetabDaLoadedSessionState(
        sessionData = sessionData,
        sessionFile = sessionFile,
        workflowData = workflowData,
        daData = daData,
        session = session,
        updateSelectInput = mismatch_capture$updateSelectInput,
        updateTextAreaInput = mismatch_capture$updateTextAreaInput,
        logInfo = function(...) invisible(NULL),
        logWarn = function(...) invisible(NULL),
        debugLog = debugLog
      )
    },
    showNotification = mismatch_capture$showNotification,
    removeNotification = mismatch_capture$removeNotification,
    logInfo = function(...) invisible(NULL),
    logError = function(...) invisible(NULL)
  )
  expect_identical(mismatch$status, "error")
  expect_match(mismatch$errorMessage, "assay_names do not match", fixed = TRUE)
})

test_that("MCI-016.2 formula and contrast preflight rejects invalid DA inputs before analysis", {
  object <- module_ci_metab_da_object(layout = "combined")
  valid_cases <- list(
    two_group = list(
      contrasts = module_ci_metab_da_contrasts("two_group"),
      formula = "~ 0 + group"
    ),
    multi_group = list(
      contrasts = module_ci_metab_da_contrasts("multi_group"),
      formula = "~ 0 + group"
    ),
    batch_aware = list(
      contrasts = module_ci_metab_da_contrasts("batch_aware"),
      formula = "~ 0 + group + batch"
    ),
    reversed = list(
      contrasts = module_ci_metab_da_contrasts("reversed"),
      formula = "~ 0 + group"
    )
  )

  for (case in valid_cases) {
    preflight <- validateMetabDesignDaPreflight(
      designMatrix = object@design_matrix,
      assayList = object@metabolite_data,
      contrastsTbl = case$contrasts,
      formulaString = case$formula,
      sampleIdCol = object@sample_id,
      groupCol = object@group_id
    )
    expect_true(preflight$valid)
    expect_true(all(case$contrasts$contrasts %in% preflight$contrast_validation$contrasts_tbl$contrasts))
  }

  duplicate <- validateMetabDesignDaPreflight(
    designMatrix = object@design_matrix,
    assayList = object@metabolite_data,
    contrastsTbl = module_ci_metab_da_contrasts("duplicate"),
    formulaString = "~ 0 + group",
    sampleIdCol = object@sample_id,
    groupCol = object@group_id
  )
  expect_false(duplicate$valid)
  expect_true(any(grepl("duplicate friendly contrast labels", duplicate$errors, fixed = TRUE)))

  invalid_term <- validateMetabDesignDaPreflight(
    designMatrix = object@design_matrix,
    assayList = object@metabolite_data,
    contrastsTbl = module_ci_metab_da_contrasts("invalid_term"),
    formulaString = "~ 0 + group",
    sampleIdCol = object@sample_id,
    groupCol = object@group_id
  )
  expect_false(invalid_term$valid)
  expect_true(any(grepl("references terms absent", invalid_term$errors, fixed = TRUE)))

  invalid_formula <- validateMetabDesignDaPreflight(
    designMatrix = object@design_matrix,
    assayList = object@metabolite_data,
    contrastsTbl = module_ci_metab_da_contrasts("two_group"),
    formulaString = "~ 0 + missing_factor",
    sampleIdCol = object@sample_id,
    groupCol = object@group_id
  )
  expect_false(invalid_formula$valid)
  expect_true(any(grepl("cannot produce a model matrix", invalid_formula$errors, fixed = TRUE)))

  no_contrasts <- resolveMetabDaAnalysisInputs(
    currentS4Object = object,
    workflowData = list(contrasts_tbl = module_ci_metab_da_contrasts("no_contrasts")),
    formulaString = "~ 0 + group"
  )
  expect_false(no_contrasts$ok)
  expect_identical(
    no_contrasts$errorMessage,
    "No contrasts defined. Please define contrasts in the design tab."
  )

  invalid_workflow <- module_ci_metab_da_workflow(
    object,
    contrasts = module_ci_metab_da_contrasts("invalid_term")
  )
  capture <- module_ci_metab_da_capture_ui()
  shell_called <- FALSE
  analysis_state <- runMetabDaAnalysisObserverEntry(
    currentS4Object = object,
    workflowData = invalid_workflow,
    formulaString = "~ 0 + group",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    daData = module_ci_metab_da_data_env(),
    session = "module-ci-metab-da",
    experimentPaths = list(),
    runAnalysisShell = function(...) {
      shell_called <<- TRUE
      stop("analysis shell must not run for invalid preflight")
    },
    showNotification = capture$showNotification
  )
  expect_identical(analysis_state$status, "error")
  expect_false(shell_called)
  expect_match(analysis_state$errorMessage, "references terms absent", fixed = TRUE)
  expect_match(capture$notifications[[1L]][[1L]], "references terms absent", fixed = TRUE)
})

test_that("MCI-016.3 DA orchestration preserves schema, provenance, edge cases, and intensity columns", {
  object <- module_ci_metab_da_object(layout = "combined")
  contrasts <- module_ci_metab_da_contrasts("multi_group")
  observed <- module_ci_metab_da_bind_orchestration_stubs("mixed")

  suppressMessages(capture.output(
    result <- runMetabolitesDA(
      theObject = object,
      contrasts_tbl = contrasts,
      formula_string = "~ 0 + group",
      da_q_val_thresh = 0.05,
      treat_lfc_cutoff = 0
    )
  ))

  module_ci_metab_da_assert_long_schema(result$da_metabolites_long)
  expect_named(result$per_assay_results, c("LCMS_Pos", "GCMS"))
  expect_identical(sort(unique(result$da_metabolites_long$assay)), c("GCMS", "LCMS_Pos"))
  expect_identical(sort(unique(result$da_metabolites_long$friendly_name)), c("KO_vs_WT", "RES_vs_WT"))
  expect_identical(sort(unique(result$da_metabolites_long$comparison)), sort(contrasts$contrasts))
  expect_true(all(module_ci_metab_da_feature_ids() %in% result$da_metabolites_long$metabolite_id))
  expect_true(all(c("intensity.S1.WT", "intensity.S3.KO", "intensity.S5.RES") %in% names(result$da_metabolites_long)))
  expect_gt(result$significant_counts$LCMS_Pos$up, 0)
  expect_gt(result$significant_counts$LCMS_Pos$down, 0)
  expect_gt(result$significant_counts$LCMS_Pos$ns, 0)
  expect_identical(observed$calls[[1L]]$sample_cols, module_ci_metab_da_samples())

  cases <- list(
    no_significant = function(res) all(res$da_metabolites_long$significant == "NS"),
    all_significant = function(res) all(res$da_metabolites_long$significant != "NS"),
    no_variance = function(res) all(res$da_metabolites_long$logFC == 0) &&
      all(res$da_metabolites_long$fdr_qvalue == 1),
    tied_pvalues = function(res) all(res$da_metabolites_long$fdr_qvalue == 0.02)
  )
  for (case_name in names(cases)) {
    module_ci_metab_da_bind_orchestration_stubs(case_name)
    suppressMessages(capture.output(
      case_result <- runMetabolitesDA(
        theObject = object,
        contrasts_tbl = contrasts,
        formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0
      )
    ))
    expect_true(cases[[case_name]](case_result), info = case_name)
  }

  module_ci_metab_da_bind_orchestration_stubs("mixed")
  missing_metadata_object <- module_ci_metab_da_object(
    layout = "lc",
    missing_metadata = TRUE
  )
  suppressMessages(capture.output(
    missing_metadata <- runMetabolitesDA(
      theObject = missing_metadata_object,
      contrasts_tbl = module_ci_metab_da_contrasts("two_group"),
      formula_string = "~ 0 + group"
    )
  ))
  expect_identical(
    missing_metadata$da_metabolites_long$metabolite_name,
    missing_metadata$da_metabolites_long$metabolite_id
  )
})

test_that("MCI-016.4 render matrix covers combined/per-assay filters, plots, tables, thresholds, and selectors", {
  results_list <- module_ci_metab_da_results_list("mixed")
  long <- results_list$da_metabolites_long

  combined <- filterMetabDaDisplayResults(
    results = long,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "Combined"
  )
  expect_setequal(unique(combined$assay), c("LCMS_Pos", "GCMS"))

  per_assay <- filterMetabDaDisplayResults(
    results = long,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "LCMS_Pos",
    significanceFilter = "up",
    daQValThresh = 0.05,
    treatLfcCutoff = 1
  )
  expect_true(nrow(per_assay) > 0)
  expect_identical(unique(per_assay$assay), "LCMS_Pos")
  expect_true(all(per_assay$logFC > 1))
  expect_true(all(per_assay$fdr_qvalue < 0.05))

  empty <- filterMetabDaDisplayResults(
    results = long,
    selectedContrast = "DOES_NOT_EXIST",
    selectedAssay = "LCMS_Pos"
  )
  expect_equal(nrow(empty), 0)
  expect_null(buildMetabDaResultsTableRenderOutput(
    daResultsList = list(da_metabolites_long = empty),
    selectedContrast = "DOES_NOT_EXIST",
    selectedAssay = "LCMS_Pos",
    requireInputs = function(...) invisible(NULL)
  ))

  captured_volcano <- NULL
  volcano <- buildMetabDaStaticVolcanoRenderOutput(
    daResultsList = results_list,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "Combined",
    daQValThresh = 0.05,
    treatLfcCutoff = 1,
    requireInputs = function(...) invisible(NULL),
    buildPlot = function(...) {
      captured_volcano <<- list(...)
      "volcano-rendered"
    }
  )
  expect_identical(volcano, "volcano-rendered")
  expect_identical(captured_volcano$selectedAssay, "Combined")
  expect_identical(captured_volcano$treatLfcCutoff, 1)

  captured_heatmap <- NULL
  da_data <- module_ci_metab_da_data_env()
  heatmap <- buildMetabDaHeatmapRenderOutput(
    daResultsList = results_list,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "LCMS_Pos",
    topN = 4,
    heatmapClustering = "row",
    treeCutMethod = "kmeans",
    nClusters = 2,
    daData = da_data,
    requireInputs = function(...) invisible(NULL),
    buildPlot = function(...) {
      captured_heatmap <<- list(...)
      list(plot = "heatmap-rendered", row_clusters = c(M_UP = 1, M_DOWN = 2))
    }
  )
  expect_type(heatmap, "list")
  expect_false(captured_heatmap$clusterCols)
  expect_true(captured_heatmap$clusterRows)
  expect_identical(captured_heatmap$topN, 4)
  expect_identical(captured_heatmap$nClusters, 2)

  table_output <- buildMetabDaResultsTableRenderOutput(
    daResultsList = results_list,
    selectedContrast = "KO_vs_WT",
    selectedAssay = "All",
    significanceFilter = "significant",
    daQValThresh = 0.05,
    maxRows = 5,
    requireInputs = function(...) invisible(NULL),
    buildDatatable = function(results) results
  )
  expect_s3_class(table_output, "data.frame")
  expect_lte(nrow(table_output), 5)
  expect_true(all(table_output$significant != "NS"))

  capture <- module_ci_metab_da_capture_ui()
  selector_state <- updateMetabDaResultsSelectorInputs(
    daResultsLong = long,
    session = "module-ci-metab-da",
    updateSelectInput = capture$updateSelectInput,
    logInfo = function(...) invisible(NULL)
  )
  expect_identical(selector_state$contrastChoices, c("KO_vs_WT", "RES_vs_WT"))
  expect_identical(selector_state$assayChoices, c("Combined", "LCMS_Pos", "GCMS"))
  expect_identical(selector_state$tableAssayChoices, c("All", "LCMS_Pos", "GCMS"))
  expect_identical(capture$select_updates[[4L]]$selected, "Combined")
  expect_identical(capture$select_updates[[6L]]$selected, "All")
})

test_that("MCI-016.5 export matrix preserves TSV/XLSX schema, assay provenance, intensity columns, and prefixes", {
  tmp_root <- tempfile("mci-016-export-")
  da_output_dir <- file.path(tmp_root, "da_metabolites")
  publication_graphs_dir <- file.path(tmp_root, "Publication")
  dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  module_ci_metab_da_local_binding(
    package_ns,
    "generateMetabDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  module_ci_metab_da_local_binding(
    package_ns,
    "generateMetabDAHeatmap",
    function(...) NULL
  )

  expect_true(outputMetabDaResultsAllContrasts(
    da_results_list = module_ci_metab_da_results_list("mixed"),
    da_output_dir = da_output_dir,
    publication_graphs_dir = publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  pos_tsv <- file.path(da_output_dir, "de_posmode_metabolites_groupKO-groupWT_long_annot.tsv")
  gc_tsv <- file.path(da_output_dir, "de_gcms_metabolites_groupKO-groupWT_long_annot.tsv")
  expect_true(file.exists(pos_tsv))
  expect_true(file.exists(gc_tsv))
  expect_true(file.exists(file.path(da_output_dir, "de_posmode_metabolites_groupKO-groupWT_long_annot.xlsx")))
  expect_true(file.exists(file.path(da_output_dir, "de_gcms_metabolites_groupKO-groupWT_long_annot.xlsx")))

  exported <- read.delim(pos_tsv, check.names = FALSE)
  expect_equal(
    names(exported)[seq_len(6)],
    c("metabolite_id", "metabolite_name", "logFC", "raw_pvalue", "fdr_qvalue", "significant")
  )
  expect_true("assay" %in% names(exported))
  expect_true(all(exported$assay == "LCMS_Pos"))
  expect_true(all(c("comparison", "friendly_name", "numerator", "denominator") %in% names(exported)))
  expect_true(any(grepl("^intensity\\.", names(exported))))
  expect_identical(unique(exported$comparison), "groupKO-groupWT")
  expect_identical(unique(exported$friendly_name), "KO_vs_WT")

  numsig <- read.delim(
    file.path(publication_graphs_dir, "NumSigDeMolecules", "metabolites_num_sig_de_molecules.tab"),
    check.names = FALSE
  )
  expect_setequal(numsig$mode, c("posmode", "gcms"))
  expect_setequal(numsig$assay, c("LCMS_Pos", "GCMS"))
})

test_that("MCI-016.6 summary/report consumer smoke accepts exported metabolomics DA schemas", {
  tmp_root <- tempfile("mci-016-summary-smoke-")
  paths <- list(
    results_dir = file.path(tmp_root, "results"),
    results_summary_dir = file.path(tmp_root, "results_summary"),
    publication_graphs_dir = file.path(tmp_root, "Publication"),
    time_dir = file.path(tmp_root, "time"),
    qc_dir = file.path(tmp_root, "qc"),
    da_output_dir = file.path(tmp_root, "da_metabolites"),
    pathway_dir = file.path(tmp_root, "pathway"),
    source_dir = file.path(tmp_root, "source"),
    feature_qc_dir = file.path(tmp_root, "feature_qc")
  )
  lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_root, recursive = TRUE, force = TRUE))

  package_ns <- asNamespace("MultiScholaR")
  module_ci_metab_da_local_binding(
    package_ns,
    "generateMetabDAVolcanoStatic",
    function(...) {
      ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
        ggplot2::geom_point()
    }
  )
  module_ci_metab_da_local_binding(
    package_ns,
    "generateMetabDAHeatmap",
    function(...) NULL
  )

  expect_true(outputMetabDaResultsAllContrasts(
    da_results_list = module_ci_metab_da_results_list("mixed"),
    da_output_dir = paths$da_output_dir,
    publication_graphs_dir = paths$publication_graphs_dir,
    da_q_val_thresh = 0.05,
    lfc_threshold = 1
  ))

  tsv_paths <- file.path(
    paths$da_output_dir,
    c(
      "de_posmode_metabolites_groupKO-groupWT_long_annot.tsv",
      "de_gcms_metabolites_groupKO-groupWT_long_annot.tsv"
    )
  )
  smoke <- module_ci_metab_da_report_schema_smoke(tsv_paths)
  expect_identical(smoke$table_count, 2L)
  expect_setequal(smoke$assays, c("GCMS", "LCMS_Pos"))
  expect_true(length(smoke$intensity_columns) > 0L)

  writeLines("module-ci metabolomics study parameters", file.path(paths$source_dir, "study_parameters.txt"))
  write.table(
    module_ci_metab_da_object()@design_matrix,
    file.path(paths$source_dir, "design_matrix.tab"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write.table(
    module_ci_metab_da_contrasts("multi_group"),
    file.path(paths$source_dir, "contrasts_tbl.tab"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  old_project_dirs <- if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
    get("project_dirs", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  had_project_dirs <- exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)
  assign("project_dirs", list(metabolomics = paths), envir = .GlobalEnv)
  withr::defer({
    if (had_project_dirs) {
      assign("project_dirs", old_project_dirs, envir = .GlobalEnv)
    } else if (exists("project_dirs", envir = .GlobalEnv, inherits = FALSE)) {
      rm("project_dirs", envir = .GlobalEnv)
    }
  })

  suppressMessages(capture.output(
    copy_state <- copyToResultsSummary(
      omic_type = "metabolomics",
      experiment_label = "module_ci",
      contrasts_tbl = module_ci_metab_da_contrasts("multi_group"),
      design_matrix = module_ci_metab_da_object()@design_matrix,
      force = TRUE,
      current_rmd = NULL
    )
  ))

  expect_true(file.exists(file.path(
    paths$results_summary_dir,
    "Publication_tables",
    "DA_results_metabolomics.xlsx"
  )))
  expect_type(copy_state, "list")
})
