# Extracted from test-metab-02aa-da-display-filter-characterization.R:2153

# prequel ----------------------------------------------------------------------
library(testthat)
repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}
loadSelectedFunctions(
  paths = file.path(repo_root, "R", c(
    "mod_metab_da_display_helpers.R",
    "mod_metab_da_observer_helpers.R",
    "mod_metab_da_registration_helpers.R",
    "mod_metab_da_server_helpers.R",
    "mod_metab_da_server.R",
    "mod_metab_da.R"
  )),
  symbols = c(
    "filterMetabDaDisplayResults",
    "summarizeMetabDaDisplayResults",
    "buildMetabDaSummaryStatsText",
    "buildMetabDaSummaryStatsRenderOutput",
    "shapeMetabDaTableDisplayResults",
    "buildMetabDaResultsDatatable",
    "buildMetabDaResultsTableRenderOutput",
    "buildMetabDaClusterSummaryText",
    "buildMetabDaClusterSummaryRenderOutput",
    "buildMetabDaHeatmapManualSaveWarning",
    "buildMetabDaHeatmapManualSaveWarningRenderOutput",
    "buildMetabDaContrastsDisplayText",
    "buildMetabDaContrastsRenderOutput",
    "buildMetabDaStatusText",
    "buildMetabDaStatusRenderOutput",
    "registerMetabDaOverviewOutputs",
    "buildMetabDaGlimmaCombinedViewInfoBanner",
    "resolveMetabDaGlimmaWidgetOutput",
    "buildMetabDaGlimmaErrorBanner",
    "buildMetabDaGlimmaRenderOutput",
    "buildMetabDaStaticVolcanoPlot",
    "buildMetabDaStaticVolcanoRenderOutput",
    "buildMetabDaHeatmapPlotOutput",
    "buildMetabDaHeatmapRenderOutput",
    "registerMetabDaVisualizationOutputs",
    "registerMetabDaResultsOutputs",
    "registerMetabDaSaveHeatmapObserver",
    "registerMetabDaLoadFilteredSessionObserver",
    "registerMetabDaRunAnalysisObserver",
    "registerMetabDaMainOutputsAndObservers",
    "createMetabDaServerState",
    "registerMetabDaServerBodyOutputs",
    "startMetabDaServerBody",
    "initializeMetabDaServerBody",
    "createMetabDaServerModuleHandler",
    "runMetabDaServerModuleCallback",
    "runMetabDaServerModuleShell",
    "runMetabDaServerModule",
    "runMetabDaServerEntry",
    "mod_metab_da_server",
    "runMetabDaSaveHeatmapObserverShell",
    "runMetabDaSaveHeatmapObserverEntry",
    "updateMetabDaResultsSelectorInputs",
    "writeMetabDaResultsArtifacts",
    "writeMetabDaResultsDownloadCsv",
    "buildMetabDaResultsDownloadHandler",
    "resolveMetabDaAnalysisInputs",
    "runMetabDaAnalysisObserverEntry",
    "runMetabDaAnalysisObserverShell",
    "resolveMetabDaLoadSessionFile",
    "runMetabDaLoadSessionObserverEntry",
    "runMetabDaLoadSessionObserverShell",
    "restoreMetabDaLoadedSessionState"
  ),
  env = environment()
)
helper_env <- environment()
if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass("MetaboliteAssayData", slots = c(args = "list"))
}
makeMetabDaObserverS4 <- function(args = list()) {
  if (all(c(
    "metabolite_data", "design_matrix", "sample_id", "group_id",
    "metabolite_id_column", "annotation_id_column"
  ) %in% methods::slotNames("MetaboliteAssayData"))) {
    return(methods::new(
      "MetaboliteAssayData",
      metabolite_data = list(
        LCMS_Pos = data.frame(
          database_identifier = "M1",
          Name = "Met One",
          S1 = 10,
          stringsAsFactors = FALSE
        ),
        LCMS_Neg = data.frame(
          database_identifier = "M2",
          Name = "Met Two",
          S1 = 20,
          stringsAsFactors = FALSE
        )
      ),
      metabolite_id_column = "database_identifier",
      annotation_id_column = "Name",
      database_identifier_type = "database_identifier",
      internal_standard_regex = "^IS_",
      design_matrix = data.frame(
        Run = "S1",
        group = "A",
        replicates = "R1",
        stringsAsFactors = FALSE
      ),
      sample_id = "Run",
      group_id = "group",
      technical_replicate_id = NA_character_,
      args = args
    ))
  }

  methods::new("MetaboliteAssayData", args = args)
}
buildMetabDaDisplayResults <- function() {
  data.frame(
    metabolite_id = c("M1", "M2", "M3", "M4"),
    metabolite_name = c("Met One", "Met Two", "Met Three", "Met Four"),
    comparison = c("groupB-groupA", "groupB-groupA", "groupC-groupA", "groupB-groupA"),
    friendly_name = c("B vs A", "B vs A", "C vs A", "B vs A"),
    assay = c("LCMS_Pos", "LCMS_Neg", "LCMS_Pos", "LCMS_Pos"),
    logFC = c(1.2, -1.1, -0.9, 0.2),
    raw_pvalue = c(0.001, 0.004, 0.02, 0.3),
    fdr_qvalue = c(0.01, 0.02, 0.03, 0.20),
    significant = c("Up", "Down", "Down", "NS"),
    stringsAsFactors = FALSE
  )
}

# test -------------------------------------------------------------------------
captured_notifications <- list()
captured_removed <- character()
captured_run <- NULL
captured_selectors <- NULL
captured_artifacts <- NULL
captured_info_logs <- character()
captured_error_logs <- character()
workflow_state <- new.env(parent = emptyenv())
workflow_state$tab_status <- list(differential_analysis = "pending", other = "keep")
omicsTestOwner(workflow = workflow_state, omic_type = "metabolomics")
da_state <- new.env(parent = emptyenv())
current_s4 <- makeMetabDaObserverS4()
contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )
results <- list(
    da_metabolites_long = buildMetabDaDisplayResults(),
    significant_counts = list(LCMS_Pos = list(up = 1, down = 2, ns = 3))
  )
wrapper <- omicsMockWrapper("runMetabReferenceDaAnalysis", list(results = results))
analysis_state <- runMetabDaAnalysisObserverShell(
    currentS4 = current_s4,
    contrastsTbl = contrasts_tbl,
    formulaString = "~ 0 + group",
    daQValThresh = 0.01,
    treatLfcCutoff = 1.5,
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    updateSelectors = function(...) {
      captured_selectors <<- list(...)
      invisible(list(updated = TRUE))
    },
    writeArtifacts = function(...) {
      captured_artifacts <<- list(...)
      invisible(list(status = "written"))
    },
    showNotification = function(...) {
      captured_notifications[[length(captured_notifications) + 1]] <<- list(...)
      invisible(NULL)
    },
    removeNotification = function(id) {
      captured_removed <<- c(captured_removed, id)
      invisible(NULL)
    },
    logInfo = function(message) {
      captured_info_logs <<- c(captured_info_logs, message)
      invisible(NULL)
    },
    logError = function(message) {
      captured_error_logs <<- c(captured_error_logs, message)
      invisible(NULL)
    }
  )
expect_identical(analysis_state$status, "success")
expect_identical(analysis_state$results, results)
expect_identical(analysis_state$selectorState, list(updated = TRUE))
expect_identical(analysis_state$diskState, list(status = "written"))
expect_identical(wrapper$calls[[1L]][[1L]], workflow_state$state_manager)
expect_identical(wrapper$calls[[1L]][[3L]], list(formula_string = "~ 0 + group", da_q_val_thresh = 0.01,
    treat_lfc_cutoff = 1.5))
expect_identical(wrapper$calls[[1L]][[4L]], contrasts_tbl)
expect_identical(da_state$da_results_list, results)
expect_true(isTRUE(da_state$analysis_complete))
expect_identical(
    workflow_state$tab_status,
    list(differential_analysis = "complete", other = "keep")
  )
expect_identical(captured_removed, "da_running")
expect_identical(captured_selectors$daResultsLong, results$da_metabolites_long)
expect_identical(captured_selectors$session, "metab-session")
expect_identical(captured_artifacts$results, results)
expect_identical(captured_artifacts$experimentPaths$da_output_dir, "/tmp/metab-da-output")
