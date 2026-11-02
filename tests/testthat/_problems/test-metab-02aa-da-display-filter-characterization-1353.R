# Extracted from test-metab-02aa-da-display-filter-characterization.R:1353

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
omicsStandInReaders()
captured_event <- NULL
captured_entry <- NULL
captured_log <- NULL
workflow_state <- new.env(parent = emptyenv())
workflow_state$contrasts_tbl <- data.frame(
    friendly_names = c("B vs A"),
    contrasts = c("groupB-groupA"),
    stringsAsFactors = FALSE
  )
da_state <- new.env(parent = emptyenv())
da_state$current_state_name <- "loaded_for_de"
workflow_state$state_manager <- list(getState = function(state_name = NULL) {
    if (identical(state_name, "loaded_for_de")) "current-s4" else NULL
  })
registration <- registerMetabDaRunAnalysisObserver(
    input = list(
      run_da_analysis = 2,
      formula_string = "~ 0 + group",
      da_q_val_thresh = 0.01,
      treat_lfc_cutoff = 1.5
    ),
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    experimentPaths = list(
      da_output_dir = "/tmp/metab-da-output",
      publication_graphs_dir = "/tmp/metab-publication-graphs"
    ),
    observeEvent = function(eventExpr, handlerExpr) {
      captured_event <<- eventExpr
      eval(substitute(handlerExpr), envir = parent.frame())
      "observer-registered"
    },
    runAnalysisObserverEntry = function(...) {
      captured_entry <<- list(...)
      invisible(list(status = "success"))
    },
    logInfo = function(message) {
      captured_log <<- message
      invisible(NULL)
    }
  )
expect_identical(registration, "observer-registered")
expect_identical(captured_event, 2)
expect_identical(captured_log, "=== RUN DA ANALYSIS BUTTON CLICKED ===")
expect_identical(captured_entry$currentS4Object, "current-s4")
