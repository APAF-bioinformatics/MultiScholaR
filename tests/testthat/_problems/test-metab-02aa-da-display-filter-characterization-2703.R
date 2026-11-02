# Extracted from test-metab-02aa-da-display-filter-characterization.R:2703

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
captured_save_state <- NULL
captured_warn_logs <- character()
workflow_state <- new.env(parent = emptyenv())
omicsTestOwner(workflow = workflow_state, omic_type = "metabolomics")
adoption <- omicsMockWrapper("ingestWorkflowOmicsLegacySession", TRUE)
da_state <- new.env(parent = emptyenv())
session_data <- list(
    current_s4_object = list(kind = "non-s4"),
    r6_current_state_name = NULL,
    contrasts_tbl = NULL,
    assay_names = NULL
  )
restore_state <- restoreMetabDaLoadedSessionState(
    sessionData = session_data,
    sessionFile = "/tmp/metab_filtered_session_data_latest.rds",
    workflowData = workflow_state,
    daData = da_state,
    session = "metab-session",
    updateSelectInput = function(...) stop("should not update selectors"),
    updateTextAreaInput = function(...) stop("should not update formula"),
    logInfo = function(...) invisible(NULL),
    logWarn = function(message) {
      captured_warn_logs <<- c(captured_warn_logs, message)
      invisible(NULL)
    }
  )
expect_identical(restore_state$stateName, "loaded_for_de")
expect_identical(restore_state$contrastChoices, character())
expect_identical(restore_state$assayChoices, character())
expect_null(restore_state$formulaValue)
expect_identical(da_state$current_state_name, "loaded_for_de")
expect_false(exists("current_s4_object", envir = da_state, inherits = FALSE))
expect_identical(adoption$calls[[1L]]$state_name, "loaded_for_de")
expect_identical(adoption$calls[[1L]][[1L]]$current_s4_object, session_data$current_s4_object)
expect_length(captured_warn_logs, 1)
expect_match(captured_warn_logs[[1]], "^Could not extract formula from S4:")
