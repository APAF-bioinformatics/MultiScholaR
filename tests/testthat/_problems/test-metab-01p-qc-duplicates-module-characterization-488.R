# Extracted from test-metab-01p-qc-duplicates-module-characterization.R:488

# prequel ----------------------------------------------------------------------
library(testthat)
repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)
skipIfMissingMetabQcDuplicatesTargetFiles <- function() {
  required_paths <- c(
    "R/mod_metab_qc_duplicates_ui_helpers.R",
    "R/mod_metab_qc_duplicates_render_helpers.R",
    "R/mod_metab_qc_duplicates_server_helpers.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only metab QC duplicate helper file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}
skipIfMissingMetabQcDuplicatesTargetFiles()
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
  paths = c(
    file.path(repo_root, "R", "mod_metab_qc_duplicates_ui_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_render_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_server_helpers.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_ui.R"),
    file.path(repo_root, "R", "mod_metab_qc_duplicates_server.R")
  ),
  symbols = c(
    "reportMetabDuplicateDetection",
    "reportMetabDuplicateDetectionError",
    "runMetabDuplicateDetectionObserverShell",
    "runMetabDuplicateResolutionObserver",
    "runMetabDuplicateResolutionObserverShell",
    "runMetabDuplicateRevertObserverShell",
    "prepareMetabDuplicateResolutionState",
    "applyMetabDuplicateResolutionState",
    "runMetabDuplicateResolutionWorkflow",
    "resolveMetabDuplicateAssayData",
    "buildMetabDuplicateResolutionSummary",
    "buildMetabDuplicateSummaryUi",
    "buildMetabDuplicateTablesUi",
    "registerMetabDuplicateTableRenderers",
    "detectMetabDuplicateFeatures",
    "revertMetabDuplicateResolution",
    "renderMetabDuplicateFilterPlot",
    "mod_metab_qc_duplicates_server"
  ),
  env = environment()
)
if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      design_matrix = data.frame(
        Sample_ID = character(),
        group = character(),
        stringsAsFactors = FALSE
      ),
      sample_id = "Sample_ID",
      group_id = "group"
    )
  )
}
newMetabDuplicateModuleObject <- function(metabolite_data, metabolite_id_column = "metabolite_id") {
  sample_ids <- unique(unlist(lapply(metabolite_data, function(assay) {
    names(assay)[vapply(assay, is.numeric, logical(1))]
  }), use.names = FALSE))
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = metabolite_data,
    metabolite_id_column = metabolite_id_column,
    design_matrix = data.frame(
      Sample_ID = sample_ids,
      group = rep("A", length(sample_ids)),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group"
  )
}

# test -------------------------------------------------------------------------
captured <- new.env(parent = emptyenv())
captured$resolution_stats_updates <- list()
captured$filter_plot_updates <- list()
resolved_assay_data <- list(
    Plasma = data.frame(
      metabolite_id = "m1",
      Sample1 = 2,
      stringsAsFactors = FALSE
    )
  )
resolved_s4 <- newMetabDuplicateModuleObject(resolved_assay_data)
expected_stats <- list(
    Plasma = list(original = 2, resolved = 1, removed = 1)
  )
workflow_data <- omicsTestOwner(omic_type = "metabolomics", config = list(config = "value"))
resolution <- omicsMockWrapper("runMetabReferenceDuplicateResolution",
    list(object = resolved_s4, value = list(stats = expected_stats)))
visible <- withVisible(
    runMetabDuplicateResolutionWorkflow(
      workflowData = workflow_data,
      omicType = "metabolomics",
      setResolutionStatsFn = function(value) {
        captured$resolution_stats_updates <- c(
          captured$resolution_stats_updates,
          list(value)
        )
        invisible(NULL)
      },
      setFilterPlotFn = function(value) {
        captured$filter_plot_updates <- c(
          captured$filter_plot_updates,
          list(value)
        )
        invisible(NULL)
      },
      applyResolutionStateFn = function(currentS4, statsList, workflowData, omicType, setResolutionStatsFn, setFilterPlotFn) {
        captured$apply_call <- list(
          metabolite_data = currentS4@metabolite_data,
          stats_list = statsList,
          workflow_data = workflowData,
          omic_type = omicType
        )
        setResolutionStatsFn(statsList)
        setFilterPlotFn("plot-token")
        list(
          stateName = "metab_duplicates_resolved",
          qcPlot = "plot-token"
        )
      },
      buildResolutionSummaryFn = function(statsList, stateName) {
        captured$summary_call <- list(
          stats_list = statsList,
          state_name = stateName
        )
        list(
          totalRemoved = 1,
          resultText = "summary-token"
        )
      }
    )
  )
