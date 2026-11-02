# Extracted from test-metab-03b-norm-module-helper-characterization.R:5501

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
  paths = c(
    file.path(repo_root, "R", "mod_metab_norm_session_helpers.R"),
    file.path(repo_root, "R", "mod_metab_norm_pipeline_helpers.R"),
    file.path(repo_root, "R", "mod_metab_norm_render_helpers.R"),
    file.path(repo_root, "R", "mod_metab_norm_binding_helpers.R"),
    file.path(repo_root, "R", "mod_metab_norm_entry_helpers.R")
  ),
  symbols = c(
    "updateMetabNormDesignDrivenChoices",
    "initializeMetabNormAssayNames",
    "runMetabNormAutoPreNormalizationQcObserverShell",
    "renderMetabNormItsdSelectionUi",
    "renderMetabNormItsdSelectionTable",
    "registerMetabNormItsdSelectionTracking",
    "runMetabNormItsdSelectionTrackingObserverShell",
    "runMetabNormItsdSelectionTableObserverShell",
    "renderMetabNormRuvQcUi",
    "renderMetabNormRuvCancorPlot",
    "renderMetabNormRuvOptimizationSummary",
    "renderMetabNormRuvResultsTable",
    "runMetabNormRuvBindingObserverShell",
    "runMetabNormQcImageBindingShell",
    "runMetabNormAssayLabelBindingShell",
    "runMetabNormNormalizationObserverWrapper",
    "runMetabNormResetNormalizationObserverWrapper",
    "runMetabNormApplyCorrelationObserverWrapper",
    "runMetabNormSkipCorrelationObserverWrapper",
    "runMetabNormExportSessionObserverWrapper",
    "buildMetabNormCorrelationFilterSummary",
    "resolveMetabNormFinalQcRenderState",
    "buildMetabNormFinalQcPcaPlot",
    "getPlotAesthetics",
    "appendMetabNormNormalizationLog",
    "renderMetabNormNormalizationLog",
    "renderMetabNormCorrelationFilterSummary",
    "renderMetabNormFinalQcPlot",
    "renderMetabNormAssayLabel",
    "renderMetabNormQcImageForAssay",
    "resolveMetabNormExportSourceDir",
    "collectMetabNormFeatureCountsPerAssay",
    "buildMetabNormExportSessionData",
    "saveMetabNormExportSessionRdsFiles",
    "saveMetabNormExportMetadataFiles",
    "saveMetabNormExportSummaryFile",
    "runMetabNormExportSessionWorkflow",
    "checkMetabNormExportSessionReady",
    "dispatchMetabNormExportSession",
    "handleMetabNormExportSessionOutcome",
    "resolveMetabNormManualItsdFeatureIds",
    "runMetabNormPreNormalizationQcStep",
    "runMetabNormItsdProgressApplyShell",
    "runMetabNormItsdNormalizationStep",
    "runMetabNormLog2ProgressApplyShell",
    "runMetabNormLog2TransformationStep",
    "runMetabNormBetweenSampleProgressApplyShell",
    "runMetabNormBetweenSampleNormalizationStep",
    "runMetabNormPostNormalizationQcStep",
    "runMetabNormRuvProgressApplyShell",
    "runMetabNormRuvOptimizationStep",
    "runMetabNormRuvCorrectionStep",
    "runMetabNormRuvQcStep",
    "buildMetabNormLabelPlot",
    "buildMetabNormTitlePlot",
    "loadMetabNormImageAsPlot",
    "generateMetabNormCompositeFromFiles",
    "generateMetabNormPreNormalizationQc",
    "runMetabNormCompositeQcFigureStep",
    "runMetabNormCompositeQcRefreshShell",
    "runMetabNormNormalizationPipelineShell",
    "runMetabNormNormalizationObserverShell",
    "runMetabNormExportSessionObserverShell",
    "runMetabNormResetNormalizationObserverShell",
    "resolveMetabNormSkipCorrelationInputObject",
    "completeMetabNormSkipCorrelationState",
    "handleMetabNormSkipCorrelationOutcome",
    "dispatchMetabNormSkipCorrelation",
    "runMetabNormSkipCorrelationObserverEntry",
    "runMetabNormSkipCorrelationObserverShell",
    "runMetabNormApplyCorrelationObserverShell",
    "runMetabNormApplyCorrelationWorkflow",
    "handleMetabNormApplyCorrelationOutcome",
    "dispatchMetabNormApplyCorrelation",
    "runMetabNormApplyCorrelationObserverEntry"
  ),
  env = environment()
)
if (!methods::isClass("MetabNormSummaryMock")) {
  methods::setClass("MetabNormSummaryMock", slots = c(design_matrix = "data.frame"))
}
if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      annotation_id_column = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      annotation_id_column = "annotation_id"
    )
  )
}
makeMetabNormTestS4 <- function(metabolite_data = list(
                                  Plasma = data.frame(
                                    metabolite_id = "M1",
                                    annotation_id = "A1",
                                    stringsAsFactors = FALSE
                                  )
                                ),
                                metabolite_id_column = "metabolite_id",
                                annotation_id_column = "annotation_id",
                                database_identifier_type = "Unknown",
                                internal_standard_regex = NA_character_,
                                design_matrix = data.frame(
                                  Sample_ID = c("Sample_1", "Sample_2"),
                                  group = c("A", "B"),
                                  stringsAsFactors = FALSE
                                ),
                                sample_id = "Sample_ID",
                                group_id = "group",
                                technical_replicate_id = NA_character_,
                                args = list()) {
  id_columns <- c(
    metabolite_id_column,
    if (!is.na(annotation_id_column) && nzchar(annotation_id_column)) annotation_id_column else character()
  )
  sample_ids <- as.character(design_matrix[[sample_id]])
  normalize_assay <- function(assay_data) {
    if (is.null(assay_data) || !is.data.frame(assay_data)) {
      assay_data <- data.frame(stringsAsFactors = FALSE)
    }
    if (!metabolite_id_column %in% names(assay_data)) {
      assay_data[[metabolite_id_column]] <- if (nrow(assay_data) > 0) {
        paste0("M", seq_len(nrow(assay_data)))
      } else {
        character()
      }
    }
    if (!is.na(annotation_id_column) &&
        nzchar(annotation_id_column) &&
        !annotation_id_column %in% names(assay_data)) {
      assay_data[[annotation_id_column]] <- assay_data[[metabolite_id_column]]
    }
    for (sample_id_value in sample_ids) {
      if (!sample_id_value %in% names(assay_data)) {
        assay_data[[sample_id_value]] <- if (nrow(assay_data) > 0) {
          seq_len(nrow(assay_data))
        } else {
          numeric()
        }
      }
    }
    assay_data[, unique(c(id_columns, setdiff(names(assay_data), id_columns))), drop = FALSE]
  }
  normalized_metabolite_data <- lapply(metabolite_data, normalize_assay)

  methods::new(
    "MetaboliteAssayData",
    metabolite_data = normalized_metabolite_data,
    metabolite_id_column = metabolite_id_column,
    annotation_id_column = annotation_id_column,
    database_identifier_type = database_identifier_type,
    internal_standard_regex = internal_standard_regex,
    design_matrix = design_matrix,
    sample_id = sample_id,
    group_id = group_id,
    technical_replicate_id = technical_replicate_id,
    args = args
  )
}

# test -------------------------------------------------------------------------
capture <- new.env(parent = emptyenv())
capture$req <- list()
capture$log <- character()
capture$notifications <- list()
capture$save_state <- NULL
workflow_data <- omicsTestOwner(omic_type = "metabolomics", config = list(method = "TIC"))
state_manager <- workflow_data$state_manager
post_filter_object <- structure(list(stage = "post_filter"), class = "MetabNormResetMock")
norm_data <- new.env(parent = emptyenv())
norm_data$post_filter_obj <- post_filter_object
norm_data$normalization_complete <- TRUE
norm_data$ruv_complete <- TRUE
norm_data$correlation_filtering_complete <- TRUE
norm_data$post_norm_obj <- "post-norm"
norm_data$ruv_corrected_obj <- "ruv-corrected"
norm_data$correlation_filtered_obj <- "corr-filtered"
norm_data$ruv_optimization_results <- list(Plasma = data.frame(best_k = 2))
visible <- withVisible(
    runMetabNormResetNormalizationObserverShell(
      workflowData = workflow_data,
      normData = norm_data,
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      showNotificationFn = function(message, type = NULL, duration = NULL) {
        capture$notifications <- c(
          capture$notifications,
          list(list(message = message, type = type, duration = duration))
        )
        invisible(NULL)
      },
      reqFn = function(value) {
        capture$req <- c(capture$req, list(value))
        invisible(value)
      }
    )
  )
expect_false(visible$visible)
expect_identical(capture$req, list(state_manager))
expect_identical(state_manager$current_state, "metab_reset")
expect_identical(state_manager$getState("metab_reset"), post_filter_object)
expect_false(norm_data$normalization_complete)
expect_false(norm_data$ruv_complete)
expect_false(norm_data$correlation_filtering_complete)
