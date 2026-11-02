# Extracted from test-metab-03b-norm-module-helper-characterization.R:4792

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
omicsStandInReaders()
capture <- new.env(parent = emptyenv())
capture$order <- character()
capture$req <- list()
initial_s4 <- structure(list(stage = "initial"), class = "MetabNormPipelineShellMock")
pre_norm_s4 <- structure(list(stage = "pre"), class = "MetabNormPipelineShellMock")
itsd_s4 <- structure(list(stage = "itsd"), class = "MetabNormPipelineShellMock")
log2_s4 <- structure(list(stage = "log2"), class = "MetabNormPipelineShellMock")
between_s4 <- structure(list(stage = "between"), class = "MetabNormPipelineShellMock")
post_norm_s4 <- structure(list(stage = "post"), class = "MetabNormPipelineShellMock")
ruv_s4 <- structure(list(stage = "ruv"), class = "MetabNormPipelineShellMock")
final_s4 <- structure(list(stage = "final"), class = "MetabNormPipelineShellMock")
workflow_data <- list(
    state_manager = list(
      getState = function() {
        capture$order <- c(capture$order, "getState")
        initial_s4
      }
    )
  )
input_values <- list(
    apply_itsd = TRUE,
    itsd_aggregation = "median",
    log_offset = 1,
    norm_method = "median",
    ruv_mode = "per_assay",
    auto_percentage_min = 0.1,
    auto_percentage_max = 0.3,
    ruv_grouping_variable = "group",
    separation_metric = "pearson",
    k_penalty_weight = 0.25,
    adaptive_k_penalty = FALSE,
    ruv_k = 2,
    ruv_percentage = 0.2
  )
experiment_paths <- list(metabolite_qc_dir = "/tmp/metab-qc")
norm_data <- new.env(parent = emptyenv())
norm_data$assay_names <- c("Plasma", "Urine")
norm_data$itsd_selections <- list(Plasma = 1:2)
composite_generator <- function(...) {
    "unused-generator"
  }
save_plot <- function(...) {
    "unused-save"
  }
warn_log <- function(message) {
    capture$warn <- c(capture$warn, message)
  }
composite_qc_state <- list(qcDir = "/tmp/metab-qc", compositeSaved = TRUE)
visible <- withVisible(
    runMetabNormNormalizationPipelineShell(
      workflowData = workflow_data,
      inputValues = input_values,
      experimentPaths = experiment_paths,
      omicType = "metabolite",
      normData = norm_data,
      getPlotAestheticsFn = function() {
        capture$order <- c(capture$order, "aesthetics")
        list(color_var = "group", shape_var = "batch")
      },
      addLogFn = function(message) {
        capture$log <- c(capture$log, message)
      },
      reqFn = function(value) {
        capture$order <- c(capture$order, "req")
        capture$req <- c(capture$req, list(value))
        invisible(value)
      },
      generateCompositeFromFilesFn = composite_generator,
      savePlotFn = save_plot,
      logWarnFn = warn_log,
      runPreNormalizationQcStepFn = function(currentS4, totalSteps, experimentPaths, groupingVariable, shapeVariable, normData, addLogFn) {
        capture$order <- c(capture$order, "pre")
        capture$pre <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = pre_norm_s4)
      },
      runItsdProgressApplyShellFn = function(currentS4, totalSteps, applyItsd, itsdAggregation, itsdSelections, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "itsd")
        capture$itsd <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          applyItsd = applyItsd,
          itsdAggregation = itsdAggregation,
          itsdSelections = itsdSelections,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = itsd_s4)
      },
      runLog2ProgressApplyShellFn = function(currentS4, totalSteps, logOffset, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "log2")
        capture$log2 <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          logOffset = logOffset,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = log2_s4)
      },
      runBetweenSampleProgressApplyShellFn = function(currentS4, totalSteps, normMethod, workflowData, normData, addLogFn) {
        capture$order <- c(capture$order, "between")
        capture$between <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          normMethod = normMethod,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = between_s4)
      },
      runPostNormalizationQcStepFn = function(currentS4, experimentPaths, groupingVariable, shapeVariable, addLogFn) {
        capture$order <- c(capture$order, "post")
        capture$post <- list(
          currentS4 = currentS4,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          addLogFn = addLogFn
        )
        list(currentS4 = post_norm_s4)
      },
      runRuvProgressApplyShellFn = function(
        currentS4,
        totalSteps,
        ruvMode,
        autoPercentageMin,
        autoPercentageMax,
        ruvGroupingVariable,
        separationMetric,
        kPenaltyWeight,
        adaptiveKPenalty,
        manualK,
        manualPercentage,
        experimentPaths,
        groupingVariable,
        shapeVariable,
        workflowData,
        normData,
        addLogFn
      ) {
        capture$order <- c(capture$order, "ruv")
        capture$ruv <- list(
          currentS4 = currentS4,
          totalSteps = totalSteps,
          ruvMode = ruvMode,
          autoPercentageMin = autoPercentageMin,
          autoPercentageMax = autoPercentageMax,
          ruvGroupingVariable = ruvGroupingVariable,
          separationMetric = separationMetric,
          kPenaltyWeight = kPenaltyWeight,
          adaptiveKPenalty = adaptiveKPenalty,
          manualK = manualK,
          manualPercentage = manualPercentage,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          shapeVariable = shapeVariable,
          workflowData = workflowData,
          normData = normData,
          addLogFn = addLogFn
        )
        list(currentS4 = ruv_s4)
      },
      runCompositeQcRefreshShellFn = function(
        currentS4,
        experimentPaths,
        assayNames,
        ruvMode,
        omicType,
        normData,
        addLogFn,
        generateCompositeFromFilesFn,
        savePlotFn,
        logWarnFn
      ) {
        capture$order <- c(capture$order, "composite")
        capture$composite <- list(
          currentS4 = currentS4,
          experimentPaths = experimentPaths,
          assayNames = assayNames,
          ruvMode = ruvMode,
          omicType = omicType,
          normData = normData,
          addLogFn = addLogFn,
          generateCompositeFromFilesFn = generateCompositeFromFilesFn,
          savePlotFn = savePlotFn,
          logWarnFn = logWarnFn
        )
        list(
          currentS4 = final_s4,
          compositeQcState = composite_qc_state,
          plotRefreshTrigger = 12
        )
      }
    )
  )
