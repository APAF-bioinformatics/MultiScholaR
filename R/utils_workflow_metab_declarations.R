# KIT-032: the metabolomics declarations, the KIT-031 twins (the step map `kit-001/step-map-metab.json`).
# KIT-073: built by the shared omics core (utils_workflow_omics_declarations.R) from the binding below: the class, the
# stored slots, the names and strings the metabolomics route records, and the adaptors that call the metabolomics
# kernels under the names the metab modules call them by (so the recorded `@args` keep their keys). The names below
# are the ones the modules, the record, the harness and the tests use.

metabolomicsOperationFields <- c("metabolite_data", "design_matrix")

#' The merge rules of a partitioned metabolomics declaration.
#' @noRd
metabolomicsPartitionMerge <- c(metabolite_data = "concatenate", design_matrix = "design")

metabRuvParameterNames <- omicsRuvParameterNames

#' The metabolomics binding of the omics declaration core.
#' @noRd
metabolomicsDeclarationBinding <- list(
    omic_type = "metabolomics", prefix = "metab", class = "MetaboliteAssayData", data_slot = "metabolite_data",
    id_slot = "metabolite_id_column", fields = metabolomicsOperationFields, merge = metabolomicsPartitionMerge,
    design_field = "design_matrix", feature_noun = "metabolite", long_table = "da_metabolites_long", design_input_prefix = "metab_assay_",
    step_scope = "the step", da_scope = "the metabolomics DA",
    intensity_parameters = c("metabolites_intensity_cutoff_percentile", "metabolites_proportion_of_samples_below_cutoff"),
    intensity_threshold = "min_metabolite_intensity_threshold",
    step_labels = c(metab_intensity_filter = "2_Intensity_Filtered", metab_duplicate_resolution = "3_Duplicates_Resolved",
        metab_qc_finalize = "4_QC_Complete", metab_itsd_normalisation = "5_ITSD_Normalised", metab_log2 = "6_Log2",
        metab_between_sample_normalisation = "7_Normalised", metab_ruv = "8_RUV_Corrected",
        metab_pearson_selection = "9_Correlation_Filtered", metab_norm_complete = "9_Normalisation_Complete",
        metab_reset = "0_Reset"),
    # The module calls the kernel by its own name, so the adaptor does too (the `@args` key `metaboliteIntensityFiltering`).
    filter_intensity = function(object, parameters) {
        metaboliteIntensityFiltering(object,
            metabolites_intensity_cutoff_percentile = parameters$metabolites_intensity_cutoff_percentile,
            metabolites_proportion_of_samples_below_cutoff = parameters$metabolites_proportion_of_samples_below_cutoff)
    },
    # As `resolveMetabDuplicateAssayData` runs it per assay; the module reads its per-assay statistics.
    resolve_duplicates = function(object) {
        resolved <- resolveMetabDuplicateAssayData(assayList = object@metabolite_data,
            metaboliteIdCol = object@metabolite_id_column, logWarnFn = function(...) NULL)
        object@metabolite_data <- resolved$resolvedAssayList
        list(object = object, value = list(stats = resolved$statsList))
    },
    optimise_ruv = function(object, mode, params) {
        runPerAssayRuvOptimizationFn <- runPerAssayRuvOptimization
        runPerAssayRuvOptimizationFn(theObject = object, ruv_mode = mode, params = params, experiment_paths = NULL)
    },
    extract_best_k = function(results) extractBestKPerAssay(results),
    extract_controls = function(results) extractCtrlPerAssay(results),
    filter_correlation = function(object, correlations, threshold) {
        filterSamplesByMetaboliteCorrelationThresholdFn <- filterSamplesByMetaboliteCorrelationThreshold
        filterSamplesByMetaboliteCorrelationThresholdFn(theObject = object, pearson_correlation_per_pair = correlations,
            min_pearson_correlation_threshold = threshold)
    },
    run_da = function(object, parameters) {
        runMetabolitesDAFn <- runMetabolitesDA
        runMetabolitesDAFn(theObject = object, contrasts_tbl = parameters$contrasts, formula_string = parameters$formula,
            da_q_val_thresh = parameters$qvalue, treat_lfc_cutoff = parameters$lfc, eBayes_trend = TRUE, eBayes_robust = TRUE)
    },
    # KIT-076: the kernel's halves and steps, for one declaration per contrast.
    da_preflight = function(object, contrasts, formula) invisible(NULL),
    da_assay_input = function(...) prepareMetabDaAssayInput(...),
    da_fit = function(...) fitMetabDaAssayModel(...),
    da_apply = function(...) applyMetabDaContrasts(...),
    da_format = function(...) formatMetabDaContrastTable(...),
    da_count = function(...) countMetabDaSignificant(...),
    # KIT-079: the kernels' per-page and per-batch halves, for the partitioned variants.
    filter_intensity_page = function(page, threshold, proportion, id_column) metaboliteIntensityFilteringHelper(assay_table = page,
        min_metabolite_intensity_threshold = threshold, metabolites_proportion_of_samples_below_cutoff = proportion, metabolite_id_column = id_column),
    resolve_duplicates_batch = function(batch, id_column, sample_columns) resolveDuplicateFeaturesByIntensity(assay_tibble = batch,
        id_col = id_column, sample_cols = sample_columns),
    # As `resolveMetabDuplicateAssayData` reports each assay.
    duplicate_value = function(before, after, without_samples) list(stats = lapply(stats::setNames(names(before), names(before)), function(assay) {
        list(original = before[[assay]], resolved = after[[assay]], removed = if (without_samples[[assay]]) 0 else before[[assay]] - after[[assay]])
    })),
    # KIT-078: what the worker reads from the matrix for the step's figures and files.
    qc_dir_role = "metabolite_qc_dir",
    qc_plots = function(...) generateMetabQcPlots(...),
    composite_figure = function(qc_dir, output_dir, assay_names, ruv_mode) runMetabNormCompositeQcFigureStep(list(metabolite_qc_dir = qc_dir), assay_names, ruv_mode,
        "metabolomics", generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles, outputDir = output_dir),
    filtering_metrics = function(object) computeMetabFilteringMetrics(object),
    da_export = function(...) outputMetabDaResultsAllContrasts(...)
)

metabolomicsStepLabels <- metabolomicsDeclarationBinding$step_labels

#' Unique feature ids per assay, by the id column when the table has it.
#' @noRd
metabolomicsFeatureCounts <- function(object) omicsFeatureCounts(metabolomicsDeclarationBinding, object)

#' The metabolomics state read for one action, and its metadata without a table (KIT-075).
#' @noRd
readMetabWorkflowState <- function(manager, reference = NULL) readWorkflowOmicsState(metabolomicsDeclarationBinding, manager, reference)

#' @noRd
readMetabWorkflowMetadata <- function(manager, reference = NULL) readWorkflowOmicsMetadata(metabolomicsDeclarationBinding, manager, reference)

#' Identify the metabolomics DA manifest.
#' @noRd
isMetabResultManifest <- function(value) {
    is.list(value) && isTRUE(value[["schema", exact = TRUE]] %in% c("metab.da_results.v1", "metab.da_results.v2"))
}

#' Every metabolomics declaration, by name.
#' @noRd
metabolomicsWorkflowOperations <- function() omicsWorkflowOperations(metabolomicsDeclarationBinding)

metabIntensityFilterOperation <- function() omicsIntensityFilterOperation(metabolomicsDeclarationBinding)
metabDuplicateResolutionOperation <- function() omicsDuplicateResolutionOperation(metabolomicsDeclarationBinding)
metabQcFinalizeOperation <- function() omicsStoreOperation(metabolomicsDeclarationBinding, "metab_qc_finalize")
metabNormCompleteOperation <- function() omicsStoreOperation(metabolomicsDeclarationBinding, "metab_norm_complete")
metabResetOperation <- function() omicsStoreOperation(metabolomicsDeclarationBinding, "metab_reset")
metabItsdNormalisationOperation <- function() omicsItsdNormalisationOperation(metabolomicsDeclarationBinding)
metabLog2Operation <- function() omicsLog2Operation(metabolomicsDeclarationBinding)
metabBetweenSampleNormalisationOperation <- function() omicsBetweenSampleNormalisationOperation(metabolomicsDeclarationBinding)
metabRuvOperation <- function() omicsRuvOperation(metabolomicsDeclarationBinding)
metabPearsonPairsOperation <- function() omicsPearsonPairsOperation(metabolomicsDeclarationBinding)
metabPearsonSelectionOperation <- function() omicsPearsonSelectionOperation(metabolomicsDeclarationBinding)
metabDaOperation <- function() omicsDaOperation(metabolomicsDeclarationBinding)
metabDaContrastOperation <- function() omicsDaContrastOperation(metabolomicsDeclarationBinding)

# The metabolomics names of the wrappers (KIT-032), one body each in the core.

#' @noRd
runMetabReferenceIntensityFilter <- function(manager, percentile, proportion, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceIntensityFilter(metabolomicsDeclarationBinding, manager, percentile, proportion, config, description,
        state_name, hydrate)
}

#' @noRd
runMetabReferenceDuplicateResolution <- function(manager, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceDuplicateResolution(metabolomicsDeclarationBinding, manager, config, description, state_name, hydrate)
}

#' @noRd
runMetabReferenceItsdNormalisation <- function(manager, aggregation, feature_ids, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceItsdNormalisation(metabolomicsDeclarationBinding, manager, aggregation, feature_ids, config,
        description, state_name, hydrate)
}

#' @noRd
runMetabReferenceLog2 <- function(manager, offset, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceLog2(metabolomicsDeclarationBinding, manager, offset, config, description, state_name, hydrate)
}

#' @noRd
runMetabReferenceBetweenSampleNormalisation <- function(manager, method, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceBetweenSampleNormalisation(metabolomicsDeclarationBinding, manager, method, config, description,
        state_name, hydrate)
}

#' @noRd
runMetabReferenceRuv <- function(manager, parameters, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceRuv(metabolomicsDeclarationBinding, manager, parameters, config, description, state_name, hydrate)
}

#' @noRd
runMetabReferenceCorrelation <- function(manager, group, threshold, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceCorrelation(metabolomicsDeclarationBinding, manager, group, threshold, config, description,
        state_name, hydrate)
}

#' @noRd
recordMetabDesignProvenance <- function(manager, object, source_dir = NULL, design_file = NULL) {
    recordWorkflowOmicsDesignProvenance(metabolomicsDeclarationBinding, manager, object, source_dir, design_file)
}

#' @noRd
runMetabReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    runWorkflowOmicsReferenceDaAnalysis(metabolomicsDeclarationBinding, manager, source, settings, contrasts)
}

#' @noRd
replayMetabolomicsStep <- function(manager, step, parameters) {
    replayWorkflowOmicsStep(metabolomicsDeclarationBinding, manager, step, parameters)
}
