# KIT-031: the lipidomics declarations, the kit's completeness test (the step map `kit-001/step-map-lipid.json`).
# KIT-073: built by the shared omics core (utils_workflow_omics_declarations.R) from the binding below: the class, the
# stored slots, the names and strings the lipid route records, and the adaptors that call the lipid kernels under the
# names the lipid modules call them by (so the recorded `@args` keep their keys). The names below are the ones the
# modules, the record, the harness and the tests use.

lipidomicsOperationFields <- c("lipid_data", "design_matrix")

#' The merge rules of a partitioned lipid declaration: the assay tables concatenate, the design is recomputed once.
#' @noRd
lipidomicsPartitionMerge <- c(lipid_data = "concatenate", design_matrix = "design")

lipidRuvParameterNames <- omicsRuvParameterNames

#' The lipidomics binding of the omics declaration core.
#' @noRd
lipidomicsDeclarationBinding <- list(
    omic_type = "lipidomics", prefix = "lipid", class = "LipidomicsAssayData", data_slot = "lipid_data",
    id_slot = "lipid_id_column", fields = lipidomicsOperationFields, merge = lipidomicsPartitionMerge,
    feature_noun = "lipid", long_table = "da_lipids_long", design_input_prefix = "lipid_assay_",
    step_scope = "the lipid step", da_scope = "the lipid DA",
    intensity_parameters = c("lipids_intensity_cutoff_percentile", "lipids_proportion_of_samples_below_cutoff"),
    intensity_threshold = "min_lipid_intensity_threshold",
    step_labels = c(lipid_intensity_filter = "2_Intensity_Filtered", lipid_duplicate_resolution = "3_Duplicates_Resolved",
        lipid_qc_finalize = "4_QC_Complete", lipid_itsd_normalisation = "5_ITSD_Normalised", lipid_log2 = "6_Log2",
        lipid_between_sample_normalisation = "7_Normalised", lipid_ruv = "8_RUV_Corrected",
        lipid_pearson_selection = "9_Correlation_Filtered", lipid_norm_complete = "9_Normalisation_Complete",
        lipid_reset = "0_Reset"),
    filter_intensity = function(object, parameters) {
        lipidIntensityFilteringFn <- lipidIntensityFiltering
        lipidIntensityFilteringFn(object,
            lipids_intensity_cutoff_percentile = parameters$lipids_intensity_cutoff_percentile,
            lipids_proportion_of_samples_below_cutoff = parameters$lipids_proportion_of_samples_below_cutoff)
    },
    # Per assay, the row with the highest mean intensity is kept; the module reads the row counts before and after.
    resolve_duplicates = function(object) {
        resolveLipidDuplicateFeaturesByIntensityFn <- resolveLipidDuplicateFeaturesByIntensity
        assays <- object@lipid_data
        names(assays) <- names(assays) %||% paste0("Assay_", seq_along(assays))
        before <- vapply(assays, nrow, integer(1))
        resolved <- lapply(assays, function(assay) {
            sample_cols <- names(assay)[vapply(assay, is.numeric, logical(1))]
            if (!length(sample_cols)) return(assay)
            resolveLipidDuplicateFeaturesByIntensityFn(assay_tibble = assay, id_col = object@lipid_id_column,
                sample_cols = sample_cols)
        })
        object@lipid_data <- resolved
        list(object = object, value = list(counts = list(before = before, after = vapply(resolved, nrow, integer(1)))))
    },
    optimise_ruv = function(object, mode, params) {
        runLipidPerAssayRuvOptimizationFn <- runLipidPerAssayRuvOptimization
        runLipidPerAssayRuvOptimizationFn(theObject = object, ruv_mode = mode, params = params, experiment_paths = NULL)
    },
    extract_best_k = function(results) extractLipidBestKPerAssay(results),
    extract_controls = function(results) extractLipidCtrlPerAssay(results),
    filter_correlation = function(object, correlations, threshold) {
        filterSamplesByLipidCorrelationThresholdFn <- filterSamplesByLipidCorrelationThreshold
        filterSamplesByLipidCorrelationThresholdFn(theObject = object, pearson_correlation_per_pair = correlations,
            min_pearson_correlation_threshold = threshold)
    },
    run_da = function(object, parameters) {
        runLipidsDAFn <- runLipidsDA
        runLipidsDAFn(theObject = object, contrasts_tbl = parameters$contrasts, formula_string = parameters$formula,
            da_q_val_thresh = parameters$qvalue, treat_lfc_cutoff = parameters$lfc, eBayes_trend = TRUE, eBayes_robust = TRUE)
    }
)

lipidomicsStepLabels <- lipidomicsDeclarationBinding$step_labels

#' Unique feature ids per assay (the QC tabs' count), by the id column when the table has it.
#' @noRd
lipidomicsFeatureCounts <- function(object) omicsFeatureCounts(lipidomicsDeclarationBinding, object)

#' The lipidomics state read for one action, and its metadata without a table (KIT-075).
#' @noRd
readLipidWorkflowState <- function(manager, reference = NULL) readWorkflowOmicsState(lipidomicsDeclarationBinding, manager, reference)

#' @noRd
readLipidWorkflowMetadata <- function(manager, reference = NULL) readWorkflowOmicsMetadata(lipidomicsDeclarationBinding, manager, reference)

#' Identify the lipidomics DA manifest.
#' @noRd
isLipidResultManifest <- function(value) {
    is.list(value) && identical(value[["schema", exact = TRUE]], "lipid.da_results.v1")
}

#' Every lipidomics declaration, by name.
#' @noRd
lipidomicsWorkflowOperations <- function() omicsWorkflowOperations(lipidomicsDeclarationBinding)

lipidIntensityFilterOperation <- function() omicsIntensityFilterOperation(lipidomicsDeclarationBinding)
lipidDuplicateResolutionOperation <- function() omicsDuplicateResolutionOperation(lipidomicsDeclarationBinding)
lipidQcFinalizeOperation <- function() omicsStoreOperation(lipidomicsDeclarationBinding, "lipid_qc_finalize")
lipidNormCompleteOperation <- function() omicsStoreOperation(lipidomicsDeclarationBinding, "lipid_norm_complete")
lipidResetOperation <- function() omicsStoreOperation(lipidomicsDeclarationBinding, "lipid_reset")
lipidItsdNormalisationOperation <- function() omicsItsdNormalisationOperation(lipidomicsDeclarationBinding)
lipidLog2Operation <- function() omicsLog2Operation(lipidomicsDeclarationBinding)
lipidBetweenSampleNormalisationOperation <- function() omicsBetweenSampleNormalisationOperation(lipidomicsDeclarationBinding)
lipidRuvOperation <- function() omicsRuvOperation(lipidomicsDeclarationBinding)
lipidPearsonPairsOperation <- function() omicsPearsonPairsOperation(lipidomicsDeclarationBinding)
lipidPearsonSelectionOperation <- function() omicsPearsonSelectionOperation(lipidomicsDeclarationBinding)
lipidDaOperation <- function() omicsDaOperation(lipidomicsDeclarationBinding)

# The lipid names of the wrappers (KIT-031), one body each in the core.

#' @noRd
runLipidReferenceOperation <- runWorkflowOmicsReferenceOperation

#' @noRd
markLipidReferenceStoreStep <- markWorkflowOmicsStoreStep

#' @noRd
commitLipidReferenceStoreStep <- commitWorkflowOmicsStoreStep

#' @noRd
runLipidReferenceIntensityFilter <- function(manager, percentile, proportion, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceIntensityFilter(lipidomicsDeclarationBinding, manager, percentile, proportion, config, description,
        state_name, hydrate)
}

#' @noRd
runLipidReferenceDuplicateResolution <- function(manager, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceDuplicateResolution(lipidomicsDeclarationBinding, manager, config, description, state_name, hydrate)
}

#' @noRd
runLipidReferenceItsdNormalisation <- function(manager, aggregation, feature_ids, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceItsdNormalisation(lipidomicsDeclarationBinding, manager, aggregation, feature_ids, config, description,
        state_name, hydrate)
}

#' @noRd
runLipidReferenceLog2 <- function(manager, offset, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceLog2(lipidomicsDeclarationBinding, manager, offset, config, description, state_name, hydrate)
}

#' @noRd
runLipidReferenceBetweenSampleNormalisation <- function(manager, method, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceBetweenSampleNormalisation(lipidomicsDeclarationBinding, manager, method, config, description,
        state_name, hydrate)
}

#' @noRd
runLipidReferenceRuv <- function(manager, parameters, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceRuv(lipidomicsDeclarationBinding, manager, parameters, config, description, state_name, hydrate)
}

#' @noRd
runLipidReferenceCorrelation <- function(manager, group, threshold, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceCorrelation(lipidomicsDeclarationBinding, manager, group, threshold, config, description,
        state_name, hydrate)
}

#' @noRd
recordLipidDesignProvenance <- function(manager, object, source_dir = NULL, design_file = NULL) {
    recordWorkflowOmicsDesignProvenance(lipidomicsDeclarationBinding, manager, object, source_dir, design_file)
}

#' @noRd
runLipidReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    runWorkflowOmicsReferenceDaAnalysis(lipidomicsDeclarationBinding, manager, source, settings, contrasts)
}

#' @noRd
replayLipidomicsStep <- function(manager, step, parameters) {
    replayWorkflowOmicsStep(lipidomicsDeclarationBinding, manager, step, parameters)
}
