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
    design_field = "design_matrix", feature_noun = "lipid", long_table = "da_lipids_long", design_input_prefix = "lipid_assay_",
    step_scope = "the lipid step", da_scope = "the lipid DA",
    intensity_parameters = c("lipids_intensity_cutoff_percentile", "lipids_proportion_of_samples_below_cutoff"),
    intensity_threshold = "min_lipid_intensity_threshold",
    step_labels = c(lipid_intensity_filter = "2_Intensity_Filtered", lipid_duplicate_resolution = "3_Duplicates_Resolved",
        lipid_qc_finalize = "4_QC_Complete", lipid_itsd_normalisation = "5_ITSD_Normalised", lipid_log2 = "6_Log2",
        lipid_between_sample_normalisation = "7_Normalised", lipid_ruv = "8_RUV_Corrected",
        lipid_pearson_selection = "9_Correlation_Filtered", lipid_norm_complete = "9_Normalisation_Complete",
        lipid_reset = "0_Reset"),
    # KIT-094: the steps the shared core registers, runs, replays and records for this omic.
    steps = omicsRouteSteps,
    # KIT-078: the declarations on the worker's `omics` role.
    worker_steps = c("lipid_intensity_filter", "lipid_duplicate_resolution", "lipid_itsd_normalisation", "lipid_log2",
        "lipid_between_sample_normalisation", "lipid_ruv", "lipid_pearson_pairs", "lipid_pearson_selection", "lipid_da_contrast"),
    # KIT-089: the intensity filter defaults to `partitioned` (user, 2026-11-06, on KIT-084's synthetic pairs: exact, peak
    # .77 and .81 of `whole`). KIT-100: log2 too (user, 2026-11-07, parity with proteomics: exact, peak .925 of `whole` in
    # KIT-084); duplicate resolution and every other step stay `whole`.
    default_variants = c(lipid_intensity_filter = "partitioned", lipid_log2 = "partitioned"),
    publications = omicsRoutePublications,
    # The worker's compute handler: the steps whose figures read the filtering metrics, the step whose product holds
    # native references, and the parameter naming a product the step depends on.
    filtering_steps = c("lipid_intensity_filter", "lipid_duplicate_resolution"),
    native_steps = "lipid_da_contrast",
    dependency_parameters = "pairs",
    replay = omicsRouteReplay,
    record_outcomes = omicsRouteRecordOutcomes,
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
    },
    # KIT-076: the kernel's halves and steps, for one declaration per contrast.
    da_preflight = function(object, contrasts, formula) {
        preflight <- validateLipidDaAnalysisPreflight(currentS4 = object, contrastsTbl = contrasts, formulaString = formula)
        if (!isTRUE(preflight$valid)) stop(paste(preflight$errors, collapse = "; "))
        invisible(NULL)
    },
    da_assay_input = function(...) prepareLipidDaAssayInput(...),
    da_fit = function(...) fitLipidDaAssayModel(...),
    da_apply = function(...) applyLipidDaContrasts(...),
    da_format = function(...) formatLipidDaContrastTable(...),
    da_count = function(...) countLipidDaSignificant(...),
    # KIT-079: the kernels' per-page and per-batch halves, for the partitioned variants.
    filter_intensity_page = function(page, threshold, proportion, id_column) lipidIntensityFilteringHelper(assay_table = page,
        min_lipid_intensity_threshold = threshold, lipids_proportion_of_samples_below_cutoff = proportion, lipid_id_column = id_column),
    resolve_duplicates_batch = function(batch, id_column, sample_columns) resolveLipidDuplicateFeaturesByIntensity(assay_tibble = batch,
        id_col = id_column, sample_cols = sample_columns),
    duplicate_value = function(before, after, without_samples) list(counts = list(before = before, after = after)),
    # KIT-078: what the worker reads from the matrix for the step's figures and files.
    qc_dir_role = "lipid_qc_dir",
    qc_plots = function(...) generateLipidQcPlots(...),
    composite_figure = function(qc_dir, output_dir, assay_names, ruv_mode) writeWorkflowOmicsNormComposite(qc_dir, assay_names, ruv_mode, "lipidomics",
        output_dir = output_dir),
    filtering_metrics = function(object) computeLipidFilteringMetrics(object),
    da_export = function(...) outputLipidDaResultsAllContrasts(...),
    # KIT-082: what the tabs show from the whole matrix.
    final_pca = function(object, color, shape) omicsNormFinalQcPcaPlot(object, color, shape),
    volcano_widget = function(...) generateLipidDAVolcanoPlotGlimma(...)
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

#' The state without its matrix, and with only the id column (KIT-082).
#' @noRd
readLipidWorkflowShell <- function(manager, reference = NULL, measures = character()) readWorkflowOmicsShell(lipidomicsDeclarationBinding, manager, reference, measures)

#' @noRd
readLipidWorkflowIds <- function(manager, reference = NULL) readWorkflowOmicsIds(lipidomicsDeclarationBinding, manager, reference)

#' Identify the lipidomics DA manifest.
#' @noRd
isLipidResultManifest <- function(value) {
    is.list(value) && isTRUE(value[["schema", exact = TRUE]] %in% c("lipid.da_results.v1", "lipid.da_results.v2"))
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
lipidDaContrastOperation <- function() omicsDaContrastOperation(lipidomicsDeclarationBinding)

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

#' The design state's object (KIT-103 C7): the `LipidomicsAssayData` of a builder result or an imported
#' design, with today's argument mapping (the builder commit and the import, each at 70af7bb4): the
#' mapping's id column, its annotation and internal-standard pattern when set (else NA), the
#' `Run`, `group` and `tech_rep_group` roles and the configuration as `args`. The shape of a
#' configuration's `design$build_object(results, columnMapping)`.
#' @param results A list with `data_cln`, `design_matrix` (with `tech_rep_group`) and `config_list`.
#' @param columnMapping The column mapping.
#' @noRd
lipidDesignObject <- function(results, columnMapping) {
    annotation <- columnMapping$annotation_col
    if (is.null(annotation) || is.na(annotation) || !nzchar(annotation)) annotation <- NA_character_
    pattern <- columnMapping$is_pattern
    if (is.null(pattern) || is.na(pattern)) pattern <- NA_character_
    createLipidomicsAssayData(lipid_data = results$data_cln, design_matrix = results$design_matrix,
        lipid_id_column = columnMapping$lipid_id_col,
        annotation_id_column = annotation, sample_id = "Run", group_id = "group",
        technical_replicate_id = "tech_rep_group", database_identifier_type = "Unknown",
        internal_standard_regex = pattern, args = results$config_list)
}

#' @noRd
runLipidReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    runWorkflowOmicsReferenceDaAnalysis(lipidomicsDeclarationBinding, manager, source, settings, contrasts)
}

#' @noRd
replayLipidomicsStep <- function(manager, step, parameters) {
    replayWorkflowOmicsStep(lipidomicsDeclarationBinding, manager, step, parameters)
}
