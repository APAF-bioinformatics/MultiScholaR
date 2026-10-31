# KIT-032: the metabolomics declarations, the KIT-031 twins on `MetaboliteAssayData` (`kit-001/step-map-metab.json`).
# Every step in process (`compute = "in_process"`, `runSavedOperation` on the `whole` variant); the kernels unchanged and
# called under the names the metab module calls them by, so the recorded `@args` keep their keys. The wrappers commit
# under the module's state names with its configuration and descriptions through the omic-neutral core of
# utils_workflow_lipid_declarations.R (`runWorkflowOmicsReferenceOperation`, `commitWorkflowOmicsStoreStep`,
# `markWorkflowOmicsStoreStep`). Nothing here reads GUI inputs, the global environment or a Shiny session.

metabolomicsOperationFields <- c("metabolite_data", "design_matrix")

#' The merge rules of a partitioned metabolomics declaration.
#' @noRd
metabolomicsPartitionMerge <- c(metabolite_data = "concatenate", design_matrix = "design")

#' One metabolomics declaration: `MetaboliteAssayData` in, both stored components hydrated.
#' @noRd
newMetabolomicsOperation <- function(name, unit, parameters, product_type, kernel, output_class = "MetaboliteAssayData",
                                     global_statistics = NULL, fixed = NULL) {
    partitioned <- unit %in% c("row", "complete_sample", "feature_group")
    newWorkflowOperation(name, "MetaboliteAssayData", fields = metabolomicsOperationFields, unit = unit, kernel = kernel,
        output_class = output_class, parameters = parameters, product_type = product_type,
        global_statistics = global_statistics, merge = if (partitioned) metabolomicsPartitionMerge else NULL, fixed = fixed)
}

#' Unique feature ids per assay, by the id column when the table has it.
#' @noRd
metabolomicsFeatureCounts <- function(object) {
    vapply(object@metabolite_data, function(assay) {
        if (object@metabolite_id_column %in% names(assay)) length(unique(assay[[object@metabolite_id_column]])) else nrow(assay)
    }, numeric(1))
}

#' The intensity filter (row, after the per-assay quantile pre-pass; the `ceiling()` removed, user decision 2026-10-28).
#' The module calls the kernel by its own name, so the declaration does too (the `@args` key `metaboliteIntensityFiltering`).
#' @noRd
metabIntensityFilterOperation <- function() {
    newMetabolomicsOperation("metab_intensity_filter", "row",
        parameters = c("metabolites_intensity_cutoff_percentile", "metabolites_proportion_of_samples_below_cutoff"),
        product_type = "metab.intensity_filter_result",
        global_statistics = c(min_metabolite_intensity_threshold =
            "the percentile of every intensity of one assay table, computed once per assay before any row"),
        kernel = function(object, parameters, store) {
            before <- metabolomicsFeatureCounts(object)
            filtered <- metaboliteIntensityFiltering(object,
                metabolites_intensity_cutoff_percentile = parameters$metabolites_intensity_cutoff_percentile,
                metabolites_proportion_of_samples_below_cutoff = parameters$metabolites_proportion_of_samples_below_cutoff)
            object <- NULL
            list(final = store$checkpoint(filtered, "metab-intensity-filtered"),
                counts = list(before = before, after = metabolomicsFeatureCounts(filtered)))
        })
}

#' Duplicate resolution (feature group by the id column), as `resolveMetabDuplicateAssayData` runs it per assay.
#' @noRd
metabDuplicateResolutionOperation <- function() {
    newMetabolomicsOperation("metab_duplicate_resolution", "feature_group", parameters = character(),
        product_type = "metab.duplicate_resolution_result",
        kernel = function(object, parameters, store) {
            resolved <- resolveMetabDuplicateAssayData(assayList = object@metabolite_data,
                metaboliteIdCol = object@metabolite_id_column, logWarnFn = function(...) NULL)
            object@metabolite_data <- resolved$resolvedAssayList
            list(final = store$checkpoint(object, "metab-duplicates-resolved"), stats = resolved$statsList)
        })
}

#' A store step: the current object under a new name. No kernel, no product.
#' @noRd
metabStoreOperation <- function(name) {
    newMetabolomicsOperation(name, "store", parameters = character(), product_type = "metab.store_step",
        kernel = function(object, parameters, store) list())
}
metabQcFinalizeOperation <- function() metabStoreOperation("metab_qc_finalize")
metabNormCompleteOperation <- function() metabStoreOperation("metab_norm_complete")
metabResetOperation <- function() metabStoreOperation("metab_reset")

#' Internal-standard normalisation (complete sample after the mean-factor pre-pass).
#' @noRd
metabItsdNormalisationOperation <- function() {
    newMetabolomicsOperation("metab_itsd_normalisation", "complete_sample",
        parameters = c("itsd_aggregation", "itsd_feature_ids"), product_type = "metab.itsd_normalisation_result",
        global_statistics = c(mean_standard_factor =
            "the mean over all samples of each assay's per-sample aggregate standard signal, computed once before any sample is scaled"),
        fixed = list(method = "ITSD"),
        kernel = function(object, parameters, store) {
            normaliseUntransformedDataFn <- normaliseUntransformedData
            normalised <- normaliseUntransformedDataFn(theObject = object, method = "ITSD",
                itsd_aggregation = parameters$itsd_aggregation, itsd_feature_ids = parameters$itsd_feature_ids)
            object <- NULL
            list(final = store$checkpoint(normalised, "metab-itsd-normalised"))
        })
}

#' Log2 transformation (row).
#' @noRd
metabLog2Operation <- function() {
    newMetabolomicsOperation("metab_log2", "row", parameters = "offset", product_type = "metab.log2_result",
        kernel = function(object, parameters, store) {
            logTransformAssaysFn <- logTransformAssays
            transformed <- logTransformAssaysFn(theObject = object, offset = parameters$offset)
            object <- NULL
            list(final = store$checkpoint(transformed, "metab-log2"))
        })
}

#' Between-sample normalisation (whole matrix; `none` the module's identity, the state still committed).
#' @noRd
metabBetweenSampleNormalisationOperation <- function() {
    newMetabolomicsOperation("metab_between_sample_normalisation", "whole_matrix", parameters = "normalisation_method",
        product_type = "metab.between_sample_normalisation_result",
        kernel = function(object, parameters, store) {
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            normalised <- if (identical(parameters$normalisation_method, "none")) object else
                normaliseBetweenSamplesFn(theObject = object, normalisation_method = parameters$normalisation_method)
            object <- NULL
            list(final = store$checkpoint(normalised, "metab-normalised"),
                applied = !identical(parameters$normalisation_method, "none"))
        })
}

metabRuvParameterNames <- c("ruv_mode", "ruv_grouping_variable", "percentage_min", "percentage_max", "separation_metric",
    "k_penalty_weight", "adaptive_k_penalty", "manual_k", "manual_percentage")

#' RUV-III (whole matrix): the per-assay optimisation, then `ruvIII_C_Varying` with the best k and controls per assay
#' (the module's optimisation and correction steps in one declaration); the optimisation results a product.
#' @noRd
metabRuvOperation <- function() {
    newMetabolomicsOperation("metab_ruv", "whole_matrix", parameters = metabRuvParameterNames,
        product_type = "metab.ruv_result",
        kernel = function(object, parameters, store) {
            runPerAssayRuvOptimizationFn <- runPerAssayRuvOptimization
            ruvIII_C_VaryingFn <- ruvIII_C_Varying
            ruv_results <- runPerAssayRuvOptimizationFn(theObject = object, ruv_mode = parameters$ruv_mode,
                params = parameters[setdiff(metabRuvParameterNames, "ruv_mode")], experiment_paths = NULL)
            best_k <- extractBestKPerAssay(ruv_results)
            controls <- extractCtrlPerAssay(ruv_results)
            corrected <- ruvIII_C_VaryingFn(theObject = object, ruv_grouping_variable = parameters$ruv_grouping_variable,
                ruv_number_k = best_k, ctrl = controls)
            object <- NULL
            final <- store$checkpoint(corrected, "metab-ruv-corrected")
            corrected <- NULL
            saved_results <- lapply(ruv_results, function(result) {
                for (field in intersect(c("cancor_plot", "best_cancor_plot"), names(result))) {
                    if (inherits(result[[field]], c("gg", "ggplot"))) result[[field]] <- stripProtNormQcPlotOwners(result[[field]])
                }
                result
            })
            list(final = final, optimisation = store$product(saved_results, "optimisation", "metab-ruv-optimisation"),
                best_k = best_k,
                control_count = lapply(controls, function(x) if (is.null(x)) NA_integer_ else if (is.logical(x)) sum(x, na.rm = TRUE) else length(x)),
                percentage = lapply(ruv_results, function(x) x$best_percentage %||% NA_real_))
        })
}

#' Pearson correlations per sample pair within the grouping variable (whole matrix; a product).
#' @noRd
metabPearsonPairsOperation <- function() {
    newMetabolomicsOperation("metab_pearson_pairs", "whole_matrix", parameters = "correlation_group",
        product_type = "metab-pearson-pairs", output_class = NULL,
        kernel = function(object, parameters, store) {
            pearsonCorForSamplePairsFn <- pearsonCorForSamplePairs
            pearsonCorForSamplePairsFn(theObject = object, correlation_group = parameters$correlation_group)
        })
}

#' The sample selection by the correlation threshold on the pairs product (whole matrix).
#' @noRd
metabPearsonSelectionOperation <- function() {
    newMetabolomicsOperation("metab_pearson_selection", "whole_matrix",
        parameters = c("pairs", "min_pearson_correlation_threshold"), product_type = "metab.pearson_selection_result",
        kernel = function(object, parameters, store) {
            filterSamplesByMetaboliteCorrelationThresholdFn <- filterSamplesByMetaboliteCorrelationThreshold
            filtered <- filterSamplesByMetaboliteCorrelationThresholdFn(theObject = object,
                pearson_correlation_per_pair = parameters$correlations,
                min_pearson_correlation_threshold = parameters$min_pearson_correlation_threshold)
            object <- NULL
            list(final = store$checkpoint(filtered, "metab-correlation-filtered"))
        })
}

#' Differential abundance (whole matrix, product only): `runMetabolitesDA` once over every contrast and assay; the
#' manifest `metab.da_results.v1` is the step's output (the GUI never builds the S4 result class).
#' @noRd
metabDaOperation <- function() {
    newMetabolomicsOperation("metab_da_output", "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrasts"),
        product_type = "metab-da-references", output_class = NULL,
        fixed = list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            runMetabolitesDAFn <- runMetabolitesDA
            result <- runMetabolitesDAFn(theObject = object, contrasts_tbl = parameters$contrasts,
                formula_string = parameters$formula, da_q_val_thresh = parameters$qvalue,
                treat_lfc_cutoff = parameters$lfc, eBayes_trend = TRUE, eBayes_robust = TRUE)
            object <- NULL
            result$theObject <- NULL
            accepted <- list(schema = "metab.da_references.v1", source = store$source, tables = list(), results = list(),
                products = list(), qvalue_warnings = result$qvalue_warnings)
            if (is.data.frame(result$da_metabolites_long))
                accepted$tables$da_metabolites_long <- store$table(as.data.frame(result$da_metabolites_long))
            for (assay in names(result$per_assay_results)) {
                table <- result$per_assay_results[[assay]]
                accepted$results[[assay]] <- list(assay = assay,
                    tables = if (is.data.frame(table)) list(da_metabolites_long = store$table(as.data.frame(table))) else list(),
                    rows = if (is.data.frame(table)) nrow(table) else 0L,
                    significant_counts = result$significant_counts[[assay]])
            }
            accepted$products$da_results <- store$product(result, "da_results", "metab-da-results")
            accepted$rows <- if (is.data.frame(result$da_metabolites_long)) nrow(result$da_metabolites_long) else 0L
            accepted$columns <- names(result$da_metabolites_long)
            accepted$model <- list(model = "limma", formula_used = parameters$formula, eBayes_trend = TRUE,
                eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))
            accepted
        })
}

#' Identify the metabolomics DA manifest.
#' @noRd
isMetabResultManifest <- function(value) {
    is.list(value) && identical(value[["schema", exact = TRUE]], "metab.da_results.v1")
}

#' Every metabolomics declaration, by name.
#' @noRd
metabolomicsWorkflowOperations <- function() {
    declarations <- list(metabIntensityFilterOperation(), metabDuplicateResolutionOperation(), metabQcFinalizeOperation(),
        metabItsdNormalisationOperation(), metabLog2Operation(), metabBetweenSampleNormalisationOperation(),
        metabRuvOperation(), metabPearsonPairsOperation(), metabPearsonSelectionOperation(), metabNormCompleteOperation(),
        metabResetOperation(), metabDaOperation())
    stats::setNames(declarations, vapply(declarations, `[[`, "", "name"))
}

metabolomicsStepLabels <- c(metab_intensity_filter = "2_Intensity_Filtered", metab_duplicate_resolution = "3_Duplicates_Resolved",
    metab_qc_finalize = "4_QC_Complete", metab_itsd_normalisation = "5_ITSD_Normalised", metab_log2 = "6_Log2",
    metab_between_sample_normalisation = "7_Normalised", metab_ruv = "8_RUV_Corrected", metab_pearson_selection = "9_Correlation_Filtered",
    metab_norm_complete = "9_Normalisation_Complete", metab_reset = "0_Reset")

# ---------------------------------------------------------------------------------------------------------------
# The wrappers: the operation on the head, its product on the checkpoint, the commit under the module's state name.
# ---------------------------------------------------------------------------------------------------------------

#' @noRd
runMetabReferenceIntensityFilter <- function(manager, percentile, proportion, config, description = sprintf(
    "Applied metabolite intensity filter (percentile: %d%%, proportion: %.2f)", as.integer(percentile), proportion),
    state_name = "metab_intensity_filtered", hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, metabIntensityFilterOperation(),
        list(metabolites_intensity_cutoff_percentile = percentile, metabolites_proportion_of_samples_below_cutoff = proportion),
        state_name, config, description, hydrate = hydrate)
}

#' @noRd
runMetabReferenceDuplicateResolution <- function(manager, config,
    description = "Resolved duplicate metabolite features by keeping highest intensity", state_name = "metab_duplicates_resolved",
    hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, metabDuplicateResolutionOperation(), list(), state_name, config, description, hydrate = hydrate)
}

#' @noRd
runMetabReferenceItsdNormalisation <- function(manager, aggregation, feature_ids, config,
    description = paste("ITSD normalization (aggregation:", aggregation, ")"), state_name = "metab_itsd_norm", hydrate = TRUE) {
    parameters <- list(itsd_aggregation = aggregation)
    parameters["itsd_feature_ids"] <- list(feature_ids)
    runWorkflowOmicsReferenceOperation(manager, metabItsdNormalisationOperation(), parameters, state_name, config, description, hydrate = hydrate)
}

#' @noRd
runMetabReferenceLog2 <- function(manager, offset, config, description = paste("Log2 transformation (offset:", offset, ")"),
    state_name = "metab_log2", hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, metabLog2Operation(), list(offset = offset), state_name, config, description, hydrate = hydrate)
}

#' @noRd
runMetabReferenceBetweenSampleNormalisation <- function(manager, method, config,
    description = paste("Between-sample normalization (method:", method, ")"), state_name = "metab_normalized", hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, metabBetweenSampleNormalisationOperation(), list(normalisation_method = method),
        state_name, config, description, hydrate = hydrate)
}

#' @param parameters The nine RUV parameters by their declared names.
#' @return The operation result with `optimisation`, the hydrated per-assay optimisation results.
#' @noRd
runMetabReferenceRuv <- function(manager, parameters, config, description = "RUV-III batch correction complete",
    state_name = "metab_ruv_corrected", hydrate = TRUE) {
    result <- runWorkflowOmicsReferenceOperation(manager, metabRuvOperation(), parameters[metabRuvParameterNames], state_name, config,
        description, hydrate = hydrate)
    result$optimisation <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
    result
}

#' The correlation filter: the pairs product, then the selection on it, on the head.
#' @return The selection's result with `correlations`, the per-assay correlation tables.
#' @noRd
runMetabReferenceCorrelation <- function(manager, group, threshold, config,
    description = paste("Correlation filtering (threshold:", threshold, ")"), state_name = "metab_correlation_filtered",
    hydrate = TRUE) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        pairs <- runSavedOperation(manager, metabPearsonPairsOperation(), source, list(correlation_group = group))
        result <- runWorkflowOmicsReferenceOperation(manager, metabPearsonSelectionOperation(),
            list(pairs = pairs$reference$entity_id, min_pearson_correlation_threshold = threshold, correlations = pairs$value),
            state_name, config, description, source = source, dependencies = list(pairs$reference), hydrate = hydrate)
        result$correlations <- pairs$value
        result
    })
}

#' The design step on the owner: `design_builder` provenance with the three design roles, and the inputs (the import's
#' source copies under `import_sources/` as `metab_assay_<n>`, the design file).
#' @noRd
recordMetabDesignProvenance <- function(manager, object, source_dir = NULL, design_file = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
    setWorkflowStepProvenance(manager, "design_builder",
        parameter_names = c("sample_id", "group_id", "technical_replicate_id"),
        values = list(sample_id = object@sample_id, group_id = object@group_id,
            technical_replicate_id = object@technical_replicate_id))
    paths <- character()
    if (!is.null(source_dir) && dir.exists(file.path(source_dir, "import_sources"))) {
        sources <- sort(list.files(file.path(source_dir, "import_sources"), full.names = TRUE))
        if (length(sources)) paths[paste0("metab_assay_", seq_along(sources))] <- sources
    }
    if (is.null(design_file) && !is.null(source_dir)) design_file <- file.path(source_dir, "design_matrix.tab")
    if (!is.null(design_file) && file.exists(design_file)) paths[["design_file"]] <- design_file
    if (length(paths)) recordDIAWorkflowInputs(manager, paths)
    invisible(NULL)
}

#' Run the metabolomics DA declaration on a state and accept its manifest on the head (`metab.da_results.v1`).
#' @return `list(manifest, results, reference)`: `results` the kernel's list for the module's display and files, with
#'   `theObject` the hydrated source.
#' @noRd
runMetabReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the metabolomics DA needs the session's state owner", class = "saved_entity_scope")
    if (!isSavedWorkflowReference(source)) cli::cli_abort("the metabolomics DA needs a saved state reference", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        spec <- readWorkflowStateMetadata(manager, source)$data
        contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
        parameters <- list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh,
            lfc = settings$treat_lfc_cutoff, contrasts = contrasts)
        run <- runSavedOperation(manager, metabDaOperation(), source, parameters)
        if (identical(run$status, "computed")) {
            snapshot <- context$snapshot
            snapshot$da_dependencies[[run$reference$entity_id]] <- collectDIANativeReferences(run$value)
            publishDIAStateHead(context, snapshot)
        }
        clearDIAStepProvenance(manager)
        manifest <- list(schema = "metab.da_results.v1", source = source, metabolite_id_column = spec$slots$metabolite_id_column,
            sample_id = spec$slots$sample_id, contrasts = contrasts, results = run$value$results, tables = run$value$tables,
            products = run$value$products, parameters = parameters[c("formula", "qvalue", "lfc")],
            qvalue_warnings = run$value$qvalue_warnings, model = run$value$model)
        reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot <- context$snapshot
        snapshot$da_results <- reference
        snapshot$da_dependencies[[reference$entity_id]] <- collectDIANativeReferences(list(source, run$value))
        publishDIAStateHead(context, snapshot)
        results <- hydrateSavedWorkflowEntity(context$store, run$value$products$da_results)$value
        results$theObject <- readWorkflowExecutionComponents(manager, source, fields = metabolomicsOperationFields)
        list(manifest = manifest, results = results, reference = reference)
    })
}

#' Replay one recorded metabolomics step through the wrapper the module used, with the recorded values.
#' @noRd
replayMetabolomicsStep <- function(manager, step, parameters) {
    operation <- step$operation
    config <- readWorkflowStateMetadata(manager)$config
    description <- step$description %||% operation
    state <- step$state
    unsupported <- function() cli::cli_abort("step {.val {state}} cannot be replayed: operation {.val {operation}}",
        class = "workflow_replay_unsupported")
    switch(operation,
        metab_intensity_filter = runMetabReferenceIntensityFilter(manager, parameters$metabolites_intensity_cutoff_percentile,
            parameters$metabolites_proportion_of_samples_below_cutoff, config, description, state_name = state, hydrate = FALSE),
        metab_duplicate_resolution = runMetabReferenceDuplicateResolution(manager, config, description, state_name = state, hydrate = FALSE),
        metab_qc_finalize = , metab_norm_complete = commitWorkflowOmicsStoreStep(manager, operation, state, config, description),
        metab_itsd_normalisation = runMetabReferenceItsdNormalisation(manager, parameters$itsd_aggregation,
            parameters$itsd_feature_ids, config, description, state_name = state, hydrate = FALSE),
        metab_log2 = runMetabReferenceLog2(manager, parameters$offset, config, description, state_name = state, hydrate = FALSE),
        metab_between_sample_normalisation = runMetabReferenceBetweenSampleNormalisation(manager, parameters$normalisation_method,
            config, description, state_name = state, hydrate = FALSE),
        metab_ruv = runMetabReferenceRuv(manager, parameters, config, description, state_name = state, hydrate = FALSE),
        metab_pearson_selection = {
            previous <- Filter(function(p) identical(p$operation, "metab_pearson_pairs"), step$preceded_by %||% list())
            if (!length(previous)) unsupported()
            pairs <- workflowReplayParameters(previous[[1L]]$parameters, workflowOperationEntry("metab_pearson_pairs")$parameters)
            runMetabReferenceCorrelation(manager, pairs$correlation_group, parameters$min_pearson_correlation_threshold, config,
                description, state_name = state, hydrate = FALSE)
        },
        unsupported())
    invisible(NULL)
}
