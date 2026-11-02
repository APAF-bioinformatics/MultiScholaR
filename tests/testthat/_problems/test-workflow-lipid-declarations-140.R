# Extracted from test-workflow-lipid-declarations.R:140

# prequel ----------------------------------------------------------------------
library(testthat)

# test -------------------------------------------------------------------------
skip_if_not(workflowOwnerCapable())
root <- withr::local_tempdir()
fixture <- lipidDeclarationsOwner(root)
owner <- fixture$owner
withr::defer(closeDIAStateAdapter(owner))
journey <- suppressWarnings(suppressMessages(lipidDeclarationsJourney(fixture)))
snapshot <- owner$.context$snapshot
history <- unlist(snapshot$state_history)
expect_identical(history, c("initial", "lipid_raw_data_s4", "lipid_intensity_filtered", "lipid_duplicates_resolved",
        "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized", "lipid_ruv_corrected", "lipid_correlation_filtered"))
operations <- vapply(snapshot$step_provenance[history], function(p) p$operation, "")
expect_identical(unname(operations), c("import", "design_builder", "lipid_intensity_filter", "lipid_duplicate_resolution",
        "lipid_qc_finalize", "lipid_itsd_normalisation", "lipid_log2", "lipid_between_sample_normalisation", "lipid_ruv",
        "lipid_pearson_selection"))
expect_identical(snapshot$step_provenance$lipid_correlation_filtered$preceded_by[[1L]]$operation, "lipid_pearson_pairs")
expect_identical(snapshot$step_provenance$lipid_intensity_filtered$values$lipids_intensity_cutoff_percentile, 1)
expect_true("itsd_feature_ids" %in% names(snapshot$step_provenance$lipid_itsd_norm$values))
expect_null(snapshot$step_provenance$lipid_itsd_norm$values$itsd_feature_ids)
for (state in history[-1L]) expect_true(isDIAComponentCheckpoint(readWorkflowStateMetadata(owner, snapshot$states[[state]])$data), info = state)
for (state in c("lipid_intensity_filtered", "lipid_ruv_corrected", "lipid_correlation_filtered")) {
        products <- readWorkflowStateMetadata(owner, snapshot$states[[state]])$data$products
        expect_true(any(grepl("_science$", names(products))), info = state)
    }
args <- readWorkflowStateMetadata(owner, snapshot$states$lipid_intensity_filtered)$data$slots$args
expect_identical(args$lipidIntensityFilteringFn$lipids_intensity_cutoff_percentile, 1)
expect_identical(unname(journey$duplicates$value$counts$before - journey$duplicates$value$counts$after), c(2L, 2L))
expect_identical(unname(unlist(journey$ruv$value$best_k)), c(1L, 1L))
expect_true(all(unlist(journey$ruv$value$control_count) > 0L))
expect_length(journey$ruv$optimisation, 2L)
record <- buildWorkflowAnalysisRecord(owner)
check <- validateWorkflowAnalysisRecord(record)
expect_true(check$valid)
expect_identical(check$incomplete, "external resources: controller-process queries only")
expect_identical(record$omic_type, "lipidomics")
expect_identical(vapply(record$steps, function(s) s$operation, ""), unname(operations))
roles <- vapply(record$inputs, function(i) i$role, "")
expect_setequal(roles, c("lipid_assay_1", "lipid_assay_2", "design_file"))
expect_true(all(nzchar(vapply(record$inputs, function(i) i$sha256, ""))))
by_state <- stats::setNames(record$steps, vapply(record$steps, `[[`, "", "state"))
named <- function(parameters) stats::setNames(parameters, vapply(parameters, `[[`, "", "name"))
intensity <- named(by_state$lipid_intensity_filtered$parameters)
expect_identical(intensity$lipids_intensity_cutoff_percentile$source, "args")
expect_identical(intensity$lipids_intensity_cutoff_percentile$provenance, "no_shipped_default")
ruv <- named(by_state$lipid_ruv_corrected$parameters)
expect_identical(decodeWorkflowRecordValue(ruv$manual_k$value), 1L)
expect_identical(by_state$lipid_ruv_corrected$outcomes$k, 2L)
expect_identical(length(by_state$lipid_ruv_corrected$outcomes$k_per_assay), 2L)
expect_true(by_state$lipid_itsd_norm$outcomes$itsd_applied)
expect_identical(by_state$lipid_itsd_norm$step_label, "5_ITSD_Normalised")
da <- record$analyses$differential_abundance
expect_identical(da$operation, "lipid_da_output")
expect_identical(da$source_state, "lipid_correlation_filtered")
expect_identical(vapply(da$results, function(r) r$contrast, ""), c("LCMS_Pos", "LCMS_Neg"))
expect_true(all(vapply(da$results, function(r) r$rows, 1L) > 0L))
expect_identical(named(da$parameters)$qvalue$provenance, "default")
expect_identical(da$model$model, "limma")
manifest <- resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results)
expect_true(isLipidResultManifest(manifest))
expect_identical(manifest$source$entity_id, snapshot$states$lipid_correlation_filtered$entity_id)
expect_null(owner$.context$pending_provenance)
direct <- suppressWarnings(suppressMessages(runLipidsDA(theObject = owner$getState("lipid_correlation_filtered"),
        contrasts_tbl = fixture$workflow$contrasts_tbl, formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0, eBayes_trend = TRUE, eBayes_robust = TRUE)))
expect_equal(as.data.frame(journey$da$results$da_lipids_long), as.data.frame(direct$da_lipids_long))
expect_equal(journey$da$results$per_assay_results, direct$per_assay_results)
expect_identical(journey$da$results$significant_counts, direct$significant_counts)
expect_s4_class(journey$da$results$theObject, "LipidomicsAssayData")
