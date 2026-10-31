# KIT-031: the lipidomics declarations, the kit's completeness test. Every step of the lipidomics route is a registered
# in-process declaration with the unit of the step map; on a real owner in a temporary project the wrappers the
# modules call commit each state as a component checkpoint with its operation's provenance, the analysis record is
# complete (no legacy step, no missing parameter, the inputs recorded), the DA manifest `lipid.da_results.v1` is the
# head's and its display list equals the kernel's direct output, the report overlay reads the lipid variables, and
# the record replays exact state by state on a copy of the store. The fixture is helper-workflow-lipid-declarations.R.
library(testthat)

test_that("KIT-031: the lipidomics declarations are registered in process with the units, merges and statistics of the step map", {
    units <- c(lipid_intensity_filter = "row", lipid_duplicate_resolution = "feature_group", lipid_qc_finalize = "store",
        lipid_itsd_normalisation = "complete_sample", lipid_log2 = "row", lipid_between_sample_normalisation = "whole_matrix",
        lipid_ruv = "whole_matrix", lipid_pearson_pairs = "whole_matrix", lipid_pearson_selection = "whole_matrix",
        lipid_norm_complete = "store", lipid_reset = "store", lipid_da_output = "whole_matrix")
    declarations <- lipidomicsWorkflowOperations()
    expect_setequal(names(declarations), names(units))
    for (name in names(units)) {
        declaration <- declarations[[name]]
        expect_identical(declaration$unit, units[[name]], info = name)
        expect_identical(declaration$input_class, "LipidomicsAssayData", info = name)
        expect_identical(declaration$fields, c("lipid_data", "design_matrix"), info = name)
        partitioned <- units[[name]] %in% c("row", "complete_sample", "feature_group")
        expect_identical(declaration$merge, if (partitioned) c(lipid_data = "concatenate", design_matrix = "design") else NULL, info = name)
        entry <- workflowActionSpec(name)
        expect_identical(entry$compute, "in_process", info = name)
        expect_identical(entry$omic_type, "lipidomics", info = name)
        expect_null(entry$handler, info = name)
        expect_identical(entry$variants, "whole", info = name)
        expect_error(assertWorkflowOperationDispatchable(name), class = "workflow_execution_unsupported")
        expect_true(all(entry$parameters %in% c(names(workflowOperationDefaults[[name]]), names(entry$declaration$fixed))), info = name)
        expect_true(all(entry$parameters %in% names(workflowParameterLabels[[name]])), info = name)
    }
    # The pre-pass statistics the step map records: the per-assay quantile and the mean standard factor.
    expect_named(declarations$lipid_intensity_filter$global_statistics, "min_lipid_intensity_threshold")
    expect_named(declarations$lipid_itsd_normalisation$global_statistics, "mean_standard_factor")
    expect_identical(declarations$lipid_itsd_normalisation$fixed, list(method = "ITSD"))
    expect_identical(declarations$lipid_da_output$fixed, list(eBayes_trend = TRUE, eBayes_robust = TRUE))
    expect_null(declarations$lipid_da_output$output_class)
    expect_null(declarations$lipid_pearson_pairs$output_class)
    expect_identical(workflowActionSpec("lipid_intensity_filter")$step_label, "2_Intensity_Filtered")
    expect_null(workflowActionSpec("lipid_da_output")$step_label)
    # The route's sections name its declarations; the proteomics map is untouched.
    expect_identical(workflowOperationForSection("lipidIntensityFilteringFn", "LIPIDOMICS"), "lipid_intensity_filter")
    expect_identical(workflowOperationForSection("ruvIII_C_Varying", "LIPIDOMICS"), "lipid_ruv")
    expect_length(workflowOperationForSection("lipidIntensityFilteringFn"), 0L)
    expect_identical(workflowOperationForSection("chooseBestProteinAccession"), "protein_accession_cleanup")
    # The DA manifest is declared, with the operation and the long table the record reads.
    output <- workflowOutputDeclaration("lipid.da_results.v1")
    expect_identical(output$operation, "lipid_da_output")
    expect_identical(output$long_table, "da_lipids_long")
    expect_identical(output$identifies, "isLipidResultManifest")
    expect_true(isLipidResultManifest(list(schema = "lipid.da_results.v1")))
    expect_false(isLipidResultManifest(list(schema = "dia.da_results.v1")))
    expect_identical(workflowOutputDeclaration("dia.da_results.v1")$long_table, "da_proteins_long")
})

test_that("KIT-031: the wrappers refuse a manager that is not the session owner and the store steps leave a double alone", {
    double <- list2env(list(getState = function() NULL, saveState = function(...) invisible(TRUE)))
    expect_error(runLipidReferenceIntensityFilter(double, 1, 0.5, list()), class = "saved_entity_scope")
    expect_error(runLipidReferenceCorrelation(double, "group", 0.5, list()), class = "saved_entity_scope")
    expect_error(runLipidReferenceDaAnalysis(double, NULL, list(), data.frame()), class = "saved_entity_scope")
    expect_null(markLipidReferenceStoreStep(double, "lipid_qc_finalize"))
    expect_null(recordLipidDesignProvenance(double, NULL))
})

test_that("KIT-031: the lipid journey on the owner commits every state with its declaration, records completely, resolves its DA manifest, overlays the report and replays exact", {
    skip_if_not(workflowOwnerCapable())
    root <- withr::local_tempdir()
    fixture <- lipidDeclarationsOwner(root)
    owner <- fixture$owner
    withr::defer(closeDIAStateAdapter(owner))
    journey <- suppressWarnings(suppressMessages(lipidDeclarationsJourney(fixture)))

    # Every commit is a component checkpoint attributed to its declaration (the design step to `design_builder`).
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
    # Each state carries its operation's product; the kernels ran as the modules call them (the `@args` keys).
    for (state in c("lipid_intensity_filtered", "lipid_ruv_corrected", "lipid_correlation_filtered")) {
        products <- readWorkflowStateMetadata(owner, snapshot$states[[state]])$data$products
        expect_true(any(grepl("_science$", names(products))), info = state)
    }
    args <- readWorkflowStateMetadata(owner, snapshot$states$lipid_intensity_filtered)$data$slots$args
    expect_identical(args$lipidIntensityFilteringFn$lipids_intensity_cutoff_percentile, 1)
    # The intensity filter counted its rows; the duplicates were resolved (two duplicate ids per assay in the fixture).
    expect_identical(unname(journey$duplicates$value$counts$before - journey$duplicates$value$counts$after), c(2L, 2L))
    expect_identical(unname(unlist(journey$ruv$value$best_k)), c(1L, 1L))
    expect_true(all(unlist(journey$ruv$value$control_count) > 0L))
    expect_length(journey$ruv$optimisation, 2L)

    # The record: valid, no legacy step, no missing declared parameter, the inputs recorded with their hashes.
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

    # The DA manifest is the head's and its display list equals the kernel's direct output on the same state.
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
    long <- materializeWorkflowTable(owner$.context, manifest$results$LCMS_Pos$tables$da_lipids_long)
    expect_identical(nrow(long), manifest$results$LCMS_Pos$rows)
    # A second run on the same state and settings is a hit: the manifest's tables are the same entities.
    again <- suppressWarnings(suppressMessages(runLipidReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), fixture$workflow$contrasts_tbl)))
    expect_identical(again$manifest$results$LCMS_Pos$tables$da_lipids_long$entity_id, manifest$results$LCMS_Pos$tables$da_lipids_long$entity_id)

    # The report overlay: the lipid template's variables take the record's values; the proteomics-only ones are left.
    json <- file.path(root, "analysis_record.json")
    writeWorkflowAnalysisRecord(owner, json)
    expect_true(validateWorkflowAnalysisRecord(jsonlite::read_json(json, simplifyVector = FALSE))$valid)
    overlay <- workflowRecordReportParameters(json, list(itsd_aggregation = "sum", normalisation_method = "none", best_k = 0,
        num_neg_ctrl = 0L, s4_de_q_val_thresh = 0.1, intensity_cutoff_percentile = 5, log_transform_offset = 2,
        ruv_applied = FALSE, itsd_applied = FALSE, qvalue_threshold = 0.01, ruv_grouping_variable = "x"), list())
    expect_identical(overlay$itsd_aggregation, "median")
    expect_true(overlay$itsd_applied)
    expect_identical(overlay$normalisation_method, "cyclicloess")
    expect_identical(overlay$best_k, 2L)
    expect_true(overlay$ruv_applied)
    expect_identical(overlay$ruv_grouping_variable, "group")
    expect_identical(overlay$s4_de_q_val_thresh, 0.05)
    expect_identical(overlay$intensity_cutoff_percentile, 1)
    expect_identical(overlay$log_transform_offset, 1)
    expect_identical(overlay$qvalue_threshold, 0.01)
    expect_true(overlay$record_present)

    # The replay on a copy of the store cut to the design state: every state hash-identical, the DA re-run.
    copy <- withr::local_tempdir()
    stopifnot(system2("cp", c("-a", file.path(fixture$source_dir, "."), copy)) == 0L)
    harness <- new.env(parent = globalenv())
    sys.source(file.path(testthat::test_path("..", ".."), "tools", "profiling", "workflow_kit.R"), harness)
    api <- function(name) get(name, asNamespace("MultiScholaR"))
    store_root <- dirname(list.files(fixture$source_dir, "^store[.]json$", recursive = TRUE, all.files = TRUE, full.names = TRUE)[[1L]])
    store_root <- dirname(store_root)
    copy_root <- file.path(copy, substring(store_root, nchar(fixture$source_dir) + 2L))
    prepared <- harness$kitPrepareReplayProject(copy_root, c("initial", "lipid_raw_data_s4"), api)
    inputs <- c(lipid_assay_1 = file.path(fixture$source_dir, "import_sources", "LCMS_Neg_export.txt"),
        lipid_assay_2 = file.path(fixture$source_dir, "import_sources", "LCMS_Pos_export.txt"),
        design_file = file.path(fixture$source_dir, "design_matrix.tab"))
    scratch <- withr::local_tempdir()
    record_read <- jsonlite::read_json(json, simplifyVector = FALSE)
    replay <- suppressWarnings(suppressMessages(harness$kitReplayJourney(record_read, prepared, store_root, inputs, fixture$paths, scratch, api)))
    expect_identical(replay$states, history)
    expect_identical(as.character(replay$mismatches), character())
    expect_true(any(vapply(replay$per_step, function(x) identical(x$operation, "lipid_da_output"), logical(1))))
    # A record whose input file changed is refused before any step runs.
    tampered <- inputs
    tampered[["design_file"]] <- file.path(fixture$source_dir, "import_sources", "LCMS_Pos_export.txt")
    expect_error(replayWorkflowAnalysisRecord(record_read, prepared, tampered, fixture$paths, scratch), class = "workflow_replay_inputs")
})
