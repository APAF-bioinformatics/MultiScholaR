test_that("physical filter and grouped consumers preserve inclusion and exclusion science", {
    harness <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools/profiling/workflow_checkpoint_dia_driver.R"), harness)
    root <- withr::local_tempdir()
    fixture <- harness$preparePerf107GuiInputs(file.path(root, "input"))
    report <- as.data.frame(vroom::vroom(fixture$files$report, show_col_types = FALSE,
        progress = FALSE))
    design <- as.data.frame(vroom::vroom(fixture$files$design, show_col_types = FALSE,
        progress = FALSE))
    design$Run <- design$sample
    design$tech_rep_group <- paste(design$group, design$replicates, sep = "_")
    design$omit <- design$Run == "b3"
    object <- PeptideQuantitativeDataDiann(report, design,
        technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "DIA")))
    object <- rollUpPrecursorToPeptide(object, core_utilisation = NA)
    parameters <- list(peptides_per_sample_cutoff = 100,
        inclusion_list = c("b3", "b3", "absent", "", NA_character_))
    for (exclude in c(FALSE, TRUE)) {
        rlang::local_bindings(config_list = list(
            removePeptidesWithOnlyOneReplicate = list(replicate_group_column = "replicates"),
            peptideMissingValueImputation = list(proportion_missing_values = 0.5,
                exclusion_column = if (exclude) "omit" else NULL)
        ), .env = .GlobalEnv)
        native_filter <- filterMinNumPeptidesPerSample(object, 100,
            core_utilisation = NA, inclusion_list = parameters$inclusion_list)
        native_replicate <- removePeptidesWithOnlyOneReplicate(native_filter,
            replicate_group_column = "replicates", core_utilisation = NA)
        native_imputed <- peptideMissingValueImputation(native_replicate,
            proportion_missing_values = 0.5, core_utilisation = NA,
            exclusion_column = if (exclude) "omit" else NULL)
        expect_identical(nrow(native_filter@peptide_data), 860L)
        expect_identical(nrow(native_replicate@peptide_data), 859L)
        expect_identical(nrow(native_imputed@peptide_data), if (exclude) 860L else 960L)
        store_root <- file.path(root, paste0("store-", exclude)); dir.create(store_root)
        manager <- newDIAStateAdapter(newWorkflowEntityStore(store_root,
            basename(store_root), "DIA"), WorkflowState$new(audit_enabled = FALSE))
        withr::defer(closeDIAStateAdapter(manager))
        manager$.context$reference_states <- TRUE
        manager$saveState("input", object, list(), "input")
        workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
        result <- computeDIAReferenceSampleFilter(manager, getDIAStateReference(manager),
            parameters, output_mode = "physical")
        finishDIAReferencePeptideStep(workflow, result$checkpoint, result$state_name,
            result$description, result$config, result$section, result$counts,
            observed_inclusion_absence = result$summary$inclusion_entries_absent_from_data)
        runDIAReferencePeptideGroups(workflow, "removePeptidesWithOnlyOneReplicate",
            list(replicate_group_column = "replicates"), "replicate_filtered", "replicates")
        runDIAReferencePeptideGroups(workflow, "peptideMissingValueImputation",
            list(proportion_missing_values = 0.5, exclusion_column = if (exclude) "omit" else NULL),
            "imputed", "imputation")
        expect_identical(diaStateOwnership(manager)$full_s4_reads, 0L)
        expect_identical(manager$.context$store$open_resources, 0L)
        for (entry in list(list("sample_filtered", native_filter, "filterMinNumPeptidesPerSample"),
            list("replicate_filtered", native_replicate, "removePeptidesWithOnlyOneReplicate"),
            list("imputed", native_imputed, "peptideMissingValueImputation"))) {
            expect_true(harness$verifyPerf107ScientificState(manager$getState(entry[[1L]]),
                entry[[2L]], entry[[3L]]))
        }
        # A successful filter receipt is insufficient if a downstream consumer
        # keeps the old design or presents the pre-imputation state as complete.
        actual_imputed <- manager$getState("imputed")
        stale_design <- actual_imputed
        stale_design@design_matrix <- object@design_matrix
        expect_error(harness$verifyPerf107ScientificState(stale_design,
            native_imputed, "peptideMissingValueImputation"), "design_matrix")
        skipped_imputation <- manager$getState("replicate_filtered")
        expect_false(identical(skipped_imputation@peptide_data,
            native_imputed@peptide_data))
        expect_error(harness$verifyPerf107ScientificState(skipped_imputation,
            native_imputed, "peptideMissingValueImputation"), "peptide_data")
        # Check the positive again so a detector that rejects every input fails.
        expect_true(harness$verifyPerf107ScientificState(actual_imputed,
            native_imputed, "peptideMissingValueImputation"))
        closeDIAStateAdapter(manager)
    }
})
