# KIT-015: the protein QC declarations on the fixtures: each compute-only entry point against the module's
# kernel path, the extracted dedup kernel against the module's inline code, IQ's per-group check, the worker
# entry point and the controller's completion for a protein operation with a peptide input.

proteinDeclarationFixtureManager <- function(object, root, state) {
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState(state, object, list(), "fixture")
    manager
}

imputedPeptideFixture <- function() {
    precursor <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- filterMinNumPeptidesPerSample(rollUpPrecursorToPeptide(precursor, core_utilisation = NA),
        peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    rlang::local_bindings(config_list = list(peptideMissingValueImputation = list(proportion_missing_values = 0.5)), .env = .GlobalEnv)
    peptideMissingValueImputation(rolled, proportion_missing_values = 0.5, core_utilisation = NA)
}

test_that("the four protein QC declarations validate with the units of the step map and the whole placement", {
    expect_identical(diaProteinIQRollupOperation()$unit, "feature_group")
    expect_identical(diaProteinAccessionCleanupOperation()$unit, "row")
    expect_identical(diaProteinIntensityFilterOperation()$unit, "whole_matrix")
    expect_identical(diaProteinDuplicateRemovalOperation()$unit, "feature_group")
    expect_null(diaProteinIntensityFilterOperation()$merge)
    expect_identical(unname(diaProteinDuplicateRemovalOperation()$merge[["counts"]]), "sum")
    for (name in c("protein_iq_rollup", "protein_accession_cleanup", "protein_intensity_filter", "protein_duplicate_removal")) {
        entry <- workflowActionSpec(name)
        expect_identical(entry$variants, "whole", info = name)
        expect_identical(entry$compute, "protein")
        expect_identical(entry$completion, "protein")
        expect_true(is.function(get(entry$handler, mode = "function")), info = name)
        expect_true(is.function(get(entry$validate, mode = "function")), info = name)
    }
    expect_identical(workflowActionSpec("protein_iq_rollup")$input, "peptide")
    expect_identical(workflowActionSpec("protein_replicate_filter")$default_variant, "rows")
})

test_that("the dedup kernel is the module's inline code: identical table and duplicates on a table with duplicates", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    table <- object@protein_quant_table
    table <- rbind(table, table[c(1L, 1L, 5L), , drop = FALSE])
    table[nrow(table), 2L] <- NA
    expected_duplicates <- table |>
        dplyr::group_by(.data[[object@protein_id_column]]) |>
        dplyr::filter(dplyr::n() > 1) |>
        dplyr::select(dplyr::all_of(object@protein_id_column)) |>
        dplyr::distinct() |>
        dplyr::pull(dplyr::all_of(object@protein_id_column))
    numeric_columns <- setdiff(colnames(table), object@protein_id_column)
    expected <- table |>
        dplyr::group_by(!!rlang::sym(object@protein_id_column)) |>
        dplyr::summarise(dplyr::across(dplyr::all_of(numeric_columns), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
        dplyr::ungroup()
    resolved <- resolveProteinDuplicateRows(table, object@protein_id_column, mean)
    expect_identical(resolved$duplicates, expected_duplicates)
    expect_identical(resolved$table, expected)
    expect_identical(sort(resolved$duplicates), sort(unique(table[[object@protein_id_column]][c(1L, 5L)])))
    # The module's in-memory path on the same object gives the same table.
    object@protein_quant_table <- table
    state <- WorkflowState$new(audit_enabled = FALSE)
    state$saveState("protein_intensity_filtered", object, list(), "fixture")
    workflow <- list2env(list(state_manager = state, qc_params = list()))
    module <- runProteinDuplicateRemovalStep(workflow, "mean", logInfoFn = function(...) NULL)
    expect_identical(module$deduplicatedS4@protein_quant_table, expected)
    expect_identical(module$duplicates, expected_duplicates)
})

test_that("the duplicate removal and the intensity filter declarations reproduce the module's kernel path on the saved fixture", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    table <- object@protein_quant_table
    object@protein_quant_table <- rbind(table, table[c(2L, 2L), , drop = FALSE])
    rlang::local_bindings(config_list = list(removeRowsWithMissingValuesPercent = list()), .env = .GlobalEnv)
    # Duplicate removal.
    manager <- proteinDeclarationFixtureManager(object, withr::local_tempdir(), "protein_intensity_filtered")
    withr::defer(closeDIAStateAdapter(manager))
    value <- runProjected(manager, diaProteinDuplicateRemovalOperation(), getDIAStateReference(manager), list(aggregation_method = "mean"))
    expect_identical(value$state_name, "duplicates_removed")
    expect_identical(value$execution$variant, "whole")
    expect_identical(value$config$num_duplicates, 1L)
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    result <- finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
    expect_identical(manager$current_state, "duplicates_removed")
    expect_true(isSavedWorkflowReference(result$stateReference))
    expect_identical(result$duplicates, value$duplicates)
    native <- resolveProteinDuplicateRows(object@protein_quant_table, object@protein_id_column, mean)$table
    expect_equal(as.data.frame(manager$getState("duplicates_removed")@protein_quant_table), as.data.frame(native), ignore_attr = TRUE)
    # Intensity filter (flexible), against the kernel on the in-memory object with the same parameter updates.
    manager2 <- proteinDeclarationFixtureManager(object, withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(manager2))
    parameters <- list(strict_mode = FALSE, min_reps_per_group = 1L, min_groups = 1L, proteins_intensity_cutoff_percentile = 5)
    validateWorkflowActionParameters("protein_intensity_filter", parameters, workflowExecutionPolicy(FALSE))
    value2 <- runProjected(manager2, diaProteinIntensityFilterOperation(), getDIAStateReference(manager2), parameters)
    expect_identical(value2$state_name, "protein_intensity_filtered")
    expect_identical(value2$qc_name, "intensity_filter")
    expect_false(isTRUE(value2$config$strict_mode))
    expect_identical(value2$config$proteins_intensity_cutoff_percentile, 5)
    workflow2 <- list2env(list(state_manager = manager2, qc_params = list(protein_qc = list())))
    result2 <- finishDIAReferenceProteinStep(workflow2, value2, publish = FALSE)
    expect_identical(manager2$current_state, "protein_intensity_filtered")
    expect_identical(names(workflow2$qc_params$protein_qc), "intensity_filter")
    kernel_object <- updateMissingValueParameters(object, min_reps_per_group = 1L, min_groups = 1L)
    kernel_object <- updateConfigParameter(kernel_object, function_name = "removeRowsWithMissingValuesPercent",
        parameter_name = "proteins_intensity_cutoff_percentile", new_value = 5)
    removeRowsWithMissingValuesPercentFn <- removeRowsWithMissingValuesPercent
    native2 <- removeRowsWithMissingValuesPercentFn(kernel_object)
    committed <- manager2$getState("protein_intensity_filtered")
    expect_equal(as.data.frame(committed@protein_quant_table), as.data.frame(native2@protein_quant_table), ignore_attr = TRUE)
    expect_identical(committed@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
        native2@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff)
    expect_true(grepl(paste("Proteins remaining:", value2$protein_count), result2$resultText, fixed = TRUE))
})

test_that("the accession cleanup declaration skips without FASTA data as the module does and records the skip", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- proteinDeclarationFixtureManager(object, withr::local_tempdir(), "protein_s4_created")
    withr::defer(closeDIAStateAdapter(manager))
    parameters <- list(delimiter = ";", aggregation_method = "mean")
    validateWorkflowActionParameters("protein_accession_cleanup", parameters, workflowExecutionPolicy(FALSE))
    expect_error(validateWorkflowActionParameters("protein_accession_cleanup", list(delimiter = ";"), workflowExecutionPolicy(FALSE)))
    value <- runProjected(manager, diaProteinAccessionCleanupOperation(), getDIAStateReference(manager), parameters)
    expect_identical(value$state_name, "protein_accession_cleaned")
    expect_false(value$config$cleanup_applied)
    expect_identical(value$description, "Skipped accession cleanup (no FASTA data)")
    expect_identical(value$qc_name, "accession_cleanup")
    expect_false(value$accession_cleanup_results$cleanup_applied)
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    result <- finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
    expect_identical(manager$current_state, "protein_accession_cleaned")
    expect_false(result$cleanupApplied)
    expect_identical(workflow$accession_cleanup_results$proteins_after, value$protein_count)
    expect_equal(as.data.frame(manager$getState("protein_accession_cleaned")@protein_quant_table),
        as.data.frame(object@protein_quant_table), ignore_attr = TRUE)
    # The module's in-memory path without FASTA data gives the same skip record.
    state <- WorkflowState$new(audit_enabled = FALSE)
    state$saveState("protein_s4_created", object, list(), "fixture")
    bus <- list2env(list(state_manager = state, qc_params = list(), aa_seq_tbl_final = NULL, fasta_metadata = NULL))
    module <- runProteinAccessionCleanupStep(bus, ";", "mean", logInfoFn = function(...) NULL, logWarnFn = function(...) NULL)
    expect_false(module$cleanupApplied)
    expect_identical(bus$accession_cleanup_results$proteins_after, value$protein_count)
})

test_that("the IQ rollup declaration reproduces the module's reference path on the imputed fixture, and IQ per protein group equals IQ on the whole table (partition check)", {
    skip_if_not_installed("iq")
    object <- imputedPeptideFixture()
    manager <- proteinDeclarationFixtureManager(object, withr::local_tempdir(), "imputed")
    withr::defer(closeDIAStateAdapter(manager))
    work <- withr::local_tempdir()
    parameters <- list(peptide_qc_dir = file.path(work, "peptide"), protein_qc_dir = file.path(work, "protein"))
    validateWorkflowActionParameters("protein_iq_rollup", parameters, workflowExecutionPolicy(FALSE))
    value <- runProjected(manager, diaProteinIQRollupOperation(), getDIAStateReference(manager), parameters)
    expect_identical(value$state_name, "protein_s4_created")
    expect_identical(value$checkpoint$class, "ProteinQuantitativeData")
    expect_identical(basename(value$config$iq_output_file), "iq_output_file.txt")
    expect_identical(value$section, "protein_rollup")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
    result <- finishDIAReferenceProteinStep(workflow, value, publish = FALSE)
    expect_identical(manager$current_state, "protein_s4_created")
    declared <- manager$getState("protein_s4_created")
    expect_s4_class(declared, "ProteinQuantitativeData")
    expect_identical(nrow(declared@protein_quant_table), value$protein_count)
    # The module's reference path on a second store of the same fixture.
    manager2 <- proteinDeclarationFixtureManager(object, withr::local_tempdir(), "imputed")
    withr::defer(closeDIAStateAdapter(manager2))
    project <- withr::local_tempdir()
    paths <- list(source_dir = project, protein_qc_dir = file.path(project, "protein_qc"), peptide_qc_dir = file.path(project, "peptide_qc"))
    for (path in paths) dir.create(path, showWarnings = FALSE)
    manager2$.context$experiment_paths <- paths
    bus <- list2env(list(state_manager = manager2, qc_params = list()))
    module <- runProteinIqRollupApplyStep(bus, paths, showNotificationFn = function(...) NULL,
        logInfoFn = function(...) NULL, logWarnFn = function(...) NULL)
    expect_true(isSavedWorkflowReference(module$stateReference))
    reference <- manager2$getState("protein_s4_created")
    expect_equal(as.data.frame(declared@protein_quant_table), as.data.frame(reference@protein_quant_table), ignore_attr = TRUE)
    expect_equal(as.data.frame(declared@protein_id_table), as.data.frame(reference@protein_id_table), ignore_attr = TRUE)
    expect_equal(as.data.frame(declared@design_matrix), as.data.frame(reference@design_matrix), ignore_attr = TRUE)
    # Partition check (ticket: IQ's unit is checked, never assumed): IQ on each protein group's rows alone equals
    # its row in the whole-table output, with the same column set when the group covers every sample.
    input <- file.path(parameters$peptide_qc_dir, "peptide_values_imputed.tsv")
    data <- vroom::vroom(input, show_col_types = FALSE, .name_repair = "minimal")
    whole <- vroom::vroom(value$config$iq_output_file, show_col_types = FALSE, .name_repair = "minimal")
    protein <- object@protein_id_column
    groups <- head(unique(data[[protein]]), 5L)
    for (group in groups) {
        part_in <- file.path(work, paste0("part-", make.names(group), ".tsv"))
        part_out <- file.path(work, paste0("part-", make.names(group), "-out.txt"))
        vroom::vroom_write(data[data[[protein]] == group, , drop = FALSE], part_in)
        iq::process_long_format(part_in, output_filename = part_out, sample_id = object@sample_id, primary_id = protein,
            secondary_id = object@peptide_sequence_column, intensity_col = "Peptide.Imputed",
            filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"), normalization = "none")
        part <- vroom::vroom(part_out, show_col_types = FALSE, .name_repair = "minimal")
        expect_identical(nrow(part), 1L, info = group)
        row <- whole[whole[[protein]] == group, , drop = FALSE]
        shared <- intersect(names(part), names(row))
        expect_true(length(shared) > 1L, info = group)
        expect_equal(unlist(part[1L, setdiff(shared, protein)]), unlist(row[1L, setdiff(shared, protein)]), tolerance = 1e-9, info = group)
    }
})

test_that("the worker entry point computes the IQ rollup from a frozen request (peptide input, protein compute role) and the controller's completion finishes it through the protein finishing step", {
    skip_if_not_installed("iq")
    object <- imputedPeptideFixture()
    store_root <- withr::local_tempdir()
    manager <- proteinDeclarationFixtureManager(object, store_root, "imputed")
    withr::defer(closeDIAStateAdapter(manager))
    store <- manager$.context$store
    reference <- getDIAStateReference(manager)
    directory <- withr::local_tempdir()
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("d", 64), .package = "MultiScholaR")
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = "protein_iq_rollup",
        input_references = list(state = reference), parameters = list(),
        token = list(request_id = "wfx-test"), owner = store$identity, project_root = dirname(store$directory),
        directory = directory, config = list(mode = "adapters"), worker_cap = 12 * 1024^3)
    path <- file.path(directory, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    closeDIAStateAdapter(manager)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"),
        expect_null(suppressMessages(workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))))))
    result <- readRDS(file.path(directory, "worker-result.rds"))
    expect_identical(result$status, "complete")
    expect_identical(result$result$operation, "protein_iq_rollup")
    expect_identical(result$result$checkpoint$class, "ProteinQuantitativeData")
    expect_true(file.exists(result$result$config$iq_output_file))
    expect_true(startsWith(result$result$config$iq_output_file, directory))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(store$directory), store$identity$project_id, store$identity$workflow_id), reopen = TRUE)
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
    owner <- list(manager = manager, accepting = TRUE)
    completed <- completeWorkflowAction(owner, "protein_iq_rollup", result$result, context = list(peptide_qc = list()))
    expect_identical(manager$current_state, "protein_s4_created")
    expect_true(isSavedWorkflowReference(completed$stateReference))
    expect_true(!is.null(manager$.context$snapshot$output_products[[completed$stateReference$entity_id]]))
    expect_true(grepl("IQ Protein Rollup", completed$resultText, fixed = TRUE))
    # The protein result reader hands the product to a module with the QC entries merged.
    bus <- shiny::reactiveValues(qc_params = list(protein_qc = list(earlier = list(a = 1))))
    value <- readManagedProteinResult(manager, list(state_reference = completed$stateReference), bus)
    expect_identical(value$stateReference$entity_id, completed$stateReference$entity_id)
    expect_true(grepl("IQ Protein Rollup", value$resultText, fixed = TRUE))
    expect_identical(names(shiny::isolate(bus$qc_params)$protein_qc), "earlier")
})

test_that("the controller's publication phase publishes every protein state's files and figures from the recorded product (KIT-015 audit PE-02)", {
    skip_if_not_installed("iq")
    imputed <- imputedPeptideFixture()
    protein <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    rlang::local_bindings(config_list = list(removeRowsWithMissingValuesPercent = list()), .env = .GlobalEnv)
    withr::local_dir(withr::local_tempdir())
    publish <- function(operation, manager, value) {
        workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list())))
        completed <- completeDIAProteinAction(workflow, value, publish = FALSE)
        root <- withr::local_tempdir()
        request_dir <- file.path(root, "request")
        dir.create(request_dir)
        owner <- list(manager = manager, control = list(generation = 1L))
        # KIT-023: the publication phase adds no global (the test's own config binding is already in place).
        globals_before <- sort(ls(.GlobalEnv, all.names = TRUE))
        seen <- new.env()
        testthat::local_mocked_bindings(workflowExportBoundary = function(phase) {
            if (identical(phase, "before_qvalue_figures")) seen$globals <- sort(ls(.GlobalEnv, all.names = TRUE))
            invisible(NULL)
        }, .package = "MultiScholaR")
        invisible(capture.output(suppressMessages(publishWorkflowProteinResult(owner, list(request_dir = request_dir, root = root),
            list(operation = operation, state_reference = completed$stateReference))), type = "output"))
        expect_identical(seen$globals, globals_before)
        list.files(file.path(request_dir, "required-output"), recursive = TRUE)
    }
    request <- withr::local_tempdir()
    m1 <- proteinDeclarationFixtureManager(imputed, withr::local_tempdir(), "imputed")
    withr::defer(closeDIAStateAdapter(m1))
    v1 <- computeDIAReferenceIQRollup(m1, getDIAStateReference(m1),
        list(peptide_qc_dir = file.path(request, "iq", "peptide"), protein_qc_dir = file.path(request, "iq", "protein")))
    f1 <- publish("protein_iq_rollup", m1, v1)
    expect_true(any(grepl("protein_qc/iq_output_file[.]txt$", f1)))
    expect_true(any(grepl("peptide_qc/peptide_values_imputed[.]tsv$", f1)))
    expect_true(any(grepl("9_protein_s4_created_combined_plots[.]png$", f1)))
    # IQ's default QC plot is not written into the working directory (audit PE-01).
    expect_false(file.exists("qc-plots.pdf"))
    m2 <- proteinDeclarationFixtureManager(protein, withr::local_tempdir(), "protein_s4_created")
    withr::defer(closeDIAStateAdapter(m2))
    v2 <- computeDIAReferenceAccessionCleanup(m2, getDIAStateReference(m2), list(delimiter = ";", aggregation_method = "mean"))
    f2 <- publish("protein_accession_cleanup", m2, v2)
    expect_true(any(grepl("accession_cleanup_results[.]RDS$", f2)))
    expect_true(any(grepl("10_protein_accession_cleaned_combined_plots[.]png$", f2)))
    m3 <- proteinDeclarationFixtureManager(protein, withr::local_tempdir(), "protein_accession_cleaned")
    withr::defer(closeDIAStateAdapter(m3))
    v3 <- computeDIAReferenceProteinIntensityFilter(m3, getDIAStateReference(m3),
        list(strict_mode = FALSE, min_reps_per_group = 1L, min_groups = 1L, proteins_intensity_cutoff_percentile = 5))
    expect_true(any(grepl("11_protein_intensity_filtered_combined_plots[.]png$", publish("protein_intensity_filter", m3, v3))))
    m4 <- proteinDeclarationFixtureManager(protein, withr::local_tempdir(), "protein_intensity_filtered")
    withr::defer(closeDIAStateAdapter(m4))
    v4 <- computeDIAReferenceProteinDuplicateRemoval(m4, getDIAStateReference(m4), list(aggregation_method = "mean"))
    expect_true(any(grepl("12_duplicates_removed_combined_plots[.]png$", publish("protein_duplicate_removal", m4, v4))))
    # The replicate value carries its field and text too (audit PE-04: one completion shape).
    m5 <- proteinDeclarationFixtureManager(protein, withr::local_tempdir(), "duplicates_removed")
    withr::defer(closeDIAStateAdapter(m5))
    v5 <- computeDIAReferenceProteinReplicateFilter(m5, getDIAStateReference(m5), list(grouping_variable = "group", parallel_cores = 1L))
    expect_identical(v5$field, "filteredS4")
    expect_true(grepl("Protein Replicate Filter Applied Successfully", v5$text, fixed = TRUE))
    expect_error(completeDIAProteinAction(list2env(list(state_manager = m5)), list(state_name = "x")), class = "workflow_execution_invalid")
    expect_true(any(grepl("remove_proteins_with_only_one_rep[.]tsv$", publish("protein_replicate_filter", m5, v5))))
})
