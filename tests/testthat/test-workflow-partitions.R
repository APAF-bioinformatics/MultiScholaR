# KIT-006: the complete-sample partitioner and the two declarations on it; KIT-009: the feature-group partitioner.

test_that("complete-group batches never split a group and keep the declared order", {
    groups <- sprintf("run_%03d", 1:32)
    visited <- list()
    consumed <- integer()
    value <- partitionBySamples(groups, 5L, stage = function(selection, index) {
        visited[[index]] <<- selection
        length(selection)
    }, between = function() consumed <<- c(consumed, -1L), consume = function(staged, index) {
        consumed <<- c(consumed, staged)
        staged * 2L
    })
    expect_identical(unlist(visited), groups)
    expect_identical(value$batch_sizes, as.integer(lengths(visited)))
    expect_identical(sum(value$batch_sizes), 32L)
    expect_length(value$batch_sizes, 5L)
    expect_identical(consumed[[1L]], -1L)
    expect_identical(consumed[-1L], value$batch_sizes)
    expect_identical(unlist(value$consumed), value$batch_sizes * 2L)
    expect_identical(partitionBySamples(character(), 3L, stage = function(selection, index) selection)$batch_sizes, 0L)
    expect_error(partitionBySamples(c("a", "a"), 2L, stage = function(selection, index) NULL), class = "saved_entity_limit")
    one <- partitionBySamples(groups, 1L, stage = function(selection, index) selection)
    expect_identical(one$staged[[1L]], groups)
})

test_that("the admission chooses the fewest batches that fit and refuses above the cap", {
    state <- c(100, 90, 80, 70, 60, 50, 40, 30)
    expect_identical(completeGroupBatchBound(state, 2L), 100 + 90 + 80 + 70)
    expect_identical(completeGroupBatchBound(state, 3L), 100 + 90 + 80)
    expect_identical(completeGroupBatchBound(state, 8L), 100)
    expect_identical(completeGroupBatchBound(numeric(), 3L), 0)
    estimate_for <- function(b) 1000 + 2 * completeGroupBatchBound(state, b)
    admitted <- chooseCompleteGroupBatches(3:8, estimate_for, 1000 + 2 * (100 + 90) + 1, "test-policy")
    expect_identical(admitted$batches, 4L)
    expect_identical(admitted$estimated_bytes, 1000 + 2 * (100 + 90))
    expect_length(admitted$estimates, 8L)
    expect_null(admitted$estimates[[2L]])
    expect_identical(admitted$estimates[[3L]], 1000 + 2 * 270)
    expect_identical(admitted$policy, "test-policy")
    expect_identical(chooseCompleteGroupBatches(3:8, estimate_for, 1e9, "p")$batches, 3L)
    expect_identical(chooseCompleteGroupBatches(6L, estimate_for, 1e9, "p")$batches, 6L)
    expect_error(chooseCompleteGroupBatches(3:8, estimate_for, 1000 + 2 * 100 - 1, "p"), class = "workflow_execution_unsupported")
    expect_error(chooseCompleteGroupBatches(3:8, estimate_for, NA_real_, "p"), class = "workflow_execution_unsupported")
})

test_that("the partitioner names no scientific table, column or operation", {
    body <- readLines(file.path(if (exists("workflow_test_root")) workflow_test_root else here::here(), "R/utils_workflow_partitions.R"))
    expect_false(any(grepl("P14|peptide|protein|precursor|proteom|rollup|DIA|dia_|\"Run\"|filterMin", body)))
})

test_that("the sample filter takes adaptive complete-sample batches: three at the normal cap, more under a reduced cap, exact either way", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    native <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    samples <- length(unique(object@peptide_data$Run))
    expect_gte(samples, 4L)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "kit006-filter", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "fixture")
    reference <- getDIAStateReference(manager)
    parameters <- list(peptides_per_sample_cutoff = 3L, inclusion_list = NULL)
    normal <- computeDIAReferenceSampleFilter(manager, reference, parameters, output_mode = "physical")
    expect_identical(normal$execution$admission$batches, 3L)
    expect_length(normal$execution$admission$count_batch_sizes, 3L)
    expect_identical(sum(normal$execution$admission$count_batch_sizes), samples)
    expect_identical(normal$execution$admission$policy, "sample-adaptive-complete-batches-and-expanded-transfer-v4")
    expect_identical(normal$execution$bulk_query$declared_unit, "complete_sample")
    estimates <- unlist(normal$execution$admission$estimates)
    expect_length(estimates, samples - 2L)
    expect_true(all(diff(estimates) <= 0))
    # A cap between the four- and three-batch estimates admits exactly four batches.
    cap <- (normal$execution$admission$estimates[[4L]] + normal$execution$admission$estimates[[3L]]) / 2
    reduced <- computeDIAReferenceSampleFilter(manager, reference, parameters, worker_cap = cap, output_mode = "physical")
    expect_identical(reduced$execution$admission$batches, 4L)
    expect_length(reduced$execution$admission$count_batch_sizes, 4L)
    for (value in list(normal, reduced)) {
        expect_identical(materializeDIATable(manager$.context, value$reference), native@peptide_data)
        expect_identical(value$summary, native@args$filterMinNumPeptidesPerSample$filter_summary)
        expect_identical(materializeDIATable(manager$.context, value$checkpoint$support$filterMinNumPeptidesPerSample$support_table),
            native@args$filterMinNumPeptidesPerSample$support_table)
    }
    expect_error(computeDIAReferenceSampleFilter(manager, reference, parameters,
        worker_cap = normal$execution$admission$estimates[[samples]] - 1, output_mode = "physical"),
        class = "workflow_execution_unsupported")
    expect_identical(getDIAStateReference(manager), reference)
    expect_identical(store$open_resources, 0L)
})

test_that("the P14 rollup and the sample filter are declarations on the complete-sample unit whose kernels run the reference entry points", {
    rollup <- diaPrecursorRollupOperation()
    filter <- diaSampleFilterOperation()
    expect_identical(rollup$unit, "complete_sample")
    expect_identical(filter$unit, "complete_sample")
    expect_identical(rollup$name, "precursor_rollup")
    expect_identical(filter$name, "minimum_peptides_per_sample")
    expect_true(all(c("global_validation", "protein_annotation_lookup") %in% names(rollup$global_statistics)))
    expect_identical(unname(filter$merge[["counts"]]), "sum")
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "kit006-declarations", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "fixture")
    reference <- getDIAStateReference(manager)
    filtered <- runProjected(manager, filter, reference, list(peptides_per_sample_cutoff = 3L, inclusion_list = NULL))
    native <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    expect_identical(materializeDIATable(manager$.context, filtered$reference), native@peptide_data)
    expect_identical(filtered$state_name, "sample_filtered")
    expect_identical(filtered$execution$admission$batches, 3L)
    rolled <- runProjected(manager, rollup, reference, list(worker_cap = 12 * 1024^3))
    expect_identical(rolled$state_name, "precursor_rollup")
    expect_identical(rolled$execution$independence_unit, "complete_sample")
    expect_identical(rolled$execution$profile$name, "P14")
    expect_gte(rolled$execution$admission$batches, 3L)
    expect_identical(sum(rolled$execution$batch_sizes), length(unique(object@peptide_data$Run)))
    expected <- rollUpPrecursorToPeptide(object)
    expect_identical(rolled$summary, expected@args$rollUpPrecursorToPeptide$rollup_summary)
    expect_identical(getDIAStateReference(manager), reference)
    expect_identical(store$open_resources, 0L)
})

# KIT-009: the feature-group partitioner and the peptide-group family on it.

peptideGroupFixtureManager <- function(object, root) {
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState("input", object, list(), "fixture")
    manager
}

test_that("feature-group batches never split a group, run the kernel once per batch and record the admission", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager <- peptideGroupFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    checkpoint <- readDIAStateMetadata(manager)$data
    scope <- newWorkflowTableReadScope(manager$.context$store, checkpoint$components$peptide_data)
    withr::defer(releaseWorkflowTableReadScope(scope))
    columns <- c(object@sample_id, object@protein_id_column, object@peptide_sequence_column)
    seen <- list()
    consumed <- integer()
    value <- partitionByFeatureGroup(scope, columns, object@protein_id_column,
        kernel = function(part, index) {
            seen[[index]] <<- unique(as.character(part[[object@protein_id_column]]))
            expect_identical(length(attr(part, "workflow_source_order")), nrow(part))
            nrow(part)
        }, consume = function(rows, index) consumed <<- c(consumed, rows),
        worker_cap = 12 * 1024^3, working_bytes_per_row = 2048, overhead_bytes = 1024^3, policy = "test-policy")
    groups <- unlist(seen)
    expect_false(anyDuplicated(groups) > 0L)
    expect_setequal(groups, unique(as.character(object@peptide_data[[object@protein_id_column]])))
    expect_identical(value$groups, length(groups))
    expect_identical(value$batches, length(seen))
    expect_identical(value$batch_sizes, consumed)
    expect_identical(sum(value$batch_sizes), nrow(object@peptide_data))
    expect_identical(value$admission$policy, "test-policy")
    expect_identical(value$admission$working_bytes_per_row, 2048)
    expect_gte(value$admission$batch_rows, max(table(object@peptide_data[[object@protein_id_column]])))
    expect_identical(value$admission$estimated_bytes,
        1024^3 + value$admission$batch_rows * (value$admission$row_bytes + 2048))
    expect_error(partitionByFeatureGroup(scope, columns, object@protein_id_column,
        kernel = function(part, index) NULL, consume = function(value, index) NULL,
        worker_cap = value$admission$estimated_bytes - 1, working_bytes_per_row = 2048, overhead_bytes = 1024^3,
        policy = "p"), class = "workflow_execution_unsupported")
    expect_identical(manager$.context$store$open_resources, 1L)
})

test_that("every peptide-group section is per-group independent on the fixture: the partitioned and grouped variants equal the kernel run per group", {
    # The fixture is precursor-level: the protein-evidence filter consumes it as production does (before the
    # rollup, with Modified.Sequence present); the other sections consume the rolled-up, sample-filtered table.
    precursor <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- filterMinNumPeptidesPerSample(rollUpPrecursorToPeptide(precursor, core_utilisation = NA),
        peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    rlang::local_bindings(config_list = list(
        filterMinNumPeptidesPerProtein = list(), peptideIntensityFiltering = list(),
        removePeptidesWithOnlyOneReplicate = list(replicate_group_column = "replicates"),
        peptideMissingValueImputation = list(proportion_missing_values = 0.5, exclusion_column = NULL)), .env = .GlobalEnv)
    sections <- list(
        list(section = "filterMinNumPeptidesPerProtein", parameters = list(num_peptides_per_protein_thresh = 2L,
            num_peptidoforms_per_protein_thresh = 1L), description = "protein evidence", object = precursor),
        list(section = "peptideIntensityFiltering", parameters = list(strict_mode = FALSE, min_reps_per_group = 1L, min_groups = 1L,
            peptides_intensity_cutoff_percentile = 5, grouping_variable = "group"), description = "intensity", object = rolled),
        list(section = "removePeptidesWithOnlyOneReplicate", parameters = list(replicate_group_column = "replicates"),
            description = "replicates", object = rolled),
        list(section = "peptideMissingValueImputation", parameters = list(proportion_missing_values = 0.5, exclusion_column = NULL),
            description = "imputation", object = rolled))
    for (entry in sections) {
        object <- entry$object
        root <- withr::local_tempdir()
        manager <- peptideGroupFixtureManager(object, root)
        reference <- getDIAStateReference(manager)
        partitioned <- computeDIAReferencePeptideGroups(manager, reference, entry$section, entry$parameters,
            variant = "partitioned", description = entry$description)
        grouped <- computeDIAReferencePeptideGroups(manager, reference, entry$section, entry$parameters,
            variant = "grouped", description = entry$description)
        expect_identical(partitioned$execution$variant, "partitioned")
        expect_identical(partitioned$execution$admission$policy, "feature-group-cursor-batches-v1")
        expect_identical(partitioned$execution$batch_sizes, grouped$execution$batch_sizes)
        expect_identical(partitioned$execution$groups, length(unique(object@peptide_data[[object@protein_id_column]])))
        expect_lte(partitioned$execution$engine_limit_bytes, 1024^3)
        data_partitioned <- materializeDIATable(manager$.context, partitioned$reference)
        data_grouped <- materializeDIATable(manager$.context, grouped$reference)
        expect_identical(data_partitioned, data_grouped)
        expect_identical(partitioned$summary, grouped$summary)
        expect_identical(names(partitioned$checkpoint$support[[entry$section]]), names(grouped$checkpoint$support[[entry$section]]))
        for (name in names(partitioned$checkpoint$support[[entry$section]]))
            expect_identical(materializeDIATable(manager$.context, partitioned$checkpoint$support[[entry$section]][[name]]),
                materializeDIATable(manager$.context, grouped$checkpoint$support[[entry$section]][[name]]))
        # Per-group independence: the kernel on each protein group alone, concatenated in source order, is the
        # same table (the imputation generates rows, so the comparison is on the retained and generated sets).
        groups <- split(seq_len(nrow(object@peptide_data)), as.character(object@peptide_data[[object@protein_id_column]]))
        per_group <- lapply(groups, function(rows) {
            part <- object
            part@peptide_data <- object@peptide_data[rows, , drop = FALSE]
            if (entry$section == "peptideIntensityFiltering") {
                quantity <- object@norm_quantity_column
                values <- object@peptide_data[[quantity]]
                threshold <- unname(stats::quantile(values[is.finite(values)], probs = 5 / 100, na.rm = TRUE)[[1L]])
                peptideIntensityFilteringHelper(part@peptide_data, object@design_matrix,
                    min_peptide_intensity_threshold = threshold, sample_id_column = object@sample_id, grouping_variable = "group",
                    protein_id_column = object@protein_id_column, peptide_sequence_column = object@peptide_sequence_column,
                    peptide_quantity_column = quantity, core_utilisation = NA, min_reps_per_group = 1L, min_groups = 1L,
                    strict_mode = FALSE, groupwise_percentage_cutoff = 1, max_groups_percentage_cutoff = 50,
                    return_filter_result = TRUE, require_complete_design_runs = FALSE)$data
            } else do.call(get(entry$section, mode = "function"), c(list(theObject = part), entry$parameters))@peptide_data
        })
        expected <- dplyr::bind_rows(per_group)
        key <- function(data) {
            columns <- intersect(c(object@sample_id, object@protein_id_column, object@peptide_sequence_column,
                "Modified.Sequence", object@norm_quantity_column, "Peptide.Imputed", "is_imputed"), names(data))
            out <- as.data.frame(data[columns])
            rownames(out) <- NULL
            out[do.call(order, unname(as.list(out[intersect(c(object@protein_id_column, object@peptide_sequence_column, object@sample_id), columns)]))), , drop = FALSE]
        }
        actual <- key(data_partitioned)
        wanted <- key(expected)
        rownames(actual) <- rownames(wanted) <- NULL
        expect_equal(actual, wanted, ignore_attr = TRUE, info = entry$section)
        expect_identical(getDIAStateReference(manager), reference)
        closeDIAStateAdapter(manager)
    }
})

test_that("the imputation is a declaration on the feature-group unit whose kernel runs the compute-only entry point", {
    declaration <- diaPeptideImputationOperation()
    expect_identical(declaration$unit, "feature_group")
    expect_identical(declaration$name, "peptide_imputation")
    expect_identical(declaration$merge[["peptide_data"]], "concatenate")
    expect_identical(declaration$merge[["design_matrix"]], "design")
    object <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    object <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    rlang::local_bindings(config_list = list(peptideMissingValueImputation = list(proportion_missing_values = 0.5)), .env = .GlobalEnv)
    manager <- peptideGroupFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager)
    native <- peptideMissingValueImputation(object, proportion_missing_values = 0.5, core_utilisation = NA)
    value <- runProjected(manager, declaration, reference, list(proportion_missing_values = 0.5, exclusion_column = NULL))
    expect_identical(value$state_name, "imputed")
    expect_identical(value$section, "peptideMissingValueImputation")
    expect_identical(value$execution$variant, "partitioned")
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
        value$section, value$counts)
    expect_identical(manager$current_state, "imputed")
    expect_identical(nrow(manager$getState("imputed")@peptide_data), nrow(native@peptide_data))
    expect_identical(manager$getState("imputed")@args$peptideMissingValueImputation$imputation_summary$imputed_rows,
        native@args$peptideMissingValueImputation$imputation_summary$imputed_rows)
})

test_that("a peptide-group section computed outside a session and without a global configuration keeps its configuration private (worker process)", {
    object <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    object <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    manager <- peptideGroupFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    if (exists("config_list", envir = globalenv(), inherits = FALSE)) {
        saved <- get("config_list", envir = globalenv())
        rm("config_list", envir = globalenv())
        withr::defer(assign("config_list", saved, envir = globalenv()))
    }
    value <- computeDIAReferencePeptideGroups(manager, getDIAStateReference(manager), "peptideMissingValueImputation",
        list(proportion_missing_values = 0.5, exclusion_column = NULL))
    expect_identical(value$state_name, "imputed")
    expect_identical(value$checkpoint$slots$args$peptideMissingValueImputation$proportion_missing_values, 0.5)
    expect_false(exists("config_list", envir = globalenv(), inherits = FALSE))
})
