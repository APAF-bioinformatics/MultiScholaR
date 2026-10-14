# KIT-006: the complete-sample partitioner and the two declarations on it; KIT-009: the feature-group partitioner;
# KIT-010: the row partitioner and wide-table sample projection.

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

# KIT-010: the row partitioner, wide-table projection and the protein replicate filter on the row unit.

proteinFixtureManager <- function(object, root, state = "duplicates_removed") {
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    manager$.context$reference_states <- TRUE
    manager$saveState(state, object, list(), "fixture")
    manager
}

wideProteinRecipe <- function() {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_samples <- 300L
    recipe$n_proteins <- 300L
    recipe$peptides_per_protein <- 1L
    recipe$forms_per_peptide <- 1L
    recipe
}

test_that("row pages cover a long table in source order, the kernel once per page, and refuse above the cap", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    # 4,608 rows: more than one 10,000-row page needs a larger table; repeat the rows to 23,040.
    big <- object
    big@peptide_data <- do.call(rbind, replicate(5L, object@peptide_data, simplify = FALSE))
    rownames(big@peptide_data) <- NULL
    manager <- peptideGroupFixtureManager(big, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    checkpoint <- readDIAStateMetadata(manager)$data
    scope <- newWorkflowTableReadScope(manager$.context$store, checkpoint$components$peptide_data, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    columns <- c(object@sample_id, object@protein_id_column, object@peptide_sequence_column)
    offsets <- integer(); first <- character()
    value <- partitionByRows(scope, columns, kernel = function(page, index, offset) {
        offsets[[index]] <<- offset
        first[[index]] <<- as.character(page[[object@protein_id_column]][[1L]])
        nrow(page)
    }, consume = function(rows, index) NULL, worker_cap = 12 * 1024^3, working_bytes_per_row = 256,
        overhead_bytes = 1e6, policy = "test-rows")
    expect_identical(value$batches, 3L)
    expect_identical(value$batch_sizes, c(10000L, 10000L, 3040L))
    expect_identical(value$rows, 23040L)
    expect_identical(offsets, c(0L, 10000L, 20000L))
    expect_identical(first, as.character(big@peptide_data[[object@protein_id_column]][c(1L, 10001L, 20001L)]))
    expect_identical(value$admission$page_rows, 10000L)
    expect_identical(value$admission$estimated_bytes, 1e6 + 10000 * (value$admission$row_bytes + 256))
    expect_error(partitionByRows(scope, columns, kernel = function(page, index, offset) NULL,
        consume = function(value, index) NULL, worker_cap = value$admission$estimated_bytes - 1,
        working_bytes_per_row = 256, overhead_bytes = 1e6, policy = "p"), class = "workflow_execution_unsupported")
})

test_that("a wide 300-sample protein table streams in row pages with every column and projects to a subset of samples (audit B24)", {
    object <- makeWorkflowCheckpointProtein(wideProteinRecipe())
    expect_identical(ncol(object@protein_quant_table), 301L)
    manager <- proteinFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager)
    checkpoint <- readDIAStateMetadata(manager)$data
    component <- checkpoint$components$protein_quant_table
    expect_identical(component$metadata$shape[[2L]], 301L)
    scope <- newWorkflowTableReadScope(manager$.context$store, component, reader = "duckdb")
    withr::defer(releaseWorkflowTableReadScope(scope))
    seen <- NULL
    value <- partitionByRows(scope, scope$metadata$columns, kernel = function(page, index, offset) {
        expect_identical(names(page), names(object@protein_quant_table))
        page
    }, consume = function(page, index) seen <<- rbind(seen, page), worker_cap = 12 * 1024^3,
        working_bytes_per_row = 48 * 301, overhead_bytes = 1e6, policy = "test-wide")
    expect_identical(value$batches, 1L)
    expect_identical(value$admission$columns, 301L)
    rownames(seen) <- NULL
    expect_equal(seen, object@protein_quant_table, ignore_attr = TRUE)
    # Sample projection through the declaration: the protein key and 150 of the 300 sample columns.
    samples <- object@design_matrix$Run[seq_len(150L)]
    projected <- readDIAExecutionComponents(manager, reference, fields = "protein_quant_table",
        columns = list(protein_quant_table = c("Protein.Group", samples)))
    expect_identical(names(projected@protein_quant_table), c("Protein.Group", samples))
    expect_identical(nrow(projected@protein_quant_table), 300L)
    expect_equal(projected@protein_quant_table[[samples[[150L]]]], object@protein_quant_table[[samples[[150L]]]])
    # The design follows the projection, so the object is valid with the selected samples only.
    expect_identical(as.character(projected@design_matrix$Run), as.character(samples))
})

test_that("the protein replicate filter on row pages equals the whole-object kernel and the native method (KIT-010)", {
    object <- makeWorkflowCheckpointProtein(wideProteinRecipe())
    # Make the decision non-trivial: blank out one group's values for a tenth of the proteins.
    design <- object@design_matrix
    runs_in_group <- design$Run[design$group == unique(design$group)[[1L]]]
    drop <- seq(1L, nrow(object@protein_quant_table), by = 10L)
    for (run in runs_in_group) object@protein_quant_table[drop, run] <- NA_real_
    native <- removeProteinsWithOnlyOneReplicate(object, NA, grouping_variable = "group")
    expect_lt(nrow(native@protein_quant_table), nrow(object@protein_quant_table))
    for (variant in c("whole", "rows")) {
        manager <- proteinFixtureManager(object, withr::local_tempdir())
        reference <- getDIAStateReference(manager)
        value <- computeDIAReferenceProteinReplicateFilter(manager, reference, list(grouping_variable = "group", parallel_cores = 4L),
            variant = variant)
        expect_identical(value$state_name, "protein_replicate_filtered")
        expect_identical(value$execution$variant, variant)
        expect_identical(value$config$parallel_cores, 4L)
        table <- materializeDIATable(manager$.context, value$reference)
        expect_equal(table, native@protein_quant_table, ignore_attr = TRUE, info = variant)
        design_out <- materializeDIATable(manager$.context, value$checkpoint$components$design_matrix)
        expect_equal(design_out, native@design_matrix, ignore_attr = TRUE, info = variant)
        expect_identical(value$checkpoint$slots$args$removeProteinsWithOnlyOneReplicate$grouping_variable, "group")
        expect_identical(value$protein_count, nrow(native@protein_quant_table))
        expect_identical(value$checkpoint$dimensions$rows, nrow(native@protein_quant_table))
        if (variant == "rows") {
            expect_identical(value$execution$admission$policy, "row-pages-v1")
            expect_identical(value$execution$batches, 1L)
        }
        expect_identical(getDIAStateReference(manager), reference)
        closeDIAStateAdapter(manager)
    }
})

test_that("the protein replicate filter is a declaration on the row unit and completes through the protein completion", {
    declaration <- diaProteinReplicateFilterOperation()
    expect_identical(declaration$unit, "row")
    expect_identical(declaration$name, "protein_replicate_filter")
    expect_identical(unname(declaration$merge[["protein_quant_table"]]), "concatenate")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager <- proteinFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager)
    value <- runProjected(manager, declaration, reference, list(grouping_variable = "group", parallel_cores = 1L))
    expect_identical(value$execution$variant, "rows")
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    result <- completeDIAProteinAction(workflow, value)
    expect_identical(manager$current_state, "protein_replicate_filtered")
    expect_true(isSavedWorkflowReference(result$stateReference))
    expect_identical(result$proteinCount, value$protein_count)
    native <- removeProteinsWithOnlyOneReplicate(object, NA, grouping_variable = "group")
    expect_equal(manager$getState("protein_replicate_filtered")@protein_quant_table, native@protein_quant_table, ignore_attr = TRUE)
    expect_true(!is.null(manager$.context$snapshot$output_products[[result$stateReference$entity_id]]))
})

test_that("the protein publication phase publishes the replicate filter's table and figures from the recorded product (KIT-010)", {
    skip_if_not_installed("withr")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    root <- withr::local_tempdir()
    store_root <- withr::local_tempdir()
    manager <- proteinFixtureManager(object, store_root)
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager)
    value <- computeDIAReferenceProteinReplicateFilter(manager, reference, list(grouping_variable = "group", parallel_cores = 1L))
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
    completed <- completeDIAProteinAction(workflow, value)
    request_dir <- file.path(root, "request")
    dir.create(request_dir)
    owner <- list(manager = manager, control = list(generation = 1L))
    paths <- list(request_dir = request_dir, root = root)
    result <- list(operation = "protein_replicate_filter", state_reference = completed$stateReference)
    # KIT-023: the controller's publication phase binds nothing into the global environment (the figure writer
    # takes its directories as an argument); observed at the figure boundary.
    globals_before <- sort(ls(.GlobalEnv, all.names = TRUE))
    seen <- new.env()
    testthat::local_mocked_bindings(workflowExportBoundary = function(phase) {
        if (identical(phase, "before_qvalue_figures")) seen$globals <- sort(ls(.GlobalEnv, all.names = TRUE))
        invisible(NULL)
    }, .package = "MultiScholaR")
    invisible(capture.output(suppressMessages(publishWorkflowProteinResult(owner, paths, result)), type = "output"))
    expect_identical(seen$globals, globals_before)
    expect_false(any(c("project_dirs", "config_list") %in% setdiff(seen$globals, globals_before)))
    published <- list.files(file.path(request_dir, "required-output"), recursive = TRUE)
    expect_true(any(grepl("remove_proteins_with_only_one_rep[.]tsv$", published)))
    expect_true(any(grepl("qc_params[.]RDS$", published)))
    expect_true(any(grepl("13_protein_replicate_filtered_combined_plots[.]png$", published)))
    expect_identical(getDIAStateReference(manager), completed$stateReference)
})

test_that("the worker entry point computes a protein operation from a frozen request (KIT-010 attempt 2: the worker's dispatch by compute role)", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    store_root <- withr::local_tempdir()
    manager <- proteinFixtureManager(object, store_root)
    withr::defer(closeDIAStateAdapter(manager))
    store <- manager$.context$store
    reference <- getDIAStateReference(manager)
    directory <- withr::local_tempdir()
    # The build identity needs the installed runtime; the request and the worker see the same mocked value.
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = "protein_replicate_filter",
        input_references = list(state = reference), parameters = list(grouping_variable = "group", parallel_cores = 4L),
        token = list(request_id = "wfx-test"), owner = store$identity, project_root = dirname(store$directory),
        directory = directory, config = list(mode = "adapters"), worker_cap = 12 * 1024^3)
    path <- file.path(directory, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::local_envvar(MULTISCHOLAR_WORKFLOW_WORKER = "1")
    closeDIAStateAdapter(manager)
    expect_null(suppressMessages(workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256")))))
    result <- readRDS(file.path(directory, "worker-result.rds"))
    expect_identical(result$status, "complete")
    expect_identical(result$result$operation, "protein_replicate_filter")
    expect_identical(result$result$execution$variant, "rows")
    expect_true(isDIAComponentCheckpoint(result$result$checkpoint))
    expect_identical(result$result$checkpoint$class, "ProteinQuantitativeData")
    # Every registered compute role is one the worker's entry point dispatches.
    for (name in workflowOperationNames()) expect_true(workflowActionSpec(name)$compute %in% c("fixture", "peptide", "protein", "in_process"), info = name)
})

test_that("the protein evidence, intensity and replicate filters are declarations on the feature-group unit whose kernels run the family's compute-only entry point (KIT-014)", {
    precursor <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    rolled <- filterMinNumPeptidesPerSample(rollUpPrecursorToPeptide(precursor, core_utilisation = NA),
        peptides_per_sample_cutoff = 3L, inclusion_list = NULL, core_utilisation = NA)
    rlang::local_bindings(config_list = list(filterMinNumPeptidesPerProtein = list(), peptideIntensityFiltering = list(),
        removePeptidesWithOnlyOneReplicate = list(replicate_group_column = "replicates")), .env = .GlobalEnv)
    entries <- list(
        list(declaration = diaPeptideProteinEvidenceOperation(), name = "peptide_protein_evidence_filter",
            section = "filterMinNumPeptidesPerProtein", state = "protein_peptide_filtered", object = precursor,
            parameters = list(num_peptides_per_protein_thresh = 2L, num_peptidoforms_per_protein_thresh = 1L),
            native = function(object) filterMinNumPeptidesPerProtein(object, num_peptides_per_protein_thresh = 2L,
                num_peptidoforms_per_protein_thresh = 1L)),
        list(declaration = diaPeptideIntensityFilterOperation(), name = "peptide_intensity_filter",
            section = "peptideIntensityFiltering", state = "intensity_filtered", object = rolled,
            parameters = list(strict_mode = FALSE, min_reps_per_group = 1L, min_groups = 1L,
                peptides_intensity_cutoff_percentile = 5, grouping_variable = "group"), native = NULL),
        list(declaration = diaPeptideReplicateFilterOperation(), name = "peptide_replicate_filter",
            section = "removePeptidesWithOnlyOneReplicate", state = "replicate_filtered", object = rolled,
            parameters = list(replicate_group_column = "replicates"),
            native = function(object) removePeptidesWithOnlyOneReplicate(object, replicate_group_column = "replicates")))
    for (entry in entries) {
        declaration <- entry$declaration
        expect_identical(declaration$unit, "feature_group", info = entry$name)
        expect_identical(declaration$name, entry$name)
        expect_identical(unname(declaration$merge[["peptide_data"]]), "concatenate")
        expect_identical(unname(declaration$merge[["design_matrix"]]), "design")
        expect_identical(unname(declaration$merge[["counts"]]), "sum")
        expect_identical(workflowActionSpec(entry$name)$default_variant, "partitioned")
        expect_identical(workflowActionSpec(entry$name)$variants, c("grouped", "partitioned"))
        manager <- peptideGroupFixtureManager(entry$object, withr::local_tempdir())
        reference <- getDIAStateReference(manager)
        value <- runProjected(manager, declaration, reference, entry$parameters)
        expect_identical(value$state_name, entry$state)
        expect_identical(value$section, entry$section)
        expect_identical(value$execution$variant, "partitioned")
        workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list())))
        result <- finishDIAReferencePeptideStep(workflow, value$checkpoint, value$state_name, value$description, value$config,
            value$section, value$counts)
        expect_identical(manager$current_state, entry$state)
        expect_true(isSavedWorkflowReference(result$stateReference))
        expect_identical(nrow(manager$getState(entry$state)@peptide_data), value$checkpoint$dimensions$rows)
        if (is.function(entry$native)) {
            native <- entry$native(entry$object)
            expect_identical(nrow(manager$getState(entry$state)@peptide_data), nrow(native@peptide_data), info = entry$name)
        }
        closeDIAStateAdapter(manager)
    }
    # The intensity declaration carries the dataset-wide quantile as its global statistic (audit A1).
    expect_true("min_peptide_intensity_threshold" %in% names(diaPeptideIntensityFilterOperation()$global_statistics))
    # The GUI entry point of the family goes through the registry (an unknown section refuses before any work).
    expect_identical(unname(diaPeptideGroupOperations[names(diaPeptideGroupStates)]),
        c("peptide_protein_evidence_filter", "peptide_intensity_filter", "peptide_replicate_filter", "peptide_imputation"))
})

test_that("the Q-value filter is a declaration whose compute-only entry point equals the module's kernel on the fixture and whose finishing step commits the state with its progress product (KIT-014)", {
    declaration <- diaQvalueFilterOperation()
    expect_identical(declaration$unit, "feature_group")
    expect_identical(declaration$name, "peptide_qvalue_filter")
    expect_setequal(names(declaration$global_statistics), c("protein_identification_evidence", "source_row_offsets"))
    expect_identical(workflowActionSpec("peptide_qvalue_filter")$variants, "streamed")
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    # The design step records the retained columns on the object (the kernel's default set carries no sample id).
    object@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids <- names(object@peptide_data)
    rlang::local_bindings(config_list = list(srlQvalueProteotypicPeptideClean = list()), .env = .GlobalEnv)
    parameters <- list(qvalue_threshold = 0.01, global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
        proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report")
    validateWorkflowActionParameters("peptide_qvalue_filter", parameters, workflowExecutionPolicy(FALSE))
    expect_error(validateWorkflowActionParameters("peptide_qvalue_filter", parameters[-1L], workflowExecutionPolicy(FALSE)))
    # The module's kernel on the in-memory owner (the path the GUI took before KIT-011).
    native_state <- WorkflowState$new(audit_enabled = FALSE)
    native_state$saveState("raw_data_s4", object, list(), "fixture")
    native_workflow <- list2env(list(state_manager = native_state, qc_params = list()))
    native <- runPeptideQvalueApplyStep(native_workflow, 0.01, 0.01, TRUE, logInfoFn = function(...) NULL,
        logWarnFn = function(...) NULL)$filteredS4
    manager <- peptideGroupFixtureManager(object, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager)
    value <- runProjected(manager, declaration, reference, parameters)
    expect_identical(value$state_name, "qvalue_filtered")
    expect_identical(value$section, "srlQvalueProteotypicPeptideClean")
    expect_identical(value$execution$variant, "streamed")
    # KIT-022 (KIT-014 audit PH-02): the admission estimate is recorded, and a cap below it refuses before any work.
    admission <- value$execution$admission
    expect_true(is.numeric(admission$estimated_bytes) && admission$estimated_bytes > 67108864)
    expect_identical(admission$worker_cap, 12 * 1024^3)
    expect_identical(admission$policy, "qvalue-streamed-page-v1")
    expect_true(admission$batch_rows >= 1 && admission$batch_rows <= 10000)
    expect_error(computeDIAReferenceQvalue(manager, reference, parameters, worker_cap = 67108864 + 1),
        class = "workflow_execution_unsupported")
    expect_identical(value$checkpoint$dimensions$rows, nrow(native@peptide_data))
    expect_identical(value$exclusion_summary$excluded_rows,
        native@args$srlQvalueProteotypicPeptideClean$biological_exclusion_summary$excluded_rows)
    expect_identical(value$qc$qvalue_threshold, 0.01)
    expect_true(grepl(paste("Proteins remaining:", value$proteins), value$text, fixed = TRUE))
    expect_null(value$checkpoint$products$qvalue_progress)
    workflow <- list2env(list(state_manager = manager, qc_params = list(peptide_qc = list(sample_filter = list(cutoff = 3L)))))
    result <- finishDIAReferenceQvalueStep(workflow, value)
    expect_identical(manager$current_state, "qvalue_filtered")
    expect_true(isSavedWorkflowReference(result$stateReference))
    record <- readDIAStateMetadata(manager)$data
    expect_true(isSavedWorkflowReference(record$products$qvalue_progress))
    expect_true(isSavedWorkflowReference(record$support$srlQvalueProteotypicPeptideClean$biological_exclusion_ledger))
    expect_identical(names(workflow$qc_params$peptide_qc), c("sample_filter", "qvalue_filter"))
    expect_identical(result$qcParams$peptide_qc$qvalue_filter$proteotypic_only, TRUE)
    committed <- manager$getState("qvalue_filtered")
    expect_identical(nrow(committed@peptide_data), nrow(native@peptide_data))
    expect_identical(sort(unique(committed@peptide_data$Protein.Ids)), sort(unique(native@peptide_data$Protein.Ids)))
    # The progress product is the one the managed result reader and the publication read for this step.
    expect_identical(workflowPeptideProgressProduct(record)$entity_id, record$products$qvalue_progress$entity_id)
})

test_that("the controller's completion finishes a Q-value result through the step's own finishing step and the worker entry point computes it from a frozen request (KIT-014 launch chain)", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids <- names(object@peptide_data)
    store_root <- withr::local_tempdir()
    manager <- peptideGroupFixtureManager(object, store_root)
    withr::defer(closeDIAStateAdapter(manager))
    store <- manager$.context$store
    reference <- getDIAStateReference(manager)
    parameters <- list(qvalue_threshold = 0.01, global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
        proteotypic_only = TRUE, confidence_provenance_mode = "diann_main_report")
    directory <- withr::local_tempdir()
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("e", 64), .package = "MultiScholaR")
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = "peptide_qvalue_filter",
        input_references = list(state = reference), parameters = parameters,
        token = list(request_id = "wfx-test"), owner = store$identity, project_root = dirname(store$directory),
        directory = directory, config = list(mode = "adapters"), worker_cap = 12 * 1024^3)
    path <- file.path(directory, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    closeDIAStateAdapter(manager)
    # No global configuration in the worker: the template keeps a private copy seeded from the saved state.
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"),
        expect_null(suppressMessages(workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))))))
    result <- readRDS(file.path(directory, "worker-result.rds"))
    expect_identical(result$status, "complete")
    expect_identical(result$result$operation, "peptide_qvalue_filter")
    expect_identical(result$result$section, "srlQvalueProteotypicPeptideClean")
    expect_identical(result$result$execution$variant, "streamed")
    expect_true(isDIAComponentCheckpoint(result$result$checkpoint))
    # The completion branch: a Q-value value finishes through its own step under the controller's acceptance lease.
    manager <- newDIAStateAdapter(newWorkflowEntityStore(dirname(store$directory), store$identity$project_id, store$identity$workflow_id), reopen = TRUE)
    withr::defer(closeDIAStateAdapter(manager))
    testthat::local_mocked_bindings(assertWorkflowController = function(...) invisible(TRUE), .package = "MultiScholaR")
    owner <- list(manager = manager, accepting = TRUE)
    completed <- completeWorkflowAction(owner, "peptide_qvalue_filter", result$result, context = list(peptide_qc = list()))
    expect_identical(manager$current_state, "qvalue_filtered")
    expect_true(isSavedWorkflowReference(completed$stateReference))
    expect_true(isSavedWorkflowReference(readDIAStateMetadata(manager)$data$products$qvalue_progress))
    # Every KIT-014 operation names handlers and validators that exist and the compute role the worker dispatches.
    for (name in c("peptide_qvalue_filter", "peptide_protein_evidence_filter", "peptide_intensity_filter", "peptide_replicate_filter")) {
        entry <- workflowActionSpec(name)
        expect_true(is.function(get(entry$handler, mode = "function")), info = name)
        expect_true(is.function(get(entry$validate, mode = "function")), info = name)
        expect_identical(entry$compute, "peptide")
        expect_identical(entry$completion, "peptide")
    }
})
