# KIT-006: the complete-sample partitioner and the two declarations on it.

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
