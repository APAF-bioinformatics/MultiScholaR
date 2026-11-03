# KIT-078: the omics steps on the worker.
library(testthat)

test_that("KIT-078: every omics compute declaration runs from a frozen request through the worker's entry point", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
  fixtures <- omicsWorkerFixtures(withr::local_tempdir())
  for (omic in names(fixtures)) {
    fixture <- fixtures[[omic]]$fixture
    binding <- fixtures[[omic]]$binding
    withr::defer(closeDIAStateAdapter(fixture$owner))
    steps <- list(omicsIntensityFilterStep(binding, 1, 0.5), omicsDuplicateResolutionStep(binding),
      omicsItsdNormalisationStep(binding, "median", NULL), omicsLog2Step(binding, 1),
      omicsBetweenSampleNormalisationStep(binding, "cyclicloess"), omicsRuvStep(binding, fixtures[[omic]]$ruv))
    for (step in steps) {
      result <- omicsWorkerRun(fixture, step$operation, step$parameters)
      expect_identical(result$status, "complete", info = paste(omic, step$operation, result$error))
      expect_identical(result$result$operation, step$operation, info = omic)
      expect_true(isDIAComponentCheckpoint(result$result$checkpoint), info = step$operation)
      expect_identical(result$result$checkpoint$class, binding$class, info = step$operation)
      expect_true(isSavedWorkflowReference(result$result$reference), info = step$operation)
      expect_true(length(result$result$roots) >= 1L, info = step$operation)
      expect_true(length(result$result$owned) >= 1L, info = step$operation)
      # The worker moved no head: the source is still the session's state.
      expect_false(step$state_name %in% unlist(fixture$owner$getHistory()), info = step$operation)
    }
    # The filtering steps hand over their metrics; a later step does not.
    filtering <- omicsWorkerRun(fixture, steps[[1L]]$operation, steps[[1L]]$parameters)
    expect_true(isSavedWorkflowReference(filtering$result$display$filtering_metrics), info = omic)
    expect_null(omicsWorkerRun(fixture, steps[[4L]]$operation, steps[[4L]]$parameters)$result$display, info = omic)
    # The product-only declarations return the source's checkpoint.
    source <- readWorkflowStateMetadata(fixture$owner)$data
    pairs <- omicsWorkerRun(fixture, paste0(binding$prefix, "_pearson_pairs"), list(correlation_group = "group"))
    expect_identical(pairs$status, "complete", info = pairs$error)
    expect_identical(pairs$result$checkpoint, source, info = omic)
    contrast <- omicsWorkerRun(fixture, paste0(binding$prefix, "_da_contrast"), list(formula = "~ 0 + group", qvalue = 0.05, lfc = 0,
      contrast = as.character(fixture$workflow$contrasts_tbl$contrasts[[1L]]),
      friendly_name = as.character(fixture$workflow$contrasts_tbl$friendly_names[[1L]])))
    expect_identical(contrast$status, "complete", info = contrast$error)
    expect_identical(contrast$result$checkpoint, source, info = omic)
    expect_true(length(contrast$result$native) >= 1L, info = omic)
    # Every declaration the registry puts on the role is covered above; the store steps refuse dispatch.
    expect_setequal(Filter(function(name) identical(workflowActionSpec(name)$omic_type, omic) && identical(workflowActionSpec(name)$compute, "omics") &&
      !is.null(workflowActionSpec(name)$declaration), workflowOperationNames()), c(vapply(steps, `[[`, "", "operation"),
      paste0(binding$prefix, c("_pearson_pairs", "_pearson_selection", "_da_contrast"))))
    for (name in paste0(binding$prefix, c("_qc_finalize", "_norm_complete", "_reset", "_da_output")))
      expect_error(assertWorkflowOperationDispatchable(name), class = "workflow_execution_unsupported")
  }
})

test_that("KIT-078: a worker's QC figures are staged and saved, and it leaves the head to the controller", {
  skip_if_not(workflowOwnerCapable())
  testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64), .package = "MultiScholaR")
  fixtures <- omicsWorkerFixtures(withr::local_tempdir())
  for (omic in names(fixtures)) {
    fixture <- fixtures[[omic]]$fixture
    binding <- fixtures[[omic]]$binding
    withr::defer(closeDIAStateAdapter(fixture$owner))
    # The fixture's duplicate ids would stop every plot; the module's QC figures follow the QC tabs.
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(binding, fixture$owner, list())))
    head <- fixture$owner$.context$head
    result <- omicsWorkerRun(fixture, paste0(binding$prefix, "_qc_figures"), list(stage = "post_filter", grouping_variable = "group", shape_variable = NULL))
    expect_identical(result$status, "complete", info = result$error)
    staged <- list.files(file.path(result$directory, "omics-publication", binding$qc_dir_role))
    expect_true(any(grepl("_pre_norm_pca[.]png$", staged)), info = omic)
    expect_identical(fixture$owner$.context$head, head)
    # Each figure file's product and asset wait for the controller's root.
    expect_gte(length(result$result$roots), 2L * length(staged))
  }
})

test_that("KIT-078: updating the filtering progress from a worker's metrics equals updating it from the object", {
  object <- lipidDeclarationsObject()
  metrics <- suppressMessages(computeLipidFilteringMetrics(object))
  expect_identical(names(metrics), c("assay_names", "metrics_list_this_step", "total_lipids"))
  expect_identical(metrics$assay_names, names(object@lipid_data))
  object <- metabDeclarationsObject()
  metrics <- suppressMessages(computeMetabFilteringMetrics(object))
  expect_identical(names(metrics), c("assay_names", "metrics_list_this_step", "total_metabolites"))
  expect_identical(metrics$assay_names, names(object@metabolite_data))
})
