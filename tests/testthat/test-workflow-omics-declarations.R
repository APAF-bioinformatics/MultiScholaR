# KIT-073: the shared omics declaration core changes nothing observable. The dump of the registry entries, the
# declaration fields and both omics' wrapper journeys equals the one frozen on the base commit before the core
# (tests/testthat/fixtures/kit073-omics-declarations-frozen.rds; docs/qa/workflow-execution/kit/kit-073/freeze-dump.R).
library(testthat)

test_that("KIT-073: the omics declarations and their journeys equal the frozen dump of the two files the core replaced", {
  skip_if_not(workflowOwnerCapable())
  frozen_file <- test_path("fixtures", "kit073-omics-declarations-frozen.rds")
  skip_if_not(file.exists(frozen_file))
  frozen <- readRDS(frozen_file)
  root <- withr::local_tempdir()
  current <- suppressWarnings(suppressMessages(omicsDeclarationsDump(root)))
  # KIT-076 adds one declaration per omic (`<prefix>_da_contrast`); the twelve frozen ones are unchanged.
  current$registry$lipidomics$lipid_da_contrast <- NULL
  current$registry$metabolomics$metab_da_contrast <- NULL
  expect_identical(names(current$registry$lipidomics), names(frozen$registry$lipidomics))
  expect_identical(names(current$registry$metabolomics), names(frozen$registry$metabolomics))
  expect_identical(current$registry, frozen$registry)
  for (omic in c("lipidomics", "metabolomics")) {
    expect_identical(current[[omic]]$history, frozen[[omic]]$history, info = omic)
    for (state in names(frozen[[omic]]$states)) {
      expect_identical(current[[omic]]$states[[state]], frozen[[omic]]$states[[state]], info = paste(omic, state))
    }
    expect_identical(current[[omic]]$values, frozen[[omic]]$values, info = omic)
    expect_identical(current[[omic]]$correlations, frozen[[omic]]$correlations, info = omic)
    # KIT-076: the manifest is v2 and `fit.eb` moved per contrast; the tables, the counts and the warnings are the frozen ones.
    products <- setdiff(names(frozen[[omic]]$da$results), "contrasts_results")
    expect_identical(current[[omic]]$da$results[products], frozen[[omic]]$da$results[products], info = omic)
  }
})

test_that("KIT-073: the core serves both omics from one binding each, and the omic names stay one body", {
  bindings <- omicsDeclarationBindings()
  expect_identical(names(bindings), c("lipidomics", "metabolomics"))
  expect_identical(vapply(bindings, `[[`, "", "prefix"), c(lipidomics = "lipid", metabolomics = "metab"))
  expect_identical(runLipidReferenceOperation, runWorkflowOmicsReferenceOperation)
  expect_identical(markLipidReferenceStoreStep, markWorkflowOmicsStoreStep)
  expect_identical(commitLipidReferenceStoreStep, commitWorkflowOmicsStoreStep)
  expect_identical(lipidomicsStepLabels, bindings$lipidomics$step_labels)
  expect_identical(metabolomicsStepLabels, bindings$metabolomics$step_labels)
  expect_identical(names(lipidomicsWorkflowOperations()), names(omicsWorkflowOperations(bindings$lipidomics)))
  expect_error(replayWorkflowOmicsStep(bindings$lipidomics, NULL, list(operation = "metab_log2", state = "x"), list()),
    class = "workflow_replay_unsupported")
})
