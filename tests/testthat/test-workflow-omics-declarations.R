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
  # KIT-078 adds the three publication operations per omic and moves the compute steps to the worker's `omics` role, with
  # their handler, validator and, where they write no section, `section = NULL`; every other field is the frozen one.
  for (omic in c("lipidomics", "metabolomics")) {
    prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
    current$registry[[omic]][paste0(prefix, c("_qc_figures", "_composite_figure", "_da_publication",
      # KIT-082: what the tabs read from the whole matrix outside their steps.
      "_da_tab_results", "_filtering_metrics", "_final_pca", "_volcano_page", "_session_export", "_integration_export"))] <- NULL
    # KIT-079: three steps gain a partitioned variant, `whole` staying the default; KIT-089: except the intensity filter.
    for (step in paste0(prefix, c("_intensity_filter", "_log2", "_duplicate_resolution"))) {
      expect_identical(current$registry[[omic]][[step]]$entry$variants, c("whole", "partitioned"), info = step)
      expect_identical(current$registry[[omic]][[step]]$entry$default_variant,
        if (identical(step, paste0(prefix, "_intensity_filter"))) "partitioned" else "whole", info = step)
      current$registry[[omic]][[step]]$entry$variants <- "whole"
      current$registry[[omic]][[step]]$entry$default_variant <- "whole"
      current$registry[[omic]][[step]]$declaration$partitioned <- NULL
    }
    for (name in names(current$registry[[omic]])) {
      if (!identical(current$registry[[omic]][[name]]$entry$compute, "omics")) next
      entry <- current$registry[[omic]][[name]]$entry
      expect_identical(c(entry$completion, entry$publication), c("omics", "omics"), info = name)
      entry$section <- NULL
      entry[c("compute", "handler", "validate", "completion", "publication")] <- frozen$registry[[omic]][[name]]$entry[c("compute", "handler", "validate", "completion", "publication")]
      entry["handler"] <- list(NULL)
      current$registry[[omic]][[name]]$entry <- entry[names(frozen$registry[[omic]][[name]]$entry)]
    }
  }
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
  # KIT-095: the transcriptomics binding is the third.
  expect_identical(names(bindings), c("lipidomics", "metabolomics", "transcriptomics"))
  expect_identical(vapply(bindings, `[[`, "", "prefix"), c(lipidomics = "lipid", metabolomics = "metab", transcriptomics = "trx"))
  expect_identical(runLipidReferenceOperation, runWorkflowOmicsReferenceOperation)
  expect_identical(markLipidReferenceStoreStep, markWorkflowOmicsStoreStep)
  expect_identical(commitLipidReferenceStoreStep, commitWorkflowOmicsStoreStep)
  expect_identical(lipidomicsStepLabels, bindings$lipidomics$step_labels)
  expect_identical(metabolomicsStepLabels, bindings$metabolomics$step_labels)
  expect_identical(names(lipidomicsWorkflowOperations()), names(omicsWorkflowOperations(bindings$lipidomics)))
  expect_error(replayWorkflowOmicsStep(bindings$lipidomics, NULL, list(operation = "metab_log2", state = "x"), list()),
    class = "workflow_replay_unsupported")
})
