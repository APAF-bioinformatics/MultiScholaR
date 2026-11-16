# KIT-128 (parked row 24 b, decision D4): a kernel error on one assay is named in the result and warned about, the
# other assays are kept, and the run stops only when every assay failed. Three sites: runLipidsDA, runMetabolitesDA and
# the per-contrast path's fit half (omicsDaAssayFit, the GUI's path). The reproduced case: a feature without a value in
# any sample of the assay stops limma's eBayes(trend = TRUE) ("covariate contains NA or infinite values").
library(testthat)

kit128FailureObject <- function(omic, assays = "neg") {
  object <- omicsDaContrastsFixture(omic)
  slot_name <- if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
  tables <- methods::slot(object, slot_name)
  for (assay in assays) tables[[assay]][1L, object@design_matrix$Run] <- NA_real_
  methods::slot(object, slot_name) <- tables
  object
}

kit128FailureRun <- function(omic, object) {
  suppressMessages(utils::capture.output(result <- omicsDaKernel(omic)(theObject = object,
    contrasts_tbl = omicsDaContrastsTable("groupB-groupA"), formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0)))
  result
}

for (omic in c("lipid", "metab")) {
  long_table <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"

  test_that(paste("KIT-128:", omic, "orchestrator names a failed assay, warns and keeps the others"), {
    healthy <- kit128FailureRun(omic, omicsDaContrastsFixture(omic))
    expect_null(healthy$failed_assays)
    expect_warning(result <- kit128FailureRun(omic, kit128FailureObject(omic)), "neg", class = "omics_da_assay_failed")
    expect_identical(names(result$failed_assays), "neg")
    expect_match(result$failed_assays$neg, "covariate contains NA or infinite values", fixed = TRUE)
    expect_identical(names(result$per_assay_results), "pos")
    expect_identical(result$per_assay_results$pos, healthy$per_assay_results$pos)
    expect_identical(result[[long_table]], healthy$per_assay_results$pos)
  })

  test_that(paste("KIT-128:", omic, "orchestrator stops when every assay failed"), {
    expect_error(suppressWarnings(kit128FailureRun(omic, kit128FailureObject(omic, c("pos", "neg")))),
      "covariate contains NA or infinite values", class = "omics_da_all_assays_failed")
  })

  test_that(paste("KIT-128:", omic, "per-contrast path records a failed fit under the assay's errors and keeps the others"), {
    halves <- paste0("fit", if (identical(omic, "lipid")) "Lipid" else "Metab", "DaAssayModel")
    original <- get(halves, envir = asNamespace("MultiScholaR"))
    run <- function(failing_prefix) {
      failing <- function(data, ...) {
        if (any(startsWith(rownames(data), failing_prefix))) stop("the fit half could not fit this assay")
        original(data, ...)
      }
      do.call(testthat::local_mocked_bindings, c(stats::setNames(list(failing), halves), list(.package = "MultiScholaR")))
      workflow <- omicsTestOwner(omicsDaContrastsFixture(omic), "da_input")
      manager <- workflow$state_manager
      wrapper <- if (identical(omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
      suppressWarnings(suppressMessages(utils::capture.output(out <- wrapper(manager, getWorkflowStateReference(manager),
        list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), omicsDaContrastsTable("groupB-groupA")))))
      list(run = out, manager = manager)
    }
    # The `neg` assay's features are F2_*.
    out <- run("F2_")
    expect_identical(out$run$manifest$results[["groupB-groupA"]]$errors, list(neg = "the fit half could not fit this assay"))
    expect_identical(out$run$results$failed_contrasts, list(`groupB-groupA` = list(neg = "the fit half could not fit this assay")))
    expect_identical(names(out$run$results$per_assay_results), "pos")
    expect_identical(nrow(out$run$results[[long_table]]), 60L)
    failures <- buildWorkflowAnalysisRecord(out$manager)$analyses$differential_abundance$failures
    expect_identical(vapply(failures, function(f) f$assay, ""), "neg")
    # Every assay failed: the run stops with the message.
    expect_error(run("F"), "the fit half could not fit this assay")
  })
}
