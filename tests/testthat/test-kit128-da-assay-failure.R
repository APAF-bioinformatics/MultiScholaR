# KIT-128 (parked row 24 b, decision D4): a kernel error on one assay is named in the result and warned about, the
# other assays are kept, and the run stops only when every assay failed. Three sites: runLipidsDA, runMetabolitesDA and
# the per-contrast path's fit half (omicsDaAssayFit, the GUI's path). KIT-196: a feature without a value in any sample is
# no longer a failure (it is reported as not tested), so the orchestrator's failure is made by a mocked contrast half.
library(testthat)

kit128FailureRun <- function(omic, object) {
  suppressMessages(utils::capture.output(result <- omicsDaKernel(omic)(theObject = object,
    contrasts_tbl = omicsDaContrastsTable("groupB-groupA"), formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
    treat_lfc_cutoff = 0)))
  result
}

# KIT-196: the orchestrator's contrast half raises for the assay whose features start with `failing_prefix`.
kit128ContrastFailureRun <- function(omic, failing_prefix) {
  half <- paste0("apply", if (identical(omic, "lipid")) "Lipid" else "Metab", "DaContrasts")
  original <- get(half, envir = asNamespace("MultiScholaR"))
  failing <- function(fit, ...) {
    if (any(startsWith(rownames(fit$coefficients), failing_prefix))) stop("the contrast half could not test this assay")
    original(fit, ...)
  }
  do.call(testthat::local_mocked_bindings, c(stats::setNames(list(failing), half), list(.package = "MultiScholaR")))
  kit128FailureRun(omic, omicsDaContrastsFixture(omic))
}

for (omic in c("lipid", "metab")) {
  long_table <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"

  test_that(paste("KIT-128:", omic, "orchestrator names a failed assay, warns and keeps the others"), {
    healthy <- kit128FailureRun(omic, omicsDaContrastsFixture(omic))
    expect_null(healthy$failed_assays)
    # The `neg` assay's features are F2_*.
    expect_warning(result <- kit128ContrastFailureRun(omic, "F2_"), "neg", class = "omics_da_assay_failed")
    expect_identical(names(result$failed_assays), "neg")
    expect_match(result$failed_assays$neg, "the contrast half could not test this assay", fixed = TRUE)
    expect_identical(names(result$per_assay_results), "pos")
    expect_identical(result$per_assay_results$pos, healthy$per_assay_results$pos)
    expect_identical(result[[long_table]], healthy$per_assay_results$pos)
  })

  test_that(paste("KIT-128:", omic, "orchestrator stops when every assay failed"), {
    expect_error(suppressWarnings(kit128ContrastFailureRun(omic, "F")),
      "the contrast half could not test this assay", class = "omics_da_all_assays_failed")
  })

  test_that(paste("KIT-128:", omic, "per-contrast path records a failed fit under the assay's errors and keeps the others"), {
    halves <- paste0("fit", if (identical(omic, "lipid")) "Lipid" else "Metab", "DaAssayModel")
    original <- get(halves, envir = asNamespace("MultiScholaR"))
    test_env <- environment()
    run <- function(failing_prefix) {
      failing <- function(data, ...) {
        if (any(startsWith(rownames(data), failing_prefix))) stop("the fit half could not fit this assay")
        original(data, ...)
      }
      do.call(testthat::local_mocked_bindings, c(stats::setNames(list(failing), halves), list(.package = "MultiScholaR")))
      # The owner stays open until the test ends (the record is read after the run).
      workflow <- omicsTestOwner(omicsDaContrastsFixture(omic), "da_input", env = test_env)
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

  test_that(paste("KIT-128:", omic, "orchestrator names the root cause of an error raised inside the formatting step"), {
    formatter <- paste0("format", if (identical(omic, "lipid")) "Lipid" else "Metab", "DaContrastTable")
    original <- get(formatter, envir = asNamespace("MultiScholaR"))
    run <- function(failing) {
      failing_formatter <- function(da_tbl, contrast_name, assay_name, ...) {
        if (assay_name %in% failing) stop("the table could not be formatted")
        original(da_tbl, contrast_name, assay_name, ...)
      }
      do.call(testthat::local_mocked_bindings, c(stats::setNames(list(failing_formatter), formatter), list(.package = "MultiScholaR")))
      kit128FailureRun(omic, omicsDaContrastsFixture(omic))
    }
    # The formatting step runs under purrr, which wraps the error; the record names the cause, not purrr's index.
    expect_warning(result <- run("neg"), "neg: the table could not be formatted", class = "omics_da_assay_failed")
    expect_identical(result$failed_assays, list(neg = "the table could not be formatted"))
    expect_identical(names(result$per_assay_results), "pos")
    # K128-5: a failed assay leaves no fit in the result either.
    expect_identical(names(result$contrasts_results), "pos")
    expect_error(run(c("pos", "neg")), "pos: the table could not be formatted", class = "omics_da_all_assays_failed")
  })
}
