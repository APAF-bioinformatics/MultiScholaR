# KIT-145 (decision M17): cyclic loess between-sample normalisation is refused below 10 features per assay in lipidomics
# and metabolomics. `limma::normalizeCyclicLoess` fits one curve per sample through that sample's measured features;
# with 2 or 3 points the fit is decided by rounding (KIT-134, Status 2026-11-18). The rule is behind
# `options(MultiScholaR.cyclic_loess_minimum = TRUE)` until the re-record switches it on: off, today's behaviour.
library(testthat)

loessMinimumSamples <- sprintf("S%02d", 1:6)

loessMinimumDesign <- function() {
  data.frame(Run = loessMinimumSamples, group = rep(c("A", "B"), each = 3L), tech_rep_group = loessMinimumSamples,
    stringsAsFactors = FALSE)
}

#' A log2-scale assay of `n_features` rows and six samples; `missing` is a two-column matrix of (row, column) set to NA.
loessMinimumAssay <- function(n_features, seed, missing = NULL) {
  values <- withr::with_seed(seed, matrix(rnorm(n_features * 6L, 15, 2), nrow = n_features) + rep(rnorm(6L, 0, 0.4), each = n_features))
  colnames(values) <- loessMinimumSamples
  if (!is.null(missing)) values[missing] <- NA_real_
  assay <- tibble::as_tibble(values)
  assay[["feature"]] <- sprintf("F%03d", seq_len(n_features))
  assay[c("feature", loessMinimumSamples)]
}

loessMinimumCases <- function() {
  list(
    metab = list(noun = "metabolites", slot = "metabolite_data", omic_type = "metabolomics", binding = metabolomicsDeclarationBinding,
      make = function(assays) methods::new("MetaboliteAssayData", metabolite_data = assays, design_matrix = loessMinimumDesign(),
        sample_id = "Run", group_id = "group", metabolite_id_column = "feature", technical_replicate_id = "tech_rep_group", args = list())),
    lipid = list(noun = "lipids", slot = "lipid_data", omic_type = "lipidomics", binding = lipidomicsDeclarationBinding,
      make = function(assays) methods::new("LipidomicsAssayData", lipid_data = assays, design_matrix = loessMinimumDesign(),
        sample_id = "Run", group_id = "group", lipid_id_column = "feature", technical_replicate_id = "tech_rep_group", args = list())))
}

loessMinimumRun <- function(object, method = "cyclicloess") {
  suppressWarnings(suppressMessages(normaliseBetweenSamples(object, normalisation_method = method)))
}

loessMinimumValues <- function(object, case, assay = 1L) {
  as.matrix(methods::slot(object, case$slot)[[assay]][loessMinimumSamples])
}

test_that("the minimum is the named constant 10 and the rule is off unless the option switches it on", {
  expect_identical(OMICS_CYCLIC_LOESS_MIN_FEATURES, 10L)
  withr::local_options(MultiScholaR.cyclic_loess_minimum = NULL)
  expect_false(omicsCyclicLoessMinimumActive())
  withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
  expect_true(omicsCyclicLoessMinimumActive())
})

for (omic in names(loessMinimumCases())) {
  test_that(paste(omic, "rule off (today's default): three features run as they do today, no refusal"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = NULL)
    object <- case$make(list(one = loessMinimumAssay(3L, 31L)))
    expect_s4_class(loessMinimumRun(object), class(object))
  })

  test_that(paste(omic, "rule on: nine features are refused, ten run"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    below <- case$make(list(one = loessMinimumAssay(OMICS_CYCLIC_LOESS_MIN_FEATURES - 1L, 32L)))
    expect_error(loessMinimumRun(below), class = "omics_cyclic_loess_too_few_features")
    at <- case$make(list(one = loessMinimumAssay(OMICS_CYCLIC_LOESS_MIN_FEATURES, 33L)))
    normalised <- loessMinimumRun(at)
    expect_s4_class(normalised, class(at))
    # It ran: the values moved, and none was lost.
    expect_false(isTRUE(all.equal(loessMinimumValues(normalised, case), loessMinimumValues(at, case))))
    expect_true(all(is.finite(loessMinimumValues(normalised, case))))
  })

  test_that(paste(omic, "rule on: the count is per sample, of the features measured in it"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    # Twelve features, but sample S04 has only nine measured: its curve would be fitted through nine points.
    thin <- case$make(list(one = loessMinimumAssay(12L, 34L, missing = cbind(1:3, 4L))))
    error <- expect_error(loessMinimumRun(thin), class = "omics_cyclic_loess_too_few_features")
    expect_match(conditionMessage(error), "only 9", fixed = TRUE)
    expect_match(conditionMessage(error), "S04", fixed = TRUE)
    # Twelve features, two missing in each sample (different features): every sample has ten measured.
    enough <- case$make(list(one = loessMinimumAssay(12L, 34L, missing = cbind(c(1:6, 7:12), rep(1:6, 2L)))))
    expect_s4_class(loessMinimumRun(enough), class(enough))
  })

  test_that(paste(omic, "rule on: in a multi-assay object the refusal names the assay, and nothing is returned"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    object <- case$make(list(LCMS_Pos = loessMinimumAssay(14L, 35L), GCMS = loessMinimumAssay(5L, 36L)))
    error <- expect_error(loessMinimumRun(object), class = "omics_cyclic_loess_too_few_features")
    expect_match(conditionMessage(error), "GCMS", fixed = TRUE)
    expect_no_match(conditionMessage(error), "LCMS_Pos", fixed = TRUE)
    expect_identical(error$assay, "GCMS")
    expect_identical(error$fewest, 5L)
    expect_identical(error$minimum, 10L)
  })

  test_that(paste(omic, "rule on: the message says what happened, why, and which methods to choose instead"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    object <- case$make(list(one = loessMinimumAssay(4L, 37L)))
    message <- conditionMessage(expect_error(loessMinimumRun(object), class = "omics_cyclic_loess_too_few_features"))
    expect_match(message, "Cyclic loess normalisation was not run", fixed = TRUE)
    expect_match(message, paste("only 4", case$noun), fixed = TRUE)
    expect_match(message, "needs at least 10", fixed = TRUE)
    expect_match(message, "fits a curve", fixed = TRUE)
    expect_match(message, "can differ from one computer to another", fixed = TRUE)
    # The alternatives are the tab's own labels for its other between-sample methods.
    offered <- names(omicGuiConfig(case$omic_type)$norm$panels$between_sample$methods)
    for (label in setdiff(offered, "Cyclic Loess")) expect_match(message, label, fixed = TRUE)
  })

  test_that(paste(omic, "rule on: the other methods are not refused on three features"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    object <- case$make(list(one = loessMinimumAssay(3L, 38L)))
    for (method in c("quantile", "scale", "none")) expect_s4_class(loessMinimumRun(object, method), class(object))
  })

  test_that(paste(omic, "rule on: the declared operation the worker runs fails with the refusal, nothing checkpointed"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    kernel <- omicsBetweenSampleNormalisationOperation(case$binding)$kernel
    checkpointed <- 0L
    store <- list(checkpoint = function(object, name) {
      checkpointed <<- checkpointed + 1L
      object
    })
    object <- case$make(list(one = loessMinimumAssay(5L, 39L)))
    expect_error(suppressWarnings(suppressMessages(kernel(object, list(normalisation_method = "cyclicloess"), store))),
      class = "omics_cyclic_loess_too_few_features")
    expect_identical(checkpointed, 0L)
  })

  # KIT-145 item 2: the test that would have caught the defect. At the minimum, the result does not depend on rounding:
  # a relative perturbation of 1e-12 of the input leaves it unchanged to 1e-9. (On 2 or 3 features the same perturbation
  # moved values by up to 0.71 log2 and gave 84 distinct results in 200 draws.)
  test_that(paste(omic, "at the minimum the result is unchanged under a 1e-12 relative perturbation of the input"), {
    case <- loessMinimumCases()[[omic]]
    withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
    assay <- loessMinimumAssay(OMICS_CYCLIC_LOESS_MIN_FEATURES, 40L)
    reference <- loessMinimumValues(loessMinimumRun(case$make(list(one = assay))), case)
    worst <- withr::with_seed(41L, max(vapply(seq_len(25L), function(draw) {
      perturbed <- assay
      perturbed[loessMinimumSamples] <- as.matrix(assay[loessMinimumSamples]) * (1 + stats::runif(nrow(assay) * 6L, -1e-12, 1e-12))
      max(abs(loessMinimumValues(loessMinimumRun(case$make(list(one = perturbed))), case) - reference))
    }, numeric(1))))
    expect_lt(worst, 1e-9)
  })
}

test_that("the tab's in-process run reports the refusal in its log and as an error notification, and commits no normalised state", {
  skip_if_not(exists("module_ci_metab_norm_object", mode = "function"), "the MCI-015 helper is not loaded")
  withr::local_options(MultiScholaR.cyclic_loess_minimum = TRUE)
  object <- module_ci_metab_norm_object(positive_only = TRUE)
  object@metabolite_data <- lapply(object@metabolite_data, function(assay) assay[seq_len(5L), , drop = FALSE])
  workflow <- module_ci_metab_norm_workflow(object)
  norm_data <- module_ci_metab_norm_data(object)
  root <- withr::local_tempdir()
  paths <- list(metabolite_qc_dir = file.path(root, "qc"), source_dir = file.path(root, "source"), export_dir = file.path(root, "export"))
  rec <- new.env(parent = emptyenv())
  rec$log <- character()
  rec$notes <- character()
  local_mocked_bindings(generateMetabQcPlots = function(...) invisible(NULL), .package = "MultiScholaR")
  before <- as.character(workflow$state_manager$getHistory())
  result <- suppressWarnings(suppressMessages(runWorkflowOmicsInProcessNormalisation(metabolomicsDeclarationBinding, workflow$state_manager,
    list(apply_itsd = FALSE, itsd_aggregation = "median", itsd_feature_ids = NULL, log_offset = 1, norm_method = "cyclicloess",
      ruv_mode = "skip", ruv_parameters = list(ruv_mode = "skip"), grouping_variable = "group", shape_variable = "batch",
      assay_names = norm_data$assay_names, qc_dir = paths[[1L]], experiment_paths = paths),
    norm_data, workflow$config_list, addLog = function(entry) rec$log <- c(rec$log, entry), advance = function(detail) NULL,
    notify = function(ui, ..., type = "default") rec$notes <- c(rec$notes, paste0("[", type, "] ", ui)),
    finish = function() rec$log <- c(rec$log, "finished"), itsd_table = NULL)))
  expect_false(result)
  expect_length(rec$notes, 1L)
  expect_match(rec$notes, "^\\[error\\] Error: .*Cyclic loess normalisation was not run")
  expect_match(rec$notes, "Quantile", fixed = TRUE)
  expect_true(any(grepl("^ERROR: .*Cyclic loess normalisation was not run", rec$log)))
  expect_false("finished" %in% rec$log)
  expect_false(isTRUE(norm_data$normalization_complete))
  # The log2 step before it is committed, as for any later failure; no normalised state follows it.
  after <- as.character(workflow$state_manager$getHistory())
  expect_false("metab_normalized" %in% after)
  expect_identical(setdiff(after, before), "metab_log2")
})
