# KIT-128: declared technical replicates enter the lipid and metab DA as a block (duplicateCorrelation, then
# lmFit(block, correlation); limma User's Guide 11.3 and 9.7), through the proteomics rule
# (resolveProtTechnicalReplicateBlock). The fixture is generated here: 2 groups, 4 biological units per group, 2 injections
# per unit, no true effect, 75% of the residual variance biological (KIT-124's C6 recipe). Size of the null-rate test:
# 20 seeds by 500 features (two assays of 250) = 10,000 tests per omic; the limit is KIT-124's rule 1, the one-sided
# 97.5% exact binomial upper bound of the share of P < 0.05 at most 0.075.
library(testthat)

kit128Design <- function(replicates = c("technical", "unique", "spanned")) {
  replicates <- match.arg(replicates)
  groups <- rep(c("A", "B"), each = 8L)
  unit <- rep(rep(seq_len(4L), each = 2L), times = 2L)
  runs <- sprintf("%s%d_%d", groups, unit, rep(seq_len(2L), times = 8L))
  data.frame(Run = runs, group = groups, replicates = unit,
    tech_rep_group = switch(replicates, technical = paste(groups, unit, sep = "_"), unique = runs, spanned = groups),
    stringsAsFactors = FALSE)
}

#' One assay's table: level + a biological term per unit (variance 0.1875) + a technical term per injection (0.0625).
kit128Assay <- function(seed, id, annotation, features = 250L) {
  design <- kit128Design()
  unit <- match(design$tech_rep_group, unique(design$tech_rep_group))
  values <- withr::with_seed(seed, {
    level <- stats::rnorm(features, 20, 2)
    biological <- matrix(stats::rnorm(features * 8L, 0, sqrt(0.1875)), nrow = features)
    level + biological[, unit] + matrix(stats::rnorm(features * 16L, 0, 0.25), nrow = features)
  })
  colnames(values) <- design$Run
  table <- data.frame(sprintf("F%d_%03d", seed, seq_len(features)), paste("Name", seq_len(features)), stringsAsFactors = FALSE)
  names(table) <- c(id, annotation)
  cbind(table, as.data.frame(values, check.names = FALSE))
}

kit128Object <- function(omic, seed, replicates = "technical", technical_replicate_id = "tech_rep_group") {
  design <- kit128Design(replicates)
  if (identical(omic, "lipid")) {
    methods::new("LipidomicsAssayData",
      lipid_data = list(pos = kit128Assay(seed, "lipid_id", "lipid_name"), neg = kit128Assay(seed + 1000L, "lipid_id", "lipid_name")),
      design_matrix = design, sample_id = "Run", group_id = "group", technical_replicate_id = technical_replicate_id,
      lipid_id_column = "lipid_id", annotation_id_column = "lipid_name")
  } else {
    methods::new("MetaboliteAssayData",
      metabolite_data = list(pos = kit128Assay(seed, "metabolite_id", "metabolite_name"),
        neg = kit128Assay(seed + 1000L, "metabolite_id", "metabolite_name")),
      design_matrix = design, sample_id = "Run", group_id = "group", technical_replicate_id = technical_replicate_id,
      metabolite_id_column = "metabolite_id", annotation_id_column = "metabolite_name")
  }
}

kit128Contrast <- function() omicsDaContrastsTable("groupB-groupA")

kit128Run <- function(omic, object) {
  suppressWarnings(suppressMessages(utils::capture.output(result <- omicsDaKernel(omic)(theObject = object,
    contrasts_tbl = kit128Contrast(), formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0))))
  result
}

#' The GUI's path: one declaration per contrast on the saved fit, on a session owner (as test-workflow-omics-da-contrasts.R).
kit128ContrastRun <- function(omic, object, env = parent.frame()) {
  workflow <- omicsTestOwner(object, "da_input", env = env)
  manager <- workflow$state_manager
  wrapper <- if (identical(omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
  run <- suppressWarnings(suppressMessages(wrapper(manager, getWorkflowStateReference(manager),
    list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), kit128Contrast())))
  list(manager = manager, run = run)
}

#' The route's unblocked test on one assay table (the kernel before KIT-128): lmFit, contrasts.fit, eBayes(trend, robust).
kit128UnblockedP <- function(table, design) {
  values <- as.matrix(table[, design$Run])
  rownames(values) <- table[[1L]]
  model <- stats::model.matrix(~ 0 + group, design)
  fit <- limma::contrasts.fit(limma::lmFit(values, model), limma::makeContrasts(contrasts = "groupB-groupA", levels = model))
  limma::eBayes(fit, trend = TRUE, robust = TRUE)$p.value[, 1L]
}

for (omic in c("lipid", "metab")) {
  long_table <- if (identical(omic, "lipid")) "da_lipids_long" else "da_metabolites_long"
  id_column <- if (identical(omic, "lipid")) "lipid_id" else "metabolite_id"
  data_slot <- if (identical(omic, "lipid")) "lipid_data" else "metabolite_data"
  fit_half <- if (identical(omic, "lipid")) fitLipidDaAssayModel else fitMetabDaAssayModel

  test_that(paste("KIT-128:", omic, "null rate at 0.05 on duplicated injections is within the limit (upper bound at most 0.075)"), {
    p <- unlist(lapply(seq_len(20L), function(seed) kit128Run(omic, kit128Object(omic, seed))[[long_table]]$raw_pvalue))
    expect_length(p, 10000L)
    bound <- stats::binom.test(sum(p < 0.05), length(p), alternative = "less", conf.level = 0.975)$conf.int[[2L]]
    expect_lte(bound, 0.075, label = sprintf("null rate %.4f, its 97.5%% upper bound %.4f", mean(p < 0.05), bound))
  })

  test_that(paste("KIT-128:", omic, "the fit records the block, and the per-contrast path equals the orchestrator"), {
    object <- kit128Object(omic, 1L)
    direct <- kit128Run(omic, object)
    fit <- direct$contrasts_results$pos$fit.eb
    # limma keeps the extra components through contrasts.fit and eBayes.
    expect_identical(fit$block_reason, "technical_replicates")
    expect_true(fit$uses_duplicate_correlation)
    expect_gt(fit$consensus_correlation, 0.5)
    out <- kit128ContrastRun(omic, object)
    expect_identical(out$run$results[[long_table]], direct[[long_table]])
    expect_identical(out$run$manifest$model$block_reason, "technical_replicates")
    expect_true(out$run$manifest$model$uses_duplicate_correlation)
    record <- buildWorkflowAnalysisRecord(out$manager)
    expect_identical(record$analyses$differential_abundance$model$block_reason, "technical_replicates")
    expect_true(validateWorkflowAnalysisRecord(record)$valid)
  })

  test_that(paste("KIT-128:", omic, "a design without technical replicates is fitted as before, with the reason recorded"), {
    for (case in list(c("unique", "no_duplicates"), c("spanned", "spanned_by_design"))) {
      object <- kit128Object(omic, 1L, replicates = case[[1L]])
      result <- kit128Run(omic, object)
      for (assay in c("pos", "neg")) {
        table <- methods::slot(object, data_slot)[[assay]]
        long <- result$per_assay_results[[assay]]
        expect_identical(long$raw_pvalue[match(table[[id_column]], long[[id_column]])],
          unname(kit128UnblockedP(table, object@design_matrix)), info = paste(case[[1L]], assay))
        expect_identical(result$contrasts_results[[assay]]$fit.eb$block_reason, case[[2L]])
        expect_false(result$contrasts_results[[assay]]$fit.eb$uses_duplicate_correlation)
      }
      manifest <- kit128ContrastRun(omic, object)$run$manifest
      expect_identical(manifest$model$block_reason, case[[2L]])
      expect_false(manifest$model$uses_duplicate_correlation)
    }
    # No declared column: today's test, reason `no_declaration`.
    undeclared <- kit128Run(omic, kit128Object(omic, 1L, technical_replicate_id = NA_character_))
    unique_blocks <- kit128Run(omic, kit128Object(omic, 1L, replicates = "unique"))
    expect_identical(undeclared[[long_table]], unique_blocks[[long_table]])
    expect_identical(undeclared$contrasts_results$pos$fit.eb$block_reason, "no_declaration")
  })

  test_that(paste("KIT-128:", omic, "the fit half without the two names is the unblocked lmFit, unchanged"), {
    object <- kit128Object(omic, 1L)
    table <- methods::slot(object, data_slot)$pos
    values <- as.matrix(table[, object@design_matrix$Run])
    rownames(values) <- table[[id_column]]
    fit <- suppressMessages(fit_half(values, object@design_matrix, "~ 0 + group", "Run"))
    plain <- limma::lmFit(values, design = fit$design)
    expect_identical(fit, plain)
    blocked <- suppressMessages(fit_half(values, object@design_matrix, "~ 0 + group", "Run", group_id = "group",
      technical_replicate_id = "tech_rep_group"))
    # Balanced and complete: the block changes the standard errors, not the group means.
    expect_equal(blocked$coefficients, plain$coefficients, tolerance = 1e-12)
    expect_true(all(blocked$stdev.unscaled > plain$stdev.unscaled))
  })

  test_that(paste("KIT-128:", omic, "the fit identity carries the block: a fit saved without it is not a hit"), {
    object <- kit128Object(omic, 1L)
    table <- methods::slot(object, data_slot)$pos
    values <- as.matrix(table[, object@design_matrix$Run])
    rownames(values) <- table[[id_column]]
    software <- "kit128"
    without <- newWorkflowOmicsDaFitIdentity(values, object@design_matrix, "~ 0 + group", "pos", software)
    with_block <- newWorkflowOmicsDaFitIdentity(values, object@design_matrix, "~ 0 + group", "pos", software,
      block = object@design_matrix$tech_rep_group)
    expect_false(identical(without, with_block))
    # A design the rule does not call technical keeps the identity it had.
    expect_identical(newWorkflowOmicsDaFitIdentity(values, object@design_matrix, "~ 0 + group", "pos", software, block = NULL), without)
  })
}
