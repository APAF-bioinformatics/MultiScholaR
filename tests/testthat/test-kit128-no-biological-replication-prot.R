# KIT-128 K128-1 (the user's rule, 2026-11-17): a design whose every group is one biological unit measured more than once
# (the rule's `spanned_by_design`) is refused in every omic, proteomics included: no biological replication.
library(testthat)

kit128ProtData <- function() {
  data <- matrix(c(10, 10.2, 10.5, 10.4, 12, 12.1, 11.8, 11.9, 15, 14.8, 15.1, 15.3, 13, 13.2, 12.9, 13.1), nrow = 2L,
    byrow = TRUE, dimnames = list(c("P1", "P2"), paste0("S", 1:8)))
  data
}

kit128ProtDesign <- function(replicates) {
  design <- data.frame(Run = paste0("S", 1:8), group = rep(c("A", "B"), each = 4L), replicates = replicates,
    stringsAsFactors = FALSE)
  rownames(design) <- design$Run
  design
}

test_that("KIT-128 K128-1: proteomics refuses a design whose every group is one unit measured four times", {
  expect_error(suppressMessages(runTestsContrasts(data = kit128ProtData(), contrast_strings = "groupB-groupA",
    design_matrix = kit128ProtDesign(c(1, 1, 1, 1, 2, 2, 2, 2)), formula_string = "~ 0 + group")),
    "groups A, B", class = "omics_da_no_biological_replication")
  # Two biological units per group, each measured twice: fitted with the block, as before.
  results <- suppressMessages(runTestsContrasts(data = kit128ProtData(), contrast_strings = "groupB-groupA",
    design_matrix = kit128ProtDesign(c(1, 1, 2, 2, 3, 3, 4, 4)), formula_string = "~ 0 + group"))
  expect_identical(nrow(results$results[[1L]]), 2L)
})

test_that("KIT-128 K128-1: the DIA_LIMPA route refuses it before dpcDE (the KIT-053 branch, dpcDE mocked)", {
  skip_if_not_installed("limpa")
  ids <- sprintf("P%02d", 1:30)
  samples <- paste0("S", 1:8)
  E <- withr::with_seed(128L, matrix(stats::rnorm(240, 20, 1), 30, 8, dimnames = list(ids, samples)))
  nobs <- matrix(2L, 30, 8, dimnames = dimnames(E))
  elist <- methods::new("EList", list(E = E, genes = data.frame(NPrec = rep(2L, 30L), PropObs = rep(1, 30L), row.names = ids),
    other = list(standard.error = matrix(0.2, 30, 8, dimnames = dimnames(E)), n.observations = nobs)))
  table <- data.frame(Protein.Ids = ids, as.data.frame(E, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(table) <- NULL
  object_with <- function(replicates) {
    ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids",
      protein_id_table = data.frame(Protein.Ids = ids, stringsAsFactors = FALSE),
      design_matrix = data.frame(Run = samples, group = rep(c("control", "treated"), each = 4L), replicates = replicates,
        stringsAsFactors = FALSE),
      sample_id = "Run", group_id = "group", technical_replicate_id = "replicates",
      args = list(limpa_dpc_quant_results = list(quantified_elist = elist, standard_errors = elist$other$standard.error)))
  }
  calls <- new.env(parent = emptyenv())
  calls$dpcDE <- 0L
  testthat::local_mocked_bindings(dpcDE = function(y, design, plot = TRUE, ...) {
    calls$dpcDE <- calls$dpcDE + 1L
    limma::lmFit(y$E, design)
  }, .package = "limpa")
  run <- function(object) suppressMessages(differentialAbundanceAnalysis(object,
    contrasts_tbl = data.frame(contrasts = "treated_vs_control=grouptreated-groupcontrol", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
  # Each group is one unit measured four times: refused, and dpcDE never ran.
  expect_error(run(object_with(rep(1L, 8L))), "groups control, treated", class = "omics_da_no_biological_replication")
  expect_identical(calls$dpcDE, 0L)
  # Two units per group, each measured twice: the branch fits with the block, as before.
  expect_identical(run(object_with(rep(c(1L, 1L, 2L, 2L), 2L)))$da_model$block_reason, "technical_replicates")
  expect_identical(calls$dpcDE, 1L)
})
