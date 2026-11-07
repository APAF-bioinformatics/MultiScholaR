# fidelity-coverage-compare: shared
library(testthat)

ProteinQuantitativeData <- get(
  "ProteinQuantitativeData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

newProteinNormObject <- function(values = data.frame(
                                   S1 = c(10, 20, NA_real_),
                                   S2 = c(12, 18, 30),
                                   S3 = c(15, 21, 33),
                                   check.names = FALSE
                                 ),
                                 args = list()) {
  ProteinQuantitativeData(
    protein_quant_table = data.frame(
      Protein.Ids = paste0("P", seq_len(nrow(values))),
      values,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = colnames(values),
      group = c("A", "A", "B"),
      replicate = c("R1", "R2", "R3"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate",
    protein_id_column = "Protein.Ids",
    args = args
  )
}

test_that("protein S4 normalisation methods preserve switch routing and RUV helper orchestration", {
    skip("quarantined Q-18 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
  local_mocked_bindings(
    imputeNipalsMatrix = function(x, ncomp) {
      x[is.na(x)] <- 99
      x
    },
    .package = "MultiScholaR"
  )

  overrides <- list(
    checkParamsObjectFunctionSimplify = function(theObject, name, default) {
      value <- theObject@args[[name]]
      if (is.null(value)) default else value
    },
    updateParamInObject = function(theObject, name) theObject,
    cleanDesignMatrix = function(theObject) theObject,
    normalizeCyclicLoess = function(x) x + 1,
    normalizeQuantiles = function(x) x + 2,
    normalizeMedianAbsValues = function(x) x + 3,
    ruv_cancorplot = function(Y, X, ctl) {
      list(Y = Y, X = X, ctl = ctl)
    },
    getRuvIIIReplicateMatrixHelper = function(design_matrix, sample_id_sym, grouping_sym) {
      matrix(
        c(1, 0, 0,
          0, 1, 0,
          0, 0, 1),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(design_matrix$Run, design_matrix$Run)
      )
    },
    RUVIII_C_Varying = function(k, Y, M, toCorrect, potentialControls) {
      matrix(
        c(
          10, 20, NA_real_,
          11, 21, NA_real_,
          12, 22, NA_real_
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(rownames(Y), colnames(Y))
      )
    }
  )

  normalise_method <- makeFunctionWithOverrides(
    methods::selectMethod("normaliseBetweenSamples", "ProteinQuantitativeData"),
    overrides
  )
  cancor_method <- makeFunctionWithOverrides(
    methods::selectMethod("ruvCancor", "ProteinQuantitativeData"),
    overrides
  )
  ruviii_method <- makeFunctionWithOverrides(
    methods::selectMethod("ruvIII_C_Varying", "ProteinQuantitativeData"),
    overrides
  )

  cyclic_obj <- newProteinNormObject(args = list(normalisation_method = "cyclicloess"))
  cyclic_norm <- normalise_method(cyclic_obj)
  expect_equal(cyclic_norm@protein_quant_table$S1[1], 11)

  quantile_obj <- newProteinNormObject(args = list(normalisation_method = "quantile"))
  quantile_norm <- normalise_method(quantile_obj)
  expect_equal(quantile_norm@protein_quant_table$S1[1], 12)

  scale_obj <- newProteinNormObject(args = list(normalisation_method = "scale"))
  scale_norm <- normalise_method(scale_obj)
  expect_equal(scale_norm@protein_quant_table$S1[1], 13)

  none_obj <- newProteinNormObject(args = list(normalisation_method = "none"))
  none_norm <- normalise_method(none_obj)
  expect_equal(none_norm@protein_quant_table$S1[1], 10)
  expect_true(is.na(none_norm@protein_quant_table$S1[3]))

  cancor_obj <- newProteinNormObject(args = list(
    ctrl = c(P1 = TRUE, P2 = TRUE, P3 = TRUE, P4 = TRUE, P5 = TRUE),
    num_components_to_impute = 2,
    ruv_grouping_variable = "group"
  ))
  cancor_plot <- cancor_method(cancor_obj)
  expect_true(is.list(cancor_plot))
  expect_equal(cancor_plot$Y["S1", "P3"], 99)
  expect_identical(cancor_plot$X, c("A", "A", "B"))
  expect_identical(unname(cancor_plot$ctl), rep(TRUE, 5))

  expect_error(
    cancor_method(newProteinNormObject(args = list(
      ctrl = rep(TRUE, 5),
      num_components_to_impute = 2,
      ruv_grouping_variable = "missing_group"
    ))),
    "is not a column in the design matrix",
    fixed = TRUE
  )
  expect_error(
    cancor_method(newProteinNormObject(args = list(
      ctrl = rep(TRUE, 5),
      num_components_to_impute = 0,
      ruv_grouping_variable = "group"
    ))),
    "value is invalid",
    fixed = TRUE
  )
  expect_error(
    cancor_method(newProteinNormObject(args = list(
      ctrl = c(TRUE, TRUE, TRUE, TRUE),
      num_components_to_impute = 2,
      ruv_grouping_variable = "group"
    ))),
    "less than 5",
    fixed = TRUE
  )

  ruv_obj <- newProteinNormObject(args = list(
    ctrl = c(P1 = TRUE, P2 = FALSE, P3 = TRUE),
    ruv_grouping_variable = "group",
    ruv_number_k = 1
  ))
  ruv_corrected <- ruviii_method(ruv_obj)
  expect_s4_class(ruv_corrected, "ProteinQuantitativeData")
  expect_identical(colnames(ruv_corrected@protein_quant_table), c("Protein.Ids", "S1", "S2", "S3"))
  expect_identical(ruv_corrected@protein_quant_table$Protein.Ids, c("P1", "P2"))
  expect_equal(ruv_corrected@protein_quant_table$S1, c(10, 20))
})

test_that("KIT-101 F-23: a between-sample normalisation on a DIA_LIMPA object reaches the quantified EList dpcDE fits; none leaves it byte-identical", {
  skip_if_not_installed("limpa")
  skip_if_not_installed("limma")
  ids <- sprintf("P%02d", 1:60)
  samples <- paste0("S", 1:6)
  set.seed(23)
  E <- matrix(stats::rnorm(360, 20, 1.5), 60, 6, dimnames = list(ids, samples))
  E <- sweep(E, 2, c(0, 0.3, -0.2, 0.5, -0.4, 0.1), "+")
  nobs <- matrix(2L, 60, 6, dimnames = dimnames(E))
  nobs[c(5, 70, 200)] <- 0L
  se <- matrix(0.2, 60, 6, dimnames = dimnames(E))
  elist <- methods::new("EList", list(E = E, genes = data.frame(NPrec = rep(2L, 60L), PropObs = rep(1, 60L), row.names = ids),
    other = list(standard.error = se, n.observations = nobs)))
  table <- data.frame(Protein.Ids = ids, as.data.frame(E, check.names = FALSE), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(table) <- NULL
  design <- data.frame(Run = samples, group = rep(c("A", "B"), each = 3L), replicates = rep(1:3, 2L), stringsAsFactors = FALSE)
  object <- ProteinQuantitativeData(protein_quant_table = table, protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = ids, stringsAsFactors = FALSE), design_matrix = design,
    sample_id = "Run", group_id = "group", technical_replicate_id = "replicates",
    args = list(limpa_dpc_quant_results = list(quantified_elist = elist, standard_errors = se, n_observations = nobs)))
  normalise <- function(method) {
    out <- NULL
    invisible(utils::capture.output(out <- suppressMessages(normaliseBetweenSamples(object, normalisation_method = method))))
    out
  }

  # The default (none): the DPC results, the EList and its E are byte-identical to the input.
  none <- normalise("none")
  expect_identical(serialize(none@args$limpa_dpc_quant_results$quantified_elist$E, NULL), serialize(E, NULL))
  expect_identical(none@args$limpa_dpc_quant_results, object@args$limpa_dpc_quant_results)

  # A chosen method: the EList's E is normalised as the protein table is, and nothing else in the results changes.
  cyclic <- normalise("cyclicloess")
  results <- cyclic@args$limpa_dpc_quant_results
  normalised_E <- results$quantified_elist$E
  expect_false(isTRUE(all.equal(normalised_E, E)))
  expect_equal(normalised_E, limma::normalizeCyclicLoess(E))
  table_after <- as.matrix(cyclic@protein_quant_table[match(ids, cyclic@protein_quant_table$Protein.Ids), samples])
  expect_equal(unname(normalised_E[ids, samples]), unname(table_after))
  expect_identical(results$quantified_elist$other, elist$other)
  expect_identical(results$quantified_elist$genes, elist$genes)
  expect_identical(results$standard_errors, se)
  expect_identical(results$n_observations, nobs)

  # The DA's dpcDE branch receives the normalised E.
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(dpcDE = function(y, design, plot = TRUE, ...) {
    captured$E <- y$E
    limma::lmFit(y$E, design)
  }, .package = "limpa")
  result <- suppressMessages(differentialAbundanceAnalysis(cyclic,
    contrasts_tbl = data.frame(contrasts = "B_vs_A=groupB-groupA", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = FALSE, eBayes_robust = FALSE,
    args_row_id = "Protein.Ids", qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue"))
  expect_identical(result$da_model$model, "dpcDE")
  expect_equal(captured$E, normalised_E[, samples])
})
