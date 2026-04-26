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
    impute.nipals = function(x, ncomp) {
      x[is.na(x)] <- 99
      x
    },
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
        dimnames = list(colnames(Y), rownames(Y))
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
