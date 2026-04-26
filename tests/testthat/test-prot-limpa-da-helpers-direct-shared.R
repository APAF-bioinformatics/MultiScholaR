# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("limpa DA conversion helper preserves error and standard formatting branches", {
  converter_missing_limma <- makeFunctionWithOverrides(
    convertDpcDAToStandardFormat,
    list(
      requireNamespace = function(package, quietly = TRUE) FALSE
    )
  )

  expect_error(
    converter_missing_limma(
      dpc_fit = structure(list(), class = "MArrayLM"),
      contrast_strings = "B-A",
      design_matrix = data.frame(A = 1, B = 1)
    ),
    "limma package is required",
    fixed = TRUE
  )

  skip_if_not_installed("limma")
  skip_if_not_installed("qvalue")

  exprs <- matrix(
    c(
      10, 11, 15, 16,
      12, 13, 18, 19,
      20, 21, 25, 26
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("P1", "P2", "P3"), c("S1", "S2", "S3", "S4"))
  )
  design <- stats::model.matrix(~ 0 + group, data = data.frame(group = c("A", "A", "B", "B")))
  colnames(design) <- c("A", "B")
  fit <- limma::lmFit(exprs, design)

  converted <- convertDpcDAToStandardFormat(
    dpc_fit = fit,
    contrast_strings = c(
      "B_vs_A=B-A",
      second = "A-B"
    ),
    design_matrix = design,
    eBayes_trend = FALSE,
    eBayes_robust = FALSE
  )

  expect_true(is.list(converted))
  expect_true(isTRUE(converted$dpc_method_used))
  expect_identical(names(converted$results), c("B_vs_A=B-A", "A-B"))
  expect_true(all(c(
    "uniprot_acc",
    "comparison",
    "logFC",
    "log_intensity",
    "raw_pvalue",
    "fdr_qvalue",
    "fdr_value_bh_adjustment"
  ) %in% colnames(converted$results[[1]])))
  expect_identical(unique(converted$results[[1]]$comparison), "B_vs_A")
  expect_identical(unique(converted$results[[2]]$comparison), "second")
})
