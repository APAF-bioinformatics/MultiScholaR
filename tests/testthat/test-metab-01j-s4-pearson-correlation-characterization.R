library(methods)
library(testthat)
library(dplyr)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")
      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedExpressions <- function(paths, matcher, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        eval(expr, envir = env)
      }
    }
  }
}

findSelectedExpression <- function(paths, matcher) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      if (matcher(expr)) {
        return(expr)
      }
    }
  }

  NULL
}

normalizeSelectorValue <- function(x) {
  if (is.character(x)) {
    return(x[[1]])
  }

  if (is.symbol(x)) {
    return(as.character(x))
  }

  as.character(x)
}

isTargetSetMethod <- function(expr, method_name) {
  if (!is.call(expr) || !identical(as.character(expr[[1]]), "setMethod")) {
    return(FALSE)
  }

  expr_parts <- as.list(expr)
  part_names <- names(expr_parts)
  method_arg <- expr_parts[[which(part_names == "f")[1]]]

  if (is.null(method_arg) && length(expr_parts) >= 2) {
    method_arg <- expr_parts[[2]]
  }

  !is.null(method_arg) && identical(normalizeSelectorValue(method_arg), method_name)
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      design_matrix = "data.frame",
      sample_id = "character",
      metabolite_id_column = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      design_matrix = data.frame(),
      sample_id = "Run",
      metabolite_id_column = "feature",
      technical_replicate_id = "TechRep",
      args = list()
    )
  )
}

if (!methods::isGeneric("pearsonCorForSamplePairs")) {
  methods::setGeneric(
    "pearsonCorForSamplePairs",
    function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA) {
      standardGeneric("pearsonCorForSamplePairs")
    }
  )
}

checkParamsObjectFunctionSimplifyAcceptNull <- function(theObject, param_name_string, default_value = NULL) {
  if (!is.null(default_value)) {
    return(default_value)
  }

  function_args <- theObject@args$pearsonCorForSamplePairs
  if (!is.null(function_args) && !is.null(function_args[[param_name_string]])) {
    return(function_args[[param_name_string]])
  }

  NULL
}

helper_paths <- c(
  file.path(repo_root, "R", "func_prot_qc_correlation_helpers.R"),
  file.path(repo_root, "R", "func_metab_qc_correlation_helpers.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedFunctions(
  paths = helper_paths,
  symbols = c("getPairsOfSamplesTable", "calculateMetabolitePairCorrelation"),
  env = environment()
)

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "pearsonCorForSamplePairs")
  },
  env = environment()
)

newPearsonObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        feature = c("M1", "M2", "M3"),
        Sample_1 = c(1, 2, 3),
        Sample_2 = c(1, 2, 3),
        Sample_3 = c(3, 2, 1),
        Sample_4 = c(1, 2, 3)
      )
    ),
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      TechRep = c("pair_a", "pair_a", "pool_qc", "pool_qc"),
      AltGroup = c("grp1", "grp1", "grp2", "grp2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    metabolite_id_column = "feature",
    technical_replicate_id = "TechRep",
    args = list()
  )
}

test_that("metabolomics S4 pearson-correlation method keeps non-pool pairs and preserves assay names", {
  results <- suppressWarnings(
    suppressMessages(
      pearsonCorForSamplePairs(
        newPearsonObject(),
        tech_rep_remove_regex = "pool"
      )
    )
  )

  expect_named(results, "LCMS_Pos")

  assay_results <- results$LCMS_Pos

  expect_s3_class(assay_results, "data.frame")
  expect_identical(assay_results$TechRep, "pair_a")
  expect_identical(assay_results$Run.x, "Sample_2")
  expect_identical(assay_results$Run.y, "Sample_1")
  expect_equal(assay_results$pearson_correlation, 1)
})

test_that("metabolomics S4 pearson-correlation method honors correlation-group override", {
  results <- suppressWarnings(
    suppressMessages(
      pearsonCorForSamplePairs(
        newPearsonObject(),
        correlation_group = "AltGroup"
      )
    )
  )

  assay_results <- results$LCMS_Pos |>
    dplyr::arrange(AltGroup)

  expect_identical(assay_results$AltGroup, c("grp1", "grp2"))
  expect_equal(assay_results$pearson_correlation, c(1, -1))
})

test_that("metabolomics S4 pearson-correlation source retains pair-building and helper-routing logic", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "pearsonCorForSamplePairs")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "tidyr::pivot_longer\\(")
  expect_match(target_text, "getPairsOfSamplesTable\\(")
  expect_match(target_text, "calculateMetabolitePairCorrelation\\(")
  expect_match(target_text, "!stringr::str_detect\\(")
})
