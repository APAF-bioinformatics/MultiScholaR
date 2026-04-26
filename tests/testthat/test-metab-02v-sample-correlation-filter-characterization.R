library(methods)
library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

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
      metabolite_id_column = "character",
      design_matrix = "data.frame",
      sample_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "feature",
      design_matrix = data.frame(),
      sample_id = "Run"
    )
  )
}

if (!methods::isGeneric("filterSamplesByMetaboliteCorrelationThreshold")) {
  methods::setGeneric(
    "filterSamplesByMetaboliteCorrelationThreshold",
    function(
      theObject,
      pearson_correlation_per_pair = NULL,
      min_pearson_correlation_threshold = 0.5
    ) {
      standardGeneric("filterSamplesByMetaboliteCorrelationThreshold")
    }
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "filterSamplesByMetaboliteCorrelationThreshold")
  },
  env = environment()
)

newCorrelationFilterObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = data.frame(
        feature = c("M1", "M2"),
        Sample_1 = c(10, 50),
        Sample_2 = c(11, 51),
        Sample_3 = c(30, 70),
        Sample_4 = c(31, 71),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "feature",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      Batch = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run"
  )
}

newCorrelationResults <- function() {
  list(
    LCMS_Pos = data.frame(
      Run.x = c("Sample_1", "Sample_3"),
      Run.y = c("Sample_2", "Sample_4"),
      pearson_correlation = c(0.95, 0.42),
      stringsAsFactors = FALSE
    )
  )
}

test_that("metabolomics S4 sample-correlation filtering removes low-correlation samples from assay and design matrix", {
  filtered <- suppressMessages(
    filterSamplesByMetaboliteCorrelationThreshold(
      newCorrelationFilterObject(),
      pearson_correlation_per_pair = newCorrelationResults(),
      min_pearson_correlation_threshold = 0.8
    )
  )

  filtered_assay <- filtered@metabolite_data$LCMS_Pos

  expect_named(filtered@metabolite_data, "LCMS_Pos")
  expect_identical(
    colnames(filtered_assay),
    c("feature", "Sample_1", "Sample_2")
  )
  expect_equal(filtered_assay$Sample_1, c(10, 50))
  expect_equal(filtered_assay$Sample_2, c(11, 51))
  expect_identical(filtered@design_matrix$Run, c("Sample_1", "Sample_2"))
})

test_that("metabolomics S4 sample-correlation filtering validates correlation inputs", {
  expect_error(
    suppressMessages(
      filterSamplesByMetaboliteCorrelationThreshold(
        newCorrelationFilterObject(),
        pearson_correlation_per_pair = 1,
        min_pearson_correlation_threshold = 0.8
      )
    ),
    "`pearson_correlation_per_pair` must be a list of correlation data frames \\(one per assay\\)\\."
  )

  expect_error(
    suppressMessages(
      filterSamplesByMetaboliteCorrelationThreshold(
        newCorrelationFilterObject(),
        pearson_correlation_per_pair = newCorrelationResults(),
        min_pearson_correlation_threshold = "0.8"
      )
    ),
    "`min_pearson_correlation_threshold` must be a numeric value\\."
  )
})

test_that("metabolomics S4 sample-correlation filtering source retains pivoted pair filtering and design pruning", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "filterSamplesByMetaboliteCorrelationThreshold")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "tidyr::pivot_longer\\(")
  expect_match(
    target_text,
    "pearson_correlation >= min_pearson_correlation_threshold"
  )
  expect_match(target_text, "samples_to_remove_total <<- c\\(")
  expect_match(
    target_text,
    "samples_to_remove_unique <- unique\\(samples_to_remove_total\\)"
  )
  expect_match(
    target_text,
    "!!rlang::sym\\(sample_id_col_name\\) %in% samples_to_remove_unique"
  )
})
