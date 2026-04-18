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

  !is.null(method_arg) && identical(normalizeSelectorValue(method_arg), method_name)
}

if (!methods::isClass("MetaboliteAssayData")) {
  methods::setClass(
    "MetaboliteAssayData",
    slots = c(
      metabolite_data = "list",
      metabolite_id_column = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      internal_standard_regex = "character",
      annotation_id_column = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      internal_standard_regex = "Internal Standard",
      annotation_id_column = "annotation",
      args = list()
    )
  )
}

if (!methods::isGeneric("normaliseUntransformedData")) {
  methods::setGeneric(
    "normaliseUntransformedData",
    function(
      theObject,
      method = "ITSD",
      itsd_feature_ids = NULL,
      itsd_pattern_columns = NULL,
      itsd_aggregation = "sum",
      remove_itsd_after_norm = TRUE,
      ...
    ) {
      standardGeneric("normaliseUntransformedData")
    }
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "normaliseUntransformedData")
  },
  env = environment()
)

newMetabItsdObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        Name = c("ITSD_1", "M1", "M2"),
        annotation = c("Internal Standard", "alpha", "beta"),
        Sample_1 = c(10, 5, 15),
        Sample_2 = c(20, 30, 60),
        qc_score = c(100, 200, 300)
      )
    ),
    metabolite_id_column = "Name",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    internal_standard_regex = "Internal Standard",
    annotation_id_column = "annotation",
    args = list()
  )
}

test_that("metabolomics S4 ITSD normalization preserves the current average-centered contract", {
  normalized <- normaliseUntransformedData(newMetabItsdObject(), method = "ITSD")
  normalized_assay <- normalized@metabolite_data$LCMS_Pos

  expect_s3_class(normalized_assay, "tbl_df")
  expect_identical(normalized_assay$Name, c("M1", "M2"))
  expect_equal(normalized_assay$Sample_1, c(7.5, 22.5))
  expect_equal(normalized_assay$Sample_2, c(22.5, 45))
  expect_identical(normalized_assay$qc_score, c(200, 300))
  expect_true(normalized@args$ITSDNormalization$applied)
  expect_identical(normalized@args$ITSDNormalization$method_type, "average_centered")
  expect_identical(normalized@args$ITSDNormalization$itsd_aggregation, "sum")
  expect_identical(normalized@args$ITSDNormalization$itsd_pattern_columns, "annotation")
  expect_identical(normalized@args$ITSDNormalization$removed_itsd, TRUE)
  expect_identical(normalized@args$ITSDNormalization$itsd_features_per_assay$LCMS_Pos, "ITSD_1")
  expect_identical(normalized@args$ITSDNormalization$itsd_counts_per_assay$LCMS_Pos, 1L)
  expect_s3_class(normalized@args$ITSDNormalization$timestamp, "POSIXct")
})

test_that("metabolomics S4 ITSD normalization rejects unsupported methods", {
  expect_error(
    normaliseUntransformedData(newMetabItsdObject(), method = "PQN"),
    "This method currently only supports method = 'ITSD'.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 ITSD normalization source retains default-pattern fallback and average-centered scaling", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "normaliseUntransformedData")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "itsd_pattern_columns <- annotation_col")
  expect_match(
    target_text,
    "\\.x \\* \\(average_norm_factor/sample_norm_factors_vec\\[dplyr::cur_column\\(\\)\\]\\)"
  )
  expect_match(
    target_text,
    "theObject@args\\$ITSDNormalization\\$itsd_features_per_assay <- as.list\\(itsd_collector\\$features_per_assay\\)"
  )
})
