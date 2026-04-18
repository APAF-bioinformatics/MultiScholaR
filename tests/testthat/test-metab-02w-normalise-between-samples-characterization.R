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
      sample_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      args = list()
    )
  )
}

if (!methods::isGeneric("normaliseBetweenSamples")) {
  methods::setGeneric(
    "normaliseBetweenSamples",
    function(theObject, normalisation_method = NULL) {
      standardGeneric("normaliseBetweenSamples")
    }
  )
}

if (!methods::isGeneric("cleanDesignMatrix")) {
  methods::setGeneric(
    "cleanDesignMatrix",
    function(theObject) {
      standardGeneric("cleanDesignMatrix")
    }
  )
}

if (!methods::hasMethod("cleanDesignMatrix", "MetaboliteAssayData")) {
  methods::setMethod(
    "cleanDesignMatrix",
    "MetaboliteAssayData",
    function(theObject) {
      assay_names <- names(theObject@metabolite_data)
      if (length(theObject@metabolite_data) >= 1) {
        assay_cols <- colnames(theObject@metabolite_data[[1]])
        keep_samples <- intersect(
          as.character(theObject@design_matrix[[theObject@sample_id]]),
          assay_cols
        )
        theObject@design_matrix <- theObject@design_matrix[
          theObject@design_matrix[[theObject@sample_id]] %in% keep_samples,
          ,
          drop = FALSE
        ]
      }

      names(theObject@metabolite_data) <- assay_names
      theObject@args$clean_design_matrix_called <- TRUE
      theObject
    }
  )
}

checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default = NULL) {
  param_value <- get(param_name_string, envir = parent.frame())
  if (!is.null(param_value)) {
    return(param_value)
  }

  object_value <- theObject@args[[param_name_string]]
  if (!is.null(object_value)) {
    return(object_value)
  }

  default
}

log_info <- function(...) {
  invisible(NULL)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "normaliseBetweenSamples")
  },
  env = environment()
)

newBetweenSampleNormObject <- function(assay_tbl, assay_names = "", args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = stats::setNames(list(assay_tbl), assay_names),
    metabolite_id_column = "Name",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3"),
      Batch = c("A", "A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    args = args
  )
}

test_that("metabolomics S4 between-sample normalization preserves assay reconstruction and design cleanup for method none", {
  assay_tbl <- tibble::tibble(
    Name = c("M1", "M2"),
    annotation = c("alpha", "beta"),
    Batch = c("Batch_1", "Batch_1"),
    Sample_1 = c(10, 20),
    Sample_2 = c("30", "40")
  )

  normalized <- suppressMessages(
    suppressWarnings(
      normaliseBetweenSamples(
        newBetweenSampleNormObject(
          assay_tbl = assay_tbl,
          assay_names = "",
          args = list(normalisation_method = "none")
        )
      )
    )
  )

  normalized_assay <- normalized@metabolite_data$Assay_1

  expect_named(normalized@metabolite_data, "Assay_1")
  expect_s3_class(normalized_assay, "tbl_df")
  expect_identical(
    colnames(normalized_assay),
    c("Name", "annotation", "Batch", "Sample_1", "Sample_2")
  )
  expect_equal(normalized_assay$Sample_1, c(10, 20))
  expect_equal(normalized_assay$Sample_2, c(30, 40))
  expect_identical(normalized@args$normalisation_method, "none")
  expect_true(isTRUE(normalized@args$clean_design_matrix_called))
  expect_identical(normalized@design_matrix$Run, c("Sample_1", "Sample_2"))
})

test_that("metabolomics S4 between-sample normalization preserves all-NA sample fallback", {
  assay_tbl <- tibble::tibble(
    Name = c("M1", "M2"),
    annotation = c("alpha", "beta"),
    Sample_1 = c(10, 20),
    Sample_2 = c(NA_real_, NA_real_)
  )

  original <- newBetweenSampleNormObject(
    assay_tbl = assay_tbl,
    assay_names = "LCMS_Pos",
    args = list(normalisation_method = "quantile")
  )

  normalized <- NULL
  suppressMessages(
    expect_warning(
      normalized <- normaliseBetweenSamples(original),
      "At least one sample column contains only NA values. Skipping normalization.",
      fixed = TRUE
    )
  )

  expect_equal(normalized@metabolite_data$LCMS_Pos, assay_tbl)
  expect_identical(normalized@args$normalisation_method, "quantile")
  expect_true(isTRUE(normalized@args$clean_design_matrix_called))
})

test_that("metabolomics S4 between-sample normalization source retains method dispatch and design cleanup", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "normaliseBetweenSamples")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(
    target_text,
    "theObject@args\\$normalisation_method <- normalisation_method_final"
  )
  expect_match(target_text, "limma::normalizeQuantiles\\(assay_matrix\\)")
  expect_match(target_text, "limma::normalizeMedianAbsValues\\(assay_matrix\\)")
  expect_match(target_text, "cleanDesignMatrix\\(theObject\\)")
})
