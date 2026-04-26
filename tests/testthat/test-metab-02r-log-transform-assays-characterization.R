library(testthat)

# fidelity-coverage-compare: shared

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

if (!methods::isGeneric("logTransformAssays")) {
  methods::setGeneric(
    "logTransformAssays",
    function(theObject, offset = 1, ...) {
      standardGeneric("logTransformAssays")
    }
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

if (!methods::hasMethod("logTransformAssays", "MetaboliteAssayData")) {
  loadSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "logTransformAssays")
    },
    env = environment()
  )
}

newMetabLogTransformObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2"),
        annotation = c("alpha", "beta"),
        Sample_1 = c(0, -3),
        Sample_2 = c(6, 14),
        qc_score = c(10, 20)
      )
    ),
    metabolite_id_column = "Name",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    args = list()
  )
}

test_that("metabolomics S4 log-transform method preserves current assay transformation contract", {
  transformed <- logTransformAssays(newMetabLogTransformObject(), offset = 2)
  transformed_assay <- transformed@metabolite_data$LCMS_Pos

  expect_s3_class(transformed_assay, "tbl_df")
  expect_equal(transformed_assay$Sample_1, c(1, 1))
  expect_equal(transformed_assay$Sample_2, log2(c(8, 16)))
  expect_identical(transformed_assay$Name, c("M1", "M2"))
  expect_identical(transformed_assay$annotation, c("alpha", "beta"))
  expect_identical(transformed_assay$qc_score, c(10, 20))
  expect_true(transformed@args$log_transformed)
  expect_identical(transformed@args$log_transform_offset, 2)
})

test_that("metabolomics S4 log-transform method rejects non-positive offsets", {
  expect_error(
    logTransformAssays(newMetabLogTransformObject(), offset = 0),
    "`offset` must be a single positive numeric value.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 log-transform source retains the current negative-value clamp and args updates", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "logTransformAssays")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(
    target_text,
    "ifelse\\(!is\\.na\\(.x\\)\\s*&\\s*\\.x\\s*<\\s*0,\\s*0,\\s*\\.x\\)"
  )
  expect_match(target_text, "theObject@args\\$log_transformed <- TRUE")
  expect_match(target_text, "theObject@args\\$log_transform_offset <- offset")
})
