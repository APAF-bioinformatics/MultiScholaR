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
    function(theObject, offset = 0, ...) {
      standardGeneric("logTransformAssays")
    }
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_transformation_methods.R"),
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

test_that("metabolomics S4 log-transform method accepts offset 0 and rejects negative offsets", {
  # KIT-172: offset 0 (plain log2) is the default; below 0 is refused with its class.
  expect_s4_class(suppressWarnings(suppressMessages(logTransformAssays(newMetabLogTransformObject(), offset = 0))),
    "MetaboliteAssayData")
  expect_error(
    logTransformAssays(newMetabLogTransformObject(), offset = -1),
    class = "omics_log2_offset_invalid"
  )
})

test_that("metabolomics S4 log-transform method keeps the negative-value clamp above offset 0 and reads non-positive values as missing at 0", {
  # KIT-172: above offset 0 the clamp of earlier releases (a negative value becomes 0, then log2(x + offset)); at
  # offset 0 a value of 0 or below becomes missing. The args record the offset either way.
  clamped <- logTransformAssays(newMetabLogTransformObject(), offset = 1)
  expect_equal(clamped@metabolite_data$LCMS_Pos$Sample_1, c(0, 0))
  expect_true(clamped@args$log_transformed)
  expect_identical(clamped@args$log_transform_offset, 1)
  plain <- suppressWarnings(logTransformAssays(newMetabLogTransformObject(), offset = 0))
  expect_true(all(is.na(plain@metabolite_data$LCMS_Pos$Sample_1)))
  expect_equal(plain@metabolite_data$LCMS_Pos$Sample_2, log2(c(6, 14)))
  expect_identical(plain@args$log_transform_offset, 0)
})
