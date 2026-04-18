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

if (!methods::isGeneric("cleanDesignMatrix")) {
  methods::setGeneric(
    "cleanDesignMatrix",
    function(theObject) {
      standardGeneric("cleanDesignMatrix")
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
    isTargetSetMethod(expr, "cleanDesignMatrix")
  },
  env = environment()
)

newCleanDesignMatrixObject <- function(assay_tbl, design_matrix) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(LCMS_Pos = assay_tbl),
    metabolite_id_column = "Name",
    design_matrix = design_matrix,
    sample_id = "Run",
    args = list()
  )
}

test_that("metabolomics S4 cleanDesignMatrix preserves assay sample ordering while dropping unmatched design rows", {
  object <- newCleanDesignMatrixObject(
    assay_tbl = tibble::tibble(
      Name = c("M1", "M2"),
      annotation = c("alpha", "beta"),
      Sample_2 = c(20, 40),
      Sample_1 = c(10, 30)
    ),
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_3", "Sample_2"),
      Batch = c("A", "B", "C"),
      stringsAsFactors = FALSE
    )
  )

  cleaned <- cleanDesignMatrix(object)

  expect_identical(cleaned@design_matrix$Run, c("Sample_2", "Sample_1"))
  expect_identical(cleaned@design_matrix$Batch, c("C", "A"))
})

test_that("metabolomics S4 cleanDesignMatrix leaves the object unchanged when no assay sample ids match", {
  object <- newCleanDesignMatrixObject(
    assay_tbl = tibble::tibble(
      Name = c("M1", "M2"),
      Sample_1 = c(10, 20),
      Sample_2 = c(30, 40)
    ),
    design_matrix = data.frame(
      Run = c("Missing_1", "Missing_2"),
      Batch = c("A", "B"),
      stringsAsFactors = FALSE
    )
  )

  cleaned <- NULL
  expect_warning(
    cleaned <- cleanDesignMatrix(object),
    "No sample columns identified in the first assay matching design matrix sample IDs.",
    fixed = TRUE
  )

  expect_identical(cleaned@design_matrix, object@design_matrix)
})

test_that("metabolomics S4 cleanDesignMatrix source retains join-and-rename cleanup steps", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "cleanDesignMatrix")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "dplyr::inner_join\\(")
  expect_match(target_text, "temp_sample_id")
  expect_match(target_text, "dplyr::rename\\(")
  expect_match(target_text, "dplyr::filter\\(")
})
