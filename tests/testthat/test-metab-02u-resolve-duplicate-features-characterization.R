library(testthat)
library(dplyr)

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
      annotation_id_column = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "database_identifier",
      annotation_id_column = "annotation_id",
      internal_standard_regex = "^ISTD",
      design_matrix = data.frame(),
      sample_id = "Run"
    )
  )
}

if (!methods::isGeneric("resolveDuplicateFeatures")) {
  methods::setGeneric(
    "resolveDuplicateFeatures",
    function(theObject, itsd_pattern_columns = NULL) {
      standardGeneric("resolveDuplicateFeatures")
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
    isTargetSetMethod(expr, "resolveDuplicateFeatures")
  },
  env = environment()
)

newResolveDuplicateFeaturesObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        database_identifier = c("ISTD_LOW", "ISTD_HIGH", "DB_DUP", "DB_DUP", "DB_KEEP"),
        metabolite = c("ISTD Mix", "ISTD Mix", "Non-ITSD low", "Non-ITSD high", "Unique"),
        annotation_id = c("ISTD_alpha", "ISTD_alpha", "feature_low", "feature_high", "feature_unique"),
        Sample_1 = c(10, 30, 5, 50, 1),
        Sample_2 = c(20, 40, 15, 60, 2)
      )
    ),
    metabolite_id_column = "database_identifier",
    annotation_id_column = "annotation_id",
    internal_standard_regex = "^ISTD",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run"
  )
}

test_that("metabolomics S4 duplicate resolution preserves ITSD renaming and non-ITSD deduplication", {
  resolved <- suppressWarnings(
    suppressMessages(resolveDuplicateFeatures(newResolveDuplicateFeaturesObject()))
  )
  resolved_assay <- resolved@metabolite_data$LCMS_Pos

  expect_named(resolved@metabolite_data, "LCMS_Pos")
  expect_equal(nrow(resolved_assay), 4)
  expect_false(".original_row_id" %in% names(resolved_assay))

  itsd_rows <- resolved_assay %>%
    filter(grepl("^ISTD Mix_[12]$", database_identifier)) %>%
    arrange(database_identifier)

  expect_equal(itsd_rows$database_identifier, c("ISTD Mix_1", "ISTD Mix_2"))
  expect_equal(itsd_rows$Sample_1, c(30, 10))
  expect_equal(itsd_rows$Sample_2, c(40, 20))

  non_itsd_rows <- resolved_assay %>%
    filter(database_identifier %in% c("DB_DUP", "DB_KEEP")) %>%
    arrange(database_identifier)

  expect_equal(non_itsd_rows$database_identifier, c("DB_DUP", "DB_KEEP"))
  expect_equal(non_itsd_rows$metabolite, c("Non-ITSD high", "Unique"))
  expect_equal(non_itsd_rows$Sample_1, c(50, 1))
  expect_equal(non_itsd_rows$Sample_2, c(60, 2))
})

test_that("metabolomics S4 duplicate-resolution source retains ITSD separation and ranked deduplication", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "resolveDuplicateFeatures")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "columns_to_check_for_itsd <- annot_col_slot_name")
  expect_match(target_text, "\\.original_row_id = dplyr::row_number\\(\\)")
  expect_match(target_text, "stringr::str_detect\\(values, itsd_regex\\)")
  expect_match(
    target_text,
    "\\.new_itsd_id = paste0\\(!!rlang::sym\\(specific_name_col\\)"
  )
  expect_match(
    target_text,
    "dplyr::slice_max\\(order_by = \\.data\\$avg_intensity,"
  )
  expect_match(target_text, "dplyr::select\\(-\\.original_row_id\\)")
})
