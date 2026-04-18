library(testthat)
library(dplyr)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) && length(expr) >= 3 && as.character(expr[[1]]) %in% c("<-", "=")
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

if (!methods::isGeneric("metaboliteIntensityFiltering")) {
  methods::setGeneric(
    "metaboliteIntensityFiltering",
    function(
      theObject,
      metabolites_intensity_cutoff_percentile = NULL,
      metabolites_proportion_of_samples_below_cutoff = NULL
    ) {
      standardGeneric("metaboliteIntensityFiltering")
    }
  )
}

checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default_value = NULL) {
  param_value <- get(param_name_string, envir = parent.frame())

  object_value <- NULL
  if (!is.null(theObject@args) &&
      !is.null(theObject@args$metaboliteIntensityFiltering) &&
      is.list(theObject@args$metaboliteIntensityFiltering)) {
    object_value <- theObject@args$metaboliteIntensityFiltering[[param_name_string]]
  }

  if (!is.null(param_value)) {
    return(param_value)
  }

  if (!is.null(object_value)) {
    return(object_value)
  }

  default_value
}

updateParamInObject <- function(theObject, param_name_string) {
  if (is.null(theObject@args$metaboliteIntensityFiltering)) {
    theObject@args$metaboliteIntensityFiltering <- list()
  }

  theObject@args$metaboliteIntensityFiltering[[param_name_string]] <- get(param_name_string, envir = parent.frame())
  theObject
}

helper_paths <- c(
  file.path(repo_root, "R", "func_metab_qc_filtering_helpers.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedFunctions(
  paths = helper_paths,
  symbols = "metaboliteIntensityFilteringHelper",
  env = environment()
)

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_qc_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "metaboliteIntensityFiltering")
  },
  env = environment()
)

newMetabFilteringObject <- function(args = list()) {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = tibble::tibble(
        Name = c("M1", "M2", "M3"),
        annotation = c("alpha", "beta", "gamma"),
        Sample_1 = c(10, 100, 60),
        Sample_2 = c(20, 200, NA)
      )
    ),
    metabolite_id_column = "Name",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    args = args
  )
}

test_that("metabolomics S4 metabolite-intensity filtering preserves config cleanup and assay filtering", {
  filtered <- metaboliteIntensityFiltering(
    newMetabFilteringObject(
      args = list(
        metaboliteIntensityFiltering = list(
          metabolites_intensity_cutoff_percentile = "50 # config percentile",
          metabolites_proportion_of_samples_below_cutoff = "0.6 # config cutoff"
        )
      )
    )
  )

  filtered_assay <- filtered@metabolite_data$LCMS_Pos

  expect_s3_class(filtered_assay, "tbl_df")
  expect_identical(filtered_assay$Name, c("M2", "M3"))
  expect_identical(filtered_assay$annotation, c("beta", "gamma"))
  expect_equal(filtered_assay$Sample_1, c(100, 60))
  expect_equal(filtered_assay$Sample_2, c(200, NA))
  expect_null(
    filtered@args$metaboliteIntensityFiltering$metabolites_intensity_cutoff_percentile
  )
  expect_null(
    filtered@args$metaboliteIntensityFiltering$metabolites_proportion_of_samples_below_cutoff
  )
})

test_that("metabolomics S4 metabolite-intensity filtering rejects non-numeric cutoff values after comment stripping", {
  expect_error(
    suppressWarnings(
      metaboliteIntensityFiltering(
        newMetabFilteringObject(
          args = list(
            metaboliteIntensityFiltering = list(
              metabolites_intensity_cutoff_percentile = "bad_value # config percentile",
              metabolites_proportion_of_samples_below_cutoff = "0.6"
            )
          )
        )
      )
    ),
    "Failed to convert cleaned metabolites_intensity_cutoff_percentile ('bad_value' from raw 'bad_value # config percentile') to numeric. Check config.ini or parameter value.",
    fixed = TRUE
  )
})

test_that("metabolomics S4 metabolite-intensity filtering source retains comment cleanup and helper routing", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "metaboliteIntensityFiltering")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(
    target_text,
    "cleaned_intensity_percentile <- trimws\\(sub\\(\"#\\.\\*\\$\", \"\","
  )
  expect_match(
    target_text,
    "ceiling\\(quantile\\(all_intensity_values,"
  )
  expect_match(target_text, "metaboliteIntensityFilteringHelper\\(")
  expect_match(
    target_text,
    "updateParamInObject\\(theObject, config_proportion_cutoff\\)"
  )
})
