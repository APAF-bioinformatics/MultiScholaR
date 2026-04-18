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

isTargetSetClass <- function(expr, class_name) {
  is.call(expr) &&
    identical(as.character(expr[[1]]), "setClass") &&
    length(expr) >= 2 &&
    identical(normalizeSelectorValue(expr[[2]]), class_name)
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
      sample_id = "character"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id",
      design_matrix = data.frame(),
      sample_id = "Run"
    )
  )
}

if (!methods::isGeneric("plotInteractiveVolcano")) {
  methods::setGeneric(
    "plotInteractiveVolcano",
    function(objectsList, anno_list = NULL) {
      standardGeneric("plotInteractiveVolcano")
    }
  )
}

left_join <- dplyr::left_join
join_by <- dplyr::join_by
sym <- rlang::sym
column_to_rownames <- tibble::column_to_rownames
decideTests <- function(...) {
  matrix(integer(), nrow = 0, ncol = 0)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_da_results.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults") ||
      isTargetSetMethod(expr, "plotInteractiveVolcano")
  },
  env = environment()
)

newDaInteractiveResultObject <- function() {
  counts_table <- data.frame(
    metabolite_id = c("M1", "M2"),
    S1 = c(10, 20),
    S2 = c(30, 40),
    stringsAsFactors = FALSE
  )

  fit_eb <- list(
    coefficients = matrix(numeric(), nrow = 2, ncol = 0),
    p.value = matrix(numeric(), nrow = 2, ncol = 0),
    t = matrix(numeric(), nrow = 2, ncol = 0),
    stdev.unscaled = matrix(numeric(), nrow = 2, ncol = 0),
    genes = data.frame(feature = c("alpha", "beta"), stringsAsFactors = FALSE)
  )

  methods::new(
    "MetabolomicsDifferentialAbundanceResults",
    theObject = methods::new(
      "MetaboliteAssayData",
      metabolite_data = list(primary = counts_table),
      metabolite_id_column = "metabolite_id",
      design_matrix = data.frame(
        Run = c("S1", "S2"),
        genotype_group = c("A", "B"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run"
    ),
    fit.eb = fit_eb,
    contrasts_results_table = list()
  )
}

test_that("metabolomics S4 interactive-volcano method preserves zero-coefficient fallback", {
  output <- NULL

  expect_warning(
    output <- plotInteractiveVolcano(list(primary = newDaInteractiveResultObject())),
    "components were missing rownames"
  )

  expect_identical(names(output), "primary")
  expect_s4_class(output$primary, "MetabolomicsDifferentialAbundanceResults")
  expect_identical(output$primary@interactive_volcano_plot, list())
})

test_that("metabolomics S4 interactive-volcano definition retains qvalue and Glimma hooks", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "plotInteractiveVolcano")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "qvalue::qvalue", fixed = TRUE)
  expect_match(target_text, "Glimma::glimmaVolcano", fixed = TRUE)
  expect_match(target_text, "display.columns = if \\(!is.null\\(anno_tbl\\)\\)")
  expect_match(target_text, "colnames\\(anno_tbl\\)")
})
