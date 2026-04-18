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
      group_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      args = list()
    )
  )
}

if (!methods::isClass("MetabolomicsDifferentialAbundanceResults")) {
  methods::setClass(
    "MetabolomicsDifferentialAbundanceResults",
    slots = c(
      theObject = "MetaboliteAssayData",
      fit.eb = "list",
      contrasts_results_table = "list"
    )
  )
}

if (!methods::isGeneric("differentialAbundanceAnalysisHelper")) {
  methods::setGeneric(
    "differentialAbundanceAnalysisHelper",
    function(
      theObject,
      contrasts_tbl = NULL,
      formula_string = NULL,
      group_id = NULL,
      da_q_val_thresh = NULL,
      treat_lfc_cutoff = NULL,
      eBayes_trend = NULL,
      eBayes_robust = NULL,
      args_group_pattern = NULL
    ) {
      standardGeneric("differentialAbundanceAnalysisHelper")
    }
  )
}

one_of <- function(x) {
  rlang::sym(x)
}

rownames_to_column <- tibble::rownames_to_column

checkParamsObjectFunctionSimplify <- function(theObject, param_name, default) {
  current_value <- get(param_name, envir = parent.frame())

  if (!is.null(current_value)) {
    return(current_value)
  }

  object_value <- theObject@args[[param_name]]
  if (!is.null(object_value)) {
    return(object_value)
  }

  default
}

updateParamInObject <- function(theObject, param_name) {
  theObject@args[[param_name]] <- get(param_name, envir = parent.frame())
  theObject
}

run_tests_calls <- list()
runTestsContrasts <- function(
  data_matrix,
  contrast_strings,
  design_matrix,
  formula_string,
  treat_lfc_cutoff,
  eBayes_trend,
  eBayes_robust
) {
  run_tests_calls[[length(run_tests_calls) + 1]] <<- list(
    data_matrix = data_matrix,
    contrast_strings = contrast_strings,
    design_matrix = design_matrix,
    formula_string = formula_string,
    treat_lfc_cutoff = treat_lfc_cutoff,
    eBayes_trend = eBayes_trend,
    eBayes_robust = eBayes_robust
  )

  list(
    fit.eb = list(status = "ok"),
    results = list(
      primary = data.frame(
        comparison = contrast_strings,
        logFC = 1.5,
        stringsAsFactors = FALSE,
        row.names = "M1"
      )
    )
  )
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_da_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

loadSelectedExpressions(
  paths = target_paths,
  matcher = function(expr) {
    isTargetSetMethod(expr, "differentialAbundanceAnalysisHelper")
  },
  env = environment()
)

newMetabDaObject <- function() {
  methods::new(
    "MetaboliteAssayData",
    metabolite_data = list(
      LCMS_Pos = data.frame(
        Name = c("M1"),
        Sample_1 = c(10),
        Sample_2 = c(20),
        stringsAsFactors = FALSE
      )
    ),
    metabolite_id_column = "Name",
    design_matrix = data.frame(
      Run = c("Sample_1", "Sample_2"),
      group = c("1A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Run",
    group_id = "group",
    args = list()
  )
}

test_that("metabolomics S4 DA helper preserves numeric-group sanitization and result restoration", {
  run_tests_calls <<- list()
  target_object <- newMetabDaObject()

  output <- differentialAbundanceAnalysisHelper(
    target_object,
    contrasts_tbl = data.frame(contrasts = "1A-B", stringsAsFactors = FALSE),
    formula_string = "~ 0 + group",
    group_id = "group",
    treat_lfc_cutoff = 2,
    eBayes_trend = "FALSE",
    eBayes_robust = "YES"
  )

  expect_s4_class(output, "MetabolomicsDifferentialAbundanceResults")
  expect_length(run_tests_calls, 1)
  expect_identical(run_tests_calls[[1]]$contrast_strings, "grp_1A-B")
  expect_identical(run_tests_calls[[1]]$design_matrix$group, c("grp_1A", "B"))
  expect_identical(run_tests_calls[[1]]$formula_string, "~ 0 + group")
  expect_identical(run_tests_calls[[1]]$treat_lfc_cutoff, 2)
  expect_identical(run_tests_calls[[1]]$eBayes_trend, FALSE)
  expect_identical(run_tests_calls[[1]]$eBayes_robust, TRUE)
  expect_identical(rownames(run_tests_calls[[1]]$design_matrix), c("Sample_1", "Sample_2"))
  expect_equal(output@contrasts_results_table[[1]]$comparison, "1A-B")
  expect_equal(output@contrasts_results_table[[1]]$Name, "M1")
  expect_identical(output@theObject@args$treat_lfc_cutoff, 2)
  expect_identical(output@theObject@args$eBayes_trend, "FALSE")
  expect_identical(output@theObject@args$eBayes_robust, "YES")
})

test_that("metabolomics S4 DA helper definition retains parameter capture and limma handoff", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "differentialAbundanceAnalysisHelper")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "group_mapping <- setNames\\(original_groups, safe_groups\\)")
  expect_match(target_text, "updateParamInObject\\(theObject, \"treat_lfc_cutoff\"\\)")
  expect_match(target_text, "runTestsContrasts\\(")
  expect_match(target_text, "new\\(\"MetabolomicsDifferentialAbundanceResults\"")
})
