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
      metabolite_id_column = "character"
    )
  )
}

if (!methods::isGeneric("plotVolcanoS4")) {
  methods::setGeneric(
    "plotVolcanoS4",
    function(
      objectsList,
      da_q_val_thresh = 0.05,
      qvalue_column = "fdr_qvalue",
      log2fc_column = "logFC"
    ) {
      standardGeneric("plotVolcanoS4")
    }
  )
}

bind_rows <- dplyr::bind_rows
mutate <- dplyr::mutate
sym <- rlang::sym
case_when <- dplyr::case_when
group_by <- dplyr::group_by
nest <- tidyr::nest
str_split_i <- stringr::str_split_i

scale_color_manual <- function(values, name, limits) {
  structure(
    list(
      values = values,
      name = name,
      limits = limits
    ),
    class = "mock_scale"
  )
}

plotOneVolcanoNoVerticalLines <- function(
  input_data,
  input_title,
  log_q_value_column,
  log_fc_column
) {
  structure(
    list(
      data = input_data,
      title = input_title,
      log_q_value_column = log_q_value_column,
      log_fc_column = log_fc_column
    ),
    class = "mock_plot"
  )
}

`+.mock_plot` <- function(e1, e2) {
  e1$scale <- e2
  e1
}

loadSelectedExpressions(
  paths = c(
    file.path(repo_root, "R", "func_metab_s4_da_results.R"),
    file.path(repo_root, "R", "func_metab_s4_objects.R")
  ),
  matcher = function(expr) {
    isTargetSetClass(expr, "MetabolomicsDifferentialAbundanceResults") ||
      isTargetSetMethod(expr, "plotVolcanoS4")
  },
  env = environment()
)

test_that("metabolomics DA volcano list-method preserves grouped plot names and significance scales", {
  contrasts_results_table <- list(
    "groupB-groupA = B vs A" = data.frame(
      metabolite_id = c("M1", "M2", "M3"),
      fdr_qvalue = c(0.001, 0.02, 0.6),
      logFC = c(1.8, -1.2, 0.1),
      stringsAsFactors = FALSE
    ),
    "groupC-groupA = C vs A" = data.frame(
      metabolite_id = c("M4", "M5"),
      fdr_qvalue = c(0.04, 0.8),
      logFC = c(0.7, -0.3),
      stringsAsFactors = FALSE
    )
  )

  result_object <- methods::new(
    "MetabolomicsDifferentialAbundanceResults",
    theObject = methods::new(
      "MetaboliteAssayData",
      metabolite_data = list(),
      metabolite_id_column = "metabolite_id"
    ),
    contrasts_results_table = contrasts_results_table
  )

  output <- plotVolcanoS4(list(primary = result_object))

  expect_identical(names(output), "primary")
  expect_identical(names(output$primary@volcano_plot), names(contrasts_results_table))

  primary_plot <- output$primary@volcano_plot[[1]]
  secondary_plot <- output$primary@volcano_plot[[2]]

  expect_s3_class(primary_plot, "mock_plot")
  expect_s3_class(secondary_plot, "mock_plot")
  expect_match(primary_plot$title, "^primary - groupB-groupA\\s*$")
  expect_match(secondary_plot$title, "^primary - groupC-groupA\\s*$")
  expect_identical(primary_plot$log_q_value_column, "lqm")
  expect_identical(primary_plot$log_fc_column, "logFC")
  expect_equal(
    as.character(primary_plot$data$label),
    c("Significant Up", "Significant Down", "Not significant")
  )
  expect_identical(
    levels(primary_plot$data$label),
    c("Significant Up", "Significant Down", "Not significant")
  )
  expect_equal(as.character(primary_plot$data$colour), c("red", "blue", "grey"))
  expect_identical(levels(primary_plot$data$colour), c("blue", "grey", "red"))
  expect_equal(
    round(primary_plot$data$lqm, 6),
    round(-log10(c(0.001, 0.02, 0.6)), 6)
  )
  expect_identical(
    primary_plot$scale$values,
    c(
      "Significant Up" = "red",
      "Significant Down" = "blue",
      "Not significant" = "grey"
    )
  )
  expect_identical(
    primary_plot$scale$limits,
    c("Significant Up", "Significant Down", "Not significant")
  )
})
