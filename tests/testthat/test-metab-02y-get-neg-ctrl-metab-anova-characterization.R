# fidelity-coverage-compare: shared
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
      group_id = "character",
      internal_standard_regex = "character",
      annotation_id_column = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      internal_standard_regex = "Internal Standard",
      annotation_id_column = "annotation",
      technical_replicate_id = "TechRep",
      args = list()
    )
  )
}

if (!methods::isGeneric("getNegCtrlMetabAnova")) {
  methods::setGeneric(
    "getNegCtrlMetabAnova",
    function(theObject,
             ruv_grouping_variable = NULL,
             percentage_as_neg_ctrl = NULL,
             num_neg_ctrl = NULL,
             ruv_qval_cutoff = NULL,
             ruv_fdr_method = NULL) {
      standardGeneric("getNegCtrlMetabAnova")
    }
  )
}

checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default_value = NULL) {
  calling_frame <- parent.frame()

  if (exists(param_name_string, envir = calling_frame, inherits = FALSE)) {
    param_value <- get(param_name_string, envir = calling_frame, inherits = FALSE)
    if (!is.null(param_value)) {
      return(param_value)
    }
  }

  object_value <- theObject@args[[param_name_string]]
  if (!is.null(object_value)) {
    return(object_value)
  }

  default_value
}

updateParamInObject <- function(theObject, param_name_string) {
  theObject
}

log_info <- function(...) {
  invisible(NULL)
}

log_warn <- function(...) {
  invisible(NULL)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

if (!methods::hasMethod("getNegCtrlMetabAnova", "MetaboliteAssayData")) {
  loadSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "getNegCtrlMetabAnova")
    },
    env = environment()
  )
}

getTargetMethodEnv <- function(generic_name, signature_name, fallback_env = environment()) {
  if (methods::hasMethod(generic_name, signature_name)) {
    return(environment(methods::getMethod(generic_name, signature_name)))
  }

  fallback_env
}

localBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

helper_capture <- new.env(parent = emptyenv())

getNegCtrlProtAnovaHelper <- function(...) {
  stop("getNegCtrlProtAnovaHelper test stub was not configured")
}

newNegCtrlMetabObject <- function(assay_tbl, assay_name = "LCMS_Pos", args = list(), design_matrix = NULL) {
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      group = c("A", "A", "B", "B"),
      ruv_group = c("Batch_1", "Batch_1", "Batch_2", "Batch_2"),
      stringsAsFactors = FALSE
    )
  }

  assay_list <- if (is.null(assay_name)) {
    list(assay_tbl)
  } else {
    stats::setNames(list(assay_tbl), assay_name)
  }

  methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list,
    metabolite_id_column = "Name",
    design_matrix = design_matrix,
    sample_id = "Run",
    group_id = "group",
    args = args
  )
}

resetHelperCapture <- function() {
  helper_names <- ls(helper_capture, all.names = TRUE)
  if (length(helper_names) > 0) {
    rm(list = helper_names, envir = helper_capture)
  }
}

test_that("metabolomics S4 getNegCtrlMetabAnova forwards filtered assay data and helper arguments", {
  resetHelperCapture()

  localBinding(
    getTargetMethodEnv("getNegCtrlMetabAnova", "MetaboliteAssayData"),
    "getNegCtrlProtAnovaHelper",
    function(data_matrix,
             design_matrix,
             grouping_variable,
             percentage_as_neg_ctrl,
             num_neg_ctrl,
             ruv_qval_cutoff,
             ruv_fdr_method) {
      helper_capture$data_matrix <- data_matrix
      helper_capture$design_matrix <- design_matrix
      helper_capture$grouping_variable <- grouping_variable
      helper_capture$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
      helper_capture$num_neg_ctrl <- num_neg_ctrl
      helper_capture$ruv_qval_cutoff <- ruv_qval_cutoff
      helper_capture$ruv_fdr_method <- ruv_fdr_method

      selected <- seq_len(nrow(data_matrix)) <= num_neg_ctrl
      stats::setNames(selected, rownames(data_matrix))
    }
  )

  assay_tbl <- tibble::tibble(
    Name = paste0("M", seq_len(5)),
    Annotation = letters[seq_len(5)],
    Sample_1 = c("10", "20", "30", "40", NA),
    Sample_2 = c(11, 21, 31, 41, NA),
    Sample_3 = c(12, 22, 32, 42, 99),
    Sample_4 = c(13, 23, 33, 43, NA)
  )

  result <- suppressMessages(
    suppressWarnings(
      getNegCtrlMetabAnova(
        newNegCtrlMetabObject(
          assay_tbl = assay_tbl,
          args = list(
            ruv_fdr_method = "BH"
          )
        ),
        ruv_grouping_variable = "ruv_group",
        percentage_as_neg_ctrl = c(LCMS_Pos = 40)
      )
    )
  )

  expect_named(result, "LCMS_Pos")
  expect_equal(sum(result$LCMS_Pos), 2)
  expect_identical(names(result$LCMS_Pos), paste0("M", seq_len(4)))

  expect_equal(dim(helper_capture$data_matrix), c(4, 4))
  expect_identical(rownames(helper_capture$data_matrix), paste0("M", seq_len(4)))
  expect_identical(colnames(helper_capture$data_matrix), paste0("Sample_", seq_len(4)))
  expect_false("group" %in% colnames(helper_capture$design_matrix))
  expect_true("ruv_group" %in% colnames(helper_capture$design_matrix))
  expect_identical(rownames(helper_capture$design_matrix), paste0("Sample_", seq_len(4)))
  expect_identical(helper_capture$grouping_variable, "ruv_group")
  expect_identical(helper_capture$percentage_as_neg_ctrl, 40)
  expect_identical(helper_capture$num_neg_ctrl, 2L)
  expect_identical(helper_capture$ruv_qval_cutoff, 0.05)
  expect_identical(helper_capture$ruv_fdr_method, "BH")
})

test_that("metabolomics S4 getNegCtrlMetabAnova falls back to current hardcoded defaults for unnamed assays", {
  resetHelperCapture()

  localBinding(
    getTargetMethodEnv("getNegCtrlMetabAnova", "MetaboliteAssayData"),
    "getNegCtrlProtAnovaHelper",
    function(data_matrix,
             design_matrix,
             grouping_variable,
             percentage_as_neg_ctrl,
             num_neg_ctrl,
             ruv_qval_cutoff,
             ruv_fdr_method) {
      helper_capture$grouping_variable <- grouping_variable
      helper_capture$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
      helper_capture$num_neg_ctrl <- num_neg_ctrl
      helper_capture$ruv_qval_cutoff <- ruv_qval_cutoff
      helper_capture$ruv_fdr_method <- ruv_fdr_method

      selected <- seq_len(nrow(data_matrix)) <= num_neg_ctrl
      stats::setNames(selected, rownames(data_matrix))
    }
  )

  assay_tbl <- tibble::tibble(
    Name = paste0("M", seq_len(4)),
    Sample_1 = c(10, 20, 30, 40),
    Sample_2 = c(11, 21, 31, 41),
    Sample_3 = c(12, 22, 32, 42),
    Sample_4 = c(13, 23, 33, 43)
  )

  result <- suppressMessages(
    suppressWarnings(
      getNegCtrlMetabAnova(
        newNegCtrlMetabObject(
          assay_tbl = assay_tbl,
          assay_name = NULL,
          args = list()
        )
      )
    )
  )

  expect_named(result, "Assay_1")
  expect_equal(sum(result$Assay_1), 0)
  expect_identical(helper_capture$grouping_variable, "group")
  expect_identical(helper_capture$percentage_as_neg_ctrl, 10)
  expect_identical(helper_capture$num_neg_ctrl, 0L)
  expect_identical(helper_capture$ruv_qval_cutoff, 0.05)
  expect_identical(helper_capture$ruv_fdr_method, "BH")
})

test_that("metabolomics S4 getNegCtrlMetabAnova source retains helper-dispatch and filtered-result cleanup", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "getNegCtrlMetabAnova")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "updateParamInObject\\(theObject, \"ruv_grouping_variable\"\\)")
  expect_match(target_text, "dplyr::select\\(.*-dplyr::any_of\\(group_id\\)\\)")
  expect_match(target_text, "getNegCtrlProtAnovaHelper\\(")
  expect_match(target_text, "final_control_list <- control_features_list\\[!sapply\\(control_features_list,\\s*is.null\\)\\]")
})
