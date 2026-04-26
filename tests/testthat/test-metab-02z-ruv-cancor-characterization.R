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
      internal_standard_regex = "character",
      annotation_id_column = "character",
      group_id = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      internal_standard_regex = "Internal Standard",
      annotation_id_column = "annotation",
      group_id = "group",
      technical_replicate_id = "TechRep",
      args = list()
    )
  )
}

if (!methods::isGeneric("ruvCancor")) {
  methods::setGeneric(
    "ruvCancor",
    function(theObject,
             ctrl = NULL,
             num_components_to_impute = NULL,
             ruv_grouping_variable = NULL) {
      standardGeneric("ruvCancor")
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

log_error <- function(...) {
  invisible(NULL)
}

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

if (!methods::hasMethod("ruvCancor", "MetaboliteAssayData")) {
  loadSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "ruvCancor")
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
  target_env <- env
  while (!exists(name, envir = target_env, inherits = FALSE) && environmentIsLocked(target_env)) {
    parent_env <- parent.env(target_env)
    if (identical(parent_env, emptyenv())) {
      break
    }
    target_env <- parent_env
  }

  had_binding <- exists(name, envir = target_env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = target_env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, target_env)

  if (was_locked) {
    unlockBinding(name, target_env)
  }
  assign(name, value, envir = target_env)
  if (was_locked) {
    lockBinding(name, target_env)
  }

  withr::defer({
    if (exists(name, envir = target_env, inherits = FALSE) && bindingIsLocked(name, target_env)) {
      unlockBinding(name, target_env)
    }
    if (had_binding) {
      assign(name, old_value, envir = target_env)
    } else if (exists(name, envir = target_env, inherits = FALSE)) {
      rm(list = name, envir = target_env)
    }
    if (was_locked && exists(name, envir = target_env, inherits = FALSE)) {
      lockBinding(name, target_env)
    }
  }, envir = .local_envir)
}

cancor_capture <- new.env(parent = emptyenv())

ruv_cancorplot <- function(...) {
  stop("ruv_cancorplot test stub was not configured")
}

newRuvCancorObject <- function(assay_tbl, assay_name = "LCMS_Pos", args = list(), design_matrix = NULL) {
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Run = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      group = c("A", "A", "B", "B"),
      ruv_group = c("Batch_1", "Batch_1", "Batch_2", "Batch_2"),
      stringsAsFactors = FALSE
    )
  }

  assay_list <- if (is.null(assay_name)) {
    stats::setNames(list(assay_tbl), "Assay_1")
  } else {
    stats::setNames(list(assay_tbl), assay_name)
  }

  sample_ids <- grep("^Sample_", colnames(assay_tbl), value = TRUE)
  valid_design_matrix <- design_matrix[design_matrix$Run %in% sample_ids, , drop = FALSE]

  if (nrow(valid_design_matrix) == 0) {
    valid_design_matrix <- as.data.frame(
      lapply(design_matrix, function(column) rep(NA, length(sample_ids))),
      stringsAsFactors = FALSE
    )
    names(valid_design_matrix) <- names(design_matrix)
    valid_design_matrix$Run <- sample_ids
  }

  object <- methods::new(
    "MetaboliteAssayData",
    metabolite_data = assay_list,
    metabolite_id_column = "Name",
    design_matrix = valid_design_matrix,
    sample_id = "Run",
    args = args
  )

  if (is.null(assay_name)) {
    names(object@metabolite_data) <- ""
  }

  object@design_matrix <- design_matrix
  object
}

resetCancorCapture <- function() {
  capture_names <- ls(cancor_capture, all.names = TRUE)
  if (length(capture_names) > 0) {
    rm(list = capture_names, envir = cancor_capture)
  }
}

test_that("metabolomics S4 ruvCancor forwards filtered data and assay-specific controls", {
  resetCancorCapture()

  localBinding(
    getTargetMethodEnv("ruvCancor", "MetaboliteAssayData"),
    "ruv_cancorplot",
    function(Y, X, ctl) {
      cancor_capture$Y <- Y
      cancor_capture$X <- X
      cancor_capture$ctl <- ctl

      list(Y = Y, X = X, ctl = ctl)
    }
  )

  assay_tbl <- tibble::tibble(
    Name = paste0("M", seq_len(6)),
    Annotation = letters[seq_len(6)],
    Sample_1 = c("10", "20", "30", "40", "50", "60"),
    Sample_2 = c(11, 21, 31, 41, 51, 61),
    Sample_3 = c(12, 22, 32, 42, 52, 62),
    Sample_4 = c(13, 23, 33, 43, 53, 63)
  )

  design_matrix <- data.frame(
    Run = c("Sample_4", "Sample_2", "Sample_1", "Sample_3", "Unused"),
    group = c("B", "A", "A", "B", "Z"),
    ruv_group = c("Batch_2", "Batch_1", "Batch_1", "Batch_2", "Batch_9"),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(
    suppressWarnings(
      ruvCancor(
        newRuvCancorObject(
          assay_tbl = assay_tbl,
          design_matrix = design_matrix
        ),
        ctrl = list(LCMS_Pos = paste0("M", 2:6)),
        ruv_grouping_variable = "ruv_group"
      )
    )
  )

  expect_named(result, "LCMS_Pos")
  expect_identical(
    rownames(cancor_capture$Y),
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(colnames(cancor_capture$Y), paste0("M", seq_len(6)))
  expect_equal(
    unname(cancor_capture$Y[, "M1"]),
    c(13, 11, 10, 12)
  )
  expect_identical(
    cancor_capture$X,
    c("Batch_2", "Batch_1", "Batch_1", "Batch_2")
  )
  expect_identical(
    cancor_capture$ctl,
    c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  )
})

test_that("metabolomics S4 ruvCancor preserves unnamed-assay global control handling", {
  resetCancorCapture()

  localBinding(
    getTargetMethodEnv("ruvCancor", "MetaboliteAssayData"),
    "ruv_cancorplot",
    function(Y, X, ctl) {
      cancor_capture$Y <- Y
      cancor_capture$X <- X
      cancor_capture$ctl <- ctl

      "stub-cancor-plot"
    }
  )

  assay_tbl <- tibble::tibble(
    Name = paste0("M", seq_len(5)),
    Sample_1 = c(1, 2, 3, 4, 5),
    Sample_2 = c(2, 3, 4, 5, 6),
    Sample_3 = c(3, 4, 5, 6, 7),
    Sample_4 = c(4, 5, 6, 7, 8)
  )

  result <- suppressMessages(
    suppressWarnings(
      ruvCancor(
        newRuvCancorObject(
          assay_tbl = assay_tbl,
          assay_name = NULL
        ),
        ctrl = paste0("M", seq_len(5)),
        ruv_grouping_variable = "ruv_group"
      )
    )
  )

  expect_identical(names(result), "")
  expect_identical(result[[1]], "stub-cancor-plot")
  expect_identical(rownames(cancor_capture$Y), paste0("Sample_", seq_len(4)))
  expect_identical(cancor_capture$X, c("Batch_1", "Batch_1", "Batch_2", "Batch_2"))
  expect_identical(cancor_capture$ctl, rep(TRUE, 5))
})

test_that("metabolomics S4 ruvCancor source retains imputation and control-selection flow", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "ruvCancor")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "mixOmics::impute.nipals\\(")
  expect_match(target_text, "ruv_cancorplot\\(")
  expect_match(target_text, "feature_names_in_assay %in%\\s*ctrl_assay")
  expect_match(target_text, "num_controls_found < 5")
})
