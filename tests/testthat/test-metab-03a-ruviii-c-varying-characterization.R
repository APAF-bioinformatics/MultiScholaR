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
      technical_replicate_id = "character",
      internal_standard_regex = "character",
      annotation_id_column = "character",
      args = "list"
    ),
    prototype = list(
      metabolite_data = list(),
      metabolite_id_column = "Name",
      design_matrix = data.frame(),
      sample_id = "Run",
      group_id = "group",
      technical_replicate_id = "TechRep",
      internal_standard_regex = "Internal Standard",
      annotation_id_column = "annotation",
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

if (!methods::isGeneric("ruvIII_C_Varying")) {
  methods::setGeneric(
    "ruvIII_C_Varying",
    function(theObject,
             ruv_grouping_variable = NULL,
             ruv_number_k = NULL,
             ctrl = NULL) {
      standardGeneric("ruvIII_C_Varying")
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

log_info <- function(...) {
  invisible(NULL)
}

log_warn <- function(...) {
  invisible(NULL)
}

log_error <- function(...) {
  invisible(NULL)
}

is_replicate_temp <- "is_replicate_temp"

target_paths <- c(
  file.path(repo_root, "R", "func_metab_s4_norm_methods.R"),
  file.path(repo_root, "R", "func_metab_s4_objects.R")
)

if (!methods::hasMethod("cleanDesignMatrix", "MetaboliteAssayData") ||
    !methods::hasMethod("ruvIII_C_Varying", "MetaboliteAssayData")) {
  loadSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "cleanDesignMatrix") ||
        isTargetSetMethod(expr, "ruvIII_C_Varying")
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

helper_capture <- new.env(parent = emptyenv())
ruv_capture <- new.env(parent = emptyenv())

resetRuvIIIHarness <- function() {
  helper_names <- ls(helper_capture, all.names = TRUE)
  if (length(helper_names) > 0) {
    rm(list = helper_names, envir = helper_capture)
  }

  ruv_names <- ls(ruv_capture, all.names = TRUE)
  if (length(ruv_names) > 0) {
    rm(list = ruv_names, envir = ruv_capture)
  }
}

getRuvIIIReplicateMatrixHelper <- function(design_matrix, sample_id_column, grouping_variable, temp_column = is_replicate_temp) {
  sample_col <- rlang::as_name(rlang::ensym(sample_id_column))
  grouping_col <- rlang::as_name(rlang::ensym(grouping_variable))

  helper_capture$design_matrix <- design_matrix
  helper_capture$sample_col <- sample_col
  helper_capture$grouping_col <- grouping_col

  replicate_matrix <- stats::model.matrix(
    ~ 0 + group_label,
    data = data.frame(group_label = design_matrix[[grouping_col]], stringsAsFactors = FALSE)
  )
  rownames(replicate_matrix) <- as.character(design_matrix[[sample_col]])
  colnames(replicate_matrix) <- sub("^group_label", "", colnames(replicate_matrix))
  replicate_matrix
}

RUVIII_C_Varying <- function(k, Y, M, toCorrect, potentialControls) {
  ruv_capture$k <- k
  ruv_capture$Y <- Y
  ruv_capture$M <- M
  ruv_capture$toCorrect <- toCorrect
  ruv_capture$potentialControls <- potentialControls

  corrected <- Y + 100
  dimnames(corrected) <- dimnames(Y)
  corrected
}

newRuvIIIObject <- function(assay_tbl, assay_name = "LCMS_Pos", args = list(), design_matrix = NULL) {
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Run = c("Sample_4", "Sample_2", "Sample_1", "Sample_3", "Unused"),
      group = c("B", "A", "A", "B", "Z"),
      ruv_group = c("Batch_2", "Batch_1", "Batch_1", "Batch_2", "Batch_9"),
      TechRep = c("TR2", "TR1", "TR1", "TR2", "TRX"),
      stringsAsFactors = FALSE
    )
  }

  assay_list <- if (is.null(assay_name)) {
    list(assay_tbl)
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
    group_id = "group",
    technical_replicate_id = "TechRep",
    args = args
  )

  object@design_matrix <- design_matrix
  object
}

test_that("metabolomics S4 ruvIII_C_Varying forwards per-assay controls and cleans reordered design rows", {
  resetRuvIIIHarness()

  method_env <- getTargetMethodEnv("ruvIII_C_Varying", "MetaboliteAssayData")
  localBinding(method_env, "getRuvIIIReplicateMatrixHelper", getRuvIIIReplicateMatrixHelper)
  localBinding(method_env, "RUVIII_C_Varying", RUVIII_C_Varying)

  assay_tbl <- tibble::tibble(
    Name = c("M1", "M2", "M3"),
    Annotation = c("alpha", "beta", "gamma"),
    Sample_1 = c(10, 20, 30),
    Sample_2 = c(11, 21, 31),
    Sample_3 = c(12, 22, 32),
    Sample_4 = c(13, 23, 33)
  )

  corrected <- suppressMessages(
    suppressWarnings(
      ruvIII_C_Varying(
        newRuvIIIObject(
          assay_tbl = assay_tbl
        ),
        ruv_grouping_variable = "ruv_group",
        ruv_number_k = list(LCMS_Pos = 2),
        ctrl = list(LCMS_Pos = c("M2", "M3"))
      )
    )
  )

  corrected_assay <- corrected@metabolite_data$LCMS_Pos

  expect_identical(helper_capture$sample_col, "Run")
  expect_identical(helper_capture$grouping_col, "ruv_group")
  expect_identical(
    helper_capture$design_matrix$Run,
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(rownames(ruv_capture$Y), c("Sample_4", "Sample_2", "Sample_1", "Sample_3"))
  expect_identical(colnames(ruv_capture$Y), c("M1", "M2", "M3"))
  expect_identical(rownames(ruv_capture$M), rownames(ruv_capture$Y))
  expect_identical(ruv_capture$k, 2L)
  expect_identical(ruv_capture$toCorrect, c("M1", "M2", "M3"))
  expect_identical(ruv_capture$potentialControls, c("M2", "M3"))
  expect_equal(
    corrected_assay$Sample_4,
    c(113, 123, 133)
  )
  expect_equal(
    corrected_assay$Sample_2,
    c(111, 121, 131)
  )
  expect_identical(
    corrected@design_matrix$Run,
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(corrected@args$ruv_grouping_variable, "ruv_group")
  expect_identical(corrected@args$ruv_number_k$LCMS_Pos, 2)
  expect_identical(corrected@args$ctrl$LCMS_Pos, c("M2", "M3"))
})

test_that("metabolomics S4 ruvIII_C_Varying source retains helper dispatch and cleanup wiring", {
  target_expr <- findSelectedExpression(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "ruvIII_C_Varying")
    }
  )

  expect_false(is.null(target_expr))

  target_text <- paste(deparse(target_expr), collapse = "\n")

  expect_match(target_text, "getRuvIIIReplicateMatrixHelper\\(")
  expect_match(target_text, "RUVIII_C_Varying\\(")
  expect_match(target_text, "potentialControls = potential_controls_names")
  expect_match(target_text, "cleanDesignMatrix\\(theObject\\)")
})
