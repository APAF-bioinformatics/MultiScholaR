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

if (!methods::isClass("LipidomicsAssayData")) {
  methods::setClass(
    "LipidomicsAssayData",
    slots = c(
      lipid_data = "list",
      lipid_id_column = "character",
      annotation_id_column = "character",
      database_identifier_type = "character",
      internal_standard_regex = "character",
      design_matrix = "data.frame",
      sample_id = "character",
      group_id = "character",
      technical_replicate_id = "character",
      args = "list"
    ),
    prototype = list(
      lipid_data = list(),
      lipid_id_column = "database_identifier",
      annotation_id_column = "lipid_identification",
      database_identifier_type = "MockDB",
      internal_standard_regex = "ISTD",
      design_matrix = data.frame(),
      sample_id = "Sample_ID",
      group_id = "group",
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

is_replicate_temp <- "is_replicate_temp"

target_paths <- c(
  file.path(repo_root, "R", "func_lipid_norm_ruv_helpers.R"),
  file.path(repo_root, "R", "func_lipid_s4_objects.R")
)

if (!methods::hasMethod("getNegCtrlMetabAnova", "LipidomicsAssayData") ||
    !methods::hasMethod("ruvCancor", "LipidomicsAssayData") ||
    !methods::hasMethod("ruvIII_C_Varying", "LipidomicsAssayData")) {
  loadSelectedExpressions(
    paths = target_paths,
    matcher = function(expr) {
      isTargetSetMethod(expr, "getNegCtrlMetabAnova") ||
        isTargetSetMethod(expr, "ruvCancor") ||
        isTargetSetMethod(expr, "ruvIII_C_Varying")
    },
    env = environment()
  )
}

if (!methods::hasMethod("cleanDesignMatrix", "LipidomicsAssayData")) {
  methods::setMethod(
    f = "cleanDesignMatrix",
    signature = "LipidomicsAssayData",
    definition = function(theObject) {
      assay_list <- theObject@lipid_data
      if (length(assay_list) == 0) {
        return(theObject)
      }

      sample_id_col <- theObject@sample_id
      design_matrix <- theObject@design_matrix
      sample_cols <- grep("^Sample_", colnames(assay_list[[1]]), value = TRUE)
      sample_index <- match(sample_cols, as.character(design_matrix[[sample_id_col]]))
      sample_index <- sample_index[!is.na(sample_index)]

      theObject@design_matrix <- design_matrix[sample_index, , drop = FALSE]
      theObject
    }
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

newLipidRuvObject <- function(assay_tbl,
                              assay_name = "Positive Mode",
                              args = list(),
                              design_matrix = NULL) {
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Sample_ID = c("Sample_1", "Sample_2", "Sample_3", "Sample_4"),
      group = c("A", "A", "B", "B"),
      ruv_group = c("Batch_1", "Batch_1", "Batch_2", "Batch_2"),
      TechRep = c("TR1", "TR1", "TR2", "TR2"),
      stringsAsFactors = FALSE
    )
  }

  assay_list <- if (is.null(assay_name)) {
    list(assay_tbl)
  } else {
    stats::setNames(list(assay_tbl), assay_name)
  }

  sample_ids <- grep("^Sample_", colnames(assay_tbl), value = TRUE)
  valid_design_matrix <- design_matrix[design_matrix$Sample_ID %in% sample_ids, , drop = FALSE]

  if (nrow(valid_design_matrix) == 0) {
    valid_design_matrix <- data.frame(
      Sample_ID = sample_ids,
      group = rep("group", length(sample_ids)),
      ruv_group = rep("batch", length(sample_ids)),
      TechRep = rep("rep", length(sample_ids)),
      stringsAsFactors = FALSE
    )
  }

  object <- methods::new(
    "LipidomicsAssayData",
    lipid_data = assay_list,
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
    database_identifier_type = "MockDB",
    internal_standard_regex = "ISTD",
    design_matrix = valid_design_matrix,
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = "TechRep",
    args = args
  )

  object@design_matrix <- design_matrix
  object
}

helper_capture <- new.env(parent = emptyenv())
cancor_capture <- new.env(parent = emptyenv())
ruv_capture <- new.env(parent = emptyenv())

resetCapture <- function(env) {
  capture_names <- ls(env, all.names = TRUE)
  if (length(capture_names) > 0) {
    rm(list = capture_names, envir = env)
  }
}

test_that("lipidomics S4 getNegCtrlMetabAnova forwards filtered assay data and helper arguments", {
  resetCapture(helper_capture)

  localBinding(
    getTargetMethodEnv("getNegCtrlMetabAnova", "LipidomicsAssayData"),
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
    database_identifier = paste0("L", seq_len(5)),
    lipid_identification = letters[seq_len(5)],
    Sample_1 = c("10", "20", "30", "40", NA),
    Sample_2 = c(11, 21, 31, 41, NA),
    Sample_3 = c(12, 22, 32, 42, 99),
    Sample_4 = c(13, 23, 33, 43, NA)
  )

  result <- suppressMessages(
    suppressWarnings(
      getNegCtrlMetabAnova(
        newLipidRuvObject(
          assay_tbl = assay_tbl,
          args = list(ruv_fdr_method = "BH")
        ),
        ruv_grouping_variable = "ruv_group",
        percentage_as_neg_ctrl = stats::setNames(40, "Positive Mode")
      )
    )
  )

  expect_named(result, "Positive Mode")
  expect_equal(sum(result$`Positive Mode`), 2)
  expect_identical(names(result$`Positive Mode`), paste0("L", seq_len(4)))

  expect_equal(dim(helper_capture$data_matrix), c(4, 4))
  expect_identical(rownames(helper_capture$data_matrix), paste0("L", seq_len(4)))
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

test_that("lipidomics S4 getNegCtrlMetabAnova falls back to current hardcoded defaults for unnamed assays", {
  resetCapture(helper_capture)

  localBinding(
    getTargetMethodEnv("getNegCtrlMetabAnova", "LipidomicsAssayData"),
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
    database_identifier = paste0("L", seq_len(4)),
    Sample_1 = c(10, 20, 30, 40),
    Sample_2 = c(11, 21, 31, 41),
    Sample_3 = c(12, 22, 32, 42),
    Sample_4 = c(13, 23, 33, 43)
  )

  result <- suppressMessages(
    suppressWarnings(
      getNegCtrlMetabAnova(
        newLipidRuvObject(
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

test_that("lipidomics S4 ruvCancor forwards filtered data and assay-specific controls", {
  resetCapture(cancor_capture)

  localBinding(
    getTargetMethodEnv("ruvCancor", "LipidomicsAssayData"),
    "ruv_cancorplot",
    function(Y, X, ctl) {
      cancor_capture$Y <- Y
      cancor_capture$X <- X
      cancor_capture$ctl <- ctl

      list(Y = Y, X = X, ctl = ctl)
    }
  )

  assay_tbl <- tibble::tibble(
    database_identifier = paste0("L", seq_len(6)),
    lipid_identification = letters[seq_len(6)],
    Sample_1 = c("10", "20", "30", "40", "50", "60"),
    Sample_2 = c(11, 21, 31, 41, 51, 61),
    Sample_3 = c(12, 22, 32, 42, 52, 62),
    Sample_4 = c(13, 23, 33, 43, 53, 63)
  )

  design_matrix <- data.frame(
    Sample_ID = c("Sample_4", "Sample_2", "Sample_1", "Sample_3", "Unused"),
    group = c("B", "A", "A", "B", "Z"),
    ruv_group = c("Batch_2", "Batch_1", "Batch_1", "Batch_2", "Batch_9"),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(
    suppressWarnings(
      ruvCancor(
        newLipidRuvObject(
          assay_tbl = assay_tbl,
          design_matrix = design_matrix
        ),
        ctrl = stats::setNames(list(paste0("L", 2:6)), "Positive Mode"),
        ruv_grouping_variable = "ruv_group"
      )
    )
  )

  expect_named(result, "Positive Mode")
  expect_identical(
    rownames(cancor_capture$Y),
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(colnames(cancor_capture$Y), paste0("L", seq_len(6)))
  expect_equal(unname(cancor_capture$Y[, "L1"]), c(13, 11, 10, 12))
  expect_identical(cancor_capture$X, c("Batch_2", "Batch_1", "Batch_1", "Batch_2"))
  expect_identical(cancor_capture$ctl, c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
})

test_that("lipidomics S4 ruvCancor falls back to default assay names and global controls", {
  resetCapture(cancor_capture)

  localBinding(
    getTargetMethodEnv("ruvCancor", "LipidomicsAssayData"),
    "ruv_cancorplot",
    function(Y, X, ctl) {
      cancor_capture$Y <- Y
      cancor_capture$X <- X
      cancor_capture$ctl <- ctl

      "stub-cancor-plot"
    }
  )

  assay_tbl <- tibble::tibble(
    database_identifier = paste0("L", seq_len(5)),
    Sample_1 = c(1, 2, 3, 4, 5),
    Sample_2 = c(2, 3, 4, 5, 6),
    Sample_3 = c(3, 4, 5, 6, 7),
    Sample_4 = c(4, 5, 6, 7, 8)
  )

  result <- suppressMessages(
    suppressWarnings(
      ruvCancor(
        newLipidRuvObject(
          assay_tbl = assay_tbl,
          assay_name = NULL
        ),
        ctrl = paste0("L", seq_len(5)),
        ruv_grouping_variable = "ruv_group"
      )
    )
  )

  expect_named(result, "Assay_1")
  expect_identical(result$Assay_1, "stub-cancor-plot")
  expect_identical(rownames(cancor_capture$Y), paste0("Sample_", seq_len(4)))
  expect_identical(cancor_capture$X, c("Batch_1", "Batch_1", "Batch_2", "Batch_2"))
  expect_identical(cancor_capture$ctl, rep(TRUE, 5))
})

test_that("lipidomics S4 ruvIII_C_Varying forwards per-assay controls and cleans reordered design rows", {
  resetCapture(helper_capture)
  resetCapture(ruv_capture)

  method_env <- getTargetMethodEnv("ruvIII_C_Varying", "LipidomicsAssayData")
  localBinding(
    method_env,
    "getRuvIIIReplicateMatrixHelper",
    function(design_matrix, sample_id_column, grouping_variable, temp_column = is_replicate_temp) {
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
  )
  localBinding(
    method_env,
    "RUVIII_C_Varying",
    function(k, Y, M, toCorrect, potentialControls) {
      ruv_capture$k <- k
      ruv_capture$Y <- Y
      ruv_capture$M <- M
      ruv_capture$toCorrect <- toCorrect
      ruv_capture$potentialControls <- potentialControls

      corrected <- Y + 100
      dimnames(corrected) <- dimnames(Y)
      corrected
    }
  )

  assay_tbl <- tibble::tibble(
    database_identifier = c("L1", "L2", "L3"),
    lipid_identification = c("alpha", "beta", "gamma"),
    Sample_1 = c(10, 20, 30),
    Sample_2 = c(11, 21, 31),
    Sample_3 = c(12, 22, 32),
    Sample_4 = c(13, 23, 33)
  )

  design_matrix <- data.frame(
    Sample_ID = c("Sample_4", "Sample_2", "Sample_1", "Sample_3", "Unused"),
    group = c("B", "A", "A", "B", "Z"),
    ruv_group = c("Batch_2", "Batch_1", "Batch_1", "Batch_2", "Batch_9"),
    TechRep = c("TR2", "TR1", "TR1", "TR2", "TRX"),
    stringsAsFactors = FALSE
  )

  corrected <- suppressMessages(
    suppressWarnings(
      ruvIII_C_Varying(
        newLipidRuvObject(
          assay_tbl = assay_tbl,
          design_matrix = design_matrix
        ),
        ruv_grouping_variable = "ruv_group",
        ruv_number_k = stats::setNames(list(2), "Positive Mode"),
        ctrl = stats::setNames(list(c("L2", "L3")), "Positive Mode")
      )
    )
  )

  corrected_assay <- corrected@lipid_data$`Positive Mode`

  expect_identical(helper_capture$sample_col, "Sample_ID")
  expect_identical(helper_capture$grouping_col, "ruv_group")
  expect_identical(
    helper_capture$design_matrix$Sample_ID,
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(rownames(ruv_capture$Y), c("Sample_4", "Sample_2", "Sample_1", "Sample_3"))
  expect_identical(colnames(ruv_capture$Y), c("L1", "L2", "L3"))
  expect_identical(rownames(ruv_capture$M), rownames(ruv_capture$Y))
  expect_identical(ruv_capture$k, 2L)
  expect_identical(ruv_capture$toCorrect, c("L1", "L2", "L3"))
  expect_identical(ruv_capture$potentialControls, c("L2", "L3"))
  expect_equal(corrected_assay$Sample_4, c(113, 123, 133))
  expect_equal(corrected_assay$Sample_2, c(111, 121, 131))
  expect_identical(
    corrected@design_matrix$Sample_ID,
    c("Sample_4", "Sample_2", "Sample_1", "Sample_3")
  )
  expect_identical(corrected@args$ruv_grouping_variable, "ruv_group")
  expect_identical(corrected@args$ruv_number_k$`Positive Mode`, 2)
  expect_identical(corrected@args$ctrl$`Positive Mode`, c("L2", "L3"))
})
