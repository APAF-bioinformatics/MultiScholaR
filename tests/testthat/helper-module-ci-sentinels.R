module_ci_sentinel_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

module_ci_sentinel_canonical_vector <- function(values, order_sensitive = TRUE) {
  values <- as.character(unlist(values, use.names = FALSE))
  values[is.na(values)] <- "<NA>"
  if (!isTRUE(order_sensitive)) {
    values <- sort(unique(values))
  }
  paste(seq_along(values), values, sep = ":", collapse = "\n")
}

module_ci_sentinel_hash <- function(value) {
  temp <- tempfile(fileext = ".txt")
  on.exit(unlink(temp), add = TRUE)

  lines <- if (is.data.frame(value)) {
    capture.output(utils::write.table(
      value,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "NA"
    ))
  } else {
    capture.output(dput(value, control = c("keepNA", "keepInteger", "showAttributes")))
  }

  writeLines(lines, temp, useBytes = TRUE)
  unname(tools::md5sum(temp))
}

module_ci_sentinel_digest <- function(values,
                                      kind,
                                      order_sensitive = TRUE,
                                      metadata = list()) {
  values <- as.character(unlist(values, use.names = FALSE))
  values[is.na(values)] <- "<NA>"
  digest_values <- if (isTRUE(order_sensitive)) values else sort(unique(values))
  payload <- module_ci_sentinel_canonical_vector(digest_values, order_sensitive = TRUE)

  structure(
    list(
      kind = kind,
      count = length(values),
      values = values,
      digest_values = digest_values,
      hash = module_ci_sentinel_hash(payload),
      order_sensitive = isTRUE(order_sensitive),
      metadata = metadata
    ),
    class = "module_ci_sentinel_digest"
  )
}

module_ci_sentinel_digest_values <- function(x) {
  if (inherits(x, "module_ci_sentinel_digest")) {
    return(x$values)
  }
  as.character(unlist(x, use.names = FALSE))
}

module_ci_sentinel_apply_expected_changes <- function(values,
                                                      expected_dropped = character(),
                                                      expected_renamed = character()) {
  values <- as.character(unlist(values, use.names = FALSE))
  expected_dropped <- as.character(unlist(expected_dropped, use.names = FALSE))
  values <- values[!values %in% expected_dropped]

  if (length(expected_renamed) > 0L) {
    rename_names <- names(expected_renamed)
    expected_renamed <- as.character(expected_renamed)
    if (is.null(rename_names) || any(!nzchar(rename_names))) {
      stop("expected_renamed must be a named character vector of old_id = new_id", call. = FALSE)
    }
    names(expected_renamed) <- rename_names
    for (old_id in rename_names) {
      values[values == old_id] <- expected_renamed[[old_id]]
    }
  }

  values
}

module_ci_sentinel_assert_identity <- function(before,
                                               after,
                                               expected_dropped = character(),
                                               expected_renamed = character(),
                                               order_sensitive = TRUE,
                                               label = "identity") {
  before_values <- module_ci_sentinel_digest_values(before)
  after_values <- module_ci_sentinel_digest_values(after)
  expected_after <- module_ci_sentinel_apply_expected_changes(
    before_values,
    expected_dropped = expected_dropped,
    expected_renamed = expected_renamed
  )

  if (isTRUE(order_sensitive)) {
    testthat::expect_true(
      identical(after_values, expected_after),
      info = sprintf(
        "%s drift or order change; expected=[%s], actual=[%s]",
        label,
        paste(expected_after, collapse = ", "),
        paste(after_values, collapse = ", ")
      )
    )
  } else {
    unexpected_missing <- setdiff(expected_after, after_values)
    unexpected_added <- setdiff(after_values, expected_after)
    testthat::expect_true(
      length(unexpected_missing) == 0L && length(unexpected_added) == 0L,
      info = sprintf(
        "%s set drift; missing=[%s], added=[%s]",
        label,
        paste(unexpected_missing, collapse = ", "),
        paste(unexpected_added, collapse = ", ")
      )
    )
  }

  invisible(module_ci_sentinel_digest(after_values, kind = label, order_sensitive = order_sensitive))
}

module_ci_sentinel_guess_sample_columns <- function(data) {
  if (!is.data.frame(data)) {
    return(character())
  }

  numeric_columns <- names(data)[vapply(data, is.numeric, logical(1))]
  pattern_columns <- grep("^(WT|KO|RES|S|Run|Sample)[._-]?[0-9]+$", names(data), value = TRUE)
  unique(c(pattern_columns, numeric_columns))
}

module_ci_sentinel_sample_ids <- function(object,
                                          sample_col = NULL,
                                          sample_columns = NULL) {
  if (!is.null(sample_columns)) {
    return(as.character(sample_columns))
  }

  if (is.list(object) && !is.data.frame(object)) {
    for (field in c("sample_ids", "sample_columns", "sample_order")) {
      if (!is.null(object[[field]])) {
        return(as.character(unlist(object[[field]], use.names = FALSE)))
      }
    }
  }

  if (methods::is(object, "S4")) {
    slot_names <- methods::slotNames(object)
    if (all(c("design_matrix", "sample_id") %in% slot_names)) {
      design <- object@design_matrix
      sample_slot <- object@sample_id
      if (is.data.frame(design) && length(sample_slot) == 1L && sample_slot %in% names(design)) {
        return(as.character(design[[sample_slot]]))
      }
    }
  }

  if (is.data.frame(object)) {
    candidate_cols <- unique(c(sample_col, "sample", "Run", "Sample", "sample_id"))
    candidate_cols <- candidate_cols[!is.na(candidate_cols) & nzchar(candidate_cols)]
    for (candidate in candidate_cols) {
      if (candidate %in% names(object)) {
        return(as.character(object[[candidate]]))
      }
    }
    return(module_ci_sentinel_guess_sample_columns(object))
  }

  character()
}

module_ci_sentinel_feature_ids <- function(object, feature_col = NULL) {
  if (is.list(object) && !is.data.frame(object) && !is.null(object$feature_ids)) {
    return(as.character(unlist(object$feature_ids, use.names = FALSE)))
  }

  if (methods::is(object, "S4")) {
    slot_names <- methods::slotNames(object)
    if (all(c("protein_quant_table", "protein_id_column") %in% slot_names)) {
      table <- object@protein_quant_table
      id_col <- object@protein_id_column
      if (is.data.frame(table) && id_col %in% names(table)) {
        return(as.character(table[[id_col]]))
      }
    }
    if (all(c("lipid_data", "lipid_id_column") %in% slot_names)) {
      lipid_data <- object@lipid_data
      id_col <- object@lipid_id_column
      return(unique(unlist(lapply(lipid_data, function(assay) {
        if (is.data.frame(assay) && id_col %in% names(assay)) {
          as.character(assay[[id_col]])
        } else {
          character()
        }
      }), use.names = FALSE)))
    }
    if (all(c("metabolite_data", "metabolite_id_column") %in% slot_names)) {
      metabolite_data <- object@metabolite_data
      id_col <- object@metabolite_id_column
      return(unique(unlist(lapply(metabolite_data, function(assay) {
        if (is.data.frame(assay) && id_col %in% names(assay)) {
          as.character(assay[[id_col]])
        } else {
          character()
        }
      }), use.names = FALSE)))
    }
  }

  if (is.data.frame(object)) {
    candidates <- unique(c(
      feature_col,
      "Protein",
      "Protein.Group",
      "protein_id",
      "Peptide",
      "Feature.Name",
      "feature_id",
      "database_identifier",
      "LipidName",
      "lipid_id",
      "lipid",
      "Metabolite"
    ))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    for (candidate in candidates) {
      if (candidate %in% names(object)) {
        return(as.character(object[[candidate]]))
      }
    }
  }

  character()
}

module_ci_sentinel_assay_names <- function(object) {
  if (exists("module_ci_extract_assay_names", mode = "function")) {
    extracted <- module_ci_extract_assay_names(object)
    if (length(extracted) > 0L) {
      return(as.character(extracted))
    }
  }

  if (is.data.frame(object)) {
    for (field in c("assay", "Assay", "assay_name")) {
      if (field %in% names(object)) {
        return(unique(as.character(object[[field]])))
      }
    }
  }
  if (is.list(object) && !is.null(object$assay_names)) {
    return(as.character(unlist(object$assay_names, use.names = FALSE)))
  }
  if (is.list(object) && !is.null(names(object)) && all(nzchar(names(object)))) {
    return(names(object))
  }

  character()
}

module_ci_sentinel_assert_assays <- function(before, after, label = "assay provenance") {
  before_assays <- module_ci_sentinel_assay_names(before)
  after_assays <- module_ci_sentinel_assay_names(after)
  testthat::expect_true(
    setequal(after_assays, before_assays),
    info = sprintf(
      "%s drift; expected=[%s], actual=[%s]",
      label,
      paste(before_assays, collapse = ", "),
      paste(after_assays, collapse = ", ")
    )
  )
  invisible(module_ci_sentinel_digest(after_assays, kind = "assay", order_sensitive = FALSE))
}

module_ci_sentinel_path_value <- function(object, path, missing = structure(list(), class = "module_ci_sentinel_missing")) {
  parts <- strsplit(path, ".", fixed = TRUE)[[1]]
  current <- object

  for (part in parts) {
    if (inherits(current, "module_ci_sentinel_missing")) {
      return(missing)
    }
    if (is.environment(current)) {
      if (!exists(part, envir = current, inherits = FALSE)) {
        return(missing)
      }
      current <- get(part, envir = current, inherits = FALSE)
    } else if (is.list(current)) {
      if (is.null(current[[part]])) {
        return(missing)
      }
      current <- current[[part]]
    } else if (methods::is(current, "S4")) {
      if (!part %in% methods::slotNames(current)) {
        return(missing)
      }
      current <- methods::slot(current, part)
    } else {
      return(missing)
    }
  }

  current
}

module_ci_sentinel_value_string <- function(value) {
  paste(capture.output(dput(value, control = c("keepNA", "keepInteger"))), collapse = "\n")
}

module_ci_sentinel_assert_parameters <- function(live_state,
                                                 serialized_state,
                                                 paths) {
  paths <- as.character(paths)
  problems <- character()
  for (path in paths) {
    live_value <- module_ci_sentinel_path_value(live_state, path)
    serialized_value <- module_ci_sentinel_path_value(serialized_state, path)
    if (inherits(live_value, "module_ci_sentinel_missing")) {
      problems <- c(problems, sprintf("live_missing:%s", path))
      next
    }
    if (inherits(serialized_value, "module_ci_sentinel_missing")) {
      problems <- c(problems, sprintf("serialized_missing:%s", path))
      next
    }
    if (!identical(module_ci_sentinel_value_string(serialized_value), module_ci_sentinel_value_string(live_value))) {
      problems <- c(problems, sprintf("drift:%s", path))
    }
  }
  testthat::expect_identical(problems, character())
  invisible(TRUE)
}

module_ci_sentinel_env_snapshot <- function(env,
                                            names = ls(env, all.names = TRUE)) {
  stats::setNames(lapply(names, function(name) {
    exists_now <- exists(name, envir = env, inherits = FALSE)
    value <- if (exists_now) get(name, envir = env, inherits = FALSE) else NULL
    list(
      exists = exists_now,
      class = if (exists_now) class(value) else character(),
      hash = if (exists_now) module_ci_sentinel_hash(value) else NA_character_
    )
  }), names)
}

module_ci_sentinel_assert_env_unchanged <- function(before,
                                                   after,
                                                   allowed_changed = character(),
                                                   allowed_added = character(),
                                                   allowed_removed = character()) {
  all_names <- union(names(before), names(after))
  unexpected <- character()

  for (name in all_names) {
    before_exists <- isTRUE(before[[name]]$exists)
    after_exists <- isTRUE(after[[name]]$exists)
    if (!before_exists && after_exists && !name %in% allowed_added) {
      unexpected <- c(unexpected, sprintf("added:%s", name))
    } else if (before_exists && !after_exists && !name %in% allowed_removed) {
      unexpected <- c(unexpected, sprintf("removed:%s", name))
    } else if (before_exists && after_exists && !identical(before[[name]]$hash, after[[name]]$hash) && !name %in% allowed_changed) {
      unexpected <- c(unexpected, sprintf("changed:%s", name))
    }
  }

  testthat::expect_identical(unexpected, character())
  invisible(TRUE)
}

module_ci_sentinel_assert_project_dirs_isolated <- function(project_dirs,
                                                           omics = names(project_dirs)) {
  omics <- as.character(omics)
  problems <- character()
  missing_omics <- setdiff(omics, names(project_dirs))
  if (length(missing_omics) > 0L) {
    problems <- c(problems, sprintf("missing_omics:%s", paste(missing_omics, collapse = ",")))
  }

  path_rows <- do.call(rbind, lapply(omics, function(omic) {
    paths <- unlist(project_dirs[[omic]], recursive = TRUE, use.names = TRUE)
    paths <- paths[nzchar(paths)]
    data.frame(omic = omic, field = names(paths), path = as.character(paths), stringsAsFactors = FALSE)
  }))

  duplicated_paths <- path_rows$path[duplicated(path_rows$path)]
  if (length(duplicated_paths) > 0L) {
    problems <- c(problems, sprintf("shared_paths:%s", paste(unique(duplicated_paths), collapse = ",")))
  }
  testthat::expect_identical(problems, character())
  invisible(path_rows)
}

module_ci_sentinel_read_artifact <- function(path, type) {
  switch(
    type,
    tsv = utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    xlsx = {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("openxlsx is required for xlsx artifact schema validation", call. = FALSE)
      }
      openxlsx::read.xlsx(path)
    },
    rds = readRDS(path),
    html = readLines(path, warn = FALSE),
    text = readLines(path, warn = FALSE),
    stop(sprintf("Unsupported artifact type: %s", type), call. = FALSE)
  )
}

module_ci_sentinel_assert_artifact_schema <- function(schema) {
  path <- schema$path
  type <- schema$type
  if (is.null(path) || !nzchar(path)) {
    stop("artifact schema path must be non-empty", call. = FALSE)
  }
  if (is.null(type) || !nzchar(type)) {
    stop("artifact schema type must be non-empty", call. = FALSE)
  }

  problems <- character()
  if (!file.exists(path)) {
    problems <- c(problems, sprintf("missing:%s", path))
    testthat::expect_identical(problems, character())
    return(invisible(NULL))
  }
  if (is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    problems <- c(problems, sprintf("empty:%s", path))
  }

  artifact <- tryCatch(
    module_ci_sentinel_read_artifact(path, type),
    error = function(e) {
      problems <<- c(problems, sprintf("read_error:%s:%s", path, conditionMessage(e)))
      NULL
    }
  )

  if (is.data.frame(artifact)) {
    required_columns <- module_ci_sentinel_coalesce(schema$required_columns, character())
    min_rows <- module_ci_sentinel_coalesce(schema$min_rows, 1L)
    if (nrow(artifact) < min_rows) {
      problems <- c(problems, sprintf("min_rows:%s", path))
    }
    missing_columns <- setdiff(required_columns, names(artifact))
    if (length(missing_columns) > 0L) {
      problems <- c(problems, sprintf("missing_columns:%s:%s", path, paste(missing_columns, collapse = ",")))
    }
  } else if (identical(type, "rds") && !is.null(artifact)) {
    expected_class <- schema$class
    if (!is.null(expected_class)) {
      if (!inherits(artifact, expected_class)) {
        problems <- c(problems, sprintf("class:%s:%s", path, expected_class))
      }
    }
  } else if (!is.null(artifact)) {
    contains <- module_ci_sentinel_coalesce(schema$contains, character())
    for (pattern in contains) {
      if (!any(grepl(pattern, artifact, fixed = TRUE))) {
        problems <- c(problems, sprintf("missing_text:%s:%s", path, pattern))
      }
    }
  }

  testthat::expect_identical(problems, character())
  invisible(artifact)
}
