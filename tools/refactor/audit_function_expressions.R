#!/usr/bin/env Rscript

# Recursively inventory every function expression in R/, including anonymous
# callbacks and functions nested inside Shiny modules.

abort <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

parse_args <- function(args) {
  result <- list(
    repo_root = ".",
    baseline_ref = "main",
    target_root = ".",
    output_dir = ".refactor-fidelity-audit/function-expressions"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--repo-root", "--baseline-ref", "--target-root", "--output-dir")) {
      abort("Unknown argument: ", key)
    }
    if (i == length(args)) {
      abort("Missing value for ", key)
    }
    result[[gsub("-", "_", sub("^--", "", key))]] <- args[[i + 1L]]
    i <- i + 2L
  }

  result
}

call_name <- function(expr) {
  if (!is.call(expr)) {
    return(NULL)
  }
  head <- expr[[1L]]
  if (is.symbol(head)) as.character(head) else NULL
}

ast_text <- function(expr) {
  paste(
    deparse(
      expr,
      width.cutoff = 500L,
      control = c("keepNA", "keepInteger", "niceNames", "hexNumeric")
    ),
    collapse = "\n"
  )
}

ast_hash <- function(expr) {
  digest::digest(ast_text(expr), algo = "md5", serialize = FALSE)
}

selector_text <- function(expr) {
  if (is.null(expr)) {
    return(NA_character_)
  }
  if (is.character(expr)) {
    return(paste(expr, collapse = ","))
  }
  ast_text(expr)
}

named_call_arg <- function(expr, arg_name, fallback_position = NULL) {
  args <- as.list(expr)[-1L]
  arg_names <- names(args) %||% rep("", length(args))
  named_index <- match(arg_name, arg_names)
  if (!is.na(named_index)) {
    return(args[[named_index]])
  }
  if (!is.null(fallback_position) && length(args) >= fallback_position) {
    return(args[[fallback_position]])
  }
  NULL
}

function_arg_index <- function(expr, named_arg) {
  args <- as.list(expr)[-1L]
  arg_names <- names(args) %||% rep("", length(args))
  named_index <- match(named_arg, arg_names)
  if (!is.na(named_index)) {
    return(named_index)
  }

  function_indices <- which(vapply(args, function(arg) {
    is.call(arg) && identical(call_name(arg), "function")
  }, logical(1L)))
  if (length(function_indices)) tail(function_indices, 1L) else NA_integer_
}

set_method_signature <- function(expr) {
  signature <- named_call_arg(expr, "signature")
  if (!is.null(signature)) {
    return(signature)
  }

  args <- as.list(expr)[-1L]
  if (length(args) < 2L || identical(function_arg_index(expr, "definition"), 2L)) {
    return(NULL)
  }
  args[[2L]]
}

assignment_name <- function(expr) {
  if (!is.call(expr) || !isTRUE(call_name(expr) %in% c("<-", "="))) {
    return(NULL)
  }
  lhs <- expr[[2L]]
  if (is.symbol(lhs)) as.character(lhs) else NULL
}

source_line <- function(expr) {
  srcref <- attr(expr, "srcref")
  if (is.null(srcref) || length(srcref) < 1L) NA_integer_ else as.integer(srcref[[1L]])
}

top_owner <- function(expr, expression_index) {
  symbol <- assignment_name(expr)
  if (!is.null(symbol)) {
    return(paste0("symbol::", symbol))
  }

  kind <- call_name(expr)
  if (isTRUE(kind %in% c("setMethod", "setGeneric", "setClass"))) {
    entity <- selector_text(named_call_arg(expr, if (kind == "setMethod") "f" else "name", 1L))
    signature <- if (kind == "setMethod") {
      selector_text(set_method_signature(expr))
    } else {
      NA_character_
    }
    return(paste(c(kind, entity, signature), collapse = "::"))
  }

  paste0("expression::", expression_index)
}

inventory_file <- function(path, relative_path, side) {
  expressions <- parse(file = path, keep.source = TRUE)
  expression_srcrefs <- attr(expressions, "srcref")
  records <- list()
  current_top_line <- NA_integer_

  add_function <- function(expr, kind, name, owner, context, function_depth) {
    formals_expr <- expr[[2L]]
    body_expr <- expr[[3L]]
    records[[length(records) + 1L]] <<- data.frame(
      side = side,
      file_path = relative_path,
      line = source_line(expr),
      top_expression_line = current_top_line,
      kind = kind,
      name = name %||% "",
      owner = owner,
      context = context,
      function_depth = function_depth,
      full_hash = ast_hash(expr),
      formals_hash = ast_hash(formals_expr),
      body_hash = ast_hash(body_expr),
      formals = ast_text(formals_expr),
      stringsAsFactors = FALSE
    )
  }

  walk <- function(node, owner, context, function_depth = 0L, hint_kind = NULL, hint_name = NULL) {
    if (!is.call(node)) {
      return(invisible(NULL))
    }

    node_call <- call_name(node)
    if (identical(node_call, "function")) {
      add_function(
        node,
        hint_kind %||% "anonymous",
        hint_name %||% "",
        owner,
        context,
        function_depth
      )

      formal_values <- as.list(node[[2L]])
      if (length(formal_values)) {
        formal_names <- names(formal_values) %||% rep("", length(formal_values))
        for (i in seq_along(formal_values)) {
          value <- formal_values[[i]]
          if (!rlang::is_missing(value) && is.language(value)) {
            walk(
              value,
              owner,
              paste0(context, "/formal:", formal_names[[i]]),
              function_depth + 1L
            )
          }
        }
      }
      walk(node[[3L]], owner, paste0(context, "/body"), function_depth + 1L)
      return(invisible(NULL))
    }

    assigned_name <- assignment_name(node)
    if (!is.null(assigned_name) && length(node) >= 3L &&
        is.call(node[[3L]]) && identical(call_name(node[[3L]]), "function")) {
      assigned_kind <- if (function_depth == 0L) "top_level_symbol" else "nested_assignment"
      walk(
        node[[3L]],
        owner,
        paste0(context, "/assignment:", assigned_name),
        function_depth,
        assigned_kind,
        assigned_name
      )
      return(invisible(NULL))
    }

    children <- as.list(node)[-1L]
    child_names <- names(children) %||% rep("", length(children))
    for (i in seq_along(children)) {
      if (rlang::is_missing(children[[i]])) {
        next
      }
      child <- children[[i]]
      if (!is.language(child)) {
        next
      }

      label <- child_names[[i]]
      if (is.na(label) || !nzchar(label)) {
        label <- as.character(i)
      }
      child_kind <- NULL
      child_name <- NULL

      if (is.call(child) && identical(call_name(child), "function")) {
        if (identical(node_call, "setMethod") &&
            identical(i, function_arg_index(node, "definition"))) {
          child_kind <- "setMethod_definition"
          child_name <- selector_text(named_call_arg(node, "f", 1L))
        } else if (identical(node_call, "setGeneric") &&
                   identical(i, function_arg_index(node, "def"))) {
          child_kind <- "setGeneric_definition"
          child_name <- selector_text(named_call_arg(node, "name", 1L))
        } else {
          child_kind <- "call_argument"
          child_name <- paste0(node_call %||% "call", "::", label)
        }
      }

      walk(
        child,
        owner,
        paste0(context, "/", node_call %||% "call", ":", label),
        function_depth,
        child_kind,
        child_name
      )
    }

    invisible(NULL)
  }

  for (i in seq_along(expressions)) {
    current_top_line <- if (length(expression_srcrefs) >= i) {
      as.integer(expression_srcrefs[[i]][[1L]])
    } else {
      NA_integer_
    }
    owner <- top_owner(expressions[[i]], i)
    walk(expressions[[i]], owner, paste0("expr:", i))
  }

  if (!length(records)) {
    return(NULL)
  }
  do.call(rbind, records)
}

inventory_root <- function(root, side) {
  r_root <- file.path(root, "R")
  files <- sort(list.files(r_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE))
  records <- lapply(files, function(path) {
    relative_path <- file.path("R", substring(path, nchar(r_root) + 2L))
    inventory_file(path, relative_path, side)
  })
  records <- Filter(Negate(is.null), records)
  if (!length(records)) {
    return(data.frame())
  }
  result <- do.call(rbind, records)
  rownames(result) <- NULL
  result$occurrence_id <- paste0(side, "-function-", seq_len(nrow(result)))
  result
}

inventory_todo_scaffolds <- function(root, baseline_inventory, target_inventory) {
  r_root <- file.path(root, "R")
  files <- sort(list.files(r_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE))
  records <- list()
  function_pattern <- "^\\s*#\\s*Function\\s+[0-9]+:\\s*([.A-Za-z][.A-Za-z0-9_]*)\\s*\\(\\)"

  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    todo_lines <- grep("^\\s*#\\s*TODO:\\s*Extract", lines)
    if (!length(todo_lines)) {
      next
    }

    todo_line <- min(todo_lines)
    if (todo_line >= length(lines)) {
      next
    }

    relative_path <- file.path("R", substring(path, nchar(r_root) + 2L))
    block_lines <- seq.int(todo_line + 1L, length(lines))
    matches <- regexec(function_pattern, lines[block_lines], perl = TRUE)
    captures <- regmatches(lines[block_lines], matches)
    matched <- lengths(captures) >= 2L
    if (!any(matched)) {
      next
    }

    matched_lines <- block_lines[matched]
    matched_names <- vapply(captures[matched], `[[`, character(1L), 2L)
    records[[length(records) + 1L]] <- data.frame(
      todo_file = relative_path,
      todo_line = matched_lines,
      name = matched_names,
      stringsAsFactors = FALSE
    )
  }

  if (!length(records)) {
    return(data.frame())
  }

  result <- do.call(rbind, records)
  baseline_counts <- table(baseline_inventory$name[nzchar(baseline_inventory$name)])
  target_counts <- table(target_inventory$name[nzchar(target_inventory$name)])
  result$baseline_named_occurrences <- unname(baseline_counts[result$name])
  result$target_named_occurrences <- unname(target_counts[result$name])
  result$baseline_named_occurrences[is.na(result$baseline_named_occurrences)] <- 0L
  result$target_named_occurrences[is.na(result$target_named_occurrences)] <- 0L
  rownames(result) <- NULL
  result
}

extract_git_ref <- function(repo_root, ref) {
  snapshot_root <- tempfile("function-parity-baseline-")
  dir.create(snapshot_root, recursive = TRUE)
  archive_path <- tempfile("function-parity-baseline-", fileext = ".tar")
  on.exit(unlink(archive_path), add = TRUE)

  output <- system2(
    "git",
    c("-C", repo_root, "archive", "--format=tar", paste0("--output=", archive_path), ref),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status") %||% 0L
  if (status != 0L) {
    abort("Could not archive baseline ref ", ref, ": ", paste(output, collapse = "\n"))
  }
  utils::untar(archive_path, exdir = snapshot_root)
  snapshot_root
}

count_by_hash <- function(values) {
  counts <- table(values)
  setNames(as.integer(counts), names(counts))
}

compare_inventory <- function(baseline, target) {
  target_full <- count_by_hash(target$full_hash)
  target_body <- count_by_hash(target$body_hash)
  baseline_full <- count_by_hash(baseline$full_hash)

  identity_key <- function(data) paste(data$kind, data$name, data$owner, sep = "::")
  baseline$identity_key <- identity_key(baseline)
  target$identity_key <- identity_key(target)
  target_identities <- unique(target$identity_key)

  target_location <- paste(
    paste0(target$file_path, ":", target$top_expression_line),
    target$kind,
    ifelse(nzchar(target$name), target$name, "<anonymous>"),
    target$owner,
    target$context,
    sep = " | "
  )
  target_by_full <- split(seq_len(nrow(target)), target$full_hash)
  target_by_body <- split(seq_len(nrow(target)), target$body_hash)
  target_by_identity <- split(seq_len(nrow(target)), target$identity_key)

  baseline$baseline_full_multiplicity <- unname(baseline_full[baseline$full_hash])
  baseline$target_full_multiplicity <- unname(target_full[baseline$full_hash])
  baseline$target_full_multiplicity[is.na(baseline$target_full_multiplicity)] <- 0L
  baseline$target_body_multiplicity <- unname(target_body[baseline$body_hash])
  baseline$target_body_multiplicity[is.na(baseline$target_body_multiplicity)] <- 0L

  baseline$status <- ifelse(
    baseline$target_full_multiplicity > 0L,
    "exact_function_ast_present",
    ifelse(
      baseline$target_body_multiplicity > 0L,
      "exact_body_ast_present_formals_drift",
      ifelse(
        baseline$identity_key %in% target_identities,
        "named_or_owned_function_drift",
        "unmatched_function_expression"
      )
    )
  )
  baseline$match_basis <- ifelse(
    baseline$status == "exact_function_ast_present",
    "full_ast",
    ifelse(
      baseline$status == "exact_body_ast_present_formals_drift",
      "body_ast",
      ifelse(
        baseline$status == "named_or_owned_function_drift",
        "identity",
        "none"
      )
    )
  )
  matching_indices <- lapply(seq_len(nrow(baseline)), function(i) {
    switch(
      baseline$match_basis[[i]],
      full_ast = target_by_full[[baseline$full_hash[[i]]]],
      body_ast = target_by_body[[baseline$body_hash[[i]]]],
      identity = target_by_identity[[baseline$identity_key[[i]]]],
      integer()
    ) %||% integer()
  })
  baseline$target_match_count <- vapply(matching_indices, length, integer(1L))
  baseline$target_match_locations <- vapply(
    matching_indices,
    function(indices) paste(target_location[indices], collapse = " || "),
    character(1L)
  )
  baseline$multiplicity_reduced <- baseline$target_full_multiplicity > 0L &
    baseline$target_full_multiplicity < baseline$baseline_full_multiplicity
  baseline
}

main <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) {
    abort("The digest package is required")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    abort("The jsonlite package is required")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    abort("The rlang package is required")
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(args$repo_root, mustWork = TRUE)
  target_root <- normalizePath(args$target_root, mustWork = TRUE)
  output_dir <- if (grepl("^/", args$output_dir)) {
    args$output_dir
  } else {
    file.path(repo_root, args$output_dir)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  baseline_root <- extract_git_ref(repo_root, args$baseline_ref)
  on.exit(unlink(baseline_root, recursive = TRUE, force = TRUE), add = TRUE)

  baseline <- inventory_root(baseline_root, "baseline")
  target <- inventory_root(target_root, "target")
  comparison <- compare_inventory(baseline, target)
  todo_scaffolds <- inventory_todo_scaffolds(target_root, baseline, target)

  unique_status <- aggregate(
    occurrence_id ~ full_hash + status,
    data = comparison,
    FUN = length
  )
  status_counts <- sort(table(comparison$status), decreasing = TRUE)
  summary <- list(
    baseline_ref = args$baseline_ref,
    target_root = target_root,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    baseline = list(
      file_count = length(unique(baseline$file_path)),
      function_occurrence_count = nrow(baseline),
      unique_function_ast_count = length(unique(baseline$full_hash)),
      unique_body_ast_count = length(unique(baseline$body_hash))
    ),
    target = list(
      file_count = length(unique(target$file_path)),
      function_occurrence_count = nrow(target),
      unique_function_ast_count = length(unique(target$full_hash)),
      unique_body_ast_count = length(unique(target$body_hash))
    ),
    todo_scaffolds = list(
      file_count = length(unique(todo_scaffolds$todo_file)),
      entry_count = nrow(todo_scaffolds),
      unique_name_count = length(unique(todo_scaffolds$name)),
      target_resolved_entry_count = sum(todo_scaffolds$target_named_occurrences > 0L),
      target_unresolved_entry_count = sum(todo_scaffolds$target_named_occurrences == 0L)
    ),
    comparison = list(
      occurrence_status_counts = as.list(status_counts),
      unique_baseline_function_asts_exactly_present = length(unique(
        comparison$full_hash[comparison$status == "exact_function_ast_present"]
      )),
      unique_baseline_function_asts_not_exactly_present = length(unique(
        comparison$full_hash[comparison$status != "exact_function_ast_present"]
      )),
      unmatched_occurrence_count = sum(comparison$status == "unmatched_function_expression"),
      multiplicity_reduction_occurrence_count = sum(comparison$multiplicity_reduced)
    )
  )

  utils::write.csv(
    comparison,
    file.path(output_dir, "baseline-function-comparison.csv"),
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    target,
    file.path(output_dir, "target-function-inventory.csv"),
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    unique_status,
    file.path(output_dir, "unique-function-status.csv"),
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    todo_scaffolds,
    file.path(output_dir, "todo-function-inventory.csv"),
    row.names = FALSE,
    na = ""
  )
  jsonlite::write_json(
    summary,
    file.path(output_dir, "summary.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  cat(jsonlite::toJSON(summary, auto_unbox = TRUE), "\n")
}

main()
