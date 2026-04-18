`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

abort <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

trim_quotes <- function(x) {
  gsub('^"|"$', "", x)
}

normalize_call_name <- function(expr) {
  if (is.null(expr) || !is.call(expr)) {
    return(NULL)
  }

  head <- expr[[1]]
  if (is.symbol(head)) {
    return(as.character(head))
  }

  NULL
}

normalize_selector_value <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.character(x)) {
    if (length(x) == 1L) {
      return(trim_quotes(x[[1]]))
    }

    quoted <- sprintf('"%s"', trim_quotes(x))
    return(sprintf("c(%s)", paste(quoted, collapse = ", ")))
  }

  if (is.symbol(x)) {
    return(as.character(x))
  }

  trim_quotes(paste(deparse(x, width.cutoff = 500L), collapse = ""))
}

is_top_level_assignment <- function(expr) {
  is.call(expr) && length(expr) >= 3 && identical(normalize_call_name(expr), "<-")
}

is_top_level_equals_assignment <- function(expr) {
  is.call(expr) && length(expr) >= 3 && identical(normalize_call_name(expr), "=")
}

extract_top_level_symbol <- function(expr) {
  if (is_top_level_assignment(expr) || is_top_level_equals_assignment(expr)) {
    lhs <- expr[[2]]
    if (is.symbol(lhs)) {
      return(as.character(lhs))
    }
  }

  NULL
}

extract_method_signature <- function(expr) {
  if (!is.call(expr) || !identical(normalize_call_name(expr), "setMethod")) {
    return(NULL)
  }

  args <- as.list(expr)[-1]
  arg_names <- names(args) %||% rep("", length(args))
  named_idx <- match("signature", arg_names)

  signature_expr <- NULL
  if (!is.na(named_idx)) {
    signature_expr <- args[[named_idx]]
  } else if (length(args) >= 2) {
    signature_expr <- args[[2]]
  }

  normalize_selector_value(signature_expr)
}

extract_call_match <- function(expr) {
  call_name <- normalize_call_name(expr)

  if (is.null(call_name) || !call_name %in% c("setMethod", "setGeneric", "setClass")) {
    return(list(kind = NULL, value = NULL, signature = NULL))
  }

  value <- if (length(expr) >= 2) normalize_selector_value(expr[[2]]) else NULL
  signature <- if (identical(call_name, "setMethod")) extract_method_signature(expr) else NULL
  list(kind = call_name, value = value, signature = signature)
}

read_source_lines <- function(path) {
  readLines(path, warn = FALSE)
}

line_from_srcref <- function(srcref, index) {
  if (is.null(srcref) || length(srcref) < index) {
    return(NA_integer_)
  }

  as.integer(srcref[[index]])
}

normalize_block_text <- function(lines, start_line, end_line) {
  if (!length(lines) || is.na(start_line) || is.na(end_line) || start_line < 1L || end_line < start_line) {
    return(NULL)
  }

  block <- lines[start_line:end_line]

  while (length(block) > 0 && identical(block[[length(block)]], "")) {
    block <- block[-length(block)]
  }

  if (!length(block)) {
    return("")
  }

  paste0(paste(block, collapse = "\n"), "\n")
}

normalize_text_for_fidelity <- function(text) {
  if (is.null(text) || identical(text, "")) {
    return("")
  }

  gsub("\\s+", " ", trimws(text))
}

hash_text <- function(text) {
  if (is.null(text)) {
    return(NULL)
  }

  temp_path <- tempfile("fidelity-hash-")
  on.exit(unlink(temp_path), add = TRUE)
  writeChar(text, temp_path, eos = NULL, useBytes = TRUE)
  unname(tools::md5sum(temp_path)[[1]])
}

ast_fingerprint <- function(expr) {
  paste(deparse(expr, width.cutoff = 500L), collapse = "\n")
}

build_inventory <- function(path) {
  exprs <- parse(file = path, keep.source = TRUE)
  srcrefs <- attr(exprs, "srcref") %||% vector("list", length(exprs))

  lapply(seq_along(exprs), function(i) {
    expr <- exprs[[i]]
    call_match <- extract_call_match(expr)

    list(
      index = i,
      expr = expr,
      srcref = srcrefs[[i]],
      symbol = extract_top_level_symbol(expr),
      call_kind = call_match$kind,
      call_value = call_match$value,
      call_signature = call_match$signature
    )
  })
}

list_r_files <- function(repo_root) {
  r_dir <- file.path(repo_root, "R")
  if (!dir.exists(r_dir)) {
    return(character())
  }

  files <- list.files(r_dir, pattern = "\\.R$", recursive = TRUE, full.names = FALSE)
  sort(file.path("R", files))
}

read_description_collate <- function(repo_root) {
  description_path <- file.path(repo_root, "DESCRIPTION")
  if (!file.exists(description_path)) {
    return(character())
  }

  dcf <- read.dcf(description_path)
  if (!"Collate" %in% colnames(dcf)) {
    return(character())
  }

  collate_text <- dcf[1, "Collate"]
  if (is.na(collate_text) || !nzchar(trimws(collate_text))) {
    return(character())
  }

  scan(text = collate_text, what = character(), quiet = TRUE)
}

read_namespace_exports <- function(repo_root) {
  namespace_path <- file.path(repo_root, "NAMESPACE")
  exports <- list(
    functions = character(),
    methods = character(),
    classes = character()
  )

  if (!file.exists(namespace_path)) {
    return(exports)
  }

  lines <- readLines(namespace_path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  lines <- lines[nzchar(trimws(lines))]

  if (!length(lines)) {
    return(exports)
  }

  exprs <- parse(text = lines, keep.source = FALSE)

  for (expr in exprs) {
    call_name <- normalize_call_name(expr)
    if (is.null(call_name) || length(expr) < 2) {
      next
    }

    value <- normalize_selector_value(expr[[2]])
    if (is.null(value)) {
      next
    }

    if (identical(call_name, "export")) {
      exports$functions <- c(exports$functions, value)
    } else if (identical(call_name, "exportMethods")) {
      exports$methods <- c(exports$methods, value)
    } else if (identical(call_name, "exportClasses")) {
      exports$classes <- c(exports$classes, value)
    }
  }

  lapply(exports, function(values) sort(unique(values)))
}

make_entity_key <- function(entity_kind, entity_name, signature_key = NULL) {
  parts <- c(entity_kind, entity_name)
  if (!is.null(signature_key) && nzchar(signature_key)) {
    parts <- c(parts, signature_key)
  }
  paste(parts, collapse = "::")
}

is_entity_exported <- function(entity_kind, entity_name, exports) {
  if (identical(entity_kind, "setMethod")) {
    return(entity_name %in% exports$methods)
  }

  if (identical(entity_kind, "setClass")) {
    return(entity_name %in% exports$classes)
  }

  entity_name %in% exports$functions
}

entity_record <- function(rec, relative_path, lines, collate_index, exports) {
  start_line <- line_from_srcref(rec$srcref, 1)
  end_line <- line_from_srcref(rec$srcref, 3)
  raw_text <- normalize_block_text(lines, start_line, end_line)
  ast_text <- ast_fingerprint(rec$expr)

  records <- list()

  if (!is.null(rec$symbol)) {
    records[[length(records) + 1L]] <- list(
      entity_key = make_entity_key("symbol", rec$symbol),
      entity_kind = "symbol",
      entity_name = rec$symbol,
      signature_key = NULL,
      file_path = relative_path,
      line_start = start_line,
      line_end = end_line,
      exported = is_entity_exported("symbol", rec$symbol, exports),
      collate_index = collate_index,
      hash_raw = hash_text(raw_text),
      hash_normalized = hash_text(normalize_text_for_fidelity(raw_text)),
      hash_ast = hash_text(ast_text)
    )
  }

  if (!is.null(rec$call_kind) && !is.null(rec$call_value)) {
    records[[length(records) + 1L]] <- list(
      entity_key = make_entity_key(rec$call_kind, rec$call_value, rec$call_signature),
      entity_kind = rec$call_kind,
      entity_name = rec$call_value,
      signature_key = rec$call_signature,
      file_path = relative_path,
      line_start = start_line,
      line_end = end_line,
      exported = is_entity_exported(rec$call_kind, rec$call_value, exports),
      collate_index = collate_index,
      hash_raw = hash_text(raw_text),
      hash_normalized = hash_text(normalize_text_for_fidelity(raw_text)),
      hash_ast = hash_text(ast_text)
    )
  }

  records
}

inventory_snapshot <- function(repo_root) {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
  collate <- read_description_collate(repo_root)
  collate_map <- setNames(seq_along(collate), file.path("R", collate))
  exports <- read_namespace_exports(repo_root)
  r_files <- list_r_files(repo_root)
  entities <- list()
  parse_failures <- character()

  for (relative_path in r_files) {
    full_path <- file.path(repo_root, relative_path)
    file_inventory <- tryCatch(
      build_inventory(full_path),
      error = function(e) e
    )

    if (inherits(file_inventory, "error")) {
      parse_failures <- c(parse_failures, sprintf("%s: %s", relative_path, conditionMessage(file_inventory)))
      next
    }

    lines <- read_source_lines(full_path)
    collate_index <- unname(collate_map[[relative_path]] %||% NA_integer_)

    for (rec in file_inventory) {
      entities <- c(entities, entity_record(rec, relative_path, lines, collate_index, exports))
    }
  }

  entity_kind_counts <- table(vapply(entities, `[[`, character(1), "entity_kind"))

  list(
    repo_root = repo_root,
    collate = lapply(seq_along(collate), function(i) list(file_path = file.path("R", collate[[i]]), collate_index = i)),
    exports = list(
      functions = I(unname(exports$functions)),
      methods = I(unname(exports$methods)),
      classes = I(unname(exports$classes))
    ),
    r_files = I(unname(r_files)),
    entities = entities,
    parse_failures = I(unname(parse_failures)),
    counts = list(
      r_files = length(r_files),
      collate = length(collate),
      exports = list(
        functions = length(exports$functions),
        methods = length(exports$methods),
        classes = length(exports$classes)
      ),
      entities = as.list(stats::setNames(as.integer(entity_kind_counts), names(entity_kind_counts))),
      parse_failures = length(parse_failures)
    )
  )
}
