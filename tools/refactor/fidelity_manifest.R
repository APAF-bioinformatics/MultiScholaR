suppressPackageStartupMessages({
  library(yaml)
})

normalize_manifest <- function(manifest) {
  if (is.null(manifest$version)) {
    abort("Manifest is missing 'version'")
  }

  if (is.null(manifest$entries) || !length(manifest$entries)) {
    abort("Manifest is missing 'entries'")
  }

  manifest
}

read_manifest <- function(path) {
  normalize_manifest(yaml::read_yaml(path))
}

find_leading_comment_start <- function(lines, expr_line1) {
  i <- expr_line1 - 1L
  if (i < 1L) {
    return(expr_line1)
  }

  saw_comment <- FALSE
  while (i >= 1L && grepl("^\\s*#", lines[[i]])) {
    saw_comment <- TRUE
    i <- i - 1L
  }

  if (saw_comment) {
    return(i + 1L)
  }

  expr_line1
}

match_selector_records <- function(selector, inventory, entry_id = NULL) {
  selector_kind <- selector$kind
  selector_signature <- if (!is.null(selector$signature)) normalize_selector_value(selector$signature) else NULL

  matches <- switch(
    selector_kind,
    symbol = Filter(function(rec) identical(rec$symbol, selector$value), inventory),
    setMethod = Filter(function(rec) {
      identical(rec$call_kind, "setMethod") &&
        identical(rec$call_value, selector$value) &&
        (is.null(selector_signature) || identical(rec$call_signature, selector_signature))
    }, inventory),
    setGeneric = Filter(function(rec) identical(rec$call_kind, "setGeneric") && identical(rec$call_value, selector$value), inventory),
    setClass = Filter(function(rec) identical(rec$call_kind, "setClass") && identical(rec$call_value, selector$value), inventory),
    expr_index = {
      idx <- as.integer(selector$value)
      if (is.na(idx) || idx < 1L || idx > length(inventory)) {
        abort("Invalid expr_index", if (!is.null(entry_id)) paste0(" for entry ", entry_id) else "", ": ", selector$value)
      }
      inventory[idx]
    },
    NULL
  )

  if (is.null(matches)) {
    abort("Unsupported selector kind", if (!is.null(entry_id)) paste0(" for entry ", entry_id) else "", ": ", selector_kind)
  }

  matches
}

resolve_expr_selector <- function(entry, selector, lines, inventory) {
  matches <- match_selector_records(selector, inventory, entry$id)

  if (!length(matches)) {
    abort("Selector did not match any block for entry ", entry$id)
  }

  if (length(matches) != 1L) {
    abort("Selector matched ", length(matches), " blocks for entry ", entry$id)
  }

  rec <- matches[[1]]
  sr <- rec$srcref
  start_line <- if (isTRUE(selector$include_leading_comments %||% TRUE)) {
    find_leading_comment_start(lines, sr[[1]])
  } else {
    sr[[1]]
  }
  end_line <- sr[[3]]

  list(
    id = entry$id,
    action = entry$action %||% "extract",
    source = entry$source,
    target = entry$target %||% NA_character_,
    group = entry$group %||% NA_character_,
    selector_kind = selector$kind,
    selector_value = selector$value,
    start_line = start_line,
    end_line = end_line,
    text = normalize_block_text(lines, start_line, end_line)
  )
}

resolve_anchor_range <- function(entry, selector, lines) {
  start_hits <- which(grepl(selector$start, lines, perl = TRUE))
  if (!length(start_hits)) {
    abort("Anchor start did not match for entry ", entry$id)
  }

  occurrence <- as.integer(selector$occurrence %||% 1L)
  if (occurrence < 1L || occurrence > length(start_hits)) {
    abort("Invalid occurrence for entry ", entry$id, ": ", occurrence)
  }

  start_line <- start_hits[[occurrence]]
  end_hits <- which(seq_along(lines) > start_line & grepl(selector$end, lines, perl = TRUE))
  if (!length(end_hits)) {
    abort("Anchor end did not match for entry ", entry$id)
  }

  end_anchor <- end_hits[[1]]
  end_line <- if (isTRUE(selector$end_inclusive %||% FALSE)) end_anchor else end_anchor - 1L

  if (end_line < start_line) {
    abort("Anchor range is empty for entry ", entry$id)
  }

  list(
    id = entry$id,
    action = entry$action %||% "extract",
    source = entry$source,
    target = entry$target %||% NA_character_,
    group = entry$group %||% NA_character_,
    selector_kind = selector$kind,
    selector_value = paste(selector$start, "->", selector$end),
    start_line = start_line,
    end_line = end_line,
    text = normalize_block_text(lines, start_line, end_line)
  )
}

resolve_entry <- function(entry, repo_root, cache) {
  source_path <- file.path(repo_root, entry$source)
  if (!file.exists(source_path)) {
    abort("Source file does not exist: ", entry$source)
  }

  if (is.null(cache[[entry$source]])) {
    cache[[entry$source]] <- list(
      lines = read_source_lines(source_path),
      inventory = build_inventory(source_path)
    )
  }

  selector <- entry$selector
  if (is.null(selector$kind)) {
    abort("Entry is missing selector.kind: ", entry$id)
  }

  if (identical(selector$kind, "anchor_range")) {
    block <- resolve_anchor_range(entry, selector, cache[[entry$source]]$lines)
  } else {
    block <- resolve_expr_selector(entry, selector, cache[[entry$source]]$lines, cache[[entry$source]]$inventory)
  }

  list(block = block, cache = cache)
}

validate_selector <- function(entry, inventory, lines) {
  selector <- entry$selector
  kind <- selector$kind

  if (identical(kind, "anchor_range")) {
    start_hits <- which(grepl(selector$start, lines, perl = TRUE))
    end_hits <- which(grepl(selector$end, lines, perl = TRUE))
    if (!length(start_hits)) {
      return("anchor start did not match")
    }
    if (!length(end_hits)) {
      return("anchor end did not match")
    }
    return(NULL)
  }

  matches <- tryCatch(
    match_selector_records(selector, inventory, entry$id),
    error = function(e) e
  )
  if (inherits(matches, "error")) {
    return(conditionMessage(matches))
  }

  if (!length(matches)) {
    return("selector matched 0 blocks")
  }

  if (length(matches) != 1L) {
    return(sprintf("selector matched %d blocks", length(matches)))
  }

  NULL
}

build_block_from_record <- function(entry, selector, lines, rec, source_path, target_path = entry$target %||% NA_character_) {
  sr <- rec$srcref
  start_line <- if (isTRUE(selector$include_leading_comments %||% TRUE)) {
    find_leading_comment_start(lines, sr[[1]])
  } else {
    sr[[1]]
  }
  end_line <- sr[[3]]

  list(
    id = entry$id,
    action = entry$action %||% "extract",
    source = source_path,
    target = target_path,
    group = entry$group %||% NA_character_,
    selector_kind = selector$kind,
    selector_value = selector$value %||% NA_character_,
    start_line = start_line,
    end_line = end_line,
    text = normalize_block_text(lines, start_line, end_line)
  )
}

block_ast_hash <- function(text) {
  if (is.null(text) || !nzchar(text)) {
    return(NULL)
  }

  exprs <- tryCatch(
    parse(text = text, keep.source = FALSE),
    error = function(e) NULL
  )
  if (is.null(exprs)) {
    return(NULL)
  }

  hash_text(paste(vapply(as.list(exprs), ast_fingerprint, character(1)), collapse = "\n\n"))
}

strip_comment_lines <- function(text) {
  if (is.null(text) || !nzchar(text)) {
    return("")
  }

  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- lines[!grepl("^\\s*#", lines)]
  paste(lines, collapse = "\n")
}

augment_block_fidelity <- function(block) {
  text <- block$text %||% ""
  block$hash_raw <- hash_text(text)
  block$normalized_text <- normalize_text_for_fidelity(text)
  block$hash_normalized <- hash_text(block$normalized_text)
  block$hash_ast <- block_ast_hash(text)
  block$noncomment_normalized <- normalize_text_for_fidelity(strip_comment_lines(text))
  block
}

resolve_target_cache <- function(target_root, target_path, cache) {
  if (!nzchar(target_path)) {
    return(list(cache = cache, state = NULL))
  }

  if (is.null(cache[[target_path]])) {
    absolute_path <- file.path(target_root, target_path)
    if (!file.exists(absolute_path)) {
      return(list(
        cache = cache,
        state = list(
          absolute_path = absolute_path,
          exists = FALSE
        )
      ))
    }

    inventory_result <- tryCatch(
      build_inventory(absolute_path),
      error = function(e) e
    )
    cache[[target_path]] <- list(
      absolute_path = absolute_path,
      exists = TRUE,
      lines = read_source_lines(absolute_path),
      inventory = if (inherits(inventory_result, "error")) NULL else inventory_result,
      inventory_error = if (inherits(inventory_result, "error")) conditionMessage(inventory_result) else NULL
    )
  }

  list(cache = cache, state = cache[[target_path]])
}

build_target_candidates <- function(entry, target_path, lines, inventory) {
  if (is.null(inventory) || !length(inventory)) {
    return(list())
  }

  lapply(seq_along(inventory), function(i) {
    rec <- inventory[[i]]
    candidate <- build_block_from_record(entry, entry$selector, lines, rec, target_path, target_path)
    candidate$candidate_index <- i
    augment_block_fidelity(candidate)
  })
}

build_file_block <- function(entry, target_path, lines) {
  block <- list(
    id = entry$id,
    action = entry$action %||% "extract",
    source = entry$source,
    target = target_path,
    group = entry$group %||% NA_character_,
    selector_kind = entry$selector$kind,
    selector_value = entry$selector$value %||% NA_character_,
    start_line = if (length(lines)) 1L else NA_integer_,
    end_line = if (length(lines)) length(lines) else NA_integer_,
    text = if (length(lines)) paste0(paste(lines, collapse = "\n"), "\n") else ""
  )
  augment_block_fidelity(block)
}

resolve_target_block <- function(entry, source_block, target_root, cache) {
  target_path <- entry$target %||% ""
  if (!nzchar(target_path)) {
    return(list(
      cache = cache,
      status = "content_missing",
      target_resolver = "missing_target_path",
      note = "manifest entry has no target path",
      block = NULL
    ))
  }

  cached <- resolve_target_cache(target_root, target_path, cache)
  cache <- cached$cache
  target_state <- cached$state

  if (is.null(target_state) || !isTRUE(target_state$exists)) {
    return(list(
      cache = cache,
      status = "content_missing",
      target_resolver = "missing_target_file",
      note = "target file does not exist",
      block = NULL
    ))
  }

  lines <- target_state$lines
  inventory <- target_state$inventory
  candidates <- build_target_candidates(entry, target_path, lines, inventory)
  file_block <- build_file_block(entry, target_path, lines)
  selector <- entry$selector
  selector_status <- NULL

  if (!is.null(inventory) && selector$kind %in% c("symbol", "setMethod", "setGeneric", "setClass")) {
    selector_matches <- tryCatch(
      match_selector_records(selector, inventory, entry$id),
      error = function(e) e
    )
    if (!inherits(selector_matches, "error")) {
      if (length(selector_matches) == 1L) {
        block <- build_block_from_record(entry, selector, lines, selector_matches[[1]], target_path, target_path)
        return(list(
          cache = cache,
          status = "resolved",
          target_resolver = "selector",
          note = NULL,
          block = augment_block_fidelity(block)
        ))
      }
      if (length(selector_matches) > 1L) {
        selector_status <- "selector_ambiguous"
      }
    } else {
      selector_status <- conditionMessage(selector_matches)
    }
  }

  raw_matches <- Filter(function(candidate) identical(candidate$hash_raw, source_block$hash_raw), candidates)
  if (length(raw_matches) == 1L) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "raw_text_search",
      note = NULL,
      block = raw_matches[[1]]
    ))
  }

  if (identical(file_block$hash_raw, source_block$hash_raw)) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "file_raw_text_search",
      note = NULL,
      block = file_block
    ))
  }

  normalized_matches <- Filter(function(candidate) identical(candidate$hash_normalized, source_block$hash_normalized), candidates)
  if (length(normalized_matches) == 1L) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "normalized_text_search",
      note = NULL,
      block = normalized_matches[[1]]
    ))
  }

  if (identical(file_block$hash_normalized, source_block$hash_normalized)) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "file_normalized_text_search",
      note = NULL,
      block = file_block
    ))
  }

  ast_matches <- Filter(function(candidate) {
    !is.null(source_block$hash_ast) &&
      !is.null(candidate$hash_ast) &&
      identical(candidate$hash_ast, source_block$hash_ast)
  }, candidates)
  if (length(ast_matches) == 1L) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "ast_search",
      note = NULL,
      block = ast_matches[[1]]
    ))
  }

  if (!is.null(source_block$hash_ast) && !is.null(file_block$hash_ast) && identical(file_block$hash_ast, source_block$hash_ast)) {
    return(list(
      cache = cache,
      status = "resolved",
      target_resolver = "file_ast_search",
      note = NULL,
      block = file_block
    ))
  }

  if (identical(selector_status, "selector_ambiguous")) {
    return(list(
      cache = cache,
      status = "selector_ambiguous",
      target_resolver = "selector_ambiguous",
      note = "selector matched multiple target blocks",
      block = NULL
    ))
  }

  note <- target_state$inventory_error
  if (is.null(note) && is.character(selector_status) && nzchar(selector_status)) {
    note <- selector_status
  }
  if (is.null(note)) {
    note <- "no matching target block found"
  }

  list(
    cache = cache,
    status = "content_missing",
    target_resolver = "content_missing",
    note = note,
    block = NULL
  )
}

classify_manifest_result <- function(entry, source_block, target_result) {
  if (identical(entry$action %||% "extract", "manual_merge")) {
    return(list(
      status = "manual_merge_expected",
      raw_exact = FALSE,
      normalized_text_exact = FALSE,
      ast_exact = FALSE,
      exception_candidate = "manual_merge_expected"
    ))
  }

  if (!identical(target_result$status, "resolved") || is.null(target_result$block)) {
    exception_type <- if (identical(target_result$status, "selector_ambiguous")) {
      "selector_ambiguous"
    } else {
      "missing_target_block"
    }
    return(list(
      status = target_result$status,
      raw_exact = FALSE,
      normalized_text_exact = FALSE,
      ast_exact = FALSE,
      exception_candidate = exception_type
    ))
  }

  target_block <- target_result$block
  raw_exact <- identical(source_block$hash_raw, target_block$hash_raw)
  normalized_text_exact <- identical(source_block$hash_normalized, target_block$hash_normalized)
  ast_exact <- !is.null(source_block$hash_ast) &&
    !is.null(target_block$hash_ast) &&
    identical(source_block$hash_ast, target_block$hash_ast)

  exception_candidate <- NULL
  status <- "semantic_drift_suspected"

  if (raw_exact) {
    status <- "raw_exact"
  } else if (normalized_text_exact) {
    status <- "normalized_only"
    exception_candidate <- "whitespace_only_drift"
  } else if (ast_exact) {
    status <- "ast_only"
    if (identical(source_block$noncomment_normalized, target_block$noncomment_normalized)) {
      exception_candidate <- "comment_only_drift"
    } else {
      exception_candidate <- "ast_only_drift"
    }
  } else {
    exception_candidate <- "semantic_mismatch"
  }

  list(
    status = status,
    raw_exact = raw_exact,
    normalized_text_exact = normalized_text_exact,
    ast_exact = ast_exact,
    exception_candidate = exception_candidate
  )
}

manifest_entry_catalog <- function(entry) {
  list(
    entry_id = entry$id,
    action = entry$action %||% "extract",
    selector_kind = entry$selector$kind %||% NA_character_,
    selector_payload = entry$selector,
    source_path = entry$source %||% NA_character_,
    target_path = entry$target %||% NA_character_,
    group_name = entry$group %||% NA_character_
  )
}

compare_manifest_entry <- function(entry, baseline_root, target_root, baseline_cache, target_cache) {
  source_result <- resolve_entry(entry, baseline_root, baseline_cache)
  baseline_cache <- source_result$cache
  source_block <- augment_block_fidelity(source_result$block)
  target_result <- resolve_target_block(entry, source_block, target_root, target_cache)
  target_cache <- target_result$cache
  classification <- classify_manifest_result(entry, source_block, target_result)
  target_block <- target_result$block

  list(
    baseline_cache = baseline_cache,
    target_cache = target_cache,
    comparison = c(
      manifest_entry_catalog(entry),
      list(
        baseline_source_resolver = paste0("selector:", entry$selector$kind),
        target_resolver = target_result$target_resolver %||% NA_character_,
        raw_exact = classification$raw_exact,
        normalized_text_exact = classification$normalized_text_exact,
        ast_exact = classification$ast_exact,
        hash_raw_baseline = source_block$hash_raw %||% NA_character_,
        hash_raw_target = if (!is.null(target_block)) target_block$hash_raw %||% NA_character_ else NA_character_,
        hash_normalized_baseline = source_block$hash_normalized %||% NA_character_,
        hash_normalized_target = if (!is.null(target_block)) target_block$hash_normalized %||% NA_character_ else NA_character_,
        hash_ast_baseline = source_block$hash_ast %||% NA_character_,
        hash_ast_target = if (!is.null(target_block)) target_block$hash_ast %||% NA_character_ else NA_character_,
        status = classification$status,
        exception_candidate = classification$exception_candidate %||% NA_character_,
        note = target_result$note %||% NA_character_,
        source_start_line = source_block$start_line %||% NA_integer_,
        source_end_line = source_block$end_line %||% NA_integer_,
        target_start_line = if (!is.null(target_block)) target_block$start_line %||% NA_integer_ else NA_integer_,
        target_end_line = if (!is.null(target_block)) target_block$end_line %||% NA_integer_ else NA_integer_
      )
    )
  )
}
