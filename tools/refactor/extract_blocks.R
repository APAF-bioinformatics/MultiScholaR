#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

abort <- function(...) {
  msg <- paste0(..., collapse = "")
  writeLines(msg, con = stderr())
  quit(status = 1)
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
    return(trim_quotes(x[[1]]))
  }

  if (is.symbol(x)) {
    return(as.character(x))
  }

  trim_quotes(paste(deparse(x), collapse = ""))
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

extract_call_match <- function(expr) {
  call_name <- normalize_call_name(expr)

  if (is.null(call_name) || !call_name %in% c("setMethod", "setGeneric", "setClass")) {
    return(list(kind = NULL, value = NULL))
  }

  value <- if (length(expr) >= 2) normalize_selector_value(expr[[2]]) else NULL
  list(kind = call_name, value = value)
}

parse_args <- function(args) {
  opts <- list(
    manifest = "tools/refactor/manifest-wave1.yml",
    repo_root = ".",
    output_root = NULL,
    dry_run = FALSE,
    write_targets = FALSE,
    preserve_existing_targets = FALSE,
    rewrite_sources = FALSE,
    emit_collate = NULL,
    verbose = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--manifest") {
      i <- i + 1
      opts$manifest <- args[[i]]
    } else if (arg == "--repo-root") {
      i <- i + 1
      opts$repo_root <- args[[i]]
    } else if (arg == "--output-root") {
      i <- i + 1
      opts$output_root <- args[[i]]
    } else if (arg == "--dry-run") {
      opts$dry_run <- TRUE
    } else if (arg == "--write-targets") {
      opts$write_targets <- TRUE
    } else if (arg == "--preserve-existing-targets") {
      opts$preserve_existing_targets <- TRUE
    } else if (arg == "--rewrite-sources") {
      opts$rewrite_sources <- TRUE
    } else if (arg == "--emit-collate") {
      i <- i + 1
      opts$emit_collate <- args[[i]]
    } else if (arg == "--verbose") {
      opts$verbose <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1
  }

  if (!opts$dry_run && !opts$write_targets && !opts$rewrite_sources && is.null(opts$emit_collate)) {
    opts$dry_run <- TRUE
  }

  if (opts$rewrite_sources && !opts$write_targets) {
    abort("--rewrite-sources requires --write-targets")
  }

  if (is.null(opts$output_root)) {
    opts$output_root <- opts$repo_root
  }

  opts
}

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

read_source_lines <- function(path) {
  readLines(path, warn = FALSE)
}

build_inventory <- function(path) {
  exprs <- parse(file = path, keep.source = TRUE)
  srcrefs <- attr(exprs, "srcref")

  lapply(seq_along(exprs), function(i) {
    expr <- exprs[[i]]
    call_match <- extract_call_match(expr)

    list(
      index = i,
      expr = expr,
      srcref = srcrefs[[i]],
      symbol = extract_top_level_symbol(expr),
      call_kind = call_match$kind,
      call_value = call_match$value
    )
  })
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

normalize_block_text <- function(lines, start_line, end_line) {
  block <- lines[start_line:end_line]

  while (length(block) > 0 && block[[length(block)]] == "") {
    block <- block[-length(block)]
  }

  paste0(paste(block, collapse = "\n"), "\n")
}

resolve_expr_selector <- function(entry, selector, lines, inventory) {
  matches <- switch(
    selector$kind,
    symbol = Filter(function(rec) identical(rec$symbol, selector$value), inventory),
    setMethod = Filter(function(rec) identical(rec$call_kind, "setMethod") && identical(rec$call_value, selector$value), inventory),
    setGeneric = Filter(function(rec) identical(rec$call_kind, "setGeneric") && identical(rec$call_value, selector$value), inventory),
    setClass = Filter(function(rec) identical(rec$call_kind, "setClass") && identical(rec$call_value, selector$value), inventory),
    expr_index = {
      idx <- as.integer(selector$value)
      if (is.na(idx) || idx < 1L || idx > length(inventory)) {
        abort("Invalid expr_index for entry ", entry$id, ": ", selector$value)
      }
      inventory[idx]
    },
    NULL
  )

  if (is.null(matches)) {
    abort("Unsupported selector kind for entry ", entry$id, ": ", selector$kind)
  }

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

group_blocks_by_target <- function(blocks) {
  out <- list()
  for (block in blocks) {
    target <- block$target
    out[[target]] <- c(out[[target]], list(block))
  }
  out
}

write_target_files <- function(grouped_blocks, output_root, preserve_existing_targets = FALSE, verbose = FALSE) {
  for (target in names(grouped_blocks)) {
    if (!nzchar(target)) {
      next
    }

    target_path <- file.path(output_root, target)
    dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)

    contents <- vapply(grouped_blocks[[target]], function(block) block$text, character(1))
    new_text <- paste(contents, collapse = "\n")

    if (isTRUE(preserve_existing_targets) && file.exists(target_path)) {
      existing_lines <- readLines(target_path, warn = FALSE)
      existing_text <- paste(existing_lines, collapse = "\n")
      pieces <- character()
      if (nzchar(existing_text)) {
        pieces <- c(pieces, existing_text)
      }
      if (nzchar(new_text)) {
        pieces <- c(pieces, new_text)
      }
      writeLines(paste(pieces, collapse = "\n"), con = target_path)
    } else {
      writeLines(new_text, con = target_path)
    }

    if (isTRUE(verbose)) {
      message("Wrote target: ", target)
    }
  }
}

merge_ranges <- function(ranges) {
  if (!length(ranges)) {
    return(ranges)
  }

  ord <- order(vapply(ranges, `[[`, integer(1), "start"), vapply(ranges, `[[`, integer(1), "end"))
  ranges <- ranges[ord]
  merged <- list(ranges[[1]])

  for (rng in ranges[-1]) {
    last <- merged[[length(merged)]]
    if (rng$start <= (last$end + 1L)) {
      last$end <- max(last$end, rng$end)
      merged[[length(merged)]] <- last
    } else {
      merged[[length(merged) + 1L]] <- rng
    }
  }

  merged
}

rewrite_source_files <- function(blocks, repo_root, verbose = FALSE) {
  blocks <- Filter(function(block) identical(block$action, "extract"), blocks)
  by_source <- split(blocks, vapply(blocks, `[[`, character(1), "source"))

  for (source in names(by_source)) {
    source_path <- file.path(repo_root, source)
    lines <- read_source_lines(source_path)
    ranges <- lapply(by_source[[source]], function(block) list(start = block$start_line, end = block$end_line))
    ranges <- merge_ranges(ranges)

    keep <- rep(TRUE, length(lines))
    for (rng in ranges) {
      keep[rng$start:rng$end] <- FALSE
    }

    rewritten <- lines[keep]
    writeLines(rewritten, con = source_path)

    if (isTRUE(verbose)) {
      message("Rewrote source: ", source)
    }
  }
}

resolve_emit_collate_path <- function(repo_root, output_root, path) {
  if (grepl("^/", path)) {
    return(path)
  }

  if (dirname(path) == ".") {
    return(file.path(output_root, path))
  }

  file.path(repo_root, path)
}

emit_collate <- function(manifest, repo_root, output_root, path) {
  targets <- manifest$collate_targets %||% unique(vapply(Filter(function(x) {
    !identical(x$action %||% "extract", "skip") && !is.null(x$target)
  }, manifest$entries), function(entry) entry$target, character(1)))

  target_path <- resolve_emit_collate_path(repo_root, output_root, path)
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(targets, con = target_path)
}

print_summary <- function(resolved, manual_merge) {
  cat("Resolved entries:\n")
  for (block in resolved) {
    cat(
      sprintf(
        "- %s [%s] %s:%d-%d -> %s\n",
        block$id,
        block$selector_kind,
        block$source,
        block$start_line,
        block$end_line,
        block$target
      )
    )
  }

  if (length(manual_merge)) {
    cat("\nManual merge entries:\n")
    for (block in manual_merge) {
      cat(
        sprintf(
          "- %s %s:%d-%d (%s)\n",
          block$id,
          block$source,
          block$start_line,
          block$end_line,
          block$selector_value
        )
      )
    }
  }
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  output_root <- normalizePath(opts$output_root, mustWork = FALSE)
  manifest_path <- file.path(repo_root, opts$manifest)

  if (!file.exists(manifest_path)) {
    abort("Manifest does not exist: ", opts$manifest)
  }

  manifest <- read_manifest(manifest_path)
  cache <- list()
  resolved <- list()

  for (entry in manifest$entries) {
    action <- entry$action %||% "extract"
    resolved_entry <- resolve_entry(entry, repo_root, cache)
    cache <- resolved_entry$cache
    resolved[[length(resolved) + 1L]] <- resolved_entry$block
  }

  manual_merge <- Filter(function(block) identical(block$action, "manual_merge"), resolved)
  extractable <- Filter(function(block) identical(block$action, "extract"), resolved)

  print_summary(resolved, manual_merge)

  if (isTRUE(opts$dry_run)) {
    cat("\nDry run only. No files written.\n")
  }

  if (isTRUE(opts$write_targets)) {
    grouped <- group_blocks_by_target(extractable)
    write_target_files(
      grouped,
      output_root,
      preserve_existing_targets = opts$preserve_existing_targets,
      verbose = opts$verbose
    )
  }

  if (isTRUE(opts$rewrite_sources)) {
    rewrite_source_files(extractable, repo_root, verbose = opts$verbose)
  }

  if (!is.null(opts$emit_collate)) {
    emit_collate(manifest, repo_root, output_root, opts$emit_collate)
    if (isTRUE(opts$verbose)) {
      message("Wrote collate list: ", opts$emit_collate)
    }
  }

  if (length(manual_merge)) {
    cat("\nUnresolved manual merges: ", length(manual_merge), "\n", sep = "")
  }
}

main()
