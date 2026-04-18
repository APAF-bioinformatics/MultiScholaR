#!/usr/bin/env Rscript

resolve_refactor_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    ofile <- frames[[idx]]$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }

  normalizePath("tools/refactor", mustWork = TRUE)
}

refactor_dir <- resolve_refactor_dir()
source(file.path(refactor_dir, "fidelity_inventory.R"))
source(file.path(refactor_dir, "fidelity_manifest.R"))

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
    cat("\nNOTE: manual_merge entries were not written. Resolve them manually.\n")
  }
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), con = stderr())
    quit(status = 1)
  }
)
