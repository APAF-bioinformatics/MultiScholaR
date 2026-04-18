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

extract_top_level_symbol <- function(expr) {
  call_name <- normalize_call_name(expr)

  if (is.call(expr) && length(expr) >= 3 && !is.null(call_name) && call_name %in% c("<-", "=")) {
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
    } else if (arg == "--verbose") {
      opts$verbose <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1
  }

  opts
}

read_manifest <- function(path) {
  yaml::read_yaml(path)
}

build_inventory <- function(path) {
  exprs <- parse(file = path, keep.source = TRUE)

  lapply(seq_along(exprs), function(i) {
    expr <- exprs[[i]]
    call_match <- extract_call_match(expr)

    list(
      index = i,
      symbol = extract_top_level_symbol(expr),
      call_kind = call_match$kind,
      call_value = call_match$value
    )
  })
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

  matches <- switch(
    kind,
    symbol = Filter(function(rec) identical(rec$symbol, selector$value), inventory),
    setMethod = Filter(function(rec) identical(rec$call_kind, "setMethod") && identical(rec$call_value, selector$value), inventory),
    setGeneric = Filter(function(rec) identical(rec$call_kind, "setGeneric") && identical(rec$call_value, selector$value), inventory),
    setClass = Filter(function(rec) identical(rec$call_kind, "setClass") && identical(rec$call_value, selector$value), inventory),
    expr_index = {
      idx <- as.integer(selector$value)
      if (is.na(idx) || idx < 1L || idx > length(inventory)) {
        return("expr_index is out of bounds")
      }
      inventory[idx]
    },
    NULL
  )

  if (is.null(matches)) {
    return(paste("unsupported selector kind", kind))
  }

  if (!length(matches)) {
    return("selector matched 0 blocks")
  }

  if (length(matches) != 1L) {
    return(sprintf("selector matched %d blocks", length(matches)))
  }

  NULL
}

parse_all_r_files <- function(repo_root) {
  files <- list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE)
  errs <- character()

  for (file in files) {
    ok <- tryCatch({
      parse(file = file, keep.source = TRUE)
      TRUE
    }, error = function(e) {
      errs <<- c(errs, sprintf("%s: %s", basename(file), conditionMessage(e)))
      FALSE
    })

    if (!ok) {
      next
    }
  }

  errs
}

check_target_duplicates <- function(repo_root, manifest) {
  targets <- unique(manifest$collate_targets %||% vapply(Filter(function(entry) {
    !identical(entry$action %||% "extract", "skip") && !is.null(entry$target)
  }, manifest$entries), function(entry) entry$target, character(1)))

  existing_targets <- Filter(function(path) file.exists(file.path(repo_root, path)), targets)
  if (!length(existing_targets)) {
    return(character())
  }

  records <- list()
  for (target in existing_targets) {
    inventory <- build_inventory(file.path(repo_root, target))
    for (rec in inventory) {
      if (!is.null(rec$symbol)) {
        records[[length(records) + 1L]] <- data.frame(
          symbol = rec$symbol,
          file = target,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(records)) {
    return(character())
  }

  df <- do.call(rbind, records)
  dup_symbols <- unique(df$symbol[duplicated(df$symbol)])

  if (!length(dup_symbols)) {
    return(character())
  }

  vapply(dup_symbols, function(sym) {
    files <- paste(df$file[df$symbol == sym], collapse = ", ")
    sprintf("%s defined in %s", sym, files)
  }, character(1))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  manifest_path <- file.path(repo_root, opts$manifest)

  if (!file.exists(manifest_path)) {
    abort("Manifest does not exist: ", opts$manifest)
  }

  manifest <- read_manifest(manifest_path)
  failures <- character()

  parse_failures <- parse_all_r_files(repo_root)
  failures <- c(failures, parse_failures)

  by_source <- split(manifest$entries, vapply(manifest$entries, function(entry) entry$source, character(1)))
  for (source in names(by_source)) {
    source_path <- file.path(repo_root, source)
    if (!file.exists(source_path)) {
      failures <- c(failures, sprintf("%s: source file does not exist", source))
      next
    }

    inventory <- build_inventory(source_path)
    lines <- readLines(source_path, warn = FALSE)

    for (entry in by_source[[source]]) {
      err <- validate_selector(entry, inventory, lines)
      if (!is.null(err)) {
        failures <- c(failures, sprintf("%s: %s", entry$id, err))
      }
    }
  }

  duplicate_failures <- check_target_duplicates(repo_root, manifest)
  failures <- c(failures, duplicate_failures)

  if (length(failures)) {
    cat("Verification failures:\n")
    for (failure in failures) {
      cat("- ", failure, "\n", sep = "")
    }
    quit(status = 1)
  }

  cat("Verification passed.\n")
}

main()
