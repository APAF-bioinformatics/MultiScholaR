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

parse_all_r_files <- function(repo_root) {
  files <- list_r_files(repo_root)
  errs <- character()

  for (relative_path in files) {
    full_path <- file.path(repo_root, relative_path)
    ok <- tryCatch({
      parse(file = full_path, keep.source = TRUE)
      TRUE
    }, error = function(e) {
      errs <<- c(errs, sprintf("%s: %s", relative_path, conditionMessage(e)))
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

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), con = stderr())
    quit(status = 1)
  }
)
