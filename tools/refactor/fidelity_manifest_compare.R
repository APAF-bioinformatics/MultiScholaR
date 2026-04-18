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
    manifest = NULL,
    baseline_root = NULL,
    target_root = NULL,
    pretty = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--manifest")) {
      i <- i + 1L
      opts$manifest <- args[[i]]
    } else if (identical(arg, "--baseline-root")) {
      i <- i + 1L
      opts$baseline_root <- args[[i]]
    } else if (identical(arg, "--target-root")) {
      i <- i + 1L
      opts$target_root <- args[[i]]
    } else if (identical(arg, "--pretty")) {
      opts$pretty <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1L
  }

  if (is.null(opts$manifest) || is.null(opts$baseline_root) || is.null(opts$target_root)) {
    abort("--manifest, --baseline-root, and --target-root are required")
  }

  opts
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  manifest_path <- normalizePath(opts$manifest, mustWork = TRUE)
  baseline_root <- normalizePath(opts$baseline_root, mustWork = TRUE)
  target_root <- normalizePath(opts$target_root, mustWork = TRUE)
  manifest <- read_manifest(manifest_path)

  baseline_cache <- list()
  target_cache <- list()
  comparisons <- list()

  for (entry in manifest$entries) {
    result <- compare_manifest_entry(entry, baseline_root, target_root, baseline_cache, target_cache)
    baseline_cache <- result$baseline_cache
    target_cache <- result$target_cache
    comparisons[[length(comparisons) + 1L]] <- result$comparison
  }

  output <- list(
    manifest_path = manifest_path,
    entry_count = length(comparisons),
    comparisons = comparisons
  )

  cat(
    jsonlite::toJSON(
      output,
      auto_unbox = TRUE,
      pretty = isTRUE(opts$pretty),
      null = "null",
      na = "null"
    )
  )
  cat("\n")
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), con = stderr())
    quit(status = 1)
  }
)
